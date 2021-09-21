# NOTE
# This script is kept here only for archiving purpose.
# It was used for a one-time analysis to assess variant QC, but is not used as a regular part of gnomAD production

from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.evaluation import add_rank
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.variant_qc import *
import argparse
import logging

logger = logging.getLogger("select_qc_set")

def get_trio_samples_to_keep(trios: hl.Table, n_keep: int) -> hl.Table:
    trios = trios.annotate(r_prob=hl.rand_unif(0, 1))
    trios = trios.order_by('r_prob')
    trios = trios.add_index()
    keep_trios = trios.filter(trios.idx < n_keep)
    return keep_trios

def unpack_trios(trios: hl.Table) -> hl.expr.SetExpression:
    probands = trios.aggregate(hl.agg.collect_as_set(trios.proband_id))
    fathers = trios.aggregate(hl.agg.collect_as_set(trios.father_id))
    mothers = trios.aggregate(hl.agg.collect_as_set(trios.mother_id))
    relateds = hl.set(probands.union(fathers).union(mothers))
    return relateds


def main(args):
    hl.init(log='/select_qc_set.log')
    data_type = 'exomes' if args.exomes else 'genomes'
    N = args.sample_size

    if args.subsample:
        # Load callset data
        mt = get_gnomad_data(data_type)  # NOTE: This has been split!

        # Load family data
        exome_trios = hl.import_table(fam_path('exomes'), no_header=True)
        exome_trios = exome_trios.rename(
            {'f0': 'fam_id', 'f1': 'proband_id', 'f2': 'father_id', 'f3': 'mother_id', 'f4': 'sex', 'f5': 'phenotype'})

        genome_trios = hl.import_table(fam_path('genomes'), no_header=True)
        genome_trios = genome_trios.rename(
            {'f0': 'fam_id', 'f1': 'proband_id', 'f2': 'father_id', 'f3': 'mother_id', 'f4': 'sex', 'f5': 'phenotype'})

        trios = exome_trios if args.exomes else genome_trios
        n_exome_trios = exome_trios.count()
        n_genome_trios = genome_trios.count()

        if n_exome_trios > n_genome_trios:
            desired_n_trios = n_genome_trios
            if args.exomes:
                keep_trios = get_trio_samples_to_keep(exome_trios, desired_n_trios)
            else:
                keep_trios = genome_trios
        else:
            desired_n_trios = n_exome_trios
            if args.exomes:
                keep_trios = exome_trios
            else:
                keep_trios = get_trio_samples_to_keep(genome_trios, desired_n_trios)

        keep_count = keep_trios.count()
        total_count = trios.count()
        logger.info('Keeping {0} trios out of {1} from {2}'.format(keep_count, total_count, data_type))

        # Annotate keep status for relateds
        relateds = unpack_trios(trios)
        keep_relateds = unpack_trios(keep_trios)
        mt = mt.annotate_cols(is_related=relateds.contains(mt.s), keep_related=mt.meta.high_quality & keep_relateds.contains(mt.s))
        n_high_quality_related = mt.aggregate_cols(hl.agg.count_where(mt.keep_related))
        logger.info('Keeping {0} high-quality individuals out of {1} individuals in trios'.format(n_high_quality_related, hl.eval(hl.len(keep_relateds))))

        # Randomly sample high-quality unrelateds up to total number desired
        n_unrelateds_to_keep = N - n_high_quality_related
        mt = mt.annotate_cols(r_prob=hl.cond(mt.is_related, 1, hl.rand_unif(0, 1)))
        sample_ht = mt.cols()
        sample_ht = sample_ht.filter(sample_ht.meta.high_quality).order_by('r_prob')
        keep_unrelateds = sample_ht.aggregate(hl.agg.collect_as_set(sample_ht.s))
        keep_unrelateds = keep_unrelateds[:n_unrelateds_to_keep]
        mt = mt.annotate_cols(keep_unrelated=keep_unrelateds.contains(mt.s))

        # Filter to keep samples
        mt = mt.filter_cols(mt.keep_unrelated | mt.keep_related)
        mt = mt.annotate_rows(n_nonref=hl.agg.count_where(mt.GT.is_non_ref()))
        mt = mt.filter_rows(mt.n_nonref > 0)
        mt.write('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.mt'.format(data_type, N), overwrite=True)

        mt = hl.read_matrix_table('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.mt'.format(data_type, N))
        v,s = mt.count()
        logger.info('{0} samples and {1} variants found in subsampled callset'.format(s, v))
        ht = mt.rows().repartition(1000, shuffle=False)
        ht.write('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.ht'.format(data_type, N), overwrite=True)

    if args.add_vqsr_rank:
        logger.info('Ranking indels and SNPs...')
        ht = hl.read_table('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.ht'.format(data_type, N))
        ht = ht.key_by('locus', 'alleles')
        ht = add_rank(ht, ht.info.VQSLOD)
        ht.write('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.vqsr.ranked.ht'.format(data_type, N), overwrite=True)

        logger.info('Adding rankings for bi-allelics and singletons...')
        ht = hl.read_table('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.vqsr.ranked.ht'.format(data_type, N))
        # ht = add_full_rankings(ht, ht.info.VQSLOD) # FIXME:  This is stale
        ht.write('gs://gnomad/variant_qc/temp/gnomad.{0}.{1}_samples.vqsr.full_rankings.ht'.format(data_type, N), overwrite=True)  # TODO: move file locations once script is finalized

    # TODO: add code to rank RF?

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--sample_size', help='Total number of samples desired in subsample', default=10000, type=int)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--subsample', help='Subsample original callset.', action='store_true')
    parser.add_argument('--add_vqsr_rank', help='Add variant rankings.', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
