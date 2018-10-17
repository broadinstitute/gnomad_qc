from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *
import datetime

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("finalize_sample_qc")
logger.setLevel(logging.INFO)


def add_release_annotations(ht: hl.Table) -> hl.Table:
    """

    :param Table ht: Table containing meta column annotations for the dataset
    :return: Table containing final 'high_quality' and 'release' sample status annotations
    :rtype: Table
    """
    ht = ht.annotate(high_quality=(hl.len(ht.hard_filters) == 0) & (hl.len(ht.pop_platform_filters) == 0))
    return ht.annotate(release=ht.high_quality & (hl.len(ht.perm_filters) == 0) & ~ht.related)


def main(args):
    subpop_nfe_ht = hl.read_table(subpop_ht_path('nfe')).select('subpop', 'known_subpop')
    subpop_eas_ht = hl.read_table(subpop_ht_path('eas')).select('subpop', 'known_subpop')
    subpop_ht = subpop_nfe_ht.union(subpop_eas_ht)

    exome_ht_a = hl.read_table(qc_ht_path('exomes', 'hard_filters'))
    exome_ht_b = hl.read_table(qc_ht_path('exomes', 'platforms'))
    exome_ht_b = exome_ht_b.rename({f'PC{i}': f'callrate_pc{i}' for i in range(1, 11)})
    exome_ht_c = hl.read_table(qc_ht_path('exomes', 'pop_platform'))
    exome_ht_c = exome_ht_c.transmute(**{f'PC{i+1}': exome_ht_c.scores[i] for i in range(10)}).drop('qc_platform')

    genome_ht_a = hl.read_table(qc_ht_path('genomes', 'hard_filters'))
    genome_ht_c = hl.read_table(qc_ht_path('genomes', 'pop_platform'))
    genome_ht_c = genome_ht_c.transmute(**{f'PC{i+1}': genome_ht_c.scores[i] for i in range(10)}).drop('qc_platform')

    super_meta_exome_ht = exome_ht_a.annotate(**exome_ht_b[exome_ht_a.key],
                                              **exome_ht_c[exome_ht_a.key],
                                              **subpop_ht[exome_ht_a.key])
    super_meta_exome_ht = add_release_annotations(super_meta_exome_ht).repartition(10)
    super_meta_exome_ht.write(metadata_exomes_ht_path(args.metadata_version_label), args.overwrite)
    super_meta_exome_ht.flatten().export(metadata_exomes_tsv_path(args.metadata_version_label))

    super_meta_genome_ht = genome_ht_a.annotate(**genome_ht_c[genome_ht_a.key],
                                                **subpop_ht[genome_ht_a.key])
    super_meta_genome_ht = add_release_annotations(super_meta_genome_ht).repartition(10)
    super_meta_genome_ht.write(metadata_genomes_ht_path(args.metadata_version_label), args.overwrite)
    super_meta_genome_ht.flatten().export(metadata_genomes_tsv_path(args.metadata_version_label))

    logger.info('{0} exome, {1} genome high-quality samples remaining after sample QC'.format(
        super_meta_exome_ht.aggregate(hl.agg.count_where(super_meta_exome_ht.high_quality)),
        super_meta_genome_ht.aggregate(hl.agg.count_where(super_meta_genome_ht.high_quality))))
    logger.info('{0} exome, {1} genome release samples remaining after sample QC'.format(
        super_meta_exome_ht.aggregate(hl.agg.count_where(super_meta_exome_ht.release)),
        super_meta_genome_ht.aggregate(hl.agg.count_where(super_meta_genome_ht.release))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--metadata_version_label', help='Label for metadata version (usually YYYY-MM-DD)', default=datetime.datetime.today().strftime('%Y-%m-%d'))
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
