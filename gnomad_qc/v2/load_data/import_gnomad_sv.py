from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources import *
import argparse
import logging

logger = logging.getLogger('import_gnomad_sv')


def import_vcf():
    mt = hl.import_vcf(gnomad_sv_vcf_path, force_bgz=True, min_partitions=300)
    meta = hl.import_table(gnomad_sv_release_samples_list_path, force=True, key='s')
    mt = mt.annotate_cols(
        release=hl.is_defined(meta[mt.col_key])
    )
    return mt


def merge_genomes_sv_samples(genomes_mt: hl.MatrixTable, sv_mt: hl.MatrixTable) -> hl.Table:
    genomes_samples = genomes_mt.cols()

    genomes_meta = get_gnomad_meta('genomes')
    genomes_samples = genomes_samples.annotate(
        **genomes_meta[genomes_samples.s].select('project_id', 'release')
    )

    genomes_samples = genomes_samples.key_by(
        sv_s=(genomes_samples.project_id + "_" + genomes_samples.s).replace('\W+', '_')
    )
    genomes_samples = genomes_samples.select(
        genomes_s=genomes_samples.s,
        in_v2_genomes=True
    )

    gnomad_sv_samples = sv_mt.cols()

    gnomad_sv_samples = gnomad_sv_samples.annotate(
        in_v2_sv=True
    )

    merged_samples = genomes_samples.join(gnomad_sv_samples, how='outer')
    merged_samples = merged_samples.annotate(
        in_v2_genomes=hl.is_defined(merged_samples.in_v2_genomes),
        in_v2_sv=hl.is_defined(merged_samples.in_v2_sv)
    ).persist()

    print("Summary of sample sharing between gnomAD genomes and gnomad SVs:")
    merged_samples.group_by(
        'in_v2_genomes',
        'in_v2_sv'
    ).aggregate(
        n=hl.agg.count()
    ).show()

    merged_samples = merged_samples.filter(
        merged_samples.in_v2_genomes & merged_samples.in_v2_sv
    )

    return merged_samples


def order_cols(mt1: hl.MatrixTable, mt2: hl.MatrixTable) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    cols = mt1.add_col_index(name='idx1').cols().select('idx1')
    cols = cols.join(mt2.add_col_index(name='idx2').cols().select('idx2'))
    return (
        mt1.choose_cols(cols.idx1.collect()),
        mt2.choose_cols(cols.idx2.collect())
    )


def merge_with_short_variants(min_af: float):
    # Load gnomAD genomes and filter to PASS / release sites and samples
    gnomad_genomes = get_gnomad_data('genomes', adj=True, release_samples=True, release_annotations=True)
    gnomad_genomes = gnomad_genomes.filter_rows(
        (hl.len(gnomad_genomes.filters) == 0) &
        (gnomad_genomes.popmax[0].AF >= min_af)
    )

    # Load gnomAD SVs and filter to PASS / release sites and samples
    gnomad_sv = hl.read_matrix_table(gnomad_sv_mt_path)
    gnomad_sv = gnomad_sv.filter_cols(gnomad_sv.release)
    gnomad_sv = gnomad_sv.filter_rows(
        hl.len(gnomad_sv.filters) == 0
    )

    # Select common samples and unify samples IDs (using SV sample ID as key)
    merged_samples = merge_genomes_sv_samples(gnomad_genomes, gnomad_sv)

    merged_samples = merged_samples.key_by('sv_s')
    gnomad_sv = gnomad_sv.annotate_cols(
        genomes_s=merged_samples[gnomad_sv.s].genomes_s
    )
    gnomad_sv = gnomad_sv.filter_cols(hl.is_defined(gnomad_sv.genomes_s))

    merged_samples = merged_samples.key_by('genomes_s')
    gnomad_genomes = gnomad_genomes.filter_cols(
        hl.is_defined(merged_samples[gnomad_genomes.col_key])
    )
    gnomad_genomes = gnomad_genomes.annotate_cols(genomes_s=gnomad_genomes.s)
    gnomad_genomes = gnomad_genomes.key_cols_by(
        s=merged_samples[gnomad_genomes.s].sv_s
    )

    # Unify global, col, row and entry schemas
    gnomad_genomes = gnomad_genomes.select_entries('GT', GQ=hl.null(hl.tint32))
    gnomad_genomes = gnomad_genomes.select_rows(
        sv_info=hl.null(gnomad_sv.info.dtype)
    )
    gnomad_genomes = gnomad_genomes.select_cols('genomes_s')
    gnomad_genomes = gnomad_genomes.select_globals()

    gnomad_sv = gnomad_sv.select_entries('GT', 'GQ')
    gnomad_sv = gnomad_sv.select_rows(sv_info=gnomad_sv.info)
    gnomad_sv = gnomad_sv.select_cols('genomes_s')
    gnomad_sv = gnomad_sv.select_globals()

    gnomad_genomes, gnomad_sv = order_cols(gnomad_genomes, gnomad_sv)
    merged_mt = gnomad_sv.union_rows(gnomad_genomes)
    merged_mt = merged_mt.annotate_globals(
        min_af=min_af
    )
    return merged_mt


def generate_hists():
    mt = hl.read_matrix_table(gnomad_sv_mt_path)
    meta = get_gnomad_meta('genomes')
    meta = meta.key_by(s=(meta.project_id + "_" + meta.s).replace('\W+', '_'))
    mt = mt.annotate_cols(
        age=meta[mt.col_key].age,
        in_v2=hl.is_defined(meta[mt.col_key])
    )

    logger.info("Found {} samples with age data.".format(
        mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.age)))
    ))

    hists = mt.select_rows(
        age_hist_het=hl.or_missing(~mt.filters.contains('MULTIALLELIC'), hl.agg.filter(mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10))),
        age_hist_hom=hl.or_missing(~mt.filters.contains('MULTIALLELIC'), hl.agg.filter(mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10))),
        gq_hist_alt=hl.or_missing(~mt.filters.contains('MULTIALLELIC'), hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.GQ, 0, 100, 20))),
        gq_hist_all=hl.or_missing(~mt.filters.contains('MULTIALLELIC'), hl.agg.hist(mt.GQ, 0, 100, 20)),
    ).rows()

    hist_expr = {}
    for hist in ['age_hist_het', 'age_hist_hom', 'gq_hist_alt', 'gq_hist_all']:
        hist_expr.update({
            f'{hist}_bin_freq': hl.delimit(hists[hist].bin_freq, delimiter="|"),
            f'{hist}_bin_edges': hl.delimit(hists[hist].bin_edges, delimiter="|"),
            f'{hist}_n_smaller': hists[hist].n_smaller,
            f'{hist}_n_larger': hists[hist].n_larger
        })

    hists = hists.select(
        **hist_expr
    )

    return hists


def main(args):
    if args.import_vcf:
        mt = import_vcf()
        mt.write(gnomad_sv_mt_path, overwrite=args.overwrite)

    if args.merge_with_short_variants:
        mt = merge_with_short_variants(args.short_variants_af)
        mt.write('gs://gnomad/projects/genomes_sv_ld/genomes_union_sv.mt', overwrite=args.overwrite)

    if args.generate_hists:
        hists_ht = generate_hists()
        hists_ht.write('gs://gnomad-public/papers/2019-sv/gnomad_sv_hists.ht', overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--import_vcf', help='Imports gnomAD SV VCF and writes it as MT.', action='store_true')
    parser.add_argument('--merge_with_short_variants', help='Creates an MT merging SVs and short variants. An AF cutoff is used for short variants (--short_variants_af)', action='store_true')
    parser.add_argument('--short_variants_af', help='Short variant AF cutoff (AF >= x in any pop) for merging SVs and short variants', default=0.005, type=float)
    parser.add_argument('--generate_hists', help='Generate age and GQ histograms', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
