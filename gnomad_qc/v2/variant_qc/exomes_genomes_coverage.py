# NOTE
# This script is kept here only for archiving purpose.
# It was used in the gnomAD v2 flagship LOF paper (https://www.nature.com/articles/s41586-020-2308-7), but is not used as a regular part of gnomAD production

from gnomad_qc.v2.resources import *
from gnomad_qc.v2.resources.variant_qc import get_ucsc_mappability
import argparse
import logging

logger = logging.getLogger("exomes_genomes_coverage")

COVERAGE_BINS = [1] + list(range(5, 31, 5)) + [50, 100]


def import_cds_from_gtf(overwrite: bool) -> hl.Table:
    """
    Creates a HT with a row for each base / gene that is in a CDS in gencode v19
    :param bool overwrite: Overwrite existing table
    :return: HT
    :rtype: Table
    """
    gtf = hl.experimental.import_gtf(
        'gs://hail-common/references/gencode/gencode.v19.annotation.gtf.bgz',
        reference_genome='GRCh37',
        skip_invalid_contigs=True, min_partitions=200
    )
    gtf = gtf.filter((gtf.feature == 'CDS') & (gtf.transcript_type == 'protein_coding') & (gtf.tag == 'basic'))
    gtf = gtf.annotate(locus=hl.range(gtf.interval.start.position, gtf.interval.end.position).map(lambda x: hl.locus(gtf.interval.start.contig, x, 'GRCh37')))
    gtf = gtf.key_by().select('gene_id', 'locus', 'gene_name').explode('locus')
    gtf = gtf.key_by('locus', 'gene_id').distinct()
    return gtf.checkpoint('gs://gnomad-tmp/gencode_grch37.gene_by_base.ht', overwrite=overwrite, _read_if_exists=not overwrite)


def compute_per_base_cds_coverage(overwrite: bool):
    """
    Creates a HT with a row for each base and gene that is in a CDS in gencode v19 with the following information:
    1) gene id
    2) gnomAD exomes coverage by platform
    3) gnomAD genomes coverage by PCR status

    :param bool overwrite: Whether to overwrite existing results
    """

    def get_haploid_coverage_struct(mt):
        return hl.cond(
            mt.locus.in_autosome_or_par(),
            hl.struct(
                mean=0.5 * hl.agg.sum(mt.mean * mt.n) / hl.agg.sum(mt.n),
                median=0.5 * hl.median(hl.agg.collect(mt.median)),
                **{
                    f'over_{i / 2}': hl.agg.sum(mt[f'over_{float(i)}']) / hl.agg.sum(mt.n)
                    for i in COVERAGE_BINS
                }
            ),
            hl.struct(
                mean=hl.agg.sum(hl.cond(mt.sex == 'male', mt.mean * mt.n, 0.5 * mt.mean * mt.n)) / hl.agg.sum(mt.n),
                median=hl.median(hl.agg.collect(hl.cond(mt.sex == 'male', mt.median, 0.5 * mt.median))),
                **{
                    f'over_{i / 2}':
                        hl.null(hl.tfloat32) if not f'over_{i / 2}' in mt.entry else
                        hl.agg.sum(
                            hl.cond(
                                mt.sex == 'male',
                                mt[f'over_{i / 2}'],
                                mt[f'over_{float(i)}']
                            )
                        ) / hl.agg.sum(mt.n)
                    for i in COVERAGE_BINS
                }
            )
        )

    def get_cov_ht(data_type: str) -> hl.Table:
        cov_mt = hl.read_matrix_table(coverage_mt_path(data_type, grouped=True))
        qc_platforms = cov_mt.aggregate_cols(hl.agg.collect_as_set(cov_mt.qc_platform))
        return cov_mt.select_rows(
            **{
                f'{data_type}_{platform}': hl.agg.filter(
                    cov_mt.qc_platform == platform,
                    get_haploid_coverage_struct(cov_mt)
                ) for platform in qc_platforms
            },
            **{f'{data_type}_all': get_haploid_coverage_struct(cov_mt)}
        ).rows()

    genomes_cov = get_cov_ht('genomes')
    exomes_cov = get_cov_ht('exomes')
    ucsc_mappability = get_ucsc_mappability()
    gtf = import_cds_from_gtf(overwrite)

    gene_by_base_cov_ht = gtf.annotate(
        **exomes_cov[gtf.locus],
        **genomes_cov[gtf.locus],
        mappability=ucsc_mappability[gtf.locus].duke_35_map
    )
    gene_by_base_cov_ht.write('gs://gnomad-public/papers/2019-flagship-lof/v1.1/summary_gene_coverage/gencode_grch37.gene_by_base.cov.ht', overwrite=overwrite)


def export_gene_coverage(overwrite: bool):
    """
    Exports gene coverage summary stats as HT and tsv.
    """
    min_good_coverage_dp = 'over_10.0'
    min_good_coverage_prop = 0.8

    cov = hl.read_table('gs://gnomad-lfran/exomes_genomes_coverage/gencode_grch37.gene_by_base.cov.ht')
    cov = hl.filter_intervals(cov, [hl.parse_locus_interval('Y')], keep=False)

    cov = cov.group_by(cov.gene_id).aggregate(
        gene_name=hl.agg.take(cov.gene_name, 1)[0],
        gene_interval=hl.interval(
            hl.locus(hl.agg.take(cov.locus.contig, 1)[0], hl.agg.min(cov.locus.position)),
            hl.locus(hl.agg.take(cov.locus.contig, 1)[0], hl.agg.max(cov.locus.position))
        ),
        mean_mappability=hl.agg.mean(cov.mappability),
        median_mappability=hl.agg.approx_quantiles(cov.mappability, 0.5),
        cov_stats=[
            hl.struct(
                data_type=x.split("_")[0],
                platform=hl.delimit(x.split("_")[1:], delimiter="_"),
                frac_well_covered_bases=hl.agg.fraction(cov[x][min_good_coverage_dp] >= min_good_coverage_prop),
                mean_well_covered_samples=hl.agg.mean(cov[x][min_good_coverage_dp])
            ) for x in list(cov.row_value) if x.startswith('exomes') or x.startswith('genomes')
        ]
    )

    cov = cov.explode('cov_stats')
    cov = cov.select(
        'gene_name',
        'gene_interval',
        'mean_mappability',
        'median_mappability',
        **cov.cov_stats

    )
    cov = cov.annotate_globals(
        min_good_coverage_prop=min_good_coverage_prop,
        min_good_coverage_dp=min_good_coverage_dp
    )

    cov = cov.checkpoint('gs://gnomad-public/papers/2019-flagship-lof/v1.1/summary_gene_coverage/gencode_grch37_gene_by_platform_coverage_summary.ht', overwrite=overwrite, _read_if_exists=not overwrite)
    if not overwrite and hl.hadoop_is_file('gs://gnomad-public/papers/2019-flagship-lof/v1.1/summary_gene_coverage/gencode_grch37_gene_by_platform_coverage_summary.tsv.gz'):
        logger.warn("gs://gnomad-public/papers/2019-flagship-lof/v1.1/summary_gene_coverage/gencode_grch37_gene_by_platform_coverage_summary.tsv.gz not exported as it already exists and --overwrite is not set.")
    else:
        cov.export('gs://gnomad-public/papers/2019-flagship-lof/v1.1/summary_gene_coverage/gencode_grch37_gene_by_platform_coverage_summary.tsv.gz')


def main(args):
    if args.compute_per_base_cds_coverage:
        print("Computing per-base CDS coverage")
        compute_per_base_cds_coverage(args.overwrite)

    if args.export_gene_coverage:
        print("Exporting gene coverage")
        export_gene_coverage(args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compute_per_base_cds_coverage', help='Computes per-base coverage for CDS regions.', action='store_true')
    parser.add_argument('--export_gene_coverage', help='Exports gene coverage in long format to tsv.', action='store_true')
    parser.add_argument('--overwrite', help='Overwrites existing results if set.', action='store_true')
    args = parser.parse_args()
    main(args)
