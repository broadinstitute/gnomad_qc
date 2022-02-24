import argparse
import logging
from typing import Tuple

import hail as hl

from gnomad.sample_qc.pipeline import annotate_sex
from gnomad.sample_qc.sex import get_ploidy_cutoffs, get_sex_expr

from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.sample_qc import hard_filtered_ac_an_af, sex

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


def compute_sex(
    vds,
    use_gnomad_methods: bool = False,
    aaf_threshold: float = 0.001,
    f_stat_cutoff: float = 0.5,
    high_cov_intervals: bool = False,
    per_platform_high_cov_intervals: bool = False,
    normalization_contig: str = 'chr20',
    x_cov: int = None,
    y_cov: int = None,
    norm_cov: int = None,
    prop_samples_x: float = None,
    prop_samples_y: float = None,
    prop_samples_norm: float = None,
    reference_only: bool = False,
    variants_only: bool = False,
) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    :param use_gnomad_methods:
    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :return: Table with inferred sex annotation
    :rtype: hl.Table
    """

    # TODO: Determine what variants to use for f-stat calculation, current implementation will use a HT created by performing a full densify to get callrate and AF, note, this might not be needed, see CCDG f-stat values
    freq_ht = hard_filtered_ac_an_af.ht()

    if high_cov_intervals:
        coverage_mt = interval_coverage(interval_name, padding).mt()
        coverage_mt = coverage_mt.filter_rows(hl.literal({"chrY", "chrX", "chr20"}).contains(coverage_mt.interval.start.contig))
        coverage_mt = coverage_mt.annotate_rows(
            **{
                f"over_{y_cov}x": hl.agg.fraction(coverage_mt.mean_dp > y_cov),
                f"over_{x_cov}x": hl.agg.fraction(coverage_mt.mean_dp > x_cov),
                f"over_{norm_cov}x": hl.agg.fraction(coverage_mt.mean_dp > norm_cov),
            }
        )
        coverage_mt = coverage_mt.filter_rows(
            ((coverage_mt.interval.start.contig == "chrX") & (coverage_mt[f"over_{x_cov}x"] > prop_samples_x))
            | ((coverage_mt.interval.start.contig == "chrY") & (coverage_mt[f"over_{y_cov}x"] > prop_samples_y))
            | ((coverage_mt.interval.start.contig == normalization_contig) & (coverage_mt[f"over_{norm_cov}x"] > prop_samples_norm))
        )
        calling_intervals_ht = coverage_mt.rows()
    else:
        calling_intervals_ht = calling_intervals(interval_name, padding).ht()

    sex_ht = annotate_sex(
        hl.vds.to_merged_sparse_mt(vds) if use_gnomad_methods else vds,
        included_intervals=calling_intervals_ht,
        normalization_contig=normalization_contig,
        sites_ht=freq_ht,
        aaf_expr="AF",
        gt_expr="LGT",
        f_stat_cutoff=f_stat_cutoff,
        aaf_threshold=aaf_threshold,
        reference_only=reference_only,
        variants_only=variants_only,
    )

    return sex_ht


def reannotate_sex(
    cov_threshold: int,
    x_ploidy_cutoffs: Tuple[float, Tuple[float, float], float],
    y_ploidy_cutoffs: Tuple[Tuple[float, float], float],
):
    """
    Rerun sex karyotyping annotations without re-computing sex imputation metrics.

    :param cov_threshold: Filtering threshold to use for chr20 coverage in `compute_hard_filters`
    :param x_ploidy_cutoffs: Tuple of X chromosome ploidy cutoffs: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for triple X)
    :param y_ploidy_cutoffs: Tuple of Y chromosome ploidy cutoffs: ((lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)
    :return: Table with sex karyotyping annotations
    :rtype: hl.Table
    """
    # Copy HT to temp location to overwrite annotation
    sex_ht = sex.ht().checkpoint("gs://gnomad-tmp/sex_ht_checkpoint.ht", overwrite=True)
    hard_filter_ht = compute_hard_filters(cov_threshold, include_sex_filter=False)

    # Copy HT to temp location because it uses sex_ht for chr20 coverage
    hard_filter_ht = hard_filter_ht.checkpoint(
        "gs://gnomad-tmp/hardfilter_checkpoint.ht", overwrite=True
    )
    new_x_ploidy_cutoffs, new_y_ploidy_cutoffs = get_ploidy_cutoffs(
        sex_ht.filter(hl.is_missing(hard_filter_ht[sex_ht.key])), f_stat_cutoff=0.5
    )
    x_ploidy_cutoffs = hl.struct(
        upper_x=x_ploidy_cutoffs[0] if x_ploidy_cutoffs[0] else new_x_ploidy_cutoffs[0],
        lower_xx=x_ploidy_cutoffs[1][0]
        if x_ploidy_cutoffs[1][0]
        else new_x_ploidy_cutoffs[1][0],
        upper_xx=x_ploidy_cutoffs[1][1]
        if x_ploidy_cutoffs[1][1]
        else new_x_ploidy_cutoffs[1][1],
        lower_xxx=x_ploidy_cutoffs[2]
        if x_ploidy_cutoffs[2]
        else new_x_ploidy_cutoffs[2],
    )
    y_ploidy_cutoffs = hl.struct(
        lower_y=y_ploidy_cutoffs[0][0]
        if y_ploidy_cutoffs[0][0]
        else new_y_ploidy_cutoffs[0][0],
        upper_y=y_ploidy_cutoffs[0][1]
        if y_ploidy_cutoffs[0][1]
        else new_y_ploidy_cutoffs[0][1],
        lower_yy=y_ploidy_cutoffs[1]
        if y_ploidy_cutoffs[1]
        else new_y_ploidy_cutoffs[1],
    )
    sex_ht = sex_ht.annotate(
        **get_sex_expr(
            sex_ht.chrX_ploidy,
            sex_ht.chrY_ploidy,
            (
                x_ploidy_cutoffs["upper_x"],
                (x_ploidy_cutoffs["lower_xx"], x_ploidy_cutoffs["upper_xx"]),
                x_ploidy_cutoffs["lower_xxx"],
            ),
            (
                (y_ploidy_cutoffs["lower_y"], y_ploidy_cutoffs["upper_y"]),
                y_ploidy_cutoffs["lower_yy"],
            ),
        )
    )
    sex_ht = sex_ht.annotate_globals(
        x_ploidy_cutoffs=x_ploidy_cutoffs, y_ploidy_cutoffs=y_ploidy_cutoffs,
    )

    return sex_ht


def main(args):
    hl.init(log="/sex_inference.log", default_reference="GRCh38")

    if args.compute_callrate_ac_af:
        vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True)
        mt = hl.vds.to_dense_mt(vds)
        ht = mt.annotate_rows(
            AN=hl.agg.count_where(hl.is_defined(mt.GT)) * 2,
            AC=hl.agg.sum(mt.GT.n_alt_alleles),
        ).rows()
        ht = ht.annotate(AF=ht.AC/ht.AN)
        ht.write(hard_filtered_ac_an_af.path)
    if args.impute_sex:
        vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True)
        # Added because without this impute_sex_chromosome_ploidy will still run even with overwrite=False
        if args.overwrite or not file_exists(sex.path):
            ht = compute_sex(
                vds,
                args.use_gnomad_methods,
                args.aaf_threshold,
                args.f_stat_cutoff,
                args.high_cov_intervals,
                args.per_platform_high_cov_intervals,
                args.normalization_contig,
                args.x_cov,
                args.y_cov,
                args.norm_cov,
                args.prop_samples_x,
                args.prop_samples_y,
                args.prop_samples_norm,
                args.reference_only,
                args.variants_only,
            )
            ht.write(sex.path, overwrite=True)
    elif args.reannotate_sex:
        reannotate_sex(
            args.min_cov,
            (args.upper_x, (args.lower_xx, args.upper_xx), args.lower_xxx),
            ((args.lower_y, args.upper_y), args.lower_yy),
        ).write(sex.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--impute-sex",
        help="Runs sex imputation. Also runs sex karyotyping annotation.",
        action="store_true",
    )
    parser.add_argument(
        "--use-gnomad-methods",
        help="Runs sex imputation. Also runs sex karyotyping annotation.",
        action="store_true",
    )
    parser.add_argument("--aaf-threshold", help="", type=float, default=0.001)
    parser.add_argument("--f-stat-cutoff", help="", type=float, default=0.5)
    parser.add_argument("--high-cov-intervals", help="", action="store_true")
    parser.add_argument("--per-platform-high-cov-intervals", help="", action="store_true")
    parser.add_argument("--x-cov", help="", type=int, default=10)
    parser.add_argument("--y-cov", help="", type=int, default=5)
    parser.add_argument("--norm-cov", help="", type=int, default=20)
    parser.add_argument("--prop-samples-x", help="", type=float, default=0.80)
    parser.add_argument("--prop-samples-y", help="", type=float, default=0.35)
    parser.add_argument("--prop-samples-norm", help="", type=float, default=0.85)

    parser.add_argument(
        "--reannotate-sex",
        help="Runs the sex karyotyping annotations again, without re-computing sex imputation metrics.",
        action="store_true",
    )
    parser.add_argument("--upper-x", help="Upper cutoff for single X", type=float)
    parser.add_argument("--lower-xx", help="Lower cutoff for double X", type=float)
    parser.add_argument("--upper-xx", help="Upper cutoff for double X", type=float)
    parser.add_argument("--lower-xxx", help="Lower cutoff for triple X", type=float)
    parser.add_argument("--lower-y", help="Lower cutoff for single Y", type=float)
    parser.add_argument("--upper-y", help="Upper cutoff for single Y", type=float)
    parser.add_argument("--lower-yy", help="Lower cutoff for double Y", type=float)

    main(parser.parse_args())
