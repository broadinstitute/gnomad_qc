import argparse
import logging

import hail as hl

from gnomad.sample_qc.pipeline import annotate_sex
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.sample_qc import (
    hard_filtered_ac_an_af,
    interval_coverage,
    platform,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


# TODO: How to incorporate use_gnomad_methods_x_ploidy and use_gnomad_methods_y_ploidy?
#  This means both need to be run if we want x_ploidy to use one method and y_ploidy to use another right?
#  Should we break this up, run and save both methods (x and y), then have another function that actually uses them to
#  determine the cutoffs? Or somehow bake this into the pipeline function? Or maybe add an option into Hail also where
#  x and y can be done in different ways (or one done and not the other)
def compute_sex(
    vds,
    variants_only_x_ploidy: bool = False,
    variants_only_y_ploidy: bool = False,
    high_cov_intervals: bool = False,
    per_platform: bool = False,
    high_cov_by_platform_all: bool = False,
    min_platform_size: bool = 100,
    normalization_contig: str = "chr20",
    x_cov: int = None,
    y_cov: int = None,
    norm_cov: int = None,
    prop_samples_x: float = None,
    prop_samples_y: float = None,
    prop_samples_norm: float = None,
    freq_ht=None,
    aaf_threshold: float = 0.001,
    f_stat_cutoff: float = 0.5,
) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Within this function there are some critical parameters to consider (more details on each are described below):
        - Filtering to intervals with high coverage for sex chromosome ploidy estimation (`high_cov_intervals`, `x_cov`,
         `y_cov`, `norm_cov`, `prop_samples_x`, `prop_samples_y`, and `prop_samples_norm`).
        - Per platform computations: either to determine per platform sex karyotype cutoffs, or if `high_cov_intervals`
         is used, per platform sex chromosome ploidy estimation and determination of sex karyotype cutoffs is performed
         using per platform high quality intervals (`x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`, `prop_samples_y`,
         `prop_samples_norm`, and `min_platform_size`).
        - Use per platform stats to determine which are high quality across all platforms (depending on platform sample
         size). Uses parameters `high_cov_by_platform_all`, `x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
         `prop_samples_y`, `prop_samples_norm`.
        - We use `annotate_sex`, which uses Hail's VDS imputation of sex ploidy `hail.vds.impute_sex_chromosome_ploidy`,
          and uses `variants_only_x_ploidy` and `variants_only_y_ploidy`.
        - Use of a `freq_ht` to filter variants for the X-chromosome heterozygosity computation. This is the f_stat
         annotation applied by Hail's `impute_sex` module. (`freq_ht`, `aaf_threshold`).

    The following are the available options for chrX and chrY relative sex ploidy estimation and sex karyotype
    cutoff determination:
        - Ploidy estimation using all intervals defined by `interval_name` and `padding`. Use options
        - Per platform ploidy estimation using all intervals defined by `interval_name` and `padding`. Sex karyotype
         cutoffs are determined by per platform ploidy distributions.
        - Ploidy estimation using only high coverage intervals across the entire dataset, defined as: chrX intervals
         having a proportion of all samples (`prop_samples_x`) with over a specified mean DP (`x_cov`), chrY intervals
         having a proportion of all samples (`prop_samples_y`) with over a specified mean DP (`y_cov`), and
         `normalization_contig` intervals having a proportion of all samples (`prop_samples_norm`) with over a
         specified mean DP (`norm_cov`). Sex karyotype cutoffs are determined by ploidy distributions of all samples.
        - Per platform ploidy estimation using only per platform high coverage intervals defined by: chrX intervals with
         the specific platform having a proportion of samples (`prop_samples_x`) with over a specified mean DP (`x_cov`),
         chrY intervals with the specific platform having a proportion of samples (`prop_samples_y`) with over a
         specified mean DP (`y_cov`), and `normalization_contig` intervals with the specific platform having a
         proportion of samples (`prop_samples_norm`) with over a specified mean DP (`norm_cov`). Sex karyotype cutoffs
         are determined by per platform ploidy distributions.
        - A combined ploidy estimation using only per platform high quality intervals. Each interval is assessed as high
         quality per platform on chrX, chrY, and `normalization_contig` as described above for "per platform ploidy
         estimation". Then this method uses only platforms that have # of samples > `min_platform_size` to determine
         intervals that have a high coverage across all platforms. Then sex ploidy estimation and sex karyotype cutoffs
         are determined using this intersection of high quality intervals across platforms.

    For each of the options described above, there is also the possibility to use a ploidy estimation that uses only
    variants within the specified calling intervals:
        - This can be defined differently for chrX and chrY using `variants_only_x_ploidy` and `variants_only_y_ploidy`.
        - `variants_only_x_ploidy` is the preferred method for chrX ploidy computation.
        - If not using only variants Hail's `hail.vds.impute_sex_chromosome_ploidy` method will only compute chromosome
         ploidy using reference block DP per calling interval. This method breaks up the reference blocks at the
         calling interval boundries, maintaining all reference block END information for the mean DP per interval
         computation. This is different from the method used on sparse MatrixTables in gnomad_methods
         `gnomad.utils.sparse_mt.impute_sex_ploidy`. In this the chromosome coverage will be computed using all
         reference blocks and variants found in the sparse MatrixTable, but it only uses specified calling intervals to
         determine the contig size and doesn't adjust break up the reference blocks in the same way the Hail method does.

    f-stat is computed on all variants unless a `freq_ht` is defined. Therefore the aaf_threshold is only used when
    `freq_ht` is supplied.

    :param use_gnomad_methods:
    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :return: Table with inferred sex annotation
    :rtype: hl.Table
    """

    def _get_filtered_coverage_mt(coverage_mt, agg_func):
        return coverage_mt.filter_rows(
            (
                (coverage_mt.interval.start.contig == "chrX")
                & agg_func(coverage_mt[f"over_{x_cov}"] > prop_samples_x)
            )
            | (
                (coverage_mt.interval.start.contig == "chrY")
                & agg_func(coverage_mt[f"over_{y_cov}"] > prop_samples_y)
            )
            | (
                (coverage_mt.interval.start.contig == "chr20")
                & agg_func(coverage_mt[f"over_{norm_cov}x"] > prop_samples_norm)
            )
        )

    def _annotate_sex(vds, calling_intervals_ht):
        if use_gnomad_methods_x_ploidy != use_gnomad_methods_y_ploidy:
            calling_intervals_x_ht = calling_intervals_ht.filter(calling_intervals_ht.interval.start.contig != "chrY")
            calling_intervals_y_ht = calling_intervals_ht.filter(calling_intervals_ht.interval.start.contig != "chrX")
        return annotate_sex(
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

    if per_platform or high_cov_by_platform_all:
        platform_ht = platform.ht()
        platforms = platform_ht.aggregate(
            hl.agg.collect_as_set(platform_ht.qc_platform)
        )

    if high_cov_intervals or high_cov_by_platform_all:
        coverage_mt = interval_coverage.mt()
        coverage_mt = coverage_mt.filter_rows(
            hl.literal({"chrY", "chrX", normalization_contig}).contains(
                coverage_mt.interval.start.contig
            )
        )

        if per_platform or high_cov_by_platform_all:
            coverage_mt = coverage_mt.annotate_cols(
                platform=platform_ht[coverage_mt.col_key].qc_platform
            )

        coverage_agg_expr = {
            f"over_{x_cov}x": hl.agg.fraction(coverage_mt.mean_dp > x_cov),
            f"over_{y_cov}x": hl.agg.fraction(coverage_mt.mean_dp > y_cov),
            f"over_{norm_cov}x": hl.agg.fraction(coverage_mt.mean_dp > norm_cov),
        }

        if per_platform or high_cov_by_platform_all:
            coverage_mt = coverage_mt.group_cols_by(coverage_mt.platform).aggregate(
                **coverage_agg_expr
            )

        if per_platform:
            per_platform_sex_hts = []
            for platform in platforms:
                coverage_platform_mt = _get_filtered_coverage_mt(
                    coverage_mt.filter_cols(coverage_mt.platform == platform),
                    hl.agg.any,
                )
                calling_intervals_ht = coverage_platform_mt.rows()
                vds = hl.vds.filter_samples(vds, coverage_platform_mt.cols())
                sex_ht = _annotate_sex(vds, calling_intervals_ht)
                sex_ht = sex_ht.annotate(platform=platform)
                per_platform_sex_hts.append(sex_ht)
            sex_ht = per_platform_sex_hts[0].union(*per_platform_sex_hts[1:])
        elif high_cov_by_platform_all:
            platform_n_ht = platform_ht.group_by(platform_ht.platform).aggregate(
                n_samples=hl.agg.count()
            )
            coverage_mt = coverage_mt.annotate_cols(
                n_samples=platform_n_ht[coverage_mt.col_key].n_samples
            )
            coverage_mt = coverage_mt.filter_cols(
                coverage_mt.n_samples >= min_platform_size
            )
            coverage_mt = _get_filtered_coverage_mt(coverage_mt, hl.agg.all)
            calling_intervals_ht = coverage_mt.rows()
            sex_ht = _annotate_sex(vds, calling_intervals_ht)
        else:
            coverage_mt = coverage_mt.annotate_rows(**coverage_agg_expr)
            coverage_mt = _get_filtered_coverage_mt(coverage_mt, lambda x: x)
            calling_intervals_ht = coverage_mt.rows()
            sex_ht = _annotate_sex(vds, calling_intervals_ht)
    else:
        calling_intervals_ht = calling_intervals(interval_name, padding).ht()
        if per_platform:
            per_platform_sex_hts = []
            for platform in platforms:
                ht = platform_ht.filter(platform_ht.platform == platform)
                vds = hl.vds.filter_samples(vds, ht)
                sex_ht = _annotate_sex(vds, calling_intervals_ht)
                sex_ht = sex_ht.annotate(platform=platform)
                per_platform_sex_hts.append(sex_ht)
            sex_ht = per_platform_sex_hts[0].union(*per_platform_sex_hts[1:])
        else:
            sex_ht = _annotate_sex(vds, calling_intervals_ht)

    return sex_ht


def main(args):
    hl.init(log="/sex_inference.log", default_reference="GRCh38")

    if args.compute_callrate_ac_af:
        vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True)
        mt = hl.vds.to_dense_mt(vds)
        n_samples = mt.count_cols()
        ht = mt.annotate_rows(
            AN=hl.agg.count_where(hl.is_defined(mt.GT)) * 2,
            AC=hl.agg.sum(mt.GT.n_alt_alleles),
        ).rows()
        ht = ht.annotate(AF=ht.AC / ht.AN, callrate=ht.AN / (n_samples * 2))
        ht.write(hard_filtered_ac_an_af.path)

    if args.impute_sex:
        vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True)
        if args.f_stat_high_callrate_common_var:
            # TODO: Determine what variants to use for f-stat calculation, current implementation will use a HT created by
            #  performing a full densify to get callrate and AF, note, this might not be needed, see CCDG f-stat values
            freq_ht = hard_filtered_ac_an_af.ht()
            freq_ht = freq_ht.filter(freq_ht.callrate > args.min_callrate)

        # Added because without this impute_sex_chromosome_ploidy will still run even with overwrite=False
        if args.overwrite or not file_exists(sex.path):
            ht = compute_sex(
                vds,
                args.use_gnomad_methods,
                freq_ht if args.f_stat_high_callrate_common_variants else None,
                args.aaf_threshold,
                args.f_stat_cutoff,
                args.high_cov_intervals,
                args.per_platform,
                args.high_cov_all_platforms,
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
    parser.add_argument(
        "--per-platform-high-cov-intervals", help="", action="store_true"
    )
    parser.add_argument("--x-cov", help="", type=int, default=10)
    parser.add_argument("--y-cov", help="", type=int, default=5)
    parser.add_argument("--norm-cov", help="", type=int, default=20)
    parser.add_argument("--prop-samples-x", help="", type=float, default=0.80)
    parser.add_argument("--prop-samples-y", help="", type=float, default=0.35)
    parser.add_argument("--prop-samples-norm", help="", type=float, default=0.85)

    args = parser.parse_args()
    main(args)

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
