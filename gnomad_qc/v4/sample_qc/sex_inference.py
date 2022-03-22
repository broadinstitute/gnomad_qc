import argparse
import logging
from typing import Callable

import hail as hl

from gnomad.sample_qc.pipeline import annotate_sex
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
    ukbb_f_stat,
)
from gnomad_qc.v4.resources.sample_qc import (
    hard_filtered_ac,
    hard_filtered_af_callrate,
    interval_coverage,
    platform,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


def compute_sex(
    vds: hl.vds.VariantDataset,
    coverage_mt: hl.MatrixTable,
    high_cov_intervals: bool = False,
    per_platform: bool = False,
    high_cov_by_platform_all: bool = False,
    platform_ht: hl.Table = None,
    min_platform_size: bool = 100,
    normalization_contig: str = "chr20",
    variants_only_x_ploidy: bool = False,
    variants_only_y_ploidy: bool = False,
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
         using per platform high coverage intervals (`x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`, `prop_samples_y`,
         `prop_samples_norm`, and `min_platform_size`).
        - Use per platform stats to determine which are high coverage across all platforms (depending on platform sample
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
        - A combined ploidy estimation using only per platform high coverage intervals. Each interval is assessed as high
         coverage per platform on chrX, chrY, and `normalization_contig` as described above for "per platform ploidy
         estimation". Then this method uses only platforms that have # of samples > `min_platform_size` to determine
         intervals that have a high coverage across all platforms. Then sex ploidy estimation and sex karyotype cutoffs
         are determined using this intersection of high coverage intervals across platforms.

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

    :param vds: Input VDS for use in sex inference
    :param coverage_mt: Input interval coverage MatrixTable
    :param high_cov_intervals: Whether to filter to high coverage intervals for the sex ploidy and karyotype inference
    :param per_platform: Whether to run the sex ploidy and karyotype inference per platform
    :param high_cov_by_platform_all: Whether to filter to high coverage intervals for the sex ploidy and karyotype
        inference. Using only intervals that are considered high coverage across all platforms
    :param platform_ht: Input platform assignment Table. This is only needed if per_platform or high_cov_by_platform_all are True
    :param min_platform_size: Required size of a platform to be considered when using `high_cov_by_platform_all`. Only
        platforms that have # of samples > 'min_platform_size' are used to determine intervals that have a
        high coverage across all platforms
    :param normalization_contig: Which autosomal chromosome to use for normalizing the coverage of chromosomes X and Y
    :param variants_only_x_ploidy: Whether to use depth of variant data within calling intervals instead of reference
        data for chrX ploidy estimation. Default will only use reference data
    :param variants_only_y_ploidy: Whether to use depth of variant data within calling intervals instead of reference
        data for chrY ploidy estimation. Default will only use reference data
    :param x_cov: Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the
        interval_coverage MatrixTable
    :param y_cov: Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the
        interval_coverage MatrixTable
    :param norm_cov: Mean coverage level used to define high coverage intervals on the normalization autosome
        (`normalization_contig`). This field must be in the interval_coverage MT
    :param prop_samples_x: Proportion samples at specified coverage `x_cov` to determine high coverage intervals on chromosome X
    :param prop_samples_y: Proportion samples at specified coverage `y_cov` to determine high coverage intervals on chromosome Y
    :param prop_samples_norm: Proportion samples at specified coverage `norm_cov` to determine high coverage intervals
        on the normalization chromosome specified by `normalization_contig`
    :param freq_ht: Table to use for f-stat allele frequency cutoff. The input VDS is filtered to sites in this Table
        prior to running Hail's `impute_sex` module, and alternate allele frequency is used from this Table with a
        `aaf_threshold` cutoff
    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff
    :return: Table with inferred sex annotation
    """

    def _get_filtered_coverage_mt(
        coverage_mt: hl.MatrixTable,
        agg_func: Callable[[hl.expr.BooleanExpression], hl.BooleanExpression],
    ) -> hl.MatrixTable:
        """
        Helper function to filter the interval coverage MatrixTable to high coverage intervals.

        High coverage intervals are determined using `agg_func`, `x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
        `prop_samples_y`, and `prop_samples_norm`.

        :param coverage_mt: Input interval coverage MatrixTable
        :param agg_func: Hail aggregation function to determine if an interval coverage meets the `cov_*`> `prop_samples_*` criteria
        :return: Filtered interval coverage MatrixTable
        """
        return coverage_mt.filter_rows(
            (
                (coverage_mt.interval.start.contig == "chrX")
                & agg_func(coverage_mt[f"over_{x_cov}x"] > prop_samples_x)
            )
            | (
                (coverage_mt.interval.start.contig == "chrY")
                & agg_func(coverage_mt[f"over_{y_cov}x"] > prop_samples_y)
            )
            | (
                (coverage_mt.interval.start.contig == "chr20")
                & agg_func(coverage_mt[f"over_{norm_cov}x"] > prop_samples_norm)
            )
        )

    def _annotate_sex(vds, calling_intervals_ht):
        """
        Helper function to perform `annotate_sex` using unchanged parameters with changes to the VDS and calling intervals.

        :param vds: Input VDS to use for sex annotation
        :param calling_intervals_ht: Calling intervals to filter to for sex annotation
        :return: Table containing sex annotation for samples in the input VDS
        """
        return annotate_sex(
            vds,
            included_intervals=calling_intervals_ht,
            normalization_contig=normalization_contig,
            sites_ht=freq_ht,
            aaf_expr="AF",
            gt_expr="LGT",
            f_stat_cutoff=f_stat_cutoff,
            aaf_threshold=aaf_threshold,
            variants_only_x_ploidy=variants_only_x_ploidy,
            variants_only_y_ploidy=variants_only_y_ploidy,
        )

    if platform_ht is not None:
        logger.info("Collecting platform list...")
        platforms = platform_ht.aggregate(
            hl.agg.collect_as_set(platform_ht.qc_platform)
        )

    if high_cov_intervals or high_cov_by_platform_all:
        logger.info(
            "Running sex ploidy and sex karyotype estimation using only high coverage intervals: %d percent of samples "
            "with greater than %dx coverage on chrX, %d percent of samples with greater than %dx coverage on chrY, and "
            "%d percent of samples with greater than %dx coverage on %s...",
            int(prop_samples_x * 100),
            x_cov,
            int(prop_samples_y * 100),
            y_cov,
            int(prop_samples_norm * 100),
            norm_cov,
            normalization_contig,
        )

        if per_platform or high_cov_by_platform_all:
            logger.info("Annotating coverage MatrixTable with platform information")
            coverage_mt = coverage_mt.annotate_cols(
                platform=platform_ht[coverage_mt.col_key].qc_platform
            )
            coverage_mt = coverage_mt.group_cols_by(coverage_mt.platform).aggregate(
                **{
                    f"over_{x_cov}x": hl.agg.fraction(coverage_mt.mean_dp > x_cov),
                    f"over_{y_cov}x": hl.agg.fraction(coverage_mt.mean_dp > y_cov),
                    f"over_{norm_cov}x": hl.agg.fraction(
                        coverage_mt.mean_dp > norm_cov
                    ),
                }
            )

        if per_platform:
            logger.info(
                "Running sex ploidy and sex karyotype estimation per platform using per platform "
                "high coverage intervals..."
            )
            per_platform_sex_hts = []
            x_ploidy_cutoffs = {}
            y_ploidy_cutoffs = {}
            for platform in platforms:
                ht = platform_ht.filter(platform_ht.qc_platform == platform)
                logger.info("Platform %s has %d samples...", platform, ht.count())
                coverage_platform_mt = _get_filtered_coverage_mt(
                    coverage_mt.filter_cols(coverage_mt.platform == platform),
                    hl.agg.any,
                )
                calling_intervals_ht = coverage_platform_mt.rows()
                sex_ht = _annotate_sex(
                    hl.vds.filter_samples(vds, ht), calling_intervals_ht
                )
                sex_ht = sex_ht.annotate(platform=platform)
                per_platform_sex_hts.append(sex_ht)
                x_ploidy_cutoffs[platform] = sex_ht.index_globals().x_ploidy_cutoffs
                y_ploidy_cutoffs[platform] = sex_ht.index_globals().y_ploidy_cutoffs
            sex_ht = per_platform_sex_hts[0].union(*per_platform_sex_hts[1:])
            sex_ht = sex_ht.annotate_globals(
                x_ploidy_cutoffs=hl.struct(**x_ploidy_cutoffs),
                y_ploidy_cutoffs=hl.struct(**y_ploidy_cutoffs),
            )
        elif high_cov_by_platform_all:
            logger.info(
                "Running sex ploidy and sex karyotype estimation using high coverage intervals across all platforms. "
                "Limited to platforms with at least %s samples...",
                min_platform_size,
            )
            platform_n_ht = platform_ht.group_by(platform_ht.qc_platform).aggregate(
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
            logger.info(
                "Running sex ploidy and sex karyotype estimation using high coverage intervals across the full sample set..."
            )
            coverage_mt = coverage_mt.annotate_rows(
                **{
                    f"over_{x_cov}x": hl.agg.fraction(coverage_mt.mean_dp > x_cov),
                    f"over_{y_cov}x": hl.agg.fraction(coverage_mt.mean_dp > y_cov),
                    f"over_{norm_cov}x": hl.agg.fraction(
                        coverage_mt.mean_dp > norm_cov
                    ),
                }
            )
            coverage_mt = _get_filtered_coverage_mt(coverage_mt, lambda x: x)
            calling_intervals_ht = coverage_mt.rows()
            sex_ht = _annotate_sex(vds, calling_intervals_ht)
    else:
        calling_intervals_ht = coverage_mt.rows()
        if per_platform:
            logger.info(
                "Running sex ploidy estimation and per platform sex karyotype estimation..."
            )
            per_platform_sex_hts = []
            x_ploidy_cutoffs = {}
            y_ploidy_cutoffs = {}
            for platform in platforms:
                ht = platform_ht.filter(platform_ht.qc_platform == platform)
                logger.info("Platform %s has %d samples...", platform, ht.count())
                sex_ht = _annotate_sex(
                    hl.vds.filter_samples(vds, ht), calling_intervals_ht
                )
                sex_ht = sex_ht.annotate(platform=platform)
                per_platform_sex_hts.append(sex_ht)
                x_ploidy_cutoffs[platform] = sex_ht.index_globals().x_ploidy_cutoffs
                y_ploidy_cutoffs[platform] = sex_ht.index_globals().y_ploidy_cutoffs
            sex_ht = per_platform_sex_hts[0].union(*per_platform_sex_hts[1:])
            sex_ht = sex_ht.annotate_globals(
                x_ploidy_cutoffs=hl.struct(**x_ploidy_cutoffs),
                y_ploidy_cutoffs=hl.struct(**y_ploidy_cutoffs),
            )
        else:
            logger.info("Running sex ploidy and sex karyotype estimation...")
            sex_ht = _annotate_sex(vds, calling_intervals_ht)

        sex_ht = sex_ht.annotate_globals(
            high_cov_intervals=high_cov_intervals,
            per_platform=per_platform,
            high_cov_by_platform_all=high_cov_by_platform_all,
            min_platform_size=min_platform_size,
            normalization_contig=normalization_contig,
            variants_only_x_ploidy=variants_only_x_ploidy,
            variants_only_y_ploidy=variants_only_y_ploidy,
            x_cov=x_cov,
            y_cov=y_cov,
            norm_cov=norm_cov,
            prop_samples_x=prop_samples_x,
            prop_samples_y=prop_samples_y,
            prop_samples_norm=prop_samples_norm,
            aaf_threshold=aaf_threshold,
            f_stat_cutoff=f_stat_cutoff,
        )

    return sex_ht


def main(args):
    hl.init(log="/sex_inference.log", default_reference="GRCh38")

    try:
        if args.compute_ac:
            vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True, test=args.test)
            mt = vds.variant_data
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)
            n_samples = mt.count_cols()
            logger.info("Number of samples: %d", n_samples)
            ht = mt.annotate_rows(
                AC=hl.agg.sum(mt.LGT.n_alt_alleles()),
            ).rows()
            ht = ht.annotate_globals(n_samples=n_samples)
            ht.write(
                get_checkpoint_path("test_ac") if args.test else hard_filtered_ac.path,
                overwrite=args.overwrite,
            )
        if args.compute_af_callrate:
            if args.test:
                ac_ht = hl.read_table(get_checkpoint_path("test_ac"))
            else:
                ac_ht = hard_filtered_ac.ht()
            vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True, test=args.test)
            mt = hl.vds.to_dense_mt(vds)
            mt = mt.filter_rows(hl.len(mt.alleles) == 2)
            n_samples = ac_ht.index_globals().n_samples
            ht = mt.annotate_rows(
                AN=hl.agg.count_where(hl.is_defined(mt.LGT)) * 2,
                AC=ac_ht[mt.row_key].AC,
            ).rows()
            ht = ht.annotate(AF=ht.AC / ht.AN, callrate=ht.AN / (n_samples * 2))
            ht.write(
                get_checkpoint_path("test_af_callrate")
                if args.test
                else hard_filtered_af_callrate.path,
                overwrite=args.overwrite,
            )
        if args.impute_sex:
            vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True, test=args.test)
            if args.f_stat_high_callrate_common_var and args.f_stat_ukbb_var:
                raise ValueError(
                    "f_stat_high_callrate_common_var and f_stat_ukbb_var can't be used together."
                )
            if args.f_stat_high_callrate_common_var:
                freq_ht = (
                    hl.read_table(get_checkpoint_path("test_af_callrate"))
                    if args.test
                    else hard_filtered_af_callrate.ht()
                )
                freq_ht = freq_ht.filter(freq_ht.callrate > args.min_callrate)
            elif args.f_stat_ukbb_var:
                freq_ht = ukbb_f_stat.ht()
            else:
                freq_ht = (
                    hl.read_table(get_checkpoint_path("test_ac"))
                    if args.test
                    else hard_filtered_ac.ht()
                )
                freq_ht = freq_ht.annotate(AF=freq_ht.AC / (freq_ht.n_samples * 2))

            # Added because without this impute_sex_chromosome_ploidy will still run even with overwrite=False
            if args.overwrite or not file_exists(
                get_checkpoint_path("sex_imputation") if args.test else sex.path
            ):
                logger.info(
                    "Loading interval coverage MatrixTable and filtering to chrX, chrY and %s...",
                    args.normalization_contig,
                )
                if file_exists(interval_coverage.path):
                    coverage_mt = interval_coverage.mt()
                elif args.test:
                    test_coverage_path = get_checkpoint_path(
                        f"test_interval_coverage.{args.calling_interval_name}.pad{args.calling_interval_padding}",
                        mt=True,
                    )
                    if file_exists(test_coverage_path):
                        coverage_mt = hl.read_matrix_table(test_coverage_path)
                    else:
                        raise FileNotFoundError(
                            f"There is no final coverage MatrixTable written and a test interval coverage MatrixTable does "
                            f"not exist for calling interval {args.calling_interval_name} and interval padding "
                            f"{args.calling_interval_padding}. Please run platform_inference.py --compute_coverage with the "
                            f"--test argument and needed --calling_interval_name/--calling_interval_padding arguments."
                        )
                else:
                    raise FileNotFoundError(
                        f"There is no final coverage MatrixTable written. Please run: "
                        f"platform_inference.py --compute_coverage to compute the interval coverage MatrixTable."
                    )

                coverage_mt = coverage_mt.filter_rows(
                    hl.literal({"chrY", "chrX", args.normalization_contig}).contains(
                        coverage_mt.interval.start.contig
                    )
                )
                if args.per_platform or args.high_cov_by_platform_all:
                    logger.info("Loading platform information...")
                    if file_exists(platform.path):
                        platform_ht = platform.ht()
                    elif args.test:
                        test_platform_path = get_checkpoint_path(
                            f"test_platform_assignment.{args.calling_interval_name}.pad{args.calling_interval_padding}"
                        )
                        if file_exists(test_platform_path):
                            platform_ht = hl.read_table(test_platform_path)
                        else:
                            raise FileNotFoundError(
                                f"There is no final platform assignment Table written and a test platform assignment Table "
                                f"does not exist for calling interval {args.calling_interval_name} and interval padding "
                                f"{args.calling_interval_padding}. Please run platform_inference.py --assign_platforms "
                                f"with the --test argument and needed --calling_interval_name/--calling_interval_padding "
                                f"arguments."
                            )
                    else:
                        raise FileNotFoundError(
                            f"There is no final platform assignment Table written. Please run: "
                            f"platform_inference.py --assign_platforms to compute the platform assignment Table."
                        )
                else:
                    platform_ht = None

                ht = compute_sex(
                    vds,
                    coverage_mt,
                    args.high_cov_intervals,
                    args.per_platform,
                    args.high_cov_by_platform_all,
                    platform_ht,
                    args.min_platform_size,
                    args.normalization_contig,
                    args.variants_only_x_ploidy,
                    args.variants_only_y_ploidy,
                    args.x_cov,
                    args.y_cov,
                    args.norm_cov,
                    args.prop_samples_x,
                    args.prop_samples_y,
                    args.prop_samples_norm,
                    freq_ht,
                    args.min_af,
                    args.f_stat_cutoff,
                )
                ht = ht.annotate_globals(
                    f_stat_high_callrate_common_var=args.f_stat_high_callrate_common_var,
                    f_stat_ukbb_var=args.f_stat_ukbb_var,
                )
                ht.write(
                    get_checkpoint_path(
                        f"sex_imputation{'.per_platform' if args.per_platform else ''}"
                        f"{'.high_cov_by_platform_all' if args.high_cov_by_platform_all else ''}"
                        f"{'.high_cov' if args.high_cov_intervals else ''}"
                        f"{'.ukbb_f_stat' if args.f_stat_ukbb_var else ''}"
                        f"{'.callrate_common_f_stat' if args.f_stat_high_callrate_common_var else ''}"
                    )
                    if args.test
                    else sex.path,
                    overwrite=True,
                )
            else:
                logger.warning("File exists and overwrite is not set!")
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("sex_inference"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Test the pipeline using the gnomAD v4 test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--compute-ac",
        help="Compute the allele count, for all bi-allelic variants in the VDS.",
        action="store_true",
    )
    parser.add_argument(
        "--compute-af-callrate",
        help=(
            "Compute the callrate and allele frequency for all variants in the VDS. NOTE: This requires a full densify!!"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--impute-sex",
        help="Runs sex ploidy and sex karyotyping imputation.",
        action="store_true",
    )
    parser.add_argument(
        "--f-stat-high-callrate-common-var",
        help=(
            "Whether to use high callrate (> min_callrate) and common variants (> aaf-threshold) for f-stat computation. "
            "By default, no callrate cutoff will be used, and allele frequency will be approximated with AC/(n_samples * 2)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--f-stat-ukbb-var",
        help=(
            "Whether to use UK Biobank high callrate and common variants (UKBB allele frequency > aaf-threshold) for "
            "f-stat computation. By default, no callrate cutoff will be used, and allele frequency will be approximated "
            "with AC/(n_samples * 2)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min_callrate", help="Minimum variant callrate.", default=0.99, type=float
    )
    parser.add_argument(
        "--min-af",
        help="Minimum variant allele frequency to retain variant in qc matrix table.",
        default=0.001,
        type=float,
    )
    parser.add_argument(
        "--f-stat-cutoff",
        help=(
            "Cutoff for f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY "
            "are above cutoff."
        ),
        type=float,
        default=0.5,
    )
    parser.add_argument(
        "--high-cov-intervals",
        help="Whether to filter to high coverage intervals for the sex ploidy and karyotype inference.",
        action="store_true",
    )
    parser.add_argument(
        "--per-platform",
        help="Whether to run the sex ploidy and karyotype inference per platform.",
        action="store_true",
    )
    parser.add_argument(
        "--high-cov-by-platform-all",
        help=(
            "Whether to filter to high coverage intervals for the sex ploidy and karyotype inference. "
            "Using only intervals that are considered high coverage across all platforms."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-platform-size",
        help=(
            "Required size of a platform to be considered in '--high-cov-by-platform-all'. Only platforms that "
            "have # of samples > 'min_platform_size' are used to determine intervals that have a high coverage across "
            "all platforms."
        ),
        type=int,
        default=100,
    )
    parser.add_argument(
        "--normalization-contig",
        help="Which autosomal chromosome to use for normalizing the coverage of chromosomes X and Y.",
        type=str,
        default="chr20",
    )
    parser.add_argument(
        "--variants-only-x-ploidy",
        help=(
            "Whether to use depth of variant data for the x ploidy estimation instead of the default behavior that "
            "will use reference blocks."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--variants-only-y-ploidy",
        help=(
            "Whether to use depth of variant data for the y ploidy estimation instead of the default behavior that "
            "will use reference blocks."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--x-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the "
            "interval_coverage MT!"
        ),
        type=int,
        default=10,
    )
    parser.add_argument(
        "--y-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the "
            "interval_coverage MT!"
        ),
        type=int,
        default=5,
    )
    parser.add_argument(
        "--norm-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on the normalization autosome. This field must "
            "be in the interval_coverage MT!"
        ),
        type=int,
        default=20,
    )
    parser.add_argument(
        "--prop-samples-x",
        help="Proportion samples at specified coverage '--x-cov' to determine high coverage intervals on chromosome X.",
        type=float,
        default=0.80,
    )
    parser.add_argument(
        "--prop-samples-y",
        help="Proportion samples at specified coverage '--y-cov' to determine high coverage intervals on chromosome Y.",
        type=float,
        default=0.35,
    )
    parser.add_argument(
        "--prop-samples-norm",
        help=(
            "Proportion samples at specified coverage '--norm-cov' to determine high coverage intervals on the "
            "normalization chromosome specified by '--normalization-contig'."
        ),
        type=float,
        default=0.85,
    )
    parser.add_argument(
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'. Only "
            "used if '--test' is set and final coverage MT and/or platform assignment HT are not already written."
        ),
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    parser.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or 50 bp. Only used if '--test' is "
            "set and final coverage MT and/or platform assignment HT are not already written."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()
    main(args)

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
