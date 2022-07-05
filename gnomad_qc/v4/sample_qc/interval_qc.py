import argparse
import logging
from typing import Callable

import hail as hl

from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.sample_qc import interval_coverage, platform

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


def compute_interval_qc(
    vds: hl.vds.VariantDataset,
    coverage_mt: hl.MatrixTable,
    per_platform: bool = False,
    high_cov_by_platform_all: bool = False,
    platform_ht: hl.Table = None,
    min_platform_size: bool = 100,
    x_cov: int = None,
    y_cov: int = None,
    autosome_cov: int = None,
    prop_samples_x: float = None,
    prop_samples_y: float = None,
    prop_samples_autosome: float = None,
) -> hl.Table:
    """


    - Use per platform stats to determine which are high coverage intervals across all platforms (depending on
    platform sample size). Uses parameters `high_cov_by_platform_all`, `x_cov`, `y_cov`, `autosome_cov`,
    `prop_samples_x`, `prop_samples_y`, `prop_samples_autosome`.

    - Addition of an annotation indication high coverage intervals across the entire dataset, defined as: chrX intervals
    having a proportion of all samples (`prop_samples_x`) with over a specified mean DP (`x_cov`), chrY intervals
    having a proportion of all samples (`prop_samples_y`) with over a specified mean DP (`y_cov`), and
    autosome intervals having a proportion of all samples (`prop_samples_autosome`) with over a
    specified mean DP (`autosome_cov`).

    - Addition of a per platform high coverage interval annotation. Defined by: chrX intervals with
    the specific platform having a proportion of samples (`prop_samples_x`) with over a specified mean DP (`x_cov`),
    chrY intervals with the specific platform having a proportion of samples (`prop_samples_y`) with over a
    specified mean DP (`y_cov`), and autosomal intervals with the specific platform having a
    proportion of samples (`prop_samples_norm`) with over a specified mean DP (`autosome_cov`). Sex karyotype cutoffs
    are determined by per platform ploidy distributions.

    :param vds: Input VDS for use in sex inference
    :param coverage_mt: Input interval coverage MatrixTable
    :param per_platform: Whether to run the sex ploidy and karyotype inference per platform
    :param high_cov_by_platform_all: Whether to filter to high coverage intervals for the sex ploidy and karyotype
        inference. Using only intervals that are considered high coverage across all platforms
    :param platform_ht: Input platform assignment Table. This is only needed if per_platform or high_cov_by_platform_all are True
    :param min_platform_size: Required size of a platform to be considered when using `high_cov_by_platform_all`. Only
        platforms that have # of samples > 'min_platform_size' are used to determine intervals that have a
        high coverage across all platforms
    :param x_cov: Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the
        interval_coverage MatrixTable
    :param y_cov: Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the
        interval_coverage MatrixTable
    :param autosome_cov: Mean coverage level used to define high coverage intervals on the autosomes. This field must
        be in the interval_coverage MT
    :param prop_samples_x: Proportion samples at specified coverage `x_cov` to determine high coverage intervals on chromosome X
    :param prop_samples_y: Proportion samples at specified coverage `y_cov` to determine high coverage intervals on chromosome Y
    :param prop_samples_autosome: Proportion samples at specified coverage `autosome_cov` to determine high coverage intervals
        on the autosomes
    :return: Table with interval QC annotations
    """

    def _get_filtered_coverage_mt(
        coverage_mt: hl.MatrixTable,
        agg_func: Callable[
            [hl.expr.BooleanExpression], hl.BooleanExpression
        ] = hl.agg.all,
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
                    coverage_mt.filter_cols(coverage_mt.platform == platform)
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
            coverage_mt = _get_filtered_coverage_mt(coverage_mt)
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
        variant_depth_only_x_ploidy=variant_depth_only_x_ploidy,
        variant_depth_only_y_ploidy=variant_depth_only_y_ploidy,
        x_cov=x_cov,
        y_cov=y_cov,
        norm_cov=norm_cov,
        prop_samples_x=prop_samples_x,
        prop_samples_y=prop_samples_y,
        prop_samples_norm=prop_samples_norm,
        f_stat_min_af=min_af,
        f_stat_cutoff=f_stat_cutoff,
    )

    return sex_ht


def main(args):
    hl.init(log="/interval_qc.log", default_reference="GRCh38")

    try:
        vds = get_gnomad_v4_vds(remove_hard_filtered_samples=True, test=args.test)
        logger.info(
            "Loading interval coverage MatrixTable...",
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

        ht = compute_interval_qc(
            vds,
            coverage_mt,
            args.per_platform,
            args.high_cov_by_platform_all,
            platform_ht,
            args.min_platform_size,
            args.x_cov,
            args.y_cov,
            args.autosome_cov,
            args.prop_samples_x,
            args.prop_samples_y,
            args.prop_samples_autosome,
        )

        ht.write(
            get_checkpoint_path("interval_qc")
            if args.test
            else interval_qc.path,
            overwrite=args.overwrite,
        )
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("interval_qc"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Test the pipeline using the gnomAD v4 test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--per-platform",
        help=(
            "Whether to add a per platform interval QC annotation indicating if the interval is considered high "
            "coverage within the samples assigned platform."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--high-cov-by-platform-all",
        help=(
            "Whether to add an annotation to the interval QC MT that indicates if the interval is considered high "
            "coverage across all platforms."
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
        "--autosome-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on the autosomes. This field must "
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
        "--prop-samples-autosomes",
        help=(
            "Proportion samples at specified coverage '--autosome-cov' to determine high coverage intervals on the "
            "autosomes."
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
