import argparse
import logging
from typing import List, Tuple, Union

from gnomad.utils.slack import slack_notifications
import hail as hl

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    calling_intervals,
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.sample_qc import (
    hard_filtered_samples,
    interval_coverage,
    interval_qc,
    interval_qc_pass,
    platform,
    sex,
    sex_chr_coverage,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def generate_sex_chr_interval_coverage_mt(
    vds: hl.vds.VariantDataset,
    calling_intervals_ht: hl.Table,
) -> hl.MatrixTable:
    """
    Create a MatrixTable of interval-by-sample coverage on sex chromosomes with intervals split at PAR regions.

    :param vds: Input VariantDataset.
    :param calling_intervals_ht: Calling interval Table.
    :return: MatrixTable with interval coverage per sample on sex chromosomes.
    """
    contigs = ["chrX", "chrY"]
    calling_intervals_ht = calling_intervals_ht.filter(
        hl.literal(contigs).contains(calling_intervals_ht.interval.start.contig)
    )
    logger.info(
        "Filtering VariantDataset to the following contigs: %s...",
        ", ".join(contigs),
    )
    vds = hl.vds.filter_chromosomes(vds, keep=contigs)
    rg = vds.reference_data.locus.dtype.reference_genome

    par_boundaries = []
    for par_interval in rg.par:
        par_boundaries.append(par_interval.start)
        par_boundaries.append(par_interval.end)

    # Segment on PAR interval boundaries
    calling_intervals = hl.segment_intervals(calling_intervals_ht, par_boundaries)

    # Annotate intervals overlapping PAR
    calling_intervals = calling_intervals.annotate(
        overlap_par=hl.any(
            lambda x: x.overlaps(calling_intervals.interval), hl.literal(rg.par)
        )
    )

    kept_contig_filter = hl.array(contigs).map(
        lambda x: hl.parse_locus_interval(x, reference_genome=rg)
    )
    vds = hl.vds.VariantDataset(
        hl.filter_intervals(vds.reference_data, kept_contig_filter),
        hl.filter_intervals(vds.variant_data, kept_contig_filter),
    )
    mt = hl.vds.interval_coverage(vds, calling_intervals)
    mt = mt.annotate_rows(overlap_par=calling_intervals[mt.row_key].overlap_par)

    return mt


def filter_to_test(
    mt: hl.MatrixTable,
    sex_mt: hl.MatrixTable,
    num_partitions: int = 10,
) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    """
    Filter `mt` to `num_partitions` partitions on chr1 and `sex_mt` to `num_partitions` partitions on chrX and all chrY.

    :param mt: Input MatrixTable to filter to specified number of partitions on chr1.
    :param sex_mt: Input MatrixTable to filter to specified number of partitions on chrX and all of chrY.
    :param num_partitions: Number of partitions to grab from mt.
    :return: Input MatrixTables filtered to `num_partitions` on chr1, chrX, and all of chrY.
    """
    logger.info(
        "Filtering to columns in both the coverage MT and the sex coverage MT for testing...",
    )
    mt = mt.semi_join_cols(sex_mt.cols())
    sex_mt = sex_mt.semi_join_cols(mt.cols())

    logger.info(
        "Filtering to %d partitions on chr1, chrX, and all of chrY for testing...",
        num_partitions,
    )
    mt = mt._filter_partitions(range(num_partitions))
    sex_mt_chrx = sex_mt._filter_partitions(range(num_partitions))
    sex_mt_chry = sex_mt.filter_rows(sex_mt.interval.start.contig == "chrY")

    return mt, sex_mt_chrx.union_rows(sex_mt_chry)


def compute_interval_qc(
    mt: hl.MatrixTable,
    platform_ht: hl.Table,
    mean_dp_thresholds: List[int] = [5, 10, 15, 20, 25],
    split_by_sex: bool = False,
) -> hl.Table:
    """
    Compute interval QC aggregate statistics per interval, per platform, and optionally split by sex karyotype.

    The following annotations must be on `mt` (interval-by-sample MT output from `hl.vds.interval_coverage`):
        - interval - Genomic interval of interest.
        - mean_dp - Mean depth of bases across the interval.
        - fraction_over_dp_threshold - Fraction of interval (in bases) above each DP threshold. Second element must be
        dp >= 1 (dp > 0).

        If `split_by_sex`:
            - sex_karyotype - StringExpression annotation with sex karyotype information including 'XX' and 'XY' values.

    The `platform_ht` must have a 'qc_platform' annotation indicating the platform each sample was assigned.

    Returns a Table with the following annotations:
        - interval_mean_dp - Mean DP of the interval across 'all' samples and optionally split by 'XX' and 'XY'.
        - fraction_over_{dp}x - for all 'dp' in `mean_dp_thresholds`, which is the fraction of samples with mean DP
            over 'dp'. Computed across 'all' samples and optionally split by 'XX' and 'XY'.
        - mean_fraction_over_dp_0 - Mean of the fraction of the interval (in bases) that is dp > 0. Computed across
            'all' samples and optionally split by 'XX' and 'XY'.
        - platform_interval_mean_dp - Same as 'interval_mean_dp', but instead of a single value for 'all', 'XX' and
            'XY' each is a dict of per platform values.
        - platform_fraction_over_{dp}x  - Same as 'fraction_over_{dp}x', but instead of a single value for 'all',
            'XX' and 'XY' each is a dict of per platform values.
        - platform_mean_fraction_over_dp_0 - Same as 'mean_fraction_over_dp_0', but instead of a single value for
            'all', 'XX' and 'XY' each is a dict of per platform values.

    :param mt: Input interval coverage MatrixTable.
    :param platform_ht: Input platform assignment Table.
    :param mean_dp_thresholds: List of mean DP thresholds to use for computing the fraction of samples with mean
        interval DP >= the threshold.
    :param split_by_sex: Whether the interval QC should be stratified by sex. If True, `mt` must be annotated with
        'sex_karyotype'.
    :return: Table with interval QC annotations.
    """

    def _get_agg_expr(expr, agg_func=hl.agg.mean, group_by=None):
        """
        Helper function to call `agg_func` on `expr` with optional stratification by sex karyotype and `group_by`.

        :param expr: Expression to pass to `agg_func`.
        :param agg_func: Function to use for aggregation of `expr`. Default is hl.agg.mean.
        :return: Aggregation expression.
        """
        if group_by is not None:
            agg_func = hl.agg.group_by(group_by, agg_func(expr))
        else:
            agg_func = agg_func(expr)
        agg_expr = {"all": agg_func}
        if split_by_sex:
            agg_expr.update(
                {
                    "XX": hl.agg.filter(mt.sex_karyotype == "XX", agg_func),
                    "XY": hl.agg.filter(mt.sex_karyotype == "XY", agg_func),
                }
            )

        return agg_expr

    num_samples = mt.count_cols()
    mt = mt.annotate_cols(platform=platform_ht[mt.col_key].qc_platform)
    mt = mt.filter_cols(hl.is_defined(mt.platform))
    num_samples_no_platform = num_samples - mt.count_cols()
    if num_samples_no_platform > 0:
        logger.warning(
            "Number of samples in MT with no platform assignment: %d",
            num_samples - mt.count_cols(),
        )

    # Note: Default hl.vds.interval_coverage will return a list for 'fraction_over_dp_threshold' where the
    # second element is dp >= 1 (dp > 0)
    agg_groups = [("", None), ("platform_", mt.platform)]
    mt = mt.annotate_rows(
        **{
            f"{prefix}interval_mean_dp": _get_agg_expr(mt.mean_dp, group_by=group_by)
            for prefix, group_by in agg_groups
        },
        **{
            f"{prefix}fraction_over_{dp}x": _get_agg_expr(
                mt.mean_dp >= dp, agg_func=hl.agg.fraction, group_by=group_by
            )
            for dp in mean_dp_thresholds
            for prefix, group_by in agg_groups
        },
        **{
            f"{prefix}mean_fraction_over_dp_0": _get_agg_expr(
                mt.fraction_over_dp_threshold[1], group_by=group_by
            )
            for prefix, group_by in agg_groups
        },
    )

    mt = mt.annotate_globals(
        mean_dp_thresholds=mean_dp_thresholds,
        platform_n_samples=mt.aggregate_cols(
            hl.agg.group_by(mt.platform, hl.agg.count())
        ),
    )

    return mt.rows()


def get_high_qual_cutoff_dict(
    autosome_par_cutoff,
    x_nonpar_cutoff,
    y_nonpar_cutoff,
    autosome_par_qc_ann,
    x_nonpar_qc_ann,
    y_nonpar_qc_ann,
    split_by_sex: bool = False,
):
    if split_by_sex:
        xx = "XX"
        xy = "XY"
        autosome_par = "all"
    else:
        xx = xy = autosome_par = "all"

    autosome_par_cutoff = (autosome_par_qc_ann, autosome_par, autosome_par_cutoff)
    x_nonpar_cutoff = (x_nonpar_qc_ann, xx, x_nonpar_cutoff)
    y_nonpar_cutoff = (y_nonpar_qc_ann, xy, y_nonpar_cutoff)
    high_qual_cutoffs = {
        "autosome_par": [autosome_par_cutoff],
        "x_non_par": [x_nonpar_cutoff],
        "y_non_par": [y_nonpar_cutoff],
    }

    if split_by_sex:
        high_qual_cutoffs["x_non_par"].append(y_nonpar_cutoff)

    return high_qual_cutoffs


def get_interval_qc_pass(
    interval_qc_ht: hl.Table,
    high_qual_cutoffs,
    per_platform: bool = False,
    all_platforms: bool = False,
    min_platform_size: int = 100,
) -> hl.Table:
    """
    Add `interval_qc_pass` annotation to indicate whether the site falls within a high coverage interval.

    :param interval_qc_ht: Input interval QC Table.
    :param high_qual_cutoffs:
    :param per_platform: Whether to make the interval QC pass annotation a DictionaryExpression with interval QC pass
        per platform.
    :param all_platforms: Whether to consider an interval as passing QC only if it passes interval QC per platform
        across all platforms (with a sample size above `min_platform_size`).
    :param min_platform_size: Required size of a platform to be considered in `all_platforms`. Only platforms that
        have # of samples > 'min_platform_size' are used to determine intervals that have a high coverage across all
        platforms.
    :return: MatrixTable or Table with samples removed
    """
    if per_platform and all_platforms:
        raise ValueError("Only one of 'per_platform' and 'all_platforms' can be True!")

    interval_start = interval_qc_ht.interval.start
    region_exprs = {
        "autosome_par": interval_start.in_autosome_or_par(),
        "x_non_par": interval_start.in_x_nonpar(),
        "y_non_par": interval_start.in_y_nonpar(),
    }

    if per_platform or all_platforms:
        platform_n_samples = (
            interval_qc_ht.index_globals().platform_n_samples.collect()[0]
        )
        ann_prefix = "platform_"
    else:
        ann_prefix = ""

    interval_qc_globals = hl.struct(
        per_platform=per_platform,
        all_platforms=all_platforms,
        high_qual_cutoffs=high_qual_cutoffs,
    )

    def _get_qc_ann_expr(expr, platform=None):
        if platform is None:
            return expr
        return expr.get(platform, None)

    def _get_pass_expr(platform=None):
        return functools.reduce(
            operator.ior,
            [
                region_exprs[region]
                & functools.reduce(
                    operator.iand,
                    [
                        _get_qc_ann_expr(
                            interval_qc_ht[f"{ann_prefix}{ann}"].get(sex, None),
                            platform,
                        )
                        > cutoff
                        for ann, sex, cutoff in ann_cutoffs
                    ],
                )
                for region, ann_cutoffs in high_qual_cutoffs.items()
            ],
        )

    if per_platform or all_platforms:
        pass_interval_qc = hl.struct(
            **{platform: _get_pass_expr(platform) for platform in platform_n_samples}
        )
        if all_platforms:
            interval_qc_globals = interval_qc_globals.annotate(
                min_platform_size=min_platform_size
            )
            platforms = [
                platform
                for platform, n_samples in platform_n_samples.items()
                if (n_samples >= min_platform_size) & (platform != "platform_-1")
            ]
            pass_interval_qc = hl.all(
                [pass_interval_qc[platform] for platform in platforms]
            )
    else:
        pass_interval_qc = _get_pass_expr()

    interval_qc_ht = interval_qc_ht.select(pass_interval_qc=pass_interval_qc)
    interval_qc_ht = interval_qc_ht.select_globals(
        pass_interval_qc_parameters=interval_qc_globals
    )

    return interval_qc_ht


def annotate_interval_qc_filter(
    t: Union[hl.MatrixTable, hl.Table],
    **kwargs,
) -> Union[hl.MatrixTable, hl.Table]:
    """
    Wrapper function around `get_interval_qc_pass` to annotate a Table/MatrixTable with 'pass_interval_qc'.

    Passes the interval QC Table resource and kwargs to `get_interval_qc_pass`.

    :param t: Input Table or MatrixTable.
    :param kwargs: Optional keyword arguments to pass to `get_interval_qc_pass`.
    :return: Input Table or MatrixTable annotated with 'pass_interval_qc'.
    """
    interval_qc_ht = interval_qc.ht()
    interval_qc_ht = get_interval_qc_pass(interval_qc_ht, **kwargs)

    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(interval_qc_pass=interval_qc_ht[t.locus].pass_interval_qc)
    else:
        t = t.annotate(interval_qc_pass=interval_qc_ht[t.locus].pass_interval_qc)

    return t


def main(args):
    hl.init(
        log="/interval_qc.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test = args.test
    calling_interval_name = args.calling_interval_name
    calling_interval_padding = args.calling_interval_padding
    overwrite = args.overwrite
    mean_dp_thresholds = args.mean_dp_thresholds

    try:
        if args.sex_chr_interval_coverage:
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            calling_intervals_ht = calling_intervals(
                calling_interval_name, calling_interval_padding
            ).ht()
            sex_coverage_mt = generate_sex_chr_interval_coverage_mt(
                vds,
                calling_intervals_ht,
            )
            sex_coverage_mt = sex_coverage_mt.annotate_globals(
                calling_interval_name=calling_interval_name,
                calling_interval_padding=calling_interval_padding,
            )
            sex_coverage_mt.write(
                get_checkpoint_path("test_sex_imputation_cov", mt=True)
                if test
                else sex_chr_coverage.path,
                overwrite=overwrite,
            )

        if args.generate_interval_qc_ht:
            platform_ht = platform.ht()
            coverage_mt = interval_coverage.mt()

            if test:
                coverage_mt, sex_coverage_mt = filter_to_test(
                    coverage_mt,
                    hl.read_matrix_table(
                        get_checkpoint_path("test_sex_imputation_cov", mt=True)
                    ),
                )
                coverage_mt = coverage_mt.checkpoint(
                    get_checkpoint_path("interval_qc_coverage", mt=True),
                    _read_if_exists=True,
                )
                sex_coverage_mt = sex_coverage_mt.checkpoint(
                    get_checkpoint_path("interval_qc_sex_coverage", mt=True),
                    _read_if_exists=True,
                )
            else:
                sex_coverage_mt = sex_chr_coverage.mt()

            logger.info("Removing hard-filtered samples from the coverage MTs...")
            coverage_mt = coverage_mt.filter_cols(
                hl.is_missing(hard_filtered_samples.ht()[coverage_mt.col_key])
            )
            sex_coverage_mt = sex_coverage_mt.filter_cols(
                hl.is_missing(hard_filtered_samples.ht()[sex_coverage_mt.col_key])
            )

            logger.info("Computing interval QC on autosomes...")
            coverage_mt = coverage_mt.filter_rows(
                coverage_mt.interval.start.in_autosome()
            )
            coverage_mt = coverage_mt.annotate_rows(overlap_par=False)
            ht = compute_interval_qc(
                coverage_mt,
                platform_ht=platform_ht,
                mean_dp_thresholds=mean_dp_thresholds,
            )

            logger.info("Filtering to XX and XY samples...")
            sex_ht = sex.ht().select("sex_karyotype")
            sex_coverage_mt = sex_coverage_mt.annotate_cols(
                **sex_ht[sex_coverage_mt.col_key]
            )
            sex_coverage_mt = sex_coverage_mt.filter_cols(
                (sex_coverage_mt.sex_karyotype == "XX")
                | (sex_coverage_mt.sex_karyotype == "XY")
            )

            logger.info(
                "Computing interval QC on sex chromosomes and joining with autosome interval QC HT..."
            )
            ht = ht.union(
                compute_interval_qc(
                    sex_coverage_mt,
                    platform_ht=platform_ht,
                    mean_dp_thresholds=mean_dp_thresholds,
                    split_by_sex=True,
                )
            )
            ht.write(
                get_checkpoint_path("interval_qc") if test else interval_qc.path,
                overwrite=overwrite,
            )
        if args.generate_interval_qc_pass_ht:
            ht = (
                hl.read_table(get_checkpoint_path("interval_qc"))
                if test
                else interval_qc.ht()
            )
            bla = {}
            if args.by_mean_fraction_over_dp_0:
                bla["autosome_par_cutoff"] = bla["x_nonpar_cutoff"] = bla[
                    "y_nonpar_cutoff"
                ] = args.mean_fraction_over_dp_0
                bla["autosome_par_qc_ann"] = bla["x_nonpar_qc_ann"] = bla[
                    "y_nonpar_qc_ann"
                ] = "mean_fraction_over_dp_0"

            if args.by_fraction_samples_over_cov:
                annotation_format = "fraction_over_{cov}x"
                bla["autosome_par_cutoff"] = bla["x_nonpar_cutoff"] = bla[
                    "y_nonpar_cutoff"
                ] = args.fraction_samples
                bla["autosome_par_qc_ann"] = bla[
                    "x_nonpar_qc_ann"
                ] = annotation_format.format(cov=args.autosome_par_xx_cov)
                bla["y_nonpar_qc_ann"] = annotation_format.format(
                    cov=args.xy_nonpar_cov
                )

            high_qual_cutoffs = get_high_qual_cutoff_dict(
                **bla,
                split_by_sex=True,
            )
            ht = get_interval_qc_pass(
                ht,
                high_qual_cutoffs,
                per_platform=args.per_platform,
                all_platforms=args.all_platforms,
                min_platform_size=args.min_platform_size,
            )

            ht.write(
                get_checkpoint_path("interval_qc_pass")
                if test
                else interval_qc_pass.path,
                overwrite=overwrite,
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
        help="Test using only 2 partitions on chr20, chrX, and chrY.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    sex_coverage_args = parser.add_argument_group(
        "Sex chromosome interval coverage",
        "Arguments used for computing interval coverage on sex chromosomes.",
    )
    sex_coverage_args.add_argument(
        "--sex-chr-interval-coverage",
        help=(
            "Create a MatrixTable of interval-by-sample coverage on sex chromosomes with intervals split at PAR "
            "regions."
        ),
        action="store_true",
    )
    sex_coverage_args.add_argument(
        "--calling-interval-name",
        help="Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'.",
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    sex_coverage_args.add_argument(
        "--calling-interval-padding",
        help="Number of base pair padding to use on the calling intervals. One of 0 or 50 bp.",
        type=int,
        choices=[0, 50],
        default=50,
    )

    interval_qc_args = parser.add_argument_group(
        "Compute aggregate interval stats for interval QC",
        "Arguments used for computing interval QC stats.",
    )
    interval_qc_args.add_argument(
        "--generate-interval-qc-ht",
        help="Compute aggregate interval stats for interval QC from coverage MatrixTables.",
        action="store_true",
    )
    interval_qc_args.add_argument(
        "--mean-dp-thresholds",
        help=(
            "List of mean DP cutoffs to determining the fraction of samples with mean coverage >= the cutoff for each "
            "interval."
        ),
        type=int,
        nargs="+",
        default=[5, 10, 15, 20, 25],
    )

    interval_qc_pass_args = parser.add_argument_group(
        "Generate interval QC pass annotation",
        "Arguments used for determining intervals that pass QC.",
    )
    interval_qc_pass_args.add_argument(
        "--generate-interval-qc-pass-ht",
        help=(
            "Create a MatrixTable of interval-by-sample coverage on sex chromosomes with intervals split at PAR "
            "regions."
        ),
        action="store_true",
    )

    interval_qc_pass_platform_opt_parser = (
        interval_qc_pass_args.add_mutually_exclusive_group(required=False)
    )
    interval_qc_pass_platform_opt_parser.add_argument(
        "--per-platform",
        help=(
            "Whether to make the interval QC pass annotation a DictionaryExpression with interval QC pass per "
            "platform."
        ),
        action="store_true",
    )
    interval_qc_pass_platform_opt_parser.add_argument(
        "--all-platforms",
        help=(
            "Whether to consider an interval as passing QC only if it passes interval QC per platform across all "
            "platforms (with a sample size above '--min-platform-size')."
        ),
        action="store_true",
    )
    interval_qc_pass_args.add_argument(
        "--min-platform-size",
        help=(
            "Required size of a platform to be considered in '--all-platforms'. Only platforms that "
            "have # of samples > 'min_platform_size' are used to determine intervals that have a high coverage across "
            "all platforms."
        ),
        type=int,
        default=100,
    )
    interval_qc_pass_method_parser = interval_qc_pass_args.add_mutually_exclusive_group(
        required=False
    )
    interval_qc_pass_method_parser.add_argument(
        "--by-mean-fraction-over-dp-0",
        help="Whether to use the mean fraction of bases over DP 0 to determine high quality intervals. Can't be set "
        "at the same time as '--by-prop-samples-over-cov'.",
        action="store_true",
    )
    interval_qc_pass_method_parser.add_argument(
        "--by-fraction-samples-over-cov",
        help=(
            "Whether to determine high quality intervals using the fraction of samples (--fraction-samples) with a "
            "mean interval coverage over a specified coverage for intervals on the the autosomes/sex chromosome "
            "PAR/chrX in XX individuals (--autosome-par-xx-cov) and intervals on non-PAR chrX in XY individuals and "
            "non-PAR chrY in XY individuals (--xy-nonpar-cov). Can't be set at the same time as "
            "'--by-mean-fraction-over-dp-0'"
        ),
        action="store_true",
    )
    interval_qc_pass_args.add_argument(
        "--mean-fraction-over-dp-0",
        help="Mean fraction of bases over DP 0 used to define high quality intervals.",
        type=float,
        default=0.99,
    )
    interval_qc_pass_args.add_argument(
        "--autosome-par-xx-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on the the autosomes, sex chromosome PAR, "
            "and chrX in XX individuals. This field must be in the interval coverage MatrixTables!"
        ),
        type=int,
        default=20,
    )
    interval_qc_pass_args.add_argument(
        "--xy-nonpar-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on non-PAR chrX in XY individuals and non-PAR "
            "chrY in XY individuals. This field must be in the interval coverage MatrixTables!"
        ),
        type=int,
        default=10,
    )
    interval_qc_pass_args.add_argument(
        "--fraction-samples",
        help=(
            "Fraction of samples with mean coverage greater than '--autosome-par-xx-cov'/'--xy-nonpar-cov' over the "
            "interval to determine high coverage intervals."
        ),
        type=float,
        default=0.85,
    )

    args = parser.parse_args()

    if args.generate_interval_qc_pass_ht and not (
        args.by_mean_fraction_over_dp_0 or args.by_fraction_samples_over_cov
    ):
        parser.error(
            "One of --by-mean-fraction-over-dp-0 or --by-fraction-samples-over-cov is required when "
            "--generate-interval-qc-pass-ht is specified."
        )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
