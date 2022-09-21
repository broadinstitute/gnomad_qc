import argparse
import logging
from typing import List, Tuple, Union

from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.slack import slack_notifications
import hail as hl

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_checkpoint_path, get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    hard_filtered_samples,
    hard_filtered_samples_no_sex,
    interval_coverage,
    interval_qc,
    platform,
    sex,
    sex_imputation_coverage,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("interval_qc")
logger.setLevel(logging.INFO)


def filter_to_test(
    mt: hl.MatrixTable,
    sex_mt: hl.MatrixTable,
    num_partitions: int = 10,
) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    """
    Filter `mt` to `num_partitions` partitions on chr1 and `sex_mt` to `num_partitions` partitions on chrX and chrY.

    :param mt: Input MatrixTable to filter to specified number of partitions on chr1.
    :param sex_mt: Input MatrixTable to filter to specified number of partitions on chrX and all of chrY.
    :param num_partitions: Number of partitions to grab from mt.
    :return: Input MatrixTables filtered to `num_partitions` on chr1, chrX, and all of chrY.
    """
    logger.info(
        "Filtering to %d partitions on chr1, chrX, and all of chrY for testing...",
        num_partitions,
    )
    mt = mt._filter_partitions(range(num_partitions))
    sex_mt_chrx = sex_mt._filter_partitions(range(num_partitions))
    sex_mt_chry = sex_mt.filter_rows(
        (sex_mt.interval.start.contig == "chrY")
    ).repartition(100)
    sex_mt_chry = sex_mt_chry._filter_partitions(range(num_partitions))

    return mt, sex_mt_chrx.union_rows(sex_mt_chry)


def compute_interval_qc(
    mt: hl.MatrixTable,
    add_per_platform: bool = False,
    platform_ht: hl.Table = None,
    mean_dp_thresholds: List[int] = [5, 10, 15, 20, 25],
    split_by_sex: bool = False,
) -> hl.Table:
    """


    :param mt: Input interval coverage MatrixTable
    :param add_per_platform: Whether to compute the aggregate interval means per platform.
    :param platform_ht: Input platform assignment Table.
    :param mean_dp_thresholds: List of mean DP thresholds to use for computing the fraction of samples with mean
        interval DP >= the threshold.
    :param bool split_by_sex: Whether the interval QC should be stratified by sex. If True, mt must be annotated with sex_karyotype.
    :return: Table with interval QC annotations
    """
    if add_per_platform and platform_ht is None:
        raise ValueError("'platform_ht' must be defined if 'add_per_platform' is True!")

    def _get_agg_expr(expr, agg_func=hl.agg.mean, group_by=None):
        """

        :param expr:
        :param agg_func:
        :return:
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

    mt = mt.select_globals(mean_dp_thresholds=mean_dp_thresholds)
    if add_per_platform:
        mt = mt.annotate_cols(platform=platform_ht[mt.col_key].qc_platform)
        platforms = platform_ht.aggregate(
            hl.agg.collect_as_set(platform_ht.qc_platform)
        )
        mt = mt.annotate_globals(platforms=platforms)

    agg_expr = {
        "interval_mean_dp": _get_agg_expr(mt.mean_dp),
        "prop_samples_by_dp": hl.struct(
            **{
                f"over_{dp}x": _get_agg_expr(mt.mean_dp >= dp, hl.agg.fraction)
                for dp in mean_dp_thresholds
            }
        ),
        "mean_fraction_over_dp_0": _get_agg_expr(mt.fraction_over_dp_threshold[1]),
    }

    add_globals = hl.struct()
    if add_per_platform:
        logger.info("Adding per platform aggregation...")
        agg_expr.update(
            {
                "platform_interval_mean_dp": _get_agg_expr(
                    mt.mean_dp,
                    group_by=mt.platform,
                ),
                "platform_prop_samples_by_dp": hl.struct(
                    **{
                        f"over_{dp}x": _get_agg_expr(
                            mt.mean_dp >= dp,
                            agg_func=hl.agg.fraction,
                            group_by=mt.platform,
                        )
                        for dp in mean_dp_thresholds
                    }
                ),
                "platform_mean_fraction_over_dp_0": _get_agg_expr(
                    mt.fraction_over_dp_threshold[1],
                    group_by=mt.platform,
                ),
            }
        )
        add_globals = hl.struct(
            platform_n_samples=mt.aggregate_cols(
                hl.agg.group_by(mt.platform, hl.agg.count())
            )
        )

    ht = mt.select_rows(**agg_expr).rows()
    ht = ht.annotate_globals(**add_globals)

    return ht


def get_interval_qc_pass(
    interval_qc_ht: hl.Table,
    per_platform: bool = False,
    all_platforms: bool = False,
    min_platform_size: int = 100,
    by_mean_fraction_over_dp_0: bool = True,
    by_prop_samples_over_cov: bool = False,
    mean_fraction_over_dp_0: float = 0.99,
    autosome_par_xx_cov: int = 20,
    xy_nonpar_cov: int = 10,
    prop_samples: float = 0.85,
) -> hl.Table:
    """
    Add `interval_qc_pass` annotation to indicate whether the site falls within a high coverage interval.

    :param t: Input MatrixTable or Table.
    :param per_platform: Whether filter to per platform high coverage intervals for the sex ploidy imputation.
    :param all_platforms: Whether to filter to high coverage intervals for the sex ploidy imputation. Use only intervals that are considered high coverage across all platforms.
    :param by_mean_fraction_over_dp_0: Whether to use the mean fraction of bases over DP 0 to determine high coverage intervals.
    :param by_prop_samples_over_cov: Whether to determine high coverage intervals using the proportion of samples with a mean interval coverage over a specified coverage for chrX (--x-cov), chrY (--y-cov), and the normalization contig (--norm-cov).
    :param mean_fraction_over_dp_0: Mean fraction of bases over DP used to define high coverage intervals. Default is 0.99.
    :param autosome_par_xx_cov: Mean coverage level used to define high coverage intervals on the the autosomes/sex chr par/female X. Default is 20.
    :param xy_nonpar_cov: Mean coverage level used to define high coverage intervals on male X and Y. This field must be in the sex interval coverage MT. Default is 10.
    :param prop_samples: Proportion of samples with mean coverage greater than `autosome_cov`/`sex_cov` over the interval to determine high coverage intervals. Default is 0.85.
    :return: MatrixTable or Table with samples removed
    """
    if (by_mean_fraction_over_dp_0 and by_prop_samples_over_cov) or (
        not by_mean_fraction_over_dp_0 and not by_prop_samples_over_cov
    ):
        raise ValueError(
            "One and only one of 'high_cov_by_mean_fraction_over_dp_0' and 'high_cov_by_prop_samples_over_cov' must be "
            "True!"
        )
    if per_platform and all_platforms:
        raise ValueError("Only one of 'per_platform' and 'all_platforms' can be True!")

    interval_qc_ht = interval_qc_ht.annotate_globals(
        per_platform=per_platform,
        all_platforms=all_platforms,
    )
    interval_start = interval_qc_ht.interval.start
    autosome_or_par = interval_start.in_autosome_or_par()
    x_non_par = interval_start.in_x_nonpar()
    y_non_par = interval_start.in_y_nonpar()

    if per_platform or all_platforms:
        platform_n_samples = (
            interval_qc_ht.index_globals().platform_n_samples.collect()[0]
        )
        ann_prefix = "platform_"
    else:
        ann_prefix = ""

    if by_mean_fraction_over_dp_0:
        add_globals = hl.struct(mean_fraction_over_dp_0=mean_fraction_over_dp_0)
        qc_expr = interval_qc_ht[f"{ann_prefix}mean_fraction_over_dp_0"]
        qc_autosome_par_expr = qc_expr["all"]
        qc_xx_expr = qc_expr.get("XX", None)
        qc_xy_expr = qc_expr.get("XY", None)
        cutoff = mean_fraction_over_dp_0
    if by_prop_samples_over_cov:
        add_globals = hl.struct(
            autosome_par_xx_cov=autosome_par_xx_cov,
            xy_nonpar_cov=xy_nonpar_cov,
            prop_samples=prop_samples,
        )
        qc_expr = interval_qc_ht[f"{ann_prefix}prop_samples_by_dp"]
        qc_autosome_par_expr = qc_expr[f"over_{autosome_par_xx_cov}x"]["all"]
        qc_xx_expr = qc_expr[f"over_{autosome_par_xx_cov}x"].get("XX", None)
        qc_xy_expr = qc_expr[f"over_{xy_nonpar_cov}x"].get("XY", None)
        cutoff = prop_samples

    def _get_pass_expr(qc_autosome_par_expr, qc_xx_expr, qc_xy_expr):
        return (
            (autosome_or_par & (qc_autosome_par_expr > cutoff))
            | (x_non_par & (qc_xx_expr > cutoff) & (qc_xy_expr > cutoff))
            | (y_non_par & (qc_xy_expr > cutoff))
        )

    if per_platform or all_platforms:
        interval_qc_ht = interval_qc_ht.select(
            pass_interval_qc=hl.struct(
                **{
                    platform: _get_pass_expr(
                        qc_autosome_par_expr[platform],
                        qc_xx_expr[platform],
                        qc_xy_expr[platform],
                    )
                    for platform in platform_n_samples
                }
            )
        )
        if all_platforms:
            add_globals = add_globals.annotate(min_platform_size=min_platform_size)
            platforms = [
                platform
                for platform, n_samples in platform_n_samples.items()
                if (n_samples >= min_platform_size) & (platform != "platform_-1")
            ]
            interval_qc_ht = interval_qc_ht.select(
                pass_interval_qc=hl.all(
                    [
                        interval_qc_ht.pass_interval_qc[platform]
                        for platform in platforms
                    ]
                )
            )
    else:
        interval_qc_ht = interval_qc_ht.select(
            pass_interval_qc=_get_pass_expr(
                qc_autosome_par_expr,
                qc_xx_expr,
                qc_xy_expr,
            )
        )
    interval_qc_ht = interval_qc_ht.annotate_globals(**add_globals)

    return interval_qc_ht


def annotate_interval_qc_filter(
    t: Union[hl.MatrixTable, hl.Table],
    interval_qc_ht,
    **kwargs,
) -> Union[hl.MatrixTable, hl.Table]:
    # interval_qc_ht = interval_qc.ht()
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

    try:
        coverage_mt = interval_coverage.mt()
        sex_coverage_mt = sex_imputation_coverage.mt()
        if args.test:
            coverage_mt, sex_coverage_mt = filter_to_test(coverage_mt, sex_coverage_mt)

        if args.add_per_platform:
            platform_ht = platform.ht()
        else:
            platform_ht = None

        if args.remove_hardfiltered_no_sex:
            logger.info(
                "Removing hard-filtered (without sex imputation) samples from the coverage MTs..."
            )
            coverage_mt = coverage_mt.filter_cols(
                hl.is_missing(hard_filtered_samples_no_sex.ht()[coverage_mt.col_key])
            )
            sex_coverage_mt = sex_coverage_mt.filter_cols(
                hl.is_missing(
                    hard_filtered_samples_no_sex.ht()[sex_coverage_mt.col_key]
                )
            )
        if args.remove_hardfiltered:
            logger.info("Removing hard-filtered samples from the coverage MTs...")
            coverage_mt = coverage_mt.filter_cols(
                hl.is_missing(hard_filtered_samples.ht()[coverage_mt.col_key])
            )
            sex_coverage_mt = sex_coverage_mt.filter_cols(
                hl.is_missing(hard_filtered_samples.ht()[sex_coverage_mt.col_key])
            )

        ht = compute_interval_qc(
            filter_to_autosomes(coverage_mt),
            add_per_platform=args.add_per_platform,
            platform_ht=platform_ht,
            mean_dp_thresholds=args.mean_dp_thresholds,
        )

        logger.info("Filtering sex imputation coverage MT to sex chromosomes...")
        sex_coverage_mt = hl.filter_intervals(
            sex_coverage_mt,
            [hl.parse_locus_interval(contig) for contig in ["chrX", "chrY"]],
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

        ht = ht.union(
            compute_interval_qc(
                sex_coverage_mt,
                add_per_platform=args.add_per_platform,
                platform_ht=platform_ht,
                mean_dp_thresholds=args.mean_dp_thresholds,
                split_by_sex=True,
            )
        )

        ht.write(
            get_checkpoint_path("interval_qc") if args.test else interval_qc.path,
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
        help="Test using only 2 partitions on chr20, chrX, and chrY.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--add-per-platform",
        help=(
            "Whether to add a per platform interval QC annotation indicating if the interval is considered high "
            "coverage within the samples assigned platform."
        ),
        action="store_true",
    )
    hard_filter_parser = parser.add_mutually_exclusive_group(required=True)
    hard_filter_parser.add_argument(
        "--remove-hardfiltered-no-sex",
        help="",
        action="store_true",
    )
    hard_filter_parser.add_argument(
        "--remove-hardfiltered",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--mean-dp-thresholds",
        help=(
            "List of mean DP cutoffs to determining the fraction of samples with mean coverage >= the cutoff for each "
            "interval."
        ),
        type=int,
        nargs="+",
        default=[5, 10, 15, 20, 25],
    )

    args = parser.parse_args()
    main(args)

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
