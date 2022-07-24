import argparse
import logging

import hail as hl

from gnomad.utils.file_utils import file_exists
from gnomad.utils.generic import filter_to_autosomes
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_checkpoint_path, get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    interval_coverage,
    interval_qc,
    platform,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


def load_platform_ht(
    test: bool = False,
    calling_interval_name: str = "intersection",
    calling_interval_padding: int = 50,
):
    """

    :param test:
    :param calling_interval_name:
    :param calling_interval_padding:
    :return:
    """
    logger.info("Loading platform information...")
    if file_exists(platform.path):
        ht = platform.ht()
    elif test:
        test_platform_path = get_checkpoint_path(
            f"test_platform_assignment.{calling_interval_name}.pad{calling_interval_padding}"
        )
        if file_exists(test_platform_path):
            ht = hl.read_table(test_platform_path)
        else:
            raise FileNotFoundError(
                f"There is no final platform assignment Table written and a test platform assignment Table "
                f"does not exist for calling interval {calling_interval_name} and interval padding "
                f"{calling_interval_padding}. Please run platform_inference.py --assign_platforms "
                f"with the --test argument and needed --calling_interval_name/--calling_interval_padding "
                f"arguments."
            )
    else:
        raise FileNotFoundError(
            f"There is no final platform assignment Table written. Please run: "
            f"platform_inference.py --assign_platforms to compute the platform assignment Table."
        )

    return ht


def load_coverage_mt(
    calling_interval_name: str = "intersection", calling_interval_padding: int = 50
):
    """

    :param calling_interval_name:
    :param calling_interval_padding:
    :return:
    """
    logger.info("Loading interval coverage MatrixTable...")
    if file_exists(interval_coverage.path):
        mt = interval_coverage.mt()
    elif args.test:
        test_coverage_path = get_checkpoint_path(
            f"test_interval_coverage.{calling_interval_name}.pad{calling_interval_padding}",
            mt=True,
        )
        if file_exists(test_coverage_path):
            mt = hl.read_matrix_table(test_coverage_path)
        else:
            raise FileNotFoundError(
                f"There is no final coverage MatrixTable written and a test interval coverage MatrixTable does "
                f"not exist for calling interval {calling_interval_name} and interval padding "
                f"{calling_interval_padding}. Please run platform_inference.py --compute_coverage with the "
                f"--test argument and needed --calling_interval_name/--calling_interval_padding arguments."
            )
    else:
        raise FileNotFoundError(
            f"There is no final coverage MatrixTable written. Please run: "
            f"platform_inference.py --compute_coverage to compute the interval coverage MatrixTable."
        )

    return mt


def compute_interval_qc(
    mt: hl.MatrixTable,
    add_per_platform: bool = False,
    platform_ht: hl.Table = None,
    mean_dp_thresholds: tuple = (5, 10, 15, 20, 25, 30),
    split_by_sex: bool = False,
) -> hl.Table:
    """


    :param mt: Input interval coverage MatrixTable
    :param add_per_platform: Whether to run the sex ploidy and karyotype inference per platform
        inference. Using only intervals that are considered high coverage across all platforms
    :param platform_ht: Input platform assignment Table. This is only needed if per_platform or high_cov_by_platform_all are True
    :param mean_dp_thresholds: Mean coverage level used to define high coverage intervals on the autosomes. This field must
        be in the interval_coverage MT
    :param bool split_by_sex: Whether the interval QC should be stratified by sex. if True, mt must be annotated with sex_karyotype.
    :return: Table with interval QC annotations
    """

    def _get_frac_agg_expr(dp: int):
        """

        :param dp:
        :return:
        """
        frac_agg_expr = hl.agg.fraction(mt.mean_dp >= dp)
        if split_by_sex:
            frac_agg_expr = hl.agg.group_by(mt.sex_karyotype, frac_agg_expr)

        return frac_agg_expr

    agg_expr = {
        "interval_mean_dp": hl.agg.mean(mt.mean_dp),
        "pct_samples_by_dp": hl.struct(
            **{
                f"over_{dp}x": hl.agg.fraction(_get_frac_agg_expr(dp))
                for dp in mean_dp_thresholds
            }
        ),
    }

    if add_per_platform:
        logger.info("Adding per platform aggregation...")
        mt = mt.annotate_cols(platform=platform_ht[mt.col_key].qc_platform)
        agg_expr["pct_samples_by_platform"] = hl.agg.group_by(
            mt.platform,
            hl.struct(
                **{f"over_{dp}x": _get_frac_agg_expr(dp) for dp in mean_dp_thresholds}
            ),
        )

    ht = mt.annotate_rows(**agg_expr)
    ht = ht.annotate_globals(mean_dp_thresholds=mean_dp_thresholds)

    return ht


def main(args):
    hl.init(log="/interval_qc.log", default_reference="GRCh38")

    try:
        coverage_mt = load_coverage_mt(
            args.calling_interval_name, args.calling_interval_padding
        )
        if args.add_per_platform:
            platform_ht = load_platform_ht(
                args.test, args.calling_interval_name, args.calling_interval_padding
            )
        else:
            platform_ht = None

        ht = compute_interval_qc(
            filter_to_autosomes(coverage_mt),
            add_per_platform=args.add_per_platform,
            platform_ht=platform_ht,
            mean_dp_thresholds=args.mean_dp_thresholds,
        )

        logger.info("Filtering to sex chromosomes...")
        coverage_mt = hl.filter_intervals(
            coverage_mt,
            [hl.parse_locus_interval(contig) for contig in ["chrX", "chrY"]],
        )

        logger.info("Filtering to XX and XY samples...")
        sex_ht = sex.ht().select("sex_karyotype")
        coverage_mt = coverage_mt.annotate_cols(**sex_ht[coverage_mt.col_key])
        coverage_mt = coverage_mt.filter_cols(
            (coverage_mt.sex_karyotype == "XX") | (coverage_mt.sex_karyotype == "XY")
        )

        ht = ht.union(
            compute_interval_qc(
                coverage_mt,
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
        help="Test the pipeline using the gnomAD v4 test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--add-per-platform",
        help=(
            "Whether to add a per platform interval QC annotation indicating if the interval is considered high "
            "coverage within the samples assigned platform."
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
        "--autosome-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on the autosomes. This field must "
            "be in the interval_coverage MT!"
        ),
        type=int,
        default=20,
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
