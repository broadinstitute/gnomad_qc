import argparse
import logging

import hail as hl

from gnomad.sample_qc.platform import (
    assign_platform_from_pcs,
    run_platform_pca,
)
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v4.resources.basics import calling_intervals, get_gnomad_v4_vds, logging_path
from gnomad_qc.v4.resources.sample_qc import (
    hard_filtered_samples,
    interval_coverage,
    platform,
    platform_pca_eigenvalues,
    platform_pca_loadings,
    platform_pca_scores,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("platform_pca")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/platform_pca.log", default_reference="GRCh38")

    try:
        if args.compute_coverage:
            vds = get_gnomad_v4_vds(remove_hard_filtered_samples=False)
            ht = calling_intervals(
                args.calling_interval_name, args.calling_interval_padding
            ).ht()
            mt = hl.vds.interval_coverage(vds, intervals=ht)
            mt.write(interval_coverage.path)

        if args.run_platform_pca:
            logger.info("Running platform PCA...")
            mt = interval_coverage.mt()

            logger.info("Removing hard filtered samples from interval coverage MT...")
            ht = hard_filtered_samples.ht()
            mt = mt.filter_cols(hl.is_defined(ht[mt.col_key]))

            logger.info("Filter interval coverage MT to autosomes...")
            mt = mt.filter_rows(mt.interval.start.in_autosome())

            logger.info(
                "Annotating interval coverage MT with 'callrate' defined as fraction over dp 0..."
            )
            # Default `hl.vds.interval_coverage` will return a list for `fraction_over_dp_threshold` where the second element is dp >= 1 (dp > 0)
            mt = mt.annotate_entries(callrate=mt.fraction_over_dp_threshold[1])

            # NOTE: added None binarization_threshold parameter to make sure we things the same way as before parameter existed
            eigenvalues, scores_ht, loadings_ht = run_platform_pca(
                mt, binarization_threshold=None
            )
            scores_ht.write(platform_pca_scores.path, overwrite=args.overwrite)
            loadings_ht.write(platform_pca_loadings.path, overwrite=args.overwrite)

        if args.assign_platforms:
            logger.info("Assigning platforms based on platform PCA clustering")
            scores_ht = hl.read_table(platform_pca_scores.ht())
            platform_ht = assign_platform_from_pcs(
                scores_ht,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=args.hdbscan_min_samples,
            )

            # Make sure hdbscan_min_samples is not None before annotating globals
            if not args.hdbscan_min_samples:
                hdbscan_min_samples = args.hdbscan_min_cluster_size
            platform_ht = platform_ht.annotate_globals(
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )
            platform_ht = platform_ht.repartition(args.n_partitions)
            platform_ht = platform_ht.checkpoint(
                platform.path, overwrite=args.overwrite
            )
            logger.info(f"Platform PCA Table count: {platform_ht.count()}")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(logging_path())


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--compute_coverage", help="", action="store_true",
    )
    parser.add_argument(
        "--calling_interval_name",
        help="",
        type=str,
        choices=["ukb", "broad", "intersection"],
    )
    parser.add_argument(
        "--calling_interval_padding", help="", type=int, choices=[0, 50],
    )
    parser.add_argument(
        "--run_platform_pca",
        help="Runs platform PCA (assumes callrate MT was computed)",
        action="store_true",
    )
    parser.add_argument(
        "--assign_platforms",
        help="Assigns platforms based on callrate PCA results using HDBSCAN",
        action="store_true",
    )
    parser.add_argument(
        "--hdbscan_min_samples",
        help="Minimum samples parameter for HDBSCAN. If not specified, --hdbscan_min_cluster_size is used.",
        type=int,
        required=False,
    )
    parser.add_argument(
        "--hdbscan_min_cluster_size",
        help="Minimum cluster size parameter for HDBSCAN.",
        type=int,
        default=100,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
