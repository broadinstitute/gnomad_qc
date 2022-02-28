import argparse
import logging

import hail as hl

from gnomad.sample_qc.platform import (
    assign_platform_from_pcs,
    run_platform_pca,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    calling_intervals,
    get_checkpoint_path,
    get_gnomad_v4_vds,
    gnomad_v4_testset,
    gnomad_v4_testset_meta,
    get_logging_path,
)
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
            if args.test:
                logger.info("Loading test VDS...")
                vds = gnomad_v4_testset.vds()
            else:
                logger.info("Loading full v4 VDS...")
                vds = get_gnomad_v4_vds(remove_hard_filtered_samples=False)

            logger.info(
                "Loading calling intervals: %s with padding of %d...",
                args.calling_interval_name,
                args.calling_interval_padding,
            )
            ht = calling_intervals(
                args.calling_interval_name, args.calling_interval_padding
            ).ht()
            mt = hl.vds.interval_coverage(vds, intervals=ht)
            mt = mt.annotate_globals(
                calling_interval_name=args.calling_interval_name,
                calling_interval_padding=args.calling_interval_padding,
            )
            mt.write(
                get_checkpoint_path(
                    f"test_interval_coverage.{args.calling_interval_name}.pad{args.calling_interval_padding}",
                    mt=True,
                )
                if args.test
                else interval_coverage.path,
                overwrite=args.overwrite,
            )

        if args.run_platform_pca:
            logger.info("Running platform PCA...")
            if args.test:
                test_coverage_path = get_checkpoint_path(
                    f"test_interval_coverage.{args.calling_interval_name}.pad{args.calling_interval_padding}",
                    mt=True,
                )
                if file_exists(test_coverage_path):
                    mt = hl.read_matrix_table(test_coverage_path)
                else:
                    raise FileNotFoundError(
                        f"The test interval coverage MatrixTable does not exist for calling interval "
                        f"{args.calling_interval_name} and interval padding {args.calling_interval_padding}. "
                        f"Please run --compute_coverage with the --test argument and needed "
                        f"--calling_interval_name/--calling_interval_padding arguments."
                    )
            else:
                mt = interval_coverage.mt()

            logger.info(
                "Removing hard filtered samples from interval coverage MatrixTable..."
            )
            if args.test:
                ht = gnomad_v4_testset_meta.ht()
                ht = ht.filter(hl.len(ht.rand_sampling_meta.hard_filters_no_sex) == 0)
            else:
                ht = hard_filtered_samples.ht()

            mt = mt.filter_cols(hl.is_defined(ht[mt.col_key]))

            logger.info("Filter interval coverage MatrixTable to autosomes...")
            mt = mt.filter_rows(mt.interval.start.in_autosome())

            logger.info(
                "Annotating interval coverage MatrixTable with 'callrate' defined as fraction over dp 0..."
            )
            # Default `hl.vds.interval_coverage` will return a list for `fraction_over_dp_threshold` where the
            # second element is dp >= 1 (dp > 0)
            mt = mt.annotate_entries(callrate=mt.fraction_over_dp_threshold[1])

            # NOTE: added None binarization_threshold parameter to be consistent with runs before this parameter existed
            eigenvalues, scores_ht, loadings_ht = run_platform_pca(
                mt, binarization_threshold=None
            )
            scores_ht = scores_ht.annotate_globals(**mt.index_globals())
            scores_ht.write(
                get_checkpoint_path(
                    f"test_platform_scores.{args.calling_interval_name}.pad{args.calling_interval_padding}"
                )
                if args.test
                else platform_pca_scores.path,
                overwrite=args.overwrite,
            )
            loadings_ht = loadings_ht.annotate_globals(**mt.index_globals())
            loadings_ht.write(
                get_checkpoint_path(
                    f"test_platform_loadings.{args.calling_interval_name}.pad{args.calling_interval_padding}"
                )
                if args.test
                else platform_pca_loadings.path,
                overwrite=args.overwrite,
            )
            eigenvalues_ht = hl.Table.parallelize(
                hl.literal(
                    [{"PC": i + 1, "eigenvalue": x} for i, x in enumerate(eigenvalues)],
                    "array<struct{PC: int, eigenvalue: float}>",
                )
            )
            eigenvalues_ht = eigenvalues_ht.annotate_globals(**mt.index_globals())
            eigenvalues_ht.write(
                get_checkpoint_path(
                    f"test_platform_eigenvalues.{args.calling_interval_name}.pad{args.calling_interval_padding}"
                )
                if args.test
                else platform_pca_eigenvalues.path,
                overwrite=args.overwrite,
            )

        if args.assign_platforms:
            logger.info("Assigning platforms based on platform PCA clustering")
            if args.test:
                test_scores_path = get_checkpoint_path(
                    f"test_platform_scores.{args.calling_interval_name}.pad{args.calling_interval_padding}"
                )
                if file_exists(test_scores_path):
                    scores_ht = hl.read_table(test_scores_path)
                else:
                    raise FileNotFoundError(
                        f"The test platform PCA Table does not exist for calling interval "
                        f"{args.calling_interval_name} and interval padding {args.calling_interval_padding}. "
                        f"Please run --compute_coverage and -- run_platform_pca with the --test argument and needed "
                        f"--calling_interval_name/--calling_interval_padding arguments."
                    )
            else:
                scores_ht = hl.read_table(platform_pca_scores.ht())

            platform_ht = assign_platform_from_pcs(
                scores_ht,
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=args.hdbscan_min_samples,
            )

            # Make sure hdbscan_min_samples is not None before annotating globals
            if not args.hdbscan_min_samples:
                hdbscan_min_samples = args.hdbscan_min_cluster_size
            else:
                hdbscan_min_samples = args.hdbscan_min_samples
            platform_ht = platform_ht.annotate_globals(
                hdbscan_min_cluster_size=args.hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
                **scores_ht.index_globals(),
            )
            platform_ht = platform_ht.checkpoint(
                get_checkpoint_path(
                    f"test_platform_assignment.{args.calling_interval_name}.pad{args.calling_interval_padding}"
                )
                if args.test
                else platform.path,
                overwrite=args.overwrite,
            )
            logger.info(f"Platform PCA Table count: {platform_ht.count()}")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("platform_pca"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite output files (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Use the v4 test dataset instead of the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--compute_coverage",
        help="Compute per interval coverage metrics using Hail's vds.interval_coverage method.",
        action="store_true",
    )
    parser.add_argument(
        "--calling_interval_name",
        help="Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'.",
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    parser.add_argument(
        "--calling_interval_padding",
        help="Number of base pair padding to use on the calling intervals. One of 0 or 50 bp.",
        type=int,
        choices=[0, 50],
        default=50,
    )
    parser.add_argument(
        "--run_platform_pca",
        help="Runs platform PCA (assumes coverage MatrixTable was computed, --compute_coverage).",
        action="store_true",
    )
    parser.add_argument(
        "--assign_platforms",
        help="Assigns platforms based on per interval fraction of bases over DP 0 PCA results using HDBSCAN.",
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
        default=50,
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
