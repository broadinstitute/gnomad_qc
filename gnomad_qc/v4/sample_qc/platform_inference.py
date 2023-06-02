"""Script to assign platforms based on per interval fraction of bases over DP 0 PCA results using HDBSCAN."""
import argparse
import logging

import hail as hl
from gnomad.sample_qc.platform import assign_platform_from_pcs, run_platform_pca
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_logging_path,
    gnomad_v4_testset_meta,
)
from gnomad_qc.v4.resources.sample_qc import (
    hard_filtered_samples_no_sex,
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
    """Assign platforms based on PCA of per interval fraction of bases over DP 0."""
    hl.init(
        log="/platform_pca.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default.
    hl._set_flags(use_new_shuffle="1")

    calling_interval_name = args.calling_interval_name
    calling_interval_padding = args.calling_interval_padding
    hdbscan_min_cluster_size = args.hdbscan_min_cluster_size
    hdbscan_min_samples = args.hdbscan_min_samples
    n_assignment_pcs = args.n_assignment_pcs

    try:
        if args.run_platform_pca:
            logger.info("Running platform PCA...")
            if args.test:
                test_coverage_path = get_checkpoint_path(
                    f"test_interval_coverage.{calling_interval_name}.pad{calling_interval_padding}",
                    mt=True,
                )
                if file_exists(test_coverage_path):
                    mt = hl.read_matrix_table(test_coverage_path)
                else:
                    raise FileNotFoundError(
                        "The test interval coverage MatrixTable does not exist for"
                        f" calling interval {calling_interval_name} and interval"
                        f" padding {calling_interval_padding}. Please run"
                        " --compute_coverage with the --test argument and needed"
                        " --calling_interval_name/--calling_interval_padding"
                        " arguments."
                    )
            else:
                mt = interval_coverage.mt()

            logger.info(
                "Removing hard filtered samples from interval coverage MatrixTable..."
            )
            if args.test:
                ht = gnomad_v4_testset_meta.ht()
                ht = ht.filter(hl.len(ht.rand_sampling_meta.hard_filters_no_sex) != 0)
            else:
                ht = hard_filtered_samples_no_sex.ht()

            mt = mt.filter_cols(hl.is_missing(ht[mt.col_key]))

            logger.info("Filter interval coverage MatrixTable to autosomes...")
            mt = mt.filter_rows(mt.interval.start.in_autosome())

            logger.info(
                "Annotating interval coverage MatrixTable with 'callrate' defined as"
                " fraction over dp 0..."
            )
            # Default `hl.vds.interval_coverage` will return a list for
            # `fraction_over_dp_threshold` where the second element is dp >= 1 (dp > 0).
            mt = mt.annotate_entries(callrate=mt.fraction_over_dp_threshold[1])

            # NOTE: added None binarization_threshold parameter to be consistent with runs before this parameter existed. # noqa
            eigenvalues, scores_ht, loadings_ht = run_platform_pca(
                mt, binarization_threshold=None, n_pcs=args.n_platform_pcs
            )
            scores_ht = scores_ht.annotate_globals(**mt.index_globals())
            scores_ht.write(
                get_checkpoint_path(
                    f"test_platform_scores.{calling_interval_name}.pad{calling_interval_padding}"
                )
                if args.test
                else platform_pca_scores.path,
                overwrite=args.overwrite,
            )
            loadings_ht = loadings_ht.annotate_globals(**mt.index_globals())
            loadings_ht.write(
                get_checkpoint_path(
                    f"test_platform_loadings.{calling_interval_name}.pad{calling_interval_padding}"
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
                    f"test_platform_eigenvalues.{calling_interval_name}.pad{calling_interval_padding}"
                )
                if args.test
                else platform_pca_eigenvalues.path,
                overwrite=args.overwrite,
            )

        if args.assign_platforms:
            # TODO: Add drop of `gq_thresholds` for future versions.
            logger.info("Assigning platforms based on platform PCA clustering")
            if args.test:
                test_scores_path = get_checkpoint_path(
                    f"test_platform_scores.{calling_interval_name}.pad{calling_interval_padding}"
                )
                if file_exists(test_scores_path):
                    scores_ht = hl.read_table(test_scores_path)
                else:
                    raise FileNotFoundError(
                        "The test platform PCA Table does not exist for calling"
                        f" interval {calling_interval_name} and interval padding"
                        f" {calling_interval_padding}. Please run hard_filters.py"
                        " --compute-coverage and platform_inference.py"
                        " --run-platform-pca with the --test argument and needed"
                        " --calling-interval-name/--calling-interval-padding"
                        " arguments."
                    )
            else:
                scores_ht = platform_pca_scores.ht()

            platform_ht = assign_platform_from_pcs(
                scores_ht.annotate(scores=scores_ht.scores[:n_assignment_pcs]),
                hdbscan_min_cluster_size=hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )

            # Make sure hdbscan_min_samples is not None before annotating globals.
            if not hdbscan_min_samples:
                hdbscan_min_samples = hdbscan_min_cluster_size
            else:
                hdbscan_min_samples = hdbscan_min_samples
            platform_ht = platform_ht.annotate_globals(
                hdbscan_min_cluster_size=hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
                n_pcs=n_assignment_pcs,
                **scores_ht.index_globals(),
            )
            platform_ht = platform_ht.checkpoint(
                get_checkpoint_path(
                    f"test_platform_assignment.{calling_interval_name}.pad{calling_interval_padding}"
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
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb',"
            " 'broad', or 'intersection'."
        ),
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    parser.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or"
            " 50 bp."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )
    parser.add_argument(
        "--run-platform-pca",
        help=(
            "Runs platform PCA (assumes coverage MatrixTable was computed,"
            " --compute-coverage)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--n-platform-pcs",
        help="Number of platform PCs to compute.",
        type=int,
        default=30,
    )
    parser.add_argument(
        "--assign-platforms",
        help=(
            "Assigns platforms based on per interval fraction of bases over DP 0 PCA"
            " results using HDBSCAN."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--n-assignment-pcs",
        help="Number of platform PCs to use for platform assignment.",
        type=int,
        default=10,
    )
    parser.add_argument(
        "--hdbscan-min-samples",
        help=(
            "Minimum samples parameter for HDBSCAN. If not specified,"
            " --hdbscan-min-cluster-size is used."
        ),
        type=int,
        required=False,
    )
    parser.add_argument(
        "--hdbscan-min-cluster-size",
        help="Minimum cluster size parameter for HDBSCAN.",
        type=int,
        default=50,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
