"""Script to assign platforms based on per interval fraction of bases over DP 0 PCA results using HDBSCAN."""
import argparse
import logging

import hail as hl
from gnomad.sample_qc.platform import assign_platform_from_pcs, run_platform_pca
from gnomad.utils.slack import slack_notifications

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_logging_path, gnomad_v4_testset_meta
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


def get_pipeline_resources(
    test: bool,
    overwrite: bool,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the platform inference pipeline.

    :param test: Whether to gather all resources for the test dataset.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        platform inference pipeline.
    """
    # Initialize platform inference pipeline resource collection.
    platform_pipeline = PipelineResourceCollection(
        pipeline_name="platform_inference",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the platform inference pipeline.
    platform_pca = PipelineStepResourceCollection(
        "--run-platform-pca",
        input_resources={
            "hard_filters.py --compute-coverage": {
                "coverage_mt": interval_coverage(test=test)
            },
            "hard_filters.py --compute-hard-filters": {
                "hard_filter_ht": hard_filtered_samples(include_sex_filter=False)
            },
        },
        output_resources={
            "platform_pca_scores_ht": platform_pca_scores(test=test),
            "platform_pca_loadings_ht": platform_pca_loadings(test=test),
            "platform_pca_eigenvalues_ht": platform_pca_eigenvalues(test=test),
        },
    )
    assign_platforms = PipelineStepResourceCollection(
        "--assign-platforms",
        input_resources={
            "--run-platform-pca": {
                "platform_pca_scores_ht": platform_pca_scores(test=test)
            }
        },
        output_resources={"platform_ht": platform(test=test)},
    )

    # Add all steps to the platform inference pipeline resource collection.
    platform_pipeline.add_steps(
        {
            "run_platform_pca": platform_pca,
            "assign_platforms": assign_platforms,
        }
    )

    return platform_pipeline


def main(args):
    """Assign platforms based on PCA of per interval fraction of bases over DP 0."""
    hl.init(
        log="/platform_pca.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default.
    hl._set_flags(use_new_shuffle="1")

    test = args.test
    overwrite = args.overwrite
    hdbscan_min_cluster_size = args.hdbscan_min_cluster_size
    hdbscan_min_samples = args.hdbscan_min_samples
    n_assignment_pcs = args.n_assignment_pcs

    platform_resources = get_pipeline_resources(test=test, overwrite=overwrite)

    try:
        if args.run_platform_pca:
            logger.info("Running platform PCA...")
            res = platform_resources.run_platform_pca
            res.check_resource_existence()
            mt = res.coverage_mt.mt()

            logger.info(
                "Removing hard filtered samples from interval coverage MatrixTable..."
            )
            if args.test:
                ht = gnomad_v4_testset_meta.ht()
                ht = ht.filter(hl.len(ht.rand_sampling_meta.hard_filters_no_sex) != 0)
            else:
                ht = res.hard_filter_ht.ht()

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
            scores_ht.write(res.platform_pca_scores_ht.path, overwrite=overwrite)
            loadings_ht = loadings_ht.annotate_globals(**mt.index_globals())
            loadings_ht.write(res.platform_pca_loadings_ht.path, overwrite=overwrite)
            eigenvalues_ht = hl.Table.parallelize(
                hl.literal(
                    [{"PC": i + 1, "eigenvalue": x} for i, x in enumerate(eigenvalues)],
                    "array<struct{PC: int, eigenvalue: float}>",
                )
            )
            eigenvalues_ht = eigenvalues_ht.annotate_globals(**mt.index_globals())
            eigenvalues_ht.write(
                res.platform_pca_eigenvalues_ht.path, overwrite=overwrite
            )

        if args.assign_platforms:
            # TODO: Add drop of `gq_thresholds` for future versions.
            logger.info("Assigning platforms based on platform PCA clustering")
            res = platform_resources.assign_platforms
            res.check_resource_existence()
            scores_ht = res.platform_pca_scores_ht.ht()

            ht = assign_platform_from_pcs(
                scores_ht.annotate(scores=scores_ht.scores[:n_assignment_pcs]),
                hdbscan_min_cluster_size=hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
            )

            # Make sure hdbscan_min_samples is not None before annotating globals.
            if not hdbscan_min_samples:
                hdbscan_min_samples = hdbscan_min_cluster_size
            else:
                hdbscan_min_samples = hdbscan_min_samples
            ht = ht.annotate_globals(
                hdbscan_min_cluster_size=hdbscan_min_cluster_size,
                hdbscan_min_samples=hdbscan_min_samples,
                n_pcs=n_assignment_pcs,
                **scores_ht.index_globals(),
            )
            ht = ht.checkpoint(res.platform_ht.path, overwrite=overwrite)
            logger.info(f"Platform PCA Table count: {ht.count()}")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("platform_pca"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
