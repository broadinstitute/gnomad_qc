"""Script to compute relatedness estimates among pairs of samples in the callset."""
import argparse
import logging
import textwrap

import hail as hl

from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.utils.file_utils import check_file_exists_raise_error
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_checkpoint_path, get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    cuking_input_path,
    cuking_output_path,
    joint_qc_meta,
    pca_related_samples_to_drop,
    pca_samples_rankings,
    joint_qc,
    relatedness,
)

from gnomad_qc.v4.sample_qc.cuKING.mt_to_cuking_inputs import mt_to_cuking_inputs
from gnomad_qc.v4.sample_qc.cuKING.cuking_outputs_to_ht import cuking_outputs_to_ht

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def main(args):
    """Compute relatedness estimates among pairs of samples in the callset."""
    test = args.test
    overwrite = args.overwrite

    if args.print_cuking_command:
        if (
            args.prepare_cuking_inputs
            or args.create_relatedness_table
            or args.compute_related_samples_to_drop
        ):
            raise ValueError(
                "--print-cuking-command can't be used simultaneously with other run "
                "modes."
            )

        logger.warning("Make sure that $PROJECT_ID is set!")
        logger.warning(
            "This printed command assumes that the cuKING directory is in the same "
            "location where the command is being run!"
        )
        print(
            textwrap.dedent(
                f"""\
                cd cuKING && \\
                ./cloud_batch_submit.py \\
                    --location=us-central1 \\
                    --project-id=$PROJECT_ID \\
                    --tag-name=$(git describe --tags) \\
                    --input-uri={cuking_input_path(test=test)} \\
                    --output-uri={cuking_output_path(test=test)} \\
                    --requester-pays-project=$PROJECT_ID \\
                    --kin-threshold={args.min_kin_cutoff} \\
                    --split-factor={args.cuking_split_factor} &&
                cd ..
                """
            )
        )
        return

    hl.init(log="/relatedness.log", default_reference="GRCh38")

    try:
        if args.prepare_inputs:
            parquet_uri = cuking_input_path(test=test)
            check_file_exists_raise_error(parquet_uri, not overwrite)
            mt_to_cuking_inputs(
                mt=joint_qc(test=test).mt(),
                parquet_uri=parquet_uri,
                overwrite=overwrite,
            )

        if args.create_cuking_relatedness_table:
            check_file_exists_raise_error(relatedness(test=test).path, not overwrite)
            ht = cuking_outputs_to_ht(parquet_uri=cuking_output_path(test=test))
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(relatedness(test=test).path, overwrite=overwrite)

        if args.compute_related_samples_to_drop:
            # compute_related_samples_to_drop uses a rank Table as a tie breaker when
            # pruning samples.
            check_file_exists_raise_error(
                [
                    pca_samples_rankings.path,
                    pca_related_samples_to_drop(test=test).path,
                ],
                not overwrite,
            )

            rank_ht = joint_qc_meta.ht()
            rank_ht = rank_ht.select(
                rank_ht.hard_filtered,
                rank_ht.releasable,
                rank_ht.chr20_mean_dp,
                present_in_v3=hl.is_defined(rank_ht.v3_meta),
            )

            # Favor v3 release samples, then v4 samples over v3 non-release samples.
            rank_ht = rank_ht.order_by(
                hl.asc(rank_ht.hard_filtered),
                hl.desc(rank_ht.present_in_v3 & rank_ht.releasable),
                hl.desc(rank_ht.releasable),
                hl.desc(rank_ht.chr20_mean_dp),
            ).add_index(name="rank")
            rank_ht = rank_ht.key_by(rank_ht.s)
            rank_ht = rank_ht.select(rank_ht.hard_filtered, rank_ht.rank)
            rank_ht = rank_ht.checkpoint(pca_samples_rankings.path, overwrite=overwrite)

            samples_to_drop = compute_related_samples_to_drop(
                relatedness(test=test).ht(),
                rank_ht,
                args.second_degree_kin_cutoff,
                filtered_samples=hl.literal(
                    rank_ht.aggregate(
                        hl.agg.filter(
                            rank_ht.hard_filtered, hl.agg.collect_as_set(rank_ht.s)
                        )
                    )
                ),
            )
            samples_to_drop = samples_to_drop.key_by(samples_to_drop.s)
            samples_to_drop.write(
                pca_related_samples_to_drop(test=test).path, overwrite=overwrite
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("relatedness"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test", help="Use a test MatrixTableResource as input.", action="store_true"
    )
    parser.add_argument(
        "--prepare-cuking-inputs",
        help=(
            "Converts the dense QC MatrixTable to a Parquet format suitable for cuKING."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--print-cuking-command",
        help="Print the command to submit a Cloud Batch job for running cuKING.",
        action="store_true",
    )
    parser.add_argument(
        "--cuking-split-factor",
        help=(
            "Split factor to use for splitting the full relatedness matrix table "
            "into equally sized submatrices that are computed independently to "
            "parallelize the relatedness computation and decrease the memory "
            "requirements. For example, to halve memory requirements, the full matrix "
            "can be split into  equally sized submatrices (i.e. a 'split factor' of 4)."
            " Only the 'upper triangular' submatrices need to be evaluated due to "
            "the symmetry of the relatedness matrix, leading to 10 shards."
        ),
        default=4,
        type=int,
    )
    parser.add_argument(
        "--min-kin-cutoff",
        help=(
            "Minimum kinship threshold for a pair of samples to be retained in the "
            "cuKING output."
        ),
        default=0.05,
        type=float,
    )
    parser.add_argument(
        "--create-cuking-relatedness-table",
        help="Convert the cuKING outputs to a standard Hail Table.",
        action="store_true",
    )
    parser.add_argument(
        "--relatedness-n-partitions",
        help="Number of desired partitions for the relatedness Table.",
        default=50,
        type=int,
    )
    parser.add_argument(
        "--compute-related-samples-to-drop",
        help="Determine the minimal set of related samples to prune.",
        action="store_true",
    )
    parser.add_argument(
        "--second-degree-kin-cutoff",
        help=(
            "Minimum kinship threshold for filtering a pair of samples with a second "
            "degree relationship when filtering related individuals. Default is 0.1 "
            "Bycroft et al. (2018) calculates a theoretical kinship of 0.08838835 "
            "for a second degree relationship cutoff, but 0.1 was decided on after "
            "evaluation of the distributions in gnomAD v3 and v4."
        ),
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
