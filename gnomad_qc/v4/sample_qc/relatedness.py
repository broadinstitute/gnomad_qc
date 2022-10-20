"""Script to compute relatedness estimates among pairs of samples in the callset."""
import argparse
import logging
import textwrap

import hail as hl
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.utils.file_utils import check_file_exists_raise_error
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import check_resource_existence, get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    get_cuking_input_path,
    get_cuking_output_path,
    get_joint_qc,
    joint_qc_meta,
    pca_related_samples_to_drop,
    pca_samples_rankings,
    relatedness,
)
from gnomad_qc.v4.sample_qc.cuKING.cuking_outputs_to_ht import cuking_outputs_to_ht
from gnomad_qc.v4.sample_qc.cuKING.mt_to_cuking_inputs import mt_to_cuking_inputs

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def main(args):
    """Compute relatedness estimates among pairs of samples in the callset."""
    test = args.test
    overwrite = args.overwrite
    min_emission_kinship = args.min_emission_kinship

    joint_qc_mt = get_joint_qc(test=test)
    cuking_input_path = get_cuking_input_path(test=test)
    cuking_output_path = get_cuking_output_path(test=test)
    cuking_relatedness_ht = relatedness(test=test)
    pc_relate_relatedness_ht = relatedness("pc_relate", test=test)

    if args.print_cuking_command:
        check_resource_existence(
            input_pipeline_step="--prepare-cuking-inputs",
            input_resources=[cuking_input_path],
        )
        logger.info(
            "Printing a command that can be used to submit a Cloud Batch job to run "
            "cuKING on the files created by --prepare-cuking-inputs."
        )
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
                    --input-uri={cuking_input_path} \\
                    --output-uri={cuking_output_path} \\
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
        if args.prepare_cuking_inputs:
            logger.info(
                "Converting joint dense QC MatrixTable to a Parquet format suitable "
                "for cuKING."
            )
            check_resource_existence(
                input_pipeline_step="generate_qc_mt.py --generate-qc-mt",
                output_pipeline_step="--prepare-cuking-inputs",
                input_resources=[joint_qc_mt],
                output_resources=[cuking_input_path],
                overwrite=overwrite,
            )
            mt_to_cuking_inputs(
                mt=joint_qc_mt.mt(),
                parquet_uri=cuking_input_path,
                overwrite=overwrite,
            )

        if args.create_cuking_relatedness_table:
            logger.info("Converting cuKING outputs to Hail Table.")
            check_resource_existence(
                input_pipeline_step="--print-cuking-command",
                output_pipeline_step="--create-cuking-relatedness-table",
                input_resources=[cuking_output_path],
                output_resources=[cuking_relatedness_ht],
                overwrite=overwrite,
            )

            ht = cuking_outputs_to_ht(parquet_uri=cuking_output_path)
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(cuking_relatedness_ht.path, overwrite=overwrite)

        if args.run_pc_relate_pca:
            logger.info("Running PCA for PC-Relate")
            check_resource_existence(
                input_pipeline_step="generate_qc_mt.py --generate-qc-mt",
                output_pipeline_step="--run-pc-relate-pca",
                input_resources=[joint_qc_mt],
                output_resources=[pc_relate_pca_scores],
                overwrite=overwrite,
            )
            eig, scores, _ = hl.hwe_normalized_pca(
                joint_qc_mt.mt().GT, k=args.n_pca_pcs, compute_loadings=False
            )
            scores.write(pc_relate_pca_scores.path, overwrite=overwrite)

        if args.create_pc_relate_relatedness_table:
            logger.info("Running PC-Relate")
            logger.warning(
                "PC-relate requires SSDs and doesn't work with preemptible workers!"
            )
            check_resource_existence(
                input_pipeline_step=(
                    "generate_qc_mt.py --generate-qc-mt and --run-pc-relate-pca"
                ),
                output_pipeline_step="--create-pc-relate-relatedness-table",
                input_resources=[joint_qc_mt, pc_relate_pca_scores],
                output_resources=[pc_relate_relatedness_ht],
                overwrite=overwrite,
            )
            mt = joint_qc_mt.mt()
            ht = hl.pc_relate(
                mt.GT,
                min_individual_maf=args.min_individual_maf,
                scores_expr=pc_relate_pca_scores.ht()[mt.col_key].scores[
                    : args.n_pc_relate_pcs
                ],
                block_size=args.block_size,
                min_kinship=min_emission_kinship,
                statistics="all",
            )
            ht = ht.annotate_globals(
                min_individual_maf=args.min_individual_maf,
                min_emission_kinship=min_emission_kinship,
            )
            ht.write(pc_relate_relatedness_ht.path, overwrite=overwrite)

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
