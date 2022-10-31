"""Script to compute relatedness estimates among pairs of samples in the callset."""
import argparse
import logging
import textwrap

import hail as hl
import networkx as nx
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.utils.file_utils import check_file_exists_raise_error
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import check_resource_existence, get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    get_cuking_input_path,
    get_cuking_output_path,
    get_joint_qc,
    ibd,
    joint_qc_meta,
    pc_relate_pca_scores,
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
    ibd_ht = ibd(test=test)
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
            "location where the command is being run and that $PROJECT_ID is set!"
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
                    --kin-threshold={min_emission_kinship} \\
                    --split-factor={args.cuking_split_factor} &&
                cd ..
                """
            )
        )
        return

    hl.init(
        log="/relatedness.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

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

        if args.run_ibd_on_cuking_pairs:
            logger.info(
                "Running IBD on cuKING pairs with kinship over %f.",
                args.ibd_min_cuking_kin,
            )
            check_resource_existence(
                input_pipeline_step=(
                    "generate_qc_mt.py --generate-qc-mt and "
                    "--create-cuking-relatedness-table"
                ),
                output_pipeline_step="--run-ibd-on-cuking-pairs",
                input_resources=[joint_qc_mt, cuking_relatedness_ht],
                output_resources=[ibd_ht],
                overwrite=overwrite,
            )
            mt = joint_qc_mt.mt()
            ht = cuking_relatedness_ht.ht()
            ht = ht.filter(ht.kin > args.ibd_min_cuking_kin)
            ht = ht.checkpoint(new_temp_file("cuking_for_ibd", extension="ht"))
            ht_i = ht.i.collect()
            ht_j = ht.j.collect()
            mt = mt.filter_cols(hl.literal(set(ht_i) | set(ht_j)).contains(mt.s))
            mt = mt.checkpoint(new_temp_file("cuking_for_ibd", extension="mt"))

            # Build a network from the list of pairs and extract the list of
            # connected components
            pair_graph = nx.Graph()
            pair_graph.add_edges_from(list(zip(ht_i, ht_j)))
            connected_samples = list(nx.connected_components(pair_graph))

            ibd_hts = []
            sample_subset = []
            while connected_samples:
                sample_subset.extend(connected_samples.pop())
                if len(sample_subset) > args.ibd_max_samples or not connected_samples:
                    print(len(ibd_hts))
                    subset_ibd_ht = hl.identity_by_descent(
                        mt.filter_cols(hl.literal(sample_subset).contains(mt.s))
                    )
                    subset_ibd_ht = subset_ibd_ht.filter(
                        hl.is_defined(
                            hl.coalesce(
                                ht[subset_ibd_ht.i, subset_ibd_ht.j],
                                ht[subset_ibd_ht.j, subset_ibd_ht.i],
                            )
                        )
                    )
                    ibd_hts.append(
                        subset_ibd_ht.checkpoint(
                            new_temp_file("relatedness_ibd_subset", extension="ht")
                        )
                    )
                    sample_subset = []
            full_ibd_ht = ibd_hts[0].union(*ibd_hts[1:])
            full_ibd_ht = full_ibd_ht.annotate_globals(
                ibd_min_cuking_kin=args.ibd_min_cuking_kin
            )
            full_ibd_ht.distinct().write(ibd_ht.path, overwrite=overwrite)

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
            ht = ht.repartition(args.relatedness_n_partitions)
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

    relatedness_estimate_args = parser.add_argument_group(
        "Common relatedness estimate arguments",
        "Arguments relevant to both cuKING and PC-relate relatedness estimates.",
    )
    relatedness_estimate_args.add_argument(
        "--min-emission-kinship",
        help=(
            "Minimum kinship threshold for emitting a pair of samples in the "
            "relatedness output."
        ),
        default=0.05,
        type=float,
    )
    relatedness_estimate_args.add_argument(
        "--relatedness-n-partitions",
        help="Number of desired partitions for the relatedness Table.",
        default=100,
        type=int,
    )

    cuking_args = parser.add_argument_group(
        "cuKING specific relatedness arguments",
        "Arguments specific to computing relatedness estimates using cuKING.",
    )
    cuking_args.add_argument(
        "--prepare-cuking-inputs",
        help=(
            "Converts the dense QC MatrixTable to a Parquet format suitable for cuKING."
        ),
        action="store_true",
    )
    print_cuking_cmd = cuking_args.add_argument_group(
        "Print cuKING Cloud Batch job submission command",
        "Arguments used to create the cuKING Cloud Batch job submission command "
        "needed to run cuKING.",
    )
    print_cuking_cmd.add_argument(
        "--print-cuking-command",
        help="Print the command to submit a Cloud Batch job for running cuKING.",
        action="store_true",
    )
    print_cuking_cmd.add_argument(
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
    cuking_args.add_argument(
        "--create-cuking-relatedness-table",
        help="Convert the cuKING outputs to a standard Hail Table.",
        action="store_true",
    )
    print_cuking_ibd = cuking_args.add_argument_group(
        "Run IBD on related pairs identified by cuKING",
        "Arguments used run IBD on related cuKING pairs.",
    )
    print_cuking_ibd.add_argument(
        "--run-ibd-on-cuking-pairs",
        help="Run IBD on related pairs identified by cuKING.",
        action="store_true",
    )
    print_cuking_ibd.add_argument(
        "--ibd-min-cuking-kin",
        help="Min cuKING kinship for pair to be included in IBD estimates.",
        default=0.16,
        type=float,
    )
    print_cuking_ibd.add_argument(
        "--ibd-max-samples",
        help="Max number of samples to include in each IBD run.",
        default=10000,
        type=int,
    )

    pc_relate_args = parser.add_argument_group(
        "PC-relate specific relatedness arguments",
        "Arguments specific to computing relatedness estimates using PC-relate.",
    )
    run_pca = pc_relate_args.add_argument_group(
        "Run PCA for PC-relate", "Arguments used to run the PCA for PC-relate."
    )
    run_pca.add_argument(
        "--run-pc-relate-pca",
        help="Run PCA to generate the scores needed for PC-relate.",
        action="store_true",
    )
    run_pca.add_argument(
        "--n-pca-pcs",
        help="Number of PCs to compute for PC-relate.",
        type=int,
        default=20,
    )
    run_pc_relate = pc_relate_args.add_argument_group(
        "Run PC-relate", "Arguments used to run PC-relate."
    )
    run_pc_relate.add_argument(
        "--create-pc-relate-relatedness-table",
        help="Run PC-relate to create the PC-relate relatedness Table.",
        action="store_true",
    )
    run_pc_relate.add_argument(
        "--n-pc-relate-pcs",
        help="Number of PCs to compute for PC-relate.",
        type=int,
        default=10,
    )
    run_pc_relate.add_argument(
        "--min-individual-maf",
        help="The minimum individual-specific minor allele frequency.",
        default=0.01,
        type=float,
    )
    run_pc_relate.add_argument(
        "--block-size",
        help="Block size parameter to use for PC-relate.",
        default=2048,
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

    if args.print_cuking_command and (
        args.prepare_cuking_inputs
        or args.create_cuking_relatedness_table
        or args.run_ibd_on_cuking_pairs
        or args.run_pc_relate_pca
        or args.create_pc_relate_relatedness_table
        or args.compute_related_samples_to_drop
    ):
        parser.error(
            "--print-cuking-command can't be used simultaneously with other run modes."
        )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
