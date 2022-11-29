"""Script to compute relatedness estimates among pairs of samples in the callset."""
import argparse
import logging
import textwrap

import hail as hl
import networkx as nx
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
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


def compute_ibd(mt, relatedness_ht, ibd_min_cuking_kin, ibd_max_cuking_ibs0):
    """
    Add summary.

    :param mt:
    :param relatedness_ht:
    :param ibd_min_cuking_kin:
    :param ibd_max_cuking_ibs0:
    :return:
    """
    ht = relatedness_ht.filter(
        (relatedness_ht.kin > ibd_min_cuking_kin)
        | (relatedness_ht.ibs0 < ibd_max_cuking_ibs0)
    )
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

    return full_ibd_ht.distinct()


def compute_rank_ht(rank_ht):
    """
    Add summary.

    :param rank_ht:
    :return:
    """
    rank_ht = rank_ht.select(
        rank_ht.releasable,
        rank_ht.chr20_mean_dp,
        v3_release=hl.is_defined(rank_ht.v3_meta.v3_release)
        & rank_ht.v3_meta.v3_release,
        present_in_v3=hl.is_defined(rank_ht.v3_meta),
    )

    # Favor v3 release samples, then v4 samples over v3 non-release samples.
    rank_ht = rank_ht.order_by(
        hl.desc(rank_ht.v3_release),
        hl.desc(rank_ht.releasable),
        hl.asc(rank_ht.present_in_v3),
        hl.desc(rank_ht.chr20_mean_dp),
    ).add_index(name="rank")
    rank_ht = rank_ht.key_by(rank_ht.s)
    rank_ht = rank_ht.select(rank_ht.hard_filtered, rank_ht.rank)
    return rank_ht


def compute_related_samples_to_drop2(
    relatedness_ht,
    rank_ht,
    second_degree_kin_cutoff,
    must_keep=None,
):
    """
    Add summary.

    :param relatedness_ht:
    :param rank_ht:
    :param second_degree_kin_cutoff:
    :param must_keep:
    :return:
    """
    relatedness_ht = relatedness_ht.key_by(
        i=relatedness_ht.i.annotate(rank=rank_ht[relatedness_ht.i.s].rank),
        j=relatedness_ht.j.annotate(rank=rank_ht[relatedness_ht.j.s].rank),
    )
    relatedness_ht = relatedness_ht.filter(
        relatedness_ht.kin > second_degree_kin_cutoff
    )
    related_pair_graph = nx.Graph()
    related_pair_graph.add_edges_from(
        list(zip(relatedness_ht.i.collect(), relatedness_ht.j.collect()))
    )

    def remove_related(g, must_keep=None):
        to_remove = []
        cant_remove = []
        connected_samples = list(nx.connected_components(g))
        for con in connected_samples:
            if len(con) > 1:
                con_idx = {i: s for i, s in enumerate(con)}
                degree_rank_list = sorted(
                    [(g.degree[s], s["rank"], i) for i, s in con_idx.items()]
                )
                last_i = len(degree_rank_list) - 1
                last_s = con_idx[degree_rank_list[last_i][2]]
                if must_keep is not None:
                    while last_s.s in must_keep and last_i > 0:
                        last_i -= 1
                        last_s = con_idx[degree_rank_list[last_i][2]]
                    if last_s.s in must_keep:
                        all_s = [con_idx[s[2]].s for s in degree_rank_list]
                        cant_remove.extend(all_s)
                        last_s = None
                        print(
                            "MUST REMOVE ONE: \n\t",
                            "\n\t ".join(
                                map(str, [con_idx[s[2]] for s in degree_rank_list])
                            ),
                        )
                if last_s is not None:
                    to_remove.append(last_s)
                    new_con = {con_idx[s[2]] for s in degree_rank_list[:-1]}
                    if len(new_con) > 1:
                        new_g = g.subgraph(new_con)
                        to_remove.extend(remove_related(new_g))

        return to_remove

    return remove_related(
        related_pair_graph,
        must_keep=set(
            joint_qc_meta.filter(joint_qc_meta.v3_meta.v3_release).s.collect()
        ),
    )


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
    pca_samples_rankings_ht = pca_samples_rankings(test=test)
    pca_related_samples_to_drop_ht = pca_related_samples_to_drop(test=test)

    if args.print_cuking_command:
        check_resource_existence(
            input_step_resources={"--prepare-cuking-inputs": [cuking_input_path]}
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
                    --service-account=cuking@$PROJECT_ID.iam.gserviceaccount.com \\
                    --write-success-file \\
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
                input_step_resources={
                    "generate_qc_mt.py --generate-qc-mt": [joint_qc_mt]
                },
                output_step_resources={"--prepare-cuking-inputs": [cuking_input_path]},
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
                input_step_resources={"--print-cuking-command": [cuking_output_path]},
                output_step_resources={
                    "--create-cuking-relatedness-table": [cuking_relatedness_ht]
                },
                overwrite=overwrite,
            )

            ht = cuking_outputs_to_ht(parquet_uri=cuking_output_path)
            ht = ht.key_by(i=hl.struct(s=ht.i), j=hl.struct(s=ht.j))
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(cuking_relatedness_ht.path, overwrite=overwrite)

        if args.run_ibd_on_cuking_pairs:
            logger.info(
                "Running IBD on cuKING pairs with kinship over %f.",
                args.ibd_min_cuking_kin,
            )
            check_resource_existence(
                input_step_resources={
                    "generate_qc_mt.py --generate-qc-mt": [joint_qc_mt],
                    "--create-cuking-relatedness-table": [cuking_relatedness_ht],
                },
                output_step_resources={"--run-ibd-on-cuking-pairs": [ibd_ht]},
                overwrite=overwrite,
            )

            compute_ibd(
                joint_qc_mt.mt(),
                cuking_relatedness_ht.ht(),
                args.ibd_min_cuking_kin,
                args.ibd_max_cuking_ibs0,
            ).write(ibd_ht.path, overwrite=overwrite)

        if args.run_pc_relate_pca:
            logger.info("Running PCA for PC-Relate")
            check_resource_existence(
                input_step_resources={
                    "generate_qc_mt.py --generate-qc-mt": [joint_qc_mt]
                },
                output_step_resources={"--run-pc-relate-pca": [pc_relate_pca_scores]},
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
                input_step_resources={
                    "generate_qc_mt.py --generate-qc-mt": [joint_qc_mt],
                    "--run-pc-relate-pca": [pc_relate_pca_scores],
                },
                output_step_resources={
                    "--create-pc-relate-relatedness-table": [pc_relate_relatedness_ht]
                },
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
            rel_method = args.relatedness_method_for_samples_to_drop
            relatedness_ht = relatedness(rel_method, test=test)

            # compute_related_samples_to_drop uses a rank Table as a tiebreaker when
            # pruning samples.
            check_resource_existence(
                input_step_resources={
                    "generate_qc_mt.py --generate-qc-mt": [joint_qc_mt],
                    f"--create-{rel_method.replace('_', '-')}-relatedness-table": [
                        relatedness_ht
                    ],
                },
                output_step_resources={
                    "--compute-related-samples-to-drop": [
                        pca_samples_rankings_ht,
                        pca_related_samples_to_drop_ht,
                    ]
                },
                overwrite=overwrite,
            )

            rank_ht = compute_rank_ht(
                joint_qc_meta.ht().semi_join(joint_qc_mt.mt().cols())
            )
            rank_ht = rank_ht.checkpoint(
                pca_samples_rankings_ht.path, overwrite=overwrite
            )

            samples_to_drop = compute_related_samples_to_drop2(
                relatedness_ht.ht(),
                rank_ht,
                args.second_degree_kin_cutoff,
            )
            # TODO: Annotate globals?
            samples_to_drop = samples_to_drop.key_by(samples_to_drop.s)
            samples_to_drop.write(
                pca_related_samples_to_drop_ht.path, overwrite=overwrite
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
    cuking_ibd = cuking_args.add_argument_group(
        "Run IBD on related pairs identified by cuKING",
        "Arguments used run IBD on related cuKING pairs.",
    )
    cuking_ibd.add_argument(
        "--run-ibd-on-cuking-pairs",
        help="Run IBD on related pairs identified by cuKING.",
        action="store_true",
    )
    cuking_ibd.add_argument(
        "--ibd-min-cuking-kin",
        help="Min cuKING kinship for pair to be included in IBD estimates.",
        default=0.16,
        type=float,
    )
    cuking_ibd.add_argument(
        "--ibd-max-cuking-ibs0",
        help=(
            "Max cuKING IBS0 for pair to be included in IBD estimates. Note: This "
            "default was determined from looking at the cuKING Kinship vs. IBS0 plot "
            "for gnomAD v3 + gnomAD v4."
        ),
        default=50,
        type=int,
    )
    cuking_ibd.add_argument(
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
        help="Number of PCs to use in PC-relate.",
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
    related_samples_to_drop = pc_relate_args.add_argument_group(
        "Compute related samples to drop",
        "Arguments used to determine related samples that should be dropped from "
        "the ancestry PCA.",
    )
    related_samples_to_drop.add_argument(
        "--compute-related-samples-to-drop",
        help="Determine the minimal set of related samples to prune for ancestry PCA.",
        action="store_true",
    )
    related_samples_to_drop.add_argument(
        "--relatedness-method-for-samples-to-drop",
        help=(
            "Which relatedness method to use when determining related samples to drop. "
            "Options are 'cuking' and 'pc_relate'. Default is 'cuking'."
        ),
        default="cuking",
        type=str,
        choices=["cuking", "pc_relate"],
    )
    related_samples_to_drop.add_argument(
        "--second-degree-kin-cutoff",
        help=(
            "Minimum kinship threshold for filtering a pair of samples with a second "
            "degree relationship when filtering related individuals. Default is "
            "0.08838835. Bycroft et al. (2018) calculates a theoretical kinship of "
            "0.08838835 for a second degree relationship cutoff. This cutoff should"
            "be determined by evaluation of the kinship distribution."
        ),
        default=0.08838835,
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
