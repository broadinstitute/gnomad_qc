"""Script to compute relatedness estimates among pairs of samples in the callset."""
import argparse
import logging
import textwrap
from typing import Dict, Optional, Tuple

import hail as hl
import networkx as nx
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    compute_related_samples_to_drop,
    get_slope_int_relationship_expr,
)
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
    get_logging_path,
)
from gnomad_qc.v4.resources.sample_qc import (
    finalized_outlier_filtering,
    get_cuking_input_path,
    get_cuking_output_path,
    get_joint_qc,
    ibd,
    joint_qc_meta,
    pc_relate_pca_scores,
    related_samples_to_drop,
    relatedness,
    sample_rankings,
)
from gnomad_qc.v4.sample_qc.cuKING.cuking_outputs_to_ht import cuking_outputs_to_ht
from gnomad_qc.v4.sample_qc.cuKING.mt_to_cuking_inputs import mt_to_cuking_inputs

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)

DATA_TYPES = {"genomes", "exomes"}


def print_cuking_command(
    cuking_input_path: str,
    cuking_output_path: str,
    min_emission_kinship: float = 0.5,
    cuking_split_factor: int = 4,
) -> None:
    """
    Add description.

    :param cuking_input_path:
    :param cuking_output_path:
    :param min_emission_kinship:
    :param cuking_split_factor:
    :return:
    """
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
                 --split-factor={cuking_split_factor} &&
             cd ..
             """
        )
    )


def compute_ibd_on_cuking_pair_subset(
    mt: hl.Table,
    relatedness_ht: hl.Table,
    ibd_min_cuking_kin: float = 0.16,
    ibd_max_cuking_ibs0: int = 50,
    ibd_max_samples: int = 10000,
) -> hl.Table:
    """
    Run Hail's identity by descent on a subset of related pairs identified by cuKING.

    The pairs that will get an identity by descent annotation are those where either
    `ibd_min_cuking_kin` or `ibd_max_cuking_ibs0` are met.

    :param mt: QC MatrixTable.
    :param relatedness_ht: cuKING relatedness Table.
    :param ibd_min_cuking_kin: Minimum cuKING kinship for pair to be included in IBD
        estimates. Default is 0.16.
    :param ibd_max_cuking_ibs0: Maximum cuKING IBS0 for pair to be included in IBD
        estimates. Default is 50. This default was determined from looking at the cuKING
        Kinship vs. IBS0 plot for gnomAD v3 + gnomAD v4.
    :param ibd_max_samples: Maximum number of samples to include in each IBD run.
    :return: Table containing identity by descent metrics on related sample pairs.
    """
    logger.info(
        "Filtering the relatedness HT to pairs with cuKING kinship greater than %f or "
        "IBS0 less than %d for IBD annotation.",
        ibd_min_cuking_kin,
        ibd_max_cuking_ibs0,
    )
    ht = relatedness_ht.filter(
        (relatedness_ht.kin > ibd_min_cuking_kin)
        | (relatedness_ht.ibs0 < ibd_max_cuking_ibs0)
    )
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)
    ht = ht.checkpoint(new_temp_file("cuking_for_ibd", extension="ht"))

    logger.info(
        "Building a network from the list of pairs and extracting the list of "
        "connected components."
    )
    pair_graph = nx.Graph()
    pair_graph.add_edges_from(list(zip(ht.i.collect(), ht.j.collect())))
    connected_samples = list(nx.connected_components(pair_graph))

    logger.info(
        "Looping through %d connected components that include a total of %d samples. "
        "hl.identity_by_descent is run on at most %d samples at a time, where all "
        "samples from a connected component are in the same IBD run.",
        len(connected_samples),
        sum([len(c) for c in connected_samples]),
        ibd_max_samples,
    )
    ibd_hts = []
    sample_subset = []
    while connected_samples:
        sample_subset.extend(connected_samples.pop())
        if len(sample_subset) > ibd_max_samples or not connected_samples:
            ibd_ht = hl.identity_by_descent(
                mt.filter_cols(hl.literal(sample_subset).contains(mt.s))
            )
            # Order the ibd output so that i and j are the same in `ht` and `ibd_ht`.
            ht_idx1 = hl.is_defined(ht[ibd_ht.i, ibd_ht.j])
            ht_idx2 = hl.is_defined(ht[ibd_ht.j, ibd_ht.i])
            ibd_ht = ibd_ht.key_by(
                i=hl.struct(s=hl.if_else(ht_idx1, ibd_ht.i, ibd_ht.j)),
                j=hl.struct(s=hl.if_else(ht_idx2, ibd_ht.j, ibd_ht.i)),
            )
            # Keep only pairs that were in the original relatedness HT.
            ibd_ht = ibd_ht.filter(hl.is_defined(ibd_ht.key))
            ibd_hts.append(
                ibd_ht.checkpoint(
                    new_temp_file("relatedness_ibd_subset", extension="ht")
                )
            )
            sample_subset = []

    full_ibd_ht = ibd_hts[0].union(*ibd_hts[1:])
    full_ibd_ht = full_ibd_ht.annotate_globals(
        ibd_min_cuking_kin=args.ibd_min_cuking_kin
    )

    return full_ibd_ht.distinct()


def finalize_relatedness_ht(
    ht: hl.Table,
    meta_ht: hl.Table,
    relatedness_method: str,
    relatedness_args: Dict[str, float],
) -> hl.Table:
    """
    Add description.

    :param ht:
    :param meta_ht:
    :param relatedness_method:
    :param relatedness_args:
    :return:
    """
    if relatedness_method == "pc-relate":
        y_expr = ht.ibd0
        ibd1_expr = ht.ibd1
        parent_child_max_ann = "parent_child_max_ibd0"
        relatedness_args.pop("duplicate_twin_ibd1_min")
        relatedness_args.pop("duplicate_twin_ibd1_max")
    else:
        y_expr = ht.ibs0 / ht.ibs2
        ibd1_expr = None
        parent_child_max_ann = "parent_child_max_ibs0_over_ibs2"

    ht = ht.annotate(
        relationship=get_slope_int_relationship_expr(
            kin_expr=ht.kin,
            y_expr=y_expr,
            **relatedness_args,
            ibd1_expr=ibd1_expr,
        )
    )
    relatedness_args[parent_child_max_ann] = relatedness_args.pop("parent_child_max_y")
    ht = ht.annotate_globals(
        relationship_inference_method=relatedness_method,
        relationship_cutoffs=hl.struct(**relatedness_args),
    )
    ht = ht.key_by(
        i=ht.i.annotate(data_type=meta_ht[ht.i.s].data_type),
        j=ht.j.annotate(data_type=meta_ht[ht.j.s].data_type),
    )
    gnomad_v3_duplicate_expr = (ht.relationship == DUPLICATE_OR_TWINS) & (
        hl.set({ht.i.data_type, ht.j.data_type}) == hl.set(DATA_TYPES)
    )
    gnomad_v3_release_expr = hl.coalesce(
        meta_ht[ht.i.s].v3_meta.v3_release | meta_ht[ht.j.s].v3_meta.v3_release,
        False,
    )
    ht = ht.annotate(
        gnomad_v3_duplicate=gnomad_v3_duplicate_expr,
        gnomad_v3_release_duplicate=gnomad_v3_duplicate_expr & gnomad_v3_release_expr,
    )

    return ht


def compute_rank_ht(ht: hl.Table, filter_ht: Optional[hl.Table] = None) -> hl.Table:
    """
    Add a rank to each sample for use when breaking maximal independent set ties.

    Favor v3 release samples, then v4 samples over v3 non-release samples, then
    higher chr20 mean DP.

    If `filter_ht` is provided, rank based on filtering 'outlier_filtered' annotation
    first.

    :param ht: Table to add rank to.
    :param filter_ht: Optional Table with 'outlier_filtered' annotation to be used in
        ranking.
    :return: Table containing sample ID and rank.
    """
    ht = ht.select(
        "releasable",
        "chr20_mean_dp",
        in_v3=hl.is_defined(ht.v3_meta),
        in_v3_release=hl.is_defined(ht.v3_meta) & ht.v3_meta.v3_release,
    )
    rank_order = []
    ht_select = ["rank"]
    if filter_ht is not None:
        ht = ht.annotate(
            filtered=hl.coalesce(filter_ht[ht.key].outlier_filtered, False)
        )
        rank_order = [ht.filtered]
        ht_select.append("filtered")

    rank_order.extend(
        [
            hl.desc(ht.in_v3_release),
            hl.asc(ht.in_v3),
            hl.desc(ht.releasable),
            hl.desc(ht.chr20_mean_dp),
        ]
    )
    ht = ht.order_by(*rank_order).add_index(name="rank")
    ht = ht.key_by(ht.s)

    return ht.select(*ht_select)


def run_compute_related_samples_to_drop(
    ht: hl.Table,
    meta_ht: hl.Table,
    filter_ht: hl.Table,
    release: bool,
) -> Tuple[hl.Table, hl.Table]:
    """
    Add summary.

    :param ht:
    :param meta_ht:
    :param filter_ht:
    :param release:
    :return:
    """
    # Compute_related_samples_to_drop uses a rank Table as a tiebreaker when
    # pruning samples.
    rank_ht = compute_rank_ht(meta_ht, filter_ht=filter_ht if release else None)
    rank_ht = rank_ht.checkpoint(new_temp_file("rank", extension="ht"))

    filtered_samples = None
    if release:
        filtered_samples = rank_ht.aggregate(
            hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s)),
            _localize=False,
        )
    second_degree_kin_cutoff = hl.eval(ht.relationship_cutoffs.second_degree_min_kin)
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)
    v3_release_samples = meta_ht.aggregate(
        hl.agg.filter(meta_ht.v3_meta.v3_release, hl.agg.collect_as_set(meta_ht.s)),
        _localize=False,
    )
    samples_to_drop_ht = compute_related_samples_to_drop(
        ht,
        rank_ht,
        second_degree_kin_cutoff,
        filtered_samples=filtered_samples,
        keep_samples=v3_release_samples,
        keep_samples_when_related=True,
    )
    samples_to_drop_ht = samples_to_drop_ht.annotate_globals(
        keep_samples=v3_release_samples,
        second_degree_kin_cutoff=second_degree_kin_cutoff,
    )

    return rank_ht, samples_to_drop_ht


def get_relatedness_resources(
    test: bool,
    release: bool,
    relatedness_method: str,
    overwrite: bool,
) -> PipelineResourceCollection:
    """
    Add description.

    :param test:
    :param release:
    :param relatedness_method:
    :param overwrite:
    :return:
    """
    joint_qc_mt = get_joint_qc(test=test)
    joint_qc_mt_input = {
        "generate_qc_mt.py --generate-qc-mt": {"joint_qc_mt": joint_qc_mt}
    }
    joint_qc_meta_input = {
        "generate_qc_mt.py --generate-qc-meta": {"joint_qc_meta": joint_qc_meta}
    }

    relatedness_pipeline = PipelineResourceCollection(
        pipeline_name="relatedness",
        pipeline_resources={"joint_qc_mt": joint_qc_mt, "joint_qc_meta": joint_qc_meta},
        overwrite=overwrite,
    )

    prepare_cuking_inputs = PipelineStepResourceCollection(
        "--prepare-cuking-inputs",
        output_resources={"cuking_input_path": get_cuking_input_path(test=test)},
        input_resources=joint_qc_mt_input,
    )
    print_cuking_command = PipelineStepResourceCollection(
        "--print-cuking-command",
        previous_pipeline_steps=[prepare_cuking_inputs],
        output_resources={"cuking_output_path": get_cuking_output_path(test=test)},
    )
    create_cuking_relatedness_table = PipelineStepResourceCollection(
        "--create-cuking-relatedness-table",
        previous_pipeline_steps=[print_cuking_command],
        output_resources={"cuking_relatedness_ht": relatedness("cuking", test=test)},
    )
    run_ibd_on_cuking_pairs = PipelineStepResourceCollection(
        "--run-ibd-on-cuking-pairs",
        previous_pipeline_steps=[create_cuking_relatedness_table],
        add_input_resources=joint_qc_mt_input,
        output_resources={"ibd_ht": ibd(test=test)},
    )
    run_pc_relate_pca = PipelineStepResourceCollection(
        "--run-pc-relate-pca",
        input_resources=joint_qc_mt_input,
        output_resources={"pc_relate_pca_scores": pc_relate_pca_scores},
    )
    create_pc_relate_relatedness_table = PipelineStepResourceCollection(
        "--create-pc-relate-relatedness-table",
        previous_pipeline_steps=[run_pc_relate_pca],
        add_input_resources=joint_qc_mt_input,
        output_resources={
            "pc_relate_relatedness_ht": relatedness("pc_relate", test=test)
        },
    )
    finalize_relatedness_ht = PipelineStepResourceCollection(
        "--finalize-relatedness-ht",
        previous_pipeline_steps=[
            create_cuking_relatedness_table
            if relatedness_method == "cuking"
            else create_pc_relate_relatedness_table
        ],
        add_input_resources=joint_qc_meta_input,
        output_resources={"final_relatedness_ht": relatedness(test=test)},
    )
    compute_related_samples_to_drop = PipelineStepResourceCollection(
        "--compute-related-samples-to-drop",
        previous_pipeline_steps=[finalize_relatedness_ht],
        add_input_resources={
            **joint_qc_mt_input,
            **joint_qc_meta_input,
            **(
                {
                    "outlier_filtering.py --create-finalized-outlier-filter": {
                        "outlier_filter_ht": finalized_outlier_filtering(test=test)
                    }
                }
                if release
                else {}
            ),
        },
        output_resources={
            "sample_rankings_ht": sample_rankings(test=test, release=release),
            "related_samples_to_drop_ht": related_samples_to_drop(
                test=test, release=release
            ),
        },
    )

    relatedness_pipeline.add_steps(
        {
            "prepare_cuking_inputs": prepare_cuking_inputs,
            "print_cuking_command": print_cuking_command,
            "create_cuking_relatedness_table": create_cuking_relatedness_table,
            "run_ibd_on_cuking_pairs": run_ibd_on_cuking_pairs,
            "run_pc_relate_pca": run_pc_relate_pca,
            "create_pc_relate_relatedness_table": create_pc_relate_relatedness_table,
            "finalize_relatedness_ht": finalize_relatedness_ht,
            "compute_related_samples_to_drop": compute_related_samples_to_drop,
        }
    )

    return relatedness_pipeline


def main(args):
    """Compute relatedness estimates among pairs of samples in the callset."""
    test = args.test
    overwrite = args.overwrite
    min_emission_kinship = args.min_emission_kinship
    release = args.release
    rel_method = args.finalize_relatedness_method
    second_degree_kin_cutoff = args.second_degree_kin_cutoff
    relatedness_resources = get_relatedness_resources(
        test, release, rel_method, overwrite
    )

    if args.print_cuking_command:
        res = relatedness_resources.print_cuking_command
        res.check_resource_existance()
        print_cuking_command(
            res.cuking_input_path,
            res.cuking_output_path,
            min_emission_kinship,
            args.cuking_split_factor,
        )
        return

    hl.init(
        log="/relatedness.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    joint_qc_mt = relatedness_resources.joint_qc_mt.mt()
    joint_qc_meta_ht = relatedness_resources.joint_qc_meta.ht()

    try:
        if args.prepare_cuking_inputs:
            logger.info(
                "Converting joint dense QC MatrixTable to a Parquet format suitable "
                "for cuKING."
            )
            res = relatedness_resources.prepare_cuking_inputs
            res.check_resource_existance()
            mt_to_cuking_inputs(
                mt=joint_qc_mt,
                parquet_uri=res.cuking_input_path,
                overwrite=overwrite,
            )

        if args.create_cuking_relatedness_table:
            logger.info("Converting cuKING outputs to Hail Table.")
            res = relatedness_resources.prepare_cuking_inputs
            res.check_resource_existance()
            ht = cuking_outputs_to_ht(parquet_uri=res.cuking_output_path)
            ht = ht.key_by(i=hl.struct(s=ht.i), j=hl.struct(s=ht.j))
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(res.cuking_relatedness_ht.path, overwrite=overwrite)

        if args.run_ibd_on_cuking_pairs:
            logger.info(
                "Running IBD on cuKING pairs with kinship over %f.",
                args.ibd_min_cuking_kin,
            )
            res = relatedness_resources.run_ibd_on_cuking_pairs
            res.check_resource_existance()
            compute_ibd_on_cuking_pair_subset(
                joint_qc_mt,
                res.cuking_relatedness_ht.ht(),
                args.ibd_min_cuking_kin,
                args.ibd_max_cuking_ibs0,
                args.ibd_max_samples,
            ).write(res.ibd_ht.path, overwrite=overwrite)

        if args.run_pc_relate_pca:
            logger.info("Running PCA for PC-Relate")
            res = relatedness_resources.run_pc_relate_pca
            res.check_resource_existance()
            eig, scores, _ = hl.hwe_normalized_pca(
                joint_qc_mt.GT, k=args.n_pca_pcs, compute_loadings=False
            )
            scores.write(res.pc_relate_pca_scores.path, overwrite=overwrite)

        if args.create_pc_relate_relatedness_table:
            logger.info("Running PC-Relate")
            logger.warning(
                "PC-relate requires SSDs and doesn't work with preemptible workers!"
            )
            res = relatedness_resources.create_pc_relate_relatedness_table
            res.check_resource_existance()

            ht = hl.pc_relate(
                joint_qc_mt.GT,
                min_individual_maf=args.min_individual_maf,
                scores_expr=pc_relate_pca_scores.ht()[joint_qc_mt.col_key].scores[
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
            ht.write(res.pc_relate_relatedness_ht.path, overwrite=overwrite)

        if args.finalize_relatedness_ht:
            res = relatedness_resources.finalize_relatedness_ht
            res.check_resource_existance()
            relatedness_args = {
                "parent_child_max_y": args.parent_child_max_ibd0_or_ibs0_over_ibs2,
                "second_degree_sibling_lower_cutoff_slope": args.second_degree_sibling_lower_cutoff_slope,
                "second_degree_sibling_lower_cutoff_intercept": args.second_degree_sibling_lower_cutoff_intercept,
                "second_degree_upper_sibling_lower_cutoff_slope": args.second_degree_upper_sibling_lower_cutoff_slope,
                "second_degree_upper_sibling_lower_cutoff_intercept": args.second_degree_upper_sibling_lower_cutoff_intercept,
                "duplicate_twin_min_kin": args.duplicate_twin_min_kin,
                "second_degree_min_kin": second_degree_kin_cutoff,
                "duplicate_twin_ibd1_min": args.duplicate_twin_ibd1_min,
                "duplicate_twin_ibd1_max": args.duplicate_twin_ibd1_max,
            }
            finalize_relatedness_ht(
                res.relatedness_ht.ht(), joint_qc_meta_ht, rel_method, relatedness_args
            ).write(res.final_relatedness_ht.path, overwrite=overwrite)

        if args.compute_related_samples_to_drop:
            res = relatedness_resources.compute_related_samples_to_drop
            res.check_resource_existance()
            rank_ht, drop_ht = run_compute_related_samples_to_drop(
                res.final_relatedness_ht.ht(),
                joint_qc_meta_ht.semi_join(joint_qc_mt.cols()),
                res.filter_ht.ht(),
                release,
            )
            rank_ht.write(res.sample_rankings_ht.path, overwrite=overwrite)
            drop_ht.write(res.related_samples_to_drop_ht.path, overwrite=overwrite)

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
        help="Minimum cuKING kinship for pair to be included in IBD estimates.",
        default=0.16,
        type=float,
    )
    cuking_ibd.add_argument(
        "--ibd-max-cuking-ibs0",
        help=(
            "Maximum cuKING IBS0 for pair to be included in IBD estimates. Note: This "
            "default was determined from looking at the cuKING Kinship vs. IBS0 plot "
            "for gnomAD v3 + gnomAD v4."
        ),
        default=50,
        type=int,
    )
    cuking_ibd.add_argument(
        "--ibd-max-samples",
        help="Maximum number of samples to include in each IBD run.",
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

    finalize = parser.add_argument_group(
        "Finalize relatedness specific arguments",
        "Arguments specific to creating the final relatedness Table including adding a "
        "'relationship' annotation for each pair. Note: The defaults provided for the "
        "slope and intercept cutoffs were determined from visualization of the "
        "cuking kinship distribution and the IBS0/IBS2 vs. kinship plot.",
    )
    finalize.add_argument(
        "--finalize-relatedness-ht",
        help="Whether to finalize the relatedness HT.",
        action="store_true",
    )
    finalize.add_argument(
        "--finalize-relatedness-method",
        help=(
            "Which relatedness method to use for finalized relatedness Table. Options "
            "are 'cuking' and 'pc_relate'. Default is 'cuking'."
        ),
        default="cuking",
        type=str,
        choices=["cuking", "pc_relate"],
    )
    finalize.add_argument(
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
    finalize.add_argument(
        "--parent-child-max-ibd0-or-ibs0-over-ibs2",
        help=(
            "Maximum value of IBD0 (if '--finalize-relatedness-method' is 'pc_relate') "
            "or IBS0/IBS2 (if '--finalize-relatedness-method' is 'cuking') for a "
            "parent-child pair."
        ),
        default=5.2e-5,
        type=float,
    )
    finalize.add_argument(
        "--second-degree-sibling-lower-cutoff-slope",
        help=(
            "Slope of the line to use as a lower cutoff for second degree relatives "
            "and siblings from parent-child pairs."
        ),
        default=-1.9e-3,
        type=float,
    )
    finalize.add_argument(
        "--second-degree-sibling-lower-cutoff-intercept",
        help=(
            "Intercept of the line to use as a lower cutoff for second degree "
            "relatives and siblings from parent-child pairs."
        ),
        default=5.8e-4,
        type=float,
    )
    finalize.add_argument(
        "--second-degree-upper-sibling-lower-cutoff-slope",
        help=(
            "Slope of the line to use as an upper cutoff for second degree relatives "
            "and a lower cutoff for siblings."
        ),
        default=-1e-2,
        type=float,
    )
    finalize.add_argument(
        "--second-degree-upper-sibling-lower-cutoff-intercept",
        help=(
            "Intercept of the line to use as an upper cutoff for second degree "
            "relatives and a lower cutoff for siblings."
        ),
        default=2.2e-3,
        type=float,
    )
    finalize.add_argument(
        "--duplicate-twin-min-kin",
        help="Minimum kinship for duplicate or twin pairs.",
        default=0.42,
        type=float,
    )
    finalize.add_argument(
        "--duplicate-twin-ibd1-min",
        help=(
            "Minimum IBD1 cutoff for duplicate or twin pairs. Only used when "
            "'--finalize-relatedness-method' is 'pc_relate'. Note: the min is used "
            "because pc_relate can output large negative values in some corner cases."
        ),
        default=-0.15,
        type=float,
    )
    finalize.add_argument(
        "--duplicate-twin-ibd1-max",
        help=(
            "Maximum IBD1 cutoff for duplicate or twin pairs. Only used when "
            "'--finalize-relatedness-method' is 'pc_relate'."
        ),
        default=0.1,
        type=float,
    )

    drop_related_samples = parser.add_argument_group(
        "Compute related samples to drop",
        "Arguments used to determine related samples that should be dropped from "
        "the ancestry PCA or release.",
    )
    drop_related_samples.add_argument(
        "--compute-related-samples-to-drop",
        help="Determine the minimal set of related samples to prune for ancestry PCA.",
        action="store_true",
    )
    drop_related_samples.add_argument(
        "--release",
        help=(
            "Whether to determine related samples to drop for the release based on "
            "outlier filtering of sample QC metrics."
        ),
        action="store_true",
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
        or args.finalize_relatedness_ht
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
