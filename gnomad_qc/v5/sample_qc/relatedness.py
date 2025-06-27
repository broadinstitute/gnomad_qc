"""Script to compute relatedness estimates among pairs of samples in the joint gnomAD v4 (exomes + genomes) + AoU callset."""

import argparse
import logging
import textwrap
from typing import Dict, Optional, Tuple

import hail as hl
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    compute_related_samples_to_drop,
    get_slope_int_relationship_expr,
)
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import project_meta
from gnomad_qc.v5.resources.sample_qc import (
    get_cuking_input_path,
    get_cuking_output_path,
    get_joint_qc,
    related_samples_to_drop,
    relatedness,
    sample_rankings,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def print_cuking_command(
    cuking_input_path: str,
    cuking_output_path: str,
    min_emission_kinship: float = 0.5,
    cuking_split_factor: int = 4,
    location: str = "us-central1",
    machine_type: str = "a2-highgpu-2g",
    accelerator_count: int = 2,
    accelerator_type: str = "nvidia-tesla-a100",
    requester_pays_project: str = "terra-vpc-sc-93ccd8d2",
) -> None:
    """
    Print the command to submit a dsub job for running cuKING.

    :param cuking_input_path: Path to the cuKING input Parquet files.
    :param cuking_output_path: Path to the cuKING output Parquet files.
    :param min_emission_kinship: Minimum kinship threshold for emitting a pair of
        samples in the relatedness output.
    :param cuking_split_factor: Split factor to use for splitting the full relatedness
        matrix table into equally sized submatrices that are computed independently to
        parallelize the relatedness computation and decrease the memory requirements.
        For example, to halve memory requirements, the full matrix can be split into
        equally sized submatrices (i.e. a 'split factor' of 4). Only the
        'upper triangular' submatrices need to be evaluated due to the symmetry of the
        relatedness matrix, leading to 10 shards. Default is 4.
    :param location: Location to run the dsub job. Default is "us-central1".
    :param machine_type: Machine type to use for the dsub job. Default is "a2-highgpu-2g".
        The default for v4 was a2-highgpu-1g, but this was not sufficient for v5.
    :param accelerator_count: Number of accelerators to use for the dsub job. Default is 2.
    :param accelerator_type: Type of accelerator to use for the dsub job.
        Default is "nvidia-tesla-a100".
    :param requester_pays_project: Project ID to use for requester pays buckets.
        Default is "terra-vpc-sc-93ccd8d2".
    :return:
    """
    logger.info(
        "Printing a command that can be used to submit a dsub job to run "
        "cuKING on the files created by --prepare-cuking-inputs."
    )
    # TOTAL_SHARDS calculates the total number of shards using the formula k(k+1)/2.
    cuking_command = f"""\
        SPLIT_FACTOR={cuking_split_factor} \\
        TOTAL_SHARDS=$((SPLIT_FACTOR * (SPLIT_FACTOR + 1) / 2)) \\
        for SHARD_INDEX in $(seq 0 $((TOTAL_SHARDS - 1))); do
        cuKING_dsub \\
        --location={location} \\
        --machine-type={machine_type} \\
        --accelerator-count={accelerator_count} \\
        --accelerator-type={accelerator_type} \\
        --command="cuking \\
        --input_uri="{cuking_input_path}" \\
        --output_uri="{cuking_output_path}/out_split_${{SHARD_INDEX}}.parquet" \\
        --requester_pays_project={requester_pays_project} \\
        --kin_threshold={min_emission_kinship} \\
        --split_factor=${{SPLIT_FACTOR}} \\
        --shard_index=${{SHARD_INDEX}}"
        done
        """
    print(textwrap.dedent(cuking_command))


def convert_cuking_output_to_ht(cuking_output_path: str) -> hl.Table:
    """
    Convert cuKING output Parquet files to a Hail Table.

    .. note ::

        We cannot import this function from the cuKING script, which is why it is duplicated here.

    :param cuking_output_path: Path to the cuKING output Parquet files.
    :return: Hail Table containing the relatedness estimates.
    """
    spark = hl.utils.java.Env.spark_session()
    # The cuking output is sharded and stored in a single directory
    df = spark.read.parquet(f"{cuking_output_path}/*.parquet")
    ht = hl.Table.from_spark(df)
    ht = ht.key_by(i=hl.struct(s=ht.i), j=hl.struct(s=ht.j))
    return ht


def finalize_relatedness_ht(
    ht: hl.Table,
    meta_ht: hl.Table,
    relatedness_args: Dict[str, float],
) -> hl.Table:
    """
    Create the finalized relatedness Table including adding a 'relationship' annotation for each pair.

    The `relatedness_args` dictionary should have the following keys:
        - 'second_degree_min_kin': Minimum kinship threshold for filtering a pair of
          samples with a second degree relationship when filtering related individuals.
          Default is 0.08838835. Bycroft et al. (2018) calculates a theoretical kinship
          of 0.08838835 for a second degree relationship cutoff. This cutoff should be
          determined by evaluation of the kinship distribution.
        - 'parent_child_max_ibd0_or_ibs0_over_ibs2': Maximum value of IBD0 (if
          `relatedness_method` is 'pc_relate') or IBS0/IBS2 (if `relatedness_method` is
          'cuking') for a parent-child pair.
        - 'second_degree_sibling_lower_cutoff_slope': Slope of the line to use as a
          lower cutoff for second degree relatives and siblings from parent-child pairs.
        - 'second_degree_sibling_lower_cutoff_intercept': Intercept of the line to use
          as a lower cutoff for second degree relatives and siblings from parent-child
          pairs.
        - 'second_degree_upper_sibling_lower_cutoff_slope': Slope of the line to use as
          an upper cutoff for second degree relatives and a lower cutoff for siblings.
        - 'second_degree_upper_sibling_lower_cutoff_intercept': Intercept of the line to
          use as an upper cutoff for second degree relatives and a lower cutoff for
          siblings.
        - 'duplicate_twin_min_kin': Minimum kinship for duplicate or twin pairs.


    The following annotations are added to `ht`:
        - 'relationship': Relationship annotation for the pair. Returned by
          `get_slope_int_relationship_expr`.
        - 'gnomad_v4_duplicate': Whether the sample is a duplicate of a sample found in
          the gnomAD v4 genomes.
        - 'gnomad_v4_release_duplicate': Whether the sample is a duplicate of a sample
          found in the gnomAD v4 release genomes.

    :param ht: Input relatedness Table.
    :param meta_ht: Input metadata Table. Used to add v4 overlap annotations.
    :param relatedness_args: Dictionary of arguments to be passed to
        `get_slope_int_relationship_expr`.
    :return: Finalized relatedness Table
    """
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
        relationship_inference_method="cuking",
        relationship_cutoffs=hl.struct(**relatedness_args),
    )
    ht = ht.key_by(
        i=ht.i.annotate(data_type=meta_ht[ht.i.s].data_type),
        j=ht.j.annotate(data_type=meta_ht[ht.j.s].data_type),
    )
    gnomad_v4_duplicate_expr = (ht.relationship == DUPLICATE_OR_TWINS) & (
        hl.set({ht.i.data_type, ht.j.data_type}) == hl.set(["exomes", "genomes"])
    )
    gnomad_v4_release_expr = hl.coalesce(
        meta_ht[ht.i.s].release | meta_ht[ht.j.s].release,
        False,
    )
    ht = ht.annotate(
        gnomad_v4_duplicate=gnomad_v4_duplicate_expr,
        gnomad_v4_release_duplicate=gnomad_v4_duplicate_expr & gnomad_v4_release_expr,
    )

    return ht


def compute_rank_ht(ht: hl.Table) -> hl.Table:
    """
    Add a rank to each sample for use when breaking maximal independent set ties.

    Favor AoU samples, then v4 release samples, then genomes, then higher mean depth
    ('chr20_mean_dp' from gnomad and 'mean_coverage' from AoU).

    :param ht: Table to add rank to.
    :return: Table containing sample ID and rank.
    """
    ht = ht.select(
        "mean_depth",
        is_genome=hl.is_defined(ht.data_type) & ht.data_type == "genomes",
        in_aou=hl.is_defined(ht.project) & ht.project == "aou",
        in_v4_release=hl.is_defined(ht.release) & ht.release,
    )
    rank_order = []
    ht_select = ["rank"]
    rank_order.extend(
        [
            hl.desc(ht.in_aou),
            hl.desc(ht.in_v4_release),
            hl.desc(ht.is_genome),
            hl.desc(ht.mean_depth),
        ]
    )
    ht = ht.order_by(*rank_order).add_index(name="rank")
    ht = ht.key_by(ht.s)

    return ht.select(*ht_select)


def run_compute_related_samples_to_drop(
    ht: hl.Table,
    meta_ht: hl.Table,
) -> Tuple[hl.Table, hl.Table]:
    """
    Determine the minimal set of related samples to prune for ancestry PCA.

    Runs `compute_related_samples_to_drop` from gnomad_methods after computing the
    sample rankings using `compute_rank_ht`.

    :param ht: Input relatedness Table.
    :param meta_ht: Metadata Table with 'project', 'data_type', 'release', and
        'mean_depth' annotations to be used in ranking and filtering of related
        individuals.
    :return: Table with sample rank and a Table with related samples to drop for PCA.
    """
    # Compute_related_samples_to_drop uses a rank Table as a tiebreaker when
    # pruning samples.
    rank_ht = compute_rank_ht(meta_ht)
    rank_ht = rank_ht.checkpoint(new_temp_file("rank", extension="ht"))

    second_degree_min_kin = hl.eval(ht.relationship_cutoffs.second_degree_min_kin)
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)

    samples_to_drop_ht = compute_related_samples_to_drop(
        ht,
        rank_ht,
        second_degree_min_kin,
    )
    samples_to_drop_ht = samples_to_drop_ht.annotate_globals(
        second_degree_min_kin=second_degree_min_kin,
    )

    return rank_ht, samples_to_drop_ht


def main(args):
    """Compute relatedness estimates among pairs of samples in the callset."""
    test = args.test
    overwrite = args.overwrite
    release = args.release
    min_emission_kinship = args.min_emission_kinship
    second_degree_min_kin = args.second_degree_min_kin
    joint_meta = project_meta.ht()

    if args.print_cuking_command:
        print_cuking_command(
            get_cuking_input_path(test=test),
            get_cuking_output_path(test=test),
            min_emission_kinship,
            args.cuking_split_factor,
        )
        return

    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/relatedness.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")

    try:
        if args.prepare_cuking_inputs:
            from gnomad_qc.v4.sample_qc.cuKING.mt_to_cuking_inputs import (
                mt_to_cuking_inputs,
            )

            # NOTE: This step is going to be run in a notebook to avoid
            # import error associated with the above import.
            # We also cannot run `mt_to_cuking_inputs.py` directly because
            # Hail needs to be initialized with a temporary directory
            # to avoid memory errors.
            joint_qc_mt = get_joint_qc(test=test).mt()
            check_resource_existence(
                output_step_resources={"cuking_input_parquet": get_cuking_input_path()}
            )
            logger.info(
                "Converting joint dense QC MatrixTable to a Parquet format suitable "
                "for cuKING..."
            )
            if test:
                logger.info("Filtering MT to the first 2 partitions for testing...")
                joint_qc_mt = joint_qc_mt._filter_partitions(range(2))

            mt_to_cuking_inputs(
                mt=joint_qc_mt,
                parquet_uri=get_cuking_input_path(test=test),
                overwrite=overwrite,
            )
        if args.create_cuking_relatedness_table:
            logger.info("Converting cuKING outputs to Hail Table...")
            check_resource_existence(
                output_step_resources={
                    "relatedness_ht": relatedness(test=test, raw=True).path
                }
            )
            ht = convert_cuking_output_to_ht(get_cuking_output_path(test=test))
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(relatedness(test=test, raw=True).path, overwrite=overwrite)

        if args.finalize_relatedness_ht:
            check_resource_existence(
                output_step_resources={
                    "final_relatedness_ht": relatedness(test=test).path
                }
            )
            relatedness_args = {
                "parent_child_max_y": args.parent_child_max_ibs0_over_ibs2,
                "second_degree_sibling_lower_cutoff_slope": (
                    args.second_degree_sibling_lower_cutoff_slope
                ),
                "second_degree_sibling_lower_cutoff_intercept": (
                    args.second_degree_sibling_lower_cutoff_intercept
                ),
                "second_degree_upper_sibling_lower_cutoff_slope": (
                    args.second_degree_upper_sibling_lower_cutoff_slope
                ),
                "second_degree_upper_sibling_lower_cutoff_intercept": (
                    args.second_degree_upper_sibling_lower_cutoff_intercept
                ),
                "duplicate_twin_min_kin": args.duplicate_twin_min_kin,
                "second_degree_min_kin": second_degree_min_kin,
            }
            finalize_relatedness_ht(
                relatedness(test=test, raw=True).ht(), joint_meta, relatedness_args
            ).write(relatedness(test=test).path, overwrite=overwrite)

        if args.compute_related_samples_to_drop:
            check_resource_existence(
                output_resources={
                    "sample_rankings_ht": sample_rankings(test=test),
                    "related_samples_to_drop_ht": related_samples_to_drop(test=test),
                }
            )
            rank_ht, drop_ht = run_compute_related_samples_to_drop(
                relatedness(test=test).ht(),
                joint_meta.select_globals().semi_join(joint_qc_mt.cols()),
            )
            rank_ht.write(sample_rankings(test=test).path, overwrite=overwrite)
            drop_ht.write(related_samples_to_drop(test=test).path, overwrite=overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("relatedness", environment=args.environment))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Filter to the first 2 partitions for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--environment",
        help="Environment where script will run.",
        default="rwb",
        type=str,
    )
    parser.add_argument(
        "--min-emission-kinship",
        help=(
            "Minimum kinship threshold for emitting a pair of samples in the "
            "relatedness output."
        ),
        default=0.05,
        type=float,
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
        help="Create a Hail Table from the cuKING output.",
        action="store_true",
    )
    cuking_args.add_argument(
        "--relatedness-n-partitions",
        help="Number of partitions to use for the relatedness Table.",
        default=100,
        type=int,
    )

    # Finalize relatedness arguments
    finalize_args = parser.add_argument_group(
        "Finalize relatedness arguments",
        "Arguments for finalizing the relatedness table with relationship annotation.",
    )
    finalize_args.add_argument(
        "--finalize-relatedness-ht",
        help="Finalize the relatedness HT by adding relationship annotations.",
        action="store_true",
    )
    finalize_args.add_argument(
        "--parent-child-max-ibs0-over-ibs2",
        help="Maximum value of IBS0/IBS2 for a parent-child pair.",
        default=5.2e-5,
        type=float,
    )
    finalize_args.add_argument(
        "--second-degree-sibling-lower-cutoff-slope",
        help="Slope of the line to use as a lower cutoff for second degree relatives and siblings from parent-child pairs.",
        default=-1.9e-3,
        type=float,
    )
    finalize_args.add_argument(
        "--second-degree-sibling-lower-cutoff-intercept",
        help="Intercept of the line to use as a lower cutoff for second degree relatives and siblings from parent-child pairs.",
        default=5.8e-4,
        type=float,
    )
    finalize_args.add_argument(
        "--second-degree-upper-sibling-lower-cutoff-slope",
        help="Slope of the line to use as an upper cutoff for second degree relatives and a lower cutoff for siblings.",
        default=-1e-2,
        type=float,
    )
    finalize_args.add_argument(
        "--second-degree-upper-sibling-lower-cutoff-intercept",
        help="Intercept of the line to use as an upper cutoff for second degree relatives and a lower cutoff for siblings.",
        default=2.2e-3,
        type=float,
    )
    finalize_args.add_argument(
        "--duplicate-twin-min-kin",
        help="Minimum kinship for duplicate or twin pairs.",
        default=0.42,
        type=float,
    )
    finalize_args.add_argument(
        "--second-degree-min-kin",
        help="Minimum kinship threshold for filtering a pair of samples with a second degree relationship when filtering related individuals.",
        default=0.08838835,
        type=float,
    )
    drop_related_samples = parser.add_argument_group(
        "Compute related samples to drop",
        "Arguments used to determine related samples that should be dropped from "
        "the ancestry PCA or release.",
    )
    drop_related_samples.add_argument(
        "--compute-related-samples-to-drop",
        help=(
            "Determine the minimal set of related samples to prune for ancestry PCA or "
            "release if '--release' is used."
        ),
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
