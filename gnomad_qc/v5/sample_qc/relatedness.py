"""Script to compute relatedness estimates among pairs of samples in the joint gnomAD v4 (exomes + genomes) + AoU callset."""

import argparse
import logging
import textwrap
from typing import Dict, Optional, Tuple, Union

import hail as hl
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    compute_related_samples_to_drop,
    get_slope_int_relationship_expr,
)
from gnomad.utils.file_utils import check_file_exists_raise_error
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import consent_samples_to_drop, project_meta
from gnomad_qc.v5.resources.sample_qc import (
    get_cuking_input_path,
    get_cuking_output_path,
    get_joint_qc,
    hard_filtered_samples,
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
        - 'parent_child_max_ibd0_or_ibs0_over_ibs2': Maximum IBS0/IBS2 for a
          parent-child pair.
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
        - 'gnomad_exomes_duplicate': Whether the sample pair is an AoU gnomAD exomes
          duplicate.
        - 'gnomad_genomes_duplicate': Whether the sample pair is an AoU gnomAD genomes
          duplicate.
        - 'released_gnomad_exomes_duplicate': Whether the sample pair contains a released
          gnomAD exomes duplicate.
        - 'released_gnomad_genomes_duplicate': Whether the sample pair contains a released
          gnomAD genomes duplicate.

    :param ht: Input relatedness Table.
    :param meta_ht: Input metadata Table. Used to add v4 overlap annotations.
    :param relatedness_args: Dictionary of arguments to be passed to
        `get_slope_int_relationship_expr`.
    :return: Finalized relatedness Table
    """
    y_expr = ht.ibs0 / ht.ibs2
    ibd1_expr = None

    ht = ht.annotate(
        ibs0_2=ht.ibs0 / ht.ibs2,
        relationship=get_slope_int_relationship_expr(
            kin_expr=ht.kin,
            y_expr=y_expr,
            **relatedness_args,
            ibd1_expr=ibd1_expr,
        ),
    )
    # The function `get_slope_int_relationship_expr` expects parent_child_max_y for
    # flexibility but revert back to parent_child_max_ibs0_over_ibs2 for global clarity.
    relatedness_args["parent_child_max_ibs0_over_ibs2"] = relatedness_args.pop(
        "parent_child_max_y"
    )
    ht = ht.annotate_globals(
        relationship_inference_method="cuking",
        relationship_cutoffs=hl.struct(**relatedness_args),
    )

    # Mark AoU samples as duplicates of gnomAD samples.
    ht = ht.key_by(
        i=ht.i.annotate(
            data_type=meta_ht[ht.i.s].data_type,
            project=meta_ht[ht.i.s].project,
            release=meta_ht[ht.i.s].release,
        ),
        j=ht.j.annotate(
            data_type=meta_ht[ht.j.s].data_type,
            project=meta_ht[ht.j.s].project,
            release=meta_ht[ht.j.s].release,
        ),
    )

    is_aou_gnomad_duplicate = (ht.relationship == DUPLICATE_OR_TWINS) & (
        hl.set({ht.i.project, ht.j.project}) == hl.set(["aou", "gnomad"])
    )

    # Get the gnomAD sample's data type and release status (if one exists in the pair).
    gnomad_data_type = hl.if_else(
        ht.i.project == "gnomad",
        ht.i.data_type,
        hl.if_else(ht.j.project == "gnomad", ht.j.data_type, hl.missing(hl.tstr)),
    )

    gnomad_release_status = hl.if_else(
        ht.i.project == "gnomad",
        ht.i.release,
        hl.if_else(ht.j.project == "gnomad", ht.j.release, hl.missing(hl.tbool)),
    )

    ht = ht.annotate(
        is_aou_gnomad_duplicate=is_aou_gnomad_duplicate,
        gnomad_data_type=gnomad_data_type,
        gnomad_exomes_aou_duplicate=is_aou_gnomad_duplicate
        & (gnomad_data_type == "exomes"),
        gnomad_genomes_aou_duplicate=is_aou_gnomad_duplicate
        & (gnomad_data_type == "genomes"),
        released_gnomad_exomes_aou_duplicate=is_aou_gnomad_duplicate
        & (gnomad_data_type == "exomes")
        & gnomad_release_status,
        released_gnomad_genomes_duplicate=is_aou_gnomad_duplicate
        & (gnomad_data_type == "genomes")
        & gnomad_release_status,
    )

    return ht


def get_consent_samples_to_drop(write_resource: bool = False) -> Union[hl.Table, None]:
    """
    Create additional samples to drop for release HailExpression resource.

    Additional samples to drop are genomes that are no longer consented to be in gnomAD (n = 897).

    .. note ::

        We are adding this function in this script (rather than earlier in the pipeline) due to the timing of confirmation
        of the consent drop.

    :param write_resource: Whether to write the consent drop samples to a HailExpression resource.
    :return: Either Table with consent drop samples or None if writing resource.
    """
    meta_ht = v4_meta(data_type="genomes").ht()
    # RP-1061: 776 samples.
    # RP-1411: 121 samples.
    consent_drop_ht = (
        meta_ht.filter(
            (meta_ht.project_meta.research_project_key == "RP-1061")
            | (meta_ht.project_meta.research_project_key == "RP-1411")
        )
        .select_globals()
        .select()
    )
    if write_resource:
        consent_drop_s = consent_drop_ht.aggregate(
            hl.agg.collect_as_set(consent_drop_ht.s)
        )
        hl.experimental.write_expression(
            consent_drop_s,
            consent_samples_to_drop.path,
        )
    return consent_drop_ht


def compute_rank_ht(ht: hl.Table, filter_ht: Optional[hl.Table] = None) -> hl.Table:
    """
    Add a rank to each sample for use when breaking maximal independent set ties.

    Favor AoU samples, then v4 release samples, then genomes, then higher mean depth
    ('chr20_mean_dp' from gnomad and 'mean_coverage' from AoU).

    If `filter_ht` is provided, rank based on filtering 'outlier_filtered' annotation
    first.

    :param ht: Table to add rank to.
    :param filter_ht: Optional Table with 'outlier_filtered' annotation to be used in
        ranking.
    :return: Table containing sample ID and rank.
    """
    ht = ht.select(
        "mean_depth",
        in_aou=ht.project == "aou",
        is_genome=ht.data_type == "genomes",
        in_v4_release=ht.release,
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
    release: bool = False,
    filter_ht: Optional[hl.Table] = None,
) -> Tuple[hl.Table, hl.Table]:
    """
    Determine the minimal set of related samples to prune for genetic ancestry PCA.

    Runs `compute_related_samples_to_drop` from gnomad_methods after computing the
    sample rankings using `compute_rank_ht`.

    .. note ::
        After some discussion, we have decided to force v4 retention (see slack thread: https://atgu.slack.com/archives/CRA2TKTV0/p1757336429532189)
        and prioritize AoU samples in the ranking because this retains 1,608 more AoU samples.

    :param ht: Input relatedness Table.
    :param meta_ht: Metadata Table with 'project', 'data_type', 'release', and
        'mean_depth' annotations to be used in ranking and filtering of related
        individuals.
    :param release: Whether to determine related samples to drop for the release based
        on outlier filtering of sample QC metrics. Also drops non-released v4 samples and consent drop samples.
        `filter_ht` must be supplied. Default is False.
    :param filter_ht: Optional Table with outlier filtering of sample QC metrics to
        use if `release` is True. Default is None.
    :return: Table with sample rank and a Table with related samples to drop for PCA.
    """
    if release and filter_ht is None:
        raise ValueError("'filter_ht' must be supplied when 'release' is True!")

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

    second_degree_min_kin = hl.eval(ht.relationship_cutoffs.second_degree_min_kin)
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)

    # Create consent drop HailExpression resource if it does not exist.
    if not check_file_exists_raise_error(
        fname=consent_samples_to_drop.path,
        error_if_not_exists=False,
    ):
        get_consent_samples_to_drop(write_resource=True)
    consent_drop_s = hl.experimental.read_expression(consent_samples_to_drop.path)

    # Get set of 806,296 v4 release samples to keep.
    v4_release_ht = meta_ht.filter((meta_ht.project == "gnomad") & meta_ht.release)
    v4_release_s = hl.literal(
        v4_release_ht.aggregate(hl.agg.collect_as_set(v4_release_ht.s))
    )
    v4_release_keep = v4_release_s.difference(consent_drop_s)
    if hl.eval(hl.len(v4_release_keep) != 806296):
        raise ValueError(
            "Expected 806,296 v4 release samples to keep, but got {}.".format(
                hl.len(v4_release_keep)
            )
        )

    samples_to_drop_ht = compute_related_samples_to_drop(
        ht,
        rank_ht,
        second_degree_min_kin,
        filtered_samples=filtered_samples,
        # Force maximal independent set to keep v4 release samples minus consent
        # drop samples.
        keep_samples=v4_release_keep,
        keep_samples_when_related=True,
    )
    samples_to_drop_ht = samples_to_drop_ht.annotate_globals(
        second_degree_min_kin=second_degree_min_kin,
    )

    if release:
        # Add samples to drop to the drop HT.
        v4_unreleased_ht = (
            meta_ht.filter((meta_ht.project == "gnomad") & ~meta_ht.release)
            .select_globals()
            .select()
        )
        consent_drop_ht = get_consent_samples_to_drop(write_resource=False)
        samples_to_drop_ht = samples_to_drop_ht.union(
            v4_unreleased_ht,
            consent_drop_ht,
        )

    return rank_ht, samples_to_drop_ht


def main(args):
    """Compute relatedness estimates among pairs of samples in the callset."""
    test = args.test
    overwrite = args.overwrite
    min_emission_kinship = args.min_emission_kinship
    second_degree_min_kin = args.second_degree_min_kin
    release = args.release

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

    joint_meta = project_meta.ht()
    joint_qc_mt = get_joint_qc().mt()

    if test:
        logger.info("Filtering MT to the first 2 partitions for testing...")
        joint_qc_mt = joint_qc_mt._filter_partitions(range(2))

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
            check_resource_existence(
                output_step_resources={"cuking_input_parquet": get_cuking_input_path()}
            )
            logger.info(
                "Converting joint dense QC MatrixTable to a Parquet format suitable "
                "for cuKING..."
            )
            mt_to_cuking_inputs(
                mt=joint_qc_mt,
                parquet_uri=get_cuking_input_path(test=test),
                overwrite=overwrite,
            )
        if args.create_cuking_relatedness_table:
            logger.info("Converting cuKING outputs to Hail Table...")
            check_resource_existence(
                output_step_resources={
                    "relatedness_ht": [relatedness(test=test, raw=True).path]
                },
                overwrite=overwrite,
            )
            ht = convert_cuking_output_to_ht(get_cuking_output_path(test=test))
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(relatedness(test=test, raw=True).path, overwrite=overwrite)

        if args.finalize_relatedness_ht:
            logger.info("Finalizing relatedness HT...")
            check_resource_existence(
                output_step_resources={
                    "final_relatedness_ht": [relatedness(test=test).path]
                },
                overwrite=overwrite,
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
            logger.info("Computing the sample rankings and related samples to drop...")
            check_resource_existence(
                output_step_resources={
                    "sample_rankings_ht": [sample_rankings(test=test).path],
                    "related_samples_to_drop_ht": [
                        (related_samples_to_drop(test=test, release=release).path)
                    ],
                },
                overwrite=overwrite,
            )
            filter_ht = None
            if release:
                # TODO: Replace this with outlier filtering HT when that exists.
                filter_ht_path = hard_filtered_samples.path
                check_resource_existence(
                    input_step_resources={
                        "filter_ht": [filter_ht_path],
                    }
                )
                filter_ht = hl.read_table(filter_ht_path)
                # TODO: Remove this when have outlier table.
                filter_ht = filter_ht.annotate(outlier_filtered=True).select(
                    "outlier_filtered"
                )

            rank_ht, drop_ht = run_compute_related_samples_to_drop(
                relatedness(test=test).ht(),
                joint_meta.select_globals().semi_join(joint_qc_mt.cols()),
                release=release,
                filter_ht=filter_ht,
            )
            rank_ht.write(sample_rankings(test=test).path, overwrite=overwrite)
            drop_ht.write(
                related_samples_to_drop(test=test, release=release).path,
                overwrite=overwrite,
            )
            # Export the related samples to drop Table to a TSV file for SV team.
            drop_ht.export(
                related_samples_to_drop(test=test, release=release).path[:-3] + ".tsv"
            )

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
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
