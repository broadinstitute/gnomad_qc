"""Script to compute relatedness estimates among pairs of samples in the joint gnomAD v4 (exomes + genomes) + AoU callset."""

import argparse
import logging
import textwrap

import hail as hl

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.sample_qc import (
    get_cuking_input_path,
    get_cuking_output_path,
    get_joint_qc,
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


def main(args):
    """Compute relatedness estimates among pairs of samples in the callset."""
    test = args.test
    overwrite = args.overwrite
    min_emission_kinship = args.min_emission_kinship

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

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("relatedness"))


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

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
