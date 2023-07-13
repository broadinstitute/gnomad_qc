r"""
This is a batch script which adds VRS IDs to a Hail Table by creating a sharded-VCF, running a vrs-annotation script on each shard, and annotating the merged results back onto the original Hail Table.

The vrs-annotation script that generates the VRS IDs needs to be run with Query On Batch. These VRS annotations can be added back to the original Table with either Query On Batch or Spark.(https://hail.is/docs/0.2/cloud/query_on_batch.html#:~:text=Hail%20Query%2Don%2DBatch%20uses,Team%20at%20our%20discussion%20forum.)
usage: python3 vrs_annotation_batch.py \
--billing-project gnomad-annot \
--working-bucket gnomad-qin \
--image us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten_0615_vrs0_8_4 \
--version test_v4.0_exomes \
--prefix v4_vds \
--partitions-for-vcf-export 20 \ # TODO: change
--downsample 0.1 \ # not necessary for a test run
--header-path gs://gnomad/v4.0/annotations/exomes/vrs-header-fix.txt \
--run-vrs \ # only works on Hail BATCH
--annotate-original \ # not necessary for v4 exomes
--overwrite \
--backend-mode batch \
"""

import argparse
import datetime
import errno
import logging
import os
import sys

import ga4gh.vrs
import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.resources.annotations import vrs_annotations as v3_vrs_annotations
from gnomad_qc.v4.resources.annotations import get_vrs as v4_vrs_annotations

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("annotate_vrs_ids")
logger.setLevel(logging.INFO)


def init_job(
    batch,
    name: str = None,
    image: str = None,
    cpu: float = None,
    memory: float = None,
    disk_size: float = None,
):
    """Initialize a hail batch job with some default parameters.

    :param batch: Batch object
    :param name: job label which will show up in the Batch web UI
    :param image: docker image name (eg. "weisburd/image-name@sha256:aa19845da5")
    :param cpu: number of CPUs (between 0.25 to 16)
    :param memory: amount of RAM in Gb (eg. 3.75)
    :param disk_size: amount of disk in Gb (eg. 50)
    :return: new job object
    """
    j = batch.new_job(name=name)
    if image:
        j.image(image)

    if cpu:
        if cpu < 0.25 or cpu > 16:
            raise ValueError(
                f"CPU arg is {cpu}. This is outside the range of 0.25 to 16 CPUs"
            )

        j.cpu(cpu)  # Batch default is 1

    if memory:
        if memory < 0.1 or memory > 60:
            raise ValueError(
                f"Memory arg is {memory}. This is outside the range of 0.1 to 60 Gb"
            )

        j.memory(f"{memory}Gi")  # Batch default is 3.75G

    if disk_size:
        if disk_size < 1 or disk_size > 1000:
            raise ValueError(
                f"Disk size arg is {disk_size}. This is outside the range of 1 to"
                " 1000 Gb"
            )

        j.storage(f"{disk_size}Gi")

    j.command(
        "set -euxo pipefail"
    )  # set bash options for easier debugging and to make command execution more robust

    return j


def init_job_with_gcloud(
    batch,
    name: str = None,
    image: str = None,
    cpu: float = None,
    memory: float = None,
    disk_size: float = None,
    mount: str = None,
):
    """
    Create job and initialize glcoud authentication and gsutil commands.

    Wraps Ben Weisburd's init_job (https://github.com/broadinstitute/tgg_methods/blob/master/tgg/batch/batch_utils.py#L160) with additional gcloud steps.
    Parameters passed through to init_job:
        :param batch: Batch object.
        :param name: Job label which will show up in the Batch web UI.
        :param image: Docker image name (eg. "us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten_0615_vrs0_8_4").
        :param cpu: Number of CPUs (between 0.25 to 16).
        :param memory: Amount of RAM in Gb (eg. 3.75).
        :param disk_size: Amount of disk in Gb (eg. 50).
        :param mount: Name of GCP Bucket to mount using cloudfuse.
        :return: New job object.
    """
    job = init_job(batch, name, image, cpu, memory, disk_size)
    job.command(f"gcloud -q auth activate-service-account --key-file=/gsa-key/key.json")
    job.command(f"curl -sSL broad.io/install-gcs-connector | python3")
    if mount:
        job.cloudfuse(mount, "/local-vrs-mount")
    return job


def main(args):
    """Generate VRS annotations for a Hail Table of variants."""
    # Initialize Hail w/chosen mode (spark or batch for QoB), setting global
    # seed for if user decides to downsample
    hl.init(
        backend=args.backend_mode,
        tmp_dir=args.tmp_dir_hail,
        gcs_requester_pays_configuration="gnomad-vrs",
        default_reference="GRCh38",
        global_seed=args.hail_rand_seed,
    )

    logger.info(f"Hail version as: {hl.version()}")

    # Example schema
    logger.info("Example VRS schema for 2 variants:")
    ht_example = hl.read_table("gs://gnomad-vrs-io-finals/ht-inputs/two-snps-schema.ht")
    ht_example.show()

    working_bucket = args.working_bucket
    version = args.version

    # Prefix to create custom named versions of outputs
    prefix = args.prefix + version

    # TODO: add v4 exomes to this dict
    input_paths_dict = {
        "v3.1.2": public_release("genomes").path,
        "v4.0_exomes": "gs://gnomad-qin/v4_annotations/v4_vds_all_variants.ht",
        "test_v4.0_exomes": (
            f"gs://{working_bucket}/v4_annotations/{prefix}_2_partitions.ht"
        ),
        "test_v3.1.2": public_release("genomes").path,
        "test_v3_1k": (
            "gs://gnomad-vrs-io-finals/ht-inputs/ht-1k-TESTING-ONLY-repartition-10p.ht"
        ),
        "test_v3_10k": (
            "gs://gnomad-vrs-io-finals/ht-inputs/ht-10k-TESTING-ONLY-repartition-50p.ht"
        ),
        "test_v3_100k": "gs://gnomad-vrs-io-finals/ht-inputs/ht-100k-TESTING-ONLY-repartition-100p.ht",
        "test_grch37": "gs://gnomad-vrs-io-finals/working-notebooks/scratch/downsample_and_downpart_with_ychr_grch37.ht",
    }

    output_paths_dict = {
        "v3.1.2": v3_vrs_annotations.path,
        "v4.0": v4_vrs_annotations.path,
        "test_v4.0_exomes": f"gs://{working_bucket}/v4_annotations/{prefix}v4_vds_2_partitions_output.ht",
        "test_v3.1.2": (
            f"gs://gnomad-vrs-io-finals/ht-outputs/{prefix}-Full-ht-release-output.ht"
        ),
        "test_v3_1k": f"gs://{working_bucket}/ht-outputs/{prefix}-Full-ht-1k-output.ht",
        "test_v3_10k": (
            f"gs://{working_bucket}/ht-outputs/{prefix}-Full-ht-10k-output.ht"
        ),
        "test_v3_100k": (
            f"gs://{working_bucket}/ht-outputs/{prefix}-Full-ht-100k-output.ht"
        ),
        "test_grch37": "'gs://gnomad-vrs-io-finals/working-notebooks/scratch/grch37_final_output.ht",
    }

    # Read in Hail Table, partition, and export to sharded VCF
    ht_original = hl.read_table(input_paths_dict[version])
    assembly = get_reference_genome(ht_original.locus).name

    # Option to downsample for testing, if you want to annotate part of a Hail Table but not all of it
    # For this, is important that we have set the Hail random seed!
    if args.downsample < 1.00:
        logger.info("Downsampling Table...")
        ht_original = ht_original.sample(args.downsample)
        ht_original = ht_original.annotate_globals(vrs_downsample=args.downsample)

    if args.run_vrs:
        if args.backend_mode == "spark":
            raise ValueError(
                'Annotation step --run-vrs can only be run with "batch" setting for'
                " backend-mode"
            )

        # Create backend and batch for coming annotation batch jobs
        backend = hb.ServiceBackend(
            billing_project=args.billing_project, bucket=working_bucket
        )

        batch_vrs = hb.Batch(name="vrs-annotation", backend=backend)

        # Check resource existence
        check_resource_existence(
            output_step_resources={
                "--run-vrs": [
                    f"gs://gnomad-vrs-io-finals/ht-outputs/annotated-checkpoint-VRS-{prefix}.ht"
                ],
            },
            overwrite=args.overwrite,
        )
        # Use 'select' to remove all non-key rows - VRS-Allele is added back to
        # original table based on just locus and allele
        ht = ht_original.select()

        logger.info("Table read in and selected")

        # Repartition the Hail Table if requested. In the following step, the Hail
        # Table is exported to a sharded VCF with one shard per partition
        if args.partitions_for_vcf_export:
            if args.partitions_for_vcf_export < ht.n_partitions():
                logger.info("Repartitioning Table by Naive Coalsece")
                ht = ht.naive_coalesce(args.partitions_for_vcf_export)
            elif args.partitions_for_vcf_export > ht.n_partitions():
                logger.info(
                    "Repartitioning Table by Repartition, NOTE this results in a"
                    " shuffle"
                )
                ht = ht.repartition(args.partitions_for_vcf_export)

        logger.info("Exporting sharded VCF")

        hl.export_vcf(
            ht,
            f"gs://{working_bucket}/vrs-temp/shards/shard-{version}.vcf.bgz",
            append_to_header=args.header_path,
            parallel="header_per_shard",
        )

        logger.info(
            "Gathering list of files in"
            f" gs://{working_bucket}/vrs-temp/shards/shard-{version}.vcf.bgz using"
            " Hail's Hadoop method"
        )

        # Create a list of all shards of VCF
        file_dict = hl.utils.hadoop_ls(
            f"gs://{working_bucket}/vrs-temp/shards/shard-{version}.vcf.bgz/part-*.bgz"
        )

        # Create a list of all file names to later annotate in parallel
        file_list = [file_item["path"].split("/")[-1] for file_item in file_dict]

        # Query-On-Batch sometimes produces duplicate shards (same number, same content, different name) for v3.1.2 release Table
        # This block removes duplicates from the list of files to be annotated
        to_exclude = []
        for file_idx in range(len(file_list)):
            file_name = file_list[file_idx].split("-")
            file_num = file_name[1]
            if file_list[file_idx - 1].split("-")[1] == file_num:
                to_exclude.append(file_list[file_idx - 1])

        file_list = sorted(list(set(file_list) - set(to_exclude)))
        logger.info(f"Number of duplicates to be excluded: {len(to_exclude)}")

        logger.info("File list created... getting ready to start Batch Jobs")
        # Define SeqRepo path to be read in, outside of the loop, to avoid reading
        # it in for each job
        seqrepo_path = "/tmp/local-seqrepo"
        if args.seqrepo_mount:
            seqrepo_path = "/local-vrs-mount/seqrepo/2018-11-26/"

        for vcf_name in file_list:
            # Setting up jobs
            new_job = init_job_with_gcloud(
                batch=batch_vrs,
                name=f"VCF_job_{vcf_name.split('.')[0].split('-')[1]}",
                image=args.image,
                mount=args.seqrepo_mount,
                disk_size=args.disk_size,
                memory=args.memory,
            )

            # Script path for annotation is not expected to change for GA4GH if
            # dependencies are installed using pip in the DockerFile
            vrs_script_path = "/usr/local/lib/python3.9/dist-packages/ga4gh/vrs/extras/vcf_annotation.py"

            # Store VCF input and output paths
            vcf_input = f"gs://{working_bucket}/vrs-temp/shards/shard-{version}.vcf.bgz/{vcf_name}"
            vcf_output = f'/temp-vcf-annotated/annotated-{vcf_name.split(".")[0]}.vcf'

            # Perform VRS annotation on vcf_input and store output in /vrs-temp
            new_job.command("mkdir /temp-vcf-annotated/")
            new_job.command(
                f"python3 {vrs_script_path} --vcf_in"
                f" {batch_vrs.read_input(vcf_input)} --vcf_out"
                f" {vcf_output} --seqrepo_root_dir {seqrepo_path} --assembly"
                f" {assembly} --vrs_attributes"
            )

            # Copy annotated shard to its appropriate place in Google Bucket
            new_job.command(
                "gsutil cp"
                f" {vcf_output} gs://{working_bucket}/vrs-temp/annotated-shards/annotated-{version}.vcf/"
            )

        # Execute all jobs in the batch_vrs Batch
        batch_vrs.run()

        logger.info(
            "Batch jobs executed, preparing to read in sharded VCF from prior step."
            " Preparing list of files first using Hail's hadoop_ls method."
        )

        annotated_file_dict = hl.utils.hadoop_ls(
            f"gs://{working_bucket}/vrs-temp/annotated-shards/annotated-{version}.vcf/*.vcf"
        )

        annotated_file_list = [
            annotated_file_item["path"] for annotated_file_item in annotated_file_dict
        ]

        logger.info("File list created. Now reading in annotated shards.")

        # Import all annotated shards
        ht_annotated = hl.import_vcf(
            annotated_file_list,
            reference_genome=assembly,
        ).make_table()
        logger.info("Annotated table constructed")

        vrs_struct = hl.struct(
            VRS_Allele_IDs=ht_annotated.info.VRS_Allele_IDs.split(","),
            VRS_Starts=ht_annotated.info.VRS_Starts.split(",").map(lambda x: hl.int(x)),
            VRS_Ends=ht_annotated.info.VRS_Ends.split(",").map(lambda x: hl.int(x)),
            VRS_States=ht_annotated.info.VRS_States.split(","),
        )

        ht_annotated = ht_annotated.annotate(vrs=vrs_struct)
        ht_annotated = ht_annotated.select(ht_annotated.vrs)

        # Checkpoint (write) resulting annotated table
        ht_annotated = ht_annotated.checkpoint(
            f"gs://gnomad-vrs-io-finals/ht-outputs/annotated-checkpoint-VRS-{prefix}.ht",
            overwrite=args.overwrite,
        )
        logger.info(
            "Annotated Hail Table checkpointed to:"
            f" gs://gnomad-vrs-io-finals/ht-outputs/annotated-checkpoint-VRS-{prefix}.ht"
        )

        # Deleting temporary files to keep clutter down, note that this does not
        # include the annotated HT
        logger.info("Preparing to delete temporary files and sharded VCFs generated")
        delete_temps = hb.Batch(name="delete_temps", backend=backend)
        d1 = init_job_with_gcloud(
            batch=delete_temps,
            name=f"del_files",
            image=args.image,
        )
        d1.command(f"gsutil -m -q rm -r gs://{working_bucket}/vrs-temp/")
        delete_temps.run()

    if args.annotate_original:
        check_resource_existence(
            input_step_resources={
                "--run-vrs": [
                    f"gs://gnomad-vrs-io-finals/ht-outputs/annotated-checkpoint-VRS-{prefix}.ht"
                ],
            },
            output_step_resources={
                "--annotate-original": [output_paths_dict[version]],
            },
            overwrite=args.overwrite,
        )
        # Output final Hail Tables with VRS annotations
        ht_annotated = hl.read_table(
            f"gs://gnomad-vrs-io-finals/ht-outputs/annotated-checkpoint-VRS-{prefix}.ht"
        )

        # Annotate final Hail Table with GA4GH version
        ht_annotated = ht_annotated.annotate_globals(
            ga4gh_vrs_version=ga4gh.vrs.__version__
        )

        if "3.1.2" in version:
            logger.info("Adding VRS IDs to original Table")
            ht_final = ht_original.annotate(
                info=ht_original.info.annotate(
                    vrs=ht_annotated[ht_original.locus, ht_original.alleles].vrs
                )
            )
            logger.info(f"Outputting final table at: {output_paths_dict[version]}")
            ht_final.write(output_paths_dict[version], overwrite=args.overwrite)

        else:
            logger.info(
                "For test datasets, final output is identical to the checkpointed"
                " annotated HT"
            )
            logger.info(f"Outputting final table at: {output_paths_dict[version]}")
            ht_annotated.write(output_paths_dict[version], overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--billing-project", help="Project to bill.", type=str)
    parser.add_argument(
        "--image",
        default=(
            "us-central1-docker.pkg.dev/broad-mpg-gnomad/ga4gh-vrs/marten_0615_vrs0_8_4"
        ),
        help="Image in a GCP Artifact Registry repository.",
        type=str,
    )
    parser.add_argument(
        "--working-bucket",
        help=(
            "Name of GCP Bucket to output intermediate files (sharded VCFs and"
            " checkpointed HTs) to. Final outputs for test versions go to working"
            " bucket, but final outputs ran on the release HT always go in"
            " gs://gnomad-vrs-io-finals/ht-outputs/ ."
        ),
        default="gnomad-vrs-io-finals",
        type=str,
    )
    parser.add_argument(
        "--version", help="Version of HT to read in. ", default="test_v3_1k", type=str
    )
    parser.add_argument(
        "--partitions-for-vcf-export",
        help=(
            "Number of partitions to use when exporting the Table to a sharded VCF"
            " (each partition is exported to a separate VCF). This value determines the"
            " breakdown of jobs that will be ran (the VRS annotation script is ran in"
            " parallel on the VCFs)."
        ),
        default=None,
        type=int,
    )
    parser.add_argument(
        "--prefix",
        help=(
            "Prefix to append to names of all saved temp files and test version output"
            " files."
        ),
        type=str,
        default="",
    )
    parser.add_argument(
        "--header-path",
        help=(
            "Full path of txt file containing lines to append to VCF headers for fields"
            " that maybe be missing when exporting the Table to VCF."
        ),
        type=str,
        default="gs://gnomad-vrs-io-finals/header-fix.txt",
    )
    parser.add_argument(
        "--downsample",
        help="Proportion to which to downsample the original Hail Table input.",
        default=1.00,
        type=float,
    )
    parser.add_argument(
        "--disk-size",
        help="Amount of disk (GB) to allocate to each Hail Batch job.",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--memory",
        help="Amount of memory (GB) to allocate to each Hail Batch job.",
        type=float,
        default=None,
    )
    parser.add_argument(
        "--seqrepo-mount",
        help=(
            "Bucket to mount and read from using Hail Batch's CloudFuse for access to"
            " seqrepo: PLEASE note this DOES have performance implications."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--overwrite",
        help=(
            "Boolean to pass to ht.write(overwrite=_____) determining whether or not to"
            " overwrite existing output for the final Table and checkpointed files."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--run-vrs",
        help=(
            "Pass argument to run VRS annotation on dataset of choice. Specifying"
            " '--run-vrs' also requires setting 'backend-mode' to 'batch', which is the"
            " default."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--annotate-original",
        help="Pass argument to add VRS annotations back to original dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--hail-rand-seed",
        help="Random seed for hail. Default is 5.",
        default=5,
        type=int,
    )
    parser.add_argument(
        "--backend-mode",
        help="Mode in which to run Hail - either 'spark' or 'batch' (for QoB)",
        choices=["spark", "batch"],
        type=str,
        default="batch",
    )
    parser.add_argument(
        "--tmp-dir-hail",
        help="Directory for temporary files to set when initializing Hail.",
        default=None,
        type=str,
    )

    args = parser.parse_args()

    main(args)
