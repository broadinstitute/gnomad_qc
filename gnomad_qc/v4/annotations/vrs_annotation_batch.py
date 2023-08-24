r"""
This is a batch script which adds VRS IDs to a Hail Table by creating sharded VCFs, running a vrs-annotation script on each shard, and merging the results into the original Hail Table.

Advise to run from a backend (eg. hailctl batch submit) to avoid losing progress
 in case of disconnection.

Example command to use `hailctl batch submit`:
hailctl batch submit \
--image-name us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vrs084 \
gnomad_qc/gnomad_qc/v4/annotations/vrs_annotation_batch.py \
-- \ # TODO: remove this when hailctl batch submit bug is fixed in hail 0.2.121
--billing-project gnomad-annot \
--working-bucket gnomad-tmp-4day \
--image us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vrs084 \
--header-path gs://gnomad/v4.0/annotations/exomes/vrs-header-fix.txt \
--run-vrs \
--annotate-original \
--overwrite \
--backend-mode batch \
--test
"""

import argparse
import logging

import hail as hl
import hailtop.batch as hb
from gnomad.utils.reference_genome import get_reference_genome

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.annotations import get_vep as v4_input_ht
from gnomad_qc.v4.resources.annotations import get_vrs as v4_vrs_annotations

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("annotate_vrs_ids")
logger.setLevel(logging.INFO)

# Define the version of ga4gh.vrs code this was run on,
# as present in the Dockerfile
# Please change this when the Dockerfile is updated
VRS_VERSION = "0.8.4"


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
    # Retry gcloud auth and gsutil commands to avoid transient failures
    job.command("""
    retry() {
              "$@" ||
                  (sleep 2 && "$@") ||
                  (sleep 5 && "$@");
    }
    retry gcloud -q auth activate-service-account --key-file=/gsa-key/key.json

    curl_and_python() {
        curl -sSL broad.io/install-gcs-connector | python3
    }
    retry curl_and_python
    """)
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
        gcs_requester_pays_configuration="broad-mpg-gnomad",
        default_reference="GRCh38",
        global_seed=args.hail_rand_seed,
        regions=["us-central1"],
    )

    logger.info(f"Hail version as: {hl.version()}")

    # Example schema
    logger.info("Example VRS schema for 2 variants:")
    ht_example = hl.read_table("gs://gnomad-vrs-io-finals/ht-inputs/two-snps-schema.ht")
    ht_example.show()

    working_bucket = args.working_bucket

    input_path = v4_input_ht().path

    # Read in Hail Table, partition, and export to sharded VCF
    ht_original = hl.read_table(input_path)
    assembly = get_reference_genome(ht_original.locus).name

    # Option to downsample for testing, if you want to annotate part
    # of a Hail Table but not all of it
    # For this, is important that we have set the Hail random seed!
    if args.test:
        logger.info("Filtering to 200 partitions for testing...")
        ht_original = ht_original._filter_partitions(range(200))

    elif args.downsample < 1.00:
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
            billing_project=args.billing_project,
            remote_tmpdir=args.tmp_dir_hail,
        )

        batch_vrs = hb.Batch(name="vrs-annotation", backend=backend)

        # Check resource existence
        check_resource_existence(
            output_step_resources={
                "--run-vrs": [
                    v4_vrs_annotations(test=args.test).path,
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
            f"gs://{working_bucket}/vrs-temp/shards/shard.vcf.bgz",
            append_to_header=args.header_path,
            parallel="header_per_shard",
        )

        logger.info(
            "Gathering list of files in"
            f" gs://{working_bucket}/vrs-temp/shards/shard.vcf.bgz using"
            " Hail's Hadoop method"
        )

        # Create a list of all shards of VCF
        # Note: This step requires using Hail 0.2.120 or later, for
        # faster listing of files in a directory with hadoop_ls.
        file_dict = hl.utils.hadoop_ls(
            f"gs://{working_bucket}/vrs-temp/shards/shard.vcf.bgz/part-*.bgz"
        )

        # Create a list of all file names to later annotate in parallel
        file_list = [file_item["path"].split("/")[-1] for file_item in file_dict]

        # Query-On-Batch sometimes produces duplicate shards
        # (same number, same content, different name)
        # This block removes duplicates from the list of files to be annotated
        to_exclude = []
        for file_idx in range(len(file_list)):
            file_name = file_list[file_idx].split("-")
            file_num = file_name[1]
            if len(file_list) > 1 and file_list[file_idx - 1].split("-")[1] == file_num:
                to_exclude.append(file_list[file_idx - 1])

        file_list = sorted(list(set(file_list) - set(to_exclude)))
        logger.info(f"Number of duplicates to be excluded: {len(to_exclude)}")

        logger.info("File list created... getting ready to start Batch Jobs")
        # Define SeqRepo path to be read in, outside the loop, to avoid reading
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
            vcf_input = (
                f"gs://{working_bucket}/vrs-temp/shards/shard.vcf.bgz/{vcf_name}"
            )
            vcf_output = f'/temp-vcf-annotated/annotated-{vcf_name.split(".")[0]}.vcf'

            # Perform VRS annotation on vcf_input and store output in /vrs-temp
            new_job.command("mkdir /temp-vcf-annotated/")
            new_job.command(
                f"python3 {vrs_script_path} --vcf_in"
                f" {batch_vrs.read_input(vcf_input)} --vcf_out"
                f" {vcf_output} --seqrepo_root_dir {seqrepo_path} --assembly"
                f" {assembly} --vrs_attributes"
            )

            # bgzip the annotated VCF
            new_job.command(f"bgzip {vcf_output}")

            # Copy annotated shard to its appropriate place in Google Bucket
            new_job.command(
                "gsutil cp"
                f" {vcf_output}.gz"
                f" gs://{working_bucket}/vrs-temp/annotated-shards/annotated.vcf/"
            )

        # Execute all jobs in the batch_vrs Batch
        batch_vrs.run()

        logger.info(
            "Batch jobs executed, preparing to read in sharded VCF from prior step."
            " Preparing list of files first using Hail's hadoop_ls method."
        )

        annotated_file_dict = hl.utils.hadoop_ls(
            f"gs://{working_bucket}/vrs-temp/annotated-shards/annotated.vcf/*"
        )

        annotated_file_list = [
            annotated_file_item["path"] for annotated_file_item in annotated_file_dict
        ]

        logger.info("File list created. Now reading in annotated shards.")

        # Import all annotated shards
        ht_annotated = hl.import_vcf(
            annotated_file_list,
            force_bgz=True,
            reference_genome=assembly,
        ).rows()
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
            v4_vrs_annotations(test=args.test).path,
            overwrite=args.overwrite,
        )
        logger.info(
            "Annotated Hail Table checkpointed to:"
            f" {v4_vrs_annotations(test=args.test).path}"
        )

    if args.annotate_original:
        output_path = v4_vrs_annotations(test=args.test, original_annotations=True)
        check_resource_existence(
            input_step_resources={
                "--run-vrs": [v4_vrs_annotations(test=args.test).path],
            },
            output_step_resources={
                "--annotate-original": [output_path],
            },
            overwrite=args.overwrite,
        )

        # Output final Hail Tables with VRS annotations
        ht_annotated = hl.read_table(v4_vrs_annotations(test=args.test).path)

        logger.info("Adding VRS IDs and GA4GH.VRS version to original Table")
        ht_final = ht_original.annotate(
            info=hl.struct(vrs=ht_annotated[ht_original.locus, ht_original.alleles].vrs)
        )

        ht_final = ht_final.annotate_globals(vrs_version=VRS_VERSION)

        logger.info(f"Outputting final table at: {output_path}")
        ht_final.write(output_path, overwrite=args.overwrite)

        logger.info(f"Done! Final table written to {output_path}.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--billing-project", help="Project to bill.", type=str)
    parser.add_argument(
        "--image",
        default="us-central1-docker.pkg.dev/broad-mpg-gnomad/images/vrs084",
        help="Image in a GCP Artifact Registry repository.",
        type=str,
    )
    parser.add_argument(
        "--working-bucket",
        help=(
            "Name of GCP Bucket to output intermediate files (sharded VCFs and"
            " checkpointed HTs) to."
        ),
        default="gnomad-tmp-4day",
        type=str,
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
        "--test",
        help="Fiter to only 200 partitions for testing purposes.",
        action="store_true",
    )
    parser.add_argument(
        "--header-path",
        help=(
            "Full path of txt file containing lines to append to VCF headers for fields"
            " that maybe be missing when exporting the Table to VCF."
        ),
        type=str,
        default=None,
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
