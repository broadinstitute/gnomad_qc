"""Script for running VQSR isolation forest on gnomAD v4 variant QC data."""
import argparse
import json
import logging
from typing import Any, Dict, List, Optional

import hail as hl
import hailtop.batch as hb

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("isolation_forest")
logger.setLevel(logging.INFO)

FEATURES = [
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_MQ",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_FS",
    "AS_SOR",
]


def extract_variant_annotations_job(
    b: hb.Batch,
    input_vcf: str,
    annotations: List[str],
    output_prefix: str,
    resource_args: str,
    extra_args: str,
    utils: Dict[str, Any],
    gcp_billing_project: str,
    gatk_jar: str = "/root/gatk.jar",
    out_bucket: Optional[str] = None,
):
    filename = f"{output_prefix}_extract_variant_annotations"
    outpath = f"{out_bucket}extract_variant_annotations/"

    j = b.new_job("GATK: ExtractVariantAnnotations")
    # couldn't find a public image with both gatk and bcftools installed
    j.image(utils["GATK_IMAGE"])
    cpu = utils.get("CPU", 1)
    mem_gb = utils.get("COMMAND_MEM_GB", 6)
    add_mem_gb = utils.get("ADDITIONAL_MEM_GB", 1)
    disk_size = utils.get("DISK_SIZE_GB", 100)

    j.cpu(cpu)
    j.memory(f"{int(mem_gb) + int(add_mem_gb)}G")
    j.storage(f"{disk_size}G")

    j.declare_resource_group(
        extract={
            "annotations_hdf5": "{root}.annot.hdf5",
            "extracted_vcf": "{root}.vcf.gz",
            "extracted_idx": "{root}.vcf.gz.tbi",
        }
    )

    cmd = f"""set -euo pipefail
        export GATK_LOCAL_JAR={gatk_jar}
        gatk --java-options "-Xmx{mem_gb}G" \\
            ExtractVariantAnnotations \\
                -V {input_vcf} \\
                -O {j.extract} \\
                -A {" -A ".join(annotations)} \\
                --gcs-project-for-requester-pays {gcp_billing_project} \\
                --use-allele-specific-annotations \\
                {resource_args} \\
                {extra_args}
    """

    j.command(cmd)

    if out_bucket:
        b.write_output(j.extract, f"{outpath}{filename}")

    return j


def train_variant_annotations_model_job(
    b: hb.Batch,
    annotations_hdf5,
    output_prefix: str,
    extra_args: str,
    utils: Dict[str, Any],
    gcp_billing_project: str,
    hyperparameters_json=None,
    gatk_jar: str = "/root/gatk.jar",
    out_bucket: Optional[str] = None,
):
    filename = f"{output_prefix}_train_variant_annotations_model"
    outpath = f"{out_bucket}train_variant_annotations_model/"

    j = b.new_job("GATK: TrainVariantAnnotationsModel")
    # couldn't find a public image with both gatk and bcftools installed
    j.image(utils["GATK_IMAGE"])
    cpu = utils.get("CPU", 1)
    mem_gb = utils.get("COMMAND_MEM_GB", 6)
    add_mem_gb = utils.get("ADDITIONAL_MEM_GB", 1)
    disk_size = utils.get("DISK_SIZE_GB", 100)

    j.cpu(cpu)
    j.memory(f"{int(mem_gb) + int(add_mem_gb)}G")
    j.storage(f"{disk_size}G")

    j.declare_resource_group(
        isolation_forest={
            "snp_scorer": "{root}.snp.scorer.pkl",
            "snp_train": "{root}.snp.trainingScores.hdf5",
            "snp_calibration": "{root}.snp.calibrationScores.hdf5",
            "indel_train": "{root}.indel.trainingScores.hdf5",
            "indel_calibration": "{root}.indel.calibrationScores.hdf5",
            "indel_scorer": "{root}.indel.scorer.pkl",
        }
    )
    hyperparameters_json_arg = ""

    if hyperparameters_json is not None:
        hyperparameters_json_arg = " --hyperparameters-json " + hyperparameters_json

    cmd = f"""set -euo pipefail
        export GATK_LOCAL_JAR={gatk_jar}
        gatk --java-options "-Xmx{mem_gb}G" \\
            TrainVariantAnnotationsModel \\
                --annotations-hdf5 {annotations_hdf5}{hyperparameters_json_arg} \\
                -O {j.isolation_forest} \\
                --model-backend PYTHON_IFOREST \\
                --gcs-project-for-requester-pays {gcp_billing_project} {extra_args}
    """

    j.command(cmd)

    if out_bucket:
        b.write_output(j.isolation_forest, f"{outpath}{filename}")

    return j


def split_intervals(
    b: hb.Batch,
    utils: Dict,
    gcp_billing_project: str,
):
    """
    Split genome into intervals to parallelize VQSR for large sample sizes
    :param b: Batch object to add jobs to
    :param utils: a dictionary containing resources (file paths and arguments) to be used to split genome
    :param gcp_billing_project: GCP billing project for requester-pays buckets
    :return: a Job object with a single output j.intervals of type ResourceGroup
    """
    j = b.new_job(f"""Make {utils['SCATTER_COUNT']} intervals""")
    j.image(utils["GATK_IMAGE"])
    java_mem = 3
    j.memory("standard")  # ~ 4G/core ~ 4G
    j.storage("16G")
    j.declare_resource_group(
        intervals={
            f"interval_{idx}": f"{{root}}/{str(idx).zfill(4)}-scattered.interval_list"
            for idx in range(utils["SCATTER_COUNT"])
        }
    )

    j.command(f"""set -e
    gatk --java-options "-Xms{java_mem}g" SplitIntervals \\
      -L {utils['CALLING_INTERVALS']} \\
      --interval-padding 150 \\
      -O {j.intervals} \\
      -scatter {utils['SCATTER_COUNT']} \\
      -R {utils['ref_fasta']} \\
      --gcs-project-for-requester-pays {gcp_billing_project} \\
      -mode INTERVAL_SUBDIVISION
      """)
    # Could save intervals to a bucket here to avoid rerunning the job
    return j


def score_variant_annotations_job(
    b: hb.Batch,
    input_vcf: str,
    annotations: List[str],
    extracted_vcf: str,
    model_prefix: str,
    output_prefix: str,
    extra_args: str,
    resource_args: str,
    utils: Dict[str, Any],
    gcp_billing_project: str,
    gatk_jar: str = "/root/gatk.jar",
    interval: Optional[str] = None,
    scatter_idx: Optional[int] = None,
    out_bucket: Optional[str] = None,
):
    filename = f"{output_prefix}_score_variant_annotations"
    outpath = f"{out_bucket}score_variant_annotations/"

    j = b.new_job("GATK: ScoreVariantAnnotations")
    # couldn't find a public image with both gatk and bcftools installed
    j.image(utils["GATK_IMAGE"])
    cpu = utils.get("CPU", 1)
    mem_gb = utils.get("COMMAND_MEM_GB", 6)
    add_mem_gb = utils.get("ADDITIONAL_MEM_GB", 1)
    disk_size = utils.get("DISK_SIZE_GB", 100)

    j.cpu(cpu)
    j.memory(f"{int(mem_gb) + int(add_mem_gb)}G")
    j.storage(f"{disk_size}G")

    j.declare_resource_group(
        output_score={
            "scored_vcf": "{root}.vcf.gz",
            "scored_vcf_idx": "{root}.vcf.gz.tbi",
            "annotations_hdf5": "{root}.annot.hdf5",
            "scores_hdf5": "{root}.scores.hdf5",
        }
    )

    cmd = f"""set -euo pipefail
            export GATK_LOCAL_JAR={gatk_jar}
            gatk --java-options "-Xmx{mem_gb}G" \\
                ScoreVariantAnnotations \\
                    -V {input_vcf} \\
                    -O {j.output_score} \\
                    -A {" -A ".join(annotations)} \\
                    --resource:extracted,extracted=true {extracted_vcf} \\
                    --model-prefix {model_prefix} \\
                    {f'-L {interval} ' if interval else ''} \\
                    --use-allele-specific-annotations \\
                    --gcs-project-for-requester-pays {gcp_billing_project} \\
                    --model-backend PYTHON_IFOREST {extra_args} \\
                    {resource_args}
        """

    j.command(cmd)

    if out_bucket:
        if scatter_idx is not None:
            b.write_output(j.output_score, f"{outpath}{filename}{scatter_idx}")
        else:
            b.write_output(j.output_score, f"{outpath}{filename}")

    return j


# VQSR workflow including scatter step
def isolation_forest_workflow(
    sites_only_vcf: str,
    output_prefix: str,
    annotations: List[str],
    utils: Dict[str, Any],
    resource_args: str,
    out_bucket: str,
    hyperparameters_json: str,
    extract_extra_args: str,
    train_extra_args: str,
    score_extra_args: str,
    batch_billing_project: str,
    gcp_billing_project: str,
    batch_suffix: str,
    gatk_jar: str = "/root/gatk.jar",
):
    """
    Wraps all the functions into a workflow

    :param sites_only_vcf: A concatenated, sites-only version of the sharded input
        VCFs; used for extracting training and calibration sets.
    :param output_prefix: Base prefix for output files. Sharded output VCFs will be
        named following the pattern:
        "{output_prefix}.{zero_based_shard_index}.score.vcf.gz".
    :param annotations: Annotations to be used for extraction, training, and scoring.
    :param resource_args: Resource arguments to be used for extraction and scoring.
        For example, "--resource:training_and_calibration_set,training=true,calibration=true gs://path-to-training-and-calibration-set ...".
        See GATK documentation for the ExtractVariantAnnotations and
        ScoreVariantAnnotations tools..
    :param hyperparameters_json: Optional JSON file specifying model hyperparameters to
        be used by TrainVariantAnnotationsModel. See GATK documentation for this tool.
    :param extract_extra_args: Optional catch-all string to provide additional
        arguments for ExtractVariantAnnotations. This can include intervals (as string
        arguments or non-localized files), variant-type modes, arguments for enabling
        positive-unlabeled learning, etc. The \"do-not-gzip-vcf-output\" argument is
        not supported by this workflow. See GATK documentation for this tool.
    :param train_extra_args: Optional catch-all string to provide additional arguments
        for TrainVariantAnnotationsModel. This can include variant-type modes, arguments
        for enabling positive-unlabeled learning, etc. See GATK documentation for this
        tool.
    :param score_extra_args: Optional catch-all string to provide additional arguments
        for ScoreVariantAnnotations. This can include intervals (as string arguments or
        non-localized files), variant-type modes, arguments for enabling
        positive-unlabeled learning and hard filtering, etc. The
        "do-not-gzip-vcf-output" argument is not supported by this workflow. See GATK
        documentation for this tool.
    :param out_bucket: path to write, plots, evaluation results, and recalibrated VCF to.
    :param batch_billing_project: Batch billing project to be used for the workflow.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :return:
    """

    tmp_vqsr_bucket = f"{out_bucket}/"

    backend = hb.ServiceBackend(
        billing_project=batch_billing_project,
        remote_tmpdir="gs://gnomad-tmp-4day/",
    )

    logger.info(
        f"Starting hail Batch with the project {batch_billing_project}, "
        f"bucket {tmp_vqsr_bucket}"
    )

    b = hb.Batch(
        f"isolation forest pipeline{batch_suffix}",
        backend=backend,
    )

    extract_variant_annotations = extract_variant_annotations_job(
        b=b,
        input_vcf=sites_only_vcf,
        annotations=annotations,
        output_prefix=output_prefix,
        resource_args=resource_args,
        extra_args=extract_extra_args,
        utils=utils,
        gcp_billing_project=gcp_billing_project,
        gatk_jar=gatk_jar,
        out_bucket=out_bucket,
    )

    train_variant_annotations_model = train_variant_annotations_model_job(
        b,
        annotations_hdf5=extract_variant_annotations.extract["annotations_hdf5"],
        output_prefix=output_prefix,
        extra_args=train_extra_args,
        utils=utils,
        gcp_billing_project=gcp_billing_project,
        hyperparameters_json=hyperparameters_json,
        gatk_jar=gatk_jar,
        out_bucket=out_bucket,
    )

    intervals_j = split_intervals(
        b=b, utils=utils, gcp_billing_project=gcp_billing_project
    )

    score_variant_annotations_jobs = [
        score_variant_annotations_job(
            b,
            sites_only_vcf,
            annotations,
            extracted_vcf=extract_variant_annotations.extract["extracted_vcf"],
            model_prefix=train_variant_annotations_model.isolation_forest,
            output_prefix=output_prefix,
            extra_args=score_extra_args,
            resource_args=resource_args,
            utils=utils,
            gcp_billing_project=gcp_billing_project,
            gatk_jar=gatk_jar,
            out_bucket=out_bucket,
            interval=intervals_j.intervals[f"interval_{idx}"],
            scatter_idx=idx,
        )
        for idx in range(utils["SCATTER_COUNT"])
    ]

    b.run()


def main(args):
    hl.init(
        backend="batch",
        tmp_dir="gs://gnomad-tmp-4day/",
        gcs_requester_pays_configuration=args.gcp_billing_project,
        regions=["us-central1"],
    )

    with open(args.resources, "r") as f:
        utils = json.load(f)

    print("billing project as: ", args.batch_billing_project)

    resource_args = f"""-resource:hapmap,training=true,calibration=true {utils['hapmap_resource_vcf']} \\
          -resource:omni,training=true,calibration=true {utils['omni_resource_vcf']} \\
          -resource:1000G,training=true,calibration=true {utils['one_thousand_genomes_resource_vcf']} \\
          -resource:dbsnp,training=true,calibration=true {utils['dbsnp_resource_vcf']} \\
          -resource:singletons,training=true,calibration=true {args.transmitted_singletons}
    """
    isolation_forest_workflow(
        sites_only_vcf=args.input_vcf,
        output_prefix=args.out_vcf_name,
        annotations=args.annotations,
        utils=utils,
        resource_args=resource_args,
        out_bucket=args.out_bucket,
        hyperparameters_json=args.hyperparameters_json,
        extract_extra_args="",
        train_extra_args="",
        score_extra_args="",
        batch_billing_project=args.batch_billing_project,
        gcp_billing_project=args.gcp_billing_project,
        batch_suffix=args.batch_suffix,
        gatk_jar="/root/gatk.jar",
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            "If the dataset should be filtered to chr22 for testing (also filtered to "
            "evaluation interval specified by --test-intervals)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--input-vcf",
        type=str,
        required=True,
        help="Input sites VCF containing AS annotations.",
    )
    parser.add_argument(
        "--out-bucket", type=str, required=True, help="Bucket to store VQSR outputs in."
    )
    parser.add_argument(
        "--out-vcf-name",
        type=str,
        required=True,
        help="Required prefix for VQSR outputs.",
    )
    parser.add_argument(
        "--resources",
        type=str,
        required=True,
        help=(
            "Path to .json file containing paths and information for the VQSR pipeline."
        ),
    )
    parser.add_argument(
        "--batch-billing-project",
        type=str,
        required=True,
        help="Hail Batch billing project.",
    )
    parser.add_argument(
        "--gcp-billing-project",
        type=str,
        required=True,
        help="Google Cloud billing project for reading requester pays buckets.",
    )
    parser.add_argument(
        "--batch-suffix",
        type=str,
        default="",
        help="String to add to end of batch name.",
    )
    parser.add_argument(
        "--annotations",
        help="Features to use in the random forests model.",
        default=FEATURES,
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--transmitted-singletons",
        type=str,
        required=False,
        help="Path to transmitted singletons or first singleton truth set VCF.",
    )
    parser.add_argument(
        "--hyperparameters-json",
        type=str,
        required=False,
        help="",
    )

    args = parser.parse_args()
    main(args)
