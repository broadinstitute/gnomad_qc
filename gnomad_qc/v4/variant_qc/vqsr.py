"""Script to run VQSR on an AS-annotated Sites VCF."""
import argparse
import json
import logging
from typing import Dict, List, Optional

import hail as hl
import hailtop.batch as hb
from hailtop.batch.job import Job

from gnomad_qc.v4.resources.annotations import get_true_positive_vcf_path, info_vcf_path
from gnomad_qc.v4.resources.basics import calling_intervals
from gnomad_qc.v4.resources.sample_qc import interval_qc_pass
from gnomad_qc.v4.resources.variant_qc import VQSR_FEATURES, get_variant_qc_result
from gnomad_qc.v4.variant_qc.import_variant_qc_vcf import (
    import_variant_qc_vcf as import_vqsr,
)

logger = logging.getLogger(__file__)
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger.setLevel(logging.INFO)


SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]


def split_intervals(
    b: hb.Batch,
    utils: Dict,
    calling_interval_path: str,
    gcp_billing_project: str,
    gatk_image: str,
    scatter_count: int = 100,
    interval_padding: int = 150,
) -> Job:
    """
    Split genome into intervals to parallelize VQSR for large sample sizes.

    :param b: Batch object to add jobs to.
    :param utils: Dictionary containing resources (file paths and arguments) to be
        used to split genome.
    :param calling_interval_path: Path to the calling interval list with no padding.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param scatter_count: Number of intervals to split the dataset into for scattering.
    :param interval_padding: Number of bases to pad each interval by.
    :return: Job object with a single output j.intervals of type ResourceGroup.
    """
    j = b.new_job(f"""Make {scatter_count} intervals""")
    j.image(gatk_image)
    java_mem = 3
    j.memory("standard")  # ~ 4G/core ~ 4G
    j.storage("16G")
    j.declare_resource_group(
        intervals={
            f"interval_{idx}": f"{{root}}/{str(idx).zfill(4)}-scattered.interval_list"
            for idx in range(scatter_count)
        }
    )

    # Modes other than INTERVAL_SUBDIVISION will produce an unpredicted number
    # of intervals. But we have to expect exactly the `scatter_count` number of
    # output files because our workflow is not dynamic.
    j.command(f"""set -e
    gatk --java-options "-Xms{java_mem}g" SplitIntervals \\
      -L {calling_interval_path} \\
      --interval-padding {interval_padding} \\
      -O {j.intervals} \\
      -scatter {scatter_count} \\
      -R {utils['ref_fasta']} \\
      --gcs-project-for-requester-pays {gcp_billing_project} \\
      -mode INTERVAL_SUBDIVISION
      """)

    return j


# SNPs
def snps_variant_recalibrator_create_model(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    use_as_annotations: bool,
    gcp_billing_project: str,
    gatk_image: str,
    features: List[str],
    singletons_resource_vcf: Optional[str] = None,
    evaluation_interval_path: Optional[str] = None,
    out_bucket: str = None,
    is_small_callset: bool = False,
    is_large_callset: bool = False,
    max_gaussians: int = 6,
) -> Job:
    """
    First step of VQSR for SNPs: run VariantRecalibrator to subsample variants and produce a file of the VQSR model.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:

        - Run the recalibrator with the following additional arguments:
          --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
        - Apply the resulting model to each genomic interval with, running the
          recalibrator with the same base parameters, plus:
          --input-model <model-file> --output-tranches-for-scatter
        - Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    :param b: Batch object to add jobs to.
    :param sites_only_vcf: Sites only VCF file to be used to build the model.
    :param utils: a dictionary containing resources (file paths) to be used to create
        the model.
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param features: List of features to be used in the SNP VQSR model.
    :param singletons_resource_vcf: Optional singletons VCF to be used in
        building the model.
    :param evaluation_interval_path: Optional path to the evaluation interval list.
    :param out_bucket: Full path to output bucket to write model and plots to.
    :param is_small_callset: Whether the dataset is small. Used to set number of CPUs
        for the job.
    :param is_large_callset: Whether the dataset is huge. Used to set number of CPUs
        for the job.
    :param max_gaussians: Maximum number of Gaussians for the positive model.
    :return: Job object with 2 outputs: j.model_file and j.snp_rscript_file.
    """
    j = b.new_job("VQSR: SNPsVariantRecalibratorCreateModel")
    j.image(gatk_image)
    j.memory("highmem")
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G.
    else:
        ncpu = 16  # ~ 8G/core ~ 128G.
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage("50G")

    downsample_factor = 75 if is_large_callset else 10

    tranche_cmdl = " ".join([f"-tranche {v}" for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = " ".join([f"-an {v}" for v in features])

    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms{java_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads={ncpu-2}" \\
          VariantRecalibrator \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          -L {evaluation_interval_path} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          {"--use-allele-specific-annotations" if use_as_annotations else ""} \\
          --sample-every-Nth-variant {downsample_factor} \\
          --output-model {j.model_file} \\
          --max-gaussians {max_gaussians} \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {utils['hapmap_resource_vcf']} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {utils['omni_resource_vcf']} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {utils['one_thousand_genomes_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {utils['dbsnp_resource_vcf']} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {singletons_resource_vcf}' if singletons_resource_vcf else ''} \\
          --rscript-file {j.snp_rscript}
          ls $(dirname {j.snp_rscript})
          ln {j.snp_rscript}.pdf {j.snp_rscript_pdf}
          ln {j.tranches}.pdf {j.tranches_pdf}
          """)

    if out_bucket:
        b.write_output(
            j.snp_rscript, f"{out_bucket}model/SNPS/snps.features.build.RScript"
        )
        b.write_output(
            j.snp_rscript_pdf, f"{out_bucket}model/SNPS/snps.features.build.pdf"
        )
        b.write_output(
            j.tranches_pdf, f"{out_bucket}model/SNPS/snps.tranches.build.pdf"
        )
        b.write_output(j.model_file, f"{out_bucket}model/SNPS/snps.model.report")

    return j


def snps_variant_recalibrator(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    out_bucket: str,
    use_as_annotations: bool,
    gcp_billing_project,
    gatk_image: str,
    features: List[str],
    interval: Optional[hb.ResourceGroup] = None,
    tranche_idx: Optional[int] = None,
    model_file: Optional[hb.ResourceFile] = None,
    singletons_resource_vcf: Optional[str] = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for SNPs: run VariantRecalibrator scattered to apply the VQSR model file to each genomic interval.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:

        - Run the recalibrator with the following additional arguments:
          --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
        - Apply the resulting model to each genomic interval with, running the
          recalibrator with the same base parameters, plus:
          --input-model <model-file> --output-tranches-for-scatter
        - Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    :param b: Batch object to add jobs to.
    :param sites_only_vcf: Sites only VCF file to be used to build the model.
    :param model_file: Model file to be applied.
    :param utils: Dictionary containing resources (file paths).
    :param out_bucket: Full path to output bucket to write model and plots to.
    :param tranche_idx: Index for the tranches file.
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param features: List of features to be used in the SNP VQSR model.
    :param singletons_resource_vcf: Optional singletons VCF to include in
        VariantRecalibrator.
    :param interval: Genomic interval to apply the model to.
    :param max_gaussians: Maximum number of Gaussians for the positive model.
    :return: Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches.
    """
    j = b.new_job("VQSR: SNPsVariantRecalibratorScattered")

    j.image(gatk_image)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes.
    j.memory(f"{mem_gb}G")
    j.cpu(4)
    j.storage("20G")

    j.declare_resource_group(recalibration={"index": "{root}.idx", "base": "{root}"})

    tranche_cmdl = " ".join([f"-tranche {v}" for v in SNP_RECALIBRATION_TRANCHE_VALUES])
    an_cmdl = " ".join([f"-an {v}" for v in features])

    cmd = f"""set -euo pipefail
        gatk --java-options "-Xms{mem_gb-1}g -XX:+UseParallelGC -XX:ParallelGCThreads=3" \\
          VariantRecalibrator \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode SNP \\
          --max-gaussians {max_gaussians} \\
          -resource:hapmap,known=false,training=true,truth=true,prior=15 {utils['hapmap_resource_vcf']} \\
          -resource:omni,known=false,training=true,truth=true,prior=12 {utils['omni_resource_vcf']} \\
          -resource:1000G,known=false,training=true,truth=false,prior=10 {utils['one_thousand_genomes_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=7 {utils['dbsnp_resource_vcf']} \\
        """

    if interval:
        cmd += f" -L {interval}"
    if use_as_annotations:
        cmd += " --use-allele-specific-annotations"
    if model_file:
        cmd += f" --input-model {model_file} --output-tranches-for-scatter"
    if singletons_resource_vcf:
        cmd += (
            " -resource:singletons,known=true,training=true,truth=true,prior=10"
            f" {singletons_resource_vcf}"
        )

    j.command(cmd)

    if out_bucket:
        if tranche_idx is not None:
            b.write_output(
                j.tranches,
                f"{out_bucket}recalibration/SNPS/snps.{tranche_idx}.tranches",
            )
        else:
            b.write_output(j.tranches, f"{out_bucket}recalibration/SNPS/snps.tranches")

    return j


# INDELs
def indels_variant_recalibrator_create_model(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    use_as_annotations: bool,
    gcp_billing_project: str,
    gatk_image: str,
    features: List[str],
    singletons_resource_vcf: str = None,
    evaluation_interval_path: Optional[str] = None,
    out_bucket: str = None,
    is_small_callset: bool = False,
    max_gaussians: int = 4,
) -> Job:
    """
    First step of VQSR for INDELs: run VariantRecalibrator to subsample variants and produce a file of the VQSR model.

    To support cohorts with more than 10,000 WGS samples, the INDEL recalibration
    process is broken down across genomic regions for parallel processing, and done in
    3 steps:

        - Run the recalibrator with the following additional arguments:
          --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
        - Apply the resulting model to each genomic interval with, running the
          recalibrator with the same base parameters, plus:
          --input-model <model-file> --output-tranches-for-scatter
        - Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value

    :param b: Batch object to add jobs to.
    :param sites_only_vcf: Sites only VCF file to be used to build the model.
    :param utils: Dictionary containing resources (file paths).
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param features: List of features to be used in the INDEL VQSR model.
    :param singletons_resource_vcf: Optional singletons VCF to be used in building the
        model.
    :param evaluation_interval_path: Optional path to the evaluation interval list.
    :param out_bucket: Full path to output bucket to write model and plots to.
    :param is_small_callset: Whether the dataset is small. Used to set number of CPUs
        for the job.
    :param max_gaussians: Maximum number of Gaussians for the positive model.
    :return: Job object with 2 outputs: j.model_file and j.indel_rscript_file. The
        latter is useful to produce the optional tranche plot.
    """
    j = b.new_job("VQSR: INDELsVariantRecalibratorCreateModel")
    j.image(gatk_image)
    j.memory("highmem")
    if is_small_callset:
        ncpu = 8  # ~ 8G/core ~ 64G
    else:
        ncpu = 16  # ~ 8G/core ~ 128G
    j.cpu(ncpu)
    java_mem = ncpu * 8 - 10
    j.storage("50G")

    downsample_factor = 75 if not is_small_callset else 10

    tranche_cmdl = " ".join(
        [f"-tranche {v}" for v in INDEL_RECALIBRATION_TRANCHE_VALUES]
    )
    an_cmdl = " ".join([f"-an {v}" for v in features])

    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms{java_mem}g -XX:+UseParallelGC -XX:ParallelGCThreads={ncpu-2}" \\
          VariantRecalibrator \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          -L {evaluation_interval_path} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          {"--use-allele-specific-annotations" if use_as_annotations else ""} \\
          --sample-every-Nth-variant {downsample_factor} \\
          --output-model {j.model_file} \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {utils['mills_resource_vcf']} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {utils['axiom_poly_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {utils['dbsnp_resource_vcf']} \\
          {f'-resource:singletons,known=true,training=true,truth=true,prior=10 {singletons_resource_vcf}' if singletons_resource_vcf else ''} \\
          --rscript-file {j.indel_rscript}
          ls $(dirname {j.indel_rscript})
          ln {j.indel_rscript}.pdf {j.indel_rscript_pdf}
        """)

    if out_bucket:
        b.write_output(
            j.indel_rscript, f"{out_bucket}model/INDELS/indels.features.build.RScript"
        )
        b.write_output(
            j.indel_rscript_pdf, f"{out_bucket}model/INDELS/indels.features.build.pdf"
        )
        b.write_output(j.model_file, f"{out_bucket}model/INDELS/indels.model.report")

    return j


def indels_variant_recalibrator(
    b: hb.Batch,
    sites_only_vcf: str,
    utils: Dict,
    out_bucket: str,
    use_as_annotations: bool,
    gcp_billing_project: str,
    gatk_image: str,
    features: List[str],
    interval: Optional[hb.ResourceGroup] = None,
    tranche_idx: Optional[int] = None,
    model_file: Optional[hb.ResourceFile] = None,
    singletons_resource_vcf: Optional[str] = None,
    max_gaussians: int = 4,
) -> Job:
    """
    Second step of VQSR for INDELs: run VariantRecalibrator scattered to apply the VQSR model file to each genomic interval.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:

        - Run the recalibrator with the following additional arguments:
          --sample-every-Nth-variant <downsample_factor> --output-model <model_file>
        - Apply the resulting model to each genomic interval with, running the
          recalibrator with the same base parameters, plus:
          --input-model <model-file> --output-tranches-for-scatter
        - Collate the resulting per-interval tranches with GatherTranches

    The --max-gaussians parameter sets the expected number of clusters in modeling.
    If a dataset gives fewer distinct clusters, e.g. as can happen for smaller data,
    then the tool will tell you there is insufficient data with a No data found error
    message. In this case, try decrementing the --max-gaussians value.

    :param b: Batch object to add jobs to.
    :param sites_only_vcf: Sites only VCF file to be used to build the model.
    :param model_file: Model file to be applied.
    :param utils: Dictionary containing resources (file paths).
    :param out_bucket: Full path to output bucket to write model and plots to.
    :param tranche_idx: Index for the tranches file.
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param features: List of features to be used in the INDEL VQSR model.
    :param singletons_resource_vcf: Optional singletons VCF to include in
        VariantRecalibrator.
    :param interval: Genomic interval to apply the model to.
    :param max_gaussians: Maximum number of Gaussians for the positive model.
    :return: Job object with 2 outputs: j.recalibration (ResourceGroup) and j.tranches.
    """
    j = b.new_job("VQSR: INDELsVariantRecalibratorScattered")

    j.image(gatk_image)
    mem_gb = 64  # ~ twice the sum of all input resources and input VCF sizes
    j.memory(f"{mem_gb}G")
    j.cpu(4)
    j.storage("20G")

    j.declare_resource_group(recalibration={"index": "{root}.idx", "base": "{root}"})

    tranche_cmdl = " ".join(
        [f"-tranche {v}" for v in INDEL_RECALIBRATION_TRANCHE_VALUES]
    )
    an_cmdl = " ".join([f"-an {v}" for v in features])

    cmd = f"""set -euo pipefail
        gatk --java-options "-Xms{mem_gb-1}g -XX:+UseParallelGC -XX:ParallelGCThreads=3" \\
          VariantRecalibrator \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          -V {sites_only_vcf} \\
          -O {j.recalibration} \\
          --tranches-file {j.tranches} \\
          --trust-all-polymorphic \\
          {tranche_cmdl} \\
          {an_cmdl} \\
          -mode INDEL \\
          --max-gaussians {max_gaussians} \\
          -resource:mills,known=false,training=true,truth=true,prior=12 {utils['mills_resource_vcf']} \\
          -resource:axiomPoly,known=false,training=true,truth=false,prior=10 {utils['axiom_poly_resource_vcf']} \\
          -resource:dbsnp,known=true,training=false,truth=false,prior=2 {utils['dbsnp_resource_vcf']} \\
        """

    if interval:
        cmd += f" -L {interval}"
    if use_as_annotations:
        cmd += " --use-allele-specific-annotations"
    if model_file:
        cmd += f" --input-model {model_file} --output-tranches-for-scatter"
    if singletons_resource_vcf:
        cmd += (
            " -resource:singletons,known=true,training=true,truth=true,prior=10"
            f" {singletons_resource_vcf}"
        )

    j.command(cmd)

    if out_bucket:
        if tranche_idx is not None:
            b.write_output(
                j.tranches,
                f"{out_bucket}recalibration/INDELS/indels.{tranche_idx}.tranches",
            )
        else:
            b.write_output(
                j.tranches, f"{out_bucket}recalibration/INDELS/indels.tranches"
            )
    return j


def gather_tranches(
    b: hb.Batch,
    tranches: List[hb.ResourceFile],
    mode: str,
    disk_size: int,
    gcp_billing_project: str,
) -> Job:
    """
    Third step of VQSR for SNPs: run GatherTranches to gather scattered per-interval tranches outputs.

    To support cohorts with more than 10,000 WGS samples, the SNP recalibration process
    is broken down across genomic regions for parallel processing, and done in 3 steps:

        - Run the recalibrator with the following additional arguments:
          ``--sample-every-Nth-variant <downsample_factor> --output-model <model_file>``
        - Apply the resulting model to each genomic interval with, running the
          recalibrator with the same base parameters, plus:
          ``--input-model <model-file> --output-tranches-for-scatter``
        - Collate the resulting per-interval tranches with GatherTranches

    :param b: Batch object to add jobs to.
    :param tranches: Index for the tranches file.
    :param mode: Recalibration mode to employ, either SNP or INDEL.
    :param disk_size: Disk size to be used for the job.
    :return: Job object with one output j.out_tranches.
    """
    j = b.new_job(f"VQSR: {mode}GatherTranches")
    # TODO: Should this switch to a GATK image?
    j.image(
        "us.gcr.io/broad-dsde-methods/gatk-for-ccdg@sha256:9e9f105ecf3534fbda91a4f2c2816ec3edf775882917813337a8d6e18092c959"
    )
    j.memory("8G")
    j.cpu(2)
    j.storage(f"{disk_size}G")

    inputs_cmdl = " ".join([f"--input {t}" for t in tranches])
    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms6g" \\
          GatherTranches \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --mode {mode} \\
          {inputs_cmdl} \\
          --output {j.out_tranches}""")

    return j


def apply_recalibration(
    b: hb.Batch,
    input_vcf: str,
    out_vcf_name: str,
    indels_recalibration: hb.ResourceGroup,
    indels_tranches: hb.ResourceFile,
    snps_recalibration: hb.ResourceGroup,
    snps_tranches: hb.ResourceFile,
    snp_hard_filter: float,
    indel_hard_filter: float,
    disk_size: int,
    use_as_annotations: bool,
    gcp_billing_project: str,
    scatter: Optional[int] = None,
    interval: Optional[hb.ResourceGroup] = None,
    out_bucket: Optional[str] = None,
    overlap_skip: bool = False,
) -> Job:
    """
    Apply a score cutoff to filter variants based on a recalibration table.

    Runs ApplyVQSR twice to apply first indel, then SNP recalibrations.

    Targets indel_filter_level and snp_filter_level sensitivities. The tool matches
    them internally to a VQSLOD score cutoff based on the model's estimated sensitivity
    to a set of true variants.

    The filter determination is not just a pass/fail process. The tool evaluates for
    each variant which "tranche", or slice of the dataset, it falls into in terms of
    sensitivity to the truthset. Variants in tranches that fall below the specified
    truth sensitivity filter level have their FILTER field annotated with the
    corresponding tranche level. This results in a callset that is filtered to the
    desired level but retains the information necessary to increase sensitivity
    if needed.

    :param b: Batch object to add jobs to.
    :param input_vcf: Sites only VCF file to be used.
    :param out_vcf_name: Output vcf filename.
    :param indels_recalibration: Input recal file (ResourceGroup) for INDELs.
    :param indels_tranches: Input tranches file (ResourceFile) for INDELs.
    :param snps_recalibration: Input recal file (ResourceGroup) for SNPs.
    :param snps_tranches: Input tranches file (ResourceFile) for SNPs.
    :param snp_hard_filter: SNP hard filter level.
    :param indel_hard_filter: INDEL hard filter level.
    :param disk_size: Disk size to be used for the job.
    :param use_as_annotations: If set, Allele-Specific variant recalibrator will be used.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param scatter: Scatter index to be used in output VCF filename if running in
        scattered mode.
    :param interval: Genomic interval to apply the model to.
    :param out_bucket: Full path to output bucket to write output(s) to.
    :return: Job object with one ResourceGroup output j.output_vcf, corresponding
        to a VCF with tranche annotated in the FILTER field.
    """
    if scatter is not None:
        filename = f"{out_vcf_name}_vqsr_recalibrated_{scatter}"
        outpath = f"{out_bucket}apply_recalibration/scatter/"
    else:
        filename = f"{out_vcf_name}_vqsr_recalibrated"
        outpath = out_bucket

    j = b.new_job("VQSR: ApplyRecalibration")
    # Custom image from Lindo Nkambule with GATK and BCFtools installed
    # TODO: follow up on what version and how maintained
    j.image("docker.io/lindonkambule/vqsr_gatk_bcftools_img:latest")
    j.memory("8G")
    j.storage(f"{disk_size}G")
    j.declare_resource_group(
        output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"},
        intermediate_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"},
    )

    j.command(f"""set -euo pipefail
        gatk --java-options "-Xms5g" \\
          ApplyVQSR \\
          -O tmp.indel.recalibrated.vcf \\
          -V {input_vcf} \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --recal-file {indels_recalibration} \\
          --tranches-file {indels_tranches} \\
          --truth-sensitivity-filter-level {indel_hard_filter} \\
          --create-output-variant-index true \\
          {f'-L {interval} ' if interval else ''} \\
          {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
          -mode INDEL

        gatk --java-options "-Xms5g" \\
          ApplyVQSR \\
          -O {j.output_vcf['vcf.gz']} \\
          -V tmp.indel.recalibrated.vcf \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --recal-file {snps_recalibration} \\
          --tranches-file {snps_tranches} \\
          --truth-sensitivity-filter-level {snp_hard_filter} \\
          --create-output-variant-index true \\
          {f'-L {interval} ' if interval else ''} \\
          {'--use-allele-specific-annotations ' if use_as_annotations else ''} \\
          -mode SNP""")

    # NOTE: An INDEL at the beginning of a chunk will overlap with the previous chunk
    #  and will cause issues when trying to merge. This makes sure the INDEL is ONLY
    #  in ONE of two consecutive chunks (not both).
    if interval and not overlap_skip:
        # Overwrite VCF with overlap issue addressed.
        j.command(f"""
                cat {interval} | grep -Ev "^@" > tmp.interval
                interval=$(cat tmp.interval | awk 'NR == 1 {{c=$1; s=$2}} END {{ print c":"s"-"$1":"$3}}')
                bcftools view -t $interval {j.output_vcf['vcf.gz']} --output-file {j.output_vcf['vcf.gz']} --output-type z
                tabix -f {j.output_vcf['vcf.gz']}
                df -h; pwd; du -sh $(dirname {j.output_vcf['vcf.gz']})
            """)

    if out_bucket:
        b.write_output(j.output_vcf, f"{outpath}{filename}")

    return j


def gather_vcfs(
    b: hb.Batch,
    input_vcfs: List[hb.ResourceGroup],
    out_vcf_name: str,
    disk_size: int,
    gcp_billing_project: str,
    gatk_image: str,
    out_bucket: str = None,
) -> Job:
    """
    Combine recalibrated VCFs into a single VCF & saves the output VCF to a bucket `out_bucket`.

    :param b: Batch object to add jobs to.
    :param input_vcfs: List of VCFs to be gathered.
    :param out_vcf_name: Output vcf filename.
    :param disk_size: Disk size to be used for the job.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param out_bucket: Full path to output bucket to write the gathered VCF to.
    :return: Job object with one ResourceGroup output j.output_vcf.
    """
    filename = f"{out_vcf_name}_vqsr_recalibrated"
    j = b.new_job("VQSR: FinalGatherVcf")
    j.image(gatk_image)
    j.memory(f"16G")
    j.storage(f"{disk_size}G")

    j.declare_resource_group(
        output_vcf={"vcf.gz": "{root}.vcf.gz", "vcf.gz.tbi": "{root}.vcf.gz.tbi"}
    )

    input_cmdl = " ".join([f'--input {v["vcf.gz"]}' for v in input_vcfs])
    j.command(f"""set -euo pipefail
        # --ignore-safety-checks makes a big performance difference so we include it in our invocation.
        # This argument disables expensive checks that the file headers contain the same set of
        # genotyped samples and that files are in order by position of first record.
        cd /io
        mkdir tmp/
        gatk --java-options "-Xms6g -Djava.io.tmpdir=`pwd`/tmp" \\
          GatherVcfsCloud \\
          --gcs-project-for-requester-pays {gcp_billing_project} \\
          --ignore-safety-checks \\
          --gather-type BLOCK \\
          {input_cmdl} \\
          --output {j.output_vcf['vcf.gz']} \\
          --tmp-dir `pwd`/tmp

        tabix {j.output_vcf['vcf.gz']}""")
    if out_bucket:
        b.write_output(j.output_vcf, f"{out_bucket}{filename}")
    return j


def make_vqsr_jobs(
    b: hb.Batch,
    sites_only_vcf: str,
    is_small_callset: bool,
    is_large_callset: bool,
    output_vcf_name: str,
    utils: Dict,
    out_bucket: str,
    intervals: Dict,
    use_as_annotations: bool,
    gcp_billing_project: str,
    gatk_image: str,
    snp_features: List[str],
    indel_features: List[str],
    snp_hard_filter: float,
    indel_hard_filter: float,
    scatter_count: int = 100,
    singleton_vcf_path: Optional[str] = None,
    evaluation_interval_path: Optional[str] = None,
    overlap_skip: bool = False,
) -> None:
    """
    Add jobs to Batch that perform the allele-specific VQSR variant QC.

    :param b: Batch object to add jobs to.
    :param sites_only_vcf: Path to a sites only VCF created using gnomAD
        default_compute_info().
    :param is_small_callset: For small callsets, we reduce the resources per job from
        our default.
    :param is_large_callset: For large callsets, we increase resources per job from
        our default and shard our VCF to launch jobs in a scattered mode.
    :param output_vcf_name: Name, without extension, to use for the output VCF file(s).
    :param utils: Dictionary containing resource files and parameters to be used in VQSR.
    :param out_bucket: Path to write, plots, evaluation results, and recalibrated VCF to.
    :param intervals: ResourceGroup object with intervals to scatter.
    :param use_as_annotations: Whether to use allele-specific annotation for VQSR.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param gatk_image: GATK docker image.
    :param snp_features: List of features to be used in the SNP VQSR model.
    :param indel_features: List of features to be used in the INDEL VQSR model.
    :param snp_hard_filter: SNP hard filter level.
    :param indel_hard_filter: INDEL hard filter level.
    :param scatter_count: Number of intervals to scatter over.
    :param singleton_vcf_path: Full path to transmitted and/or sibling singletons VCF
        file and its index.
    :param evaluation_interval_path: Optional full path to evaluation intervals file.
    :return: None.
    """
    # Scale resources for jobs as appropriate.
    if is_small_callset:
        small_disk = 50
    elif not is_large_callset:
        small_disk = 100
    else:
        small_disk = 200

    if is_small_callset:
        huge_disk = 200
    elif not is_large_callset:
        huge_disk = 500
    else:
        huge_disk = 2000

    snp_max_gaussians = 6
    indel_max_gaussians = 4

    # If it is a large callset, run in scatter mode.
    if is_large_callset:
        # 1. Run SNP recalibrator in a scattered mode.
        # Check if model already exists:
        if hl.hadoop_exists(f"{out_bucket}model/SNPS/snps.model.report"):
            print(
                "Found existing model for SNPs:"
                f" {out_bucket}model/SNPS/snps.model.report"
            )
            snps_model_file = b.read_input(f"{out_bucket}model/SNPS/snps.model.report")
        else:
            snps_model_file = snps_variant_recalibrator_create_model(
                b=b,
                sites_only_vcf=sites_only_vcf,
                singletons_resource_vcf=singleton_vcf_path,
                evaluation_interval_path=evaluation_interval_path,
                utils=utils,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                gatk_image=gatk_image,
                features=snp_features,
                out_bucket=out_bucket,
                is_small_callset=is_small_callset,
                is_large_callset=is_large_callset,
                max_gaussians=snp_max_gaussians,
            ).model_file

        snps_recalibrator_jobs = [
            snps_variant_recalibrator(
                b=b,
                sites_only_vcf=sites_only_vcf,
                utils=utils,
                out_bucket=out_bucket,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                gatk_image=gatk_image,
                features=snp_features,
                interval=intervals[f"interval_{idx}"],
                tranche_idx=idx,
                model_file=snps_model_file,
                singletons_resource_vcf=singleton_vcf_path,
                max_gaussians=snp_max_gaussians,
            )
            for idx in range(scatter_count)
        ]

        snps_recalibrations = [j.recalibration for j in snps_recalibrator_jobs]
        snps_tranches = [j.tranches for j in snps_recalibrator_jobs]
        snps_gathered_tranches = gather_tranches(
            b=b,
            tranches=snps_tranches,
            mode="SNP",
            disk_size=small_disk,
            gcp_billing_project=gcp_billing_project,
        ).out_tranches

        # 2. Run INDEL recalibrator in a scattered mode.
        if hl.hadoop_exists(f"{out_bucket}model/INDELS/indels.model.report"):
            print(
                "Found existing model for INDELs:"
                f" {out_bucket}model/INDELS/indels.model.report"
            )
            indels_model_file = b.read_input(
                f"{out_bucket}model/INDELS/indels.model.report"
            )
        else:
            indels_model_file = indels_variant_recalibrator_create_model(
                b=b,
                sites_only_vcf=sites_only_vcf,
                singletons_resource_vcf=singleton_vcf_path,
                evaluation_interval_path=evaluation_interval_path,
                utils=utils,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                gatk_image=gatk_image,
                features=indel_features,
                out_bucket=out_bucket,
                is_small_callset=is_small_callset,
                max_gaussians=indel_max_gaussians,
            ).model_file

        indels_recalibrator_jobs = [
            indels_variant_recalibrator(
                b=b,
                sites_only_vcf=sites_only_vcf,
                utils=utils,
                out_bucket=out_bucket,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                gatk_image=gatk_image,
                features=indel_features,
                interval=intervals[f"interval_{idx}"],
                tranche_idx=idx,
                model_file=indels_model_file,
                singletons_resource_vcf=singleton_vcf_path,
                max_gaussians=indel_max_gaussians,
            )
            for idx in range(scatter_count)
        ]

        indels_recalibrations = [j.recalibration for j in indels_recalibrator_jobs]
        indels_tranches = [j.tranches for j in indels_recalibrator_jobs]
        indels_gathered_tranches = gather_tranches(
            b=b,
            tranches=indels_tranches,
            mode="INDEL",
            disk_size=small_disk,
            gcp_billing_project=gcp_billing_project,
        ).out_tranches

        # 3. Apply recalibration.
        # Doesn't require too much storage in scatter mode (<500MB on gnomad VCF
        # for each scatter), use small_disk.
        scattered_vcfs = [
            apply_recalibration(
                b=b,
                input_vcf=sites_only_vcf,
                out_vcf_name=output_vcf_name,
                indels_recalibration=indels_recalibrations[idx],
                indels_tranches=indels_gathered_tranches,
                snps_recalibration=snps_recalibrations[idx],
                snps_tranches=snps_gathered_tranches,
                snp_hard_filter=snp_hard_filter,
                indel_hard_filter=indel_hard_filter,
                disk_size=small_disk,
                use_as_annotations=use_as_annotations,
                gcp_billing_project=gcp_billing_project,
                scatter=idx,
                interval=intervals[f"interval_{idx}"],
                out_bucket=out_bucket,
                overlap_skip=overlap_skip,
            ).output_vcf
            for idx in range(scatter_count)
        ]

        # 4. Gather VCFs.
        if not overlap_skip:
            gather_vcfs(
                b=b,
                input_vcfs=scattered_vcfs,
                out_vcf_name=output_vcf_name,
                disk_size=huge_disk,
                out_bucket=out_bucket,
                gcp_billing_project=gcp_billing_project,
                gatk_image=gatk_image,
            )
        else:
            logger.info("Read VCF from shards. Remember to check for duplicates.")

    # If not is_large_callset, no need to run as scattered.
    else:
        snps_recalibrator_job = snps_variant_recalibrator(
            b=b,
            sites_only_vcf=sites_only_vcf,
            utils=utils,
            out_bucket=out_bucket,
            use_as_annotations=use_as_annotations,
            gcp_billing_project=gcp_billing_project,
            gatk_image=gatk_image,
            features=snp_features,
            singletons_resource_vcf=singleton_vcf_path,
            max_gaussians=snp_max_gaussians,
        )
        snps_recalibration = snps_recalibrator_job.recalibration
        snps_tranches = snps_recalibrator_job.tranches

        indels_variant_recalibrator_job = indels_variant_recalibrator(
            b=b,
            sites_only_vcf=sites_only_vcf,
            singletons_resource_vcf=singleton_vcf_path,
            utils=utils,
            out_bucket=out_bucket,
            use_as_annotations=use_as_annotations,
            gcp_billing_project=gcp_billing_project,
            gatk_image=gatk_image,
            features=indel_features,
            max_gaussians=indel_max_gaussians,
        )
        indels_recalibration = indels_variant_recalibrator_job.recalibration
        indels_tranches = indels_variant_recalibrator_job.tranches

        apply_recalibration(
            b=b,
            input_vcf=sites_only_vcf,
            out_vcf_name=output_vcf_name,
            indels_recalibration=indels_recalibration,
            indels_tranches=indels_tranches,
            snps_recalibration=snps_recalibration,
            snps_tranches=snps_tranches,
            snp_hard_filter=snp_hard_filter,
            indel_hard_filter=indel_hard_filter,
            disk_size=huge_disk,
            use_as_annotations=use_as_annotations,
            gcp_billing_project=gcp_billing_project,
            out_bucket=out_bucket,
        )

        return


def main(args):
    """Run VQSR variant qc workflow."""
    hl.init(
        backend="batch",
        tmp_dir="gs://gnomad-tmp-4day",
        gcs_requester_pays_configuration=args.gcp_billing_project,
        default_reference="GRCh38",
        regions=["us-central1"],
    )

    test = args.test
    scatter_count = args.scatter_count
    n_partitions = args.n_partitions
    transmitted_singletons = args.transmitted_singletons
    sibling_singletons = args.sibling_singletons
    snp_hard_filter = args.snp_hard_filter
    indel_hard_filter = args.indel_hard_filter
    true_positive_type = None
    overlap_skip = args.overlap_skip
    if transmitted_singletons and sibling_singletons:
        true_positive_type = "transmitted_singleton.sibling_singleton"
    elif transmitted_singletons:
        true_positive_type = "transmitted_singleton"
    elif sibling_singletons:
        true_positive_type = "sibling_singleton"

    logger.info(
        "Batch billing project as: %s, and GCP billing project as: %s}",
        args.batch_billing_project,
        args.gcp_billing_project,
    )

    # Run entire vqsr workflow, including generating and applying models.
    tmp_vqsr_bucket = f"{args.out_bucket}/"

    backend = hb.ServiceBackend(
        billing_project=args.batch_billing_project,
        remote_tmpdir="gs://gnomad-tmp-4day/",
    )

    with open(args.resources, "r") as f:
        utils = json.load(f)

    logger.info(
        f"Starting hail Batch with the project {args.batch_billing_project}, "
        f"bucket {tmp_vqsr_bucket}"
    )

    calling_interval_ht = calling_intervals(
        interval_name="union",
        calling_interval_padding=0,
    ).ht()

    if args.interval_qc_filter:
        evaluation_interval_ht = interval_qc_pass(all_platforms=True).ht()
        evaluation_interval_ht = evaluation_interval_ht.filter(
            evaluation_interval_ht.pass_interval_qc
        )
    elif args.calling_interval_filter:
        evaluation_interval_ht = calling_intervals(
            interval_name="intersection",
            calling_interval_padding=50,
        ).ht()
    else:
        evaluation_interval_ht = None

    if test:
        n_partitions = 10
        scatter_count = 10
        snp_hard_filter = 90.0
        indel_hard_filter = 90.0
        calling_interval_ht = calling_interval_ht.filter(
            calling_interval_ht.interval.start.contig == "chr22"
        )
        if evaluation_interval_ht is not None:
            evaluation_interval_ht = evaluation_interval_ht.filter(
                evaluation_interval_ht.interval.start.contig == "chr22"
            )

    # Write out intervals to a temp text file.
    calling_interval_path = hl.utils.new_temp_file("calling_intervals", "intervals")
    int_expr = calling_interval_ht.interval
    calling_interval_ht = calling_interval_ht.key_by(
        interval=int_expr.start.contig
        + ":"
        + hl.str(int_expr.start.position)
        + "-"
        + hl.str(int_expr.end.position)
    ).select()
    calling_interval_ht.export(calling_interval_path, header=False)

    if evaluation_interval_ht is not None:
        evaluation_interval_path = hl.utils.new_temp_file(
            "evaluation_intervals", "intervals"
        )
        int_expr = evaluation_interval_ht.interval
        evaluation_interval_ht = evaluation_interval_ht.key_by(
            interval=int_expr.start.contig
            + ":"
            + hl.str(int_expr.start.position)
            + "-"
            + hl.str(int_expr.end.position)
        ).select()
        evaluation_interval_ht.export(evaluation_interval_path, header=False)
    else:
        evaluation_interval_path = calling_interval_path

    singleton_vcf_path = get_true_positive_vcf_path(
        adj=args.adj,
        true_positive_type=true_positive_type,
    )

    b = hb.Batch(
        f"VQSR pipeline {args.batch_suffix}",
        backend=backend,
    )

    # Split passed intervals into value specified in resources json.
    # These are used for scattered jobs, if required.
    intervals_j = split_intervals(
        b=b,
        utils=utils,
        calling_interval_path=calling_interval_path,
        gcp_billing_project=args.gcp_billing_project,
        gatk_image=args.gatk_image,
        scatter_count=scatter_count,
        interval_padding=150,
    )

    # Converts run_mode into the two booleans used later in the code.
    is_small_callset = False
    is_large_callset = False
    if args.run_mode == "small":
        is_small_callset = True
    elif args.run_mode == "large":
        is_large_callset = True

    # Configure all VQSR jobs.
    make_vqsr_jobs(
        b=b,
        sites_only_vcf=info_vcf_path(info_method=args.compute_info_method),
        is_small_callset=is_small_callset,
        is_large_callset=is_large_callset,
        output_vcf_name=args.out_vcf_name,
        scatter_count=scatter_count,
        utils=utils,
        out_bucket=tmp_vqsr_bucket,
        intervals=intervals_j.intervals,
        use_as_annotations=True,
        gcp_billing_project=args.gcp_billing_project,
        gatk_image=args.gatk_image,
        snp_features=args.snp_features,
        indel_features=args.indel_features,
        snp_hard_filter=snp_hard_filter,
        indel_hard_filter=indel_hard_filter,
        singleton_vcf_path=singleton_vcf_path,
        evaluation_interval_path=evaluation_interval_path,
        overlap_skip=overlap_skip,
    )
    # Run all jobs, as loaded into the Batch b.
    b.run()

    logger.info("VQSR Batch jobs executed successfully")

    if args.load_vqsr:
        if is_large_callset and overlap_skip:
            outpath = f"{tmp_vqsr_bucket}apply_recalibration/scatter/{args.out_vcf_name}_vqsr_recalibrated_*.vcf.gz"
        else:
            outpath = f"{tmp_vqsr_bucket}{args.out_vcf_name}_vqsr_recalibrated.vcf.gz"

        hts = import_vqsr(
            outpath,
            args.model_id,
            n_partitions,
            args.header_path,
            args.array_elements_required,
            is_split=False,
            deduplicate_check=True,
        )

        for ht, split in zip(hts, [True, False]):
            ht = ht.annotate_globals(
                transmitted_singletons=args.transmitted_singletons,
                sibling_singletons=args.sibling_singletons,
                adj=args.adj,
                interval_qc_filter=args.interval_qc_filter,
                calling_interval_filter=args.calling_interval_filter,
                compute_info_method=args.compute_info_method,
                indel_features=args.indel_features,
                snp_features=args.snp_features,
            )
            ht.checkpoint(
                get_variant_qc_result(args.model_id, split=split).path,
                overwrite=args.overwrite,
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Whether to overwrite data already present in the output Table.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="If the dataset should be filtered to chr22 for testing.",
        action="store_true",
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
        "--model-id",
        help="Model ID for the variant QC result HT.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--header-path",
        help=(
            "Optional path to a header file to use for importing the variant QC result"
            " VCF."
        ),
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of desired partitions for output Table.",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--array-elements-required",
        action="store_true",
        help="Pass if you would like array elements required in import_vcf to be true.",
    )
    parser.add_argument(
        "--gatk-image",
        default="us.gcr.io/broad-gatk/gatk:4.2.6.1",
        help="GATK Image.",
        type=str,
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
        "--scatter-count",
        type=int,
        default=100,
        help="Number of intervals to scatter across.",
    )
    parser.add_argument(
        "--snp-hard-filter",
        type=float,
        default=99.7,
        help="Hard filter cutoff for SNPs.",
    )
    parser.add_argument(
        "--indel-hard-filter",
        type=float,
        default=99.0,
        help="Hard filter cutoff for INDELs.",
    )
    parser.add_argument(
        "--run-mode",
        type=str,
        default="standard",
        choices=["small", "standard", "large"],
        help=(
            "Option to pass so that the mode/resources fit the size of the database."
            " This affects the size of the clusters and, if --large is set, will run in"
            " a scattered mode (one job for each partition)."
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
        "--compute-info-method",
        help=(
            "Method of computing the INFO score to use for the variant QC features. "
            "Default is 'AS'."
        ),
        default="quasi",
        type=str,
        choices=["AS", "quasi", "set_long_AS_missing"],
    )
    parser.add_argument(
        "--adj",
        help="Use adj genotypes for transmitted/sibling singletons.",
        action="store_true",
    )
    parser.add_argument(
        "--transmitted-singletons",
        help="Include transmitted singletons in training.",
        action="store_true",
    )
    parser.add_argument(
        "--sibling-singletons",
        help="Include sibling singletons in training.",
        action="store_true",
    )
    parser.add_argument(
        "--interval-qc-filter",
        help="Whether interval QC should be applied for training.",
        action="store_true",
    )
    parser.add_argument(
        "--calling-interval-filter",
        help=(
            "Whether to filter to the intersection of Broad/DSP calling intervals with "
            "50 bp of padding for training."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--snp-features",
        help="Features to use in the SNP VQSR model.",
        default=VQSR_FEATURES["exomes"]["snv"],
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--indel-features",
        help="Features to use in the indel VQSR model.",
        default=VQSR_FEATURES["exomes"]["indel"],
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--overlap-skip",
        help="Skip code to address overlaps, output sharded VCF.",
        action="store_true",
    )
    parser.add_argument(
        "--batch-suffix",
        type=str,
        default="",
        help="String to add to end of batch name.",
    )
    parser.add_argument(
        "--load-vqsr",
        help="Run import_variant_qc_vcf() to load VQSR VCFs as HT",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
