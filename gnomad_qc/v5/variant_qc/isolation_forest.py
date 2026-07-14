"""Script to run the GATK isolation forest variant QC model on an AS-annotated sites VCF.

SNPs and INDELs are modeled separately (GATK ``--mode``): for each variant type the
workflow runs ExtractVariantAnnotations -> TrainVariantAnnotationsModel (PYTHON_IFOREST)
-> ScoreVariantAnnotations (scattered over intervals). The two per-mode scored VCFs are
merged into a single result HT in the ``--load-iforest`` step.
"""

import argparse
import logging
from typing import List, Optional

import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.utils.file_utils import file_exists
from hailtop.batch.job import Job

from gnomad_qc.v5.resources.annotations import get_true_positive_vcf_path, info_vcf_path
from gnomad_qc.v5.resources.basics import _init_hail
from gnomad_qc.v5.resources.constants import BATCH_TMP_BUCKET
from gnomad_qc.v5.resources.variant_qc import (
    IF_FEATURES,
    get_iforest_run_prefix,
    get_variant_qc_result,
)
from gnomad_qc.v5.variant_qc.import_variant_qc_vcf import import_variant_qc_vcf

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("isolation_forest")
logger.setLevel(logging.INFO)

DEFAULT_GATK_IMAGE = "us.gcr.io/broad-gatk/gatk:4.6.1.0"
"""GATK image. 4.6.1.0 ships the conda env (scikit-learn, dill) needed for PYTHON_IFOREST."""

CALLING_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
"""Contigs to score (autosomes and sex chromosomes)."""

# Public GATK resource VCFs (broad-references); used as labeled training/calibration
# sets. SNP and INDEL modes use different truth resources.
REFERENCE_RESOURCES = {
    "ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
    "dbsnp": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz",
    "hapmap": "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz",
    "omni": "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
    "one_thousand_genomes": "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    "mills": "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "axiom_poly": "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
}


def _resource_args(mode: str, singletons_vcf: Optional[str]) -> str:
    """
    Build the GATK ``--resource`` args (labeled training/calibration sets) for a mode.

    :param mode: Variant type mode, 'SNP' or 'INDEL'.
    :param singletons_vcf: Optional true-positive (transmitted/sibling singletons) VCF.
    :return: Space-joined GATK resource argument string.
    """
    if mode == "SNP":
        resources = ["hapmap", "omni", "one_thousand_genomes", "dbsnp"]
    else:
        resources = ["mills", "axiom_poly", "dbsnp"]

    args = [
        f"--resource:{r},training=true,calibration=true {REFERENCE_RESOURCES[r]}"
        for r in resources
    ]
    if singletons_vcf:
        args.append(
            f"--resource:singletons,training=true,calibration=true {singletons_vcf}"
        )
    return " ".join(args)


def _annotation_args(features: List[str]) -> str:
    """Build the GATK ``-A`` annotation args from a feature list."""
    return " ".join(f"-A {f}" for f in features)


def extract_variant_annotations_job(
    b: hb.Batch,
    mode: str,
    sites_only_vcf: str,
    features: List[str],
    resource_args: str,
    out_root: str,
    gatk_image: str,
    gcp_billing_project: str,
) -> Job:
    """
    Run GATK ExtractVariantAnnotations for one variant type mode.

    :param b: Batch to add the job to.
    :param mode: Variant type mode, 'SNP' or 'INDEL'.
    :param sites_only_vcf: AS-annotated sites-only input VCF.
    :param features: Features to extract for this mode.
    :param resource_args: GATK ``--resource`` args for the labeled sets.
    :param out_root: Output prefix (files written as ``{out_root}.*``).
    :param gatk_image: GATK docker image.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :return: Job with output ResourceGroup ``extract``.
    """
    j = b.new_job(f"GATK: ExtractVariantAnnotations {mode}")
    j.image(gatk_image)
    j.cpu(4)
    j.memory("standard")
    j.storage("100G")
    j.declare_resource_group(
        extract={
            "annot.hdf5": "{root}.annot.hdf5",
            "vcf.gz": "{root}.vcf.gz",
            "vcf.gz.tbi": "{root}.vcf.gz.tbi",
        }
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options "-Xmx6G" ExtractVariantAnnotations \\
            -V {sites_only_vcf} \\
            -O {j.extract} \\
            --mode {mode} \\
            {_annotation_args(features)} \\
            --use-allele-specific-annotations \\
            --gcs-project-for-requester-pays {gcp_billing_project} \\
            {resource_args}
        """
    )
    b.write_output(j.extract, out_root)
    return j


def train_variant_annotations_model_job(
    b: hb.Batch,
    mode: str,
    annotations_hdf5: hb.ResourceFile,
    out_root: str,
    gatk_image: str,
    gcp_billing_project: str,
    hyperparameters_json: Optional[str] = None,
) -> Job:
    """
    Run GATK TrainVariantAnnotationsModel (PYTHON_IFOREST) for one variant type mode.

    :param b: Batch to add the job to.
    :param mode: Variant type mode, 'SNP' or 'INDEL'.
    :param annotations_hdf5: Extracted annotations HDF5 from the extract step.
    :param out_root: Output prefix (files written as ``{out_root}.{mode}.*``).
    :param gatk_image: GATK docker image.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param hyperparameters_json: Optional model hyperparameters JSON path.
    :return: Job with output ResourceGroup ``model``.
    """
    # GATK writes model files with a lowercase variant-type infix (e.g.
    # .snp.scorer.pkl).
    m = mode.lower()
    j = b.new_job(f"GATK: TrainVariantAnnotationsModel {mode}")
    j.image(gatk_image)
    j.cpu(4)
    j.memory("standard")
    j.storage("100G")
    j.declare_resource_group(
        model={
            f"{m}.scorer.pkl": f"{{root}}.{m}.scorer.pkl",
            f"{m}.trainingScores.hdf5": f"{{root}}.{m}.trainingScores.hdf5",
            f"{m}.calibrationScores.hdf5": f"{{root}}.{m}.calibrationScores.hdf5",
        }
    )
    hp_arg = (
        f" --hyperparameters-json {hyperparameters_json}"
        if hyperparameters_json
        else ""
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options "-Xmx6G" TrainVariantAnnotationsModel \\
            --annotations-hdf5 {annotations_hdf5} \\
            -O {j.model} \\
            --mode {mode} \\
            --model-backend PYTHON_IFOREST \\
            --gcs-project-for-requester-pays {gcp_billing_project}{hp_arg}
        """
    )
    b.write_output(j.model, out_root)
    return j


def split_intervals_job(
    b: hb.Batch,
    calling_intervals_arg: str,
    exclude_intervals: str,
    scatter_count: int,
    gatk_image: str,
    gcp_billing_project: str,
) -> Job:
    """
    Split the calling contigs into ``scatter_count`` intervals, excluding telomeres/centromeres.

    :param b: Batch to add the job to.
    :param calling_intervals_arg: GATK ``-L`` args covering the calling contigs.
    :param exclude_intervals: Path to an intervals file to exclude (``-XL``).
    :param scatter_count: Number of interval shards to produce.
    :param gatk_image: GATK docker image.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :return: Job with output ResourceGroup ``intervals``.
    """
    j = b.new_job(f"GATK: SplitIntervals ({scatter_count})")
    j.image(gatk_image)
    j.memory("standard")
    j.storage("16G")
    j.declare_resource_group(
        intervals={
            f"interval_{idx}": f"{{root}}/{str(idx).zfill(4)}-scattered.interval_list"
            for idx in range(scatter_count)
        }
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options "-Xms3g" SplitIntervals \\
            -R {REFERENCE_RESOURCES['ref_fasta']} \\
            {calling_intervals_arg} \\
            -XL {exclude_intervals} \\
            -O {j.intervals} \\
            -scatter {scatter_count} \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            -mode INTERVAL_SUBDIVISION \\
            --gcs-project-for-requester-pays {gcp_billing_project}
        """
    )
    return j


def score_variant_annotations_job(
    b: hb.Batch,
    mode: str,
    sites_only_vcf: str,
    features: List[str],
    extracted_vcf: hb.ResourceFile,
    model: hb.ResourceGroup,
    resource_args: str,
    interval: hb.ResourceFile,
    out_root: str,
    gatk_image: str,
    gcp_billing_project: str,
) -> Job:
    """
    Run GATK ScoreVariantAnnotations for one variant type mode over one interval.

    :param b: Batch to add the job to.
    :param mode: Variant type mode, 'SNP' or 'INDEL'.
    :param sites_only_vcf: AS-annotated sites-only input VCF.
    :param features: Features for this mode (must match the trained model).
    :param extracted_vcf: Extracted training/calibration VCF from the extract step.
    :param model: Trained model ResourceGroup (its root is the ``--model-prefix``).
    :param resource_args: GATK ``--resource`` args for the labeled sets.
    :param interval: Interval file to restrict scoring to.
    :param out_root: Output prefix (files written as ``{out_root}.*``).
    :param gatk_image: GATK docker image.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :return: Job with output ResourceGroup ``output_score``.
    """
    j = b.new_job(f"GATK: ScoreVariantAnnotations {mode}")
    j.image(gatk_image)
    j.cpu(4)
    j.memory("standard")
    j.storage("100G")
    j.declare_resource_group(
        output_score={
            "vcf.gz": "{root}.vcf.gz",
            "vcf.gz.tbi": "{root}.vcf.gz.tbi",
            "annot.hdf5": "{root}.annot.hdf5",
            "scores.hdf5": "{root}.scores.hdf5",
        }
    )
    j.command(
        f"""set -euo pipefail
        gatk --java-options "-Xmx6G" ScoreVariantAnnotations \\
            -V {sites_only_vcf} \\
            -O {j.output_score} \\
            --mode {mode} \\
            {_annotation_args(features)} \\
            --resource:extracted,extracted=true {extracted_vcf} \\
            --model-prefix {model} \\
            --model-backend PYTHON_IFOREST \\
            --use-allele-specific-annotations \\
            -L {interval} \\
            --gcs-project-for-requester-pays {gcp_billing_project} \\
            {resource_args}
        """
    )
    b.write_output(j.output_score, out_root)
    return j


def isolation_forest_workflow(
    b: hb.Batch,
    sites_only_vcf: str,
    run_prefix: str,
    out_vcf_name: str,
    singletons_vcf: Optional[str],
    calling_intervals_arg: str,
    exclude_intervals: str,
    scatter_count: int,
    gatk_image: str,
    gcp_billing_project: str,
    hyperparameters_json: Optional[str] = None,
    overwrite: bool = False,
) -> None:
    """
    Add the per-mode extract/train/score jobs (SNP and INDEL) to the Batch.

    Extract and train outputs are reused if they already exist (unless ``overwrite``).

    :param b: Batch to add jobs to.
    :param sites_only_vcf: AS-annotated sites-only input VCF.
    :param run_prefix: GCS prefix for this run's GATK outputs.
    :param out_vcf_name: Base name for scored VCF shards.
    :param singletons_vcf: Optional true-positive singletons VCF.
    :param calling_intervals_arg: GATK ``-L`` args for the calling contigs.
    :param exclude_intervals: Path to the telomere/centromere exclusion intervals file.
    :param scatter_count: Number of interval shards for scoring.
    :param gatk_image: GATK docker image.
    :param gcp_billing_project: GCP billing project for requester-pays buckets.
    :param hyperparameters_json: Optional model hyperparameters JSON path.
    :param overwrite: Whether to rerun extract/train even if outputs exist.
    :return: None.
    """
    intervals = split_intervals_job(
        b=b,
        calling_intervals_arg=calling_intervals_arg,
        exclude_intervals=exclude_intervals,
        scatter_count=scatter_count,
        gatk_image=gatk_image,
        gcp_billing_project=gcp_billing_project,
    ).intervals

    # GATK requires an uppercase --mode; IF_FEATURES is keyed 'snv'/'indel'.
    for mode, feature_key in [("SNP", "snv"), ("INDEL", "indel")]:
        m = mode.lower()
        features = IF_FEATURES[feature_key]
        resource_args = _resource_args(mode, singletons_vcf)
        extract_root = f"{run_prefix}/extract/{m}/extract"
        model_root = f"{run_prefix}/model/{m}/model"

        # Reuse extract outputs if present.
        if not overwrite and file_exists(f"{extract_root}.annot.hdf5"):
            logger.info("Reusing existing %s extract outputs.", mode)
            annotations_hdf5 = b.read_input(f"{extract_root}.annot.hdf5")
            extracted_vcf = b.read_input_group(
                **{
                    "vcf.gz": f"{extract_root}.vcf.gz",
                    "vcf.gz.tbi": f"{extract_root}.vcf.gz.tbi",
                }
            )["vcf.gz"]
        else:
            extract_j = extract_variant_annotations_job(
                b=b,
                mode=mode,
                sites_only_vcf=sites_only_vcf,
                features=features,
                resource_args=resource_args,
                out_root=extract_root,
                gatk_image=gatk_image,
                gcp_billing_project=gcp_billing_project,
            )
            annotations_hdf5 = extract_j.extract["annot.hdf5"]
            extracted_vcf = extract_j.extract["vcf.gz"]

        # Reuse trained model if present.
        if not overwrite and file_exists(f"{model_root}.{m}.scorer.pkl"):
            logger.info("Reusing existing %s model.", mode)
            model = b.read_input_group(
                **{
                    f"{m}.scorer.pkl": f"{model_root}.{m}.scorer.pkl",
                    f"{m}.trainingScores.hdf5": f"{model_root}.{m}.trainingScores.hdf5",
                    f"{m}.calibrationScores.hdf5": f"{model_root}.{m}.calibrationScores.hdf5",
                }
            )
        else:
            model = train_variant_annotations_model_job(
                b=b,
                mode=mode,
                annotations_hdf5=annotations_hdf5,
                out_root=model_root,
                gatk_image=gatk_image,
                gcp_billing_project=gcp_billing_project,
                hyperparameters_json=hyperparameters_json,
            ).model

        for idx in range(scatter_count):
            score_variant_annotations_job(
                b=b,
                mode=mode,
                sites_only_vcf=sites_only_vcf,
                features=features,
                extracted_vcf=extracted_vcf,
                model=model,
                resource_args=resource_args,
                interval=intervals[f"interval_{idx}"],
                out_root=f"{run_prefix}/score/{m}/{out_vcf_name}{idx}",
                gatk_image=gatk_image,
                gcp_billing_project=gcp_billing_project,
            )


def export_exclusion_intervals() -> str:
    """
    Export telomere/centromere intervals to a GATK-readable intervals file.

    :return: Path to the exported ``chr:start-end`` intervals file.
    """
    ht = telomeres_and_centromeres.ht()
    interval = ht.interval
    ht = ht.key_by(
        interval=interval.start.contig
        + ":"
        + hl.str(interval.start.position)
        + "-"
        + hl.str(interval.end.position)
    ).select()
    path = hl.utils.new_temp_file("telomeres_centromeres", "intervals")
    ht.export(path, header=False)
    return path


def merge_iforest_result(
    run_prefix: str,
    out_vcf_name: str,
    model_id: str,
    scatter_count: int,
    n_partitions: int,
    header_path: Optional[str],
    array_elements_required: bool,
) -> hl.Table:
    """
    Merge the per-mode scored VCFs into a single variant QC result HT.

    The SNP-pass VCF contributes SNP alleles and the INDEL-pass VCF contributes all
    non-SNP alleles, so the union covers every scored allele exactly once.

    Only the current run's ``scatter_count`` shards are read (by explicit index), so a
    rerun with a smaller scatter count does not pick up stale shards from a prior run.

    :param run_prefix: GCS prefix for this run's GATK outputs.
    :param out_vcf_name: Base name used for scored VCF shards.
    :param model_id: Isolation forest model ID (starts with ``if_``).
    :param scatter_count: Number of score shards per mode (must match the scoring run).
    :param n_partitions: Number of partitions for the imported HTs.
    :param header_path: Optional VCF header file for import.
    :param array_elements_required: Value passed to hl.import_vcf.
    :return: Merged variant QC result HT.
    """
    snp_vcfs = [
        f"{run_prefix}/score/snp/{out_vcf_name}{idx}.vcf.gz"
        for idx in range(scatter_count)
    ]
    indel_vcfs = [
        f"{run_prefix}/score/indel/{out_vcf_name}{idx}.vcf.gz"
        for idx in range(scatter_count)
    ]
    snp_ht = import_variant_qc_vcf(
        snp_vcfs,
        model_id,
        n_partitions,
        header_path,
        array_elements_required,
        deduplicate_check=True,
    )[0]
    indel_ht = import_variant_qc_vcf(
        indel_vcfs,
        model_id,
        n_partitions,
        header_path,
        array_elements_required,
        deduplicate_check=True,
    )[0]

    snp_ht = snp_ht.filter(hl.is_snp(snp_ht.alleles[0], snp_ht.alleles[1]))
    indel_ht = indel_ht.filter(~hl.is_snp(indel_ht.alleles[0], indel_ht.alleles[1]))
    ht = snp_ht.union(indel_ht)
    return ht.annotate_globals(
        model_id=model_id,
        snp_features=IF_FEATURES["snv"],
        indel_features=IF_FEATURES["indel"],
    )


def main(args):
    """Run the isolation forest variant QC workflow."""
    environment = args.environment
    _init_hail("isolation_forest", environment)

    test = args.test
    model_id = args.model_id
    scatter_count = 10 if test else args.scatter_count
    contigs = ["chr22"] if test else CALLING_CONTIGS
    calling_intervals_arg = " ".join(f"-L {c}" for c in contigs)

    true_positive_type = None
    if args.transmitted_singletons and args.sibling_singletons:
        true_positive_type = "transmitted_singleton.sibling_singleton"
    elif args.transmitted_singletons:
        true_positive_type = "transmitted_singleton"
    elif args.sibling_singletons:
        true_positive_type = "sibling_singleton"
    # --test-on-v4 swaps inputs to v4 resources to validate the Batch pipeline while v5
    # data is unavailable; outputs still route via --environment.
    if args.test_on_v4:
        from gnomad_qc.v4.resources.annotations import (
            get_true_positive_vcf_path as v4_get_true_positive_vcf_path,
        )
        from gnomad_qc.v4.resources.annotations import info_vcf_path as v4_info_vcf_path

        sites_only_vcf = v4_info_vcf_path(info_method="quasi")
        singletons_vcf = (
            v4_get_true_positive_vcf_path(
                adj=args.adj, true_positive_type=true_positive_type
            )
            if true_positive_type
            else None
        )
        run_prefix = f"gs://{BATCH_TMP_BUCKET}/if_v4_test/{model_id}"
        result_ht_path = f"{run_prefix}/gnomad.exomes.v4.0.{model_id}.result.ht"
    else:
        sites_only_vcf = info_vcf_path(test=test, environment=environment)
        singletons_vcf = (
            get_true_positive_vcf_path(
                adj=args.adj,
                true_positive_type=true_positive_type,
                test=test,
                environment=environment,
            )
            if true_positive_type
            else None
        )
        run_prefix = get_iforest_run_prefix(
            model_id, test=test, environment=environment
        )
        result_ht_path = get_variant_qc_result(
            model_id, test=test, split=True, environment=environment
        ).path

    if not args.load_only:
        exclude_intervals = export_exclusion_intervals()
        backend = hb.ServiceBackend(
            billing_project=args.batch_billing_project,
            remote_tmpdir=f"gs://{BATCH_TMP_BUCKET}/",
        )
        b = hb.Batch(f"isolation forest {model_id}{args.batch_suffix}", backend=backend)
        isolation_forest_workflow(
            b=b,
            sites_only_vcf=sites_only_vcf,
            run_prefix=run_prefix,
            out_vcf_name=args.out_vcf_name,
            singletons_vcf=singletons_vcf,
            calling_intervals_arg=calling_intervals_arg,
            exclude_intervals=exclude_intervals,
            scatter_count=scatter_count,
            gatk_image=args.gatk_image,
            gcp_billing_project=args.gcp_billing_project,
            hyperparameters_json=args.hyperparameters_json,
            overwrite=args.overwrite,
        )
        b.run()

    if args.load_iforest:
        ht = merge_iforest_result(
            run_prefix=run_prefix,
            out_vcf_name=args.out_vcf_name,
            model_id=model_id,
            scatter_count=scatter_count,
            n_partitions=args.n_partitions,
            header_path=args.header_path,
            array_elements_required=args.array_elements_required,
        )
        ht.write(result_ht_path, overwrite=args.overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Whether to overwrite existing outputs.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Filter to chr22 and use a small scatter count for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--test-on-v4",
        help="Use gnomAD v4 sites/true-positive VCFs as inputs (v5 data unavailable).",
        action="store_true",
    )
    parser.add_argument(
        "--environment",
        help="Compute environment.",
        default="batch",
        type=str,
        choices=["rwb", "batch", "dataproc"],
    )
    parser.add_argument(
        "--model-id",
        help="Model ID for the isolation forest run. Must start with 'if_'.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--out-vcf-name",
        help="Base name for scored VCF shards.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--batch-billing-project",
        help="Hail Batch billing project.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--gcp-billing-project",
        help="GCP billing project for requester-pays buckets.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--gatk-image",
        help="GATK docker image.",
        default=DEFAULT_GATK_IMAGE,
        type=str,
    )
    parser.add_argument(
        "--scatter-count",
        help="Number of intervals to scatter scoring across.",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--transmitted-singletons",
        help="Include transmitted singletons as a training/calibration set.",
        action="store_true",
    )
    parser.add_argument(
        "--sibling-singletons",
        help="Include sibling singletons as a training/calibration set.",
        action="store_true",
    )
    parser.add_argument(
        "--adj",
        help="Use adj genotypes for the true-positive singletons VCF.",
        action="store_true",
    )
    parser.add_argument(
        "--hyperparameters-json",
        help="Optional GATK model hyperparameters JSON path.",
        type=str,
    )
    parser.add_argument(
        "--batch-suffix",
        help="String to append to the Batch name.",
        default="",
        type=str,
    )
    parser.add_argument(
        "--load-iforest",
        help="Merge the scored VCFs into a variant QC result HT.",
        action="store_true",
    )
    parser.add_argument(
        "--load-only",
        help="Skip the GATK Batch and only run the load step.",
        action="store_true",
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of partitions for the result HT.",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--header-path",
        help="Optional header file to use when importing the scored VCFs.",
    )
    parser.add_argument(
        "--array-elements-required",
        help="Set array_elements_required=True when importing the scored VCFs.",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
