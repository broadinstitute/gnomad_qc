"""Script containing variant QC related resources for v5."""

from gnomad.resources.resource_utils import (
    TableResource,
    VariantDatasetResource,
    VersionedTableResource,
)

from gnomad_qc.v5.resources.basics import (
    _ALL_ENVIRONMENTS,
    _get_base_bucket,
    _validate_environment,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.constants import (
    CURRENT_VARIANT_QC_RESULT_VERSION,
    CURRENT_VARIANT_QC_VERSION,
    VARIANT_QC_RESULT_VERSIONS,
)

# Isolation forest features, keyed by variant type. GATK trains a separate model per
# --mode, so SNPs/INDELs need their own lists. AS_MQ is dropped for indels (v4 genomes
# VQSR convention); AS_pab_max is kept for both (tree-based, unlike Gaussian VQSR).
IF_FEATURES = {
    "snv": [
        "AS_QD",
        "AS_pab_max",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_FS",
        "AS_SOR",
        "AS_MQ",
    ],
    "indel": [
        "AS_QD",
        "AS_pab_max",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_FS",
        "AS_SOR",
    ],
}
"""Features used in the isolation forest model, keyed by variant type."""


def _variant_qc_root(
    version: str = CURRENT_VARIANT_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    environment: str = "dataproc",
) -> str:
    """
    Return path to variant QC root folder.

    :param version: Version of variant QC path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param data_type: Whether to return 'exomes' or 'genomes' data. Default is genomes.
    :param environment: Environment to use. Default is "dataproc". Must be one of
        "rwb", "batch", or "dataproc".
    :return: Root to variant QC path.
    """
    _validate_environment(environment, _ALL_ENVIRONMENTS)
    path_suffix = f"variant_qc/{data_type}"

    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )

    return f"gs://{_get_base_bucket(environment)}/v{version}/{path_suffix}"


######################################################################
# Truth sample resources
######################################################################

truth_samples_gvcf_paths = f"{_variant_qc_root(version='5.0', environment='batch')}/aou/truth_samples/truth_samples_gvcf_paths.tsv"
"""
Path to a single-column TSV listing the GCS path to each truth-sample gVCF (one per line).

The truth samples are 8 Genomes-in-a-Bottle (GiaB) samples sequenced with the same protocol as the AoU v8 data.
Their gVCFs are stored here: gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux/qc/control_samples/

The sample IDs are intentionally not stored in this repo. The combiner therefore reads the gVCF paths from this
manifest by known object paths rather than globbing the bucket.
"""

truth_samples_vds = VariantDatasetResource(
    f"{_variant_qc_root(version='5.0', test=False, environment='batch')}/aou/truth_samples/truth_samples.vds",
)
"""
VDS containing 8 GiaB samples.

This resource does not need to be remade for future versions.
"""


def get_truth_samples_combiner_plan(test: bool = False) -> str:
    """
    Return the path to the truth-samples VDS combiner plan (combiner ``save_path``).

    The plan lets a failed or interrupted combiner run be resumed. It is written to the
    writable batch variant QC tree (or a temp path for tests).

    Like `truth_samples_vds`, this resource does not need to be remade for future versions.

    :param test: Whether to return a temporary test path. Default is False.
    :return: Path to the combiner plan JSON.
    """
    if test:
        return f"{qc_temp_prefix(environment='batch', days=4)}truth_samples_combiner_plan.json"

    return f"{_variant_qc_root(version='5.0', environment='batch')}/aou/truth_samples/truth_samples_combiner_plan.json"


######################################################################
# Variant QC model / result resources
######################################################################


def _validate_model_id(model_id: str) -> str:
    """
    Validate a variant QC model ID and return its model type prefix.

    :param model_id: Model ID. Must start with ``rf_``, ``vqsr_``, or ``if_``.
    :return: Model type ('rf', 'vqsr', or 'if').
    """
    model_type = model_id.split("_")[0]
    if model_type not in ["rf", "vqsr", "if"]:
        raise ValueError(
            f"Model ID must start with 'rf_', 'vqsr_', or 'if_', but got {model_id}"
        )
    return model_type


def get_variant_qc_result(
    model_id: str,
    test: bool = False,
    split: bool = True,
    data_type: str = "genomes",
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Get the results of variant QC filtering for a given run.

    :param model_id: Model ID of variant QC run to load. Must start with ``rf_``,
        ``vqsr_``, or ``if_``.
    :param test: Whether to use a tmp path for variant QC tests.
    :param split: Whether to return the split or unsplit variant QC result.
    :param data_type: Whether to return 'exomes' or 'genomes' data. Default is genomes.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb",
        "batch", or "dataproc".
    :return: VersionedTableResource for variant QC results.
    """
    model_type = _validate_model_id(model_id)
    return VersionedTableResource(
        CURRENT_VARIANT_QC_RESULT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test, data_type=data_type, environment=environment)}/{model_type}/models/{model_id}/gnomad.{data_type}.v{version}.{model_type}_result{'' if split else '.unsplit'}.ht"
            )
            for version in VARIANT_QC_RESULT_VERSIONS
        },
    )


def get_iforest_run_prefix(
    model_id: str,
    version: str = CURRENT_VARIANT_QC_RESULT_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    environment: str = "batch",
) -> str:
    """
    Get the GCS prefix for an isolation forest run's GATK working/output files.

    The isolation forest workflow writes its GATK intermediate outputs (extracted
    annotations, trained model, scored VCFs) under this prefix, alongside the model's
    result HT from :func:`get_variant_qc_result`.

    :param model_id: Model ID of the isolation forest run. Must start with ``if_``.
    :param version: Version of the variant QC result path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param data_type: Whether to return 'exomes' or 'genomes' data. Default is genomes.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb",
        "batch", or "dataproc".
    :return: GCS prefix (no trailing slash) for the isolation forest run's files.
    """
    if _validate_model_id(model_id) != "if":
        raise ValueError(f"model_id must start with 'if_', but got {model_id}")
    root = _variant_qc_root(
        version, test=test, data_type=data_type, environment=environment
    )
    return f"{root}/if/models/{model_id}"
