"""Script containing variant QC related resources for v5."""

from gnomad.resources.grch38 import (
    na12878_giab,
    na12878_giab_hc_intervals,
    syndip,
    syndip_hc_intervals,
)
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

VQSR_FEATURES = {
    "genomes": {
        "snv": [
            "AS_QD",
            "AS_MQRankSum",
            "AS_ReadPosRankSum",
            "AS_FS",
            "AS_SOR",
            "AS_MQ",
        ],
        "indel": [
            "AS_QD",
            "AS_MQRankSum",
            "AS_ReadPosRankSum",
            "AS_FS",
            "AS_SOR",
        ],
    },
}
"""List of features used in the VQSR model."""

SYNDIP = "CHMI_CHMI3_Nex1"
"""
String representation for syndip truth sample.
"""

NA12878 = "ASC-4Set-1573S_NA12878@1075619236"
"""
String representation for NA12878 truth sample.
"""
# TODO: Is this present ? Likely not but check.
UKB_NA12878 = "Coriell_NA12878_NA12878"
"""
String representation for the UKB Regeneron generated NA12878 truth sample.
"""

TRUTH_SAMPLES = {
    "syndip": {"s": SYNDIP, "truth_mt": syndip, "hc_intervals": syndip_hc_intervals},
    "NA12878": {
        "s": NA12878,
        "truth_mt": na12878_giab,
        "hc_intervals": na12878_giab_hc_intervals,
    },
    "UKB_NA12878": {
        "s": UKB_NA12878,
        "truth_mt": na12878_giab,
        "hc_intervals": na12878_giab_hc_intervals,
    },
}
"""
Dictionary containing necessary information for truth samples

Current truth samples available are syndip and NA12878. Available data for each are the
following:

    - s: Sample name in the callset
    - truth_mt: Truth sample MatrixTable resource
    - hc_intervals: High confidence interval Table resource in truth sample
"""


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
# Variant QC result resources
######################################################################


def get_variant_qc_result(
    model_id: str,
    test: bool = False,
    split: bool = True,
) -> VersionedTableResource:
    r"""
    Get the results of variant QC filtering for a given run.

    :param model_id: Model ID of variant QC run to load. Must start with 'rf\_',
        'vqsr\_', or 'if\_'.
    :param test: Whether to use a tmp path for variant QC tests.
    :param split: Whether to return the split (biallelic) result. Set to False for the
        multi-allelic result. Default is True.
    :return: VersionedTableResource for variant QC results.
    """
    model_type = model_id.split("_")[0]
    if model_type not in ["rf", "vqsr", "if"]:
        raise ValueError(
            f"Model ID must start with 'rf_', 'vqsr_', or 'if_', but found {model_id}"
        )
    return VersionedTableResource(
        CURRENT_VARIANT_QC_RESULT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test)}/{model_type}/models/{model_id}/aou.genomes.v{version}.{model_type}_result{'' if split else '.unsplit'}.ht"
            )
            for version in VARIANT_QC_RESULT_VERSIONS
        },
    )
