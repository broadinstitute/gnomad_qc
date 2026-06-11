"""Script containing variant QC related resources for v5."""

from gnomad.resources.resource_utils import VariantDatasetResource

from gnomad_qc.v5.resources.basics import (
    _ALL_ENVIRONMENTS,
    _get_base_bucket,
    _validate_environment,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.constants import AOU_BUCKET, CURRENT_VARIANT_QC_VERSION


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


truth_samples_aou_gvcfs_bucket = (
    f"gs://{AOU_BUCKET}/wgs/short_read/snpindel/aux/qc/control_samples/"
)
"""
AoU bucket containing gVCFs from 8 Genomes-in-a-Bottle (GiaB) samples sequenced with the same protocol as the AoU v8 data.

Stored for reference but not used directly.
"""

truth_samples_gvcf_paths = f"{_variant_qc_root(environment='batch')}/aou/truth_samples/truth_samples_gvcf_paths.tsv"
"""
Path to a single-column TSV listing the GCS path to each truth-sample gVCF (one per line).

The sample IDs are intentionally not stored in this repo, and the truth-sample bucket cannot be listed
The combiner therefore reads the gVCF paths from this manifest by known object paths rather than globbing the bucket.å
"""

truth_samples_vds = VariantDatasetResource(
    f"{_variant_qc_root(test=False, environment='batch')}/aou/truth_samples/truth_samples.vds",
)
"""
VDS containing 8 GiaB samples.
"""


def get_truth_samples_combiner_plan(test: bool = False) -> str:
    """
    Return the path to the truth-samples VDS combiner plan (combiner ``save_path``).

    The plan lets a failed or interrupted combiner run be resumed. It is written to the
    writable batch variant QC tree (or a temp path for tests).

    :param test: Whether to return a temporary test path. Default is False.
    :return: Path to the combiner plan JSON.
    """
    if test:
        return f"{qc_temp_prefix(environment='batch', days=4)}truth_samples_combiner_plan.json"

    return f"{_variant_qc_root(environment='batch')}/aou/truth_samples/truth_samples_combiner_plan.json"
