"""Script containing variant QC related resources for v5."""

from gnomad_qc.v5.resources.basics import (
    _ALL_ENVIRONMENTS,
    _get_base_bucket,
    _validate_environment,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.constants import CURRENT_VARIANT_QC_VERSION


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
