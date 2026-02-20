"""Script containing variant QC related resources."""

from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import (
    BATCH_BUCKET,
    CURRENT_VERSION,
    GNOMAD_BUCKET,
    WORKSPACE_BUCKET,
)


def _variant_qc_root(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    environment: str = "rwb",
) -> str:
    """
    Return path to variant QC root folder.

    :param version: Version of variant QC path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param data_type: Data type, e.g. "exomes" or "genomes". Default is "genomes".
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: Root to variant QC path.
    """
    path_suffix = f"variant_qc/{data_type}"

    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )

    if environment == "rwb":
        base_bucket = WORKSPACE_BUCKET
    elif environment == "batch":
        base_bucket = BATCH_BUCKET
    else:
        base_bucket = GNOMAD_BUCKET
    return f"gs://{base_bucket}/v{version}/{path_suffix}"
