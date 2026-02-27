"""Script containing variant QC related resources.

.. note::

    Sample QC was completed in the Researcher Workbench (RWB), while allele number,
    frequency, trio stats, and variant QC were run in Hail Batch. The default
    environment for variant QC resources is therefore "batch".
"""

from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import CURRENT_VERSION, _get_base_bucket


def _variant_qc_root(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    environment: str = "batch",
) -> str:
    """
    Return path to variant QC root folder.

    :param version: Version of variant QC path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param data_type: Data type, e.g. "exomes" or "genomes". Default is "genomes".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb",
        "batch", or "dataproc".
    :return: Root to variant QC path.
    """
    path_suffix = f"variant_qc/{data_type}"

    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )

    return f"gs://{_get_base_bucket(environment)}/v{version}/{path_suffix}"
