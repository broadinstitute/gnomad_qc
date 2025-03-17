"""Script containing generic resources."""

from gnomad.resources.resource_utils import VariantDatasetResource

from gnomad_qc.v5.resources.constants import CURRENT_VERSION

# v5 DRAGEN TGP test VDS.
dragen_tgp_vds = VariantDatasetResource(
    "gs://gnomad/v5.0/testing/genomes/dragen_tgp_v5.0_test.vds"
)


def qc_temp_prefix(version: str = CURRENT_VERSION) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return.
    :return: Path to bucket with temporary QC data.
    """
    return f"gs://gnomad-tmp/gnomad.genomes.v{version}.qc_data/"


def get_checkpoint_path(
    name: str, version: str = CURRENT_VERSION, mt: bool = False
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable
    :param version: Version of annotation path to return
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    """
    return f'{qc_temp_prefix(version)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(name: str, version: str = CURRENT_VERSION) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :param version: Version of annotation path to return.
    :return: Output log path.
    """
    return f"{qc_temp_prefix(version)}{name}.log"
