"""Script containing generic resources."""

from gnomad.resources.resource_utils import VariantDatasetResource

# v5 DRAGEN TGP test VDS.
dragen_tgp_vds = VariantDatasetResource(
    "gs://gnomad/v5.0/testing/genomes/dragen_tgp_v5.0_test.vds"
)


# TODO: Change default to CURRENT_VERSION.
def qc_temp_prefix(version: str = "5.0") -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return.
    :return: Path to bucket with temporary QC data.
    """
    return f"gs://gnomad-tmp/gnomad.genomes.v{version}.qc_data/"


# TODO: Change default to CURRENT_VERSION.
def get_logging_path(name: str, version: str = "5.0") -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :param version: Version of annotation path to return.
    :return: Output log path.
    """
    return f"{qc_temp_prefix(version)}{name}.log"
