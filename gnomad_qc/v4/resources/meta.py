"""Script containing metadata related resources."""
from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from gnomad_qc.v4.resources.constants import CURRENT_VERSION, VERSIONS


def _meta_root_path(version: str = CURRENT_VERSION, data_type: str = "exomes") -> str:
    """
    Retrieve the path to the root metadata directory.

    :param version: gnomAD version.
    :param data_type: Data type ("exomes" or "genomes"). Default is "exomes".
    :return: String representation of the path to the root metadata directory.
    """
    return f"gs://gnomad/v{version}/metadata/{data_type}"


def meta_tsv_path(version: str = CURRENT_VERSION, data_type: str = "exomes") -> str:
    """
    Get the path to the finalized sample metadata information after sample QC.

    :param version: gnomAD version.
    :param data_type: Data type ("exomes" or "genomes"). Default is "exomes".
    :return: String path to the finalized metadata TSV.
    """
    return (
        f"{_meta_root_path(version, data_type)}/gnomad.{data_type}."
        f"v{version}.metadata.tsv.gz"
    )


def meta(data_type: str = "exomes") -> VersionedTableResource:
    """
    Get the gnomAD v4 sample QC meta VersionedTableResource.

    :param data_type: Data type ("exomes" or "genomes"). Default is "exomes".
    :return: gnomAD v4 sample meta VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_meta_root_path(version, data_type)}/gnomad.{data_type}.v{version}.sample_qc_metadata.ht"
                )
            )
            for version in VERSIONS
        },
    )


# The following variables only apply to "exomes" data.
_project_meta_versions = {
    "4.0": TableResource(
        path="gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.project_meta.ht"
    )
}

_picard_metric_versions = {
    "4.0": TableResource(
        path="gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.picard_metrics.ht"
    )
}

_gatk_versions = {
    "4.0": TableResource(
        path="gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.gatk_versions.ht"
    )
}

project_meta = VersionedTableResource(CURRENT_VERSION, _project_meta_versions)
picard_metrics = VersionedTableResource(CURRENT_VERSION, _picard_metric_versions)
gatk_versions = VersionedTableResource(CURRENT_VERSION, _gatk_versions)
