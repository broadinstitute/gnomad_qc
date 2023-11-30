"""Script containing metadata related resources."""
from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.constants import (
    CURRENT_PROJECT_META_VERSION,
    CURRENT_SAMPLE_QC_VERSION,
)


def _meta_root_path(
    version: str = CURRENT_PROJECT_META_VERSION, data_type: str = "exomes"
) -> str:
    """
    Retrieve the path to the root metadata directory.

    :param version: gnomAD version.
    :param data_type: Data type ("exomes" or "genomes"). Default is "exomes".
    :return: String representation of the path to the root metadata directory.
    """
    return f"gs://gnomad/v{version}/metadata/{data_type}"


def meta_tsv_path(
    version: str = CURRENT_SAMPLE_QC_VERSION, data_type: str = "exomes"
) -> str:
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


def meta(
    version: str = CURRENT_SAMPLE_QC_VERSION,
    data_type: str = "exomes",
) -> TableResource:
    """
    Get the gnomAD v4 sample QC meta VersionedTableResource.

    Function will check that metadata exists for the requested gnomAD version and data type
    and will throw an error if it doesn't exist.

    :param version: gnomAD version.
    :param data_type: Data type ("exomes" or "genomes"). Default is "exomes".
    :return: gnomAD v4 sample meta VersionedTableResource.
    """
    # Check if the meta Table exists for specified version and data type combination
    check_resource_existence(
        input_step_resources={
            f"{data_type}_{version}": [
                f"{_meta_root_path(version, data_type)}/gnomad.{data_type}.v{version}.sample_qc_metadata.ht"
            ],
        },
        overwrite=False,
    )

    return TableResource(
        path=(
            f"{_meta_root_path(version, data_type)}/gnomad.{data_type}.v{version}.sample_qc_metadata.ht"
        )
    )


# This isn't used in `meta()` function anymore but is kept here to quickly get the
# currently existing versions of sample metadata.
_meta_versions = {
    "exomes": {
        "4.0": TableResource(
            path="gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.sample_qc_metadata.ht"
        ),
    },
    "genomes": {
        "4.0": TableResource(
            path="gs://gnomad/v4.0/metadata/genomes/gnomad.genomes.v4.0.sample_qc_metadata.ht"
        ),
    },
}

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

project_meta = VersionedTableResource(
    CURRENT_PROJECT_META_VERSION, _project_meta_versions
)
picard_metrics = VersionedTableResource(
    CURRENT_PROJECT_META_VERSION, _picard_metric_versions
)
gatk_versions = VersionedTableResource(CURRENT_PROJECT_META_VERSION, _gatk_versions)
