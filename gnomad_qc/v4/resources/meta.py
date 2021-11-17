import hail as hl
from gnomad.resources.resource_utils import (
    PedigreeResource,
    TableResource,
    VersionedPedigreeResource,
    VersionedTableResource,
)

from gnomad_qc.v4.resources.constants import (
    CURRENT_VERSION,
    CURRENT_RELEASE,
    RELEASES,
)


# Samples metadata
def _meta_root_path(version: str = CURRENT_RELEASE) -> str:
    """
    Retrieve the path to the root metadata directory.

    :param version: gnomAD release version
    :return: String representation of the path to the root metadata directory
    """
    return f"gs://gnomad/v{version}/metadata/exomes"


def meta_tsv_path(
    version: str = CURRENT_RELEASE, meta_version: str = CURRENT_VERSION
) -> str:
    """
    Get the path to the finalized sample metadata information after sample QC.

    :param version: gnomAD release version
    :param meta_version: Metadata version to return
    :return: String path to the finalized metadata
    """
    return (
        f"{_meta_root_path(version)}/gnomad_exomes_v{version}_metadata_v{meta_version}.tsv.gz"
    )


_meta_versions = {
    "4": TableResource(
        path="gs://gnomad/v4/metadata/exomes/gnomad_v4_sample_qc_metadata.ht"
    ),
}

_project_meta_versions = {
    "4": TableResource(path="gs://gnomad/v4/metadata/exomes/gnomad_v4_project_meta.ht")
}

_pedigree_versions = {
    f"4{x}": PedigreeResource(
        f"gs://gnomad/v4/metadata/exomes/gnomad_exomes_v4{x}.fam", delimiter="\t",
    ),
    for x in ["", "_raw"]
}


_trios_versions = {
    f"4{x}": PedigreeResource(
        f"gs://gnomad/v4/metadata/exomes/gnomad_exomes_v4_trios{x}.fam", delimiter="\t",
        for x in ["", "_raw"]
    ),
}


meta = VersionedTableResource(CURRENT_VERSION, _meta_versions)
project_meta = VersionedTableResource(CURRENT_VERSION, _project_meta_versions)
pedigree = VersionedPedigreeResource("4", _pedigree_versions)
trios = VersionedPedigreeResource("4", _trios_versions)
ped_mendel_errors = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            path=f"{_meta_root_path(release)}/gnomad_exomes_v{release}_ped_chr20_mendel_errors.ht"
        )
        for release in RELEASES
    },
)
