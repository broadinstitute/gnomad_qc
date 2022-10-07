import hail as hl
from gnomad.resources.resource_utils import (
    PedigreeResource,
    TableResource,
    VersionedPedigreeResource,
    VersionedTableResource,
)

from gnomad_qc.v4.resources.constants import CURRENT_VERSION, VERSIONS


# Samples metadata
def _meta_root_path(version: str = CURRENT_VERSION) -> str:
    """
    Retrieve the path to the root metadata directory.

    :param version: gnomAD version
    :return: String representation of the path to the root metadata directory
    """
    return f"gs://gnomad/v{version}/metadata/exomes"


def meta_tsv_path(version: str = CURRENT_VERSION) -> str:
    """
    Get the path to the finalized sample metadata information after sample QC.

    :param version: gnomAD version
    :return: String path to the finalized metadata
    """
    return f"{_meta_root_path(version)}/gnomad.exomes.v{version}.metadata.tsv.gz"


_meta_versions = {
    "4.0": TableResource(
        path="gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.sample_qc_metadata.ht"
    ),
}

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

_pedigree_versions = {
    f"4{x}": PedigreeResource(
        f"gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0{x}.fam",
        delimiter="\t",
    )
    for x in ["", "_raw"]
}


_trios_versions = {
    f"4{x}": PedigreeResource(
        f"gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.trios{x}.fam",
        delimiter="\t",
    )
    for x in ["", "_raw"]
}


meta = VersionedTableResource(CURRENT_VERSION, _meta_versions)
project_meta = VersionedTableResource(CURRENT_VERSION, _project_meta_versions)
picard_metrics = VersionedTableResource(CURRENT_VERSION, _picard_metric_versions)
pedigree = VersionedPedigreeResource("4", _pedigree_versions)
trios = VersionedPedigreeResource("4", _trios_versions)
ped_mendel_errors = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            path=f"{_meta_root_path(version)}/gnomad.exomes.v{version}.ped_chr20_mendel_errors.ht"
        )
        for version in VERSIONS
    },
)
