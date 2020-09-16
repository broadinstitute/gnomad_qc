import hail as hl
from gnomad.resources.resource_utils import (
    TableResource,
    PedigreeResource,
    VersionedPedigreeResource,
)
from gnomad.utils.file_utils import file_exists
from gnomad_qc.v3.resources import (
    CURRENT_PROJECT_META_VERSION,
    CURRENT_META_VERSION,
    CURRENT_RELEASE,
)


# Samples metadata
def get_meta_root(version: str = CURRENT_RELEASE) -> str:
    """
    Retrieves the path to the root metadata directory

    :param version: gnomAD release version
    :return: String representation of the path to the root metadata directory
    """
    return f"gs://gnomad/metadata/genomes_v{version}"


def meta(
    version: str = CURRENT_RELEASE, meta_version: str = CURRENT_META_VERSION
) -> TableResource:
    """
    Gets the TableResource for the finalized sample metadata information after sample QC

    :param version: gnomAD release version
    :param meta_version: metadata version to return
    :return: Table containing the finalized metadata
    """
    return TableResource(
        f"{get_meta_root(version)}/gnomad_v{version}_metadata_{meta_version}.ht"
    )


def meta_tsv_path(
    version: str = CURRENT_RELEASE, meta_version: str = CURRENT_META_VERSION
) -> str:
    """
    Gets the path to the finalized sample metadata information after sample QC

    :param version: gnomAD release version
    :param meta_version: metadata version to return
    :return: String path to the finalized metadata
    """
    return f"{get_meta_root(version)}/gnomad_v{version}_metadata_{meta_version}.tsv.gz"


def project_meta(
    version: str = CURRENT_RELEASE,
    proj_meta_version: str = CURRENT_PROJECT_META_VERSION,
    min_partitions: int = 100,
) -> TableResource:
    """
    Retrieves the project level matadata

    :param version: gnomAD release version
    :param proj_meta_version:
    :param min_partitions:
    :return: Table containing project level metadata
    """
    meta_path = f"{get_meta_root(version)}/{proj_meta_version}_v{version}_project_meta.txt"
    if file_exists(meta_path):
        return TableResource(
            import_func=hl.import_table,
            import_args={
                "paths": meta_path,
                "impute": True,
                "key": "s",
                "min_partitions": min_partitions,
            },
        )
    else:
        return TableResource(f"{get_meta_root(version)}/{proj_meta_version}_v{version}_project_meta.ht")


def pedigree(version: str = CURRENT_RELEASE) -> VersionedPedigreeResource:
    """

    :param version: gnomAD release version
    :return:
    """
    return VersionedPedigreeResource(
        "final",  # TODO: Make sure "final" is the best label once the family scripts are in
        {
            "raw": PedigreeResource(
                f"{get_meta_root(version)}/gnomad_v{version}_raw.fam", delimiter="\t"
            ),
            "final": PedigreeResource(
                f"{get_meta_root(version)}/gnomad_v{version}.fam", delimiter="\t"
            ),
        },
    )


def trios(version: str = CURRENT_RELEASE) -> VersionedPedigreeResource:
    """

    :param version: gnomAD release version
    :return:
    """
    return VersionedPedigreeResource(  # TODO: Should this be merged with Pedigree into a single resource?
        "final",  # TODO: Make sure "final" is the best label once the family scripts are in
        {
            "raw": PedigreeResource(
                f"{get_meta_root(version)}/gnomad_v{version}_trios_raw.fam"
            ),
            "final": PedigreeResource(
                f"{get_meta_root(version)}/gnomad_v{version}_trios.fam"
            ),
        },
    )


def ped_mendel_errors(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: gnomAD release version
    :return:
    """
    return TableResource(
        f"{get_meta_root(version)}/gnomad_v{version}_ped_chr20_mendel_errors.ht"
    )
