from typing import Optional

from gnomad.resources.resource_utils import (
    MatrixTableResource,
    TableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
)

from gnomad_qc.v3.resources.constants import CURRENT_RELEASE, RELEASES


def qual_hists_json_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Fetch filepath for qual histograms JSON.

    :param release_version: Release version. Defualts to CURRENT RELEASE
    :return: File path for histogram JSON
    :rtype: str
    """
    return f"gs://gnomad/release/{release_version}/json/gnomad.genomes.r{release_version}.json"


# TODO: Remove if not used after all python files are in
# internal_ht_path = 'gs://gnomad/release/3.0/ht/gnomad.genomes.r3.0.nested.no_subsets.sites.ht'


def release_ht_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables
    
    :param data_type: 'exomes' or 'genomes'
    :param release_version: release version
    :param public: Whether to return the desired
    :return: File path for desired Hail Table
    :rtype: str
    """
    if public:
        return f"gs://gnomad-public/release/{release_version}/ht/{data_type}/gnomad.{data_type}.r{release_version}.sites.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/gnomad.{data_type}.r{release_version}.sites.ht"


def release_subset(
    subset: str, dense: bool = False, data_type: str = "genomes",
) -> VersionedMatrixTableResource:
    """
    Get a subset release MatrixTableResource.
    
    :param subset: One of the possible release subsets (e.g., hgdp_1kg)
    :param dense: If True, return the dense MT; if False, return the sparse MT
    :param data_type: 'exomes' or 'genomes'
    :return: MatrixTableResource for specific subset
    """

    return VersionedMatrixTableResource(
        CURRENT_RELEASE,
        {
            release: MatrixTableResource(
                f"gs://gnomad/release/{release}/mt/gnomad.{data_type}.v{release}.{subset}_subset{f'_dense' if dense else '_sparse'}.mt"
            )
            for release in RELEASES
            if release != "3"
        },
    )


def release_subset_annotations(
    subset: str, data_type: str = "genomes", sample: bool = True,
) -> VersionedTableResource:
    """
    Get the subset release sample or variant TableResource.
    
    :param subset: One of the possible release subsets (e.g., hgdp_1kg)
    :param data_type: 'exomes' or 'genomes'
    :param sample: If true, will return the sample annotations, otherwise will return the variant annotations
    :return: Table resource with sample/variant annotations for the subset
    """
    return VersionedTableResource(
        CURRENT_RELEASE,
        {
            release: TableResource(
                f"gs://gnomad/release/{release}/ht/gnomad.{data_type}.v{release}.{subset}_subset{f'_sample_meta' if sample else '_variant_annotations'}.ht"
            )
            for release in RELEASES
            if release != "3"
        },
    )


def release_subset_sample_tsv(subset: str, release: str = CURRENT_RELEASE, data_type: str = "genomes") -> str:
    """
    Get the path to the subset release sample annotation text file.

    :param subset: One of the possible release subsets (e.g., hgdp_1kg)
    :param release: Version of annotation tsv path to return
    :param data_type: 'exomes' or 'genomes'
    :return: Path to file
    """
    return f"gs://gnomad/release/{release}/tsv/gnomad.{data_type}.v{release}.{subset}_subset_sample_meta.tsv.bgz"
