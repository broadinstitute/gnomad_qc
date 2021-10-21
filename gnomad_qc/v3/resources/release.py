from typing import Optional

from gnomad.resources.resource_utils import (
    DataException,
    MatrixTableResource,
    TableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
)

from gnomad_qc.v3.resources.basics import qc_temp_prefix
from gnomad_qc.v3.resources.constants import (
    CURRENT_RELEASE,
    CURRENT_HGDP_TGP_RELEASE,
    HGDP_TGP_RELEASES,
    RELEASES,
)


def annotation_hists_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Returns path to file containing ANNOTATIONS_HISTS dictionary.
    Dictionary contains histogram values for each metric.
    For example, "InbreedingCoeff": [-0.25, 0.25, 50].

    :return: Path to file with annotations histograms
    :rtype: str
    """
    return f"gs://gnomad/release/{release_version}/json/annotation_hists.json"


def qual_hists_json_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Fetch filepath for qual histograms JSON.

    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: File path for histogram JSON
    :rtype: str
    """
    return f"gs://gnomad/release/{release_version}/json/gnomad.genomes.r{release_version}.json"


def release_ht_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
    het_nonref_patch: bool = False,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables.

    :param data_type: 'exomes' or 'genomes'
    :param release_version: release version
    :param public: Determines whether release sites Table is read from public or private bucket. Defaults to private
    :param het_nonref_patch: Whether this is frequency information for only variants that need the het nonref patch applied
    :return: File path for desired Hail Table
    :rtype: str
    """
    version_prefix = "r" if release_version.startswith("3.0") else "v"
    if public:
        if het_nonref_patch:
            DataException("The patch HT will not be made public")
        return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.{version_prefix}{release_version}.sites.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/{data_type}/gnomad.{data_type}.{version_prefix}{release_version}{'.patch' if het_nonref_patch else ''}.sites.ht"


def release_sites(
    public: bool = False, het_nonref_patch: bool = False
) -> VersionedTableResource:
    """
    Retrieve versioned resource for sites-only release Table.

    :param public: Determines whether release sites Table is read from public or private bucket. Defaults to private
    :param het_nonref_patch: Whether this is frequency information for only variants that need the het nonref patch applied
    :return: Sites-only release Table
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_ht_path(
                    release_version=release,
                    public=public,
                    het_nonref_patch=het_nonref_patch,
                )
            )
            for release in RELEASES
        },
    )


def release_vcf_path(
    release_version: Optional[str] = None,
    hgdp_tgp_subset: bool = False,
    contig: Optional[str] = None,
) -> str:
    """
    Fetch bucket for release (sites-only) VCFs.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used unless
        hgdp_tgp_subset is True in which case CURRENT_HGDP_TGP_RELEASE is used.
    :param hgdp_tgp_subset: Whether to get path for HGDP + 1KG/TGP VCF. Defaults to the full callset (metrics on all samples) sites VCF path
    :param contig: String containing the name of the desired reference contig. Defaults to the full (all contigs) sites VCF path
        sites VCF path
    :return: Filepath for the desired VCF
    """

    if release_version is None:
        release_version = (
            CURRENT_HGDP_TGP_RELEASE if hgdp_tgp_subset else CURRENT_RELEASE
        )

    if hgdp_tgp_subset:
        if release_version not in HGDP_TGP_RELEASES:
            raise DataException(
                f"{release_version} is not one of the available releases for the HGP + 1KG/TGP subset: {HGDP_TGP_RELEASES}"
            )
        subset = "hgdp_tgp"
    else:
        subset = "sites"

    if contig:
        return f"gs://gnomad/release/{release_version}/vcf/genomes/gnomad.genomes.v{release_version}.{subset}.{contig}.vcf.bgz"
    else:
        # if contig is None, return path to sharded vcf bucket
        # NOTE: need to add .bgz or else hail will not bgzip shards
        return f"gs://gnomad/release/{release_version}/vcf/genomes/gnomad.genomes.v{release_version}.{subset}.vcf.bgz"


def release_header_path(
    release_version: Optional[str] = None, hgdp_tgp_subset: bool = False
) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used unless
        hgdp_tgp_subset is True in which case CURRENT_HGDP_TGP_RELEASE is used.
    :param hgdp_tgp_subset: Whether to return the header for the HGDP + 1KG/TGP subset. Default will return the header
        path for the full release.
    :return: Filepath for header dictionary pickle
    """
    subset = ""
    if release_version is None:
        release_version = (
            CURRENT_HGDP_TGP_RELEASE if hgdp_tgp_subset else CURRENT_RELEASE
        )

    if hgdp_tgp_subset:
        if release_version not in HGDP_TGP_RELEASES:
            raise DataException(
                f"{release_version} is not one of the available releases for the HGP + 1KG/TGP subset: {HGDP_TGP_RELEASES}"
            )
        subset = "_hgdp_tgp"

    return f"gs://gnomad/release/{release_version}/vcf/genomes/gnomad.genomes.v{release_version}_header_dict{subset}.pickle"


def append_to_vcf_header_path(
    subset: str, release_version: str = CURRENT_RELEASE
) -> str:
    """
    Fetch path to TSV file containing extra fields to append to VCF header.

    Extra fields are VEP and dbSNP versions.

    :param subset: One of the possible release subsets (e.g., hgdp_1kg)
    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: Filepath for extra fields TSV file
    """
    if release_version == "3.0":
        raise DataException(
            "Extra fields to append to VCF header TSV only exists for releases after v3!"
        )
    return f"gs://gnomad/release/{release_version}/vcf/genomes/extra_fields_for_header{f'_{subset}' if subset else ''}.tsv"


def hgdp_tgp_subset(
    dense: bool = False, test: bool = False
) -> VersionedMatrixTableResource:
    """
    Get the HGDP + 1KG/TGP subset release MatrixTableResource.

    :param dense: If True, return the dense MT; if False, return the sparse MT
    :param test: If true, will return the annotation resource for testing purposes
    :return: MatrixTableResource for specific subset
    """

    return VersionedMatrixTableResource(
        default_version=CURRENT_HGDP_TGP_RELEASE,
        versions={
            release: MatrixTableResource(
                f"{qc_temp_prefix(version=release) if test else f'gs://gnomad/release/{release}/mt/'}/gnomad.genomes.v{release}.hgdp_1kg_subset{f'_dense' if dense else '_sparse'}.mt"
            )
            for release in HGDP_TGP_RELEASES
        },
    )


def hgdp_tgp_subset_annotations(
    sample: bool = True, test: bool = False
) -> VersionedTableResource:
    """
    Get the HGDP + 1KG/TGP subset release sample or variant TableResource.

    :param sample: If true, will return the sample annotations, otherwise will return the variant annotations
    :param test: If true, will return the annotation resource for testing purposes
    :return: Table resource with sample/variant annotations for the subset
    """
    return VersionedTableResource(
        default_version=CURRENT_HGDP_TGP_RELEASE,
        versions={
            release: TableResource(
                f"{qc_temp_prefix(version=release) if test else f'gs://gnomad/release/{release}/ht/'}gnomad.genomes.v{release}.hgdp_1kg_subset{f'_sample_meta' if sample else '_variant_annotations'}.ht"
            )
            for release in HGDP_TGP_RELEASES
            if release != "3"
        },
    )


def hgdp_tgp_subset_sample_tsv(
    release: str = CURRENT_HGDP_TGP_RELEASE, test: bool = False
) -> str:
    """
    Get the path to the HGDP + 1KG/TGP subset release sample annotation text file.

    :param release: Version of annotation tsv path to return
    :param test: If true, will return the sample tsv path for testing purposes
    :return: Path to file
    """
    return f"{qc_temp_prefix(version=release) if test else f'gs://gnomad/release/{release}/tsv/'}gnomad.genomes.v{release}.hgdp_1kg_subset_sample_meta.tsv.bgz"
