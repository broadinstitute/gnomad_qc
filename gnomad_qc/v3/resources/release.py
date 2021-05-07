from typing import Optional

from gnomad.resources.resource_utils import (
    DataException,
    MatrixTableResource,
    VersionedMatrixTableResource,
)
from gnomad_qc.v3.resources.constants import (
    CURRENT_RELEASE,
    CURRENT_HGDP_TGP_RELEASE,
    HGDP_TGP_RELEASES,
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
    """Fetch filepath for qual histograms JSON

    :param release_version: Release version. Defaults to CURRENT RELEASE
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


def release_var_hist_path(data_source: str, freeze: int) -> str:
    """
    Fetch bucket for release variant histograms (json files).
    
    :param str data_source: One of 'regeneron' or 'broad'
    :param int freeze: One of the data freezes
    :return: Filepath for release jsons
    :rtype: str
    """
    return f"{get_release_path(data_source, freeze)}/json/{data_source}.freeze_{freeze}.json"


def release_header_path(
    subset: Optional[str] = None, release_version: str = CURRENT_RELEASE
) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param subset: Name of the subset (eg: hgdp_tgp).
    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: Filepath for header dictionary pickle
    """
    if subset:
        subset = f"_{subset}"

    return f"gs://gnomad/release/{release_version}/vcf/genomes/gnomad.genomes.v{release_version}_header_dict{subset}.pickle"


def release_vcf_path(
    release_version: str = CURRENT_RELEASE, contig: str = None, subset: str = None
) -> str:
    """
    Fetch bucket for release (sites-only) VCFs.

    :param release_version: Release version. Defaults to CURRENT RELEASE
    :param contig: String containing the name of the desired reference contig
    :param subset: Subset being written out to VCF
    :return: Filepath for the desired VCF
    """
    if subset is None:
        subset = "sites"

    if contig:
        return f"gs://gnomad/release/{release_version}/vcf/genomes/gnomad.genomes.v{release_version}.{subset}.{contig}.vcf.bgz"
    else:
        # if contig is None, return path to sharded vcf bucket
        # NOTE: need to add .bgz or else hail will not bgzip shards
        return f"gs://gnomad/release/{release_version}/vcf/genomes/gnomad.genomes.v{release_version}.sites.vcf.bgz"


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
    if release_version not in {"3.1", "3.1.1"}:
        raise DataException(
            "Extra fields to append to VCF header TSV only exists for 3.1 and 3.1.1!"
        )
    return f"gs://gnomad/release/{release_version}/vcf/genomes/extra_fields_for_header{f'_{subset}' if subset else ''}.tsv"


def release_subset(
    subset: str, dense: bool = False, data_type: str = "genomes",
) -> VersionedMatrixTableResource:
    """
    Get the subset release MatrixTableResource.

    :param subset: One of the possible release subsets (e.g., hgdp_1kg)
    :param dense: If True, will return the dense MatrixTableResource, otherwise will return the sparse MatrixTableResource
    :param data_type: 'exomes' or 'genomes'
    :return: MatrixTableResource for specific subset
    """

    return VersionedMatrixTableResource(
        CURRENT_HGDP_TGP_RELEASE,
        {
            release: MatrixTableResource(
                f"gs://gnomad/release/{release}/mt/genomes/gnomad.{data_type}.v{release}.{subset}_subset{f'_dense' if dense else '_sparse'}.mt"
            )
            for release in HGDP_TGP_RELEASES
            if release != "3"
        },
    )
