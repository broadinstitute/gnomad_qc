"""Script containing release related resources."""
from typing import Optional

from gnomad.resources.grch38.gnomad import coverage, public_release
from gnomad.resources.resource_utils import TableResource, VersionedTableResource
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v4.resources.constants import CURRENT_RELEASE, RELEASES


def annotation_hists_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Return path to file containing ANNOTATIONS_HISTS dictionary.

    Dictionary contains histogram values for each metric.
    For example, "InbreedingCoeff": [-0.25, 0.25, 50].

    :return: Path to file with annotation histograms
    """
    return f"gs://gnomad/release/{release_version}/json/annotation_hists.json"


def qual_hists_json_path(release_version: str = CURRENT_RELEASE) -> str:
    """
    Fetch filepath for qual histograms JSON.

    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: File path for histogram JSON
    """
    return f"gs://gnomad/release/{release_version}/json/gnomad.exomes.v{release_version}.json"


def release_ht_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables.

    :param data_type: 'exomes' or 'genomes'
    :param release_version: Release version
    :param public: Determines whether release sites Table is read from public or private bucket. Defaults to private
    :return: File path for desired Hail Table
    """
    if public:
        if file_exists(public_release(data_type).versions[release_version].path):
            return public_release(data_type).versions[release_version].path
        else:
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.ht"


def release_sites(public: bool = False) -> VersionedTableResource:
    """
    Retrieve versioned resource for sites-only release Table.

    :param public: Determines whether release sites Table is read from public or private bucket. Defaults to private
    :return: Sites-only release Table
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_ht_path(release_version=release, public=public)
            )
            for release in RELEASES
        },
    )


def release_vcf_path(
    release_version: Optional[str] = None,
    contig: Optional[str] = None,
) -> str:
    """
    Fetch bucket for release (sites-only) VCFs.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used.
    :param contig: String containing the name of the desired reference contig. Defaults to the full (all contigs) sites VCF path
        sites VCF path
    :return: Filepath for the desired VCF
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    if contig:
        return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}.sites.{contig}.vcf.bgz"
    else:
        # If contig is None, return path to sharded vcf bucket.
        # NOTE: need to add .bgz or else hail will not bgzip shards.
        return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}.sites.vcf.bgz"


def release_header_path(release_version: Optional[str] = None) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used
    :return: Filepath for header dictionary pickle
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}_header_dict.pickle"


def append_to_vcf_header_path(
    subset: str, release_version: str = CURRENT_RELEASE
) -> str:
    """
    Fetch path to TSV file containing extra fields to append to VCF header.

    Extra fields are VEP and dbSNP versions.

    :param subset: One of the possible release subsets
    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: Filepath for extra fields TSV file
    """
    return (
        f"gs://gnomad/release/{release_version}/vcf/exomes/extra_fields_for_header{f'_{subset}' if subset else ''}.tsv"
    )


def release_coverage_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
) -> str:
    """
    Fetch filepath for coverage release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version.
    :param public: Determines whether release coverage Table is read from public or
        private bucket. Default is public.
    :return: File path for desired coverage Hail Table.
    """
    if public:
        if file_exists(coverage(data_type).versions[release_version].path):
            return coverage(data_type).versions[release_version].path
        else:
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.coverage.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.coverage.ht"


def release_coverage(public: bool = False) -> VersionedTableResource:
    """
    Retrieve versioned resource for coverage release Table.

    :param public: Determines whether release coverage Table is read from public or
        private bucket. Default is private.
    :return: Coverage release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_coverage_path(release_version=release, public=public)
            )
            for release in RELEASES
        },
    )
