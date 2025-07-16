"""Script containing release related resources."""

import logging
from typing import Optional

from gnomad.resources.grch38.gnomad import all_sites_an, coverage, public_release
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v4.resources.release import FREQUENCY_README
from gnomad_qc.v5.resources.constants import (
    ALL_SITES_AN_RELEASES,
    COVERAGE_RELEASES,
    CURRENT_ALL_SITES_AN_RELEASE,
    CURRENT_COVERAGE_RELEASE,
    CURRENT_RELEASE,
    RELEASES,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("release_resources")
logger.setLevel(logging.INFO)


def _release_root(
    version: str = CURRENT_RELEASE,
    test: bool = False,
    data_type: str = "genomes",
    extension: str = "ht",
) -> str:
    """
    Get root path to the release files.

    :param version: Version of release path to return.
    :param test: Whether to use a tmp path for testing.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
        Default is "exomes".
    :param extension: File extension of release file. Default is "ht".
    :return: Root path of the release files.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/release/{extension}/{data_type}"
        if test
        else f"gs://gnomad/release/{version}/{extension}/{data_type}"
    )


def annotation_hists_params_path(
    release_version: str = CURRENT_RELEASE,
    data_type: str = "genomes",
) -> str:
    """
    Return path to file containing dictionary of parameters for site metric histograms.

    The keys of the dictionary are the names of the site quality metrics
    while the values are: [lower bound, upper bound, number  of bins].
    For example, "InbreedingCoeff": [-0.25, 0.25, 50].

    :param release_version: Release version. Defaults to CURRENT RELEASE.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
        Default is "genomes".
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Path to file with annotation histograms
    """
    return f"{_release_root(version=release_version, data_type=data_type, extension='json')}/gnomad.{data_type}.v{release_version}_annotation_hist_params.json"


def qual_hists_json_path(
    release_version: str = CURRENT_RELEASE,
    data_type: str = "genomes",
    test: bool = False,
) -> str:
    """
    Fetch filepath for qual histograms JSON.

    :param release_version: Release version. Defaults to CURRENT RELEASE
    :param data_type: Data type 'exomes' or 'genomes'. Default is 'genomes'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: File path for histogram JSON
    """
    return f"{_release_root(release_version, test, data_type, extension='json')}/gnomad.{data_type}.v{release_version}_qual_hists.json"


def release_ht_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
    test: bool = False,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables.

    :param data_type: Data type of release resource to return. Should be one of
        'exomes', 'genomes' or 'joint'. Default is 'genomes'.
    :param release_version: Release version. Default is CURRENT_RELEASE.
    :param public: Whether release sites Table path returned is from public instead of private
        bucket. Default is True.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: File path for desired release Hail Table.
    """
    if public:
        if file_exists(public_release(data_type).versions[release_version].path):
            return public_release(data_type).versions[release_version].path
        else:
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.ht"
    else:
        return f"{_release_root(version=release_version, test=test, data_type=data_type)}/gnomad.{data_type}.v{release_version}.sites.ht"


def release_sites(
    data_type: str = "genomes", public: bool = False, test: bool = False
) -> VersionedTableResource:
    """
    Retrieve versioned resource for sites-only release Table.

    :param data_type: Data type of release resource to return. Should be one of
        'exomes', 'genomes', or 'joint'. Default is 'genomes'.
    :param public: Whether release sites Table path returned is from public (True) or private (False)
        bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Sites-only release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_ht_path(
                    data_type=data_type,
                    release_version=release,
                    public=public,
                    test=test,
                )
            )
            for release in RELEASES
        },
    )


def release_vcf_path(
    release_version: Optional[str] = None,
    test: bool = False,
    data_type: str = "genomes",
    contig: Optional[str] = None,
) -> str:
    """
    Fetch bucket for release (sites-only) VCFs.

    :param release_version: Release version. When no release_version is supplied
        CURRENT_RELEASE is used.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param data_type: Data type of release resource to return. Should be one of
        'exomes' or 'genomes'. Default is 'genomes'.
    :param contig: String containing the name of the desired reference contig. Default
        is the full (all contigs) sites VCF path.
    :return: Filepath for the desired VCF.
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    if contig:
        return f"{_release_root(version=release_version, test=test, data_type=data_type, extension='vcf')}/gnomad.{data_type}.v{release_version}.sites.{contig}.vcf.bgz"
    else:
        # If contig is None, return path to sharded vcf bucket.
        # NOTE: need to add .bgz or else hail will not bgzip shards.
        return f"{_release_root(version=release_version, test=test, data_type=data_type, extension='vcf')}/gnomad.{data_type}.v{release_version}.sites.vcf.bgz"


def release_header_path(
    release_version: Optional[str] = None,
    data_type: str = "genomes",
    test: bool = False,
) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param release_version: Release version. When no release_version is supplied
        CURRENT_RELEASE is used
    :param data_type: Data type of release resource to return. Should be one of
        'exomes' or 'genomes'. Default is 'genomes'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Filepath for header dictionary pickle.
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    return f"{_release_root(version=release_version, test=test, data_type=data_type, extension='vcf')}/gnomad.{data_type}.v{release_version}_header_dict.pickle"


def append_to_vcf_header_path(
    subset: str = None,
    release_version: str = CURRENT_RELEASE,
    data_type: str = "genomes",
) -> str:
    """
    Fetch path to TSV file containing extra fields to append to VCF header.

    Extra fields are VEP and dbSNP versions.

    :param subset: One of the possible release subsets.
    :param release_version: Release version. Defaults to CURRENT RELEASE.
    :param data_type: Data type of release resource to return. Should be one of
        'exomes' or 'genomes'. Default is 'genomes'.
    :return: Filepath for extra fields TSV file.
    """
    return f"gs://gnomad/release/{release_version}/vcf/{data_type}/extra_fields_for_header{f'_{subset}' if subset else ''}.tsv"


def release_coverage_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = False,
    test: bool = False,
    raw: bool = True,
    coverage_type: str = "coverage",
    data_set: str = "joint",
) -> str:
    """
    Fetch filepath for all sites coverage or allele number release Table.

    .. note ::

        If `data_set` is 'aou' and `raw` is True, v5 genomes coverage returns Table with coverage, all sites AN, and qual hists.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param release_version: Release version.
    :param public: Determines whether release coverage Table is read from public (True) or
        private (False) bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param raw: Whether to return path to raw Table. Default is True. Only applies to Table in private bucket.
    :param coverage_type: 'coverage' or 'allele_number'. Default is 'coverage'.
    :param data_set: Dataset identifier. Must be one of "aou" or "joint". Default is "joint".
    :return: File path for desired coverage Hail Table.
    """
    assert coverage_type in [
        "coverage",
        "allele_number",
    ], "coverage_type must be either 'coverage' or 'allele_number'"
    assert data_set in ["aou", "joint"], "data_set must be either 'aou' or 'joint'"

    if data_set == "aou" and coverage_type == "allele_number":
        raise ValueError(
            "allele_number is not supported for AoU genomes; this information is merged into the coverage file."
        )

    if public:
        if test:
            raise ValueError("Cannot use test=True with public=True!")
        if raw:
            raise ValueError("Cannot use raw=True with public=True!")
        try:
            if coverage_type == "coverage":
                cov = coverage(data_type)
            else:
                cov = all_sites_an(data_type)
            if release_version in cov.versions:
                path = cov.versions[release_version].path
            else:
                path = None
        except DataException:
            path = None
        if path is None:
            raise ValueError(
                f"No public {coverage_type} Table found for data_type {data_type} and release {release_version}."
            )
    else:
        return f"{_release_root(release_version, test=test, data_type=data_type)}/{'aou' if data_set == 'aou' else 'gnomad'}.{data_type}.v{release_version}.{coverage_type}{'.raw' if raw else ''}.ht"


def release_coverage_tsv_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_COVERAGE_RELEASE["genomes"],
    test: bool = False,
) -> str:
    """
    Fetch path to coverage TSV file.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param release_version: Release version. Default is CURRENT_COVERAGE_RELEASE["genomes"].
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Coverage TSV path.
    """
    return f"{_release_root(release_version, test=test, data_type=data_type, extension='tsv')}/gnomad.{data_type}.v{release_version}.coverage.all.tsv.bgz"


def release_all_sites_an_tsv_path(
    data_type: str = "genomes",
    release_version: str = None,
    test: bool = False,
) -> str:
    """
    Fetch path to all sites AN TSV file.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version. Default is
        CURRENT_ALL_SITES_AN_RELEASE[data_type].
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: All sites AN TSV path.
    """
    release_version = (
        release_version
        if release_version is not None
        else CURRENT_ALL_SITES_AN_RELEASE[data_type]
    )
    return f"{_release_root(release_version, test=test, data_type=data_type, extension='tsv')}/gnomad.{data_type}.v{release_version}.allele_number.tsv.bgz"


def release_coverage(
    data_type: str = "genomes",
    public: bool = False,
    test: bool = False,
    raw: bool = True,
    data_set: str = "joint",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for coverage release Table.

    .. note ::

        If `data_set` is 'aou' and `raw` is True, v5 genomes coverage returns Table with coverage, all sites AN, and qual hists.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param public: Determines whether release coverage Table is read from public (True) or
        private (False) bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param raw: Whether to return path to raw Table. Default is True. Only applies to Table in private bucket.
    :param data_set: Dataset identifier. Must be one of "aou" or "joint". Default is "joint".
    :return: Coverage release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_COVERAGE_RELEASE[data_type],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    data_type=data_type,
                    release_version=release,
                    public=public,
                    test=test,
                    raw=raw,
                    data_set=data_set,
                )
            )
            for release in COVERAGE_RELEASES[data_type]
        },
    )


def release_all_sites_an(
    data_type: str = "genomes",
    public: bool = False,
    test: bool = False,
    raw: bool = True,
    data_set: str = "joint",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for all sites allele number release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param public: Determines whether release allele number Table is read from public or
        private bucket. Default is private.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param raw: Whether to return path to raw Table. Default is True. Only applies to Table in private bucket.
    :param data_set: Dataset identifier. Must be one of "aou" or "joint". Default is "joint".
    :return: All sites allele number release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_ALL_SITES_AN_RELEASE[data_type],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    data_type=data_type,
                    release_version=release,
                    public=public,
                    test=test,
                    raw=raw,
                    coverage_type="allele_number",
                    data_set=data_set,
                )
            )
            for release in ALL_SITES_AN_RELEASES[data_type]
        },
    )


def included_datasets_json_path(
    data_type: str = "genomes",
    test: bool = False,
    release_version: str = CURRENT_RELEASE,
) -> str:
    """
    Fetch filepath for the JSON containing all datasets used in the release.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param release_version: Release version. Defaults to CURRENT RELEASE.
    :return: File path for release versions included datasets JSON.
    """
    return f"{_release_root(release_version, test=test, data_type=data_type, extension='json')}/gnomad.{data_type}.v{release_version}.included_datasets.json"


def validated_release_ht(
    test: bool = False,
    data_type: str = "genomes",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for validated sites-only release Table.

    :param test: Whether to use a tmp path for testing. Default is False.
    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :return: Validated release Table
    """
    return VersionedTableResource(
        CURRENT_RELEASE,
        {
            version: TableResource(
                path=(
                    f"{_release_root(version, data_type=data_type, test=test)}/gnomad.{data_type}.v{version}.validated_release.ht"
                )
            )
            for version in RELEASES
        },
    )


def get_freq_array_readme(data_type: str = "genomes") -> str:
    """
    Fetch README for freq array for the specified data_type.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :return: README for freq array.
    """
    # Add information about downsamplings only to exomes
    # (genomes release does not contain downsampling)
    if data_type == "exomes":
        return FREQUENCY_README.format(
            "\ndownsampling_group_gen-anc, e.g."
            " “200_eas_adj”,\ndownsampling_group_gen-anc, e.g. “non_ukb_218035_eas_adj”"
        )
    else:
        return FREQUENCY_README.format("")
