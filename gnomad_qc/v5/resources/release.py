"""Script containing release related resources."""

import logging

from gnomad.resources.grch38.gnomad import all_sites_an, coverage
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import (
    ALL_SITES_AN_RELEASES,
    COVERAGE_RELEASES,
    CURRENT_ALL_SITES_AN_RELEASE,
    CURRENT_COVERAGE_RELEASE,
    CURRENT_RELEASE,
    GNOMAD_BUCKET,
    WORKSPACE_BUCKET,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("release_resources")
logger.setLevel(logging.INFO)


def _release_root(
    version: str = CURRENT_RELEASE,
    test: bool = False,
    data_type: str = "genomes",
    extension: str = "ht",
    environment: str = "rwb",
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
    path_suffix = f"release/{extension}/{data_type}"
    base_bucket = WORKSPACE_BUCKET if environment == "rwb" else GNOMAD_BUCKET
    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )
    return f"gs://{base_bucket}/v{version}/{path_suffix}"


def release_coverage_path(
    release_version: str = CURRENT_RELEASE,
    public: bool = False,
    test: bool = False,
    coverage_type: str = "coverage",
    environment: str = "rwb",
) -> str:
    """
    Fetch filepath for v5 (AoU + gnomAD v4 genomes) all sites coverage or allele number release Table.

    :param release_version: Release version.
    :param public: Determines whether release coverage Table is read from public (True) or
        private (False) bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param coverage_type: 'coverage' or 'allele_number'. Default is 'coverage'.
    :param environment: Environment to use. Default is "rwb".
    :return: File path for desired coverage Hail Table.
    """
    assert coverage_type in [
        "coverage",
        "allele_number",
    ], "coverage_type must be either 'coverage' or 'allele_number'"

    if public:
        if test:
            raise ValueError("Cannot use test=True with public=True!")
        try:
            if coverage_type == "coverage":
                cov = coverage("genomes")
            else:
                cov = all_sites_an("genomes")
            if release_version in cov.versions:
                path = cov.versions[release_version].path
            else:
                path = None
        except DataException:
            path = None
        if path is None:
            raise ValueError(
                f"No public {coverage_type} Table found for genomes and release {release_version}."
            )
        return path
    else:
        return f"{_release_root(release_version, test=test, environment=environment)}/gnomad.genomes.v{release_version}.{coverage_type}.ht"


def release_coverage_tsv_path(
    release_version: str = CURRENT_COVERAGE_RELEASE["genomes"],
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch path to coverage TSV file.

    :param release_version: Release version. Default is CURRENT_COVERAGE_RELEASE["genomes"].
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Coverage TSV path.
    """
    return f"{_release_root(release_version, test=test, extension='tsv', environment=environment)}/gnomad.genomes.v{release_version}.coverage.all.tsv.bgz"


def release_all_sites_an_tsv_path(
    release_version: str = None,
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch path to all sites AN TSV file.

    :param release_version: Release version. Default is None.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: All sites AN TSV path.
    """
    release_version = (
        release_version
        if release_version is not None
        else CURRENT_ALL_SITES_AN_RELEASE["genomes"]
    )
    return f"{_release_root(release_version, test=test, extension='tsv', environment=environment)}/gnomad.genomes.v{release_version}.allele_number.tsv.bgz"


def release_coverage(
    public: bool = False,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for coverage release Table.

    :param public: Determines whether release coverage Table is read from public (True) or
        private (False) bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Coverage release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_COVERAGE_RELEASE["genomes"],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    release_version=release,
                    public=public,
                    test=test,
                    environment=environment,
                )
            )
            for release in COVERAGE_RELEASES["genomes"]
        },
    )


def release_all_sites_an(
    public: bool = False,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for all sites allele number release Table.

    :param public: Determines whether release allele number Table is read from public or
        private bucket. Default is private.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: All sites allele number release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_ALL_SITES_AN_RELEASE["genomes"],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    release_version=release,
                    public=public,
                    test=test,
                    coverage_type="allele_number",
                    environment=environment,
                )
            )
            for release in ALL_SITES_AN_RELEASES["genomes"]
        },
    )
