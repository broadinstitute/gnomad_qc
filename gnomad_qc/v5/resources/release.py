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
from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import (
    ALL_SITES_AN_RELEASES,
    COVERAGE_RELEASES,
    CURRENT_ALL_SITES_AN_RELEASE,
    CURRENT_COVERAGE_RELEASE,
    CURRENT_RELEASE,
    GNOMAD_BUCKET,
    RELEASES,
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
    data_type: str = "genomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = False,
    test: bool = False,
    raw: bool = True,
    coverage_type: str = "coverage",
    data_set: str = "joint",
    environment: str = "rwb",
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
    :param environment: Environment to use. Default is "rwb".
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
        return path
    else:
        return f"{_release_root(release_version, test=test, data_type=data_type, environment=environment)}/{'aou' if data_set == 'aou' else 'gnomad'}.{data_type}.v{release_version}.{coverage_type}{'.raw' if raw else ''}.ht"


def release_coverage_tsv_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_COVERAGE_RELEASE["genomes"],
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch path to coverage TSV file.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param release_version: Release version. Default is CURRENT_COVERAGE_RELEASE["genomes"].
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Coverage TSV path.
    """
    return f"{_release_root(release_version, test=test, data_type=data_type, extension='tsv', environment=environment)}/gnomad.{data_type}.v{release_version}.coverage.all.tsv.bgz"


def release_all_sites_an_tsv_path(
    data_type: str = "genomes",
    release_version: str = None,
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch path to all sites AN TSV file.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version. Default is
        CURRENT_ALL_SITES_AN_RELEASE[data_type].
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: All sites AN TSV path.
    """
    release_version = (
        release_version
        if release_version is not None
        else CURRENT_ALL_SITES_AN_RELEASE[data_type]
    )
    return f"{_release_root(release_version, test=test, data_type=data_type, extension='tsv', environment=environment)}/gnomad.{data_type}.v{release_version}.allele_number.tsv.bgz"


def release_coverage(
    data_type: str = "genomes",
    public: bool = False,
    test: bool = False,
    raw: bool = True,
    data_set: str = "joint",
    environment: str = "rwb",
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
    :param environment: Environment to use. Default is "rwb".
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
                    environment=environment,
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
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for all sites allele number release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'genomes'.
    :param public: Determines whether release allele number Table is read from public or
        private bucket. Default is private.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param raw: Whether to return path to raw Table. Default is True. Only applies to Table in private bucket.
    :param data_set: Dataset identifier. Must be one of "aou" or "joint". Default is "joint".
    :param environment: Environment to use. Default is "rwb".
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
                    environment=environment,
                )
            )
            for release in ALL_SITES_AN_RELEASES[data_type]
        },
    )
