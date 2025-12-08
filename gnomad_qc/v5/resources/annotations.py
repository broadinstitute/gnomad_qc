"""Script containing annotation related resources."""

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import (
    ANNOTATION_VERSIONS,
    CURRENT_ANNOTATION_VERSION,
    GNOMAD_BUCKET,
    WORKSPACE_BUCKET,
)


def _annotations_root(
    version: str = CURRENT_ANNOTATION_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "aou",
    environment: str = "rwb",
) -> str:
    """
    Get root path to the variant annotation files.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "genomes".
    :param data_set: Data set of annotation resource. Default is "aou".
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb", "batch", or "dataproc".
    :return: Root path of the variant annotation files.
    """
    path_suffix = f"annotations/{data_type}/{data_set}"

    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )

    base_bucket = WORKSPACE_BUCKET if environment == "rwb" else GNOMAD_BUCKET
    return f"gs://{base_bucket}/v{version}/{path_suffix}"


######################################################################
# Variant QC annotation resources
######################################################################


def get_trio_stats(
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Get gnomAD v5 (AoU genomes only) trio stats VersionedTableResource.

    :param test: Whether to use a temporary path for testing.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb", "batch", or "dataproc".
    :return: AoU trio stats VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}."
                "trio_stats.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_sib_stats(
    test: bool = False, environment: str = "rwb"
) -> VersionedTableResource:
    """
    Get the gnomAD v5 (AoU genomes only) sibling stats VersionedTableResource.

    :param test: Whether to use a tmp path for testing.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb", "batch", or "dataproc".
    :return: AoU sibling stats VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.sib_stats.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


######################################################################
# Frequency, coverage, and AN annotation resources
######################################################################


def get_aou_downsampling(test: bool = False) -> VersionedTableResource:
    """
    Get the downsampling annotation table.

    v5 downsamplings only applies to the AoU dataset.

    :param test: Whether to use a tmp path for tests. Default is False.
    :return: Hail Table containing downsampling annotations.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.genomes.v{version}.downsampling.aou.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def group_membership(
    test: bool = False,
    data_set: str = "aou",
) -> VersionedTableResource:
    """
    Get the group membership Table for coverage, AN, quality histograms, and frequency calculations.

    :param test: Whether to use a tmp path for tests. Default is False.
    :param data_set: Data set of annotation resource. Default is "aou".
    :return: Hail Table containing group membership annotations.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, data_set=data_set)}/gnomad.genomes.v{version}.group_membership.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def qual_hists(test: bool = False) -> VersionedTableResource:
    """
    Get the quality histograms annotation table.

    :param test: Whether to use a tmp path for tests. Default is False.
    :return: Hail Table containing quality histogram annotations.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.genomes.v{version}.qual_hists.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def coverage_and_an_path(
    test: bool = False,
    data_set: str = "aou",
) -> VersionedTableResource:
    """
    Fetch filepath for all sites coverage or allele number Table.

    .. note ::

        If `data_set` is 'gnomAD', the returned table only contains coverage and AN for consent drop samples.

    :param test: Whether to use a tmp path for testing. Default is False.
    :param data_set: Dataset identifier. Must be one of "aou" or "gnomad". Default is "aou".
    :return: Coverage and allele number Hail Table.
    """
    assert data_set in ["aou", "gnomad"], "data_set must be either 'aou' or 'gnomad'"

    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, data_set=data_set)}/{'aou' if data_set == 'aou' else 'gnomad'}.genomes.v{version}.coverage_and_an.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )
