"""Script containing annotation related resources."""

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from gnomad_qc.v5.resources.basics import (
    _ALL_ENVIRONMENTS,
    _SAMPLE_DATA_ENVIRONMENTS,
    _get_base_bucket,
    _validate_environment,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.constants import (
    ANNOTATION_VERSIONS,
    CURRENT_ANNOTATION_VERSION,
    CURRENT_VEP_ANNOTATION_VERSION,
    DEFAULT_VEP_VERSION,
    VEP_ANNOTATION_VERSIONS,
)

SAMPLE_ANNOTATION_DEFAULT_ENVIRONMENT = "rwb"
VARIANT_ANNOTATION_DEFAULT_ENVIRONMENT = "batch"


def _annotations_root(
    version: str = CURRENT_ANNOTATION_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "aou",
    environment: str = "batch",
    read_only: bool = False,
) -> str:
    """
    Get root path to the variant annotation files.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "genomes".
    :param data_set: Data set of annotation resource. Default is "aou".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb", "batch", or "dataproc".
    :param read_only: Whether the path is for read-only resources. When True
        and environment is "batch", the read-only bucket is used.
        Default is False.
    :return: Root path of the variant annotation files.
    """
    _validate_environment(environment, _ALL_ENVIRONMENTS)
    path_suffix = f"annotations/{data_type}/{data_set}"

    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )

    return f"gs://{_get_base_bucket(environment, read_only=read_only)}/v{version}/{path_suffix}"


######################################################################
# Variant QC annotation resources
######################################################################


def get_trio_stats(
    test: bool = False,
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Get gnomAD v5 (AoU genomes only) trio stats VersionedTableResource.

    :param test: Whether to use a temporary path for testing.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: AoU trio stats VersionedTableResource.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/aou.genomes.v{version}."
                "trio_stats.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_sib_stats(
    test: bool = False, environment: str = "batch"
) -> VersionedTableResource:
    """
    Get the gnomAD v5 (AoU genomes only) sibling stats VersionedTableResource.

    :param test: Whether to use a tmp path for testing.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: AoU sibling stats VersionedTableResource.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/aou.genomes.v{version}.sib_stats.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


######################################################################
# Frequency, coverage, and AN annotation resources
######################################################################


def get_aou_downsampling(
    test: bool = False, environment: str = "batch"
) -> VersionedTableResource:
    """
    Get the downsampling annotation table.

    v5 downsamplings only applies to the AoU dataset.

    :param test: Whether to use a tmp path for tests. Default is False.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Hail Table containing downsampling annotations.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.downsampling.aou.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def group_membership(
    test: bool = False,
    data_set: str = "aou",
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Get the group membership Table for coverage, AN, quality histograms, and frequency calculations.

    :param test: Whether to use a tmp path for tests. Default is False.
    :param data_set: Data set of annotation resource. Default is "aou".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Hail Table containing group membership annotations.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, data_set=data_set, environment=environment)}/gnomad.genomes.v{version}.group_membership.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def qual_hists(
    test: bool = False, environment: str = "batch"
) -> VersionedTableResource:
    """
    Get the quality histograms annotation table.

    :param test: Whether to use a tmp path for tests. Default is False.
    :param environment: Environment to use for quality histograms. Default is "batch".
        Must be one of "rwb" or "batch".
    :return: Hail Table containing quality histogram annotations.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.qual_hists.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def coverage_and_an_path(
    test: bool = False,
    data_set: str = "aou",
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Fetch filepath for all sites coverage or allele number Table.

    .. note ::

        If `data_set` is 'gnomAD', the returned table only contains coverage and AN for consent drop samples.

    :param test: Whether to use a tmp path for testing. Default is False.
    :param data_set: Dataset identifier. Must be one of "aou" or "gnomad". Default is "aou".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb",
        "batch", or "dataproc".
    :return: Coverage and allele number Hail Table.
    """
    _validate_environment(environment, _ALL_ENVIRONMENTS)
    assert data_set in ["aou", "gnomad"], "data_set must be either 'aou' or 'gnomad'"

    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, data_set=data_set, environment=environment)}/{'aou' if data_set == 'aou' else 'gnomad'}.genomes.v{version}.coverage_and_an.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_freq(
    version: str = CURRENT_ANNOTATION_VERSION,
    data_type: str = "genomes",
    test: bool = False,
    data_set: str = "aou",
    environment: str = "batch",
) -> TableResource:
    """
    Get the frequency annotation Table for v5.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource ("genomes" or "exomes").
    :param test: Whether to use a tmp path for testing.
    :param data_set: Data set of annotation resource. Default is "aou".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb", "batch", or "dataproc".
    :return: Hail Table containing frequency annotations.
    """
    _validate_environment(environment, _ALL_ENVIRONMENTS)
    assert data_set in [
        "aou",
        "gnomad",
        "merged",
    ], "data_set must be either 'aou', 'gnomad', or 'merged'"
    return TableResource(
        f"{_annotations_root(version, test, data_type, data_set, environment)}/{data_set}.genomes.v{version}.frequencies.ht"
    )


######################################################################
# Variant QC annotation resources
######################################################################


def get_info_ht(
    test: bool = False, environment: str = "batch"
) -> VersionedTableResource:
    """
    Get the gnomAD v5 (AoU genomes only) info VersionedTableResource.

    :param test: Whether to use a tmp path for testing.
    :param environment: Environment to use. Default is "batch". Must be one of
        "rwb" or "batch".
    :return: Info VersionedTableResource.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.info.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def info_vcf_path(
    version: str = CURRENT_ANNOTATION_VERSION,
    test: bool = False,
    environment: str = "batch",
) -> str:
    """
    Path to sites VCF (input information for running VQSR).

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for testing.
    :param environment: Environment to use. Must be one of "rwb" or "batch". Default is "batch".
    :return: String for the path to the info VCF.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.info.vcf.bgz"


def get_aou_vcf_header(environment: str = "batch") -> str:
    """
    Get path to AoU annotation sites-only VCF header.

    This is needed for proper import of the sites-only VCF as the QUALapprox
    annotation is stated in the previous header as an int but is actually a float.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Path to the VCF header file.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    return (
        f"{_annotations_root(version='5.0', environment=environment, read_only=True)}"
        "/aou_annotation_sites_only_header.vcf"
    )


def get_aou_annotated_sites_only_vcf(environment: str = "batch") -> str:
    """
    Get path to AoU sites-only VCF with annotations needed for variant QC.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Path to the annotated sites-only VCF.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(environment, read_only=True)
    if environment == "batch":
        return f"gs://{bucket}/aou_sites_vcf/v8/echo_full_gnomad_annotated.sites-only.vcf.gz"
    return f"gs://{bucket}/echo_full_gnomad_annotated.sites-only.vcf.gz"


######################################################################
# VEP resources
######################################################################


def get_vep(
    test: bool = False,
    vep_version: str = DEFAULT_VEP_VERSION,
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Get the gnomAD v5 VEP annotation VersionedTableResource.

    :param test: Whether to use a tmp path for analysis of the test Table instead of
        the full v5 Table.
    :param vep_version: VEP version to use (e.g., "105", "115"). Default is "105".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb", "batch", or "dataproc".
    :return: gnomAD v5 VEP VersionedTableResource.
    """
    _validate_environment(environment, _ALL_ENVIRONMENTS)
    vep_version_postfix = "" if vep_version == DEFAULT_VEP_VERSION else vep_version
    return VersionedTableResource(
        CURRENT_VEP_ANNOTATION_VERSION[vep_version],
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.vep{vep_version_postfix}.ht"
                )
            )
            for version in VEP_ANNOTATION_VERSIONS[vep_version]
        },
    )


def validate_vep_path(
    test: bool = False,
    vep_version: str = DEFAULT_VEP_VERSION,
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Get the gnomAD v5 VEP annotation VersionedTableResource for validation counts.

    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v5 VDS.
    :param vep_version: VEP version to use (e.g., "105", "115"). Default is "105".
    :param environment: Environment to use. Default is "batch". Must be one of "rwb", "batch", or "dataproc".
    :return: gnomAD v5 VEP VersionedTableResource containing validity check.
    """
    _validate_environment(environment, _ALL_ENVIRONMENTS)
    vep_version_postfix = "" if vep_version == DEFAULT_VEP_VERSION else vep_version
    return VersionedTableResource(
        CURRENT_VEP_ANNOTATION_VERSION[vep_version],
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test=test, environment=environment)}/gnomad.genomes.v{version}.vep{vep_version_postfix}.validate.ht"
                )
            )
            for version in VEP_ANNOTATION_VERSIONS[vep_version]
        },
    )
