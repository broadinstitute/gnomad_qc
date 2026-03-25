"""Script containing metadata related resources."""

from gnomad.resources.resource_utils import (
    ExpressionResource,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v5.resources.basics import (
    _SAMPLE_DATA_ENVIRONMENTS,
    _get_base_bucket,
    _validate_environment,
)
from gnomad_qc.v5.resources.constants import (
    CURRENT_PROJECT_META_VERSION,
    CURRENT_SAMPLE_QC_VERSION,
)


def get_project_meta(environment: str = "batch") -> VersionedTableResource:
    """
    Get the VersionedTableResource for per-sample project-level metadata.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: VersionedTableResource for project metadata.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return VersionedTableResource(
        CURRENT_PROJECT_META_VERSION,
        {
            "5.0": TableResource(
                path=f"gs://{bucket}/v5.0/metadata/gnomad.v5.0.project_meta.ht"
            )
        },
    )


def get_sample_id_collisions(environment: str = "batch") -> TableResource:
    """
    Get the TableResource for sample IDs that collide between AoU and gnomAD v4.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: TableResource of sample ID collisions.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return TableResource(
        path=f"gs://{bucket}/v5.0/metadata/gnomad.v5.0.sample_id_collisions.ht"
    )


def get_low_quality_samples(environment: str = "batch") -> ExpressionResource:
    """
    Get the ExpressionResource for AoU-flagged low-quality sample IDs.

    SetExpression containing IDs of 3 samples with an unspecified data quality issue.

    For more information, see Known Issue #1 in the AoU QC document:
    https://support.researchallofus.org/hc/en-us/articles/29390274413716-All-of-Us-Genomic-Quality-Report.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: ExpressionResource of low-quality sample IDs.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return ExpressionResource(
        path=f"gs://{bucket}/v5.0/metadata/gnomad.v5.0.low_quality_samples.he",
    )


def get_failing_metrics_samples(environment: str = "batch") -> ExpressionResource:
    """
    Get the ExpressionResource for samples failing AoU genomic QC metrics.

    SetExpression containing IDs of 4030 samples failing coverage hard filters and
    1490 samples with non-XX/XY sex ploidies.

    For more information about samples failing coverage hard filters, see
    docstring of `get_aou_failing_genomic_metrics_samples`.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: ExpressionResource of failing-metrics sample IDs.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return ExpressionResource(
        path=f"gs://{bucket}/v5.0/metadata/gnomad.v5.0.failing_genomic_metrics_samples.he",
    )


def get_samples_to_exclude_resource(environment: str = "batch") -> ExpressionResource:
    """
    Get the ExpressionResource for the combined set of samples to exclude.

    SetExpression containing IDs of 5514 samples to exclude from v5 analysis.

    Contains samples that should not have been included in the AoU v8 release
    (3 samples with unspecified quality issues and 4030 samples failing coverage hard
    filters) and 1490 samples with non-XX/XY sex ploidies.

    The total number of samples to exclude is 5514, not 5523 because 9 samples both
    fail coverage filters and have non-XX/XY sex ploidies.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: ExpressionResource of sample IDs to exclude.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return ExpressionResource(
        path=f"gs://{bucket}/v5.0/metadata/gnomad.v5.0.samples_to_exclude.he",
    )


def get_consent_samples_to_drop(environment: str = "batch") -> TableResource:
    """
    Get the TableResource for consent-withdrawn samples.

    Table containing IDs of 897 samples that are no longer consented to be in gnomAD.

    Samples are from the following projects:
    - RP-1061: 776 samples.
    - RP-1411: 121 samples.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: TableResource of consent-withdrawn sample IDs.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return TableResource(
        path=f"gs://{bucket}/v5.0/metadata/gnomad.v5.0.consent_samples_to_drop.ht",
    )


# Backward-compatible aliases — resolve rwb environment.
project_meta = get_project_meta(environment="rwb")
sample_id_collisions = get_sample_id_collisions(environment="rwb")
low_quality_samples = get_low_quality_samples(environment="rwb")
failing_metrics_samples = get_failing_metrics_samples(environment="rwb")
samples_to_exclude = get_samples_to_exclude_resource(environment="rwb")
consent_samples_to_drop = get_consent_samples_to_drop(environment="rwb")


def _meta_root_path(
    version: str = CURRENT_PROJECT_META_VERSION,
    environment: str = "batch",
) -> str:
    """
    Retrieve the path to the root metadata directory.

    :param version: gnomAD version.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: String representation of the path to the root metadata directory.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    bucket = _get_base_bucket(
        environment, read_only=True if environment == "batch" else False
    )
    return f"gs://{bucket}/v{version}/metadata/genomes"


def meta(
    version: str = CURRENT_SAMPLE_QC_VERSION,
    data_type: str = "genomes",
    environment: str = "batch",
) -> VersionedTableResource:
    """
    Get the v5 sample QC meta VersionedTableResource.

    .. note::

        Exome data is not currently supported in this function.
        The v4 sample QC meta uses a different structure, so this function
        does not pull or duplicate that data. If exome data are needed, please
        use the v4 resource directly.

    :param version: Sample QC version.
    :param data_type: Data type. Default is "genomes". If "exomes" is supplied, a
        warning will be raised suggesting the use of v4 sample QC metadata.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Sample QC meta VersionedTableResource.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    if data_type == "exomes":
        raise ValueError(
            "Exome sample QC metadata is not supported in v5. "
            "The v4 sample QC meta has a different structure and should be "
            "imported directly from the v4 resource using 'from gnomad_qc.v4.resources.meta import meta as v4_meta'."
        )

    if data_type != "genomes":
        raise ValueError(
            f"Unsupported data_type: {data_type}. Only 'genomes' is supported."
        )

    return VersionedTableResource(
        default_version=CURRENT_SAMPLE_QC_VERSION,
        versions={
            CURRENT_SAMPLE_QC_VERSION: TableResource(
                path=f"{_meta_root_path(version, environment)}/gnomad.genomes.v{version}.sample_qc_metadata.ht"
            )
        },
    )
