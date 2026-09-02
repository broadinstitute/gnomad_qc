"""Script containing metadata related resources."""

import json
import logging
from typing import Optional

import hail as hl
import hailtop.fs as hfs
from gnomad.resources.resource_utils import (
    ExpressionResource,
    TableResource,
    VersionedTableResource,
)
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v5.resources.basics import (
    _SAMPLE_DATA_ENVIRONMENTS,
    _check_resource_existence,
    _get_base_bucket,
    _validate_environment,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.constants import (
    CURRENT_PROJECT_META_VERSION,
    CURRENT_SAMPLE_QC_VERSION,
)

logger = logging.getLogger("meta_resources")
logger.setLevel(logging.INFO)


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


def aou_sample_artifact_path(
    name: str,
    test: bool = False,
    environment: str = "batch",
    version: str = CURRENT_PROJECT_META_VERSION,
) -> str:
    """
    Return the permanent path to a run-invariant AoU sample-level artifact.

    Artifacts are small precomputed files read at VDS load (or by pipeline
    workers) instead of rescanning the ~330-partition sample metadata tables:
    the VDS-loading JSONs written by `write_aou_vds_sample_jsons`
    (``sample_id_collisions.json``, ``high_quality_samples.json``,
    ``release_samples.json``).

    Full-scope artifacts are permanent; unlike `_meta_root_path`, their root
    resolves to the writable bucket in the "batch" environment (artifacts are
    both written and read from batch). Test-scope artifacts are disposable and
    live in the 4-day temp bucket, per repo convention.

    :param name: Artifact file name (with extension).
    :param test: If True, return the test-scoped path (10-sample test VDS runs),
        in the 4-day temp bucket.
    :param environment: Environment to use. Default is "batch". Must be one of
        "rwb" or "batch".
    :param version: gnomAD version.
    :return: GCS path to the artifact.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    if test:
        # Test artifacts are disposable: 4-day temp bucket, per repo convention.
        return (
            f"{qc_temp_prefix(version=version, environment=environment, days=4)}"
            f"aou_sample_artifacts_test/{name}"
        )
    return (
        f"gs://{_get_base_bucket(environment)}/v{version}/metadata/genomes/"
        f"aou_sample_artifacts_full/{name}"
    )


def load_aou_sample_artifact_json(
    name: str,
    test: bool = False,
    environment: str = "batch",
) -> Optional[list]:
    """
    Return the parsed JSON sample artifact, or None when it has not been written.

    :param name: Artifact file name (with extension).
    :param test: If True, read the test-scoped path.
    :param environment: Environment to use. Default is "batch". Must be one of
        "rwb" or "batch".
    :return: Parsed JSON value, or None if the file is absent.
    """
    path = aou_sample_artifact_path(name, test=test, environment=environment)
    if not file_exists(path):
        return None
    with hfs.open(path) as f:
        return json.load(f)


def write_aou_vds_sample_jsons(
    environment: str = "batch",
    test: bool = False,
    overwrite: bool = False,
) -> None:
    """
    Write the permanent sample artifact JSONs used to load the AoU VDS.

    Writes three artifacts so that VDS loads read small files instead of each
    rescanning the ~330-partition sample metadata tables:

    - ``sample_id_collisions.json``: pre-prefix sample IDs that collide with
      gnomAD samples (from the sample-collisions Table).
    - ``high_quality_samples.json``: post-prefix high-quality sample IDs,
      collected with the exact ``meta_ht.high_quality`` filter that
      ``get_aou_vds(high_quality_only=True)`` applies, so the two cannot drift.
    - ``release_samples.json``: post-prefix release sample IDs, collected with
      the exact ``meta_ht.release`` filter that ``get_aou_vds(release_only=True)``
      applies.

    `get_aou_vds` prefers these files and falls back to the original scans for
    any that are absent.

    :param environment: Environment to use. Default is "batch". Must be one of
        "rwb" or "batch".
    :param test: If True, write the test-scoped paths (4-day temp bucket).
    :param overwrite: Whether to replace existing artifacts; without it, a rerun
        fails instead of silently rewriting them. Default is False.
    """
    paths = {
        name: aou_sample_artifact_path(name, test=test, environment=environment)
        for name in (
            "sample_id_collisions.json",
            "high_quality_samples.json",
            "release_samples.json",
        )
    }
    _check_resource_existence(
        environment=environment,
        output_step_resources={
            "--write-vds-sample-artifact-jsons": list(paths.values())
        },
        overwrite=overwrite,
    )

    sc_ht = get_sample_id_collisions(environment=environment).ht()
    collisions = sorted(sc_ht.aggregate(hl.agg.collect_as_set(sc_ht.s)))
    with hfs.open(paths["sample_id_collisions.json"], "w") as f:
        json.dump(collisions, f)
    logger.info(
        "Wrote %d collision sample IDs: %s",
        len(collisions),
        paths["sample_id_collisions.json"],
    )

    meta_ht = meta(data_type="genomes", environment=environment).ht()
    hq_samples = meta_ht.filter(meta_ht.high_quality).s.collect()
    with hfs.open(paths["high_quality_samples.json"], "w") as f:
        json.dump(hq_samples, f)
    logger.info(
        "Wrote %d high-quality sample IDs: %s",
        len(hq_samples),
        paths["high_quality_samples.json"],
    )

    release_samples = meta_ht.filter(meta_ht.release).s.collect()
    with hfs.open(paths["release_samples.json"], "w") as f:
        json.dump(release_samples, f)
    logger.info(
        "Wrote %d release sample IDs: %s",
        len(release_samples),
        paths["release_samples.json"],
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
