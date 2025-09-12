"""Script containing metadata related resources."""

from gnomad.resources.resource_utils import (
    ExpressionResource,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v5.resources.constants import (
    CURRENT_PROJECT_META_VERSION,
    WORKSPACE_BUCKET,
)

_project_meta_versions = {
    "5.0": TableResource(
        path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.project_meta.ht"
    )
}

project_meta = VersionedTableResource(
    CURRENT_PROJECT_META_VERSION, _project_meta_versions
)

sample_id_collisions = TableResource(
    path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.sample_id_collisions.ht"
)

low_quality_samples = ExpressionResource(
    path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.low_quality_samples.he",
)
"""
SetExpression containing IDs of 3 samples with an unspecified data quality issue.

For more information, see Known Issue #1 in the AoU QC document:
https://support.researchallofus.org/hc/en-us/articles/29390274413716-All-of-Us-Genomic-Quality-Report.
"""

failing_metrics_samples = ExpressionResource(
    path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.failing_genomic_metrics_samples.he",
)
"""
SetExpression containing IDs of 4030 samples failing coverage hard filters and 1490 samples with non-XX/XY sex ploidies.

For more information about samples failing coverage hard filters, see
docstring of `get_aou_failing_genomic_metrics_samples`.
"""

samples_to_exclude = ExpressionResource(
    path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.samples_to_exclude.he",
)
"""
SetExpression containing IDs of 5514 samples to exclude from v5 analysis.

Contains samples that should not have been included in the AoU v8 release
(3 samples with unspecified quality issues and 4030 samples failing coverage hard filters)
and 1490 samples with non-XX/XY sex ploidies.

The total number of samples to exclude is 5514, not 5523 because 9 samples both fail coverage filters
and have non-XX/XY sex ploidies.
"""

consent_samples_to_drop = TableResource(
    path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.consent_samples_to_drop.ht",
)
"""
Table containing IDs of 897 samples that are no longer consented to be in gnomAD.

Samples are from the following projects:
- RP-1061: 776 samples.
- RP-1411: 121 samples.
"""
