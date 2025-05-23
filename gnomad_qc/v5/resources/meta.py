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

samples_to_exclude = ExpressionResource(
    path=f"gs://{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.samples_to_exclude.he",
)
"""
SetExpression containing IDs of samples to exclude from v5 analysis.

Contains samples that should not have been included in the AoU v8 release
(3 samples with unspecified quality issues and 4030 samples failing coverage hard filters)
and 1490 samples with non-XX/XY sex ploidies.

For more information about samples failing coverage hard filters, see
docstring of `get_aou_failing_genomic_metrics_samples`.
"""
