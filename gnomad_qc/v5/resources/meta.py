"""Script containing metadata related resources."""

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from gnomad_qc.v5.resources.constants import (
    CURRENT_PROJECT_META_VERSION,
    WORKSPACE_BUCKET,
)

_project_meta_versions = {
    "5.0": TableResource(
        path=f"{WORKSPACE_BUCKET}/v5.0/metadata/gnomad.v5.0.project_meta.ht"
    )
}

project_meta = VersionedTableResource(
    CURRENT_PROJECT_META_VERSION, _project_meta_versions
)
