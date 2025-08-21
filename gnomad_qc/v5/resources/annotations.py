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
    data_set: str = "rwb",
) -> str:
    """
    Get root path to the variant annotation files.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "genomes".
    :param environment: Environment of annotation resource. Default is "rwb".
    :return: Root path of the variant annotation files.
    """
    path_suffix = f"sample_qc/{data_type}/{data_set}"

    if test:
        environment = "rwb" if data_set == "aou" else "dataproc"
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )

    base_bucket = WORKSPACE_BUCKET if data_set == "aou" else GNOMAD_BUCKET
    return f"gs://{base_bucket}/v{version}/{path_suffix}"


def get_downsampling(test: bool = False, subset: str = "aou") -> VersionedTableResource:
    """
    Get the downsampling annotation table.

    .. note::
        v5 downsamplings only applies to the AoU dataset.

    :param test: Whether to use a tmp path for tests. Default is False.
    :param subset: Subset to return downsampling Table for. Default is "aou".
    :return: Hail Table containing subset or overall dataset downsampling annotations.
    """
    if subset != "aou":
        raise ValueError("v5 downsamplings only applies to the AoU dataset.")

    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.genomes.v{version}.downsampling{f'.{subset}' if subset else ''}.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )
