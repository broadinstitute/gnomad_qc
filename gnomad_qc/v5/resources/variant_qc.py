"""Script containing variant QC related resources for v5."""

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from gnomad_qc.v5.resources.constants import (
    CURRENT_VARIANT_QC_RESULT_VERSION,
    CURRENT_VARIANT_QC_VERSION,
    GNOMAD_BUCKET,
    VARIANT_QC_RESULT_VERSIONS,
    VARIANT_QC_VERSIONS,
    WORKSPACE_BUCKET,
)


def _variant_qc_root(
    version: str = CURRENT_VARIANT_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    environment: str = "rwb",
) -> str:
    """
    Return path to variant QC root folder.

    :param version: Version of variant QC path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param data_type: Whether to return 'exomes' or 'genomes' data. Default is genomes.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: Root to variant QC path.
    """
    if test:
        return f"gs://gnomad-tmp/gnomad_v{version}_testing/variant_qc/{data_type}"

    base_bucket = WORKSPACE_BUCKET if environment == "rwb" else GNOMAD_BUCKET
    return f"gs://{base_bucket}/v{version}/variant_qc/{data_type}"


def get_rf_run_path(
    version: str = CURRENT_VARIANT_QC_VERSION,
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Return the path to the json file containing the RF runs list.

    :param version: Version of RF path to return.
    :param test: Whether to return the test RF runs list.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: Path to json file.
    """
    return f"{_variant_qc_root(version, test=test, environment=environment)}/rf/gnomad.genomes.v{version}.rf_runs.json"


def get_rf_model_path(
    model_id: str,
    version: str = CURRENT_VARIANT_QC_VERSION,
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Get the path to the RF model for a given run.

    :param model_id: RF run to load.
    :param version: Version of model path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: Path to the RF model.
    """
    return f"{_variant_qc_root(version, test=test, environment=environment)}/rf/models/{model_id}/gnomad.genomes.v{version}.rf.model"


def get_rf_training(
    model_id: str,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Get the training data for a given run.

    :param model_id: RF run to load.
    :param test: Whether to use a tmp path for variant QC tests.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: VersionedTableResource for RF training data.
    """
    return VersionedTableResource(
        CURRENT_VARIANT_QC_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test, environment=environment)}/rf/models/{model_id}/gnomad.genomes.v{version}.training.ht"
            )
            for version in VARIANT_QC_VERSIONS
        },
    )


def get_variant_qc_result(
    model_id: str,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    r"""
    Get the results of variant QC filtering for a given run.

    :param model_id: Model ID of variant QC run to load. Must start with 'rf\_'.
    :param test: Whether to use a tmp path for variant QC tests.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: VersionedTableResource for variant QC results.
    """
    model_type = model_id.split("_")[0]
    if model_type != "rf":
        raise ValueError(f"Model ID must start with 'rf_', but found {model_id}")
    return VersionedTableResource(
        CURRENT_VARIANT_QC_RESULT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test, environment=environment)}/{model_type}/models/{model_id}/gnomad.genomes.v{version}.{model_type}_result.ht"
            )
            for version in VARIANT_QC_RESULT_VERSIONS
        },
    )


def final_filter(
    data_type: str = "genomes",
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Get finalized variant QC filtering Table.

    :param data_type: Whether to return 'exomes' or 'genomes' data. Default is genomes.
    :param test: Whether to use a tmp path for variant QC tests. Default is False.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb",
        "batch", or "dataproc".
    :return: VersionedTableResource for final variant QC data.
    """
    return VersionedTableResource(
        CURRENT_VARIANT_QC_RESULT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test, data_type=data_type, environment=environment)}/gnomad.{data_type}.v{version}.final_filter.ht"
            )
            for version in VARIANT_QC_RESULT_VERSIONS
        },
    )
