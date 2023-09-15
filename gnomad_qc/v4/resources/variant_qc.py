"""Script containing variant QC related resources."""
from typing import Optional, Union

from gnomad.resources.grch38 import (
    na12878_giab,
    na12878_giab_hc_intervals,
    syndip,
    syndip_hc_intervals,
)
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    TableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
)

from gnomad_qc.v4.resources.constants import CURRENT_VERSION, VERSIONS

SYNDIP = "CHMI_CHMI3_Nex1"
"""
String representation for syndip truth sample.
"""

NA12878 = "ASC-4Set-1573S_NA12878@1075619236"
"""
String representation for NA12878 truth sample.
"""

UKB_NA12878 = "Coriell_NA12878_NA12878"
"""
String representation for the UKB Regeneron generated NA12878 truth sample.
"""

TRUTH_SAMPLES = {
    "syndip": {"s": SYNDIP, "truth_mt": syndip, "hc_intervals": syndip_hc_intervals},
    "NA12878": {
        "s": NA12878,
        "truth_mt": na12878_giab,
        "hc_intervals": na12878_giab_hc_intervals,
    },
    "UKB_NA12878": {
        "s": UKB_NA12878,
        "truth_mt": na12878_giab,
        "hc_intervals": na12878_giab_hc_intervals,
    },
}
"""
Dictionary containing necessary information for truth samples

Current truth samples available are syndip and NA12878. Available data for each are the following:
    - s: Sample name in the callset
    - truth_mt: Truth sample MatrixTable resource
    - hc_intervals: High confidence interval Table resource in truth sample
"""


def _variant_qc_root(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Return path to variant QC root folder.

    :param version: Version of variant QC path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :return: Root to variant QC path.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/variant_qc/exomes"
        if test
        else f"gs://gnomad/v{version}/variant_qc/exomes"
    )


def get_callset_truth_data(
    truth_sample: str, mt: bool = True
) -> Union[VersionedMatrixTableResource, VersionedTableResource]:
    """
    Get resources for the truth sample data that is subset from the full callset.

    If `mt` this will return the truth sample MatrixTable (subset from callset); otherwise it returns the
    merged truth sample Table that includes both the truth data and the data from the callset

    :param str truth_sample: Name of the truth sample
    :param bool mt: Whether path is for a MatrixTable, default is True
    :return: Path to callset truth sample MT
    """
    if mt:
        return VersionedMatrixTableResource(
            CURRENT_VERSION,
            {
                version: MatrixTableResource(
                    f"{_variant_qc_root(version)}/truth_samples/gnomad.exomes.v{version}.{truth_sample}.mt"
                )
                for version in VERSIONS
            },
        )
    else:
        return VersionedTableResource(
            CURRENT_VERSION,
            {
                version: TableResource(
                    f"{_variant_qc_root(version)}/truth_samples/gnomad.exomes.v{version}.{truth_sample}.ht"
                )
                for version in VERSIONS
            },
        )


def get_score_bins(model_id: str, aggregated: bool) -> VersionedTableResource:
    """
    Return the path to a Table containing RF or VQSR scores and annotated with a bin based on rank of the metric scores.

    :param model_id: RF or VQSR model ID for which to return score data.
    :param bool aggregated: Whether to get the aggregated data.
         If True, will return the path to Table grouped by bin that contains aggregated variant counts per bin.
    :return: Path to desired hail Table
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version)}/score_bins/gnomad.exomes.v{version}.{model_id}.{'aggregated' if aggregated else 'bins'}.ht"
            )
            for version in VERSIONS
        },
    )


def get_binned_concordance(model_id: str, truth_sample: str) -> VersionedTableResource:
    """
    Return the path to a truth sample concordance Table.

    This Table contains concordance information (TP, FP, FN) between a truth sample within the callset and the
    sample's truth data, grouped by bins of a metric (RF or VQSR scores).

    :param model_id: RF or VQSR model ID for which to return score data.
    :param truth_sample: Which truth sample concordance to analyze (e.g., "NA12878" or "syndip")
    :return: Path to binned truth data concordance Hail Table
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version)}/binned_concordance/gnomad.exomes.v{version}.{truth_sample}.{model_id}.binned_concordance.ht"
            )
            for version in VERSIONS
        },
    )


def get_rf_run_path(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Return the path to the json file containing the RF runs list.

    :param version: Version of RF path to return.
    :param test: Whether to return the test RF runs list.
    :return: Path to json file.
    """
    return (
        f"{_variant_qc_root(version, test=test)}/rf/gnomad.exomes.v{version}.rf_runs.json"
    )


def get_rf_model_path(
    model_id: str, version: str = CURRENT_VERSION, test: bool = False
) -> str:
    """
    Get the path to the RF model for a given run.

    :param model_id: RF run to load.
    :param version: Version of model path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :return: Path to the RF model.
    """
    return (
        f"{_variant_qc_root(version, test=test)}/rf/models/{model_id}/gnomad.exomes.v{version}.rf.model"
    )


def get_rf_training(model_id: str, test: bool = False) -> VersionedTableResource:
    """
    Get the training data for a given run.

    :param model_id: RF run to load.
    :param test: Whether to use a tmp path for variant QC tests.
    :return: VersionedTableResource for RF training data.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test)}/rf/models/{model_id}/gnomad.exomes.v{version}.training.ht"
            )
            for version in VERSIONS
        },
    )


def get_rf_result(
    model_id: Optional[str] = None, test: bool = False
) -> VersionedTableResource:
    """
    Get the results of RF filtering for a given run.

    :param model_id: RF run to load.
    :param test: Whether to use a tmp path for variant QC tests.
    :return: VersionedTableResource for RF filtered data.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test)}/rf/models/{model_id}/gnomad.exomes.v{version}.rf_result.ht"
            )
            for version in VERSIONS
        },
    )


def final_filter(test: bool = False) -> VersionedTableResource:
    """
    Get finalized variant QC filtering Table.

    :param test: Whether to use a tmp path for variant QC tests.
    :return: VersionedTableResource for final variant QC data.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test)}/gnomad.exomes.v{version}.final_filter.ht"
            )
            for version in VERSIONS
        },
    )
