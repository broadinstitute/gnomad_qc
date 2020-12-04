from typing import Optional, Union

import hail as hl

from gnomad.resources.grch38 import na12878_giab, na12878_giab_hc_intervals, syndip, syndip_hc_intervals
from gnomad.resources.resource_utils import (MatrixTableResource,
                                             TableResource,
                                             VersionedTableResource)
from gnomad_qc.v3.resources.constants import CURRENT_RELEASE, RELEASES

VARIANT_QC_ROOT = 'gs://gnomad/variant_qc/genomes_v3.1'

SYNDIP = "CHMI_CHMI3_WGS2"
"""
String representation for syndip truth sample
"""

NA12878 = "NWD714003"
"""
String representation for NA12878 truth sample
"""

TRUTH_SAMPLES = {
    "syndip": {
        "s": SYNDIP,
        "truth_mt": syndip,
        "hc_intervals": syndip_hc_intervals,
    },
    "NA12878": {
        "s": NA12878,
        "truth_mt": na12878_giab,
        "hc_intervals": na12878_giab_hc_intervals,
    },
}
"""
Dictionary containing necessary information for truth samples (syndip, NA12878)
"""


def get_variant_qc_root(version: str = CURRENT_RELEASE) -> str:
    """
    Return path to variant QC root folder

    :param version: Version of variant QC path to return
    :return: Root to sample QC path
    """
    return f"gs://gnomad/variant_qc/genomes_v{version}"


def get_callset_truth_data(truth_sample: str, mt: bool = True) -> Union[MatrixTableResource, TableResource]:
    """
    Get resources for the truth sample data that is subset from the full callset
    If `mt` this will return the truth sample MatrixTable (subset from callset) otherwise it returns the
    merged truth sample Table that includes both the truth data and the data from the callset
    :param str truth_sample: Name of the truth sample
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Path to callset truth sample MT
    :rtype: str
    """
    if mt:
        return MatrixTableResource(f"{VARIANT_QC_ROOT}/truth_samples/{truth_sample}.mt")
    else:
        return TableResource(f"{VARIANT_QC_ROOT}/truth_samples/{truth_sample}.ht")


def get_truth_sample_data(
        truth_sample: str = None,
        data_type: str = None,
) -> Union[str, hl.Table, hl.MatrixTable]:
    """
    Returns relevant data for truth samples.
    Current truth samples available are syndip and NA12878.
    `data_type` should be one of the following:
    - s: sample name in the callset
    - truth_mt: truth sample MatrixTable
    - hc_intervals: high confidence interval Table in truth sample
    :param str truth_sample: Name of the truth sample. One of the truth samples with information in `truth_sample_dict`
    :param str data_type: Truth sample data type. One of 's', 'truth_mt', or 'hc_intervals'
    :return: Sample name, Table, or MatrixTable of requested truth sample data
    :rtype: Union[str, hl.Table, hl.MatrixTable]
    """

    if not truth_sample or not data_type:
        raise DataException("Must specify both truth_sample and desired data_type")

    if truth_sample not in TRUTH_SAMPLES:
        raise DataException("This truth sample is not present")

    truth_samples_info = TRUTH_SAMPLES[truth_sample]
    if data_type not in truth_samples_info:
        raise DataException(
            f"This data type is not present for truth sample: {truth_sample}"
        )

    return truth_samples_info[data_type]


def get_transmitted_singleton_vcf_path(confidence: str) -> str:
    return f'{VARIANT_QC_ROOT}/transmitted_singletons_{confidence}.vcf.bgz'


def get_score_quantile_bins(model_id: str, aggregated: bool) -> TableResource:
    return TableResource('{}/{}.{}.ht'.format(
        VARIANT_QC_ROOT,
        model_id,
        'binned' if aggregated else 'rank'
    ))


def get_binned_concordance(model_id: str, truth_sample: str) -> TableResource:
    return TableResource(f'{VARIANT_QC_ROOT}/models/{model_id}/{truth_sample}_{model_id}_binned_concordance.ht')


def get_rf_annotated(adj: bool = False) -> TableResource:
    """
    Returns the path to the RF-ready annotated HT
    :param bool adj: Whether to load 'adj' or 'raw'
    :return: Table with RF annotations
    """
    return TableResource(f'{VARIANT_QC_ROOT}/annotated_rf{".adj" if adj else ""}.ht')


def rf_run_hash_path():
    """
    Returns the path to the json file containing the RF runs list.
    :return: Path to json file
    :rtype: str
    """

    return f"{VARIANT_QC_ROOT}/rf_runs.json"


def get_final_rf():
    """
    :return:
    :rtype: str
    """

    return TableResource(f"{VARIANT_QC_ROOT}/rf_final.ht")


def get_rf(
    data: str = "rf_result",
    run_hash: Optional[str] = None,
) -> Union[str, TableResource]:
    """
    Gets the path to the desired RF data.
    Data can take the following values:
        - 'training': path to the training data for a given run
        - 'model': path to pyspark pipeline RF model
        - 'rf_result' (default): path to HT containing result of RF filtering
    :param str data: One of 'training', 'model' or 'rf_result' (default)
    :param str run_hash: Hash of RF run to load
    :return: Path to desired RF data
    """

    if data == "model":
        return f"{VARIANT_QC_ROOT}/models/{run_hash}/{data}.model"
    else:
        return TableResource(f"{VARIANT_QC_ROOT}/models/{run_hash}/{data}.ht")


def get_checkpoint_path(name: str = None, mt: bool = False) -> str:
    """
    Creates a checkpoint path for Table or MatrixTable
    :param str name: Name of intermediate Table/MatrixTable
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    :rtype: str
    """
    return f'gs://gnomad-tmp/{name}.{"mt" if mt else "ht"}'


def get_filtering_model(model_id: str, split: bool = True, finalized: bool = True) -> TableResource:
    """
       Gets the specified filtering annotation resource.
       :param model_id: Filtering model id
       :param split: Split or multi-allelic version of the filtering file
       :return: Filtering annotation file
       """
    path = f"{VARIANT_QC_ROOT}/{model_id}{'.finalized' if finalized else ''}{'.split' if split else ''}.ht"
    return TableResource(path)


final_filter = VersionedTableResource(
    CURRENT_RELEASE,
    {release: TableResource(f"{get_variant_qc_root(release)}/filter_final.ht") for release in RELEASES}
)

class DataException(Exception):
    pass
