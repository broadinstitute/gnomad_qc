"""Script containing variant QC related resources."""

from typing import Union

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

from gnomad_qc.v5.resources.constants import (
    CURRENT_VARIANT_QC_RESULT_VERSION,
    CURRENT_VARIANT_QC_VERSION,
    VARIANT_QC_RESULT_VERSIONS,
    VARIANT_QC_VERSIONS,
)

VQSR_FEATURES = {
    "genomes": {
        "snv": [
            "AS_QD",
            "AS_MQRankSum",
            "AS_ReadPosRankSum",
            "AS_FS",
            "AS_SOR",
            "AS_MQ",
        ],
        "indel": [
            "AS_QD",
            "AS_MQRankSum",
            "AS_ReadPosRankSum",
            "AS_FS",
            "AS_SOR",
        ],
    },
}
"""List of features used in the VQSR model."""

SYNDIP = "CHMI_CHMI3_Nex1"
"""
String representation for syndip truth sample.
"""

NA12878 = "ASC-4Set-1573S_NA12878@1075619236"
"""
String representation for NA12878 truth sample.
"""
# TODO: Is this present ? Likely not but check.
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

Current truth samples available are syndip and NA12878. Available data for each are the
following:

    - s: Sample name in the callset
    - truth_mt: Truth sample MatrixTable resource
    - hc_intervals: High confidence interval Table resource in truth sample
"""


def _variant_qc_root(
    version: str = CURRENT_VARIANT_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
) -> str:
    """
    Return path to variant QC root folder.

    :param version: Version of variant QC path to return.
    :param test: Whether to use a tmp path for variant QC tests.
    :param data_type: Whether to return 'exomes' or 'genomes' data. Default is exomes.
    :return: Root to variant QC path.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/variant_qc/{data_type}"
        if test
        else f"gs://gnomad/v{version}/variant_qc/{data_type}"
    )


def get_variant_qc_result(
    model_id: str, test: bool = False, split: bool = True
) -> VersionedTableResource:
    r"""
    Get the results of variant QC filtering for a given run.

    :param model_id: Model ID of variant QC run to load. Must start with 'rf\_',
        'vqsr\_', or 'if\_'.
    :param test: Whether to use a tmp path for variant QC tests.
    :param split: Whether to return the split or unsplit variant QC result.
    :return: VersionedTableResource for variant QC results.
    """
    model_type = model_id.split("_")[0]
    if model_type not in ["rf", "vqsr", "if"]:
        raise ValueError(
            f"Model ID must start with 'rf_', 'vqsr_', or 'if_', but found {model_id}"
        )
    return VersionedTableResource(
        CURRENT_VARIANT_QC_RESULT_VERSION,
        {
            version: TableResource(
                f"{_variant_qc_root(version, test=test)}/{model_type}/models/{model_id}/gnomad.genomes.v{version}.{model_type}_result{'' if split else '.unsplit'}.ht"
            )
            for version in VARIANT_QC_RESULT_VERSIONS
        },
    )
