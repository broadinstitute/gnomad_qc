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


def get_callset_truth_data(
    truth_sample: str, mt: bool = True, test: bool = False
) -> Union[VersionedMatrixTableResource, VersionedTableResource]:
    """
    Get resources for the truth sample data that is subset from the full callset.

    If `mt` this will return the truth sample MatrixTable (subset from callset);
    otherwise it returns the merged truth sample Table that includes both the truth
    data and the data from the callset.

    :param str truth_sample: Name of the truth sample.
    :param bool mt: Whether path is for a MatrixTable, default is True.
    :param test: Whether to use a tmp path for variant QC tests.
    :return: Path to callset truth sample MT.
    """
    if mt:
        return VersionedMatrixTableResource(
            CURRENT_VARIANT_QC_VERSION,
            {
                version: MatrixTableResource(
                    f"{_variant_qc_root(version, test=test)}/truth_samples/gnomad.genomes.v{version}.{truth_sample}.mt"
                )
                for version in VARIANT_QC_VERSIONS
            },
        )
    else:
        return VersionedTableResource(
            CURRENT_VARIANT_QC_VERSION,
            {
                version: TableResource(
                    f"{_variant_qc_root(version, test=test)}/truth_samples/gnomad.exomes.v{version}.{truth_sample}.ht"
                )
                for version in VARIANT_QC_VERSIONS
            },
        )
