"""Script containing sample QC related resources."""

from typing import Optional

from gnomad.resources.resource_utils import (
    MatrixTableResource,
    TableResource,
    VersionedMatrixTableResource,
    VersionedTableResource,
)

from gnomad_qc.v5.resources.constants import (
    CURRENT_SAMPLE_QC_VERSION,
    SAMPLE_QC_VERSIONS,
)


def get_sample_qc_root(
    version: str = CURRENT_SAMPLE_QC_VERSION,
    test: bool = False,
    data_type="genomes",
    data_set="aou",
) -> str:
    """
    Return path to sample QC root folder.

    :param version: Version of sample QC path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used in sample QC, e.g. "aou" or "hgdp_tgp".
    :return: Root to sample QC path.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/sample_qc/{data_type}/{data_set}"
        if test
        else f"gs://gnomad/v{version}/sample_qc/{data_type}/{data_set}"
    )


######################################################################
# Ancestry inference resources
######################################################################


def _get_ancestry_pca_ht_path(
    part: str,
    version: str = CURRENT_SAMPLE_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "hgdp_tgp",
) -> str:
    """
    Get path to files related to ancestry PCA.

    :param part: String indicating the type of PCA file to return (loadings,
        eigenvalues, or scores).
    :param version: Version of sample QC path to return.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used in sample QC, e.g. "aou" or "hgdp_tgp".
    :return: Path to requested ancestry PCA file.
    """
    return f"{get_sample_qc_root(version, test, data_type, data_set)}/ancestry_inference/gnomad.{data_type}.{data_set}.v{version}.pca_{part}.ht"


def ancestry_pca_loadings(
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "hgdp_tgp",
) -> VersionedTableResource:
    """
    Get the ancestry PCA loadings VersionedTableResource.

    :param test: Whether to use a temp path.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used in sample QC, e.g. "aou" or "hgdp_tgp".
    :return: Ancestry PCA loadings
    """
    return VersionedTableResource(
        CURRENT_SAMPLE_QC_VERSION,
        {
            version: TableResource(
                _get_ancestry_pca_ht_path(
                    "loadings",
                    version,
                    test,
                    data_type,
                    data_set,
                )
            )
            for version in SAMPLE_QC_VERSIONS
        },
    )


def ancestry_pca_scores(
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "hgdp_tgp",
) -> VersionedTableResource:
    """
    Get the ancestry PCA loadings VersionedTableResource.

    :param test: Whether to use a temp path.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used in sample QC, e.g. "aou" or "hgdp_tgp".
    :return: Ancestry PCA loadings
    """
    return VersionedTableResource(
        CURRENT_SAMPLE_QC_VERSION,
        {
            version: TableResource(
                _get_ancestry_pca_ht_path(
                    "scores",
                    version,
                    test,
                    data_type,
                    data_set,
                )
            )
            for version in SAMPLE_QC_VERSIONS
        },
    )


def ancestry_pca_eigenvalues(
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "aou",
) -> VersionedTableResource:
    """
    Get the ancestry PCA loadings VersionedTableResource.

    :param test: Whether to use a temp path.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used in sample QC, e.g. "aou" or "hgdp_tgp".
    :return: Ancestry PCA loadings
    """
    return VersionedTableResource(
        CURRENT_SAMPLE_QC_VERSION,
        {
            version: TableResource(
                _get_ancestry_pca_ht_path(
                    "eigenvalues",
                    version,
                    test,
                    data_type,
                    data_set,
                )
            )
            for version in SAMPLE_QC_VERSIONS
        },
    )


def pop_rf_path(
    version: str = CURRENT_SAMPLE_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "aou",
) -> str:
    """
    Path to RF model used for inferring sample populations.

    :param version: gnomAD Version.
    :param test: Whether the RF assignment was from a test dataset.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used to infer the population, e.g. "aou" or "hgdp_tgp".
    :return: String path to sample pop RF model.
    """
    return f"{get_sample_qc_root(version, test, data_type, data_set)}/ancestry_inference/gnomad.{data_type}.{data_set}.v{version}.pop.RF_fit.pickle"


def get_pop_ht(
    version: str = CURRENT_SAMPLE_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "aou",
):
    """
    Get the TableResource of samples' inferred population for the indicated gnomAD version.

    :param version: Version of pop TableResource to return.
    :param test: Whether to use the test version of the pop TableResource.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used to infer the population, e.g. "aou" or "hgdp_tgp".
    :return: TableResource of sample pops.
    """
    return TableResource(
        f"{get_sample_qc_root(version, test, data_type, data_set)}/ancestry_inference/gnomad.{data_type}.{data_set}.v{version}.pop.ht"
    )


def get_pop_pr_ht(
    version: str = CURRENT_SAMPLE_QC_VERSION,
    test: bool = False,
    data_type: str = "genomes",
    data_set: str = "aou",
):
    """
    Get the TableResource of ancestry inference precision and recall values.

    :param version: Version of pop PR TableResource to return.
    :param test: Whether to use the test version of the pop PR TableResource.
    :param data_type: Data type used in sample QC, e.g. "genomes".
    :param data_set: Data set used to infer the population, e.g. "aou" or "hgdp_tgp".
    :return: TableResource of ancestry inference PR values.
    """
    return TableResource(
        f"{get_sample_qc_root(version, test, data_type, data_set)}/ancestry_inference/gnomad.{data_type}.{data_set}.v{version}..pop_pr.ht"
    )


######################################################################
# Ancestry MT resources
######################################################################

hgdp_tgp_unrelateds_without_outliers_mt = MatrixTableResource(
    "gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/pca_results/unrelateds_without_outliers.mt"
)


def get_union_dense_mt(test: bool = False) -> VersionedMatrixTableResource:
    """
    Get the dense MatrixTableResource at final joint v3 and v4 QC sites.

    :param test: Whether to use a tmp path for a test resource.
    :return: MatrixTableResource of QC sites.
    """
    return VersionedMatrixTableResource(
        CURRENT_SAMPLE_QC_VERSION,
        {
            version: MatrixTableResource(
                f"{get_sample_qc_root(version, test, data_set='union')}/qc_mt/gnomad.union.v{version}.dense.mt"
            )
            for version in SAMPLE_QC_VERSIONS
        },
    )
