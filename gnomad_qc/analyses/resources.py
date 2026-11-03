"""Script containing gnomAD analysis resources."""

import logging
from typing import Optional, Union

from gnomad.resources.resource_utils import TableResource

logger = logging.getLogger("analysis_resources")
logger.setLevel(logging.INFO)


def _analysis_root(
    name: str,
    version: Optional[str] = None,
    test: bool = False,
    tmp_bucket: str = "gs://gnomad-tmp",
) -> str:
    """
    Get root path to the analysis resources.

    :param name: Name of the analysis.
    :param version: Version of the analysis path to return. Default is None, which
        returns no additional directory nesting.
    :param test: Whether to use a tmp path for analysis of the test dataset instead of
        the full dataset.
    :param tmp_bucket: Tmp bucket to write test output to. Default is
        "gs://gnomad-tmp"; pass a lifecycle-suffixed bucket (e.g.
        "gs://gnomad-tmp-30day") to keep test output around longer.
    :return: Root path of analysis resources.
    """
    version = f"/{version}" if version else ""
    return (
        f"{tmp_bucket}/gnomad_{name}_testing{version}"
        if test
        else f"gs://gnomad/analyses/{name}{version}"
    )


def get_pext(
    name: str = "base_level",
    gtex_version: str = "v10",
    suffix: str = "ht",
    test: bool = False,
) -> Union[TableResource, str]:
    """
    Get pext annotation file.

    :param name: Name of the pext annotation file. Default is "base_level".
        Must be one of ["annotation_level", "base_level", "exomes", "genomes",
        "browser"].
    :param gtex_version: GTEx version. Default is "v10".
    :param suffix: Suffix of the pext annotation file. Default is "ht".
    :param test: Whether to use a tmp path for analysis of the test dataset instead of
        the full dataset. Default is False.
    :return: Pext annotation path or TableResource.
    """
    version = f"gtex_{gtex_version}"
    path = (
        f"{_analysis_root('pext', version=version, test=test)}/gnomad.pext."
        f"{version}.{name}.{suffix}"
    )
    if suffix == "ht":
        return TableResource(path=path)
    else:
        return path


def get_hom_depletion(
    data_type: str = "exomes",
    version: Optional[str] = None,
    test: bool = False,
) -> TableResource:
    """
    Get homozygote depletion analysis Table.

    :param data_type: Dataset the check was run on, one of "exomes" or "genomes".
        Default is "exomes".
    :param version: Version of the analysis path to return. Default is None, which
        returns no additional directory nesting.
    :param test: Whether to use a tmp path for analysis of the test dataset instead of
        the full dataset. Default is False.
    :return: Hom depletion TableResource.
    """
    root = _analysis_root(
        "hom_depletion",
        version=version,
        test=test,
        tmp_bucket="gs://gnomad-tmp-30day",
    )
    return TableResource(path=f"{root}/gnomad.hom_depletion.{data_type}.ht")
