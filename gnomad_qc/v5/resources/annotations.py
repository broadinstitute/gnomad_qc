"""Script containing annotation related resources."""

from typing import Optional

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
    data_set: str = "aou",
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


def get_freq(
    version: str = CURRENT_ANNOTATION_VERSION,
    data_type: str = "genomes",
    test: bool = False,
    subset: Optional[str] = None,
    finalized: bool = True,
) -> TableResource:
    """
    Get the frequency annotation Table for v5.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource ("genomes" or "exomes").
    :param test: Whether to use a tmp path for testing.
    :param subset: Optional subset ("gnomad", "aou") for separate frequency tables.
    :param finalized: Whether to return the finalized frequency table. Default is True.
    :return: Hail Table containing frequency annotations.
    """
    ht_name = f"gnomad.{data_type}.v{version}.frequencies"

    if subset:
        ht_name += f".{subset}"
        if test:
            ht_name += ".test"
        # Subset-specific frequencies go in temp directory during processing
        if not finalized:
            ht_path = f"{_annotations_root(version, test, data_type)}/temp/{ht_name}.ht"
        else:
            ht_path = f"{_annotations_root(version, test, data_type)}/{ht_name}.ht"
    else:
        # Main merged frequency table
        if finalized:
            ht_name += ".final"
        if test:
            ht_name += ".test"
        ht_path = f"{_annotations_root(version, test, data_type)}/{ht_name}.ht"

    return TableResource(ht_path)


def get_age_hist(
    version: str = CURRENT_ANNOTATION_VERSION,
    data_type: str = "genomes",
    test: bool = False,
    subset: Optional[str] = None,
) -> TableResource:
    """
    Get the age histogram Table aligned with frequency groups.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource ("genomes" or "exomes").
    :param test: Whether to use a tmp path for testing.
    :param subset: Optional subset ("gnomad", "aou") for separate age histogram tables.
    :return: Hail Table containing age histograms.
    """
    ht_name = f"gnomad.{data_type}.v{version}.age_hist"

    if subset:
        ht_name += f".{subset}"

    if test:
        ht_name += ".test"

    ht_path = f"{_annotations_root(version, test, data_type)}/{ht_name}.ht"
    return TableResource(ht_path)


def get_consent_ans(
    version: str = CURRENT_ANNOTATION_VERSION,
    data_type: str = "genomes",
    test: bool = False,
) -> TableResource:
    """
    Get the per-variant AN for consent-withdrawn sample set.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource ("genomes" or "exomes").
    :param test: Whether to use a tmp path for testing.
    :return: Hail Table containing consent withdrawal ANs.
    """
    ht_path = f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.consent_ans.ht"
    return TableResource(ht_path)


def get_group_membership(
    version: str = CURRENT_ANNOTATION_VERSION,
    data_type: str = "genomes",
    subset: Optional[str] = None,
    test: bool = False,
) -> TableResource:
    """
    Get the group membership Table used for frequency calculations.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource ("genomes" or "exomes").
    :param subset: Optional subset ("aou" for All of Us data). If None, returns gnomAD group membership.
    :param test: Whether to use a tmp path for testing.
    :return: Hail Table containing group membership annotations.
    """
    subset_suffix = f".{subset}" if subset else ""
    ht_path = f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.group_membership{subset_suffix}.ht"
    return TableResource(ht_path)
