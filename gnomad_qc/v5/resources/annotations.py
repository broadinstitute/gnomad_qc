"""v5 annotation resources (freq, age hists, group membership)."""

from typing import Optional

from gnomad.resources.resource_utils import TableResource, VersionedTableResource

from .constants import ANNOTATION_VERSIONS, CURRENT_ANNOTATION_VERSION


def _annotations_root(
    version: str = CURRENT_ANNOTATION_VERSION,
    test: bool = False,
    data_type: str = "genomes",
) -> str:
    """
    Get root path to the v5 annotation files.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for testing instead of production.
    :param data_type: Data type of annotation resource ("genomes" or "exomes").
    :return: Root path of the annotation files.
    """
    return (
        f"gs://gnomad-v5-test/annotations/{data_type}/v{version}"
        if test
        else f"gs://gnomad/v5/annotations/{data_type}/v{version}"
    )


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
