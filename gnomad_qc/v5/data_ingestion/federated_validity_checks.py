"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl
from gnomad.assessment.validity_checks import (
    check_array_struct_missingness,
    check_missingness_of_struct,
    compute_missingness,
    flatten_missingness_struct,
    summarize_variants,
)
from gnomad.resources.grch38.gnomad import public_release

from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.assessment import get_indexed_array_for_missingness_path

from typing import Dict, List

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("federated_validity_checks")
logger.setLevel(logging.INFO)


def read_file_to_dict(file_path):
    """
    Read in a file of key:value pairs and convert it to a dictionary.

    :param file_path: Path to file to process. The file should have two columns, the first column is the keys and the second column is the values. File should not contain a header.
    """
    ht = hl.import_table(file_path, no_header=True)
    file_dict = dict(hl.tuple([ht.f0, ht.f1]).collect())
    logger.info("Input file converted to the following dictionary: %s.", file_dict)
    return file_dict


def check_missingness(
    ht,
    missingness_threshold: float = 0.5,
    struct_annotations: List[str] = ["grpmax", "fafmax", "histograms"],
    indexed_array_annotations: Dict[str, str] = {
        "faf": "faf_index_dict",
        "freq": "freq_index_dict",
    },
) -> None:
    """
    Check for fraction of missing data in the Table.

    :param ht: Input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.50.
    :param struct_annotations: List of struct annotations to check for missingness. Default is ['grpmax', 'fafmax', 'histograms'].
    :param indexed_array_annotations: Dictionary of indexed array struct annotations to check for missingness, where keys are the names of the annotations, and values are the names of the globals containing the mapping of group name to index for that key. Default is {'faf':'faf_index_dict', 'freq':'freq_index_dict'}.
    :return: None
    """

    logger.info("Checking for missingness within struct annotations...")
    logger.info("Struct annotations being checked: %s.", struct_annotations)
    # Determine missingness of each struct annotation.
    metric_missingness = {}
    for metric in struct_annotations:
        metric_missingness.update(check_missingness_of_struct(ht[metric], metric))

    missingness_struct = ht.aggregate(hl.struct(**metric_missingness))
    missingness_dict = flatten_missingness_struct(missingness_struct)

    logger.info("Checking for missingness within indexed array annotations...")
    logger.info(
        "Indexed array annotations being checked: %s.", indexed_array_annotations
    )
    # Determine missingness of each indexed array annotation.
    missingness_dict.update(
        check_array_struct_missingness(ht, indexed_array_annotations)
    )

    # Report whether or not each metric pass or fails the missingness check
    # based on the missingness_threshold.
    for field, missingness in missingness_dict.items():
        if missingness > missingness_threshold:
            logger.info(
                "FAILED missingness check for %s: %.2f%% missing",
                field,
                100 * missingness,
            )
        else:
            logger.info(
                "Passed missingness check for %s: %.2f%% missing",
                field,
                100 * missingness,
            )

    logger.info("Checking for missingness of info and non-info fields...")
    # Gather info and non-info metrics (or if doesn't exist, set to an empty list)
    # and substract missingness dict.
    info_metrics = (
        set(ht.row.info) - missingness_dict.keys() if "info" in ht.row else set()
    )
    non_info_metrics = set(ht.row) - {"info"} - missingness_dict.keys()
    n_sites = ht.count()
    logger.info("Info metrics are %s", info_metrics)
    logger.info("Non-info metrics are %s", non_info_metrics)
    compute_missingness(
        ht, info_metrics, non_info_metrics, n_sites, missingness_threshold
    )


def validate_federated_data(
    ht,
    missingness_threshold: float = 0.50,
    struct_annotations_for_missingness: List[str] = ["grpmax", "fafmax", "histograms"],
    indexed_array_annotations_for_missingness_dict: Dict[str, str] = {
        "faf": "faf_index_dict",
        "freq": "freq_index_dict",
    },
) -> None:
    """
    Perform validity checks on federated data.

    :param ht: Input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.50.
    :param struct_annotations_for_missingness: List of struct annotations to check for missingness. Default is ['grpmax', 'fafmax', 'histograms'].
    :param indexed_array_annotations_for_missingness_dict: Dictionary of indexed array struct annotations to check for missingness, where keys are the names of the annotations, and values are the names of the globals containing the mapping of group name to index for that key. Default is {'faf':'faf_index_dict', 'freq':'freq_index_dict'}.
    :return: None
    """

    # Summarize variant and check that all contig exist.
    expected_contigs = [
        i
        for i in hl.get_reference("GRCh38").contigs
        if i in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    ]
    logger.info("Summarizing variants and checking contigs...")
    summarize_variants(ht, expected_contigs=expected_contigs)

    # Check for missingness.
    logger.info("Checking for missingness...")
    check_missingness(
        ht,
        missingness_threshold,
        struct_annotations=struct_annotations_for_missingness,
        indexed_array_annotations=indexed_array_annotations_for_missingness_dict,
    )

    # TODO: consider adding check_global_and_row_annot_lengths, check for raw and adj.


def main(args):
    """Perform validity checks for federated data."""
    hl.init(
        log="/federated_validity_checks.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")
    test_n_partitions = args.test_n_partitions
    missingness_threshold = args.missingness_threshold
    struct_annotations_for_missingness = args.struct_annotations_for_missingness
    indexed_array_annotations_for_missingness_dict = read_file_to_dict(
        args.indexed_array_annotations_for_missingness_dict
    )

    # indexed_array_annotations_for_missingness_dict=convert_argparse_to_dict(args.indexed_array_annotations_for_missingness_dict)

    try:
        # TODO: Add resources to intake federated data once obtained.
        ht = public_release(data_type="exomes").ht()

        if test_n_partitions:
            logger.info("Filtering to %d partitions.", test_n_partitions)
            ht = ht._filter_partitions(range(test_n_partitions))

        validate_federated_data(
            ht,
            missingness_threshold,
            struct_annotations_for_missingness,
            indexed_array_annotations_for_missingness_dict,
        )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("federated_validity_checks"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the input for testing purposes. Defaults"
            "to 2 if passed without a value."
        ),
        nargs="?",
        const=2,
        type=int,
    )

    parser.add_argument(
        "--missingness-threshold",
        help=(
            "Upper cutoff for allowed amount of missingness. Missingness above this value will be flagged as 'FAILED'. Default is 0.50."
        ),
        type=float,
        default=0.50,
    )
    parser.add_argument(
        "--struct-annotations-for-missingness",
        help="List of Struct annotation for which to check for missingness. Provided values should be separated by a space, example: 'grpmax fafmax histograms'.",
        type=str,
        nargs="+",
        default=["grpmax", "fafmax", "histograms"],
    )
    parser.add_argument(
        "--indexed-array-annotations-for-missingness-dict",
        help="File containing indexed array annotations for which to check for missingness. This file should not contain a header. "
        "The first column should contain the keys (names of the array annotations), and the second column should contain the values (names of the globals containing "
        "the mapping of group name to index for that key).",
        type=str,
        nargs="+",
        default=get_indexed_array_for_missingness_path(),
    )

    args = parser.parse_args()
    main(args)
