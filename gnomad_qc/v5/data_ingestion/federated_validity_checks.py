"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import itertools
import json
import logging
from collections import defaultdict
from typing import Dict, List, Union

import hail as hl
from gnomad.assessment.validity_checks import (
    check_array_struct_missingness,
    check_missingness_of_struct,
    compute_missingness,
    compute_and_check_summations,
    flatten_missingness_struct,
    summarize_variants,
    unfurl_array_annotations,
)
from gnomad.resources.grch38.gnomad import public_release

from gnomad_qc.v5.resources.assessment import get_indexed_array_for_missingness_path
from gnomad_qc.v5.resources.basics import get_logging_path

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("federated_validity_checks")
logger.setLevel(logging.INFO)


### DELETE

from typing import Tuple


def filter_meta_array(
    meta_expr: hl.expr.ArrayExpression,
    keys_to_keep: List[str] = None,
    keys_to_exclude: List[str] = None,
    key_value_pairs_to_keep: Dict[str, List[str]] = None,
    key_value_pairs_to_exclude: Dict[str, List[str]] = None,
    keep_combine_operator: str = "and",
    exclude_combine_operator: str = "and",
    combine_operator: str = "and",
    exact_match: bool = False,
) -> hl.expr.ArrayExpression:
    """
    Filter a metadata array expression based on keys and key-value pairs to keep/exclude.

    If `exact_match` is True, the filtering will only be applied to items with exactly
    the specified keys in `keys_to_keep` (and the keys in `key_value_pairs_to_keep`
    if provided). When `key_value_pairs_to_keep` is also provided, the keys in
    `key_value_pairs_to_keep` must also be present in the metadata item. This
    parameter is only relevant when `keys_to_keep` is provided, `combine_operator`
    is "and", and `exact_match` is True.

    :param meta_expr: Metadata array expression to filter.
    :param keys_to_keep: List of keys to keep.
    :param keys_to_exclude: List of keys to exclude.
    :param key_value_pairs_to_keep: Dictionary of key-value pairs to keep.
    :param key_value_pairs_to_exclude: Dictionary of key-value pairs to exclude.
    :param keep_combine_operator: Whether to use "and" or "or" to combine the filtering
        criteria for keys/key-value pairs to keep.
    :param exclude_combine_operator: Whether to use "and" or "or" to combine the
        filtering criteria for keys/key-value pairs to exclude.
    :param combine_operator: Whether to use "and" or "or" to combine the keep and
        exclude filtering criteria.
    :param exact_match: Whether to apply the filtering only to items with exactly the
        specified keys.
    :return: The filtered metadata array expression.
    """
    keys_to_keep = keys_to_keep or {}
    key_value_pairs_to_keep = key_value_pairs_to_keep or {}
    keys_to_exclude = keys_to_exclude or {}
    key_value_pairs_to_exclude = key_value_pairs_to_exclude or {}

    combine_operator_map = {"and": hl.all, "or": hl.any}
    for o in [keep_combine_operator, exclude_combine_operator, combine_operator]:
        if o not in combine_operator_map:
            raise ValueError(
                "The combine operators must be one of 'and' or 'or', but found" f" {o}!"
            )

    # Assign operators to their respective values in the combine_operator_map dict.
    keep_combine_operator = combine_operator_map[keep_combine_operator]
    exclude_combine_operator = combine_operator_map[exclude_combine_operator]
    combine_operator = combine_operator_map[combine_operator]

    def _get_filter(m: hl.DictExpression) -> hl.expr.BooleanExpression:
        """
        Get the filter to apply to the metadata item.

        :param m: Metadata item.
        :return: Filter to apply to the metadata item.
        """
        # If keys_to_keep is provided, filter to only metadata items with the specified
        # keys. If exact_match is True, filter to only metadata items with the exact
        # keys specified in keys_to_keep, where any keys in key_value_pairs_to_keep
        # are also present.
        if exact_match:
            keep_filter = [
                hl.set(set(keys_to_keep) | set(key_value_pairs_to_keep.keys()))
                == hl.set(m.keys())
            ]
        else:
            keep_filter = [m.contains(k) for k in keys_to_keep]

        # If key_value_pairs_to_keep is provided, filter to only metadata items with the
        # specified key-value pairs.
        keep_filter += [
            hl.literal(v if isinstance(v, list) else [v]).contains(m.get(k, ""))
            for k, v in key_value_pairs_to_keep.items()
        ]

        # If keys_to_exclude is provided, filter to only metadata items without the
        # specified keys and if key_value_pairs_to_exclude is provided, filter to only
        # metadata items without the specified key-value pairs.
        exclude_filter = [~m.contains(k) for k in keys_to_exclude] + [
            ~hl.literal(v if isinstance(v, list) else [v]).contains(m.get(k, ""))
            for k, v in key_value_pairs_to_exclude.items()
        ]

        filters = []
        if keep_filter:
            filters.append(keep_combine_operator(keep_filter))
        if exclude_filter:
            filters.append(exclude_combine_operator(exclude_filter))

        return combine_operator(filters)

    return meta_expr.filter(lambda m: _get_filter(m))


## DELETE


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


def generate_dict_for_sum_comparisons_from_arrays(
    ht: hl.Table,
    meta_expr: hl.expr.ArrayExpression,
    meta_indexed_expr: Dict[str, hl.expr.ArrayExpression],
    meta_indexed_dict: Dict[str, int],
    annotations_to_sum: List[str] = ["AC", "AN"],
    primary_groupings: Dict[str, str] = {"group": "adj"},
    secondary_groupings: List[str] = ["gen_anc", "sex"],
) -> Dict[str, Dict[str, Union[str, List[str]]]]:
    """
    Generate dictionary depicting sum expectations for indexed metadata expressions.

    Dictionary lists will contain all combinations of secondary groupings, and the annotations for combinations of primary groupings to which these
    sums are expected to match.

    Example: If defaults are used and meta_expr=ht.freq_meta, meta_indexed_expr={'freq': ht.freq}, and meta_indexed_dict=ht.freq_index_dict,
    the final dictionary will contain items such as:

    {'AC_group_adj_gen_anc': {'values_to_sum': ['AC_afr_adj',
               'AC_amr_adj',
               'AC_asj_adj',
               'AC_eas_adj',
               'AC_fin_adj',
               'AC_mid_adj',
               'AC_nfe_adj',
               'AC_remaining_adj',
               'AC_sas_adj'],
              'expected_total': 'AC_adj'},
              AN_group_adj_gen_anc_sex': {'values_to_sum': ['AN_afr_XX_adj',
               'AN_afr_XY_adj',
               'AN_amr_XX_adj',
               'AN_amr_XY_adj',
               'AN_asj_XX_adj',
               'AN_asj_XY_adj',
               'AN_eas_XX_adj',
               'AN_eas_XY_adj',
               'AN_fin_XX_adj',
               'AN_fin_XY_adj',
               'AN_mid_XX_adj',
               'AN_mid_XY_adj',
               'AN_nfe_XX_adj',
               'AN_nfe_XY_adj',
               'AN_remaining_XX_adj',
               'AN_remaining_XY_adj',
               'AN_sas_XX_adj',
               'AN_sas_XY_adj'],
              'expected_total': 'AN_adj'}}

    where keys are the names to use to represent the sum values, 'values_to_sum' contains the list of secondary grouping annotations to sum,
    and 'expected_total' contains the primary grouping annotation to which those sums are expected to equal. Annotation names are based on the
    output of unfurl_array_annotations, where the metric name precedes the annotation name.


    Example of meta_expr in table:
        [{'group': 'adj'},
         {'group': 'raw'},
         {'gen_anc': 'afr', 'group': 'adj'},...]

    Example of meta_indexed_expr in table:
        {'afr_adj': 2,
         'adj': 0,
         'raw': 1,...}

    :param ht: Table containing meta information and array expressions to be summed.
    :param meta_expr: Metadata expression that contains the values of the elements in
        `meta_indexed_expr`. The most often used expression is `freq_meta` to index into
        a 'freq' array (example: ht.freq_meta).
    :param meta_indexed_expr: Dictionary where the keys are the expression name
        and the values are the expressions indexed by the `meta_expr` such as a 'freq'
        array (example: 'freq': ht.freq).
    :param meta_indexed_dict: Dictionary depicting the unfurled meta_expr names as keys, and their index in `meta_expr` as values. (example: ht.freq_index_dict).
    :param annotations_to_sum: List of annotation fields withing `meta_expr` to sum. Default is ['AC', 'AN']).
    :param primary_groupings: Dictionary containing primary grouping keys and the corresponding values to filter to. Default is {'group': 'adj'}.
    :param secondary_groupings: List of secondary grouping keys to filter to. All relevant combinations of these values will be summed. For example, if the secondary grouping is
            ['gen_anc', 'sex'], all values within the 'gen_anc' split will be summed and compared to the totals of the primary groupings, as will all values
                within the 'sex' split, and all values split by both 'gen_anc' and 'sex' if that combination is present in the table. Default is ['gen_anc', 'sex'].
    :return: Dictionary with annotation names as keys, and expressions of annotations to compute the sums for those annotations as values.
    """
    comparison_groups = defaultdict(lambda: {"values_to_sum": []})
    # Reverse meta_indexed_dict so that the indices are the keys and the annotation names are the values.
    reversed_index_dict = {
        value: key for key, value in hl.eval(meta_indexed_dict).items()
    }

    for n_primary_combination in range(1, len(primary_groupings) + 1):
        # Iterate through each possible combination of the primary groupings.
        for primary_combination in itertools.combinations(
            primary_groupings.items(), n_primary_combination
        ):
            if primary_combination:
                # Generate name and dictionary for each of the primary grouping combinations.
                primary_grouping_dict = dict(primary_combination)
                primary_grouping_name = "_".join(
                    [f"{key}_{value}" for key, value in primary_grouping_dict.items()]
                )

                # Filter Table metadata to only data relevant to the primary grouping.
                meta = filter_meta_array(
                    meta_expr=meta_expr,
                    key_value_pairs_to_keep=primary_grouping_dict,
                    exact_match=True,
                )

                if hl.eval(hl.len(meta)) == 0:
                    logger.warning(
                        f"Skipping primary grouping combination '{primary_grouping_name}' due to no matching entries."
                    )
                    continue

                for group_meta in hl.eval(meta):
                    # Obtain the index of the specified grouping in the meta_expr array.
                    index_in_meta = hl.eval(meta_expr.index(group_meta))
                    # Obtain the name of the annotation that corresponds to the index value in meta_indexed_dict.
                    primary_group_annotation = reversed_index_dict[index_in_meta]

                for n_secondary_combination in range(1, len(secondary_groupings) + 1):
                    # Iterate through each of the possible combinations of the secondary groupings.
                    for secondary_combination in itertools.combinations(
                        secondary_groupings, n_secondary_combination
                    ):
                        # Generate name for the primary + secondary grouping combination.
                        grouping_combination_name = (
                            f"{primary_grouping_name}_{'_'.join(secondary_combination)}"
                        )

                        # Filter arrays by the grouping combinations.
                        logger.info(
                            f"Applying filtering criteria: %s and %s to find fields to sum for '%s' annotation",
                            primary_grouping_dict,
                            secondary_combination,
                            grouping_combination_name,
                        )
                        meta = filter_meta_array(
                            meta_expr=meta_expr,
                            key_value_pairs_to_keep=primary_grouping_dict,
                            keys_to_keep=secondary_combination,
                            exact_match=True,
                        )

                        if hl.eval(hl.len(meta)) == 0:
                            logger.warning(
                                f"Skipping secondary grouping combintation '%s' due to no matching entries.",
                                grouping_combination_name,
                            )
                        else:
                            for m in hl.eval(meta):
                                # Obtain the index of the specified grouping in the meta_expr array.
                                index_in_meta = hl.eval(meta_expr.index(m))
                                # Obtain the name of the annotation that corresponds to the index value in meta_indexed_dict.
                                secondary_key_name = reversed_index_dict[index_in_meta]

                                # For each annotation to sum, add to the comparison_groups dict a list of the names of the fields to be
                                # summed within 'values_to_sum' dict.
                                # Note that the formatting is based on the table structure after running unfurl_array_annotations, where the metric name
                                # will precede the annotation name.
                                for annotation in annotations_to_sum:
                                    secondary_metric_name = (
                                        f"{annotation}_{secondary_key_name}"
                                    )
                                    comparison_groups[
                                        f"{annotation}_{grouping_combination_name}"
                                    ]["values_to_sum"].append(secondary_metric_name)

                        # For each annotation to sum, add to the comparison_groups dict the name of the annotation to which the sums should total within the 'expected_total' dict.
                        for annotation in annotations_to_sum:
                            comparison_groups[
                                f"{annotation}_{grouping_combination_name}"
                            ][
                                "expected_total"
                            ] = f"{annotation}_{primary_group_annotation}"
    return comparison_groups


def validate_federated_data(
    ht,
    meta_expr_for_summations: hl.expr.ArrayExpression,
    meta_indexed_expr_for_summations: Dict[str, hl.expr.ArrayExpression],
    meta_indexed_dict_for_summations: Dict[str, int],
    annotations_to_sum: List[str] = ["AC", "AN"],
    primary_groupings_for_summations: Dict[str, str] = {"group": "adj"},
    secondary_groupings_for_summations: List[str] = ["gen_anc", "sex"],
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
    :param meta_expr_for_summations: Metadata expression that contains the values of the elements in
        `meta_indexed_expr`. The most often used expression is `freq_meta` to index into
        a 'freq' array (example: ht.freq_meta).
    :param meta_indexed_expr_for_summations: Dictionary where the keys are the expression name
        and the values are the expressions indexed by the `meta_expr` such as a 'freq'
        array (example: 'freq': ht.freq}).
    :param meta_indexed_dict_for_summations: Dictionary depicting the unfurled meta_expr names as keys, and their index in `meta_expr` as values. (example: ht.freq_index_dict).
    :param annotations_to_sum: List of annotation fields withing `meta_expr` to sum. Default is ['AC', 'AN']).
    :param primary_groupings_for_summations: Dictionary containing primary grouping keys and the corresponding values to filter to. Default is {'group': 'adj'}.
    :param secondary_groupings_for_summations: List of secondary grouping keys to filter to. All relevant combinations of these values will be summed. For example, if the secondary grouping is
            ['gen_anc', 'sex'], all values within the 'gen_anc' split will be summed and compared to the totals of the primary groupings, as will all values
                within the 'sex' split, and all values split by both 'gen_anc' and 'sex' if that combination is present in the table. Default is ['gen_anc', 'sex'].
    :return: None
    """
    # Summarize variants and check that all contig exist.
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

    # Check that subset total sum to expected totals.
    logger.info("Checking summations...")
    comparison_groups = generate_dict_for_sum_comparisons_from_arrays(
        ht=ht,
        meta_expr=meta_expr_for_summations,
        meta_indexed_expr=meta_indexed_expr_for_summations,
        meta_indexed_dict=meta_indexed_dict_for_summations,
        annotations_to_sum=annotations_to_sum,
        primary_groupings=primary_groupings_for_summations,
        secondary_groupings=secondary_groupings_for_summations,
    )

    summation_mismatch_results = compute_and_check_summations(ht, comparison_groups)
    for summation, mismatch_count in summation_mismatch_results.items():
        if mismatch_count == 0:
            logger.info(
                "PASSED summation check for %s: mismatch for %d rows",
                summation,
                mismatch_count,
            )
        else:
            logger.info(
                "FAILED summation check for %s: mismatch for %d rows",
                summation,
                mismatch_count,
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
    config_path = args.config_path

    try:
        # TODO: Add resources to intake federated data once obtained.
        ht = public_release(data_type="exomes").ht()

        if test_n_partitions:
            logger.info("Filtering to %d partitions.", test_n_partitions)
            ht = ht._filter_partitions(range(test_n_partitions))

        # Read in parameters from config file.
        with hl.hadoop_open(config_path, "r") as f:
            config = json.load(f)

        missingness_threshold = config["missingness_threshold"]
        indexed_array_annotations = config["indexed_array_annotations"]
        struct_annotations_for_missingness = config[
            "struct_annotations_for_missingness"
        ]
        indexed_array_annotations_for_missingness_dict = config[
            "indexed_array_annotations_for_missingness_dict"
        ]
        # TODO: Allow more than one array for summations
        meta_expr_for_summations = eval(config["meta_expr_for_summations"])
        meta_indexed_expr_for_summations = eval(
            config["meta_indexed_expr_for_summations"]
        )
        meta_indexed_dict_for_summations = eval(
            config["meta_indexed_dict_for_summations"]
        )
        annotations_to_sum = config["annotations_to_sum"]
        primary_groupings_for_summations = config["primary_groupings_for_summations"]
        secondary_groupings_for_summations = config[
            "secondary_groupings_for_summations"
        ]

        # Create row annotations for each element of the indexed arrays and their structs.
        annotations = unfurl_array_annotations(ht, indexed_array_annotations)
        ht = ht.annotate(**annotations)

        validate_federated_data(
            ht=ht,
            missingness_threshold=missingness_threshold,
            struct_annotations_for_missingness=struct_annotations_for_missingness,
            indexed_array_annotations_for_missingness_dict=indexed_array_annotations_for_missingness_dict,
            meta_expr_for_summations=meta_expr_for_summations,
            meta_indexed_expr_for_summations=meta_indexed_expr_for_summations,
            meta_indexed_dict_for_summations=meta_indexed_dict_for_summations,
            annotations_to_sum=annotations_to_sum,
            primary_groupings_for_summations=primary_groupings_for_summations,
            secondary_groupings_for_summations=secondary_groupings_for_summations,
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
        "--config-path",
        help=(
            "Path to JSON config file for defining parameters. Paramters to define are as follows:"
            "missingness_threshold: Float defining upper cutoff for allowed amount of missingness. Missingness above this value will be flagged as 'FAILED'."
            "indexed_array_annotations: Dictionary of indexed array annotations which will be unfurled. Example: {'faf': 'faf_index_dict', 'freq': 'freq_index_dict'}."
            "struct_annotations_for_missingness: List of struct annotations to check for missingness."
            "indexed_array_annotations_for_missingness_dict: Dictionary of indexed array struct annotations to check for missingness, where keys"
            "are the names of the annotations, and values are the names of the globals containing the mapping of group name to index for that key."
            "Example: {'faf':'faf_index_dict', 'freq':'freq_index_dict'}."
            "meta_expr_for_summations: Metadata expression that contains the values of the elements in `meta_indexed_expr`. The most often used expression"
            "is `freq_meta` to index into a 'freq' array. Example: ht.freq_meta."
            "meta_indexed_expr_for_summations: Dictionary where the keys are the expression name and the values are the expressions indexed by the `meta_expr`"
            "such as a 'freq' array. Example: 'freq': ht.freq."
            "meta_indexed_dict_for_summations: Dictionary depicting the unfurled meta_expr names as keys, and their index in `meta_expr` as values. Example: ht.freq_index_dict."
            "annotations_to_sum: List of annotation fields within `meta_expr` to sum. Example: ['AC', 'AN']."
            "primary_groupings_for_summations: Dictionary containing primary grouping keys and the corresponding values to filter to for the summation check. Example: {'group': 'adj'}."
            "secondary_groupings_for_summations: List of secondary grouping keys to filter to for the summation check. All relevant combinations of these values will be summed. Example: ['gen_anc', 'sex']."
        ),
        type=str,
    )

    args = parser.parse_args()
    main(args)
