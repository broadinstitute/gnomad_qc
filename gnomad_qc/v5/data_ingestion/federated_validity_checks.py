"""Script to perform validity checks on input federated data or final release files."""

import argparse
import importlib
import inspect
import json
import logging
import re
from collections import defaultdict
from copy import deepcopy
from io import StringIO
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from bs4 import BeautifulSoup
from gnomad.assessment.parse_validity_logs import generate_html_report, parse_log_file
from gnomad.assessment.validity_checks import (
    check_global_and_row_annot_lengths,
    check_globals_for_retired_terms,
    check_missingness_of_struct,
    check_raw_and_adj_callstats,
    check_sex_chr_metrics,
    compare_subset_freqs,
    compute_missingness,
    flatten_missingness_struct,
    sum_group_callstats,
    summarize_variant_filters,
    summarize_variants,
    unfurl_array_annotations,
)
from gnomad.resources.resource_utils import VersionedTableResource
from gnomad.utils.filtering import remove_fields_from_constant
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.vcf import ALLELE_TYPE_FIELDS, REGION_FLAG_FIELDS
from jsonschema import validate
from jsonschema.exceptions import ValidationError

from gnomad_qc.v5.configs.validity_inputs_schema import schema
from gnomad_qc.v5.resources.basics import get_logging_path

for handler in logging.root.handlers[:]:
    logging.root.removeHandler(handler)


# Configure the root logger for terminal output.
logging.basicConfig(
    format="%(levelname)s (%(module)s.%(funcName)s %(lineno)d): %(message)s",
    level=logging.INFO,
)

# Create an in-memory stream for logging.
log_stream = StringIO()

logger = logging.getLogger("gnomad.assessment.validity_checks")
logger.setLevel(logging.INFO)

# Create a stream handler for in-memory logging.
memory_handler = logging.StreamHandler(log_stream)
formatter = logging.Formatter(
    "%(levelname)s (%(module)s.%(funcName)s %(lineno)d): %(message)s"
)
memory_handler.setFormatter(formatter)
logger.addHandler(memory_handler)

# Remove original alleles for containing non-releasable alleles.
ALLELE_TYPE_FIELDS = deepcopy(ALLELE_TYPE_FIELDS)
ALLELE_TYPE_FIELDS = remove_fields_from_constant(
    ALLELE_TYPE_FIELDS, ["original_alleles"]
)

ALLELE_TYPE_FIELDS = {
    "exomes": ALLELE_TYPE_FIELDS,
    "genomes": remove_fields_from_constant(ALLELE_TYPE_FIELDS, ["has_star"]),
}

# Drop decoy, still doesn't exist on 38.
REGION_FLAG_FIELDS = deepcopy(REGION_FLAG_FIELDS)
REGION_FLAG_FIELDS = remove_fields_from_constant(
    REGION_FLAG_FIELDS, ["decoy", "nonpar"]
)
REGION_FLAG_FIELDS = {
    "exomes": (
        REGION_FLAG_FIELDS
        + [
            "fail_interval_qc",
            "outside_ukb_capture_region",
            "outside_broad_capture_region",
        ]
    ),
    "genomes": REGION_FLAG_FIELDS,
}


def get_table_kind(lines, header_index) -> str:
    """Determine whether a markdown table corresponds to "global" or "row" fields by scanning upward from the table header line.

    :param lines: The full list of lines from the markdown document.
    :param header_index: The index of the table header line (the line with column names).
    :return: String 'global' if the nearest preceding section marker indicates global fields, 'row' if it indicates row fields, or 'None' if neither is found.
    """
    for j in range(header_index - 1, -1, -1):
        prev_line = lines[j].lower().strip()
        if not prev_line:
            continue
        if "global" in prev_line:
            return "global"
        if "row" in prev_line:
            return "row"
        # Stop if hit a heading.
        if prev_line.startswith("#"):
            break
    return None


def hail_type_from_string(type_str: str) -> Any:
    """
    Convert a type string from the markdown to a Hail type.

    This function expects flattened fields (no nested dicts or nested structs).
    Complex nested types may not be fully supported.

    :param type_str: Type string from markdown text, such as `int32`.
    :return: Hail type represented by the type string.
    """
    type_str = type_str.strip().strip("`")

    if type_str in {"int64"}:
        return hl.tint64
    elif type_str in {"int32"}:
        return hl.tint32
    elif type_str in {"float64"}:
        return hl.tfloat64
    elif type_str in {"float32"}:
        return hl.tfloat32
    elif type_str in {"str", "string"}:
        return hl.tstr
    elif type_str in {"bool", "boolean"}:
        return hl.tbool
    elif type_str in {"locus<GRCh38>"}:
        return hl.tlocus("GRCh38")

    # Handle arrays.
    array_match = re.match(r"array<(.+)>", type_str)
    if array_match:
        inner_type = hail_type_from_string(array_match.group(1))
        return hl.tarray(inner_type)

    # Handle sets.
    set_match = re.match(r"set<(.+)>", type_str)
    if set_match:
        inner_type = hail_type_from_string(set_match.group(1))
        return hl.tset(inner_type)

    # Handle dictionaries.
    dict_match = re.match(r"dict<(.+),\s*(.+)>", type_str)
    if dict_match:
        key_type = hail_type_from_string(dict_match.group(1))
        value_type = hail_type_from_string(dict_match.group(2))
        return hl.tdict(key_type, value_type)

    # Handle structs (if type is struct, it will always be a parent type for
    # which we can skip obtaining the type).
    if type_str.startswith("struct"):
        return hl.tstruct()

    raise ValueError(f"Unrecognized type string: {type_str}")


def is_concrete_type(htype) -> bool:
    """Determine whether a Hail type represents a "concrete" field that should be added to field_types, as opposed to an empty container (such as an empty array or struct).

    Empty structs and arrays are not concrete types. Arrays are interpreted recursively.

    :param htype: Hail type to check (ex: hl.tint32, hl.tarray(hl.tfloat64), hl.tstruct()).
    :return: Bool of whether or not the hail type is considered "concrete".
    """
    if isinstance(htype, hl.tstruct):
        return len(htype.fields) > 0
    elif isinstance(htype, hl.tarray):
        return is_concrete_type(htype.element_type)
    else:
        return True


def parse_field_necessity_from_md(
    md_text: str,
) -> Tuple[Dict[str, str], Dict[str, Dict[str, Any]]]:
    """Create dictionary of field necessity from parsing markdown text.

    :param md_text: Markdown text to parse.
    :return: Dictionary of field names and their necessity, and dictionary split into 'global_field_types' and 'row_field_types' keys, containing field names and their types.
    """
    field_necessities = {}
    field_types = {"global_field_types": {}, "row_field_types": {}}
    lines = md_text.splitlines()
    in_table = False

    parent_types = {}

    for i, line in enumerate(lines):
        # Detect table header to distinguish between global and row fields.
        if line.strip().startswith("| Field") and "Field Necessity" in line:
            in_table = True
            table_kind = get_table_kind(lines, i)
            continue

        # Skip alignment row.
        if in_table and re.match(r"^\|[-| ]+\|$", line):
            continue
        # Process table rows.
        elif in_table and line.strip().startswith("|"):
            parts = [c.strip() for c in line.strip().split("|") if c.strip()]
            field_raw = parts[0]
            field_type = parts[1]
            necessity = parts[-1]

            # Strip HTML tags and extra formatting.
            field = BeautifulSoup(field_raw, "html.parser").get_text()
            field = re.sub(r"[*`]", "", field).strip()

            # Convert the types in string format to hail types.
            field_type = hail_type_from_string(field_type)

            field_parts = field.split(".")
            # Keep track of parent fields in the 'parent_types' dictionary.
            if len(field_parts) == 1:
                parent_types[field] = field_type

            if len(field_parts) > 1:
                parent = field_parts[0]
                parent_type = parent_types[parent]
                # If the parent type of a given field is an array, wrap the field type
                # within tarray.
                if isinstance(parent_type, hl.tarray):
                    field_type = hl.tarray(field_type)

            # Skip recording the field type if the field is a parent type with further
            # nodes (will be array or struct with no defined inner types).
            if is_concrete_type(field_type):
                if table_kind == "global":
                    field_types["global_field_types"][field] = field_type
                elif table_kind == "row":
                    field_types["row_field_types"][field] = field_type

            # Normalize necessity label.
            necessity_clean = necessity.strip().lower()
            if "required" in necessity_clean:
                field_necessities[field] = "required"
            elif "optional" in necessity_clean:
                field_necessities[field] = "optional"
            elif "not needed" in necessity_clean:
                field_necessities[field] = "not_needed"
        # End of table.
        elif in_table and not line.strip():
            in_table = False
        elif in_table and not line.strip().startswith("|"):
            in_table = False

    return field_necessities, field_types


def log_field_validation_results(
    field_issues: Dict[str, Dict[str, List[str]]],
    fields_validated: Dict[str, Dict[str, List[str]]],
    type_issues: List[str],
    types_validated: List[str],
) -> None:
    """
    Log the results of field existence and type validation.

    :param field_issues: Nested dictionary mapping necessity ("required", "optional") and annotation_kind ("row", "global") to list of missing field names.
    :param fields_validated:  Nested dictionary mapping necessity ("required", "optional") and annotation_kind ("row", "global") to list of fields successfully found.
    :param type_issues: List of strings describing fields with incorrect or mismatched types.
    :param types_validated: List of strings describing successful type validations.
    :return: None
    """
    # Log missing fields.
    for necessity, annos in field_issues.items():
        for annotation_type, fields in annos.items():
            if not fields:
                continue

            if necessity == "required":
                msg = "FAILED FIELD VALIDATIONS: missing %s fields: %s" % (
                    annotation_type,
                    ", ".join(sorted(fields)),
                )
                logger.error(msg)

            elif necessity == "optional":
                msg = "MISSING optional %s fields: %s" % (
                    annotation_type,
                    ", ".join(sorted(fields)),
                )
                logger.warning(msg)

            else:
                raise ValueError("necessity must be one of 'required' or 'optional")

    # Log found/validated fields.
    for necessity, annos in fields_validated.items():
        for annotation_type, fields in annos.items():
            if fields:
                logger.info(
                    "Found %s %s fields: %s",
                    necessity,
                    annotation_type,
                    ", ".join(sorted(fields)),
                )

    # Log type validations.
    if type_issues:
        logger.error(
            "Type issues: %s",
            " | ".join(sorted(type_issues)),
        )

    # Log type validations.
    if types_validated:
        logger.info(
            "Validated types: %s",
            " | ".join(sorted(types_validated)),
        )


def validate_config(config: Dict[str, Any], schema: Dict[str, Any]) -> None:
    """Validate JSON config inputs.

    :param config: JSON configuration for parameter inputs.
    :param schema: JSON schema to use for validation.
    :return: None.
    """
    # Validate config file against schema.
    try:
        validate(instance=config, schema=schema)
        logger.info("JSON is valid.")
    except ValidationError as e:
        raise ValueError(f"JSON validation error: {e.message}")


def validate_config_fields_in_ht(ht: hl.Table, config: Dict[str, Any]) -> None:
    """Check that necessary fields defined in the JSON config are present in the Hail Table.

    :param ht: Hail Table.
    :param config: JSON configuration for parameter inputs.
    :return: None.
    """
    # Define array struct annotations (must include frequency information as 'freq' annotation within 'freq_fields', and
    # if filtering allele frequency is present it should be provided as 'faf'
    # annotation within 'faf_fields').
    array_struct_annotations = [config["freq_fields"]["freq"]]

    if config.get("faf_fields"):
        array_struct_annotations.append(config["faf_fields"]["faf"])

    # Check that all necessary fields are present in the Table.
    missing_fields = {}

    # Check that specified global annotations are present.
    global_fields = [
        config["freq_fields"]["freq_meta"],
        config["freq_fields"]["freq_meta_sample_count"],
    ]

    if config.get("faf_fields"):
        global_fields.append(config["faf_fields"]["faf_meta"])

    missing_global_fields = [i for i in global_fields if i not in ht.globals]
    missing_fields["globals"] = missing_global_fields

    # Check that specified row annotations are present.
    row_fields = array_struct_annotations + config.get(
        "struct_annotations_to_skip_missingness"
    )

    missing_row_fields = [i for i in row_fields if i not in ht.row]
    missing_fields["rows"] = missing_row_fields

    # Check that specified info annotations are present when configured.
    if config.get("check_mono_and_only_het"):
        info_annotations = ["monoallelic", "only_het"]
        info_fields = list(ht.info.dtype)
        missing_info_fields = [f for f in info_annotations if f not in info_fields]
        missing_fields["missing_info_fields"] = missing_info_fields
    else:
        missing_fields["missing_info_fields"] = []

    # Check that freq_annotations_to_sum values are present in the 'freq' struct.
    freq_fields = list(ht.freq.dtype.element_type)
    freq_annotations = config["freq_annotations_to_sum"] + [config["nhomalt_metric"]]
    missing_freq_fields = [i for i in freq_annotations if i not in freq_fields]
    missing_fields["missing_freq_fields"] = missing_freq_fields

    # Check that sort_order values are present as keys within freq_meta_expr.
    freq_meta_list = hl.eval(ht[config["freq_fields"]["freq_meta"]])
    freq_meta_keys = set(key for item in freq_meta_list for key in item.keys())

    missing_sort_order_keys = [
        i for i in config["sort_order"] if i not in freq_meta_keys
    ]
    missing_fields["missing_sort_order_keys"] = missing_sort_order_keys

    # Check that specified subsets are present as values within the
    # freq_meta_expr subset key.
    subset_values = {i["subset"] for i in freq_meta_list if "subset" in i}
    subsets = [subset for subset in config.get("subsets", []) if subset]
    missing_subsets = set(subsets) - subset_values
    missing_fields["missing_subsets"] = missing_subsets

    if any(missing_fields.values()):
        error_message = "Validation failed. Missing fields:\n" + "\n".join(
            f"{key}: {value}" for key, value in missing_fields.items() if value
        )
        raise ValueError(error_message)
    else:
        logger.info("Validated presence of config fields in the Table.")


def validate_required_fields(
    ht: hl.Table,
    field_types: Dict[str, Dict[str, Any]],
    field_necessities: Dict[str, str],
    validate_all_fields: bool = False,
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """
    Validate that the table contains the required global and row fields and that their values are of the expected types.

    .. note::

        Required fields can be nested (e.g., 'info.QD' indicates that the 'QD' field is nested within the 'info' struct).

    :param ht: Table to validate.
    :param field_types: Nested dictionary of both global and row fields and their expected types. There are two keys: "global_field_types" and "row_field_types", respectively containing the global and row fields as keys and their expected types as values.
    :param field_necessities: Flat dictionary with annotation fields as keys and values field necessity("required" or "optional") as values.
    :param validate_all_fields: Whether to validate all fields or only the required/optional ones.
    :return: Tuple of fields checked and whether or not they passed validation checks.
    """
    field_issues = defaultdict(lambda: defaultdict(list))
    type_issues = []
    fields_validated = defaultdict(lambda: defaultdict(list))
    types_validated = []

    default_necessity = "not_needed"

    # If validate_all_fields is True, set all field necessities to "required".
    if validate_all_fields:
        field_necessities = {k: "required" for k in field_necessities}
        default_necessity = "required"

    def _check_field_exists_and_type(
        root_expr: hl.expr.Expression,
        field_path: str,
        expected_type: Any,
        annotation_kind: str = "row",
    ) -> None:
        """
        Check that the field exists and is of the expected type.

        :param root_expr: Root expression to check.
        :param field_path: Path to the field to check, where a period indicates a nested field.
        :param expected_type: Expected type of the field.
        :param annotation_kind: Kind of annotation to check ("row" or "global").
        :return: None.
        """
        parts = field_path.split(".")
        current_field = root_expr

        field_necessity = field_necessities.get(field_path, default_necessity)

        # Unless specified, only check optional and required fields.
        if field_necessity == "not_needed":
            return

        for i, part in enumerate(parts):
            dtype = current_field.dtype

            if isinstance(dtype, hl.tstruct):
                if part not in dtype.fields:
                    field_issues[field_necessity][annotation_kind].append(field_path)
                    return
                current_field = current_field[part]

            elif isinstance(dtype, hl.tarray) and isinstance(
                dtype.element_type, hl.tstruct
            ):
                if part not in dtype.element_type.fields:
                    field_issues[field_necessity][annotation_kind].append(
                        f"{field_path} (available in array struct: {list(dtype.element_type.fields.keys())})"
                    )
                    return
                current_field = current_field.map(lambda x: x[part])

            else:
                logger.info(
                    "Unsupported type while traversing %s: %s",
                    ".".join(parts[:i]),
                    dtype,
                )
                return

        # If we make it here, the field exists fully.
        fields_validated[field_necessity][annotation_kind].append(field_path)

        # Check that the field type matches expectation.
        if isinstance(expected_type, hl.tarray) and isinstance(
            current_field.dtype, hl.tarray
        ):
            if expected_type.element_type != current_field.dtype.element_type:
                type_issues.append(
                    f"{annotation_kind.capitalize()} field '{field_path}' is an array with incorrect element type: "
                    f"expected {expected_type.element_type}, found {current_field.dtype.element_type}"
                )
                return

            types_validated.append(
                f"{annotation_kind.capitalize()} field '{field_path}' IS an array with correct element type {expected_type.element_type}"
            )
        else:
            if current_field.dtype != expected_type:
                type_issues.append(
                    f"{annotation_kind.capitalize()} field '{field_path}' is not of type {expected_type}, found {current_field.dtype}"
                )
                return

            types_validated.append(
                f"{annotation_kind.capitalize()} field '{field_path}' IS of expected type {expected_type}"
            )

    # Validate global fields.
    for field, expected_type in field_types["global_field_types"].items():
        _check_field_exists_and_type(ht.globals, field, expected_type, "global")

    # Validate row fields.
    for field, expected_type in field_types["row_field_types"].items():
        _check_field_exists_and_type(ht.row, field, expected_type, "row")

    return field_issues, type_issues, fields_validated, types_validated


def check_fields_not_in_requirements(
    ht: hl.Table, field_types: Dict[str, Dict[str, Any]]
) -> None:
    """
    Warn about fields in HT missing from requirements.

    :param ht: Hail Table.
    :param field_types: Nested dictionary of both global and row fields and their expected types. There should be two keys: "global_field_types" and "row_field_types".
    :return: None.
    """

    def _flatten_dtype(dtype: hl.expr.types.HailType, prefix: str = "") -> List[str]:
        """Recursively extract nested names from a Hail DataType."""
        names = []

        # Handle structs.
        if isinstance(dtype, hl.tstruct):
            for field, field_dtype in dtype.items():
                name = f"{prefix}.{field}" if prefix else field
                # Check if this field itself is a struct or container
                names.extend(_flatten_dtype(field_dtype, name))
        # Handle arrays and sets.
        elif isinstance(dtype, (hl.tarray, hl.tset)):
            names.extend(_flatten_dtype(dtype.element_type, prefix))
        # Handle dicts.
        elif isinstance(dtype, hl.tdict):
            names.extend(_flatten_dtype(dtype.value_type, prefix))
        else:
            if prefix:
                names.append(prefix)

        return names

    # Define the mapping between HT components and the requirements dict.
    tasks = [
        ("Global", ht.globals.dtype, "global_field_types"),
        ("Row", ht.row.dtype, "row_field_types"),
    ]

    for label, dtype, req_key in tasks:
        table_fields = set(_flatten_dtype(dtype))
        required_fields = set(field_types.get(req_key, {}).keys())

        unexpected = table_fields - required_fields

        if unexpected:
            logger.warning(
                "%s fields present in Table but missing from requirements: %s",
                label,
                ", ".join(sorted(unexpected)),
            )


def filter_to_test_partitions(
    ht: hl.Table,
    test_n_partitions: int = 2,
) -> hl.Table:
    """
    Filter the Table to a specified number of partitions on autosomes and sex chromosomes for testing purposes.

    :param ht: Input Table.
    :param test_n_partitions: Number of partitions to filter to. Default is 2.
    :return: Filtered Table with only the specified number of partitions.
    """
    test_ht = ht._filter_partitions(range(test_n_partitions))
    x_ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval("chrX")]
    )._filter_partitions(range(test_n_partitions))

    y_ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval("chrY")]
    )._filter_partitions(range(test_n_partitions))

    ht = test_ht.union(x_ht, y_ht)

    return ht


def check_missingness(
    ht: hl.Table,
    missingness_threshold: float = 0.5,
    structs_to_not_traverse: Optional[Tuple[str]] = ("vep",),
) -> None:
    """
    Check for and report the fraction of missing data in row annotations.

    For struct annotations, missingness is checked recursively unless the
    annotation name is included in `structs_to_not_traverse`, in which case
    only top-level missingness of the struct itself is checked.

    :param ht: Input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of
        missingness. Default is 0.50.
    :param structs_to_not_traverse: Optional tuple of top-level struct row
        annotations that should be treated as a single field rather than
        recursively traversed. Default is ("vep",).
    :return: None
    """
    n_sites = ht.count()

    logger.info(
        "Missingness threshold (upper cutoff for allowed missingness): %.2f",
        missingness_threshold,
    )

    metric_missingness = {}
    struct_annotations_checked = []
    non_struct_annotations_checked = []
    non_traversed_struct_annotations = []

    for field, dtype in ht.row.dtype.items():
        field_expr = ht[field]

        if isinstance(dtype, hl.tstruct):
            if field in structs_to_not_traverse:
                non_traversed_struct_annotations.append(field)
                metric_missingness[field] = hl.agg.sum(hl.is_missing(field_expr))
            else:
                struct_annotations_checked.append(field)
                metric_missingness.update(
                    check_missingness_of_struct(field_expr, field)
                )
        else:
            non_struct_annotations_checked.append(field)
            metric_missingness[field] = hl.agg.sum(hl.is_missing(field_expr))

    logger.info(
        "Struct annotations being recursively checked: %s.",
        struct_annotations_checked,
    )
    logger.info(
        "Struct annotations checked only at the top level: %s.",
        non_traversed_struct_annotations,
    )
    logger.info(
        "Non-struct annotations being checked: %s.",
        non_struct_annotations_checked,
    )
    logger.info(
        "Checking missingness for %d annotations.",
        len(metric_missingness),
    )

    output = flatten_missingness_struct(ht.aggregate(hl.struct(**metric_missingness)))

    n_fail = 0
    for field, n_missing in output.items():
        frac_missing = n_missing / n_sites

        if frac_missing > missingness_threshold:
            logger.info(
                "FAILED missingness check for %s: %d sites or %.2f%% missing",
                field,
                n_missing,
                100 * frac_missing,
            )
            n_fail += 1
        else:
            logger.info(
                "Passed missingness check for %s: %d sites or %.2f%% missing",
                field,
                n_missing,
                100 * frac_missing,
            )

    logger.warning("%d missingness checks failed.", n_fail)


def check_missingness_old(
    ht: hl.Table,
    missingness_threshold: float = 0.5,
    struct_annotations_to_skip: Optional[List[str]] = None,
) -> None:
    """
    Check for and report the fraction of missing data in the Table.

    :param ht: Input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.50.
    :param struct_annotations_to_skip: Optional list of top-level struct row
        annotations to skip when automatically selecting struct annotations for
        missingness checks. Default is ['info'].
    :return: None
    """

    # By default, skip `info` because top-level info missingness is checked
    # separately via compute_missingness.
    if struct_annotations_to_skip is None:
        struct_annotations_to_skip = ["info"]

    struct_annotations = [
        field
        for field, dtype in ht.row.dtype.items()
        if isinstance(dtype, hl.tstruct) and field not in struct_annotations_to_skip
    ]

    logger.info("Checking for missingness within struct annotations...")
    logger.info("Struct annotations being checked: %s.", struct_annotations)
    # Determine missingness of each struct annotation.
    metric_missingness = {}
    for metric in struct_annotations:
        metric_missingness.update(check_missingness_of_struct(ht[metric], metric))

    missingness_struct = ht.aggregate(hl.struct(**metric_missingness))
    missingness_dict = flatten_missingness_struct(missingness_struct)

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


def run_row_to_globals_length_check(
    ht: hl.Table,
    config: Dict[str, Any],
    check_all_rows: bool = True,
) -> None:
    """
    Build the row_to_globals_check mapping from config and run check_global_and_row_annot_lengths.

    :param ht: Hail table to check.
    :param config: Configuration dictionary containing freq_fields and optional faf_fields.
    :param check_all_rows: Whether to check all rows. If False, only checks first rows. Default is True.
    :return: None
    """
    row_to_globals_check = {
        config["freq_fields"]["freq"]: [
            config["freq_fields"]["freq_meta"],
            config["freq_fields"]["freq_meta_sample_count"],
        ]
    }
    if config["freq_fields"].get("freq_index_dict"):
        row_to_globals_check[config["freq_fields"]["freq"]].append(
            config["freq_fields"]["freq_index_dict"]
        )
    if config.get("faf_fields"):
        row_to_globals_check[config["faf_fields"]["faf"]] = [
            config["faf_fields"]["faf_meta"],
        ]
        if config["faf_fields"].get("faf_index_dict"):
            row_to_globals_check[config["faf_fields"]["faf"]].append(
                config["faf_fields"]["faf_index_dict"]
            )

    check_global_and_row_annot_lengths(
        t=ht, row_to_globals_check=row_to_globals_check, check_all_rows=check_all_rows
    )


def add_info_annotations(
    ht: hl.Table, region_flag_fields: List[str], allele_type_fields: List[str]
) -> hl.Table:
    """
    Add select annotations to `info` if present in the Table.

    :param ht: Table to annotate.
    :param region_flag_fields: List of region flag fields to check for and add to info if present in the Table.
    :param allele_type_fields: List of allele type fields to check for and add to info if present in the Table.
    :return: Annotated Table with new `info` field.
    """
    info_dict = {}
    missing_region_flags = []

    if "region_flags" in ht.row:
        for field in region_flag_fields:
            if field in ht["region_flags"]:
                info_dict[field] = ht["region_flags"][field]
            else:
                missing_region_flags.append(field)

    if missing_region_flags:
        logger.warning("Missing region_flag fields: %s", missing_region_flags)

    missing_allele_info = []
    if "allele_info" in ht.row:
        for field in allele_type_fields:
            if field in ht["allele_info"]:
                info_dict[field] = ht["allele_info"][field]
            else:
                missing_allele_info.append(field)

    if missing_allele_info:
        logger.warning("Missing allele type fields: %s", missing_allele_info)

    if "monoallelic" in ht.row:
        info_dict["monoallelic"] = ht["monoallelic"]

    if "only_het" in ht.row:
        info_dict["only_het"] = ht["only_het"]

    ht = ht.annotate(info=ht.info.annotate(**info_dict))

    return ht


def validate_federated_data(
    ht: hl.Table,
    freq_meta_expr: hl.expr.ArrayExpression,
    missingness_threshold: float = 0.50,
    struct_annotations_to_skip_missingness: Optional[List[str]] = None,
    freq_annotations_to_sum: List[str] = ["AC", "AN", "homozygote_count"],
    sort_order: List[str] = ["subset", "downsampling", "gen_anc", "sex", "group"],
    nhomalt_metric: str = "nhomalt",
    verbose: bool = False,
    subsets: List[str] = None,
    variant_filter_field: str = "AS_VQSR",
    problematic_regions: List[str] = ["lcr", "non_par", "segdup"],
    site_gt_check_expr: Dict[str, hl.expr.BooleanExpression] = None,
) -> None:
    """
    Perform validity checks on federated data.

    :param ht: Input Table.
    :param freq_meta_expr: Metadata expression that contains the values of the elements in
        `meta_indexed_expr`. The most often used expression is `freq_meta` to index into
        a 'freq' array (example: ht.freq_meta).
    :param freq_annotations_to_sum: List of annotation fields within `meta_expr` to sum. Default is ['AC', 'AN', 'homozygote_count'].
    :param sort_order: Order in which groupings are unfurled into flattened annotations. Default is ["subset", "downsampling", gen_anc", "sex", "group"].
    :param nhomalt_metric: Name of metric denoting homozygous alternate count. Default is "nhomalt".
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False, show only top values of annotations that fail checks. Default is False.
    :param subsets: List of sample subsets.
    :param variant_filter_field: String of variant filtration used in the filters annotation on `ht` (e.g. RF, VQSR, AS_VQSR). Default is "AS_VQSR".
    :param problematic_regions: List of regions considered problematic to run filter check in. Default is ["lcr", "non_par", "segdup"].
    :param site_gt_check_expr: Optional dictionary of strings and boolean expressions typically used to log how many monoallelic or 100% heterozygous sites are in the Table.
    :return: None
    """
    # Summarize variants and check that all contigs exist.
    expected_contigs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    logger.info("Summarizing variants and checking contigs...")
    summarize_variants(ht, expected_contigs=expected_contigs)

    logger.info("Summarizing variant filters...")
    summarize_variant_filters(
        t=ht,
        variant_filter_field=variant_filter_field,
        problematic_regions=problematic_regions,
        single_filter_count=True,
        site_gt_check_expr=site_gt_check_expr,
    )

    # Check for missingness.
    logger.info("Checking for missingness...")
    check_missingness(
        ht,
        missingness_threshold,
        structs_to_not_traverse=struct_annotations_to_skip_missingness,
    )

    # Check that subset totals sum to expected totals.
    logger.info("Checking summations...")
    has_gen_anc = any("gen_anc" in entry for entry in hl.eval(freq_meta_expr))
    has_pop = any("pop" in entry for entry in hl.eval(freq_meta_expr))

    if has_gen_anc and has_pop:
        raise ValueError(
            "Both 'gen_anc' and 'pop' labels found within freq_meta_expr! Only one label can be used."
        )
    elif has_gen_anc:
        gen_anc_label_name = "gen_anc"
    elif has_pop:
        gen_anc_label_name = "pop"
    else:
        raise ValueError(
            "Neither 'gen_anc' nor 'pop' labels found within freq_meta_expr! One label must be used."
        )

    sum_group_callstats(
        t=ht,
        sexes={i["sex"] for i in hl.eval(freq_meta_expr) if "sex" in i},
        subsets=subsets,
        gen_anc_groups={
            i[gen_anc_label_name]
            for i in hl.eval(freq_meta_expr)
            if gen_anc_label_name in i
        },
        groups=["adj"],
        verbose=verbose,
        sort_order=sort_order,
        metric_first_field=True,
        metrics=freq_annotations_to_sum,
        gen_anc_label_name=gen_anc_label_name,
    )

    logger.info("Checking sex chromosomes metrics...")
    info_metrics = list(ht.row.info)
    contigs = ht.aggregate(hl.agg.collect_as_set(ht.locus.contig))
    check_sex_chr_metrics(
        t=ht,
        info_metrics=info_metrics,
        contigs=contigs,
        verbose=verbose,
        nhomalt_metric=nhomalt_metric,
    )

    logger.info("Checking raw and adj callstats...")
    check_raw_and_adj_callstats(
        t=ht,
        subsets=subsets,
        verbose=verbose,
        metric_first_field=True,
        nhomalt_metric=nhomalt_metric,
    )

    if subsets:
        logger.info("Comparing subset frequencies...")
        compare_subset_freqs(
            t=ht,
            subsets=subsets,
            verbose=verbose,
            metric_first_field=True,
            metrics=freq_annotations_to_sum,
        )


def create_logtest_ht(exclude_xnonpar_y: bool = False) -> hl.Table:
    """
    Create a test Hail Table with nested struct annotations to test log output.

    :param exclude_xnonpar_y: If True, exclude chrX non-pseudoautosomal region and chrY variants when making test data. Default is False.
    :return: Table to use for testing log output.
    """
    data = [
        {
            "locus": hl.locus("chr1", 100000, reference_genome="GRCh38"),
            "alleles": ["A", "T"],
            "info": hl.struct(),
            "freq": [
                {"AC": 5, "AF": 0.1, "AN": 20, "homozygote_count": 3},
                {"AC": 10, "AF": 0.05, "AN": 5, "homozygote_count": None},
                {"AC": 15, "AF": 0.07, "AN": 10, "homozygote_count": 2},
                {"AC": 20, "AF": 0.09, "AN": 15, "homozygote_count": 1},
                {"AC": 25, "AF": 0.11, "AN": 18, "homozygote_count": 2},
                {"AC": 30, "AF": 0.13, "AN": 22, "homozygote_count": 4},
                {"AC": 30, "AF": 0.13, "AN": 22, "homozygote_count": 4},
                {"AC": 20, "AF": 0.10, "AN": 20, "homozygote_count": 2},
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.set(["AS_VQSR"]),
            "region_flags": hl.struct(lcr=False, non_par=False, segdup=True),
            "allele_info": hl.struct(allele_type="del", n_alt_alleles=1),
        },
        {
            "locus": hl.locus("chr1", 200000, reference_genome="GRCh38"),
            "alleles": ["C", "G"],
            "info": hl.struct(),
            "freq": [
                {"AC": 6, "AF": 0.08, "AN": 60, "homozygote_count": 4},
                {"AC": 8, "AF": 0.50, "AN": 90, "homozygote_count": 5},
                {"AC": 12, "AF": 0.15, "AN": 50, "homozygote_count": 2},
                {"AC": 18, "AF": 0.18, "AN": 70, "homozygote_count": 3},
                {"AC": 22, "AF": 0.20, "AN": 85, "homozygote_count": 6},
                {"AC": 28, "AF": 0.25, "AN": 95, "homozygote_count": 7},
                {"AC": 24, "AF": 0.20, "AN": 90, "homozygote_count": 1},
                {"AC": 20, "AF": 0.15, "AN": 80, "homozygote_count": 0},
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.set(["AC0"]),
            "region_flags": hl.struct(lcr=False, non_par=False, segdup=False),
            "allele_info": hl.struct(allele_type="snv", n_alt_alleles=1),
        },
        {
            "locus": hl.locus("chr1", 300000, reference_genome="GRCh38"),
            "alleles": ["G", "T"],
            "info": hl.struct(),
            "freq": [
                {"AC": 65, "AF": 0.18, "AN": 200, "homozygote_count": 10},
                {"AC": 88, "AF": 0.20, "AN": 220, "homozygote_count": 12},
                {"AC": 75, "AF": 0.17, "AN": 180, "homozygote_count": 8},
                {"AC": 95, "AF": 0.22, "AN": 250, "homozygote_count": 15},
                {"AC": 100, "AF": 0.24, "AN": 275, "homozygote_count": 18},
                {"AC": 110, "AF": 0.28, "AN": 300, "homozygote_count": 20},
                {"AC": 110, "AF": 0.28, "AN": 300, "homozygote_count": 20},
                {"AC": 100, "AF": 0.20, "AN": 200, "homozygote_count": 2},
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.empty_set(hl.tstr),
            "region_flags": hl.struct(lcr=True, non_par=True, segdup=True),
            "allele_info": hl.struct(allele_type="snv", n_alt_alleles=1),
        },
        {
            "locus": hl.locus("chrX", 400000, reference_genome="GRCh38"),
            "alleles": ["T", "C"],
            "info": hl.struct(),
            "freq": [
                {"AC": 8, "AF": 0.08, "AN": 30, "homozygote_count": 1},
                {"AC": 14, "AF": 0.12, "AN": 40, "homozygote_count": 2},
                {"AC": 22, "AF": 0.14, "AN": 50, "homozygote_count": 4},
                {"AC": 30, "AF": 0.18, "AN": 60, "homozygote_count": 5},
                {"AC": 40, "AF": 0.20, "AN": 75, "homozygote_count": 7},
                {"AC": 50, "AF": 0.25, "AN": 85, "homozygote_count": 9},
                {"AC": 40, "AF": 0.20, "AN": 80, "homozygote_count": 2},
                {"AC": 30, "AF": 0.10, "AN": 70, "homozygote_count": 2},
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.set(["AS_VQSR", "AC0"]),
            "region_flags": hl.struct(lcr=True, non_par=True, segdup=False),
            "allele_info": hl.struct(allele_type="snv", n_alt_alleles=1),
        },
    ]

    if not exclude_xnonpar_y:
        chry_variant = [
            {
                "locus": hl.locus("chrY", 500000, reference_genome="GRCh38"),
                "alleles": ["G", "A"],
                "info": hl.struct(),
                "freq": [
                    {"AC": 12, "AF": 0.10, "AN": 40, "homozygote_count": 1},
                    {"AC": 20, "AF": 0.15, "AN": 50, "homozygote_count": None},
                    {"AC": 35, "AF": 0.18, "AN": 65, "homozygote_count": 3},
                    {"AC": 42, "AF": 0.22, "AN": 70, "homozygote_count": 6},
                    {"AC": 55, "AF": 0.27, "AN": 90, "homozygote_count": None},
                    {"AC": 65, "AF": 0.33, "AN": 100, "homozygote_count": 10},
                    {"AC": 45, "AF": 0.23, "AN": 50, "homozygote_count": 1},
                    {"AC": 40, "AF": 0.20, "AN": 44, "homozygote_count": 2},
                ],
                "faf": [
                    hl.struct(faf95=0.0012, faf99=0.0025),
                    hl.struct(faf95=0.0010, faf99=0.0020),
                ],
                "filters": hl.set(["AS_VQSR", "AC0"]),
                "region_flags": hl.struct(lcr=True, non_par=False, segdup=False),
                "allele_info": hl.struct(allele_type="del", n_alt_alleles=1),
            }
        ]

        chrx_nonpar_variant = [
            {
                "locus": hl.locus("chrX", 22234567, reference_genome="GRCh38"),
                "alleles": ["C", "T"],
                "info": hl.struct(),
                "freq": [
                    {"AC": 15, "AF": 0.12, "AN": 50, "homozygote_count": 2},
                    {"AC": 25, "AF": 0.18, "AN": 60, "homozygote_count": None},
                    {"AC": 40, "AF": 0.22, "AN": 80, "homozygote_count": 5},
                    {"AC": 55, "AF": 0.27, "AN": 100, "homozygote_count": 8},
                    {"AC": 68, "AF": 0.32, "AN": 120, "homozygote_count": None},
                    {"AC": 80, "AF": 0.38, "AN": 140, "homozygote_count": 12},
                    {"AC": 60, "AF": 0.28, "AN": 120, "homozygote_count": 2},
                    {"AC": 50, "AF": 0.20, "AN": 10, "homozygote_count": 5},
                ],
                "faf": [
                    hl.struct(faf95=0.0013, faf99=0.0027),
                    hl.struct(faf95=0.0011, faf99=0.0023),
                ],
                "filters": hl.set(["AC0"]),
                "region_flags": hl.struct(lcr=True, non_par=True, segdup=False),
                "allele_info": hl.struct(allele_type="del", n_alt_alleles=1),
            }
        ]

        data.extend(chry_variant + chrx_nonpar_variant)

    ht = hl.Table.parallelize(
        data,
        hl.tstruct(
            locus=hl.tlocus(reference_genome="GRCh38"),
            alleles=hl.tarray(hl.tstr),
            info=hl.tstruct(),
            freq=hl.tarray(
                hl.tstruct(
                    AC=hl.tint32,
                    AF=hl.tfloat64,
                    AN=hl.tint32,
                    homozygote_count=hl.tint64,
                )
            ),
            faf=hl.tarray(hl.tstruct(faf95=hl.tfloat64, faf99=hl.tfloat64)),
            filters=hl.tset(hl.tstr),
            region_flags=hl.tstruct(lcr=hl.tbool, non_par=hl.tbool, segdup=hl.tbool),
            allele_info=hl.tstruct(allele_type=hl.tstr, n_alt_alleles=hl.tint32),
        ),
    )

    ht = ht.annotate(
        region_flags=ht.region_flags.annotate(
            fail_interval_qc=False,
            outside_ukb_capture_region=False,
            outside_broad_capture_region=False,
        ),
        allele_info=ht.allele_info.annotate(variant_type="snv", was_mixed=False),
    )

    # Define global annotation for freq_index_dict.

    freq_meta = [
        {"group": "adj"},
        {"group": "raw"},
        {"gen_anc": "afr", "group": "adj"},
        {"gen_anc": "amr", "group": "adj"},
        {"sex": "XX", "group": "adj"},
        {"sex": "XY", "group": "adj"},
        {"subset": "non_ukb", "group": "adj"},
    ]

    faf_meta = [{"group": "adj"}, {"group": "raw"}]

    freq_meta_sample_count = [10, 20, 3, 4, 8, 12, 14]
    faf_meta_sample_count = [10, 20]

    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        faf_meta=faf_meta,
        freq_meta_sample_count=freq_meta_sample_count,
        faf_meta_sample_count=faf_meta_sample_count,
        extra_global_field="extra_global_field",
    )

    # Add in retired terms to globals.
    ht = ht.annotate_globals(
        test_meta=[{"group": "adj"}, {"group": "adj", "pop": "oth"}]
    )

    # Add grpmax and fafmax annotations.
    grpmax = hl.struct(
        gnomad=hl.struct(
            AC=hl.or_missing(hl.rand_bool(0.8), hl.max(ht.freq.AC)),
            AF=hl.or_missing(hl.rand_bool(0.2), hl.max(ht.freq.AF)),
            AN=hl.or_missing(hl.rand_bool(0.7), hl.max(ht.freq.AN)),
            homozygote_count=hl.or_missing(
                hl.rand_bool(0.3),
                hl.max(ht.freq.map(lambda x: hl.or_else(x.homozygote_count, 0))),
            ),
            gen_anc=hl.or_missing(
                hl.rand_bool(0.5),
                hl.dict(
                    hl.zip(
                        ht.freq.AC, ht.freq.map(lambda x: x.get("gen_anc", "unknown"))
                    )
                ).get(hl.max(ht.freq.AC)),
            ),
        )
    )
    fafmax = hl.struct(
        gnomad=hl.struct(
            faf95_max=hl.or_missing(hl.rand_bool(0.10), hl.max(ht.faf.faf95)),
            faf95_max_gen_anc=hl.or_missing(
                hl.rand_bool(0.75),
                hl.dict(
                    hl.zip(
                        ht.faf.faf95,
                        ht.freq_meta.map(lambda x: x.get("gen_anc", "unknown")),
                    )
                ).get(hl.max(ht.faf.faf95)),
            ),
            faf99_max=hl.or_missing(hl.rand_bool(0.3), hl.max(ht.faf.faf99)),
            faf99_max_gen_anc=hl.or_missing(
                hl.rand_bool(0.5),
                hl.dict(
                    hl.zip(
                        ht.faf.faf99,
                        ht.freq_meta.map(lambda x: x.get("gen_anc", "unknown")),
                    )
                ).get(hl.max(ht.faf.faf99)),
            ),
        )
    )

    ht = ht.annotate(grpmax=grpmax, fafmax=fafmax)
    # Add monoallelic and only_het annotations.
    ht = ht.annotate(monoallelic=hl.rand_bool(0.50), only_het=hl.rand_bool(0.10))
    ht = ht.key_by("locus", "alleles")

    return ht


def load_gnomad_data(
    gnomad_input_file: str,
    version: str = "5.0",
    data_type: str = "genomes",
    test: bool = False,
    data_set: Optional[str] = None,
    public_release: Optional[bool] = None,
    environment: Optional[str] = None,
) -> hl.Table:
    """
    Load gnomAD data based on specified input file and parameters.

    :param gnomad_input_file: Name of resource to load, either "freq" or "release_sites".
    :param version: Version to load. For example "4.0", "4.1", "5.0". Default is "5.0".
    :param data_type: Type of gnomAD data to load, either "exomes" or "genomes".
    :param test: If True, load test version of the data. Default is False.
    :param data_set: Data set of annotation resource. One of "aou", "gnomad", or "merged". Default is None.
    :param public_release: Whether or not to use the public version of the release. Default is None.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb", "batch", or
        "dataproc". Default is None.
    :return: Hail Table of the specified gnomAD data.
    """
    # Extract the first digit before any dot, ignoring a leading 'v'.
    major_v = version.split(".")[0]

    # Define module mapping based on major version.
    module_mapping = {
        "4": {
            "freq": ("gnomad_qc.v4.resources.annotations", "get_freq"),
            "release_sites": ("gnomad_qc.v4.resources.release", "release_sites"),
        },
        "5": {
            "freq": ("gnomad_qc.v5.resources.annotations", "get_freq"),
            "release_sites": ("gnomad_qc.v5.resources.release", "release_sites"),
        },
    }

    if major_v not in module_mapping:
        raise ValueError(f"Major version {major_v} not supported.")

    if gnomad_input_file not in module_mapping[major_v]:
        raise ValueError(f"Input '{gnomad_input_file}' not found for v{major_v}")

    module_path, function_name = module_mapping[major_v][gnomad_input_file]

    # Import the module and get the function to call.
    module = importlib.import_module(module_path)
    resource_func = getattr(module, function_name)

    logger.info("Loading %s version %s (%s)...", gnomad_input_file, major_v, data_type)

    # Collect all possible params for the function.
    all_params = {
        "data_type": data_type,
        "test": test,
        "version": version,
        "data_set": data_set,
        "public": public_release,
        "environment": environment,
    }

    # Filter to only the parameter that function can accept.
    sig_params = inspect.signature(resource_func).parameters
    valid_args = {
        k: v for k, v in all_params.items() if k in sig_params and v is not None
    }

    logger.info("Using valid parameters %s for function %s", valid_args, function_name)

    # Log which file and params are being used.
    arg_preview = ", ".join([f"{k}={v}" for k, v in valid_args.items()])
    logger.info(f"Calling {module_path}.{function_name}({arg_preview})")

    resource = resource_func(**valid_args)

    # Some resources (e.g. v4 release_sites) return a VersionedTableResource and do
    # not accept a version argument in their function signature. Select the requested
    # version explicitly instead of relying on the resource default.
    if isinstance(resource, VersionedTableResource):
        if version not in resource.versions:
            available_versions = ", ".join(sorted(resource.versions.keys()))
            raise ValueError(
                f"Requested version '{version}' is not available for "
                f"{gnomad_input_file}. Available versions: {available_versions}"
            )

        logger.info(
            "Using resource version '%s' for %s.",
            version,
            gnomad_input_file,
        )
        return resource.versions[version].ht()

    return resource.ht()


def main(args):
    """Perform validity checks for federated data."""
    hl.init(
        log="/federated_validity_checks.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")
    test_n_partitions = args.test_n_partitions
    config_path = args.config_path
    verbose = args.verbose
    output_base = args.output_base

    if args.exclude_xnonpar_y_in_logtest and not args.use_logtest_ht:
        raise ValueError(
            "exclude_xnonpar_y_in_logtest can only be used with use_logtest_ht."
        )

    try:
        # Read in config file and validate.
        with hl.hadoop_open(config_path, "r") as f:
            config = json.load(f)

        validate_config(config, schema)

        data_type = config["data_type"]
        allele_type_fields = ALLELE_TYPE_FIELDS[data_type]
        region_flag_fields = REGION_FLAG_FIELDS[data_type]

        # Read in field necessity markdown file.
        # When submitting hail dataproc job, include "--files field_requirements.md".
        try:
            with open("field_requirements.md", "r") as f:
                md_text = f.read()
        except FileNotFoundError:
            raise FileNotFoundError(
                "Missing required file 'field_requirements.md'.\n"
                "If running a Hail Dataproc job, be sure to include it with the --files argument:\n"
                "  hailctl dataproc submit <cluster-name> --files field_requirements.md  federated_validity_checks.py..."
            )

        field_necessities, field_types = parse_field_necessity_from_md(md_text)

        if args.use_logtest_ht:
            logger.info("Using logtest ht...")
            ht = create_logtest_ht(args.exclude_xnonpar_y_in_logtest)

        else:
            # Load data from the specified gnomAD resource function.
            ht = load_gnomad_data(
                gnomad_input_file=args.gnomad_input_file,
                version=args.gnomad_version,
                data_type=data_type,
                test=args.gnomad_test,
                data_set=args.gnomad_data_set,
                public_release=args.gnomad_public_release,
                environment=args.gnomad_environment,
            )
            output_base = f"{output_base}/{data_type}/{args.gnomad_input_file}"

            # Check that fields specified in the config are present in the Table.
            validate_config_fields_in_ht(ht=ht, config=config)

            # Confirm Table is using build GRCh38.
            build = get_reference_genome(ht.locus).name
            if build != "GRCh38":
                raise ValueError(f"Reference genome is {build}, not GRCh38!")

            if test_n_partitions:
                logger.info(
                    "Filtering to %d partitions and sex chromosomes...",
                    test_n_partitions,
                )
                ht = filter_to_test_partitions(ht, test_n_partitions)

        logger.info("Check that row and global annotations lengths match...")
        run_row_to_globals_length_check(
            ht=ht,
            config=config,
            check_all_rows=not args.check_only_first_rows_to_globals,
        )
        check_globals_for_retired_terms(ht)

        # Create row annotations for each element of the indexed arrays and their
        # structs.
        array_struct_annotations = {
            config["freq_fields"]["freq"]: config["freq_fields"]["freq_meta"],
        }

        if config.get("faf_fields"):
            array_struct_annotations[config["faf_fields"]["faf"]] = config[
                "faf_fields"
            ]["faf_meta"]

        logger.info("Validate required fields...")
        field_issues, type_issues, fields_validated, types_validated = (
            validate_required_fields(
                ht=ht,
                field_types=field_types,
                field_necessities=field_necessities,
                validate_all_fields=args.validate_all_fields,
            )
        )

        log_field_validation_results(
            field_issues, fields_validated, type_issues, types_validated
        )

        check_fields_not_in_requirements(ht, field_types)

        # TODO: Add in lof per person check.
        logger.info("Unfurl array annotations...")
        annotations = unfurl_array_annotations(
            ht=ht,
            array_meta_dicts=array_struct_annotations,
            sorted_keys=config["sort_order"],
        )
        ht = ht.annotate(info=ht.info.annotate(**annotations))

        logger.info("Creating info annotations...")
        ht = add_info_annotations(ht, region_flag_fields, allele_type_fields)

        # If config specifies to check for monoallelic and only heterozygous sites,
        # create the site_gt_check_expr to pass to validate_federated_data.
        if config.get("check_mono_and_only_het"):
            site_gt_check_expr = {
                "monoallelic": ht.info.monoallelic,
                "only_het": ht.info.only_het,
            }
        else:
            site_gt_check_expr = None

        region_flags = [f for f in region_flag_fields if f in ht.info]

        validate_federated_data(
            ht=ht,
            missingness_threshold=args.missingness_threshold,
            struct_annotations_to_skip_missingness=config.get(
                "struct_annotations_to_skip_missingness"
            ),
            freq_meta_expr=ht[config["freq_fields"]["freq_meta"]],
            freq_annotations_to_sum=config["freq_annotations_to_sum"],
            sort_order=config["sort_order"],
            nhomalt_metric=config["nhomalt_metric"],
            verbose=verbose,
            subsets=config["subsets"],
            variant_filter_field=config["variant_filter_field"],
            problematic_regions=region_flags,
            site_gt_check_expr=site_gt_check_expr,
        )

        memory_handler.flush()
        log_output = log_stream.getvalue()

        # TODO: Create resource functions when know organization of federated data.
        log_file = output_base + ".log"
        output_file = output_base + ".html"

        # Write parsed log to html file.
        with hl.hadoop_open(log_file, "w") as f:
            f.write(log_output)

        parsed_logs = parse_log_file(log_file)
        logger.info("Writing html file to %s...", output_file)
        generate_html_report(parsed_logs, output_file)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("federated_validity_checks"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    # Create a mutually exclusive group for --test-n-partitions and --use-test-ht.
    test_group = parser.add_mutually_exclusive_group()

    test_group.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the input (as well as sex chromosomes) for testing purposes. Defaults"
            " to 2 if passed without a value. Cannot be used if --use-logtest-ht is set."
        ),
        nargs="?",
        const=2,
        type=int,
    )
    test_group.add_argument(
        "--use-logtest-ht",
        help="Use a pre-defined Hail Table for testing of log output rather than loading data. Cannot be used if --test-n-partitions is set.",
        action="store_true",
    )
    parser.add_argument(
        "--exclude-xnonpar-y-in-logtest",
        help="Exclude chrX non-pseudoautosomal region and chrY variants when using the logtest data.",
        action="store_true",
    )
    parser.add_argument(
        "--check-only-first-rows-to-globals",
        help="Check only the first row when checking that the lengths of row annotations match the lengths of associated global annotations.",
        action="store_true",
    )
    parser.add_argument(
        "--config-path",
        help=(
            "Path to JSON config file for defining parameters. Parameters to define are as follows:\n"
            " - struct_annotations_to_skip_missingness: Optional list of top-level struct annotations to skip during missingness checks.\n"
            " - freq_fields: Dictionary containing the names of frequency-related fields:\n"
            "      * freq: Name of annotation containing the array of frequency metric objects\n"
            "        corresponding to each frequency metadata group.\n"
            "      * freq_meta: Name of annotation containing allele frequency metadata, an\n"
            "        ordered list containing the frequency aggregation group for each element\n"
            "        of the freq array row annotation, with at least the following groups:\n"
            "        group (adj/raw), gen_anc (inferred genetic ancestry group), and sex\n"
            "        (sex karyotype).\n"
            "      * freq_meta_sample_count: Name of annotation containing sample count per\n"
            "        sample grouping defined in the freq_meta global annotation.\n"
            " - faf_fields: Dictionary containing the names of filtering allele frequency (FAF) related fields:\n"
            "      * faf: Name of annotation containing structs of FAF information.\n"
            "      * faf_meta: Name of annotation for FAF metadata, an ordered list\n"
            "        containing the frequency aggregation group for each element of the faf\n"
            "        arrays, with at least the following groups: group (adj/raw), gen_anc\n"
            "        (inferred genetic ancestry group), and sex (sex karyotype).\n"
            " - freq_annotations_to_sum: List of annotation fields within `freq_meta` to sum. Example: ['AC', 'AN', 'homozygote_count'].\n"
            " - sort_order: Order in which groupings are unfurled into flattened annotations. Default is ['gen_anc', 'sex', 'group'].\n"
            " - nhomalt_metric: Name of metric denoting homozygous alternate count.\n"
            " - subsets: List of sample subsets to include for the subset validity check.\n"
            " - variant_filter_field: String of variant filtration used in the filters annotation of the Hail Table (e.g. 'RF', 'VQSR', 'AS_VQSR').\n"
            " - data_type: Data type to run checks on. One of 'exomes' or 'genomes'.\n"
            " - check_mono_and_only_het: Whether to run the check for monoallelic and 100 percent heterozygous sites in the Table('monoallelic' and 'only_het' annotations must be present)."
        ),
        type=str,
    )
    parser.add_argument(
        "--validate-all-fields",
        help="Validate all fields, regardless of necessity status.",
        action="store_true",
    )
    parser.add_argument(
        "--verbose",
        help="Log successes in addition to failures during validation",
        action="store_true",
    )
    parser.add_argument(
        "--output-base",
        help=(
            "Base path for output files. Will be used to create a log file and an html file."
        ),
        type=str,
        default="gs://gnomad-tmp/federated_validity_checks/federated_validity_checks",
    )
    parser.add_argument(
        "--missingness-threshold",
        help="Float defining upper cutoff for allowed amount of missingness. Missingness above this value will be flagged as 'FAILED'.",
        type=float,
        default=0.50,
    )
    # Create a group for gnomAD input arguments.
    gnomad_group = parser.add_argument_group("gnomad", "gnomAD input options")
    gnomad_group.add_argument(
        "--gnomad-input-file",
        help="Source to load gnomAD data from. 'freq' loads from get_freq and 'release_sites' loads from release_sites. Default is None.",
        choices=["freq", "release_sites"],
        type=str,
        default=None,
    )
    gnomad_group.add_argument(
        "--gnomad-version",
        help="Version of gnomAD resources to use. Default is None.",
        choices=["4.0", "4.1", "4.1.1", "5.0"],
        default=None,
        type=str,
    )
    gnomad_group.add_argument(
        "--gnomad-test",
        help="Load test dataset (smaller subset for testing).",
        action="store_true",
    )
    gnomad_group.add_argument(
        "--gnomad-data-set",
        help="Data set of annotation resource to load, if applicable. One of 'aou', 'gnomad', or 'merged'. Default is None.",
        choices=["aou", "gnomad", "merged"],
        type=str,
        default=None,
    )
    gnomad_group.add_argument(
        "--gnomad-public-release",
        help="Whether or not to use the public version of the release when loading data. Only applicable when loading 'release_sites'.",
        action="store_true",
    )
    gnomad_group.add_argument(
        "--gnomad-environment",
        help=(
            "Environment to use when loading gnomAD data. Must be one of 'rwb', 'batch', or 'dataproc'. Default is None."
        ),
        choices=["rwb", "batch", "dataproc"],
        type=str,
        default=None,
    )
    args = parser.parse_args()
    main(args)
