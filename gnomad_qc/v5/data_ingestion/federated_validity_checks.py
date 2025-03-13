"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import json
import logging
from io import StringIO
from typing import Any, Dict, List

import hail as hl

from gnomad_qc.assessment.parse_validity_logs import (
    generate_html_report,
    parse_log_file,
)

from gnomad.assessment.validity_checks import (
    check_missingness_of_struct,
    check_raw_and_adj_callstats,
    check_sex_chr_metrics,
    compute_missingness,
    flatten_missingness_struct,
    sum_group_callstats,
    summarize_variants,
    unfurl_array_annotations,
)

from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.reference_genome import get_reference_genome
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
        raise ValueError(f"JSON validation error: %s, {e.message}")


def validate_ht_fields(ht: hl.Table, config: Dict[str, Any]) -> None:
    """Check that necessary fields defined in the JSON config are present in the Hail Table.

    :param ht: Hail Table.
    :param config: JSON configuration for parameter inputs.
    :return: None.
    """
    indexed_array_annotations = config["indexed_array_annotations"]

    # Check that all neccesaary fields are present in the Table.
    missing_fields = {}

    # Check that specified global annotations are present.
    global_fields = [i for i in indexed_array_annotations.values()]
    missing_global_fields = [i for i in global_fields if i not in ht.globals]
    missing_fields["globals"] = missing_global_fields

    # Check that specified row annotations are present.
    row_fields = [i for i in indexed_array_annotations.keys()] + config[
        "struct_annotations_for_missingness"
    ]
    missing_row_fields = [i for i in row_fields if i not in ht.row]
    missing_fields["rows"] = missing_row_fields

    # Check that freq_annotations_to_sum values are present in the 'freq' struct.
    freq_field_names = list(ht.freq.dtype.element_type)
    freq_annotations = config["freq_annotations_to_sum"] + [config["nhomalt_metric"]]
    missing_freq_fields = [i for i in freq_annotations if i not in freq_field_names]
    missing_fields["missing_freq_fields"] = missing_freq_fields

    # Check that sort_order values are present as keys within freq_meta_expr.
    freq_meta_keys = set(
        key for item in hl.eval(eval(config["freq_meta_expr"])) for key in item.keys()
    )
    missing_sort_order_keys = [
        i for i in config["freq_sort_order"] if i not in freq_meta_keys
    ]
    missing_fields["missing_sort_order_keys"] = missing_sort_order_keys

    if any(missing_fields.values()):
        error_message = "Validation failed. Missing fields:\n" + "\n".join(
            f"{key}: {value}" for key, value in missing_fields.items() if value
        )
        raise ValueError(error_message)
    else:
        logger.info("Validated presence of config fields in the Table.")


def check_missingness(
    ht: hl.Table,
    missingness_threshold: float = 0.5,
    struct_annotations: List[str] = ["grpmax", "fafmax", "histograms"],
) -> None:
    """
    Check for and report the fraction of missing data in the Table.

    :param ht: Input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.50.
    :param struct_annotations: List of struct annotations to check for missingness. Default is ['grpmax', 'fafmax', 'histograms'].
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
    ht: hl.Table,
    freq_meta_expr: hl.expr.ArrayExpression,
    missingness_threshold: float = 0.50,
    struct_annotations_for_missingness: List[str] = ["grpmax", "fafmax", "histograms"],
    freq_annotations_to_sum: List[str] = ["AC", "AN", "homozygote_count"],
    freq_sort_order: List[str] = ["gen_anc", "sex", "group"],
    nhomalt_metric: str = "nhomalt",
    verbose: bool = False,
) -> None:
    """
    Perform validity checks on federated data.

    :param ht: Input Table.
    :param freq_meta_expr: Metadata expression that contains the values of the elements in
        `meta_indexed_expr`. The most often used expression is `freq_meta` to index into
        a 'freq' array (example: ht.freq_meta).
    :param freq_annotations_to_sum: List of annotation fields within `meta_expr` to sum. Default is ['AC', 'AN', 'homozygote_count'].
    :param freq_sort_order: Order in which groupings are unfurled into flattened annotations. Default is ["gen_anc", "sex", "group"].
    :param nhomalt_metric: Name of metric denoting homozygous alternate count. Default is "nhomalt".
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False, show only top values of annotations that fail checks. Default is False.
    :return: None
    """
    # Summarize variants and check that all contigs exist.
    expected_contigs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    logger.info("Summarizing variants and checking contigs...")
    summarize_variants(ht, expected_contigs=expected_contigs)

    # Check for missingness.
    logger.info("Checking for missingness...")
    check_missingness(
        ht,
        missingness_threshold,
        struct_annotations=struct_annotations_for_missingness,
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
        subsets=[""],
        pops={
            i[gen_anc_label_name]
            for i in hl.eval(freq_meta_expr)
            if gen_anc_label_name in i
        },
        groups=["adj"],
        verbose=verbose,
        sort_order=freq_sort_order,
        delimiter="_",
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
        delimiter="_",
        nhomalt_metric=nhomalt_metric,
    )

    logger.info("Checking raw and adj callstats...")
    check_raw_and_adj_callstats(
        t=ht,
        subsets=[""],
        verbose=verbose,
        delimiter="_",
        metric_first_field=True,
        nhomalt_metric=nhomalt_metric,
    )
    # TODO: consider adding check_global_and_row_annot_lengths, check for raw and adj.


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
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.empty_set(hl.tstr),
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
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.empty_set(hl.tstr),
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
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.empty_set(hl.tstr),
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
            ],
            "faf": [
                hl.struct(faf95=0.001, faf99=0.002),
                hl.struct(faf95=0.0009, faf99=0.0018),
            ],
            "filters": hl.empty_set(hl.tstr),
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
                ],
                "faf": [
                    hl.struct(faf95=0.0012, faf99=0.0025),
                    hl.struct(faf95=0.0010, faf99=0.0020),
                ],
                "filters": hl.empty_set(hl.tstr),
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
                ],
                "faf": [
                    hl.struct(faf95=0.0013, faf99=0.0027),
                    hl.struct(faf95=0.0011, faf99=0.0023),
                ],
                "filters": hl.empty_set(hl.tstr),
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
                    homozygote_count=hl.tint32,
                )
            ),
            faf=hl.tarray(hl.tstruct(faf95=hl.tfloat64, faf99=hl.tfloat64)),
            filters=hl.tset(hl.tstr),
        ),
    )

    # Define global annotation for freq_index_dict.
    freq_index_dict = {
        "adj": 0,
        "raw": 1,
        "afr_adj": 2,
        "amr_adj": 3,
        "adj_XX": 4,
        "adj_XY": 5,
    }
    faf_index_dict = {"adj": 0, "raw": 1}

    freq_meta = [
        {"group": "adj"},
        {"group": "raw"},
        {"gen_anc": "afr", "group": "adj"},
        {"gen_anc": "amr", "group": "adj"},
        {"sex": "XX", "group": "adj"},
        {"sex": "XY", "group": "adj"},
    ]

    ht = ht.annotate_globals(
        freq_index_dict=freq_index_dict,
        faf_index_dict=faf_index_dict,
        freq_meta=freq_meta,
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
    ht = ht.key_by("locus", "alleles")

    return ht


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

    if args.exclude_xnonpar_y_in_logtest and not args.use_logtest_ht:
        raise ValueError(
            "exclude_xnonpar_y_in_logtest can only be used with use_logtest_ht."
        )

    try:
        # Read in config file and validate.
        with hl.hadoop_open(config_path, "r") as f:
            config = json.load(f)

        validate_config(config, schema)

        # TODO: Add resources to intake federated data once obtained.
        ht = public_release(data_type="exomes").ht()

        # Check that fields specified in the config are present in the Table.
        validate_ht_fields(ht=ht, config=config)

        # Confirm Table is using build GRCh38.
        build = get_reference_genome(ht.locus).name
        if build != "GRCh38":
            raise ValueError(f"Reference genome is {build}, not GRCh38!")

        # Filter to test partitions if specified.
        if test_n_partitions:
            logger.info(
                "Filtering to %d partitions and sex chromosomes...", test_n_partitions
            )
            test_ht = ht._filter_partitions(range(test_n_partitions))

            x_ht = hl.filter_intervals(
                ht, [hl.parse_locus_interval("chrX")]
            )._filter_partitions(range(test_n_partitions))

            y_ht = hl.filter_intervals(
                ht, [hl.parse_locus_interval("chrY")]
            )._filter_partitions(range(test_n_partitions))

            ht = test_ht.union(x_ht, y_ht)

        if args.use_logtest_ht:
            logger.info("Using logtest ht...")
            ht = create_logtest_ht(args.exclude_xnonpar_y_in_logtest)

        # Create row annotations for each element of the indexed arrays and their
        # structs.
        annotations = unfurl_array_annotations(ht, config["indexed_array_annotations"])
        ht = ht.annotate(info=ht.info.annotate(**annotations))

        validate_federated_data(
            ht=ht,
            missingness_threshold=config["missingness_threshold"],
            struct_annotations_for_missingness=config[
                "struct_annotations_for_missingness"
            ],
            freq_meta_expr=eval(config["freq_meta_expr"]),
            freq_annotations_to_sum=config["freq_annotations_to_sum"],
            freq_sort_order=config["freq_sort_order"],
            nhomalt_metric=config["nhomalt_metric"],
            verbose=verbose,
        )

        handler.flush()
        log_output = log_stream.getvalue()

        # TODO: Create resource functions when know organization of federated data.
        log_file = args.output_base + ".log"
        output_file = args.output_base + ".html"

        # Write parsed log to html file.
        with hl.hadoop_open(log_file, "w") as f:
            f.write(log_output)

        parsed_logs = parse_log_file(log_file)
        generate_html_report(parsed_logs, output_file)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("federated_validity_checks"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    # Create a mutually exclusive group for --test-n-partitions and --use-test-ht.
    test_group = parser.add_mutually_exclusive_group()

    test_group.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the input (as well as sex chromosomes) for testing purposes. Defaults"
            "to 2 if passed without a value. Cannot be used if --use-logtest-ht is set."
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
        "--config-path",
        help=(
            "Path to JSON config file for defining parameters. Paramters to define are as follows:"
            "missingness_threshold: Float defining upper cutoff for allowed amount of missingness. Missingness above this value will be flagged as 'FAILED'."
            "indexed_array_annotations: Dictionary of indexed array annotations which will be unfurled. Example: {'faf': 'faf_index_dict', 'freq': 'freq_index_dict'}."
            "struct_annotations_for_missingness: List of struct annotations to check for missingness."
            "freq_meta_expr: Metadata expression that contains the values of the elements in `meta_indexed_expr`. The most often used expression"
            "is `freq_meta` to index into a 'freq' array. Example: ht.freq_meta."
            "freq_annotations_to_sum: List of annotation fields within `freq_meta_expr` to sum. Example: ['AC', 'AN', 'homozygote_count']."
            "freq_sort_order: Order in which groupings are unfurled into flattened annotations. Default is ['gen_anc', 'sex', 'group']."
        ),
        required=True,
        type=str,
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

    args = parser.parse_args()
    main(args)
