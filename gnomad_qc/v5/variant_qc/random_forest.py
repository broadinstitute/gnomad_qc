"""Script for running random forest model on gnomAD v5 variant QC data."""

import argparse
import json
import logging
import sys
import uuid
from typing import Any, Dict, List, Tuple, Union

import hail as hl
from gnomad.resources.grch38.reference_data import (
    get_truth_ht,
    telomeres_and_centromeres,
)
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    get_rf_runs,
    get_run_data,
    load_model,
    pretty_print_runs,
    save_model,
)

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v5.resources.annotations import (
    get_info_ht,
    get_sib_stats,
    get_trio_stats,
)
from gnomad_qc.v5.resources.basics import get_logging_path, qc_temp_prefix
from gnomad_qc.v5.resources.variant_qc import (
    get_rf_model_path,
    get_rf_run_path,
    get_rf_training,
    get_variant_qc_result,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_random_forest")
logger.setLevel(logging.INFO)

INFO_FEATURES = [
    "AS_MQRankSum",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_SOR",
]
"""List of info features to be used for variant QC."""

NON_INFO_FEATURES = [
    "variant_type",
    "allele_type",
    "n_alt_alleles",
    "has_star",
]
"""List of features to be used for variant QC that are not in the info field."""

FEATURES = INFO_FEATURES + NON_INFO_FEATURES
"""Combined list of all features used in the RF model."""

TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]
"""List of truth datasets to be used for variant QC."""

LABEL_COL = "rf_label"
PREDICTION_COL = "rf_prediction"
PROBABILITY_COL = "rf_probability"
TRAIN_COL = "rf_train"


def get_variant_qc_annotations_ht(
    test: bool = False,
    environment: str = "rwb",
) -> hl.Table:
    """
    Get the variant QC annotations Table for RF training.

    Combines info HT with trio stats, sibling stats, and truth data to create the full
    set of annotations needed for RF training.

    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "rwb".
    :return: Table with variant QC annotations for RF training.
    """
    logger.info("Loading info HT...")
    info_ht = get_info_ht(test=test, environment=environment).ht()

    logger.info("Loading trio stats HT...")
    trio_stats_ht = get_trio_stats(test=test, environment=environment).ht()

    logger.info("Loading sibling stats HT...")
    sib_stats_ht = get_sib_stats(test=test, environment=environment).ht()

    logger.info("Loading truth data HT...")
    truth_data_ht = get_truth_ht()

    # Select relevant fields from trio stats.
    # gnomad library's generate_trio_stats returns n_transmitted_{group} and
    # ac_children_{group} where group is 'raw' or 'adj'.
    trio_stats_ht = trio_stats_ht.select(
        "n_transmitted_raw",
        "n_transmitted_adj",
        "ac_children_raw",
        "ac_children_adj",
    )

    # Join all annotations into a single table.
    ht = info_ht.annotate(
        trio_stats=trio_stats_ht[info_ht.key],
        sib_stats=sib_stats_ht[info_ht.key],
        truth_data=truth_data_ht[info_ht.key],
    )

    # Extract info annotations.
    ht = ht.transmute(
        AS_MQRankSum=ht.info.AS_MQRankSum,
        AS_QD=ht.info.AS_QD,
        AS_ReadPosRankSum=ht.info.AS_ReadPosRankSum,
        AS_SOR=ht.info.AS_SOR,
        AS_FS=ht.info.AS_FS,
        AS_MQ=ht.info.AS_MQ,
    )

    # Compute transmitted singletons from trio stats.
    # A transmitted singleton is when exactly one allele was transmitted (n_transmitted == 1)
    # and the allele count in children is 2 (one from each of two children, or two from one).
    # We use AC from the info HT for the AC filter.
    ac_raw = ht.AC_info.AC_high_quality_raw
    ac_adj = ht.AC_info.AC_high_quality

    ht = ht.annotate(
        # Transmitted singletons.
        transmitted_singleton_raw=hl.or_else(
            (ht.trio_stats.n_transmitted_raw == 1) & (ac_raw == 2), False
        ),
        transmitted_singleton_adj=hl.or_else(
            (ht.trio_stats.n_transmitted_adj == 1) & (ac_adj == 2), False
        ),
        # Sibling singletons - gnomad library's generate_sib_stats returns
        # n_sib_shared_variants_{group}.
        sibling_singleton_raw=hl.or_else(
            (ht.sib_stats.n_sib_shared_variants_raw == 1) & (ac_raw == 2), False
        ),
        sibling_singleton_adj=hl.or_else(
            (ht.sib_stats.n_sib_shared_variants_adj == 1) & (ac_adj == 2), False
        ),
        # Truth data annotations.
        **{td: hl.or_else(ht.truth_data[td], False) for td in TRUTH_DATA},
        # Hard filters - using QD, FS, and MQ thresholds consistent with v4.
        fail_hard_filters=(ht.AS_QD < 2) | (ht.AS_FS > 60) | (ht.AS_MQ < 30),
    )

    # Annotate with allele info needed for features.
    ht = ht.annotate(
        allele_type=hl.if_else(
            hl.is_snp(ht.alleles[0], ht.alleles[1]),
            "snv",
            hl.if_else(
                hl.is_insertion(ht.alleles[0], ht.alleles[1]),
                "ins",
                hl.if_else(
                    hl.is_deletion(ht.alleles[0], ht.alleles[1]),
                    "del",
                    "complex",
                ),
            ),
        ),
        variant_type=hl.if_else(
            hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv", "indel"
        ),
        has_star=ht.alleles.any(lambda a: a == "*"),
        n_alt_alleles=hl.len(ht.alleles) - 1,
    )

    # Drop intermediate structs.
    ht = ht.drop("trio_stats", "sib_stats", "truth_data", "AC_info")

    return ht


def train_rf(
    ht: hl.Table,
    test: bool = False,
    features: List[str] = None,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    transmitted_singletons: bool = False,
    sibling_singletons: bool = False,
    adj: bool = False,
    filter_centromere_telomere: bool = False,
    test_intervals: Union[str, List[str]] = "chr20",
) -> Tuple[hl.Table, Any]:
    """
    Train random forest model using `train_rf_model`.

    :param ht: Table containing annotations needed for RF training.
    :param test: Whether to filter the input Table to chr20 and `test_intervals` for
        test purposes. Default is False.
    :param features: List of features to use in the random forests model. When no
        `features` list is provided, the global FEATURES is used.
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all
        training examples are used. Default is 1.0.
    :param num_trees: Number of trees in the RF model. Default is 500.
    :param max_depth: Maximum tree depth in the RF model. Default is 5.
    :param transmitted_singletons: Whether to use transmitted singletons for training.
        Default is False.
    :param sibling_singletons: Whether to use sibling singletons for training. Default
        is False.
    :param adj: Whether to use adj genotypes for transmitted/sibling singletons instead
        of raw. Default is False and raw is used.
    :param filter_centromere_telomere: Filter centromeres and telomeres before training.
        Default is False.
    :param test_intervals: Specified interval(s) will be held out for testing and
        evaluation only. Default is "chr20".
    :return: Input `ht` annotated with training information and the RF model.
    """
    if features is None:
        features = FEATURES

    if test_intervals:
        if isinstance(test_intervals, str):
            test_intervals = [test_intervals]
        test_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in test_intervals
        ]

    if test:
        logger.info("Filtering to chr22 and evaluation intervals for testing...")
        chr22_interval = [hl.parse_locus_interval("chr22", reference_genome="GRCh38")]
        ht = hl.filter_intervals(ht, chr22_interval + test_intervals)

    logger.info("Annotating true positives and false positives in HT...")
    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills
    transmit_expr = ht.transmitted_singleton_raw
    sibling_expr = ht.sibling_singleton_raw

    if adj:
        transmit_expr = ht.transmitted_singleton_adj
        sibling_expr = ht.sibling_singleton_adj

    if transmitted_singletons:
        tp_expr |= transmit_expr
    if sibling_singletons:
        tp_expr |= sibling_expr

    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    rf_ht = ht
    if filter_centromere_telomere:
        logger.info("Filtering centromeres and telomeres from HT...")
        rf_ht = rf_ht.filter(
            ~hl.is_defined(telomeres_and_centromeres.ht()[rf_ht.locus])
        )

    rf_ht, rf_model = train_rf_model(
        rf_ht,
        rf_features=features,
        tp_expr=rf_ht.tp,
        fp_expr=rf_ht.fp,
        fp_to_tp=fp_to_tp,
        num_trees=num_trees,
        max_depth=max_depth,
        test_expr=hl.literal(test_intervals).any(
            lambda interval: interval.contains(rf_ht.locus)
        ),
    )

    rf_ht = rf_ht.annotate_globals(
        fp_to_tp=fp_to_tp,
        num_trees=num_trees,
        max_depth=max_depth,
        transmitted_singletons=transmitted_singletons,
        sibling_singletons=sibling_singletons,
        adj=adj,
        filter_centromere_telomere=filter_centromere_telomere,
        test_intervals=test_intervals,
    )

    ht = ht.select("tp", "fp").join(rf_ht, how="left")

    return ht, rf_model


def add_model_to_run_list(
    ht: hl.Table,
    model_id: str,
    rf_runs: Dict[str, Any],
    rf_run_path: str,
) -> None:
    """
    Add RF model run to RF run list.

    :param ht: Table containing RF model run information as globals.
    :param model_id: ID of RF model run.
    :param rf_runs: Dictionary containing current RF run information.
    :param rf_run_path: Path to RF run list.
    :return: None
    """
    logger.info("Adding run to RF run list")
    ht = ht.annotate_globals(test_intervals=hl.str(ht.test_intervals))
    ht_globals = hl.eval(ht.globals)
    input_args = [
        "transmitted_singletons",
        "sibling_singletons",
        "adj",
        "filter_centromere_telomere",
    ]
    rf_output = ["test_intervals", "features_importance", "test_results"]
    rf_runs[model_id] = get_run_data(
        input_args={k: ht_globals[k] for k in input_args},
        **{k: ht_globals[k] for k in rf_output},
    )

    with hl.hadoop_open(rf_run_path, "w") as f:
        json.dump(rf_runs, f)


def get_variant_qc_resources(
    test: bool,
    overwrite: bool,
    model_id: str = None,
    environment: str = "rwb",
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :param model_id: Model ID to use for RF model. If not provided, a new model ID will
        be generated.
    :param environment: Environment to use. Default is "rwb".
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC pipeline.
    """
    # If no model ID is supplied, generate one and make sure it doesn't already exist.
    rf_run_path = get_rf_run_path(test=test, environment=environment)
    rf_runs = get_rf_runs(rf_run_path)
    if model_id is None:
        model_id = f"rf_{str(uuid.uuid4())[:8]}"
        while model_id in rf_runs:
            model_id = f"rf_{str(uuid.uuid4())[:8]}"

    # Initialize variant QC pipeline resource collection.
    variant_qc_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc",
        overwrite=overwrite,
        pipeline_resources={
            "RF models": {
                "rf_run_path": rf_run_path,
                "model_id": model_id,
            },
        },
    )
    # Create resource collection for each step of the variant QC pipeline.
    train_random_forest = PipelineStepResourceCollection(
        "--train-rf",
        output_resources={
            "rf_training_ht": get_rf_training(
                model_id=model_id, test=test, environment=environment
            ),
            "rf_model_path": get_rf_model_path(
                model_id=model_id, test=test, environment=environment
            ),
        },
    )
    apply_random_forest = PipelineStepResourceCollection(
        "--apply-rf",
        output_resources={
            "rf_result_ht": get_variant_qc_result(
                model_id=model_id, test=test, environment=environment
            )
        },
        pipeline_input_steps=[train_random_forest],
    )

    # Add all steps to the variant QC pipeline resource collection.
    variant_qc_pipeline.add_steps(
        {
            "train_rf": train_random_forest,
            "apply_rf": apply_random_forest,
        }
    )

    return variant_qc_pipeline


def main(args):
    """Run random forest variant QC pipeline."""
    hl.init(
        default_reference="GRCh38",
        log="/variant_qc_random_forest.log",
        tmp_dir=qc_temp_prefix(environment=args.environment),
    )

    overwrite = args.overwrite
    test = args.test
    environment = args.environment

    variant_qc_resources = get_variant_qc_resources(
        test=test,
        overwrite=overwrite,
        model_id=args.model_id,
        environment=environment,
    )
    rf_run_path = variant_qc_resources.rf_run_path
    model_id = variant_qc_resources.model_id

    try:
        if args.list_rf_runs:
            logger.info("RF runs:")
            pretty_print_runs(get_rf_runs(rf_run_path))

        if args.train_rf:
            res = variant_qc_resources.train_rf
            res.check_resource_existence()

            logger.info("Loading variant QC annotations...")
            vqc_annotation_ht = get_variant_qc_annotations_ht(
                test=test, environment=environment
            )

            ht, rf_model = train_rf(
                vqc_annotation_ht,
                test=test,
                features=args.features,
                fp_to_tp=args.fp_to_tp,
                num_trees=args.num_trees,
                max_depth=args.max_depth,
                transmitted_singletons=args.transmitted_singletons,
                sibling_singletons=args.sibling_singletons,
                adj=args.adj,
                filter_centromere_telomere=args.filter_centromere_telomere,
                test_intervals=args.test_intervals,
            )
            ht = ht.checkpoint(res.rf_training_ht.path, overwrite=overwrite)

            logger.info("Adding run to RF run list")
            add_model_to_run_list(ht, model_id, get_rf_runs(rf_run_path), rf_run_path)

            logger.info("Saving RF model")
            save_model(rf_model, res.rf_model_path, overwrite=overwrite)

        if args.apply_rf:
            res = variant_qc_resources.apply_rf
            res.check_resource_existence()
            logger.info("Applying RF model %s...", model_id)
            rf_ht = res.rf_training_ht.ht()
            rf_features = hl.eval(rf_ht.features)

            logger.info("Loading variant QC annotations...")
            vqc_annotation_ht = get_variant_qc_annotations_ht(
                test=test, environment=environment
            )
            ht = rf_ht.annotate(**vqc_annotation_ht[rf_ht.key].select(*rf_features))
            ht = apply_rf_model(
                ht,
                rf_model=load_model(res.rf_model_path),
                features=rf_features,
            )

            logger.info("Finished applying RF model...")
            summary_cols = ["tp", "fp", TRAIN_COL, LABEL_COL, PREDICTION_COL]
            ht = ht.select(*summary_cols, PROBABILITY_COL)
            ht = ht.annotate_globals(rf_model_id=model_id)
            ht = ht.checkpoint(res.rf_result_ht.path, overwrite=overwrite)
            ht.group_by(*summary_cols).aggregate(n=hl.agg.count()).show(-1)

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            get_logging_path("variant_qc_random_forest", environment=environment)
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            "If the dataset should be filtered to chr22 for testing (also filtered to "
            "evaluation interval specified by --test-intervals)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--model-id",
        help=(
            "Model ID. Created by --train-rf and only needed for --apply-rf without"
            " running --train-rf."
        ),
        required=False,
    )
    parser.add_argument(
        "--environment",
        help="Environment to run in.",
        choices=["rwb", "batch", "dataproc"],
        default="rwb",
    )

    actions = parser.add_argument_group("Actions")
    actions.add_argument(
        "--list-rf-runs",
        help=(
            "Lists all previous RF runs, along with their model ID, parameters and"
            " testing results."
        ),
        action="store_true",
    )
    actions.add_argument("--train-rf", help="Trains RF model.", action="store_true")
    actions.add_argument(
        "--apply-rf", help="Applies RF model to the data.", action="store_true"
    )

    rf_params = parser.add_argument_group("Random Forest parameters")
    rf_params.add_argument(
        "--features",
        help="Features to use in the random forests model.",
        default=FEATURES,
        type=str,
        nargs="+",
    )
    rf_params.add_argument(
        "--fp-to-tp",
        help=(
            "Ratio of FPs to TPs for training the RF model. If 0, all training examples"
            " are used. Default is 1.0."
        ),
        default=1.0,
        type=float,
    )
    rf_params.add_argument(
        "--test-intervals",
        help=(
            "The specified interval(s) will be held out for testing and evaluation"
            ' only. Default is "chr20".'
        ),
        nargs="+",
        type=str,
        default="chr20",
    )
    rf_params.add_argument(
        "--num-trees",
        help="Number of trees in the RF model. Default is 500.",
        default=500,
        type=int,
    )
    rf_params.add_argument(
        "--max-depth",
        help="Maximum tree depth in the RF model. Default is 5.",
        default=5,
        type=int,
    )

    training_params = parser.add_argument_group("Training data parameters")
    training_params.add_argument(
        "--adj",
        help="Use adj genotypes for transmitted/sibling singletons.",
        action="store_true",
    )
    training_params.add_argument(
        "--transmitted-singletons",
        help="Include transmitted singletons in training.",
        action="store_true",
    )
    training_params.add_argument(
        "--sibling-singletons",
        help="Include sibling singletons in training.",
        action="store_true",
    )
    training_params.add_argument(
        "--filter-centromere-telomere",
        help="Train RF without centromeres and telomeres.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if not args.model_id and not args.train_rf and args.apply_rf:
        sys.exit(
            "Error: --model_id is required when running --apply-rf without running"
            " --train-rf too."
        )

    if args.model_id and args.train_rf:
        sys.exit(
            "Error: --model_id and --train-rf are mutually exclusive. --train-rf will"
            " generate a run model ID."
        )

    main(args)
