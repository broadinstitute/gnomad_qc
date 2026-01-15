"""Script for running random forest model on gnomAD v5 variant QC data."""

import argparse
import json
import logging
import uuid
from typing import Any, Dict, List, Tuple, Union

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    get_rf_runs,
    get_run_data,
    load_model,
    pretty_print_runs,
    save_model,
)

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.annotations.generate_variant_qc_annotations import (
    INFO_FEATURES,
    NON_INFO_FEATURES,
)
from gnomad_qc.v5.resources.annotations import get_variant_qc_annotations
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

FEATURES = INFO_FEATURES + NON_INFO_FEATURES
"""Combined list of all features used in the RF model."""

RF_COLS = {
    "label": "rf_label",
    "prediction": "rf_prediction",
    "probability": "rf_probability",
    "train": "rf_train",
}
"""RF column names for training labels, predictions, probabilities, and train/test split."""


def load_variant_qc_annotations_ht(
    test: bool = False,
    environment: str = "rwb",
) -> hl.Table:
    """
    Load the pre-computed variant QC annotations Table for RF training.

    This Table is created by the generate_variant_qc_annotations.py script and contains
    all annotations needed for RF training:
    - Info features (in `info` struct): AS_MQRankSum, AS_pab_max, AS_MQ, AS_QD,
      AS_ReadPosRankSum, AS_SOR, AS_FS
    - Non-info features: variant_type, allele_type, n_alt_alleles, was_mixed, has_star
    - Truth data: hapmap, omni, mills, kgp_phase1_hc
    - Training labels: transmitted_singleton_raw, transmitted_singleton_adj,
      sibling_singleton_raw, sibling_singleton_adj, fail_hard_filters

    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "rwb".
    :return: Table with variant QC annotations for RF training.
    """
    logger.info("Loading variant QC annotations HT...")
    ht = get_variant_qc_annotations(test=test, environment=environment).ht()

    # Extract info features from the info struct for RF training.
    ht = ht.transmute(**{f: ht.info[f] for f in INFO_FEATURES})

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
    :param test: Whether to filter the input Table to chr22 and `test_intervals` for
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


def get_or_create_model_id(
    model_id: str = None,
    test: bool = False,
    environment: str = "rwb",
) -> Tuple[str, str]:
    """
    Get or create a unique model ID for the RF run.

    :param model_id: Existing model ID to use. If None, generates a new unique ID.
    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "rwb".
    :return: Tuple of (model_id, rf_run_path).
    """
    rf_run_path = get_rf_run_path(test=test, environment=environment)
    rf_runs = get_rf_runs(rf_run_path)

    if model_id is None:
        model_id = f"rf_{str(uuid.uuid4())[:8]}"
        while model_id in rf_runs:
            model_id = f"rf_{str(uuid.uuid4())[:8]}"

    return model_id, rf_run_path


def _initialize_hail(args) -> None:
    """
    Initialize Hail with appropriate configuration for the environment.

    :param args: Parsed command-line arguments.
    """
    environment = args.environment
    tmp_dir_days = args.tmp_dir_days

    if environment == "batch":
        batch_kwargs = {
            "backend": "batch",
            "log": get_logging_path("v5_random_forest", environment="batch"),
            "tmp_dir": (
                f"{qc_temp_prefix(environment='batch', days=tmp_dir_days)}random_forest"
            ),
            "gcs_requester_pays_configuration": args.gcp_billing_project,
            "regions": ["us-central1"],
        }
        # Add optional batch configuration parameters.
        for param in [
            "app_name",
            "driver_cores",
            "driver_memory",
            "worker_cores",
            "worker_memory",
        ]:
            value = getattr(args, param, None)
            if value is not None:
                batch_kwargs[param] = value

        hl.init(**batch_kwargs)
    else:
        hl.init(
            log=get_logging_path(
                "v5_random_forest",
                environment=environment,
            ),
            tmp_dir=f"{qc_temp_prefix(environment=environment, days=tmp_dir_days)}random_forest",
        )
    hl.default_reference("GRCh38")


def main(args):
    """Run random forest variant QC pipeline."""
    _initialize_hail(args)

    overwrite = args.overwrite
    test = args.test
    environment = args.environment

    model_id, rf_run_path = get_or_create_model_id(
        model_id=args.model_id,
        test=test,
        environment=environment,
    )

    # Get resource paths.
    vqc_annotations_ht_resource = get_variant_qc_annotations(
        test=test, environment=environment
    )
    rf_training_ht_resource = get_rf_training(
        model_id=model_id, test=test, environment=environment
    )
    rf_model_path = get_rf_model_path(
        model_id=model_id, test=test, environment=environment
    )
    rf_result_ht_resource = get_variant_qc_result(
        model_id=model_id, test=test, environment=environment
    )

    try:
        if args.list_rf_runs:
            logger.info("RF runs:")
            pretty_print_runs(get_rf_runs(rf_run_path))

        if args.train_rf:
            # Check input exists.
            check_resource_existence(
                input_resources={
                    "variant_qc_annotations_ht": vqc_annotations_ht_resource
                },
                output_resources={"rf_training_ht": rf_training_ht_resource},
                overwrite=overwrite,
            )

            logger.info("Loading variant QC annotations...")
            vqc_annotation_ht = load_variant_qc_annotations_ht(
                test=test, environment=environment
            )

            logger.info("Training RF model...")
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
            ht = ht.checkpoint(rf_training_ht_resource.path, overwrite=overwrite)

            logger.info("Adding run to RF run list")
            add_model_to_run_list(ht, model_id, get_rf_runs(rf_run_path), rf_run_path)

            logger.info("Saving RF model")
            save_model(rf_model, rf_model_path, overwrite=overwrite)

        if args.apply_rf:
            # Check inputs exist.
            check_resource_existence(
                input_resources={"rf_training_ht": rf_training_ht_resource},
                output_resources={"rf_result_ht": rf_result_ht_resource},
                overwrite=overwrite,
            )

            logger.info("Applying RF model %s...", model_id)
            rf_ht = rf_training_ht_resource.ht()
            rf_features = hl.eval(rf_ht.features)

            logger.info("Loading variant QC annotations...")
            vqc_annotation_ht = load_variant_qc_annotations_ht(
                test=test, environment=environment
            )
            ht = rf_ht.annotate(**vqc_annotation_ht[rf_ht.key].select(*rf_features))
            ht = apply_rf_model(
                ht,
                rf_model=load_model(rf_model_path),
                features=rf_features,
            )

            logger.info("Finished applying RF model...")
            summary_cols = [
                "tp",
                "fp",
                RF_COLS["train"],
                RF_COLS["label"],
                RF_COLS["prediction"],
            ]
            ht = ht.select(*summary_cols, RF_COLS["probability"])
            ht = ht.annotate_globals(rf_model_id=model_id)
            ht = ht.checkpoint(rf_result_ht_resource.path, overwrite=overwrite)
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
    parser.add_argument(
        "--tmp-dir-days",
        help="Number of days to keep temporary files.",
        default=4,
        type=int,
    )
    batch_args = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
    )

    batch_args.add_argument(
        "--gcp-billing-project",
        help="GCP billing project to use for requester pays buckets (required for --environment=batch).",
        default=None,
    )
    batch_args.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_args.add_argument(
        "--driver-cores",
        type=int,
        default=None,
        help="Number of cores for driver node.",
    )
    batch_args.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_args.add_argument(
        "--worker-cores",
        type=int,
        default=None,
        help="Number of cores for worker nodes.",
    )
    batch_args.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
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

    batch_args = [
        "app_name",
        "driver_cores",
        "driver_memory",
        "worker_cores",
        "worker_memory",
    ]
    provided_batch_args = [arg for arg in batch_args if getattr(args, arg) is not None]
    if provided_batch_args and args.environment != "batch":
        parser.error(
            f"Batch configuration arguments ({', '.join('--' + a.replace('_', '-') for a in provided_batch_args)}) "
            f"require --environment=batch"
        )

    if args.environment == "batch" and args.gcp_billing_project is None:
        parser.error("--gcp-billing-project is required when --environment=batch")

    if not args.model_id and not args.train_rf and args.apply_rf:
        parser.error(
            "Error: --model-id is required when running --apply-rf without running"
            " --train-rf too."
        )

    if args.model_id and args.train_rf:
        parser.error(
            "Error: --model-id and --train-rf are mutually exclusive. --train-rf will"
            " generate a run model ID."
        )

    main(args)
