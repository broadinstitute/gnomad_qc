"""Script for running random forest model on gnomAD v4 variant QC data."""
import argparse
import json
import logging
import sys
import uuid
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.utils.slack import slack_notifications
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
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_variant_qc_annotations
from gnomad_qc.v4.resources.sample_qc import interval_qc_pass
from gnomad_qc.v4.resources.variant_qc import (
    get_rf_model_path,
    get_rf_result,
    get_rf_run_path,
    get_rf_training,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_random_forest")
logger.setLevel(logging.INFO)

FEATURES = [
    "allele_type",
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_SOR",
    "n_alt_alleles",
    "variant_type",
]
LABEL_COL = "rf_label"
PREDICTION_COL = "rf_prediction"
TRAIN_COL = "rf_train"


def train_rf(
    ht: hl.Table,
    test: bool = False,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    transmitted_singletons: bool = False,
    sibling_singletons: bool = False,
    adj: bool = False,
    filter_centromere_telomere: bool = False,
    test_intervals: Union[str, List[str]] = "chr20",
    interval_qc_pass_ht: Optional[hl.Table] = None,
) -> Tuple[hl.Table, Any]:
    """
    Train random forest model using `train_rf_model`.

    :param ht: Table containing annotations needed for RF training.
    :param test: Whether to filter the input Table to chr20 and `test_intervals` for
        test purposes. Default is False.
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
    :param interval_qc_pass_ht: Optional interval QC pass Table that contains an
        'interval_qc_pass' annotation indicating whether the interval passes
        high-quality criteria. This annotation is used to filter the Table before
        running training the RF model. Default is None.
    :return: Input `ht` annotated with training information and the RF model.
    """
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
    if interval_qc_pass_ht is not None:
        logger.info("Filtering to intervals that pass interval QC...")
        rf_ht = rf_ht.filter(interval_qc_pass_ht[rf_ht.locus].pass_interval_qc)

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
        interval_qc_filter=interval_qc_pass_ht is not None,
        test_intervals=test_intervals,
    )

    return rf_ht, rf_model


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
        "compute_info_method",
        "transmitted_singletons",
        "sibling_singletons",
        "adj",
        "filter_centromere_telomere",
        "interval_qc_filter",
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
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :param model_id: Model ID to use for RF model. If not provided, a new model ID will
        be generated.
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC pipeline.
    """
    # If no model ID is supplied, generate one and make sure it doesn't already exist.
    rf_run_path = get_rf_run_path(test=test)
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
                "rf_runs": rf_runs,
                "rf_run_path": rf_run_path,
                "model_id": model_id,
            },
            "annotations/generate_variant_qc_annotations.py --create-variant-qc-annotation-ht": {
                "vqc_annotation_ht": get_variant_qc_annotations()
            },
        },
    )
    # Create resource collection for each step of the variant QC pipeline.
    train_random_forest = PipelineStepResourceCollection(
        "--train-rf",
        input_resources={
            "interval_qc.py --generate-interval-qc-pass-ht": {
                "interval_qc_pass_ht": interval_qc_pass(all_platforms=True)
            },
        },
        output_resources={
            "rf_training_ht": get_rf_training(model_id=model_id, test=test),
            "rf_model_path": get_rf_model_path(model_id=model_id, test=test),
        },
    )
    apply_random_forest = PipelineStepResourceCollection(
        "--apply-rf",
        output_resources={"rf_result_ht": get_rf_result(model_id=model_id, test=test)},
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
        tmp_dir="gs://gnomad-tmp-4day",
    )

    overwrite = args.overwrite
    test = args.test
    compute_info_method = args.compute_info_method

    variant_qc_resources = get_variant_qc_resources(
        test=test,
        overwrite=overwrite,
        model_id=args.model_id,
    )
    rf_runs = variant_qc_resources.rf_runs
    rf_run_path = variant_qc_resources.rf_run_path
    model_id = variant_qc_resources.model_id
    vqc_annotation_ht = variant_qc_resources.vqc_annotation_ht.ht()

    if args.list_rf_runs:
        logger.info(f"RF runs:")
        pretty_print_runs(rf_runs)

    if args.train_rf:
        res = variant_qc_resources.train_rf
        res.check_resource_existence()
        ht = vqc_annotation_ht.annotate(
            **vqc_annotation_ht[f"{compute_info_method}_info"]
        )
        if args.interval_qc_filter:
            interval_qc_pass_ht = res.interval_qc_pass_ht.ht()
        else:
            interval_qc_pass_ht = None
        ht, rf_model = train_rf(
            ht,
            test=test,
            fp_to_tp=args.fp_to_tp,
            num_trees=args.num_trees,
            max_depth=args.max_depth,
            transmitted_singletons=args.transmitted_singletons,
            sibling_singletons=args.sibling_singletons,
            adj=args.adj,
            filter_centromere_telomere=args.filter_centromere_telomere,
            test_intervals=args.test_intervals,
            interval_qc_pass_ht=interval_qc_pass_ht,
        )
        ht = ht.annotate_globals(compute_info_method=compute_info_method)
        ht = ht.checkpoint(res.rf_training_ht.path, overwrite=overwrite)

        logger.info("Adding run to RF run list")
        add_model_to_run_list(ht, model_id, rf_runs, rf_run_path)

        logger.info("Saving RF model")
        save_model(rf_model, res.rf_model_path, overwrite=overwrite)

    if args.apply_rf:
        res = variant_qc_resources.apply_rf
        res.check_resource_existence()
        logger.info(f"Applying RF model {model_id}...")
        rf_ht = res.rf_training_ht.ht()
        ht = vqc_annotation_ht.annotate(
            **vqc_annotation_ht[f"{hl.eval(rf_ht.compute_info_method)}_info"]
        )
        ht = apply_rf_model(
            ht,
            rf_model=load_model(res.rf_model_path),
            features=hl.eval(rf_ht.features),
            label=LABEL_COL,
        )

        logger.info("Finished applying RF model...")
        ht = ht.annotate_globals(rf_model_id=model_id)
        ht = ht.checkpoint(res.rf_result_ht.path, overwrite=overwrite)

        ht = ht.annotate(tp=rf_ht[ht.key].tp, fp=rf_ht[ht.key].fp)
        summary_cols = ["tp", "fp", TRAIN_COL, LABEL_COL, PREDICTION_COL]
        ht.group_by(*summary_cols).aggregate(n=hl.agg.count()).show(-1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
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
        "--compute-info-method",
        help=(
            "Method of computing the INFO score to use for the variant QC features. "
            "Default is 'AS'."
        ),
        default="AS",
        type=str,
        choices=["AS", "quasi", "set_long_AS_missing"],
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
        help="Maxmimum tree depth in the RF model. Default is 5.",
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
    training_params.add_argument(
        "--interval-qc-filter",
        help="Whether interval QC should be applied for RF training.",
        action="store_true",
    )

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

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
