# noqa: D100

import argparse
import json
import logging
import sys
import uuid
from typing import List, Optional, Union

import hail as hl
from gnomad.resources.grch38.reference_data import (
    telomeres_and_centromeres,
)
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.pipeline import train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    get_rf_runs,
    get_run_data,
    load_model,
    median_impute_features,
    pretty_print_runs,
    save_model,
)

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import (
    get_vqsr_filters,
    get_variant_qc_annotations # since get_info and others are being rolled into one
)

from gnomad_qc.v4.resources.basics import get_checkpoint_path
from gnomad_qc.v4.resources.variant_qc import (
    get_rf_model_path,
    get_rf_result,
    get_rf_training,
    rf_run_path,
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
INFO_FEATURES = [
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_SOR",
]
LABEL_COL = "rf_label"
PREDICTION_COL = "rf_prediction"
TRAIN_COL = "rf_train"
TRUTH_DATA = ["hapmap", "omni", "mills", "kgp_phase1_hc"]

def train_rf(
    ht: hl.Table,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    no_transmitted_singletons: bool = False,
    no_sibling_singletons: bool = False,
    filter_centromere_telomere: bool = False,
    test_intervals: Union[str, List[str]] = "chr20",
):
    """
    Train random forest model using `train_rf_model`.

    :param ht: Table containing annotations needed for RF training, built with ???
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used.
    :param num_trees: Number of trees in the RF model.
    :param max_depth: Maxmimum tree depth in the RF model.
    :param no_transmitted_singletons: Do not use transmitted singletons for training.
    :param no_sibling_singletons: Do not use sibling singletons for training. 
    :param filter_centromere_telomere: Filter centromeres and telomeres before training.
    :param test_intervals: Specified interval(s) will be held out for testing and evaluation only. (default to "chr20")
    :return: `ht` annotated with training information and the RF model
    """
    features = FEATURES
    test_intervals = test_intervals

    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills
    if not no_transmitted_singletons:
        tp_expr = tp_expr | ht.transmitted_singleton_adj 
    if not no_sibling_singletons:
        tp_expr = tp_expr | ht.sibling_singleton_adj # TODO: make sure this is acceptable Hail logic

    if test_intervals:
        if isinstance(test_intervals, str):
            test_intervals = [test_intervals]
        test_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in test_intervals
        ]

    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    if filter_centromere_telomere:
        logger.info("Filtering centromeres and telomeres from HT...")
        rf_ht = ht.filter(~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus]))
    else:
        rf_ht = ht

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

    logger.info("Joining original RF Table with training information")
    ht = ht.join(rf_ht, how="left")

    return ht, rf_model


def main(args):  # noqa: D103
    hl.init(log="/variant_qc_random_forest.log")

    if args.list_rf_runs:
        logger.info(f"RF runs:")
        pretty_print_runs(get_rf_runs(rf_run_path()))

    if args.train_rf:
        ht = get_variant_qc_annotations(test=args.test_dataset).ht()
        model_id = f"rf_{str(uuid.uuid4())[:8]}"
        rf_runs = get_rf_runs(rf_run_path())
        while model_id in rf_runs:
            model_id = f"rf_{str(uuid.uuid4())[:8]}"

        ht, rf_model = train_rf(ht,
            fp_to_tp=args.fp_to_tp,
            num_trees=args.num_trees,
            max_depth=args.max_depth,
            no_transmitted_singletons=args.no_transmitted_singletons,
            no_sibling_singletons=args.no_sibling_singletons,
            filter_centromere_telomere=args.filter_centromere_telomere,
            test_intervals=args.test_intervals)

        ht = ht.checkpoint(
            get_rf_training(model_id=model_id).path,
            overwrite=args.overwrite,
        )

        logger.info("Adding run to RF run list")
        rf_runs[model_id] = get_run_data(
            input_args={
                "transmitted_singletons": (
                    not args.no_transmitted_singletons
                ),
                # TODO: this will also require editing get_run_data maybe to accept No Sibling Singletons? 
                "sibling_singletons": (
                    not args.no_sibling_singletons
                ),
                "adj": args.adj,
                "filter_centromere_telomere": args.filter_centromere_telomere,
            },
            test_intervals=args.test_intervals,
            features_importance=hl.eval(ht.features_importance),
            test_results=hl.eval(ht.test_results),
        )

        with hl.hadoop_open(rf_run_path(), "w") as f:
            json.dump(rf_runs, f)

        logger.info("Saving RF model")
        save_model(
            rf_model,
            get_rf_model_path(model_id=model_id),
            overwrite=args.overwrite,
        )

    else:
        model_id = args.model_id

    if args.apply_rf:
        logger.info(f"Applying RF model {model_id}...")
        rf_model = load_model(get_rf_model_path(model_id=model_id))
        ht = get_rf_training(model_id=model_id).ht()
        features = hl.eval(ht.features)
        ht = apply_rf_model(ht, rf_model, features, label=LABEL_COL)

        logger.info("Finished applying RF model")
        ht = ht.annotate_globals(rf_model_id=model_id)
        ht = ht.checkpoint(
            get_rf_result(model_id=model_id).path,
            overwrite=args.overwrite,
        )

        ht_summary = ht.group_by(
            "tp", "fp", TRAIN_COL, LABEL_COL, PREDICTION_COL
        ).aggregate(n=hl.agg.count())
        ht_summary.show(n=20)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset. (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--test-dataset",
        help="If test of Variant QC table should be used. (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--model_id",
        help=(
            "Model ID. Created by --train_rf and only needed for --apply_rf without"
            " running --train_rf."
        ),
        required=False,
    )

    actions = parser.add_argument_group("Actions")
    actions.add_argument(
        "--list_rf_runs",
        help=(
            "Lists all previous RF runs, along with their model ID, parameters and"
            " testing results."
        ),
        action="store_true",
    )
    actions.add_argument(
        "--annotate_for_rf",
        help="Creates an annotated HT with features for RF.",
        action="store_true",
    )
    actions.add_argument("--train_rf", help="Trains RF model.", action="store_true")
    actions.add_argument(
        "--apply_rf", help="Applies RF model to the data.", action="store_true"
    )

    annotate_params = parser.add_argument_group("Annotate features parameters")
    annotate_params.add_argument(
        "--impute_features",
        help="If set, feature imputation is performed.",
        action="store_true",
    )
    annotate_params.add_argument(
        "--n_partitions",
        help="Desired number of partitions for annotated RF Table.",
        type=int,
        default=5000,
    )

    rf_params = parser.add_argument_group("Random Forest parameters")
    rf_params.add_argument(
        "--fp_to_tp",
        help=(
            "Ratio of FPs to TPs for training the RF model. If 0, all training examples"
            " are used. (default=1.0)"
        ),
        default=1.0,
        type=float,
    )
    rf_params.add_argument(
        "--test_intervals",
        help=(
            "The specified interval(s) will be held out for testing and evaluation"
            ' only. (default to "chr20")'
        ),
        nargs="+",
        type=str,
        default="chr20",
    )
    rf_params.add_argument(
        "--num_trees",
        help="Number of trees in the RF model. (default=500)",
        default=500,
        type=int,
    )
    rf_params.add_argument(
        "--max_depth",
        help="Maxmimum tree depth in the RF model. (default=5)",
        default=5,
        type=int,
    )

    training_params = parser.add_argument_group("Training data parameters")
    training_params.add_argument(
        "--adj", help="Use adj genotypes.", action="store_true"
    )
    training_params.add_argument(
        "--vqsr_model_id",
        help=(
            "If a VQSR model ID is provided the VQSR training annotations will be used"
            " for training."
        ),
        default="vqsr_alleleSpecificTrans",
        choices=["vqsr_classic", "vqsr_alleleSpecific", "vqsr_alleleSpecificTrans"],
        type=str,
    )
    training_params.add_argument(
        "--no_transmitted_singletons",
        help="Do not use transmitted singletons for training.",
        action="store_true",
    )
    training_params.add_argument(
        "--no_sibling_singletons",
        help="Do not use sibling singletons for training.",
        action="store_true",
    )
    training_params.add_argument(
        "--filter_centromere_telomere",
        help="Train RF without centromeres and telomeres.",
        action="store_true",
    )

    args = parser.parse_args()

    if not args.model_id and not args.train_rf and args.apply_rf:
        sys.exit(
            "Error: --model_id is required when running --apply_rf without running"
            " --train_rf too."
        )

    if args.model_id and args.train_rf:
        sys.exit(
            "Error: --model_id and --train_rf are mutually exclusive. --train_rf will"
            " generate a run model ID."
        )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
