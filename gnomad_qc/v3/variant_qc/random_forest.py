import argparse
import json
import logging
import sys
from typing import List, Optional, Union
import uuid

import hail as hl

from gnomad.resources.grch38.reference_data import (
    get_truth_ht,
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
from gnomad_qc.v3.resources.basics import get_checkpoint_path
from gnomad_qc.v3.resources.annotations import (
    allele_data,
    fam_stats,
    get_freq,
    get_info,
    get_vqsr_filters,
    qc_ac,
)
from gnomad_qc.v3.resources.variant_qc import (
    get_rf_annotations,
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
    "InbreedingCoeff",
    "n_alt_alleles",
    "variant_type",
]
INBREEDING_COEFF_HARD_CUTOFF = -0.3
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


def create_rf_ht(
    impute_features: bool = True,
    adj: bool = False,
    n_partitions: int = 5000,
    checkpoint_path: Optional[str] = None,
) -> hl.Table:
    """
    Creates a Table with all necessary annotations for the random forest model.

    Annotations that are included:

        Features for RF:
            - InbreedingCoeff
            - variant_type
            - allele_type
            - n_alt_alleles
            - has_star
            - AS_QD
            - AS_pab_max
            - AS_MQRankSum
            - AS_SOR
            - AS_ReadPosRankSum

        Training sites (bool):
            - transmitted_singleton
            - fail_hard_filters - (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30)

    :param bool impute_features: Whether to impute features using feature medians (this is done by variant type)
    :param str adj: Whether to use adj genotypes
    :param int n_partitions: Number of partitions to use for final annotated table
    :param str checkpoint_path: Optional checkpoint path for the Table before median imputation and/or aggregate summary
    :return: Hail Table ready for RF
    :rtype: Table
    """

    group = "adj" if adj else "raw"

    ht = get_info(split=True).ht()
    ht = ht.transmute(**ht.info)
    ht = ht.select("lowqual", "AS_lowqual", "FS", "MQ", "QD", *INFO_FEATURES)

    inbreeding_ht = get_freq().ht()
    inbreeding_ht = inbreeding_ht.select(
        InbreedingCoeff=hl.if_else(
            hl.is_nan(inbreeding_ht.InbreedingCoeff),
            hl.null(hl.tfloat32),
            inbreeding_ht.InbreedingCoeff,
        )
    )
    trio_stats_ht = fam_stats.ht()
    trio_stats_ht = trio_stats_ht.select(
        f"n_transmitted_{group}", f"ac_children_{group}"
    )

    truth_data_ht = get_truth_ht()
    allele_data_ht = allele_data.ht()
    allele_counts_ht = qc_ac.ht()

    logger.info("Annotating Table with all columns from multiple annotation Tables")
    ht = ht.annotate(
        **inbreeding_ht[ht.key],
        **trio_stats_ht[ht.key],
        **truth_data_ht[ht.key],
        **allele_data_ht[ht.key].allele_data,
        **allele_counts_ht[ht.key],
    )
    # Filter to only variants found in high quality samples and are not lowqual
    ht = ht.filter((ht[f"ac_qc_samples_{group}"] > 0) & ~ht.AS_lowqual)
    ht = ht.select(
        "a_index",
        "was_split",
        *FEATURES,
        *TRUTH_DATA,
        **{
            "transmitted_singleton": (ht[f"n_transmitted_{group}"] == 1)
            & (ht[f"ac_qc_samples_{group}"] == 2),
            "fail_hard_filters": (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
        },
        singleton=ht.ac_release_samples_raw == 1,
        ac_raw=ht.ac_qc_samples_raw,
        ac=ht.ac_release_samples_adj,
        ac_qc_samples_unrelated_raw=ht.ac_qc_samples_unrelated_raw,
    )

    ht = ht.repartition(n_partitions, shuffle=False)
    if checkpoint_path:
        ht = ht.checkpoint(checkpoint_path, overwrite=True)

    if impute_features:
        ht = median_impute_features(ht, {"variant_type": ht.variant_type})

    summary = ht.group_by("omni", "mills", "transmitted_singleton",).aggregate(
        n=hl.agg.count()
    )
    logger.info("Summary of truth data annotations:")
    summary.show(20)

    return ht


def train_rf(
    ht: hl.Table,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    no_transmitted_singletons: bool = False,
    no_inbreeding_coeff: bool = False,
    vqsr_training: bool = False,
    vqsr_model_id: str = False,
    filter_centromere_telomere: bool = False,
    test_intervals: Union[str, List[str]] = "chr20",
):
    """
    Train random forest model using `train_rf_model`

    :param ht: Table containing annotations needed for RF training, built with `create_rf_ht`
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used.
    :param num_trees: Number of trees in the RF model.
    :param max_depth: Maxmimum tree depth in the RF model.
    :param no_transmitted_singletons: Do not use transmitted singletons for training.
    :param no_inbreeding_coeff: Do not use inbreeding coefficient as a feature for training.
    :param vqsr_training: Use VQSR training sites to train the RF.
    :param vqsr_model_id: VQSR model to use for vqsr_training. `vqsr_training` must be True for this parameter to be used.
    :param filter_centromere_telomere: Filter centromeres and telomeres before training.
    :param test_intervals: Specified interval(s) will be held out for testing and evaluation only. (default to "chr20")
    :return: `ht` annotated with training information and the RF model
    """
    features = FEATURES
    test_intervals = test_intervals
    if no_inbreeding_coeff:
        logger.info("Removing InbreedingCoeff from list of features...")
        features.remove("InbreedingCoeff")

    if vqsr_training:
        logger.info("Using VQSR training sites for RF training...")
        vqsr_ht = get_vqsr_filters(vqsr_model_id, split=True).ht()
        ht = ht.annotate(
            vqsr_POSITIVE_TRAIN_SITE=vqsr_ht[ht.key].info.POSITIVE_TRAIN_SITE,
            vqsr_NEGATIVE_TRAIN_SITE=vqsr_ht[ht.key].info.NEGATIVE_TRAIN_SITE,
        )
        tp_expr = ht.vqsr_POSITIVE_TRAIN_SITE
        fp_expr = ht.vqsr_NEGATIVE_TRAIN_SITE
    else:
        fp_expr = ht.fail_hard_filters
        tp_expr = ht.omni | ht.mills
        if not no_transmitted_singletons:
            tp_expr = tp_expr | ht.transmitted_singleton

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


def main(args):
    hl.init(log="/variant_qc_random_forest.log")

    if args.list_rf_runs:
        logger.info(f"RF runs:")
        pretty_print_runs(get_rf_runs(rf_run_path()))

    if args.annotate_for_rf:
        ht = create_rf_ht(
            impute_features=args.impute_features,
            adj=args.adj,
            n_partitions=args.n_partitions,
            checkpoint_path=get_checkpoint_path("rf_annotation"),
        )
        ht.write(
            get_rf_annotations(args.adj).path, overwrite=args.overwrite,
        )
        logger.info(f"Completed annotation wrangling for random forests model training")

    if args.train_rf:
        model_id = f"rf_{str(uuid.uuid4())[:8]}"
        rf_runs = get_rf_runs(rf_run_path())
        while model_id in rf_runs:
            model_id = f"rf_{str(uuid.uuid4())[:8]}"

        ht, rf_model = train_rf(
            get_rf_annotations(args.adj).ht(),
            fp_to_tp=args.fp_to_tp,
            num_trees=args.num_trees,
            max_depth=args.max_depth,
            no_transmitted_singletons=args.no_transmitted_singletons,
            no_inbreeding_coeff=args.no_inbreeding_coeff,
            vqsr_training=args.vqsr_training,
            vqsr_model_id=args.vqsr_model_id,
            filter_centromere_telomere=args.filter_centromere_telomere,
            test_intervals=args.test_intervals,
        )

        ht = ht.checkpoint(
            get_rf_training(model_id=model_id).path, overwrite=args.overwrite,
        )

        logger.info("Adding run to RF run list")
        rf_runs[model_id] = get_run_data(
            input_args={
                "transmitted_singletons": None
                if args.vqsr_training
                else not args.no_transmitted_singletons,
                "adj": args.adj,
                "vqsr_training": args.vqsr_training,
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
            rf_model, get_rf_model_path(model_id=model_id), overwrite=args.overwrite,
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
            get_rf_result(model_id=model_id).path, overwrite=args.overwrite,
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
        "--model_id",
        help="Model ID. Created by --train_rf and only needed for --apply_rf without running --train_rf.",
        required=False,
    )

    actions = parser.add_argument_group("Actions")
    actions.add_argument(
        "--list_rf_runs",
        help="Lists all previous RF runs, along with their model ID, parameters and testing results.",
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
        help="Ratio of FPs to TPs for training the RF model. If 0, all training examples are used. (default=1.0)",
        default=1.0,
        type=float,
    )
    rf_params.add_argument(
        "--test_intervals",
        help='The specified interval(s) will be held out for testing and evaluation only. (default to "chr20")',
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
        "--vqsr_training", help="Use VQSR training examples.", action="store_true"
    )
    training_params.add_argument(
        "--vqsr_model_id",
        help="If a VQSR model ID is provided the VQSR training annotations will be used for training.",
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
        "--no_inbreeding_coeff",
        help="Train RF without inbreeding coefficient as a feature.",
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
            "Error: --model_id is required when running --apply_rf without running --train_rf too."
        )

    if args.model_id and args.train_rf:
        sys.exit(
            "Error: --model_id and --train_rf are mutually exclusive. --train_rf will generate a run model ID."
        )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
