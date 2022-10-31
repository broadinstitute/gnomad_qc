# noqa: D100

import argparse
import logging

import hail as hl
import pandas as pd
from gnomad.utils.file_utils import check_file_exists_raise_error
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v3.resources.basics import get_checkpoint_path, get_logging_path
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import subpop_outliers
from gnomad_qc.v3.sample_qc.sample_qc import ancestry_pca_scores, assign_pops
from gnomad_qc.v3.sample_qc.subpop_analysis import CURATED_SUBPOPS

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subpop_analysis")
logger.setLevel(logging.INFO)


def main(args):  # noqa: D103
    pop = args.pop
    include_unreleasable_samples = args.include_unreleasable_samples
    high_quality = args.high_quality

    try:
        check_file_exists_raise_error(
            ancestry_pca_scores(include_unreleasable_samples, high_quality, pop).path,
            error_if_not_exists=True,
            error_if_not_exists_msg=(
                "PCs have not yet been computed for the supplied parameters. Please run"
                " subpop_analysis.py with the desired parameters for"
                " `include_unreleasable_samples`, `high_quality`, and `pop`."
            ),
        )

        # Read in metadata
        meta_ht = meta.ht()

        # Initiate lists for storing results
        merged = []

        if args.remove_outliers:
            check_file_exists_raise_error(
                subpop_outliers(pop).path,
                error_if_not_exists=True,
                error_if_not_exists_msg=(
                    "The --remove-outliers option was used, but a Table of outlier"
                    f" samples does not exist for population {pop} at"
                    f" {subpop_outliers(pop).path}"
                ),
            )
            outliers_ht = subpop_outliers(pop).ht()
        else:
            outliers_ht = None

        # Run assign pops for each test combination of values in min_prob_list and
        # max_prop_mislabeled_list
        for min_prob in args.min_prob_list:
            for (
                max_proportion_mislabeled_training_samples
            ) in args.max_prop_mislabeled_list:
                logger.info(
                    "Running RF for min_prob: %s\nand max_proportion_mislabeled: %s ",
                    min_prob,
                    max_proportion_mislabeled_training_samples,
                )
                joint_pca_ht, joint_pca_fit = assign_pops(
                    min_prob=min_prob,
                    max_proportion_mislabeled_training_samples=max_proportion_mislabeled_training_samples,
                    include_unreleasable_samples=False,
                    pcs=args.pcs,
                    withhold_prop=args.withhold_prop,
                    pop=pop,
                    curated_subpops=CURATED_SUBPOPS[pop],
                    additional_samples_to_drop=outliers_ht,
                    high_quality=high_quality,
                    missing_label="Other",
                )

                # Filter to just evaluation samples
                ht = joint_pca_ht.filter(joint_pca_ht.evaluation_sample)

                # Calculate number of mislabeled samples
                # NOTE: `training_pop` is the known label for the sample, but could have been used for either training or evaluation # noqa
                total = ht.count()
                mislabelled = ht.aggregate(
                    hl.agg.count_where(ht.training_pop != ht.subpop)
                )
                correct = 1 - mislabelled / total

                # Calculate TP/FP/FN
                # Definitions:
                # TP = Samples where `training_pop` equals assigned `subpop`
                # FP = Samples where `training_pop` does not equal assigned `subpop` and `subpop` is not `Other`
                # FN = Samples where `training_pop` does not equal assigned `subpop` and
                # `subpop` is `Other`

                ht = ht.annotate(
                    result_status=hl.case()
                    .when(ht.training_pop == ht.subpop, "TP")
                    .when(
                        (ht.training_pop != ht.subpop) & (ht.subpop != "Other"),
                        "FP",
                    )
                    .when(
                        (ht.training_pop != ht.subpop) & (ht.subpop == "Other"),
                        "FN",
                    )
                    .default(hl.missing(hl.tstr))
                )

                tp = ht.aggregate(hl.agg.count_where(ht.result_status == "TP"))
                fp = ht.aggregate(hl.agg.count_where(ht.result_status == "FP"))
                fn = ht.aggregate(hl.agg.count_where(ht.result_status == "FN"))

                merged.append(
                    {
                        "min_prob": min_prob,
                        "max_proportion_mislabeled_training_samples": max_proportion_mislabeled_training_samples,
                        "percent_correct": correct * 100,
                        "TP": tp,
                        "FP": fp,
                        "FN": fn,
                        "error_rate": hl.eval(ht.assign_pops_from_pc_params.error_rate),
                    }
                )

        # Convert results to Table
        results_ht = hl.Table.parallelize(
            hl.literal(
                merged,
                "array<struct{min_prob: float,"
                " max_proportion_mislabeled_training_samples: float, percent_correct:"
                " float, TP: int, FP: int, FN: int, error_rate: float}>",
            )
        )

        # Calculate precision and recall
        # Precision: TP/(TP+FP)
        # Recall: TP/(TP+FN)
        results_ht = results_ht.annotate(
            precision=results_ht.TP / (results_ht.TP + results_ht.FP),
            recall=results_ht.TP / (results_ht.TP + results_ht.FN),
        )

        # Write results
        file_path = get_checkpoint_path(
            "subpop_analysis/parameter_testing_{}{}{}".format(
                pop,
                "_with_unreleasable_samples" if include_unreleasable_samples else "",
                ".high_quality" if high_quality else "",
            )
        )

        results_ht = results_ht.checkpoint(file_path, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("subpop_parameter_testing"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "This script generates a Table with precision and recall results to test"
            " the performance of subpop assignment"
        )
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to",
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--pop",
        help=(
            "Population to test (need to first run subpop_analysis.py with this `pop`)"
        ),
        type=str,
    )

    parser.add_argument(
        "--min-prob-list",
        help="List of minimum RF probabilities to test for subpop assignment",
        type=float,
        nargs="+",
        default=[
            0.10,
            0.20,
            0.30,
            0.40,
            0.50,
            0.60,
            0.70,
            0.75,
            0.80,
            0.85,
            0.90,
            0.95,
        ],
    )
    parser.add_argument(
        "--max-prop-mislabeled-list",
        help=(
            "List of maximum number of training samples that can be mislabelled to test"
            " for subpop assignment"
        ),
        type=float,
        nargs="+",
        default=[1],
    )
    parser.add_argument(
        "--pcs",
        help="List of PCs to use in the subpop RF assignment",
        type=int,
        nargs="+",
        default=list(range(1, 20)),
    )
    parser.add_argument(
        "--withhold-prop",
        help=(
            "Proportion of training pop samples to withhold from training, all samples"
            " will be kept if this flag is not used. NOTE: 20 percent of the training"
            " samples are used for evaluation, withholding samples will put aside a"
            " certain proportion of samples before splitting the remainder into"
            " training and evaluation. Setting a value for `withhold-prop` is"
            " recommended if values are supplied to `max-prop-mislabeled-list` in order"
            " to produce a reasonable PR curve when using iterations."
        ),
        type=float,
        default=None,
    )
    parser.add_argument(
        "--include-unreleasable-samples",
        help="Includes unreleasable samples for computing PCA",
        action="store_true",
    )
    parser.add_argument(
        "--high-quality",
        help="Filter to only high-quality samples when computing the PCA",
        action="store_true",
    )
    parser.add_argument(
        "--remove-outliers",
        help=(
            "Whether to remove outliers when training the random forest model. Outliers"
            " are manually determined after visualizing the PC plots"
        ),
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        from gnomad_qc.slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
