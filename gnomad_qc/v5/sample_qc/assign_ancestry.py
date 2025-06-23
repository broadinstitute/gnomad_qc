"""Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels or spike ins."""

import argparse
import json
import logging
import pickle
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.sample_qc import hgdp_tgp_pop_outliers
from gnomad_qc.v4.sample_qc.assign_ancestry import (
    V3_SPIKE_PROJECTS,
    V4_POP_SPIKE_DICT,
    write_pca_results,
)
from gnomad_qc.v4.resources.sample_qc import (
    joint_qc_meta as v4_joint_qc_meta,
    related_samples_to_drop,  # TODO: remove when switch to v5.
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET

from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_checkpoint_path,
    get_logging_path,
)
from gnomad_qc.v5.resources.meta import project_meta, sample_id_collisions

from gnomad_qc.v5.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_joint_qc,
    get_pop_ht,
    get_pop_pr_ht,
    per_pop_min_rf_probs_json_path,
    pop_rf_path,
    # related_samples_to_drop, #TODO: switch from v4 related_samples_to_drop to v5 related_samples_to_drop once ready.
)

# TODO: Switch from using pop to gen_anc?

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("ancestry_assignment")
logger.setLevel(logging.INFO)


def run_pca(
    related_samples_to_drop: hl.Table,
    meta_ht: hl.Table,
    include_unreleasable_samples: hl.bool = False,
    n_pcs: int = 30,
    test: hl.bool = False,
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Run population PCA using `run_pca_with_relateds`.

    :param related_samples_to_drop: Table of related samples to drop from PCA run.
    :param meta_ht: Table of meta data containing 'releasable' field, used to drop samples that are not releasable unless include_unreleasable_samples is True.
    :param include_unreleasable_samples: Should unreleasable samples be included in the
        PCA.
    :param n_pcs: Number of PCs to compute.
    :param test: Subset QC MT to small test dataset.
    :return: Eigenvalues, scores and loadings from PCA.
    """
    logger.info("Running population PCA")
    qc_mt = get_joint_qc(test=test).mt()

    if args.test_on_chr20:
        logger.info("Filtering QC MT to chromosome 20...")
        qc_mt = hl.filter_intervals(
            qc_mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
        )

    samples_to_drop = related_samples_to_drop.select()

    if not include_unreleasable_samples:
        logger.info("Excluding unreleasable samples for PCA.")
        samples_to_drop = samples_to_drop.union(
            qc_mt.filter_cols(~meta_ht[qc_mt.col_key].releasable).cols().select()
        )
    else:
        logger.info("Including unreleasable samples for PCA.")

    return run_pca_with_relateds(qc_mt, samples_to_drop, n_pcs=n_pcs)


def main(args):
    """Assign global ancestry labels to samples."""
    hl.init(
        spark_conf={"spark.memory.offHeap.enabled": "false"},
        log="/home/jupyter/workspaces/gnomadproduction/assign_ancestry.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")
    try:

        include_unreleasable_samples = args.include_unreleasable_samples
        overwrite = args.overwrite
        test = args.test

        if test and args.test_on_chr20:
            raise ValueError("Both test and test_on_chr20 cannot be set to True.")

        # Use tmp path if either test dataset or test on chr20 is specified.
        use_tmp_path = args.test_on_chr20 | test
        include_v2_known_in_training = args.include_v2_known_in_training

        if args.run_pca:
            pop_eigenvalues, pop_scores_ht, pop_loadings_ht = run_pca(
                meta_ht=project_meta.ht(),
                related_samples_to_drop=related_samples_to_drop(release=False).ht(),
                include_unreleasable_samples=include_unreleasable_samples,
                n_pcs=args.n_pcs,
                test=test,
            )

            write_pca_results(
                pop_eigenvalues,
                pop_scores_ht,
                pop_loadings_ht,
                overwrite,
                include_unreleasable_samples,
                use_tmp_path,
            )
        if args.assign_pops:
            # Rename sample collision in v4 joint qc meta.
            v4_joint_qc_meta = v4_joint_qc_meta.ht()
            v4_joint_qc_meta = add_project_prefix_to_sample_collisions(
                t=v4_joint_qc_meta,
                sample_collisions=sample_id_collisions.ht(),
                project="gnomad",
            )
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("genetic_ancestry_assignment", environment="rwb"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Run script on test dataset.", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="Overwrite output files.", action="store_true"
    )
    parser.add_argument("--run-pca", help="Compute ancestry PCA", action="store_true")
    parser.add_argument(
        "--n-pcs",
        help="Number of PCs to compute for ancestry PCA. Defaults to 30.",
        default=30,
        type=int,
    )
    parser.add_argument(
        "--include-unreleasable-samples",
        help="Include unreleasable samples when computing PCA.",
        action="store_true",
    )
    parser.add_argument(
        "--test-on-chr20",
        help="Filter the QC Matrix Table to chromosome 20 before running PCA.",
        action="store_true",
    )
    parser.add_argument(
        "--assign-pops", help="Assigns pops from PCA.", action="store_true"
    )
    parser.add_argument(
        "--pop-pcs",
        help=(
            "List of PCs to use for ancestry assignment. The values provided should be"
            " 1-based. If a single integer is passed, the script assumes this"
            " represents the total PCs to use e.g. --pop-pcs=6 will use PCs"
            " 1,2,3,4,5,and 6. Defaults to 20 PCs."
        ),
        default=[20],
        type=int,
        nargs="+",
    )
    parser.add_argument(
        "--min-pop-prob",
        help="Minimum RF prob for pop assignment. Defaults to 0.75.",
        default=0.75,
        type=float,
    )
    parser.add_argument(
        "--include-v2-known-in-training",
        help=(
            "Whether to train RF classifier using v2 known pop labels. Default is"
            " False."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--v4-population-spike",
        help="List of v4 populations to spike into the RF training populations.",
        type=str,
        nargs="+",
        choices=["Arab", "Bedouin", "Persian", "Qatari"],
    )
    parser.add_argument(
        "--v3-population-spike",
        help="List of v3 populations to spike into the RF training populations.",
        type=str,
        nargs="+",
        choices=["asj", "ami", "afr", "amr", "eas", "sas", "fin", "nfe"],
    )
    parser.add_argument(
        "--compute-precision-recall",
        help=(
            "Compute precision and recall for the RF model using evaluation samples. "
            "This is computed for all evaluation samples as well as per population."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--number-pr-points",
        help=(
            "Number of min prob cutoffs to compute PR metrics for. e.g. 100 will "
            "compute PR metrics for min prob of 0 to 1 in increments of 0.01. Default "
            "is 100."
        ),
        default=100,
        type=int,
    )
    parser.add_argument(
        "--apply-per-pop-min-rf-probs",
        help=(
            "Apply per ancestry group minimum RF probabilities for finalized pop "
            "assignment instead of using '--min-pop-prob' for all samples. There must "
            "be a JSON file located in the path defined by the "
            "'per_pop_min_rf_probs_json_path' resource, or "
            "'--infer-per-pop-min-rf-probs' must be used."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--infer-per-pop-min-rf-probs",
        help=(
            "Whether to infer per ancestry group minimum RF probabilities and write "
            "them out to 'per_pop_min_rf_probs_json_path' before determining the "
            "finalized pop assignment."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-recall",
        help=(
            "Minimum recall value to choose per ancestry group minimum RF "
            "probabilities. This cutoff is applied first, and if the chosen cutoff "
            "results in a precision lower than '--min-precision', the minimum RF "
            "probabilities with the highest recall that meets '--min-precision' is "
            "used. Default is 0.99."
        ),
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--min-precision",
        help=(
            "Minimum precision value to choose per ancestry group minimum RF "
            "probabilities. This cutoff is applied if the chosen minimum RF "
            "probabilities cutoff using '--min-recall' results in a precision lower "
            "than this value. The minimum RF probabilities with the highest recall "
            "that meets '--min-precision' is used. Default is 0.99."
        ),
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--set-ami-exomes-to-remaining",
        help=(
            "Whether to change the ancestry group for any exomes inferred as 'ami' to"
            " 'remaining'. Should be used in cases where only a few exomes were"
            " inferred as amish to avoid having ancestry groups with only a few"
            " samples."
        ),
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
