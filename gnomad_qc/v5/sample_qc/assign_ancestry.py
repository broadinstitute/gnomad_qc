"""Script to assign global ancestry labels to samples using HGDP or TGP labels and spike ins from known v3 and v4 genetic ancestry labels."""

import argparse
import json
import logging
import pickle
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from hail.utils.misc import new_temp_file

from gnomad_qc.v3.resources.sample_qc import hgdp_tgp_pop_outliers
from gnomad_qc.v4.resources.sample_qc import joint_qc_meta as v4_joint_qc_meta
from gnomad_qc.v4.resources.sample_qc import (
    related_samples_to_drop,  # TODO: remove when switch to v5.
)
from gnomad_qc.v4.sample_qc.assign_ancestry import V3_SPIKE_PROJECTS, V4_POP_SPIKE_DICT
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_checkpoint_path,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import project_meta, sample_id_collisions
from gnomad_qc.v5.resources.sample_qc import (  # related_samples_to_drop, #TODO: switch from v4 related_samples_to_drop to v5 related_samples_to_drop once ready.
    gen_anc_rf_path,
    genetic_ancestry_pca_eigenvalues,
    genetic_ancestry_pca_loadings,
    genetic_ancestry_pca_scores,
    get_gen_anc_ht,
    get_gen_anc_pr_ht,
    get_joint_qc,
    per_gen_anc_min_rf_probs_json_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("genetic_ancestry_assignment")
logger.setLevel(logging.INFO)


def run_pca(
    qc_mt: hl.MatrixTable,
    meta_ht: hl.Table,
    related_samples_to_drop: hl.Table,
    include_unreleasable_samples: hl.bool = False,
    n_pcs: int = 30,
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Run genetic ancestry PCA using `run_pca_with_relateds`.

    :param qc_mt: QC Matrix Table to run PCA on.
    :param related_samples_to_drop: Table of related samples to drop from PCA run.
    :param meta_ht: Table of meta data containing 'releasable' field, used to drop samples that are not releasable unless `include_unreleasable_samples` is True.
    :param include_unreleasable_samples: Should unreleasable samples be included in the
        PCA.
    :param n_pcs: Number of PCs to compute.
    :return: Eigenvalues, scores and loadings from PCA.
    """
    logger.info("Running population PCA")
    samples_to_drop = related_samples_to_drop.select()

    if not include_unreleasable_samples:
        logger.info("Excluding unreleasable samples for PCA.")
        samples_to_drop = samples_to_drop.union(
            qc_mt.filter_cols(~meta_ht[qc_mt.col_key].releasable).cols().select()
        )
    else:
        logger.info("Including unreleasable samples for PCA.")

    return run_pca_with_relateds(qc_mt, samples_to_drop, n_pcs=n_pcs)


def write_pca_results(
    gen_anc_pca_eigenvalues: List[float],
    gen_anc_pca_scores_ht: hl.Table,
    gen_anc_pca_loadings_ht: hl.Table,
    overwrite: hl.bool = False,
    included_unreleasables: hl.bool = False,
    test: hl.bool = False,
):
    """
    Write out the eigenvalue hail Table, scores hail Table, and loadings hail Table returned by `run_pca()`.

    :param gen_anc_pca_eigenvalues: List of eigenvalues returned by run_pca.
    :param gen_anc_pca_scores_ht: Table of scores returned by run_pca.
    :param gen_anc_pca_loadings_ht: Table of loadings returned by run_pca.
    :param overwrite: Whether to overwrite an existing file.
    :param included_unreleasables: Whether run_pca included unreleasable samples.
    :param test: Whether the test QC MT was used in the PCA.
    :return: None.
    """
    gen_anc_pca_eigenvalues_ht = hl.Table.parallelize(
        hl.literal(
            [
                {"PC": i + 1, "eigenvalue": x}
                for i, x in enumerate(gen_anc_pca_eigenvalues)
            ],
            "array<struct{PC: int, eigenvalue: float}>",
        )
    )
    gen_anc_pca_eigenvalues_ht.write(
        genetic_ancestry_pca_eigenvalues(included_unreleasables, test).path,
        overwrite=overwrite,
    )
    gen_anc_pca_scores_ht.write(
        genetic_ancestry_pca_scores(included_unreleasables, test).path,
        overwrite=overwrite,
    )
    gen_anc_pca_loadings_ht.write(
        genetic_ancestry_pca_loadings(included_unreleasables, test).path,
        overwrite=overwrite,
    )


def assign_pops(
    min_prob: float,
    include_unreleasable_samples: bool = False,
    pcs: List[int] = list(range(1, 21)),
    missing_label: str = "remaining",
    test: bool = False,
    overwrite: bool = False,
    include_v2_known_in_training: bool = False,
    v4_population_spike: Optional[List[str]] = None,
    v3_population_spike: Optional[List[str]] = None,
) -> Tuple[hl.Table, Any]:
    """
    Use a random forest model to assign global population labels based on the results from `run_pca`.

    Training data is the known label for HGDP and 1KG samples and all v2 samples with
    known pops unless specificied to restrict only to 1KG and HGDP samples. Can also
    specify a list of pops with known v3/v4 labels to include
    (v3_population_spike/v4_population_spike) for training. Pops supplied for v4 are
    specified by race/ethnicity and converted to a ancestry group using
    V4_POP_SPIKE_DICT. The method assigns a population label to all samples in the
    dataset.

    :param min_prob: Minimum RF probability for pop assignment.
    :param include_unreleasable_samples: Whether unreleasable samples were included in
        PCA.
    :param pcs: List of PCs to use in the RF.
    :param missing_label: Label for samples for which the assignment probability is
        smaller than `min_prob`.
    :param test: Whether running assigment on a test dataset.
    :param overwrite: Whether to overwrite existing files.
    :param include_v2_known_in_training: Whether to train RF classifier using v2 known
        pop labels. Default is False.
    :param v4_population_spike: Optional List of v4 populations to spike into the RF.
        Must be in v4_pop_spike dictionary. Defaults to None.
    :param v3_population_spike: Optional List of v3 populations to spike into the RF.
        Must be in v4_pop_spike dictionary. Defaults to None.
    :return: Table of pop assignments and the RF model.
    """
    logger.info("Prepping HT for RF...")
    pop_pca_scores_ht = prep_ht_for_rf(
        include_unreleasable_samples,
        test,
        include_v2_known_in_training,
        v4_population_spike,
        v3_population_spike,
    )
    pop_field = "pop"
    logger.info(
        "Running RF with PCs %s using %d training examples: %s",
        pcs,
        pop_pca_scores_ht.aggregate(
            hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
        ),
        pop_pca_scores_ht.aggregate(hl.agg.counter(pop_pca_scores_ht.training_pop)),
    )
    # Run the pop RF.
    pop_ht, pops_rf_model = assign_population_pcs(
        pop_pca_scores_ht,
        pc_cols=pcs,
        known_col="training_pop",
        output_col=pop_field,
        min_prob=min_prob,
        missing_label=missing_label,
    )

    pop_ht = pop_ht.annotate_globals(
        min_prob=min_prob,
        include_unreleasable_samples=include_unreleasable_samples,
        pcs=pcs,
        include_v2_known_in_training=include_v2_known_in_training,
    )

    if v3_population_spike:
        pop_ht = pop_ht.annotate_globals(
            v3_population_spike=v3_population_spike,
        )
    if v4_population_spike:
        pop_ht = pop_ht.annotate_globals(
            v4_population_spike=v4_population_spike,
        )

    return pop_ht, pops_rf_model


def main(args):
    """Assign genetic ancestry group labels to samples."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/assign_ancestry.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")
    try:

        include_unreleasable_samples = args.include_unreleasable_samples
        overwrite = args.overwrite
        test = args.test
        test_on_chr20 = args.test_on_chr20

        if test and test_on_chr20:
            raise ValueError("Both test and test_on_chr20 cannot be set to True.")

        # Use tmp path if either test dataset or test on chr20 is specified.
        use_tmp_path = test_on_chr20 or test
        include_v2_known_in_training = args.include_v2_known_in_training

        if args.run_pca:
            qc_mt = get_joint_qc(test=test).mt()

            if test_on_chr20:
                logger.info("Filtering QC MT to chromosome 20...")
                qc_mt = hl.filter_intervals(
                    qc_mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
                )

            gen_anc_eigenvalues, gen_anc_scores_ht, gen_anc_loadings_ht = run_pca(
                qc_mt=qc_mt,
                meta_ht=project_meta.ht(),
                related_samples_to_drop=related_samples_to_drop(release=False).ht(),
                include_unreleasable_samples=include_unreleasable_samples,
                n_pcs=args.n_pcs,
            )

            write_pca_results(
                gen_anc_eigenvalues,
                gen_anc_scores_ht,
                gen_anc_loadings_ht,
                overwrite,
                include_unreleasable_samples,
                use_tmp_path,
            )
        if args.assign_gen_anc:
            # Rename sample collision in v4 joint qc meta.
            v4_joint_qc_meta = v4_joint_qc_meta.ht()
            v4_joint_qc_meta = add_project_prefix_to_sample_collisions(
                t=v4_joint_qc_meta,
                sample_collisions=sample_id_collisions.ht(),
                project="gnomad",
            )

            pop_pcs = args.pop_pcs
            pop_pcs = list(range(1, pop_pcs[0] + 1)) if len(pop_pcs) == 1 else pop_pcs
            logger.info("Using following PCs: %s", pop_pcs)
            pop_ht, pops_rf_model = assign_pops(
                args.min_pop_prob,
                include_unreleasable_samples,
                pcs=pop_pcs,
                test=test,
                overwrite=overwrite,
                include_v2_known_in_training=include_v2_known_in_training,
                v4_population_spike=args.v4_population_spike,
                v3_population_spike=args.v3_population_spike,
            )

            logger.info("Writing pop ht...")
            pop_ht.write(get_pop_ht(test=test).path, overwrite=overwrite)

            with hl.hadoop_open(
                pop_rf_path(test=test),
                "wb",
            ) as out:
                pickle.dump(pops_rf_model, out)
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
    parser.add_argument(
        "--run-pca", help="Compute genetic ancestry PCA", action="store_true"
    )
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
        "--assign-gen-anc",
        help="Assigns genetic ancestry groups from PCA.",
        action="store_true",
    )
    parser.add_argument(
        "--gen-anc-pcs",
        help=(
            "List of PCs to use for genetic ancestry assignment. The values provided should be"
            " 1-based. If a single integer is passed, the script assumes this"
            " represents the total PCs to use e.g. --gen-anc-pcs=6 will use PCs"
            " 1,2,3,4,5,and 6. Defaults to 20 PCs."
        ),
        default=[20],
        type=int,
        nargs="+",
    )
    parser.add_argument(
        "--min-gen-anc-prob",
        help="Minimum RF prob for genetic ancestry assignment. Defaults to 0.75.",
        default=0.75,
        type=float,
    )
    parser.add_argument(
        "--include-v2-known-in-training",
        help=(
            "Whether to train RF classifier using v2 known genetic ancestry labels. Default is"
            " False."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--v4-spike",
        help="List of v4 exomes gen anc groups to spike into the RF training gen anc groups.",
        type=str,
        nargs="+",
        choices=["Arab", "Bedouin", "Persian", "Qatari"],
    )
    parser.add_argument(
        "--v3-spike",
        help="List of v3 genomes gen anc groups to spike into the RF training gen anc groups.",
        type=str,
        nargs="+",
        choices=["asj", "ami", "afr", "amr", "eas", "sas", "fin", "nfe"],
    )
    parser.add_argument(
        "--compute-precision-recall",
        help=(
            "Compute precision and recall for the RF model using evaluation samples. "
            "This is computed for all evaluation samples as well as per gen anc group."
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
        "--apply-per-grp-min-rf-probs",
        help=(
            "Apply per genetic ancestry group minimum RF probabilities for finalized group "
            "assignment instead of using '--min-grp-prob' for all samples. There must "
            "be a JSON file located in the path defined by the "
            "'per_grp_min_rf_probs_json_path' resource, or "
            "'--infer-per-grp-min-rf-probs' must be used."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--infer-per-grp-min-rf-probs",
        help=(
            "Whether to infer per genetic ancestry group minimum RF probabilities and write "
            "them out to 'per_grp_min_rf_probs_json_path' before determining the "
            "finalized group assignment."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-recall",
        help=(
            "Minimum recall value to choose per genetic ancestry group minimum RF "
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
            "Minimum precision value to choose per genetic ancestry group minimum RF "
            "probabilities. This cutoff is applied if the chosen minimum RF "
            "probabilities cutoff using '--min-recall' results in a precision lower "
            "than this value. The minimum RF probabilities with the highest recall "
            "that meets '--min-precision' is used. Default is 0.99."
        ),
        default=0.99,
        type=float,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
