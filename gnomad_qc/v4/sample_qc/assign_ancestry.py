"""Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels."""
import argparse
import json
import logging
import pickle
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.v3.resources.sample_qc import hgdp_tgp_pop_outliers
from gnomad_qc.v4.resources.basics import get_checkpoint_path
from gnomad_qc.v4.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_joint_qc,
    get_pop_ht,
    get_pop_pr_ht,
    joint_qc_meta,
    per_pop_min_rf_probs_json_path,
    pop_rf_path,
    related_samples_to_drop,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("ancestry_assignment")
logger.setLevel(logging.INFO)


V4_POP_SPIKE_DICT = {
    "Arab": "mid",
    "Bedouin": "mid",
    "Persian": "mid",
    "Qatari": "mid",
}
"""
Dictionary with potential pops to use for training (with v4 race/ethnicity as key and corresponding pop as value).
"""

V3_SPIKE_PROJECTS = {
    "asj": ["Jewish_Genome_Project"],
    "ami": ["NHLBI_WholeGenome_Sequencing"],
    "afr": ["TOPMED_Tishkoff_Cardiometabolics_Phase4"],
    "amr": [
        "PAGE: Global Reference Panel",
        "PAGE: Multiethnic Cohort (MEC)",
        "CostaRica",
    ],
    "eas": ["osaka"],
    "sas": ["CCDG_PROMIS", "TOPMED_Saleheen_PROMIS_Phase4"],
    "fin": [
        "G4L Initiative Stanley Center",
        "WGSPD3_Palotie_FinnishBP_THL_WGS",
        "WGSPD3_Palotie_Finnish_WGS",
        "WGSPD3_Palotie_Finnish_WGS_December2018",
    ],
    "nfe": [
        "Estonia_University of Tartu_Whole Genome Sequencing",
        "CCDG_Atrial_Fibrillation_Munich",
        "CCDG_Atrial_Fibrillation_Norway",
        "CCDG_Atrial_Fibrillation_Sweden",
        "Estonia_University of Tartu_Whole Genome Sequencing",
        "WGSPD",
    ],
}
"""
Dictionary with v3 pops as keys and approved cohorts to use for training for those pops as values. Decisions were made based on results of analysis to determine which v3 samples/cohorts to use as training samples (per sample calculating the mean distance to all samples in a given population, and an alternative calculating per sample the mean distance limited to only HGDP/1KG samples in each population).
Projects that were excluding based on these analyses are:
afr: NHLBI_WholeGenome_Sequencing
ami: Pedigree-Based Whole Genome Sequencing of Affective and Psychotic Disorders
amr: NHLBI_WholeGenome_Sequencing
amr: PAGE: Women''s Health Initiative (WHI)
"""


def run_pca(
    related_samples_to_drop: hl.Table,
    include_unreleasable_samples: hl.bool = False,
    n_pcs: int = 30,
    test: hl.bool = False,
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Run population PCA using `run_pca_with_relateds`.

    :param related_samples_to_drop: Table of related samples to drop from PCA run.
    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA.
    :param n_pcs: Number of PCs to compute.
    :param test: Subset QC MT to small test dataset.
    :return: Eigenvalues, scores and loadings from PCA.
    """
    logger.info("Running population PCA")
    qc_mt = get_joint_qc(test=test).mt()
    joint_meta = joint_qc_meta.ht()
    samples_to_drop = related_samples_to_drop.select()

    if not include_unreleasable_samples:
        logger.info("Excluding unreleasable samples for PCA.")
        samples_to_drop = samples_to_drop.union(
            qc_mt.filter_cols(~joint_meta[qc_mt.col_key].releasable).cols().select()
        )
    else:
        logger.info("Including unreleasable samples for PCA.")

    return run_pca_with_relateds(qc_mt, samples_to_drop, n_pcs=n_pcs)


def prep_ht_for_rf(
    include_unreleasable_samples: bool = False,
    test: bool = False,
    include_v2_known_in_training: bool = False,
    v4_population_spike: Optional[List[str]] = None,
    v3_population_spike: Optional[List[str]] = None,
) -> hl.Table:
    """
    Prepare the PCA scores hail Table for the random forest population assignment runs.

    Either train the RF with only HGDP and TGP, or HGDP and TGP and all v2 known labels. Can also specify list of pops with known v3/v4 labels to include (v3_population_spike/v4_population_spike) for training. Pops supplied for v4 are specified by race/ethnicity and converted to a ancestry group using V4_POP_SPIKE_DICT.

    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA.
    :param test: Whether RF should run on the test QC MT.
    :param include_v2_known_in_training: Whether to train RF classifier using v2 known pop labels. Default is False.
    :param v4_population_spike: Optional List of populations to spike into training. Must be in v4_population_spike dictionary. Default is None.
    :param v3_population_spike: Optional List of populations to spike into training. Must be in v3_population_spike dictionary. Default is None.
    :return Table with input for the random forest.
    """
    pop_pca_scores_ht = ancestry_pca_scores(include_unreleasable_samples, test).ht()
    joint_meta_ht = joint_qc_meta.ht()

    # Collect sample names of hgdp/tgp outliers to remove (these are outliers
    # found by Alicia Martin's group during pop-specific PCA analyses as well
    # as one duplicate sample)
    hgdp_tgp_outliers = hl.literal(hgdp_tgp_pop_outliers.ht().s.collect())

    joint_meta_ht = joint_meta_ht.annotate(
        hgdp_or_tgp=hl.or_else(
            joint_meta_ht.v3_meta.v3_subsets.hgdp
            | joint_meta_ht.v3_meta.v3_subsets.tgp,
            False,
        )
    )
    # Note: Excluding v3_project_pop="oth" samples from training. These are samples from
    #  Oceania and there are only a few known Oceania samples, and in past inference
    #  analyses no new samples are inferred as belonging to this group.
    training_pop = hl.or_missing(
        joint_meta_ht.hgdp_or_tgp
        & (joint_meta_ht.v3_meta.v3_project_pop != "oth")
        & ~hgdp_tgp_outliers.contains(joint_meta_ht.s),
        joint_meta_ht.v3_meta.v3_project_pop,
    )

    if include_v2_known_in_training:
        training_pop = hl.coalesce(
            training_pop,
            hl.or_missing(
                joint_meta_ht.data_type == "exomes", joint_meta_ht.v2_meta.v2_known_pop
            ),
        )

    if v4_population_spike:
        logger.info(
            "Spiking v4 pops, %s, into the RF training data...", v4_population_spike
        )

        v4_spike_err = set(v4_population_spike) - set(V4_POP_SPIKE_DICT.keys())
        if len(v4_spike_err) > 0:
            raise ValueError(f"Pops: {v4_spike_err}, are not in V4_POP_SPIKE_DICT")

        v4_population_spike = {r: V4_POP_SPIKE_DICT[r] for r in v4_population_spike}

        training_pop = hl.coalesce(
            training_pop,
            hl.literal(v4_population_spike).get(joint_meta_ht.v4_race_ethnicity),
        )

    if v3_population_spike:
        logger.info(
            "Spiking v3 pops, %s, into the RF training data", v3_population_spike
        )
        v3_spike_err = set(v3_population_spike) - set(V3_SPIKE_PROJECTS.keys())
        if len(v3_spike_err) > 0:
            raise ValueError(f"Pops: {v3_spike_err}, are not in V3_SPIKE_PROJECTS")

        # Filter to only pre-determined list of v3 cohorts for the v3 spike-ins
        training_pop = hl.coalesce(
            training_pop,
            hl.or_missing(
                hl.literal(V3_SPIKE_PROJECTS)
                .get(joint_meta_ht.v3_meta.v3_project_pop, hl.empty_array(hl.tstr))
                .contains(joint_meta_ht.v3_meta.v3_research_project),
                joint_meta_ht.v3_meta.v3_project_pop,
            ),
        )

    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        **joint_meta_ht.select(
            training_pop=training_pop,
            hgdp_or_tgp=joint_meta_ht.hgdp_or_tgp,
        )[pop_pca_scores_ht.key]
    )

    return pop_pca_scores_ht


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

    Training data is the known label for HGDP and 1KG samples and all v2 samples with known pops unless specificied to restrict only to 1KG and HGDP samples. Can also specify a list of pops with known v3/v4 labels to include (v3_population_spike/v4_population_spike) for training. Pops supplied for v4 are specified by race/ethnicity and converted to a ancestry group using V4_POP_SPIKE_DICT. The method assigns
    a population label to all samples in the dataset.

    :param min_prob: Minimum RF probability for pop assignment.
    :param include_unreleasable_samples: Whether unreleasable samples were included in PCA.
    :param pcs: List of PCs to use in the RF.
    :param missing_label: Label for samples for which the assignment probability is smaller than `min_prob`.
    :param test: Whether running assigment on a test dataset.
    :param overwrite: Whether to overwrite existing files.
    :param include_v2_known_in_training: Whether to train RF classifier using v2 known pop labels. Default is False.
    :param v4_population_spike: Optional List of v4 populations to spike into the RF. Must be in v4_pop_spike dictionary. Defaults to None.
    :param v3_population_spike: Optional List of v3 populations to spike into the RF. Must be in v4_pop_spike dictionary. Defaults to None.
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
    pop_ht = pop_ht.checkpoint(
        get_checkpoint_path(f"assign_pops_rf_iter_1_pop_lots_of_loggers_{pcs[-1]}"),
        overwrite=overwrite,
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


def write_pca_results(
    pop_pca_eigenvalues: List[float],
    pop_pca_scores_ht: hl.Table,
    pop_pca_loadings_ht: hl.Table,
    overwrite: hl.bool = False,
    included_unreleasables: hl.bool = False,
    test: hl.bool = False,
):
    """
    Write out the eigenvalue hail Table, scores hail Table, and loadings hail Table returned by run_pca().

    :param pop_pca_eigenvalues: List of eigenvalues returned by run_pca.
    :param pop_pca_scores_ht: Table of scores returned by run_pca.
    :param pop_pca_loadings_ht: Table of loadings returned by run_pca.
    :param overwrite: Whether to overwrite an existing file.
    :param included_unreleasables: Whether run_pca included unreleasable samples.
    :param test: Whether the test QC MT was used in the PCA.
    :return: None
    """
    pop_pca_eigenvalues_ht = hl.Table.parallelize(
        hl.literal(
            [{"PC": i + 1, "eigenvalue": x} for i, x in enumerate(pop_pca_eigenvalues)],
            "array<struct{PC: int, eigenvalue: float}>",
        )
    )
    pop_pca_eigenvalues_ht.write(
        ancestry_pca_eigenvalues(included_unreleasables, test).path,
        overwrite=overwrite,
    )
    pop_pca_scores_ht.write(
        ancestry_pca_scores(included_unreleasables, test).path,
        overwrite=overwrite,
    )
    pop_pca_loadings_ht.write(
        ancestry_pca_loadings(included_unreleasables, test).path,
        overwrite=overwrite,
    )


def get_most_likely_pop_expr(
    ht: hl.Table,
) -> Tuple[hl.expr.StructExpression, List[Tuple[str, str]]]:
    """
    Get StructExpression with 'pop' and 'prob' for the most likely population based on RF probabilities.

    :param ht: Input population inference Table with random forest probabilities.
    :return: Struct Expression with 'pop' and 'prob' for the highest RF probability.
    """
    # Get list of all population RF probability annotations.
    prob_ann = [(x.split("prob_")[-1], x) for x in ht.row if x.startswith("prob_")]

    # Sort a list of all population RF probabilities and keep the highest as the most
    # likely pop.
    most_likely_pop_expr = hl.rbind(
        hl.sorted(
            [hl.struct(pop=pop, prob=ht[ann]) for pop, ann in prob_ann],
            key=lambda el: el["prob"],
            reverse=True,
        ),
        lambda pop_probs: pop_probs[0],
    )

    return most_likely_pop_expr, prob_ann


def compute_precision_recall(ht: hl.Table, num_pr_points: int = 1000) -> hl.Table:
    """
    Create Table with false positives (FP), true positives (TP), false negatives (FN), precision, and recall.

    Includes population specific calculations.

    :param ht: Input population inference Table with random forest probabilities.
    :param num_pr_points: Number of min prob cutoffs to compute PR metrics for.
    :return: Table with FP, TP, FN, precision, and recall.
    """
    # Use only RF evaluation samples to compute metrics.
    ht = ht.filter(ht.evaluation_sample)
    most_likely_pop, prob_ann = get_most_likely_pop_expr(ht)
    ht = ht.annotate(most_likely_pop=most_likely_pop)

    # Add a list of min_prob_cutoffs from 0 to 1 in increments of 1/num_pr_points and
    # explode so each sample has one row for each min_prob_cutoff.
    ht = ht.annotate(
        min_prob_cutoff=hl.literal(list(range(num_pr_points))) / float(num_pr_points)
    )
    ht = ht.explode(ht.min_prob_cutoff)

    # For each 'pop', filter to any sample that has 'pop' as the known population or
    # most likely population based on RF probabilities.
    pops = [None] + [pop for pop, _ in prob_ann]

    pr_agg = hl.struct()
    for pop in pops:
        # Group by min_prob_cutoffs and aggregate to get TP, FP, and FN per
        # min_prob_cutoffs.
        meet_cutoff = ht.most_likely_pop.prob >= ht.min_prob_cutoff
        same_training_most_likely_pop = ht.training_pop == ht.most_likely_pop.pop
        if pop:
            is_most_likely_pop = ht.most_likely_pop.pop == pop
            is_training_pop = ht.training_pop == pop
            tp = hl.agg.count_where(meet_cutoff & is_most_likely_pop & is_training_pop)
            fp = hl.agg.count_where(meet_cutoff & is_most_likely_pop & ~is_training_pop)
            fn = hl.agg.count_where(is_training_pop & ~meet_cutoff)
        else:
            tp = hl.agg.count_where(meet_cutoff & same_training_most_likely_pop)
            fp = hl.agg.count_where(meet_cutoff & ~same_training_most_likely_pop)
            fn = hl.agg.count_where(~meet_cutoff)

        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        agg = hl.struct(TP=tp, FP=fp, FN=fn, precision=precision, recall=recall)
        if pop:
            agg = hl.struct(**{pop: agg})
        pr_agg = pr_agg.annotate(**agg)

    ht = ht.group_by("min_prob_cutoff").aggregate(**pr_agg)
    ht = ht.annotate_globals(pops=[pop for pop, _ in prob_ann])

    return ht


def infer_per_pop_min_rf_probs(
    ht: hl.Table, min_recall: float = 0.99, min_precision: float = 0.99
) -> Dict[str, Dict[str, float]]:
    """
    Infer per ancestry group minimum RF probabilities from precision and recall values.

    Minimum recall (`min_recall`) is used to choose per ancestry group minimum RF
    probabilities. This `min_recall` cutoff is applied first, and if the chosen minimum
    RF probabilities cutoff results in a precision lower than `min_precision`, the
    minimum RF probabilities with the highest recall that meets `min_precision` is
    used.

    :param ht: Precision recall Table returned by `compute_precision_recall`.
    :param min_recall: Minimum recall value to choose per ancestry group minimum RF
        probabilities. Default is 0.99.
    :param min_precision: Minimum precision value to choose per ancestry group minimum
        RF probabilities. Default is 0.99.
    :return: Dictionary of per pop min probability cutoffs `min_prob_cutoff`.
    """
    # Get list of all populations.
    pops = hl.eval(ht.index_globals().pops)
    min_prob_per_pop = {}
    for pop in pops:
        pr_vals = hl.tuple(
            [ht.min_prob_cutoff, ht[pop].precision, ht[pop].recall]
        ).collect()
        min_probs, precisions, recalls = map(list, zip(*pr_vals))

        try:
            idx = next(
                x for x, r in reversed(list(enumerate(recalls))) if r >= min_recall
            )
        except StopIteration:
            raise ValueError(
                f"Recall never reaches {min_recall} for {pop}. Recall values are: "
                f"{recalls}"
            )
        precision_cutoff = precisions[idx]
        if precision_cutoff < min_precision:
            try:
                idx = next(x for x, p in enumerate(precisions) if p >= min_precision)
            except StopIteration:
                raise ValueError(
                    f"Precision never reaches {min_precision} for {pop}. Precision "
                    f"values are: {precisions}"
                )
        min_prob_per_pop[pop] = {
            "precision_cutoff": precisions[idx],
            "recall_cutoff": recalls[idx],
            "min_prob_cutoff": min_probs[idx],
        }

    return min_prob_per_pop


def assign_pop_with_per_pop_probs(
    pop_ht: hl.Table,
    min_prob_cutoffs: Dict[str, float],
    missing_label: str = "remaining",
) -> hl.Table:
    """
    Assign samples to populations based on population-specific minimum RF probabilities.

    :param pop_ht: Table containing results of population inference.
    :param min_prob_cutoffs: Dictionary with population as key, and minimum RF probability required to assign a sample to that population as value.
    :param missing_label: Label for samples for which the assignment probability is smaller than required minimum probability.
    :return: Table with 'pop' annotation based on supplied per pop min probabilities.
    """
    min_prob_cutoffs = hl.literal(min_prob_cutoffs)
    pop_prob, _ = get_most_likely_pop_expr(pop_ht)
    pop = pop_prob.pop
    pop_ht = pop_ht.annotate(
        pop=hl.if_else(pop_prob.prob >= min_prob_cutoffs.get(pop), pop, missing_label)
    )
    pop_ht = pop_ht.annotate_globals(min_prob_cutoffs=min_prob_cutoffs)

    return pop_ht


def main(args):
    """Assign global ancestry labels to samples."""
    hl.init(
        log="/assign_ancestry.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    try:
        include_unreleasable_samples = args.include_unreleasable_samples
        overwrite = args.overwrite
        test = args.test
        include_v2_known_in_training = args.include_v2_known_in_training

        if args.run_pca:
            pop_eigenvalues, pop_scores_ht, pop_loadings_ht = run_pca(
                related_samples_to_drop(release=False).ht(),
                include_unreleasable_samples,
                args.n_pcs,
                test,
            )

            write_pca_results(
                pop_eigenvalues,
                pop_scores_ht,
                pop_loadings_ht,
                overwrite,
                include_unreleasable_samples,
                test,
            )

        if args.assign_pops:
            pop_pcs = args.pop_pcs
            pop_pcs = (
                list(range(1, pop_pcs + 1)) if isinstance(pop_pcs, int) else pop_pcs
            )
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

        if args.compute_precision_recall:
            ht = compute_precision_recall(
                get_pop_ht(test=test).ht(), num_pr_points=args.number_pr_points
            )
            ht.write(get_pop_pr_ht(test=test).path, overwrite=overwrite)

        if args.apply_per_pop_min_rf_probs:
            pop = get_pop_ht(test=test)
            if args.infer_per_pop_min_rf_probs:
                min_probs = infer_per_pop_min_rf_probs(
                    get_pop_pr_ht(test=test).ht(),
                    min_recall=args.min_recall,
                    min_precision=args.min_precision,
                )
                logger.info(
                    "Per ancestry group min prob cutoff inference: %s", min_probs
                )
                min_probs = {
                    pop: cutoff["min_prob_cutoff"] for pop, cutoff in min_probs.items()
                }
                with hl.hadoop_open(per_pop_min_rf_probs_json_path(), "w") as d:
                    d.write(json.dumps(min_probs))

            with hl.hadoop_open(per_pop_min_rf_probs_json_path(), "r") as d:
                min_probs = json.load(d)
            logger.info(
                "Using the following min prob cutoff per ancestry group: %s", min_probs
            )
            pop_ht = pop.ht().checkpoint(new_temp_file("pop_ht", extension="ht"))
            pop_ht = assign_pop_with_per_pop_probs(pop_ht, min_probs)
            pop_ht.write(pop.path, overwrite=overwrite)

    finally:
        hl.copy_log(
            f"gs://gnomad-tmp-4day/ancestry_assignment/ancestry_assignment.rf.log"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )
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
        default=20,
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
            "Number of min prob cutoffs to compute PR metrics for. e.g. 1000 will "
            "compute PR metrics for min prob of 0 to 1 in increments of 0.001. Default "
            "is 1000."
        ),
        default=1000,
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

    args = parser.parse_args()

    if args.slack_channel:
        from gnomad_qc.slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
