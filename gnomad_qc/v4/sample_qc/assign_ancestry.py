"""Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels."""
import argparse
import json
import logging
import pickle
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds

from gnomad.utils.slack import slack_notifications
from gnomad_qc.v3.resources.meta import meta as v3_meta
from gnomad_qc.v3.resources.sample_qc import hgdp_tgp_pop_outliers
from gnomad_qc.v4.resources.basics import get_checkpoint_path
from gnomad_qc.v4.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_joint_qc,
    get_pop_ht,
    joint_qc_meta,
    pca_related_samples_to_drop,
    pop_rf_path,
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
    seed: int = 24,
    test: bool = False,
    include_v2_known_in_training: bool = False,
    v4_population_spike: Optional[List[str]] = None,
    v3_population_spike: Optional[List[str]] = None,
) -> hl.Table:
    """
    Prepare the PCA scores hail Table for the random forest population assignment runs.


    Either train the RF with only HGDP and TGP, or HGDP and TGP and all v2 known labels. Can also specify list of pops with known v3/v4 labels to include (v3_population_spike/v4_population_spike) for training. Pops supplied for v4 are specified by race/ethnicity and converted to a ancestry group using V4_POP_SPIKE_DICT.


    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA.
    :param seed: Random seed, defaults to 24.
    :param test: Whether RF should run on the test QC MT.
    :param include_v2_known_in_training: Whether to train RF classifier using v2 known pop labels. Default is False.
    :param v4_population_spike: Optional List of populations to spike into training. Must be in v4_population_spike dictionary. Default is None.
    :param v3_population_spike: Optional List of populations to spike into training. Must be in v3_population_spike dictionary. Default is None.
    :return Table with input for the random forest.
    """
    pop_pca_scores_ht = ancestry_pca_scores(include_unreleasable_samples, test).ht()
    global joint_qc_meta
    global v3_meta

    joint_meta = joint_qc_meta.ht()[pop_pca_scores_ht.key]

    # Collect sample names of hgdp/tgp outliers to remove (these are outliers found by Alicia Martin's group during pop-specific PCA analyses as well as one duplicate sample)
    hgdp_tgp_outliers = hl.literal(hgdp_tgp_pop_outliers.ht().s.collect())

    training_pop = hl.or_missing(
        (joint_meta.v3_meta.v3_subsets.hgdp | joint_meta.v3_meta.v3_subsets.tgp)
        & (
            joint_meta.v3_meta.v3_project_pop != "oth"
        )  # Not using v3_project_pop="oth" samples as these are the samples from Oceania (there are only a few known Oceania samples and in past inference analyses no new samples are inferred as belonging to this group)
        & ~hgdp_tgp_outliers.contains(pop_pca_scores_ht.s),
        joint_meta.v3_meta.v3_project_pop,
    )

    if include_v2_known_in_training:
        training_pop = hl.if_else(
            hl.is_defined(training_pop),
            training_pop,
            hl.if_else(
                (hl.is_defined(joint_meta.v2_meta.v2_known_pop))
                & (hl.is_missing(joint_meta.v3_meta.v3_subsets.hgdp)),
                joint_meta.v2_meta.v2_known_pop,
                hl.missing(hl.tstr),
                missing_false=True,
            ),
            missing_false=True,
        )

    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        training_pop=training_pop,
        hgdp_or_tgp=joint_meta.v3_meta.v3_subsets.hgdp
        | joint_meta.v3_meta.v3_subsets.tgp,
    )

    if v4_population_spike:
        logger.info(
            "Spiking v4 pops, %s, into the RF training data...", v4_population_spike
        )

        pop_spiking = hl.dict(
            [
                (pop, V4_POP_SPIKE_DICT[pop])
                if pop in V4_POP_SPIKE_DICT
                else ValueError(f"Supplied pop: {pop} is not in V4_POP_SPIKE_DICT")
                for pop in v4_population_spike
            ]
        )

        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.case()
            .when(
                hl.is_defined(pop_pca_scores_ht.training_pop),
                pop_pca_scores_ht.training_pop,
            )
            .when(
                pop_spiking.contains(
                    joint_qc_meta.ht()[pop_pca_scores_ht.key].v4_race_ethnicity
                ),
                pop_spiking[
                    (joint_qc_meta.ht()[pop_pca_scores_ht.key].v4_race_ethnicity)
                ],
            )
            .or_missing()
        )

    if v3_population_spike:
        logger.info(
            "Spiking v3 pops, %s, into the RF training data", v3_population_spike
        )
        pop_spiking = hl.dict(
            [
                (pop, V3_SPIKE_PROJECTS[pop])
                if pop in V3_SPIKE_PROJECTS
                else ValueError(f"Supplied pop: {pop} is not in V3_SPIKE_PROJECTS")
                for pop in v3_population_spike
            ]
        )

        joint_qc_meta = joint_qc_meta.ht()
        v3_meta = v3_meta.ht()

        # TODO: Edit code once the joint QC is updated
        # Annotate research project for eas samples from Osaka
        v3_meta = v3_meta.annotate(
            project_meta=v3_meta.project_meta.annotate(
                research_project=hl.if_else(
                    (v3_meta.project_meta.project_pop == "eas")
                    & (v3_meta.s.startswith("JPNW")),
                    "osaka",
                    v3_meta.project_meta.research_project,
                )
            )
        )
        joint_qc_meta = joint_qc_meta.annotate(v3=v3_meta[joint_qc_meta.key])

        # Filter to only pre-determined list of v3 cohorts for the v3 spike-ins
        joint_qc_meta = joint_qc_meta.filter(
            hl.literal(V3_SPIKE_PROJECTS).contains(joint_qc_meta.v3_meta.v3_project_pop)
        )
        joint_qc_meta = joint_qc_meta.filter(
            hl.literal(V3_SPIKE_PROJECTS)[
                joint_qc_meta.v3_meta.v3_project_pop
            ].contains(joint_qc_meta.v3.project_meta.research_project)
        )

        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.case()
            .when(
                hl.is_defined(pop_pca_scores_ht.training_pop),
                pop_pca_scores_ht.training_pop,
            )
            .when(
                pop_spiking.contains(
                    joint_qc_meta[pop_pca_scores_ht.key].v3_meta.v3_project_pop
                ),
                (joint_qc_meta[pop_pca_scores_ht.key].v3_meta.v3_project_pop),
            )
            .or_missing()
        )

    # Keep track of original training population labels, this is useful if
    # samples are withheld to create PR curves
    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        original_training_pop=pop_pca_scores_ht.training_pop,
    )

    return pop_pca_scores_ht


def assign_pops(
    min_prob: float,
    include_unreleasable_samples: bool = False,
    pcs: List[int] = list(range(1, 21)),
    missing_label: str = "remaining",
    seed: int = 24,
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
    :param seed: Random seed, defaults to 24.
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
        seed,
        test,
        include_v2_known_in_training,
        v4_population_spike,
        v3_population_spike,
    )
    pop_field = "pop"
    logger.info(
        "Running RF for the first time with PCs %s using %d training examples: %s",
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

    pop_ht = pop_ht.annotate(
        original_training_pop=pop_pca_scores_ht[pop_ht.key].original_training_pop
    )
    pop_ht = pop_ht.annotate_globals(
        min_prob=min_prob,
        include_unreleasable_samples=include_unreleasable_samples,
        pcs=pcs,
        include_v2_known_in_training=include_v2_known_in_training,
        v3_population_spike=v3_population_spike,
        v4_population_spike=v4_population_spike,
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


def assign_pop_with_per_pop_probs(
    pop_ht: hl.Table,
    min_prob_decisions: Dict[str, float],
    missing_label: str = "remaining",
) -> hl.Table:
    """
    Assign samples to populations based on population-specific minimum RF probabilities.

    :param pop_ht: Table containing results of population inference.
    :param min_prob_decisions: Dictionary with population as key, and minimum RF probability required to assign a sample to that population as value.
    :param missing_label: Label for samples for which the assignment probability is smaller than required minimum probability.
    :return: Table with 'pop' annotation based on supplied per pop min probabilities rather than 'min-pop-prob'..
    """
    pop_ht = pop_ht.annotate(
        most_likely_pop=hl.rbind(
            hl.sorted(
                [
                    hl.struct(pop=x[-3:], prob=pop_ht[x])
                    for x in pop_ht.row
                    if x.startswith("prob_")
                ],
                key=lambda el: el[1],
                reverse=True,
            ),
            lambda pop_probs: hl.struct(
                pop=pop_probs[0][0],
                prob=pop_probs[0][1],
                second_prob=pop_probs[1][1],
                prob_diff=pop_probs[0][1] - pop_probs[1][1],
            ),
        )
    )
    pop_ht = pop_ht.annotate(
        pop=hl.if_else(
            pop_ht.most_likely_pop.prob
            > hl.literal(min_prob_decisions).get(pop_ht.most_likely_pop.pop),
            pop_ht.most_likely_pop.pop,
            missing_label,
        )
    )
    pop_ht = pop_ht.annotate_globals(min_prob_decisions=min_prob_decisions)

    return pop_ht


def main(args):
    """Assign global ancestry labels to samples."""
    try:
        include_unreleasable_samples = args.include_unreleasable_samples
        overwrite = args.overwrite
        test = args.test
        include_v2_known_in_training = args.include_v2_known_in_training

        if args.run_pca:
            pop_eigenvalues, pop_scores_ht, pop_loadings_ht = run_pca(
                pca_related_samples_to_drop().ht(),
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

            if args.min_prob_decisions_json:
                with hl.hadoop_open(args.min_prob_decisions_json, "r") as d:
                    min_prob_decisions = json.load(d)
                logger.info(
                    "Using min prob decisions per ancestry group: %s",
                    min_prob_decisions,
                )

                pop_ht = assign_pop_with_per_pop_probs(pop_ht, min_prob_decisions)

            logger.info("Writing pop ht...")
            pop_ht = pop_ht.checkpoint(
                get_pop_ht(
                    test=test,
                ).path,
                overwrite=overwrite,
                _read_if_exists=not overwrite,
            )

            with hl.hadoop_open(
                pop_rf_path(test=test),
                "wb",
            ) as out:
                pickle.dump(pops_rf_model, out)
    finally:
        hl.copy_log(
            f"gs://gnomad-tmp-4day/ancestry_assignment/ancestry_assignment.rf.log"
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
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
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )
    parser.add_argument(
        "--test", help="Run script on test dataset.", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="Overwrite output files.", action="store_true"
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
        "--min-prob-decisions-json",
        help=(
            "Optional path to JSON file containing decisions on minimum RF prob for pop"
            " assignment to use per each ancestry group instead of 'min-pop-prob'."
        ),
        type=str,
    )

    args = parser.parse_args()

    if args.slack_channel:
        from gnomad_qc.slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
