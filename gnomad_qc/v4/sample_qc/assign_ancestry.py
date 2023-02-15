"""Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels."""
import argparse
import logging
import pickle
from typing import Any, List, Tuple

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
Dictionary with v3 pops as keys and approved cohorts to use for training for those pops as values.
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


def calculate_mislabeled_training(pop_ht: hl.Table, pop_field: str) -> [int, float]:
    """
    Calculate the number and proportion of mislabeled training samples.

    :param pop_ht: Table with assigned pops/subpops that is returned by `assign_population_pcs`.
    :param pop_field: Name of field in the Table containing the assigned pop/subpop.
    :return: The number and proportion of mislabeled training samples.
    """
    logger.info("Calculating mislabeled training samples")
    n_mislabeled_samples = pop_ht.aggregate(
        hl.agg.count_where(pop_ht.training_pop != pop_ht[pop_field])
    )

    logger.info("Counting where training pop is defined")
    defined_training_pops = pop_ht.aggregate(
        hl.agg.count_where(hl.is_defined(pop_ht.training_pop))
    )

    prop_mislabeled_samples = n_mislabeled_samples / defined_training_pops

    return n_mislabeled_samples, prop_mislabeled_samples


def prep_ht_for_rf(
    include_unreleasable_samples: bool = False,
    withhold_prop: hl.float = None,
    seed: int = 24,
    test: bool = False,
    only_train_on_hgdp_tgp: bool = False,
    v4_population_spike: List[str] = None,
    v3_population_spike: List[str] = None,
) -> hl.Table:
    """
    Prepare the PCA scores hail Table for the random forest population assignment runs.

    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA.
    :param withhold_prop: Proportion of training pop samples to withhold from training. All samples are kept when set as `None`.
    :param seed: Random seed, defaults to 24.
    :param test: Whether RF should run on the test QC MT.
    :param only_train_on_hgdp_tgp: Whether to train RF classifier using only the HGDP and 1KG populations. Default is False.
    :param v4_population_spike: List of populations to spike into training. Must be in v4_population_spike dictionary. Default is None.
    :param v3_population_spike: List of populations to spike into training. Must be in v3_population_spike dictionary. Default is None.
    :return Table with input for the random forest.
    """
    pop_pca_scores_ht = ancestry_pca_scores(include_unreleasable_samples, test).ht()
    global joint_qc_meta
    global v3_meta

    joint_meta = joint_qc_meta.ht()[pop_pca_scores_ht.key]

    hgdp_tgp_outliers = hgdp_tgp_pop_outliers.ht()
    hgdp_tgp_outliers = hgdp_tgp_outliers.s.collect()

    # TODO: Add code for subpopulations
    # Either train the RF with only HGDP and TGP, or HGDP and TGP and all v2 known labels
    if only_train_on_hgdp_tgp:
        logger.info("Training using HGDP and 1KG samples only...")
        training_pop = hl.or_missing(
            (joint_meta.v3_meta.v3_subsets.hgdp | joint_meta.v3_meta.v3_subsets.tgp)
            & (joint_meta.v3_meta.v3_project_pop != "oth")
            & (~hl.literal(hgdp_tgp_outliers).contains(pop_pca_scores_ht.s)),
            joint_meta.v3_meta.v3_project_pop,
        )
    else:
        logger.info(
            "Training using HGDP and 1KG samples and all v2 samples with known pops..."
        )
        training_pop = hl.if_else(
            (joint_meta.v3_meta.v3_subsets.hgdp | joint_meta.v3_meta.v3_subsets.tgp)
            & (joint_meta.v3_meta.v3_project_pop != "oth")
            & (~hl.literal(hgdp_tgp_outliers).contains(pop_pca_scores_ht.s)),
            joint_meta.v3_meta.v3_project_pop,
            hl.if_else(
                (hl.is_defined(joint_meta.v2_meta.v2_known_pop))
                & (hl.is_missing(joint_meta.v3_meta.v3_subsets.hgdp)),
                joint_meta.v2_meta.v2_known_pop,
                hl.null(hl.tstr),
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
            "Spiking v4 pops, %s, into the RF training data", v4_population_spike
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
            hl.literal(V3_SPIKE_PROJECTS)[joint_qc_meta.v3_meta.v3_project_pop].contains(
                joint_qc_meta.v3.project_meta.research_project
            )
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
                pop_spiking[
                    (joint_qc_meta[pop_pca_scores_ht.key].v3_meta.v3_project_pop)
                ],
            )
            .or_missing()
        )

    # Keep track of original training population labels, this is useful if
    # samples are withheld to create PR curves
    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        original_training_pop=pop_pca_scores_ht.training_pop,
    )

    # Use the withhold proportion to create PR curves for when the RF removes
    # samples because it will remove samples that are misclassified
    if withhold_prop:
        logger.info("Withholding proportion is %s", withhold_prop)
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                hl.is_defined(pop_pca_scores_ht.training_pop)
                & hl.rand_bool(1.0 - withhold_prop, seed=seed),
                pop_pca_scores_ht.training_pop,
            )
        )

        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            withheld_sample=hl.is_defined(pop_pca_scores_ht.original_training_pop)
            & (~hl.is_defined(pop_pca_scores_ht.training_pop))
        )
    return pop_pca_scores_ht


def assign_pops(
    min_prob: float,
    include_unreleasable_samples: bool = False,
    max_number_mislabeled_training_samples: int = None,
    max_proportion_mislabeled_training_samples: float = None,
    pcs: List[int] = list(range(1, 21)),
    withhold_prop: float = None,
    missing_label: str = "remaining",
    seed: int = 24,
    test: bool = False,
    overwrite: bool = False,
    only_train_on_hgdp_tgp: bool = False,
    v4_population_spike: List[str] = None,
    v3_population_spike: List[str] = None,
) -> Tuple[hl.Table, Any]:
    """
    Use a random forest model to assign global population labels based on the results from `run_pca`.

    Training data is the v3 project metadata field `project_pop` when defined, otherwise, the v3 inferred population
    with the exception of `oth`. Training data is restricted to 1KG and HGDP samples, if specified. The method assigns
    a population label to all samples in the dataset. If a maximum number or proportion of mislabeled samples is set and the number
    of truth outliers exceeds this threshold, the outliers are removed from the training data, and the random
    forest runs until the number of truth outliers is less than or equal to `max_number_mislabeled_training_samples` or
    `max_proportion_mislabeled_training_samples`. Only one of `max_number_mislabeled_training_samples`
    or `max_proportion_mislabeled_training_samples` can be set.
    :param min_prob: Minimum RF probability for pop assignment.
    :param include_unreleasable_samples: Whether unreleasable samples were included in PCA.
    :param max_number_mislabeled_training_samples: If set, run the population assignment until the number of mislabeled training samples is less than this number threshold.
    :param max_proportion_mislabeled_training_samples: If set, run the population assignment until the number of mislabeled training samples is less than this proportion threshold.
    :param pcs: List of PCs to use in the RF.
    :param withhold_prop: Proportion of training pop samples to withhold from training. Keep all samples if `None`.
    :param missing_label: Label for samples for which the assignment probability is smaller than `min_prob`.
    :param seed: Random seed, defaults to 24.
    :param test: Whether running assigment on a test dataset.
    :param overwrite: Whether to overwrite existing files.
    :param only_train_on_hgdp_tgp: Whether to train the RF classifier using only the HGDP and 1KG populations. Defaults to False.
    :param v4_population_spike: List of v4 populations to spike into the RF. Must be in v4_pop_spike dictionary. Defaults to None.
    :param v3_population_spike: List of v3 populations to spike into the RF. Must be in v4_pop_spike dictionary. Defaults to None.
    :return: Table of pop assignments and the RF model.
    """
    logger.info("Assigning global population labels")

    logger.info("Checking passed args...")
    if (
        max_number_mislabeled_training_samples is not None
        and max_proportion_mislabeled_training_samples is not None
    ):
        raise ValueError(
            "Only one of max_number_mislabeled_training_samples or"
            " max_proportion_mislabeled_training_samples can be set"
        )
    elif max_proportion_mislabeled_training_samples is not None:
        max_mislabeled = max_proportion_mislabeled_training_samples
    elif max_number_mislabeled_training_samples is not None:
        max_mislabeled = max_number_mislabeled_training_samples
    else:
        max_mislabeled = None

    logger.info("Prepping HT for RF.")
    pop_pca_scores_ht = prep_ht_for_rf(
        include_unreleasable_samples,
        withhold_prop,
        seed,
        test,
        only_train_on_hgdp_tgp,
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
    # Run the pop RF for the first time
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
    # Calculate number and proportion of mislabeled samples
    n_mislabeled_samples, prop_mislabeled_samples = calculate_mislabeled_training(
        pop_ht, pop_field
    )

    mislabeled = (
        n_mislabeled_samples
        if max_number_mislabeled_training_samples
        # TODO: Run analysis on what makes sense as input here, utilize withhold
        # arg for PR curves
        else prop_mislabeled_samples
    )
    pop_assignment_iter = 1
    logger.info(
        "The maximum number of samples with conflicting training and assigned labels"
        " allowed is %s.",
        str(max_mislabeled),
    )
    # Rerun the RF until the number of mislabeled samples (known pop !=
    # assigned pop) is below our max mislabeled threshold. The basis of this
    # decision should be rooted in the reliability of the samples' provided
    # population label.
    if max_mislabeled:
        while mislabeled > max_mislabeled:
            pop_assignment_iter += 1
            logger.info(
                "Found %i(%d%%) training samples labeled differently from their known"
                " pop. Re-running assignment without them.",
                n_mislabeled_samples,
                round(prop_mislabeled_samples * 100, 2),
            )

            pop_ht = pop_ht[pop_pca_scores_ht.key]

            # Remove mislabled samples from RF training unless they are from 1KG or HGDP
            logger.info(
                "Reannotating training pop after removing misclassified samples unless"
                " part of HGDP or 1KG"
            )
            pop_pca_scores_ht = pop_pca_scores_ht.annotate(
                training_pop=hl.or_missing(
                    (
                        (pop_ht.training_pop == pop_ht[pop_field])
                        | (pop_pca_scores_ht.hgdp_or_tgp)
                    ),
                    pop_pca_scores_ht.training_pop,
                ),
            )
            logger.info("Checkpointing reannotated scores ht")
            pop_pca_scores_ht = pop_pca_scores_ht.checkpoint(
                get_checkpoint_path(
                    f"assign_pops_rf_iter_{pop_assignment_iter}{'_test' if test else ''}_lots_of_loggers"
                ),
                overwrite=overwrite,
            )

            logger.info(
                "Counted %d training examples",
                pop_pca_scores_ht.aggregate(
                    hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
                ),
            )
            logger.info("Rerunning gnomAD method's assign population_pcs")
            pop_ht, pops_rf_model = assign_population_pcs(
                pop_pca_scores_ht,
                pc_cols=pcs,
                known_col="training_pop",
                output_col=pop_field,
                min_prob=min_prob,
                missing_label=missing_label,
            )
            logger.info("Completed gnomAD method's assign population_pcs")
            pop_ht = pop_ht.checkpoint(
                get_checkpoint_path(
                    f"assign_pops_rf_iter_{pop_assignment_iter}{'_test' if test else ''}pop_ht_lots_of_loggers_pcs{pcs[-1]}"
                ),
                overwrite=overwrite,
            )
            (
                n_mislabeled_samples,
                prop_mislabeled_samples,
            ) = calculate_mislabeled_training(pop_ht, pop_field)

            mislabeled = (
                n_mislabeled_samples
                if max_number_mislabeled_training_samples
                # TODO: Run analysis on what makes sense as input here, utilize withhold
                # arg for PR curves
                else prop_mislabeled_samples
            )

    pop_ht = pop_ht.annotate(
        original_training_pop=pop_pca_scores_ht[pop_ht.key].original_training_pop
    )
    pop_ht = pop_ht.annotate_globals(
        min_prob=min_prob,
        include_unreleasable_samples=include_unreleasable_samples,
        pop_assignment_iterations=pop_assignment_iter,
        pcs=pcs,
        only_train_on_hgdp_tgp=only_train_on_hgdp_tgp,
        v3_population_spike=v3_population_spike,
        v4_population_spike=v4_population_spike,
    )

    if max_mislabeled:
        pop_ht = pop_ht.annotate_globals(
            max_mislabeled=max_mislabeled,
        )

    if withhold_prop:
        pop_ht = pop_ht.annotate_globals(withhold_prop=withhold_prop)
        pop_ht = pop_ht.annotate(
            withheld_sample=pop_pca_scores_ht[pop_ht.key].withheld_sample
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


def main(args):
    """Assign global ancestry labels to samples."""
    try:
        include_unreleasable_samples = args.include_unreleasable_samples
        overwrite = args.overwrite
        test = args.test
        only_train_on_hgdp_tgp = args.only_train_on_hgdp_tgp

        if (
            args.max_number_mislabeled_training_samples is not None
            and args.max_proportion_mislabeled_training_samples is not None
        ):
            raise ValueError(
                "Only one of max_number_mislabeled_training_samples or"
                " max_proportion_mislabeled_training_samples can be set"
            )
        elif args.max_proportion_mislabeled_training_samples is not None:
            max_mislabeled = args.max_proportion_mislabeled_training_samples
        elif args.max_number_mislabeled_training_samples is not None:
            max_mislabeled = args.max_number_mislabeled_training_samples
        else:
            max_mislabeled = None

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
                max_number_mislabeled_training_samples=args.max_number_mislabeled_training_samples,
                max_proportion_mislabeled_training_samples=args.max_proportion_mislabeled_training_samples,
                pcs=pop_pcs,
                withhold_prop=args.withhold_prop,
                test=test,
                overwrite=overwrite,
                only_train_on_hgdp_tgp=only_train_on_hgdp_tgp,
                v4_population_spike=args.v4_population_spike,
                v3_population_spike=args.v3_population_spike,
            )
            logger.info("Writing pop ht...")
            pop_ht = pop_ht.checkpoint(
                get_pop_ht(
                    test=test, only_train_on_hgdp_tgp=only_train_on_hgdp_tgp
                ).path,
                overwrite=overwrite,
                _read_if_exists=not overwrite,
            )

            with hl.hadoop_open(
                pop_rf_path(test=test, only_train_on_hgdp_tgp=only_train_on_hgdp_tgp),
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
        "--withhold-prop",
        help=(
            "Proportion of training samples to withhold from pop assignment RF"
            " training. All samples with known labels will be used for training if this"
            " flag is not used."
        ),
        type=float,
    )
    mislabel_parser = parser.add_mutually_exclusive_group(required=False)
    mislabel_parser.add_argument(
        "--max-number-mislabeled-training-samples",
        help=(
            "If set, run the population assignment until number of training samples"
            " that are mislabelled is under this number threshold. The basis of this"
            " decision should be rooted in the reliability of the samples' provided"
            " population label. Can't be used if"
            " `max-proportion-mislabeled-training-samples` is already set."
        ),
        type=int,
        default=None,
    )
    mislabel_parser.add_argument(
        "--max-proportion-mislabeled-training-samples",
        help=(
            "If set, run the population assignment until number of training samples"
            " that are mislabelled is under this proportion threshold. The basis of"
            " this decision should be rooted in the reliability of the samples'"
            " provided population label. Can't be used if"
            " `max-number-mislabeled-training-samples` is already set."
        ),
        type=float,
        default=None,
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
        "--only-train-on-hgdp-tgp",
        help="Whether to train RF classifier using only the HGDP and TGP populations.",
        action="store_true",
    )
    parser.add_argument(
        "--v4-population-spike",
        help="List of v4 populations to spike into the RF training populations.",
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--v3-population-spike",
        help="List of v3 populations to spike into the RF training populations.",
        type=str,
        nargs="+",
    )

    args = parser.parse_args()

    if args.slack_channel:
        from gnomad_qc.slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
