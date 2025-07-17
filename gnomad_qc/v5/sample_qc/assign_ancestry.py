"""Script to assign genetic ancestry group labels to samples using HGDP or TGP labels and spike ins from known v3 and v4 genetic ancestry labels."""

import argparse
import json
import logging
import pickle
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.sample_qc import hgdp_tgp_pop_outliers
from gnomad_qc.v4.resources.sample_qc import joint_qc_meta as v4_joint_qc_meta
from gnomad_qc.v4.sample_qc.assign_ancestry import V3_SPIKE_PROJECTS, V4_POP_SPIKE_DICT
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import project_meta, sample_id_collisions
from gnomad_qc.v5.resources.sample_qc import (
    gen_anc_rf_path,
    genetic_ancestry_pca_eigenvalues,
    genetic_ancestry_pca_loadings,
    genetic_ancestry_pca_scores,
    get_gen_anc_ht,
    get_gen_anc_pr_ht,
    get_joint_qc,
    per_grp_min_rf_probs_json_path,
    related_samples_to_drop,
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
    :param meta_ht: Table of meta data containing 'releasable' field, used to drop samples that are not releasable unless `include_unreleasable_samples` is True.
    :param related_samples_to_drop: Table of related samples to drop from PCA run.
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


def prep_ht_for_rf(
    gen_anc_pca_scores_ht: hl.Table,
    joint_meta_ht: hl.Table,
    include_v2_known_in_training: bool = False,
    v4_gen_anc_spike: Optional[List[str]] = None,
    v3_gen_anc_spike: Optional[List[str]] = None,
    include_v3_oceania: bool = False,
) -> hl.Table:
    """
    Prepare the PCA scores hail Table for the random forest genetic ancestry group assignment runs.

    Either train the RF with only HGDP and TGP, or HGDP and TGP and all v2 known labels.

    Can also specify list of genetic ancestry groups with known v3/v4 labels to include
    (v3_gen_anc_spike/v4_gen_anc_spike) for training. Known groups supplied for v3/v4 are
    specified by project name or self-reported ancestry and mapped to a genetic ancestry group using
    `V3_SPIKE_PROJECTS`/`V4_POP_SPIKE_DICT`.

    :param gen_anc_pca_scores_ht: Table of PCA scores to use for RF training.
    :param joint_meta_ht: Table of joint QC meta data.
    :param include_v2_known_in_training: Whether to train RF classifier using v2 known
        genetic ancestry labels. Default is False.
    :param v4_gen_anc_spike: Optional List of genetic ancestry groups to spike into training.
        Must be in V4_POP_SPIKE_DICT dictionary. Default is None.
    :param v3_gen_anc_spike: Optional List of genetic ancestry groups to spike into training.
        Must be in V3_SPIKE_PROJECTS dictionary. Default is None.
    :param include_v3_oceania: Whether to include v3 Oceania samples in training. Default is False.
    :return: Table with input for the random forest.
    """
    # Collect sample names of hgdp/tgp outliers to remove (these are outliers
    # found by Alicia Martin's group during group-specific PCA analyses as well
    # as one duplicate sample).
    hgdp_tgp_outliers_ht = hgdp_tgp_pop_outliers.ht()
    hgdp_tgp_outliers_ht = add_project_prefix_to_sample_collisions(
        t=hgdp_tgp_outliers_ht,
        sample_collisions=sample_id_collisions.ht(),
        project="gnomad",
    )
    hgdp_tgp_outliers = hl.literal(hgdp_tgp_outliers_ht.s.collect())

    joint_meta_ht = joint_meta_ht.annotate(
        hgdp_or_tgp=hl.or_else(
            joint_meta_ht.v3_meta.v3_subsets.hgdp
            | joint_meta_ht.v3_meta.v3_subsets.tgp,
            False,
        )
    )

    if include_v3_oceania:
        logger.info("Including v3 Oceania samples in training...")
        training_gen_anc = hl.or_missing(
            joint_meta_ht.hgdp_or_tgp & ~hgdp_tgp_outliers.contains(joint_meta_ht.s),
            hl.if_else(
                joint_meta_ht.v3_meta.v3_project_pop != "oth",
                joint_meta_ht.v3_meta.v3_project_pop,
                "oce",
            ),
        )
    else:
        logger.info("Excluding v3 Oceania samples from training...")
        # Note: Excluding v3_project_pop="oth" samples from training. These are samples from
        #  Oceania and there are only a few known Oceania samples, and in past inference
        #  analyses no new samples are inferred as belonging to this group.
        training_gen_anc = hl.or_missing(
            joint_meta_ht.hgdp_or_tgp
            & (joint_meta_ht.v3_meta.v3_project_pop != "oth")
            & ~hgdp_tgp_outliers.contains(joint_meta_ht.s),
            joint_meta_ht.v3_meta.v3_project_pop,
        )

    if include_v2_known_in_training:
        training_gen_anc = hl.coalesce(
            training_gen_anc,
            hl.or_missing(
                joint_meta_ht.v5_meta.data_type == "exomes",
                joint_meta_ht.v2_meta.v2_known_pop,
            ),
        )

    if v4_gen_anc_spike:
        logger.info(
            "Spiking v4 genetic ancestry groups, %s, into the RF training data...",
            v4_gen_anc_spike,
        )

        v4_spike_err = set(v4_gen_anc_spike) - set(V4_POP_SPIKE_DICT.keys())
        if len(v4_spike_err) > 0:
            raise ValueError(f"Pops: {v4_spike_err}, are not in V4_POP_SPIKE_DICT")

        v4_gen_anc_spike = {r: V4_POP_SPIKE_DICT[r] for r in v4_gen_anc_spike}

        training_gen_anc = hl.coalesce(
            training_gen_anc,
            hl.literal(v4_gen_anc_spike).get(joint_meta_ht.v4_meta.v4_race_ethnicity),
        )

    if v3_gen_anc_spike:
        logger.info(
            "Spiking v3 genetic ancestry groups, %s, into the RF training data",
            v3_gen_anc_spike,
        )
        v3_spike_err = set(v3_gen_anc_spike) - set(V3_SPIKE_PROJECTS.keys())
        if len(v3_spike_err) > 0:
            raise ValueError(f"Pops: {v3_spike_err}, are not in V3_SPIKE_PROJECTS")

        # Filter to only pre-determined list of v3 cohorts for the v3 spike-ins.
        training_gen_anc = hl.coalesce(
            training_gen_anc,
            hl.or_missing(
                hl.literal(V3_SPIKE_PROJECTS)
                .get(joint_meta_ht.v3_meta.v3_project_pop, hl.empty_array(hl.tstr))
                .contains(joint_meta_ht.v3_meta.v3_research_project),
                joint_meta_ht.v3_meta.v3_project_pop,
            ),
        )

    gen_anc_pca_scores_ht = gen_anc_pca_scores_ht.annotate(
        **joint_meta_ht.select(
            training_gen_anc=training_gen_anc,
            hgdp_or_tgp=joint_meta_ht.hgdp_or_tgp,
        )[gen_anc_pca_scores_ht.key]
    )

    return gen_anc_pca_scores_ht


def assign_gen_anc(
    gen_anc_pca_scores_ht: hl.Table,
    joint_meta_ht: hl.Table,
    min_prob: float,
    include_unreleasable_samples: bool = False,
    pcs: List[int] = list(range(1, 21)),
    missing_label: str = "remaining",
    include_v2_known_in_training: bool = False,
    v4_gen_anc_spike: Optional[List[str]] = None,
    v3_gen_anc_spike: Optional[List[str]] = None,
    n_partitions: int = 100,
    include_v3_oceania: bool = False,
) -> Tuple[hl.Table, Any]:
    """
    Use a random forest model to assign global genetic ancestry group labels based on the results from `run_pca`.

    Training data is the known label for HGDP and 1KG samples. All v2 samples with
    known genetic ancestry can be included if specificied. Can also
    specify a list of genetic ancestry groups with known v3/v4 labels to include
    (`v3_gen_anc_spike`/`v4_gen_anc_spike`) for training (see docstring of `prep_ht_for_rf`
    for more details on these spike groups). The method assigns a genetic ancestry group
    label to all samples in the dataset.

    :param gen_anc_pca_scores_ht: Table of PCA scores to use for RF training.
    :param joint_meta_ht: Table of joint QC meta data.
    :param min_prob: Minimum RF probability for genetic ancestry assignment.
    :param include_unreleasable_samples: Whether unreleasable samples were included in
        PCA.
    :param pcs: List of PCs to use in the RF.
    :param missing_label: Label for samples for which the assignment probability is
        smaller than `min_prob`.
    :param include_v2_known_in_training: Whether to train RF classifier using v2 known
        genetic ancestry labels. Default is False.
    :param v4_gen_anc_spike: Optional List of v4 genetic ancestry groups to spike into the RF.
        Must be in v4_gen_anc_spike dictionary. Defaults to None.
    :param v3_gen_anc_spike: Optional List of v3 genetic ancestry groups to spike into the RF.
        Must be in v3_gen_anc_spike dictionary. Defaults to None.
    :param n_partitions: Number of partitions to repartition the genetic ancestry group inference table to. Default is 100.
    :param include_v3_oceania: Whether to include v3 Oceania samples in training. Default is False.
    :return: Table of genetic ancestry group assignments and the RF model.
    """
    logger.info("Prepping HT for RF...")
    gen_anc_pca_scores_ht = prep_ht_for_rf(
        gen_anc_pca_scores_ht,
        joint_meta_ht,
        include_v2_known_in_training,
        v4_gen_anc_spike,
        v3_gen_anc_spike,
        include_v3_oceania,
    )
    gen_anc_field = "gen_anc"
    logger.info(
        "Running RF with PCs %s using %d training examples: %s",
        pcs,
        gen_anc_pca_scores_ht.aggregate(
            hl.agg.count_where(hl.is_defined(gen_anc_pca_scores_ht.training_gen_anc))
        ),
        gen_anc_pca_scores_ht.aggregate(
            hl.agg.counter(gen_anc_pca_scores_ht.training_gen_anc)
        ),
    )
    # Run the gen anc RF.
    gen_anc_ht, gen_anc_rf_model = assign_population_pcs(
        gen_anc_pca_scores_ht,
        pc_cols=pcs,
        known_col="training_gen_anc",
        output_col=gen_anc_field,
        min_prob=min_prob,
        missing_label=missing_label,
        n_partitions=n_partitions,
    )

    gen_anc_ht = gen_anc_ht.annotate_globals(
        min_prob=min_prob,
        include_unreleasable_samples=include_unreleasable_samples,
        pcs=pcs,
        include_v2_known_in_training=include_v2_known_in_training,
    )

    if v3_gen_anc_spike:
        gen_anc_ht = gen_anc_ht.annotate_globals(
            v3_gen_anc_spike=v3_gen_anc_spike,
        )
    if v4_gen_anc_spike:
        gen_anc_ht = gen_anc_ht.annotate_globals(
            v4_gen_anc_spike=v4_gen_anc_spike,
        )

    return gen_anc_ht, gen_anc_rf_model


def get_most_likely_gen_anc_expr(
    ht: hl.Table,
) -> Tuple[hl.expr.StructExpression, List[Tuple[str, str]]]:
    """
    Get StructExpression with 'gen_anc' and 'prob' for the most likely genetic ancestry group based on RF probabilities.

    :param ht: Input genetic ancestry group inference Table with random forest probabilities.
    :return: Struct Expression with 'gen_anc' and 'prob' for the highest RF probability.
    """
    # Get list of all genetic ancestry group RF probability annotations.
    prob_ann = [(x.split("prob_")[-1], x) for x in ht.row if x.startswith("prob_")]

    # Sort a list of all genetic ancestry group RF probabilities and keep the highest as the most
    # likely gen_anc.
    most_likely_gen_anc_expr = hl.rbind(
        hl.sorted(
            [hl.struct(gen_anc=gen_anc, prob=ht[ann]) for gen_anc, ann in prob_ann],
            key=lambda el: el["prob"],
            reverse=True,
        ),
        lambda gen_anc_probs: gen_anc_probs[0],
    )

    return most_likely_gen_anc_expr, prob_ann


def compute_precision_recall(ht: hl.Table, num_pr_points: int = 100) -> hl.Table:
    """
    Create Table with false positives (FP), true positives (TP), false negatives (FN), precision, and recall.

    Includes genetic ancestry group specific calculations.

    :param ht: Input genetic ancestry group inference Table with random forest probabilities.
    :param num_pr_points: Number of min prob cutoffs to compute PR metrics for.
    :return: Table with FP, TP, FN, precision, and recall.
    """
    # Use only RF evaluation samples to compute metrics.
    ht = ht.filter(ht.evaluation_sample)
    most_likely_gen_anc, prob_ann = get_most_likely_gen_anc_expr(ht)
    ht = ht.annotate(most_likely_gen_anc=most_likely_gen_anc)

    # Add a list of min_prob_cutoffs from 0 to 1 in increments of 1/num_pr_points and
    # explode so each sample has one row for each min_prob_cutoff.
    ht = ht.annotate(
        min_prob_cutoff=hl.literal(list(range(num_pr_points))) / float(num_pr_points)
    )
    ht = ht.explode(ht.min_prob_cutoff)

    # For each 'gen_anc', filter to any sample that has 'gen_anc' as the known genetic ancestry group or
    # most likely genetic ancestry group based on RF probabilities.
    gen_ancs = [None] + [gen_anc for gen_anc, _ in prob_ann]

    pr_agg = hl.struct()
    for gen_anc in gen_ancs:
        # Group by min_prob_cutoffs and aggregate to get TP, FP, and FN per
        # min_prob_cutoffs.
        meet_cutoff = ht.most_likely_gen_anc.prob >= ht.min_prob_cutoff
        same_training_most_likely_gen_anc = (
            ht.training_gen_anc == ht.most_likely_gen_anc.gen_anc
        )
        if gen_anc:
            is_most_likely_gen_anc = ht.most_likely_gen_anc.gen_anc == gen_anc
            is_training_gen_anc = ht.training_gen_anc == gen_anc
            tp = hl.agg.count_where(
                meet_cutoff & is_most_likely_gen_anc & is_training_gen_anc
            )
            fp = hl.agg.count_where(
                meet_cutoff & is_most_likely_gen_anc & ~is_training_gen_anc
            )
            fn = hl.agg.count_where(is_training_gen_anc & ~meet_cutoff)
        else:
            tp = hl.agg.count_where(meet_cutoff & same_training_most_likely_gen_anc)
            fp = hl.agg.count_where(meet_cutoff & ~same_training_most_likely_gen_anc)
            fn = hl.agg.count_where(~meet_cutoff)

        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        agg = hl.struct(TP=tp, FP=fp, FN=fn, precision=precision, recall=recall)
        if gen_anc:
            agg = hl.struct(**{gen_anc: agg})
        pr_agg = pr_agg.annotate(**agg)

    ht = ht.group_by("min_prob_cutoff").aggregate(**pr_agg)
    ht = ht.annotate_globals(gen_ancs=[gen_anc for gen_anc, _ in prob_ann])

    return ht


def infer_per_grp_min_rf_probs(
    ht: hl.Table, min_recall: float = 0.99, min_precision: float = 0.99
) -> Dict[str, Dict[str, float]]:
    """
    Infer per genetic ancestry group minimum RF probabilities from precision and recall values.

    Minimum recall (`min_recall`) is used to choose per genetic ancestry group minimum RF
    probabilities. This `min_recall` cutoff is applied first, and if the chosen minimum
    RF probabilities cutoff results in a precision lower than `min_precision`, the
    minimum RF probabilities with the highest recall that meets `min_precision` is
    used.

    :param ht: Precision recall Table returned by `compute_precision_recall`.
    :param min_recall: Minimum recall value to choose per genetic ancestry group minimum RF
        probabilities. Default is 0.99.
    :param min_precision: Minimum precision value to choose per genetic ancestry group minimum
        RF probabilities. Default is 0.99.
    :return: Dictionary of per genetic ancestry group min probability cutoffs `min_prob_cutoff`.
    """
    # Get list of all genetic ancestry groups.
    gen_ancs = hl.eval(ht.index_globals().gen_ancs)
    min_prob_per_gen_anc = {}
    for gen_anc in gen_ancs:
        pr_vals = hl.tuple(
            [ht.min_prob_cutoff, ht[gen_anc].precision, ht[gen_anc].recall]
        ).collect()
        min_probs, precisions, recalls = map(list, zip(*pr_vals))

        try:
            idx = next(
                x for x, r in reversed(list(enumerate(recalls))) if r >= min_recall
            )
        except StopIteration:
            raise ValueError(
                f"Recall never reaches {min_recall} for {gen_anc}. Recall values are: "
                f"{recalls}"
            )
        precision_cutoff = precisions[idx]
        if precision_cutoff < min_precision:
            try:
                idx = next(x for x, p in enumerate(precisions) if p >= min_precision)
            except StopIteration:
                raise ValueError(
                    f"Precision never reaches {min_precision} for {gen_anc}. Precision "
                    f"values are: {precisions}"
                )
        min_prob_per_gen_anc[gen_anc] = {
            "precision_cutoff": precisions[idx],
            "recall_cutoff": recalls[idx],
            "min_prob_cutoff": min_probs[idx],
        }

    return min_prob_per_gen_anc


def assign_gen_anc_with_per_gen_anc_probs(
    gen_anc_ht: hl.Table,
    min_prob_cutoffs: Dict[str, float],
    missing_label: str = "remaining",
) -> hl.Table:
    """
    Assign samples to genetic ancestry groups based on genetic ancestry group-specific minimum RF probabilities.

    :param gen_anc_ht: Table containing results of genetic ancestry group inference.
    :param min_prob_cutoffs: Dictionary with genetic ancestry group as key, and minimum RF
        probability required to assign a sample to that genetic ancestry group as value.
    :param missing_label: Label for samples for which the assignment probability is
        smaller than required minimum probability.
    :return: Table with 'gen_anc' annotation based on supplied per genetic ancestry group min probabilities.
    """
    min_prob_cutoffs = hl.literal(min_prob_cutoffs)
    gen_anc_prob, _ = get_most_likely_gen_anc_expr(gen_anc_ht)
    gen_anc = gen_anc_prob.gen_anc
    gen_anc_ht = gen_anc_ht.annotate(
        gen_anc=hl.if_else(
            gen_anc_prob.prob >= min_prob_cutoffs.get(gen_anc), gen_anc, missing_label
        )
    )
    gen_anc_ht = gen_anc_ht.annotate_globals(min_prob_cutoffs=min_prob_cutoffs)

    return gen_anc_ht


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
            check_resource_existence(
                output_step_resources={
                    "genetic_ancestry_pca_eigenvalues": genetic_ancestry_pca_eigenvalues(
                        include_unreleasable_samples, use_tmp_path
                    ).path,
                    "genetic_ancestry_pca_scores": genetic_ancestry_pca_scores(
                        include_unreleasable_samples, use_tmp_path
                    ).path,
                    "genetic_ancestry_pca_loadings": genetic_ancestry_pca_loadings(
                        include_unreleasable_samples, use_tmp_path
                    ).path,
                }
            )
            qc_mt = get_joint_qc(test=test).mt()

            if test_on_chr20:
                logger.info("Filtering QC MT to chromosome 20...")
                qc_mt = hl.filter_intervals(
                    qc_mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
                )

            gen_anc_eigenvalues, gen_anc_scores_ht, gen_anc_loadings_ht = run_pca(
                qc_mt=qc_mt,
                meta_ht=project_meta.ht(),
                related_samples_to_drop=related_samples_to_drop().ht(),
                # TODO: add in 'release' arg if it gets added to resource.
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
            check_resource_existence(
                output_step_resources={
                    "gen_anc_ht": get_gen_anc_ht(test=use_tmp_path).path
                }
            )

            gen_anc_pca_scores_ht = genetic_ancestry_pca_scores(
                include_unreleasable_samples=include_unreleasable_samples,
                test=use_tmp_path,
            ).ht()

            # Rename sample collision in v4 joint qc meta.
            v4_joint_qc_meta_ht = v4_joint_qc_meta.ht()
            v4_joint_qc_meta_ht = add_project_prefix_to_sample_collisions(
                t=v4_joint_qc_meta_ht,
                sample_collisions=sample_id_collisions.ht(),
                project="gnomad",
            )

            # Combine v4 and v5 project meta tables.
            meta_ht = project_meta.ht()
            meta_ht = meta_ht.select(v5_meta=hl.struct(**meta_ht.row_value))
            meta_ht = meta_ht.annotate(
                v4_meta=v4_joint_qc_meta_ht[meta_ht.s].select(
                    "v4_race_ethnicity",
                ),
                **v4_joint_qc_meta_ht[meta_ht.s].select("v2_meta", "v3_meta"),
            )

            gen_anc_pcs = args.gen_anc_pcs
            gen_anc_pcs = (
                list(range(1, gen_anc_pcs[0] + 1))
                if len(gen_anc_pcs) == 1
                else gen_anc_pcs
            )
            logger.info("Using following PCs: %s", gen_anc_pcs)
            gen_anc_ht, gen_anc_rf_model = assign_gen_anc(
                gen_anc_pca_scores_ht=gen_anc_pca_scores_ht,
                joint_meta_ht=meta_ht,
                min_prob=args.min_gen_anc_prob,
                include_unreleasable_samples=include_unreleasable_samples,
                pcs=gen_anc_pcs,
                include_v2_known_in_training=include_v2_known_in_training,
                v4_gen_anc_spike=args.v4_spike,
                v3_gen_anc_spike=args.v3_spike,
                n_partitions=args.gen_anc_partitions,
                include_v3_oceania=args.include_v3_oceania,
            )

            logger.info("Writing genetic ancestry group ht...")
            gen_anc_ht.write(
                get_gen_anc_ht(test=use_tmp_path).path, overwrite=overwrite
            )

            with hl.hadoop_open(
                gen_anc_rf_path(test=use_tmp_path),
                "wb",
            ) as out:
                pickle.dump(gen_anc_rf_model, out)

        if args.compute_precision_recall:
            check_resource_existence(
                output_step_resources={
                    "gen_anc_pr_ht": get_gen_anc_pr_ht(test=use_tmp_path).path
                }
            )

            ht = compute_precision_recall(
                get_gen_anc_ht(test=use_tmp_path).ht(),
                num_pr_points=args.number_pr_points,
            )
            ht.write(get_gen_anc_pr_ht(test=use_tmp_path).path, overwrite=overwrite)

        if args.apply_per_grp_min_rf_probs:
            check_resource_existence(
                output_step_resources={
                    "gen_anc_ht": get_gen_anc_ht(test=use_tmp_path).path
                }
            )

            gen_anc = get_gen_anc_ht(test=use_tmp_path)
            if args.infer_per_grp_min_rf_probs:
                min_probs = infer_per_grp_min_rf_probs(
                    get_gen_anc_pr_ht(test=use_tmp_path).ht(),
                    min_recall=args.min_recall,
                    min_precision=args.min_precision,
                )
                logger.info(
                    "Per genetic ancestry group min prob cutoff inference: %s",
                    min_probs,
                )
                min_probs = {
                    gen_anc: cutoff["min_prob_cutoff"]
                    for gen_anc, cutoff in min_probs.items()
                }
                with hl.hadoop_open(per_grp_min_rf_probs_json_path(), "w") as d:
                    d.write(json.dumps(min_probs))

            with hl.hadoop_open(per_grp_min_rf_probs_json_path(), "r") as d:
                min_probs = json.load(d)
            logger.info(
                "Using the following min prob cutoff per ancestry group: %s", min_probs
            )
            gen_anc_ht = gen_anc.ht().checkpoint(
                new_temp_file("gen_anc_ht", extension="ht")
            )
            gen_anc_ht = assign_gen_anc_with_per_gen_anc_probs(gen_anc_ht, min_probs)
            gen_anc_ht.write(gen_anc.path, overwrite=overwrite)
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
        "--include-v3-oceania",
        help="Whether to include v3 Oceania samples in training. Default is False.",
        action="store_true",
    )
    parser.add_argument(
        "--v4-spike",
        help="List of v4 exomes groups to spike into the RF training gen anc groups.",
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
        "--gen-anc-partitions",
        help="Number of partitions to repartition the genetic ancestry group inference table to. Defaults to 100.",
        default=100,
        type=int,
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
