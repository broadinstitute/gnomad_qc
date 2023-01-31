"""Add description."""
import argparse
import logging
from typing import List, Optional

import hail as hl
from gnomad.assessment.validity_checks import compare_row_counts
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SIBLINGS,
    UNRELATED,
)
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import all_ukb_samples_to_remove, get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta, project_meta
from gnomad_qc.v4.resources.sample_qc import (
    contamination,
    finalized_outlier_filtering,
    get_pop_ht,
    get_sample_qc,
    hard_filtered_samples,
    hard_filtered_samples_no_sex,
    joint_qc_meta,
    pca_related_samples_to_drop,
    platform,
    relatedness,
    release_related_samples_to_drop,
    sample_chr20_mean_dp,
    sample_qc_mt_callrate,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_metadata")
logger.setLevel(logging.INFO)


def reformat_sex_imputation_ht() -> hl.Table:
    """

    :return:
    """
    ht = sex.ht()
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    ht = ht.transmute(impute_sex_stats=hl.struct(**{x: ht[x] for x in impute_stats}))
    ht = ht.select(sex_imputation=hl.struct(**ht.row.drop("s")))
    ht = ht.select_globals(sex_imputation_ploidy_cutoffs=ht.globals)

    return ht


def reformat_hard_filters_ht() -> hl.Table:
    """

    :return:
    """
    ht = hard_filtered_samples.ht()

    ht_ex = ht.explode(ht.hard_filters)
    hard_filters = ht_ex.aggregate(hl.agg.collect_as_set(ht_ex.hard_filters))
    ht = ht.select(
        **{
            v: hl.if_else(
                hl.is_defined(ht.hard_filters),
                ht.hard_filters.contains(v),
                False,
            )
            for v in hard_filters
        },
        hard_filters=ht.hard_filters,
        hard_filtered=hl.is_defined(ht.hard_filters) & (hl.len(ht.hard_filters) > 0),
    )

    return ht


def reformat_relatedness_ht(ht, method: str = "cuking") -> hl.Table:
    rel_ht = get_relatedness_set_ht(relatedness(method).ht())

    # Annotating meta HT with related filter booleans.
    # Any sample that is hard filtered will have missing values for these bools.
    # Any sample that was filtered for relatedness will have True for
    #  sample_filters.related.
    # If a filtered related sample had a relationship with a higher degree than
    #  second-degree (duplicate, parent-child, sibling), that filter will also be True.
    rel_dict = {
        "related": hl.null(hl.tbool),
        "duplicate_or_twin": DUPLICATE_OR_TWINS,
        "parent_child": PARENT_CHILD,
        "sibling": SIBLINGS,
    }
    rel_set_expr_dict = {
        "release": hl.is_defined(pca_related_samples_to_drop.ht()[ht.key]),
        "all_samples": hl.is_defined(release_related_samples_to_drop.ht()[ht.key]),
    }
    ht = ht.select(
        "sample_filters",
        # relatedness_inference_relationships=rel_ht[ht.key].relationships,
    )
    ht = ht.annotate(
        sample_filters=ht.sample_filters.annotate(
            **{
                f"{k}_{rel}": hl.or_missing(
                    ~ht.sample_filters.hard_filtered,
                    v,  # & ht.relatedness_inference_relationships.contains(rel_val),
                )
                for rel, rel_val in rel_dict.items()
                for k, v in rel_set_expr_dict.items()
            }
        )
    )

    logger.info("Adding relatedness globals (cutoffs)")
    ht = ht.annotate_globals(relatedness_inference_cutoffs=ht.index_globals())

    return ht


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Parse relatedness Table to get every relationship (except UNRELATED) per sample.

    Return Table keyed by sample with all sample relationships in a set.

    :param Table relatedness_ht: Table with inferred relationship information output by pc_relate.
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationships annotated as a set.
    :rtype: hl.Table
    """
    relatedness_ht = relatedness_ht.filter(relatedness_ht.relationship != UNRELATED)
    relatedness_ht = relatedness_ht.select("relationship", s=relatedness_ht.i.s).union(
        relatedness_ht.select("relationship", s=relatedness_ht.j.s)
    )
    relatedness_ht = relatedness_ht.group_by(relatedness_ht.s).aggregate(
        relationships=hl.agg.collect_as_set(relatedness_ht.relationship)
    )
    return relatedness_ht


def get_relationship_filter_expr(
    hard_filtered_expr: hl.expr.BooleanExpression,
    relationship: str,
    relationship_set: hl.expr.SetExpression,
) -> hl.expr.builders.CaseBuilder:
    """
    Return case statement to populate relatedness filters in sample_filters struct.

    :param hl.expr.BooleanExpression hard_filtered_expr: Boolean for whether sample was hard filtered.
    :param str relationship: Relationship to check for. One of DUPLICATE_OR_TWINS, PARENT_CHILD, or SIBLINGS.
    :param hl.expr.SetExpression relationship_set: Set containing all possible relationship strings for sample.
    :return: Case statement used to population sample_filters related filter field.
    :rtype: hl.expr.builders.CaseBuilder
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.null(hl.tbool))
        .when(hl.is_defined(relationship_set), relationship_set.contains(relationship))
        .default(False)
    )


def name_me(
    left_ht: hl.Table,
    right_ht: hl.Table,
    logger_str=None,
    ann_expr=None,
    ann_label=None,
    global_expr=None,
    global_label=None,
    left_missing_approved: Optional[List[str]] = None,
    right_missing_approved: Optional[List[str]] = None,
    sample_count_match: bool = True,
) -> hl.Table:
    """
    Join left and right tables using specified keys and join types and returns result.

    Also prints warning if sample counts are not the same.

    :param Table left_ht: Left Table to be joined
    :param Table right_ht: Right Table to be joined
    :param bool sample_count_match: Are the sample counts expected to match in the tables
    :return: Table with annotations
    :rtype: Table
    """
    if logger_str is not None:
        logger.info("\n\nAnnotating with the %s Table.", logger_str)

    if sample_count_match:
        if not compare_row_counts(left_ht, right_ht):
            logger.warning("Sample counts in left and right tables do not match!")

            in_left_not_right = left_ht.anti_join(right_ht)
            if right_missing_approved:
                in_left_not_right = in_left_not_right.filter(
                    hl.literal(set(right_missing_approved)).contains(
                        in_left_not_right.s
                    ),
                    keep=False,
                )
            if in_left_not_right.count() != 0:
                logger.warning(
                    f"The following {in_left_not_right.count()} samples are found in "
                    "the left HT, but are not found in the right HT or in the "
                    "approved missing list:"
                )
                in_left_not_right.select().show(n=-1)
            elif right_missing_approved:
                logger.info(
                    "All samples found in the left HT, but not found in the right HT, "
                    "are in the approved missing list."
                )

            in_right_not_left = right_ht.anti_join(left_ht)
            if left_missing_approved:
                in_right_not_left = in_right_not_left.filter(
                    hl.literal(set(left_missing_approved)).contains(
                        in_right_not_left.s
                    ),
                    keep=False,
                )
            if in_right_not_left.count() != 0:
                logger.warning(
                    f"The following {in_right_not_left.count()} samples are found in "
                    "the right HT, but are not found in left HT or in the approved "
                    "missing list:"
                )
                in_right_not_left.select().show(n=-1)
            elif left_missing_approved:
                logger.info(
                    "All samples found in the right HT, but not found in the left HT, "
                    "are in the approved missing list."
                )
        else:
            logger.info("Sample counts match.")
    else:
        logger.info("No sample count check requested.")

    if ann_label is None:
        if ann_expr is None:
            ann_expr = right_ht[left_ht.key]
        else:
            ann_expr = ann_expr.annotate(**right_ht[left_ht.key])
    else:
        if ann_expr is None:
            ann_expr = hl.struct()
        ann_expr = ann_expr.annotate(**{ann_label: right_ht[left_ht.key]})

    if global_expr is None:
        global_expr = hl.struct()

    if global_label is None:
        global_expr = global_expr.annotate(**right_ht.index_globals())
    else:
        global_expr = global_expr.annotate(**{global_label: right_ht.index_globals()})

    return ann_expr, global_expr


def main(args):
    """Add description."""
    hl.init(log="/sample_metadata.log", default_reference="GRCh38")

    # Get list of UKB samples that should be removed.
    removed_ukb_samples = hl.import_table(
        all_ukb_samples_to_remove, no_header=True
    ).f0.collect()
    # Get list of hard filtered samples before sex imputation.
    hf_samples_no_sex = hard_filtered_samples_no_sex.ht().s.collect()
    # Get list of hard filtered samples with sex imputation.
    hf_samples = hard_filtered_samples.ht().s.collect()
    # Get list of v3 samples (expected in relatedness and pop).
    v3_samples = joint_qc_meta.ht().s.collect()

    logger.info("Loading the VDS columns to begin creation of the meta HT.")
    ht = get_gnomad_v4_vds(remove_hard_filtered_samples=False).variant_data.cols()

    # Note: 71 samples are found in the right HT, but are not found in left HT.
    #  They overlap with the UKB withheld samples indicating they were not removed
    #  when the metadata HT was created. Is this expected?
    ann_expr, global_expr = name_me(
        ht,
        project_meta.ht(),
        logger_str="project meta",
        left_missing_approved=removed_ukb_samples,
    )

    # Note: the withdrawn UKB list was updated after the sample QC HT creation, so
    #  the sample QC HT has 5 samples more in it than the final sample list.
    ann_expr, global_expr = name_me(
        ht,
        get_sample_qc("bi_allelic").ht(),
        ann_expr=ann_expr,
        ann_label="sample_qc",
        global_expr=global_expr,
        global_label="sample_qc_parameters",
        logger_str="sample QC",
        left_missing_approved=removed_ukb_samples,
    )

    # TODO: Number of PCs will be 9 because that is what was used, should we modify to
    #  have all 30 PCs, or add all 30 PCs to another annotation, or only keep the 9
    #  since that is all that was used?
    ann_expr, global_expr = name_me(
        ht,
        platform.ht().drop("gq_thresholds"),
        ann_expr=ann_expr,
        ann_label="platform_inference",
        global_expr=global_expr,
        global_label="platform_inference_parameters",
        logger_str="platform assignment",
        right_missing_approved=hf_samples_no_sex,
    )

    # TODO: Keep or drop is_female?
    # TODO: Should it be a struct with raw and adj or id this OK?
    #  chrx_frac_hom_alt: float64,
    #  chrx_frac_hom_alt_adj: float64,
    ann_expr, global_expr = name_me(
        ht,
        reformat_sex_imputation_ht(),
        ann_expr=ann_expr,
        global_expr=global_expr,
        logger_str="sex imputation",
        right_missing_approved=hf_samples_no_sex,
    )

    hard_filter_metrics_expr, hard_filter_global_expr = name_me(
        ht,
        contamination.ht(),
        global_label="contamination_approximation_parameters",
        logger_str="contamination approximation",
    )

    hard_filter_metrics_expr, hard_filter_global_expr = name_me(
        ht,
        sample_chr20_mean_dp.ht().drop("gq_thresholds"),
        ann_expr=hard_filter_metrics_expr,
        global_expr=hard_filter_global_expr,
        global_label="chr20_mean_dp_parameters",
        logger_str="chr20 sample mean DP",
    )

    # TODO: Should it be a struct with raw and adj or id this OK
    #  sample_qc_mt_callrate: float64,
    #  sample_qc_mt_callrate_adj: float64
    hard_filter_metrics_expr, hard_filter_global_expr = name_me(
        ht,
        sample_qc_mt_callrate.ht(),
        ann_expr=hard_filter_metrics_expr,
        global_expr=hard_filter_global_expr,
        global_label="sample_qc_mt_callrate_parameters",
        logger_str="sample QC MT callrate",
    )

    ann_expr = ann_expr.annotate(hard_filter_metrics=hard_filter_metrics_expr)
    global_expr = global_expr.annotate(hard_filter_parameters=hard_filter_global_expr)

    # TODO: How to handle PCs, different number in tables than used?
    # TODO: How to specify pop PC table to use?
    ann_expr, global_expr = name_me(
        ht,
        # get_pop_ht().ht(),
        hl.read_table(
            "gs://gnomad/v4.0/sample_qc/joint/ancestry_inference/gnomad.joint.v4.0"
            ".hgdp_tgp_training.pop.rf_w_16pcs_0.75_pop_prob.ht"
        ),
        ann_expr=ann_expr,
        ann_label="population_inference",
        global_expr=global_expr,
        global_label="population_inference_parameters",
        logger_str="population inference",
        left_missing_approved=v3_samples,
        right_missing_approved=hf_samples,
    )

    # TODO: rerun to get the callrate cutoff global annotation on the Table. Confirm
    #  results are identical otherwise
    sample_filters_expr, global_expr = name_me(
        ht,
        reformat_hard_filters_ht(),
        logger_str="hard filters",
        ann_expr=ann_expr,
        global_expr=global_expr,
        sample_count_match=False,
    )

    # TODO: Add relationship expression code.
    # TODO: Where to compute the release related samples to drop?
    # logger.info(logging_statement.format("PCA related samples to drop HT"))
    # ann_ht = reformat_relatedness_ht(ht)
    # ann_expr = ann_expr.annotate(**name_me(ht, ann_ht))
    # global_expr = global_expr.annotate(**ann_ht.index_globals())

    sample_filters_expr, global_expr = name_me(
        ht,
        finalized_outlier_filtering().ht(),
        logger_str="outlier detection",
        ann_expr=sample_filters_expr,
        global_expr=global_expr,
        global_label="outlier_detection_parameters",
        right_missing_approved=hf_samples,
        left_missing_approved=removed_ukb_samples,
    )
    ann_expr = ann_expr.annotate(sample_filters=sample_filters_expr)

    ht = ht.annotate(**ann_expr)
    ht = ht.annotate_globals(**global_expr)

    logger.info("Annotating high_quality field")
    ht = ht.annotate(
        high_quality=~ht.sample_filters.hard_filtered
        & ~ht.sample_filters.outlier_filtered
    )

    logger.info("Annotating releasable field")
    ht = ht.annotate(
        release=ht.project_meta.releasable
        & ht.high_quality
        # & ~ht.sample_filters.release_related
    )

    # ht = ht.checkpoint(meta.path, overwrite=args.overwrite)

    # logger.info(f"Release sample counts:{ht.aggregate(hl.agg.count_where(ht.release))}")
    ht.describe()
    # ht.summarize()
    # logger.info(f"Final sample count: {ht.count()}")

    # TODO: Add more nearest neighbor info?
    # TODO: Add trio info?
    # TODO: joint that has v3 info?


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )

    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
