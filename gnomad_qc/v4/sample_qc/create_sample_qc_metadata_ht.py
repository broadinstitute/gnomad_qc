"""Add description."""
import argparse
import logging

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
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta, project_meta
from gnomad_qc.v4.resources.sample_qc import (
    finalized_outlier_filtering,
    get_pop_ht,
    get_sample_qc,
    hard_filtered_samples,
    pca_related_samples_to_drop,
    relatedness,
    release_related_samples_to_drop,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_qc")
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
    ht = ht.transmute(
        sample_filters=hl.struct(
            **{
                v: hl.if_else(
                    hl.is_defined(ht.hard_filters),
                    ht.hard_filters.contains(v),
                    False,
                )
                for v in hard_filters
            },
            hard_filters=ht.hard_filters,
            hard_filtered=hl.is_defined(ht.hard_filters)
            & (hl.len(ht.hard_filters) > 0),
        )
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


def reformat_outlier_filtering_ht() -> hl.Table:
    """

    :return:
    """
    ht = finalized_outlier_filtering().ht()
    ht = ht.transmute(
        sample_filters=ht.sample_filters.annotate(
            **{x: ht[x] for x in ht.row if x.startswith("fail_")},
            qc_metrics_filters=ht.qc_metrics_filters,
        )
    )
    ht = ht.select_globals(outlier_detection_metrics=ht.index_globals())

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
    # left_key: str,
    right_ht: hl.Table,
    # right_key: str,
    sample_count_match: bool = True,
) -> hl.Table:
    """
    Join left and right tables using specified keys and join types and returns result.

    Also prints warning if sample counts are not the same.

    :param Table left_ht: Left Table to be joined
    :param str left_key: Key of left Table
    :param Table right_ht: Right Table to be joined
    :param str right_key: Key of right Table
    :param str join_type: Type of join
    :param bool sample_count_match: Are the sample counts expected to match in the tables
    :return: Table with annotations
    :rtype: Table
    """
    if sample_count_match and not compare_row_counts(left_ht, right_ht):
        logger.warning("Sample counts in left and right tables do not match!")

        in_left_not_right = left_ht.anti_join(right_ht)
        if in_left_not_right.count() != 0:
            logger.warning(
                f"The following {in_left_not_right.count()} samples are found in the"
                " left HT, but are not found in the right HT"
            )
            in_left_not_right.select().show(n=-1)

        in_right_not_left = right_ht.anti_join(left_ht)
        if in_right_not_left.count() != 0:
            logger.warning(
                f"The following {in_right_not_left.count()} samples are found in the"
                " right HT, but are not found in left HT"
            )
            in_right_not_left.select().show(n=-1)

    return right_ht[left_ht.key]


def main(args):
    """Add description."""
    hl.init(log="/hail.log", default_reference="GRCh38")
    logging_statement = "Reading in the {}"

    logger.info(
        "Loading project metadata file with subset, age, and releasable information to "
        "begin creation of the meta HT"
    )
    ht = get_gnomad_v4_vds(remove_hard_filtered_samples=False).variant_data.cols()

    logger.info(logging_statement.format("project meta HT"))
    ann_expr = name_me(ht, project_meta.ht())

    # TODO: Add platform
    logger.info(logging_statement.format("sex HT"))
    ann_ht = reformat_sex_imputation_ht()
    ann_expr = ann_expr.annotate(**name_me(ht, ann_ht))
    global_expr = ann_ht.index_globals()

    logger.info(logging_statement.format("sample QC HT"))
    ann_ht = get_sample_qc("bi_allelic").ht()
    ann_expr = ann_expr.annotate(sample_qc=name_me(ht, ann_ht))
    global_expr = global_expr.annotate(**ann_ht.index_globals())

    logger.info(logging_statement.format("population PCA HT"))
    # ann_ht = get_pop_ht().ht()
    ann_ht = hl.read_table(
        "gs://gnomad/v4.0/sample_qc/joint/ancestry_inference/gnomad.joint.v4.0.hgdp_tgp_training.pop.rf_w_16pcs_0.75_pop_prob.ht"
    )
    ann_expr = ann_expr.annotate(
        population_inference=name_me(ht, ann_ht, sample_count_match=False)
    )
    global_expr = global_expr.annotate(
        population_inference_pca_metrics=ann_ht.index_globals()
    )

    logger.info(logging_statement.format("hard filters HT"))
    ann_ht = reformat_hard_filters_ht()
    ann_expr = ann_expr.annotate(**name_me(ht, ann_ht, sample_count_match=False))
    global_expr = global_expr.annotate(**ann_ht.index_globals())

    # logger.info(logging_statement.format("PCA related samples to drop HT"))
    # TODO: Start back up here
    # ann_ht = reformat_relatedness_ht(ht)
    # ann_expr = ann_expr.annotate(**name_me(ht, ann_ht))
    # global_expr = global_expr.annotate(**ann_ht.index_globals())

    # TODO: what to add for outlier filtering?
    logger.info(logging_statement.format("outlier HT"))
    ann_ht = reformat_outlier_filtering_ht()
    ann_expr = ann_expr.annotate(**name_me(ht, ann_ht))
    global_expr = global_expr.annotate(**ann_ht.index_globals())

    ht = ht.annotate(**ann_expr)
    ht = ht.annotate_gobals(**global_expr)

    logger.info("Annotating high_quality field")
    ht = ht.annotate(
        high_quality=~ht.sample_filters.hard_filtered
        & (hl.len(ht.sample_filters.qc_metrics_filters) == 0)
    )

    logger.info("Annotating releasable field")
    ht = ht.annotate(
        release=ht.project_meta.releasable
        & ht.high_quality
        & ~ht.sample_filters.release_related
    )

    ht = ht.checkpoint(
        meta.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite
    )

    logger.info(f"Release sample counts:{ht.aggregate(hl.agg.count_where(ht.release))}")
    ht.describe()
    ht.summarize()
    logger.info(f"Final sample count: {ht.count()}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )

    parser.add_argument(
        "--regressed_metrics_outlier",
        help="Should metadata HT use regression outlier model.",
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
