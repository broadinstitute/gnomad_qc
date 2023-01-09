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

from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta, meta_tsv_path, project_meta
from gnomad_qc.v3.resources.sample_qc import (
    get_sample_qc,
    hard_filtered_samples,
    pca_related_samples_to_drop,
    picard_metrics,
    pop,
    regressed_metrics,
    relatedness,
    release_related_samples_to_drop,
    sex,
    stratified_metrics,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_qc")
logger.setLevel(logging.INFO)


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


def join_tables(
    left_ht: hl.Table,
    left_key: str,
    right_ht: hl.Table,
    right_key: str,
    join_type: str,
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

        if join_type != "outer":
            logger.warning(
                "Join type is not an outer join so some samples will be filtered!"
            )

    return left_ht.key_by(left_key).join(right_ht.key_by(right_key), how=join_type)


def generate_metadata(regressed_metrics_outlier: bool = True) -> hl.Table:
    """
    Pull all sample QC information together in one Table.

    :param regressed_metrics_outlier: Should the outlier table be from the residuals instead of pop stratified
    :return: Table annotated with all metadata
    :rtype: hl.Table
    """
    logging_statement = "Reading in {} and joining with meta HT"

    logger.info(
        "Loading metadata file with subset, age, and releasable information to begin"
        " creation of the meta HT"
    )
    left_ht = get_gnomad_v3_mt(remove_hard_filtered_samples=False).cols()

    right_ht = project_meta.ht()
    right_ht = right_ht.select(
        project_meta=hl.struct(**right_ht.row.drop(*(SUBSETS + ["s"]))),
        subsets=hl.struct(**{x: right_ht[x] for x in SUBSETS}),
    )
    left_ht = join_tables(left_ht, "s", right_ht.select_globals(), "s", "right")

    logger.info(logging_statement.format("picard metric HT"))
    right_ht = picard_metrics.ht()
    right_ht = right_ht.select("bam_metrics")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left", sample_count_match=False)

    logger.info(logging_statement.format("sex HT"))
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    right_ht = sex.ht()
    right_ht = right_ht.transmute(
        impute_sex_stats=hl.struct(**{x: right_ht[x] for x in impute_stats})
    )

    # Create struct for join
    right_ht = right_ht.select(sex_imputation=hl.struct(**right_ht.row.drop("s")))
    right_ht = right_ht.select_globals(sex_imputation_ploidy_cutoffs=right_ht.globals)
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info(logging_statement.format("sample QC HT"))
    right_ht = get_sample_qc("bi_allelic").ht()

    # Remove annotations that cannot be computed from the sparse format
    not_in_sparse = ["n_called", "n_not_called", "n_filtered", "call_rate"]
    right_ht = right_ht.annotate(
        **{x: right_ht[x].drop(*not_in_sparse) for x in right_ht.row_value}
    )
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info(logging_statement.format("population PCA HT"))
    right_ht = pop.ht()
    right_ht = right_ht.select_globals(
        population_inference_pca_metrics=right_ht.globals
    )
    right_ht = right_ht.select(population_inference=hl.struct(**right_ht.row.drop("s")))
    left_ht = join_tables(
        left_ht, "s", right_ht, "s", "outer", sample_count_match=False
    )

    logger.info(
        "Reading hard filters HT, renaming hard filters struct to sample_filters, and"
        " joining with meta HT"
    )
    right_ht = hard_filtered_samples.ht()
    left_ht = join_tables(
        left_ht, "s", right_ht, "s", "outer", sample_count_match=False
    )

    # Change sample_filters to a struct
    ex_right_ht = right_ht.explode(right_ht.hard_filters)
    hard_filters = ex_right_ht.aggregate(
        hl.agg.collect_as_set(ex_right_ht.hard_filters)
    )
    left_ht = left_ht.transmute(
        sample_filters=hl.struct(
            **{
                v: hl.if_else(
                    hl.is_defined(left_ht.hard_filters),
                    left_ht.hard_filters.contains(v),
                    False,
                )
                for v in hard_filters
            },
            hard_filters=left_ht.hard_filters,
            hard_filtered=hl.is_defined(left_ht.hard_filters)
            & (hl.len(left_ht.hard_filters) > 0),
        )
    )

    logger.info(
        "Reading in PCA related samples to drop HT and preparing to annotate meta HT's"
        " sample_filter struct with relatedness booleans"
    )
    related_samples_to_drop_ht = pca_related_samples_to_drop.ht()
    release_related_samples_to_drop_ht = release_related_samples_to_drop.ht()
    relatedness_ht = get_relatedness_set_ht(relatedness.ht())
    related_samples_to_drop_ht = related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[related_samples_to_drop_ht.s].relationships
    )
    release_related_samples_to_drop_ht = release_related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[release_related_samples_to_drop_ht.s].relationships
    )

    # Annotating meta HT with related filter booleans
    # Any sample that is hard filtered will have missing values for these bools
    # Any sample that was filtered for relatedness will have True for sample_filters.related
    # If a filtered related sample had a relationship with a higher degree than second-degree (duplicate, parent-child, sibling),
    # that filter will also be True
    release_else_expr = release_related_samples_to_drop_ht[left_ht.key].relationships
    all_else_expr = related_samples_to_drop_ht[left_ht.key].relationships
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            release_related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(release_related_samples_to_drop_ht[left_ht.key]),
            ),
            release_duplicate=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                DUPLICATE_OR_TWINS,
                release_else_expr,
            ),
            release_parent_child=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                PARENT_CHILD,
                release_else_expr,
            ),
            release_sibling=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                SIBLINGS,
                release_else_expr,
            ),
            all_samples_related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(related_samples_to_drop_ht[left_ht.key]),
            ),
            all_samples_duplicate=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                DUPLICATE_OR_TWINS,
                all_else_expr,
            ),
            all_samples_parent_child=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                PARENT_CHILD,
                all_else_expr,
            ),
            all_samples_sibling=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                SIBLINGS,
                all_else_expr,
            ),
        )
    )
    left_ht = left_ht.annotate(
        relatedness_inference=hl.struct(
            relationships=relatedness_ht[left_ht.s].relationships,
        )
    )

    logger.info("Adding relatedness globals (cutoffs)")
    left_ht = left_ht.annotate_globals(
        relatedness_inference_cutoffs=hl.struct(**relatedness_ht.index_globals())
    )

    logger.info(logging_statement.format("outlier HT"))
    if regressed_metrics_outlier:
        right_ht = regressed_metrics.ht()
    else:
        right_ht = stratified_metrics.ht()

    right_ht = right_ht.select_globals(
        outlier_detection_metrics=hl.struct(
            **right_ht.index_globals(), used_regressed_metrics=regressed_metrics_outlier
        )
    )

    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
    left_ht = left_ht.transmute(
        sample_filters=left_ht.sample_filters.annotate(
            **{x: left_ht[x] for x in left_ht.row if x.startswith("fail_")},
            qc_metrics_filters=left_ht.qc_metrics_filters,
        )
    )
    if regressed_metrics_outlier:
        left_ht = left_ht.transmute(
            sample_qc=left_ht.sample_qc.annotate(
                **{x: left_ht[x] for x in left_ht.row if x.endswith("_residual")},
            )
        )

    logger.info("Annotating high_quality field")
    left_ht = left_ht.annotate(
        high_quality=~left_ht.sample_filters.hard_filtered
        & (hl.len(left_ht.sample_filters.qc_metrics_filters) == 0)
    )

    logger.info("Annotating releasable field")
    left_ht = left_ht.annotate(
        release=left_ht.project_meta.releasable
        & left_ht.high_quality
        & ~left_ht.project_meta.exclude
        & ~left_ht.sample_filters.release_related
    ).persist()

    logger.info(
        "Release sample counts:"
        f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.release)))}"
    )
    left_ht.describe()
    left_ht.summarize()
    logger.info(f"Final count: {left_ht.count()}")
    logger.info("Complete")

    return left_ht


def main(args):
    """Add description."""
    hl.init(log="/hail.log", default_reference="GRCh38")

    meta_ht = generate_metadata(args.regressed_metrics_outlier)
    meta_ht.checkpoint(
        meta.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite
    )
    n_pcs = meta_ht.aggregate(
        hl.agg.min(hl.len(meta_ht.population_inference.pca_scores))
    )

    meta_ht = meta_ht.annotate(
        population_inference=meta_ht.population_inference.transmute(
            **{
                f"PC{i + 1}": meta_ht.population_inference.pca_scores[i]
                for i in range(n_pcs)
            }
        ),
        hard_filters=hl.or_missing(
            hl.len(meta_ht.sample_filters.hard_filters) > 0,
            hl.delimit(meta_ht.sample_filters.hard_filters),
        ),
        qc_metrics_filters=hl.or_missing(
            hl.len(meta_ht.sample_filters.qc_metrics_filters) > 0,
            hl.delimit(meta_ht.sample_filters.qc_metrics_filters),
        ),
    )
    meta_ht.flatten().export(meta_tsv_path())


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

    main(parser.parse_args())
