"""Script to merge the output of all sample QC modules into a single Table."""

import argparse
import logging

from typing import Dict

import hail as hl

from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SECOND_DEGREE_RELATIVES,
    SIBLINGS,
    UNRELATED,
)

from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_logging_path,
    get_samples_to_exclude,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import project_meta, sample_id_collisions
from gnomad_qc.v5.resources.genetic_ancestry import (
    get_gen_anc_ht,
)
from gnomad_qc.v5.resources.sample_qc import (
    finalized_outlier_filtering,
    get_gen_anc_ht,
    hard_filtered_samples,
    related_samples_to_drop,
    relatedness,
)
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v4.resources.variant_qc import TRUTH_SAMPLES

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_metadata")
logger.setLevel(logging.INFO)


def restructure_meta_fields(meta_ht: hl.Table) -> hl.Table:
    """
    Restructure the meta Table to group related fields into nested structs.

    .. note::

        Creates:
        - `sex_imputation`: struct with sex karyotype and mean depth
        - `project_meta`: struct with project-level metadata
        - `metrics`: struct with sequencing and QC metrics

    :param meta_ht: Table with metadata.
    :return: Annotated Table with new structured fields.
    """
    meta_ht = meta_ht.transmute(
        sex_imputation=hl.struct(
            sex_karyotype=meta_ht.sex_karyotype,
            mean_dp=meta_ht.mean_depth,
        )
    )

    meta_ht = meta_ht.transmute(
        project_meta=hl.struct(
            age=meta_ht.age,
            releasable=meta_ht.releasable,
            sample_source=meta_ht.sample_source,
            data_type=meta_ht.data_type,
            project=meta_ht.project,
        ),
        metrics=hl.struct(
            contamination_freemix=meta_ht.contamination_freemix,
            contamination_charr=meta_ht.contamination_charr,
            chimeras_rate=meta_ht.chimeras_rate,
            bases_over_20x_coverage=meta_ht.bases_over_20x_coverage,
            median_insert_size=meta_ht.median_insert_size,
            callrate=meta_ht.callrate,
        ),
    )

    return meta_ht


def prepare_meta_with_hard_filters(
    meta_ht: hl.Table,
    samples_to_exclude: hl.expr.SetExpression,
    aou_hard_filters_ht: hl.Table,
) -> hl.Table:
    """
    Add hard filter annotations to the metadata Table

    :param meta_ht: Table with metadata.
    :param samples_to_exclude: Table with samples to exclude.
    :return: Annotated meta Table with hard filter fields
    """
    # Build exclusion list.
    exclusion_ht = hl.Table.parallelize(
        [{"s": x} for x in hl.eval(samples_to_exclude)]
    ).key_by("s")

    exclusion_ht = add_project_prefix_to_sample_collisions(
        t=exclusion_ht, sample_collisions=sample_id_collisions.ht(), project="aou"
    )

    samples_to_exclude = hl.literal(exclusion_ht.s.collect())

    # Obtain AoU hard filters and set hard_filter to sample_qc_metrics if any hard filters were applied.
    aou_hard_filters_ht = aou_hard_filters_ht.annotate(
        hard_filters=hl.if_else(
            hl.len(aou_hard_filters_ht.sample_qc_metric_hard_filters) > 0,
            {"sample_qc_metrics"},
            hl.empty_set(
                aou_hard_filters_ht.sample_qc_metric_hard_filters.dtype.element_type
            ),
        )
    )

    # Add sample exclusion to hard filters if sample is in exclusion list.
    aou_hard_filters_ht = aou_hard_filters_ht.annotate(
        hard_filters=hl.if_else(
            samples_to_exclude.contains(aou_hard_filters_ht.s),
            aou_hard_filters_ht.hard_filters.union({"sample_exclusion"}),
            aou_hard_filters_ht.hard_filters,
        )
    )

    # Add AoU hard filter to the meta Table.
    meta_ht = meta_ht.annotate(
        hard_filters=hl.if_else(
            meta_ht.project == "gnomad",
            meta_ht.hard_filters,
            aou_hard_filters_ht[meta_ht.key].hard_filters,
        )
    )

    meta_ht = meta_ht.annotate(
        hard_filtered=hl.or_else(hl.len(meta_ht.hard_filters) > 0, False)
    )

    return meta_ht


def annotate_genetic_ancestry(
    meta_ht: hl.Table, v4_meta_ht: hl.Table, aou_gen_anc_ht: hl.Table
) -> hl.Table:
    """
    Annotate metadata Table with genetic ancestry information.

    .. note::

        - For 'gnomad' project, uses v4 genetic ancestry inference.
        - For 'aou' project, uses AoU projected genetic ancestry.

    :param meta_ht: Table with metadata.
    :param v4_meta_ht: Table with v4 genetic ancestry inference.
    :param aou_gen_anc_ht: Table with AoU projected genetic ancestry.
    :return: Table with genetic ancestry inference annotated.
    """

    # Prefix collisions for``.
    v4_meta_ht = add_project_prefix_to_sample_collisions(
        t=v4_meta_ht,
        sample_collisions=sample_id_collisions.ht(),
        project="gnomad",
    )

    # Nest all row fields (except key) under gen`etic_ancestry_inference for AoU data.
    fields_to_nest = [f for f in aou_gen_anc_ht.row if f != "s"]
    aou_gen_anc_ht = aou_gen_anc_ht.transmute(
        genetic_ancestry_inference=hl.struct(
            **{f: aou_gen_anc_ht[f] for f in fields_to_nest}
        )
    )

    # Drop unnecessary fields from v4 meta that are not in the AoU metadata.
    v4_meta_ht = v4_meta_ht.transmute(
        genetic_ancestry_inference=v4_meta_ht.population_inference.drop(
            "prob_oth", "training_pop_all", "training_pop"
        )
    )

    # Rename 'pop' to 'gen_anc' while maintainf order of annotations.
    v4_meta_ht = v4_meta_ht.annotate(
        genetic_ancestry_inference=hl.struct(
            gen_anc=v4_meta_ht.genetic_ancestry_inference.pop,
            **{
                k: v
                for k, v in v4_meta_ht.genetic_ancestry_inference.items()
                if k != "pop"
            },
        )
    )

    # Align AoU fields with v4 field order.
    v4_fields = [f for f in v4_meta_ht.genetic_ancestry_inference.dtype.fields]
    aou_gen_anc_ht = aou_gen_anc_ht.annotate(
        genetic_ancestry_inference=hl.struct(
            **{f: aou_gen_anc_ht["genetic_ancestry_inference"][f] for f in v4_fields}
        )
    )

    # Annotate genetic ancestry inference by project.
    meta_ht = meta_ht.annotate(
        genetic_ancestry_inference=hl.case()
        .when(
            meta_ht.project_meta.project == "gnomad",
            v4_meta_ht[meta_ht.s].genetic_ancestry_inference,
        )
        .when(
            meta_ht.project_meta.project == "aou",
            aou_gen_anc_ht[meta_ht.s].genetic_ancestry_inference,
        )
        .or_missing()
    )

    meta_ht = meta_ht.annotate(
        genetic_ancestry_inference=meta_ht.genetic_ancestry_inference.annotate(
            gen_anc=meta_ht.gen_anc
        )
    )

    return meta_ht


def add_sample_filter_annotations(
    meta_ht: hl.Table, outlier_filters_ht: hl.Table
) -> hl.Table:
    """
    Annotate a metadata Hail Table with sample filter information, including outlier and hard filters, and derive a 'high_quality' flag.

    :param meta_ht: Table with metadata.
    :param outlier_filters_ht: Table with outlier filters.
    :return: Table with sample filter annotations.
    """
    # Add outlier filters based on project
    meta_ht = meta_ht.annotate(
        outlier_filters=hl.if_else(
            meta_ht.project_meta.project == "gnomad",
            meta_ht.outlier_filters,
            outlier_filters_ht[meta_ht.key].qc_metrics_filters,
        )
    )

    # Create sample_filters struct containing information on outlier filters and hard filters.
    meta_ht = meta_ht.transmute(
        sample_filters=hl.struct(
            outlier_filters=meta_ht.outlier_filters,
            hard_filters=meta_ht.hard_filters,
        )
    )

    # Add bool flags for each filter.
    meta_ht = meta_ht.annotate(
        sample_filters=meta_ht.sample_filters.annotate(
            hard_filtered=hl.len(meta_ht.sample_filters.hard_filters) > 0,
            outlier_filtered=hl.len(meta_ht.sample_filters.outlier_filters) > 0,
        )
    )

    # Define high quality flag (sample is neither hard nor outlier filtered).
    meta_ht = meta_ht.annotate(
        high_quality=~meta_ht.sample_filters.hard_filtered
        & ~meta_ht.sample_filters.outlier_filtered
    )

    return meta_ht


def add_relatedness_inference(meta_ht: hl.Table, relatedness_ht: hl.Table) -> hl.Table:
    """
    Add relationship inference and filters to metadata Table.

    :param meta_ht: Hail Table containing sample metadata, hard filter flag, and outlier filter flag.
    :param relatedness_ht: Table containing relatedness inference results.
    :return: Hail Table annotated with relatedness filter fields.
    """

    # Obtain relationships from relatedness inference Table.
    relatedness_inference_ht = annotate_relationships(
        relatedness_ht, outlier_filters_ht
    )

    # Obtain relatedness filters.
    relatedness_filters_ht = annotate_relatedness_filters(
        ht=meta_ht,
        relationship_ht=relatedness_inference_ht,
        hard_filtered_expr=meta_ht.sample_filters.hard_filtered,
        outlier_filtered_expr=meta_ht.sample_filters.outlier_filtered,
    )

    # Annotate metadata Table with relatedness inference and filters.
    meta_ht = meta_ht.annotate(
        relatedness_inference=hl.struct(
            release_relatedness_inference=relatedness_filters_ht.release_relatedness_inference,
            relationships=relatedness_inference_ht[meta_ht.key].relationships,
        )
    )

    return meta_ht


def get_relatedness_dict_ht(
    ht: hl.Table,
    filter_exprs: Dict[str, hl.expr.BooleanExpression] = None,
) -> hl.Table:
    """
    Parse relatedness Table to get every relationship (except UNRELATED) per sample.

    Return Table keyed by sample with all sample relationships in dictionary where the
    key is the relationship and the value is a set of all samples with that
    relationship to the given sample.

    :param ht: Table with inferred relationship information. Keyed by sample pair (i, j).
    :param filter_exprs: Optional dictionary of filter expressions to apply to `ht`
        before creating the 'relationships' annotations. Keyed by the postfix to add to
        'relationships' as the annotation label, and with boolean expressions as the
        values. By default, no additional filtering is applied, and a single
        'relationships' annotation is created.
    :return: Table keyed by sample (s) with all relationships annotated as a dict.
    """
    if filter_exprs is None:
        filter_exprs = {"": True}

    # Add annotations for the relationship of each item in filter_expr.
    ht = ht.annotate(
        **{
            f"_relationship{label}": hl.or_missing(expr, ht.relationship)
            for label, expr in filter_exprs.items()
        }
    )
    relationship_labels = [f"_relationship{label}" for label in filter_exprs]

    # Filter to only pairs that are related (second-degree or closer).
    ht = ht.filter(ht.relationship != UNRELATED)

    # Build a Table of relationships duplicating the info for each pair. Pair (i, j)
    # will have two rows, one for (s: i.s, pair: j.s) and one for (s: j.s, pair: i.s).
    ht = ht.select(*relationship_labels, s=ht.i.s, pair=ht.j.s).union(
        ht.select(*relationship_labels, s=ht.j.s, pair=ht.i.s)
    )

    # Group the Table by the sample name (s) and aggregate to get a relationships
    # dictionary per sample. The dictionary is built by grouping the relationship
    # annotation (becomes the key) and collecting the set of all samples (value) with
    # each relationship. Each item in filter_expr has its own annotation.
    ht = ht.group_by(ht.s).aggregate(
        **{
            f"relationships{label}": hl.agg.filter(
                hl.is_defined(ht[f"_relationship{label}"]),
                hl.agg.group_by(
                    ht[f"_relationship{label}"], hl.agg.collect_as_set(ht.pair)
                ),
            )
            for label in filter_exprs
        }
    )

    return ht


def get_relationship_filter_expr(
    hard_filtered_expr: hl.expr.BooleanExpression,
    related_drop_expr: hl.expr.BooleanExpression,
    relationship_set: hl.expr.SetExpression,
    relationship: str,
) -> hl.expr.builders.CaseBuilder:
    """
    Return case statement to populate relatedness filters in sample_filters struct.

    :param hard_filtered_expr: Boolean for whether sample was hard filtered.
    :param related_drop_expr: Boolean for whether sample was filtered due to
        relatedness.
    :param relationship_set: Set containing all possible relationship strings for
        sample.
    :param relationship: Relationship to check for. One of DUPLICATE_OR_TWINS,
        PARENT_CHILD, SIBLINGS, or SECOND_DEGREE_RELATIVES.
    :return: Case statement used to populate sample_filters related filter field.
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.missing(hl.tbool))
        .when(relationship == SECOND_DEGREE_RELATIVES, related_drop_expr)
        .when(
            hl.is_defined(relationship_set) & related_drop_expr,
            relationship_set.contains(relationship),
        )
        .default(False)
    )


def annotate_relationships(ht: hl.Table, outlier_filter_ht: hl.Table) -> hl.Table:
    """
    Get relatedness relationship annotations for the combined meta Table.

    Table has the following annotations:
        - relationships: A dictionary of all relationships (except UNRELATED) the
          sample has with other samples in the dataset. The key is the relationship
          and the value is a set of all samples with that relationship to the given
          sample.
        - gnomad_v3_duplicate: Sample is in the gnomAD v3.1 sample set that passed hard
          filtering.
        - gnomad_v3_release_duplicate: Sample is in the gnomAD v3.1 release.

    :param ht: Sample QC filter Table to add relatedness filter annotations to.
    :param outlier_filter_ht: Table with 'outlier_filtered' annotation indicating if a
        sample was filtered during outlier detection on sample QC metrics.
    :return: Table with related filters added and Table with relationship and gnomad v3
        overlap information.
    """
    relatedness_inference_parameters = ht.index_globals()

    logger.info("Aggregating sample relationship information...")
    # Filter to pairs passing hard filtering (all pairs in the
    # relatedness Table pass hard filtering) for 'relationships' annotation and to only
    # pairs passing both hard-filtering and outlier filtering for
    # 'relationships_high_quality' annotation.
    filter_expr = {
        "": True,
        "_high_quality": (
            ~outlier_filter_ht[ht.i.s].outlier_filtered
            & ~outlier_filter_ht[ht.j.s].outlier_filtered
        ),
    }
    rel_dict_ht = get_relatedness_dict_ht(ht, filter_expr)

    # Use the outlier HT samples as a base HT to annotate the relatedness info on
    # because it includes all samples that pass hard filters, which is what was used
    # in relatedness inference. rel_dict_ht only includes samples with relatedness info,
    # and we want to make sure all samples that went through relatedness inference have
    # empty relationship dictionaries, and are False for v3 duplicate bools.
    ht = outlier_filter_ht.select()
    ht = ht.annotate(
        **{
            r: hl.bind(
                lambda x: hl.or_else(
                    x, hl.empty_dict(x.dtype.key_type, x.dtype.value_type)
                ),
                rel_dict_ht[ht.key][r],
            )
            for r in rel_dict_ht.row_value
        },
    )
    ht = ht.select_globals(**relatedness_inference_parameters)

    return ht


def annotate_relatedness_filters(
    ht: hl.Table,
    relationship_ht: hl.Table,
    hard_filtered_expr: hl.expr.BooleanExpression,
    outlier_filtered_expr: hl.expr.BooleanExpression,
) -> hl.Table:
    """
    Get relatedness filtering Table for the combined meta Table.

    Add the following related filter boolean annotations to the input `ht` under a
    `relatedness_filters` struct:

        - related: Whether the sample was filtered for second-degree (or closer)
          relatedness in the ancestry inference PCA.
        - duplicate_or_twin: Whether the filtered sample has a duplicate or twin among
          all samples that are not hard-filtered.
        - parent_child: Whether the filtered sample has a parent or child among all
          samples that are not hard-filtered.
        - sibling: Whether the filtered sample has a sibling among all samples that are
          not hard-filtered.
        - Any sample in `ht` that is hard-filtered will have a missing value for these
          annotations.

    These related filter annotations are also provided for release filtered samples
    added to the input `ht` under a `release_relatedness_filters` struct:

        - related: Whether the release filtered sample was filtered for
          second-degree (or closer) relatedness in the final release.
        - duplicate_or_twin: Whether the release filtered sample has a
          duplicate or twin among all samples that are not hard-filtered or
          outlier-filtered.
        - parent_child: Whether the release filtered sample has a parent or
          child among all samples that are not hard-filtered or outlier-filtered.
        - sibling: Whether the release filtered sample has a sibling among all
          samples that are not hard-filtered or outlier-filtered.
        - Any sample in `ht` that is hard-filtered or outlier-filtered will have a
          missing value for these annotations.

    :param ht: Sample QC filter Table to add relatedness filter annotations to.
    :param relationship_ht: Table with relationships annotations.
    :param hard_filtered_expr: Boolean Expression indicating whether the sample was
        hard filtered.
    :param outlier_filtered_expr: Boolean Expression indicating whether the sample was
        outlier filtered.
    :return: Table with related filters added and Table with relationship and gnomad v3
        overlap information.
    """
    rel_dict = {
        "related": SECOND_DEGREE_RELATIVES,
        "duplicate_or_twin": DUPLICATE_OR_TWINS,
        "parent_child": PARENT_CHILD,
        "sibling": SIBLINGS,
    }
    relationships = relationship_ht[ht.key]
    relatedness_filters_ht = ht.annotate(
        relatedness_filters=hl.struct(
            **{
                rel: get_relationship_filter_expr(
                    hard_filtered_expr,
                    hl.is_defined(related_samples_to_drop(release=False).ht()[ht.key]),
                    relationships.relationships,
                    rel_val,
                )
                for rel, rel_val in rel_dict.items()
            }
        ),
        release_relatedness_filters=hl.struct(
            **{
                rel: get_relationship_filter_expr(
                    hard_filtered_expr | outlier_filtered_expr,
                    hl.is_defined(related_samples_to_drop(release=True).ht()[ht.key]),
                    relationships.relationships_high_quality,
                    rel_val,
                )
                for rel, rel_val in rel_dict.items()
            }
        ),
    )

    return relatedness_filters_ht


def main(args):
    """Merge the output of all sample QC modules into a single Table."""
    hl.init(
        spark_conf={"spark.memory.offHeap.enabled": "false"},
        log="/home/jupyter/workspaces/gnomadproduction/create_sample_meta.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )

    try:  # Add AoU hard filters to the meta Table.
        meta_ht = prepare_meta_with_hard_filters(
            meta_ht=project_meta.ht(),
            samples_to_exclude=get_samples_to_exclude(),
            aou_hard_filters_ht=hard_filtered_samples.ht(),
        )
        # Restructure metadata fields.
        meta_ht = restructure_meta_fields(meta_ht=meta_ht)

        # Annotate genetic ancestry inference.
        meta_ht = annotate_genetic_ancestry(
            meta_ht=meta_ht,
            v4_meta_ht=v4_meta("4.0", "genomes").ht(),
            aou_gen_anc_ht=get_gen_anc_ht(projection_only=True).ht(),
        )
        # Add sample filter annotations.
        meta_ht = add_sample_filter_annotations(
            meta_ht=meta_ht, outlier_filters_ht=finalized_outlier_filtering().ht()
        )

        logger.info("\n\nAnnotating high_quality field and releasable field.")
        # Excluding samples in the ELGH2 project from the high_quality and release
        # samples because we identified that they do not have the full set of 'AS'
        # annotations in 'gvcf_info' so we need to exclude them from variant QC and release.
        hq_expr = (
            ~ht.sample_filters.hard_filtered
            & ~ht.sample_filters.outlier_filtered
            & ~ht.sample_filters.elgh2_project
        )

        ht = ht.annotate(
            high_quality=hq_expr,
            release=(
                ht.project_meta.releasable
                & hq_expr
                & ~ht.sample_filters.release_relatedness_filters.related
                & ~ht.sample_filters.control
            ),
        )

        # Add relatedness inference and filters to the metadata Table.
        meta_ht = add_relatedness_inference(
            meta_ht=meta_ht, relatedness_ht=relatedness().ht()
        )

        # ht = ht.checkpoint(meta().path, overwrite=args.overwrite)

        logger.info("Total sample count: %s", ht.count())

        # TODO: Add just count genomes
        logger.info(
            "Release sample count: %s", ht.aggregate(hl.agg.count_where(ht.release))
        )
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("create_sample_meta", environment="rwb"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
