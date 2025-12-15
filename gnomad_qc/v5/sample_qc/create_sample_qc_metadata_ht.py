"""Script to merge the output of all sample QC modules into a single Table."""

import argparse
import logging

import hail as hl
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SECOND_DEGREE_RELATIVES,
    SIBLINGS,
)

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v4.resources.sample_qc import hgdp_tgp_meta_updated
from gnomad_qc.v4.sample_qc.create_sample_qc_metadata_ht import (
    get_relatedness_dict_ht,
    get_relationship_filter_expr,
)
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_logging_path,
    get_samples_to_exclude,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import meta, project_meta, sample_id_collisions
from gnomad_qc.v5.resources.sample_qc import (
    finalized_outlier_filtering,
    get_gen_anc_ht,
    hard_filtered_samples,
    related_samples_to_drop,
    relatedness,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_metadata")
logger.setLevel(logging.INFO)


def restructure_meta_fields(meta_ht: hl.Table) -> hl.Table:
    """
    Restructure the meta Table to group related fields into nested structs.

    .. note::

        Creates:
        - `project_meta`: struct with project-level metadata
        - `metrics`: struct with sequencing and QC metrics

    :param meta_ht: Table with metadata.
    :return: Annotated Table with new structured fields.
    """
    # Note: 'mean_dp' in gnomAD is derived from chr20_mean_dp, whereas for AoU
    # it comes from mean_coverage.
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
            mean_dp=meta_ht.mean_depth,
        ),
    ).drop("callrate", "gen_anc")

    return meta_ht


def annotate_hard_filters(
    meta_ht: hl.Table,
    samples_to_exclude: hl.expr.SetExpression,
    aou_hard_filters_ht: hl.Table,
) -> hl.Table:
    """
    Add hard-filter annotations to the metadata Table.

    :param meta_ht: Table with metadata.
    :param samples_to_exclude: Expression with samples to exclude.
    :param aou_hard_filters_ht: Table with AoU hard filter information.
    :return: Annotated meta Table with hard-filter fields.
    """
    # Add AoU hard filters and samples to exclude.
    # Build hard-filters annotation by project.
    # Note: Verified that none of these excluded sample IDs overlap with
    # gnomAD sample IDs, so no need to call
    # 'add_project_prefix_to_sample_collisions' function.
    is_excluded = samples_to_exclude.contains(meta_ht.s)
    aou_qc_filters = aou_hard_filters_ht[meta_ht.s].sample_qc_metric_hard_filters
    has_qc_filters = hl.is_defined(aou_qc_filters) & (hl.len(aou_qc_filters) > 0)
    meta_ht = meta_ht.annotate(
        hard_filters=(
            hl.case()
            .when(
                meta_ht.project_meta.project == "gnomad",
                hl.or_else(meta_ht.hard_filters, hl.empty_set(hl.tstr)),
            )
            .when(
                meta_ht.project_meta.project == "aou",
                hl.empty_set(hl.tstr)
                .union(
                    hl.if_else(has_qc_filters, aou_qc_filters, hl.empty_set(hl.tstr))
                )
                .union(
                    hl.if_else(
                        is_excluded, hl.set(["sample_exclusion"]), hl.empty_set(hl.tstr)
                    )
                ),
            )
            .default(hl.empty_set(hl.tstr))
        )
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
    :param v4_meta_ht: gnomAD v4 metadata Table containing inferred genetic ancestry.
    :param aou_gen_anc_ht: Table with AoU projected genetic ancestry.
    :return: Table annotated with genetic ancestry inference outputs.
    """
    # Nest genetic ancestry inference fields under
    # 'genetic_ancestry_inference' and drop unnecessary fields.
    v4_meta_ht = v4_meta_ht.annotate(
        genetic_ancestry_inference=hl.struct(
            gen_anc=v4_meta_ht.population_inference.pop,
            **{
                k: v
                for k, v in v4_meta_ht.population_inference.items()
                if k not in ["pop", "prob_oth", "training_pop_all", "training_pop"]
            },
        )
    )

    # Align AoU fields with v4 fields.
    v4_fields = [f for f in v4_meta_ht.genetic_ancestry_inference.dtype.fields]
    fields_exprs = {f: aou_gen_anc_ht[f] for f in v4_fields}
    aou_gen_anc_ht = aou_gen_anc_ht.annotate(
        genetic_ancestry_inference=hl.struct(**fields_exprs)
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

    return meta_ht


def update_hgdp_tgp_outlier_annotation(
    meta_ht: hl.Table, v4_meta_ht: hl.Table
) -> hl.Table:
    """
    Update `outlier_filters` in ht.sample_filters based on the annotations in hgdp_tgp_meta_updated.

    .. note::

        For 'gnomad' samples, we cannot rely solely on the length of 'outlier_filters' because some v3-filtered samples were
        rescued in v4 and still retain populated 'outlier_filters'. The 'hgdp_tgp_meta_updated' contains an updated 'outlier' bool. For samples
    present in "hgdp_tgp_meta_updated" and marked as either "hgdp" or "tgp" in the v4 meta subsets, if the updated outlier annotation is True,
    outlier_filters will be populated with hgdp_tgp_outlier. If the updated outlier annotation is False, the outlier_filters annotation will be
    updated to an empty set.

    :param meta_ht: Table with metadata.
    :param v4_meta_ht: Table with v4 metadata.
    :return: Annotated Table with updated outlier_filters for the hgdp_tgp samples.
    """
    # Confirmed do not need to run sample collision function on
    # hgdp_tgp_meta_updated.ht().
    hg_ht = hgdp_tgp_meta_updated.ht()

    # Filter to HGDP or TGP subset samples only in the v4 meta.
    subset_ht = v4_meta_ht.filter(v4_meta_ht.subsets.hgdp | v4_meta_ht.subsets.tgp)

    # Annotate subset table with updated outlier info.
    subset_ht = subset_ht.annotate(
        updated_outlier=hg_ht[
            subset_ht.project_meta.sample_id
        ].hgdp_tgp_meta.subcontinental_pca.outlier
    )

    # Annotate main metadata with the updated outlier filters
    meta_ht = meta_ht.annotate(
        outlier_filters=hl.if_else(
            hl.is_defined(subset_ht[meta_ht.s].updated_outlier),
            hl.if_else(
                subset_ht[meta_ht.s].updated_outlier,
                meta_ht.outlier_filters.union({"hgdp_tgp_outlier"}),
                hl.empty_set(hl.tstr),
            ),
            meta_ht.outlier_filters,
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
    # Add outlier filters based on project.
    meta_ht = meta_ht.annotate(
        outlier_filters=hl.or_else(
            hl.if_else(
                meta_ht.project_meta.project == "gnomad",
                meta_ht.outlier_filters,
                outlier_filters_ht[meta_ht.key].qc_metrics_filters,
            ),
            hl.empty_set(hl.tstr),
        )
    )

    # Create sample_filters struct containing information on outlier filters
    # and hard filters.
    meta_ht = meta_ht.transmute(
        sample_filters=hl.struct(
            outlier_filters=meta_ht.outlier_filters,
            hard_filters=meta_ht.hard_filters,
        )
    )

    # Add bool flags for each filter.
    meta_ht = meta_ht.annotate(
        sample_filters=meta_ht.sample_filters.annotate(
            outlier_filtered=((hl.len(meta_ht.sample_filters.outlier_filters) > 0)),
            hard_filtered=hl.len(meta_ht.sample_filters.hard_filters) > 0,
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

    :param meta_ht: Hail Table containing sample metadata, hard-filter flag, and outlier-filter flag.
    :param relatedness_ht: Table containing relatedness inference results.
    :return: Hail Table annotated with relatedness filter fields.
    """
    # Obtain relationships from relatedness inference Table.
    relatedness_inference_ht = annotate_relationships(relatedness_ht, meta_ht)

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
            relatedness_filters=relatedness_filters_ht[
                meta_ht.key
            ].relatedness_inference.relatedness_filters,
            release_relatedness_filters=relatedness_filters_ht[
                meta_ht.key
            ].relatedness_inference.release_relatedness_filters,
            relationships=relatedness_inference_ht[meta_ht.key].relationships,
            released_gnomad_exomes_aou_duplicate=relatedness_inference_ht[
                meta_ht.key
            ].released_gnomad_exomes_aou_duplicate,
            released_gnomad_genomes_duplicate=relatedness_inference_ht[
                meta_ht.key
            ].released_gnomad_genomes_duplicate,
        )
    )

    return meta_ht


def annotate_relationships(relatedness_ht: hl.Table, meta_ht: hl.Table) -> hl.Table:
    """
    Get relatedness relationship annotations for the combined meta Table.

    Table gets the annotations 'relationships' and 'relationships_high_quality': dictionaries of all relationships (except UNRELATED) the
          sample has with other samples in the dataset. The key is the relationship
          and the value is a set of all samples with that relationship to the given
          sample.

    :param relatedness_ht: Table with relatedness information.
    :param meta_ht: Table with 'outlier_filtered' annotation nested under 'sample_filters' indicating if a
        sample was filtered during outlier detection on sample QC metrics.
    :return: Table with relationships added.
    """
    relatedness_inference_parameters = relatedness_ht.index_globals()

    logger.info("Aggregating sample relationship information...")
    # Filter to pairs passing hard filtering (all pairs in the
    # relatedness Table pass hard filtering) for 'relationships' annotation and to only
    # pairs passing both hard-filtering and outlier filtering for
    # 'relationships_high_quality' annotation.
    filter_expr = {
        "": True,
        "_high_quality": (
            ~meta_ht[relatedness_ht.i.s].sample_filters.outlier_filtered
            & ~meta_ht[relatedness_ht.j.s].sample_filters.outlier_filtered
        ),
    }
    rel_dict_ht = get_relatedness_dict_ht(relatedness_ht, filter_expr)

    # Generate duplicate samples lists.
    exome_dups = relatedness_ht.filter(
        relatedness_ht.released_gnomad_exomes_aou_duplicate
    )
    genome_dups = relatedness_ht.filter(
        relatedness_ht.released_gnomad_genomes_duplicate
    )

    exome_dup_samples = hl.literal(
        exome_dups.aggregate(
            hl.agg.collect_as_set(exome_dups.i.s).union(
                hl.agg.collect_as_set(exome_dups.j.s)
            )
        )
    )

    genome_dup_samples = hl.literal(
        genome_dups.aggregate(
            hl.agg.collect_as_set(genome_dups.i.s).union(
                hl.agg.collect_as_set(genome_dups.j.s)
            )
        )
    )

    # Use the meta HT samples as a base HT to annotate the relatedness info on
    # because it includes all samples and defines those that pass hard filters, which is what was used
    # in relatedness inference. rel_dict_ht only includes samples with relatedness info,
    # and we want to make sure all samples that went through relatedness inference have
    # empty relationship dictionaries.
    ht = meta_ht.filter(~meta_ht.sample_filters.hard_filtered)
    ht = ht.select()

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

    # Annotate duplicate flags by membership
    ht = ht.annotate(
        released_gnomad_exomes_aou_duplicate=exome_dup_samples.contains(ht.s),
        released_gnomad_genomes_duplicate=genome_dup_samples.contains(ht.s),
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
    `relatedness_filters` and `release_relatedness_filters` struct:

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
        hard-filtered.
    :param outlier_filtered_expr: Boolean Expression indicating whether the sample was
        outlier-filtered.
    :return: Table of relatedness filters.
    """
    rel_dict = {
        "related": SECOND_DEGREE_RELATIVES,
        "duplicate_or_twin": DUPLICATE_OR_TWINS,
        "parent_child": PARENT_CHILD,
        "sibling": SIBLINGS,
    }
    relationships = relationship_ht[ht.key]
    relatedness_filters_ht = ht.annotate(
        relatedness_inference=hl.struct(
            relatedness_filters=hl.struct(
                **{
                    rel: get_relationship_filter_expr(
                        hard_filtered_expr,
                        hl.is_defined(
                            related_samples_to_drop(release=False).ht()[ht.key]
                        ),
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
                        hl.is_defined(
                            related_samples_to_drop(release=True).ht()[ht.key]
                        ),
                        relationships.relationships_high_quality,
                        rel_val,
                    )
                    for rel, rel_val in rel_dict.items()
                }
            ),
        )
    )

    return relatedness_filters_ht


def main(args):
    """Merge the output of all sample QC modules into a single Table."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/create_sample_meta.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )

    try:
        check_resource_existence(
            output_step_resources={
                "sample_qc_meta": [meta().path],
            },
            overwrite=args.overwrite,
        )

        # Restructure metadata fields.
        meta_ht = restructure_meta_fields(meta_ht=project_meta.ht())

        # Add AoU hard filters to the meta Table.
        meta_ht = annotate_hard_filters(
            meta_ht=meta_ht,
            samples_to_exclude=get_samples_to_exclude(),
            aou_hard_filters_ht=hard_filtered_samples.ht(),
        )

        # Add "gnomad" prefix to gnomAD v4 meta Table.
        v4_meta_ht = v4_meta("4.0", "genomes").ht()
        v4_meta_ht = add_project_prefix_to_sample_collisions(
            t=v4_meta_ht,
            sample_collisions=sample_id_collisions.ht(),
            project="gnomad",
        )

        # Annotate genetic ancestry inference.
        meta_ht = annotate_genetic_ancestry(
            meta_ht=meta_ht,
            v4_meta_ht=v4_meta_ht,
            aou_gen_anc_ht=get_gen_anc_ht(projection_only=True).ht(),
        )

        meta_ht = update_hgdp_tgp_outlier_annotation(
            meta_ht=meta_ht, v4_meta_ht=v4_meta_ht
        )

        # Add sample outlier and filter annotations.
        meta_ht = add_sample_filter_annotations(
            meta_ht=meta_ht, outlier_filters_ht=finalized_outlier_filtering().ht()
        )

        # Add relatedness inference and filters to the metadata Table.
        meta_ht = add_relatedness_inference(
            meta_ht=meta_ht,
            relatedness_ht=relatedness().ht(),
        )

        # Move high_quality to the end of the Table.
        meta_ht = meta_ht.select(
            *[f for f in meta_ht.row if f not in ["s", "high_quality"]],
            "high_quality",
        )

        logger.info("Annotating release field...")
        # Drop 'release' and re-annotate so that it will appear at the end.
        meta_ht = meta_ht.drop("release")
        meta_ht = meta_ht.annotate(
            release=(
                meta_ht.project_meta.releasable
                & meta_ht.high_quality
                & ~meta_ht.relatedness_inference.release_relatedness_filters.related
            ),
        )

        meta_ht = meta_ht.filter(meta_ht.project_meta.data_type == "genomes")
        meta_ht = meta_ht.checkpoint(meta().path, overwrite=args.overwrite)

        logger.info("Total genome sample count: %s", meta_ht.count())

        logger.info(
            "Total genome release sample count: %s",
            meta_ht.aggregate(hl.agg.count_where(meta_ht.release)),
        )
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("create_sample_meta", environment="rwb"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
