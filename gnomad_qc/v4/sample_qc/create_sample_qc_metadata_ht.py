"""Script to merge the output of all sample QC modules into a single Table."""
import argparse
import logging
from typing import Dict, List, Optional

import hail as hl
from gnomad.assessment.validity_checks import compare_row_counts
from gnomad.sample_qc.relatedness import (
    DUPLICATE_OR_TWINS,
    PARENT_CHILD,
    SECOND_DEGREE_RELATIVES,
    SIBLINGS,
    UNRELATED,
)
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

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
    platform,
    related_samples_to_drop,
    relatedness,
    sample_chr20_mean_dp,
    sample_qc_mt_callrate,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_metadata")
logger.setLevel(logging.INFO)


# TODO: Add annotation documentation to globals
# TODO: How to handle PCs in platform and population Tables? For example in the
#  platform Table the number of PCs will be 9 because that is what was used, should we
#  modify to have all 30 PCs, or add all 30 PCs to another annotation, or only keep the
#  9 since that is all that was used?
# TODO: Keep or drop is_female from sex HT?
# TODO: Add more nearest neighbor info?
# TODO: Add trio info?
# TODO: joint that has v3 info?
def get_sex_imputation_ht() -> hl.Table:
    """
    Load and reformat sex imputation Table for annotation on the combined meta Table.

    :return: Reformatted sex imputation Table.
    """
    ht = sex.ht()
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    ht = ht.transmute(impute_sex_stats=hl.struct(**{x: ht[x] for x in impute_stats}))

    return ht


def get_hard_filter_metric_ht(ht) -> hl.Table:
    """
    Combine sample contamination, chr20 mean DP, and QC MT callrate into a single Table.

    :param ht: Input Table to add annotations to.
    :return: Table with hard filter metric annotations added.
    """
    # TODO: Add drop of `gq_thresholds` to the sample_chr20_mean_dp code.
    hard_filter_metric_hts = {
        "contamination_approximation": contamination.ht(),
        "chr20_sample_mean_dp": sample_chr20_mean_dp.ht().drop("gq_thresholds"),
        "sample_qc_mt_callrate": sample_qc_mt_callrate.ht(),
    }
    for ann in hard_filter_metric_hts:
        ht = add_annotations(
            {"ht": ht}, {ann: hard_filter_metric_hts[ann]}, ann_top_level=True
        )

    return ht


def get_sample_filter_ht(
    ht: hl.Table, relationship_ht: hl.Table, ukb_remove: List[str], hf_s: List[str]
) -> hl.Table:
    """
    Combine sample filters into a single Table to be added to the metadata Table.

    Includes hard-filters, sample QC outlier filters, and relatedness filters.

    :param ht: Input Table to add annotations to.
    :param relationship_ht: Table with relationships annotations.
    :param ukb_remove: List of UKB samples that have been removed from the v4 VDS.
    :param hf_s: List of hard-filtered samples.
    :return: Table with hard filter metric annotations added.
    """
    ht = add_annotations(
        {"ht": ht},
        {"hard_filters": get_hard_filters_ht()},
        ann_top_level=True,
        global_top_level=True,
        sample_count_match=False,
    )

    ht = add_annotations(
        {"ht": ht},
        {"outlier_detection": finalized_outlier_filtering().ht()},
        ann_top_level=True,
        ht_missing=ukb_remove,
        ann_ht_missing=hf_s,
    )

    # Checkpoint to a temp location to prevent class too large error.
    ht = ht.checkpoint(new_temp_file("quality_filters", extension="ht"), overwrite=True)

    ht = add_annotations(
        {"ht": ht},
        {
            "relatedness_filters": annotate_relatedness_filters(
                ht.select_globals(), relationship_ht
            )
        },
        ann_top_level=True,
        global_top_level=True,
    )

    return ht


def get_hard_filters_ht() -> hl.Table:
    """
    Load and reformat hard-filters Table for annotation on the combined meta Table.

    :return: Reformatted hard-filters Table.
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
            f"relationships{label}": hl.agg.group_by(
                ht[f"_relationship{label}"], hl.agg.collect_as_set(ht.pair)
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

    logger.info("Creating gnomAD v3/v4 overlap annotation Table...")
    v3_duplicate_ht = ht.filter(ht.gnomad_v3_duplicate)
    v3_duplicate_ht = v3_duplicate_ht.key_by(
        s=hl.if_else(
            v3_duplicate_ht.i.data_type == "exomes",
            v3_duplicate_ht.i.s,
            v3_duplicate_ht.j.s,
        )
    )
    v3_duplicate_ht = v3_duplicate_ht.select(
        "gnomad_v3_duplicate", "gnomad_v3_release_duplicate"
    )

    logger.info("Aggregating sample relationship information...")
    # TODO: should we add v3 relationships to the relationships set, or have a
    #  different annotation for that?
    # Use a filter to only exome-exome pairs passing hard filtering (all pairs in the
    # relatedness Table pass hard filtering) for 'relationships' annotation and to only
    # exome-exome pairs passing both hard-filtering and outlier filtering for
    # 'relationships_high_quality' annotation.
    exome_filter_expr = (ht.i.data_type == "exomes") & (ht.j.data_type == "exomes")
    filter_expr = {
        "": exome_filter_expr,
        "_high_quality": exome_filter_expr
        & ~outlier_filter_ht[ht.i.s].outlier_filtered
        & ~outlier_filter_ht[ht.j.s].outlier_filtered,
    }
    ht = get_relatedness_dict_ht(ht, filter_expr)
    ht = ht.annotate(**v3_duplicate_ht[ht.key])
    ht = ht.select_globals(**relatedness_inference_parameters)

    ht = ht.checkpoint(new_temp_file("relationships", extension="ht"), overwrite=True)

    return ht


def annotate_relatedness_filters(ht: hl.Table, relationship_ht: hl.Table) -> hl.Table:
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
    added to the input `ht` under a `release_relatedness_filters` struct::
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
    sample_filters_ht = ht.annotate(
        relatedness_filters=hl.struct(
            **{
                rel: get_relationship_filter_expr(
                    ht.hard_filtered,
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
                    ht.outlier_filtered,
                    hl.is_defined(related_samples_to_drop(release=True).ht()[ht.key]),
                    relationships.relationships_high_quality,
                    rel_val,
                )
                for rel, rel_val in rel_dict.items()
            }
        ),
    )

    return sample_filters_ht


def add_annotations(
    ht: Dict[str, hl.Table],
    ann_ht: Dict[str, hl.Table],
    ann_top_level: bool = False,
    global_top_level: bool = False,
    ht_missing: Optional[List[str]] = None,
    ann_ht_missing: Optional[List[str]] = None,
    sample_count_match: bool = True,
) -> hl.Table:
    """
    Annotate `ht` with contents of `ann_ht` and optionally check that sample counts match.

    :param ht: Table to annotate.
    :param ann_ht: Table with annotations to add to `ht`.
    :param ann_top_level: Whether to add all annotations on `ann_ht` to the top level
        of `ht` instead of grouping them under a new annotation, `label`.
    :param global_top_level: Whether to add all global annotations on `ann_ht` to the
        top level instead of grouping them under a new annotation, `label`_parameters.
    :param ht_missing: Optional list of approved samples missing from `ht`, but
        present in `ann_ht`.
    :param ann_ht_missing: Optional list of approved samples missing from `ann_ht`, but
        present in `ht`.
    :param sample_count_match: Check whether the sample counts match in the two input
        tables. Default is True.
    :return: Table with additional annotations.
    """
    if len(ht) != 1 or len(ann_ht) != 1:
        raise ValueError("Both input dicts `ht` and `ann_ht` must be of length 1!")

    ht_label, ht = list(ht.items())[0]
    ann_label, ann_ht = list(ann_ht.items())[0]
    logger.info("\n\nAnnotating with the %s Table.", ann_label)

    def _sample_check(
        ht1: hl.Table,
        ht2: hl.Table,
        approved_missing: List[str],
        ht1_label: str,
        ht2_label: str,
    ) -> None:
        """
        Report samples found in `ht1` but not in `ht2` or the `approved_missing` list.

        :param ht1: Input Table.
        :param ht2: Input Table to compare to samples in `ht1`.
        :param approved_missing: List of approved samples that are missing from `ht2`.
        :param ht1_label: Label to use as reference to `ht1` in logging message.
        :param ht2_label: Label to use as reference to `ht2` in logging message.
        :return: None.
        """
        missing = ht1.anti_join(ht2)
        if approved_missing:
            missing = missing.filter(
                ~hl.literal(set(approved_missing)).contains(missing.s)
            )
        if missing.count() != 0:
            logger.warning(
                f"The following {missing.count()} samples are found in the {ht1_label} "
                f"Table, but are not found in the {ht2_label} Table or in the approved "
                "missing list:"
            )
            missing.select().show(n=-1)
        elif approved_missing:
            logger.info(
                f"All samples found in the {ht1_label} Table, but not found in the "
                f"{ht2_label} Table, are in the approved missing list."
            )

    if sample_count_match:
        if not compare_row_counts(ht, ann_ht):
            logger.warning("Sample counts in Tables do not match!")
            _sample_check(ht, ann_ht, ann_ht_missing, ht_label, ann_label)
            _sample_check(ann_ht, ht, ht_missing, ann_label, ht_label)
        else:
            logger.info("Sample counts match.")
    else:
        logger.info("No sample count check requested.")

    ann = ann_ht[ht.key]
    global_ann = ann_ht.index_globals()
    if not ann_top_level:
        ann = {ann_label: ann_ht[ht.key]}
    if not global_top_level:
        global_ann = {f"{ann_label}_parameters": ann_ht.index_globals()}

    ht = ht.annotate(**ann)
    ht = ht.annotate_globals(**global_ann)
    return ht


def main(args):
    """Merge the output of all sample QC modules into a single Table."""
    hl.init(
        log="/sample_metadata.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    # Get list of UKB samples that should be removed.
    ukb_remove = hl.import_table(all_ukb_samples_to_remove, no_header=True).f0.collect()
    # Get list of hard filtered samples before sex imputation.
    hf_no_sex_s = hard_filtered_samples_no_sex.ht().s.collect()
    # Get list of hard filtered samples with sex imputation.
    hf_s = hard_filtered_samples.ht().s.collect()
    # Get list of v3 samples (expected in relatedness and pop).
    v3_s = joint_qc_meta.ht().s.collect()

    logger.info("Loading the VDS columns to begin creation of the meta HT.")
    vds = get_gnomad_v4_vds(remove_hard_filtered_samples=False)
    vds_sample_ht = vds.variant_data.cols().select().select_globals()
    relatedness_inference_ht = annotate_relationships(
        relatedness().ht(), finalized_outlier_filtering().ht()
    )

    # Note: 71 samples are found in the right HT, but are not found in left HT.
    #  They overlap with the UKB withheld samples indicating they were not removed
    #  when the metadata HT was created.
    ht = add_annotations(
        {"vds_sample_ht": vds_sample_ht},
        {"project_meta": project_meta.ht()},
        ann_top_level=True,
        global_top_level=True,
        ht_missing=ukb_remove,
    )
    # Note: the withdrawn UKB list was updated after the sample QC HT creation, so
    #  the sample QC HT has 5 samples more in it than the final sample list.
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {"sample_qc": get_sample_qc("bi_allelic").ht()},
        ht_missing=ukb_remove,
    )
    # TODO: Add drop of `gq_thresholds` to the platform_inference code.
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {"platform_inference": platform.ht().drop("gq_thresholds")},
        ann_ht_missing=hf_no_sex_s,
    )
    # TODO: Add drop of `is_female` to the sex_inference code.
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {"sex_imputation": get_sex_imputation_ht().drop("is_female")},
        ann_ht_missing=hf_no_sex_s,
    )
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {"hard_filter_metrics": get_hard_filter_metric_ht(vds_sample_ht)},
    )
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {"population_inference": get_pop_ht().ht()},
        ht_missing=v3_s,
        ann_ht_missing=hf_s,
    )
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {"relatedness_inference": relatedness_inference_ht},
        sample_count_match=False,
    )
    # Checkpoint to a temp location to prevent class too large error.
    ht = ht.checkpoint(new_temp_file("sample_qc_meta", extension="ht"), overwrite=True)
    ht = add_annotations(
        {"vds_sample_ht": ht},
        {
            "sample_filters": get_sample_filter_ht(
                vds_sample_ht, relatedness_inference_ht, ukb_remove, hf_s
            )
        },
    )

    logger.info("\n\nAnnotating high_quality field and releasable field.")
    high_quality_expr = (
        ~ht.sample_filters.hard_filtered & ~ht.sample_filters.outlier_filtered
    )
    ht = ht.annotate(
        high_quality=high_quality_expr,
        release=ht.project_meta.releasable
        & high_quality_expr
        & ~ht.sample_filters.release_relatedness_filters.related,
    )

    ht = ht.checkpoint(meta.path, overwrite=args.overwrite)
    logger.info(
        "Release sample count: %s", ht.aggregate(hl.agg.count_where(ht.release))
    )
    ht.describe()
    logger.info("Final sample count: %s", ht.count())


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
