"""Script to merge the output of all sample QC modules into a single Table."""
import argparse
import logging
from datetime import datetime
from functools import reduce
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
from gnomad_qc.v4.resources.meta import gatk_versions, meta, project_meta
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


# TODO: Add annotation documentation to globals. gnomad_production issue #898.
# TODO: How to handle PCs in platform and population Tables? For example in the
#  platform Table the number of PCs will be 9 because that is what was used, should we
#  modify to have all 30 PCs, or add all 30 PCs to another annotation, or only keep the
#  9 since that is all that was used? gnomad_production issue #899.
# TODO: Add more nearest neighbor info? gnomad_production issue #900.
# TODO: Add trio info? gnomad_production issue #901.
# TODO: Should we have a joint HT that has v3 info? Including adding v3 relationships
#  to the relationships set, or have different annotation for that. gnomad_production
#  issue #902.
# TODO: Add GATK version resource and an annotation for it. gnomad_production
#  issue #903.
# TODO: Add an annotation indicating a sample is a test sample like CHM.
#  gnomad_production issue #905.


def get_project_meta() -> hl.Table:
    """Load project-specific metadata Table and add GATK version."""
    fixed_homalt_ver = hl.literal({"4.1.4.1", "4.1.8.0"})
    ht = project_meta.ht()

    # Add an annotation at the project_meta level indicating the sample belongs to UKB.
    project_meta_expr = ht.project_meta.annotate(
        UKB_sample=ht.project_meta.project == "UKB"
    )
    # Add GATK version to project metadata.
    project_meta_expr = project_meta_expr.annotate(
        gatk_version=hl.or_else(
            gatk_versions.ht()[ht.key].gatk_version,
            hl.or_missing(project_meta_expr.UKB_sample, "4.0.10.1"),
        )
    )
    # Add an annotation indicating whether the sample was processed with the fixed
    # homalt model.
    project_meta_expr = project_meta_expr.annotate(
        fixed_homalt_model=fixed_homalt_ver.contains(project_meta_expr.gatk_version),
    )
    ht = ht.annotate(project_meta=project_meta_expr)

    return ht


def get_sex_imputation_ht() -> hl.Table:
    """
    Load and reformat sex imputation Table for annotation on the combined meta Table.

    :return: Reformatted sex imputation Table.
    """
    ht = sex.ht()
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    ht = ht.transmute(impute_sex_stats=hl.struct(**{x: ht[x] for x in impute_stats}))

    return ht


def get_hard_filters_ht(ht: hl.Table) -> hl.Table:
    """
    Load and reformat hard-filters Table for annotation on the combined meta Table.

    :param: Input HT to add hard-filter annotations to.
    :return: Reformatted hard-filters Table.
    """
    # Get hard filtered samples before sex imputation.
    hf_no_sex_s = hl.literal(hard_filtered_samples_no_sex.ht().s.collect())

    hf_ht = hard_filtered_samples.ht()
    hf_ht = hf_ht.annotate(hf_no_sex=hf_no_sex_s.contains(hf_ht.s))

    hf_explode_ht = hf_ht.explode(hf_ht.hard_filters)
    filters = hf_explode_ht.aggregate(
        hl.agg.group_by(
            hl.if_else(hf_explode_ht.hf_no_sex, "hf_no_sex", "hf_sex"),
            hl.agg.collect_as_set(hf_explode_ht.hard_filters),
        )
    )

    hf = hf_ht[ht.key]
    hf_expr = hf.hard_filters

    # If the sample doesn't have chimera data, it can't be filtered based on chimeras,
    # so it should have a missing value instead of False.
    no_chimera_expr = hl.is_missing(project_meta.ht()[ht.key].bam_metrics.chimeras_rate)
    ann_expr = {
        "chimera": hl.if_else(
            no_chimera_expr,
            hl.missing(hl.tbool),
            hl.or_else(hf_expr.contains("chimera"), False),
        )
    }
    for f in filters["hf_no_sex"] - {"chimera"}:
        ann_expr[f] = hl.or_else(hf_expr.contains(f), False)

    # If the sample was filtered before sex imputation, it should have a missing value
    # for the sex imputation filters.
    for f in filters["hf_sex"]:
        ann_expr[f] = hl.if_else(
            hf.hf_no_sex,
            hl.missing(hl.tbool),
            hl.or_else(hf_expr.contains(f), False),
            missing_false=True,
        )

    ann_expr["hard_filters"] = hl.or_else(
        hf_expr, hl.empty_set(hf_expr.dtype.element_type)
    )
    ann_expr["hard_filtered"] = hl.or_else(hl.len(hf_expr) > 0, False)

    ht = ht.annotate(**ann_expr)
    ht = ht.annotate_globals(**hf_ht.index_globals())
    ht = ht.checkpoint(new_temp_file("hard_filters", extension="ht"), overwrite=True)

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

    logger.info("Creating gnomAD v3/v4 overlap annotation Table...")
    v3_duplicate_ht = ht.filter(ht.gnomad_v3_duplicate)
    v3_duplicate_ht = v3_duplicate_ht.key_by(
        s=hl.if_else(
            v3_duplicate_ht.i.data_type == "exomes",
            v3_duplicate_ht.i.s,
            v3_duplicate_ht.j.s,
        )
    )

    logger.info("Aggregating sample relationship information...")
    # Filter to only exome-exome pairs passing hard filtering (all pairs in the
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
        **{
            r: hl.or_else(v3_duplicate_ht[ht.key][r], False)
            for r in ["gnomad_v3_duplicate", "gnomad_v3_release_duplicate"]
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


def add_annotations(
    base_ht: hl.Table,
    ann_ht: hl.Table,
    ann_label: str,
    ann_top_level: bool = False,
    global_top_level: bool = False,
    base_ht_missing: Optional[List[str]] = None,
    ann_ht_missing: Optional[List[str]] = None,
    sample_count_match: bool = True,
) -> hl.Table:
    """
    Annotate `base_ht` with contents of `ann_ht` and optionally check that sample counts match.

    :param base_ht: Table to annotate.
    :param ann_ht: Table with annotations to add to `base_ht`.
    :param ann_label: Label to use for Struct annotation of `ann_ht` on `base_ht` if
        `ann_top_level` is True. Also used and for logging message describing
        annotations being added.
    :param ann_top_level: Whether to add all annotations on `ann_ht` to the top level
        of `base_ht` instead of grouping them under a new annotation, `ann_label`.
    :param global_top_level: Whether to add all global annotations on `ann_ht` to the
        top level instead of grouping them under a new annotation, `ann_label`_parameters.
    :param base_ht_missing: Optional list of approved samples missing from `base_ht`,
        but present in `ann_ht`.
    :param ann_ht_missing: Optional list of approved samples missing from `ann_ht`, but
        present in `base_ht`.
    :param sample_count_match: Check whether the sample counts match in the two input
        tables. Default is True.
    :return: Table with additional annotations.
    """
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
        if not compare_row_counts(base_ht, ann_ht):
            logger.warning("Sample counts in Tables do not match!")
            _sample_check(base_ht, ann_ht, ann_ht_missing, "input", ann_label)
            _sample_check(ann_ht, base_ht, base_ht_missing, ann_label, "input")
        else:
            logger.info("Sample counts match.")
    else:
        logger.info("No sample count check requested.")

    ann = ann_ht[base_ht.key]
    global_ann = ann_ht.index_globals()
    if not ann_top_level:
        ann = {ann_label: ann_ht[base_ht.key]}
    if not global_top_level:
        global_ann = {f"{ann_label}_parameters": ann_ht.index_globals()}

    ht = base_ht.annotate(**ann)
    ht = ht.annotate_globals(**global_ann)

    return ht


def get_sample_filter_ht(base_ht: hl.Table, relationship_ht: hl.Table) -> hl.Table:
    """
    Combine sample filters into a single Table to be added to the metadata Table.

    Includes hard-filters, sample QC outlier filters, and relatedness filters.

    :param base_ht: Input Table to add annotations to.
    :param relationship_ht: Table with relationships annotations.
    :return: Table with hard filter metric annotations added.
    """
    logger.info("Combining sample filters Tables for 'sample_filters' struct.")
    # Get list of UKB samples that should be removed.
    ukb_remove = hl.import_table(all_ukb_samples_to_remove, no_header=True).f0.collect()
    # Get list of hard filtered samples with sex imputation.
    hf_s = hard_filtered_samples.ht().s.collect()
    # Get list of unreleasable samples, the outlier filtering Table only includes
    # releasable samples.
    meta_ht = project_meta.ht()
    unreleasable_s = meta_ht.filter(~meta_ht.project_meta.releasable).s.collect()

    hard_filters_ht = get_hard_filters_ht(base_ht)
    outlier_filters_ht = finalized_outlier_filtering().ht()
    relatedness_filters_ht = annotate_relatedness_filters(
        ht=base_ht,
        relationship_ht=relationship_ht,
        hard_filtered_expr=hard_filters_ht[base_ht.key].hard_filtered,
        outlier_filtered_expr=outlier_filters_ht[base_ht.key].outlier_filtered,
    )
    sample_filters = {
        "hard_filters": {
            "ann_ht": hard_filters_ht,
            "global_top_level": True,
        },
        "outlier_detection": {
            "ann_ht": outlier_filters_ht,
            "base_ht_missing": ukb_remove,
            "ann_ht_missing": hf_s + unreleasable_s,
        },
        "relatedness_filters": {
            "ann_ht": relatedness_filters_ht,
            "global_top_level": True,
        },
    }

    sample_filters = [
        {**ann_params, **{"ann_label": ann, "ann_top_level": True}}
        for ann, ann_params in sample_filters.items()
    ]
    sample_filters_ht = reduce(
        lambda ht, ann_params: add_annotations(ht, **ann_params),
        [base_ht] + sample_filters,
    )
    sample_filters_ht = sample_filters_ht.checkpoint(
        new_temp_file("sample_filters", extension="ht"), overwrite=True
    )

    return sample_filters_ht


def get_hard_filter_metric_ht(base_ht: hl.Table) -> hl.Table:
    """
    Combine sample contamination, chr20 mean DP, and QC MT callrate into a single Table.

    :param base_ht: Input Table to add annotations to.
    :return: Table with hard filter metric annotations added.
    """
    logger.info("Combining hard-filter metric Tables for 'hard_filter_metrics' struct.")

    # NOTE: Forgot to drop the `gq_thresholds` in the sample_chr20_mean_dp code.
    # NOTE: Bi-allelic sample QC was used for hard-filtering instead of the under
    # three alt alleles sample QC metrics which were used for outlier detection
    # because we realized the large sample size significantly decreases the number of
    # bi-allelic variants.
    hard_filter_metrics = {
        "contamination_approximation": contamination.ht(),
        "chr20_sample_mean_dp": sample_chr20_mean_dp.ht().drop("gq_thresholds"),
        "sample_qc_mt_callrate": sample_qc_mt_callrate.ht(),
        "bi_allelic_sample_qc": get_sample_qc("bi_allelic").ht(),
    }
    hard_filter_metrics = [
        {
            "ann_ht": ann_ht,
            "ann_label": ann,
            "ann_top_level": False if ann == "bi_allelic_sample_qc" else True,
        }
        for ann, ann_ht in hard_filter_metrics.items()
    ]
    hard_filter_metrics_ht = reduce(
        lambda ht, ann_params: add_annotations(ht, **ann_params),
        [base_ht] + hard_filter_metrics,
    )
    hard_filter_metrics_ht = hard_filter_metrics_ht.checkpoint(
        new_temp_file("hard_filter_metrics", extension="ht"), overwrite=True
    )

    return hard_filter_metrics_ht


def get_sample_qc_meta_ht(base_ht: hl.Table) -> hl.Table:
    """
    Combine all sample QC metadata Tables into a single Table.

    :param base_ht: Input Table to add all sample QC annotations to.
    :return: Table with all sample QC annotation Tables added to `base_ht`.
    """
    # Get list of UKB samples that should be removed.
    ukb_remove = hl.import_table(all_ukb_samples_to_remove, no_header=True).f0.collect()
    # Get list of hard filtered samples before sex imputation.
    hf_no_sex_s = hard_filtered_samples_no_sex.ht().s.collect()
    # Get list of hard filtered samples with sex imputation.
    hf_s = hard_filtered_samples.ht().s.collect()
    # Get list of v3 samples (expected in relatedness and pop).
    v3_s = joint_qc_meta.ht().s.collect()

    hard_filter_metrics_ht = get_hard_filter_metric_ht(base_ht)
    relatedness_inference_ht = annotate_relationships(
        relatedness().ht(), finalized_outlier_filtering().ht()
    )
    sample_filters_ht = get_sample_filter_ht(base_ht, relatedness_inference_ht)

    sample_qc_meta = {
        "project_meta": {
            "ann_ht": get_project_meta(),
            "ann_top_level": True,
            "global_top_level": True,
            # Note: 71 samples are found in the project meta HT, but are not found in
            #  the loaded VDS. They overlap with the UKB withheld samples indicating
            #  they were not removed when the project metadata HT was created.
            "base_ht_missing": ukb_remove,
        },
        "sample_qc": {
            "ann_ht": get_sample_qc("under_three_alt_alleles").ht(),
            # Note: the withdrawn UKB list was updated after the sample QC HT creation,
            #  so the sample QC HT has 5 samples more in it than the loaded VDS.
            "base_ht_missing": ukb_remove,
        },
        "platform_inference": {
            # Note: Forgot to drop `gq_thresholds` in the platform_inference code.
            "ann_ht": platform.ht().drop("gq_thresholds"),
            "ann_ht_missing": hf_no_sex_s,
        },
        "sex_imputation": {
            # Note: Forgot to drop `is_female` in the sex_inference code.
            "ann_ht": get_sex_imputation_ht().drop("is_female"),
            "ann_ht_missing": hf_no_sex_s,
        },
        "hard_filter_metrics": {"ann_ht": hard_filter_metrics_ht},
        "population_inference": {
            "ann_ht": get_pop_ht().ht(),
            "base_ht_missing": v3_s,
            "ann_ht_missing": hf_s,
        },
        "relatedness_inference": {
            "ann_ht": relatedness_inference_ht,
            "sample_count_match": False,
        },
        "sample_filters": {"ann_ht": sample_filters_ht},
    }

    sample_qc_meta = [
        {**ann_params, **{"ann_label": ann}}
        for ann, ann_params in sample_qc_meta.items()
    ]
    sample_qc_meta_ht = reduce(
        lambda ht, ann_params: add_annotations(ht, **ann_params),
        [base_ht] + sample_qc_meta,
    )

    return sample_qc_meta_ht


def main(args):
    """Merge the output of all sample QC modules into a single Table."""
    hl.init(
        log="/sample_metadata.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    logger.info("Loading the VDS columns to begin creation of the meta HT.")
    vds = get_gnomad_v4_vds(remove_hard_filtered_samples=False)
    ht = get_sample_qc_meta_ht(vds.variant_data.cols().select().select_globals())

    logger.info("\n\nAnnotating high_quality field and releasable field.")
    hq_expr = ~ht.sample_filters.hard_filtered & ~ht.sample_filters.outlier_filtered

    ht = ht.annotate(
        high_quality=hq_expr,
        release=(
            ht.project_meta.releasable
            & hq_expr
            & ~ht.sample_filters.release_relatedness_filters.related
        ),
    )
    ht = ht.annotate_globals(date=datetime.now().isoformat())
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
