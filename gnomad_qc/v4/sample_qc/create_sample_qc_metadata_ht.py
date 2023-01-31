"""Script to merge the output of all sample QC modules into a single Table."""
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


# TODO: Add documentation to globals
def get_sex_imputation_ht() -> hl.Table:
    """
    Load and reformat sex imputation Table for annotation on the combined meta Table.

    :return: Reformatted sex imputation Table.
    """
    ht = sex.ht()
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    ht = ht.transmute(impute_sex_stats=hl.struct(**{x: ht[x] for x in impute_stats}))

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


def get_relatedness_ht(ht, method: str = "cuking") -> hl.Table:
    """
    Load and reformat relatedness Table for annotation on the combined meta Table.

    :param ht:
    :param method:
    :return: Reformatted relatedness Table.
    """
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

    :param relatedness_ht: Table with inferred relationship information output by pc_relate.
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationships annotated as a set.
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

    :param hard_filtered_expr: Boolean for whether sample was hard filtered.
    :param relationship: Relationship to check for. One of DUPLICATE_OR_TWINS, PARENT_CHILD, or SIBLINGS.
    :param relationship_set: Set containing all possible relationship strings for sample.
    :return: Case statement used to population sample_filters related filter field.
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.null(hl.tbool))
        .when(hl.is_defined(relationship_set), relationship_set.contains(relationship))
        .default(False)
    )


def name_me(
    ht: hl.Table,
    ann_ht: hl.Table,
    label: str,
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
    :param label: Label to use for new annotation, global annotation prefix, and log
        output. `label` is modified to lowercase and spaces are replaced by underscores
        after printing logger and before use as annotation label.
    :param ann_top_level: Whether to add all annotations on `ann_ht` to the top level
        of `ht` instead of grouping them under a new annotation, `label`.
    :param global_top_level: Whether to add all global annotations on `ann_ht` to the
        top level instead of grouping them under a new annotation, `label`_parameters.
    :param ht_missing:
    :param ann_ht_missing:
    :param sample_count_match: Check whether the sample counts match in the two input
        tables. Default is True.
    :return: Table with additional annotations.
    """
    logger.info("\n\nAnnotating with the %s Table.", label)
    label = label.lower().replace(" ", "_")

    if sample_count_match:
        if not compare_row_counts(ht, ann_ht):
            logger.warning("Sample counts in left and right tables do not match!")

            in_left_not_right = ht.anti_join(ann_ht)
            if ann_ht_missing:
                in_left_not_right = in_left_not_right.filter(
                    hl.literal(set(ann_ht_missing)).contains(in_left_not_right.s),
                    keep=False,
                )
            if in_left_not_right.count() != 0:
                logger.warning(
                    f"The following {in_left_not_right.count()} samples are found in "
                    "the left HT, but are not found in the right HT or in the "
                    "approved missing list:"
                )
                in_left_not_right.select().show(n=-1)
            elif ann_ht_missing:
                logger.info(
                    "All samples found in the left HT, but not found in the right HT, "
                    "are in the approved missing list."
                )

            in_right_not_left = ht.anti_join(ann_ht)
            if ht_missing:
                in_right_not_left = in_right_not_left.filter(
                    hl.literal(set(ht_missing)).contains(in_right_not_left.s),
                    keep=False,
                )
            if in_right_not_left.count() != 0:
                logger.warning(
                    f"The following {in_right_not_left.count()} samples are found in "
                    "the right HT, but are not found in left HT or in the approved "
                    "missing list:"
                )
                in_right_not_left.select().show(n=-1)
            elif ht_missing:
                logger.info(
                    "All samples found in the right HT, but not found in the left HT, "
                    "are in the approved missing list."
                )
        else:
            logger.info("Sample counts match.")
    else:
        logger.info("No sample count check requested.")

    if ann_top_level:
        ht = ht.annotate(**ann_ht[ht.key])
    else:
        ht = ht.annotate(**{label: ann_ht[ht.key]})

    if global_top_level:
        ht = ht.annotate_globals(**ann_ht.index_globals())
    else:
        ht = ht.annotate_globals(**{f"{label}_parameter": ann_ht.index_globals()})

    return ht


def main(args):
    """Merge the output of all sample QC modules into a single Table."""
    hl.init(
        log="/sample_metadata.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

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
    ht = (
        get_gnomad_v4_vds(remove_hard_filtered_samples=False)
        .variant_data.cols()
        .select()
        .select_globals()
    )

    # Note: 71 samples are found in the right HT, but are not found in left HT.
    #  They overlap with the UKB withheld samples indicating they were not removed
    #  when the metadata HT was created. Is this expected?
    ht = name_me(ht, project_meta.ht(), "project meta", True, True, removed_ukb_samples)

    # Note: the withdrawn UKB list was updated after the sample QC HT creation, so
    #  the sample QC HT has 5 samples more in it than the final sample list.
    ht = name_me(
        ht,
        get_sample_qc("bi_allelic").ht(),
        "sample QC",
        ht_missing=removed_ukb_samples,
    )

    # TODO: Number of PCs will be 9 because that is what was used, should we modify to
    #  have all 30 PCs, or add all 30 PCs to another annotation, or only keep the 9
    #  since that is all that was used?
    ht = name_me(
        ht,
        platform.ht().drop("gq_thresholds"),
        "platform inference",
        ann_ht_missing=hf_samples_no_sex,
    )

    # TODO: Keep or drop is_female?
    # TODO: Should it be a struct with raw and adj or id this OK?
    #  chrx_frac_hom_alt: float64,
    #  chrx_frac_hom_alt_adj: float64,
    ht = name_me(
        ht, get_sex_imputation_ht(), "sex imputation", ht_missing=hf_samples_no_sex
    )

    hard_filter_metric_ht = name_me(
        ht.select().select_globals(),
        contamination.ht(),
        "contamination approximation",
        no_ann_label=True,
    )
    hard_filter_metric_ht = name_me(
        hard_filter_metric_ht,
        sample_chr20_mean_dp.ht().drop("gq_thresholds"),
        "chr20 sample mean DP",
        no_ann_label=True,
    )

    # TODO: Should it be a struct with raw and adj or id this OK
    #  sample_qc_mt_callrate: float64,
    #  sample_qc_mt_callrate_adj: float64
    hard_filter_metric_ht = name_me(
        hard_filter_metric_ht,
        sample_qc_mt_callrate.ht(),
        "sample QC MT callrate",
        no_ann_label=True,
    )

    ht = ht.annotate(hard_filter_metrics=hard_filter_metric_ht[ht.key])
    ht = ht.annotate_globals(
        hard_filter_parameters=hard_filter_metric_ht.index_globals()
    )

    # Adding checkpoint to prevent class too large error
    ht = ht.checkpoint(new_temp_file("sample_qc_meta", extension="ht"), overwrite=True)

    # TODO: How to handle PCs, different number in tables than used?
    # TODO: How to specify pop PC table to use?
    ht = name_me(
        ht,
        # get_pop_ht().ht(),
        hl.read_table(
            "gs://gnomad/v4.0/sample_qc/joint/ancestry_inference/gnomad.joint.v4.0"
            ".hgdp_tgp_training.pop.rf_w_16pcs_0.75_pop_prob.ht"
        ),
        "population inference",
        ht_missing=v3_samples,
        ann_ht_missing=hf_samples,
    )

    # TODO: rerun to get the callrate cutoff global annotation on the Table. Confirm
    #  results are identical otherwise
    sample_filters_ht = name_me(
        ht.select().select_globals(),
        get_hard_filters_ht(),
        "hard filters",
        no_ann_label=True,
        no_global_label=True,
        sample_count_match=False,
    )

    # TODO: Add relationship expression code.
    # TODO: Where to compute the release related samples to drop?
    # logger.info(logging_statement.format("PCA related samples to drop HT"))
    # ann_ht = reformat_relatedness_ht(ht)
    # ann_expr = ann_expr.annotate(**name_me(ht, ann_ht))
    # global_expr = global_expr.annotate(**ann_ht.index_globals())

    sample_filters_ht = name_me(
        sample_filters_ht,
        finalized_outlier_filtering().ht(),
        "outlier detection",
        no_ann_label=True,
        ht_missing=hf_samples,
        ann_ht_missing=removed_ukb_samples,
    )
    ht = ht.annotate(sample_filters=sample_filters_ht[ht.key])
    ht = ht.annotate_globals(**sample_filters_ht.index_globals())

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

    ht = ht.checkpoint(meta.path, overwrite=args.overwrite)

    logger.info(
        "Release sample count: %s", ht.aggregate(hl.agg.count_where(ht.release))
    )
    ht.describe()
    # ht.summarize()
    logger.info("Final sample count: %s", ht.count())

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
