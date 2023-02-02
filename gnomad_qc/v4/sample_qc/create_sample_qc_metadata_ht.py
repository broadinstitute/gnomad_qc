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


# TODO: Add annotation documentation to globals
# TODO: How to handle PCs in platform and population Tables? For example in the
#  platform Table the number of PCs will be 9 because that is what was used, should we
#  modify to have all 30 PCs, or add all 30 PCs to another annotation, or only keep the
#  9 since that is all that was used?
# TODO: Keep or drop is_female from sex HT?
# TODO: Add relatedness relationship expression code.
# TODO: Modify for any changes to get_pop_ht, ideally we change so that there is a
#  single final pop HT that keeps track of the parameters used.
# TODO: Where to compute the release related samples to drop?
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
    hard_filter_metric_hts = {
        "contamination approximation": contamination.ht(),
        "chr20 sample mean DP": sample_chr20_mean_dp.ht().drop("gq_thresholds"),
        "sample QC MT callrate": sample_qc_mt_callrate.ht(),
    }
    for label, ann_ht in hard_filter_metric_hts.items():
        ht = add_annotations(ht, ann_ht, label, ann_top_level=True)

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


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Parse relatedness Table to get every relationship (except UNRELATED) per sample.

    Return Table keyed by sample with all sample relationships in a set.

    :param relatedness_ht: Table with inferred relationship information. Keyed by
        sample pair (i, j).
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
    :param relationship: Relationship to check for. One of DUPLICATE_OR_TWINS,
        PARENT_CHILD, or SIBLINGS.
    :param relationship_set: Set containing all possible relationship strings for
        sample.
    :return: Case statement used to population sample_filters related filter field.
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.null(hl.tbool))
        .when(hl.is_defined(relationship_set), relationship_set.contains(relationship))
        .default(False)
    )


def annotate_relatedness(ht, method: str = "cuking") -> hl.Table:
    """
    Annotate Table with relatedness information for the combined meta Table.

    Create related filter boolean annotations. Where any sample that is hard filtered
    will have a missing value for these annotations. Any sample filtered for
    second-degree relatedness will have True for the 'related' annotation.

    These related filter annotations are provided for all samples and release only
    samples.

    :param method: Relatedness method to load for the `related` annotation.
    :return: Reformatted relatedness Table.
    """
    # TODO: We don't have a 'relationship' annotation on this Table yet so this might
    #  change, but this is close to what was in v3.1 and will be similar.
    rel_ht = get_relatedness_set_ht(relatedness(method).ht())
    ht = ht.select(relatedness_inference_relationships=rel_ht[ht.key].relationships)
    rel_dict = {
        "related": hl.null(hl.tbool),
        "duplicate_or_twin": DUPLICATE_OR_TWINS,
        "parent_child": PARENT_CHILD,
        "sibling": SIBLINGS,
    }
    rel_set_expr_dict = {
        "release": hl.is_defined(release_related_samples_to_drop.ht()[ht.key]),
        "all_samples": hl.is_defined(pca_related_samples_to_drop().ht()[ht.key]),
    }
    ht = ht.annotate(
        **{
            f"{k}_{rel}": get_relationship_filter_expr(
                ht.hard_filtered,
                rel_val,
                v,
            )
            for rel, rel_val in rel_dict.items()
            for k, v in rel_set_expr_dict.items()
        }
    )

    logger.info("Adding relatedness globals (cutoffs)")
    ht = ht.annotate_globals(relatedness_inference_cutoffs=rel_ht.index_globals())

    return ht


def add_annotations(
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
    :param ht_missing: Optional list of approved samples missing from `ht`, but
        present in `ann_ht`.
    :param ann_ht_missing: Optional list of approved samples missing from `ann_ht`, but
        present in `ht`.
    :param sample_count_match: Check whether the sample counts match in the two input
        tables. Default is True.
    :return: Table with additional annotations.
    """
    logger.info("\n\nAnnotating with the %s Table.", label)
    label = label.lower().replace(" ", "_")

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
            _sample_check(ht, ann_ht, ann_ht_missing, "'ht'", "'ann_ht'")
            _sample_check(ann_ht, ht, ht_missing, "'ann_ht'", "'ht'")
        else:
            logger.info("Sample counts match.")
    else:
        logger.info("No sample count check requested.")

    ann = ann_ht[ht.key]
    global_ann = ann_ht.index_globals()
    if not ann_top_level:
        ann = {label: ann_ht[ht.key]}
    if not global_top_level:
        global_ann = {f"{label}_parameter": ann_ht.index_globals()}

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

    # Note: 71 samples are found in the right HT, but are not found in left HT.
    #  They overlap with the UKB withheld samples indicating they were not removed
    #  when the metadata HT was created.
    ht = add_annotations(
        vds_sample_ht,
        project_meta.ht(),
        "project meta",
        ann_top_level=True,
        global_top_level=True,
        ht_missing=ukb_remove,
    )

    # Note: the withdrawn UKB list was updated after the sample QC HT creation, so
    #  the sample QC HT has 5 samples more in it than the final sample list.
    ann_ht = get_sample_qc("bi_allelic").ht()
    ht = add_annotations(ht, ann_ht, "sample QC", ht_missing=ukb_remove)

    ann_ht = platform.ht().drop("gq_thresholds")
    ht = add_annotations(ht, ann_ht, "platform inference", ann_ht_missing=hf_no_sex_s)

    ann_ht = get_sex_imputation_ht()
    ht = add_annotations(ht, ann_ht, "sex imputation", ann_ht_missing=hf_no_sex_s)

    ann_ht = get_hard_filter_metric_ht(vds_sample_ht)
    ht = add_annotations(ht, ann_ht, "hard filter metrics")

    # Checkpoint to a temp location to prevent class too large error.
    ht = ht.checkpoint(new_temp_file("sample_qc_meta", extension="ht"), overwrite=True)

    ann_ht = get_pop_ht(only_train_on_hgdp_tgp=args.use_pop_only_train_on_hgdp_tgp).ht()
    ht = add_annotations(
        ht, ann_ht, "population inference", ht_missing=v3_s, ann_ht_missing=hf_s
    )

    sample_filters_ht = add_annotations(
        vds_sample_ht,
        get_hard_filters_ht(),
        "hard filters",
        ann_top_level=True,
        global_top_level=True,
        sample_count_match=False,
    )

    sample_filters_ht = annotate_relatedness(sample_filters_ht, args.relatedness_method)

    sample_filters_ht = add_annotations(
        sample_filters_ht,
        finalized_outlier_filtering().ht(),
        "outlier detection",
        ann_top_level=True,
        ht_missing=hf_s,
        ann_ht_missing=ukb_remove,
    )
    ht = ht.annotate(sample_filters=sample_filters_ht[ht.key])
    ht = ht.transmute(
        relatedness_inference_relationships=ht.sample_filters.relatedness_inference_relationships
    )
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
        & ~ht.sample_filters.release_related
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
        "--relatedness-method",
        help=(
            "Which relatedness method to use when determining related samples. "
            "Options are 'cuking' and 'pc_relate'. Default is 'cuking'."
        ),
        default="cuking",
        type=str,
        choices=["cuking", "pc_relate"],
    )
    parser.add_argument(
        "--use-pop-only-train-on-hgdp-tgp",
        help=(
            "Whether to use the population inference table created from an RF "
            "classifier trained using only the HGDP and 1KG (tgp) populations."
        ),
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
