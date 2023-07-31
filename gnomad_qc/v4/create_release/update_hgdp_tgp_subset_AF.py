"""Script to update the AFs of the HGDP + 1KG subset for v4."""
import argparse
import logging

import hail as hl

from gnomad_qc.v4.resources.sample_qc import (
    hgdp_recomputed_freemix,
    hgdp_tgp_meta,
    hgdp_tgp_pop_outliers,
    hgdp_tgp_related_samples_to_drop,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("update HGDP + 1KG subset")
logger.setLevel(logging.INFO)


def add_updated_sample_qc_annotations(ht: hl.Table) -> hl.Table:
    """Add updated sample QC annotations to the HGDP + 1KG subset.

    .. Note::
    The following annotations are updated based on the latest sample QC:
    - `sample_filters.hard_filtered_updated` to apply the recomputed freemix filter for HGDP samples;
    - `sample_filters.relatedness_inference.related_updated` to apply the updated relatedness inference implemented by Alicia Martin's group;
    - `sample_filters.pop_outlier_updated` to apply the updated pop outlier filter implemented by Alicia Martin's group.

    :param ht: Table with the HGDP + 1KG subset metadata from the last release.
    :return: Table with updated sample QC annotations.
    """
    ht = hgdp_tgp_meta.ht()

    contamination_ht = hgdp_recomputed_freemix.ht()
    samples_contaminated = contamination_ht.filter(
        contamination_ht.recomputed_contam_rate > 0.05
    ).s.collect(_localize=False)
    samples_hard_filtered = ht.filter(ht.gnomad_sample_filters.hard_filtered).s.collect(
        _localize=False
    )

    ht = ht.annotate(
        gnomad_sample_filters=ht.gnomad_sample_filters.annotate(
            hard_filtered_updated=(
                hl.case()
                .when(samples_contaminated.contains(ht.s), True)
                .when(samples_hard_filtered.contains(ht.s), True)
                .default(False)
            )
        )
    )
    num_hard_filtered_old = ht.filter(ht.gnomad_sample_filters.hard_filtered).count()
    num_hard_filtered_new = ht.filter(
        ht.gnomad_sample_filters.hard_filtered_updated
    ).count()

    logger.info("%d sample(s) hard filtered", num_hard_filtered_new)
    logger.info(
        "%d more sample(s) hard filtered than last release",
        num_hard_filtered_new - num_hard_filtered_old,
    )

    relatedness_ht = hgdp_tgp_related_samples_to_drop.ht()
    samples_related = relatedness_ht.node.s.collect(_localize=False)
    ht = ht.annotate(
        relatedness_inference=ht.relatedness_inference.annotate(
            related_updated=hl.if_else(samples_related.contains(ht.s), True, False)
        )
    )

    logger.info(
        "%d sample(s) related",
        ht.aggregate(hl.agg.count_where(ht.relatedness_inference.related_updated)),
    )
    logger.info(
        "%d sample(s) were assigned differently for the relatedness than last release",
        ht.aggregate(
            hl.agg.count_where(
                ht.relatedness_inference.related_updated
                != ht.relatedness_inference.related
            )
        ),
    )

    pop_outliers_ht = hgdp_tgp_pop_outliers.ht()
    samples_pop_outliers = pop_outliers_ht.s.collect(_localize=False)
    ht = ht.annotate(
        hgdp_tgp_meta=ht.hgdp_tgp_meta.annotate(
            subcontinental_pca=ht.hgdp_tgp_meta.subcontinental_pca.annotate(
                outlier_updated=hl.if_else(
                    samples_pop_outliers.contains(ht.s), True, False
                )
            )
        )
    )

    logger.info(
        "%d sample(s) are pop outliers",
        ht.aggregate(
            hl.agg.count_where(ht.hgdp_tgp_meta.subcontinental_pca.outlier_updated)
        ),
    )
    logger.info(
        "%d sample(s) were assigned differently than last release",
        ht.aggregate(
            hl.agg.count_where(
                ht.hgdp_tgp_meta.subcontinental_pca.outlier_updated
                != ht.hgdp_tgp_meta.subcontinental_pca.outlier
            )
        ),
    )

    return ht


def compare_old_new_filters(ht: hl.Table) -> None:
    """Compare the old and new sample filters.

    :param ht: Table with the HGDP + 1KG subset metadata with the old and new sample filters.
    :return: Table with the old and new sample filters.
    """
    # old filters used in v3.1.2 release sites HT
    old_filters = ht.filter(
        ~ht.gnomad_sample_filters.hard_filtered
        & ~ht.gnomad_sample_filters.release_related
        & (hl.len(ht.gnomad_sample_filters.qc_metrics_filters) == 0)
    ).select()

    # new filters to be used in v4 release sites HT
    new_filters = ht.filter(
        ~ht.gnomad_sample_filters.hard_filtered_updated
        & ~ht.hgdp_tgp_meta.subcontinental_pca.outlier_updated
        & ~ht.relatedness_inference.related_updated
        & ~ht.gnomad_sample_filters.release_related
    ).select()

    logger.info(
        "%d samples were kept after the old filters and %d were kept after the new"
        " filters",
        old_filters.count(),
        new_filters.count(),
    )

    logger.info(
        "%d sample(s) were kept by both filters",
        old_filters.semi_join(new_filters).count(),
    )
    logger.info(
        "%d sample(s) were kept by the old filter but not the new filter",
        old_filters.anti_join(new_filters).count(),
    )
    logger.info(
        "%d sample(s) were kept by the new filter but not the old filter",
        new_filters.anti_join(old_filters).count(),
    )
