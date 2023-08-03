"""Script to update the AFs of the HGDP + 1KG subset for v4."""
import argparse
import logging

import hail as hl
from gnomad.utils.annotations import (
    annotate_freq,
    merge_freq_arrays,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.release import make_freq_index_dict
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.release import (
    hgdp_tgp_subset,
    hgdp_tgp_subset_annotations,
    release_ht_path,
)
from gnomad_qc.v4.resources.release import hgdp_tgp_updated_AF
from gnomad_qc.v4.resources.sample_qc import (
    hgdp_recomputed_freemix,
    hgdp_tgp_pop_outliers,
    hgdp_tgp_populations_updated,
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

    populations_ht = hgdp_tgp_populations_updated.ht()
    ht = ht.annotate(
        hgdp_tgp_meta=ht.hgdp_tgp_meta.annotate(
            population_updated=populations_ht[ht.s].population
        )
    )

    logger.info(
        "%d samples have different population labels than last release",
        ht.aggregate(
            hl.agg.count_where(
                ht.hgdp_tgp_meta.population_updated != ht.hgdp_tgp_meta.population
            )
        ),
    )

    return ht


def _get_filtered_samples(ht: hl.Table) -> tuple[hl.Table, hl.Table]:
    """Compare the old and new sample filters.

    .. Note::
    We will use the new sample filters for the v4 release sites HT. To prepare for the two sets of samples, we will use anti-join, while taking into account the samples in the two populations `Han` & `Papuan`.

    :param ht: Table with the HGDP + 1KG subset metadata with the old and new sample filters.
    :return: Table with the old and new sample filters.
    """
    # old filters used in v3.1.2 release sites HT
    samples_in_v31_release = ht.filter(
        ~ht.gnomad_sample_filters.hard_filtered
        & ~ht.gnomad_sample_filters.release_related
        & (hl.len(ht.gnomad_sample_filters.qc_metrics_filters) == 0)
    )

    samples_in_old_pops = samples_in_v31_release.filter(
        (samples_in_v31_release.hgdp_tgp_meta.population == "Papuan")
        | (samples_in_v31_release.hgdp_tgp_meta.population == "Han")
    ).select()

    # new filters to be used in v4 release sites HT
    # TODO: Julia will double check the relatedness filters, if we should use
    # two filters or one
    samples_in_v4_release = ht.filter(
        ~ht.gnomad_sample_filters.hard_filtered_updated
        & ~ht.hgdp_tgp_meta.subcontinental_pca.outlier_updated
        & ~ht.relatedness_inference.related_updated
        & ~ht.gnomad_sample_filters.release_related
    )

    samples_in_new_pops = samples_in_v4_release.filter(
        (samples_in_v4_release.hgdp_tgp_meta.population_updated == "PapuanHighlands")
        | (samples_in_v4_release.hgdp_tgp_meta.population_updated == "PapuanSepik")
        | (samples_in_v4_release.hgdp_tgp_meta.population_updated == "Han")
        | (samples_in_v4_release.hgdp_tgp_meta.population_updated == "NorthernHan")
    ).select()

    samples_to_add = samples_in_v4_release.anti_join(samples_in_v31_release)
    logger.info(
        "%d samples will be added after the new filters",
        samples_to_add.count(),
    )

    samples_to_subtract = samples_in_v31_release.anti_join(samples_in_v4_release)
    logger.info(
        "%d sample(s) will be removed after the new filters",
        samples_to_subtract.count(),
    )

    samples_to_add = samples_to_add.select().join(samples_in_new_pops, how="outer")
    samples_to_subtract = samples_to_subtract.select().join(
        samples_in_old_pops, how="outer"
    )
    logger.info("%d samples in the `sum` set", samples_to_add.count())
    logger.info("%d samples in the `subtract` set", samples_to_subtract.count())

    return samples_to_add, samples_to_subtract


def calculate_AFs_for_selected_samples(
    mt: hl.MatrixTable, samples_ht: hl.Table
) -> hl.Table:
    """Calculate the AFs for selected samples.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with the samples to be added or subtracted.
    :return: Table with the AFs for the selected samples.
    """
    mt = mt.filter_cols(samples_ht[mt.col_key].s)

    logger.info("Generating frequency data...")
    mt = annotate_freq(
        mt,
        sex_expr=mt.meta.sex_imputation.sex_karyotype,
        pop_expr=mt.meta.hgdp_tgp_meta.population.lower(),
    )
    # make population labels lowercase to match the constants in release HT.

    mt = mt.annotate_globals(
        freq_index_dict=make_freq_index_dict(freq_meta=hl.eval(mt.freq_meta))
    )

    mt = mt.annotate_rows(freq=set_female_y_metrics_to_na_expr(mt))
    ht = mt.select_rows(mt.freq).rows()
    return ht


def update_hgdp_pop_labels(ht: hl.Table, old_pop: hl.str, new_pop: hl.str) -> hl.Table:
    """Update the population labels for HGDP samples in the release HT.

    :param ht: release HT with old population labels.
    :param old_pop: Old population label.
    :param new_pop: New population label.
    :return: release HT with updated population labels.
    """
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: x.replace(old_pop, new_pop))
        ),
        freq_index_dict=hl.dict(
            hl.zip(
                ht.freq_index_dict.keys().map(lambda k: k.replace(old_pop, new_pop)),
                ht.freq_index_dict.values(),
            )
        ),
    )

    return ht


def main(args):
    """Script to update update Allele Frequencies for HGDP + 1KG subset for v4."""
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = test_n_partitions or test_gene

    hl.init(
        log="/update_AFs_HGDP+TGP.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    mt = hgdp_tgp_subset(dense=True).mt()

    if test_gene or test_n_partitions:
        if test_gene:
            logger.info("Filtering to DRD2 in MT for testing purposes...")
            test_interval = [
                hl.parse_locus_interval(
                    "chr11:113409605-113475691", reference_genome="GRCh38"
                )
            ]

            mt = hl.filter_intervals(mt, test_interval)

        elif test_n_partitions:
            logger.info(
                "Filtering to %s partitions for testing purposes...", test_n_partitions
            )

            mt = mt._filter_partitions(range(test_n_partitions))

    logger.info("Loading HGDP_TGP subset meta HT...")
    meta_ht = hgdp_tgp_subset_annotations(sample=True).ht()

    logger.info("Adding updated sample QC annotations to meta HT...")
    meta_ht = add_updated_sample_qc_annotations(meta_ht)
    samples_to_add, samples_to_subtract = _get_filtered_samples(meta_ht)

    logger.info("Calculating AFs for selected samples...")
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    af_added_samples = calculate_AFs_for_selected_samples(mt, samples_to_add)
    af_subtracted_samples = calculate_AFs_for_selected_samples(mt, samples_to_subtract)

    logger.info("Loading the release HT...")
    ht = hl.read_table(release_ht_path(release_version="3.1.2"))

    logger.info("Updating the HGDP pop labels...")
    ht = update_hgdp_pop_labels(ht, old_pop="miaozu", new_pop="miao")
    ht = update_hgdp_pop_labels(ht, old_pop="yizu", new_pop="yi")
    ht = update_hgdp_pop_labels(ht, old_pop="bantusafrica", new_pop="bantusouthafrica")

    logger.info("Merging AFs for subtracted samples...")
    [freq, freq_meta] = merge_freq_arrays(
        [ht.freq, af_subtracted_samples.freq],
        [ht.freq_meta, af_subtracted_samples.freq_meta],
        operation="diff",
    )

    ht = ht.annotate(freq=freq)
    ht = ht.annotate_globals(freq_meta=freq_meta)

    logger.info("Merging AFs for added samples...")
    [freq, freq_meta] = merge_freq_arrays(
        [ht.freq, af_added_samples.freq], [ht.freq_meta, af_added_samples.freq_meta]
    )

    ht = ht.annotate(freq=freq)
    ht = ht.annotate_globals(freq_meta=freq_meta)

    logger.info("Writing out the AF HT...")
    if test:
        ht.write(hgdp_tgp_updated_AF(test=True).path, overwrite=True)
    else:
        ht.write(hgdp_tgp_updated_AF().path, overwrite=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script updates AFs for HGDP + 1KG subset for v4 release HT."
    )
    parser.add_argument(
        "--test_n_partitions",
        help="Test on a subset of partitions",
        type=int,
    )
    parser.add_argument(
        "--test_gene",
        help="Test on a subset of variants in DRD2 gene",
        action="store_true",
    )

    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
