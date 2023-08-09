"""Script to update the AFs of the HGDP + 1KG subset for v4 release HT."""
import argparse
import logging
from typing import List, Optional

import hail as hl
from gnomad.utils.annotations import (
    annotate_freq,
    merge_freq_arrays,
    missing_callstats_expr,
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
    hgdp_tgp_meta_updated,
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
        "%d sample(s) were assigned differently as subcontinental outlier than last"
        " release",
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
            population_updated=populations_ht[ht.s].population,
            latitude_updated=populations_ht[ht.s].latitude,
            longitude_updated=populations_ht[ht.s].longitude,
        ),
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
    # TODO: Julia will double check the relatedness filters
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

    samples_to_add = samples_to_add.select_globals()
    samples_to_subtract = samples_to_subtract.select_globals()

    samples_to_add = samples_to_add.select().join(samples_in_new_pops, how="outer")
    samples_to_subtract = samples_to_subtract.select().join(
        samples_in_old_pops, how="outer"
    )
    logger.info(
        "%d samples in the `sum` set, including the splitted Papuan & Han",
        samples_to_add.count(),
    )
    logger.info(
        "%d samples in the `subtract` set, including the unsplitted Papuan & Han",
        samples_to_subtract.count(),
    )
    return samples_to_add, samples_to_subtract


def calculate_AFs_for_selected_samples(
    mt: hl.MatrixTable,
    samples_ht: hl.Table,
    subsets: Optional[List[str]] = None,
    pop_old: hl.bool = False,
) -> hl.Table:
    """Calculate the AFs for selected samples.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with the samples to be added or subtracted.
    :param pop_old: If True, use the old population labels. Otherwise, use the updated population labels.
    :param subset: Subsets to be used for the AFs calculation: `hgdp` or `tgp`, or both.
    :return: Table with the AFs for the selected samples.
    """
    mt = mt.filter_cols(hl.is_defined(samples_ht[mt.col_key]))

    if pop_old:
        mt = annotate_freq(
            mt,
            sex_expr=mt.gnomad_sex_imputation.sex_karyotype,
            pop_expr=mt.hgdp_tgp_meta.population.lower(),
        )
    else:
        mt = annotate_freq(
            mt,
            sex_expr=mt.gnomad_sex_imputation.sex_karyotype,
            pop_expr=mt.hgdp_tgp_meta.population_updated.lower(),
        )
    # make population labels lowercase to match the constants in release HT.

    if subsets:
        freq_meta = [
            {**x, **{"subset": "|".join(subsets)}} for x in hl.eval(mt.freq_meta)
        ]
        mt = mt.annotate_globals(freq_meta=freq_meta)
        mt = mt.annotate_globals(
            freq_index_dict=make_freq_index_dict(
                freq_meta=freq_meta,
                subsets=["|".join(subsets)],
            )
        )
    else:
        mt = mt.annotate_globals(
            freq_index_dict=make_freq_index_dict(freq_meta=hl.eval(mt.freq_meta))
        )

    mt = mt.annotate_rows(freq=set_female_y_metrics_to_na_expr(mt))
    ht = mt.rows()
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


def _concatenate_subset_frequencies(
    global_freq_ht: hl.Table, subset_freq_hts: List[hl.Table]
) -> hl.Table:
    """Concatenate subset frequencies into a single Table.

    .. note::
    This is the same as we did gnomad_qc/v3/create_release/create_release_sites_ht.py#L291-L310, except excluding the pops & downsampling options.

    :param global_freq_ht: Global frequency HT.
    :param subset_freq_hts: List of subset frequency HTs.
    :return: Concatenated frequency HT.
    """
    freq_ht = hl.Table.multi_way_zip_join(
        [global_freq_ht.select("freq").select_globals("freq_meta")] + subset_freq_hts,
        data_field_name="freq",
        global_field_name="freq_meta",
    )
    freq_ht = freq_ht.transmute(freq=freq_ht.freq.flatmap(lambda x: x.freq))
    freq_ht = freq_ht.transmute_globals(
        freq_meta=freq_ht.freq_meta.flatmap(lambda x: x.freq_meta)
    )

    # Create frequency index dictionary on concatenated array (i.e., including
    # all subsets)
    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict(
            freq_meta=hl.eval(freq_ht.freq_meta),
        )
    )
    return freq_ht


def _pre_process_subset_freq(subset_ht: hl.Table, global_ht: hl.Table) -> hl.Table:
    """Pre-process frequency HT.

    :param subset_ht: Subset frequency HT.
    :param global_ht: Global frequency HT.
    :return: Pre-processed frequency HT.
    """
    # Fill in missing freq structs
    ht = subset_ht.join(global_ht.select().select_globals(), how="right")
    ht = ht.annotate(
        freq=hl.if_else(
            hl.is_missing(ht.freq),
            hl.map(lambda x: missing_callstats_expr(), hl.range(hl.len(ht.freq_meta))),
            ht.freq,
        )
    )
    ht = ht.select("freq").select_globals("freq_meta")
    return ht


def main(args):
    """Script to update update Allele Frequencies for HGDP + 1KG subset for v4."""
    test = args.test

    hl.init(
        log="/update_AFs_HGDP+TGP.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    logger.info("Loading the HGDP + 1KG subset dense MT...")
    mt = hl.read_matrix_table(hgdp_tgp_subset(dense=True))

    logger.info("Loading the release HT...")
    ht = hl.read_table(release_ht_path(release_version="3.1.2"))

    logger.info("Updating the HGDP pop labels...")
    ht = update_hgdp_pop_labels(ht, old_pop="miaozu", new_pop="miao")
    ht = update_hgdp_pop_labels(ht, old_pop="yizu", new_pop="yi")
    ht = update_hgdp_pop_labels(ht, old_pop="bantusafrica", new_pop="bantusouthafrica")

    logger.info("Selecting `freq` and `freq_meta` from the release HT...")
    ht = ht.select("freq").select_globals("freq_meta")

    if test:
        logger.info("Filtering to 10kb in DRD2 in MT for testing purposes...")
        test_interval = [
            hl.parse_locus_interval(
                "chr11:113425000-113435000", reference_genome="GRCh38"
            )
        ]
        mt = hl.filter_intervals(mt, test_interval)
        ht = hl.filter_intervals(ht, test_interval)

    logger.info("Loading HGDP_TGP subset meta HT...")
    meta_ht = hgdp_tgp_subset_annotations(sample=True).ht()

    logger.info("Adding updated sample QC annotations to meta HT...")
    meta_ht = add_updated_sample_qc_annotations(meta_ht)
    meta_ht = meta_ht.checkpoint(hgdp_tgp_meta_updated.path, _read_if_exists=True)

    logger.info("Filtering samples in meta HT that will be added and subtracted...")
    samples_to_add, samples_to_subtract = _get_filtered_samples(meta_ht)

    logger.info("Annotating MT with the new meta HT...")
    mt = mt.select_globals()
    mt = mt.annotate_cols(**meta_ht[mt.col_key])

    logger.info("Calculating AFs for added samples...")
    af_added_samples_all = calculate_AFs_for_selected_samples(
        mt, samples_to_add, pop_old=False
    )
    af_added_samples_hgdp = calculate_AFs_for_selected_samples(
        mt.filter_cols(mt.hgdp_tgp_meta.project == "HGDP"),
        samples_to_add,
        subsets=["hgdp"],
        pop_old=False,
    )
    af_added_samples_tgp = calculate_AFs_for_selected_samples(
        mt.filter_cols(mt.hgdp_tgp_meta.project == "1000 Genomes"),
        samples_to_add,
        subsets=["tgp"],
        pop_old=False,
    )

    logger.info("Concatenating AFs for added samples...")
    subset_freq_hts = [
        af_added_samples_hgdp.select("freq").select_globals("freq_meta"),
        af_added_samples_tgp.select("freq").select_globals("freq_meta"),
    ]
    freq_ht_added = _concatenate_subset_frequencies(
        af_added_samples_all, subset_freq_hts
    )

    logger.info("Calculating AFs for subtracted samples...")
    af_subtracted_samples_all = calculate_AFs_for_selected_samples(
        mt, samples_to_subtract, pop_old=True
    )
    af_subtracted_samples_hgdp = calculate_AFs_for_selected_samples(
        mt.filter_cols(mt.hgdp_tgp_meta.project == "HGDP"),
        samples_to_subtract,
        subsets=["hgdp"],
        pop_old=True,
    )
    af_subtracted_samples_tgp = calculate_AFs_for_selected_samples(
        mt.filter_cols(mt.hgdp_tgp_meta.project == "1000 Genomes"),
        samples_to_subtract,
        subsets=["tgp"],
        pop_old=True,
    )

    logger.info("Concatenating AFs for subtracted samples...")
    subset_freq_hts = [
        af_subtracted_samples_hgdp.select("freq").select_globals("freq_meta"),
        af_subtracted_samples_tgp.select("freq").select_globals("freq_meta"),
    ]
    freq_ht_subtracted = _concatenate_subset_frequencies(
        af_subtracted_samples_all, subset_freq_hts
    )

    logger.info("Annotating HT with AFs for added and subtracted samples...")
    ht = ht.annotate(
        freq_added_samples=freq_ht_added[ht.key].freq,
        freq_subtracted_samples=freq_ht_subtracted[ht.key].freq,
    )
    ht = ht.annotate_globals(
        freq_meta_added_samples=freq_ht_added.index_globals().freq_meta,
        freq_meta_subtracted_samples=freq_ht_subtracted.index_globals().freq_meta,
    )

    logger.info("Merging AFs for subtracted samples first...")
    [freq, freq_meta] = merge_freq_arrays(
        [ht.freq, ht.freq_subtracted_samples],
        [ht.freq_meta, ht.freq_meta_subtracted_samples],
        operation="diff",
        set_negatives_to_zero=True,
    )

    ht = ht.annotate(freq1=freq)
    ht = ht.annotate_globals(freq_meta1=freq_meta)

    logger.info("Merging AFs for added samples...")
    [freq, freq_meta] = merge_freq_arrays(
        [ht.freq1, ht.freq_added_samples],
        [ht.freq_meta1, ht.freq_meta_added_samples],
        set_negatives_to_zero=True,
    )

    ht = ht.annotate(freq2=freq)
    ht = ht.annotate_globals(freq_meta2=freq_meta)

    logger.info("Writing out the AF HT...")
    if test:
        ht.write(hgdp_tgp_updated_AF(test=True).path, overwrite=args.overwrite)
    else:
        ht.write(hgdp_tgp_updated_AF().path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script updates AFs for HGDP + 1KG subset for v4 release HT."
    )
    parser.add_argument(
        "--test",
        help="Test on a subset of variants in DRD2 gene",
        action="store_true",
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
