"""Script to update the call stats of the HGDP + 1KG subset for v4 release HT."""
import argparse
import logging
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.utils.annotations import (
    annotate_freq,
    merge_freq_arrays,
    set_female_y_metrics_to_na_expr,
    update_structured_annotations,
)
from gnomad.utils.filtering import filter_freq_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.release import (
    hgdp_tgp_subset,
    hgdp_tgp_subset_annotations,
    release_ht_path,
)
from gnomad_qc.v4.resources.release import hgdp_tgp_updated_callstats
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
    """
    Add updated sample QC annotations to the HGDP + 1KG subset.

    .. note::

        The following annotations are updated based on the latest sample QC:
            - `sample_filters.hard_filtered`: to apply the recomputed freemix
              filter for HGDP samples.
            - `sample_filters.relatedness_inference.related`: to apply the
              updated relatedness inference implemented by Alicia Martin's group.
            - `sample_filters.pop_outlier`: to apply the updated pop outlier
              filter implemented by Alicia Martin's group.

    :param ht: Table with the HGDP + 1KG subset metadata from the last release.
    :return: Table with updated sample QC annotations.
    """
    # Load all updated sample QC Tables.
    contamination_ht = hgdp_recomputed_freemix.ht()
    relatedness_ht = hgdp_tgp_related_samples_to_drop.ht()
    relatedness_ht = relatedness_ht.key_by(s=relatedness_ht.node.s)
    pop_outliers_ht = hgdp_tgp_pop_outliers.ht()
    populations_ht = hgdp_tgp_populations_updated.ht()

    # TODO: rerun the pc_relate to get the updated relatedness HT,
    #  because Alicia's group didn't checkpoint the results.
    # TODO: get the global & subcontinental PCs for the updated HGDP + 1KG subset:
    #  hgdp_tgp_meta.global_pca_scores,
    #  hgdp_tgp_meta.subcontinental_pca.pca_scores,
    #  hgdp_tgp_meta.subcontinental_pca.pca_scores_outliers_removed,
    #  could update the meta HT after MVP release

    # Update annotations impacted by the recomputed freemix.
    freemix = hl.coalesce(
        contamination_ht[ht.s].recomputed_contam_rate * 100, ht.bam_metrics.freemix
    )
    hard_filters = ht.gnomad_sample_filters.hard_filters | hl.if_else(
        freemix > 5.0, hl.set({"contamination"}), hl.empty_set(hl.tstr)
    )
    hard_filtered = hl.len(hard_filters) > 0
    # Update annotations impacted by the updated relatedness inference.
    related = hl.is_defined(relatedness_ht[ht.key])
    # Update annotations impacted by the updated HGDP + 1KG metadata.
    outlier = hl.is_defined(pop_outliers_ht[ht.key])
    populations = populations_ht[ht.key]
    # TODO: Will need to update based on gnomad_sample_filters.release_related decision.
    gnomad_release = (
        ~hard_filtered
        & ~outlier
        & ~related  # & ~ht.gnomad_sample_filters.release_related
    )

    sample_annotations_to_update = {
        "bam_metrics": {"freemix": freemix},
        "gnomad_sample_filters": {
            "hard_filters": hard_filters,
            "hard_filtered": hard_filtered,
        },
        "gnomad_high_quality": ht.gnomad_high_quality & ~hard_filtered,
        "gnomad_release": gnomad_release,
        "relatedness_inference": {"related": related},
        "hgdp_tgp_meta": {
            "subcontinental_pca": {"outlier": outlier},
            "population": populations.population,
            "latitude": populations.latitude,
            "longitude": populations.longitude,
        },
        "high_quality": ~hard_filtered & ~outlier,
    }
    updated_ht = update_structured_annotations(
        ht,
        sample_annotations_to_update,
        annotation_update_label="sample_annotations_updated",
    )
    updated_ht = updated_ht.checkpoint(
        hl.utils.new_temp_file("hgdp_tgp_meta_update", extension="ht"), overwrite=True
    )
    updated_counts = updated_ht.aggregate(
        hl.struct(
            n_hard_filtered=hl.agg.count_where(
                updated_ht.gnomad_sample_filters.hard_filtered
            ),
            n_related=hl.agg.count_where(updated_ht.relatedness_inference.related),
            n_outlier=hl.agg.count_where(
                updated_ht.hgdp_tgp_meta.subcontinental_pca.outlier
            ),
            n_diff=hl.agg.explode(
                lambda x: hl.agg.counter(x), updated_ht.sample_annotations_updated
            ),
        )
    )
    logger.info(
        """%d sample(s) hard filtered
        %d more sample(s) hard filtered than last release
        %d sample(s) related
        %d sample(s) were assigned differently for the relatedness than last release
        %d sample(s) are pop outliers
        %d sample(s) were assigned differently as subcontinental outlier than last release
        %d samples have different population labels than last release""",
        updated_counts["n_hard_filtered"],
        updated_counts["n_diff"]["gnomad_sample_filters.hard_filtered"],
        updated_counts["n_related"],
        updated_counts["n_diff"]["relatedness_inference.related"],
        updated_counts["n_outlier"],
        updated_counts["n_diff"]["hgdp_tgp_meta.subcontinental_pca.outlier"],
        updated_counts["n_diff"]["hgdp_tgp_meta.population"],
    )
    return updated_ht


def get_filtered_samples(ht: hl.Table) -> Tuple[hl.Table, hl.Table]:
    """
    Compare the old and new sample filters.

    .. note::

        We will use the new sample filters for the v4 release sites HT. To prepare for
        the two sets of samples, we will use anti-join, while taking into account the
        samples in the two populations `Han` & `Papuan`.

    :param ht: Table with the HGDP + 1KG subset metadata with the old and new sample filters.
    :return: Table with the old and new sample filters.
    """
    # Filter to all samples with a change in the gnomad_release status.
    release_diff_ht = ht.filter(
        ht.sample_annotations_updated.contains("gnomad_release")
    )
    release_add_ht = release_diff_ht.filter(release_diff_ht.gnomad_release)
    release_subtract_ht = release_diff_ht.filter(~release_diff_ht.gnomad_release)
    logger.info(
        """
        %d samples will be added after the new filters
        %d sample(s) will be removed after the new filters""",
        release_add_ht.count(),
        release_subtract_ht.count(),
    )
    # Filter to all samples with a population name split.
    split_pops = hl.literal(["PapuanHighlands", "PapuanSepik", "Han", "NorthernHan"])
    pop_diff_ht = ht.filter(split_pops.contains(ht.hgdp_tgp_meta.population))
    release_add_ht = release_add_ht.union(pop_diff_ht)
    logger.info(
        "%d samples in the `sum` set, including the splitted Papuan & Han",
        release_add_ht.count(),
    )
    return release_add_ht, release_subtract_ht


def calculate_callstats_for_selected_samples(
    mt: hl.MatrixTable,
    samples_ht: hl.Table,
    subsets: Optional[List[str]] = None,
) -> hl.Table:
    """
    Calculate the call stats for samples in `samples_ht`.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with the samples to be added or subtracted.
    :param subsets: Subsets to be used for the AFs calculation: `hgdp` or `tgp`, or both.
    :return: Table with the call stats for the selected samples.
    """
    logger.info("Filtering MT to sample subset and annotating with metadata...")
    mt = mt.filter_cols(hl.is_defined(samples_ht[mt.col_key]))
    mt = mt.annotate_cols(**samples_ht[mt.col_key])

    ht = annotate_freq(
        mt,
        sex_expr=mt.gnomad_sex_imputation.sex_karyotype,
        pop_expr=(
            mt.hgdp_tgp_meta.population.lower()
            if subsets
            else mt.gnomad_population_inference.pop
        ),
        annotate_mt=False,
    )
    # make population labels lowercase to match the constants in release HT.
    # TODO: Julia, could you double check the issue that I messaged you about?
    # I got extra groupings if I use pop_expr for the combined HGDP + 1KG
    # subset.

    if subsets:
        freq_meta = [
            {**x, **{"subset": "|".join(subsets)}} for x in hl.eval(ht.freq_meta)
        ]
        ht = ht.annotate_globals(freq_meta=freq_meta)
        # ht = ht.annotate_globals(
        #     freq_index_dict=make_freq_index_dict(
        #         freq_meta=freq_meta,
        #         subsets=["|".join(subsets)],
        #     )
        # )
    # else:
    #     ht = ht.annotate_globals(
    #         freq_index_dict=make_freq_index_dict_from_meta(freq_meta=hl.eval(ht.freq_meta))
    #     )

    # ht = ht.annotate(freq=set_female_y_metrics_to_na_expr(ht))
    ht = ht.checkpoint(
        hl.utils.new_temp_file("hgdp_tgp_subset_freq", extension="ht"), overwrite=True
    )
    return ht


def update_hgdp_pop_labels(ht: hl.Table, pop_map: Dict[str, Any]) -> hl.Table:
    """
    Update the population labels for HGDP samples in the release HT.

    :param ht: release HT with old population labels.
    :param pop_map: dictionary with old and new population labels.
    :return: release HT with updated population labels.
    """
    pop_map = hl.literal(pop_map)
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: hl.or_else(pop_map.get(x), x))
        ),
        freq_index_dict=hl.dict(
            hl.zip(
                ht.freq_index_dict.keys().map(lambda k: hl.or_else(pop_map.get(k), k)),
                ht.freq_index_dict.values(),
            )
        ),
    )
    return ht


def concatenate_subset_frequencies(
    global_freq_ht: hl.Table, subset_freq_hts: List[hl.Table]
) -> hl.Table:
    """
    Concatenate subset frequencies into a single Table.

    .. note::

        This is the same as we did in
        gnomad_qc/v3/create_release/create_release_sites_ht.py#L291-L310, except
        excluding the pops & downsampling options.

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

    return freq_ht


def calculate_concatenate_callstats(
    mt: hl.MatrixTable, samples_ht: hl.Table
) -> hl.Table:
    """
    Calculate the call stats for samples in `samples_ht`.

    .. note::
        This function is used to calculate callstats for all samples provided in
        the samples_ht, also samples that belong to HGDP and TGP project separately,
        then concatenate the callstats into a single Table.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with the samples.
    :return: Table with the call stats for the selected samples.
    """
    logger.info("Calculating AFs for selected samples...")
    freq_ht_all = calculate_callstats_for_selected_samples(
        mt,
        samples_ht,
    )
    freq_ht_hgdp = calculate_callstats_for_selected_samples(
        mt,
        samples_ht.filter(samples_ht.hgdp_tgp_meta.project == "HGDP"),
        subsets=["hgdp"],
    )
    freq_ht_tgp = calculate_callstats_for_selected_samples(
        mt,
        samples_ht.filter(samples_ht.hgdp_tgp_meta.project == "1000 Genomes"),
        subsets=["tgp"],
    )

    logger.info("Concatenating AFs for selected samples...")
    subset_freq_hts = [
        freq_ht_hgdp.select("freq").select_globals("freq_meta"),
        freq_ht_tgp.select("freq").select_globals("freq_meta"),
    ]
    freq_ht = concatenate_subset_frequencies(freq_ht_all, subset_freq_hts)

    return freq_ht


def remove_pops_from_freq_meta(ht: hl.Table, pops_to_remove: List[str]) -> hl.Table:
    """
    Remove populations from the freq_meta.

    :param ht: Table with call stats.
    :param pops_to_remove: List of populations to be removed.
    :return: Table with the populations removed.
    """
    pops_to_remove = hl.literal(pops_to_remove)
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.filter(lambda x: ~pops_to_remove.contains(x))
        ),
        freq_index_dict=hl.dict(
            hl.zip(
                ht.freq_index_dict.keys().filter(lambda k: ~pops_to_remove.contains(k)),
                ht.freq_index_dict.values(),
            )
        ),
    )
    return ht


def main(args):
    """Script to update update call stats for HGDP + 1KG subset for v4."""
    test = args.test

    hl.init(
        log="/update_hgdp_tgp_af.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    logger.info("Loading the HGDP + 1KG subset dense MT...")
    mt = hl.read_matrix_table(hgdp_tgp_subset(dense=True)).select_globals()

    logger.info("Loading the release HT...")
    ht = hl.read_table(release_ht_path(release_version="3.1.2"))

    logger.info("Updating the HGDP pop labels...")
    pop_map = {
        "bantusafrica": "bantusouthafrica",
        "biakaPygmy": "biaka",
        "italian": "bergamoitalian",
        "mbutiPygmy": "mbuti",
        "melanesian": "bougainville",
        "mongola": "mongolian",
        "miaozu": "miao",
        "yizu": "yi",
    }
    ht = update_hgdp_pop_labels(ht, pop_map)

    logger.info("Selecting `freq` and `freq_meta` from the release HT...")
    ht = ht.select("freq").select_globals("freq_meta")

    logger.info("Removing `Han` and `Papuan` populations from freq and freq_meta...")
    pops_to_remove = {"pop": ["pop", "papuan"]}
    freq, freq_meta = filter_freq_by_meta(
        ht.freq, ht.freq_meta, pops_to_remove, keep=False, operator="or"
    )
    ht = ht.annotate(freq=freq)
    ht = ht.annotate_globals(freq_meta=freq_meta)

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

    if args.update_annotations:
        logger.info("Adding updated sample QC annotations to meta HT...")
        meta_ht = add_updated_sample_qc_annotations(meta_ht)
        # TODO: temporarily using _read_if_exists, until we have new fields to be
        # updated.
        meta_ht = meta_ht.checkpoint(hgdp_tgp_meta_updated.path, _read_if_exists=True)

    logger.info("Filtering samples in meta HT that will be added and subtracted...")
    samples_to_add, samples_to_subtract = get_filtered_samples(meta_ht)

    logger.info(
        "Calculating and concatenating callstats for samples to be added and samples to"
        " be subtracted..."
    )
    freq_ht_added = calculate_concatenate_callstats(mt, samples_to_add)
    freq_ht_added = freq_ht_added.checkpoint(
        hgdp_tgp_updated_callstats(test=test, subset="added").path,
        overwrite=args.overwrite,
    )

    freq_ht_subtracted = calculate_concatenate_callstats(mt, samples_to_subtract)
    freq_ht_subtracted = freq_ht_subtracted.checkpoint(
        hgdp_tgp_updated_callstats(test=test, subset="subtracted").path,
        overwrite=args.overwrite,
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
    freq, freq_meta = merge_freq_arrays(
        [ht.freq, ht.freq_subtracted_samples],
        [ht.freq_meta, ht.freq_meta_subtracted_samples],
        operation="diff",
        set_negatives_to_zero=True,
    )

    # TODO: temporarily not overwriting freq or freq_meta, so we can examine the output
    ht = ht.annotate(freq1=freq)
    ht = ht.annotate_globals(freq_meta1=freq_meta)

    logger.info("Merging AFs for added samples...")
    freq, freq_meta = merge_freq_arrays(
        [ht.freq1, ht.freq_added_samples],
        [ht.freq_meta1, ht.freq_meta_added_samples],
        operation="sum",
    )

    ht = ht.annotate(freq2=freq)
    ht = ht.annotate_globals(freq_meta2=freq_meta)
    ht = ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(ht.freq_meta2)
    )
    ht = ht.annotate(freq2=set_female_y_metrics_to_na_expr(ht))

    logger.info("Writing out the AF HT...")
    if test:
        ht.write(
            hgdp_tgp_updated_callstats(test=True, subset="final").path,
            overwrite=args.overwrite,
        )
    else:
        ht.write(
            hgdp_tgp_updated_callstats(subset="final").path, overwrite=args.overwrite
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script updates AFs for HGDP + 1KG subset for v4 release HT."
    )
    parser.add_argument(
        "--test",
        help="Test on a subset of variants in DRD2 gene",
        action="store_true",
    )
    parser.add_argument(
        "--update-annotations", help="Update sample QC annotations", action="store_true"
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
