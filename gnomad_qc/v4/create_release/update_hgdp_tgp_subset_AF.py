"""Script to update the call stats of the HGDP + 1KG subset for v4 release HT."""
import argparse
import logging
from typing import List, Optional, Tuple

import hail as hl
from gnomad.utils.annotations import (
    annotate_freq,
    merge_freq_arrays,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.release import make_freq_index_dict
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

    def _update_annotations(expr, sample_annotations):
        if isinstance(sample_annotations, dict):
            updated = {}
            updated_flag = {}
            for ann, updated_expr in sample_annotations.items():
                bla1, bla2 = _update_annotations(expr[ann], updated_expr)
                updated_flag.update(
                    {ann + ("." + k if k else ""): v for k, v in bla1.items()}
                )
                updated[ann] = bla2
            if isinstance(expr, hl.Table):
                updated_flag = add_filters_expr(filters=updated_flag)
                return expr.annotate(**updated, sample_annotations_updated=updated_flag)
            return updated_flag, expr.annotate(**updated)
        else:
            return {"": sample_annotations != expr}, sample_annotations

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
    updated_ht = _update_annotations(ht, sample_annotations_to_update)
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


def _get_filtered_samples(ht: hl.Table) -> Tuple[hl.Table, hl.Table]:
    """
    Compare the old and new sample filters.

    .. note::

        We will use the new sample filters for the v4 release sites HT. To prepare for
        the two sets of samples, we will use anti-join, while taking into account the
        samples in the two populations `Han` & `Papuan`.

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
    :param subsets: Subsets to be used for the AFs calculation: `hgdp` or `tgp`, or both.
    :return: Table with the AFs for the selected samples.
    """
    mt = mt.filter_cols(hl.is_defined(samples_ht[mt.col_key]))

    if pop_old:
        mt = annotate_freq(
            mt,
            sex_expr=mt.gnomad_sex_imputation.sex_karyotype,
            pop_expr=mt.hgdp_tgp_meta.population.lower() if subsets else None,
        )
    else:
        mt = annotate_freq(
            mt,
            sex_expr=mt.gnomad_sex_imputation.sex_karyotype,
            pop_expr=mt.hgdp_tgp_meta.population_updated.lower() if subsets else None,
        )
    # make population labels lowercase to match the constants in release HT.
    # TODO: Julia, could you double check the issue that I messaged you about?
    # I got extra groupings if I use pop_expr for the combined HGDP + 1KG
    # subset.

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
    # TODO: warnings at this step
    # Hail: WARN: Name collision: field 'cols' already in object dict.
    #   This field must be referenced with __getitem__ syntax: obj['cols']
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

    # TODO: temporarily not overwriting freq or freq_meta, so we can examine the output
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

    # TODO: writing step is quite slow, should be optimized
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
