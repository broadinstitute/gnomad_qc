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
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq
from gnomad_qc.v3.resources.release import (
    hgdp_tgp_subset,
    hgdp_tgp_subset_annotations,
    release_sites,
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

POP_MAP = {
    "bantusafrica": "bantusouthafrica",
    "biakaPygmy": "biaka",
    "italian": "bergamoitalian",
    "mbutiPygmy": "mbuti",
    "melanesian": "bougainville",
    "mongola": "mongolian",
    "miaozu": "miao",
    "yizu": "yi",
}


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
            "gnomad_labeled_subpop": populations.population.lower(),
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


def get_filtered_samples(ht: hl.Table) -> Tuple[hl.Table, hl.Table, hl.Table]:
    """
    Compare the old and new sample filters.

    .. note::

        We will use the new sample filters for the v4 release sites HT. To prepare for
        the two sets of samples, we will use anti-join, while taking into account the
        samples in the two populations `Han` & `Papuan`.

    :param ht: Table with the HGDP + 1KG subset metadata with the old and new sample
        filters.
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

    # Filter to all samples with a pop name split and no change in release status.
    split_pops = hl.literal(["PapuanHighlands", "PapuanSepik", "Han", "NorthernHan"])
    pop_diff_ht = ht.filter(
        split_pops.contains(ht.hgdp_tgp_meta.population)
        & ~ht.sample_annotations_updated.contains("gnomad_release")
    )
    logger.info("%d samples in the split Papuan & Han set", pop_diff_ht.count())

    return release_add_ht, release_subtract_ht, pop_diff_ht


def filter_to_test(t):
    """
    Filter to 10kb in DRD2 in MT for testing purposes.

    :param t: MT or HT to filter
    """
    logger.info("Filtering to 10kb in DRD2 in MT for testing purposes...")
    test_interval = [
        hl.parse_locus_interval("chr11:113425000-113435000", reference_genome="GRCh38")
    ]

    return hl.filter_intervals(t, test_interval)


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


def calculate_callstats_for_selected_samples(
    mt: hl.MatrixTable,
    samples_ht: hl.Table,
    subsets: Optional[List[str]] = None,
) -> hl.Table:
    """
    Calculate the call stats for samples in `samples_ht`.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with selected samples and their metadata
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

    if subsets:
        freq_meta = [
            {**x, **{"subset": "|".join(subsets)}} for x in hl.eval(ht.freq_meta)
        ]
        ht = ht.annotate_globals(freq_meta=freq_meta)

    ht = ht.checkpoint(
        hl.utils.new_temp_file("hgdp_tgp_subset_freq", extension="ht"), overwrite=True
    )

    return ht


def concatenate_subset_frequencies(subset_freq_hts: List[hl.Table]) -> hl.Table:
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
        subset_freq_hts,
        data_field_name="freq",
        global_field_name="freq_meta",
    )
    freq_ht = freq_ht.transmute(freq=freq_ht.freq.flatmap(lambda x: x.freq))
    freq_ht = freq_ht.transmute_globals(
        freq_meta=freq_ht.freq_meta.flatmap(lambda x: x.freq_meta),
        freq_meta_sample_count=freq_ht.freq_meta.flatmap(
            lambda x: x.freq_meta_sample_count
        ),
    )

    return freq_ht


def calculate_concatenate_callstats(
    mt: hl.MatrixTable, samples_ht: hl.Table, compute_freq_all: bool = True
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
    projects = []
    subset_freq_hts = []
    logger.info("Calculating AFs for selected samples...")
    if compute_freq_all:
        projects.append(("All", None))

    projects.extend([("HGDP", ["hgdp"]), ("1000 Genomes", ["tgp"])])
    project_expr = samples_ht.hgdp_tgp_meta.project
    project_n = samples_ht.aggregate(hl.agg.counter(project_expr))

    for project, subset in projects:
        if project == "All":
            subset_ht = samples_ht
        else:
            subset_ht = samples_ht.filter(project_expr == project)
        if project == "All" or project_n.get(project, 0) > 0:
            freq_ht = calculate_callstats_for_selected_samples(mt, subset_ht, subset)
            freq_ht = freq_ht.select("freq").select_globals(
                "freq_meta", "freq_meta_sample_count"
            )
            subset_freq_hts.append(freq_ht)

    logger.info("Concatenating AFs for selected samples...")
    freq_ht = concatenate_subset_frequencies(subset_freq_hts)

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


def get_v4_genomes_release_resources(
    test: bool, overwrite: bool
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed to create the gnomAD v4 genomes release.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        gnomAD v4 genomes release pipeline.
    """
    # Initialize gnomAD v4 genomes release pipeline resource collection.
    v4_genome_release_pipeline = PipelineResourceCollection(
        pipeline_name="gnomad_v4_genomes_release",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the v4 genomes release pipeline.
    update_annotations = PipelineStepResourceCollection(
        "--update-annotations",
        input_resources={
            "Released HGDP + 1KG sample metadata": {
                "meta_ht": hgdp_tgp_subset_annotations(sample=True).versions["3.1.2"]
            }
        },
        output_resources={"updated_meta_ht": hgdp_tgp_meta_updated},
    )

    get_callstats_for_updated_samples = PipelineStepResourceCollection(
        "--get-callstats-for-updated-samples",
        pipeline_input_steps=[update_annotations],
        add_input_resources={
            "gnomAD v3.1.2 HGDP + 1KG subset dense MT": {
                "dense_mt": hgdp_tgp_subset(dense=True, public=True).versions["3.1.2"],
            }
        },
        output_resources={
            "freq_added_ht": hgdp_tgp_updated_callstats(test=test, subset="added"),
            "freq_subtracted_ht": hgdp_tgp_updated_callstats(
                test=test, subset="subtracted"
            ),
            "freq_pop_diff_ht": hgdp_tgp_updated_callstats(
                test=test, subset="pop_diff"
            ),
        },
    )
    update_release_callstats = PipelineStepResourceCollection(
        "--update-release-callstats",
        pipeline_input_steps=[get_callstats_for_updated_samples],
        add_input_resources={
            "gnomAD v3.1.2 release sites HT": {
                "sites_ht": release_sites(public=True).versions["3.1.2"]
            }
        },
        output_resources={
            "v4_freq_ht": hgdp_tgp_updated_callstats(test=test, subset="final"),
        },
    )

    # Add all steps to the gnomAD v4 genomes release pipeline resource collection.
    v4_genome_release_pipeline.add_steps(
        {
            "update_annotations": update_annotations,
            "get_callstats_for_updated_samples": get_callstats_for_updated_samples,
            "update_release_callstats": update_release_callstats,
        }
    )

    return v4_genome_release_pipeline


def get_sample_count_from_v3_subsets(ht: hl.Table) -> hl.Table:
    """
    Get the freq sample counts from the v3 subsets.

    :param ht: Table of the 3.1.2 release sites HT, which only contains the freq
       sample counts for non subset groupings.
    :return: Table with the freq sample counts from the v3 subsets added.
    """
    SUBSETS = [
        "non_v2",
        "non_topmed",
        "non_cancer",
        "controls_and_biobanks",
        "non_neuro",
        "tgp",
        "hgdp",
    ]

    sc = ht.freq_sample_count.collect(_localize=False)[0]

    for subset in SUBSETS:
        freq_subset_ht = hl.read_table(
            get_freq(version="3.1", subset=subset, het_nonref_patch=True).path
        )
        sc_subset = freq_subset_ht.freq_sample_count.collect(_localize=False)[0]
        sc = sc.extend(sc_subset)

    ht = ht.annotate_globals(freq_sample_count=sc)

    return ht


def main(args):
    """
    Script to update call stats for HGDP + 1KG subset for v4.

    .. note::
    This code is specifically designed for the update HGDP + 1KG subset in
    v4 release HT. There are a few major changes compared to the v3 release:
        - the following new sample filters are applied: hard filters,
        pop PC outliers,relatedness within the subset and relatedness to
        the rest of the release
        - the new pop labels
        - the new split of the `Han` and `Papuan` samples
    In order to avoid re-calculating the callstats for the whole subset / whole
    release, we will calculate the callstats for the samples that will be added
    and subtracted, then merge the callstats with the old callstats in the
    release HT.
    """
    test = args.test
    overwrite = args.overwrite

    hl.init(
        log="/update_hgdp_tgp_af.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    v4_genome_release_resources = get_v4_genomes_release_resources(
        test=test, overwrite=overwrite
    )

    if args.update_annotations:
        res = v4_genome_release_resources.update_annotations
        res.check_resource_existence()
        logger.info("Adding updated sample QC annotations to meta HT...")
        meta_ht = res.meta_ht.ht()
        meta_ht = add_updated_sample_qc_annotations(meta_ht)
        # TODO: temporarily using _read_if_exists, until we have new fields to
        #  be updated.
        meta_ht.write(hgdp_tgp_meta_updated.path, overwrite=True)

    if args.get_callstats_for_updated_samples:
        res = v4_genome_release_resources.get_callstats_for_updated_samples
        res.check_resource_existence()

        logger.info("Loading HGDP + 1KG subset dense MT...")
        mt = res.dense_mt.mt().select_globals()

        if test:
            mt = filter_to_test(mt)

        logger.info("Filtering samples in meta HT that will be added and subtracted...")
        meta_ht = res.updated_meta_ht.ht()
        add_ht, subtract_ht, pop_diff_ht = get_filtered_samples(meta_ht)

        logger.info(
            "Calculating and concatenating callstats for samples to be added and "
            "samples to be subtracted..."
        )
        freq_added_ht = calculate_concatenate_callstats(mt, add_ht)
        freq_added_ht.write(res.freq_added_ht.path, overwrite=overwrite)
        freq_subtracted_ht = calculate_concatenate_callstats(mt, subtract_ht)
        freq_subtracted_ht.write(res.freq_subtracted_ht.path, overwrite=overwrite)
        freq_pop_diff_ht = calculate_concatenate_callstats(
            mt, pop_diff_ht, compute_freq_all=False
        )
        freq_meta_expr, freq_expr = filter_arrays_by_meta(
            freq_pop_diff_ht.freq_meta,
            {
                "freq": freq_pop_diff_ht.freq,
                "freq_meta_sample_count": (
                    freq_pop_diff_ht.index_globals().freq_meta_sample_count
                ),
            },
            ["pop"],
        )
        freq_pop_diff_ht = freq_pop_diff_ht.annotate(freq=freq_expr["freq"])
        freq_pop_diff_ht = freq_pop_diff_ht.annotate_globals(
            freq_meta=freq_meta_expr,
            freq_meta_sample_count=freq_expr["freq_meta_sample_count"],
        )
        freq_pop_diff_ht.write(res.freq_pop_diff_ht.path, overwrite=overwrite)

    if args.update_release_callstats:
        res = v4_genome_release_resources.update_release_callstats
        res.check_resource_existence()

        logger.info("Loading v3.1.2 full sites HT...")
        ht = res.sites_ht.ht()

        if test:
            ht = filter_to_test(ht)

        logger.info(
            "Removing `Han` and `Papuan` populations from freq and freq_meta..."
        )
        pops_to_remove = {"pop": ["han", "papuan"]}
        freq_meta_expr, freq_expr = filter_arrays_by_meta(
            ht.freq_meta,
            {"freq": ht.freq},
            pops_to_remove,
            keep=False,
            combine_operator="or",
        )
        ht = ht.annotate(freq=freq_expr["freq"])
        ht = ht.annotate_globals(freq_meta=freq_meta_expr)

        logger.info("Updating HGDP pop labels...")
        ht = update_hgdp_pop_labels(ht, POP_MAP)

        # TODO: We won't need this anymore because we need to update more annotations
        #  than just this.
        logger.info("Selecting `freq` and `freq_meta` from the release HT...")
        ht = ht.select("freq").select_globals("freq_meta")

        logger.info("Annotating HT with AFs for added and subtracted samples...")
        freq_hts = [
            ht,
            res.freq_pop_diff_ht.ht().select("freq").select_globals("freq_meta"),
            res.freq_added_ht.ht().select("freq").select_globals("freq_meta"),
            res.freq_subtracted_ht.ht().select("freq").select_globals("freq_meta"),
        ]
        ht = hl.Table.multi_way_zip_join(freq_hts, "ann_array", "global_array")
        ht = ht.checkpoint(
            "gs://gnomad-tmp/julia/hgdp_tgp_update/update_hgdp_tgp_af.ht",
            overwrite=True,
        )

        logger.info("Merging AFs from diff pop samples and added samples...")
        freq_expr, freq_meta_expr = merge_freq_arrays(
            [ht.ann_array[i].freq for i in range(3)],
            [ht.global_array[i].freq_meta for i in range(3)],
            operation="sum",
        )

        logger.info("Merging AFs from subtracted samples...")
        freq_expr, freq_meta_expr = merge_freq_arrays(
            [freq_expr, ht.ann_array[3].freq],
            [freq_meta_expr, ht.global_array[3].freq_meta],
            operation="diff",
            set_negatives_to_zero=False,
        )
        ht = ht.select(freq=freq_expr)
        logger.info("Making freq_index_dict from freq_meta...")
        ht = ht.select_globals(
            freq_meta=freq_meta_expr,
            freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta_expr)),
        )

        logger.info(
            "Set Y-variant frequency callstats for female-specific "
            "metrics to missing structs..."
        )
        ht = ht.annotate(freq=set_female_y_metrics_to_na_expr(ht))

        logger.info("Writing out the AF HT...")
        ht.write(res.v4_freq_ht.path, overwrite=overwrite)


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

    parser.add_argument(
        "--update-annotations", help="Update sample QC annotations", action="store_true"
    )
    parser.add_argument(
        "--get-callstats-for-updated-samples", help="", action="store_true"
    )
    parser.add_argument("--update-release-callstats", help="", action="store_true")

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
