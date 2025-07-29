"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies for samples that will be removed from gnomAD v4
in v5 due to relatedness filtering or ancestry changes. It uses a differential approach
to avoid densifying the full dataset.
"""

import argparse
import logging
from typing import Optional, Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq,
    build_freq_stratification_list,
    generate_freq_group_membership_array,
    merge_freq_arrays,
    merge_histograms,
)
from gnomad.utils.release import make_freq_index_dict_from_meta

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
    check_resource_existence,
)
from gnomad_qc.v4.resources.annotations import get_freq as get_v4_freq
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v5.resources.annotations import get_age_hist_ht, get_freq_ht
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import V5_DOWNSAMPLINGS, WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import project_meta, sample_id_collisions
from gnomad_qc.v5.resources.sample_qc import get_gen_anc_ht, related_samples_to_drop

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def get_freq_resources(
    overwrite: bool = False,
    test: Optional[bool] = False,
) -> PipelineResourceCollection:
    """
    Get frequency calculation resources.

    :param overwrite: Whether to overwrite existing results.
    :param test: Whether to run in test mode.
    :return: PipelineResourceCollection for frequency calculations.
    """
    freq_pipeline = PipelineResourceCollection(
        "v5_frequency",
        overwrite=overwrite,
        test=test,
    )

    # Add frequency calculation step
    freq_step = PipelineStepResourceCollection(
        "--run-frequency-calculations",
        input_resources={
            "v4_freq_ht": get_v4_freq(data_type="exomes"),
            "v4_meta_ht": v4_meta(data_type="joint"),
            "project_meta": project_meta,
            "related_samples_to_drop": related_samples_to_drop(test=test),
        },
        output_resources={
            "freq_ht": get_freq_ht(test=test),
            "age_hist_ht": get_age_hist_ht(test=test),
        },
    )

    freq_pipeline.add_steps({"freq_calculation": freq_step})

    return freq_pipeline


def identify_samples_to_remove_for_relatedness(test: bool = False) -> hl.Table:
    """
    Identify samples that will be removed from v5 due to relatedness filtering.

    :param test: Whether to run in test mode.
    :return: Table with samples to remove for relatedness.
    """
    logger.info("Identifying samples to remove for relatedness...")

    # Get samples to drop due to relatedness
    samples_to_drop = related_samples_to_drop(test=test, release=True).ht()

    # Load v5 metadata
    v5_meta_ht = project_meta.ht()

    # Join with metadata to get sample information
    meta_with_drops = v5_meta_ht.annotate(
        will_be_dropped=samples_to_drop[v5_meta_ht.key].s.is_defined()
    )

    return meta_with_drops


def identify_samples_with_ancestry_changes(test: bool = False) -> hl.Table:
    """
    Identify samples that have changed ancestry between v4 and v5.

    :param test: Whether to run in test mode.
    :return: Table with samples that changed ancestry.
    """
    logger.info("Identifying samples with ancestry changes...")

    # Load v4 and v5 metadata
    v4_meta_ht = v4_meta(data_type="joint").ht()
    v5_gen_anc_ht = get_gen_anc_ht(test=test).ht()

    # Compare ancestry assignments
    ancestry_changes = v5_gen_anc_ht.annotate(
        v4_gen_anc=v4_meta_ht[v5_gen_anc_ht.key].gen_anc,
        ancestry_changed=hl.if_else(
            v4_meta_ht[v5_gen_anc_ht.key].pop != v5_gen_anc_ht.gen_anc, True, False
        ),
    ).filter(ancestry_changes.ancestry_changed)

    return ancestry_changes


def calculate_frequency_for_samples(
    samples_ht: hl.Table,
    test: bool = False,
    use_v5_ancestry: bool = False,
) -> hl.Table:
    """
    Calculate frequency statistics for a set of samples using gnomAD v4 VDS.

    This function loads the gnomAD v4 VDS, resolves sample ID collisions, filters to the
    specified samples, and calculates frequency statistics using annotate_freq. The function
    handles both v4 and v5 ancestry data based on the use_v5_ancestry parameter.

    :param samples_ht: Table with samples to process. Must have 's' field as sample ID.
    :param test: Whether to run in test mode. If True, uses test datasets.
    :param use_v5_ancestry: Whether to use v5 ancestry data instead of v4. If True,
        loads v5 genetic ancestry HT and annotates samples with v5 ancestry for
        stratification. If False, uses v4 ancestry from metadata.
    :return: Frequency statistics Table with freq, freq_meta, and freq_meta_sample_count
        annotations. Compatible with v4 frequency HT structure for differential analysis.
    """
    logger.info("Calculating frequency stats for samples...")

    sample_ids = samples_ht.s.collect()
    vds = get_gnomad_v4_vds(release_only=True, annotate_meta=True, test=test)

    logger.info("Resolving sample ID collisions...")
    sample_collisions = sample_id_collisions.ht()

    # Access the sparse MTs individually and rename collisions
    reference_mt = add_project_prefix_to_sample_collisions(
        t=vds.reference_data,
        sample_collisions=sample_collisions,
        project="gnomad",
    )
    variant_mt = add_project_prefix_to_sample_collisions(
        t=vds.variant_data,
        sample_collisions=sample_collisions,
        project="gnomad",
    )

    vds = hl.vds.VariantDataset(reference_data=reference_mt, variant_data=variant_mt)

    vds_filtered = hl.vds.filter_samples(vds, sample_ids, keep=True)

    mt = hl.vds.to_dense_mt(vds_filtered)

    if use_v5_ancestry:
        v5_gen_anc_ht = get_gen_anc_ht(test=test).ht()
        mt = mt.annotate_cols(gen_anc=v5_gen_anc_ht[mt.col_key].gen_anc)
        gen_anc_expr = mt.gen_anc
    else:
        gen_anc_expr = mt.meta.pop

    freq_ht = annotate_freq(
        mt,
        sex_expr=mt.meta.sex_karyotype,
        gen_anc_expr=gen_anc_expr,
        additional_strata_expr=[
            {"data_set": "gnomad"},
            {"data_set": "gnomad", "gen_anc": gen_anc_expr},
        ],
        # Use v4 downsamplings for compatibility with v4 frequency HT when doing differential analysis
        # annotate_freq will automatically filter out downsamplings larger than
        # available sample counts
        downsamplings=hl.literal(DOWNSAMPLINGS["v4"]),
    )

    return freq_ht


def calculate_age_histograms_for_samples(
    samples_ht: hl.Table,
    test: bool = False,
) -> hl.Table:
    """
    Calculate age histograms for a set of samples using gnomAD v4 VDS.

    This function loads the gnomAD v4 VDS, resolves sample ID collisions, filters to the
    specified samples, and calculates age histograms using age_hists_expr. The function
    generates age distributions for heterozygous and homozygous variant carriers.

    :param samples_ht: Table with samples to process. Must have 's' field as sample ID.
    :param test: Whether to run in test mode. If True, uses test datasets.
    :return: Age histogram Table with age_hist_het and age_hist_hom annotations.
        Compatible with v4 age histogram HT structure for differential analysis.
    """
    logger.info("Calculating age histograms for samples...")

    sample_ids = samples_ht.s.collect()

    vds = get_gnomad_v4_vds(release_only=True, annotate_meta=True, test=test)

    logger.info("Resolving sample ID collisions...")
    sample_collisions = sample_id_collisions.ht()

    # Access the sparse MTs individually and rename collisions
    reference_mt = add_project_prefix_to_sample_collisions(
        t=vds.reference_data,
        sample_collisions=sample_collisions,
        project="gnomad",
    )
    variant_mt = add_project_prefix_to_sample_collisions(
        t=vds.variant_data,
        sample_collisions=sample_collisions,
        project="gnomad",
    )

    vds = hl.vds.VariantDataset(reference_data=reference_mt, variant_data=variant_mt)

    vds_filtered = hl.vds.filter_samples(vds, sample_ids, keep=True)

    mt = hl.vds.to_dense_mt(vds_filtered)

    age_hist_expr = age_hists_expr(
        mt.GT,
        mt.meta.age,
        mt.meta.sex_karyotype,
    )

    age_hist_ht = mt.aggregate_rows(
        age_hist_het=hl.agg.explode(age_hist_expr.age_hist_het),
        age_hist_hom=hl.agg.explode(age_hist_expr.age_hist_hom),
    )

    return age_hist_ht


def subtract_frequency_stats(
    base_freq_ht: hl.Table,
    removed_freq_stats: hl.Table,
) -> hl.Table:
    """
    Subtract frequency statistics of removed samples from base frequency HT.

    This function performs differential analysis by subtracting the frequency statistics
    of removed samples from the base frequency HT using merge_freq_arrays with
    operation="diff". This is used to update v4 frequency HT for samples that will
    be removed in v5.

    :param base_freq_ht: Base frequency HT (v4) containing freq, freq_meta, and
        freq_meta_sample_count annotations.
    :param removed_freq_stats: Frequency statistics for removed samples with same
        structure as base_freq_ht.
    :return: Modified frequency HT with removed samples' stats subtracted. Maintains
        same structure as input with updated freq, freq_meta, freq_meta_sample_count,
        and freq_index_dict annotations.
    """
    logger.info(
        "Subtracting removed samples' frequency statistics from base frequency HT..."
    )

    # Join the removed frequency stats with the base frequency HT to use
    # merge_freq_arrays
    joined_ht = base_freq_ht.annotate(
        removed_freq=removed_freq_stats[base_freq_ht.key].freq
    )

    merged_freq, merged_meta, sample_counts = merge_freq_arrays(
        [joined_ht.freq, joined_ht.removed_freq],
        [base_freq_ht.freq_meta, removed_freq_stats.freq_meta],
        operation="diff",
        count_arrays={
            "counts": [
                base_freq_ht.freq_meta_sample_count,
                removed_freq_stats.freq_meta_sample_count,
            ]
        },
    )

    # Update the frequency HT with the merged results
    modified_ht = joined_ht.select(freq=merged_freq).select_globals(
        freq_meta=merged_meta,
        freq_meta_sample_count=sample_counts["counts"],
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(merged_meta)),
    )

    return modified_ht


def subtract_age_histograms(
    base_age_hist_ht: hl.Table,
    removed_age_hist_stats: hl.Table,
) -> hl.Table:
    """
    Subtract age histogram statistics of removed samples from base age histogram HT.

    :param base_age_hist_ht: Base age histogram HT (v4).
    :param removed_age_hist_stats: Age histogram statistics for removed samples.
    :return: Modified age histogram HT with removed samples' stats subtracted.
    """
    logger.info(
        "Subtracting removed samples' age histogram statistics from base age histogram HT..."
    )

    joined_ht = base_age_hist_ht.annotate(
        removed_age_hist_het=removed_age_hist_stats[base_age_hist_ht.key].age_hist_het,
        removed_age_hist_hom=removed_age_hist_stats[base_age_hist_ht.key].age_hist_hom,
    )

    modified_ht = joined_ht.annotate(
        age_hist_het=hl.if_else(
            hl.is_defined(joined_ht.removed_age_hist_het),
            merge_histograms(
                [joined_ht.age_hist_het, joined_ht.removed_age_hist_het],
                operation="diff",
            ),
            joined_ht.age_hist_het,
        ),
        age_hist_hom=hl.if_else(
            hl.is_defined(joined_ht.removed_age_hist_hom),
            merge_histograms(
                [joined_ht.age_hist_hom, joined_ht.removed_age_hist_hom],
                operation="diff",
            ),
            joined_ht.age_hist_hom,
        ),
    )

    return modified_ht.select("age_hist_het", "age_hist_hom")


def add_frequency_stats(
    base_freq_ht: hl.Table,
    added_freq_stats: hl.Table,
) -> hl.Table:
    """
    Add frequency statistics of samples to base frequency HT.

    :param base_freq_ht: Base frequency HT.
    :param added_freq_stats: Frequency statistics for samples to add.
    :return: Modified frequency HT with added samples' stats.
    """
    logger.info("Adding samples' frequency statistics to base frequency HT...")

    # Join the added frequency stats with the base frequency HT to use merge_freq_arrays
    joined_ht = base_freq_ht.annotate(
        added_freq=added_freq_stats[base_freq_ht.key].freq
    )

    merged_freq, merged_meta, sample_counts = merge_freq_arrays(
        [joined_ht.freq, joined_ht.added_freq],
        [base_freq_ht.freq_meta, added_freq_stats.freq_meta],
        operation="sum",
        count_arrays={
            "counts": [
                base_freq_ht.freq_meta_sample_count,
                added_freq_stats.freq_meta_sample_count,
            ]
        },
    )

    modified_ht = joined_ht.select(freq=merged_freq).select_globals(
        freq_meta=merged_meta,
        freq_meta_sample_count=sample_counts["counts"],
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(merged_meta)),
    )

    return modified_ht


def add_age_histograms(
    base_age_hist_ht: hl.Table,
    added_age_hist_stats: hl.Table,
) -> hl.Table:
    """
    Add age histogram statistics of samples to base age histogram HT.

    :param base_age_hist_ht: Base age histogram HT.
    :param added_age_hist_stats: Age histogram statistics for samples to add.
    :return: Modified age histogram HT with added samples' stats.
    """
    logger.info("Adding samples' age histogram statistics to base age histogram HT...")

    # Join the added age histogram stats with the base age histogram HT to use
    # merge_histograms
    joined_ht = base_age_hist_ht.annotate(
        added_age_hist_het=added_age_hist_stats[base_age_hist_ht.key].age_hist_het,
        added_age_hist_hom=added_age_hist_stats[base_age_hist_ht.key].age_hist_hom,
    )

    modified_ht = joined_ht.annotate(
        age_hist_het=hl.if_else(
            hl.is_defined(joined_ht.added_age_hist_het),
            merge_histograms(
                [joined_ht.age_hist_het, joined_ht.added_age_hist_het],
                operation="sum",
            ),
            joined_ht.age_hist_het,
        ),
        age_hist_hom=hl.if_else(
            hl.is_defined(joined_ht.added_age_hist_hom),
            merge_histograms(
                [joined_ht.age_hist_hom, joined_ht.added_age_hist_hom],
                operation="sum",
            ),
            joined_ht.age_hist_hom,
        ),
    )

    return modified_ht.select("age_hist_het", "age_hist_hom")


def process_gnomad_dataset(
    test: bool = False,
    overwrite: bool = False,
) -> Tuple[hl.Table, hl.Table]:
    """
    Process gnomAD dataset to update v4 frequency HT for v5 changes.

    This function performs differential analysis on the gnomAD dataset by:
    1. Loading v4 frequency HT as the base
    2. Identifying samples to remove due to relatedness
    3. Identifying samples with ancestry changes between v4 and v5
    4. Subtracting frequency stats for removed samples
    5. Handling ancestry changes by subtracting with v4 ancestry and adding with v5 ancestry
    6. Writing updated frequency and age histogram HTs

    The function uses v4 downsamplings for compatibility with v4 frequency HT structure.

    :param test: Whether to run in test mode. If True, filters v4 freq HT to first 2 partitions.
    :param overwrite: Whether to overwrite existing output files.
    :return: Tuple of (updated frequency HT, updated age histogram HT) for gnomAD dataset.
        Both HTs are compatible with v4 structure for differential analysis.
    """
    logger.info("Processing gnomAD dataset...")

    v4_freq_ht = get_v4_freq(data_type="exomes").ht()

    if test:
        logger.info("Filtering v4 freq HT to the first 2 partitions for testing...")
        v4_freq_ht = v4_freq_ht._filter_partitions(range(2))

    logger.info(
        "Removing gnomAD samples due to relatedness and updating frequency and histogram HT..."
    )
    samples_to_remove_relatedness = identify_samples_to_remove_for_relatedness(
        test=test
    )
    logger.info("Calculating frequency stats for removed gnomAD samples...")
    removed_freq_stats = calculate_frequency_for_samples(
        samples_to_remove_relatedness, test=test
    )
    removed_age_hist_stats = calculate_age_histograms_for_samples(
        samples_to_remove_relatedness, test=test
    )

    freq_ht_after_relatedness = subtract_frequency_stats(v4_freq_ht, removed_freq_stats)
    age_hist_ht_after_relatedness = subtract_age_histograms(
        v4_freq_ht, removed_age_hist_stats
    )

    logger.info("Handling ancestry changes...")
    samples_with_ancestry_changes = identify_samples_with_ancestry_changes(test=test)

    logger.info(
        "Calculating frequency stats for samples with ancestry changes using v4 ancestry..."
    )
    ancestry_changed_freq_stats_v4 = calculate_frequency_for_samples(
        samples_with_ancestry_changes, test=test, use_v5_ancestry=False
    )
    ancestry_changed_age_hist_stats_v4 = calculate_age_histograms_for_samples(
        samples_with_ancestry_changes, test=test
    )

    freq_ht_after_ancestry_removal = subtract_frequency_stats(
        freq_ht_after_relatedness, ancestry_changed_freq_stats_v4
    )
    age_hist_ht_after_ancestry_removal = subtract_age_histograms(
        age_hist_ht_after_relatedness, ancestry_changed_age_hist_stats_v4
    )

    logger.info(
        "Calculating frequency stats for samples with ancestry changes using v5 ancestry..."
    )
    ancestry_changed_freq_stats_v5 = calculate_frequency_for_samples(
        samples_with_ancestry_changes, test=test, use_v5_ancestry=True
    )
    ancestry_changed_age_hist_stats_v5 = calculate_age_histograms_for_samples(
        samples_with_ancestry_changes, test=test
    )

    freq_ht = add_frequency_stats(
        freq_ht_after_ancestry_removal, ancestry_changed_freq_stats_v5
    )

    age_hist_ht = add_age_histograms(
        age_hist_ht_after_ancestry_removal, ancestry_changed_age_hist_stats_v5
    )

    return freq_ht, age_hist_ht


def process_aou_dataset(
    test: bool = False,
    overwrite: bool = False,
) -> Tuple[hl.Table, hl.Table]:
    """
    Process All of Us dataset for frequency calculations and age histograms.

    :param test: Whether to run in test mode.
    :param overwrite: Whether to overwrite existing results.
    :return: Tuple of (Frequency Table for AoU dataset, Age histogram Table for AoU dataset).
    """
    logger.info("Processing All of Us dataset...")

    check_resource_existence(
        get_aou_vds(remove_hard_filtered_samples=True, test=test), "AoU VDS"
    )

    aou_vds = get_aou_vds(remove_hard_filtered_samples=True, test=test)

    logger.info("Resolving sample ID collisions...")
    sample_collisions = sample_id_collisions.ht()

    aou_reference_mt = add_project_prefix_to_sample_collisions(
        t=aou_vds.reference_data,
        sample_collisions=sample_collisions,
        project="aou",
    )
    aou_variant_mt = add_project_prefix_to_sample_collisions(
        t=aou_vds.variant_data,
        sample_collisions=sample_collisions,
        project="aou",
    )

    aou_vds = hl.vds.VariantDataset(
        reference_data=aou_reference_mt, variant_data=aou_variant_mt
    )

    aou_meta_ht = project_meta.ht()
    aou_relatedness_ht = related_samples_to_drop(test=test, release=True).ht()

    aou_samples_to_remove = aou_meta_ht.annotate(
        will_be_dropped=aou_relatedness_ht[aou_meta_ht.key].s.is_defined()
    )

    aou_removed_sample_ids = aou_samples_to_remove.filter(
        aou_samples_to_remove.will_be_dropped
    ).s.collect()

    aou_vds_filtered = hl.vds.filter_samples(aou_vds, aou_removed_sample_ids, keep=True)

    aou_mt = hl.vds.to_dense_mt(aou_vds_filtered)

    aou_freq_ht = annotate_freq(
        aou_mt,
        sex_expr=aou_mt.meta.sex_karyotype,
        gen_anc_expr=aou_mt.meta.pop,
        additional_strata_expr=[
            {"data_set": "aou"},
            {"data_set": "aou", "gen_anc": aou_mt.meta.pop},
        ],
        # Use v4 downsamplings for compatibility with v4 frequency HT when doing differential analysis
        # annotate_freq will automatically filter out downsamplings larger than
        # available sample counts
        downsamplings=hl.literal(DOWNSAMPLINGS["v4"]),
    )

    aou_age_hist_expr = age_hists_expr(
        aou_mt.GT,
        aou_mt.meta.age,
        aou_mt.meta.sex_karyotype,
    )

    aou_age_hist_ht = aou_mt.aggregate_rows(
        age_hist_het=hl.agg.explode(aou_age_hist_expr.age_hist_het),
        age_hist_hom=hl.agg.explode(aou_age_hist_expr.age_hist_hom),
    )

    return aou_freq_ht, aou_age_hist_ht


def merge_gnomad_and_aou_frequencies(
    gnomad_freq_ht: hl.Table,
    aou_freq_ht: hl.Table,
    gnomad_age_hist_ht: hl.Table,
    aou_age_hist_ht: hl.Table,
    test: bool = False,
) -> Tuple[hl.Table, hl.Table]:
    """
    Merge frequency data and age histograms from gnomAD and All of Us datasets.

    :param gnomad_freq_ht: Frequency Table for gnomAD.
    :param aou_freq_ht: Frequency Table for AoU.
    :param gnomad_age_hist_ht: Age histogram Table for gnomAD.
    :param aou_age_hist_ht: Age histogram Table for AoU.
    :param test: Whether to run in test mode.
    :return: Tuple of (Merged frequency Table, Merged age histogram Table).
    """
    logger.info("Merging frequency data and age histograms from both datasets...")

    joined_freq_ht = gnomad_freq_ht.annotate(
        aou_freq=aou_freq_ht[gnomad_freq_ht.key].freq
    )

    merged_freq, merged_meta, sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.aou_freq],
        [gnomad_freq_ht.freq_meta, aou_freq_ht.freq_meta],
        operation="sum",
        count_arrays={
            "counts": [
                gnomad_freq_ht.freq_meta_sample_count,
                aou_freq_ht.freq_meta_sample_count,
            ]
        },
    )

    merged_freq_ht = joined_freq_ht.select(freq=merged_freq).select_globals(
        freq_meta=merged_meta,
        freq_meta_sample_count=sample_counts["counts"],
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(merged_meta)),
    )

    joined_age_hist_ht = gnomad_age_hist_ht.annotate(
        aou_age_hist_het=aou_age_hist_ht[gnomad_age_hist_ht.key].age_hist_het,
        aou_age_hist_hom=aou_age_hist_ht[gnomad_age_hist_ht.key].age_hist_hom,
    )

    merged_age_hist_ht = joined_age_hist_ht.annotate(
        age_hist_het=hl.if_else(
            hl.is_defined(joined_age_hist_ht.aou_age_hist_het),
            merge_histograms(
                [joined_age_hist_ht.age_hist_het, joined_age_hist_ht.aou_age_hist_het],
                operation="sum",
            ),
            joined_age_hist_ht.age_hist_het,
        ),
        age_hist_hom=hl.if_else(
            hl.is_defined(joined_age_hist_ht.aou_age_hist_hom),
            merge_histograms(
                [joined_age_hist_ht.age_hist_hom, joined_age_hist_ht.aou_age_hist_hom],
                operation="sum",
            ),
            joined_age_hist_ht.age_hist_hom,
        ),
    )

    return merged_age_hist_ht.select("age_hist_het", "age_hist_hom"), merged_freq_ht


def main(args):
    """Generate v5 frequency data."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/generate_frequency.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")
    hl._set_flags(use_ssa_logs="1")

    test = args.test
    overwrite = args.overwrite

    try:
        if args.process_gnomad:
            logger.info("Processing gnomAD dataset...")

            gnomad_freq_resource = get_freq_ht(test=test, data_set="gnomad")
            gnomad_age_hist_resource = get_age_hist_ht(test=test, data_set="gnomad")

            check_resource_existence(
                gnomad_freq_resource, "gnomAD frequency HT", overwrite=overwrite
            )
            check_resource_existence(
                gnomad_age_hist_resource, "gnomAD age histogram HT", overwrite=overwrite
            )

            gnomad_freq_ht, gnomad_age_hist_ht = process_gnomad_dataset(
                test=test, overwrite=overwrite
            )

            logger.info(
                f"Writing gnomAD frequency HT to {gnomad_freq_resource.path}..."
            )
            gnomad_freq_ht.write(gnomad_freq_resource.path, overwrite=overwrite)

            logger.info(
                f"Writing gnomAD age histogram HT to {gnomad_age_hist_resource.path}..."
            )
            gnomad_age_hist_ht.write(gnomad_age_hist_resource.path, overwrite=overwrite)

        if args.process_aou:
            logger.info("Processing All of Us dataset...")

            aou_freq_resource = get_freq_ht(test=test, data_set="aou")
            aou_age_hist_resource = get_age_hist_ht(test=test, data_set="aou")

            check_resource_existence(
                aou_freq_resource, "AoU frequency HT", overwrite=overwrite
            )
            check_resource_existence(
                aou_age_hist_resource, "AoU age histogram HT", overwrite=overwrite
            )

            aou_freq_ht, aou_age_hist_ht = process_aou_dataset(
                test=test, overwrite=overwrite
            )

            logger.info(f"Writing AoU frequency HT to {aou_freq_resource.path}...")
            aou_freq_ht.write(aou_freq_resource.path, overwrite=overwrite)

            logger.info(
                f"Writing AoU age histogram HT to {aou_age_hist_resource.path}..."
            )
            aou_age_hist_ht.write(aou_age_hist_resource.path, overwrite=overwrite)

        if args.merge_datasets:
            logger.info(
                "Merging frequency data and age histograms from both datasets..."
            )

            merged_freq_resource = get_freq_ht(test=test)
            merged_age_hist_resource = get_age_hist_ht(test=test)

            check_resource_existence(
                merged_freq_resource, "merged frequency HT", overwrite=overwrite
            )
            check_resource_existence(
                merged_age_hist_resource, "merged age histogram HT", overwrite=overwrite
            )

            check_resource_existence(
                get_freq_ht(test=test, data_set="gnomad"), "gnomAD frequency HT"
            )
            check_resource_existence(
                get_freq_ht(test=test, data_set="aou"), "AoU frequency HT"
            )
            check_resource_existence(
                get_age_hist_ht(test=test, data_set="gnomad"),
                "gnomAD age histogram HT",
            )
            check_resource_existence(
                get_age_hist_ht(test=test, data_set="aou"),
                "AoU age histogram HT",
            )

            gnomad_freq_ht = get_freq_ht(test=test, data_set="gnomad").ht()
            aou_freq_ht = get_freq_ht(test=test, data_set="aou").ht()
            gnomad_age_hist_ht = get_age_hist_ht(test=test, data_set="gnomad").ht()
            aou_age_hist_ht = get_age_hist_ht(test=test, data_set="aou").ht()

            merged_freq_ht, merged_age_hist_ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq_ht,
                aou_freq_ht,
                gnomad_age_hist_ht,
                aou_age_hist_ht,
                test=test,
            )

            logger.info(
                f"Writing merged frequency HT to {merged_freq_resource.path}..."
            )
            merged_freq_ht.write(merged_freq_resource.path, overwrite=overwrite)

            logger.info(
                f"Writing merged age histogram HT to {merged_age_hist_resource.path}..."
            )
            merged_age_hist_ht.write(merged_age_hist_resource.path, overwrite=overwrite)

    finally:
        hl.copy_log(get_logging_path("v5_frequency"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Filter to the first 2 partitions for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--process-gnomad",
        help="Process gnomAD dataset for frequency calculations.",
        action="store_true",
    )
    parser.add_argument(
        "--process-aou",
        help="Process All of Us dataset for frequency calculations.",
        action="store_true",
    )
    parser.add_argument(
        "--merge-datasets",
        help="Merge frequency data from both gnomAD and AoU datasets.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
