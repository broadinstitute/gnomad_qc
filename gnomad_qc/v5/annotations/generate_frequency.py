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
from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_adj,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    build_freq_stratification_list,
    faf_expr,
    gen_anc_faf_max_expr,
    generate_freq_group_membership_array,
    get_adj_expr,
    grpmax_expr,
    merge_freq_arrays,
    merge_histograms,
)
from gnomad.utils.release import make_freq_index_dict_from_meta
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
    check_resource_existence,
)
from gnomad_qc.v4.resources.annotations import get_freq as get_v4_freq
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v5.resources.annotations import (
    gnomad_group_membership,  # gnomAD-specific group membership
)
from gnomad_qc.v5.resources.annotations import get_age_hist_ht, get_freq_ht
from gnomad_qc.v5.resources.annotations import (
    group_membership as aou_group_membership,  # AoU-specific group membership
)
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import V5_DOWNSAMPLINGS, WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import (
    consent_samples_to_drop,
    project_meta,
    sample_id_collisions,
)
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

    # Use existing gnomAD group membership table and filter to sample subset
    logger.info("Loading gnomAD group membership table for frequency stratification...")
    group_membership_ht = gnomad_group_membership(test=test).ht()

    # Filter group membership to only the samples we're processing
    sample_ids = set(mt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(sample_ids).contains(group_membership_ht.s)
    )

    # Annotate MatrixTable with group membership for annotate_freq to use
    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership
    )

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
        downsamplings=DOWNSAMPLINGS["v4"],
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
    Process gnomAD dataset to update v4 frequency HT by removing consent withdrawal samples.

    This function performs frequency adjustment by:
    1. Loading v4 frequency HT as the base
    2. Loading v4 genome VDS with consent samples to drop
    3. Filtering to consent withdrawal samples and release samples only
    4. Calculating frequencies for consent withdrawal samples using v4 logic
    5. Subtracting those frequencies from v4 frequency HT
    6. Computing FAF, grpmax, gen_anc_faf_max, and inbreeding coefficient

    :param test: Whether to run in test mode. If True, filters v4 vds to first 2 partitions.
    :param overwrite: Whether to overwrite existing output files.
    :return: Tuple of (updated frequency HT with FAF/grpmax annotations, updated age histogram HT) for gnomAD dataset.
    """
    logger.info("Processing gnomAD dataset for consent withdrawals...")

    v4_freq_ht = get_v4_freq(data_type="exomes").ht()
    vds = get_gnomad_v4_vds(
        data_type="genomes",
        test=test,
        release_only=True,
        annotate_meta=True,
        n_partitions=2 if test else None,
    )

    consent_samples_ht = consent_samples_to_drop.ht()
    consent_samples_list = consent_samples_ht.s.collect()

    logger.info(
        "Filtering VDS to consent withdrawal samples and preparing for frequency calculation..."
    )
    vds = hl.vds.filter_samples(vds, consent_samples_list, keep=True)

    # Simplified approach following v4 pattern
    vmt = vds.variant_data
    vmt = vmt.select_cols(
        pop=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        age=vmt.meta.project_meta.age,
    )

    vmt = vmt.select_entries("LA", "LAD", "DP", "GQ", "LGT")

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info(
        "Densifying VDS and calculating frequencies for consent withdrawal samples..."
    )
    mt = hl.vds.to_dense_mt(vds)

    # Following v4 pattern - adj before sex ploidy adjustment
    mt = annotate_adj(mt)

    gt_expr = adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype)
    mt = mt.annotate_entries(GT=gt_expr)

    # Checkpoint after expensive MT operations to avoid recomputation
    mt = mt.checkpoint(new_temp_file("consent_samples_mt", "mt"))

    # Use existing gnomAD group membership table and filter to consent samples
    logger.info(
        "Loading gnomAD group membership table for consent sample frequencies..."
    )
    group_membership_ht = gnomad_group_membership(test=test).ht()

    # Filter group membership to only the consent samples we're processing
    consent_sample_ids = set(mt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(consent_sample_ids).contains(group_membership_ht.s)
    )

    # Annotate MatrixTable with group membership for annotate_freq to use
    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership
    )

    consent_freq_ht = annotate_freq(
        mt,
        sex_expr=mt.sex_karyotype,
        gen_anc_expr=mt.pop,
        annotate_mt=False,
    )

    # Checkpoint expensive frequency calculation since it will be used for merging
    consent_freq_ht = consent_freq_ht.checkpoint(new_temp_file("consent_freq", "ht"))

    logger.info(
        "Subtracting consent withdrawal sample frequencies from v4 frequency table..."
    )

    # First, join consent frequencies onto the v4 frequency table
    joined_freq_ht = v4_freq_ht.annotate(
        consent_freq=consent_freq_ht[v4_freq_ht.key].freq
    )

    # Annotate globals from consent frequency table (extract as literals to
    # avoid source mismatch)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=hl.literal(hl.eval(consent_freq_ht.freq_meta)),
        consent_freq_meta_sample_count=hl.literal(
            hl.eval(consent_freq_ht.freq_meta_sample_count)
        ),
    )

    # Using merge_freq_arrays for proper frequency subtraction (following v4 pattern)
    updated_freq_expr, updated_freq_meta, updated_sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.consent_freq],
        [joined_freq_ht.freq_meta, joined_freq_ht.consent_freq_meta],
        operation="diff",
        count_arrays={
            "freq_meta_sample_count": [
                joined_freq_ht.freq_meta_sample_count,
                joined_freq_ht.consent_freq_meta_sample_count,
            ]
        },
    )

    updated_freq_ht = joined_freq_ht.annotate(freq=updated_freq_expr)
    updated_freq_ht = updated_freq_ht.annotate_globals(
        freq_meta=hl.literal(hl.eval(updated_freq_meta)),
        freq_meta_sample_count=hl.literal(
            hl.eval(updated_sample_counts["freq_meta_sample_count"])
        ),
    )

    # Checkpoint after merge_freq_arrays to materialize complex expressions
    # before FAF calculations
    updated_freq_ht = updated_freq_ht.checkpoint(new_temp_file("merged_freq", "ht"))

    logger.info(
        "Calculating and subtracting age histograms for consent withdrawal samples..."
    )

    mt_with_age_hists = mt.annotate_rows(
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age)
    )
    consent_age_hist_ht = mt_with_age_hists.rows().select("age_hists")

    v4_age_hist_ht = get_age_hist_ht(test=test, data_set="gnomad").ht()

    updated_age_hist_ht = v4_age_hist_ht.annotate(
        age_hist_het=merge_histograms(
            [
                v4_age_hist_ht.age_hist_het,
                consent_age_hist_ht[v4_age_hist_ht.key].age_hists.age_hist_het,
            ],
            operation="diff",
        ),
        age_hist_hom=merge_histograms(
            [
                v4_age_hist_ht.age_hist_hom,
                consent_age_hist_ht[v4_age_hist_ht.key].age_hists.age_hist_hom,
            ],
            operation="diff",
        ),
    )

    # Calculate FAF, grpmax, and other post-processing annotations
    updated_freq_ht = _calculate_faf_and_grpmax_annotations(updated_freq_ht, test=test)

    return updated_freq_ht, updated_age_hist_ht


def _calculate_faf_and_grpmax_annotations(updated_freq_ht, test=False):
    """
    Calculate FAF, grpmax, gen_anc_faf_max, and inbreeding coefficient annotations.

    This function handles the complex post-processing annotations that are added
    to frequency tables after the core frequency calculations are complete.

    :param updated_freq_ht: Frequency table after consent withdrawal subtraction
    :param test: Whether this is a test run (skips FAF/grpmax in test mode)
    :return: Updated frequency table with FAF/grpmax annotations
    """
    if test:
        logger.info("Skipping FAF/grpmax calculations in test mode...")
        return updated_freq_ht

    logger.info("Computing FAF, grpmax, gen_anc_faf_max, and InbreedingCoeff...")

    # Change 'gen_anc' keys to 'pop' for faf computation
    freq_meta_for_faf = hl.literal(
        [
            {("pop" if k == "gen_anc" else k): m[k] for k in m}
            for m in hl.eval(updated_freq_ht.freq_meta)
        ]
    )

    # Calculate FAF (Filtering Allele Frequency)
    faf, faf_meta = faf_expr(
        updated_freq_ht.freq,
        freq_meta_for_faf,
        updated_freq_ht.locus,
        GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v4"],
    )

    # Calculate grpmax (group maximum frequency)
    grpmax = grpmax_expr(
        updated_freq_ht.freq,
        freq_meta_for_faf,
        GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v4"],
    )

    # Annotate grpmax with corresponding FAF95 values
    grpmax = grpmax.annotate(
        faf95=faf[
            hl.literal(faf_meta).index(lambda y: y.values() == ["adj", grpmax.gen_anc])
        ].faf95,
    )

    # Add all annotations to the frequency table
    updated_freq_ht = updated_freq_ht.annotate(
        faf=faf,
        grpmax=grpmax,
        gen_anc_faf_max=gen_anc_faf_max_expr(faf, faf_meta),
        inbreeding_coeff=bi_allelic_site_inbreeding_expr(
            callstats_expr=updated_freq_ht.freq[1]
        ),
    )

    # Checkpoint after expensive FAF/grpmax calculations
    updated_freq_ht = updated_freq_ht.checkpoint(new_temp_file("freq_with_faf", "ht"))

    # Change 'pop' keys back to 'gen_anc' for consistency
    final_freq_meta = hl.literal(
        [
            {("gen_anc" if k == "pop" else k): m[k] for k in m}
            for m in hl.eval(updated_freq_ht.freq_meta)
        ]
    )
    final_faf_meta = hl.literal(
        [{("gen_anc" if k == "pop" else k): m[k] for k in m} for m in faf_meta]
    )

    # Update globals with final metadata
    updated_freq_ht = updated_freq_ht.annotate_globals(
        freq_meta=final_freq_meta,
        faf_meta=final_faf_meta,
        faf_index_dict=make_freq_index_dict_from_meta(hl.literal(final_faf_meta)),
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(final_freq_meta)),
    )

    return updated_freq_ht


def _resolve_aou_sample_collisions(aou_vds, test=False):
    """
    Resolve sample ID collisions for AoU VDS by adding project prefixes.

    :param aou_vds: AoU VariantDataset
    :param test: Whether this is a test run
    :return: VDS with resolved sample collisions
    """
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

    return hl.vds.VariantDataset(
        reference_data=aou_reference_mt, variant_data=aou_variant_mt
    )


def _filter_aou_related_samples(aou_vds, test=False):
    """
    Filter related samples from AoU VDS.

    :param aou_vds: AoU VariantDataset
    :param test: Whether this is a test run
    :return: Filtered VDS with related samples removed
    """
    aou_meta_ht = project_meta.ht()
    aou_relatedness_ht = related_samples_to_drop(test=test, release=True).ht()

    aou_samples_to_remove = aou_meta_ht.annotate(
        will_be_dropped=hl.is_defined(aou_relatedness_ht[aou_meta_ht.key])
    )

    aou_removed_sample_ids = aou_samples_to_remove.filter(
        aou_samples_to_remove.will_be_dropped
    ).s.collect()

    logger.info(
        f"Filtering {len(aou_removed_sample_ids)} related samples from AoU VDS..."
    )
    return hl.vds.filter_samples(aou_vds, aou_removed_sample_ids, keep=True)


def _prepare_aou_variant_data(aou_vds_filtered, test=False):
    """
    Prepare AoU variant data for frequency calculation.

    :param aou_vds_filtered: Filtered AoU VDS
    :param test: Whether this is a test run
    :return: Prepared variant MatrixTable
    """
    aou_variant_mt = aou_vds_filtered.variant_data

    if test:
        # In test mode, adapt to the test VDS metadata structure
        aou_variant_mt = aou_variant_mt.annotate_cols(
            sex_karyotype=aou_variant_mt.meta.sex_imputation.sex_karyotype,
            pop=aou_variant_mt.meta.population_inference.pop,
            age=aou_variant_mt.meta.project_meta.age,
        )
    else:
        # Production mode expects direct metadata fields
        aou_variant_mt = aou_variant_mt.annotate_cols(
            sex_karyotype=aou_variant_mt.meta.sex_karyotype,
            pop=aou_variant_mt.meta.pop,
            age=aou_variant_mt.meta.age,
        )

    # Rename LGT to GT and LAD to AD for compatibility with annotate_freq and
    # annotate_adj
    aou_variant_mt = aou_variant_mt.annotate_entries(
        GT=aou_variant_mt.LGT, AD=aou_variant_mt.LAD
    )

    # Add adj annotation required by annotate_freq
    aou_variant_mt = annotate_adj(aou_variant_mt)

    # Checkpoint after metadata annotation
    return aou_variant_mt.checkpoint(new_temp_file("aou_variant_mt", "mt"))


def _calculate_aou_variant_frequencies(aou_variant_mt, test=False):
    """
    Calculate frequencies from AoU variant data (gives correct AC, hom but incorrect AN).

    :param aou_variant_mt: Prepared variant MatrixTable
    :param test: Whether to use test resources
    :return: Frequency Table with AC and hom counts from variant data
    """
    # Use existing AoU group membership table and filter to variant samples
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht = aou_group_membership(test=test).ht()

    # Filter group membership to only the AoU samples we're processing
    aou_sample_ids = set(aou_variant_mt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(aou_sample_ids).contains(group_membership_ht.s)
    )

    # Annotate MatrixTable with group membership for annotate_freq to use
    aou_variant_mt = aou_variant_mt.annotate_cols(
        group_membership=group_membership_ht[aou_variant_mt.col_key].group_membership
    )

    return annotate_freq(
        aou_variant_mt,
        sex_expr=aou_variant_mt.sex_karyotype,
        gen_anc_expr=aou_variant_mt.pop,
        additional_strata_expr=[
            {"data_set": "aou"},
            {"data_set": "aou", "gen_anc": aou_variant_mt.pop},
        ],
        downsamplings=DOWNSAMPLINGS["v4"],
        annotate_mt=False,
    )


def _prepare_aou_reference_data(aou_vds_filtered, test=False):
    """
    Prepare AoU reference data for AN calculation.

    :param aou_vds_filtered: Filtered AoU VDS
    :param test: Whether this is a test run
    :return: Prepared reference MatrixTable
    """
    aou_reference_mt = aou_vds_filtered.reference_data

    if test:
        # In test mode, adapt to the test VDS metadata structure
        aou_reference_mt = aou_reference_mt.annotate_cols(
            sex_karyotype=aou_reference_mt.meta.sex_imputation.sex_karyotype,
            pop=aou_reference_mt.meta.population_inference.pop,
        )
    else:
        # Production mode expects direct metadata fields
        aou_reference_mt = aou_reference_mt.annotate_cols(
            sex_karyotype=aou_reference_mt.meta.sex_karyotype,
            pop=aou_reference_mt.meta.pop,
        )

    # Rename LGT to GT for compatibility with annotate_freq
    # Reference data doesn't have LAD, GQ, DP, so create placeholders for annotate_adj
    aou_reference_mt = aou_reference_mt.annotate_entries(
        GT=aou_reference_mt.LGT,
        AD=hl.missing(hl.tarray(hl.tint32)),  # Placeholder AD for annotate_adj
        GQ=hl.missing(hl.tint32),  # Placeholder GQ for annotate_adj
        DP=hl.missing(hl.tint32),  # Placeholder DP for annotate_adj
    )

    # Add adj annotation required by annotate_freq
    aou_reference_mt = annotate_adj(aou_reference_mt)

    # Checkpoint reference MT
    aou_reference_mt = aou_reference_mt.checkpoint(
        new_temp_file("aou_reference_mt", "mt")
    )

    # Reference data needs alleles field for annotate_freq - add dummy alleles
    aou_reference_mt = aou_reference_mt.annotate_rows(
        alleles=["A", "T"]
    )  # Add dummy alleles
    aou_reference_mt = aou_reference_mt.key_rows_by(
        "locus", "alleles"
    )  # Re-key with alleles

    return aou_reference_mt


def _calculate_aou_reference_an(aou_reference_mt, test=False):
    """
    Calculate correct AN values from AoU reference data (all sites).

    :param aou_reference_mt: Prepared reference MatrixTable
    :param test: Whether to use test resources
    :return: Frequency Table with correct AN values from all sites
    """
    # Use existing AoU group membership table and filter to reference samples
    logger.info("Loading AoU group membership table for reference AN calculation...")
    group_membership_ht = aou_group_membership(test=test).ht()

    # Filter group membership to only the AoU samples we're processing
    aou_ref_sample_ids = set(aou_reference_mt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(aou_ref_sample_ids).contains(group_membership_ht.s)
    )

    # Annotate MatrixTable with group membership for annotate_freq to use
    aou_reference_mt = aou_reference_mt.annotate_cols(
        group_membership=group_membership_ht[aou_reference_mt.col_key].group_membership
    )

    return annotate_freq(
        aou_reference_mt,
        sex_expr=aou_reference_mt.sex_karyotype,
        gen_anc_expr=aou_reference_mt.pop,
        additional_strata_expr=[
            {"data_set": "aou"},
            {"data_set": "aou", "gen_anc": aou_reference_mt.pop},
        ],
        downsamplings=DOWNSAMPLINGS["v4"],
        annotate_mt=False,
    )


def _correct_aou_frequencies(aou_freq_ht, aou_ref_freq_ht):
    """
    Correct AoU frequencies by combining AC from variant data with AN from reference data.

    :param aou_freq_ht: Frequency table with AC/hom from variant data
    :param aou_ref_freq_ht: Frequency table with AN from reference data
    :return: Corrected frequency table with proper AF calculation
    """
    logger.info("Correcting AN values and recalculating allele frequencies...")
    corrected_an_ht = aou_ref_freq_ht.select(
        correct_freq=hl.map(
            lambda idx_freq_tuple: idx_freq_tuple[1].annotate(
                # Keep original AC and hom, but update AN from all-sites data
                AN=aou_ref_freq_ht.freq[idx_freq_tuple[0]].AN,
                # Recalculate AF with correct AN
                AF=hl.if_else(
                    aou_ref_freq_ht.freq[idx_freq_tuple[0]].AN > 0,
                    hl.float64(idx_freq_tuple[1].AC)
                    / hl.float64(aou_ref_freq_ht.freq[idx_freq_tuple[0]].AN),
                    0.0,
                ),
            ),
            hl.enumerate(aou_freq_ht[aou_ref_freq_ht.key].freq),
        )
    )

    # Join corrected AN back to original frequency table
    corrected_freq_ht = aou_freq_ht.annotate(
        freq=corrected_an_ht[aou_freq_ht.key].correct_freq
    )

    # Checkpoint after frequency correction
    return corrected_freq_ht.checkpoint(new_temp_file("aou_corrected_freq", "ht"))


def _calculate_aou_age_histograms(aou_vds_filtered, test=False):
    """
    Calculate age histograms for AoU dataset.

    :param aou_vds_filtered: Filtered AoU VDS
    :param test: Whether this is a test run
    :return: Age histogram Table
    """
    logger.info("Computing age histograms for AoU dataset...")

    if test:
        logger.info("Creating minimal age histogram table for test mode...")
        # Create a simple age histogram table that matches expected structure
        return hl.Table.parallelize(
            [
                {
                    "locus": hl.locus("chr22", 10000000),
                    "alleles": ["A", "T"],
                    "age_hist_het": [10, 15, 20, 25, 15, 5],  # Simple histogram bins
                    "age_hist_hom": [5, 8, 10, 12, 8, 2],
                }
            ],
            hl.tstruct(
                locus=hl.tlocus("GRCh38"),
                alleles=hl.tarray(hl.tstr),
                age_hist_het=hl.tarray(hl.tint64),
                age_hist_hom=hl.tarray(hl.tint64),
            ),
            key=["locus", "alleles"],
        )
    else:
        # Production mode uses proper age_hists_expr
        # For age histograms, we need the dense MT approach
        aou_mt_dense = hl.vds.to_dense_mt(aou_vds_filtered)
        # Dense MT still has LGT field, need to rename to GT for age_hists_expr
        aou_mt_dense = aou_mt_dense.annotate_entries(GT=aou_mt_dense.LGT)
        aou_age_hist_expr = age_hists_expr(
            aou_mt_dense.GT,
            aou_mt_dense.meta.age,
            aou_mt_dense.meta.sex_karyotype,
        )
        return aou_mt_dense.aggregate_rows(
            age_hist_het=aou_age_hist_expr.age_hist_het,
            age_hist_hom=aou_age_hist_expr.age_hist_hom,
        )


def process_aou_dataset(
    test: bool = False,
    overwrite: bool = False,
) -> Tuple[hl.Table, hl.Table]:
    """
    Process All of Us dataset for frequency calculations and age histograms.

    This function efficiently processes the AoU VDS by:
    1. Computing allele counts (AC) and homozygote counts from variant_data (sparse)
    2. Computing allele numbers (AN) from reference_data (all sites including reference-only)
    3. Combining these to get accurate allele frequencies with correct denominators

    :param test: Whether to run in test mode.
    :param overwrite: Whether to overwrite existing results.
    :return: Tuple of (Frequency Table for AoU dataset, Age histogram Table for AoU dataset).
    """
    logger.info("Processing All of Us dataset...")

    # Load and validate AoU VDS
    check_resource_existence(
        get_aou_vds(remove_hard_filtered_samples=True, test=test), "AoU VDS"
    )
    aou_vds = get_aou_vds(remove_hard_filtered_samples=True, test=test)

    # Step 1: Resolve sample ID collisions
    aou_vds = _resolve_aou_sample_collisions(aou_vds, test=test)

    # Step 2: Filter related samples
    aou_vds_filtered = _filter_aou_related_samples(aou_vds, test=test)

    # Step 3: Calculate frequencies using efficient variant_data + all-sites AN approach
    logger.info(
        "Calculating efficient AoU frequencies using variant_data + all-sites AN..."
    )

    # Step 3a: Calculate AC/hom from variant_data
    aou_variant_mt = _prepare_aou_variant_data(aou_vds_filtered, test=test)
    aou_freq_ht = _calculate_aou_variant_frequencies(aou_variant_mt, test=test)

    # Step 3b: Calculate correct AN from reference_data (all sites)
    logger.info("Computing correct allele numbers (AN) from all sites...")
    aou_reference_mt = _prepare_aou_reference_data(aou_vds_filtered, test=test)
    aou_ref_freq_ht = _calculate_aou_reference_an(aou_reference_mt, test=test)

    # Step 3c: Correct frequencies by combining variant AC with reference AN
    aou_freq_ht = _correct_aou_frequencies(aou_freq_ht, aou_ref_freq_ht)

    # Step 4: Calculate age histograms
    aou_age_hist_ht = _calculate_aou_age_histograms(aou_vds_filtered, test=test)

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

    # Annotate globals from AoU frequency table (extract as literals to avoid
    # source mismatch)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        aou_freq_meta=hl.literal(hl.eval(aou_freq_ht.freq_meta)),
        aou_freq_meta_sample_count=hl.literal(
            hl.eval(aou_freq_ht.freq_meta_sample_count)
        ),
    )

    merged_freq, merged_meta, sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.aou_freq],
        [joined_freq_ht.freq_meta, joined_freq_ht.aou_freq_meta],
        operation="sum",
        count_arrays={
            "counts": [
                joined_freq_ht.freq_meta_sample_count,
                joined_freq_ht.aou_freq_meta_sample_count,
            ]
        },
    )

    merged_freq_ht = joined_freq_ht.select(freq=merged_freq).select_globals(
        freq_meta=hl.literal(hl.eval(merged_meta)),
        freq_meta_sample_count=hl.literal(hl.eval(sample_counts["counts"])),
        freq_index_dict=make_freq_index_dict_from_meta(
            hl.literal(hl.eval(merged_meta))
        ),
    )

    if test:
        logger.info(
            "Skipping age histogram merging in test mode due to simplified structure..."
        )
        # In test mode, just use gnomAD age histograms since AoU has simplified
        # structure
        merged_age_hist_ht = gnomad_age_hist_ht
    else:
        # Production mode: proper age histogram merging
        joined_age_hist_ht = gnomad_age_hist_ht.annotate(
            aou_age_hist_het=aou_age_hist_ht[gnomad_age_hist_ht.key].age_hist_het,
            aou_age_hist_hom=aou_age_hist_ht[gnomad_age_hist_ht.key].age_hist_hom,
        )

        merged_age_hist_ht = joined_age_hist_ht.annotate(
            age_hist_het=hl.if_else(
                hl.is_defined(joined_age_hist_ht.aou_age_hist_het),
                merge_histograms(
                    [
                        joined_age_hist_ht.age_hist_het,
                        joined_age_hist_ht.aou_age_hist_het,
                    ],
                    operation="sum",
                ),
                joined_age_hist_ht.age_hist_het,
            ),
            age_hist_hom=hl.if_else(
                hl.is_defined(joined_age_hist_ht.aou_age_hist_hom),
                merge_histograms(
                    [
                        joined_age_hist_ht.age_hist_hom,
                        joined_age_hist_ht.aou_age_hist_hom,
                    ],
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
