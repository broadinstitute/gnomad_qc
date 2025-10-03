"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies for samples that will be removed from gnomAD v4
in v5 due to relatedness filtering or ancestry changes. It uses a differential approach
to avoid densifying the full dataset.
"""

import argparse
import logging
from typing import Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    agg_by_strata,
    annotate_adj,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    gen_anc_faf_max_expr,
    get_adj_expr,
    grpmax_expr,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.release import make_freq_index_dict_from_meta
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.annotations import get_freq as get_v4_freq
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds
from gnomad_qc.v5.resources.annotations import (
    get_consent_ans,
    get_freq,
    get_group_membership,
)
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import (
    consent_samples_to_drop,
    project_meta,
    sample_id_collisions,
)
from gnomad_qc.v5.resources.sample_qc import related_samples_to_drop

# Constants
TEST_PARTITIONS = 2  # Number of partitions to use in test mode
TMP_DIR_DAYS = 4  # Number of days for temp directory retention

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def high_ab_het(
    entry: hl.StructExpression, col: hl.StructExpression
) -> hl.Int32Expression:
    """
    Determine if a call is considered a high allele balance heterozygous call.

    High allele balance heterozygous calls were introduced in certain GATK versions.
    Track how many calls appear at each site to correct them to homozygous
    alternate calls downstream in frequency calculations and histograms.

    Assumes the following annotations are present in `entry` struct:
        - GT
        - adj
        - _high_ab_het_ref

    Assumes the following annotations are present in `col` struct:
        - fixed_homalt_model

    :param entry: Entry struct.
    :param col: Column struct.
    :return: 1 if high allele balance heterozygous call, else 0.
    """
    return hl.int(
        entry.GT.is_het_ref()
        & entry.adj
        & ~col.fixed_homalt_model
        & entry._high_ab_het_ref
    )


def mt_hists_fields(mt: hl.MatrixTable) -> hl.StructExpression:
    """
    Annotate quality metrics histograms and age histograms onto MatrixTable.

    :param mt: Input MatrixTable.
    :return: Struct with quality metrics histograms and age histograms.
    """
    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    high_ab_gt_expr = hl.if_else(high_ab_het(mt, mt) == 1, hl.call(1, 1), mt.GT)

    return hl.struct(
        **qual_hist_expr(
            gt_expr=mt.GT,
            gq_expr=mt.GQ,
            dp_expr=mt.DP,
            adj_expr=mt.adj,
            ab_expr=mt._het_ab,
            split_adj_and_raw=True,
        ),
        high_ab_het_adjusted_ab_hists=qual_hist_expr(
            gt_expr=high_ab_gt_expr,
            adj_expr=mt.adj,
            ab_expr=mt._het_ab,
        ),
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        high_ab_het_adjusted_age_hists=age_hists_expr(mt.adj, high_ab_gt_expr, mt.age),
    )


def _prepare_consent_vds(
    v4_freq_ht: hl.Table, test: bool = False, runtime_test: bool = False
) -> hl.MatrixTable:
    """
    Load and prepare VDS for consent withdrawal sample processing.

    :param v4_freq_ht: v4 frequency table for AF annotation.
    :param test: Whether running in test mode.
    :param runtime_test: Whether to use test VDS.
    :return: Prepared MatrixTable with consent samples, split multiallelics, and annotations.
    """
    logger.info("Loading and preparing VDS for consent withdrawal samples...")

    vds = get_gnomad_v4_genomes_vds(
        test=runtime_test,
        remove_hard_filtered_samples=False,
        release_only=True,
        annotate_meta=True,
        filter_partitions=list(range(TEST_PARTITIONS)) if test else None,
    )

    consent_samples_ht = consent_samples_to_drop.ht()
    consent_samples_list = consent_samples_ht.s.collect()

    logger.info("Filtering VDS to consent withdrawal samples...")
    vds = hl.vds.filter_samples(vds, consent_samples_list, keep=True)

    # Prepare variant data with metadata
    vmt = vds.variant_data
    vmt = vmt.select_cols(
        pop=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        age=vmt.meta.project_meta.age,
    )
    # For genomes, fixed_homalt_model is always False since we apply v3-style correction to all samples
    # (following v3 and v4 genomes approach - no GATK version-based differentiation)
    vmt = vmt.annotate_cols(fixed_homalt_model=hl.bool(False))
    vmt = vmt.annotate_globals(
        age_distribution=vmt.aggregate_cols(hl.agg.hist(vmt.age, 30, 80, 10))
    )

    # Select required entries and annotate non_ref hets before split_multi
    logger.info("Selecting entries and annotating non_ref hets pre-split...")
    vmt = vmt.select_entries(
        "LA", "LAD", "DP", "GQ", "LGT", _het_non_ref=vmt.LGT.is_het_non_ref()
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = vds.checkpoint(new_temp_file("consent_samples_vds", "vds"))

    logger.info("Splitting multiallelics in gnomAD sample withdrawal VDS...")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    # Annotate with v4 frequencies and apply quality adjustments
    vmt = vds.variant_data
    vmt = vmt.annotate_rows(v4_af=v4_freq_ht[vmt.row_key].freq[0].AF)

    # This follows the v3/v4 genomes workflow for adj and sex adjusted genotypes.
    # The correct thing to do here is to adjust sex ploidy before determining the adj
    # annotation because haploid GTs have different adj filtering criteria, but the
    # option to adjust ploidy after adj is included for consistency with v3.1, where we
    # added the adj annotation before adjusting for sex ploidy.
    logger.info("Computing sex adjusted genotypes and quality annotations...")
    vmt = vmt.annotate_entries(
        adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD),
    )
    ab_cutoff = 0.9
    ab_expr = vmt.AD[1] / vmt.DP
    vmt = vmt.select_entries(
        "AD",
        "DP",
        "GQ",
        "_het_non_ref",
        "adj",
        GT=adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype),
        _het_ab=ab_expr,
        _high_ab_het_ref=(ab_expr > ab_cutoff) & ~vmt._het_non_ref,
    )

    # We set use_v3_1_correction to True to mimic the v4 genomes approach.
    logger.info("Applying v4 genomes hom alt depletion fix...")
    gt_with_depletion = hom_alt_depletion_fix(
        vmt.GT,
        het_non_ref_expr=vmt._het_non_ref,
        af_expr=vmt.v4_af,
        ab_expr=vmt.AD[1] / vmt.DP,
        use_v3_1_correction=True,
    )
    vmt = vmt.annotate_entries(GT=gt_with_depletion)
    vmt.entries().show()

    return vmt.checkpoint(new_temp_file("consent_samples_vmt", "mt"))


def _calculate_consent_frequencies(vmt: hl.MatrixTable, test: bool = False) -> hl.Table:
    """
    Calculate frequencies for consent withdrawal samples.

    :param vmt: Prepared MatrixTable with consent samples.
    :param test: Whether running in test mode.
    :return: Frequency table for consent samples.
    """
    logger.info("Loading group membership and calculating consent frequencies...")

    # Load and filter group membership
    group_membership_ht = get_group_membership(test=test).ht()
    consent_sample_ids = set(vmt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(consent_sample_ids).contains(group_membership_ht.s)
    )

    # Annotate MatrixTable with group membership
    vmt = vmt.annotate_cols(
        group_membership=group_membership_ht[vmt.col_key].group_membership,
    )
    vmt = vmt.annotate_rows(hists_fields=mt_hists_fields(vmt))

    # Load consent ANs
    # consent_ans_ht = get_consent_ans(test=test).ht()

    logger.info(
        "Calculating AC and hom alt counts per group membership using agg_by_strata..."
    )
    logger.info("Using updated agg_by_strata implementation - version check")
    # Use efficient agg_by_strata approach with localize_entries
    consent_freq_ht = agg_by_strata(
        vmt.select_entries(
            "GT",
            "adj",  # Required by agg_by_strata for quality filtering
            n_alt_alleles=vmt.GT.n_alt_alleles(),
            is_hom_var=vmt.GT.is_hom_var(),
        ),
        {
            "AC": (lambda t: t.n_alt_alleles, hl.agg.sum),
            "homozygote_count": (lambda t: t.is_hom_var, hl.agg.count_where),
        },
        group_membership_ht=group_membership_ht,
    )

    # Now combine with consent ANs to create proper frequency structs
    # consent_ans_ht = get_consent_ans(test=test).ht()
    consent_freq_ht = consent_freq_ht.annotate(
        freq=hl.range(hl.len(consent_freq_ht.AC)).map(
            lambda i: hl.struct(
                AC=hl.int32(consent_freq_ht.AC[i]),
                AF=hl.if_else(
                    consent_freq_ht.AC[i] > 0,
                    consent_freq_ht.AC[i]
                    / hl.float32(866 * 2),  # consent_ans_ht[consent_freq_ht.key].AN[i],
                    0.0,
                ),
                AN=hl.int(866 * 2),  # consent_ans_ht[consent_freq_ht.key].AN[i],
                homozygote_count=hl.int32(consent_freq_ht.homozygote_count[i]),
            )
        )
    )
    # Create freq_meta from group membership table (agg_by_strata doesn't create this)

    group_membership_globals = group_membership_ht.index_globals()
    consent_freq_ht = consent_freq_ht.annotate_globals(
        consent_freq_meta=group_membership_globals.freq_meta,
        consent_freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )

    return consent_freq_ht.checkpoint(new_temp_file("consent_freq", "ht"))


def _subtract_consent_frequencies_and_histograms(
    v4_freq_ht: hl.Table,
    consent_freq_ht: hl.Table,
    vmt: hl.MatrixTable,
    test: bool = False,
) -> hl.Table:
    """
    Subtract consent withdrawal frequencies and age histograms from v4 frequency table.

    :param v4_freq_ht: v4 frequency table (contains both freq and histograms.age_hists).
    :param consent_freq_ht: Consent withdrawal frequency table.
    :param vmt: MatrixTable with consent samples for age histogram calculation.
    :param test: Whether running in test mode.
    :return: Updated frequency table with consent frequencies and age histograms subtracted.
    """
    logger.info(
        "Subtracting consent withdrawal frequencies and age histograms from v4 frequency table..."
    )

    # Calculate consent age histograms
    logger.info("Calculating age histograms for consent withdrawal samples...")
    vmt_adj = vmt.filter_entries(vmt.adj)
    vmt_with_age_hists = vmt_adj.annotate_rows(
        age_hists=age_hists_expr(
            hl.call(vmt_adj.GT), vmt_adj.age, vmt_adj.sex_karyotype
        )
    )
    consent_age_hist_ht = vmt_with_age_hists.rows().select("age_hists")

    # Join all data together
    joined_freq_ht = v4_freq_ht.annotate(
        consent_freq=consent_freq_ht[v4_freq_ht.key].freq,
        consent_age_hists=consent_age_hist_ht[v4_freq_ht.key].age_hists,
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=consent_freq_ht.index_globals().consent_freq_meta,
        consent_freq_meta_sample_count=consent_freq_ht.index_globals().consent_freq_meta_sample_count,
    )

    # Subtract consent frequencies
    logger.info("Subtracting consent frequencies...")
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
        set_negatives_to_zero=True,
        # TODO: Remove this once we have the consent AN file, the AN place holder
        # is occassionally more than site ANs in the freq HT
    )

    # Subtract consent age histograms
    logger.info("Subtracting consent age histograms...")
    updated_age_hist_het = merge_histograms(
        [
            joined_freq_ht.histograms.age_hists.age_hist_het,
            joined_freq_ht.consent_age_hists.age_hist_het,
        ],
        operation="diff",
    )
    updated_age_hist_hom = merge_histograms(
        [
            joined_freq_ht.histograms.age_hists.age_hist_hom,
            joined_freq_ht.consent_age_hists.age_hist_hom,
        ],
        operation="diff",
    )

    # Update the frequency table with all changes
    joined_freq_ht = joined_freq_ht.annotate(
        freq=updated_freq_expr,
        histograms=joined_freq_ht.histograms.annotate(
            age_hists=joined_freq_ht.histograms.age_hists.annotate(
                age_hist_het=updated_age_hist_het,
                age_hist_hom=updated_age_hist_hom,
            )
        ),
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        freq_meta=hl.literal(hl.eval(updated_freq_meta)),
        freq_meta_sample_count=hl.literal(
            hl.eval(updated_sample_counts["freq_meta_sample_count"])
        ),
    )

    # Clean up temporary fields
    joined_freq_ht = joined_freq_ht.drop("consent_freq", "consent_age_hists")

    joined_freq_ht.select("freq", "histograms").show()
    return joined_freq_ht.checkpoint(new_temp_file("merged_freq_and_hists", "ht"))


def process_gnomad_dataset(
    data_test: bool = False,
    runtime_test: bool = False,
    overwrite: bool = False,
) -> hl.Table:
    """
    Process gnomAD dataset to update v4 frequency HT by removing consent withdrawal samples.

    This function performs frequency adjustment by:
    1. Loading v4 frequency HT (contains both frequencies and age histograms)
    2. Loading consent withdrawal VDS
    3. Filtering to sites present in BOTH consent VDS AND v4 frequency table
    4. Calculating frequencies and age histograms for consent withdrawal samples
    5. Subtracting both frequencies and age histograms from v4 frequency HT
    6. Only overwriting fields that were actually updated in the final output
    7. Computing FAF, grpmax, gen_anc_faf_max, and inbreeding coefficient

    :param data_test: Whether to run in data test mode. If True, filters v4 vds to first 2 partitions.
    :param runtime_test: Whether to run in runtime test mode. If True, filters v4 test vds to first 2 partitions.
    :param overwrite: Whether to overwrite existing output files.
    :return: Updated frequency HT with FAF/grpmax annotations and updated age histograms for gnomAD dataset.
    """
    test = runtime_test or data_test

    logger.info("Processing gnomAD dataset for consent withdrawals...")

    # Load v4 frequency table (contains both frequencies and age histograms)
    logger.info("Loading v4 frequency table...")
    v4_freq_ht = get_v4_freq(data_type="genomes").ht()

    # Prepare consent VDS (without filtering yet)
    logger.info("Loading consent VDS...")
    vmt = _prepare_consent_vds(v4_freq_ht, test=test, runtime_test=runtime_test)

    # Filter v4 frequency table to sites present in consent VDS
    logger.info("Filtering v4 frequency table to sites present in consent VDS...")
    consent_sites = vmt.rows().key_by("locus", "alleles").select().distinct()
    v4_freq_ht_filtered = v4_freq_ht.semi_join(consent_sites)

    # Calculate frequencies for consent samples (vmt already contains only
    # relevant sites)
    consent_freq_ht = _calculate_consent_frequencies(vmt, test=test)

    # Subtract consent frequencies and age histograms from v4 frequency table
    updated_freq_ht = _subtract_consent_frequencies_and_histograms(
        v4_freq_ht_filtered, consent_freq_ht, vmt, test=test
    )

    # Calculate FAF, grpmax, and other post-processing annotations
    updated_freq_ht = _calculate_faf_and_grpmax_annotations(updated_freq_ht, test=test)

    # Only overwrite fields that were actually updated (merge back with
    # original full table)
    logger.info("Preparing final frequency table with only updated fields...")
    final_freq_ht = _merge_updated_frequency_fields(v4_freq_ht, updated_freq_ht)

    return final_freq_ht


def _calculate_faf_and_grpmax_annotations(
    updated_freq_ht: hl.Table, test: bool = False
) -> hl.Table:
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


def _merge_updated_frequency_fields(
    original_freq_ht: hl.Table, updated_freq_ht: hl.Table
) -> hl.Table:
    """
    Merge frequency tables, only overwriting fields that were actually updated.

    For sites that exist in updated_freq_ht, use the updated values.
    For sites that don't exist in updated_freq_ht, keep original values.

    :param original_freq_ht: Original v4 frequency table.
    :param updated_freq_ht: Updated frequency table with consent withdrawals subtracted.
    :return: Final frequency table with selective field updates.
    """
    logger.info("Merging frequency tables with selective field updates...")

    # For sites in updated_freq_ht, use updated values; otherwise keep original
    final_freq_ht = original_freq_ht.annotate(
        **{
            field: hl.if_else(
                hl.is_defined(updated_freq_ht[original_freq_ht.key]),
                updated_freq_ht[original_freq_ht.key][field],
                original_freq_ht[field],
            )
            for field in [
                "freq",
                "faf",
                "grpmax",
                "gen_anc_faf_max",
                "inbreeding_coeff",
                "histograms",
            ]
            if field in updated_freq_ht.row
        }
    )

    # Update globals from updated table if they exist
    updated_globals = {}
    for global_field in ["freq_meta", "faf_meta", "freq_index_dict", "faf_index_dict"]:
        if global_field in updated_freq_ht.globals:
            updated_globals[global_field] = updated_freq_ht.index_globals()[
                global_field
            ]

    if updated_globals:
        final_freq_ht = final_freq_ht.annotate_globals(**updated_globals)

    return final_freq_ht


def _resolve_aou_sample_collisions(
    aou_vds: hl.vds.VariantDataset, test: bool = False
) -> hl.vds.VariantDataset:
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


def _filter_aou_related_samples(
    aou_vds: hl.vds.VariantDataset, test: bool = False
) -> hl.vds.VariantDataset:
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
    if len(aou_removed_sample_ids) == 0:
        logger.info("No related samples found to remove")
        return aou_vds
    return hl.vds.filter_samples(aou_vds, aou_removed_sample_ids, keep=False)


def _prepare_aou_variant_data(
    aou_vds_filtered: hl.vds.VariantDataset, test: bool = False
) -> hl.MatrixTable:
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
            fixed_homalt_model=hl.bool(
                False
            ),  # Always False for genomes (v3-style correction applied to all)
        )
    else:
        # Production mode expects direct metadata fields
        aou_variant_mt = aou_variant_mt.annotate_cols(
            sex_karyotype=aou_variant_mt.meta.sex_karyotype,
            pop=aou_variant_mt.meta.pop,
            age=aou_variant_mt.meta.age,
            fixed_homalt_model=hl.bool(
                False
            ),  # Always False for genomes (v3-style correction applied to all)
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


def _calculate_aou_variant_frequencies(
    aou_variant_mt: hl.MatrixTable, test: bool = False
) -> hl.Table:
    """
    Calculate complete frequency struct for AoU variant data using imported AN values.

    :param aou_variant_mt: Prepared variant MatrixTable
    :param test: Whether to use test resources
    :return: Frequency Table with complete freq struct including imported AN
    """
    # Use existing AoU group membership table and filter to variant samples
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht = get_group_membership(subset="aou", test=test).ht()

    # Filter group membership to only the AoU samples we're processing
    aou_sample_ids = set(aou_variant_mt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(aou_sample_ids).contains(group_membership_ht.s)
    )

    # Load AN values from consent_ans (calculated by another script)
    logger.info("Loading AN values from consent_ans...")
    consent_ans_ht = get_consent_ans(test=test).ht()

    logger.info("Calculating AoU variant frequencies using efficient agg_by_strata...")
    # Use efficient agg_by_strata approach for AoU variant data
    aou_variant_freq_ht = agg_by_strata(
        aou_variant_mt.select_entries(
            "GT",
            "adj",  # Required by agg_by_strata for quality filtering
            n_alt_alleles=aou_variant_mt.GT.n_alt_alleles(),
            is_hom_var=aou_variant_mt.GT.is_hom_var(),
        ),
        {
            "AC": (lambda t: t.n_alt_alleles, hl.agg.sum),
            "homozygote_count": (lambda t: t.is_hom_var, hl.agg.count_where),
        },
        group_membership_ht=group_membership_ht,
    )

    # Join with consent AN values and build complete freq struct
    logger.info("Building complete frequency struct with imported AN values...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        consent_an=consent_ans_ht[aou_variant_freq_ht.key].AN
    )

    # Build complete freq struct using AC from variant data and AN from consent_ans
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        freq=hl.map(
            lambda freq_data, an_data: hl.struct(
                AC=freq_data.AC,
                AF=hl.if_else(an_data > 0, freq_data.AC / an_data, 0.0),
                AN=an_data,
                homozygote_count=freq_data.homozygote_count,
            ),
            aou_variant_freq_ht.freq,
            aou_variant_freq_ht.consent_an,
        )
    ).drop("consent_an")

    # Create freq_meta from group membership table (agg_by_strata doesn't create this)
    group_membership_globals = group_membership_ht.index_globals()
    aou_variant_freq_ht = aou_variant_freq_ht.annotate_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )
    return aou_variant_freq_ht


def _calculate_aou_age_histograms(
    aou_vds_filtered: hl.vds.VariantDataset, test: bool = False
) -> hl.Table:
    """
    Calculate age histograms for AoU dataset.

    :param aou_vds_filtered: Filtered AoU VDS with resolved collisions and related samples removed
    :param test: Whether this is a test run (returns simplified histogram for compatibility)
    :return: Age histogram Table with age_hist_het and age_hist_hom fields
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
        # Production mode uses sparse VDS variant data directly
        logger.info("Calculating age histograms using sparse VDS variant data...")
        aou_variant_mt = aou_vds_filtered.variant_data

        # Ensure proper field access and add adj filtering
        aou_variant_mt = aou_variant_mt.annotate_entries(GT=aou_variant_mt.LGT)
        aou_variant_mt = annotate_adj(aou_variant_mt)

        # Filter to high-quality calls for age histogram calculation
        aou_variant_mt = aou_variant_mt.filter_entries(aou_variant_mt.adj)

    aou_age_hist_expr = age_hists_expr(
        aou_variant_mt.GT,
        aou_variant_mt.meta.age,
        aou_variant_mt.meta.sex_karyotype,
    )

    logger.info("Aggregating age histograms across all variants...")
    return aou_variant_mt.aggregate_rows(
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
    1. Computing complete frequency struct using imported AN from consent_ans
    2. Age histograms are calculated within the frequency calculation

    :param test: Whether to run in test mode.
    :param overwrite: Whether to overwrite existing results.
    :return: Tuple of (Frequency Table for AoU dataset, Age histogram Table for AoU dataset).
    """
    logger.info("Processing All of Us dataset...")

    # Load and validate AoU VDS
    # Note: No need to check AoU VDS existence as it's a function call that
    # will fail if missing
    aou_vds = get_aou_vds(remove_hard_filtered_samples=True, test=test)

    # Step 1: Resolve sample ID collisions
    aou_vds = _resolve_aou_sample_collisions(aou_vds, test=test)

    # Step 2: Filter related samples
    aou_vds_filtered = _filter_aou_related_samples(aou_vds, test=test)

    # Step 3: Calculate complete frequencies using imported AN values
    logger.info("Calculating AoU frequencies with imported AN values...")
    aou_variant_mt = _prepare_aou_variant_data(aou_vds_filtered, test=test)
    aou_freq_ht = _calculate_aou_variant_frequencies(aou_variant_mt, test=test)

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
        aou_freq_meta=aou_freq_ht.freq_meta,
        aou_freq_meta_sample_count=aou_freq_ht.freq_meta_sample_count,
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
        freq_meta=merged_meta,
        freq_meta_sample_count=sample_counts["counts"],
        freq_index_dict=make_freq_index_dict_from_meta(merged_meta),
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

    return merged_freq_ht, merged_age_hist_ht.select("age_hist_het", "age_hist_hom")


def main(args):
    """Generate v5 frequency data."""
    environment = args.environment
    data_test = args.data_test
    runtime_test = args.runtime_test
    test = data_test or runtime_test
    overwrite = args.overwrite

    hl.init(
        log=get_logging_path("v5_frequency_generation", environment=environment),
        tmp_dir=(
            f"gs://{WORKSPACE_BUCKET}/tmp/{TMP_DIR_DAYS}_day"
            if environment == "rwb"
            else "gs://gnomad-tmp-4day"
        ),
    )
    hl.default_reference("GRCh38")

    try:
        logger.info("Running generate_frequency.py...")
        if args.process_gnomad:
            logger.info("Processing gnomAD dataset...")

            gnomad_freq = get_freq(test=test, data_type="genomes", subset="gnomad")

            check_resource_existence(
                output_step_resources={"process-gnomad": [gnomad_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = process_gnomad_dataset(
                data_test=data_test, runtime_test=runtime_test, overwrite=overwrite
            )

            logger.info(
                f"Writing gnomAD frequency HT (with embedded age histograms) to {gnomad_freq.path}..."
            )
            gnomad_freq_ht.write(gnomad_freq.path, overwrite=overwrite)

        if args.process_aou:
            logger.info("Processing All of Us dataset...")

            aou_freq = get_freq(test=test, data_type="genomes", subset="aou")
            # Note: AoU still uses separate age histogram tables for now
            aou_age_hist = get_freq(
                test=test, data_type="genomes", subset="aou"
            )  # Placeholder - will be separate age hist table

            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq, aou_age_hist]},
                overwrite=overwrite,
            )

            aou_freq_ht, aou_age_hist_ht = process_aou_dataset(
                test=test, overwrite=overwrite
            )

            logger.info(f"Writing AoU frequency HT to {aou_freq.path}...")
            aou_freq_ht.write(aou_freq.path, overwrite=overwrite)

            logger.info(f"Writing AoU age histogram HT to {aou_age_hist.path}...")
            aou_age_hist_ht.write(aou_age_hist.path, overwrite=overwrite)

        if args.merge_datasets:
            logger.info(
                "Merging frequency data and age histograms from both datasets..."
            )

            merged_freq = get_freq(test=test, data_type="genomes")
            # Note: merged output will be a single frequency table with embedded age
            # histograms

            check_resource_existence(
                output_step_resources={"merge-datasets": [merged_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_input = get_freq(
                test=test, data_type="genomes", subset="gnomad"
            )
            aou_freq_input = get_freq(test=test, data_type="genomes", subset="aou")
            # Note: gnomAD freq table now contains embedded age histograms
            aou_age_hist_input = get_freq(
                test=test, data_type="genomes", subset="aou"
            )  # Placeholder for AoU age hist

            check_resource_existence(
                input_step_resources={
                    "process-gnomad": [gnomad_freq_input],
                    "process-aou": [aou_freq_input, aou_age_hist_input],
                }
            )

            gnomad_freq_ht = get_freq(test=test, data_type="genomes", subset="gnomad")
            aou_freq_ht = get_freq(test=test, data_type="genomes", subset="aou")
            # Extract age histograms from gnomAD frequency table for merging
            gnomad_age_hist_ht = gnomad_freq_ht.ht().select(
                age_hist_het=gnomad_freq_ht.ht().histograms.age_hists.age_hist_het,
                age_hist_hom=gnomad_freq_ht.ht().histograms.age_hists.age_hist_hom,
            )
            aou_age_hist_ht = get_freq(
                test=test, data_type="genomes", subset="aou"
            )  # Placeholder for AoU age hist

            merged_freq_ht, merged_age_hist_ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq_ht,
                aou_freq_ht,
                gnomad_age_hist_ht,
                aou_age_hist_ht,
                test=test,
            )

            logger.info(f"Writing merged frequency HT to {merged_freq.path}...")
            merged_freq_ht.write(merged_freq.path, overwrite=overwrite)

            # Note: Age histograms are now embedded in the merged frequency table
            logger.info("Age histograms are embedded in the merged frequency table.")

    finally:
        hl.copy_log(get_logging_path("v5_frequency", environment=environment))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--data-test",
        help=f"Filter to the first {TEST_PARTITIONS} partitions of full VDS for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--runtime-test",
        help="Load test dataset and filter to test partitions.",
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
    parser.add_argument(
        "--environment",
        help="Environment to run in.",
        choices=["rwb", "dataproc"],
        default="rwb",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()
    main(args)
