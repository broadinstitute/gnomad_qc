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
from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    agg_by_strata,
    annotate_adj,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    build_freq_stratification_list,
    compute_freq_by_strata,
    faf_expr,
    gen_anc_faf_max_expr,
    generate_freq_group_membership_array,
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
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v5.resources.annotations import (
    get_age_hist,
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


def generate_freq_ht(
    mt: hl.MatrixTable,
    ds_ht: hl.Table,
    meta_ht: hl.Table,
) -> hl.Table:
    """
    Generate frequency Table.

    Assumes all necessary annotations are present:
        `mt` annotations:
            - GT
            - adj
            - _high_ab_het_ref
            - fixed_homalt_model
            - fixed_homalt_model

        `ds_ht` annotations:
            - downsampling
            - downsamplings
            - ds_pop_counts

        `meta_ht` annotations:
            - pop
            - sex_karyotype
            - gatk_version
            - age
            - ukb_sample

    :param mt: Input MatrixTable.
    :param ds_ht: Table with downsampling annotations.
    :param meta_ht: Table with sample metadata annotations.
    :return: Hail Table with frequency annotations.
    """
    meta_ht = meta_ht.semi_join(mt.cols())
    additional_strata_expr = [
        {"gatk_version": meta_ht.gatk_version},
        {"gatk_version": meta_ht.gatk_version, "pop": meta_ht.pop},
    ]
    logger.info("Building frequency stratification list...")
    strata_expr = build_freq_stratification_list(
        sex_expr=meta_ht.sex_karyotype,
        gen_anc_expr=meta_ht.pop,
        additional_strata_expr=additional_strata_expr,
        downsampling_expr=ds_ht[meta_ht.key].downsampling,
    )
    group_membership_ht = generate_freq_group_membership_array(
        meta_ht,
        strata_expr,
        downsamplings=hl.eval(ds_ht.downsamplings),
        ds_gen_anc_counts=hl.eval(ds_ht.ds_pop_counts),
    )
    group_membership = group_membership_ht[mt.col_key].group_membership

    logger.info("Annotating frequencies and counting high AB het calls...")
    freq_ht = compute_freq_by_strata(
        mt.annotate_cols(group_membership=group_membership),
        entry_agg_funcs={"high_ab_hets_by_group": (high_ab_het, hl.agg.sum)},
        select_fields=["hists_fields"],
    )
    # Note: To use "multi_way_zip_join" need globals to be the same but an if_else based
    #  on strata doesn't work because hail looks for the annotation
    #  "non_ukb_downsamplings" regardless of the conditional value and throws an error
    #  if it doesn't exist.
    group_membership_globals = group_membership_ht.index_globals()
    freq_ht = freq_ht.annotate_globals(**group_membership_globals)

    return freq_ht


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

    logger.info("Computing sex adjusted genotypes and quality annotations...")
    ab_expr = vmt.AD[1] / vmt.DP
    ab_cutoff = 0.9
    gt_expr = adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype)

    vmt = vmt.select_entries(
        "AD",
        "DP",
        "GQ",
        "_het_non_ref",
        GT=gt_expr,
        adj=get_adj_expr(gt_expr, vmt.GQ, vmt.DP, vmt.AD),
        _het_ab=ab_expr,
        _high_ab_het_ref=(ab_expr > ab_cutoff) & ~vmt._het_non_ref,
    )

    logger.info("Applying v4 genomes hom alt depletion fix...")
    gt_with_depletion = hom_alt_depletion_fix(
        vmt.GT,
        het_non_ref_expr=vmt._het_non_ref,
        af_expr=vmt.v4_af,
        ab_expr=vmt.AD[1] / vmt.DP,
        use_v3_1_correction=True,
    )
    vmt = vmt.annotate_entries(GT=gt_with_depletion)

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
    logger.info(
        f"Consent frequency metadata: {hl.eval(consent_freq_ht.consent_freq_meta)}"
    )

    return consent_freq_ht.checkpoint(new_temp_file("consent_freq", "ht"))


def _subtract_consent_frequencies(
    v4_freq_ht: hl.Table, consent_freq_ht: hl.Table
) -> hl.Table:
    """
    Subtract consent withdrawal frequencies from v4 frequency table.

    :param v4_freq_ht: v4 frequency table.
    :param consent_freq_ht: Consent withdrawal frequency table.
    :return: Updated frequency table with consent frequencies subtracted.
    """
    logger.info("Subtracting consent withdrawal frequencies from v4 frequency table...")

    # Join consent frequencies onto v4 frequency table
    joined_freq_ht = v4_freq_ht.annotate(
        consent_freq=consent_freq_ht[v4_freq_ht.key].freq
    )

    # Add consent globals first (following v4 pattern)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=consent_freq_ht.index_globals().consent_freq_meta,
        consent_freq_meta_sample_count=consent_freq_ht.index_globals().consent_freq_meta_sample_count,
    )

    logger.info("Subtracting consent frequencies from v4 frequency table...")
    # Use merge_freq_arrays for proper frequency subtraction
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

    # Apply changes - use hl.literal() to break source references from merge_freq_arrays
    joined_freq_ht = joined_freq_ht.annotate(joined_freq=updated_freq_expr)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        joined_freq_meta=hl.literal(hl.eval(updated_freq_meta)),
        joined_freq_meta_sample_count=hl.literal(
            hl.eval(updated_sample_counts["freq_meta_sample_count"])
        ),
    )
    joined_freq_ht.show()
    return joined_freq_ht.checkpoint(new_temp_file("merged_freq", "ht"))


def _calculate_and_subtract_age_histograms(
    vmt: hl.MatrixTable, test: bool = False
) -> hl.Table:
    """
    Calculate and subtract age histograms for consent withdrawal samples.

    :param vmt: MatrixTable with consent samples.
    :param test: Whether running in test mode.
    :return: Updated age histogram table.
    """
    logger.info(
        "Calculating and subtracting age histograms for consent withdrawal samples..."
    )

    # Calculate age histograms using adjusted genotypes
    vmt_adj = vmt.filter_entries(vmt.adj)
    vmt_with_age_hists = vmt_adj.annotate_rows(
        age_hists=age_hists_expr(vmt_adj.GT, vmt_adj.age, vmt_adj.sex_karyotype)
    )
    consent_age_hist_ht = vmt_with_age_hists.rows().select("age_hists")

    # Load v4 age histograms and subtract consent histograms
    v4_age_hist_ht = get_age_hist(test=test, data_type="genomes").ht()

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

    return updated_age_hist_ht


def process_gnomad_dataset(
    data_test: bool = False,
    runtime_test: bool = False,
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

    :param data_test: Whether to run in data test mode. If True, filters v4 vds to first 2 partitions.
    :param runtime_test: Whether to run in runtime test mode. If True, filters v4 test vds to first 2 partitions.
    :param overwrite: Whether to overwrite existing output files.
    :return: Tuple of (updated frequency HT with FAF/grpmax annotations, updated age histogram HT) for gnomAD dataset.
    """
    test = runtime_test or data_test

    logger.info("Processing gnomAD dataset for consent withdrawals...")
    v4_freq_ht = get_v4_freq(data_type="genomes").ht()

    # Prepare VDS with consent samples
    vmt = _prepare_consent_vds(v4_freq_ht, test=test, runtime_test=runtime_test)

    # Calculate frequencies for consent samples
    consent_freq_ht = _calculate_consent_frequencies(vmt, test=test)

    # Subtract consent frequencies from v4 frequencies
    updated_freq_ht = _subtract_consent_frequencies(v4_freq_ht, consent_freq_ht)

    # Calculate and subtract age histograms
    updated_age_hist_ht = _calculate_and_subtract_age_histograms(vmt, test=test)

    # Calculate FAF, grpmax, and other post-processing annotations
    updated_freq_ht = _calculate_faf_and_grpmax_annotations(updated_freq_ht, test=test)

    return updated_freq_ht, updated_age_hist_ht


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
    Calculate frequencies from AoU variant data (gives correct AC, hom but incorrect AN).

    :param aou_variant_mt: Prepared variant MatrixTable
    :param test: Whether to use test resources
    :return: Frequency Table with AC and hom counts from variant data
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
    # Create freq_meta from group membership table (agg_by_strata doesn't create this)
    group_membership_globals = group_membership_ht.index_globals()
    aou_variant_freq_ht = aou_variant_freq_ht.annotate_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )
    logger.info(
        f"AoU variant frequency metadata: {hl.eval(aou_variant_freq_ht.freq_meta)}"
    )

    return aou_variant_freq_ht


def _prepare_aou_reference_data(
    aou_vds_filtered: hl.vds.VariantDataset, test: bool = False
) -> hl.MatrixTable:
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
            fixed_homalt_model=hl.bool(
                False
            ),  # Always False for genomes (v3-style correction applied to all)
        )
    else:
        # Production mode expects direct metadata fields
        aou_reference_mt = aou_reference_mt.annotate_cols(
            sex_karyotype=aou_reference_mt.meta.sex_karyotype,
            pop=aou_reference_mt.meta.pop,
            fixed_homalt_model=hl.bool(
                False
            ),  # Always False for genomes (v3-style correction applied to all)
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


def _calculate_aou_reference_an(
    aou_reference_mt: hl.MatrixTable, test: bool = False
) -> hl.Table:
    """
    Calculate correct AN values from AoU reference data (all sites).

    :param aou_reference_mt: Prepared reference MatrixTable
    :param test: Whether to use test resources
    :return: Frequency Table with correct AN values from all sites
    """
    # Use existing AoU group membership table and filter to reference samples
    logger.info("Loading AoU group membership table for reference AN calculation...")
    group_membership_ht = get_group_membership(subset="aou", test=test).ht()

    # Filter group membership to only the AoU samples we're processing
    aou_ref_sample_ids = set(aou_reference_mt.s.collect())
    group_membership_ht = group_membership_ht.filter(
        hl.literal(aou_ref_sample_ids).contains(group_membership_ht.s)
    )

    logger.info("Calculating AoU reference AN using efficient agg_by_strata...")
    # Use efficient agg_by_strata approach for AoU reference AN calculation
    aou_reference_an_ht = agg_by_strata(
        aou_reference_mt.select_entries(
            "GT",
            "adj",  # Required by agg_by_strata for quality filtering
            ploidy=aou_reference_mt.GT.ploidy,
        ),
        {
            "AN": (lambda t: t.ploidy, hl.agg.sum),
        },
        group_membership_ht=group_membership_ht,
    )
    # Create freq_meta from group membership table (agg_by_strata doesn't create this)
    group_membership_globals = group_membership_ht.index_globals()
    aou_reference_an_ht = aou_reference_an_ht.annotate_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )
    logger.info(f"AoU reference AN metadata: {hl.eval(aou_reference_an_ht.freq_meta)}")

    return aou_reference_an_ht


def _correct_aou_frequencies(
    aou_freq_ht: hl.Table, aou_ref_freq_ht: hl.Table
) -> hl.Table:
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
        # Production mode uses proper age_hists_expr with densified MT
        logger.info("Densifying AoU VDS for age histogram calculation...")
        aou_mt_dense = hl.vds.to_dense_mt(aou_vds_filtered)

        # Ensure proper field access and add adj filtering
        aou_mt_dense = aou_mt_dense.annotate_entries(GT=aou_mt_dense.LGT)
        aou_mt_dense = annotate_adj(aou_mt_dense)

        # Filter to high-quality calls for age histogram calculation
        aou_mt_dense = aou_mt_dense.filter_entries(aou_mt_dense.adj)

        aou_age_hist_expr = age_hists_expr(
            aou_mt_dense.GT,
            aou_mt_dense.meta.age,
            aou_mt_dense.meta.sex_karyotype,
        )

        logger.info("Aggregating age histograms across all variants...")
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
    # Note: No need to check AoU VDS existence as it's a function call that
    # will fail if missing
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
            gnomad_age_hist = get_age_hist(
                test=test, data_type="genomes", subset="gnomad"
            )

            check_resource_existence(
                output_step_resources={
                    "process-gnomad": [gnomad_freq, gnomad_age_hist]
                },
                overwrite=overwrite,
            )

            gnomad_freq_ht, gnomad_age_hist_ht = process_gnomad_dataset(
                data_test=data_test, runtime_test=runtime_test, overwrite=overwrite
            )

            logger.info(f"Writing gnomAD frequency HT to {gnomad_freq.path}...")
            gnomad_freq_ht.write(gnomad_freq.path, overwrite=overwrite)

            logger.info(f"Writing gnomAD age histogram HT to {gnomad_age_hist.path}...")
            gnomad_age_hist_ht.write(gnomad_age_hist.path, overwrite=overwrite)

        if args.process_aou:
            logger.info("Processing All of Us dataset...")

            aou_freq = get_freq(test=test, data_type="genomes", subset="aou")
            aou_age_hist = get_age_hist(test=test, data_type="genomes", subset="aou")

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
            merged_age_hist = get_age_hist(test=test, data_type="genomes")

            check_resource_existence(
                output_step_resources={
                    "merge-datasets": [merged_freq, merged_age_hist]
                },
                overwrite=overwrite,
            )

            gnomad_freq_input = get_freq(
                test=test, data_type="genomes", subset="gnomad"
            )
            aou_freq_input = get_freq(test=test, data_type="genomes", subset="aou")
            gnomad_age_hist_input = get_age_hist(
                test=test, data_type="genomes", subset="gnomad"
            )
            aou_age_hist_input = get_age_hist(
                test=test, data_type="genomes", subset="aou"
            )

            check_resource_existence(
                input_step_resources={
                    "process-gnomad": [gnomad_freq_input, gnomad_age_hist_input],
                    "process-aou": [aou_freq_input, aou_age_hist_input],
                }
            )

            gnomad_freq_ht = get_freq(test=test, data_type="genomes", subset="gnomad")
            aou_freq_ht = get_freq(test=test, data_type="genomes", subset="aou")
            gnomad_age_hist_ht = get_age_hist(
                test=test, data_type="genomes", subset="gnomad"
            )
            aou_age_hist_ht = get_age_hist(test=test, data_type="genomes", subset="aou")

            merged_freq_ht, merged_age_hist_ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq_ht,
                aou_freq_ht,
                gnomad_age_hist_ht,
                aou_age_hist_ht,
                test=test,
            )

            logger.info(f"Writing merged frequency HT to {merged_freq.path}...")
            merged_freq_ht.write(merged_freq.path, overwrite=overwrite)

            logger.info(f"Writing merged age histogram HT to {merged_age_hist.path}...")
            merged_age_hist_ht.write(merged_age_hist.path, overwrite=overwrite)

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
