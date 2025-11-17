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
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj as annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    get_consent_ans,
    get_freq,
    get_group_membership,
)
from gnomad_qc.v5.resources.basics import (
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import consent_samples_to_drop

# Constants
TEST_PARTITIONS = 2  # Number of partitions to use in test mode
TMP_DIR_DAYS = 4  # Number of days for temp directory retention

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def _prepare_consent_vds(
    v4_freq_ht: hl.Table, test: bool = False, runtime_test: bool = False
) -> hl.vds.VariantDataset:
    """
    Load and prepare VDS for consent withdrawal sample processing.

    :param v4_freq_ht: v4 frequency table for AF annotation.
    :param test: Whether running in test mode.
    :param runtime_test: Whether to use test VDS.
    :return: Prepared VDS with consent samples, split multiallelics, and annotations.
    """
    logger.info("Loading and preparing VDS for consent withdrawal samples...")

    vds = get_gnomad_v5_genomes_vds(
        test=runtime_test,
        release_only=True,
        consent_drop_only=True,
        annotate_meta=True,
        filter_partitions=list(range(TEST_PARTITIONS)) if test else None,
    )

    logger.info(
        "VDS has been filtered to %s consent withdrawal samples...",
        vds.variant_data.count_cols(),
    )

    # Prepare variant data with metadata
    vmt = vds.variant_data
    vmt.describe()
    if test:
        vmt = vmt.select_cols(
            gen_anc=vmt.meta.population_inference.pop,
            sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
            age=vmt.meta.project_meta.age,
        )
    else:
        vmt = vmt.select_cols(
            gen_anc=vmt.meta.gen_anc,
            sex_karyotype=vmt.meta.sex_karyotype,
            age=vmt.meta.age,
        )
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

    # This follows the v3/v4 genomes workflow for adj and sex adjusted genotypes which
    # were added before the hom alt depletion fix.
    # The correct order is to do the hom alt fix before adjusting sex ploidy and before
    # determining the adj annotation because haploid GTs have different adj filtering
    # criteria, but the option to adjust ploidy after adj is included for consistency
    # with v3.1, where we added the adj annotation before adjusting for sex ploidy.
    logger.info("Computing sex adjusted genotypes and quality annotations...")
    vmt = vmt.annotate_entries(
        adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD),
    )
    vmt = vmt.select_entries(
        "AD",
        "DP",
        "GQ",
        "_het_non_ref",
        "adj",
        GT=adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype),
    )

    # We set use_v3_1_correction to True to mimic the v4 genomes approach.
    logger.info("Applying v4 genomes hom alt depletion fix...")
    vmt = vmt.annotate_entries(
        GT=hom_alt_depletion_fix(
            vmt.GT,
            het_non_ref_expr=vmt._het_non_ref,
            af_expr=vmt.v4_af,
            ab_expr=vmt.AD[1] / vmt.DP,
            use_v3_1_correction=True,
        )
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    return vds.checkpoint(new_temp_file("consent_samples_vds_prepared", "vds"))


def _calculate_consent_frequencies_and_age_histograms(
    vds: hl.vds.VariantDataset, test: bool = False
) -> hl.Table:
    """
    Calculate frequencies and age histograms for consent withdrawal samples.

    :param vds: Prepared VDS with consent samples.
    :param test: Whether running in test mode.
    :return: Table with freq and age_hists annotations for consent samples.
    """
    logger.info("Densifying VDS for frequency calculations...")

    # Densify the VDS to MatrixTable
    # The metadata (pop, sex_karyotype, age) was already extracted to column
    # level in _prepare_consent_vds
    mt = hl.vds.to_dense_mt(vds)

    # Try to load group membership resource, fall back to generating from
    # scratch for testing
    try:
        logger.info("Loading group membership resource...")
        group_membership_ht = get_group_membership(subset="gnomad").ht()
        consent_sample_ids = set(mt.s.collect())
        group_membership_ht = group_membership_ht.filter(
            hl.literal(consent_sample_ids).contains(group_membership_ht.s)
        )
        logger.info("Successfully loaded group membership resource")
    except Exception:
        logger.info(
            "Group membership resource not found, generating from scratch for testing..."
        )

        # Create meta table from MatrixTable columns (only consent samples)
        # Note: mt.cols() returns a Table keyed by the column key (typically 's')
        # Following the v4 pattern from compute_coverage.py:get_genomes_group_membership_ht
        # We need to get the cols table first, then reference fields from it
        cols_ht = mt.cols()
        n_cols = cols_ht.count()
        logger.info(f"MatrixTable has {n_cols} columns (samples)")

        meta_ht = cols_ht.select(
            gen_anc=cols_ht.gen_anc,
            sex_karyotype=cols_ht.sex_karyotype,
        )
        n_meta_rows = meta_ht.count()
        logger.info(f"meta_ht has {n_meta_rows} rows after select")

        if n_meta_rows == 0:
            raise ValueError(
                "meta_ht is empty after select. Check that gen_anc and sex_karyotype "
                "fields exist on MatrixTable columns."
            )

        # Build frequency stratification list using expressions from meta_ht
        # Following the v4 pattern - build_freq_stratification_list takes expressions
        # bound to meta_ht and returns a list of dictionaries
        logger.info("Building frequency stratification list...")
        strata_expr = build_freq_stratification_list(
            sex_expr=meta_ht.sex_karyotype,
            gen_anc_expr=meta_ht.gen_anc,
        )

        # Generate group membership array
        # Both meta_ht and strata_expr expressions must reference the same table
        logger.info("Generating group membership array...")
        group_membership_ht = generate_freq_group_membership_array(
            meta_ht,
            strata_expr,
        )

        # Verify group_membership_ht has rows
        n_samples_in_group_membership = group_membership_ht.count()
        logger.info(
            f"Group membership table has {n_samples_in_group_membership} samples"
        )
        if n_samples_in_group_membership == 0:
            raise ValueError(
                "Generated group_membership_ht is empty. This likely means meta_ht "
                "has no samples or all samples were filtered out."
            )

    # Annotate MatrixTable with group membership
    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership,
    )

    # Verify group_membership was annotated correctly
    n_samples_with_group_membership = mt.aggregate_cols(
        hl.agg.count_where(hl.is_defined(mt.group_membership))
    )
    logger.info(
        f"Annotated {n_samples_with_group_membership} samples with group_membership"
    )
    if n_samples_with_group_membership == 0:
        raise ValueError(
            "No samples have group_membership annotation. Check that group_membership_ht "
            "keys match MatrixTable column keys."
        )

    logger.info(
        "Calculating frequencies and age histograms using compute_freq_by_strata..."
    )

    # Annotate hists_fields on MatrixTable rows before calling compute_freq_by_strata
    # This is required when using select_fields=["hists_fields"]
    # Following v4 pattern but simplified - we only need age_hists
    logger.info("Annotating hists_fields on MatrixTable rows...")
    mt = mt.annotate_rows(
        hists_fields=hl.struct(
            age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        )
    )

    # Use compute_freq_by_strata to calculate frequencies and age histograms
    # This follows the v4 approach and automatically calculates AC, AN, AF, homozygote_count
    # and includes age histograms when select_fields=["hists_fields"] is specified
    consent_freq_ht = compute_freq_by_strata(
        mt,
        select_fields=["hists_fields"],
    )

    # Extract age_hists from hists_fields struct
    consent_freq_ht = consent_freq_ht.annotate(
        age_hists=consent_freq_ht.hists_fields.age_hists,
    ).drop("hists_fields")

    # Annotate globals from group membership table
    group_membership_globals = group_membership_ht.index_globals()
    consent_freq_ht = consent_freq_ht.annotate_globals(
        consent_freq_meta=group_membership_globals.freq_meta,
        consent_freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )

    return consent_freq_ht.checkpoint(new_temp_file("consent_freq_and_hists", "ht"))


def _subtract_consent_frequencies_and_age_histograms(
    v4_freq_ht: hl.Table,
    consent_freq_ht: hl.Table,
    test: bool = False,
) -> hl.Table:
    """
    Subtract consent withdrawal frequencies and age histograms from v4 frequency table.

    :param v4_freq_ht: v4 frequency table (contains both freq and histograms.age_hists).
    :param consent_freq_ht: Consent withdrawal table with freq and age_hists annotations.
    :param test: Whether running in test mode.
    :return: Updated frequency table with consent frequencies and age histograms subtracted.
    """
    logger.info(
        "Subtracting consent withdrawal frequencies and age histograms from v4 frequency table..."
    )
    if test:
        v4_freq_ht = v4_freq_ht.filter(hl.is_defined(consent_freq_ht[v4_freq_ht.key]))

    joined_freq_ht = v4_freq_ht.annotate(
        consent_freq=consent_freq_ht[v4_freq_ht.key].freq,
        consent_age_hists=consent_freq_ht[v4_freq_ht.key].age_hists,
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=consent_freq_ht.index_globals().consent_freq_meta,
        consent_freq_meta_sample_count=consent_freq_ht.index_globals().consent_freq_meta_sample_count,
    )

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
    )

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

    return joined_freq_ht.checkpoint(new_temp_file("merged_freq_and_hists", "ht"))


def _calculate_faf_and_grpmax_annotations(
    updated_freq_ht: hl.Table,
) -> hl.Table:
    """
    Calculate FAF, grpmax, gen_anc_faf_max, and inbreeding coefficient annotations.

    This function handles the complex post-processing annotations that are added
    to frequency tables after the core frequency calculations are complete.

    :param updated_freq_ht: Frequency table after consent withdrawal subtraction
    :return: Updated frequency table with FAF/grpmax annotations
    """
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


def process_gnomad_dataset(
    data_test: bool = False,
    runtime_test: bool = False,
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

    :param data_test: Whether to run in data test mode. If True, filters full v4 vds to first 2 partitions.
    :param runtime_test: Whether to run in runtime test mode. If True, filters v4 test vds to first 2 partitions.
    :return: Updated frequency HT with updated frequencies and age histograms for gnomAD dataset.
    """
    # Test is used for formatting file paths
    test = runtime_test or data_test

    logger.info("Processing gnomAD dataset for consent withdrawals...")

    # Load v4 frequency table (contains both frequencies and age histograms)
    logger.info("Loading v4 frequency table...")
    v4_freq_ht = get_v4_freq(data_type="genomes").ht()

    # Prepare consent VDS (without filtering yet)
    logger.info("Loading consent VDS...")
    vds = _prepare_consent_vds(v4_freq_ht, test=test, runtime_test=runtime_test)

    # Calculate frequencies and age histograms for consent samples
    consent_freq_ht = _calculate_consent_frequencies_and_age_histograms(vds, test=test)

    logger.info(
        "Subtracting consent frequencies and age histograms from v4 frequency table..."
    )

    updated_freq_ht = _subtract_consent_frequencies_and_age_histograms(
        v4_freq_ht, consent_freq_ht, test=test
    )

    # Only overwrite fields that were actually updated (merge back with
    # original full table)
    # Note: FAF/grpmax annotations will be calculated on the final merged dataset
    final_freq_ht = _merge_updated_frequency_fields(v4_freq_ht, updated_freq_ht)

    return final_freq_ht


def _merge_updated_frequency_fields(
    original_freq_ht: hl.Table, updated_freq_ht: hl.Table
) -> hl.Table:
    """
    Merge frequency tables, only overwriting fields that were actually updated.

    For sites that exist in updated_freq_ht, use the updated values.
    For sites that don't exist in updated_freq_ht, keep original values.

    Note: FAF/grpmax annotations are not calculated during consent withdrawal processing
    and will be calculated later on the final merged dataset.

    :param original_freq_ht: Original v4 frequency table.
    :param updated_freq_ht: Updated frequency table with consent withdrawals subtracted.
    :return: Final frequency table with selective field updates.
    """
    logger.info("Merging frequency tables with selective field updates...")

    # For sites in updated_freq_ht, use updated values; otherwise keep original
    # Only updates fields that are present in updated_freq_ht (freq and histograms)
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


def mt_hists_fields(mt: hl.MatrixTable) -> hl.StructExpression:
    """
    Annotate quality metrics histograms and age histograms onto MatrixTable.

    :param mt: Input MatrixTable.
    :return: Struct with quality metrics histograms and age histograms.
    """
    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )

    return hl.struct(
        **qual_hist_expr(
            gt_expr=mt.GT,
            gq_expr=mt.GQ,
            dp_expr=mt.DP,
            adj_expr=mt.adj,
            ab_expr=mt._het_ab,
            split_adj_and_raw=True,
        ),
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
    )


def _prepare_aou_vds(aou_vds: hl.vds.VariantDataset) -> hl.MatrixTable:
    """
    Prepare AoU VDS for frequency calculations.

    :param aou_vds: AoU VariantDataset
    :return: Prepared AoU VariantDataset
    """
    aou_vmt = aou_vds.variant_data
    # Use existing AoU group membership table and filter to variant samples
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht = get_group_membership(subset="aou").ht()

    # Ploidy is already adjusted in the AoU VDS because of DRAGEN, do not need
    # to adjust it here
    aou_vmt = aou_vmt.annotate_cols(
        sex_karyotype=aou_vmt.meta.sex_karyotype,
        gen_anc=aou_vmt.meta.gen_anc,
        age=aou_vmt.meta.age,
        group_membership=group_membership_ht[aou_vmt.col_key].group_membership,
    )
    aou_vmt = aou_vmt.annotate_globals(
        freq_meta=group_membership_ht.index_globals().freq_meta,
        freq_meta_sample_count=group_membership_ht.index_globals().freq_meta_sample_count,
    )
    # Add adj annotation required by annotate_freq
    aou_vmt = annotate_adj_no_dp(aou_vmt)

    # Create freq_meta from group membership table (agg_by_strata doesn't create this)
    group_membership_globals = group_membership_ht.index_globals()
    aou_vmt = aou_vmt.annotate_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )

    # Rename LGT to GT and LAD to AD for compatibility with annotate_freq
    aou_vds = hl.vds.split_multi(hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt))

    return aou_vds


def _calculate_aou_frequencies_and_hists_using_all_sites_ans(
    aou_variant_mt: hl.MatrixTable, test: bool = False
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data.

    :param aou_variant_mt: Prepared variant MatrixTable
    :param test: Whether to use test resources
    :return: Table with freq and age_hists annotations
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    aou_variant_mt = aou_variant_mt.annotate_rows(
        hists_fields=mt_hists_fields(aou_variant_mt)
    )
    logger.info("Annotating frequencies and age histograms...")

    logger.info("Using all sites ANs for frequency calculations...")
    group_membership_ht = get_group_membership(subset="aou").ht()
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

    # Load AN values from consent_ans (calculated by another script but used
    # same group membership HT so same strata order)
    logger.info("Loading AN values from consent_ans...")
    consent_ans_ht = get_consent_ans(test=test).ht()
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        consent_an=consent_ans_ht[aou_variant_freq_ht.key].AN
    )

    logger.info("Building complete frequency struct with imported AN values...")
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
        ),
    ).drop("consent_an")

    return aou_variant_freq_ht


def _calculate_aou_frequencies_and_hists_using_densify(
    aou_vds: hl.vds.VariantDataset, test: bool = False
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using densify.

    :param aou_vds: Prepared AoU VariantDataset
    :param test: Whether to use test resources
    :return: Table with freq and age_hists annotations
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    aou_mt = hl.vds.to_dense_mt(aou_vds)
    aou_mt = aou_mt.annotate_rows(hists_fields=mt_hists_fields(aou_mt))
    aou_freq_ht = compute_freq_by_strata(aou_mt, select_fields=["hists_fields"])
    return aou_freq_ht


def process_aou_dataset(
    test: bool = False, use_all_sites_ans: bool = False
) -> hl.Table:
    """
    Process All of Us dataset for frequency calculations and age histograms.

    This function efficiently processes the AoU VDS by:
    1. Computing complete frequency struct using imported AN from consent_ans
    2. Age histograms are calculated within the frequency calculation

    :param test: Whether to run in test mode.
    :param use_all_sites_ans: Whether to use all sites ANs for frequency calculations.
    :return: Table with freq and age_hists annotations for AoU dataset.
    """
    logger.info("Processing All of Us dataset...")
    aou_vds = get_aou_vds(release_only=True, test=test)
    aou_vds = _prepare_aou_vds(aou_vds)

    # Calculate frequencies and age histograms together
    logger.info("Calculating AoU frequencies and age histograms...")
    if use_all_sites_ans:
        logger.info("Using all sites ANs for frequency calculations...")
        aou_freq_ht = _calculate_aou_frequencies_and_hists_using_all_sites_ans(
            aou_vds.variant_data, test=test
        )
    else:
        logger.info("Using densify for frequency calculations...")
        aou_freq_ht = _calculate_aou_frequencies_and_hists_using_densify(
            aou_vds, test=test
        )

    aou_freq_ht = aou_freq_ht.transmute(**aou_freq_ht.hists_fields)

    return aou_freq_ht


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
    use_all_sites_ans = args.use_all_sites_ans
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
                data_test=data_test, runtime_test=runtime_test
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

            aou_freq_ht = process_aou_dataset(
                test=test, use_all_sites_ans=use_all_sites_ans
            )

            logger.info(f"Writing AoU frequency HT to {aou_freq.path}...")
            aou_freq_ht.write(aou_freq.path, overwrite=overwrite)

            # Write age histograms separately if needed
            logger.info(f"Writing AoU age histogram HT to {aou_age_hist.path}...")
            aou_freq_ht.select("age_hists").write(
                aou_age_hist.path, overwrite=overwrite
            )

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
            # Extract age histograms from frequency tables for merging
            gnomad_age_hist_ht = gnomad_freq_ht.ht().select(
                age_hist_het=gnomad_freq_ht.ht().histograms.age_hists.age_hist_het,
                age_hist_hom=gnomad_freq_ht.ht().histograms.age_hists.age_hist_hom,
            )
            # Extract age histograms from AoU frequency table
            aou_age_hist_ht = aou_freq_ht.ht().select(
                age_hist_het=aou_freq_ht.ht().age_hists.age_hist_het,
                age_hist_hom=aou_freq_ht.ht().age_hists.age_hist_hom,
            )

            merged_freq_ht, merged_age_hist_ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq_ht,
                aou_freq_ht,
                gnomad_age_hist_ht,
                aou_age_hist_ht,
                test=test,
            )

            # Calculate FAF, grpmax, and other post-processing annotations on merged
            # dataset
            logger.info(
                "Calculating FAF, grpmax, and other annotations on merged dataset..."
            )
            merged_freq_ht = _calculate_faf_and_grpmax_annotations(merged_freq_ht)

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
        "--use-all-sites-ans",
        help="Use all sites ANs in frequency calculations to avoid a densify.",
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
