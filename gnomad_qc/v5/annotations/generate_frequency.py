"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies and histograms for:
1. gnomAD dataset - updating v4 frequencies by subtracting consent withdrawal samples
2. AoU dataset - using either pre-computed allele numbers or a densify approach

Processing Workflow:
--------------------
gnomAD (--process-gnomad):
1. Load v4 frequency table (contains frequencies and age histograms)
2. Prepare consent withdrawal VDS (split multiallelics, annotate metadata)
3. Calculate frequencies and age histograms for consent samples
4. Subtract from v4 frequencies to get updated gnomAD v5 frequencies

AoU (--process-aou):
1. Load AoU VDS with metadata
2. Prepare VDS (annotate group membership, adjust for ploidy, split multi-allelics)
3. Calculate frequencies using either: All sites ANs (efficient, requires pre-computed
AN values) or Densify approach (standard, more resource intensive)
4. Generate age histograms during frequency calculation

Usage Examples:
---------------
# Process AoU dataset using all-sites ANs.
python generate_frequency.py --process-aou --use-all-sites-ans --environment rwb

# Process AoU on batch/QoB with custom resources.
python generate_frequency.py --process-aou --environment batch --app-name "aou_freq" --driver-cores 8 --worker-memory highmem

# Process gnomAD consent withdrawals
python generate_frequency.py --process-gnomad --environment dataproc

# Run gnomAD in test mode
python generate_frequency.py --process-gnomad --test --test-partitions 2
"""

import argparse
import logging

import hail as hl
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    agg_by_strata,
    compute_freq_by_strata,
    get_adj_expr,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_freq,
    group_membership,
)
from gnomad_qc.v5.resources.basics import (
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.meta import meta

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def mt_hist_fields(mt: hl.MatrixTable) -> hl.StructExpression:
    """
    Annotate allele balance quality metrics histograms and age histograms onto MatrixTable.

    :param mt: Input MatrixTable.
    :return: Struct with allele balance, quality metrics histograms, and age histograms.
    """
    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    qual_hists = qual_hist_expr(
        gt_expr=mt.GT,
        gq_expr=mt.GQ,
        dp_expr=hl.sum(mt.AD),
        adj_expr=mt.adj,
        ab_expr=(mt.AD[1] / hl.sum(mt.AD)),
        split_adj_and_raw=True,
    )
    return hl.struct(
        qual_hists=qual_hists,
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
    )


def _prepare_aou_vds(
    aou_vds: hl.vds.VariantDataset,
    use_all_sites_ans: bool = False,
    test: bool = False,
    environment: str = "rwb",
) -> hl.vds.VariantDataset:
    """
    Prepare AoU VDS for frequency calculations.

    :param aou_vds: AoU VariantDataset.
    :param use_all_sites_ans: Whether to use all sites ANs for frequency calculations.
    :param test: Whether running in test mode.
    :param environment: Environment being used. Default is "rwb". Must be one of "rwb", "batch", or "dataproc".
    :return: Prepared AoU VariantDataset.
    """
    aou_vmt = aou_vds.variant_data
    # Use existing AoU group membership table and filter to variant samples.
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht = group_membership(
        test=test, data_set="aou", environment=environment
    ).ht()

    logger.info("Selecting cols for frequency stratification...")
    aou_vmt = aou_vmt.select_cols(
        sex_karyotype=aou_vmt.meta.sex_karyotype,
        gen_anc=aou_vmt.meta.genetic_ancestry_inference.gen_anc,
        age=aou_vmt.meta.project_meta.age,
        group_membership=group_membership_ht[aou_vmt.col_key].group_membership,
    )
    # NOTE: At the time of writing this code, it was not yet decided if we would use
    # the all sites AN for the frequency calculcation, a cost-savings approach to
    # avoid a densify, or if we would densify within the frequency script to get
    # call stats. The split timing depends on whether we use all-sites ANs:
    # - With all-sites ANs: adjust ploidy -> adj -> split (preserves LGT for adj)
    # - Without: split -> adjust ploidy -> adj (consistent with v4 workflow)
    # Please see project notes for more details:
    # https://docs.google.com/document/d/109NuIQlZ3a8uw__iBMK7FpoO2bFFGECzEDnBpB1xW8E/edit?tab=t.0#heading=h.neezqe6t09ul
    if use_all_sites_ans:
        aou_vmt = aou_vmt.select_entries(
            LGT=adjusted_sex_ploidy_expr(
                aou_vmt.locus, aou_vmt.LGT, aou_vmt.sex_karyotype
            ),
            GQ=aou_vmt.GQ,
            LAD=aou_vmt.LAD,
            LA=aou_vmt.LA,
        )
        aou_vmt = annotate_adj_no_dp(aou_vmt)
        aou_vds = hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt)
        aou_vds = hl.vds.split_multi(aou_vds, filter_changed_loci=True)
        aou_vmt = aou_vds.variant_data
    else:
        aou_vds = hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt)
        aou_vds = hl.vds.split_multi(aou_vds, filter_changed_loci=True)
        aou_vmt = aou_vds.variant_data
        aou_vmt = aou_vmt.select_entries(
            GT=adjusted_sex_ploidy_expr(
                aou_vmt.locus, aou_vmt.GT, aou_vmt.sex_karyotype
            ),
            GQ=aou_vmt.GQ,
            AD=aou_vmt.AD,
        )
        aou_vmt = annotate_adj_no_dp(aou_vmt)

    logger.info("Annotating globals...")
    group_membership_globals = group_membership_ht.index_globals()
    aou_vmt = aou_vmt.select_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
        age_distribution=aou_vmt.aggregate_cols(hl.agg.hist(aou_vmt.age, 30, 80, 10)),
        downsamplings=group_membership_globals.downsamplings,
    )

    aou_vds = hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt)

    return aou_vds


def _calculate_aou_frequencies_and_hists_using_all_sites_ans(
    aou_variant_mt: hl.MatrixTable, test: bool = False, environment: str = "rwb"
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using all sites ANs.

    :param aou_variant_mt: Prepared variant MatrixTable.
    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb", "batch", or "dataproc".
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    all_sites_an_ht = coverage_and_an_path(test=test, environment=environment).ht()
    aou_variant_mt = aou_variant_mt.annotate_rows(
        hist_fields=mt_hist_fields(aou_variant_mt)
    )

    logger.info("Annotating frequencies with all sites ANs...")
    group_membership_ht = group_membership(
        test=test, data_set="aou", environment=environment
    ).ht()
    aou_variant_freq_ht = agg_by_strata(
        aou_variant_mt.select_entries(
            "GT",
            "adj",
            n_alt_alleles=aou_variant_mt.GT.n_alt_alleles(),
            is_hom_var=aou_variant_mt.GT.is_hom_var(),
        ),
        {
            "AC": (lambda t: t.n_alt_alleles, hl.agg.sum),
            "homozygote_count": (lambda t: t.is_hom_var, hl.agg.count_where),
        },
        group_membership_ht=group_membership_ht,
        select_fields=["hist_fields"],
    )

    # Load AN values from all sites ANs table (calculated by another script but used
    # same group membership HT so same strata order).
    logger.info("Annotating AN values from all sites ANs...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        all_sites_an=all_sites_an_ht[aou_variant_freq_ht.locus].AN
    )

    logger.info("Building complete frequency struct with imported AN values...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        freq=hl.map(
            lambda AC, hom_alt, AN: hl.struct(
                AC=hl.int32(AC),
                AF=hl.if_else(AN > 0, AC / AN, hl.missing(hl.tfloat64)),
                AN=hl.int32(AN),
                homozygote_count=hl.int32(hom_alt),
            ),
            aou_variant_freq_ht.AC,
            aou_variant_freq_ht.homozygote_count,
            aou_variant_freq_ht.all_sites_an,
        ),
    ).drop("all_sites_an")

    # Nest histograms to match gnomAD structure.
    # Note: hists_fields.qual_hists already contains raw_qual_hists and qual_hists
    # as nested fields due to split_adj_and_raw=True in qual_hist_expr.
    aou_variant_freq_ht = aou_variant_freq_ht.select(
        freq=aou_variant_freq_ht.freq,
        histograms=hl.struct(
            qual_hists=aou_variant_freq_ht.hist_fields.qual_hists.qual_hists,
            raw_qual_hists=aou_variant_freq_ht.hist_fields.qual_hists.raw_qual_hists,
            age_hists=aou_variant_freq_ht.hist_fields.age_hists,
        ),
    )

    return aou_variant_freq_ht


def _calculate_aou_frequencies_and_hists_using_densify(
    aou_vds: hl.vds.VariantDataset,
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using densify.

    :param aou_vds: Prepared AoU VariantDataset.
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    aou_mt = hl.vds.to_dense_mt(aou_vds)
    aou_mt = aou_mt.annotate_rows(hist_fields=mt_hist_fields(aou_mt))
    # hl.agg.call_stats is used within compute_freq_by_strata which returns int32s for
    # AC, AN, and homozygote_count so do not need to convert to int32 here.
    aou_freq_ht = compute_freq_by_strata(aou_mt, select_fields=["hist_fields"])

    # Nest histograms to match gnomAD structure.
    # Note: hist_fields.qual_hists already contains raw_qual_hists and qual_hists
    # as nested fields due to split_adj_and_raw=True in qual_hist_expr.
    aou_freq_ht = aou_freq_ht.transmute(
        histograms=hl.struct(
            qual_hists=aou_freq_ht.hist_fields.qual_hists.qual_hists,
            raw_qual_hists=aou_freq_ht.hist_fields.qual_hists.raw_qual_hists,
            age_hists=aou_freq_ht.hist_fields.age_hists,
        )
    )

    return aou_freq_ht


def process_aou_dataset(
    test: bool = False, use_all_sites_ans: bool = False, environment: str = "rwb"
) -> hl.Table:
    """
    Process All of Us dataset for frequency calculations and age histograms.

    This function efficiently processes the AoU VDS by:
    1. Computing complete frequency struct (uses imported AN from AoU all site ANs if requested)
    2. Generating age histograms within the frequency calculation

    :param test: Whether to run in test mode.
    :param use_all_sites_ans: Whether to use all sites ANs for frequency calculations.
    :param environment: Environment to use. Default is "rwb". Must be one of "rwb", "batch", or "dataproc".
    :return: Table with freq and age_hists annotations for AoU dataset.
    """
    aou_vds = get_aou_vds(annotate_meta=True, release_only=True, test=test)
    aou_vds = _prepare_aou_vds(
        aou_vds, use_all_sites_ans=use_all_sites_ans, test=test, environment=environment
    )

    # Calculate frequencies and age histograms together
    logger.info("Calculating AoU frequencies and age histograms...")
    if use_all_sites_ans:
        logger.info("Using all sites ANs for frequency calculations...")
        aou_freq_ht = _calculate_aou_frequencies_and_hists_using_all_sites_ans(
            aou_vds.variant_data, test=test, environment=environment
        )
    else:
        logger.info("Using densify for frequency calculations...")
        aou_freq_ht = _calculate_aou_frequencies_and_hists_using_densify(aou_vds)
    aou_freq_ht = select_final_dataset_fields(aou_freq_ht, dataset="aou")

    return aou_freq_ht


def _prepare_consent_vds(
    v4_ht: hl.Table,
    test: bool = False,
    test_partitions: int = 2,
) -> hl.vds.VariantDataset:
    """
    Load and prepare VDS for consent withdrawal sample processing.

    :param v4_ht: v4 release table for AF annotation.
    :param test: Whether running in test mode.
    :param test_partitions: Number of partitions to use in test mode. Default is 2.
    :return: Prepared VDS with consent samples, split multiallelics, and annotations.
    """
    logger.info("Loading and preparing VDS for consent withdrawal samples...")

    vds = get_gnomad_v5_genomes_vds(
        release_only=True,
        consent_drop_only=True,
        annotate_meta=True,
        filter_partitions=list(range(test_partitions)) if test else None,
    )

    logger.info(
        "VDS has been filtered to %s consent withdrawal samples...",
        vds.variant_data.count_cols(),
    )

    vmt = vds.variant_data
    vmt = vmt.select_cols(
        gen_anc=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        age=vmt.meta.project_meta.age,
    )

    logger.info("Selecting entries and annotating non_ref hets pre-split...")
    vmt = vmt.select_entries(
        "LA", "LAD", "DP", "GQ", "LGT", _het_non_ref=vmt.LGT.is_het_non_ref()
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = vds.checkpoint(new_temp_file("consent_samples_vds", "vds"))

    logger.info("Splitting multiallelics in gnomAD sample withdrawal VDS...")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    # Annotate with v4 frequencies for hom alt depletion fix
    vmt = vds.variant_data
    vmt = vmt.annotate_rows(v4_af=v4_ht[vmt.row_key].freq[0].AF)

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
    logger.info("Annotating age distribution...")
    vmt = vmt.annotate_globals(
        age_distribution=vmt.aggregate_cols(hl.agg.hist(vmt.age, 30, 80, 10))
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    return vds.checkpoint(new_temp_file("consent_samples_vds_prepared", "vds"))


def _calculate_consent_frequencies_and_age_histograms(
    vds: hl.vds.VariantDataset,
    test: bool = False,
) -> hl.Table:
    """
    Calculate frequencies and age histograms for consent withdrawal samples.

    :param vds: Prepared VDS with consent samples.
    :param test: Whether running in test mode.
    :return: Table with freq and age_hists annotations for consent samples.
    """
    logger.info("Densifying VDS for frequency calculations...")
    mt = hl.vds.to_dense_mt(vds)
    # Group membership table is already filtered to consent drop samples and is in GCS.
    group_membership_ht = group_membership(
        test=test, data_set="gnomad", environment="dataproc"
    ).ht()

    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership,
    )
    mt = mt.annotate_globals(
        freq_meta=group_membership_ht.index_globals().freq_meta,
        freq_meta_sample_count=group_membership_ht.index_globals().freq_meta_sample_count,
    )

    logger.info(
        "Calculating frequencies and age histograms using compute_freq_by_strata..."
    )

    # Annotate hists_fields on MatrixTable rows before calling compute_freq_by_strata so
    # to keep the age_hists annotation on the frequency table.
    logger.info("Annotating hists_fields on MatrixTable rows...")
    mt = mt.annotate_rows(
        hists_fields=hl.struct(
            age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        )
    )

    logger.info(
        "Computing frequencies for consent samples using compute_freq_by_strata..."
    )
    consent_freq_ht = compute_freq_by_strata(
        mt,
        select_fields=["hists_fields"],
    )

    consent_freq_ht = consent_freq_ht.transmute(
        age_hists=consent_freq_ht.hists_fields.age_hists,
    )

    return consent_freq_ht.checkpoint(new_temp_file("consent_freq_and_hists", "ht"))


def _subtract_consent_frequencies_and_age_histograms(
    v4_ht: hl.Table,
    consent_freq_ht: hl.Table,
) -> hl.Table:
    """
    Subtract consent withdrawal frequencies and age histograms from v4 frequency table.

    :param v4_ht: v4 release table (contains both freq and histograms.age_hists).
    :param consent_freq_ht: Consent withdrawal table with freq and age_hists annotations.
    :return: Updated frequency table with consent frequencies and age histograms subtracted.
    """
    logger.info(
        "Subtracting consent withdrawal frequencies and age histograms from v4 release table..."
    )

    joined_freq_ht = v4_ht.annotate(
        consent_freq=consent_freq_ht[v4_ht.key].freq,
        consent_age_hists=consent_freq_ht[v4_ht.key].age_hists,
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=consent_freq_ht.index_globals().freq_meta,
        consent_freq_meta_sample_count=consent_freq_ht.index_globals().freq_meta_sample_count,
    )

    logger.info("Subtracting consent frequencies...")
    updated_freq_expr, updated_freq_meta, updated_sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.consent_freq],
        [
            joined_freq_ht.index_globals().freq_meta,
            joined_freq_ht.index_globals().consent_freq_meta,
        ],
        operation="diff",
        count_arrays={
            "freq_meta_sample_count": [
                joined_freq_ht.index_globals().freq_meta_sample_count,
                joined_freq_ht.index_globals().consent_freq_meta_sample_count,
            ],
        },
    )
    # Update the frequency table with freq changes.
    joined_freq_ht = joined_freq_ht.annotate(freq=updated_freq_expr)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        freq_meta=updated_freq_meta,
        freq_meta_sample_count=updated_sample_counts["freq_meta_sample_count"],
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

    # Update the frequency table with age hist changes.
    joined_freq_ht = joined_freq_ht.annotate(
        histograms=joined_freq_ht.histograms.annotate(
            age_hists=joined_freq_ht.histograms.age_hists.annotate(
                age_hist_het=updated_age_hist_het,
                age_hist_hom=updated_age_hist_hom,
            )
        ),
    )

    return joined_freq_ht.checkpoint(new_temp_file("merged_freq_and_hists", "ht"))


def select_final_dataset_fields(ht: hl.Table, dataset: str = "gnomad") -> hl.Table:
    """
    Create final freq Table with only desired annotations.

    :param ht: Hail Table containing all annotations.
    :param dataset: Dataset to select final fields, either "gnomad" or "aou".
    :return: Hail Table with final annotations.
    """
    if dataset not in ["gnomad", "aou"]:
        raise ValueError(f"Invalid dataset: {dataset}")

    final_globals = ["freq_meta", "freq_meta_sample_count", "age_distribution"]
    final_fields = ["freq", "histograms"]

    if dataset == "aou":
        # AoU has one extra 'downsamplings' global field that is not present in gnomAD.
        final_globals.append("downsamplings")

    # Convert all int64 annotations in the freq struct to int32s for merging type
    # compatibility.
    ht = ht.annotate(
        freq=ht.freq.map(
            lambda x: x.annotate(
                **{k: hl.int32(v) for k, v in x.items() if v.dtype == hl.tint64}
            )
        )
    )

    return ht.select(*final_fields).select_globals(*final_globals)


def _fix_v4_global_age_distribution(freq_ht: hl.Table) -> hl.Table:
    """
    Fix the age distribution global annotation in the frequency table.

    :param freq_ht: Frequency table to annotate with the age distribution.
    :return: Frequency table with the age distribution global annotation fixed.
    """
    # Use v5 meta as age is already fixed in the v5 project metadata as are the consent
    # withdrawal samples' releasable field.
    meta_ht = meta().ht()
    meta_ht = meta_ht.filter(
        (meta_ht.release) & (meta_ht.project_meta.project == "gnomad")
    )
    meta_ht = meta_ht.annotate_globals(
        age_distribution=meta_ht.aggregate(
            hl.agg.hist(meta_ht.project_meta.age, 30, 80, 10)
        )
    )
    freq_ht = freq_ht.annotate_globals(
        age_distribution=meta_ht.index_globals().age_distribution
    )

    return freq_ht


def process_gnomad_dataset(
    test: bool = False,
    test_partitions: int = 2,
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

    :param test: Whether to run in test mode. If True, filters full v4 vds to first N partitions (N controlled by test_partitions).
    :param test_partitions: Number of partitions to use in test mode. Default is 2.
    :return: Updated frequency HT with updated frequencies and age histograms for gnomAD dataset.
    """
    v4_ht = release_sites(data_type="genomes").ht()

    vds = _prepare_consent_vds(
        v4_ht,
        test=test,
        test_partitions=test_partitions,
    )

    logger.info("Calculating frequencies and age histograms for consent samples...")
    consent_freq_ht = _calculate_consent_frequencies_and_age_histograms(vds, test)

    if test:
        v4_ht = v4_ht.filter(hl.is_defined(consent_freq_ht[v4_ht.key]))
        v4_ht = v4_ht.naive_coalesce(100).checkpoint(
            new_temp_file("v4_ht_filtered_test", "ht")
        )

    logger.info("Subtracting consent frequencies and age histograms from v4...")
    updated_freq_ht = _subtract_consent_frequencies_and_age_histograms(
        v4_ht, consent_freq_ht
    )

    logger.info("Merging updated frequency fields...")
    freq_ht = _merge_updated_frequency_fields(v4_ht, updated_freq_ht)

    logger.info("Reannotating gnomAD's age distribution global annotation...")
    freq_ht = _fix_v4_global_age_distribution(freq_ht)

    # Select only the fields that were updated as FAF/grpmax/inbreeding_coeff annotations
    # will be calculated on the final merged dataset.
    logger.info("Selecting final dataset fields...")
    final_freq_ht = select_final_dataset_fields(freq_ht, dataset="gnomad")

    return final_freq_ht


def _merge_updated_frequency_fields(
    v4_release_ht: hl.Table, updated_freq_ht: hl.Table
) -> hl.Table:
    """
    Merge frequency tables, only overwriting fields that were actually updated.

    For sites that exist in updated_freq_ht, use the updated values.
    For sites that don't exist in updated_freq_ht, keep original values.

    Note: FAF/grpmax/inbreeding_coeff annotations are not calculated during consent
    withdrawal processing and will be calculated later on the final merged dataset.

    :param v4_release_ht: Original v4 release table.
    :param updated_freq_ht: Updated frequency table with consent withdrawals subtracted.
    :return: Final frequency table with selective field updates.
    """
    logger.info("Merging frequency tables with selective field updates...")

    # Bring in updated values with a single lookup.
    updated_row = updated_freq_ht[v4_release_ht.key]

    # Update freq and age_hists in a single annotate to avoid source mismatch:
    # - freq: use updated if present, otherwise keep original
    # - histograms.age_hists: update only age_hists, preserving qual_hists and raw_qual_hists
    final_freq_ht = v4_release_ht.annotate(
        freq=hl.coalesce(updated_row.freq, v4_release_ht.freq),
        histograms=v4_release_ht.histograms.annotate(
            age_hists=hl.coalesce(
                updated_row.histograms.age_hists,
                v4_release_ht.histograms.age_hists,
            )
        ),
    )

    # Update globals from updated table.
    updated_globals = {}
    for global_field in ["freq_meta", "freq_meta_sample_count"]:
        if global_field in updated_freq_ht.globals:
            updated_globals[global_field] = updated_freq_ht.index_globals()[
                global_field
            ]

    if updated_globals:
        final_freq_ht = final_freq_ht.annotate_globals(**updated_globals)

    return final_freq_ht


def _initialize_hail(args) -> None:
    """
    Initialize Hail with appropriate configuration for the environment.

    :param args: Parsed command-line arguments.
    """
    environment = args.environment
    tmp_dir_days = args.tmp_dir_days

    if environment == "batch":
        batch_kwargs = {
            "backend": "batch",
            "log": get_logging_path("v5_frequency_generation", environment="batch"),
            "tmp_dir": (
                f"{qc_temp_prefix(environment='dataproc', days=tmp_dir_days)}frequency_generation"
            ),
            "gcs_requester_pays_configuration": args.gcp_billing_project,
            "regions": ["us-central1"],
        }
        # Add optional batch configuration parameters.
        for param in [
            "app_name",
            "driver_cores",
            "driver_memory",
            "worker_cores",
            "worker_memory",
        ]:
            value = getattr(args, param, None)
            if value is not None:
                batch_kwargs[param] = value

        hl.init(**batch_kwargs)
    else:
        hl.init(
            log=get_logging_path(
                "v5_frequency_generation",
                environment=environment,
                tmp_dir_days=tmp_dir_days,
            ),
            tmp_dir=f"{qc_temp_prefix(environment=environment, days=tmp_dir_days)}frequency_generation",
        )
    hl.default_reference("GRCh38")


def main(args):
    """Generate v5 frequency data."""
    environment = args.environment
    use_all_sites_ans = args.use_all_sites_ans
    test = args.test
    test_partitions = args.test_partitions
    overwrite = args.overwrite
    tmp_dir_days = args.tmp_dir_days

    _initialize_hail(args)

    try:
        logger.info("Running generate_frequency.py...")

        if args.process_gnomad:
            logger.info("Processing gnomAD dataset...")

            gnomad_freq = get_freq(
                test=test,
                data_type="genomes",
                data_set="gnomad",
                environment=environment,
            )

            check_resource_existence(
                output_step_resources={"process-gnomad": [gnomad_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = process_gnomad_dataset(
                test=test,
                test_partitions=test_partitions,
            )

            logger.info(
                "Writing gnomAD frequency HT (with embedded age histograms) to %s...",
                gnomad_freq.path,
            )
            gnomad_freq_ht.write(gnomad_freq.path, overwrite=overwrite)

        if args.process_aou:
            logger.info("Processing All of Us dataset...")
            aou_freq = get_freq(
                test=test,
                data_type="genomes",
                data_set="aou",
                environment=environment,
            )

            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq]},
                overwrite=overwrite,
            )

            aou_freq_ht = process_aou_dataset(
                test=test,
                use_all_sites_ans=use_all_sites_ans,
                environment=environment,
            )

            logger.info("Writing AoU frequency HT to %s...", aou_freq.path)
            aou_freq_ht.write(aou_freq.path, overwrite=overwrite)

    finally:
        hl.copy_log(
            get_logging_path(
                "v5_frequency_run", environment=environment, tmp_dir_days=tmp_dir_days
            )
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate frequency data for gnomAD v5."
    )

    # General arguments.
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )

    # Test/debug arguments.
    test_group = parser.add_argument_group("testing options")
    test_group.add_argument(
        "--test",
        help="Filter to the first N partitions of full VDS for testing (N controlled by --test-partitions).",
        action="store_true",
    )
    test_group.add_argument(
        "--test-partitions",
        type=int,
        default=2,
        help="Number of partitions to use in test mode. Default is 2.",
    )

    # Processing step arguments.
    processing_group = parser.add_argument_group("processing steps")
    processing_group.add_argument(
        "--process-gnomad",
        help="Process gnomAD dataset for frequency calculations.",
        action="store_true",
    )
    processing_group.add_argument(
        "--process-aou",
        help="Process All of Us dataset for frequency calculations.",
        action="store_true",
    )
    processing_group.add_argument(
        "--use-all-sites-ans",
        help="Use all sites ANs in frequency calculations to avoid a densify.",
        action="store_true",
    )

    # Environment configuration.
    env_group = parser.add_argument_group("environment configuration")
    env_group.add_argument(
        "--environment",
        help="Environment to run in.",
        choices=["rwb", "batch", "dataproc"],
        default="rwb",
    )
    env_group.add_argument(
        "--tmp-dir-days",
        type=int,
        default=4,
        help="Number of days for temp directory retention. Default is 4.",
    )
    env_group.add_argument(
        "--gcp-billing-project",
        type=str,
        default="broad-mpg-gnomad",
        help="Google Cloud billing project for reading requester pays buckets.",
    )

    # Batch-specific configuration.
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
    )
    batch_group.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_group.add_argument(
        "--driver-cores",
        type=int,
        default=None,
        help="Number of cores for driver node.",
    )
    batch_group.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_group.add_argument(
        "--worker-cores",
        type=int,
        default=None,
        help="Number of cores for worker nodes.",
    )
    batch_group.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    batch_args = [
        "app_name",
        "driver_cores",
        "driver_memory",
        "worker_cores",
        "worker_memory",
    ]
    provided_batch_args = [arg for arg in batch_args if getattr(args, arg) is not None]
    if provided_batch_args and args.environment != "batch":
        parser.error(
            f"Batch configuration arguments ({', '.join('--' + a.replace('_', '-') for a in provided_batch_args)}) "
            f"require --environment=batch"
        )

    main(args)
