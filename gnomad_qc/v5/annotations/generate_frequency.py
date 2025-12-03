"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies for samples that will be removed from gnomAD v4
in v5 due to consent withdrawal.

Processing Workflow:
----------------------
1. process_gnomad_dataset: Update v4 frequency table by subtracting consent withdrawal samples
   - Load v4 frequency table (contains frequencies and age histograms)
   - Prepare consent withdrawal VDS (split multiallelics, annotate metadata)
   - Calculate frequencies and age histograms for consent samples
   - Subtract from v4 frequencies to get updated gnomAD v5 frequencies


Usage Examples:
---------------
# Process gnomAD consent withdrawals
python generate_frequency.py --process-gnomad --environment dataproc

# Run in test mode
python generate_frequency.py --process-gnomad --test --test-partitions 2
"""

import argparse
import logging

import hail as hl
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    compute_freq_by_strata,
    get_adj_expr,
    merge_freq_arrays,
    merge_histograms,
)
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.annotations import get_freq as get_v4_freq
from gnomad_qc.v5.resources.annotations import get_freq, group_membership
from gnomad_qc.v5.resources.basics import (
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def _prepare_consent_vds(
    v4_freq_ht: hl.Table,
    test: bool = False,
    test_partitions: int = 2,
) -> hl.vds.VariantDataset:
    """
    Load and prepare VDS for consent withdrawal sample processing.

    :param v4_freq_ht: v4 frequency table for AF annotation.
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
    vmt.describe()
    vmt = vmt.select_cols(
        gen_anc=vmt.meta.gen_anc,
        sex_karyotype=vmt.meta.sex_karyotype,
        age=vmt.meta.age,
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
    # Group membership table is already filtered to consent drop samples.
    group_membership_ht = group_membership(test=test, data_set="gnomad").ht()

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

    joined_freq_ht = v4_freq_ht.annotate(
        consent_freq=consent_freq_ht[v4_freq_ht.key].freq,
        consent_age_hists=consent_freq_ht[v4_freq_ht.key].age_hists,
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
            ]
        },
    )
    # Update the frequency table with all changes.
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

    # Update the frequency table with all changes
    joined_freq_ht = joined_freq_ht.annotate(
        histograms=joined_freq_ht.histograms.annotate(
            age_hists=joined_freq_ht.histograms.age_hists.annotate(
                age_hist_het=updated_age_hist_het,
                age_hist_hom=updated_age_hist_hom,
            )
        ),
    )

    return joined_freq_ht.checkpoint(new_temp_file("merged_freq_and_hists", "ht"))


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
    # Test is used for formatting file paths.
    v4_freq_ht = get_v4_freq(data_type="genomes").ht()

    vds = _prepare_consent_vds(
        v4_freq_ht,
        test=test,
        test_partitions=test_partitions,
    )

    logger.info("Calculating frequencies and age histograms for consent samples...")
    consent_freq_ht = _calculate_consent_frequencies_and_age_histograms(vds, test)

    if test:
        v4_freq_ht = v4_freq_ht.filter(hl.is_defined(consent_freq_ht[v4_freq_ht.key]))

    logger.info("Subtracting consent frequencies and age histograms from v4...")
    updated_freq_ht = _subtract_consent_frequencies_and_age_histograms(
        v4_freq_ht, consent_freq_ht, test
    )

    # Only overwrite fields that were actually updated (merge back with original table)
    # Note: FAF/grpmax annotations will be calculated on the final merged dataset.
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
        # Add optional batch configuration parameters
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
            default_reference="GRCh38",  # Remove if running in 0.2.134 or later
        )
    # hl.default_reference("GRCh38") # Testing in hail 0.2.130 so commenting out for now


def main(args):
    """Generate v5 frequency data."""
    environment = args.environment
    test = args.test
    test_partitions = args.test_partitions
    overwrite = args.overwrite
    tmp_dir_days = args.tmp_dir_days

    _initialize_hail(args)

    try:
        logger.info("Running generate_frequency.py...")
        if args.process_gnomad:
            logger.info("Processing gnomAD dataset...")

            gnomad_freq = get_freq(test=test, data_type="genomes", data_set="gnomad")

            check_resource_existence(
                output_step_resources={"process-gnomad": [gnomad_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = process_gnomad_dataset(
                test=test,
                test_partitions=test_partitions,
            )

            logger.info(
                f"Writing gnomAD frequency HT (with embedded age histograms) to {gnomad_freq.path}..."
            )
            gnomad_freq_ht.write(gnomad_freq.path, overwrite=overwrite)

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

    # General arguments
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )

    # Test/debug arguments
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

    # Processing step arguments
    processing_group = parser.add_argument_group("processing steps")
    processing_group.add_argument(
        "--process-gnomad",
        help="Process gnomAD dataset for frequency calculations.",
        action="store_true",
    )

    # Environment configuration
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

    # Batch-specific configuration
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch)",
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
    main(args)
