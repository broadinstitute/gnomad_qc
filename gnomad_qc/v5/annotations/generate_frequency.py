"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies and age histograms for the AoU dataset
using either pre-computed allele numbers or a densify approach.

Processing Workflow:
--------------------
1. Load AoU VDS with metadata
2. Prepare VDS (annotate group membership, adjust for ploidy, split multi-allelics)
3. Calculate frequencies using either:
   - All sites ANs (efficient, requires pre-computed AN values)
   - Densify approach (standard, more resource intensive)
4. Generate age histograms during frequency calculation

Usage Examples:
---------------
# Process AoU dataset using all-sites ANs.
python generate_frequency.py --process-aou --use-all-sites-ans --environment rwb

# Process AoU on batch/QoB with custom resources.
python generate_frequency.py --process-aou --environment batch --app-name "aou_freq" --driver-cores 8 --worker-memory highmem
"""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import (
    age_hists_expr,
    agg_by_strata,
    compute_freq_by_strata,
    qual_hist_expr,
)

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj as annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_freq,
    group_membership,
)
from gnomad_qc.v5.resources.basics import get_aou_vds, get_logging_path, qc_temp_prefix

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def mt_hists_fields(mt: hl.MatrixTable) -> hl.StructExpression:
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
    aou_vds: hl.vds.VariantDataset, test: bool = False
) -> hl.vds.VariantDataset:
    """
    Prepare AoU VDS for frequency calculations.

    :param aou_vds: AoU VariantDataset.
    :param test: Whether running in test mode.
    :return: Prepared AoU VariantDataset.
    """
    logger.info(f"Using test mode: {test}")
    aou_vmt = aou_vds.variant_data
    # Use existing AoU group membership table and filter to variant samples.
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht = group_membership(test=test, data_set="aou").ht()
    group_membership_globals = group_membership_ht.index_globals()
    # Ploidy is already adjusted in the AoU VDS because of DRAGEN, do not need
    # to adjust it here.
    aou_vmt = aou_vmt.annotate_cols(
        sex_karyotype=aou_vmt.meta.sex_karyotype,
        gen_anc=aou_vmt.meta.genetic_ancestry_inference.gen_anc,
        age=aou_vmt.meta.project_meta.age,
        group_membership=group_membership_ht[aou_vmt.col_key].group_membership,
    )
    aou_vmt = aou_vmt.annotate_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
    )

    # Add adj annotation required by annotate_freq.
    aou_vmt = annotate_adj_no_dp(aou_vmt)

    # Rename LGT to GT and LAD to AD for compatibility with annotate_freq.
    aou_vds = hl.vds.split_multi(hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt))

    return aou_vds


def _calculate_aou_frequencies_and_hists_using_all_sites_ans(
    aou_variant_mt: hl.MatrixTable, test: bool = False
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using all sites ANs.

    :param aou_variant_mt: Prepared variant MatrixTable.
    :param test: Whether to use test resources.
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    all_sites_an_ht = coverage_and_an_path(test=test).ht()
    aou_variant_mt = aou_variant_mt.annotate_rows(
        hists_fields=mt_hists_fields(aou_variant_mt)
    )
    # Add all of the qual hists from the consent_ans table into the
    # aou_variant_freq_ht qual_hists struct
    aou_variant_mt = aou_variant_mt.annotate_rows(
        qual_hists=aou_variant_mt.qual_hists.annotate(
            **all_sites_an_ht[aou_variant_mt.row_key].qual_hists,
        )
    )
    logger.info("Annotating frequencies with all sites ANs...")
    group_membership_ht = group_membership(test=test, data_set="aou").ht()
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
        select_fields=["hists_fields"],
    )

    # Load AN values from consent_ans (calculated by another script but used
    # same group membership HT so same strata order)
    logger.info("Loading AN values from consent_ans...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        all_sites_an=all_sites_an_ht[aou_variant_freq_ht.key].AN
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
            aou_variant_freq_ht.all_sites_an,
        ),
    ).drop("all_sites_an")

    # Nest histograms to match gnomAD structure
    # Note: hists_fields.qual_hists already contains raw_qual_hists and qual_hists
    # as nested fields due to split_adj_and_raw=True in qual_hist_expr
    aou_variant_freq_ht = aou_variant_freq_ht.select(
        freq=aou_variant_freq_ht.freq,
        histograms=hl.struct(
            qual_hists=aou_variant_freq_ht.hists_fields.qual_hists.qual_hists,
            raw_qual_hists=aou_variant_freq_ht.hists_fields.qual_hists.raw_qual_hists,
            age_hists=aou_variant_freq_ht.hists_fields.age_hists,
        ),
    )

    return aou_variant_freq_ht


def _calculate_aou_frequencies_and_hists_using_densify(
    aou_vds: hl.vds.VariantDataset,
    test: bool = False,
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using densify.

    :param aou_vds: Prepared AoU VariantDataset.
    :param test: Whether to use test resources.
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    aou_mt = hl.vds.to_dense_mt(aou_vds)
    aou_mt = aou_mt.annotate_rows(hists_fields=mt_hists_fields(aou_mt))
    aou_freq_ht = compute_freq_by_strata(aou_mt, select_fields=["hists_fields"])

    # Nest histograms to match gnomAD structure
    # Note: hists_fields.qual_hists already contains raw_qual_hists and qual_hists
    # as nested fields due to split_adj_and_raw=True in qual_hist_expr
    aou_freq_ht = aou_freq_ht.transmute(
        histograms=hl.struct(
            qual_hists=aou_freq_ht.hists_fields.qual_hists.qual_hists,
            raw_qual_hists=aou_freq_ht.hists_fields.qual_hists.raw_qual_hists,
            age_hists=aou_freq_ht.hists_fields.age_hists,
        )
    )

    return aou_freq_ht


def process_aou_dataset(
    test: bool = False, use_all_sites_ans: bool = False
) -> hl.Table:
    """
    Process All of Us dataset for frequency calculations and age histograms.

    This function efficiently processes the AoU VDS by:
    1. Computing complete frequency struct using imported AN from consent_ans
    2. Generating age histograms within the frequency calculation

    :param test: Whether to run in test mode.
    :param use_all_sites_ans: Whether to use all sites ANs for frequency calculations.
    :return: Table with freq and age_hists annotations for AoU dataset.
    """
    logger.info("Processing All of Us dataset...")
    logger.info(f"Using test mode: {test}")
    aou_vds = get_aou_vds(annotate_meta=True, release_only=True, test=test)
    aou_vds = _prepare_aou_vds(aou_vds, test=test)

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

    return aou_freq_ht


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
            default_reference="GRCh38",
        )
    # hl.default_reference("GRCh38") Commenting out for now to test in hail 0.2.130


def main(args):
    """Generate v5 frequency data."""
    environment = args.environment
    data_test = args.data_test
    runtime_test = args.runtime_test
    use_all_sites_ans = args.use_all_sites_ans
    test = data_test or runtime_test
    overwrite = args.overwrite
    tmp_dir_days = args.tmp_dir_days

    _initialize_hail(args)

    try:
        logger.info("Running generate_frequency.py...")
        if args.process_aou:
            logger.info("Processing All of Us dataset...")
            logger.info(f"Using test mode: {test}")

            aou_freq = get_freq(test=test, data_type="genomes", data_set="aou")

            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq]},
                overwrite=overwrite,
            )

            aou_freq_ht = process_aou_dataset(
                test=test, use_all_sites_ans=use_all_sites_ans
            )

            logger.info(f"Writing AoU frequency HT to {aou_freq.path}...")
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

    # General arguments
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )

    # Test/debug arguments
    test_group = parser.add_argument_group("testing options")
    test_group.add_argument(
        "--data-test",
        help="Filter to the first N partitions of full VDS for testing (N controlled by --test-partitions).",
        action="store_true",
    )
    test_group.add_argument(
        "--runtime-test",
        help="Load test dataset and filter to test partitions.",
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
        "--process-aou",
        help="Process All of Us dataset for frequency calculations.",
        action="store_true",
    )
    processing_group.add_argument(
        "--use-all-sites-ans",
        help="Use all sites ANs in frequency calculations to avoid a densify.",
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
