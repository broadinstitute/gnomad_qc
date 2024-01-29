"""Produce per-exome stats, per-genome stats, and variant type counts for v4."""
import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import TableResource, VersionedTableResource
from gnomad.utils.slack import slack_notifications

import gnomad_qc.v4.resources.meta as meta
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import calling_intervals, get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import (
    CURRENT_RELEASE,
    RELEASES,
    _release_root,
    release_sites,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("release")
logger.setLevel(logging.INFO)


def create_per_sample_counts_resource(
    data_type: str = "exomes",
    test: bool = False,
    overwrite: bool = False,
    filter_variant_qc: bool = False,
    ukb_regions: bool = False,
    broad_regions: bool = False,
    suffix: str = None,
) -> None:
    """
    Write out Hail's sample_qc output for desired samples and variants.
    Useful for finding the amount of variants (total or specific types) per sample.
    Also prints distribution of total variants called per sample.

    :param data_type: String of either "exomes" or "genomes" for the sample type.
    :param test: Boolean for if you would like to use a small test set.
    :param overwrite: Boolean to overwrite checkpoint files if requested.
    :param filter_variant_qc: Bool to only count variants that pass all variant qc filters (i.e. VQSR, AC0).
    :param filter_calling_intervals: Bool to filter only to variants in calling intervals.
    :param singletons_by_pop: Bool to calculate the mean number of singletons by population group.
    :param calling: String of how to filter to calling intervals.
    :param suffix: String of arguments to append to name.
    """

    # Read in release data and metadata and VDS.
    ht_release_data = release_sites(data_type).ht()
    if test:
        logger.info("Test: filtering to only chr22 sites and 467 test samples in VDS")
        ht_release_data = ht_release_data.filter(
            ht_release_data.locus.contig == "chr22"
        )

    if filter_variant_qc:
        logger.info("Filter to those which pass all filters")
        ht_release_data = ht_release_data.filter(hl.len(ht_release_data.filters) == 0)

    meta_ht = meta(data_type=data_type).ht()
    vds = get_gnomad_v4_vds(test=test, release_only=True)

    # Filter to only variants and samples in release data.
    logger.info("Filter to only variants and samples in gnomAD release.")
    vds = hl.vds.filter_variants(
        vds,
        vds.variant_data.filter_rows(
            ~hl.is_missing(ht_release_data[vds.variant_data.row_key])
        )
        .rows()
        .select(),
    )

    vds = hl.vds.filter_samples(
        vds,
        vds.reference_data.filter_cols(meta_ht[vds.reference_data.s].release)
        .cols()
        .select(),
    )

    # Filter to variants only in requested calling intervals.
    if ukb_regions:
        logger.info("Filter to only variants in UKB capture regions.")
        filtered_ukb = (
            vds.variant_data.filter_rows(
                ht_release_data[
                    vds.variant_data.row_key
                ].region_flags.outside_ukb_capture_region
            )
            .rows()
            .select()
        )
        vds = hl.vds.filter_variants(vds, filtered_ukb)

    if broad_regions:
        logger.info("Filter to only variants in Broad capture regions.")
        filtered_broad = (
            vds.variant_data.filter_rows(
                ht_release_data[
                    vds.variant_data.row_key
                ].region_flags.outside_broad_capture_region
            )
            .rows()
            .select()
        )
        vds = hl.vds.filter_variants(vds, filtered_broad)

    # Perform Hail's VDS Sample QC module.
    # Resulting table is of most interest for further analysis.
    # Write said table out to determined path.
    logger.info("Running Hail's VDS Sample QC method on filtered VDS.")
    sample_qc_ht = hl.vds.sample_qc(vds)
    sample_qc_ht = sample_qc_ht.checkpoint(
        get_per_sample_counts(test=test, data_type=data_type, suffix=suffix).path,
        overwrite=overwrite,
    )


def aggregate_and_stratify(
    data_type: str = "exomes",
    test: bool = False,
    singletons_by_pop: bool = False,
    suffix: str = None,
) -> None:
    """
    Method to read in Sample QC table and output relevant summary stats

    :param data_type: String to read in either "exomes" or "genomes".
    :param test: Bool of whether or not to read in a test output.
    :param singletons_by_pop: Bool to calculate number of singletons by v4 pop.
    :param suffix: String of suffix of Sample QC table to read in.
    """
    logger.info("Reading Sample QC HT in from resource")
    sample_qc_ht = hl.read_table(
        get_per_sample_counts(test=test, data_type=data_type, suffix=suffix).path
    )
    logger.info(
        "via:"
        f" {get_per_sample_counts(test=test,data_type=data_type,suffix=suffix).path}"
    )

    # Column 'n_non_ref' is a per-sample metric of the number of variants called.
    # Report some stats of interest from this table.
    called_per_distribution = sample_qc_ht.aggregate(
        hl.struct(
            mean=hl.agg.mean(sample_qc_ht.n_non_ref),
            quantiles=hl.agg.approx_quantiles(
                sample_qc_ht.n_non_ref, qs=[0, 0.25, 0.50, 0.75, 1.0]
            ),
        )
    )

    logger.info(
        f"Quantiles called per {data_type} by 0th, 25th, 50th, 75th, and 100th"
        f" quantiles \n: {called_per_distribution.quantiles}"
        f" and a mean of {called_per_distribution.mean}"
    )

    # Calculate number of singletons by inputed ancestry.
    if singletons_by_pop:
        logger.info("Reading in Meta HT.")
        meta_ht = meta(data_type=data_type).ht()

        # Add population inference, as keyed by 's'.
        sample_qc_ht = sample_qc_ht.annotate(
            pop=meta_ht[sample_qc_ht.s].population_inference.pop
        )

        # Aggregate and print outputs.
        dictionary_results = sample_qc_ht.aggregate(
            hl.agg.group_by(
                sample_qc_ht.pop,
                hl.struct(
                    mean=hl.agg.mean(sample_qc_ht.n_singleton),
                    quantiles=hl.agg.approx_quantiles(
                        sample_qc_ht.n_singleton,
                        qs=[0, 0.25, 0.50, 0.75, 1.0],
                    ),
                ),
            )
        )

        logger.info(dictionary_results)

        for pop_record in dictionary_results.items():
            pop_label = pop_record[0]
            mean = pop_record[1].mean
            quantiles = pop_record[1].quantiles
            logger.info(
                f"For pop {pop_label}, mean # singletons called is {mean}"
                " with distribution by 0th, 25th, 50th, 75th, and 100th"
                f" percentiles as: {quantiles}"
            )


def variant_types():
    # what is needed: stats from your notebook for:
    # Indels and Snps per data type
    # How many (number and proportion) pass VQSR and pass All Filters
    pass


def vep_consequences():
    # what is needed: stats from notebook for:
    # VEP consequence breakdown , in those two different returns
    # let users specify: pass all filters, pass VQSR, or all variants
    pass


def main(args):
    """Collect summary stats of gnomAD v4 data."""
    hl.init(
        log="/per_exome.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    data_type = args.data_type
    test = args.test
    overwrite = args.overwrite
    filter_variant_qc = args.filter_variant_qc
    ukb_regions = args.ukb_regions
    broad_regions = args.broad_regions
    singletons_by_pop = args.singletons_by_pop

    suffix = (
        f"{'ukb_calling' if ukb_regions else ''}.{'broad_calling' if broad_regions else ''}.{'vqc' if filter_variant_qc else ''}"
    )

    if args.create_per_sample_counts_resource:
        create_per_sample_counts_resource(
            data_type=data_type,
            test=test,
            overwrite=overwrite,
            filter_variant_qc=filter_variant_qc,
            ukb_regions=ukb_regions,
            broad_regions=broad_regions,
            suffix=suffix,
        )

    if args.aggregate_and_stratify:
        aggregate_and_stratify(
            data_type=data_type,
            test=test,
            singletons_by_pop=singletons_by_pop,
            suffix=suffix,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite data (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--test",
        help="Test on a small number of variants",
        action="store_true",
    )
    parser.add_argument(
        "--data-type",
        default="exomes",
        choices=["exomes", "genomes"],
        help="Data type (exomes or genomes) to produce summary stats for.",
    )
    parser.add_argument(
        "--create-per-sample-counts-resource",
        help="Run code and output to create_per_sample_counts_resource",
        action="store_true",
    )
    parser.add_argument(
        "--filter-variant-qc",
        help=(
            "Only calculate per-sample variants for those variants that pass all"
            " variant qc filters."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--ukb-regions",
        help="Filter to variants only in UKB capture regions",
        action="store_true",
    )
    parser.add_argument(
        "--broad-regions",
        help="Filter to variants only in Broad capture regions",
        action="store_true",
    )
    parser.add_argument(
        "--aggregate-and-stratify",
        help="Run code to aggregate and stratify (exising?) sample qc table",
        action="store_true",
    )
    parser.add_argument(
        "--singletons-by-pop",
        help="Output statistics for number of singletons by population group",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
