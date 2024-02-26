"""Produce per-exome stats, per-genome stats, and variant type counts for v4."""
import argparse
import logging

import hail as hl
from hail.utils.misc import new_temp_file

from gnomad.resources.resource_utils import TableResource, VersionedTableResource
from gnomad.utils import vep
from gnomad.utils.slack import slack_notifications

import gnomad_qc.v4.resources.meta as meta
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import (
    CURRENT_RELEASE,
    RELEASES,
    _release_root,
    get_per_sample_counts,
    release_sites,
)
from gnomad_qc.v4.resources.temp_hail_methods import (
    vmt_sample_qc,
    vmt_sample_qc_variant_annotations,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("release")
logger.setLevel(logging.INFO)

from typing import Optional, Sequence


def create_per_sample_counts_resource(
    data_type: str = "exomes",
    test: bool = False,
    overwrite: bool = False,
    filter_variant_qc: bool = False,
    filter_ukb_regions: bool = False,
    filter_broad_regions: bool = False,
    agg_variant_qc: bool = False,
    agg_ukb_regions: bool = False,
    agg_broad_regions: bool = False,
    agg_lof: bool = False,
    agg_rare_variants: bool = False,
    suffix: str = None,
) -> None:
    """
    Write out Hail's sample_qc output for desired samples and variants.
    Useful for finding the amount of variants (total or specific types) per sample.
    Also prints distribution of total variants called per sample.

    :param data_type: String of either "exomes" or "genomes" for the sample type.
    :param test: Boolean for if you would like to use a small test set.
    :param overwrite: Boolean to overwrite checkpoint files if requested.
    :param filter_variant_qc: Filter to only variants which pass qc before any aggregation.
    :param filter_ukb_regions: Filter to only variants in UKB regions before any aggregation.
    :param filter_broad_regions: Filter to only variants in Broad regions before any aggregation.
    :param agg_variant_qc: Stratify by variants which pass variant qc.
    :param agg_ukb_regions: Stratify by variants in UKB regions.
    :param agg_broad_regions: Stratify by variants in Broad regions.
    :param agg_lof: Stratify by variants which are loss-of-function.
    :param agg_rare_variants: Stratify by variants which have adj AF <0.1%.
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
        logger.info("Filter to those which pass all variant qc filters")
        ht_release_data = ht_release_data.filter(hl.len(ht_release_data.filters) == 0)

    meta_ht = meta(data_type=data_type).ht()
    vds = get_gnomad_v4_vds(test=test, release_only=True, split=True)

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
    if filter_ukb_regions:
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

    if filter_broad_regions:
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

    # Create variant matrix table from VDS and annotate with all relevant info
    vmt = vds.variant_data
    logger.info("VDS created and checkpointing, without annotations")
    vmt = vmt.checkpoint(
        new_temp_file(f'vmt_{"TEST" if test else ""}_mt', extension="mt"),
        overwrite=overwrite
    )

    # Add extra Allele Count and A-Type, according to Hail standards, to help their computation
    ac_and_atype = vmt_sample_qc_variant_annotations(
        global_gt=vmt.GT, alleles=vmt.alleles
    )
    vmt = vmt.annotate_rows(**ac_and_atype)

    # NOTE: not size efficient, but we are never checkpointing THIS I don't believe
    vmt = vmt.annotate_rows(**ht_release_data[vmt.row_key])

    if agg_lof:
        logger.info("Looking at variants that are VEP LOF")
        vmt_vep = vep.filter_vep_transcript_csqs(
            vmt,
            synonymous=False,
            filter_empty_csq=True,
            canonical=False,
            mane_select=True,
        )

        ht_vep = vep.get_most_severe_consequence_for_summary(vmt_vep.rows())

        vmt = vmt.annotate_rows(
            **ht_vep[vmt.row_key].select(
                "most_severe_csq", "protein_coding", "lof", "no_lof_flags"
            )
        )

    # vmt = vmt.annotate_rows(all_true=True)
    # Odd Hail behavior: all annotations need to be defined before an argument_dictionary is constructed.
    # Or else it throws an error about differing sources.
    argument_dictionary = {"all_variants": ~hl.is_missing(vmt.locus)}

    if agg_variant_qc:
        argument_dictionary["pass_all_vqc"] = hl.len(vmt.filters) == 0
    if agg_ukb_regions:
        argument_dictionary["inside_ukb_regions"] = (
            ~vmt.region_flags.outside_ukb_capture_region
        )
    if agg_broad_regions:
        argument_dictionary["inside_broad_regions"] = (
            ~vmt.region_flags.outside_broad_capture_region
        )
    if agg_rare_variants:
        argument_dictionary["rare_variants"] = vmt.freq[0].AF < 0.001
    if agg_lof:
        argument_dictionary["is_lof"] = ~hl.is_missing(vmt.lof)

    # Run Hail's stratified vmt_sample_qc for selected stratifications.
    qc = hl.struct(
        **{
            ann: hl.agg.filter(
                expr,
                vmt_sample_qc(
                    global_gt=vmt.GT,
                    gq=vmt.GQ,
                    variant_ac=vmt.variant_ac,
                    variant_atypes=vmt.variant_atypes,
                    dp=vmt.DP,
                ),
            )
            for ann, expr in argument_dictionary.items()
        }
    )
    # Annotate stratified sample qc onto columns, then take and output columns
    vmt = vmt.annotate_cols(aggregated_vmt_sample_qc=qc)
    sample_qc_ht = vmt.cols()

    sample_qc_ht = sample_qc_ht.checkpoint(
        get_per_sample_counts(test=test, data_type=data_type, suffix=suffix).path,
        overwrite=overwrite,
    )


def aggregate_and_stratify(
    data_type: str = "exomes",
    test: bool = False,
    agg_variant_qc: bool = False,
    agg_ukb_regions: bool = False,
    agg_broad_regions: bool = False,
    agg_lof: bool = False,
    agg_rare_variants: bool = False,
    counts_by_pop: bool = False,
    suffix: str = None,
) -> None:
    """
    Method to read in Sample QC table and output relevant summary stats

    :param data_type: String to read in either "exomes" or "genomes".
    :param test: Bool of whether or not to read in a test output.
    :param agg_variant_qc: Stratify by variants which pass variant qc.
    :param agg_ukb_regions: Stratify by variants in UKB regions.
    :param agg_broad_regions: Stratify by variants in Broad regions.
    :param agg_lof: Stratify by variants which are loss-of-function.
    :param agg_rare_variants: Stratify by variants which have adj AF <0.1%.
    :param counts_by_pop: Bool to calculate number of singletons, n_het, and n_hom by v4 pop.
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

    # Create list of calls to stratify by
    stratifications = ["all_variants"]

    if agg_variant_qc:
        stratifications.append("pass_all_vqc")
    if agg_ukb_regions:
        stratifications.append("inside_ukb_regions")
    if agg_broad_regions:
        stratifications.append("inside_broad_regions")
    if agg_rare_variants:
        stratifications.append("rare_variants")
    if agg_lof:
        stratifications.append("is_lof")

    # Column 'n_non_ref' is a per-sample metric of the number of variants called.
    # Report some stats of interest from this table.

    if counts_by_pop:
        logger.info("Reading in Meta HT.")
        meta_ht = meta(data_type=data_type).ht()

    for strat in stratifications:
        logger.info(f"STRATIFICATION: {strat}")
        called_per_distribution = sample_qc_ht.aggregate(
            hl.struct(
                mean=hl.agg.mean(
                    sample_qc_ht.aggregated_vmt_sample_qc[strat].n_non_ref
                ),
                quantiles=hl.agg.approx_quantiles(
                    sample_qc_ht.aggregated_vmt_sample_qc[strat].n_non_ref,
                    qs=[0, 0.25, 0.50, 0.75, 1.0],
                ),
            )
        )

        logger.info(
            f"For the stratification of {strat} in all ancestries:"
            f"\n\tmean non-ref calls of {called_per_distribution.mean:.2f}"
            f"\n\tQuantiles called per {data_type} by 0th, 25th, 50th, 75th, and 100th:"
            f"\n\t{called_per_distribution.quantiles}"
        )

        # Calculate number of singletons by inputed ancestry by stratification
        # We don't care about singletons by pop by the stratifiactions,
        # But we do care about n_het and n_hom
        if counts_by_pop:
            # Add population inference, as keyed by 's'.
            sample_qc_ht = sample_qc_ht.annotate(
                pop=meta_ht[sample_qc_ht.s].population_inference.pop
            )

            # Aggregate and print outputs.
            dictionary_results = sample_qc_ht.aggregate(
                hl.agg.group_by(
                    sample_qc_ht.pop,
                    hl.struct(
                        mean_singletons=hl.agg.mean(
                            sample_qc_ht.aggregated_vmt_sample_qc[strat].n_singleton
                        ),
                        quantiles_singletons=hl.agg.approx_quantiles(
                            sample_qc_ht.aggregated_vmt_sample_qc[strat].n_singleton,
                            qs=[0, 0.25, 0.50, 0.75, 1.0],
                        ),
                        mean_het=hl.agg.mean(
                            sample_qc_ht.aggregated_vmt_sample_qc[strat].n_het
                        ),
                        quantiles_het=hl.agg.approx_quantiles(
                            sample_qc_ht.aggregated_vmt_sample_qc[strat].n_het,
                            qs=[0, 0.25, 0.50, 0.75, 1.0],
                        ),
                        mean_hom=hl.agg.mean(
                            sample_qc_ht.aggregated_vmt_sample_qc[strat].n_hom_var
                        ),
                        quantiles_hom=hl.agg.approx_quantiles(
                            sample_qc_ht.aggregated_vmt_sample_qc[strat].n_hom_var,
                            qs=[0, 0.25, 0.50, 0.75, 1.0],
                        ),
                    ),
                )
            )

            for pop_record in dictionary_results.items():
                pop_label = pop_record[0]
                mean_singletons = pop_record[1].mean_singletons
                quantiles_singletons = pop_record[1].quantiles_singletons
                mean_het = pop_record[1].mean_het
                mean_hom = pop_record[1].mean_hom
                logger.info(
                    f"For pop {pop_label} in stratification {strat}:\n\tmean #"
                    f" singletons called is {mean_singletons:.2f} \n\twith distribution"
                    " by 0th, 25th, 50th, 75th, and 100th percentiles as:"
                    f" {quantiles_singletons} \n\twith mean (n_het ,n_hom) of:"
                    f" ({mean_het:.2f} , {mean_hom:.2f})"
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
    filter_ukb_regions = args.filter_ukb_regions
    filter_broad_regions = args.filter_broad_regions
    agg_variant_qc = args.agg_variant_qc
    agg_ukb_regions = args.agg_ukb_regions
    agg_broad_regions = args.agg_broad_regions
    agg_lof = args.agg_lof
    agg_rare_variants = args.agg_rare_variants
    counts_by_pop = args.counts_by_pop

    suffix = (
        f"{'ukb_calling.' if filter_ukb_regions else ''}{'broad_calling.' if filter_broad_regions else ''}{'vqc.' if filter_variant_qc else ''}{args.custom_suffix if args.custom_suffix else ''}"
    )

    if args.create_per_sample_counts_resource:
        create_per_sample_counts_resource(
            data_type=data_type,
            test=test,
            overwrite=overwrite,
            filter_variant_qc=filter_variant_qc,
            filter_ukb_regions=filter_ukb_regions,
            filter_broad_regions=filter_broad_regions,
            agg_variant_qc=agg_variant_qc,
            agg_ukb_regions=agg_ukb_regions,
            agg_broad_regions=agg_broad_regions,
            agg_lof=agg_lof,
            agg_rare_variants=agg_rare_variants,
            suffix=suffix,
        )

    if args.aggregate_and_stratify:
        aggregate_and_stratify(
            data_type=data_type,
            test=test,
            agg_variant_qc=agg_variant_qc,
            agg_ukb_regions=agg_ukb_regions,
            agg_broad_regions=agg_broad_regions,
            agg_lof=agg_lof,
            agg_rare_variants=agg_rare_variants,
            counts_by_pop=counts_by_pop,
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
        "--filter-ukb-regions",
        help="Filter to variants only in UKB capture regions",
        action="store_true",
    )
    parser.add_argument(
        "--filter-broad-regions",
        help="Filter to variants only in Broad capture regions",
        action="store_true",
    )
    parser.add_argument(
        "--agg-broad-regions",
        help="Stratify by variants in Broad regions.",
        action="store_true",
    )
    parser.add_argument(
        "--agg-ukb-regions",
        help="Stratify by variants in UKB regions.",
        action="store_true",
    )
    parser.add_argument(
        "--agg-rare-variants",
        help="Stratify by variants which have adj AF <0.1%.",
        action="store_true",
    )
    parser.add_argument(
        "--agg-lof",
        help="Stratify by variants which are loss-of-function.",
        action="store_true",
    )
    parser.add_argument(
        "--agg-variant-qc",
        help="Stratify by variants which pass variant qc.",
        action="store_true",
    )
    parser.add_argument(
        "--aggregate-and-stratify",
        help="Run code to aggregate and stratify (exising?) sample qc table",
        action="store_true",
    )
    parser.add_argument(
        "--counts-by-pop",
        help=(
            "Output statistics for number of singletons, n_het, and n_hom by population"
            " group"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--custom-suffix",
        type=str,
        default=None,
        help="Custom string to append to output names",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
