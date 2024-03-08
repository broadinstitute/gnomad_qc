"""Script to get per-sample stats by variant type in v4."""
import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import TableResource, VersionedTableResource
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)
from hail.utils.misc import new_temp_file

import gnomad_qc.v4.resources.meta as meta
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import get_per_sample_counts, release_sites
from gnomad_qc.v4.resources.temp_hail_methods import (
    vmt_sample_qc,
    vmt_sample_qc_variant_annotations,
)

# from hail.vds.sample_qc import vmt_sample_qc, vmt_sample_qc_variant_annotations


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("release")
logger.setLevel(logging.INFO)


def create_per_sample_counts(
    data_type: str = "exomes",
    test: bool = False,
    overwrite: bool = False,
    pass_filters: bool = False,
    ukb_regions: bool = False,
    broad_regions: bool = False,
    by_csqs: bool = False,
    rare_variants: bool = False,
) -> hl.Table:
    """
    Write out Hail's sample_qc output for desired samples and variants.

    Useful for finding the amount of variants (total or specific types) per sample.
    Also prints distribution of total variants called per sample.

    :param data_type: String of either "exomes" or "genomes" for the sample type.
    :param test: Boolean for if you would like to use a small test set.
    :param overwrite: Boolean to overwrite checkpoint files if requested.
    :param pass_filters: Filter to variants which pass all variant qc filters.
    :param ukb_regions: Filter to variants in UKB regions.
    :param broad_regions: Filter to variants in Broad regions.
    :param by_csqs: Stratify by variant consequence.
    :param rare_variants: Stratify by variants which have adj AF <0.1%.
    """
    vds = get_gnomad_v4_vds(test=test, release_only=True, split=True)

    # Read in release data and select relevant annotations.
    ht = release_sites(data_type=data_type).ht()
    ht = ht.select("freq", "filters", "region_flags", "vep")

    if test:
        logger.info("Test: filtering to variants on chr22")
        ht = ht.filter(ht.locus.contig == "chr22")

    # logger.info("Filter VDS to variants in gnomAD release.")
    # vds = hl.vds.filter_variants(vds, ht.select())

    # Create variant matrix table from VDS and annotate with all relevant info.
    vmt = vds.variant_data
    logger.info("Variant MT created and checkpointing, without annotations")
    vmt = vmt.checkpoint(
        new_temp_file(f'vmt_{"TEST" if test else ""}_mt', extension="mt"),
        overwrite=overwrite,
    )

    # Add extra Allele Count and Allele Type annotations to variant Matrix
    # Table, according to Hail standards, to help their computation.
    ac_and_atype = vmt_sample_qc_variant_annotations(
        global_gt=vmt.GT, alleles=vmt.alleles
    )
    vmt = vmt.annotate_rows(**ac_and_atype)

    if by_csqs:
        ht = filter_vep_transcript_csqs(
            ht,
            synonymous=False,
            filter_empty_csq=True,
            canonical=False,
            mane_select=True,
        )

        ht = get_most_severe_consequence_for_summary(ht)

        ht = ht.select("freq", "filters", "region_flags", "most_severe_csq", "lof")

    # NOTE: not size efficient, but we are never checkpointing THIS I don't believe
    vmt = vmt.annotate_rows(**ht[vmt.row_key])

    arg_dict = {"all_variants": ~hl.is_missing(vmt.locus)}

    if pass_filters:
        arg_dict.update({"pass_filters": hl.len(vmt.filters) == 0})
    if ukb_regions:
        arg_dict.update({"ukb_regions": ~vmt.region_flags.outside_ukb_capture_region})
    if broad_regions:
        arg_dict.update(
            {"broad_regions": ~vmt.region_flags.outside_broad_capture_region}
        )
    if ukb_regions & broad_regions:
        arg_dict.update(
            {
                "ukb_and_broad_regions": (
                    ~vmt.region_flags.outside_ukb_capture_region
                ) & (~vmt.region_flags.outside_broad_capture_region)
            }
        )
    if rare_variants:
        arg_dict.update({"rare_variants": vmt.freq[0].AF < 0.001})
    if by_csqs:
        arg_dict.update({"lof": ~hl.is_missing(vmt.lof)})
        arg_dict.update({"missense": vmt.most_severe_csq == "missense_variant"})
        arg_dict.update({"synonymous": vmt.most_severe_csq == "synonymous_variant"})

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
            for ann, expr in arg_dict.items()
        }
    )
    # Annotate stratified sample qc onto columns, then take and output columns
    vmt = vmt.annotate_cols(aggregated_vmt_sample_qc=qc)
    sample_qc_ht = vmt.cols()

    return sample_qc_ht


def aggregate_and_stratify(
    data_type: str = "exomes",
    test: bool = False,
    overwrite: bool = False,
    agg_variant_qc: bool = False,
    agg_ukb_regions: bool = False,
    agg_broad_regions: bool = False,
    agg_lof: bool = False,
    agg_rare_variants: bool = False,
    counts_by_pop: bool = False,
    suffix: str = None,
    write_counts: bool = False,
) -> None:
    """
    Read in per-sample variant counts Table and output relevant summary stats.

    :param data_type: String to read in either "exomes" or "genomes".
    :param test: Bool of whether or not to read in a test output.
    :param Overwrite: Write over existing file.
    :param agg_variant_qc: Stratify by variants which pass variant qc.
    :param agg_ukb_regions: Stratify by variants in UKB regions.
    :param agg_broad_regions: Stratify by variants in Broad regions.
    :param agg_lof: Stratify by variants which are loss-of-function.
    :param agg_rare_variants: Stratify by variants which have adj AF <0.1%.
    :param counts_by_pop: Bool to calculate number of singletons, n_het, and n_hom by v4 pop.
    :param suffix: String of suffix of Sample QC table to read in.
    :param write_counts: Write out Hail dict of results.
    """
    logger.info("Reading Sample QC HT in from resource")
    sample_qc_ht = hl.read_table(
        get_per_sample_counts(test=test, data_type=data_type, suffix=suffix).path
    )
    logger.info(
        "via:"
        f" {get_per_sample_counts(test=test, data_type=data_type, suffix=suffix).path}"
    )

    # Define internal method to return all desired counts in a struct.
    # This code would be repeated in code otherwise.
    def _count_struct(input_ht, input_strat):
        ret_struct = hl.struct(
            mean_nonref=hl.agg.mean(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_non_ref
            ),
            quantiles_nonref=hl.agg.approx_quantiles(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_non_ref,
                qs=[0, 0.25, 0.50, 0.75, 1.0],
            ),
            mean_singletons=hl.agg.mean(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_singleton
            ),
            quantiles_singletons=hl.agg.approx_quantiles(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_singleton,
                qs=[0, 0.25, 0.50, 0.75, 1.0],
            ),
            mean_het=hl.agg.mean(input_ht.aggregated_vmt_sample_qc[input_strat].n_het),
            quantiles_het=hl.agg.approx_quantiles(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_het,
                qs=[0, 0.25, 0.50, 0.75, 1.0],
            ),
            mean_hom=hl.agg.mean(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_hom_var
            ),
            quantiles_hom=hl.agg.approx_quantiles(
                input_ht.aggregated_vmt_sample_qc[input_strat].n_hom_var,
                qs=[0, 0.25, 0.50, 0.75, 1.0],
            ),
        )

        return ret_struct

    # Print stats about an aggregation
    # MUST be generated by _count_struct or have identical formatting
    def _logger_struct(
        input_struct: hl.struct, input_strat: str, input_pop: str = "ALL"
    ) -> None:
        logger.info(
            f"For pop {input_pop} in stratification {input_strat}:"
            f"\n\tmean nonref as {input_struct.mean_nonref}"
            f"\n\tmean singletons as {input_struct.mean_singletons}"
            f"\n\tmean het as {input_struct.mean_het}"
            f"\n\tand mean hom as {input_struct.mean_hom}"
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
        logger.info("Reading in meta HT.")
        meta_ht = meta(data_type=data_type).ht()

    # Dictionary to write out if user desires.
    returned_counts = {}

    for strat in stratifications:
        logger.info(f"STRATIFICATION: {strat}")
        called_per_distribution = sample_qc_ht.aggregate(
            _count_struct(sample_qc_ht, strat)
        )

        _logger_struct(input_struct=called_per_distribution, input_strat=strat)

        returned_counts[f"{strat}_all_samples"] = called_per_distribution

        # Calculate number of singletons by inputed ancestry by stratification
        if counts_by_pop:
            # Add population inference, as keyed by 's'.
            sample_qc_ht = sample_qc_ht.annotate(
                pop=meta_ht[sample_qc_ht.s].population_inference.pop
            )

            # Aggregate and print outputs.
            dictionary_results = sample_qc_ht.aggregate(
                hl.agg.group_by(sample_qc_ht.pop, _count_struct(sample_qc_ht, strat))
            )

            for pop_record in dictionary_results.items():
                pop_label = pop_record[0]
                pop_struct = pop_record[1]
                _logger_struct(
                    input_struct=pop_struct, input_strat=strat, input_pop=pop_label
                )

                returned_counts[f"{strat}_{pop_label}_specific"] = pop_struct


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
    write_counts = args.write_counts

    if args.create_per_sample_counts_resource:
        per_sample_ht = create_per_sample_counts(
            data_type=data_type,
            test=test,
            pass_filters=filter_variant_qc,
            ukb_regions=filter_ukb_regions,
            broad_regions=filter_broad_regions,
            by_csqs=True,
            rare_variants=agg_rare_variants,
        )

        per_sample_ht = per_sample_ht.checkpoint(
            get_per_sample_counts(test=test, data_type=data_type, suffix=suffix).path,
            overwrite=overwrite,
        )

    if args.aggregate_and_stratify:
        aggregate_and_stratify(
            data_type=data_type,
            test=test,
            overwrite=overwrite,
            agg_variant_qc=agg_variant_qc,
            agg_ukb_regions=agg_ukb_regions,
            agg_broad_regions=agg_broad_regions,
            agg_lof=agg_lof,
            agg_rare_variants=agg_rare_variants,
            counts_by_pop=counts_by_pop,
            suffix=suffix,
            write_counts=write_counts,
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
    parser.add_argument(
        "--write-counts",
        action="store_true",
        help="Write counts to dictionary in tmp dir",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
