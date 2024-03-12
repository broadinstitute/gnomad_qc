"""Script to get per-sample stats by variant type in v4."""
import argparse
import logging
from typing import List

import hail as hl
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources import meta
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import get_per_sample_counts, release_sites
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


def create_per_sample_counts(
    data_type: str = "exomes",
    test: bool = False,
    overwrite: bool = False,
    pass_filters: bool = False,
    ukb_regions: bool = False,
    broad_regions: bool = False,
    by_csqs: bool = False,
    rare_variants: bool = False,
    vep_canonical: bool = True,
    vep_mane: bool = False,
) -> hl.Table:
    """
    Write out Hail's sample_qc output for desired samples and variants.

    Useful for finding the amount of variants (total or specific types) per sample.
    Also prints distribution of total variants called per sample.

    :param data_type: String of either "exomes" or "genomes" for the sample type.
    :param test: Boolean for if you would like to use a small test set.
    :param overwrite: Boolean to overwrite checkpoint files if requested.
    :param pass_filters: Stratify by variants which pass all variant qc filters.
    :param ukb_regions: Stratify by variants in UKB regions.
    :param broad_regions: Stratify by variants in Broad regions.
    :param by_csqs: Stratify by variant consequence.
    :param rare_variants: Stratify by variants which have adj AF <0.1%.
    :param vep_canonical: Stratify by canonical transcripts. If trying to get
        variants in all transcripts, set it to False. Default is True.
    :param vep_mane: Stratify by MANE transcripts. Default is False.
    """
    vds = get_gnomad_v4_vds(test=test, release_only=True, split=True)

    # Read in release data and select relevant annotations.
    ht = release_sites(data_type=data_type).ht()
    ht = ht.select("freq", "filters", "region_flags", "vep")

    if test:
        logger.info("Test: filtering to variants on chr22")
        ht = ht.filter(ht.locus.contig == "chr22")

    logger.info("Filter VDS to variants in gnomAD release.")
    vds = hl.vds.filter_variants(vds, ht.select())

    # Create variant matrix table from VDS and annotate with all relevant info.
    vmt = vds.variant_data
    logger.info("Variant MT created and checkpointing, without annotations")
    vmt = vmt.checkpoint(
        new_temp_file(f'vmt_{"TEST" if test else ""}', extension="mt"),
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
            canonical=vep_canonical,
            mane_select=vep_mane,
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
    if ukb_regions and broad_regions:
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


def compute_stratified_agg_stats(
    ht: hl.Table,
    data_type: str = "exomes",
    qc_metrics: List[str] = ["n_het", "n_hom_var", "n_non_ref", "n_singleton"],
    strats: List[str] = ["all_variants"],
    by_ancestry: bool = False,
) -> hl.Table:
    """
    Compute aggregate statistics for stratified sample QC metrics.

    :param ht: Table containing sample QC metrics.
    :param data_type: String of either "exomes" or "genomes". Default is "exomes".
    :param qc_metrics: List of sample QC metrics to compute aggregate statistics for.
    :param strats: List of strata to compute aggregate statistics for.
    :param by_ancestry: Boolean indicating whether to stratify by ancestry.
    :return: Table containing aggregate statistics for stratified sample QC metrics.
    """

    def _get_metric_stats(ht, strat, metric):
        metric_expr = ht.aggregated_vmt_sample_qc[strat][metric]
        metric_mean = hl.agg.mean(metric_expr)
        metric_quantiles = hl.agg.approx_quantiles(
            metric_expr, [0.0, 0.25, 0.5, 0.75, 1.0]
        )

        return hl.struct(mean=metric_mean, quantiles=metric_quantiles)

    def _get_agg_expr(ht):
        return hl.struct(
            **{
                strat: hl.struct(
                    **{
                        _metric: _get_metric_stats(ht, strat, _metric)
                        for _metric in qc_metrics
                    }
                )
                for strat in strats
            }
        )

    ht = ht.annotate_globals(all_samples=ht.aggregate(_get_agg_expr(ht)))

    if by_ancestry:
        meta_ht = meta(data_type=data_type).ht()
        ht = ht.annotate(gen_anc=meta_ht[ht.s].population_inference.pop)
        by_ancestry_agg = ht.aggregate(hl.agg.group_by(ht.gen_anc, _get_agg_expr(ht)))
        # Remove 'None' key, otherwise it can't be annotated to globals
        by_ancestry_agg.pop(None)
        ht = ht.annotate_globals(**by_ancestry_agg)

    return ht


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
    pass_filters = args.pass_filters
    ukb_regions = args.filter_ukb_regions
    broad_regions = args.filter_broad_regions
    by_csqs = args.by_csqs
    rare_variants = args.rare_variants

    if args.create_per_sample_counts:
        per_sample_ht = create_per_sample_counts(
            data_type=data_type,
            test=test,
            pass_filters=pass_filters,
            ukb_regions=ukb_regions,
            broad_regions=broad_regions,
            by_csqs=by_csqs,
            rare_variants=rare_variants,
            vep_canonical=args.vep_canonical,
            vep_mane=args.vep_mane,
        )

        per_sample_ht.checkpoint(
            get_per_sample_counts(test=test, data_type=data_type).path,
            overwrite=overwrite,
        )

    if args.stratify_and_aggregate:
        per_sample_ht = get_per_sample_counts(test=test, data_type=data_type).ht()
        strats = []
        stratify_dictionary = {
            "all_variants": True,
            "pass_filters": pass_filters,
            "ukb_regions": ukb_regions,
            "broad_regions": broad_regions,
            "ukb_and_broad_regions": ukb_regions and broad_regions,
            "rare_variants": rare_variants,
            "lof": by_csqs,
            "missense": by_csqs,
            "synonymous": by_csqs,
        }

        for key, item in stratify_dictionary.items():
            strats.append(key)
        stratified_agg_ht = compute_stratified_agg_stats(
            per_sample_ht,
            data_type=data_type,
            strats=strats,
            by_ancestry=args.by_ancestry,
        )

        stratified_agg_ht.checkpoint(
            get_per_sample_counts(
                test=test,
                data_type=data_type,
                suffix="aggregated_",
            ).path,
            overwrite=overwrite,
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
        "--create-per-sample-counts",
        help="Create per-sample counts resource",
        action="store_true",
    )
    parser.add_argument(
        "--stratify-and-aggregate",
        help="Run code to aggregate and stratify (exising?) sample qc table",
        action="store_true",
    )
    parser.add_argument(
        "--pass-filters",
        help=(
            "Calculate per-sample variants for those variants that pass all"
            " variant qc filters."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--filter-ukb-regions",
        help="Stratify by variants in UKB capture regions",
        action="store_true",
    )
    parser.add_argument(
        "--filter-broad-regions",
        help="Stratify by variants in Broad capture regions",
        action="store_true",
    )
    parser.add_argument(
        "--rare-variants",
        help="Stratify by variants with adj AF <0.1%.",
        action="store_true",
    )
    parser.add_argument(
        "--by-csqs",
        help="Stratify by variants consequences.",
        action="store_true",
    )
    parser.add_argument(
        "--by-ancestry",
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
        "--vep-canonical",
        help="Canonical argument for filter_vep_transcript_csqs",
        action="store_true",
    )
    parser.add_argument(
        "--vep-mane",
        help="Mane_select argument for filter_vep_transcript_csqs",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
