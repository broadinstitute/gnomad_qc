"""
Script to get per-sample variant counts and aggregate sample statistics.

The following per-sample variant counts can be calculated:

    - Total number of variants
    - Number of variants that pass all variant qc filters
    - Number of variants in UK Biobank capture regions
    - Number of variants in Broad capture regions
    - Number of rare variants (adj AF <0.1%)
    - Number of loss-of-function variants
    - Number of missense variants
    - Number of synonymous variants
    - Number of singletons
    - Number of heterozygous variants
    - Number of homozygous variants
    - Number of non-reference variants
    - Number of variants in UK Biobank and Broad capture regions

The following aggregate sample stats of all of the above per-sample counts can be
computed:

        - Mean
        - Quantiles (0.0, 0.25, 0.5, 0.75, 1.0)

Aggregate statistics can also be computed by ancestry.
"""
# TODO: Maybe move to a folder called assessment and rename to
#  calculate_per_sample_stats.py, also add a resource file for assessments.
import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources import meta
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import get_per_sample_counts, release_sites
from gnomad_qc.v4.resources.temp_hail_methods import (
    vmt_sample_qc,
    vmt_sample_qc_variant_annotations,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("per_sample_stats")
logger.setLevel(logging.INFO)


def create_per_sample_counts_ht(
    mt: hl.MatrixTable,
    annotation_ht: hl.Table,
    pass_filters: bool = False,
    ukb_capture: bool = False,
    broad_capture: bool = False,
    by_csqs: bool = False,
    rare_variants: bool = False,
    vep_canonical: bool = True,
    vep_mane: bool = False,
) -> hl.Table:
    """
    Create Table of Hail's sample_qc output broken down by requested variant groupings.

    Useful for finding the number of variants (total or specific types) per sample.

    The following per-sample variant counts are calculated:

        - Total number of variants
        - Number of variants that pass all variant qc filters
        - Number of variants in UK Biobank capture regions
        - Number of variants in Broad capture regions
        - Number of rare variants (adj AF <0.1%)
        - Number of loss-of-function variants
        - Number of missense variants
        - Number of synonymous variants
        - Number of singletons
        - Number of heterozygous variants
        - Number of homozygous variants
        - Number of non-reference variants
        - Number of variants in UK Biobank and Broad capture regions

    Additionally, counts can be computed for the following variant filtering criteria:

        - Variants that pass all variant qc filters
        - Variants in UK Biobank capture regions
        - Variants in Broad capture regions
        - Variants in UK Biobank and Broad capture regions
        - Rare variants (adj AF <0.1%)
        - Loss-of-function variants
        - Variant consequences:

            - Loss-of-function
            - Missense
            - Synonymous

        - Variants in all transcripts
        - Variants in MANE transcripts
        - Variants in canonical transcripts

    :param mt: Input MatrixTable containing variant data. Must have multi-allelic sites
        split.
    :param annotation_ht: Table containing variant annotations. The following
        annotations are required: 'freq', 'filters', and 'region_flags'. If `by_csqs` is
        True, 'vep' is also required.
    :param pass_filters: Include count of variants that pass all variant QC filters.
    :param ukb_capture: Include count of variants that are in UKB capture intervals.
    :param broad_capture: Include count of variants that are in Broad capture intervals
    :param by_csqs: Include count of variants by variant consequence: loss-of-function,
        missense, and synonymous.
    :param rare_variants: Include count of rare variants, defined as those which have
        adj AF <0.1%.
    :param vep_canonical: If `by_csqs` is True, filter to only canonical transcripts.
        If trying count variants in all transcripts, set it to False. Default is True.
    :param vep_mane: If `by_csqs` is True, filter to only MANE transcripts. Default is
        False.
    :return: Table containing per-sample variant counts.
    """
    logger.info("Filtering input MT to variants in the supplied annotation HT...")
    mt = mt.semi_join_rows(annotation_ht)

    # Add extra Allele Count and Allele Type annotations to variant MatrixTable,
    # according to Hail standards, to help their computation.
    mt = mt.annotate_rows(
        **vmt_sample_qc_variant_annotations(global_gt=mt.GT, alleles=mt.alleles)
    )

    keep_annotations = ["freq", "filters", "region_flags"]
    if by_csqs:
        annotation_ht = filter_vep_transcript_csqs(
            annotation_ht,
            synonymous=False,
            canonical=vep_canonical,
            mane_select=vep_mane,
        )
        annotation_ht = get_most_severe_consequence_for_summary(annotation_ht)
        keep_annotations.extend(["most_severe_csq", "lof"])

    # Annotate the MT with the needed annotations.
    annotation_ht = annotation_ht.select(*keep_annotations).checkpoint(
        hl.utils.new_temp_file("annotation_ht", "ht")
    )
    mt = mt.annotate_rows(**annotation_ht[mt.row_key])

    filter_expr = {"all_variants": True}

    if pass_filters:
        filter_expr["pass_filters"] = hl.len(mt.filters) == 0
    if ukb_capture:
        filter_expr["ukb_capture"] = ~mt.region_flags.outside_ukb_capture_region
    if broad_capture:
        filter_expr["broad_capture"] = ~mt.region_flags.outside_broad_capture_region
    if ukb_capture and broad_capture:
        filter_expr["ukb_broad_capture_intersect"] = (
            filter_expr["ukb_capture"] & filter_expr["broad_capture"]
        )
        filter_expr["ukb_broad_capture_union"] = (
            filter_expr["ukb_capture"] | filter_expr["broad_capture"]
        )
    if rare_variants:
        # TODO: Maybe make this a parameter?
        filter_expr["rare_variants"] = mt.freq[0].AF < 0.001
    if by_csqs:
        filter_expr["lof"] = ~hl.is_missing(mt.lof)
        filter_expr["missense"] = mt.most_severe_csq == "missense_variant"
        filter_expr["synonymous"] = mt.most_severe_csq == "synonymous_variant"

    # Run Hail's 'vmt_sample_qc' for all requested filter groups.
    ht = mt.select_cols(
        _sample_qc=hl.struct(
            **{
                ann: hl.agg.filter(
                    expr,
                    vmt_sample_qc(
                        global_gt=mt.GT,
                        gq=mt.GQ,
                        variant_ac=mt.variant_ac,
                        variant_atypes=mt.variant_atypes,
                        dp=mt.DP,
                    ),
                )
                for ann, expr in filter_expr.items()
            }
        )
    ).cols()

    return ht.select(**ht._sample_qc)


def compute_agg_sample_stats(
    ht: hl.Table,
    meta_ht: Optional[hl.Table] = None,
    by_ancestry: bool = False,
) -> hl.Struct:
    """
    Compute aggregate statistics for per-sample QC metrics.

    :param ht: Table containing sample QC metrics.
    :param meta_ht: Optional Table containing sample metadata. Required if
        `by_ancestry` is True.
    :param by_ancestry: Boolean indicating whether to stratify by ancestry.
    :return: Struct of aggregate statistics for per-sample QC metrics.
    """
    if by_ancestry and meta_ht is None:
        raise ValueError(
            "If `by_ancestry` is True, a Table containing sample metadata is required."
        )

    if by_ancestry:
        ht = ht.annotate(gen_anc=meta_ht[ht.s].population_inference.pop)

    agg_expr = {
        strat: hl.struct(
            **{
                metric: hl.struct(
                    mean=hl.agg.mean(ht[strat][metric]),
                    quantiles=hl.agg.approx_quantiles(
                        ht[strat][metric], [0.0, 0.25, 0.5, 0.75, 1.0]
                    ),
                )
                for metric in ht[strat]
                if isinstance(ht[strat][metric], hl.expr.NumericExpression)
            }
        )
        for strat in ht.row_value
        if isinstance(ht[strat], hl.expr.StructExpression)
    }

    if by_ancestry:
        agg_expr = hl.struct(
            all_samples=agg_expr,
            # Remove missing gen_anc from stratified stats.
            stratified_by_gen_anc=hl.agg.filter(
                hl.is_defined(ht.gen_anc), hl.agg.group_by(ht.gen_anc, agg_expr)
            ),
        )

    return ht.aggregate(agg_expr)


def main(args):
    """Collect per-sample variant counts and aggregate sample statistics."""
    hl.init(
        log="/per_sample_stats",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    data_type = args.data_type
    test = args.test
    overwrite = args.overwrite

    if data_type == "exomes":
        logger.info("Calculating per-sample variant statistics for exomes...")
        mt = get_gnomad_v4_vds(
            test=test, release_only=True, split=True, chrom="chr22" if test else None
        ).variant_data
    else:
        logger.info("Calculating per-sample variant statistics for genomes...")
        mt = get_gnomad_v4_genomes_vds(
            test=test, release_only=True, split=True, chrom="chr22" if test else None
        ).variant_data

    ht = create_per_sample_counts_ht(
        mt,
        release_sites(data_type=data_type).ht(),
        pass_filters=args.pass_filters,
        ukb_capture=args.filter_ukb_capture_intervals,
        broad_capture=args.filter_broad_capture_intervals,
        by_csqs=args.by_csqs,
        rare_variants=args.rare_variants,
        vep_canonical=args.vep_canonical,
        vep_mane=args.vep_mane,
    ).checkpoint(hl.utils.new_temp_file("per_sample_counts", "ht"))

    if args.add_aggregate_sample_stats_global:
        logger.info("Computing aggregate sample statistics...")
        sample_qc_agg_stats = compute_agg_sample_stats(
            ht,
            meta_ht=meta(data_type=data_type).ht(),
            by_ancestry=args.by_ancestry,
        )
        logger.info("Aggregate sample statistics: %s", sample_qc_agg_stats)
        ht = ht.annotate_globals(sample_qc_agg_stats=sample_qc_agg_stats)

    ht.write(
        get_per_sample_counts(
            test=test, data_type=data_type, suffix=args.custom_suffix
        ).path,
        overwrite=overwrite,
    )


if __name__ == "__main__":
    # TODO: Could make some of these arguments --skip-<filter> instead of --<filter> so
    #  by default they all run.
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
        "--add-aggregate-sample-stats-global",
        help=(
            "Compute aggregate sample statistics from the per-sample counts and add "
            "them as globals to the output Table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--pass-filters",
        help=(
            "Calculate per-sample variants for those variants that pass all"
            " variant QC filters."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--filter-ukb-capture-intervals",
        help="Include counts of variants filtered to UK Biobank capture regions.",
        action="store_true",
    )
    parser.add_argument(
        "--filter-broad-capture-intervals",
        help="Include counts of variants filtered to Broad capture regions.",
        action="store_true",
    )
    parser.add_argument(
        "--rare-variants",
        help="Include counts of rare variants (adj AF <0.1%).",
        action="store_true",
    )
    parser.add_argument(
        "--by-csqs",
        help=(
            "Include counts of variants by consequence type: loss-of-function, "
            "missense, and synonymous."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--vep-canonical",
        help="Whether to filter to only canonical transcripts. when using --by-csqs.",
        action="store_true",
    )
    parser.add_argument(
        "--vep-mane",
        help="Whether to filter to only MANE transcripts. when using --by-csqs.",
        action="store_true",
    )
    parser.add_argument(
        "--by-ancestry",
        help=(
            "Output statistics for number of singletons, n_het, and n_hom by inferred "
            "genetic ancestry group. Only relevant when using "
            "--compute-aggregate-sample-stats"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--custom-suffix",
        type=str,
        default=None,
        help="Custom string to append to output names.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
