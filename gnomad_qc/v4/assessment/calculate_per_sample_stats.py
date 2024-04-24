"""
Script to get per-sample variant counts and aggregate sample statistics.

The following per-sample variant counts (including heterozygous, homozygous, non-ref, singletons etc.) can be calculated:

    - Total number of variants
    - Number of variants that pass all variant qc filters
    - Number of variants in UK Biobank capture regions
    - Number of variants in Broad capture regions
    - Number of variants in the intersect of UK Biobank and Broad capture regions
    - Number of variants in the union of UK Biobank and Broad capture regions
    - Number of rare variants (adj AF <0.1%)
    - Number of loss-of-function variants
    - Number of missense variants
    - Number of synonymous variants

The following aggregate sample stats of all of the above per-sample counts can be
computed:

        - Mean
        - Quantiles (0.0, 0.25, 0.5, 0.75, 1.0)

Aggregated statistics can also be computed by ancestry.
"""
import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    CSQ_CODING,
    CSQ_NON_CODING,
    LOF_CSQ_SET,
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)
from hail.vds.sample_qc import vmt_sample_qc, vmt_sample_qc_variant_annotations

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources import meta
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import get_per_sample_counts, release_sites

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
    rare_variants_af: float = 0.001,
) -> hl.Table:
    """
    Create Table of Hail's sample_qc output broken down by requested variant groupings.

    Useful for finding the number of variants per sample, either all variants, or
    variants fall into specific capture regions, or variants that are rare
    (adj AF <0.1%), or variants categorized by predicted consequences in all, canonical
    or mane transcripts.

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
    :param rare_variants_af: The allele frequency threshold to use for rare variants.
        Default is 0.001.
    :return: Table containing per-sample variant counts.
    """
    logger.info("Filtering out low confidence regions...")
    mt = filter_low_conf_regions(mt, filter_decoy=False)

    logger.info("Filtering input MT to variants in the supplied annotation HT...")
    mt = mt.semi_join_rows(annotation_ht)

    # Add extra Allele Count and Allele Type annotations to variant MatrixTable,
    # according to Hail standards, to help their computation.
    variant_ac, variant_types = vmt_sample_qc_variant_annotations(
        global_gt=mt.GT, alleles=mt.alleles
    )

    mt = mt.annotate_rows(
        variant_ac=variant_ac,
        variant_atypes=variant_types,
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
        keep_annotations.extend(["most_severe_csq", "lof", "no_lof_flags"])

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
    if all(["ukb_capture", "broad_capture", "pass_filters"]):
        filter_expr["ukb_broad_capture_union_pass_filters"] = (
            filter_expr["ukb_broad_capture_union"] & filter_expr["pass_filters"]
        )
    if rare_variants:
        filter_expr["rare_variants"] = mt.freq[0].AF < rare_variants_af
    if by_csqs:
        filter_expr["coding"] = filter_expr["pass_filters"] & CSQ_CODING.contains(
            mt.most_severe_csq
        )
        filter_expr["non_coding"] = filter_expr[
            "pass_filters"
        ] & CSQ_NON_CODING.contains(mt.most_severe_csq)
        filter_expr["lof"] = (
            filter_expr["pass_filters"]
            & hl.if_else(
                hl.any(lambda csq: mt.most_severe_csq == csq, LOF_CSQ_SET), True, False
            )
            & (mt.lof == "HC")
            & mt.no_lof_flags
        )
        filter_expr["missense"] = filter_expr["pass_filters"] & (
            mt.most_severe_csq == "missense_variant"
        )
        filter_expr["synonymous"] = filter_expr["pass_filters"] & (
            mt.most_severe_csq == "synonymous_variant"
        )
        filter_expr["intronic"] = filter_expr["pass_filters"] & (
            mt.most_severe_csq == "intron_variant"
        )
        filter_expr["intergenic"] = filter_expr["pass_filters"] & (
            mt.most_severe_csq == "intergenic_variant"
        )

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
    by_subset: bool = False,
) -> hl.Table:
    """
    Compute aggregate statistics for per-sample QC metrics.

    :param ht: Table containing sample QC metrics.
    :param meta_ht: Optional Table containing sample metadata. Required if
        `by_ancestry` is True.
    :param by_ancestry: Boolean indicating whether to stratify by ancestry.
    :param by_subset: Boolean indicating whether to stratify by subset. This is only
         working on "exomes" data.
    :return: Struct of aggregate statistics for per-sample QC metrics.
    """
    if by_ancestry and meta_ht is None:
        raise ValueError(
            "If `by_ancestry` is True, a Table containing sample metadata is required."
        )

    subset = (
        ["gnomad"]
        if not by_subset
        else [
            "gnomad",
            hl.if_else(meta_ht[ht.s].project_meta.ukb_sample, "ukb", "non-ukb"),
        ]
    )
    gen_anc = (
        ["global"]
        if not by_ancestry
        else ["global", meta_ht[ht.s].population_inference.pop]
    )

    all_strats = [
        strat
        for strat in ht.row_value
        if isinstance(ht[strat], hl.expr.StructExpression)
    ]

    ht = ht.transmute(
        subset=subset,
        gen_anc=gen_anc,
        stats_array=[(strat, ht[strat]) for strat in all_strats],
    )

    ht = ht.explode("_stats_array").explode("gen_anc").explode("subset")

    ht = ht.group_by("subset", "gen_anc", variant_filter=ht.stats_array[0]).aggregate(
        **{
            metric: hl.struct(
                mean=hl.agg.mean(ht.stats_array[1][metric]),
                quantiles=hl.agg.approx_quantiles(
                    ht.stats_array[1][metric], [0.0, 0.25, 0.5, 0.75, 1.0]
                ),
            )
            for metric in ht.stats_array[1]
            if isinstance(ht.stats_array[1][metric], hl.expr.NumericExpression)
        }
    )

    return ht


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
    per_sample_res = get_per_sample_counts(
        test=test, data_type=data_type, suffix=args.custom_suffix
    )
    per_sample_agg_res = get_per_sample_counts(
        test=test,
        data_type=data_type,
        suffix=args.custom_suffix,
        aggregated=True,
        by_ancestry=args.by_ancestry,
        by_subset=args.by_subset,
    )
    if data_type != "exomes" and args.by_subset:
        raise ValueError("Stratifying by subset is only working on exomes data type.")

    if args.create_per_sample_counts_ht:
        chrom = "chr22" if test else None
        if data_type == "exomes":
            logger.info("Calculating per-sample variant statistics for exomes...")
            mt = get_gnomad_v4_vds(
                test=test, release_only=True, split=True, chrom=chrom
            ).variant_data
        else:
            logger.info("Calculating per-sample variant statistics for genomes...")
            mt = get_gnomad_v4_genomes_vds(
                test=test, release_only=True, split=True, chrom=chrom
            ).variant_data

        release_ht = release_sites(data_type=data_type).ht()
        if test:
            release_ht = hl.filter_intervals(
                release_ht, [hl.parse_locus_interval("chr22")]
            )

        create_per_sample_counts_ht(
            mt,
            release_ht,
            pass_filters=not args.skip_pass_filters,
            ukb_capture=not args.skip_filter_ukb_capture_intervals,
            broad_capture=not args.skip_filter_broad_capture_intervals,
            by_csqs=not args.skip_by_csqs,
            rare_variants=not args.skip_rare_variants,
            vep_canonical=args.vep_canonical,
            vep_mane=args.vep_mane,
            rare_variants_af=args.rare_variants_af,
        ).write(per_sample_res.path, overwrite=overwrite)

    if args.aggregate_sample_stats:
        logger.info("Computing aggregate sample statistics...")
        ht = per_sample_res.ht().checkpoint(
            hl.utils.new_temp_file("per_sample_counts", "ht")
        )
        ht = compute_agg_sample_stats(
            ht,
            meta(data_type=data_type).ht(),
            by_ancestry=args.by_ancestry,
            by_subset=args.by_subset,
        )
        ht.write(per_sample_agg_res.path, overwrite=overwrite)


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
        "--create-per-sample-counts-ht",
        help="Create per-sample variant counts Table.",
        action="store_true",
    )
    parser.add_argument(
        "--aggregate-sample-stats",
        help=(
            "Compute aggregate sample statistics from the per-sample counts and add "
            "them as globals to the output Table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--skip-pass-filters",
        help=(
            "Whether to skip calculating the number of per-sample variants for those "
            "variants that pass all variant QC filters."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--skip-filter-ukb-capture-intervals",
        help=(
            "Whether to skip calculating the number of variants filtered to UK Biobank "
            "capture regions."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--skip-filter-broad-capture-intervals",
        help=(
            "Whether to skip calculating the number of variants filtered to Broad "
            "capture regions."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--skip-rare-variants",
        help="Whether to skip calculating the number of rare variants (adj AF <0.1%).",
        action="store_true",
    )
    parser.add_argument(
        "--rare-variants-af",
        type=float,
        default=0.001,
        help="The allele frequency threshold to use for rare variants.",
    )
    parser.add_argument(
        "--skip-by-csqs",
        help=(
            "Whether to skip calculating the number of variants by consequence type:"
            " missense, synonymous, and loss-of-function (splice_acceptor_variant,"
            " splice_donor_variant, stop_gained, frameshift_variant)."
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
        "--by-subset",
        help=(
            "Get aggregate statistics for the whole dataset and for ukb and non-ukb "
            "subsets. This is only working on exomes data type."
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
