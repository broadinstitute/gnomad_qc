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
from typing import Dict, Optional

import hail as hl
from gnomad.assessment.summary_stats import (
    get_summary_stats_csq_filter_expr,
    get_summary_stats_variant_filter_expr,
)
from gnomad.utils.annotations import annotate_with_ht
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    CSQ_CODING,
    CSQ_NON_CODING,
    LOF_CSQ_SET,
    LOFTEE_LABELS,
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)
from hail.vds.sample_qc import vmt_sample_qc, vmt_sample_qc_variant_annotations

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources import meta
from gnomad_qc.v4.resources.basics import (
    get_gnomad_v4_genomes_vds,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.release import get_per_sample_counts, release_sites

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("per_sample_stats")
logger.setLevel(logging.INFO)


def get_capture_filter_exprs(
    ht: hl.Table,
    pass_filter_expr: Optional[hl.expr.BooleanExpression] = None,
    ukb_capture: bool = False,
    broad_capture: bool = False,
) -> Dict[str, hl.expr.BooleanExpression]:
    """
    Get filter expressions for UK Biobank and Broad capture regions.

    :param ht: Table containing variant annotations. The following annotations are
        required: 'region_flags'.
    :param pass_filter_expr: Expression for variants that pass all variant QC filters.
    :param ukb_capture: Include count of variants that are in UKB capture intervals.
    :param broad_capture: Include count of variants that are in Broad capture intervals.
    :return: Dictionary of filter expressions for UK Biobank and Broad capture regions.
    """
    filter_expr = {}
    log_list = []
    if ukb_capture:
        log_list.append("variants in UK Biobank capture regions")
        filter_expr["ukb_capture"] = ~ht.region_flags.outside_ukb_capture_region
    if broad_capture:
        log_list.append("variants in Broad capture regions")
        filter_expr["broad_capture"] = ~ht.region_flags.outside_broad_capture_region
    if ukb_capture and broad_capture:
        log_list.append("variants in the intersect of UKB and Broad capture regions")
        filter_expr["ukb_broad_capture_intersect"] = (
            filter_expr["ukb_capture"] & filter_expr["broad_capture"]
        )
        log_list.append("variants in the union of UKB and Broad capture regions")
        filter_expr["ukb_broad_capture_union"] = (
            filter_expr["ukb_capture"] | filter_expr["broad_capture"]
        )
        if pass_filter_expr is not None:
            log_list.append(
                "PASS variants in the intersect of UKB and Broad capture regions"
            )
            filter_expr["ukb_broad_capture_union_pass_filters"] = (
                filter_expr["ukb_broad_capture_union"] & pass_filter_expr
            )

    logger.info("Adding filtering for:\n\t%s...", "\n\t".join(log_list))

    return filter_expr


def get_summary_stats_filter_groups_ht(
    ht: hl.Table,
    pass_filters: bool = False,
    ukb_capture: bool = False,
    broad_capture: bool = False,
    by_csqs: bool = False,
    vep_canonical: bool = True,
    vep_mane: bool = False,
    rare_variants_afs: Optional[list[float]] = None,
) -> hl.Table:
    """
    Create Table of filter groups for summary stats.

    :param ht: Table containing variant annotations. The following annotations are
        required: 'freq', 'filters', and 'region_flags'. If `by_csqs` is True, 'vep' is
        also required.
    :param pass_filters: Include count of variants that pass all variant QC filters.
    :param ukb_capture: Include count of variants that are in UKB capture intervals.
    :param broad_capture: Include count of variants that are in Broad capture intervals
    :param by_csqs: Include count of variants by variant consequence: loss-of-function,
        missense, and synonymous.
    :param vep_canonical: If `by_csqs` is True, filter to only canonical transcripts.
        If trying count variants in all transcripts, set it to False. Default is True.
    :param vep_mane: If `by_csqs` is True, filter to only MANE transcripts. Default is
        False.
    :param rare_variants_afs: The allele frequency thresholds to use for rare variants.
    :return: Table containing filter groups for summary stats.
    """
    # Filter to only canonical or MANE transcripts if requested and get the most severe
    # consequence for each variant.
    if by_csqs:
        ht = filter_vep_transcript_csqs(
            ht,
            synonymous=False,
            canonical=vep_canonical,
            mane_select=vep_mane,
        )
        ht = get_most_severe_consequence_for_summary(ht)

    # Create filter expressions for the requested variant groupings.
    filter_exprs = get_summary_stats_variant_filter_expr(
        ht,
        filter_lcr=True,
        filter_expr=ht.filters if pass_filters else None,
        freq_expr=ht.freq[0].AF,
        max_af=rare_variants_afs,
    )

    # Create filter expressions for the requested variant groupings.
    pass_expr = filter_exprs.get("pass_filters")
    filter_exprs.update(
        get_capture_filter_exprs(ht, pass_expr, ukb_capture, broad_capture)
    )

    # Create filter expressions for the requested consequence types.
    if by_csqs:
        additional_csqs = [
            "missense_variant",
            "synonymous_variant",
            "intron_variant",
            "intergenic_variant",
        ]
        csq_filter_expr = get_summary_stats_csq_filter_expr(
            ht,
            lof_csq_set=LOF_CSQ_SET,
            lof_label_set=LOFTEE_LABELS,
            lof_no_flags=True,
            lof_any_flags=True,
            lof_loftee_combinations=True,
            additional_csq_sets={
                "coding": set(CSQ_CODING),
                "non_coding": set(CSQ_NON_CODING),
            },
            additional_csqs=additional_csqs,
        )

        # For all csq breakdowns, filter to only variants that pass all filters.
        csq_filter_expr = {k: v & pass_expr for k, v in csq_filter_expr.items()}
        filter_exprs.update(csq_filter_expr)

    # Remove 'pass_filters' filter expression if not requested as a summary stats group.
    if not pass_filters:
        filter_exprs.pop("pass_filters")

    no_lcr_expr = filter_exprs.pop("no_lcr")
    filter_groups = list(filter_exprs.keys())
    ht = ht.select(
        _no_lcr=no_lcr_expr, filter_groups=[filter_exprs[g] for g in filter_groups]
    )
    ht = ht.annotate_globals(filter_groups_meta=filter_groups)
    logger.info("Filter groups for summary stats: %s", filter_groups)

    # Filter to only variants that are not in low confidence regions.
    ht = ht.filter(ht._no_lcr).drop("_no_lcr")
    ht = ht.checkpoint(hl.utils.new_temp_file("stats_annotation", "ht"))

    return ht


def create_per_sample_counts_ht(
    mt: hl.MatrixTable, filter_group_ht: hl.Table
) -> hl.Table:
    """
    Create Table of Hail's sample_qc output broken down by requested variant groupings.

    Useful for finding the number of variants per sample, either all variants, or
    variants fall into specific capture regions, or variants that are rare
    (adj AF <0.1%), or variants categorized by predicted consequences in all, canonical
    or mane transcripts.

    :param mt: Input MatrixTable containing variant data. Must have multi-allelic sites
        split.
    :param filter_group_ht: Table containing filter groups for summary stats.
    :return: Table containing per-sample variant counts.
    """
    # Add extra Allele Count and Allele Type annotations to variant MatrixTable,
    # according to Hail standards, to help their computation.
    variant_ac, variant_types = vmt_sample_qc_variant_annotations(
        global_gt=mt.GT, alleles=mt.alleles
    )
    mt = mt.annotate_rows(variant_ac=variant_ac, variant_atypes=variant_types)

    # Annotate the MT with the needed annotations.
    mt = annotate_with_ht(mt, filter_group_ht, filter_missing=True)

    # Run Hail's 'vmt_sample_qc' for all requested filter groups.
    qc_expr = vmt_sample_qc(
        global_gt=mt.GT,
        gq=mt.GQ,
        variant_ac=mt.variant_ac,
        variant_atypes=mt.variant_atypes,
        dp=mt.DP,
    )
    filter_groups = hl.eval(filter_group_ht.filter_groups_meta)
    ht = mt.select_cols(
        _qc=hl.agg.array_agg(lambda f: hl.agg.filter(f, qc_expr), mt.filter_groups)
    ).cols()
    ht = ht.checkpoint(hl.utils.new_temp_file("per_sample_counts", "ht"))

    # Add 'n_indel' to the output Table.
    qc_expr = ht._qc.map(lambda x: x.annotate(n_indel=x.n_insertion + x.n_deletion))

    # Convert each element of the summary stats array to a row annotation on the Table.
    ht = ht.select(**{f: qc_expr[i] for i, f in enumerate(filter_groups)})

    return ht


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
    if meta_ht is None and by_ancestry:
        raise ValueError("If `by_ancestry` is True, `meta_ht` is required.")
    if meta_ht is None and by_subset:
        raise ValueError("If `by_subset` is True, `meta_ht` is required.")

    subset = ["gnomad"]
    gen_anc = ["global"]
    if meta_ht is not None:
        meta_s = meta_ht[ht.s]
        subset_expr = hl.if_else(meta_s.project_meta.ukb_sample, "ukb", "non-ukb")
        subset += [subset_expr] if by_subset else []
        gen_anc += [meta_s.population_inference.pop] if by_ancestry else []

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

    ht = ht.explode("stats_array").explode("gen_anc").explode("subset")

    ht = ht.group_by("subset", "gen_anc", variant_filter=ht.stats_array[0]).aggregate(
        **{
            metric: hl.struct(
                mean=hl.agg.mean(ht.stats_array[1][metric]),
                min=hl.agg.min(ht.stats_array[1][metric]),
                max=hl.agg.max(ht.stats_array[1][metric]),
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
    ukb_capture_intervals = not args.skip_filter_ukb_capture_intervals
    broad_capture_intervals = not args.skip_filter_broad_capture_intervals
    rare_variants_afs = args.rare_variants_afs if not args.skip_rare_variants else None
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

    try:
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
                ukb_capture_intervals = False
                broad_capture_intervals = False

            release_ht = release_sites(data_type=data_type).ht()
            if test:
                release_ht = hl.filter_intervals(
                    release_ht, [hl.parse_locus_interval("chr22")]
                )

            filter_groups_ht = get_summary_stats_filter_groups_ht(
                release_ht,
                pass_filters=not args.skip_pass_filters,
                ukb_capture=ukb_capture_intervals,
                broad_capture=broad_capture_intervals,
                by_csqs=not args.skip_by_csqs,
                vep_canonical=args.vep_canonical,
                vep_mane=args.vep_mane,
                rare_variants_afs=rare_variants_afs,
            )
            create_per_sample_counts_ht(mt, filter_groups_ht).write(
                per_sample_res.path, overwrite=overwrite
            )

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
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("per_sample_stats"))


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
        "--rare-variants-afs",
        type=float,
        default=[0.0001, 0.001, 0.01],
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
