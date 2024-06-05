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
from copy import deepcopy
from typing import Dict, Optional

import hail as hl
from gnomad.assessment.summary_stats import (
    get_summary_stats_csq_filter_expr,
    get_summary_stats_filter_group_meta,
    get_summary_stats_variant_filter_expr,
)
from gnomad.utils.annotations import annotate_with_ht
from gnomad.utils.filtering import filter_to_adj, filter_to_autosomes
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    CSQ_CODING,
    CSQ_NON_CODING,
    LOF_CSQ_SET,
    LOFTEE_LABELS,
    filter_vep_transcript_csqs,
    get_most_severe_consequence_for_summary,
)
from hail.genetics.allele_type import AlleleType
from hail.vds.sample_qc import vmt_sample_qc, vmt_sample_qc_variant_annotations

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources import meta
from gnomad_qc.v4.resources.assessment import (
    get_per_sample_counts,
    get_summary_stats_filtering_groups,
)
from gnomad_qc.v4.resources.basics import (
    get_gnomad_v4_genomes_vds,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("per_sample_stats")
logger.setLevel(logging.INFO)

# V4 has no OS in LOFTEE_LABELS.
LOFTEE_LABELS = deepcopy(LOFTEE_LABELS)
LOFTEE_LABELS.remove("OS")

SUM_STAT_FILTERS = {
    "variant_qc": ["none", "pass"],  # Quality control status of the variant.
    "capture": [  # Capture methods used.
        "ukb",
        "broad",
        "ukb_broad_intersect",
        "ukb_broad_union",
    ],
    "max_af": [0.0001, 0.001, 0.01],  # Maximum allele frequency thresholds.
    "csq_set": ["lof", "coding", "non_coding"],  # Consequence sets.
    "lof_csq": deepcopy(LOF_CSQ_SET),  # Loss-of-function consequence set.
    "csq": [  # Additional consequence types.
        "missense_variant",
        "synonymous_variant",
        "intron_variant",
        "intergenic_variant",
    ],
}
"""
Dictionary of default filter settings for summary stats.
"""

COMMON_FILTERS = {"variant_qc": ["pass"], "capture": ["ukb_broad_intersect"]}
"""
Dictionary of common filter settings to use for most summary stats.
"""

COMMON_FILTER_COMBOS = [["variant_qc"], ["variant_qc", "capture"]]
"""
List of common variant filter combinations to use for summary stats.
"""

LOF_FILTERS_FOR_COMBO = {
    "lof_csq_set": ["lof"],  # Loss-of-function consequence set.
    "loftee_label": deepcopy(LOFTEE_LABELS),  # LOFTEE loss-of-function labels.
    "loftee_HC": ["HC"],  # High-confidence LOFTEE label.
    "loftee_flags": ["no_flags", "with_flags"],  # High-confidence LOFTEE flag options.
}
"""
Dictionary of an additional filter group to use for loss-of-function filter
combinations.
"""

LOF_FILTER_COMBOS = [
    ["lof_csq", "loftee_label"],
    ["lof_csq_set", "loftee_label"],
    ["lof_csq_set", "loftee_HC", "loftee_flags"],
    ["lof_csq", "loftee_HC", "loftee_flags"],
]
"""
List of loss-of-function consequence combinations to use for summary stats.
"""

MAP_FILTER_FIELD_TO_META = {
    "lof_csq": "csq",
    "loftee_HC": "loftee_label",
    "lof_csq_set": "csq_set",
}
"""
Dictionary to rename keys in `SUM_STAT_FILTERS`, `COMMON_FILTERS`, or
`LOF_FILTERS_FOR_COMBO` to final metadata keys.
"""

ALLELE_TYPE_MAP = {
    "insertion": AlleleType.INSERTION,
    "deletion": AlleleType.DELETION,
    "transition": AlleleType.TRANSITION,
    "transversion": AlleleType.TRANSVERSION,
    "star": AlleleType.STAR,
}


def get_capture_filter_exprs(
    ht: hl.Table,
    ukb_capture: bool = False,
    broad_capture: bool = False,
) -> Dict[str, hl.expr.BooleanExpression]:
    """
    Get filter expressions for UK Biobank and Broad capture regions.

    :param ht: Table containing variant annotations. The following annotations are
        required: 'region_flags'.
    :param ukb_capture: Expression for variants that are in UKB capture intervals.
    :param broad_capture: Expression for variants that are in Broad capture intervals.
    :return: Dictionary of filter expressions for UK Biobank and Broad capture regions.
    """
    filter_expr = {}
    log_list = []
    if ukb_capture:
        log_list.append("variants in UK Biobank capture regions")
        filter_expr["capture_ukb"] = ~ht.region_flags.outside_ukb_capture_region

    if broad_capture:
        log_list.append("variants in Broad capture regions")
        filter_expr["capture_broad"] = ~ht.region_flags.outside_broad_capture_region

    if ukb_capture and broad_capture:
        log_list.append("variants in the intersect of UKB and Broad capture regions")
        filter_expr["capture_ukb_broad_intersect"] = (
            filter_expr["capture_ukb"] & filter_expr["capture_broad"]
        )

        log_list.append("variants in the union of UKB and Broad capture regions")
        filter_expr["capture_ukb_broad_union"] = (
            filter_expr["capture_ukb"] | filter_expr["capture_broad"]
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
    Create Table annotated with an array of booleans indicating whether a variant belongs to certain filter groups.

    A 'filter_groups' annotation is added to the Table containing an ArrayExpression of
    BooleanExpressions for each requested filter group.

    A 'filter_group_meta' global annotation is added to the Table containing an array
    of dictionaries detailing the filters used in each filter group.

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
    :return: Table containing an ArrayExpression of filter groups for summary stats.
    """
    csq_filter_expr = {}
    if by_csqs:
        # Filter to only canonical or MANE transcripts if requested and get the most
        # severe consequence for each variant.
        ht = filter_vep_transcript_csqs(
            ht,
            synonymous=False,
            canonical=vep_canonical,
            mane_select=vep_mane,
            filter_empty_csq=False,
        )
        ht = get_most_severe_consequence_for_summary(ht)

        # Create filter expressions for the requested consequence types.
        csq_filter_expr.update(
            get_summary_stats_csq_filter_expr(
                ht,
                lof_csq_set=LOF_CSQ_SET,
                lof_label_set=LOFTEE_LABELS,
                lof_no_flags=True,
                lof_any_flags=True,
                additional_csq_sets={
                    "coding": set(CSQ_CODING),
                    "non_coding": set(CSQ_NON_CODING),
                },
                additional_csqs=set(SUM_STAT_FILTERS["csq"]),
            )
        )

    # Create filter expressions for LCR, variant QC filters, and rare variant AFs if
    # requested.
    filter_exprs = {
        "all_variants": hl.literal(True),
        **get_capture_filter_exprs(ht, ukb_capture, broad_capture),
        **get_summary_stats_variant_filter_expr(
            ht,
            filter_lcr=True,
            filter_expr=ht.filters if pass_filters else None,
            freq_expr=ht.freq[0].AF,
            max_af=rare_variants_afs,
        ),
        **csq_filter_expr,
    }

    # Create the metadata for all requested filter groups.
    ss_filters = deepcopy(SUM_STAT_FILTERS)
    ss_filters["max_af"] = rare_variants_afs
    filter_group_meta = get_summary_stats_filter_group_meta(
        ss_filters,
        common_filter_combos=COMMON_FILTER_COMBOS,
        common_filter_override=COMMON_FILTERS,
        lof_filter_combos=LOF_FILTER_COMBOS,
        lof_filter_override=LOF_FILTERS_FOR_COMBO,
        filter_key_rename=MAP_FILTER_FIELD_TO_META,
    )

    # Create filter expressions for each filter group by combining the filter
    # expressions for each filter in the filter group metadata.
    filter_groups_expr = []
    final_meta = []
    for filter_group in filter_group_meta:
        # Initialize filter expression for the filter group with True to allow for
        # a filtering group that has no filters, e.g. all variants.
        filter_expr = hl.literal(True)
        filter_group_requested = True
        for k, v in filter_group.items():
            # Rename "loftee_flags" to "loftee" to match the filter expression keys.
            k = k.replace("loftee_flags", "loftee")
            # Determine the correct key for filter_expr, it can be a combination of
            # the key and value, or just the key followed by using the value to get the
            # filter expression from a struct.
            f_expr = filter_exprs.get(f"{k}_{v}")
            f_struct = filter_exprs.get(k)
            # If the filter group is in the combinations, but not filter_exprs, then
            # the filter group was not in the requested list.
            if f_expr is None and f_struct is None:
                filter_group_requested = False
                break
            filter_expr &= f_struct[v] if f_expr is None else f_expr

        if filter_group_requested:
            filter_groups_expr.append(filter_expr)
            final_meta.append(filter_group)
        else:
            logger.warning(
                "Filter group %s was not requested and will not be included in the "
                "summary stats.",
                filter_group,
            )

    # Remove 'no_lcr' filter expression from filter groups and annotate the Table with
    # the no_lcr filter and an array of the filter groups.
    ht = ht.select(
        _no_lcr=filter_exprs["no_lcr"],
        filter_groups=filter_groups_expr,
        variant_ac=[hl.missing(hl.tint32), ht.freq[0].AC],
        variant_af=ht.freq[0].AF,
    )
    ht = ht.select_globals(filter_group_meta=final_meta)
    logger.info("Filter groups for summary stats: %s", filter_group_meta)

    # Filter to only variants that are not in low confidence regions.
    return ht.filter(ht._no_lcr).drop("_no_lcr")


def create_per_sample_counts_ht(
    mt: hl.MatrixTable, filter_group_ht: hl.Table, autosomes_only: bool = False
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
    :param autosomes_only: Whether to restrict analysis to autosomes only. Default is
        False.
    :return: Table containing per-sample variant counts.
    """
    # Add Allele Type annotations to variant MatrixTable.
    _, variant_types = vmt_sample_qc_variant_annotations(
        global_gt=mt.GT, alleles=mt.alleles
    )
    mt = mt.annotate_rows(variant_atypes=variant_types)

    # Annotate the MT with the needed annotations.
    mt = annotate_with_ht(mt, filter_group_ht, filter_missing=True)
    if autosomes_only:
        ab_cutoff = 0.9
        mt = filter_to_autosomes(mt)
        mt = mt.annotate_entries(
            GT=hl.if_else(
                (mt.variant_af > 0.01)
                & ((mt.AD[1] / mt.DP) > ab_cutoff)
                & mt.GT.is_het()
                & ~mt._het_non_ref,
                hl.call(1, 1),
                mt.GT,
            )
        )

    mt = mt.filter_entries(mt.GT.is_non_ref())
    mt = filter_to_adj(mt)

    mt = mt.select_entries("GT", "GQ", "DP")

    # Run Hail's 'vmt_sample_qc' for all requested filter groups.
    qc_expr = vmt_sample_qc(
        global_gt=mt.GT,
        gq=mt.GQ,
        variant_ac=mt.variant_ac,
        variant_atypes=mt.variant_atypes,
        dp=mt.DP,
    )
    ht = mt.select_cols(
        summary_stats=hl.agg.array_agg(
            lambda f: hl.agg.filter(f, qc_expr), mt.filter_groups
        )
    ).cols()
    ht = ht.annotate_globals(
        summary_stats_meta=filter_group_ht.index_globals().filter_group_meta
    )
    ht = ht.checkpoint(hl.utils.new_temp_file("per_sample_counts", "ht"))

    # Add 'n_indel' to the output Table.
    ht = ht.annotate(
        summary_stats=ht.summary_stats.map(
            lambda x: x.annotate(n_indel=x.n_insertion + x.n_deletion)
        )
    )

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

    ht = ht.transmute(
        subset=subset,
        gen_anc=gen_anc,
        summary_stats=hl.zip(ht.summary_stats_meta, ht.summary_stats),
    )

    ht = ht.explode("summary_stats").explode("gen_anc").explode("subset")

    ht = ht.group_by("subset", "gen_anc", variant_filter=ht.summary_stats[0]).aggregate(
        **{
            metric: hl.struct(
                mean=hl.agg.mean(ht.summary_stats[1][metric]),
                min=hl.agg.min(ht.summary_stats[1][metric]),
                max=hl.agg.max(ht.summary_stats[1][metric]),
                quantiles=hl.agg.approx_quantiles(
                    ht.summary_stats[1][metric], [0.0, 0.25, 0.5, 0.75, 1.0]
                ),
            )
            for metric in ht.summary_stats[1]
            if isinstance(ht.summary_stats[1][metric], hl.expr.NumericExpression)
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
    hl._set_flags(use_ssa_logs="1", no_whole_stage_codegen="1")

    data_type = args.data_type
    autosomes_only = args.autosomes_only_stats
    test_dataset = args.test_dataset
    test_partitions = (
        list(range(args.test_n_partitions)) if args.test_n_partitions else None
    )

    # The following four exome partitions have given us trouble with out-of-memory
    # errors in the past. We are testing on these partitions to ensure that the
    # per-sample stats script can handle them before running the entire dataset.
    if data_type == "genomes" and args.test_difficult_partitions:
        raise ValueError(
            "Difficult partitions for testing have only been chosen for exomes data."
        )

    test_partitions = (
        [20180, 40916, 41229, 46085]
        if args.test_difficult_partitions
        else test_partitions
    )
    test = test_partitions or test_dataset
    overwrite = args.overwrite
    ukb_capture_intervals = (
        False if data_type == "genomes" else not args.skip_filter_ukb_capture_intervals
    )
    broad_capture_intervals = (
        False
        if data_type == "genomes"
        else not args.skip_filter_broad_capture_intervals
    )
    rare_variants_afs = args.rare_variants_afs if not args.skip_rare_variants else None
    per_sample_res = get_per_sample_counts(
        test=test,
        data_type=data_type,
        suffix=args.custom_suffix,
        autosomes_only=autosomes_only,
    )
    per_sample_agg_res = get_per_sample_counts(
        test=test,
        data_type=data_type,
        suffix=args.custom_suffix,
        autosomes_only=autosomes_only,
        aggregated=True,
        by_ancestry=args.by_ancestry,
        by_subset=args.by_subset,
    )
    if data_type != "exomes" and args.by_subset:
        raise ValueError("Stratifying by subset is only working on exomes data type.")

    try:
        if args.create_filter_group_ht:
            logger.info("Creating Table of filter groups for summary stats...")
            release_ht = release_sites(data_type=data_type).ht()
            filtering_groups_res = get_summary_stats_filtering_groups(
                data_type=data_type, test=test
            )
            if args.test_n_partitions:
                release_ht = release_ht._filter_partitions(test_partitions)

            get_summary_stats_filter_groups_ht(
                release_ht,
                pass_filters=not args.skip_pass_filters,
                ukb_capture=ukb_capture_intervals,
                broad_capture=broad_capture_intervals,
                by_csqs=not args.skip_by_csqs,
                vep_canonical=args.vep_canonical,
                vep_mane=args.vep_mane,
                rare_variants_afs=rare_variants_afs,
            ).write(filtering_groups_res.path, overwrite=overwrite)

        if args.create_per_sample_counts_ht:
            logger.info(
                "Calculating per-sample variant statistics for %s...", data_type
            )
            if test:
                logger.warning(
                    "This test requires that the full filtering groups Table has been "
                    "created. Please run --create-filter-group-ht without "
                    "--test-dataset or --test-n-partitions if it doesn't already exist."
                )
            filter_groups_ht = get_summary_stats_filtering_groups(
                data_type, test=test, autosomes_only=autosomes_only
            ).ht()
            vds_load_func = (
                get_gnomad_v4_vds
                if data_type == "exomes"
                else get_gnomad_v4_genomes_vds
            )
            mt = vds_load_func(
                test=test_dataset,
                split=True,
                release_only=True,
                filter_partitions=test_partitions,
                filter_variant_ht=filter_groups_ht,
                entries_to_keep=["GT", "GQ", "DP", "AD"],
                annotate_het_non_ref=True,
            ).variant_data

            create_per_sample_counts_ht(mt, filter_groups_ht, autosomes_only).write(
                per_sample_res.path, overwrite=overwrite
            )

        if args.aggregate_sample_stats:
            logger.info("Computing aggregate sample statistics...")
            compute_agg_sample_stats(
                per_sample_res.ht(),
                meta(data_type=data_type).ht(),
                by_ancestry=args.by_ancestry,
                by_subset=args.by_subset,
            ).write(per_sample_agg_res.path, overwrite=overwrite)
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("per_sample_stats.original.non_ukb.big_executors"))


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
        "--test-n-partitions",
        type=int,
        help="Number of partitions to use for testing.",
    )
    parser.add_argument(
        "--test-difficult-partitions",
        help="Whether to test on a set of 4 difficult exome partitions.",
        action="store_true",
    )
    parser.add_argument(
        "--test-dataset",
        help="Test on the test dataset instead of the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--data-type",
        default="exomes",
        choices=["exomes", "genomes"],
        help="Data type (exomes or genomes) to produce summary stats for.",
    )
    parser.add_argument(
        "--create-filter-group-ht",
        help="Create Table of filter groups for summary stats.",
        action="store_true",
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
        "--autosomes-only-stats",
        help="Whether to restrict per-sample summary stats to autosomes only.",
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
        default=SUM_STAT_FILTERS["max_af"],
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
