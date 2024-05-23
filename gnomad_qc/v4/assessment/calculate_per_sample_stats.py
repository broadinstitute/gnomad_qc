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
import itertools
import logging
from copy import deepcopy
from typing import Dict, List, Optional, Union

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
    "loftee_label": deepcopy(LOFTEE_LABELS),  # LOFTEE loss-of-function labels.
    "loftee_flags": ["no_flags", "with_flags"],  # High-confidence LOFTEE flag options.
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

LOF_FILTERS_FOR_COMBO = {"loftee_HC": ["HC"], "lof_csq_set": ["lof"]}
"""
Dictionary of an additional filter group to use for loss-of-function filter
combinations.
"""

LOF_FILTER_COMBOS = [
    ["lof_csq", "loftee_label"],
    ["lof_csq_set", "loftee_HC"],
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


def generate_filter_combinations(
    combos: List[Union[List[str], Dict[str, List[str]]]],
    combo_options: Optional[Dict[str, List[str]]] = None,
) -> List[Dict[str, str]]:
    """
    Generate list of all possible filter combinations from a list of filter options.

    Example input:

    .. code-block:: python

        [
            {'pass_filters': [False, True]},
            {'pass_filters': [False, True], 'capture': ['ukb', 'broad']}
        ]

    Example output:

    .. code-block:: python

        [
            {'pass_filters': False},
            {'pass_filters': True},
            {'pass_filters': False, 'capture': 'ukb'},
            {'pass_filters': False, 'capture': 'broad'},
            {'pass_filters': True, 'capture': 'ukb'},
            {'pass_filters': True, 'capture': 'broad'},
        ]

    :param combos: List of filter groups and their options.
    :param combo_options: Dictionary of filter groups and their options that can be
        supplied if `combos` is a list of lists.
    :return: List of all possible filter combinations for each filter group.
    """
    if isinstance(combos[0], list):
        if combo_options is None:
            raise ValueError(
                "If `combos` is a list of lists, `combo_options` must be provided."
            )
        combos = [{k: combo_options[k] for k in combo} for combo in combos]

    # Create combinations of filter options.

    def _expand_combinations(filter_dict: Dict[str, List[str]]) -> List[Dict[str, str]]:
        """
        Expand filter combinations.

        :param filter_dict: Dictionary of filter options.
        :return: List of dicts of expanded filter options.
        """
        keys, values = zip(*filter_dict.items())
        return [dict(zip(keys, combo)) for combo in itertools.product(*values)]

    # Flatten list of combinations.
    expanded_meta = sum([_expand_combinations(sublist) for sublist in combos], [])

    return expanded_meta


def get_filter_group_meta(
    all_sum_stat_filters: Dict[str, List[str]] = SUM_STAT_FILTERS,
    common_filter_combos: List[List[str]] = COMMON_FILTER_COMBOS,
    common_combo_override: Dict[str, List[str]] = COMMON_FILTERS,
    lof_combos: List[List[str]] = LOF_FILTER_COMBOS,
    lof_combo_override: Dict[str, List[str]] = LOF_FILTERS_FOR_COMBO,
    filter_group_key_rename: Dict[str, str] = MAP_FILTER_FIELD_TO_META,
) -> List[Dict[str, str]]:
    """
    Generate list of filter combinations for summary stats.

    This function combines various filter settings for summary statistics and generates
    all possible filter combinations. It ensures that the generated combinations include
    both common filters and specific loss-of-function (LOF) filters.

    .. note::

        - The "variant_qc" filter group is removed if the value is "none", which can
          lead to a filter group of {} (no filters).
        - The `filter_group_key_rename` parameter can be used to rename keys in the
          `all_sum_stat_filters`, `common_combo_override`, or `lof_combo_override`
          after creating all combinations.

    Example:
        Given the following input:

        .. code-block:: python

            all_sum_stat_filters = {
                "variant_qc": ["none", "pass"],
                "capture": ["1", "2"],
                "max_af": [0.01],
                "lof_csq": ["stop_gained"],
            }
            common_filter_combos = [["variant_qc"], ["variant_qc", "capture"]]
            common_combo_override = {"variant_qc": ["pass"], "capture": ["1"]}
            lof_combos = [
                ["loftee_HC", "loftee_flags"], ["lof_csq", "loftee_HC", "loftee_flags"]
            ]
            lof_combo_override = {"loftee_HC": ["HC"], "loftee_flags": ["with_flags"]}
            filter_group_key_rename = {"lof_csq": "csq","loftee_HC": "loftee_labels"}

        The function will generate the following filter combinations:

        .. code-block:: python

            [
               # Combinations of all common filter keys and their possible values.
                {},
                {'capture': '1'},
                {'capture': '2'},
                {'variant_qc': 'pass'},
                {'variant_qc': 'pass', 'capture': '1'},
                {'variant_qc': 'pass', 'capture': '2'},

                # Combinations of all requested common filter combinations with all
                # possible other filter keys and values.
                {'variant_qc': 'pass', 'max_af': '0.01'},
                {'variant_qc': 'pass', 'csq': 'stop_gained'},
                {'variant_qc': 'pass', 'capture': '1', 'max_af': '0.01'},
                {'variant_qc': 'pass', 'capture': '1', 'csq': 'stop_gained'},

                # Combinations of all requested common filter combinations with all
                # requested LOF filter combination keys and their requested values.
                {
                    'variant_qc': 'pass', 'loftee_labels': 'HC',
                    'loftee_flags': 'with_flags'
                },
                {
                    'variant_qc': 'pass', 'csq': 'stop_gained', 'loftee_labels': 'HC',
                    'loftee_flags': 'with_flags'
                },
                 {
                    'variant_qc': 'pass', 'capture': '1', 'loftee_labels': 'HC',
                    'loftee_flags': 'with_flags'
                },
                {
                    'variant_qc': 'pass', 'capture': '1', 'csq': 'stop_gained',
                    'loftee_labels': 'HC', 'loftee_flags': 'with_flags'
                }
            ]

    :param all_sum_stat_filters: Dictionary of all possible filter types.
    :param common_filter_combos: List of lists of common filter keys to use for creating
        common filter combinations.
    :param common_combo_override: Dictionary of filter groups and their options to
        override the values in `all_sum_stat_filters` for use with values in
        `common_filter_combos`.
    :param lof_combos: List of loss-of-function keys in all_sum_stat_filters to use for
        creating filter combinations.
    :param lof_combo_override: Dictionary of filter groups and their options to override
        the values in `all_sum_stat_filters` for use with values in `lof_combos`.
    :param filter_group_key_rename: Dictionary to rename keys in `all_sum_stat_filters`,
        `common_combo_override`, or `lof_combo_override` to final metadata keys.
    :return: Dictionary of filter field to metadata.
    """
    all_sum_stat_filters = deepcopy(all_sum_stat_filters)

    # Initialize list to store filter metadata combinations.
    filter_combinations = generate_filter_combinations(
        [{f: all_sum_stat_filters[f] for f in combo} for combo in common_filter_combos]
    )

    # Update the common filter combinations with the common filter override, and remove
    # them from the all_sum_stat_filters and into a common filter dictionary.
    all_sum_stat_filters.update(common_combo_override)
    common_filters = {
        k: all_sum_stat_filters.pop(k) for k in set(sum(common_filter_combos, []))
    }

    # Add combinations of common filters with all other filters.
    filter_combinations.extend(
        generate_filter_combinations(
            [c + [f] for c in common_filter_combos for f in all_sum_stat_filters],
            {**all_sum_stat_filters, **common_filters},
        )
    )

    # Add combinations of common filters with LOF specific filters.
    all_sum_stat_filters.update(lof_combo_override)
    filter_combinations.extend(
        generate_filter_combinations(
            [c + f for c in common_filter_combos for f in lof_combos],
            {**all_sum_stat_filters, **common_filters},
        )
    )
    filter_combinations = [
        {
            filter_group_key_rename.get(str(k), str(k)): str(v)
            for k, v in filter_group.items()
            if not (k == "variant_qc" and v == "none")
        }
        for filter_group in filter_combinations
    ]

    return filter_combinations


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
    Create Table annotated with an array of booleans indicating which filter groups a variant belongs to.

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
    filter_group_meta = get_filter_group_meta(all_sum_stat_filters=ss_filters)

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
    ht = ht.select(_no_lcr=filter_exprs["no_lcr"], filter_groups=filter_groups_expr)

    ht = ht.select_globals(filter_group_meta=final_meta)
    logger.info("Filter groups for summary stats: %s", filter_group_meta)

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
