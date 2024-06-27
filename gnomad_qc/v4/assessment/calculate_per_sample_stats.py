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
from typing import Dict, Optional, Tuple

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
from hail.methods.qc import _qc_allele_type
from hail.utils import new_temp_file
from hail.utils.misc import divide_null
from hail.vds.sample_qc import vmt_sample_qc, vmt_sample_qc_variant_annotations

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources import meta
from gnomad_qc.v4.resources.assessment import (
    get_per_sample_counts,
    get_summary_stats_filtering_groups,
)
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
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
        variant_atypes=[_qc_allele_type(ht.alleles[0], ht.alleles[1])],
    )
    ht = ht.select_globals(filter_group_meta=final_meta)
    logger.info("Filter groups for summary stats: %s", filter_group_meta)

    # Filter to only variants that are not in low confidence regions.
    return ht.filter(ht._no_lcr).drop("_no_lcr")


def prepare_mt_for_sample_counts(
    mt: hl.MatrixTable,
    filter_group_ht: hl.Table,
    autosomes_only: bool = False,
) -> hl.MatrixTable:
    """
    Prepare MatrixTable for computing per-sample stats counts.

    This function prepares the MatrixTable for computing per-sample stats counts. It
    does the following:

        - Filters to non-ref genotypes (likely unnecessary if the MT is already a
          variant data MT).
        - Annotates the rows with the filter group metadata from `filter_group_ht`.
        - If `autosomes_only` is True, filters to autosomes only, performs the high
          AB het -> hom alt adjustment of GT, and adds a 'high_ab_het_ref' annotation.
        - Filters to `adj` genotypes and selects only the necessary entries for
          downstream processing.

    :param mt: Input MatrixTable containing variant data. Must have multi-allelic sites
        split.
    :param filter_group_ht: Table containing filter group metadata from
        `get_summary_stats_filter_groups_ht`.
    :param autosomes_only: Boolean indicating whether to filter to autosomes only. When
        True, the high AB het -> hom alt adjustment of GT is performed.
    :return: MatrixTable prepared for computing per-sample stats counts.
    """
    # Annotate the MT with all annotations on the filter group HT.
    mt = annotate_with_ht(mt, filter_group_ht)

    # Filter to non-ref genotypes, which is likely unnecessary if the MT is already a
    # variant data MT.
    mt = mt.filter_entries(mt.GT.is_non_ref())

    # For now, only add a high_ab_het_ref annotation (or correction) if autosomes_only
    # is True.
    # TODO: Modify this when we have fixed the script to correctly handle sex
    #  chromosomes.
    if autosomes_only:
        logger.info(
            "Filtering to autosomes only and performing high AB het -> hom alt "
            "adjustment of GT..."
        )
        mt = filter_to_autosomes(mt)
        base_args = [mt.GT, mt._het_non_ref, mt.variant_af, mt.AD[1] / mt.DP]
        mt = mt.annotate_entries(
            GT=hom_alt_depletion_fix(*base_args, return_bool=False),
            high_ab_het_ref=hom_alt_depletion_fix(*base_args, return_bool=True),
        )

    # Filter to adj genotypes and select only the necessary entries to be localized.
    mt = filter_to_adj(mt)
    mt = mt.select_entries("GT", "GQ", "DP", "high_ab_het_ref")

    # Add the filter_group_meta global in filter_group_ht to the MT.
    mt = mt.annotate_globals(
        filter_group_meta=filter_group_ht.index_globals().filter_group_meta
    )

    return mt


def create_per_sample_counts_ht(
    mt: hl.MatrixTable,
    gq_bins: Tuple[int] = (60,),
    dp_bins: Tuple[int] = (20, 30),
) -> hl.Table:
    """
    Create Table of Hail's sample_qc output broken down by requested variant groupings.

    Useful for finding the number of variants per sample, either all variants, or
    variants fall into specific capture regions, or variants that are rare
    (adj AF <0.1%), or variants categorized by predicted consequences in all, canonical
    or mane transcripts.

    :param mt: Input MatrixTable containing variant data. Must have multi-allelic sites
        split.
    :param gq_bins: Tuple of GQ bins to use for filtering. Default is (60,).
    :param dp_bins: Tuple of DP bins to use for filtering. Default is (20, 30).
    :return: Table containing per-sample variant counts.
    """
    # Run Hail's 'vmt_sample_qc' for all requested filter groups.
    qc_expr = vmt_sample_qc(
        global_gt=mt.GT,
        gq=mt.GQ,
        variant_ac=mt.variant_ac,
        variant_atypes=mt.variant_atypes,
        dp=mt.DP,
        gq_bins=gq_bins,
        dp_bins=dp_bins,
    ).annotate(n_high_ab_het_ref=hl.agg.count_where(mt.high_ab_het_ref))
    ht = mt.select_cols(
        summary_stats=hl.agg.array_agg(
            lambda f: hl.agg.filter(f, qc_expr), mt.filter_groups
        )
    ).cols()
    ht = ht.select_globals(summary_stats_meta=ht.filter_group_meta)
    ht = ht.checkpoint(hl.utils.new_temp_file("per_sample_counts", "ht"))

    # Add 'n_indel' and ' n_non_ref_alleles' to the output Table and rename the DP and
    # GQ fields.
    ht = ht.annotate(
        summary_stats=ht.summary_stats.map(
            lambda x: x.annotate(
                n_non_ref_alleles=x.n_non_ref + x.n_hom_var,
                n_indel=x.n_insertion + x.n_deletion,
                **{
                    f"n_over_{k}_{b}": x[f"bases_over_{k}_threshold"][i]
                    for k, bins in {"gq": gq_bins, "dp": dp_bins}.items()
                    for i, b in enumerate(bins)
                },
            ).drop(*[f"bases_over_{k}_threshold" for k in {"gq", "dp"}])
        )
    )

    return ht


def create_intermediate_mt_for_sample_counts(
    mt: hl.MatrixTable,
    gq_bins: Tuple[int] = (60,),
    dp_bins: Tuple[int] = (20, 30),
) -> hl.Table:
    """
    Create intermediate MatrixTable for computing per-sample stats counts.

    This function creates an intermediate Table for use in
    `create_per_sample_counts_ht`. It does the following:

        - Converts the MT to a HT with a `sample_idx_by_stat` row annotation that is
          a struct of arrays of sample indices for each genotype level stat. The stats
          are: non_ref, het, hom_var, high_ab_het_ref, over_gq_{gq}, and over_dp_{dp}.
        - Changes the `filter_groups` annotation from an array of booleans to an array
          of structs where `group_filter` is the original filter group boolean value,
          mapping to the filter_group_meta, and adds boolean annotations to the struct
          indicating whether the variant is included in each of the following variant
          level stats: singleton, singleton_ti, singleton_tv, insertion, deletion,
          transition, and transversion.
    This structure creates a large intermediate Table that allows for memory efficient
    computation of per-sample stats counts. Smaller datasets do not require this
    intermediate step and can compute per-sample stats counts directly from the MT and
    with the use of Hail's `vmt_sample_qc` function.

    :param mt: Input MatrixTable containing variant data. Must have multi-allelic sites
        split.
    :param gq_bins: Tuple of GQ bins to use for filtering. Default is (60,).
    :param dp_bins: Tuple of DP bins to use for filtering. Default is (20, 30).
    :return: Intermediate Table for computing per-sample stats counts.
    """
    # Convert MT to HT with a row annotation that is an array of all samples entries
    # for that variant.
    ht = mt.localize_entries("_entries", "_cols")

    # Add an index for each sample to the entries array to map entries back to samples.
    entry = hl.map(
        lambda i, e: (i, e.GT, e.GQ, e.DP, e.high_ab_het_ref),
        hl.range(hl.len(ht._cols)),
        ht._entries,
    ).filter(lambda x: hl.is_defined(x[1]))

    # Filter the entry array to genotype should be counted for each genotype level stat.
    stat_expr = hl.struct(
        non_ref=entry.filter(lambda x: True),
        het=entry.filter(lambda x: x[1].is_het()),
        hom_var=entry.filter(lambda x: x[1].is_hom_var()),
        high_ab_het_ref=entry.filter(lambda x: x[4]),
        **{f"over_gq_{gq}": entry.filter(lambda x: x[2] >= gq) for gq in gq_bins},
        **{f"over_dp_{dp}": entry.filter(lambda x: x[3] >= dp) for dp in dp_bins},
    )

    # Annotate the HT with the filter group metadata and the sample index by stat.
    # The filter groups annotation contains an array of structs with a group_filter
    # boolean indicating whether the variant belongs to the filter group and boolean
    # expressions for each variant level stat, i.e. singleton, indicating whether the
    # variant should be counted for that stat.
    # The sample_idx_by_stat annotation contains a struct with arrays of sample indices
    # for each genotype level stat.
    ac1 = ht.variant_ac[1] == 1
    allele_type_expr = hl.struct(
        **{k: ht.variant_atypes[0] == v for k, v in ALLELE_TYPE_MAP.items()}
    )
    ht = ht.select(
        filter_groups=ht.filter_groups.map(
            lambda x: hl.struct(
                group_filter=x,
                singleton=x & ac1,
                singleton_ti=x & allele_type_expr.transition & ac1,
                singleton_tv=x & allele_type_expr.transversion & ac1,
                **{k: x & allele_type_expr[k] for k in ALLELE_TYPE_MAP.keys()},
            )
        ),
        sample_idx_by_stat=stat_expr.annotate(
            **{k: stat_expr[k].map(lambda x: x[0]) for k in stat_expr}
        ),
    )

    # Annotate the HT with globals for the sample IDs and filter group metadata.
    ht = ht.select_globals(
        samples=ht._cols.map(lambda x: x.s),
        filter_group_meta=ht.filter_group_meta,
    )

    return ht


def create_per_sample_counts_from_intermediate_ht(ht: hl.Table) -> hl.Table:
    """
    Create Table of sample QC metrics broken down by requested variant groupings.

    .. warning::

        This function is memory intensive and should be run on a cluster with
        n1-highmem-8 workers and big-executors enabled. It also requires shuffling,
        so a cluster with all workers is recommended.

    Takes an intermediate Table (returned by `create_intermediate_mt_for_sample_counts`)
    with an array of variant level stats by filter groups and an array of sample
    indices for each genotype level stat.

    This function restructures the Table a few times to get around memory errors that
    occur when trying to aggregate the Table directly.

    The following steps are taken:

        - The Table is grouped by filter group (including all variant level stats) and
          aggregated to get the count of variants for each sample in each filter group.
          Leads to a Table where the number of rows is approximately the number of
          filter groups times the possible permutations of variant level stats. If
          this number is too small, it might be better to use
          `create_per_sample_counts_ht`.
        - Sample IDs are mapped to the sample indices and the Table is exploded so that
          each row is counts for a sample and a filter group. The number of rows in the
          Table is approximately the number of samples times the number of rows in the
          previous Table.
        - Group the Table by sample to get a struct of the count of variants by sample
          QC metric for each sample in each filter group.
        - Add sample QC stats that are combinations of other stats.

    :param ht: Input Table containing variant data. Must have multi-allelic sites split.
    :return: Table containing per-sample variant counts.
    """
    # Get the number of filter groups, samples, and variants in the Table to determine
    # a reasonable number of partitions to use for the aggregated filter group HT, and
    # the number of random groups needed to have approximately the desired partitions.
    n_filter_groups = len(hl.eval(ht.filter_group_meta))
    n_samples = hl.eval(hl.len(ht.samples))
    n_variants = hl.eval(ht.count())
    # 2^3 = 8 because there are about 3 variant level stats: singletons (singleton_ti
    # and singleton_tv will never overlap, and one of them will always also be a
    # singleton), transition/transversion, insertion/deletion.
    n_filter_permutations = n_filter_groups * 8
    # Split the Table into approximately 10,000 variants per partition.
    n_mem_partitions = max(int(n_variants / 10000), n_filter_permutations)
    # Compute the number of random groups needed to have approximately n_mem_partitions.
    n_rand_groups = max(int(n_mem_partitions / n_filter_permutations), 1)
    logger.info(
        "Aggregating counts for each sample in each filter group...\n"
        "\tNumber of samples: %d\n"
        "\tNumber of variants: %d\n"
        "\tNumber of filter groups: %d\n"
        "\tApproximate number of filter permutations: %d\n\n"
        "\tCalculated number of partitions to improve memory usage: %d\n"
        "\tNumber of random groups to have approximately the above partitions: %d",
        n_samples,
        n_variants,
        n_filter_groups,
        n_filter_permutations,
        n_mem_partitions,
        n_rand_groups,
    )

    # Group by filter groups and the random group to get the count of variants matching
    # the filter group and random group for each sample.
    ht = ht.annotate(_rand_group=hl.rand_int32(n_rand_groups))
    ht = ht.group_by("filter_groups", "_rand_group").aggregate(
        counts=hl.struct(
            **{
                k: hl.agg.explode(lambda x: hl.agg.counter(x), v)
                for k, v in ht.sample_idx_by_stat.items()
            }
        )
    )
    tmp_path = new_temp_file("agg_filter_groups", "ht")
    ht = ht.write(tmp_path)
    ht = hl.read_table(tmp_path, _n_partitions=n_mem_partitions)
    logger.info("Reading in as %s partitions.", ht.n_partitions())

    # Map the sample indices to the sample IDs and fill in samples with no counts for
    # the variant level stat with 0, then explode the stat counts array. This
    # transforms the Table so there is one row for each sample and filter group. Since
    # we need aggregate counts for each sample by filter group, and the number of
    # samples is much larger than the number of filter groups, this allows for
    # repartitioning the Table to have more partitions than the number of filter groups,
    # helping the final aggregation get around memory errors.
    ht = (
        ht.group_by("filter_groups")
        .aggregate(
            counts=hl.zip(
                ht.samples,
                hl.agg.array_agg(
                    lambda i: hl.struct(
                        **{k: hl.agg.sum(v.get(i, 0)) for k, v in ht.counts.items()}
                    ),
                    hl.range(n_samples),
                ),
            )
        )
        .explode("counts")
    )
    ht = ht.annotate(s=ht.counts[0], counts=ht.counts[1])

    # After the explode, the number of rows is much larger, at approximately the number
    # of samples times the number of filter groups, so it's important to repartition the
    # Table. This is done as a repartition on read, and the key_by is needed because
    # otherwise the repartitioning will not work, and will only be as many partitions
    # as the number of filter_groups.
    n_partitions = max(int((n_filter_permutations * n_samples) / 50000), 50)
    tmp_path = new_temp_file("stat_counts_explode_sample", "ht")
    ht.key_by("s").write(tmp_path)
    ht = hl.read_table(tmp_path, _n_partitions=n_partitions)
    logger.info("Reading in as %s partitions.", ht.n_partitions())

    # Group by sample to get a struct of the count of variants for each sample QC stat
    # for each sample in each filter group. Uses 'filter_groups' annotation to filter to
    # variants that belong to each filter group, and the 'counts' annotation to get
    # the count of variants for each genotype level stat. The 'non_ref_alleles' genotype
    # level count is used with all variant level stat filter to get the count of
    # variants for each.
    ht = ht.annotate(
        counts=ht.counts.annotate(non_ref_alleles=ht.counts.non_ref + ht.counts.hom_var)
    )
    variant_filters = [s for s in ht.filter_groups[0] if s != "group_filter"]
    agg_func = lambda x, f, s: hl.agg.sum(hl.int(x[f]) * ht.counts[s])
    ht = (
        ht.group_by("s")
        .aggregate(
            summary_stats=hl.agg.array_agg(
                lambda x: hl.struct(
                    **{f"n_{s}": agg_func(x, "group_filter", s) for s in ht.counts},
                    **{
                        f"n_{f}": agg_func(x, f, "non_ref_alleles")
                        for f in variant_filters
                    },
                ),
                ht.filter_groups,
            )
        )
        .checkpoint(new_temp_file("per_sample_counts", "ht"))
    )

    # Add sample QC stats that are combinations of other stats.
    ratios = {
        "r_ti_tv": ("n_transition", "n_transversion"),
        "r_ti_tv_singleton": ("n_singleton_ti", "n_singleton_tv"),
        "r_het_hom_var": ("n_het", "n_hom_var"),
        "r_insertion_deletion": ("n_insertion", "n_deletion"),
    }.items()
    ht = ht.annotate(
        summary_stats=ht.summary_stats.map(
            lambda x: x.annotate(
                n_indel=x.n_insertion + x.n_deletion,
                n_snp=x.n_transition + x.n_transversion,
                **{k: divide_null(hl.float64(x[n]), x[d]) for k, (n, d) in ratios},
            )
        )
    )
    ht = ht.select_globals(summary_stats_meta=ht.filter_group_meta)

    n_partitions = max(int(n_samples * 0.001), min(50, n_samples))
    logger.info("Naive coalescing to %s partitions.", n_partitions)
    ht = ht.naive_coalesce(n_partitions)

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

    ht = ht.annotate(subset=subset, gen_anc=gen_anc)
    ht = ht.explode("gen_anc").explode("subset")

    ht = ht.group_by("subset", "gen_anc").aggregate(
        summary_stats=hl.agg.array_agg(
            lambda x: hl.struct(
                **{
                    m: hl.struct(
                        mean=hl.agg.mean(x[m]),
                        min=hl.agg.min(x[m]),
                        max=hl.agg.max(x[m]),
                        quantiles=hl.agg.approx_quantiles(
                            x[m], [0.0, 0.25, 0.5, 0.75, 1.0]
                        ),
                    )
                    for m in x
                }
            ),
            ht.summary_stats,
        )
    )
    ht = ht.annotate(
        summary_stats=hl.zip(ht.summary_stats_meta, ht.summary_stats)
    ).select_globals()
    ht = ht.explode("summary_stats")
    ht = ht.annotate(
        filter_group_meta=ht.summary_stats[0], summary_stats=ht.summary_stats[1]
    )
    ht = ht.key_by(*ht.key, "filter_group_meta")

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
    create_filter_group = args.create_filter_group_ht
    create_per_sample_counts = args.create_per_sample_counts_ht
    autosomes_only = args.autosomes_only_stats
    test_dataset = args.test_dataset
    test_gene = args.test_gene
    test_difficult_partitions = args.test_difficult_partitions
    test_partitions = (
        list(range(args.test_n_partitions)) if args.test_n_partitions else None
    )

    if create_filter_group:
        err_msg = ""
        if test_partitions and create_per_sample_counts:
            err_msg = (
                "Cannot test on a custom number of partitions (--test-n-partitions) "
                "when using both --create-filter-group-ht and "
                "--create-per-sample-counts-ht because there is not an overlap in "
                "partitions."
            )
        if test_difficult_partitions:
            err_msg = (
                "Difficult partitions (--test-difficult-partitions) are only chosen "
                "for testing per-sample counts, we recommend using only --test-gene to "
                "test that --create-filter-group-ht is working as expected."
            )
        if test_dataset and not create_filter_group:
            err_msg = (
                "The use of the test VDS (--test-dataset) is only relevant when also"
                " testing per-sample counts (--create-per-sample-counts-ht), we "
                "recommend using only --test-gene to test that --create-filter-group-ht"
                " is working as expected."
            )
        if err_msg:
            raise ValueError(
                err_msg
                + " For a quick test of both --create-filter-group-ht and "
                "--create-per-sample-counts-ht arguments of the pipeline at the "
                "same time use '--test-gene --test-dataset'. The full run of the "
                "--create-filter-group-ht step is relatively quick and after tests "
                "are complete on the filter group creation, the full run of "
                "--create-filter-group-ht should be done to further test "
                "--create-per-sample-counts-ht for memory errors."
            )

    # The following four exome partitions have given us trouble with out-of-memory
    # errors in the past. We are testing on these partitions to ensure that the
    # per-sample stats script can handle them before running the entire dataset.
    if test_difficult_partitions:
        if data_type == "genomes":
            raise ValueError(
                "Difficult partitions for testing have only been chosen for exomes."
            )
        if test_dataset or test_partitions:
            raise ValueError(
                "Cannot test on difficult partitions (--test-difficult-partitions) and "
                "test dataset (--test-dataset) or a custom number of partitions "
                "(--test-n-partitions) at the same time. Difficult partitions only "
                "apply to the full exomes VDS."
            )
        logger.info(
            "Testing on difficult exome partitions to make sure the tests can pass "
            "without memory errors..."
        )
        test_partitions = [20180, 40916, 41229, 46085]

    if test_gene and test_partitions:
        raise ValueError(
            "Cannot use --test-gene and --test-n-partitions or "
            "--test-difficult-partitions at the same time."
        )
    filter_intervals = ["chr1:55039447-55064852"] if test_gene else None

    test = test_partitions or test_dataset or test_gene
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
    temp_intermediate_ht_path = get_checkpoint_path(
        "per_sample_summary_stats_intermediate" + (".test" if test else "")
    )
    if data_type != "exomes" and args.by_subset:
        raise ValueError("Stratifying by subset is only working on exomes data type.")

    try:
        if create_filter_group:
            logger.info(
                "Creating Table of filter groups for %s summary stats...", data_type
            )
            release_ht = release_sites(data_type=data_type).ht()
            filtering_groups_res = get_summary_stats_filtering_groups(
                data_type=data_type, test=test
            )
            if test_partitions:
                release_ht = release_ht._filter_partitions(test_partitions)
            if test_gene:
                release_ht = hl.filter_intervals(
                    release_ht,
                    [
                        hl.parse_locus_interval(x, reference_genome="GRCh38")
                        for x in filter_intervals
                    ],
                )

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

        if args.create_intermediate_mt_for_sample_counts or (
            create_per_sample_counts and not args.use_intermediate_mt_for_sample_counts
        ):
            if test and not test_gene:
                logger.warning(
                    "This test requires that the full filtering groups Table has been "
                    "created. Please run --create-filter-group-ht without "
                    "--test-gene or --test-n-partitions if it doesn't already exist."
                )

            filter_groups_ht = get_summary_stats_filtering_groups(
                data_type, test=test_gene
            ).ht()

            vds_load_func = (
                get_gnomad_v4_vds
                if data_type == "exomes"
                else get_gnomad_v4_genomes_vds
            )
            mt = prepare_mt_for_sample_counts(
                vds_load_func(
                    test=test_dataset,
                    split=True,
                    release_only=True,
                    filter_partitions=test_partitions,
                    filter_variant_ht=filter_groups_ht,
                    filter_intervals=filter_intervals,
                    split_reference_blocks=False,
                    entries_to_keep=["GT", "GQ", "DP", "AD"],
                    annotate_het_non_ref=True,
                ).variant_data,
                filter_groups_ht,
                autosomes_only,
            )

            if args.create_intermediate_mt_for_sample_counts:
                logger.info(
                    "Creating intermediate MatrixTable for per-sample counts..."
                )
                create_intermediate_mt_for_sample_counts(mt).write(
                    temp_intermediate_ht_path, overwrite=overwrite
                )

        if create_per_sample_counts:
            logger.info(
                "Calculating per-sample variant statistics for %s...", data_type
            )
            if args.use_intermediate_mt_for_sample_counts:
                logger.info("Using intermediate MatrixTable for per-sample counts...")
                ht = hl.read_table(temp_intermediate_ht_path)
                if args.test_n_partitions:
                    ht = ht._filter_partitions(test_partitions)
                if test_gene:
                    ht = hl.filter_intervals(
                        ht,
                        [
                            hl.parse_locus_interval(x, reference_genome="GRCh38")
                            for x in filter_intervals
                        ],
                    )
                create_per_sample_counts_from_intermediate_ht(ht).write(
                    per_sample_res.path, overwrite=overwrite
                )
            else:
                create_per_sample_counts_ht(mt).write(
                    per_sample_res.path, overwrite=overwrite
                )

        if args.aggregate_sample_stats:
            logger.info("Computing aggregate sample statistics...")
            if test:
                logger.warning(
                    "Using whatever per-sample counts testing Table that was most "
                    "recently created."
                )
            compute_agg_sample_stats(
                per_sample_res.ht(),
                meta(data_type=data_type).ht(),
                by_ancestry=args.by_ancestry,
                by_subset=args.by_subset,
            ).write(per_sample_agg_res.path, overwrite=overwrite)
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
        "--test-gene",
        help=(
            "Filter Tables/VDS to only the PCSK9 gene for testing. Recommended for a "
            "quick test (if used with --test-dataset) that the full pipeline is "
            "working. This is not recommended for testing potential memory errors in "
            "--create-per-sample-counts-ht."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-n-partitions",
        type=int,
        help=(
            "Number of partitions to use for testing --create-per-sample-counts-ht. "
            "Recommended as a quick test of --create-per-sample-counts-ht "
            "if combined with --test-dataset. It can also be useful as a test for "
            "memory errors in --create-per-sample-counts-ht when not combined with "
            "--test-dataset."
        ),
    )
    parser.add_argument(
        "--test-difficult-partitions",
        help=(
            "Whether to test on a set of 4 exome partitions that have been "
            "particularly difficult and caused memory errors in other parts of the "
            "QC workflow. Recommended to test for memory errors in "
            "--create-per-sample-counts-ht."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-dataset",
        help=(
            "Use the test VDS instead of the full VDS. Recommended for testing "
            "--create-per-sample-counts-ht. Use in combination with --test-gene or"
            "--test-n-partitions for a quick test of --create-per-sample-counts-ht "
            "or alone as a test that there are no memory errors when aggregating "
            "stats for full exomes."
        ),
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
        "--create-intermediate-mt-for-sample-counts",
        help=(
            "Create intermediate MatrixTable for per-sample variant counts. The output "
            "MatrixTable will be written to a temporary location!"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--create-per-sample-counts-ht",
        help="Create per-sample variant counts Table.",
        action="store_true",
    )
    parser.add_argument(
        "--use-intermediate-mt-for-sample-counts",
        help=(
            "Use intermediate MatrixTable for per-sample variant counts. This is only "
            "relevant when using --create-per-sample-counts-ht. Note that the this "
            "step is memory intensive and should be run on a cluster with n1-highmem-8 "
            "workers and big-executors enabled. It also requires shuffling, so a "
            "cluster with all workers is recommended."
        ),
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
