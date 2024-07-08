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
from typing import Callable, Dict, List, Optional, Tuple

import hail as hl
from gnomad.assessment.summary_stats import (
    get_summary_stats_csq_filter_expr,
    get_summary_stats_filter_group_meta,
    get_summary_stats_variant_filter_expr,
)
from gnomad.resources.resource_utils import TableResource
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import annotate_with_ht, get_adj_expr
from gnomad.utils.filtering import filter_to_adj
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
from hail.vds.sample_qc import vmt_sample_qc

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
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


def process_test_args(
    data_type: str,
    test_dataset: bool = False,
    test_gene: bool = False,
    test_partitions: Optional[List[int]] = None,
    test_difficult_partitions: bool = False,
    create_filter_group: bool = False,
    create_intermediate: bool = False,
    create_per_sample_counts: bool = False,
    use_intermediate: bool = False,
    sex_chr_only: bool = False,
    autosomes_only: bool = False,
) -> Tuple[Optional[List[hl.tinterval]], Optional[List[int]]]:
    """
    Process test arguments for the per-sample stats pipeline.

    Raises a ValueError if the test argument combination is invalid or returns the
    intervals and partitions to filter to for testing.

    The test argument combination is invalid if:

        - `create_filter_group` is True, and:

            - `test_partitions` and `create_per_sample_counts` are both True because the
              partitions in the release HT and the raw VDS are different.
            - `test_difficult_partitions` is True because difficult partitions are only
              relevant for testing per-sample counts.
            - `test_dataset` is True and `create_per_sample_counts` is False because the
              test dataset is only relevant when testing per-sample counts.

        - `test_difficult_partitions` is True, and:

            - `data_type` is "genomes" because difficult partitions are only chosen for
              the raw exome VDS.
            - `test_dataset` or `test_partitions` is True because difficult partitions
              only apply to the full exomes VDS.

        - `test_gene` is True, and `test_partitions` or `test_difficult_partitions` is
          True because there is no (or likely no) overlap in partitions between the
          gene test intervals and the test partitions or difficult partitions.

    :param data_type: Data type of the dataset.
    :param test_dataset: Boolean indicating whether to use the test dataset instead of
        the full dataset. Default is False.
    :param test_gene: Boolean indicating whether to filter the dataset to PCSK9 and/or
        TRPC5 and ZFY (depends on values of `sex_chr_only` and `autosomes_only` because
        PCSK9 is on chr1 and TRPC5 is on chrX and ZFY is on chrY). Default is False.
    :param test_partitions: Optional list of partitions to test on. Default is None.
    :param test_difficult_partitions: Boolean indicating whether to test on difficult
        partitions. Default is False.
    :param create_filter_group: Boolean indicating whether the filter group HT is being
        created. Default is False.
    :param create_intermediate: Boolean indicating whether the intermediate MT is being
        created. Default is False.
    :param create_per_sample_counts: Boolean indicating whether the per-sample counts HT
        is being created. Default is False.
    :param use_intermediate: Boolean indicating whether the intermediate MT is being
        used. Default is False.
    :param sex_chr_only: Boolean indicating whether the dataset is being filtered to
        sex chromosomes only. Default is False.
    :param autosomes_only: Boolean indicating whether the dataset is being filtered to
        autosomes only. Default is False.
    :return: Tuple of filter intervals and partitions to filter to for testing.
    """
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
        if test_dataset and not create_per_sample_counts:
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

    if (
        create_intermediate
        or (create_per_sample_counts and not use_intermediate)
        and not test_gene
    ):
        logger.warning(
            "This test requires that the full filtering groups Table has been "
            "created. Please run --create-filter-group-ht without "
            "--test-gene or --test-n-partitions if it doesn't already exist."
        )

    filter_intervals = None
    if test_gene:
        filter_intervals = []
        if not sex_chr_only:
            # Filter to 10kb in PCSK9.
            filter_intervals.append("chr1:55039447-55064852")
        if not autosomes_only:
            # Filter to TRPC5 on chrX and ZFY on chrY.
            filter_intervals.extend(
                ["chrX:111776000-111786000", "chrY:2935281-2982506"]
            )
        filter_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in filter_intervals
        ]

    return filter_intervals, test_partitions


def get_capture_filter_exprs(
    ht: hl.Table,
    ukb_capture: bool = True,
    broad_capture: bool = True,
) -> Dict[str, hl.expr.BooleanExpression]:
    """
    Get filter expressions for UK Biobank and Broad capture regions.

    :param ht: Table containing variant annotations. The following annotations are
        required: 'region_flags'.
    :param ukb_capture: Expression for variants that are in UKB capture intervals.
        Default is True.
    :param broad_capture: Expression for variants that are in Broad capture intervals.
        Default is True.
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
    capture_regions: bool = False,
    vep_canonical: bool = True,
    vep_mane: bool = False,
    rare_variants_afs: Optional[list[float]] = None,
) -> hl.Table:
    """
    Create Table annotated with an array of booleans indicating whether a variant belongs to certain filter groups.

    A 'filter_groups' annotation is added to the Table containing an ArrayExpression of
    BooleanExpressions for the following filter groups:

        - All variants
        - Variants that pass all variant QC filters
        - Variants in UK Biobank capture regions
        - Variants in Broad capture regions
        - Variants in the intersect of UK Biobank and Broad capture regions
        - Variants in the union of UK Biobank and Broad capture regions
        - Variant consequence: loss-of-function, missense, and synonymous
        - Rare variants with allele frequencies less than the thresholds in
          `rare_variants_afs`

    A 'filter_group_meta' global annotation is added to the Table containing an array
    of dictionaries detailing the filters used in each filter group.

    :param ht: Table containing variant annotations. The following annotations are
        required: 'freq', 'filters', 'region_flags', and 'vep'.
    :param capture_regions: Whether to include filtering groups for variants in UK
        Biobank and Broad capture regions. Default is False.
    :param vep_canonical: Whether to filter to only canonical transcripts. If trying
        count variants in all transcripts, set it to False. Default is True.
    :param vep_mane: Whether to filter to only MANE transcripts. Default is False.
    :param rare_variants_afs: Optional list of allele frequency thresholds to use for
        rare variant filtering.
    :return: Table containing an ArrayExpression of filter groups for summary stats.
    """
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

    # Create filter expressions for the consequence types.
    csq_filter_expr = get_summary_stats_csq_filter_expr(
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

    # Create filter expressions for LCR, variant QC filters, and rare variant AFs if
    # requested.
    filter_exprs = {
        "all_variants": hl.literal(True),
        **(get_capture_filter_exprs(ht) if capture_regions else {}),
        **get_summary_stats_variant_filter_expr(
            ht,
            filter_lcr=True,
            filter_expr=ht.filters,
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
        sex_chr_nonpar_group=(
            hl.case()
            .when(ht.locus.in_autosome_or_par(), "autosome_or_par")
            .when(ht.locus.in_x_nonpar(), "x_nonpar")
            .when(ht.locus.in_y_nonpar(), "y_nonpar")
            .or_missing()
        ),
        filter_groups=filter_groups_expr,
        variant_ac=[hl.missing(hl.tint32), ht.freq[0].AC],
        variant_af=ht.freq[0].AF,
        variant_atypes=[_qc_allele_type(ht.alleles[0], ht.alleles[1])],
    )
    ht = ht.select_globals(filter_group_meta=final_meta)
    logger.info("Filter groups for summary stats: %s", filter_group_meta)

    # Filter to only variants that are not in low confidence regions.
    return ht.filter(ht._no_lcr).drop("_no_lcr")


def load_mt_for_sample_counts(
    data_type: str,
    filter_group_ht: hl.Table,
    adjust_ploidy: bool = False,
    **kwargs,
) -> hl.MatrixTable:
    """
    Load VDS variant data and prepare MatrixTable for computing per-sample stats counts.

    This function loads the VDS for the requested data type and prepares the variant
    data MatrixTable for computing per-sample stats counts. It does the following:

        - Filters to non-ref genotypes (likely unnecessary if the MT is already a
          variant data MT).
        - Annotates the rows with the filter group metadata from `filter_group_ht`.
        - Performs the high AB het -> hom alt adjustment of GT, and adds a
          'high_ab_het_ref' annotation.
        - Filters to `adj` genotypes and selects only the necessary entries for
          downstream processing.

    :param data_type: Data type of the MatrixTable to return. Either 'genomes' or
        'exomes'.
    :param filter_group_ht: Table containing filter group metadata from
        `get_summary_stats_filter_groups_ht`.
    :param adjust_ploidy: Whether to adjust ploidy. If `data_type` is 'exomes', ploidy
        is adjusted before determining the adj annotation. If `data_type` is 'genomes',
        the ploidy adjustment is done after determining the adj annotation. This
        difference is only added for consistency with v3.1, where we added the adj
        annotation before adjusting ploidy, but the correct thing to do is to adjust
        ploidy before determining the adj annotation. Default is False.
    :param kwargs: Additional keyword arguments to pass to the VDS loading function.
    :return: MatrixTable prepared for computing per-sample stats counts.
    """
    # Load the VDS for the requested data type and filter to the variants in the
    # filter group HT.
    vds_load_func = (
        get_gnomad_v4_vds if data_type == "exomes" else get_gnomad_v4_genomes_vds
    )
    mt = vds_load_func(
        split=True,
        release_only=True,
        split_reference_blocks=False,
        entries_to_keep=["GT", "GQ", "DP", "AD"],
        annotate_het_non_ref=True,
        filter_variant_ht=filter_group_ht,
        **kwargs,
    ).variant_data

    # Annotate the MT with all annotations on the filter group HT.
    mt = annotate_with_ht(mt, filter_group_ht)

    # Filter to non-ref genotypes, which is likely unnecessary if the MT is already a
    # variant data MT.
    mt = mt.filter_entries(mt.GT.is_non_ref())

    # NOTE: The correct thing to do here is to adjust ploidy before determining the adj
    # annotation because haploid GTs have different adj filtering criteria, but the
    # option to adjust ploidy after adj is included for consistency with v3.1, where we
    # added the adj annotation before adjusting for sex ploidy.
    gt_expr = mt.GT
    adj_expr = get_adj_expr(gt_expr, mt.GQ, mt.DP, mt.AD)
    if adjust_ploidy:
        logger.info("Adjusting ploidy on sex chromosomes...")
        gt_expr = adjusted_sex_ploidy_expr(
            mt.locus, gt_expr, mt.meta.sex_imputation.sex_karyotype
        )
        if data_type == "exomes":
            adj_expr = get_adj_expr(gt_expr, mt.GQ, mt.DP, mt.AD)

    logger.info("Performing high AB het -> hom alt adjustment of GT...")
    base_args = [gt_expr, mt._het_non_ref, mt.variant_af, mt.AD[1] / mt.DP]
    mt = mt.annotate_entries(
        GT=hom_alt_depletion_fix(*base_args, return_bool=False),
        adj=adj_expr,
        high_ab_het_ref=hom_alt_depletion_fix(*base_args, return_bool=True),
    )

    # Filter to adj genotypes and select only the necessary entries to be localized.
    mt = filter_to_adj(mt)
    mt = mt.select_entries("GT", "GQ", "DP", "high_ab_het_ref").select_cols()

    # Add the filter_group_meta global in filter_group_ht to the MT.
    mt = mt.annotate_globals(
        filter_group_meta=filter_group_ht.index_globals().filter_group_meta
    )

    return mt


def annotate_per_sample_stat_combinations(
    ht: hl.Table,
    sums: Dict[str, List[str]] = {
        "n_indel": ["n_insertion", "n_deletion"],
        "n_snp": ["n_transition", "n_transversion"],
    },
    diffs: Dict[str, Tuple[str, str]] = {},
    ratios: Dict[str, Tuple[str, str]] = {
        "r_ti_tv": ("n_transition", "n_transversion"),
        "r_ti_tv_singleton": ("n_singleton_ti", "n_singleton_tv"),
        "r_het_hom_var": ("n_het", "n_hom_var"),
        "r_insertion_deletion": ("n_insertion", "n_deletion"),
    },
    additional_stat_combos: Optional[Dict[str, Callable]] = {},
) -> hl.Table:
    """
    Annotate the per-sample stats Table with ratios of other per-sample stats.

    :param ht: Input Table containing per-sample stats.
    :param sums: Dictionary of per-sample stats to sum. The key is the name of the sum
        and the value is a List of the stats to sum. Default is to sum the number of
        insertions and deletions to get the total number of indels, and the number of
        transitions and transversions to get the total number of transitions.
    :param diffs: Dictionary of per-sample stats to subtract. The key is the name of the
        difference and the value is a tuple of the stats to subtract, where the second
        stat is subtracted from the first. Default is {}.
    :param ratios: Dictionary of ratios to compute. The key is the name of the ratio
        and the value is a tuple of the numerator and denominator per-sample stats.
        Default is to compute the transition/transversion ratio, the singleton
        transition/transversion ratio, the heterozygous/homozygous variant ratio, and
        the insertion/deletion ratio.
    :param additional_stat_combos: Optional dictionary of additional per-sample stat
        combinations to compute. The key is the name of the stat and the value is a
        function that takes the per-sample stats struct and returns the computed stat.
        Default is {}.
    :return: Table containing per-sample stats annotated with the requested ratios.
    """
    return ht.annotate(
        summary_stats=ht.summary_stats.map(
            lambda x: x.annotate(
                **{k: hl.sum([x[s] for s in v]) for k, v in sums.items()},
                **{k: x[d1] - x[d2] for k, (d1, d2) in diffs.items()},
                **{
                    k: divide_null(hl.float64(x[n]), x[d])
                    for k, (n, d) in ratios.items()
                },
                **{k: v(x) for k, v in additional_stat_combos.items()},
            )
        )
    )


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
    )

    # Add aggregations for the number of hemizygous variants and high AB het ref
    # genotypes and rename the DP and GQ fields.
    qc_expr = qc_expr.annotate(
        n_hemi_var=hl.agg.count_where(mt.GT.is_haploid()),
        n_high_ab_het_ref=hl.agg.count_where(mt.high_ab_het_ref),
        **{
            f"n_over_{k}_{b}": qc_expr[f"bases_over_{k}_threshold"][i]
            for k, bins in {"gq": gq_bins, "dp": dp_bins}.items()
            for i, b in enumerate(bins)
        },
    ).drop(*[f"bases_over_{k}_threshold" for k in {"gq", "dp"}])

    ht = (
        mt.select_cols(
            _ss=hl.agg.group_by(
                mt.sex_chr_nonpar_group,
                hl.agg.array_agg(lambda f: hl.agg.filter(f, qc_expr), mt.filter_groups),
            ).items()
        )
        .cols()
        .explode("_ss")
    )
    ht = ht.transmute(sex_chr_nonpar_group=ht._ss[0], summary_stats=ht._ss[1])

    # Add 'n_indel' and ' n_non_ref_alleles' to the output Table.
    ht = annotate_per_sample_stat_combinations(
        ht, diffs={"n_hom_var": ("n_hom_var", "n_hemi_var")}, sums={}, ratios={}
    )
    ht = annotate_per_sample_stat_combinations(
        ht,
        sums={
            "n_indel": ["n_insertion", "n_deletion"],
            "n_non_ref_alleles": ["n_non_ref", "n_hom_var"],
        },
    )

    return ht.select_globals(summary_stats_meta=ht.filter_group_meta)


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
        hom_var=entry.filter(lambda x: x[1].is_hom_var() & ~x[1].is_haploid()),
        hemi_var=entry.filter(lambda x: x[1].is_haploid()),
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
        "sex_chr_nonpar_group",
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
        n1-highmem-8 workers. It also requires shuffling, so a cluster with all
        workers is recommended.

    Takes an intermediate Table (returned by `create_intermediate_mt_for_sample_counts`)
    with an array of variant level stats by filter groups and an array of sample
    indices for each genotype level stat.

    This function restructures the Table a few times to get around memory errors that
    occur when trying to aggregate the Table directly.

    The following steps are taken:

        - The Table is grouped by filter group (including all variant level stats) and
          aggregated to get the count of variants for each sample in each filter group.
          Leads to a Table where the number of rows is approximately the number of
          filter groups times the possible permutations of variant level stats.

            - Since this number can sometimes be too few partitions to get around the
              memory errors, we approximate a reasonable number of partitions as the
              number required to split the table into approximately 10,000 variants per
              partition.
            - Then we artificially increase the number of aggregation groups by
              assigning a random group number to each variant, so that when the Table
              is grouped by both the filter group and the random group instead
              of only the filter group, the number of rows will be approximately the
              number of desired partitions.
            - Therefore, the number of random groups is the ideal number of partitions
              divided by the number of possible filter groups.

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
    ht = ht.group_by("sex_chr_nonpar_group", "filter_groups", "_rand_group").aggregate(
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
        ht.group_by("sex_chr_nonpar_group", "filter_groups")
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
    # Table. In this case, we try to have about 50,000 rows per partition, which seems
    # to work well, but only going as low as 50 partitions. The repartition is done as a
    # repartition on read, and the key_by is needed because otherwise the
    # repartitioning will not work, and will only be as many partitions as the number
    # of filter_groups.
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
        ht.group_by("s", "sex_chr_nonpar_group")
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

    # Add 'n_indel', 'n_snp' and sample QC stats that are combinations of other stats.
    ht = annotate_per_sample_stat_combinations(ht)

    n_partitions = max(int(n_samples * 0.001), min(50, n_samples))
    logger.info("Naive coalescing to %s partitions.", n_partitions)
    ht = ht.naive_coalesce(n_partitions)

    return ht.select_globals(summary_stats_meta=ht.filter_group_meta)


def combine_autosome_and_sex_chr_stats(
    autosome_ht: hl.Table,
    sex_chr_ht: hl.Table,
) -> hl.Table:
    """
    Combine autosomal and sex chromosome per-sample stats Tables.

    :param autosome_ht: Table containing per-sample stats for autosomes.
    :param sex_chr_ht: Table containing per-sample stats for sex chromosomes.
    :return: Table containing combined per-sample stats.
    """
    ht = autosome_ht.union(sex_chr_ht)
    ht = ht.group_by("s", "sex_chr_nonpar_group").aggregate(
        summary_stats=hl.agg.array_agg(
            lambda x: hl.struct(
                **{
                    k: hl.agg.sum(hl.or_else(v, 0))
                    for k, v in x.items()
                    if k.startswith("n_")
                },
            ),
            ht.summary_stats,
        )
    )
    return annotate_per_sample_stat_combinations(ht, sums={})


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
        if by_subset:
            subset_expr = hl.if_else(meta_s.project_meta.ukb_sample, "ukb", "non-ukb")
            subset += [subset_expr] if by_subset else []
        gen_anc += [meta_s.population_inference.pop] if by_ancestry else []

    ht = (
        ht.annotate(
            subset=subset,
            gen_anc=gen_anc,
            summary_stats=hl.zip(ht.summary_stats_meta, ht.summary_stats),
        )
        .select_globals()
        .explode("gen_anc")
        .explode("subset")
        .explode("summary_stats")
    )
    ht = ht.annotate(
        filter_group_meta=ht.summary_stats[0], summary_stats=ht.summary_stats[1]
    ).checkpoint(new_temp_file("agg_sample_stats.explode", "ht"))
    grouped_ht = ht.group_by(
        "subset", "gen_anc", "sex_chr_nonpar_group", "filter_group_meta"
    )
    ht1 = grouped_ht.aggregate(
        **{
            k: hl.struct(mean=hl.agg.mean(v), min=hl.agg.min(v), max=hl.agg.max(v))
            for k, v in ht.summary_stats.items()
        }
    ).checkpoint(new_temp_file("agg_sample_stats.first_agg", "ht"))
    ht2 = grouped_ht.aggregate(
        **{
            k: hl.agg.approx_quantiles(v, [0.0, 0.25, 0.5, 0.75, 1.0])
            for k, v in ht.summary_stats.items()
        }
    ).checkpoint(new_temp_file("agg_sample_stats.second_agg", "ht"))
    ht2_keyed = ht2[ht1.key]
    ht = ht1.annotate(
        **{k: ht1[k].annotate(quantiles=ht2_keyed[k]) for k in ht1.row_value}
    )

    return ht


def get_pipeline_resources(
    data_type: str,
    test: bool,
    test_gene: bool,
    autosomes_only: bool,
    sex_chr_only: bool,
    use_intermediate_mt_for_sample_counts: bool,
    by_ancestry: bool,
    by_subset: bool,
    overwrite: bool,
    custom_suffix: Optional[str] = None,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the per-sample stats pipeline.

    :param data_type: Data type of the dataset.
    :param test: Whether to gather all resources for testing.
    :param test_gene: Whether the test is being performed on specific gene(s).
    :param autosomes_only: Whether to gather resources for autosomes only.
    :param sex_chr_only: Whether to gather resources for sex chromosomes only.
    :param use_intermediate_mt_for_sample_counts: Whether to use an intermediate MT for
        computing per-sample counts.
    :param by_ancestry: Whether to return resource for aggregate stats stratified by
        ancestry.
    :param by_subset: Whether to return resource for aggregate stats stratified by
        subset.
    :param overwrite: Whether to overwrite resources if they exist.
    :param custom_suffix: Optional custom suffix to add to the resource names.
    :return: PipelineResourceCollection containing resources for all steps of the
        per-sample stats pipeline.
    """
    # Initialize pipeline resource collection.
    per_sample_pipeline = PipelineResourceCollection(
        pipeline_name="calculate_per_sample_stats",
        overwrite=overwrite,
    )

    res_base_args = {"test": test, "data_type": data_type, "suffix": custom_suffix}

    # Create resource collection for each step of the pipeline.
    create_filter_group = PipelineStepResourceCollection(
        "--create-filter-group-ht",
        input_resources={
            "v4 release HT": {"release_ht": release_sites(data_type=data_type)}
        },
        output_resources={
            "filter_groups_ht": get_summary_stats_filtering_groups(
                data_type=data_type, test=test
            )
        },
    )

    # Uses `test_gene` instead of just `test` because the filter groups HT will not
    # have the same partitions as the raw VDS, so unless the test is on a specific
    # gene(s), the test will not work as expected.
    filter_groups_ht = get_summary_stats_filtering_groups(data_type, test=test_gene)
    create_intermediate = PipelineStepResourceCollection(
        "--create-intermediate-mt-for-sample-counts",
        input_resources={
            "--create-filter-group-ht": {"filter_groups_ht": filter_groups_ht}
        },
        output_resources={
            "temp_intermediate_ht": TableResource(
                get_checkpoint_path(
                    "per_sample_summary_stats_intermediate"
                    + (".test" if test else "")
                    + (".autosomes" if autosomes_only else "")
                    + (".sex_chr" if sex_chr_only else "")
                )
            )
        },
    )
    create_per_sample_counts = PipelineStepResourceCollection(
        "--create-per-sample-counts-ht",
        pipeline_input_steps=(
            [create_intermediate] if use_intermediate_mt_for_sample_counts else []
        ),
        add_input_resources={
            "--create-filter-group-ht": {"filter_groups_ht": filter_groups_ht}
        },
        output_resources={
            "per_sample_ht": get_per_sample_counts(
                **res_base_args,
                autosomes=autosomes_only,
                sex_chr=sex_chr_only,
            )
        },
    )
    combine_stats = PipelineStepResourceCollection(
        "--combine-autosome-and-sex-chr-stats",
        input_resources={
            "--create-per-sample-counts-ht --autosomes-only-stats": {
                "autosomes_ht": get_per_sample_counts(**res_base_args, autosomes=True)
            },
            "--create-per-sample-counts-ht --sex-chr-only-stats": {
                "sex_chr_ht": get_per_sample_counts(**res_base_args, sex_chr=True)
            },
        },
        output_resources={"per_sample_ht": get_per_sample_counts(**res_base_args)},
    )
    aggregate_stats = PipelineStepResourceCollection(
        "--aggregate-sample-stats",
        pipeline_input_steps=[create_per_sample_counts],
        add_input_resources={
            "sample_qc/create_sample_qc_metadata_ht.py": {
                "meta_ht": meta(data_type=data_type)
            }
        },
        output_resources={
            "per_sample_agg_ht": get_per_sample_counts(
                **res_base_args,
                autosomes=autosomes_only,
                sex_chr=sex_chr_only,
                aggregated=True,
                by_ancestry=by_ancestry,
                by_subset=by_subset,
            )
        },
    )
    # Add all steps to the pipeline resource collection.
    per_sample_pipeline.add_steps(
        {
            "create_filter_group": create_filter_group,
            "create_intermediate": create_intermediate,
            "create_per_sample_counts": create_per_sample_counts,
            "combine_stats": combine_stats,
            "aggregate_stats": aggregate_stats,
        }
    )

    return per_sample_pipeline


def main(args):
    """Collect per-sample variant counts and aggregate sample statistics."""
    hl.init(
        log="/per_sample_stats",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    overwrite = args.overwrite
    data_type = args.data_type
    create_filter_group = args.create_filter_group_ht
    create_per_sample_counts = args.create_per_sample_counts_ht
    create_intermediate = args.create_intermediate_mt_for_sample_counts
    use_intermediate = args.use_intermediate_mt_for_sample_counts
    combine_chr_stats = args.combine_autosome_and_sex_chr_stats
    aggregate_stats = args.aggregate_sample_stats
    autosomes_only = args.autosomes_only_stats
    sex_chr_only = args.sex_chr_only_stats
    exomes = data_type == "exomes"
    rare_variants_afs = args.rare_variants_afs

    if autosomes_only and sex_chr_only:
        raise ValueError(
            "Cannot use --autosomes-only-stats and --sex-chr-only-stats at the same "
            "time."
        )
    if not exomes and args.by_subset:
        raise ValueError("Stratifying by subset is only working on exomes data type.")

    # Handle test arguments.
    test_dataset = args.test_dataset
    test_gene = args.test_gene
    test_difficult_partitions = args.test_difficult_partitions
    test_partitions = (
        list(range(args.test_n_partitions)) if args.test_n_partitions else None
    )
    test = test_partitions or test_dataset or test_gene
    filter_intervals = None
    if test:
        filter_intervals, test_partitions = process_test_args(
            data_type,
            test_dataset=test_dataset,
            test_gene=test_gene,
            test_partitions=test_partitions,
            test_difficult_partitions=test_difficult_partitions,
            create_filter_group=create_filter_group,
            create_intermediate=create_intermediate,
            create_per_sample_counts=create_per_sample_counts,
            use_intermediate=use_intermediate,
            sex_chr_only=sex_chr_only,
            autosomes_only=autosomes_only,
        )

    # Get per-sample stats pipeline resources.
    per_sample_stats_resources = get_pipeline_resources(
        data_type,
        test=test,
        test_gene=test_gene,
        autosomes_only=autosomes_only,
        sex_chr_only=sex_chr_only,
        use_intermediate_mt_for_sample_counts=use_intermediate,
        by_ancestry=args.by_ancestry,
        by_subset=args.by_subset,
        overwrite=overwrite,
        custom_suffix=args.custom_suffix,
    )

    # Define VDS load function and parameters.
    vds_load_params = {
        "test": test_dataset,
        "filter_partitions": test_partitions,
        "autosomes_only": autosomes_only,
        "sex_chr_only": sex_chr_only,
        "filter_intervals": filter_intervals,
        "annotate_meta": not autosomes_only,
    }

    logger.info("Starting per-sample stats pipeline for %s...", data_type)
    try:
        if create_filter_group:
            logger.info("Creating Table of filter groups for %s summary stats...")
            res = per_sample_stats_resources.create_filter_group
            res.check_resource_existence()

            ht = res.release_ht.ht()
            ht = ht._filter_partitions(test_partitions) if test_partitions else ht
            ht = hl.filter_intervals(ht, filter_intervals) if test_gene else ht

            get_summary_stats_filter_groups_ht(
                ht,
                capture_regions=exomes,
                vep_canonical=args.vep_canonical,
                vep_mane=args.vep_mane,
                rare_variants_afs=rare_variants_afs,
            ).write(res.filter_groups_ht.path, overwrite=overwrite)

        if create_intermediate:
            logger.info("Creating intermediate MatrixTable for per-sample counts...")
            res = per_sample_stats_resources.create_intermediate
            res.check_resource_existence()

            create_intermediate_mt_for_sample_counts(
                load_mt_for_sample_counts(
                    data_type,
                    res.filter_groups_ht.ht(),
                    adjust_ploidy=not autosomes_only,
                    **vds_load_params,
                )
            ).write(res.temp_intermediate_ht.path, overwrite=overwrite)

        if create_per_sample_counts:
            logger.info("Calculating per-sample variant statistics for %s...")
            res = per_sample_stats_resources.create_per_sample_counts
            res.check_resource_existence()

            if use_intermediate:
                logger.info("Using intermediate MatrixTable for per-sample counts...")
                ht = res.temp_intermediate_ht.ht()
                ht = ht._filter_partitions(test_partitions) if test_partitions else ht
                ht = hl.filter_intervals(ht, filter_intervals) if test_gene else ht
                ht = create_per_sample_counts_from_intermediate_ht(ht)
            else:
                ht = create_per_sample_counts_ht(
                    load_mt_for_sample_counts(
                        data_type,
                        res.filter_groups_ht.ht(),
                        adjust_ploidy=not autosomes_only,
                        **vds_load_params,
                    )
                )

            ht.write(res.per_sample_ht.path, overwrite=overwrite)

        if test and (combine_chr_stats or aggregate_stats):
            logger.warning(
                "Testing: Using whatever per-sample counts testing Table(s) was most "
                "recently created."
            )

        if combine_chr_stats:
            logger.info(
                "Combining per-sample stats for autosomes with per-sample stats for "
                "sex chromosomes..."
            )
            res = per_sample_stats_resources.combine_stats
            res.check_resource_existence()

            combine_autosome_and_sex_chr_stats(
                res.autosomes_ht.ht(), res.sex_chr_ht.ht()
            ).write(res.per_sample_ht.path, overwrite=overwrite)

        if aggregate_stats:
            logger.info("Computing aggregate sample statistics...")
            res = per_sample_stats_resources.aggregate_stats
            res.check_resource_existence()

            compute_agg_sample_stats(
                res.per_sample_ht.ht(),
                res.meta_ht.ht(),
                by_ancestry=args.by_ancestry,
                by_subset=args.by_subset,
            ).write(res.per_sample_agg_ht.path, overwrite=overwrite)
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
            "Filter Tables/VDS to only the PCSK9 and/or TRPC5 and ZFY genes for "
            "testing. Recommended for a quick test (if used with --test-dataset) that "
            "the full pipeline is working. This is not recommended for testing "
            "potential memory errors in --create-per-sample-counts-ht."
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
            "workers. It also requires shuffling, so a cluster with all workers is "
            "recommended."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--combine-autosome-and-sex-chr-stats",
        help=(
            "Combine per-sample counts for autosomes with per-sample counts for sex "
            "chromosomes."
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
        "--sex-chr-only-stats",
        help="Whether to restrict per-sample summary stats to sex chromosomes only.",
        action="store_true",
    )
    parser.add_argument(
        "--rare-variants-afs",
        type=float,
        default=SUM_STAT_FILTERS["max_af"],
        help="The allele frequency threshold to use for rare variants.",
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
