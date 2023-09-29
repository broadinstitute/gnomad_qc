"""
Script to generate the frequency data annotations across v4 exomes.

This script first splits the v4 VDS into multiples VDSs based which are then densified
and annotated with frequency data and histograms. The VDSs are then merged back together
in a hail Table. Next the script corrects for the high AB heterozygous GATK artifact in
existing annotations when given the AF threshold using a high AB het array. The script
then computes the inbreeding coefficient using the raw call stats. Finally, it computes
the filtering allele frequency and grpmax with the AB-adjusted frequencies.
"""
import argparse
import logging
from copy import deepcopy
from typing import Dict, List, Optional, Tuple, Union

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_downsamplings,
    bi_allelic_site_inbreeding_expr,
    build_freq_stratification_list,
    compute_freq_by_strata,
    faf_expr,
    generate_freq_group_membership_array,
    get_adj_expr,
    merge_freq_arrays,
    merge_histograms,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta, split_vds_by_strata
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SORT_ORDER

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_downsampling, get_freq, get_split_vds
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds, get_logging_path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency")
logger.setLevel(logging.INFO)

AGE_HISTS = [
    "age_hist_het",
    "age_hist_hom",
]
"""Age histograms to compute and keep on the frequency Table."""
QUAL_HISTS = [
    "gq_hist_all",
    "dp_hist_all",
    "gq_hist_alt",
    "dp_hist_alt",
    "ab_hist_alt",
]
"""Quality histograms to compute and keep on the frequency Table."""
FREQ_HIGH_AB_HET_ROW_FIELDS = [
    "high_ab_hets_by_group",
    "high_ab_het_adjusted_ab_hists",
    "high_ab_het_adjusted_age_hists",
]
"""
List of top level row and global annotations relating to the high allele balance
heterozygote correction that we want on the frequency HT before deciding on the AF
cutoff.
"""
FREQ_ROW_FIELDS = [
    "freq",
    "qual_hists",
    "raw_qual_hists",
    "age_hists",
]
"""
List of top level row and global annotations with no high allele balance heterozygote
correction that we want on the frequency HT.
"""
ALL_FREQ_ROW_FIELDS = FREQ_ROW_FIELDS + FREQ_HIGH_AB_HET_ROW_FIELDS
"""
List of final top level row and global annotations created from dense data that we
want on the frequency HT before deciding on the AF cutoff.
"""
FREQ_GLOBAL_FIELDS = [
    "downsamplings",
    "freq_meta",
    "age_distribution",
    "freq_index_dict",
    "freq_meta_sample_count",
]
"""
List of final global annotations created from dense data that we want on the frequency
HT before deciding on the AF cutoff.
"""
SUBSET_DICT = {"gnomad": 0, "non_ukb": 1}
"""
Dictionary for accessing the annotations with subset specific annotations such as
age_hists, popmax, and faf.
"""


def get_freq_resources(
    overwrite: bool = False,
    test: Optional[bool] = False,
    chrom: Optional[str] = None,
    ukb: bool = False,
    non_ukb: bool = False,
) -> PipelineResourceCollection:
    """
    Get frequency resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :param chrom: Chromosome used in freq calculations.
    :param ukb: Whether to get frequency resources for UKB subset.
    :param non_ukb: Whether to get frequency resources for non UKB subset.
    :return: Frequency resources.
    """
    freq_output_resources = []
    if ukb:
        freq_output_resources.append("ukb")
    if non_ukb:
        freq_output_resources.append("non_ukb")
    freq_pipeline = PipelineResourceCollection(
        pipeline_name="frequency",
        overwrite=overwrite,
    )
    write_split_vds_and_downsampling_ht = PipelineStepResourceCollection(
        "--write-split-vds-and-downsampling-ht",
        output_resources={
            "split_vds": get_split_vds(test=test),
            "ds_ht": get_downsampling(test=test),
        },
    )
    freq_output_resources = {
        f"{s}_freq_ht": get_freq(
            test=test,
            hom_alt_adjusted=False,
            chrom=chrom,
            intermediate_subset=s,
            finalized=False,
        )
        for s in freq_output_resources
    }
    if non_ukb:
        freq_output_resources["non_ukb_ds_ht"] = get_downsampling(
            test=test, subset="non_ukb"
        )
    run_freq_and_dense_annotations = PipelineStepResourceCollection(
        "--run-freq-and-dense-annotations",
        pipeline_input_steps=[write_split_vds_and_downsampling_ht],
        output_resources=freq_output_resources,
    )
    combine_freq = PipelineStepResourceCollection(
        "--combine-freq-hts",
        pipeline_input_steps=[run_freq_and_dense_annotations],
        output_resources={
            "freq_ht": get_freq(
                test=test,
                hom_alt_adjusted=False,
                chrom=chrom,
                finalized=False,
            )
        },
    )
    correct_for_high_ab_hets = PipelineStepResourceCollection(
        "--correct-for-high-ab-hets",
        pipeline_input_steps=[combine_freq],
        output_resources={
            "corrected_freq_ht": get_freq(
                test=test, hom_alt_adjusted=True, chrom=chrom, finalized=False
            )
        },
    )
    finalize_freq_ht = PipelineStepResourceCollection(
        "--finalize-freq-ht",
        pipeline_input_steps=[correct_for_high_ab_hets],
        output_resources={"final_freq_ht": get_freq(test=test, finalized=True)},
    )
    freq_pipeline.add_steps(
        {
            "write_split_vds_and_downsampling_ht": write_split_vds_and_downsampling_ht,
            "run_freq_and_dense_annotations": run_freq_and_dense_annotations,
            "combine_freq": combine_freq,
            "correct_for_high_ab_hets": correct_for_high_ab_hets,
            "finalize_freq_ht": finalize_freq_ht,
        }
    )
    return freq_pipeline


def get_vds_for_freq(
    use_test_dataset: hl.bool = False,
    test_gene: hl.bool = False,
    test_n_partitions: Optional[hl.int] = None,
    chrom: Optional[hl.int] = None,
) -> hl.vds.VariantDataset:
    """
    Prepare VDS for frequency calculation by filtering to release samples and only adding necessary annotations.

    :param use_test_dataset: Whether to use test dataset.
    :param test_gene: Whether to filter to DRD2 for testing purposes.
    :param test_n_partitions: Number of partitions to use for testing.
    :param chrom: Chromosome to filter to.
    :return: Hail VDS with only necessary annotations.
    """
    logger.info(
        "Reading the %s gnomAD v4 VDS...", "test" if use_test_dataset else "full"
    )
    if test_n_partitions:
        test_partitions = range(test_n_partitions)
    else:
        test_partitions = None

    test = use_test_dataset or test_gene or test_n_partitions

    vds = get_gnomad_v4_vds(
        test=use_test_dataset,
        release_only=True,
        filter_partitions=test_partitions,
        chrom=chrom,
        annotate_meta=True,
    )

    if test_gene:
        logger.info("Filtering to DRD2 in VDS for testing purposes...")
        test_interval = [
            hl.parse_locus_interval(
                "chr11:113409605-113475691", reference_genome="GRCh38"
            )
        ]
        vds = hl.vds.filter_intervals(vds, test_interval, split_reference_blocks=True)

    logger.info("Annotating VDS with only necessary sample metadata...")
    rmt = vds.reference_data
    vmt = vds.variant_data
    project_meta_expr = vmt.meta.project_meta
    vmt = vmt.select_cols(
        pop=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        fixed_homalt_model=project_meta_expr.fixed_homalt_model,
        gatk_version=project_meta_expr.gatk_version,
        age=project_meta_expr.age,
        ukb_sample=hl.if_else(project_meta_expr.ukb_sample, "ukb", "non_ukb"),
    )

    logger.info("Dropping excessively multi-allelic site at chr19:5787204...")
    vmt = vmt.filter_rows(
        vmt.locus != hl.parse_locus("chr19:5787204", reference_genome="GRCh38")
    )

    logger.info("Getting age distribution of all samples in the callset...")
    vmt = vmt.annotate_globals(
        age_distribution=vmt.aggregate_cols(hl.agg.hist(vmt.age, 30, 80, 10))
    )

    logger.info("Annotating non_ref hets pre-split...")
    vmt = vmt.annotate_entries(_het_non_ref=vmt.LGT.is_het_non_ref())

    logger.info("Selecting only required fields to reduce memory usage...")
    vmt = vmt.select_entries("LA", "LAD", "DP", "GQ", "LGT", "_het_non_ref")

    logger.info("Spltting mutliallelics in VDS...")
    vds = hl.vds.VariantDataset(rmt, vmt)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    return vds


def annotate_freq_index_dict(ht: hl.Table) -> hl.Table:
    """
    Create frequency index dictionary.

    The keys are the strata over which frequency aggregations where calculated and
    the values are the strata's index in the frequency array.

    :param ht: Input Table.
    :return: Table with 'freq_index_dict' global field.
    """
    logger.info("Making freq index dict...")
    # Add additional strata to the sort order, keeping group, i.e. adj, at the end.
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version"]
    ht = ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=ht.freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
        )
    )
    return ht


def filter_freq_arrays_for_non_ukb_subset(
    ht: hl.Table,
    items_to_filter: Union[List[str], Dict[str, List[str]]],
    keep: bool = True,
    combine_operator: str = "and",
    annotations: Union[List[str], Tuple[str]] = ("freq", "high_ab_hets_by_group"),
    remove_subset_from_meta: bool = False,
) -> hl.Table:
    """
    Filter frequency arrays by metadata.

    Filter 'annotations' and `freq_meta` array fields to only `items_to_filter` by
    using the 'freq_meta' array field values.

    If `remove_subset_from_meta` is True, update 'freq_meta' dicts by removing items
    with the key "subset" removed. If False, update 'freq_meta' dicts to include a
    "subset" key with "non_ukb" value.

    Also rename the 'downsamplings' global field to "non_ukb_downsamplings" so we can
    merge them later without losing the non-UKB downsampling information.

    :param ht: Input Table.
    :param items_to_filter: Items to filter by.
    :param keep: Whether to keep or remove items. Default is True.
    :param combine_operator: Operator ("and" or "or") to use when combining items in
        'items_to_filter'. Default is "and".
    :param annotations: Annotations in 'ht' to filter by `items_to_filter`.
    :param remove_subset_from_meta: Whether to remove the "subset" key from 'freq_meta'
        or add "subset" key with "non_ukb" value. Default is False.
    :return: Table with filtered 'annotations' and 'freq_meta' array fields.
    """
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            **{a: ht[a] for a in annotations},
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=items_to_filter,
        keep=keep,
        combine_operator=combine_operator,
    )
    if remove_subset_from_meta:
        logger.info("Dropping non_ukb subset key from freq_meta...")
        freq_meta = freq_meta.map(
            lambda d: hl.dict(d.items().filter(lambda x: x[0] != "subset"))
        )
    else:
        logger.info("Adding non_ukb subset key from freq_meta...")
        freq_meta = freq_meta.map(
            lambda d: hl.dict(d.items().append(("subset", "non_ukb")))
        )

    ht = ht.annotate(**{a: array_exprs[a] for a in annotations})
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )

    return ht


def high_ab_het(
    entry: hl.StructExpression, col: hl.StructExpression
) -> hl.Int32Expression:
    """
    Determine if a call is considered a high allele balance heterozygous call.

    High allele balance heterozygous calls were introduced in certain GATK versions.
    Track how many calls appear at each site to correct them to homozygous
    alternate calls downstream in frequency calculations and histograms.

    Assumes the following annotations are present in `entry` struct:
        - GT
        - adj
        - _high_ab_het_ref

    Assumes the following annotations are present in `col` struct:
        - fixed_homalt_model

    :param entry: Entry struct.
    :param col: Column struct.
    :return: 1 if high allele balance heterozygous call, else 0.
    """
    return hl.int(
        entry.GT.is_het_ref()
        & entry.adj
        & ~col.fixed_homalt_model
        & entry._high_ab_het_ref
    )


def generate_freq_ht(
    mt: hl.MatrixTable,
    ds_ht: hl.Table,
    meta_ht: hl.Table,
    non_ukb_ds_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Generate frequency Table.

    Assumes all necessary annotations are present:
        `mt` annotations:
            - GT
            - adj
            - _high_ab_het_ref
            - fixed_homalt_model
            - fixed_homalt_model

        `ds_ht` annotations:
            - downsampling
            - downsamplings
            - ds_pop_counts

        `meta_ht` annotations:
            - pop
            - sex_karyotype
            - gatk_version
            - age
            - ukb_sample

    :param mt: Input MatrixTable.
    :param ds_ht: Table with downsampling annotations.
    :param meta_ht: Table with sample metadata annotations.
    :param non_ukb_ds_ht: Optional Table with non-UKB downsampling annotations.
    :return: Hail Table with frequency annotations.
    """
    meta_ht = meta_ht.semi_join(mt.cols())
    additional_strata_expr = [
        {"gatk_version": meta_ht.gatk_version},
        {"gatk_version": meta_ht.gatk_version, "pop": meta_ht.pop},
    ]
    logger.info("Building frequency stratification list...")
    strata_expr = build_freq_stratification_list(
        sex_expr=meta_ht.sex_karyotype,
        pop_expr=meta_ht.pop,
        additional_strata_expr=additional_strata_expr,
        downsampling_expr=ds_ht[meta_ht.key].downsampling,
    )
    group_membership_ht = generate_freq_group_membership_array(
        meta_ht,
        strata_expr,
        downsamplings=hl.eval(ds_ht.downsamplings),
        ds_pop_counts=hl.eval(ds_ht.ds_pop_counts),
    )
    group_membership = group_membership_ht[mt.col_key].group_membership
    group_membership_globals = group_membership_ht.index_globals()
    if non_ukb_ds_ht is not None:
        logger.info("Building non-ukb downsampling stratification list...")
        non_ukb_strata_expr = build_freq_stratification_list(
            pop_expr=meta_ht.pop,
            downsampling_expr=non_ukb_ds_ht[meta_ht.key].downsampling,
        )
        non_ukb_strata_expr = [e for e in non_ukb_strata_expr if "downsampling" in e]

        logger.info("Generating non-ukb downsampling group_membership array....")
        non_ukb_group_membership_ht = generate_freq_group_membership_array(
            meta_ht,
            non_ukb_strata_expr,
            downsamplings=hl.eval(non_ukb_ds_ht.downsamplings),
            ds_pop_counts=hl.eval(non_ukb_ds_ht.ds_pop_counts),
        )
        # Remove the first entry because it is the index for the full subset.
        group_membership = group_membership.extend(
            non_ukb_group_membership_ht[mt.col_key].group_membership[1:]
        )
        non_ukb_globals = non_ukb_group_membership_ht.index_globals()
        non_ukb_globals = non_ukb_globals.annotate(
            freq_meta=non_ukb_globals.freq_meta.map(
                lambda d: hl.dict(d.items().append(("subset", "non_ukb")))
            )
        )
        # Remove the first two because they are adj and raw for the full subset.
        group_membership_globals = group_membership_globals.annotate(
            **{
                g: group_membership_globals[g].extend(non_ukb_globals[g][2:])
                for g in ["freq_meta", "freq_meta_sample_count"]
            },
            non_ukb_downsamplings=non_ukb_globals.downsamplings,
            non_ukb_ds_pop_counts=non_ukb_globals.ds_pop_counts,
        )
    else:
        group_membership_globals = group_membership_globals.annotate(
            non_ukb_downsamplings=hl.missing(hl.tarray(hl.tint)),
            non_ukb_ds_pop_counts=hl.missing(hl.tdict(hl.tstr, hl.tint)),
        )

    logger.info("Annotating frequencies and counting high AB het calls...")
    freq_ht = compute_freq_by_strata(
        mt.annotate_cols(group_membership=group_membership),
        entry_agg_funcs={"high_ab_hets_by_group": (high_ab_het, hl.agg.sum)},
        select_fields=["hists_fields"],
    )
    # Note: To use "multi_way_zip_join" need globals to be the same but an if_else based
    #  on strata doesn't work because hail looks for the annotation
    #  "non_ukb_downsamplings" regardless of the conditional value and throws an error
    #  if it doesn't exist.
    freq_ht = freq_ht.annotate_globals(**group_membership_globals)

    return freq_ht


def densify_and_prep_vds_for_freq(
    vds: hl.vds.VariantDataset,
    ab_cutoff: float = 0.9,
) -> hl.MatrixTable:
    """
    Densify VDS and select necessary annotations for frequency and histogram calculations.

    Select entry annotations required for downstream work. 'DP', 'GQ', and '_het_ab' are
    required for histograms. 'GT' and 'adj' are required for frequency calculations.
    '_high_ab_het_ref' is required for high AB call corrections in frequency and
    histogram annotations.

    Assumes all necessary annotations are present:
        - adj
        - _het_non_ref
        - GT
        - GQ
        - AD
        - DP
        - sex_karyotype

    :param vds: Input VDS.
    :param ab_cutoff: Allele balance cutoff to use for high AB het annotation.
    :return: Dense MatrixTable with only necessary entry annotations.
    """
    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Computing sex adjusted genotypes...")
    ab_expr = mt.AD[1] / mt.DP
    gt_expr = adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype)

    mt = mt.select_entries(
        "DP",
        "GQ",
        GT=gt_expr,
        adj=get_adj_expr(gt_expr, mt.GQ, mt.DP, mt.AD),
        _het_ab=ab_expr,
        _high_ab_het_ref=(ab_expr > ab_cutoff) & ~mt._het_non_ref,
    )
    return mt


def get_downsampling_ht(mt: hl.MatrixTable, non_ukb: bool = False) -> hl.Table:
    """
    Get Table with downsampling groups for all samples or the non-UKB subset.

    :param mt: Input MatrixTable.
    :param non_ukb: Whether to get downsampling groups for the non-UKB subset. Default
        is False.
    :return: Table with downsampling groups.
    """
    logger.info(
        "Determining downsampling groups for %s...",
        "the non-UKB subset" if non_ukb else "all samples",
    )
    meta_ht = mt.cols()
    downsamplings = DOWNSAMPLINGS["v4"]
    if non_ukb:
        downsamplings = downsamplings[:-1]
    ds_ht = annotate_downsamplings(meta_ht, downsamplings, pop_expr=meta_ht.pop)

    return ds_ht


def update_non_ukb_freq_ht(freq_ht: hl.Table) -> hl.Table:
    """
    Update non-UKB subset frequencies to be ready for combining with the frequencies of other samples.

    Duplicates frequency info for all groups except the non-UKB specific downsamplings
    and adds "subset" annotation to 'freq_meta' for one of the duplicates.

    This allows for the non-UKB subset frequencies to be merged with the frequencies of
    the other samples to provide full dataset frequencies while keeping the non-UKB
    subset specific frequency information.

    :param freq_ht: Non-UKB frequency Table.
    :return: Restructured non-UKB frequency Table.
    """
    annotations = ["freq", "high_ab_hets_by_group"]
    global_annotations = ["freq_meta", "freq_meta_sample_count"]

    logger.info("Filtering to non_ukb subset downsamplings...")
    non_ukb_ds_ht = filter_freq_arrays_for_non_ukb_subset(
        freq_ht,
        items_to_filter=["downsampling", "subset"],
    )

    # Filter to only non_ukb group, pop, and sex strata so can add subset-specific
    # freqs to main array.
    # NOTE: This is duplicated data here, but it's necessary to merge split vds strata
    #  properly and still retain the subset freq data.
    logger.info("Filtering to non_ukb subset strata...")
    non_ukb_ht = filter_freq_arrays_for_non_ukb_subset(
        freq_ht,
        items_to_filter=["downsampling", "gatk_version"],
        keep=False,
        combine_operator="or",
    )

    logger.info(
        "Filtering non_ukb freqs to main freq strata to be combined with ukb..."
    )
    freq_ht = filter_freq_arrays_for_non_ukb_subset(
        freq_ht,
        items_to_filter=["downsampling", "subset"],
        remove_subset_from_meta=True,
        keep=False,
    )

    # There are no overlap of strata groups here so can do basic flatmap on the
    # arrays.
    freqs = hl.array(
        [freq_ht.row_value, non_ukb_ht[freq_ht.key], non_ukb_ds_ht[freq_ht.key]]
    )
    freq_ht = freq_ht.annotate(
        **{a: hl.flatmap(lambda x: x[a], freqs) for a in annotations}
    )
    freq_globals = hl.array(
        [
            x.select_globals(*global_annotations).index_globals()
            for x in [freq_ht, non_ukb_ht, non_ukb_ds_ht]
        ]
    )
    freq_ht = freq_ht.annotate_globals(
        **{a: hl.flatmap(lambda x: x[a], freq_globals) for a in global_annotations}
    )

    return freq_ht


def mt_hists_fields(mt: hl.MatrixTable) -> hl.StructExpression:
    """
    Annotate quality metrics histograms and age histograms onto MatrixTable.

    :param mt: Input MatrixTable.
    :return: Struct with quality metrics histograms and age histograms.
    """
    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    high_ab_gt_expr = hl.if_else(high_ab_het(mt, mt) == 1, hl.call(1, 1), mt.GT)

    return hl.struct(
        **qual_hist_expr(
            gt_expr=mt.GT,
            gq_expr=mt.GQ,
            dp_expr=mt.DP,
            adj_expr=mt.adj,
            ab_expr=mt._het_ab,
            split_adj_and_raw=True,
        ),
        high_ab_het_adjusted_ab_hists=qual_hist_expr(
            gt_expr=high_ab_gt_expr,
            adj_expr=mt.adj,
            ab_expr=mt._het_ab,
        ),
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        high_ab_het_adjusted_age_hists=age_hists_expr(mt.adj, high_ab_gt_expr, mt.age),
    )


def combine_freq_hts(
    freq_hts: Dict[str, hl.Table],
    row_annotations: List[str],
    global_annotations: List[str],
    age_hists: List[str] = AGE_HISTS,
    qual_hists: List[str] = QUAL_HISTS,
) -> hl.Table:
    """
    Combine frequency HTs into a single HT.

    :param freq_hts: Dictionary of frequency HTs.
    :param row_annotations: List of annotations to put onto one hail Table.
    :param global_annotations: List of global annotations to put onto one hail Table.
    :param age_hists: List of age histogram annotations to merge.
    :param qual_hists: List of quality histogram annotations to merge.
    :return: HT with all freq_hts annotations.
    """
    ht_i = range(len(freq_hts))
    freq_ht = hl.Table.multi_way_zip_join(
        list(freq_hts.values()), "ann_array", "global_array"
    )

    # Combine freq arrays and high ab het counts by group arrays into single
    # annotations.
    logger.info(
        "Merging frequency arrays, metadata, and high ab het counts by group array..."
    )
    a_array = freq_ht.ann_array
    g_array = freq_ht.global_array
    comb_freq, comb_freq_meta, count_arrays_dict = merge_freq_arrays(
        farrays=[a_array[i].freq for i in ht_i],
        fmeta=[g_array[i].freq_meta for i in ht_i],
        count_arrays={
            "high_ab_hets": [a_array[i].high_ab_hets_by_group for i in ht_i],
            "freq_meta_sample_count": [g_array[i].freq_meta_sample_count for i in ht_i],
        },
    )

    logger.info("Annotating merged freq array and high ab hets...")
    freq_ht = freq_ht.annotate(
        freq=comb_freq,
        high_ab_hets_by_group=count_arrays_dict["high_ab_hets"],
    )

    # Merge all histograms into single annotations.
    logger.info("Merging all histograms...")
    hist_structs = {
        "qual_hists": qual_hists,
        "raw_qual_hists": qual_hists,
        "high_ab_het_adjusted_ab_hists": ["ab_hist_alt", "ab_hist_alt_adj"],
        "age_hists": age_hists,
        "high_ab_het_adjusted_age_hists": age_hists,
    }
    hists_expr = {
        hist_struct: hl.struct(
            **{
                h: merge_histograms(
                    [freq_ht.ann_array[i][hist_struct][h] for i in ht_i]
                )
                for h in hists
            }
        )
        for hist_struct, hists in hist_structs.items()
    }
    freq_ht = freq_ht.annotate(**hists_expr)

    freq_ht = freq_ht.annotate_globals(
        downsamplings={
            "gnomad": freq_ht.global_array[0].downsamplings,
            "non_ukb": freq_hts["non_ukb"].index_globals().non_ukb_downsamplings,
        },
        age_distribution=freq_ht.global_array[0].age_distribution,
        freq_meta=comb_freq_meta,
        freq_meta_sample_count=hl.eval(count_arrays_dict["freq_meta_sample_count"]),
    )
    freq_ht = annotate_freq_index_dict(freq_ht)
    logger.info("Setting Y metrics to NA for XX groups...")
    freq_ht = freq_ht.annotate(freq=set_female_y_metrics_to_na_expr(freq_ht))
    freq_ht = freq_ht.select(*row_annotations)
    freq_ht = freq_ht.select_globals(*global_annotations)

    logger.info("Combined frequency HT schema...")
    freq_ht.describe()

    return freq_ht


def correct_for_high_ab_hets(ht: hl.Table, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct for high allele balance heterozygous calls in call statistics and histograms.

    High allele balance GTs were being called heterozygous instead of homozygous
    alternate in GATK until version 4.1.4.1. This corrects for those calls by adjusting
    them to homozygote alternate within our call statistics and histograms when a
    site's allele frequency is greater than the passed `af_threshold`. Raw data is not
    adjusted.

    :param ht: Table with frequency and histogram annotations for correction as well as
        high AB het annotations.
    :param af_threshold: Allele frequency threshold for high AB adjustment. Default is
        0.01.
    :return: Table with corrected call statistics and histograms.
    """
    # Correct call statistics by passing through 'freq' and 'high_ab_hets_by_group'
    # arrays to adjust each annotation accordingly.
    call_stats_expr = hl.map(
        lambda f, g: hl.struct(
            AC=hl.int32(f.AC + g),
            AF=hl.if_else(f.AN > 0, (f.AC + g) / f.AN, hl.missing(hl.tfloat64)),
            AN=f.AN,
            homozygote_count=f.homozygote_count + g,
        ),
        ht.freq,
        ht.high_ab_hets_by_group,
    )

    # Add already computed 'ab_adjusted' histograms to the 'qual_hist_expr'.
    qual_hist_expr = {
        f"ab_adjusted_{x}": ht[x].annotate(
            ab_hist_alt=ht.high_ab_het_adjusted_ab_hists[
                f"ab_hist_alt{'' if x.startswith('raw') else '_adj'}"
            ]
        )
        for x in FREQ_ROW_FIELDS
        if "qual_hist" in x
    }

    # If a sites AF is greater than the 'af_threshold', add high AB het adjusted
    # annotations, otherwise use original annotations.
    no_ab_adjusted_expr = {f"ab_adjusted_{x}": ht[x] for x in FREQ_ROW_FIELDS}
    ht = ht.select(
        "high_ab_hets_by_group",
        *FREQ_ROW_FIELDS,
        **hl.if_else(
            ht.freq[0].AF > af_threshold,
            hl.struct(
                ab_adjusted_freq=call_stats_expr,
                **qual_hist_expr,
                ab_adjusted_age_hists=ht.high_ab_het_adjusted_age_hists,
            ),
            hl.struct(**no_ab_adjusted_expr),
        ),
    )
    ht = ht.annotate_globals(af_threshold_for_freq_adjustment=af_threshold)
    return ht


def generate_faf_grpmax(ht: hl.Table) -> hl.Table:
    """
    Compute filtering allele frequencies ('faf') and 'grpmax' with the AB-adjusted frequencies.

    :param ht: Hail Table containing 'freq', 'ab_adjusted_freq', 'high_ab_het'
        annotations.
    :return: Hail Table with 'faf' and 'grpmax' annotations.
    """
    logger.info(
        "Filtering frequencies to just 'non_ukb' subset entries for 'faf' "
        "calculations..."
    )
    # Remove the key "subset" from each dict in 'non_ukb_freq_meta' list so 'faf' will
    # pull the correct indices, i.e. `faf_expr` requires specific lengths to consider
    # each element in the freq list.
    non_ukb_ht = filter_freq_arrays_for_non_ukb_subset(
        ht,
        items_to_filter={"subset": ["non_ukb"]},
        combine_operator="or",
        annotations=("ab_adjusted_freq",),
        remove_subset_from_meta=True,
    )
    freq_metas = {
        "gnomad": (ht.ab_adjusted_freq, ht.index_globals().freq_meta),
        "non_ukb": (
            non_ukb_ht[ht.key].ab_adjusted_freq,
            non_ukb_ht.index_globals().freq_meta,
        ),
    }
    faf_exprs = []
    faf_meta_exprs = []
    grpmax_exprs = {}
    for dataset, (freq, meta) in freq_metas.items():
        faf, faf_meta = faf_expr(freq, meta, ht.locus, POPS_TO_REMOVE_FOR_POPMAX)
        grpmax = pop_max_expr(freq, meta, POPS_TO_REMOVE_FOR_POPMAX)
        grpmax = grpmax.annotate(
            gen_anc=grpmax.pop,
            faf95=faf[
                hl.literal(faf_meta).index(lambda y: y.values() == ["adj", grpmax.pop])
            ].faf95,
        ).drop("pop")
        # Add subset back to non_ukb faf meta.
        if dataset == "non_ukb":
            faf_meta = [{**x, **{"subset": "non_ukb"}} for x in faf_meta]
        faf_exprs.append(faf)
        faf_meta_exprs.append(faf_meta)
        grpmax_exprs[dataset] = grpmax

    logger.info("Annotating 'faf' and 'grpmax'...")
    ht = ht.annotate(
        faf=hl.flatten(faf_exprs),
        grpmax=hl.struct(**grpmax_exprs),
    )
    faf_meta_exprs = hl.flatten(faf_meta_exprs)
    ht = ht.annotate_globals(
        faf_meta=faf_meta_exprs,
        faf_index_dict=make_freq_index_dict_from_meta(faf_meta_exprs),
    )

    return ht


def compute_inbreeding_coeff(ht: hl.Table) -> hl.Table:
    """
    Compute inbreeding coefficient using raw call stats.

    :param ht: Hail Table containing 'freq' array with struct entries of 'AC', 'AN', and
        'homozygote_count'.
    :return: Hail Table with inbreeding coefficient annotation.
    """
    ht = ht.annotate(
        InbreedingCoeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1])
    )
    return ht


def create_final_freq_ht(ht: hl.Table) -> hl.Table:
    """
    Create final freq Table with only desired annotations.

    :param ht: Hail Table containing all annotations.
    :return: Hail Table with final annotations.
    """
    logger.info("Dropping gatk_version from freq ht array annotations...")
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            "freq": ht.ab_adjusted_freq,
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=["gatk_version"],
        keep=False,
    )

    ht = ht.select(
        freq=array_exprs["freq"],
        faf=ht.faf,
        grpmax=ht.grpmax,
        inbreeding_coeff=ht.InbreedingCoeff,
        histograms=hl.struct(
            **{
                "qual_hists": ht.ab_adjusted_qual_hists,
                "raw_qual_hists": ht.ab_adjusted_raw_qual_hists,
                "age_hists": ht.ab_adjusted_age_hists,
            },
        ),
    )

    def _replace_pop_key_w_gen_anc(
        meta: hl.expr.ArrayExpression,
    ) -> hl.expr.ArrayExpression:
        """
        Replace 'pop' key with 'gen_anc' key in meta globals.

        :param meta: Array of dicts.
        :return: Array of dicts with 'pop' key replaced with 'gen_anc' key.
        """
        return meta.map(
            lambda d: hl.dict(
                d.items().map(lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x))
            )
        )

    g = ht.index_globals()
    ht = ht.select_globals(
        downsamplings=g.downsamplings,
        freq_meta=_replace_pop_key_w_gen_anc(freq_meta),
        freq_index_dict=make_freq_index_dict_from_meta(freq_meta),
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
        faf_meta=_replace_pop_key_w_gen_anc(g.faf_meta),
        faf_index_dict=g.faf_index_dict,
        age_distribution=g.age_distribution,
    )

    return ht


def main(args):
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    overwrite = args.overwrite
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = use_test_dataset or test_n_partitions or test_gene
    chrom = args.chrom
    ab_cutoff = args.ab_cutoff
    af_threshold = args.af_threshold
    ukb_only = args.ukb_only
    non_ukb_only = args.non_ukb_only

    ukb = True
    non_ukb = True
    if ukb_only and non_ukb_only:
        raise ValueError("Cannot specify both --ukb-only and --non-ukb-only.")
    elif ukb_only:
        non_ukb = False
    elif non_ukb_only:
        ukb = False

    hl.init(
        log="/generate_frequency_data.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")
    resources = get_freq_resources(overwrite, test, chrom, ukb, non_ukb)

    try:
        if args.write_split_vds_and_downsampling_ht:
            logger.info(
                "Getting multi-allelic split VDS with adj and _het_AD entry"
                " annotations..."
            )
            res = resources.write_split_vds_and_downsampling_ht
            res.check_resource_existence()

            vds = get_vds_for_freq(
                use_test_dataset, test_gene, test_n_partitions, chrom
            )
            vds = vds.checkpoint(res.split_vds.path, overwrite=overwrite)

            ds_ht = get_downsampling_ht(vds.variant_data)
            ds_ht.write(res.ds_ht.path, overwrite=overwrite)

        if args.run_freq_and_dense_annotations:
            logger.info("Running dense dependent steps...")
            res = resources.run_freq_and_dense_annotations
            res.check_resource_existence()

            vds = res.split_vds.vds()
            ds_ht = res.ds_ht.ht()
            meta_ht = vds.variant_data.cols()

            logger.info(
                "Splitting VDS by ukb_sample annotation to reduce data size for"
                " densification..."
            )
            vds_dict = split_vds_by_strata(vds, strata_expr=vds.variant_data.ukb_sample)
            for strata, vds in vds_dict.items():
                if ukb_only and strata == "non_ukb" or non_ukb_only and strata == "ukb":
                    continue

                if strata == "non_ukb":
                    non_ukb_ds_ht = get_downsampling_ht(vds.variant_data, non_ukb=True)
                    non_ukb_ds_ht = non_ukb_ds_ht.checkpoint(
                        res.non_ukb_ds_ht.path, overwrite=overwrite
                    )
                else:
                    non_ukb_ds_ht = None

                logger.info("Generating frequency and histograms for %s VDS...", strata)
                mt = densify_and_prep_vds_for_freq(vds, ab_cutoff=ab_cutoff)
                mt = mt.annotate_rows(hists_fields=mt_hists_fields(mt))
                freq_ht = generate_freq_ht(
                    mt, ds_ht, meta_ht, non_ukb_ds_ht=non_ukb_ds_ht
                )
                freq_ht = freq_ht.transmute(**freq_ht.hists_fields)
                freq_ht.write(
                    getattr(res, f"{strata}_freq_ht").path, overwrite=overwrite
                )

        if args.combine_freq_hts:
            logger.info("Combining frequency HTs...")
            res = resources.combine_freq
            res.check_resource_existence()
            freq_ht = combine_freq_hts(
                {
                    "ukb": res.ukb_freq_ht.ht(),
                    "non_ukb": update_non_ukb_freq_ht(res.non_ukb_freq_ht.ht()),
                },
                ALL_FREQ_ROW_FIELDS,
                FREQ_GLOBAL_FIELDS,
            )
            freq_ht.write(res.freq_ht.path, overwrite=overwrite)

        if args.correct_for_high_ab_hets:
            logger.info(
                "Adjusting annotations impacted by high AB het -> hom alt adjustment..."
            )
            res = resources.correct_for_high_ab_hets
            res.check_resource_existence()
            ht = res.freq_ht.ht()

            logger.info("Correcting call stats, qual AB hists, and age hists...")
            ht = correct_for_high_ab_hets(ht, af_threshold=af_threshold)

            logger.info("Computing FAF & grpmax...")
            ht = generate_faf_grpmax(ht)

            logger.info("Calculating InbreedingCoeff...")
            ht = compute_inbreeding_coeff(ht)

            logger.info("High AB het corrected frequency HT schema...")
            ht.describe()

            logger.info("Writing corrected frequency Table...")
            ht.write(res.corrected_freq_ht.path, overwrite=overwrite)

        if args.finalize_freq_ht:
            logger.info("Writing final frequency Table...")
            res = resources.finalize_freq_ht
            res.check_resource_existence()
            ht = create_final_freq_ht(res.corrected_freq_ht.ht())

            logger.info("Final frequency HT schema...")
            ht.describe()
            ht.write(res.final_freq_ht.path, overwrite=overwrite)
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("frequency_data"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--use-test-dataset",
        help="Runs a test on the gnomad test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--test-gene",
        help="Runs a test on the DRD2 gene in the gnomad test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the VDS as input for testing purposes. Defaults"
            "to 2 if passed without a value."
        ),
        nargs="?",
        const=2,
        type=int,
    )
    parser.add_argument(
        "--chrom",
        help="If passed, script will only run on passed chromosome.",
        type=str,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--write-split-vds-and-downsampling-ht",
        help="Write split VDS and downsampling HT.",
        action="store_true",
    )
    parser.add_argument(
        "--run-freq-and-dense-annotations",
        help=(
            "Calculate frequencies, histograms, and high AB sites per sample grouping."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--ukb-only",
        help="Only run frequency and histogram calculations for UKB samples.",
        action="store_true",
    )
    parser.add_argument(
        "--non-ukb-only",
        help="Only run frequency and histogram calculations for non-UKB samples.",
        action="store_true",
    )
    parser.add_argument(
        "--combine-freq-hts",
        help=(
            "Combine frequency and histogram Tables for UKB and non-UKB samples into a"
            " single Table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--correct-for-high-ab-hets",
        help=(
            "Correct each frequency entry to account for homozygous alternate depletion"
            " present in GATK versions released prior to 4.1.4.1 and run chosen"
            " downstream annotations."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--ab-cutoff",
        help=(
            "Allele balance threshold to use when adjusting heterozygous calls to "
            "homozygous alternate calls at sites for samples that used GATK versions"
            " released prior to 4.1.4.1."
        ),
        type=float,
        default=0.9,
    )
    parser.add_argument(
        "--af-threshold",
        help=(
            "Threshold at which to adjust site group frequencies at sites for"
            " homozygous alternate depletion present in GATK versions released prior to"
            " 4.1.4.1."
        ),
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--finalize-freq-ht",
        help=(
            "Finalize frequency Table by dropping unnecessary fields and renaming"
            " remaining fields."
        ),
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
