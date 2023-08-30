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
from typing import Dict, List, Optional

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_downsamplings,
    annotate_freq,
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
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SORT_ORDER
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_freq
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
# TODO: Add documentation.
QUAL_HISTS = [
    "gq_hist_all",
    "dp_hist_all",
    "gq_hist_alt",
    "dp_hist_alt",
    "ab_hist_alt",
]
# TODO: Add documentation.
FREQ_HIGH_AB_HET_ROW_FIELDS = [
    "high_ab_hets_by_group",
    "high_ab_het_adjusted_ab_hists",
    "high_ab_het_adjusted_age_hists",
]
FREQ_ROW_FIELDS = [
    "freq",
    "qual_hists",
    "raw_qual_hists",
    "age_hists",
]
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
    "age_hist_index_dict",
]
"""
List of final global annotations created from dense data that we want on the frequency
HT before deciding on the AF cutoff.
"""

# Dictionary for accessing the annotations with subset specific annotations
# such as age_hists, popmax, and faf.
SUBSET_DICT = {"gnomad": 0, "non_ukb": 1}


def get_freq_resources(
    overwrite: bool = False, test: Optional[bool] = False, chrom: Optional[str] = None
) -> PipelineResourceCollection:
    """
    Get frequency resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :param chrom: Chromosome used in freq calculations.
    :return: Frequency resources.
    """
    freq_pipeline = PipelineResourceCollection(
        pipeline_name="frequency",
        overwrite=overwrite,
    )
    run_freq_and_dense_annotations = PipelineStepResourceCollection(
        "--run-freq-and-dense-annotations",
        output_resources={
            "freq_and_dense_annotations": get_freq(
                test=test, hom_alt_adjusted=False, chrom=chrom
            ),
        },
    )
    correct_for_high_ab_hets = PipelineStepResourceCollection(
        "--correct-for-high-ab-hets",
        pipeline_input_steps=[run_freq_and_dense_annotations],
        output_resources={
            "freq_ht": get_freq(test=test, hom_alt_adjusted=True, chrom=chrom),
        },
    )
    freq_pipeline.add_steps(
        {
            "run_freq_and_dense_annotations": run_freq_and_dense_annotations,
            "correct_for_high_ab_hets": correct_for_high_ab_hets,
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

    # TODO: Need to decide which method to use for adj annotation.
    # logger.info(
    #    "Computing adj and _het_AD as part of reducing fields to reduce memory"
    #    " usage during dense dependent steps..."
    # )
    # vds = annotate_adj_and_select_fields(vds)

    return vds


def is_haploid(
    locus_expr: hl.expr.LocusExpression,
    karyotype_expr: hl.expr.StringExpression,
) -> hl.expr.CallExpression:
    """
    Add description.

    :param locus_expr: Hail locus expression.
    :param karyotype_expr:
    """
    xy = karyotype_expr == "XY"
    xx = karyotype_expr == "XX"
    x_nonpar = locus_expr.in_x_nonpar()
    y_par = locus_expr.in_y_par()
    y_nonpar = locus_expr.in_y_nonpar()

    return hl.or_missing(~(xx & (y_par | y_nonpar)), xy & (x_nonpar | y_nonpar))


def annotate_adj_and_select_fields(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Annotate adj, _het_ad, and select fields to reduce memory usage.

    :param vds: Hail VDS to annotate adj onto variant data.
    :return: Hail VDS with adj annotation.
    """
    rmt = vds.reference_data
    vmt = vds.variant_data

    logger.info("Computing sex adjusted genotypes...")
    rmt_sex_expr = vmt.cols()[rmt.col_key].sex_karyotype
    rmt.filter_entries(
        (rmt.locus.in_y_par() | rmt.locus.in_y_nonpar()) & (rmt_sex_expr == "XX"),
        keep=False,
    )
    rmt = rmt.annotate_entries(
        adj=(rmt.GQ >= 20)
        & hl.if_else(
            ~rmt.locus.in_autosome() & is_haploid(rmt.locus, rmt_sex_expr),
            rmt.DP >= 5,
            rmt.DP >= 10,
        )
    )

    vmt_gt_expr = hl.if_else(
        vmt.locus.in_autosome(),
        vmt.GT,
        adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype),
    )
    ab_expr = vmt.AD[1] / vmt.DP
    # TODO: Add ab cutoff parameter if using this method.
    vmt = vmt.select_entries(
        "DP",
        "GQ",
        GT=vmt_gt_expr,
        adj=get_adj_expr(vmt_gt_expr, vmt.GQ, vmt.DP, vmt.AD),
        _het_ab=ab_expr,
        # Skip adjusting genotypes if sample originally had a het nonref genotype.
        _high_ab_het_ref=(ab_expr > 0.9) & ~vmt._het_non_ref,
    )

    return hl.vds.VariantDataset(rmt, vmt)


def annotate_freq_index_dict(ht: hl.Table) -> hl.Table:
    """
    Create frequency index dictionary.

    The keys are the strata over which frequency aggregations where calculated and
    the values are the strata's index in the frequency array.

    :param ht:
    :return:
    """
    logger.info("Making freq index dict...")
    # Add additional strata to the sort order, keeping group, i.e. adj, at the end.
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version", "ukb_sample"]
    ht = ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=ht.freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
        )
    )
    return ht


def generate_freq_and_hists_ht(
    vds: hl.vds.VariantDataset,
    ds_ht: hl.Table,
    split_strata: Optional[str] = None,
    ab_cutoff: float = 0.9,
) -> hl.Table:
    """
    Generate frequency and histogram annotations.

    Assumes all necessary annotations are present:
        - adj
        - _het_ad
        - _het_non_ref
        - GT
        - GQ
        - s
        - pop
        - sex_karyotype
        - fixed_homalt_model
        - gatk_version
        - age
        - ukb_sample

    :param vds: Input VDS.
    :param ds_ht: Table with downsampling annotations.
    :param split_strata: Strata used when splitting VDS. Defaults to None.
    :param ab_cutoff: Allele balance cutoff to use for high AB het annotation.
    :return: Hail Table with frequency and histogram annotations.
    """
    final_rows_anns = {}

    ht = vds.variant_data.cols()
    additional_strata_expr = [
        {"gatk_version": ht.gatk_version},
        {"gatk_version": ht.gatk_version, "pop": ht.pop},
        {"ukb_sample": ht.ukb_sample},
    ]
    logger.info("Building frequency stratification list...")
    strata_expr = build_freq_stratification_list(
        sex_expr=ht.sex_karyotype,
        pop_expr=ht.pop,
        additional_strata_expr=additional_strata_expr,
        downsampling_expr=ds_ht[ht.key].downsampling,
    )
    logger.info("Generating group_membership array....")
    ht = generate_freq_group_membership_array(
        ht,
        strata_expr,
        downsamplings=hl.eval(ds_ht.downsamplings),
        ds_pop_counts=hl.eval(ds_ht.ds_pop_counts),
    )

    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Computing sex adjusted genotypes...")
    ab_expr = mt.AD[1] / mt.DP
    gt_expr = adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype)
    # Select entries required for downstream work. DP, GQ, and _het_ab are
    # required for histograms. GT and adj are required for frequency calculations.
    # _high_ab_het_ref is required for high AB call corrections in frequency and
    # histogram annotaitons.
    mt = mt.select_entries(
        "DP",
        "GQ",
        GT=gt_expr,
        adj=get_adj_expr(gt_expr, mt.GQ, mt.DP, mt.AD),
        _het_ab=ab_expr,
        _high_ab_het_ref=(ab_expr > ab_cutoff) & ~mt._het_non_ref,
    )

    def _high_ab_het(entry, col):
        """
        Determine if a call is considered a high allele balance heterozygous call.

        High allele balance heterozygous calls were introduced in certain GATK versions.
        Track how many calls appear at each site to correct them to homozygous
        alternate calls downstream in frequency calculations and histograms.

        :param entry:
        :param col:
        :return:
        """
        return hl.int(
            entry.GT.is_het_ref()
            & entry.adj
            & ~col.fixed_homalt_model
            & entry._high_ab_het_ref
        )

    logger.info("Annotating frequencies and counting high AB het calls...")
    freq_ht = compute_freq_by_strata(
        mt.annotate_cols(group_membership=ht[mt.col_key].group_membership),
        entry_agg_funcs={"high_ab_hets_by_group": (_high_ab_het, hl.agg.sum)},
    )
    freq_ht = freq_ht.annotate_globals(**ht.index_globals())

    # TODO: Add in non-ukb specific downsamplings and freq calculations. Its another full pass
    #  of the data but this reduces the strata size by 132 groups if we do it here instead of in one go above
    # I'm guessing well hit the OOM error again if we try adding them all above the non-ukb ht freq
    # calculation.
    if split_strata == "non_ukb":
        logger.info("Downsampling the non_ukb VDS for separate freq calcuation...")
        non_ukb_ds_freq_ht = annotate_freq(
            mt,
            pop_expr=mt.pop,
            downsamplings=DOWNSAMPLINGS["v4"][:-1],
            annotate_mt=False,
            entry_agg_funcs={"high_ab_hets_by_group": (_high_ab_het, hl.agg.sum)},
        )

        # Filter freq array field to only downsampling indices by using the freq_meta
        # array field values.
        def _select_non_ukb_entries(arr_expr):
            return (
                hl.zip(arr_expr, non_ukb_ds_freq_ht.freq_meta)
                .filter(lambda i: i[1].keys().contains("downsampling"))
                .map(lambda k: k[0])
            )

        non_ukb_ds_freq_ht = non_ukb_ds_freq_ht.annotate(
            freq=_select_non_ukb_entries(non_ukb_ds_freq_ht.freq),
            high_ab_hets_by_group=_select_non_ukb_entries(
                non_ukb_ds_freq_ht.high_ab_hets_by_group
            ),
        )
        # Filter freq_meta array field to dict entries with a "downsampling" key and
        # update it to have "subset" key with "non_ukb" value. We don't want any other
        # strata that are returned from annotate_freq, e.g. 'pop'. Also rename the
        # downsamplings global field to "non_ukb_downsamplings" so we can merge two
        # later and not lose the non_ukb downsampling information.
        non_ukb_ds_freq_ht = non_ukb_ds_freq_ht.annotate_globals(
            freq_meta=non_ukb_ds_freq_ht.freq_meta.filter(
                lambda i: i.keys().contains("downsampling")
            ).map(
                lambda d: hl.dict(
                    hl.zip(d.keys(), d.values()).append(("subset", "non_ukb"))
                )
            ),
            non_ukb_downsamplings=non_ukb_ds_freq_ht.downsamplings,
            freq_meta_sample_count=_select_non_ukb_entries(
                non_ukb_ds_freq_ht.freq_meta_sample_count
            ),
        )
        # Etiher filter main freq_ht freq here or after checkpoint to have only
        # pop and sex groups, add subset key and then merge
        non_ukb_ds_freq_ht = non_ukb_ds_freq_ht.checkpoint(
            new_temp_file("freq_non_ukb_ds", extension="ht")
        )

        # Grab/duplicate non-ukb subset specific freq strata for ancestry group and pop
        # and add "subset" key so downstream merge of non-ukb and ukb will be summed
        # properly but also contain subset freq breakdown.
        non_ukb_freq_meta, non_ukb_exprs = filter_arrays_by_meta(
            freq_ht.freq_meta,
            {
                "freq": freq_ht.freq,
                "freq_meta_sample_count": (
                    freq_ht.index_globals().freq_meta_sample_count
                ),
                "high_ab_hets_by_group": freq_ht.high_ab_hets_by_group,
            },
            items_to_filter=["downsampling", "gatk_version", "ukb_sample"],
            keep=False,
            combine_operator="or",
        )
        non_ukb_freq_meta = non_ukb_freq_meta.map(
            lambda d: hl.dict(
                hl.zip(d.keys(), d.values()).append(("subset", "non_ukb"))
            )
        )

        # There are no overlap of strata groups here so can do basic extend on the
        # arrays
        freq_ht = freq_ht.annotate(
            freq=freq_ht.freq.extend(non_ukb_exprs["freq"]).extend(
                non_ukb_ds_freq_ht[freq_ht.key].freq
            ),
            high_ab_hets_by_group=freq_ht.high_ab_hets_by_group.extend(
                non_ukb_exprs["high_ab_hets_by_group"].extend(
                    non_ukb_ds_freq_ht[freq_ht.key].high_ab_hets_by_group
                )
            ),
        )
        non_ukb_globals = non_ukb_ds_freq_ht.index_globals()
        freq_ht = freq_ht.annotate_globals(
            freq_meta=freq_ht.freq_meta.extend(non_ukb_freq_meta).extend(
                non_ukb_globals.freq_meta
            ),
            freq_meta_sample_count=freq_ht.freq_meta_sample_count.extend(
                non_ukb_exprs["freq_meta_sample_count"]
            ).extend(non_ukb_globals.freq_meta_sample_count),
            non_ukb_downsamplings=non_ukb_globals.non_ukb_downsamplings,
        )
    # Note:To use "multi_way_zip_join" need globals to be the same but an if_else based
    # on strata doesnt work because hail looks for the annotation "non_ukb_downsamplings"
    # regardless of the conditional value and throws an error if it doesnt exist.
    else:
        freq_ht = freq_ht.annotate_globals(
            non_ukb_downsamplings=hl.missing(hl.tarray(hl.tint))
        )

    logger.info("Setting Y metrics to NA for XX groups...")
    freq_ht = annotate_freq_index_dict(freq_ht)
    freq_ht = freq_ht.annotate(freq=set_female_y_metrics_to_na_expr(freq_ht))

    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    high_ab_gt_expr = hl.if_else(_high_ab_het(mt, mt) == 1, hl.call(1, 1), mt.GT)
    mt = mt.select_rows(
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

    hists = mt.rows()[freq_ht.key]
    final_rows_anns.update({r: hists[r] for r in mt.row_value})

    freq_ht = freq_ht.annotate(**final_rows_anns)

    return freq_ht


# TODO: add to params.
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
    n_hts_range = range(len(freq_hts))
    freq_ht = hl.Table.multi_way_zip_join(
        list(freq_hts.values()), "ann_array", "global_array"
    )

    # Combine freq arrays and high ab het counts by group arrays into single
    # annotations.
    logger.info(
        "Merging frequency arrays, metadata, and high ab het counts by group array..."
    )
    comb_freq, comb_freq_meta, count_arrays_dict = merge_freq_arrays(
        farrays=[freq_ht.ann_array[i].freq for i in n_hts_range],
        fmeta=[freq_ht.global_array[i].freq_meta for i in n_hts_range],
        count_arrays={
            "high_ab_hets": [
                freq_ht.ann_array[i].high_ab_hets_by_group for i in n_hts_range
            ],
            "freq_meta_sample_count": [
                freq_ht.global_array[i].freq_meta_sample_count for i in n_hts_range
            ],
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
                    [freq_ht.ann_array[i][hist_struct][h] for i in n_hts_range]
                )
                for h in hists
            }
        )
        for hist_struct, hists in hist_structs.items()
    }
    freq_ht = freq_ht.annotate(**hists_expr)

    # Add in non-ukb subset age hists
    logger.info("Adding non_ukb subset's age histograms...")
    freq_ht = freq_ht.transmute(
        age_hists=hl.array(
            [freq_ht.age_hists, freq_hts["non_ukb"][freq_ht.key].age_hists]
        ),
        high_ab_het_adjusted_age_hists=hl.array(
            [
                freq_ht.high_ab_het_adjusted_age_hists,
                freq_hts["non_ukb"][freq_ht.key].high_ab_het_adjusted_age_hists,
            ]
        ),
    )
    freq_ht = freq_ht.annotate_globals(age_hist_index_dict=SUBSET_DICT)

    freq_ht = freq_ht.annotate_globals(
        downsamplings={
            "global": freq_ht.global_array[0].downsamplings,
            "non_ukb": freq_hts["non_ukb"].index_globals().non_ukb_downsamplings,
        },
        age_distribution=freq_ht.global_array[0].age_distribution,
        freq_meta=comb_freq_meta,
        freq_meta_sample_count=hl.eval(count_arrays_dict["freq_meta_sample_count"]),
    )
    freq_ht = annotate_freq_index_dict(freq_ht)
    freq_ht = freq_ht.select(*row_annotations)
    freq_ht = freq_ht.select_globals(*global_annotations)

    logger.info("Final frequency HT schema...")
    freq_ht.describe()

    return freq_ht


def correct_for_high_ab_hets(ht: hl.Table, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct for high allele balance heterozygous calls in call statistics and histograms.

    High allele balance GTs were being called heterozygous instead of homozygous
    alternate in GATK until version 4.1.4.1. This corrects for those calls by adjusting
    them  to homozygote alternate within our call statistics and histograms when a
    site's allele frequency is greater than the passed af_threshold. Raw data is not
    adjusted.

    :param ht: Table with frequency and histogram annotations for correction as well as high AB het annotations.
    :param af_threshold: Allele frequency threshold for high AB adjustment. Default is 0.01.
    :return:
    """
    # Correct call statistics by passing through freq and high_ab_hets_by_group arrays
    # to adjust each annotation accordingly
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

    # Add already computed ab_adjusted histograms to the qual_hist_expr
    qual_hist_expr = {
        f"ab_adjusted_{x}": ht[x].annotate(
            ab_hist_alt=ht.high_ab_het_adjusted_ab_hists[
                f"ab_hist_alt{'' if x.startswith('raw') else '_adj'}"
            ]
        )
        for x in FREQ_ROW_FIELDS
        if "qual_hist" in x
    }

    # If a sites AF is greater than the af_threshold, add high AB het adjusted annotations,
    # otherwise use original annotations.
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

    return ht


def generate_faf_grpmax(ht: hl.Table) -> hl.Table:
    """
    Compute filtering allele frequencies and grpmax with the AB-adjusted frequencies.

    :param ht: Hail Table containing freq, ab_adjusted_freq, high_ab_het annotations.
    :return: Hail Table with faf & grpmax annotations.
    """
    # TODO: Clean this up, lots of repetivie code, should be able to iterate
    # over some list with map
    logger.info(
        "Filtering frequencines to just non_ukb subset entries for faf calculations"
    )
    non_ukb_freq_meta, non_ukb_freq = filter_arrays_by_meta(
        ht.freq_meta,
        {"ab_adjusted_freq": ht.ab_adjusted_freq},
        items_to_filter={"subset": ["non_ukb"]},
        keep=True,
        combine_operator="or",
    )

    # Remove the key "subset" from each dict in non_ukb_freq_meta list so faf will pull
    # the correct indices, i.e. faf_expr requires specific lengths to consider each
    # element in the freq list
    logger.info("Dropping subset key for non_ukb for faf calculations...")
    non_ukb_freq_meta = non_ukb_freq_meta.map(
        lambda d: hl.dict(
            hl.zip(d.keys(), d.values()).filter(lambda x: x[0] != "subset")
        )
    )
    # Create list of tuples to iterate through for faf and grpmax calculations
    freqs_and_metas = [
        (ht.ab_adjusted_freq, ht.freq_meta),
        (non_ukb_freq["ab_adjusted_freq"], non_ukb_freq_meta),
    ]
    fafs = [
        faf_expr(x[0], x[1], ht.locus, POPS_TO_REMOVE_FOR_POPMAX)
        for x in freqs_and_metas
    ]

    logger.info("Annotating faf and grpmax...")
    ht = ht.annotate(
        faf=[faf[0] for faf in fafs],
        grpmax=[
            pop_max_expr(x[0], x[1], POPS_TO_REMOVE_FOR_POPMAX) for x in freqs_and_metas
        ],
    )
    ht = ht.annotate_globals(faf_meta=[faf[1] for faf in fafs])
    ht = ht.annotate_globals(
        faf_index_dict=[
            make_faf_index_dict(hl.eval(x), label_delimiter="-")
            for x in hl.eval(ht.faf_meta)
        ]
    )
    logger.info("Adding in faf95 for grpmax...")
    ht = ht.annotate(
        grpmax=hl.enumerate(ht.grpmax).map(
            lambda x: x[1].annotate(
                faf95=ht.faf[x[0]][
                    ht.faf_meta[x[0]].index(
                        lambda x: x.values() == ["adj", ht.grpmax[1].pop]
                    )
                ].faf95
            )
        )
    )
    return ht


def compute_inbreeding_coeff(ht: hl.Table) -> hl.Table:
    """
    Compute inbreeding coefficient using raw call stats.

    :param ht: Hail Table containing freq array with struct entries of AC, AN, and homozygote_count.
    :return: Hail Table with inbreeding coefficient annotation.
    """
    ht = ht.annotate(
        InbreedingCoeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1])
    )
    return ht


def main(args):
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = use_test_dataset or test_n_partitions or test_gene
    chrom = args.chrom
    ab_cutoff = args.ab_cutoff
    af_threshold = args.af_threshold

    hl.init(
        log="/generate_frequency_data.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    resources = get_freq_resources(args.overwrite, test, chrom)

    try:
        if args.run_freq_and_dense_annotations:
            logger.info("Running dense dependent steps...")
            res = resources.run_freq_and_dense_annotations
            res.check_resource_existence()

            logger.info(
                "Getting multi-allelic split VDS with adj and _het_AD entry"
                " annotations..."
            )
            vds = get_vds_for_freq(
                use_test_dataset, test_gene, test_n_partitions, chrom
            )

            logger.info("Determining downsampling groups...")
            ds_ht = (
                annotate_downsamplings(
                    vds.variant_data, DOWNSAMPLINGS["v4"], pop_expr=vds.variant_data.pop
                )
                .cols()
                .checkpoint(new_temp_file("downsamplings", extension="ht"))
            )

            if args.split_vds_by_annotation:
                logger.info(
                    "Splitting VDS by ukb_sample annotation to reduce data size for"
                    " densification..."
                )
                vds_dict = split_vds_by_strata(
                    vds, strata_expr=vds.variant_data.ukb_sample
                )
                freq_ht_dict = {}
                for strata, vds in vds_dict.items():
                    logger.info(
                        "Generating frequency and histograms for %s VDS...", strata
                    )
                    ht = generate_freq_and_hists_ht(
                        vds, ds_ht, split_strata=strata, ab_cutoff=ab_cutoff
                    )
                    # TODO: Change to use a fixed location in gnomad_tmp instead of
                    #  new_temp_file so we can easily rerun only failed ones if needed?
                    # TODO: Actually, do we want to parallelize this in some way?
                    ht = ht.checkpoint(new_temp_file(f"freq_{strata}", extension="ht"))
                    freq_ht_dict[strata] = ht

                freq_ht = combine_freq_hts(
                    freq_ht_dict,
                    ALL_FREQ_ROW_FIELDS,
                    FREQ_GLOBAL_FIELDS,
                )
            else:
                freq_ht = generate_freq_and_hists_ht(vds, ds_ht, ab_cutoff=ab_cutoff)

            freq_ht.write(res.freq_and_dense_annotations.path, overwrite=args.overwrite)

        if args.correct_for_high_ab_hets:
            logger.info(
                "Adjusting annotations impacted by high AB het -> hom alt adjustment..."
            )
            res = resources.correct_for_high_ab_hets
            res.check_resource_existence()
            ht = res.freq_and_dense_annotations.ht()

            logger.info(
                "Correcting call stats, qual AB histograms, and age histograms..."
            )
            ht = correct_for_high_ab_hets(ht, af_threshold=af_threshold)

            logger.info("Computing FAF & grpmax...")
            ht = generate_faf_grpmax(ht)

            logger.info("Calculating InbreedingCoeff...")
            ht = compute_inbreeding_coeff(ht)

            # TODO: I think we should add a finalize option that does what you describe
            #  below.
            # TODO: Leaving in know while we test but need to drop fields we do not want
            # -- 'age_high_ab_his', all annotations from the split VDSs, only keep combinged,
            # rename ab_adjusted_freq to just freq and decide if we want to store uncorrect?
            # Probably,just rename it, also remove age bins and gatk versions from freq fields?
            # Also, change captialization of hists depending on decision from DP slack
            # thread
            logger.info("Writing frequency table...")
            ht.describe()
            ht.write(res.freq_ht.path, overwrite=args.overwrite)
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
        "--run-freq-and-dense-annotations",
        help=(
            "Calculate frequencies, histograms, and high AB sites per sample grouping."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--split-vds-by-annotation",
        help=(
            "Split VDS by annotation to reduce data size for densification."
            " Defaults to False."
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
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
