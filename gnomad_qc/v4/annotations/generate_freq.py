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
from typing import List, Optional

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_downsamplings,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    merge_freq_arrays,
    merge_histograms,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.filtering import split_vds_by_strata
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
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

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
FREQ_ROW_FIELDS = [
    "freq",
    "high_ab_hets_by_group",
    "qual_hists",
    "raw_qual_hists",
    "age_hists",
]
"""
List of final top level row and global annotations created from dense data that we
want on the frequency HT before deciding on the AF cutoff.
"""

FREQ_GLOBAL_FIELDS = [
    "downsamplings",
    "freq_meta",
    "age_distribution",
    "freq_index_dict",
]
"""
List of final global annotations created from dense data that we want on the frequency
HT before deciding on the AF cutoff.
"""


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

    # TODO: Add comment and logger.
    vmt = vmt.annotate_globals(
        age_distribution=vmt.aggregate_cols(hl.agg.hist(vmt.age, 30, 80, 10))
    )

    # Downsamplings are done outside the above function as we need to annotate
    # globals and rows.
    logger.info("Annotating downsampling groups...")
    vmt = annotate_downsamplings(vmt, DOWNSAMPLINGS["v4"], pop_expr=vmt.pop)

    logger.info("Annotating non_ref hets pre-split...")
    vmt = vmt.annotate_entries(_het_non_ref=vmt.LGT.is_het_non_ref())

    logger.info("Selecting only required fields to reduce memory usage...")
    vmt = vmt.select_entries("LA", "LAD", "DP", "GQ", "LGT", "_het_non_ref")

    logger.info("Spltting mutliallelics in VDS...")
    vds = hl.vds.VariantDataset(rmt, vmt)
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

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
    vmt = vmt.select_entries(
        "_het_non_ref",
        "DP",
        "GQ",
        GT=vmt_gt_expr,
        adj=get_adj_expr(vmt_gt_expr, vmt.GQ, vmt.DP, vmt.AD),
        _het_ad=vmt.AD[1],
    )

    return hl.vds.VariantDataset(rmt, vmt)


def generate_freq_and_hists_ht(
    vds: hl.vds.VariantDataset,
    ab_cutoff: float = 0.9,
    idx: Optional[int] = None,
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
        - sample_age_bin
        - ukb_sample
        - downsampling
        - downsamplings

    :param vds: Input VDS.
    :param ab_cutoff: Allele balance cutoff to use for high AB het annotation.
    :param idx: Optional index to append to temp file name.
    :return: Hail Table with frequency and histogram annotations.
    """
    final_rows_anns = {}
    final_globals_anns = {}

    logger.info("Densifying VDS # %s...", idx)
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Computing sex adjusted genotypes...")
    mt = mt.transmute_entries(
        GT=adjusted_sex_ploidy_expr(mt.locus, mt.GT, mt.sex_karyotype),
    )
    mt = mt.annotate_entries(adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD))

    logger.info("Annotating frequencies and counting high AB het calls...")
    additional_strata_expr = [
        {"gatk_version": mt.gatk_version},
        {"gatk_version": mt.gatk_version, "pop": mt.pop},
        # TODO: confirm this is how we want to name this.
        {"ukb_sample": mt.ukb_sample},
    ]

    def _needs_high_ab_het_fix(entry, col):
        return hl.int(
            entry.GT.is_het_ref()
            # & (entry._het_ad / entry.DP > ab_cutoff)
            & (entry.AD[1] / entry.DP > ab_cutoff)
            & entry.adj
            & ~col.fixed_homalt_model
            & ~entry._het_non_ref
        )  # Skip adjusting genotypes if sample originally had a het nonref genotype.

    freq_ht = annotate_freq(
        mt,
        sex_expr=mt.sex_karyotype,
        pop_expr=mt.pop,
        downsamplings=hl.eval(mt.downsamplings),
        downsampling_expr=mt.downsampling,
        ds_pop_counts=hl.eval(mt.ds_pop_counts),
        additional_strata_expr=additional_strata_expr,
        entry_agg_funcs={"high_ab_hets_by_group": (_needs_high_ab_het_fix, hl.agg.sum)},
        annotate_mt=False,
    )

    logger.info("Making freq index dict...")
    # Add additional strata to the sort order, keeping group, i.e. adj, at the end.
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version", "ukb_sample"]

    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=freq_ht.freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
        )
    )
    logger.info("Setting Y metrics to NA for XX groups...")
    freq_ht = freq_ht.annotate(freq=set_female_y_metrics_to_na_expr(freq_ht))

    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    mt = mt.select_rows(
        **qual_hist_expr(
            gt_expr=mt.GT,
            gq_expr=mt.GQ,
            dp_expr=mt.DP,
            adj_expr=mt.adj,
            # ab_expr=mt._het_ad / mt.DP,
            ab_expr=mt.AD[1] / mt.DP,
            split_adj_and_raw=True,
        ),
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        high_ab_hets_corrected_age_hists=age_hists_expr(
            mt.adj,
            hl.if_else(_needs_high_ab_het_fix(mt, mt) == 1, hl.call(1, 1), mt.GT),
            mt.age,
        ),
    )

    hists = mt.rows()[freq_ht.key]
    final_rows_anns.update({r: hists[r] for r in mt.row_value})

    freq_ht = freq_ht.annotate(**final_rows_anns)
    freq_ht = freq_ht.annotate_globals(**final_globals_anns)

    return freq_ht


# TODO: add to params.
def combine_freq_hts(
    freq_hts: hl.DictExpression,
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
    :return: HT with all freq_hts annotations.
    """
    n_hts_range = range(len(freq_hts))
    freq_ht = hl.Table.multi_way_zip_join(freq_hts, "ann_array", "global_array")

    # Combine freq arrays and high ab het counts by group arrays into single
    # annotations.
    logger.info(
        "Merging frequency arrays, metadata, and high ab het counts by group array..."
    )
    comb_freq, comb_freq_meta, comb_high_ab_hets = merge_freq_arrays(
        farrays=[freq_ht.ann_array[i].freq for i in n_hts_range],
        fmeta=[freq_ht.global_array[i].freq_meta for i in n_hts_range],
        count_arrays=[freq_ht.ann_array[i].high_ab_hets_by_group for i in n_hts_range],
    )
    freq_ht = freq_ht.annotate(
        freq=comb_freq,
        high_ab_hets_by_group=comb_high_ab_hets,
    )

    # Merge all histograms into single annotations.
    logger.info("Merging all histograms...")
    hist_structs = {
        "qual_hists": qual_hists,
        "raw_qual_hists": qual_hists,
        "age_hists": age_hists,
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

    logger.info("Making freq index dict...")
    # Add our additional strata to the sort order, keeping group, i.e. adj, at the end.
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version", "ukb_sample"]
    freq_meta = hl.eval(comb_freq_meta)
    freq_ht = freq_ht.annotate_globals(
        downsamplings=freq_ht.global_array[0].downsamplings,
        age_distribution=freq_ht.global_array[0].age_distribution,
        freq_meta=freq_meta,
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
        ),
    )
    freq_ht = freq_ht.select(*row_annotations)
    freq_ht = freq_ht.select_globals(*global_annotations)

    logger.info("Final frequency HT schema...")
    freq_ht.describe()

    return freq_ht


# Functions to correct frequencies and hists for high ab hets.
def create_high_ab_age_hists_expr(ht: hl.Table, age_group_key="sample_age_bin"):
    """
    Create histograms of high ab counts using age bins to account for high AB hets becoming hom alts.

    :param ht: Hail Table containing age hists, AB annotation.
    :param age_group_key: Age group key to use for age histogram.
    :return: Hail struct containing age histogram of high ab counts.
    """
    non_range_entries = hl.set(["n_larger", "n_smaller"])
    age_bins_indices = hl.sorted(
        hl.enumerate(ht["freq_meta"], index_first=False)
        .filter(lambda x: x[0].contains(age_group_key))
        .map(lambda x: (x[0][age_group_key], x[1]))
    )

    age_bins_indices_dict = hl.dict(age_bins_indices)
    age_bin_indices_no_edges = age_bins_indices.filter(
        lambda x: ~non_range_entries.contains(x[0])
    )
    return hl.struct(
        bin_freq=hl.starmap(
            lambda x, y: ht.high_ab_hets_by_group[y],
            age_bin_indices_no_edges,
        ),
        n_smaller=ht.high_ab_hets_by_group[age_bins_indices_dict["n_smaller"]],
        n_larger=ht.high_ab_hets_by_group[age_bins_indices_dict["n_larger"]],
    )


def generate_faf_grpmax(ht: hl.Table) -> hl.Table:
    """
    Compute filtering allele frequencies and grpmax with the AB-adjusted frequencies.

    :param ht: Hail Table containing freq, ab_adjusted_freq, high_ab_het annotations.
    :return: Hail Table with faf & grpmax annotations.
    """
    faf, faf_meta = faf_expr(
        ht.ab_adjusted_freq, ht.freq_meta, ht.locus, POPS_TO_REMOVE_FOR_POPMAX
    )
    ht = ht.annotate(
        faf=faf,
        grpmax=pop_max_expr(
            ht.ab_adjusted_freq, ht.freq_meta, POPS_TO_REMOVE_FOR_POPMAX
        ),
    )
    ht = ht.annotate_globals(
        faf_meta=faf_meta,
        faf_index_dict=make_faf_index_dict(faf_meta, label_delimiter="-"),
    )
    ht = ht.annotate(
        grpmax=ht.grpmax.annotate(
            faf95=ht.faf[
                ht.faf_meta.index(lambda x: x.values() == ["adj", ht.grpmax.pop])
            ].faf95
        )
    )
    return ht


def correct_call_stats(ht: hl.Table, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct frequencies at sites with an AF greater than the af_threshold.

    :param ht: Hail Table containing freq and high_ab_het annotations.
    :param af_threshold: AF threshold at which to correct frequency. Default is 0.01.
    :return: Hail Table with adjusted frequencies.
    """
    ht = ht.annotate(
        ab_adjusted_freq=hl.if_else(
            ht.freq[0].AF > af_threshold,
            hl.map(
                lambda f, g: hl.struct(
                    AC=hl.int32(f.AC + g),
                    AN=f.AN,
                    homozygote_count=f.homozygote_count + g,
                    AF=hl.if_else(f.AN > 0, (f.AC + g) / f.AN, hl.missing(hl.tfloat64)),
                ),
                ht.freq,
                ht.high_ab_hets_by_group,
            ),
            ht.freq,
        )
    )

    return ht


def correct_qual_hists(ht: hl.Table) -> hl.Table:  # add ab_threshold as arg
    """
    Correct quality metrics histograms.

    Correct by accessing the qual_hist and raw_qual_hist structs and removing
    all counts from the ab_hist_alt array where bin_edges exceed 0.9 AB.

    :param ht: Hail Table containing qual hists, AB annotation.
    :return: Hail Table
    """

    def _correct_ab_hist_alt(ab_hist_alt):
        return hl.struct(
            bin_edges=ab_hist_alt.bin_edges,
            bin_freq=hl.map(
                lambda edge, freq: hl.if_else(edge >= 0.9, 0, freq),
                ab_hist_alt.bin_edges[:-1],
                ab_hist_alt.bin_freq,
            ),
            n_smaller=ab_hist_alt.n_smaller,
            n_larger=0,
        )

    qual_hists = ["qual_hists", "raw_qual_hists"]
    ht = ht.annotate(
        **{
            x: ht[x].annotate(ab_hist_alt=_correct_ab_hist_alt(ht[x].ab_hist_alt))
            for x in qual_hists
        }
    )
    return ht


def correct_age_hists(ht: hl.Table) -> hl.Table:
    """
    Correct age histograms.

    Correct by subtracting age_high_ab_hists from age_hist_het and adding
    age_high_ab_hists to age_hist_hom to account for high AB hets becoming hom alts.

    :param ht: Hail Table containing age hists and hist of AB counts by age annotation.
    :return: Hail Table
    """
    ht = ht.annotate(age_high_ab_hist=create_high_ab_age_hists_expr(ht))

    return ht.annotate(
        age_hist_het=hl.struct(
            bin_edges=ht.age_hist_het.bin_edges,
            bin_freq=hl.map(
                lambda x, y: x - y,
                ht.age_hist_het.bin_freq,
                ht.age_high_ab_hist.bin_freq,
            ),
            n_smaller=ht.age_hist_het.n_smaller - ht.age_high_ab_hist.n_smaller,
            n_larger=ht.age_hist_het.n_larger - ht.age_high_ab_hist.n_larger,
        ),
        age_hist_hom=hl.struct(
            bin_edges=ht.age_hist_hom.bin_edges,
            bin_freq=hl.map(
                lambda x, y: x + y,
                ht.age_hist_hom.bin_freq,
                ht.age_high_ab_hist.bin_freq,
            ),
            n_smaller=ht.age_hist_hom.n_smaller + ht.age_high_ab_hist.n_smaller,
            n_larger=ht.age_hist_hom.n_larger + ht.age_high_ab_hist.n_larger,
        ),
    )


# TODO: add automatic copy of log file.
def main(args):
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = use_test_dataset or test_n_partitions or test_gene
    chrom = args.chrom
    ab_cutoff = args.ab_cutoff
    af_threshold = args.af_threshold
    correct_for_high_ab_hets = args.correct_for_high_ab_hets

    hl.init(
        log="/generate_frequency_data.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    resources = get_freq_resources(args.overwrite, test, chrom)

    if args.run_freq_and_dense_annotations:
        logger.info("Running dense dependent steps...")
        res = resources.run_freq_and_dense_annotations
        res.check_resource_existence()

        logger.info(
            "Getting multi-allelic split VDS with adj and _het_AD entry annotations..."
        )
        vds = get_vds_for_freq(use_test_dataset, test_gene, test_n_partitions, chrom)
        if args.split_vds_by_annotation:
            logger.info(
                "Splitting VDS by ukb_sample annotation to reduce data size for"
                " densification..."
            )
            vds_list = split_vds_by_strata(vds, strata_expr=vds.variant_data.ukb_sample)
            freq_hts = []
            for idx, vds in enumerate(vds_list, start=1):
                freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff, idx=idx)
                # TODO: Change to use a fixed location in gnomad_tmp instead of
                #  new_temp_file so we can easily rerun only failed ones if needed?
                # TODO: Actually, do we want to parallelize this in some way?
                freq_hts.append(
                    freq_ht.checkpoint(
                        new_temp_file(f"freq_ht_{idx}", extension="ht"),
                        overwrite=args.overwrite,
                        _read_if_exists=False,
                    )
                )
            freq_ht = combine_freq_hts(freq_hts, FREQ_ROW_FIELDS, FREQ_GLOBAL_FIELDS)
        else:
            freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff)

        freq_ht.write(res.freq_and_dense_annotations.path, overwrite=args.overwrite)

    if correct_for_high_ab_hets:
        logger.info(
            "Adjusting annotations impacted by high AB het -> hom alt adjustment..."
        )
        res = resources.correct_for_high_ab_hets
        res.check_resource_existence()
        ht = res.freq_and_dense_annotations.ht()

        logger.info("Correcting call stats...")
        ht = correct_call_stats(ht, af_threshold)

        logger.info("Correcting qual AB histograms...")
        ht = correct_qual_hists(ht)

        logger.info("Correcting age histograms...")
        ht = correct_age_hists(ht)

        logger.info("computing FAF & grpmax...")
        ht = generate_faf_grpmax(ht)

        logger.info("Calculating InbreedingCoeff...")
        ht = ht.annotate(
            InbreedingCoeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1])
        )

        # TODO: Leaving in know while we test but need to drop fields we do not want
        # -- 'age_high_ab_his', all annotations from the split VDSs, only keep combinged,
        # rename ab_adjusted_freq to just freq and decide if we want to store uncorrect?
        # Probably,just rename it, also remove age bins and gatk versions from freq fields?
        # Also, change captialization of hists depending on decision from DP slack
        # thread
        logger.info("Writing frequency table...")
        ht.describe()
        ht.write(res.freq_ht.path, overwrite=args.overwrite)


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
