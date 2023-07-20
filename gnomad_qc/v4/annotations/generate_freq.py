"""
Script to generate the frequency data annotations across v4 exomes.

This script first splits the v4 VDS into multiples VDSs based which are then densified
and annotated with frequency data and histograms. The VDSs are then merged back together
in a hail Table. Next the script corrects for the high AB heterzygous GATK artifact in
existing annotations when given the AF threshold using a high AB het array. The script
then computes the inbreeding coefficient using the raw call stats. Finally, it computes
the filtering allele frequency and grpmax with the AB-adjusted frequencies.
"""
import argparse
import itertools
import logging
from copy import deepcopy
from typing import Dict, List, Optional, Union

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    merge_freq_arrays,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
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
from gnomad_qc.v4.resources.meta import meta

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency")
logger.setLevel(logging.INFO)


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
    :param chr: Chromosome to filter to.
    :return: Hail VDS with only necessary annotations.
    """
    logger.info(
        "Reading the %s gnomAD v4 VDS...", "test" if use_test_dataset else "full"
    )
    if use_test_dataset:
        vds = get_gnomad_v4_vds(test=True)
    else:
        vds = hl.vds.read_vds("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds")

    meta_ht = meta.ht()
    vds = hl.vds.filter_samples(vds, meta_ht.filter(meta_ht.release))

    if test_gene or test_n_partitions or chrom:
        if test_gene:
            logger.info("Filtering to DRD2 in VDS for testing purposes...")
            test_interval = [
                hl.parse_locus_interval(
                    "chr11:113409605-113475691", reference_genome="GRCh38"
                )
            ]
            vds = hl.vds.filter_intervals(
                vds, test_interval, split_reference_blocks=True
            )

        elif test_n_partitions:
            logger.info(
                "Filtering to %s partitions for testing purposes...", test_n_partitions
            )
            test_vmt = vds.variant_data._filter_partitions(range(test_n_partitions))
            test_rmt = vds.reference_data._filter_partitions(range(test_n_partitions))
            vds = hl.vds.VariantDataset(test_rmt, test_vmt)
        else:
            logger.info("Filtering to chromosome %s...", chrom)
            vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    (
        variants,
        samples,
    ) = vds.variant_data.count()  # TODO: Remove this logger after tests?
    logger.info("VDS has %s variants in %s samples...", variants, samples)

    logger.info("Annotating VDS with only necessary sample metadata...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.annotate_cols(
            pop=meta_ht[vds.variant_data.col_key].population_inference.pop,
            sex_karyotype=meta_ht[
                vds.variant_data.col_key
            ].sex_imputation.sex_karyotype,
            fixed_homalt_model=meta_ht[
                vds.variant_data.col_key
            ].project_meta.fixed_homalt_model,
            gatk_version=meta_ht[vds.variant_data.col_key].project_meta.gatk_version,
            age=meta_ht[vds.variant_data.col_key].project_meta.age,
            sample_age_bin=get_sample_age_bin(
                meta_ht[vds.variant_data.col_key].project_meta.age
            ),
            ukb_sample=hl.if_else(
                meta_ht[vds.variant_data.col_key].project_meta.ukb_sample,
                "ukb",
                "non_ukb",
            ),
        ),
    )

    # Downsamplings are done outside the above function as we need to annotate
    # globals and rows
    logger.info("Annotating downsampling groups...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        annotate_downsamplings(vds.variant_data, pop_expr=vds.variant_data.pop),
    )

    logger.info("Annotating non_ref hets pre-split...")
    vds = annotate_non_ref_het(vds)

    logger.info("Selecting only required fields to reduce memory usage...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.select_entries("LA", "LAD", "DP", "GQ", "LGT", "_het_non_ref"),
    )

    return vds


def annotate_non_ref_het(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Annotate non-ref heterozygous calls to prevent incorrect call adjustment post-split.

    :param vds: Hail VDS to annotate non-ref het onto variant data.
    :return: Hail VDS with _non_ref_het annotation.
    """
    vds.variant_data = vds.variant_data.annotate_entries(
        _het_non_ref=vds.variant_data.LGT.is_het_non_ref()
    )
    return vds


def annotate_adj_and_select_fields(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Annotate adj, _het_ad, and select fields to reduce memory usage.

    :param vds: Hail VDS to annotate adj onto variant data.
    :return: Hail VDS with adj annotation.
    """
    rmt = vds.reference_data
    vmt = vds.variant_data

    rmt = rmt.annotate_entries(
        adj=(rmt.DP >= 10) & (rmt.GQ >= 20)
    )  # TODO: confirm this doesn't mess up haploids
    vmt = vmt.select_entries(
        "_het_non_ref",
        "DP",
        "GQ",
        "GT",
        adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD),
        _het_ad=vmt.AD[1],
    )
    return hl.vds.VariantDataset(rmt, vmt)


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
                    AF=hl.if_else(f.AN > 0, (f.AC + g) / f.AN, hl.missing(hl.tfloat64)),
                    AN=f.AN,
                    homozygote_count=f.homozygote_count + g,
                ),
                ht.freq,
                ht.high_ab_hets_by_group_membership,
            ),
            ht.freq,
        )
    )

    return ht


def create_high_ab_age_hists(ht: hl.Table, age_group_key="sample_age_bin") -> hl.Table:
    """
    Create age histograms of high ab counts to account for high AB hets becoming hom alts.

    :param ht: Hail Table containing age hists, AB annotation.
    :param age_group_key: Age group key to use for age histogram.
    :return: Hail Table
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
    ht = ht.annotate(
        age_high_ab_hist=hl.struct(
            bin_freq=hl.starmap(
                lambda x, y: ht.high_ab_hets_by_group_membership[y],
                age_bin_indices_no_edges,
            ),
            n_smaller=ht.high_ab_hets_by_group_membership[
                age_bins_indices_dict["n_smaller"]
            ],
            n_larger=ht.high_ab_hets_by_group_membership[
                age_bins_indices_dict["n_larger"]
            ],
        )
    )


def get_sample_age_bin(
    sample_age: hl.expr.Int32Expression,
    lower_bound: int = 30,
    upper_bound: int = 80,
    number_of_bins: int = 10,
) -> hl.expr.StringExpression:
    """
    Get the age bin for a sample.

    :param sample_age: Sample age.
    :return: Sample age bin.
    """
    bin_size = (upper_bound - lower_bound) / number_of_bins
    lower_bin = hl.int(
        hl.floor((sample_age - lower_bound) / bin_size) * bin_size + lower_bound
    )
    upper_bin = hl.int(lower_bin + bin_size)

    bin_label = hl.if_else(
        sample_age < lower_bound,
        "n_smaller",
        hl.if_else(
            sample_age >= upper_bound,
            "n_larger",
            hl.str(lower_bin) + "-" + hl.str(upper_bin),
        ),
    )

    return hl.or_missing(hl.is_defined(sample_age), bin_label)


def compute_age_hist(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute age histograms for each variant.

    :param mt: Input MT with age annotation.
    :return: MatrixTable with age histogram annotations.
    """
    mt = mt.annotate_rows(**age_hists_expr(mt.adj, mt.GT, mt.age))

    # Compute callset-wide age histogram global
    mt = mt.annotate_globals(
        age_distribution=mt.aggregate_cols(hl.agg.hist(mt.age, 30, 80, 10))
    )
    return mt


def correct_age_hists(ht: hl.Table) -> hl.Table:
    """
    Correct age histograms.

    Correct by subtracting age_high_ab_hists from age_hist_het and adding
    age_high_ab_hists to age_hist_hom to account for high AB hets becoming hom alts.

    :param ht: Hail Table containing age hists and hist of AB counts by age annotation.
    :return: Hail Table
    """
    return ht.annotate(
        age_hist_het=hl.struct(
            bin_freq=hl.map(
                lambda x, y: x.bin_freq - y.bin_freq,
                ht.age_hist_het.bin_freq,
                ht.age_high_ab_hist.bin_freq,
            ),
            n_smaller=ht.age_hist_het.n_smaller - ht.age_high_ab_hist.n_smaller,
            n_larger=ht.age_hist_het.n_larger - ht.age_high_ab_hist.n_larger,
        ),
        age_hist_hom=hl.struct(
            bin_freq=hl.map(
                lambda x, y: x.bin_freq + y.bin_freq,
                ht.age_hist_hom.bin_freq,
                ht.age_high_ab_hist.bin_freq,
            ),
            n_smaller=ht.age_hist_hom.n_smaller + ht.age_high_ab_hist.n_smaller,
            n_larger=ht.age_hist_hom.n_larger + ht.age_high_ab_hist.n_larger,
        ),
    )


def compute_qual_hists(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate quality metrics histograms.

    :param mt: Input MT.
    :return: MatrixTable with qual histogram annotations.
    """
    mt = mt.annotate_rows(
        qual_hists=qual_hist_expr(
            gt_expr=mt.GT,
            gq_expr=mt.GQ,
            dp_expr=mt.DP,
            adj_expr=mt.adj,
            ab_expr=mt._het_ad / mt.DP,
        )
    )
    mt = mt.annotate_rows(
        qual_hists=hl.Struct(
            **{
                i.replace("_adj", ""): mt.qual_hists[i]
                for i in mt.qual_hists
                if "_adj" in i
            }
        ),
        raw_qual_hists=hl.Struct(
            **{i: mt.qual_hists[i] for i in mt.qual_hists if "_adj" not in i}
        ),
    )
    return mt


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

    ht = ht.annotate(
        qual_hist=ht.qual_hist.annotate(
            ab_hist_alt=_correct_ab_hist_alt(ht.qual_hists.ab_hist_alt)
        ),
        raw_qual_hist=ht.raw_qual_hist.annotate(
            ab_hist_alt=_correct_ab_hist_alt(ht.raw_qual_hists.ab_hist_alt)
        ),
    )
    return ht


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


def annotate_downsamplings(
    mt: hl.MatrixTable,
    pop_expr: Optional[hl.expr.StringExpression] = None,
    downsamplings=DOWNSAMPLINGS["v4"],
) -> hl.MatrixTable:
    """
    Annotate downsampling groups.

    :param mt: Input MT.
    :param pop_expr: Optional expression for population group. When provided, population sample sizes are added as values to downsamplings.
    :return: MatrixTable with downsampling annotations.
    """
    if pop_expr is not None:
        # We need to add all pop counts to the downsamplings list
        mt = mt.annotate_cols(pop=pop_expr)
        pop_counts = mt.aggregate_cols(hl.agg.counter(mt.pop))
        downsamplings = [x for x in downsamplings if x <= sum(pop_counts.values())]
        downsamplings = sorted(set(downsamplings + list(pop_counts.values())))

    logger.info("Found %i downsamplings: %s", len(downsamplings), downsamplings)
    downsampling_ht = mt.cols()
    downsampling_ht = downsampling_ht.annotate(r=hl.rand_unif(0, 1))
    downsampling_ht = downsampling_ht.order_by(downsampling_ht.r)
    scan_expr = {"global_idx": hl.scan.count()}

    # Index by pop if pop_expr is provided. This was pulled out of the existing
    # annotate_freq code. We need to do this here so that we can annotate the
    # pop_idx field in the downsampling_expr.
    if pop_expr is not None:
        scan_expr["pop_idx"] = hl.scan.counter(downsampling_ht.pop).get(
            downsampling_ht.pop, 0
        )

    downsampling_ht = downsampling_ht.annotate(**scan_expr)
    downsampling_ht = downsampling_ht.key_by("s").select(*scan_expr)
    mt = mt.annotate_cols(downsampling=downsampling_ht[mt.s])
    mt = mt.annotate_globals(
        downsamplings=downsamplings,
        ds_pop_counts=pop_counts if pop_expr is not None else None,
    )

    return mt


def split_vds(
    vds: hl.vds.VariantDataset, split_annotation: hl.expr.Expression
) -> hl.ArrayExpression:
    """
    Split a VDS into a list of VDSs based on the split_annotation.

    :param vds: Input VDS.
    :param split_annotation: Annotation on VDS variant_data MT split on.
    :return: List of VDSs.
    """
    vds_list = []
    for split in vds.variant_data.aggregate_cols(
        hl.agg.collect_as_set(split_annotation)
    ):
        samples = vds.variant_data.filter_cols(split_annotation == split).s.collect()
        vds_list.append(hl.vds.filter_samples(vds, samples))
    return vds_list


# The function below is an updated copy of Tim's freq updates found here:
# https://github.com/broadinstitute/gnomad_methods/pull/537. I've merged the
# additional_strata_expr and added the downsampling_expr functionality
# that supports external downsampling annotations into this function
def annotate_freq_and_high_ab_hets(
    mt: hl.MatrixTable,
    sex_expr: Optional[hl.expr.StringExpression] = None,
    pop_expr: Optional[hl.expr.StringExpression] = None,
    additional_strata_expr: Optional[
        Union[
            List[Dict[str, hl.expr.StringExpression]],
            Dict[str, hl.expr.StringExpression],
        ]
    ] = None,
    downsamplings: Optional[List[int]] = None,
    ds_pop_counts: Optional[Dict[str, int]] = None,
    downsampling_expr: Optional[hl.expr.StructExpression] = None,
    ab_cutoff: float = 0.9,
) -> hl.Table:
    """
    Annotate `mt` with stratified allele frequencies.

    The output Matrix table will include:
        - row annotation `freq` containing the stratified allele frequencies
        - global annotation `freq_meta` with metadata
        - global annotation `freq_sample_count` with sample count information

    .. note::

        Currently this only supports bi-allelic sites.

        The input `mt` needs to have the following entry fields:
          - GT: a CallExpression containing the genotype
          - adj: a BooleanExpression containing whether the genotype is of high quality or not.

        All expressions arguments need to be expression on the input `mt`.

    .. rubric:: `freq` row annotation

    The `freq` row annotation is an Array of Struct, with each Struct containing the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data, and the metadata about these annotations is
    stored in the globals.

    .. rubric:: Global `freq_meta` metadata annotation

    The global annotation `freq_meta` is added to the input `mt`. It is a list of dict.
    Each element of the list contains metadata on a frequency stratification and the index in the list corresponds
    to the index of that frequency stratification in the `freq` row annotation.

    .. rubric:: Global `freq_sample_count` annotation

    The global annotation `freq_sample_count` is added to the input `mt`. This is a sample count per sample grouping
    defined in the `freq_meta` global annotation.

    .. rubric:: The `downsamplings` parameter

    If the `downsamplings` parameter is used, frequencies will be computed for all samples and by population
    (if `pop_expr` is specified) by downsampling the number of samples without replacement to each of the numbers specified in the
    `downsamplings` array, provided that there are enough samples in the dataset.
    In addition, if `pop_expr` is specified, a downsampling to each of the exact number of samples present in each population is added.
    Note that samples are randomly sampled only once, meaning that the lower downsamplings are subsets of the higher ones.

    .. rubric:: The `additional_strata_expr` parameter

    If the `additional_strata_expr` parameter is used, frequencies will be computed for each of the strata dictionaries across all
    values. For example, if `additional_strata_expr` is set to `[{'platform': mt.platform}, {'platform':mt.platform, 'pop': mt.pop},
    {'age_bin': mt.age_bin}]`, then frequencies will be computed for each of the values of `mt.platform`, each of the combined values
    of `mt.platform` and `mt.pop`, and each of the values of `mt.age_bin`.

    :param mt: Input MatrixTable
    :param sex_expr: When specified, frequencies are stratified by sex. If `pop_expr` is also specified, then a pop/sex stratifiction is added.
    :param pop_expr: When specified, frequencies are stratified by population. If `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param additional_strata_expr: When specified, frequencies are stratified by the given additional strata. This can e.g. be used to stratify by platform, platform-pop, platform-pop-sex.
    :param downsamplings: When specified, frequencies are computed by downsampling the data to the number of samples given in the list. Note that if `pop_expr` is specified, downsamplings by population is also computed.
    :param ds_pop_counts: When specified, frequencies are computed by downsampling the data to the number of samples per pop in the dict. The key is the population and the value is the number of samples.
    :param downsampling_expr: When specified, frequencies are computed using the downsampling indices in the provided StructExpression. Note that if `pop_idx` is specified within the struct, downsamplings by population is also computed.
    :param ab_cutoff: Allele balance cutoff for high AB het filtering. Default is 0.9.
    :return: MatrixTable with `freq` annotation
    """
    if downsampling_expr is not None and downsamplings is None:
        raise NotImplementedError(
            "annotate_freq requires `downsamplings` when using `downsampling_expr`"
        )

    if (
        downsampling_expr.get("pop_idx") is not None
        and pop_expr is None
        and ds_pop_counts is None
    ):
        raise NotImplementedError(
            "annotate_freq requires `pop_expr` and `ds_pop_counts` when using"
            " `downsampling_expr` with pop_idx"
        )

    if downsampling_expr is not None and downsampling_expr.get("global_idx") is None:
        raise NotImplementedError(
            "annotate_freq requires `downsampling_expr` with key 'global_idx'"
        )

    if downsampling_expr is not None and (
        downsampling_expr.get("pop_idx") is None and pop_expr is not None
    ):
        raise NotImplementedError(
            "annotate_freq requires `downsampling_expr` with key 'pop_idx' when using"
            " `pop_expr`"
        )

    if additional_strata_expr is None:
        additional_strata_expr = [{}]

    if isinstance(additional_strata_expr, dict):
        additional_strata_expr = [additional_strata_expr]

    _freq_meta_expr = hl.struct(
        **{k: v for d in additional_strata_expr for k, v in d.items()}
    )
    if sex_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(sex=sex_expr)
    if pop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(pop=pop_expr)
    if downsampling_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(downsampling=downsampling_expr)

    # Annotate cols with provided cuts
    mt = mt.annotate_cols(_freq_meta=_freq_meta_expr)

    # Get counters for sex, pop and if set additional strata
    cut_dict = {
        cut: hl.agg.filter(
            hl.is_defined(mt._freq_meta[cut]), hl.agg.counter(mt._freq_meta[cut])
        )
        for cut in mt._freq_meta
    }
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    sample_group_filters = []

    # Create downsamplings if needed
    if downsampling_expr is not None:
        # Create downsampled sample groups
        sample_group_filters.extend(
            [
                (
                    {"downsampling": str(ds), "pop": "global"},
                    mt.downsampling.global_idx < ds,
                )
                for ds in hl.eval(downsamplings)
            ]
        )
        sample_group_filters.extend(
            [
                (
                    {"downsampling": str(ds), "pop": pop},
                    (mt._freq_meta.downsampling.pop_idx < ds)
                    & (mt._freq_meta.pop == pop),
                )
                for ds in hl.eval(downsamplings)
                for pop, pop_count in hl.eval(ds_pop_counts).items()
                if (ds <= pop_count) & (pop is not None)
            ]
        )
    # Add additional groupings to strata, e.g. strata-pop, strata-sex, strata-pop-sex
    additional_strata_filters = []
    for additional_strata in additional_strata_expr:
        additional_strata_values = [
            cut_data.get(strata, {}) for strata in additional_strata
        ]
        additional_strata_combinations = itertools.product(*additional_strata_values)

        additional_strata_filters.extend(
            [
                (
                    {
                        strata: str(value)
                        for strata, value in zip(additional_strata, combination)
                    },
                    hl.all(
                        list(
                            mt._freq_meta[strata] == value
                            for strata, value in zip(additional_strata, combination)
                        )
                    ),
                )
                for combination in additional_strata_combinations
            ]
        )

    # Add all desired strata, starting with the full set and ending with
    # downsamplings (if any)
    sample_group_filters = (
        [({}, True)]
        + [({"pop": pop}, mt._freq_meta.pop == pop) for pop in cut_data.get("pop", {})]
        + [({"sex": sex}, mt._freq_meta.sex == sex) for sex in cut_data.get("sex", {})]
        + [
            (
                {"pop": pop, "sex": sex},
                (mt._freq_meta.sex == sex) & (mt._freq_meta.pop == pop),
            )
            for sex in cut_data.get("sex", {})
            for pop in cut_data.get("pop", {})
        ]
        + additional_strata_filters
        + sample_group_filters
    )

    freq_sample_count = mt.aggregate_cols(
        [hl.agg.count_where(x[1]) for x in sample_group_filters]
    )

    logger.info("number of filters: %i", len(sample_group_filters))
    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta and sample count information
    freq_meta_expr = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]
    freq_meta_expr.insert(1, {"group": "raw"})
    freq_sample_count.insert(1, freq_sample_count[0])
    mt = mt.annotate_globals(
        freq_meta=freq_meta_expr,
        freq_sample_count=freq_sample_count,
    )

    mt = mt.annotate_globals(
        sample_group_filters_range_array=hl.range(len(sample_group_filters))
    )

    ht = mt.localize_entries("entries", "cols")
    ht = ht.annotate_globals(
        indices_by_category=hl.range(hl.len(ht.cols)).aggregate(
            lambda i: hl.agg.array_agg(
                lambda j: hl.agg.filter(
                    ht.cols[i].group_membership[j], hl.agg.collect(hl.int(i))
                ),
                ht.sample_group_filters_range_array,
            ),
        )
    )

    ht = ht.annotate(
        adj_array=ht.entries.map(lambda e: e.adj),
        gt_array=ht.entries.map(lambda e: e.GT),
    )

    def needs_high_ab_het_fix(entry, col):
        return (
            ((entry._het_ad / entry.DP) > ab_cutoff)
            & entry.adj
            & ~col.fixed_homalt_model
            & ~entry._het_non_ref
        )  # Skip adjusting genotypes if sample originally had a het nonref genotype

    ht = ht.annotate(
        high_ab_hets_by_group_membership=ht.indices_by_category.map(
            lambda sample_indices: hl.len(
                sample_indices.filter(
                    lambda i: ht.gt_array[i].is_het_ref()
                    & needs_high_ab_het_fix(ht.entries[i], ht.cols[i])
                )
            )
        )
    )

    tmp_path = hl.utils.new_temp_file("just_gt_adj", "ht")
    ht = ht.select(
        "gt_array", "adj_array", "high_ab_hets_by_group_membership"
    ).checkpoint(tmp_path)

    freq_expr = ht.indices_by_category.map(
        lambda sample_indices: sample_indices.aggregate(
            lambda i: hl.agg.filter(
                ht.adj_array[i], hl.agg.call_stats(ht.gt_array[i], ht.alleles)
            ),
        )
    )

    freq_expr = (
        freq_expr[:1]
        .extend([ht.gt_array.aggregate(lambda gt: hl.agg.call_stats(gt, ht.alleles))])
        .extend(freq_expr[1:])
    )

    # Select non-ref allele (assumes bi-allelic)
    freq_expr = freq_expr.map(
        lambda cs: cs.annotate(
            AC=cs.AC[1],
            AF=cs.AF[
                1
            ],  # TODO: This is NA in case AC and AN are 0 -- should we set it to 0?
            homozygote_count=cs.homozygote_count[1],
        )
    )

    ht = ht.annotate(freq=freq_expr)
    return ht.select("freq", "high_ab_hets_by_group_membership")


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

    :param mt: Input MT.
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

    logger.info("Annotating frequencies and counting high AB het calls...")
    additional_strata_expr = [
        {"gatk_version": mt.gatk_version},
        {
            "gatk_version": mt.gatk_version,
            "pop": mt.pop,
        },
        {
            "sample_age_bin": mt.sample_age_bin
        },  # We need this to get high AB hets for age histogram correction, dont care about the frequency of it, should we drop it?
        {"ukb_sample": mt.ukb_sample},  # confirm this is how we want to name this
    ]

    freq_ht = annotate_freq_and_high_ab_hets(
        mt,
        sex_expr=mt.sex_karyotype,
        pop_expr=mt.pop,
        downsamplings=mt.downsamplings,
        downsampling_expr=mt.downsampling,
        ds_pop_counts=mt.ds_pop_counts,
        additional_strata_expr=additional_strata_expr,
        ab_cutoff=ab_cutoff,
    )

    # Note: Tim's approach doesnt account for 'raw'in the sample_group_filters when
    # calculcating high AB hets per group. I think this is because the high AB hets
    # are dependent on the adj annotation so he assumed we didnt need it. However,
    # to properly match the high AB hets to the correct group using the same index,
    # we need to account for 'raw' in the high_ab_hets array. The number of high AB
    # hets in the 'raw' group matches the number in the adj group, so we can just
    # add insert the 'raw' group into the high_ab_hets array at the first index using
    # the value of 'adj', which is the zero index in the high_ab_hets array. Yuck.
    logger.info("Inserting raw group into high_ab_hets array...")
    freq_ht = freq_ht.annotate(
        high_ab_hets_by_group_membership=hl.array(
            [freq_ht.high_ab_hets_by_group_membership[0]]
        )
        .append(freq_ht.high_ab_hets_by_group_membership[0])
        .extend(freq_ht.high_ab_hets_by_group_membership[1:])
    )

    logger.info("Making freq index dict...")
    # Add our additional strata to the sort order, keeping group, i.e. adj, at the end
    sort_order = deepcopy(SORT_ORDER)
    sort_order[-1:-1] = ["gatk_version", "ukb_sample", "sample_age_bin"]

    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=freq_ht.freq_meta,
            label_delimiter="_",
            sort_order=sort_order,
            # TODO: Check if we actually want to see age_bin, I dont thin we do
        )
    )
    logger.info("Setting Y metrics to NA for XX groups...")
    freq_ht = freq_ht.annotate(freq=set_female_y_metrics_to_na_expr(freq_ht))

    logger.info("Annotating quality metrics histograms...")
    mt = compute_qual_hists(mt)
    final_rows_anns.update(
        {
            "qual_hists": mt.rows()[freq_ht.key].qual_hists,
            "raw_qual_hists": mt.rows()[freq_ht.key].raw_qual_hists,
        }
    )
    logger.info("Computing age histograms for each variant...")
    mt = compute_age_hist(mt)  # globla age distribution is a global
    final_rows_anns.update(
        {
            "age_hist_het": mt.rows()[freq_ht.key].age_hist_het,
            "age_hist_hom": mt.rows()[freq_ht.key].age_hist_hom,
        }
    )
    final_globals_anns.update({"age_distribution": mt.index_globals().age_distribution})
    freq_ht = freq_ht.annotate(**final_rows_anns)
    freq_ht = freq_ht.annotate_globals(**final_globals_anns)
    freq_ht.describe()
    freq_ht = freq_ht.checkpoint(
        new_temp_file(f"freq_ht{idx}", extension="ht"),
        overwrite=args.overwrite,
    )

    return freq_ht


def merge_histograms(ht: hl.Table, indices: List[int]) -> hl.Table:
    """
    Merge histogram annotations.

    This function merges all split histogram annotations by
    summing the arrays in an element-wise fashion across the like histograms. Histograms
    that should be merged are named similarly with the index as a suffix. It keeps one
    bin_edge annotation but merges the bin_freq, n_smaller, and n_larger annotations by
    summing them.
    :param ht: Hail Table with histogram annotations.
    :param indices: List of indices to merge.
    :return: Hail Table with merged histogram annotations.
    """
    age_hists = ["age_hist_het", "age_hist_hom"]
    qual_hists = [
        "gq_hist_all",
        "dp_hist_all",
        "gq_hist_alt",
        "dp_hist_alt",
        "ab_hist_alt",
    ]
    hist_structs = {"qual_hists": qual_hists, "raw_qual_hists": qual_hists}

    def _hist_merge(arrays: List[hl.expr.StructExpression]):
        """
        Merge histograms.

        :param arrays: List of histogram structs to merge.
        :return: Merged histogram struct.
        """
        return hl.fold(
            lambda i, j: hl.struct(
                **{
                    "bin_edges": (
                        i.bin_edges
                    ),  # Bin edges are the same for all histograms
                    "bin_freq": hl.zip(i.bin_freq, j.bin_freq).map(
                        lambda x: x[0] + x[1]
                    ),
                    "n_smaller": i.n_smaller + j.n_smaller,
                    "n_larger": i.n_larger + j.n_larger,
                }
            ),
            arrays[0].select("bin_edges", "bin_freq", "n_smaller", "n_larger"),
            arrays[1:],
        )

    ht = ht.annotate(
        **{
            age_hist: _hist_merge([ht[f"{age_hist}_{idx}"] for idx in indices])
            for age_hist in age_hists
        },
        **{
            hist_struct: hl.struct(
                **{
                    hist: _hist_merge(
                        [ht[f"{hist_struct}_{idx}"][hist] for idx in indices]
                    )
                    for hist in hists
                }
            )
            for hist_struct, hists in hist_structs.items()
        },
    )
    ht = ht.annotate_globals(
        age_distribution=_hist_merge(
            [ht.index_globals()[f"age_distribution_{i}"] for i in indices]
        )
    )

    return ht


def combine_freq_hts(
    freq_hts: hl.DictExpression,
    row_annotations: List[str],
    globals_annotations: List[str],
) -> hl.Table:
    """
    Combine frequency HTs into a single HT.

    :param freq_hts: Dictionary of frequency HTs.
    :param row_annotations: List of annotations to put onto one hail Table.
    :param globals_annotations: List of global annotations to put onto one hail Table.
    :return: HT with all freq_hts annotations.
    """
    # Create new HT with just variants and downsamplings global annotation to join onto
    freq_ht = freq_hts[1].select().select_globals("downsamplings")

    # Annotate all hts' annotations with a idx suffix to the new HT. We don't remove
    # variants when splitting the VDS so each table has all rows.
    logger.info("Annotating frequency HT with split HT dense dependent annotations")
    freq_ht = freq_ht.annotate(
        **{
            f"{ann}_{idx}": freq_hts[idx][freq_ht.key][ann]
            for idx in freq_hts.keys()
            for ann in row_annotations
        }
    )
    freq_ht = freq_ht.annotate_globals(
        **{
            f"{ann}_{idx}": freq_hts[idx].index_globals()[ann]
            for idx in freq_hts.keys()
            for ann in globals_annotations
        }
    )

    # Combine freq arrays and high ab het counts by group arrays into single annotations
    logger.info(
        "Merging frequency arrays, metadata, and high ab het counts by group array..."
    )
    comb_freq, comb_freq_meta, comb_high_ab_hets = merge_freq_arrays(
        farrays=[freq_ht[f"freq_{idx}"] for idx in freq_hts.keys()],
        fmeta=[freq_ht[f"freq_meta_{idx}"] for idx in freq_hts.keys()],
        count_arrays=[
            freq_ht[f"high_ab_hets_by_group_membership_{idx}"]
            for idx in freq_hts.keys()
        ],
    )
    # TODO: We want to drop all _idx annotations here but keeping them around for now for testing
    # This could also live in the merge function, passed as a boolean parameter
    freq_ht = freq_ht.annotate(
        freq=comb_freq,
        high_ab_hets_by_group_membership=comb_high_ab_hets,
    )

    # Merge all histograms into single annotations
    logger.info("Merging all histograms...")
    freq_ht = merge_histograms(freq_ht, freq_hts.keys())

    freq_ht = freq_ht.annotate_globals(
        freq_meta=hl.eval(comb_freq_meta),
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=comb_freq_meta, label_delimiter="_"
        ),
    )

    return freq_ht


def main(args):  # noqa: D103
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

        logger.info("Getting VDS with fields required for splitting VDS...")
        vds = get_vds_for_freq(use_test_dataset, test_gene, test_n_partitions, chrom)

        logger.info("Spltting mutliallelics in VDS...")
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

        logger.info(
            "Computing adj and _het_AD as part of reducing fields to reduce memory"
            " usage during dense dependent steps..."
        )
        vds = annotate_adj_and_select_fields(vds)

        if args.split_vds_by_annotation:
            logger.info(
                "Splitting VDS by sample annotation to reduce data size for"
                " densification..."
            )
            vds_list = split_vds(vds, split_annotation=vds.variant_data.ukb_sample)
            freq_hts = {}

            # Lists of final top level row and global annotations created from dense
            # data that we want on the frequency HT before deciding on the AF cutoff.
            final_rows_anns = [
                "freq",
                "high_ab_hets_by_group_membership",
                "qual_hists",
                "raw_qual_hists",
                "age_hist_het",
                "age_hist_hom",
            ]
            final_globals_anns = ["freq_meta", "age_distribution"]

            # Make a dict of freq hts with idx as key which will be used when joining
            # annotations. This is more straight forward to me than using joins since
            # hail keeps the existing annotation on the left table and the right table
            # gets a suffix added to the same annotation. This way, every table's
            # annotation will end up with a suffix in the merge.
            for idx, vds in enumerate(vds_list, start=1):
                freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff, idx=idx)
                freq_hts[idx] = freq_ht

            freq_ht = combine_freq_hts(freq_hts, final_rows_anns, final_globals_anns)
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
        ht = compute_inbreeding_coeff(ht)

        # TODO: Drop Age bins from freq fields?
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
