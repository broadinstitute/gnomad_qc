"""Script to generate the frequency data annotations across v4 exomes."""  # TODO: Expand on script description once nearing completion.
import argparse
import logging
from typing import Optional, Union

import hail as hl
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq_and_high_ab_hets,
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
from gnomad_qc.v4.resources.annotations import (
    get_freq,  # TODO: Need to generate this resource for raw callstats; TODO: add the args used in this script to the resource on this new branch, use old branch as reference
)
from gnomad_qc.v4.resources.meta import meta

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency")
logger.setLevel(logging.INFO)


def get_freq_resources(
    overwrite: bool, test: bool, chrom: str
) -> PipelineResourceCollection:
    """
    Get frequency resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :return: Frequency resources.
    """
    freq_pipeline = PipelineResourceCollection(
        pipeline_name="frequency",
        overwrite=overwrite,
    )
    run_freq_and_dense_annotations = PipelineStepResourceCollection(
        "--run-freq-and-dense-annotations",
        output_resources={
            "run_freq_and_dense_annotations": get_freq(
                test=test, hom_alt_adjustment=False, chr=chrom
            ),
        },
    )
    correct_for_high_ab_hets = PipelineStepResourceCollection(
        "--correct-for-high-ab-hets",
        pipeline_input_steps=[run_freq_and_dense_annotations],
        output_resources={
            "freq_ht": get_freq(test=test, hom_alt_adjustment=True, chr=chrom),
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
    test_gene: hl.bool, test_n_partitions: hl.int, chrom: hl.int
) -> hl.vds.VariantDataset:
    """
    Prepare VDS for frequency calculation by filtering to release samples and only adding necessary annotations.

    :param test_gene: Whether to filter to DRD2 for testing purposes.
    :param test_n_partitions: Number of partitions to use for testing.
    :param chr: Chromosome to filter to.
    :return: Hail VDS with only necessary annotations.
    """
    logger.info("Reading gnomAD v4 VDS...")
    vds = vds = hl.vds.read_vds("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds")
    meta_ht = meta.ht()
    vds = hl.vds.filter_samples(vds, meta_ht.filter(meta_ht.release))

    if test_gene:
        logger.info("Filtering to DRD2 in VDS for testing purposes...")
        test_interval = [
            hl.parse_locus_interval(
                "chr11:113409605-113475691", reference_genome="GRCh38"
            )
        ]
        vds = hl.vds.filter_intervals(vds, test_interval, split_reference_blocks=True)
        variants, samples = vds.variant_data.count()
        logger.info(
            "Test VDS has %s variants in DRD2 in %s samples...", variants, samples
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

    logger.info("Annotating only necessary sample metadata...")
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

    # Downsamplings are down outside the above function as we need to store
    # global and row annotations
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
    mt = mt.annotate_rows(**age_hists_expr(mt.adj, mt.GT, mt.meta.project_meta.age))

    # Compute callset-wide age histogram global
    mt = mt.annotate_globals(
        age_distribution=mt.aggregate_cols(
            hl.agg.hist(mt.meta.project_meta.age, 30, 80, 10)
        )
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


def generate_faf_popmax(ht: hl.Table) -> hl.Table:
    """
    Compute filtering allele frequencies and popmax with the AB-adjusted frequencies.

    :param ht: Hail Table containing freq, ab_adjusted_freq, high_ab_het annotations.
    :return: Hail Table with faf & popmax annotations.
    """
    faf, faf_meta = faf_expr(
        ht.ab_adjusted_freq, ht.freq_meta, ht.locus, POPS_TO_REMOVE_FOR_POPMAX
    )
    ht = ht.annotate(
        faf=faf,
        popmax=pop_max_expr(  # TODO: Update popmax to grpmax
            ht.ab_adjusted_freq, ht.freq_meta, POPS_TO_REMOVE_FOR_POPMAX
        ),
    )
    ht = ht.annotate_globals(
        faf_meta=faf_meta,
        faf_index_dict=make_faf_index_dict(faf_meta, label_delimiter="-"),
    )
    ht = ht.annotate(
        popmax=ht.popmax.annotate(  # TODO: Update popmax to grpmax
            faf95=ht.faf[
                ht.faf_meta.index(lambda x: x.values() == ["adj", ht.popmax.pop])
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
        mt = mt.annotate_cols(pop=pop_expr)
        pop_sizes = list(mt.aggregate_cols(hl.agg.counter(mt.pop)).values())
        downsamplings = [x for x in downsamplings if x <= sum(pop_sizes)]
        downsamplings = sorted(set(downsamplings + pop_sizes))

    logger.info("Found %i downsamplings: %s", len(downsamplings), downsamplings)
    downsampling_ht = mt.cols()
    downsampling_ht = downsampling_ht.annotate(r=hl.rand_unif(0, 1))
    downsampling_ht = downsampling_ht.order_by(downsampling_ht.r)
    scan_expr = {"global_idx": hl.scan.count()}

    if pop_expr is not None:
        scan_expr["pop_idx"] = hl.scan.counter(downsampling_ht.pop).get(
            downsampling_ht.pop, 0
        )

    downsampling_ht = downsampling_ht.annotate(**scan_expr)
    downsampling_ht = downsampling_ht.key_by("s").select(*scan_expr)
    mt = mt.annotate_cols(downsampling=downsampling_ht[mt.s])
    mt = mt.annotate_globals(downsamplings=downsamplings)

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
        vds_list.append(
            hl.vds.VariantDataset(
                vds.reference_data,
                vds.variant_data.filter_cols(split_annotation == split),
            )
        )

    return vds_list


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
        {"sample_age_bin": mt.sample_age_bin},
        {"ukb_sample": mt.ukb_sample},  # confirm this is how we want to name this
    ]

    freq_ht = annotate_freq_and_high_ab_hets(  # TODO: Update this function to use _het_ad when passed or maybe create AB ann?
        mt,
        sex_expr=mt.sex_karyotype,
        pop_expr=mt.pop,
        downsamplings=mt.downsamplings,
        downsampling_expr=mt.downsampling,
        additional_strata_expr=additional_strata_expr,
        ab_cutoff=ab_cutoff,
    )

    logger.info("Making freq index dict...")
    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=hl.eval(freq_ht.freq_meta),
            label_delimiter="_",
            sort_order=SORT_ORDER.insert(
                -1, "gatk_version"
            ),  # TODO: Check if we actually want to see age_bin, I dont thin we do
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
    final_globals_anns.update({"age_distribution": mt[freq_ht.key].age_distribution})
    freq_ht = freq_ht.annotate(**final_rows_anns)
    freq_ht = freq_ht.annotate_globals(**final_globals_anns)
    freq_ht = freq_ht.checkpoint(
        new_temp_file(f"freq_ht{idx}", extension="ht"),
        overwrite=args.overwrite,
    )
    return freq_ht


def combine_freq_hts(freq_hts: hl.DictExpression, hists=[]) -> hl.Table:
    """
    Combine frequency HTs into a single HT.

    :param freq_hts: Dictionary of frequency HTs.
    :param hists: List of histograms to combine.
    :return: Combined frequency HT.
    """
    freq_ht = hl.fold(
        lambda x, y: x.union(y),
        freq_hts.values()[0],
        freq_hts.values()[1:],
    )
    freq, freq_meta = merge_freq_arrays(freq_ht.freq, freq_hts[1].freq_meta)
    freq_ht = freq_ht.annotate(
        freq=freq,
        high_ab_hets_by_group_membership=hl.agg.array_agg(  # TODO: This may need to be change
            lambda x: hl.agg.sum(x), freq_hts.values().map(lambda x: x.high_ab_hets)
        ),
    )
    freq_ht = freq_ht.annotate_globals(
        freq_meta=freq_meta,
        freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta=freq_meta, label_delimiter="_"
        ),
    )
    # TODO Merge all histograms here

    return freq_ht


def main(args):  # noqa: D103
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    test_gene = args.test_gene
    test_n_partitions = args.test_n_partitions
    test = test_gene or test_n_partitions
    chrom = args.chrom
    ab_cutoff = args.ab_cutoff
    af_threshold = args.af_threshold
    correct_for_high_ab_hets = args.correct_for_high_ab_hets

    hl.init(
        log="/generate_frequency_data.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    if args.run_dense_dependent_steps:
        logger.info("Running dense dependent steps...")
        resources = get_freq_resources(args.overwrite, test, chrom)
        res = resources.run_freq_and_dense_annotations

        logger.info("Getting VDS with fields required for splitting VDS...")
        vds = get_vds_for_freq(test_gene, test_n_partitions, chrom)

        logger.info("Spltting mutliallelics in VDS...")
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

        logger.info(
            "Computing adj and _het_AD as part of reducing fields to reduce memory"
            " usage during dense dependent steps..."
        )
        vds = annotate_adj_and_select_fields(vds)

        logger.info(
            "Splitting VDS by sample annotation to reduce data size for"
            " densification..."
        )
        if args.split_vds_by_annotation:
            vds_list = split_vds(vds, split_annotation=vds.variant_data.ukb_sample)
            freq_hts = {}
            for idx, vds in enumerate(vds_list, start=1):
                freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff, idx=idx)
                freq_hts[idx] = freq_ht
            freq_ht = combine_freq_hts(
                freq_hts,
                ["qual_hists", "raw_qual_hists", "age_hist_het", "age_hist_hom"],
            )
        else:
            freq_ht = generate_freq_and_hists_ht(vds, ab_cutoff=ab_cutoff)

        freq_ht.write(res.freq_high_ab_ht.path, overwrite=args.overwrite)

    if correct_for_high_ab_hets:
        logger.info(
            "Adjusting annotations impacted by high AB het -> hom alt adjustment..."
        )
        res = resources.correct_call_stats
        ht = res.freq_high_ab_ht.ht()

        logger.info("Correcting call stats...")
        ht = correct_call_stats(ht, af_threshold)

        logger.info("Correcting qual AB histograms...")
        ht = correct_qual_hists(ht)

        logger.info("Correcting age histograms...")
        ht = correct_age_hists(ht)

        logger.info("computing FAF & popmax...")
        ht = generate_faf_popmax(ht)

        logger.info("Calculating InbreedingCoeff...")
        ht = compute_inbreeding_coeff(ht)

        # TODO: Drop Age bins from freq fields
        logger.info("Writing frequency table...")
        ht.describe()
        ht.write(res.freq_ht.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test-dataset",
        help="Runs a test on two partitions of the MT.",
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
