"""Updates the freq HT with the the correct AN."""

import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.resources.grch38.gnomad import release_all_sites_an
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    annotate_adj,
    build_freq_stratification_list,
    generate_freq_group_membership_array,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.sparse_mt import (
    compute_stats_per_ref_site,
    densify_all_reference_sites,
    get_allele_number_agg_func,
)

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.annotations.generate_freq import (
    compute_inbreeding_coeff,
    correct_for_high_ab_hets,
    create_final_freq_ht,
    generate_faf_grpmax,
)
from gnomad_qc.v4.resources.annotations import get_downsampling, get_freq
from gnomad_qc.v4.resources.basics import (
    calling_intervals,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.meta import meta

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_exomes_AN_fix")
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


def get_freq_an_resources(
    overwrite: bool = False,
    test: Optional[bool] = False,
    chrom: Optional[str] = None,
    calling_interval_name: Optional[str] = None,
    calling_interval_padding: Optional[int] = None,
) -> PipelineResourceCollection:
    """
    Get update frequency HT AN resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :param chrom: Chromosome used in freq calculations.
    :param calling_interval_name: Name of calling intervals to use.
    :param calling_interval_padding: Padding to use for calling intervals.

    :return: Update frequency HT AN resources.
    """
    update_freq_an_pipeline = PipelineResourceCollection(
        pipeline_name="frequency_an_fix",
        overwrite=overwrite,
    )

    regenerate_ref_an_and_gq_dp_hists = PipelineStepResourceCollection(
        "--regenerate-ref-an-and-hists",
        input_resources={
            "vds": get_gnomad_v4_vds(
                release_only=True,
                test=test,
                filter_partitions=range(2) if test else None,
                annotate_meta=True,
            ),
            "interval list": {
                "interval_ht": calling_intervals(
                    interval_name=calling_interval_name,
                    calling_interval_padding=calling_interval_padding,
                )
            },
            "metadata Table": {"meta_ht": meta()},
            "downsampling Tables": {
                "ds_ht": get_downsampling(),
                "non_ukb_ds_ht": get_downsampling(subset="non_ukb"),
            },
        },
        output_resources={
            "ref_an_ht": release_all_sites_an(public=False, test=test),
            "ref_gq_dp_hists": hl.utils.new_temp_file(
                "ref_hists", "ht"
            ),  # NOTE: Placeholder, until decide if this writes
        },
    )

    update_freq_an = PipelineStepResourceCollection(
        "--update-freq-an",
        pipeline_input_steps={
            "annotations/generate_freq.py --combine-freq-hts": {
                "freq_ht": get_freq(
                    version="4.0",
                    test=test,
                    hom_alt_adjusted=False,
                    chrom=chrom,
                    finalized=False,
                )
            },
            "annotations/compute_coverage.py --compute-allele-number-ht": {
                "an_ht": release_all_sites_an(public=False, test=test)
            },
        },
        output_resources={
            "updated_an_freq_ht": get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=False,
                chrom=chrom,
                finalized=False,
            )
        },
    )
    update_af_dependent_annotations = PipelineStepResourceCollection(
        "--update-af-dependent-annotations",
        pipeline_input_steps=[update_freq_an],
        output_resources={
            "updated_freq_ht": get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=True,
                chrom=chrom,
                finalized=False,
            )
        },
    )
    finalize_freq_ht = PipelineStepResourceCollection(
        "--finalize-freq-ht",
        pipeline_input_steps=[update_af_dependent_annotations],
        output_resources={
            "final_freq_ht": get_freq(version="4.1", test=test, finalized=True)
        },
    )

    update_freq_an_pipeline.add_steps(
        {
            "update_freq_an": update_freq_an,
            "update_af_dependent_annotations": update_af_dependent_annotations,
            "finalize_freq_ht": finalize_freq_ht,
        }
    )
    return update_freq_an_pipeline


def get_group_membership_ht(
    meta_ht: hl.Table,
    ds_ht: hl.Table,
    non_ukb_ds_ht: hl.Table,
) -> hl.Table:
    """
    Get group membership HT for all sites allele number stratification.

    :param meta_ht: Metadata HT.
    :param ds_ht: Full frequency downsampling HT.
    :param non_ukb_ds_ht: Non-UKB frequency downsampling HT.
    :return: Group membership HT.
    """
    # Filter to release samples.
    meta_ht = meta_ht.filter(meta_ht.release)

    # Filter to non-UKB samples.
    non_ukb_meta_ht = meta_ht.filter(~meta_ht.project_meta.ukb_sample)

    # Create group membership HT.
    ht = generate_freq_group_membership_array(
        meta_ht,
        build_freq_stratification_list(
            sex_expr=meta_ht.sex_imputation.sex_karyotype,
            pop_expr=meta_ht.population_inference.pop,
            downsampling_expr=ds_ht[meta_ht.key].downsampling,
        ),
        downsamplings=hl.eval(ds_ht.downsamplings),
        ds_pop_counts=hl.eval(ds_ht.ds_pop_counts),
    )

    # Add non-UKB group membership to HT.
    non_ukb_ht = generate_freq_group_membership_array(
        non_ukb_meta_ht,
        build_freq_stratification_list(
            sex_expr=non_ukb_meta_ht.sex_imputation.sex_karyotype,
            pop_expr=non_ukb_meta_ht.population_inference.pop,
            downsampling_expr=non_ukb_ds_ht[non_ukb_meta_ht.key].downsampling,
        ),
        downsamplings=hl.eval(non_ukb_ds_ht.downsamplings),
        ds_pop_counts=hl.eval(non_ukb_ds_ht.ds_pop_counts),
    )
    n_non_ukb = hl.eval(hl.len(non_ukb_ht.freq_meta))
    ht = ht.annotate(
        group_membership=ht.group_membership.extend(
            hl.coalesce(
                non_ukb_ht[ht.key].group_membership,
                hl.range(n_non_ukb).map(lambda x: False),
            )
        )
    )

    # Add non-UKB group metadata and sample counts to HT.
    non_ukb_globals = non_ukb_ht.index_globals()
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.extend(
            non_ukb_globals.freq_meta.map(
                lambda d: hl.dict(d.items().append(("subset", "non_ukb")))
            )
        ).map(
            lambda d: hl.dict(
                d.items().map(lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x))
            )
        ),
        freq_meta_sample_count=ht.freq_meta_sample_count.extend(
            non_ukb_globals.freq_meta_sample_count
        ),
    )

    return ht


def adjust_interval_padding(ht: hl.Table, padding: int) -> hl.Table:
    """
    Adjust interval padding in HT.

    .. warning::

        This function can lead to overlapping intervals, so it is not recommended for
        most applications. For example, it can be used to filter a variant list to all
        variants within the returned interval list, but would not work for getting an
        aggregate statistic for each interval if the desired output is independent
        statistics.

    :param ht: HT to adjust.
    :param padding: Padding to use.
    :return: HT with adjusted interval padding.
    """
    return ht.key_by(
        interval=hl.locus_interval(
            ht.interval.start.contig,
            ht.interval.start.position - padding,
            ht.interval.end.position + padding,
            reference_genome=ht.interval.start.dtype.reference_genome,
        )
    )


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


def vds_annotate_adj(vds: hl.vds.VariantDataset, freq_ht) -> hl.vds.VariantDataset:
    """
    Annotate adj, _het_ad, and select fields to reduce memory usage.

    :param vds: Hail VDS to annotate adj onto variant data.
    :return: Hail VDS with adj annotation.
    """
    # freq_ht = freq_ht.key_by("locus")
    freq_ht = hl.Table(
        hl.ir.TableKeyBy(
            freq_ht._tir, ["locus"], is_sorted=True
        )  # Prevents hail from running sort on HT which is already sorted.
    )

    rmt = vds.reference_data
    vmt = vds.variant_data

    logger.info("Computing sex adjusted genotypes...")
    rmt_sex_expr = vmt.cols()[rmt.col_key].meta.sex_imputation.sex_karyotype
    rmt.filter_entries(
        (rmt.locus.in_y_par() | rmt.locus.in_y_nonpar()) & (rmt_sex_expr == "XX"),
        keep=False,
    )
    rmt = rmt.annotate_entries(
        adj=(
            (rmt.GQ >= 20)
            & hl.if_else(
                ~rmt.locus.in_autosome() & is_haploid(rmt.locus, rmt_sex_expr),
                rmt.DP >= 5,
                rmt.DP >= 10,
            )
        )
    )

    vmt_gt_expr = hl.if_else(
        vmt.locus.in_autosome(),
        vmt.LGT,
        adjusted_sex_ploidy_expr(
            vmt.locus, vmt.LGT, vmt.meta.sex_imputation.sex_karyotype
        ),
    )
    vmt = vmt.annotate_entries(
        LGT=vmt_gt_expr,
        adj=(
            (vmt.GQ >= 20)
            & hl.if_else(
                vmt_gt_expr.is_haploid(),
                vmt.DP >= 5,
                vmt.DP >= 10,
            )
        ),
        adj_ab=hl.is_defined(vmt.LA[1]) & ((vmt.LAD[1] / vmt.DP) >= 0.2),
    )

    vmt = vmt.annotate_rows(in_freq=hl.is_defined(freq_ht[vmt.locus]))

    return hl.vds.VariantDataset(rmt, vmt)


def update_freq_an(freq_ht: hl.Table, an_ht: hl.Table) -> hl.Table:
    """
    Update frequency HT AN field with the correct AN from the AN HT.

    This module updates the freq_ht AN field which is an array of ints with the correct
    AN from the AN HT integer array annotation. It uses the an_ht dictionary and the
    freq_ht dictionary to match the indexing of array annotations so each group is
    correctly updated. It then recomputes AF, AC/AN.

    :param freq_ht: Frequency HT to update.
    :param an_ht: AN HT to update from.
    :return: Updated frequency HT.
    """
    freq_meta = freq_ht.index_globals().freq_meta
    an_meta = an_ht.index_globals().strata_meta

    # Set all freq_ht's ANs to 0 so can use the merge_freq function to update
    # the ANs (expects int64s)
    freq_ht = freq_ht.annotate(
        freq=freq_ht.freq.map(
            lambda x: x.select(
                AC=x.AC, AN=hl.int64(0), homozygote_count=x.homozygote_count, AF=x.AF
            )
        )
    )

    # Create freq array annotation of structs with AC and nhomozygotes set to 0 to
    # use merge_freq function (expects int64s)
    an_ht = an_ht.transmute(
        freq=an_ht.AN.map(
            lambda x: hl.struct(
                AC=0, AN=x, homozygote_count=0, AF=hl.missing(hl.tfloat64)
            )
        )
    )

    # Add an_ht to freq_ht to call merge_freq_arrays
    freq_ht = freq_ht.annotate(an_freq=an_ht[freq_ht.locus].freq)
    freq_ht = freq_ht.annotate_globals(an_strata_meta=an_meta)
    freq, freq_meta = merge_freq_arrays(
        [freq_ht.freq, freq_ht.an_freq], [freq_ht.freq_meta, freq_ht.an_strata_meta]
    )
    freq_ht = freq_ht.annotate(freq=freq)
    freq_ht = freq_ht.annotate_globals(freq_meta=freq_meta)
    return freq_ht


def drop_gatk_groupings(ht: hl.Table) -> hl.Table:
    """
    Drop GATK groupings from the frequency HT.

    :param ht: Frequency HT to drop GATK groupings from.
    :return: Frequency HT with GATK groupings dropped.
    """
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            "freq": ht.freq,
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=["gatk_version"],
        keep=False,
    )
    ht = ht.annotate(freq=array_exprs["freq"])
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )
    return ht


def compute_allele_number_per_ref_site_with_adj(
    vds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    interval_ht: hl.Table,
    group_membership_ht: hl.Table,
) -> hl.Table:
    """
    Compute the allele number per reference site including per allele adj.

    :param vds: Input VariantDataset.
    :param reference_ht: Table of reference sites.
    :return: Table of allele number per reference site.
    """

    def _transform_adj(t: hl.MatrixTable) -> hl.expr.BooleanExpression:
        return hl.or_else(t.adj, True)

    def _transform_fail_adj_an(t: hl.MatrixTable) -> hl.expr.Expression:
        return hl.or_missing(~t.adj_ab, [t.LA[1], t.LGT.ploidy])

    def _get_fail_adj_an(adj_expr: hl.expr.CallExpression) -> hl.expr.Expression:
        # Get the source Table for the CallExpression to grab alleles.
        t = adj_expr._indices.source
        fail_adj_an_expr = hl.or_missing(
            t.in_freq, hl.agg.group_by(adj_expr[0], hl.agg.sum(adj_expr[1]))
        )

        return fail_adj_an_expr

    def _get_hists(qual_expr) -> hl.expr.Expression:
        t = qual_expr._indices.source
        return hl.or_missing(
            t.in_freq,
            qual_hist_expr(
                gq_expr=qual_expr[0],
                dp_expr=qual_expr[1],
                adj_expr=qual_expr[2] == 1,
                split_adj_and_raw=True,
            ),
        )

    entry_agg_funcs = {
        "AN_fail_adj": (_transform_fail_adj_an, _get_fail_adj_an),
        "AN": get_allele_number_agg_func("LGT"),
        "qual_hists": (lambda t: [t.GQ, t.DP, t.adj], _get_hists),
    }

    ht = compute_stats_per_ref_site(
        vds,
        reference_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["LA", "adj_ab", "GQ", "DP"],
        row_keep_fields=["in_freq"],
        entry_agg_group_membership={"qual_hists": [{"group": "raw"}]},
    )

    ht = ht.checkpoint(
        "gs://gnomad-tmp/julia/test_all_sites_an_and_qual_hists.intermediate.rerun_2_2_24.ht",
        # _read_if_exists=True,
        overwrite=True,
    )
    vmt = vds.variant_data
    vmt = vmt.select_entries("LA", "LGT", "LAD", "DP", "GQ")
    vmt = hl.experimental.sparse_split_multi(vmt)
    # vmt = annotate_adj(vmt)
    # TODO: adjust sex ploidy
    vmt = vmt.annotate_entries(
        adj_ab=(
            hl.case()
            .when(~vmt.GT.is_het(), True)
            .when(vmt.GT.is_het_ref(), vmt.AD[vmt.GT[1]] / vmt.DP >= 0.2)
            .default(
                (vmt.AD[vmt.GT[0]] / vmt.DP >= 0.2)
                & (vmt.AD[vmt.GT[1]] / vmt.DP >= 0.2)
            )
        )
    )
    vmt = vmt.filter_entries(~vmt.adj_ab)
    vmt_hists_ht = (
        vmt.annotate_rows(
            **qual_hist_expr(
                gq_expr=vmt.GQ,
                dp_expr=vmt.DP,
            )
        )
        .rows()
        .select("gq_hist_all", "dp_hist_all")
    )

    freq_correction_ht = vds.variant_data.rows()
    freq_correction_ht = freq_correction_ht.annotate_globals(**ht.index_globals())
    freq_correction = ht[freq_correction_ht.locus]

    freq_correction_ht = freq_correction_ht.select(
        AN=hl.enumerate(freq_correction_ht.alleles[1:]).map(
            lambda i: (
                i[1],
                hl.map(
                    lambda an, an_fail_adj, m: hl.if_else(
                        m.get("group") == "adj", an - an_fail_adj.get(i[0] + 1, 0), an
                    ),
                    freq_correction.AN,
                    freq_correction.AN_fail_adj,
                    freq_correction_ht.strata_meta,
                ),
            )
        ),
        qual_hists=freq_correction.qual_hists[0],
    )
    freq_correction_ht = freq_correction_ht.explode("AN")
    freq_correction_ht = freq_correction_ht.key_by(
        **hl.min_rep(
            freq_correction_ht.locus,
            [freq_correction_ht.alleles[0], freq_correction_ht.AN[0]],
        )
    )
    freq_correction_ht = freq_correction_ht.annotate(AN=freq_correction.AN[1])

    return ht.select("AN"), freq_correction_ht, vmt_hists_ht


def main(args):
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    overwrite = args.overwrite
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = use_test_dataset or test_n_partitions or test_gene
    chrom = args.chrom
    af_threshold = args.af_threshold

    hl.init(
        log="/frequency_an_update.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")
    resources = get_freq_an_resources(overwrite, test, chrom)

    try:
        # Fill this in, can we do this on VDS? run DP and GQ hists over sparse MT
        # then merge hists? Should be fine for variant data histograms but will
        # not work for ref data because of ref blocks. Investigating if we can do
        #  this in the AN script for ref data and then filter to variant sites and then
        #  merge with variant histogram data.
        if args.regenerate_ref_an_and_gq_dp_hists:
            logger.info("Regenerating reference AN and GQ/DP histograms...")
            ref_ht = vep_context.versions["105"].ht()
            ref_ht = ref_ht._filter_partitions(range(5))

            # Retain only 'locus' annotation from context Table
            ref_ht = ref_ht.key_by("locus").select().distinct()

            # Read in VDS.
            vds = get_gnomad_v4_vds(
                release_only=True,
                # test=False,
                test=True,
                filter_partitions=range(2),
                annotate_meta=True,
            )

            group_membership_ht = get_group_membership_ht(
                meta().ht(),
                get_downsampling().ht(),
                get_downsampling(subset="non_ukb").ht(),
            )

            interval_ht = adjust_interval_padding(
                calling_intervals(
                    interval_name="union", calling_interval_padding=0
                ).ht(),
                150,
            )

            # Read in freq and filter to only those.
            freq_ht = get_freq("4.0").ht()
            vds = vds_annotate_adj(vds, freq_ht)

            ht, freq_correction_ht, vmt_hists_ht = (
                compute_allele_number_per_ref_site_with_adj(
                    vds,
                    ref_ht,
                    interval_ht,
                    group_membership_ht,
                )
            )

            pass
        if args.update_freq_an:
            logger.info("Updating AN in freq HT...")
            res = resources.update_freq_an
            res.check_resource_existence()
            freq_ht = res.freq_ht.ht()
            freq_ht = drop_gatk_groupings(freq_ht)
            an_ht = res.an_ht.ht()
            ht = update_freq_an(freq_ht, an_ht)
            ht.write(res.updated_an_freq_ht.path, overwrite=overwrite)

        if args.update_af_dependent_annotations:
            logger.info("Updatings AF dependent annotations...")
            res = resources.update_af_dependent_annotations
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
        "--run-freq-and-dense-annotations",
        help=(
            "Calculate frequencies, histograms, and high AB sites per sample grouping."
        ),
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

    main(args)
