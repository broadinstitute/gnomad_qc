"""Updates the v4.0 freq HT with the correct AN for the v4.1 release."""

import argparse
import logging
from typing import Optional, Tuple

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    agg_by_strata,
    get_gq_dp_adj_expr,
    get_het_ab_adj_expr,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.sparse_mt import (
    compute_stats_per_ref_site,
    get_allele_number_agg_func,
)

from gnomad_qc.v4.annotations.compute_coverage import (
    adjust_interval_padding,
    get_group_membership_ht,
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
from gnomad_qc.v4.resources.release import release_all_sites_an

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_4.1_exomes_AN_fix")
logger.setLevel(logging.INFO)


def vds_annotate_adj(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Annotate adj, _het_ad, and select fields to reduce memory usage.

    :param vds: Hail VDS to annotate adj onto variant data.
    :return: Hail VDS with adj annotation.
    """
    rmt = vds.reference_data
    vmt = vds.variant_data

    logger.info(
        "Filtering chrY reference data MT entries for XX sex karyotype samples..."
    )

    # Annotating the cols and rows with annotations used in the filter_entries to
    # prevent the need to recompute these annotations for every entry in the
    # filter_entries.
    rmt = rmt.annotate_cols(
        sex_karyotype=vmt.cols()[rmt.col_key].meta.sex_imputation.sex_karyotype
    )
    rmt = rmt.annotate_rows(
        in_y_par=rmt.locus.in_y_par(),
        in_y_nonpar=rmt.locus.in_y_nonpar(),
        in_non_par=~rmt.locus.in_autosome_or_par(),
    )
    rmt.filter_entries(
        (rmt.in_y_par | rmt.in_y_nonpar) & (rmt.sex_karyotype == "XX"),
        keep=False,
    )

    logger.info("Annotating reference data MT with DP and GQ adj...")
    rmt = rmt.annotate_entries(
        adj=get_gq_dp_adj_expr(
            rmt.GQ, rmt.DP, locus_expr=rmt.locus, karyotype_expr=rmt.sex_karyotype
        )
    )

    logger.info("Adjusting sex ploidy for variant data MT and annotating adj info...")
    vmt = vmt.annotate_rows(in_non_par=~vmt.locus.in_autosome_or_par())
    vmt = vmt.annotate_entries(
        LGT=hl.if_else(
            vmt.in_non_par,
            adjusted_sex_ploidy_expr(
                vmt.locus, vmt.LGT, vmt.meta.sex_imputation.sex_karyotype
            ),
            vmt.LGT,
        )
    )
    vmt = vmt.annotate_entries(
        adj=get_gq_dp_adj_expr(vmt.GQ, vmt.DP, gt_expr=vmt.LGT),
        fail_adj_ab=~get_het_ab_adj_expr(vmt.LGT, vmt.DP, vmt.LAD),
    )

    return hl.vds.VariantDataset(rmt, vmt)


def compute_an_and_hists_het_fail_adj_ab(
    vds: hl.vds.VariantDataset,
    group_membership_ht: hl.Table,
) -> hl.Table:
    """
    Compute allele number and histograms for het fail adj ab.

    :param vds: Input VariantDataset.
    :param group_membership_ht: Table of samples group memberships.
    :return: Table of allele number and histograms for het fail adj ab.
    """
    vds = vds_annotate_adj(vds)
    vmt = vds.variant_data

    # Only need to keep rows where there is at least one non-ref genotype that fails
    # the adj allele balance filter. This saves on computation by not keeping rows
    # that will be filtered out later or don't have any genotypes that contribute
    # to the aggregate counts.
    vmt = vmt.filter_rows(hl.agg.any(vmt.adj & vmt.fail_adj_ab & vmt.LGT.is_non_ref()))

    logger.info(
        "Filtering variant data MT entries to those passing GQ and DP adj thresholds, "
        "but failing the het allele balance adj threshold..."
    )
    # Keep only the entries where the genotype passes the DP and GQ adj filters,
    # but fails the adj allele balance filter.
    # This optimization avoids splitting multi-allelic entries that do not contribute
    # to the aggregate counts.
    vmt = vmt.filter_entries(vmt.adj & vmt.fail_adj_ab)
    vmt = vmt.select_entries("LA", "LGT", "DP", "GQ", "LAD")

    logger.info("Splitting multi-allelic sites and filtering to sites in freq HT...")
    vmt = hl.experimental.sparse_split_multi(vmt)

    # Now that the multi-allelic sites are split, we can filter genotypes that
    # are non ref (we don't want to count the ref), and correctly handle the
    #  allele balance filtering of het non-ref genotypes.
    vmt = vmt.filter_entries(
        vmt.GT.is_non_ref()
        & ~get_het_ab_adj_expr(
            vmt.GT, vmt.DP, vmt.AD
        )  # TODO: Why reannotate this? its already fail_adj_ab?
    )

    logger.info(
        "Annotating variant data MT with DP and GQ qual hists for het fail adj allele "
        "balance..."
    )
    vmt = vmt.annotate_rows(
        qual_hists_het_fail_adj_ab=qual_hist_expr(gq_expr=vmt.GQ, dp_expr=vmt.DP)
    )

    logger.info("Computing allele number for het fail adj allele balance...")
    # Set all the group membership fields to "raw" since we already handled the adj
    # filtering above. This avoids the built-in adj filtering in `agg_by_strata`. This is
    # technically the raw - adj group.
    group_membership_ht = group_membership_ht.annotate_globals(
        freq_meta=group_membership_ht.freq_meta.map(
            lambda x: hl.dict(
                x.items().map(
                    lambda m: hl.if_else(m[0] == "group", ("group", "raw"), m)
                )
            )
        )
    )
    ht = agg_by_strata(
        vmt,
        {"AN_het_fail_adj_ab": get_allele_number_agg_func("GT")},
        select_fields=["qual_hists_het_fail_adj_ab"],
        group_membership_ht=group_membership_ht,
    )

    return ht


def compute_allele_number_per_ref_site_with_adj(
    vds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    interval_ht: hl.Table,
    group_membership_ht: hl.Table,
    freq_ht: hl.Table,
    het_fail_adj_ab_ht: hl.Table,
) -> Tuple[hl.Table, hl.Table]:
    """
    Compute allele number per reference site and histograms for frequency correction.

    :param vds: Input VariantDataset.
    :param reference_ht: Table of reference sites.
    :param interval_ht: Table of intervals.
    :param group_membership_ht: Table of samples group memberships.
    :param freq_ht: Table of frequency data.
    :param het_fail_adj_ab_ht: Table of variant data histograms for het fail adj ab.
    :return: Table of allele number per reference site and Table of AN and GQ/DP hists
        for frequency correction.
    """

    def _get_hists(qual_expr) -> hl.expr.Expression:
        return qual_hist_expr(
            gq_expr=qual_expr[0],
            dp_expr=qual_expr[1],
            adj_expr=qual_expr[2] == 1,
            split_adj_and_raw=True,
        )

    vds = vds_annotate_adj(vds)
    entry_agg_funcs = {
        "AN": get_allele_number_agg_func("LGT"),
        "qual_hists": (lambda t: [t.GQ, t.DP, t.adj], _get_hists),
    }

    logger.info("Computing allele number and histograms per reference site...")
    ht = compute_stats_per_ref_site(
        vds,
        reference_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["GQ", "DP"],
        entry_agg_group_membership={"qual_hists": [{"group": "raw"}]},
    )
    ht = ht.checkpoint(hl.utils.new_temp_file("an_hist_ref_sites", "ht"))

    logger.info(
        "Adjusting allele number and histograms for het fail adj allele balance..."
    )
    global_annotations = ht.index_globals()
    freq_correction = ht[freq_ht.locus]
    het_fail_adj_ab = het_fail_adj_ab_ht[freq_ht.key]
    qual_hists = freq_correction.qual_hists[0].qual_hists

    # Convert the het fail adj allele balance histograms to negative values so that
    # they can be subtracted from the reference raw histograms using merge_histograms
    # to create adj histograms.
    sub_hists = het_fail_adj_ab.qual_hists_het_fail_adj_ab
    sub_hists = sub_hists.annotate(
        **{
            k: v.annotate(
                bin_freq=v.bin_freq.map(lambda x: -x),
                n_smaller=-v.n_smaller,
                n_larger=-v.n_larger,
            )
            for k, v in sub_hists.items()
        }
    )

    freq_correction_ht = freq_ht.select(
        AN=hl.coalesce(
            hl.map(
                lambda an, an_fail_adj, m: hl.if_else(
                    m.get("group") == "adj", an - hl.coalesce(an_fail_adj, 0), an
                ),
                freq_correction.AN,
                het_fail_adj_ab.AN_het_fail_adj_ab,
                global_annotations.strata_meta,
            ),
            freq_correction.AN,
        ),
        qual_hists=freq_correction.qual_hists[0].annotate(
            qual_hists=hl.struct(
                **{
                    k: merge_histograms([v, sub_hists[k]])
                    for k, v in qual_hists.items()
                }
            )
        ),
    )
    freq_correction_ht = freq_correction_ht.annotate_globals(**global_annotations)

    return ht.select("AN", "qual_hists"), freq_correction_ht


def update_freq_an_and_hists(
    freq_ht: hl.Table,
    freq_correction_ht: hl.Table,
) -> hl.Table:
    """
    Update frequency HT AN field and histograms with the correct AN from the AN HT.

    This module updates the freq_ht AN field which is an array of ints with the correct
    AN from the AN HT integer array annotation. It uses the freq_correction_ht dictionary and the
    freq_ht dictionary to match the indexing of array annotations so each group is
    correctly updated. It then recomputes AF, AC/AN.

    :param freq_ht: Frequency HT to update.
    :param an_ht: AN HT to update from.
    :return: Updated frequency HT.
    """
    freq_meta = freq_ht.index_globals().freq_meta
    freq_correction_meta = freq_correction_ht.index_globals().strata_meta

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
    freq_correction_ht = freq_correction_ht.transmute(
        freq=freq_correction_ht.AN.map(
            lambda x: hl.struct(
                AC=0, AN=x, homozygote_count=0, AF=hl.missing(hl.tfloat64)
            )
        )
    )
    freq_ht = freq_ht.annotate(
        freq_correction_freq=freq_correction_ht[freq_ht.locus].freq
    )
    freq_ht = freq_ht.annotate_globals(freq_correction_strata_meta=freq_correction_meta)

    logger.info("Merging freq arrays from v4 and freq correction HT...")
    freq, freq_meta = merge_freq_arrays(
        [freq_ht.freq, freq_ht.freq_correction_freq],
        [freq_ht.freq_meta, freq_ht.freq_correction_strata_meta],
    )
    freq_ht = freq_ht.annotate(freq=freq)
    freq_ht = freq_ht.annotate_globals(freq_meta=freq_meta)

    logger.info("Updating freq_ht with corrected GQ and DP histograms...")
    corrected_index = freq_correction_ht[freq_ht.key]
    freq_ht = freq_ht.annotate(
        qual_hists=freq_ht.qual_hists.annotate(
            gq_hist_all=corrected_index.qual_hists.gq_hist_all,
            dp_hist_all=corrected_index.qual_hists.dp_hist_all,
        ),
        raw_qual_hists=freq_ht.raw_qual_hists.annotate(
            gq_hist_all=corrected_index.raw_qual_hists.gq_hist_all,
            dp_hist_all=corrected_index.raw_qual_hists.dp_hist_all,
        ),
    )
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
            "high_ab_hets_by_group": ht.high_ab_hets_by_group,
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=["gatk_version"],
        keep=False,
    )
    ht = ht.annotate(
        freq=array_exprs["freq"],
        high_ab_hets_by_group=array_exprs["high_ab_hets_by_group"],
    )
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )
    return ht


def main(args):
    """Script to generate all sites AN and v4.0 exomes frequency fix."""
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
    hl._set_flags(
        use_ssa_logs="1"
    )  # TODO: Do we want to try the new restart run option out?

    try:
        # Fill this in, can we do this on VDS? run DP and GQ hists over sparse MT
        # then merge hists? Should be fine for variant data histograms but will
        # not work for ref data because of ref blocks. Investigating if we can do
        #  this in the AN script for ref data and then filter to variant sites and then
        #  merge with variant histogram data.
        if args.regenerate_ref_an_hists:
            logger.info("Regenerating reference AN and GQ/DP histograms...")

            ref_ht = vep_context.versions["105"].ht()
            if test:
                ref_ht = ref_ht._filter_partitions(range(5))

            ref_ht = ref_ht.key_by("locus").select().distinct()

            vds = get_gnomad_v4_vds(
                release_only=True,
                test=test,
                filter_partitions=range(2) if test else None,
                annotate_meta=True,
            )
            freq_ht = get_freq("4.0").ht()
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

            het_fail_adj_ab_ht = compute_an_and_hists_het_fail_adj_ab(
                vds,
                group_membership_ht,
            )
            het_fail_adj_ab_ht = het_fail_adj_ab_ht.checkpoint(
                "gs://gnomad/v4.1/temp/frequency_fix/het_fail_adj_ab_ht.ht",
                overwrite=overwrite,
            )

            an_ht, freq_correction_ht = compute_allele_number_per_ref_site_with_adj(
                vds,
                ref_ht,
                interval_ht,
                group_membership_ht,
                freq_ht,
                het_fail_adj_ab_ht,
            )

            an_ht.write(
                release_all_sites_an(public=False, test=test).path, overwrite=overwrite
            )
            freq_correction_ht.write(
                "gs://gnomad/v4.1/temp/frequency_fix/freq_an_correction.ht",
                overwrite=overwrite,
            )

        if args.update_freq_an_and_hists:
            logger.info("Updating AN in freq HT...")

            freq_correction_ht = hl.read_table(
                "gs://gnomad/v4.1/temp/frequency_fix/freq_an_correction.ht"
            )
            freq_ht = get_freq(
                version="4.0",
                test=test,
                hom_alt_adjusted=False,
                chrom=chrom,
                finalized=False,
            )
            freq_ht = drop_gatk_groupings(freq_ht)
            ht = update_freq_an_and_hists(freq_ht, freq_correction_ht)
            ht.write(
                get_freq(
                    version="4.1",
                    test=test,
                    hom_alt_adjusted=False,
                    chrom=chrom,
                    finalized=False,
                ).path,
                overwrite=overwrite,
            )

        if args.update_af_dependent_anns:
            logger.info("Updatings AF dependent annotations...")

            ht = get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=False,
                chrom=chrom,
                finalized=False,
            ).ht()

            logger.info("Correcting call stats, qual AB hists, and age hists...")
            ht = correct_for_high_ab_hets(ht, af_threshold=af_threshold)

            logger.info("Computing FAF & grpmax...")
            ht = generate_faf_grpmax(ht)

            logger.info("Calculating InbreedingCoeff...")
            ht = compute_inbreeding_coeff(ht)

            logger.info("High AB het corrected frequency HT schema...")
            ht.describe()

            logger.info("Writing corrected frequency Table...")
            ht.write(
                get_freq(
                    version="4.1",
                    test=test,
                    hom_alt_adjusted=True,
                    chrom=chrom,
                    finalized=False,
                ).path,
                overwrite=overwrite,
            )

        if args.finalize_freq_ht:
            logger.info("Writing final frequency Table...")
            freq_ht = get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=True,
                chrom=chrom,
                finalized=False,
            ).ht()
            freq_ht = create_final_freq_ht(freq_ht)

            logger.info("Final frequency HT schema...")
            ht.describe()
            ht.write(
                get_freq(version="4.1", test=test, finalized=True).path,
                overwrite=overwrite,
            )
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("frequency_fix"))


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
        "--regenerate-ref-an-hists",
        help="Calculate reference allele number and GQ/DP histograms for all sites.",
        action="store_true",
    )
    parser.add_argument(
        "--update-freq-an-and-hists",
        help=(
            "Update frequency HT AN field and GQ and DP hists with the correct AN and"
            " hists from the AN HT."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--update-af-dependent-anns",
        help=(
            "Update AF dependent annotations in frequency HT with the new AF based on"
            " updated AN."
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
