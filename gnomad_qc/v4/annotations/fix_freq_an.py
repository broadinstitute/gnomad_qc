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
    get_is_haploid_expr,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.sparse_mt import compute_stats_per_ref_site

from gnomad_qc.v4.annotations.compute_coverage import (
    adjust_interval_padding,
    get_group_membership_ht,
)
from gnomad_qc.v4.resources.annotations import (
    get_all_sites_an_and_qual_hists,
    get_downsampling,
    get_freq,
)
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


def vds_annotate_adj(
    vds: hl.vds.VariantDataset, freq_ht: Optional[hl.Table] = None
) -> hl.vds.VariantDataset:
    """
    Annotate adj, _het_ad, and select fields to reduce memory usage.

    :param vds: Hail VDS to annotate adj onto variant data.
    :param freq_ht: Optional Hail Table containing frequency information.
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
        ),
        ploidy=hl.if_else(
            rmt.in_non_par
            & get_is_haploid_expr(
                locus_expr=rmt.locus, karyotype_expr=rmt.sex_karyotype
            ),
            1,
            2,
        ),
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
        ploidy=hl.int32(vmt.LGT.ploidy),
    )

    if freq_ht is not None:
        logger.info("Annotating variant data MT with in_freq field...")
        freq_ht = hl.Table(
            hl.ir.TableKeyBy(
                freq_ht._tir, ["locus"], is_sorted=True
            )  # Prevents hail from running sort on HT which is already sorted.
        )
        vmt = vmt.annotate_rows(in_freq=hl.is_defined(freq_ht[vmt.locus]))

    return hl.vds.VariantDataset(rmt, vmt)


def compute_an_and_hists_het_fail_adj_ab(
    vds: hl.vds.VariantDataset,
    group_membership_ht: hl.Table,
    freq_ht: hl.Table,
) -> hl.Table:
    """
    Compute allele number and histograms for het fail adj ab.

    :param vds: Input VariantDataset.
    :param group_membership_ht: Table of samples group memberships.
    :param freq_ht: Table of frequency data.
    :return: Table of allele number and histograms for het fail adj ab.
    """
    vmt = vds.variant_data

    # Only need to keep rows where the locus is in the freq HT and there is at least
    # one non-ref genotype that fails the adj allele balance filter.
    # This saves on computation by not keeping rows that will be filtered out later or
    # don't have any genotypes that contribute to the aggregate counts.
    vmt = vmt.filter_rows(
        vmt.in_freq & hl.agg.any(vmt.adj & vmt.fail_adj_ab & vmt.LGT.is_non_ref())
    )

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

    logger.info("Splitting multi-allelic sites...")
    vmt = hl.experimental.sparse_split_multi(vmt)

    # Now that the multi-allelic sites are split, we can filter genotypes that
    # are non ref (we don't want to count the ref), and correctly handle the
    # allele balance filtering of het non-ref genotypes.
    vmt = vmt.filter_rows(hl.is_defined(freq_ht[vmt.row_key]))
    vmt = vmt.filter_entries(
        vmt.GT.is_non_ref() & ~get_het_ab_adj_expr(vmt.GT, vmt.DP, vmt.AD)
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
    # filtering above, keeping only entries that fail the AB adj condition. This
    # avoids the built-in adj filtering in `agg_by_strata`.
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
        {"AN_het_fail_adj_ab": (lambda t: t.ploidy, hl.agg.sum)},
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
    :return: Table of AN and GQ/DP hists per reference site and Table of AN and GQ/DP hists
        for frequency correction.
    """

    def _get_hists(qual_expr) -> hl.expr.Expression:
        return qual_hist_expr(
            gq_expr=qual_expr[0],
            dp_expr=qual_expr[1],
            adj_expr=qual_expr[2] == 1,
            split_adj_and_raw=True,
        )

    entry_agg_funcs = {
        "AN": (lambda t: t.ploidy, hl.agg.sum),
        "qual_hists": (lambda t: [t.GQ, t.DP, t.adj], _get_hists),
    }

    logger.info("Computing allele number and histograms per reference site...")
    # Below we use just the raw group for qual hist computations because qual hists
    # has its own built-in adj filtering when adj is passed as an argument and will
    # produce both adj and raw histograms.
    ht = compute_stats_per_ref_site(
        vds,
        reference_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["GQ", "DP", "ploidy"],
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
    # they can be subtracted from the partial adj (GQ/DP pass) histograms using
    # merge_histograms to create complete adj (GQ, DP, and AB pass) histograms.
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

    # Get the sum of het GTs that fail the adj allele balance filter for all alleles at
    # each locus.
    het_fail_adj_ab_ht = het_fail_adj_ab_ht.group_by("locus").aggregate(
        AN_het_fail_adj_ab=hl.agg.array_sum(het_fail_adj_ab_ht.AN_het_fail_adj_ab),
        qual_hists_het_fail_adj_ab=hl.struct(
            **{
                k: hl.struct(
                    bin_edges=hl.agg.take(v.bin_edges, 1)[0],
                    bin_freq=hl.agg.array_sum(v.bin_freq),
                    n_smaller=hl.agg.sum(v.n_smaller),
                    n_larger=hl.agg.sum(v.n_larger),
                )
                for k, v in het_fail_adj_ab_ht.qual_hists_het_fail_adj_ab.items()
            }
        ),
    )

    # Annotate the all sites AN and qual hists Table with the number of hets that fail
    # the adj allele balance filter at the locus and add the AN corrected for the
    # allele balance filter to the all sites AN and qual hists Table.
    het_fail_adj_ab = het_fail_adj_ab_ht[ht.locus]
    ht = ht.select(
        "AN",
        "qual_hists",
        AN_het_fail_adj_ab=het_fail_adj_ab.AN_het_fail_adj_ab,
        qual_hists_het_fail_adj_ab=het_fail_adj_ab.qual_hists_het_fail_adj_ab,
        AN_adj_ab=hl.map(
            lambda an, an_fail_adj, m: hl.if_else(
                m.get("group") == "adj", an - hl.coalesce(an_fail_adj, 0), an
            ),
            ht.AN,
            het_fail_adj_ab.AN_het_fail_adj_ab,
            global_annotations.strata_meta,
        ),
    )

    return ht, freq_correction_ht


def main(args):
    """Script to generate all sites AN and v4.0 exomes frequency fix."""
    overwrite = args.overwrite
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test = use_test_dataset or test_n_partitions

    hl.init(
        log="/frequency_an_update.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    hl._set_flags(use_ssa_logs="1")

    try:
        if args.regenerate_ref_an_hists:
            logger.info("Regenerating reference AN and GQ/DP histograms...")

            ref_ht = vep_context.versions["105"].ht()
            ref_ht = ref_ht.key_by("locus").select().distinct()

            vds = get_gnomad_v4_vds(
                release_only=True,
                test=use_test_dataset,
                filter_partitions=range(2) if test_n_partitions else None,
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

            # Filter out interval HT to only include intervals that have at least one
            # locus in the VDS variant data. This saves on computation by reducing the
            # number of ref loci we compute AN and hists on during the test
            if test:
                tmp_interval_ht = interval_ht.annotate(in_interval=interval_ht.interval)
                tmp_ht = vds.variant_data.rows()
                tmp_ht = tmp_ht.annotate(
                    in_interval=tmp_interval_ht[tmp_ht.locus].in_interval
                )
                test_intervals = hl.literal(
                    tmp_ht.aggregate(hl.agg.collect_as_set(tmp_ht.in_interval))
                )
                interval_ht = interval_ht.filter(
                    test_intervals.contains(interval_ht.interval)
                )

            het_fail_adj_ab_ht = compute_an_and_hists_het_fail_adj_ab(
                vds_annotate_adj(vds, freq_ht),
                group_membership_ht,
                freq_ht,
            )
            het_fail_adj_ab_ht = het_fail_adj_ab_ht.checkpoint(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/het_fail_adj_ab_ht.ht",
                overwrite=overwrite,
            )

            an_ht, freq_correction_ht = compute_allele_number_per_ref_site_with_adj(
                vds_annotate_adj(vds),
                ref_ht,
                interval_ht,
                group_membership_ht,
                freq_ht,
                het_fail_adj_ab_ht,
            )
            an_ht = an_ht.checkpoint(
                get_all_sites_an_and_qual_hists(test=test).path, overwrite=overwrite
            )
            an_ht = an_ht.select(AN=an_ht.AN_adj_ab)
            an_ht.write(
                release_all_sites_an(public=False, test=test).path, overwrite=overwrite
            )
            freq_correction_ht.write(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/freq_correction.ht",
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
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--regenerate-ref-an-hists",
        help="Calculate reference allele number and GQ/DP histograms for all sites.",
        action="store_true",
    )

    args = parser.parse_args()

    main(args)
