"""
Updates the v4.0 freq HT with the correct AN for the v4.1 release.

There is an error in the v4.0 freq HT where the AN can be incorrect for variants that
are only present in one of the UKB or non-UKB subsets. This is because the frequency
data was computed on each subset separately and then combined. In order to filter the
dataset to only samples that are in the UKB or non-UKB subsets, the
`hail.vds.filter_samples` was used. Unfortunately, the behavior of this function was
not the expected behavior due to the following line of code:

    .. code-block::

        variant_data = variant_data.filter_rows(hl.agg.count() > 0)

This line of code filters any row that has no genotypes in the sample subset. When the
VariantDataset is densified prior to the frequency calculations this results in no
reference data being filled in for that row, unless a row is part of multiallelic site that is not
exclusive to the sample subset, and therefore a missing AN for that row. When combining
the frequency data for the UKB and non-UKB subsets, the AN for these rows will
only have the AN for the subset that has non-ref genotypes for that row.

This script performs the following steps:

    - Generates allele number and GQ/DP histograms for all sites in the exome calling
      interval.

        - For raw genotypes, the AN is just the sum of the ploidy (more specifics on
          this below) for all genotypes.
        - For adj genotypes, the AN is the sum of the ploidy for all genotypes that pass
          the GQ, DP, and allele balance adj filters. Since this is for a site, rather
          than a specific variant, this number will not always be the same as the AN for
          a specific variant at that site. This AN can be smaller if it is a
          multi-allelic site where there are non-ref genotypes for another variant that
          fails the adj filters.
        - For raw genotypes on non-PAR regions of the X or Y chromosome, the AN is the
          sum of the ploidy after adjusting for sex karyotype, BUT NOT taking into
          account the genotype when adjusting the ploidy. For instance, a het genotype
          for an XY sample will still have a ploidy of 1, even though the genotype is a
          het and is therefore set to missing for the frequency calculations.

    - Generates allele number and GQ/DP histograms for the frequency correction.

"""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    agg_by_strata,
    annotate_and_index_source_mt_for_sex_ploidy,
    get_gq_dp_adj_expr,
    get_het_ab_adj_expr,
    get_is_haploid_expr,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.sparse_mt import compute_stats_per_ref_site

from gnomad_qc.v4.annotations.compute_coverage import (
    adjust_interval_padding,
    get_exomes_group_membership_ht,
)
from gnomad_qc.v4.annotations.generate_freq import (
    compute_inbreeding_coeff,
    correct_for_high_ab_hets,
    create_final_freq_ht,
    generate_faf_grpmax,
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


def prep_vds_for_all_sites_stats(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Prepare VDS for all sites stats.

    Adds the following annotations to the VDS:

        - 'ploidy': The ploidy of the genotype adjusted for sex karyotype.
        - 'adj': The adj filter for GQ and DP.

    The following steps are performed on the reference data MT:

        - Segments the reference data MT at the PAR boundaries on chrX and chrY.

            - This is done to ensure that entries are annotated with the correct
              'ploidy' field when adjusting for sex karyotype.
            - If the END field of an entry spans a PAR boundary, the entry will be
              split into two separate entries so that a locus 'in_non_par' annotation
              can be used to correctly annotate if the full reference block defined by
              an entry is non-PAR or not.

        - Annotates each entry with the ploidy of the reference block adjusting for
          sex karyotype.

            - chrY non-PAR XX entries are set to missing ploidy.
            - chrX and chrY non-PAR XY entries are set to ploidy 1.
            - All other entries are set to ploidy 2.

        - Annotates each entry with the adj filter for GQ and DP taking into account
          sex karyotype for DP.

    The following steps are performed on the variant data MT:

        - Annotates each entry with the ploidy of the genotype adjusting for sex
          karyotype.

            - chrY non-PAR XX entries are set to missing ploidy.
            - chrX and chrY non-PAR XY entries are set to ploidy 1 if LGT is defined.

                - This deviates from our standard ploidy adjustment where if the
                  genotype is a het, the ploidy is set to missing. This is because we
                  are adjusting ploidy before splitting multi-allelic sites and this
                  missing annotation will be propagated to all split entries including
                  the reference genotypes causing discrepancies in adjusting ploidy
                  before splitting and after splitting. Since we are using this VDS
                  downstream to compute a reference AN, we want to keep this as a
                  ploidy of 1 even if the genotype is a het.
            - All other entries are set to the ploidy of the genotype.

        - Annotates each entry with the adj filter for GQ and DP taking into account
          sex karyotype for DP.

          - The sex karyotype and locus are used instead of the genotype to determine
            if the genotype is haploid or not. This is for the same reason as described
            above for the ploidy adjustment of non-PAR XY het genotypes.

    :param vds: Hail VDS to annotate adj onto variant data.
    :return: Hail VDS with adj annotation.
    """
    rmt = vds.reference_data
    vmt = vds.variant_data
    vmt = vmt.annotate_cols(sex_karyotype=vmt.meta.sex_imputation.sex_karyotype)

    # In order to correctly annotate the ploidy on the reference data MT, there needs to
    # be an END field where the start and end position are both in PAR or both in
    # non-PAR regions. We can make sure this is the case by segmenting the reference
    # data MT at the PAR boundaries.
    # This can be done by building an array of intervals that cover the entire reference
    # genome, then segment the intervals at the PAR boundaries and pass that interval
    # list to hl.vds.segment_reference_blocks. This will segment the reference data MT
    # so any sample with an END field that spans a PAR boundary will be split into two
    # separate entries.
    rg = rmt.locus.dtype.reference_genome
    rg_intervals = hl.array(
        [
            hl.parse_locus_interval(f"[{contig}:1-{rg.lengths[contig]}]")
            for contig in rg.contigs[:24]
        ]
    )
    rg_intervals = hl.Table.parallelize(
        rg_intervals.map(lambda x: hl.struct(interval=x)),
        schema=hl.tstruct(interval=rg_intervals.dtype.element_type),
        key="interval",
    )

    par_boundaries = []
    for par_interval in rg.par:
        par_boundaries.append(par_interval.start)
        par_boundaries.append(par_interval.end)

    rg_intervals = hl.segment_intervals(rg_intervals, par_boundaries).checkpoint(
        hl.utils.new_temp_file("segment_par_boundaries", "ht")
    )
    rmt = hl.vds.segment_reference_blocks(rmt, rg_intervals)

    # Annotating the cols and rows with annotations used in the filter_entries to
    # prevent the need to recompute these annotations for every entry in the
    # filter_entries.
    rmt = rmt.annotate_cols(sex_karyotype=vmt.cols()[rmt.col_key].sex_karyotype)
    rmt = rmt.annotate_rows(in_non_par=~rmt.locus.in_autosome_or_par())

    logger.info("Annotating reference data MT with GQ and DP adj as well as ploidy...")
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

    logger.info(
        "Adjusting sex ploidy for variant data MT and annotating with DP/GQ adj..."
    )
    # An optimization that annotates the locus's source matrix table with the
    # fields in the case statements below, so they are not re-computed for every entry.
    c_idx, r_idx = annotate_and_index_source_mt_for_sex_ploidy(
        locus_expr=vmt.locus, karyotype_expr=vmt.sex_karyotype
    )
    ploidy_expr = (
        hl.case(missing_false=True)
        .when(~r_idx.in_non_par, vmt.LGT.ploidy)
        .when(c_idx.xx & (r_idx.y_par | r_idx.y_nonpar), hl.missing(hl.tint32))
        # For XY non-PAR het genotypes, the reference ploidy is 1, so we set it to one
        # even if the genotype is a het. This is to keep the reference ploidy consistent
        # with the adjusted ploidy after splitting multi-allelic sites.
        .when(c_idx.xy & r_idx.in_non_par & hl.is_defined(vmt.LGT.ploidy), 1)
        .default(vmt.LGT.ploidy)
    )

    # Use the locus and karyotype to annotate the GQ and DP adj correctly for sex
    # ploidy.
    vmt = vmt.annotate_entries(
        adj=get_gq_dp_adj_expr(
            vmt.GQ, vmt.DP, locus_expr=vmt.locus, karyotype_expr=vmt.sex_karyotype
        ),
        ploidy=ploidy_expr,
    )

    return hl.vds.VariantDataset(rmt, vmt)


def compute_an_and_hists_het_fail_adj_ab(
    vds: hl.vds.VariantDataset,
    group_membership_ht: hl.Table,
    freq_ht: hl.Table,
) -> hl.Table:
    """
    Compute allele number and histograms for het fail adj ab and non-PAR XY het.

    This module computes the allele number and quality histograms for only the
    genotypes that are het and fail the adj allele balance or are non-PAR XY het
    genotypes.

    :param vds: Input VariantDataset.
    :param group_membership_ht: Table of samples group memberships.
    :param freq_ht: Table of frequency data.
    :return: Table of allele number and histograms for het fail adj ab.
    """
    vmt = vds.variant_data

    logger.info(
        "Filtering variant data MT entries to those passing GQ and DP adj thresholds, "
        "but failing the het allele balance adj threshold..."
    )
    # Use the locus and karyotype to determine the GQ and DP adj correctly for sex
    # ploidy on unsplit data. The adj allele balance filter is not dependent on sex
    # ploidy, so we can annotate that before adjusting the ploidy of the genotype.
    vmt = vmt.select_cols(
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        xy=vmt.meta.sex_imputation.sex_karyotype == "XY",
    )
    vmt = vmt.annotate_rows(in_non_par=~vmt.locus.in_autosome_or_par())
    vmt = vmt.annotate_entries(
        adj=get_gq_dp_adj_expr(
            vmt.GQ, vmt.DP, locus_expr=vmt.locus, karyotype_expr=vmt.sex_karyotype
        ),
        fail_het_ab_adj=~get_het_ab_adj_expr(vmt.LGT, vmt.DP, vmt.LAD),
        is_nonpar_xy_het=vmt.in_non_par & vmt.xy & vmt.LGT.is_het(),
    )

    # Keep only the entries where the genotype passes the DP and GQ adj filters,
    # but fails the adj allele balance filter. This optimization avoids splitting
    # multi-allelic entries that do not contribute to the aggregate counts.
    vmt = vmt.filter_entries((vmt.fail_het_ab_adj & vmt.adj) | vmt.is_nonpar_xy_het)

    # Only need to keep rows where the locus is in the freq HT and there is at least
    # one non-ref genotype. This saves on computation by not keeping rows that will be
    # filtered out later or don't have any genotypes that contribute to the aggregate
    # counts.
    locus_freq_ht = hl.Table(
        hl.ir.TableKeyBy(
            freq_ht._tir, ["locus"], is_sorted=True
        )  # Prevents hail from running sort on HT which is already sorted.
    )
    vmt = vmt.filter_rows(
        hl.is_defined(locus_freq_ht[vmt.locus]) & hl.agg.any(vmt.LGT.is_non_ref())
    )
    vmt = vmt.select_entries("LA", "LGT", "DP", "GQ", "LAD", "adj")

    logger.info("Splitting multi-allelic sites...")
    vmt = hl.experimental.sparse_split_multi(vmt)

    # Now that the multi-allelic sites are split, we can filter genotypes that are non
    # ref (we don't want to count the ref), and correctly handle the allele balance
    # filtering of het non-ref genotypes.
    vmt = vmt.filter_rows(hl.is_defined(freq_ht[vmt.row_key]))

    # Redo the 'is_nonpar_xy_het' and 'fail_het_ab_adj' annotations after splitting
    # multi-allelic sites to correctly handle the allele specific adj allele balance
    # and non-PAR XY het genotypes.
    gt_expr = adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype)
    vmt = vmt.annotate_entries(
        is_nonpar_xy_het=vmt.in_non_par & vmt.xy & vmt.GT.is_het(),
        GT=gt_expr,
        fail_het_ab_adj=~get_het_ab_adj_expr(gt_expr, vmt.DP, vmt.AD),
    )
    vmt = vmt.filter_entries(
        (vmt.GT.is_non_ref() & vmt.fail_het_ab_adj) | vmt.is_nonpar_xy_het
    )

    logger.info(
        "Annotating variant data MT with DP and GQ qual hists for het fail adj allele "
        "balance..."
    )
    adj_hist_include_expr = vmt.adj & (vmt.fail_het_ab_adj | vmt.is_nonpar_xy_het)
    vmt = vmt.annotate_rows(
        qual_hists_het_fail_adj_ab_or_is_nonpar_xy_het=qual_hist_expr(
            gq_expr=hl.or_missing(adj_hist_include_expr, vmt.GQ),
            dp_expr=hl.or_missing(adj_hist_include_expr, vmt.DP),
        ),
        # NOTE: This is not needed for the adjustment of the raw qual hists. They are
        # not impacted by the non-PAR XY het genotypes being set to None because the
        # GQ and DP are maintained. Adj qual hists are impacted by the non-PAR XY het
        # genotypes being set to None because the adj annotation is dependent on the
        # GT, and therefore the GQ and DP are not included in the histogram
        # calculations.
        raw_qual_hists_is_nonpar_xy_het=qual_hist_expr(
            gq_expr=hl.or_missing(vmt.is_nonpar_xy_het, vmt.GQ),
            dp_expr=hl.or_missing(vmt.is_nonpar_xy_het, vmt.DP),
        ),
    )

    logger.info("Computing allele number for het fail adj allele balance...")
    an_nonpar_xy_het_meta = [{"group": "raw"}, {"group": "raw", "subset": "non_ukb"}]
    ht = agg_by_strata(
        vmt.select_entries(
            "adj",
            # Adjust the ploidy for non-PAR XY het genotypes to 1.
            ploidy=hl.if_else(vmt.is_nonpar_xy_het, 1, vmt.GT.ploidy),
            is_nonpar_xy_het=vmt.is_nonpar_xy_het,
        ),
        {
            "AN_het_fail_adj_ab_or_nonpar_xy_het": (lambda t: t.ploidy, hl.agg.sum),
            "AN_nonpar_xy_het": (lambda t: t.is_nonpar_xy_het, hl.agg.count_where),
        },
        select_fields=[
            "qual_hists_het_fail_adj_ab_or_is_nonpar_xy_het",
            "raw_qual_hists_is_nonpar_xy_het",
        ],
        group_membership_ht=group_membership_ht,
        entry_agg_group_membership={"AN_nonpar_xy_het": an_nonpar_xy_het_meta},
    )
    ht = ht.naive_coalesce(1000).checkpoint(
        hl.utils.new_temp_file("an_qual_hist_adjust", "ht")
    )

    # Adjust the raw groups using the count of 'non_par_xy_hets' rather than the sum of
    # the ploidies of the failed AB adj hets.
    freq_meta = group_membership_ht.index_globals().freq_meta
    an_nonpar_xy_het_meta = hl.array(an_nonpar_xy_het_meta)
    ht = ht.select(
        AN_adjust=hl.map(
            lambda m, an: hl.if_else(
                an_nonpar_xy_het_meta.contains(m),
                ht.AN_nonpar_xy_het[an_nonpar_xy_het_meta.index(m)],
                an,
            ),
            freq_meta,
            ht.AN_het_fail_adj_ab_or_nonpar_xy_het,
        ),
        qual_hists_adjust=hl.struct(
            qual_hists=ht.qual_hists_het_fail_adj_ab_or_is_nonpar_xy_het,
            raw_qual_hists=ht.raw_qual_hists_is_nonpar_xy_het,
        ),
    )

    return ht


def compute_allele_number_per_ref_site(
    vds: hl.vds.VariantDataset,
    reference_ht: hl.Table,
    interval_ht: hl.Table,
    group_membership_ht: hl.Table,
) -> hl.Table:
    """
    Compute allele number per reference site and histograms for frequency correction.

    :param vds: Input VariantDataset.
    :param reference_ht: Table of reference sites.
    :param interval_ht: Table of intervals.
    :param group_membership_ht: Table of samples group memberships.
    :return: Table of AN and GQ/DP hists per reference site.
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
    return ht


def get_all_sites_nonref_adjustment(
    vds: hl.vds.VariantDataset,
    group_membership_ht: hl.Table,
) -> hl.Table:
    """
    Compute allele number and histograms for het non ref adjustments.

    The adjustment numbers computed in `compute_an_and_hists_het_fail_adj_ab` are
    correct for per allele AN and histogram adjustments for the frequency data, but
    they are not correct for the all sites AN and histogram adjustments. This is
    because when computing the sum of all per allele adjustment values, the
    heterozygous non-ref genotypes can be counted 2 times if both alleles are adjusted.

    This function will compute the number of het non-ref genotypes that are being
    counted twice so the all sites AN and histograms can be adjusted correctly.

    :param vds: Input VariantDataset.
    :param group_membership_ht: Table of samples group memberships.
    :return: Table of allele number and histograms for het fail adj ab.
    """
    vmt = vds.variant_data
    vmt = vmt.select_cols(
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        xy=vmt.meta.sex_imputation.sex_karyotype == "XY",
    )
    vmt = vmt.annotate_rows(in_non_par=~vmt.locus.in_autosome_or_par())

    # Filter to only het non-ref genotypes that are either failing the adj allele
    # balance for both alleles or are non-PAR XY het genotypes.
    vmt = vmt.filter_entries(
        vmt.LGT.is_het_non_ref()
        & (
            (vmt.in_non_par & vmt.xy)
            | (
                (vmt.LAD[vmt.LGT[0]] / vmt.DP < 0.2)
                & (vmt.LAD[vmt.LGT[1]] / vmt.DP < 0.2)
            )
        )
    )
    vmt = vmt.filter_rows(hl.agg.any(hl.is_defined(vmt.LGT)))

    # Add the DQ/DP adj filter to the entries.
    vmt = vmt.annotate_entries(
        adj=get_gq_dp_adj_expr(
            vmt.GQ, vmt.DP, locus_expr=vmt.locus, karyotype_expr=vmt.sex_karyotype
        ),
    )

    # Calculate adj histograms for the het non-ref genotypes.
    vmt = vmt.annotate_rows(
        qual_hists_het_nonref=qual_hist_expr(
            gq_expr=hl.or_missing(vmt.adj, vmt.GQ),
            dp_expr=hl.or_missing(vmt.adj, vmt.DP),
        )
    )

    vmt = vmt.select_entries(
        "adj",
        is_nonpar_xy_het=vmt.in_non_par & vmt.xy & vmt.LGT.is_het_non_ref(),
        ploidy=adjusted_sex_ploidy_expr(vmt.locus, vmt.LGT, vmt.sex_karyotype).ploidy,
    )

    # Compute the number of het non-ref genotypes that are being counted twice.
    # XY het genotypes in non-PAR regions are set to ploidy 1, so we can use the
    # ploidy to determine AN.
    # The raw AN adjustment is only the count of XY het genotypes in non-PAR regions
    # since it doesn't include adj genotypes.
    an_nonref_raw_meta = [{"group": "raw"}, {"group": "raw", "subset": "non_ukb"}]
    ht = agg_by_strata(
        vmt.annotate_entries(ploidy=hl.if_else(vmt.is_nonpar_xy_het, 1, vmt.ploidy)),
        {
            "AN_het_nonref_adj": (lambda t: t.ploidy, hl.agg.sum),
            "AN_het_nonref_raw": (lambda t: t.is_nonpar_xy_het, hl.agg.count_where),
        },
        select_fields=["qual_hists_het_nonref"],
        group_membership_ht=group_membership_ht,
        entry_agg_group_membership={"AN_het_nonref_raw": an_nonref_raw_meta},
    )
    ht = ht.key_by("locus").drop("alleles")
    ht = ht.checkpoint(hl.utils.new_temp_file("an_qual_hist_het_nonref_adjust", "ht"))

    # Modify the raw groups to use 'AN_het_nonref_raw' instead of 'AN_het_nonref_adj'.
    freq_meta = group_membership_ht.index_globals().freq_meta
    an_nonref_raw_meta = hl.array(an_nonref_raw_meta)
    ht = ht.select(
        "qual_hists_het_nonref",
        AN_het_nonref=hl.map(
            lambda m, an: hl.if_else(
                an_nonref_raw_meta.contains(m),
                ht.AN_het_nonref_raw[an_nonref_raw_meta.index(m)],
                an,
            ),
            freq_meta,
            ht.AN_het_nonref_adj,
        ),
    )

    return ht


def get_sub_hist_expr(
    hist_expr: hl.expr.StructExpression,
    sub_hist_expr: hl.expr.StructExpression,
) -> hl.expr.StructExpression:
    """
    Subtract sub histograms from histograms.

    :param hist_expr: Histograms to subtract from.
    :param sub_hist_expr: Histograms to subtract.
    :return: Subtracted histograms.
    """
    sub_hist_expr = sub_hist_expr.annotate(
        **{
            k: v.annotate(
                bin_freq=v.bin_freq.map(lambda x: -x),
                n_smaller=-v.n_smaller,
                n_larger=-v.n_larger,
            )
            for k, v in sub_hist_expr.items()
        }
    )
    hist_expr = hl.if_else(
        hl.is_defined(sub_hist_expr),
        hl.struct(
            **{k: merge_histograms([v, sub_hist_expr[k]]) for k, v in hist_expr.items()}
        ),
        hist_expr,
    )

    return hist_expr


def adjust_per_site_an_and_hists_for_frequency(
    ht: hl.Table,
    freq_ht: hl.Table,
    het_fail_adj_ab_ht: hl.Table,
) -> hl.Table:
    """
    Adjust allele number and histograms for frequency correction.

    The all sites values on `ht` are only the reference AN and quality histograms
    at the site level (keyed by locus) and don't take into account the differences in
    the alleles. This function will adjust the all sites AN and quality histograms
    by the `het_fail_adj_ab_ht` (computed by `compute_an_and_hists_het_fail_adj_ab`)
    so they can be used for the variant frequency correction.

    The reasons a variant's AN and quality histograms can be different from the all
    sites values are:

        - The reference AN adj filter doesn't take into account the allele balance adj
          filter, so the variant AN can be smaller than the reference AN if there are
          samples with het genotypes that fail the adj allele balance filter.
        - The reference AN doesn't remove the non-PAR het genotypes in XY samples, so
          the variant AN can be smaller than the reference AN if there are XY samples
          with a het genotype in a non-PAR region on chrX or chrY.

    :param ht: Table of AN and GQ/DP hists per reference site.
    :param freq_ht: Table of frequency data.
    :param het_fail_adj_ab_ht: Table of variant data histograms for het fail adj ab.
    :return: Table of adjusted AN and GQ/DP hists for frequency correction.
    """
    logger.info(
        "Adjusting allele number and histograms for het fail adj allele balance..."
    )
    global_annotations = ht.index_globals()
    freq_correction = ht[freq_ht.locus]
    het_fail_adj_ab = het_fail_adj_ab_ht[freq_ht.key]
    qual_hists = freq_correction.qual_hists[0].qual_hists
    het_fail_adj_ab_qual_hists = het_fail_adj_ab.qual_hists_adjust.qual_hists

    # Convert the het fail adj allele balance histograms to negative values so that
    # they can be subtracted from the partial adj (GQ/DP pass) histograms using
    # merge_histograms to create complete adj (GQ, DP, and AB pass) histograms.
    freq_correction_ht = freq_ht.select(
        AN=hl.coalesce(
            hl.map(
                lambda an, an_fail_adj, m: an - hl.coalesce(an_fail_adj, 0),
                freq_correction.AN,
                het_fail_adj_ab.AN_adjust,
                global_annotations.strata_meta,
            ),
            freq_correction.AN,
        ),
        qual_hists=freq_correction.qual_hists[0].annotate(
            qual_hists=get_sub_hist_expr(qual_hists, het_fail_adj_ab_qual_hists)
        ),
    )
    freq_correction_ht = freq_correction_ht.annotate_globals(**global_annotations)

    return freq_correction_ht


def adjust_all_sites_an_and_hists(
    ht: hl.Table,
    het_fail_adj_ab_ht: hl.Table,
    het_nonref_ht: hl.Table,
) -> hl.Table:
    """
    Adjust all sites allele number and histograms Table.

    :param ht: Table of AN and GQ/DP hists per reference site.
    :param het_fail_adj_ab_ht: Table of variant data histograms for het fail adj ab.
    :param het_nonref_ht: Table of AN and GQ/DP hists for only het non-ref calls.
    :return: Adjusted all sites AN and GQ/DP hists Table.
    """
    # Get the sum of het GTs that fail the adj allele balance filter for all alleles at
    # each locus.
    het_fail_adj_ab_ht = het_fail_adj_ab_ht.group_by("locus").aggregate(
        AN_adjust=hl.agg.array_sum(het_fail_adj_ab_ht.AN_adjust),
        qual_hists_adjust=hl.struct(
            **{
                k: hl.struct(
                    bin_edges=hl.agg.take(v.bin_edges, 1)[0],
                    bin_freq=hl.agg.array_sum(v.bin_freq),
                    n_smaller=hl.agg.sum(v.n_smaller),
                    n_larger=hl.agg.sum(v.n_larger),
                )
                for k, v in het_fail_adj_ab_ht.qual_hists_adjust.qual_hists.items()
            }
        ),
    )

    # Annotate the all sites AN and qual hists Table with the number of hets that fail
    # the adj allele balance filter and het non-refs at the locus and add the AN
    # corrected for the allele balance filter to the all sites AN and qual hists Table.
    het_fail_adj_ab = het_fail_adj_ab_ht[ht.locus]
    het_nonref = het_nonref_ht[ht.locus]
    ht = ht.select(
        AN_pre_adjust=ht.AN,
        qual_hists_pre_adjust=ht.qual_hists[0],
        AN_adjust=het_fail_adj_ab.AN_adjust,
        AN_het_nonref=het_nonref.AN_het_nonref,
        qual_hists_adjust=het_fail_adj_ab.qual_hists_adjust,
        qual_hists_het_nonref=het_nonref.qual_hists_het_nonref,
    )
    # Remove duplicate counts of the het non-ref sites from the adjustment histograms.
    ht = ht.annotate(
        _tmp_hist_adjust=get_sub_hist_expr(
            ht.qual_hists_adjust, ht.qual_hists_het_nonref
        )
    )
    # When present, remove duplicate counts of the het non-ref sites from the
    # adjustment AN and then adjust the AN based on result. Adjust qual histograms with
    # the previously het non-ref adjusted histograms.
    ht = ht.annotate(
        AN=hl.if_else(
            hl.is_defined(ht.AN_adjust),
            hl.enumerate(ht.AN_pre_adjust).map(
                lambda x: (
                    x[1]
                    - (
                        ht.AN_adjust[x[0]]
                        - hl.if_else(
                            hl.is_defined(ht.AN_het_nonref), ht.AN_het_nonref[x[0]], 0
                        )
                    )
                )
            ),
            ht.AN_pre_adjust,
        ),
        qual_hists=ht.qual_hists_pre_adjust.annotate(
            qual_hists=get_sub_hist_expr(
                ht.qual_hists_pre_adjust.qual_hists, ht._tmp_hist_adjust
            )
        ),
    ).drop("_tmp_hist_adjust")

    return ht


def update_freq_an_and_hists(
    freq_ht: hl.Table,
    freq_correction_ht: hl.Table,
) -> hl.Table:
    """
    Update frequency HT AN field and histograms with the correct AN from the AN HT.

    This module updates the freq_ht AN field which is an array of ints with the correct
    AN from the AN HT integer array annotation. It uses the freq_correction_ht
    dictionary and the freq_ht dictionary to match the indexing of array annotations so
    each group is correctly updated. It then recomputes AF, AC/AN.

    :param freq_ht: Frequency HT to update.
    :param freq_correction_ht: HT to update AN and histograms from.
    :return: Updated frequency HT.
    """
    freq_correction_meta = freq_correction_ht.index_globals().strata_meta

    # Set all freq_ht's ANs to 0 so can use the merge_freq function to update
    # the ANs (expects int32s).
    freq_ht = freq_ht.annotate(
        freq=freq_ht.freq.map(
            lambda x: x.select(
                AC=x.AC, AN=hl.int32(0), homozygote_count=x.homozygote_count, AF=x.AF
            )
        )
    )

    # Create freq array annotation of structs with AC and homozygote_count set to 0 to
    # use merge_freq function (expects int32s).
    freq_correction_ht = freq_correction_ht.transmute(
        freq=freq_correction_ht.AN.map(
            lambda x: hl.struct(
                AC=0, AN=hl.int32(x), homozygote_count=0, AF=hl.missing(hl.tfloat64)
            )
        )
    )
    freq_ht = freq_ht.annotate(
        freq_correction_freq=freq_correction_ht[freq_ht.key].freq
    )
    freq_ht = freq_ht.annotate_globals(freq_correction_strata_meta=freq_correction_meta)

    logger.info("Merging freq arrays from v4 and freq correction HT...")
    freq, freq_meta = merge_freq_arrays(
        [freq_ht.freq, freq_ht.freq_correction_freq],
        [freq_ht.freq_meta, freq_ht.freq_correction_strata_meta],
    )
    freq_ht = freq_ht.annotate(freq=freq)
    freq_ht = freq_ht.annotate_globals(freq_meta=freq_meta)
    freq_ht = freq_ht.drop("freq_correction_freq", "freq_correction_strata_meta")

    logger.info("Updating freq_ht with corrected GQ and DP histograms...")
    corrected_index = freq_correction_ht[freq_ht.key].qual_hists
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
    test_chr22_chrx_chry = args.test_chr22_chrx_chry
    test_n_partitions = args.test_n_partitions
    test = use_test_dataset or test_n_partitions or test_chr22_chrx_chry
    af_threshold = args.af_threshold

    hl.init(
        log="/frequency_an_update.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    hl._set_flags(use_ssa_logs="1")

    freq_ht = get_freq("4.0").ht()
    try:
        if (
            args.regenerate_ref_an_hists
            or args.calculate_var_an_hists_adjustment
            or args.calculate_all_sites_het_nonref_adjustment
        ):
            logger.info("Regenerating reference AN and GQ/DP histograms...")

            ref_ht = vep_context.versions["105"].ht()
            ref_ht = ref_ht.key_by("locus").select().distinct()

            if test_chr22_chrx_chry:
                chrom = ["chr22", "chrX", "chrY"]
                ref_ht = hl.filter_intervals(
                    ref_ht, [hl.parse_locus_interval(c) for c in chrom]
                )

            vds = get_gnomad_v4_vds(
                release_only=True,
                test=use_test_dataset or test_chr22_chrx_chry,
                filter_partitions=range(2) if test_n_partitions else None,
                annotate_meta=True,
                chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
            )

            group_membership_ht = get_exomes_group_membership_ht(
                meta().ht(),
                get_downsampling().ht(),
                get_downsampling(subset="non_ukb").ht(),
            ).checkpoint(hl.utils.new_temp_file("group_membership", "ht"))

            interval_ht = adjust_interval_padding(
                calling_intervals(
                    interval_name="union", calling_interval_padding=0
                ).ht(),
                150,
            ).checkpoint(hl.utils.new_temp_file("interval", "ht"))

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

            if args.regenerate_ref_an_hists:
                compute_allele_number_per_ref_site(
                    prep_vds_for_all_sites_stats(vds),
                    ref_ht,
                    interval_ht,
                    group_membership_ht,
                ).naive_coalesce(10000).write(
                    f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/all_sites_an_before_adjustment.ht",
                    overwrite=overwrite,
                )

            if args.calculate_var_an_hists_adjustment:
                logger.info(
                    "Calculating adjustments for all sites AN and GQ/DP histograms..."
                )
                compute_an_and_hists_het_fail_adj_ab(
                    vds,
                    group_membership_ht,
                    freq_ht,
                ).write(
                    f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/het_fail_adj_ab_ht.ht",
                    overwrite=overwrite,
                )

            if args.calculate_all_sites_het_nonref_adjustment:
                logger.info("Calculating all sites het non ref adjustment...")
                get_all_sites_nonref_adjustment(
                    vds,
                    group_membership_ht,
                ).write(
                    f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/het_nonref_adjust.ht",
                    overwrite=overwrite,
                )

        if args.adjust_all_sites_an_and_hists:
            an_ht = hl.read_table(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/all_sites_an_before_adjustment.ht",
            )
            het_fail_adj_ab_ht = hl.read_table(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/het_fail_adj_ab_ht.ht",
            )
            het_nonref_ht = hl.read_table(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/het_nonref_adjust.ht",
            )
            an_ht = adjust_all_sites_an_and_hists(
                an_ht, het_fail_adj_ab_ht, het_nonref_ht
            ).checkpoint(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/all_sites_an_after_adjustment.ht",
                overwrite=overwrite,
            )
            an_ht = an_ht.select("AN", "qual_hists").checkpoint(
                get_all_sites_an_and_qual_hists(test=test).path,
                overwrite=overwrite,
            )
            an_ht.select("AN").write(
                release_all_sites_an(public=False, test=test).path,
                overwrite=overwrite,
            )

        if args.adjust_freq_an_hists:
            an_ht = hl.read_table(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/all_sites_an_before_adjustment.ht",
            )
            het_fail_adj_ab_ht = hl.read_table(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/het_fail_adj_ab_ht.ht",
            )

            freq_correction_ht = adjust_per_site_an_and_hists_for_frequency(
                an_ht, freq_ht, het_fail_adj_ab_ht
            )
            freq_correction_ht.write(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/freq_correction.ht",
                overwrite=overwrite,
            )
        if args.update_freq_an_and_hists:
            logger.info("Updating AN in freq HT...")

            freq_correction_ht = hl.read_table(
                f"gs://gnomad{'-tmp-4day' if test else ''}/v4.1/temp/frequency_fix/freq_correction.ht"
            )
            # The new correction HT has used "gen_anc" as the pop key from the
            # start so replacing it with "pop" here to match the freq HT. Updating the
            # freq_ht to "gen_anc" is the appropriate fix but it requires far more code
            # updates to the imported functions of generate_freq.
            freq_correction_ht = freq_correction_ht.annotate_globals(
                strata_meta=freq_correction_ht.strata_meta.map(
                    lambda d: hl.dict(
                        d.items().map(
                            lambda x: hl.if_else(x[0] == "gen_anc", ("pop", x[1]), x)
                        )
                    )
                )
            )
            freq_ht = get_freq(
                version="4.0",
                hom_alt_adjusted=False,
                finalized=False,
            ).ht()

            if test:
                freq_ht = freq_ht.filter(hl.is_defined(freq_correction_ht[freq_ht.key]))

            freq_ht = drop_gatk_groupings(freq_ht)
            ht = update_freq_an_and_hists(freq_ht, freq_correction_ht)

            ht.write(
                get_freq(
                    version="4.1",
                    test=test,
                    hom_alt_adjusted=False,
                    finalized=False,
                ).path,
                overwrite=overwrite,
            )

        if args.update_af_dependent_anns:
            logger.info("Updating AF dependent annotations...")

            ht = get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=False,
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
                finalized=False,
            ).ht()
            freq_ht = create_final_freq_ht(freq_ht)

            logger.info("Final frequency HT schema...")
            freq_ht.describe()
            freq_ht.write(
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
        "--test-chr22-chrx-chry",
        help=(
            "Whether to run a test using only the chr22, chrX, and chrY chromosomes of"
            " the VDS test dataset."
        ),
        action="store_true",
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
        "--calculate-var-an-hists-adjustment",
        help=(
            "Calculate adjustments that need to be made to the all sites AN and GQ/DP "
            "histograms in order to be used for the multi-allelic split frequency "
            "Table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--calculate-all-sites-het-nonref-adjustment",
        help=(
            "Calculate the number of het non-ref genotypes that are being counted "
            "twice so the all sites AN and histograms can be adjusted correctly."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--adjust-all-sites-an-and-hists",
        help="Adjust the all sites AN and GQ/DP histograms for the frequency Table.",
        action="store_true",
    )
    parser.add_argument(
        "--adjust-freq-an-hists",
        help="Adjust the all sites AN and GQ/DP histograms for the frequency Table.",
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
