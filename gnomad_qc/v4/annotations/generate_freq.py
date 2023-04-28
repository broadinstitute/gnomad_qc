"""Script to generate the frequency data annotations across v4 exomes."""  # TODO: Expand on script description once nearing completion.
import argparse
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: subsets will be changed to UKB, non-UKB, non-TopMed
)
from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    POPS,
    POPS_TO_REMOVE_FOR_POPMAX,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import (  # TODO: Need to generate this resource for raw callstats
    get_freq,
)
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    qc_temp_prefix,
)
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v4.resources.variant_qc import (  # TODO: Confirm we want this in frequency calculations
    SYNDIP,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency")
logger.setLevel(logging.INFO)


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


def annotate_high_ab_hets_by_group_membership(
    mt: hl.MatrixTable,
    ab_cutoff: hl.float = 0.9,
) -> hl.MatrixTable:
    """
    Annotate high AB hets by group membership.

    :param mt: MatrixTable to annotate high AB hets onto.
    :param ab_cutoff: Allele balance cutoff for hom alt depletion fix. Defaults to 0.9
    :return: _description_
    """
    logger.info("Annotating number of high AB het sites in each freq group...")
    mt = mt.annotate_rows(
        high_ab_hets_by_group_membership=hl.agg.array_agg(
            lambda i: hl.agg.filter(
                mt.group_membership[i]
                & needs_high_ab_het_fix_expr(
                    mt,
                    ab_cutoff,
                ),
                hl.agg.count(),
            ),
            hl.range(hl.len(mt.group_membership)),
        )
    )
    return mt


def needs_high_ab_het_fix_expr(
    mt: hl.MatrixTable,
    ab_cutoff: hl.float = 0.9,
) -> hl.MatrixTable:
    """
    Annotate the MatrixTable rows with the count of high AB het genotypes per frequency sample grouping.

    This arrary allows the AC to be corrected with the homalt depletion fix.

    :param mt: Hail MatrixTable to annotate.
    :param ab_cutoff: Alelle balance threshold at which the GT changes to homalt. Default is 0.9.
    :return mt:
    """
    return (
        (mt.AD[1] / mt.DP > ab_cutoff)
        & mt.adj
        & mt.GT.is_het_ref()
        & ~mt.meta.project_meta.fixed_homalt_model
        & ~mt._het_non_ref  # Skip adjusting genotypes if sample originally had a het nonref genotype
    )


def correct_call_stats(mt: hl.MatrixTable, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct frequencies at sites with an AF greater than the af_threshold.

    :param ht: Hail Table containing freq and high_ab_het annotations.
    :param af_threshold: AF threshold at which to correct frequency. Default is 0.01.
    :return: Hail Table with adjusted frequencies.
    """
    mt = mt.annotate_rows(
        ab_adjusted_freq=hl.if_else(
            mt.freq[0].AF > af_threshold,
            mt.freq,
            hl.map(
                lambda f, g: hl.struct(
                    AC=hl.int32(f.AC + g),
                    AF=hl.if_else(f.AN > 0, (f.AC + g) / f.AN, hl.missing(hl.tfloat64)),
                    AN=f.AN,
                    homozygote_count=f.homozygote_count + g,
                ),
                mt.freq,
                mt.high_ab_hets_by_group_membership,
            ),
        )
    )

    return mt


# TODO: not sure if this step will be very expensive, may need to combine with functions
def set_high_ab_het_to_hom_alt(
    mt: hl.MatrixTable,
    ab_cutoff: hl.float = 0.9,
    af_cutoff: float = 0.01,
) -> hl.MatrixTable:
    """
    Adjust MT genotypes with temporary fix for the depletion of homozygous alternate genotypes.

    More details about the problem can be found on the gnomAD blog:
    https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#tweaks-and-updates

    :param mt: Input MT that needs hom alt genotype fix
    :param af_cutoff: Allele frequency cutoff for variants that need the hom alt fix. Default is 0.01
    :param ab_cutoff: Allele balance cutoff to determine which genotypes need the hom alt fix. Default is 0.9
    :return: MatrixTable with genotypes adjusted for the hom alt depletion fix
    """
    return mt.annotate_entries(
        GT=hl.if_else(
            needs_high_ab_het_fix_expr(mt, ab_cutoff)
            & (mt.ab_adjusted_freq[0].AF > af_cutoff),
            hl.call(1, 1),
            mt.GT,
        )
    )


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


def annotate_quality_metrics_hist(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate quality metrics histograms.

    :param mt: Input MT with age annotation.
    :return: MatrixTable with age histogram annotations.
    """
    mt = mt.annotate_rows(qual_hists=qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD, mt.adj))
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


def generate_faf_popmax(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute filtering allele frequencies and popmax with the AB-adjusted frequencies.

    :param mt: Hail Table containing freq, ab_adjusted_freq, high_ab_het annotations.
    :return: Hail MatirxTable with faf & popmax annotations.
    """
    faf, faf_meta = faf_expr(
        mt.ab_adjusted_freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX
    )
    mt = mt.annotate_rows(
        faf=faf,
        popmax=pop_max_expr(
            mt.ab_adjusted_freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX
        ),
    )
    mt = mt.annotate_globals(
        faf_meta=faf_meta,
        faf_index_dict=make_faf_index_dict(faf_meta, label_delimiter="-"),
    )
    mt = mt.annotate_rows(
        popmax=mt.popmax.annotate(
            faf95=mt.faf[
                mt.faf_meta.index(lambda x: x.values() == ["adj", mt.popmax.pop])
            ].faf95
        )
    )
    return mt


def main(args):  # noqa: D103
    # TODO: Determine if splitting subset freq from whole callset agg
    subsets = args.subsets
    test = args.test
    chrom = args.chrom
    af_threshold = args.af_threshold
    adjust_freqs = args.adjust_freqs

    hl.init(
        log=f"/generate_frequency_data{'.' + '_'.join(subsets) if subsets else ''}.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # TODO: Update to release = True once outlier filtering is complete,
    # possibly sample_meta=True if added
    vds = get_gnomad_v4_vds(test=test)
    meta_ht = meta.ht()
    final_anns = []

    logger.info("Adding metadata to VDS variant data cols...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.annotate_cols(meta=meta_ht[vds.variant_data.col_key]),
    )

    if test or chrom:
        if test:
            logger.info("Filtering to DRD2 in test VDS for testing purposes...")
            test_interval = [
                hl.parse_locus_interval("chr11:113409605-113475691")
            ]  # NOTE: test on gene DRD2 for now, will update to chr22 when original PR goes in
            vds = hl.vds.filter_intervals(vds, test_interval)
            variants, samples = vds.variant_data.count()
            logger.info(
                "Test VDS has %s variants in DRD2 in %s samples...", variants, samples
            )
        else:
            logger.info("Filtering to chromosome %s...")
            vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    logger.info("Annotating non_ref hets pre-split...")
    vds = annotate_non_ref_het(vds)

    # if args.subsets:
    #     vds = hl.vds.filter_samples(
    #         vds,
    #     )  # TODO: Where is subset information coming meta.project_meta.ukb is defined for non-ukb and ukd, how do we determine non-topmed v3.1 meta did (topmed=ht.s.startswith("NWD"))

    logger.info("Splitting VDS....")  # TODO: Move this to get_gnomad_v4_vds
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)  # TODO: Move this to get_gnomad_v4_vds

    logger.info("Computing adj and sex adjusted genotypes...")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus, mt.GT, mt.meta.sex_imputation.sex_karyotype
        ),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    if args.get_freq_and_high_ab:
        logger.info("Annotating frequencies...")
        mt = annotate_freq(
            mt,
            sex_expr=mt.meta.sex_imputation.sex_karyotype,
            pop_expr=mt.meta.population_inference.pop,
            # downsamplings=DOWNSAMPLINGS["v4"],
            additional_strata_expr={"gatk_version": mt.meta.project_meta.gatk_version},
            additional_strata_grouping_expr={"pop": mt.meta.population_inference.pop},
        )
        mt = annotate_high_ab_hets_by_group_membership(mt)
        final_anns.extend(["freq", "high_ab_hets_by_group_membership"])

    if adjust_freqs:
        logger.info("Adjusting frequencies by accounting for high AB hets...")
        mt = correct_call_stats(mt, af_threshold)
        final_anns.extend(
            ["ab_adjusted_freq"]
        )  # NOTE: Do we want to keep original freqs? If not, overwrite the freq ann

    if (
        args.set_high_ab_het_to_hom_alt
    ):  # NOTE: MW: We want to avoid this if we can but to avoid it will need to change age hists and inbreeding coefficient methods
        logger.info(
            "Setting het genotypes at sites with >1% AF (using adjusted frequencies)"
            " and > 0.9 AB to homalt..."  # TODO: Update AF threshold once we analyze
        )
        # using the adjusted allele frequencies to fix the het to hom alt
        mt = set_high_ab_het_to_hom_alt(mt)

    if args.calculate_inbreeding_coeff:
        logger.info("Calculating InbreedingCoeff...")
        # NOTE: This is not the ideal location to calculate this, but added here to avoid another densify # noqa
        mt = mt.annotate_rows(InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT))
        final_anns.extend(["InbreedingCoeff"])

    if args.calculate_hists:
        logger.info("Computing age histograms for each variant...")
        mt = compute_age_hist(mt)
        final_anns.extend(["age_hist_het", "age_hist_hom"])

        logger.info("Annotating quality metrics histograms...")
        mt = annotate_quality_metrics_hist(mt)
        final_anns.extend(["qual_hists", "raw_qual_hists"])

    if args.faf_popmax:
        logger.info("computing FAF & popmax...")
        mt = generate_faf_popmax(mt)
        final_anns.extend(["faf", "popmax"])

    logger.info("Writing frequency table...")
    mt.describe()
    logger.info(f"{final_anns} are the final annotations")
    ht = mt.rows()
    ht = ht.select_rows(*final_anns)
    ht = ht.write(
        get_freq(test=test, hom_alt_adjustment=adjust_freqs, chr=chrom).path,
        overwrite=True,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the MT.", action="store_true"
    )
    parser.add_argument(
        "--chrom",
        help="If passed, script will only run on passed chromosome.",
        type=int,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--subsets",
        help="Subsets to run frequency calculation on.",
        choices=SUBSETS["v4"],
    )
    parser.add_argument(
        "--get-freq-and-high-ab",
        help="Calculate frequencies and high AB sites per frequency grouping.",
        action="store_true",
    )
    parser.add_argument(
        "--adjust-freqs",
        help=(
            "Adjust each frequency entry to account for homozygous alternate depletion"
            " present in GATK versions released prior to 4.1.4.1."
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
        "--set-high-ab-het-to-hom-alt",
        help="Set high AB hets to hom alt.",
        action="store_true",
    )
    parser.add_argument(
        "--calculate-inbreeding-coeff",
        help="Calculate inbreeding coefficient.",
        action="store_true",
    )
    parser.add_argument(
        "--calculate-hists",
        help="Calculate age histograms and quality metrics histograms.",
        action="store_true",
    )
    parser.add_argument(
        "--faf-popmax",
        help=(
            "compute filtering allele frequency and population that has the highest AF "
            "based on the adjusted allele frequency"
        ),
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
