"""Script to generate the frequency data annotations across v4 exomes."""  # TODO: Expand on script description once nearing completion.
import argparse
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: subsets will be changed to UKB, non-UKB, non-TopMed
)
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq_and_high_ab_hets,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict
from gnomad.utils.slack import slack_notifications

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
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


def get_freq_resources(
    overwrite: bool, test: bool, chr: str
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
    get_freq_and_high_ab = PipelineStepResourceCollection(
        "--get-freq-and-high-ab",
        output_resources={
            "freq_high_ab_ht": get_freq(test=test, hom_alt_adjustment=False, chr=chr),
        },
    )
    correct_call_stats = PipelineStepResourceCollection(
        "--correct-call-stats",
        pipeline_input_steps=[get_freq_and_high_ab],
        output_resources={
            "freq_ht": get_freq(test=test, hom_alt_adjustment=True, chr=chr),
        },
    )
    freq_pipeline.add_steps(
        {
            "get_freq_and_high_ab": get_freq_and_high_ab,
            "correct_call_stats": correct_call_stats,
        }
    )
    return freq_pipeline


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

    :param mt: Input MT.
    :return: MatrixTable with qual histogram annotations.
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
        popmax=pop_max_expr(
            ht.ab_adjusted_freq, ht.freq_meta, POPS_TO_REMOVE_FOR_POPMAX
        ),
    )
    ht = ht.annotate_globals(
        faf_meta=faf_meta,
        faf_index_dict=make_faf_index_dict(faf_meta, label_delimiter="-"),
    )
    ht = ht.annotate(
        popmax=ht.popmax.annotate(
            faf95=ht.faf[
                ht.faf_meta.index(lambda x: x.values() == ["adj", ht.popmax.pop])
            ].faf95
        )
    )
    return ht


def main(args):  # noqa: D103
    subsets = args.subsets
    test_dataset = args.test_dataset
    test_n_partitions = args.test_n_partitions
    test = test_dataset or test_n_partitions
    chrom = args.chrom

    ab_cutoff = args.ab_cutoff
    af_threshold = args.af_threshold

    correct_callstats = args.adjust_callstats

    hl.init(
        log=f"/generate_frequency_data{'.' + '_'.join(subsets) if subsets else ''}.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # TODO: Determine if splitting subset freq from whole callset agg
    resources = get_freq_resources(args.overwrite, test, chrom)

    vds = get_gnomad_v4_vds(
        test=test_dataset, release_only=True, chrom=chrom, n_partitions=10000
    )
    meta_ht = meta.ht()
    final_anns = {}

    logger.info("Adding metadata to VDS variant data cols...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.annotate_cols(meta=meta_ht[vds.variant_data.col_key]),
    )

    if test or chrom:
        if test_dataset:
            logger.info("Filtering to DRD2 in test VDS for testing purposes...")
            test_interval = [
                hl.parse_locus_interval("chr11:113409605-113475691")
            ]  # NOTE: test on gene DRD2 for now, will update to chr22 when original PR goes in
            vds = hl.vds.filter_intervals(vds, test_interval)
            variants, samples = vds.variant_data.count()
            logger.info(
                "Test VDS has %s variants in DRD2 in %s samples...", variants, samples
            )
        elif test_n_partitions:
            test_vd = vds.variant_data._filter_partitions(range(test_n_partitions))
            test_rd = vds.reference_data._filter_partitions(range(test_n_partitions))
            vds = hl.vds.VariantDataset(test_rd, test_vd)
        else:
            logger.info("Filtering to chromosome %s...", chrom)
            vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    logger.info("Annotating non_ref hets pre-split...")
    vds = annotate_non_ref_het(vds)

    if args.subsets:
        vds = hl.vds.filter_samples(
            vds,
            vds.variant_data.meta.subsets.contains(
                hl.literal(subsets)
            ),  # TODO: Fix this logic and all subset logic
        )
    logger.info("Splitting VDS....")  # TODO: Move this to get_gnomad_v4_vds
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.select_entries("_het_non_ref", "AD", "DP", "GQ", "GT"),
    )

    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)  # TODO: Move this to get_gnomad_v4_vds

    logger.info("Computing adj and sex adjusted genotypes...")
    mt = mt.select_entries(
        "_het_non_ref",
        "AD",
        "DP",
        GT=adjusted_sex_ploidy_expr(
            mt.locus, mt.GT, mt.meta.sex_imputation.sex_karyotype
        ),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    mt = mt.checkpoint(
        "gs://gnomad-tmp/julia/dense_mt_chr20_10k_partitions.mt", overwrite=True
    )

    if args.get_freq_and_high_ab:
        logger.info("Annotating frequencies and counting high AB het calls...")
        res = resources.get_freq_and_high_ab
        freq_ht = annotate_freq_and_high_ab_hets(
            mt,
            sex_expr=mt.meta.sex_imputation.sex_karyotype,
            pop_expr=mt.meta.population_inference.pop,
            downsamplings=DOWNSAMPLINGS["v4"],
            additional_strata_expr={"gatk_version": mt.meta.project_meta.gatk_version},
            additional_strata_grouping_expr={"pop": mt.meta.population_inference.pop},
            ab_cutoff=ab_cutoff,
        )
        final_anns = {}
        if args.calculate_qual_hists:
            logger.info("Annotating quality metrics histograms...")
            mt = annotate_quality_metrics_hist(mt)
            final_anns.update(
                {
                    "qual_hists": mt[freq_ht.key].qual_hists,
                    "raw_qual_hists": mt[freq_ht.key].raw_qual_hists,
                }
            )
        if args.calculate_age_hists:
            logger.info("Computing age histograms for each variant...")
            mt = compute_age_hist(mt)
            final_anns.update(
                {"age_hist_het": mt.age_hist_het, "age_hist_hom": mt.age_hist_hom}
            )

        freq_ht = freq_ht.annotate(**final_anns)
        freq_ht.write(res.freq_high_ab_ht.path, overwrite=args.overwrite)

    if correct_callstats:  # TOOD: Update with read freq so script can start here
        logger.info("Adjusting frequencies by accounting for high AB hets...")
        res = resources.correct_call_stats
        freq_ht = correct_call_stats(freq_ht, af_threshold)
        final_anns = {"ab_adjusted_freq": freq_ht[mt.row_key].ab_adjusted_freq}

        if args.calculate_inbreeding_coeff:
            logger.info("Calculating InbreedingCoeff...")
            # NOTE: This is not the ideal location to calculate this, but added here to avoid another densify # noqa
            mt = mt.annotate_rows(
                InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT)
            )
            final_anns["InbreedingCoeff"] = mt.InbreedingCoeff

        if args.correct_qual_hists:
            # TODO: Add step to correct qual hists (ab based), should subtract high AB
            # het count from above 0.9 bin
            pass

        if args.faf_popmax:
            logger.info("computing FAF & popmax...")
            mt = generate_faf_popmax(mt)
            final_anns.update({"faf": mt.faf, "popmax": mt.popmax})

        logger.info("Writing frequency table...")
        mt.describe()
        logger.info(f"{final_anns.keys()} are the final annotations")
        ht = mt.select_rows(**final_anns).rows()
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
        "--correct-call-stats",
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
