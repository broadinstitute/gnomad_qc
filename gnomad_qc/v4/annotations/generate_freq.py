"""Script to generate the frequency data annotations across v4 exomes."""  # TODO: Expand on script description once nearing completion.
import argparse
import logging
from typing import (  # TODO: REMOVE WITH ANNOTATE_FREQ
    Any,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: subsets will be changed to UKB, non-UKB,
)
from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    POPS,
    POPS_TO_REMOVE_FOR_POPMAX,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_adj,
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
    gatk_expr: hl.expr.StringExpression,
    ab_cutoff: hl.float = 0.9,
    gatk_versions_to_fix: hl.expr.SetExpression = hl.set(["4.0.10.1", "4.1.0.0"]),
) -> hl.MatrixTable:
    """
    Annotate high AB hets by group membership.

    :param mt: MatrixTable to annotate high AB hets onto.
    :param gatk_expr: GATK version expression.
    :param ab_cutoff: Allele balance cutoff for hom alt depletion fix. Defaults to 0.9
    :param gatk_versions_to_fix: GATK versions that need the hom alt depletion fix. Defaults to hl.set(["4.0.10.1", "4.1.0.0"])
    :return: _description_
    """
    logger.info("Annotating number of high AB het sites in each freq group...")
    mt = mt.annotate_rows(
        high_ab_hets_by_group_membership=hl.agg.array_agg(
            lambda i: hl.agg.filter(
                mt.group_membership[i]
                & needs_high_ab_het_fix_expr(
                    mt, gatk_expr, ab_cutoff, gatk_versions_to_fix
                ),
                hl.agg.count(),
            ),
            hl.range(hl.len(mt.group_membership)),
        )
    )
    return mt


def needs_high_ab_het_fix_expr(
    mt: hl.MatrixTable,
    gatk_expr: hl.expr.StringExpression,
    ab_cutoff: hl.float = 0.9,
    gatk_versions_to_fix: hl.expr.SetExpression = hl.set(["4.0.10.1", "4.1.0.0"]),
) -> hl.MatrixTable:
    """
    Annotate the MatrixTable rows with the count of high AB het genotypes per frequency sample grouping.

    This arrary allows the AC to be corrected with the homalt depletion fix.

    :param mt: Hail MatrixTable to annotate.
    :param gatk_expr: Expression indicating which version of GATK a sample's gVCF was called with.
    :param ab_cutoff: Alelle balance threshold at which the GT changes to homalt. Default is 0.9.
    :param gatk_versions_to_fix: Set of GATK versions impacted by hom alt depletion. Default is {"4.0.10.1", "4.1.0.0"}.
    :return mt:
    """
    return (
        (mt.AD[1] / mt.DP > ab_cutoff)
        & mt.adj
        & mt.GT.is_het_ref()
        & (gatk_versions_to_fix.contains(gatk_expr))
        & ~mt._het_non_ref  # Skip adjusting genotypes if sample originally had a het nonref genotype
    )


def annotate_gatk_version(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate MatrixTable's samples with HaplotypeCaller's GATK version. GATK hom alt depletion fix has been in since GATK version 4.1.4.1.

    :param mt: Hail MatrixTable to annotate.
    :return mt:
    """
    meta_ht = meta.ht()
    gatk_ht = hl.read_table(
        "gs://gnomad-mwilson/v4/homalt_depletion/v4_sample_htc_versions.ht"
    ).key_by(
        "s"
    )  # TODO: Discuss if this should be a resource or annotation on metadata, misssing GATK version info for UKB samples
    mt = mt.annotate_cols(
        gatk_version=hl.case()
        .when(
            hl.is_defined(gatk_ht[mt.col_key].gatk_version),
            gatk_ht[mt.col_key].gatk_version,
        )
        .when(hl.is_defined(meta_ht[mt.col_key].project_meta.ukb_meta.ukb_batch), "ukb")
        .or_missing()
    )
    return mt


def subtract_high_ab_hets_from_ac(ht: hl.Table, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct frequencies at sites with an AF greater than the af_threshold.

    :param ht: Hail Table containing freq and high_ab_het annotations.
    :param af_threshold: AF threshold at which to correct frequency. Default is 0.01.
    :return: Hail Table with adjusted frequencies.
    """
    ht = ht.annotate(
        ab_adjusted_freq=hl.if_else(
            ht.freq[0].AF < af_threshold,
            ht.freq,
            hl.map(
                lambda f, g: hl.struct(
                    AC=hl.int32(f.AC - g),
                    AF=hl.if_else(f.AN > 0, (f.AC - g) / f.AN, hl.missing(hl.tfloat64)),
                    AN=f.AN,
                    homozygote_count=f.homozygote_count,
                ),
                ht.freq,
                ht.high_ab_hets_by_group_membership,
            ),
        )
    )

    return ht


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
            vds = hl.vds.filter_chromosomes(vds, keep=f"chr{chrom}")

    logger.info("Annotating non_ref hets pre-split...")
    vds = annotate_non_ref_het(vds)

    # if args.subsets:
    #     vds = hl.vds.filter_samples(
    #         vds,
    #     )  # TODO: Where is subset information coming meta.project_meta.ukb is defined for non-ukb and ukd, how do we determine non-topmed v3.1 meta did (topmed=ht.s.startswith("NWD"))

    logger.info("Splitting VDS....")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Annotating GATK version for frequency groupings....")
    mt = annotate_gatk_version(mt)

    logger.info("Computing adj and sex adjusted genotypes...")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus, mt.GT, mt.meta.sex_imputation.sex_karyotype
        ),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    final_fields = []
    if args.get_freq_and_high_ab:
        logger.info("Annotating frequencies...")
        mt = annotate_freq(
            mt,
            sex_expr=mt.meta.sex_imputation.sex_karyotype,
            pop_expr=mt.meta.population_inference.pop,
            # downsamplings=DOWNSAMPLINGS["v4"],
            additional_strata_expr={"gatk_version": mt.gatk_version},
            additional_strata_grouping_expr={"pop": mt.meta.population_inference.pop},
        )
        mt = annotate_high_ab_hets_by_group_membership(mt, gatk_expr=mt.gatk_version)
        ht = mt.rows()
        final_fields = ["freq", "high_ab_hets_by_group_membership"]

    if adjust_freqs:
        # TODO: This should apply fix given the AF threshold and adjusts the
        # frequencies using the list created above
        ht = subtract_high_ab_hets_from_ac(ht, af_threshold)
        final_fields.append("ab_adjusted_freq")

    logger.info("Writing frequency table...")
    ht = ht.select(*final_fields)
    ht = ht.write(
        get_freq(hom_alt_adjustment=adjust_freqs, test=test).path,
        overwrite=args.overwrite,
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

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)