"""Script to generate the frequency data annotations across v4 exomes."""  # TODO: Expand on script description once nearing completion.
import argparse
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: subsets will be changed to UKB, non-UKB,
)
from gnomad.resources.grch38.gnomad import (  # TODO: Update these values in gnomad_methods for v4
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
from gnomad_qc.v3.utils import (  # TODO: Move this to gnomad_methods? Note impacts pre-4.1.4.1 GATK generated gVCFs
    hom_alt_depletion_fix,
)
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


def annotate_homalt_gt_change(
    mt: hl.MatrixTable,
    ab_cutoff: hl.tfloat = 0.9,
) -> hl.MatrixTable:
    """
    Annotate the MatrixTable rows with the count of genotypes which would change if the homalt depletion fix was implemented.

    :param mt: Hail MatrixTable to annotate.
    :param ab_cutoff: Alelle balance threshold at which the GT changes to homalt. Default is 0.9.
    :return mt:
    """
    logger.info("Annotating high AB hets...")
    mt = mt.annotate_entries(
        het_high_ab=(
            mt.GT.is_het()
            # Skip adjusting genotypes if sample originally had a het nonref genotype
            & ~mt.GT.is_het_non_ref()
            & (mt.AD[1] / mt.DP > ab_cutoff)
        )
    )
    return mt


def annotate_gatk_version(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate MatrixTable's samples with HaplotypeCaller's GATK version. GATK hom alt depletion fix has been in since GATK version 4.1.4.1.

    :param mt: Hail MatrixTable to annotate.
    :return mt:
    """
    logger.info("Annotating MT with GATK version if it exists....")
    meta_ht = meta.ht()
    gatk_ht = hl.read_table(
        "gs://gnomad-mwilson/v4/homalt_depletion/v4_sample_htc_versions.ht"
    ).key_by(
        "s"
    )  # TODO: Discuss if this should be a resource or annotation on metadata, misssing GATK version info for UKB samples
    mt = mt.annotate_cols(
        gatk_version=hl.case()
        .when(hl.is_defined(mt.gatk_version), mt.gatk_version)
        .when(hl.is_defined(meta_ht[mt.col_key].project_meta.ukb_meta.ukb_batch), "ukb")
        .or_missing()
    )
    logger.info(mt.aggregate_cols(hl.agg.counter(mt.gatk_version)))  # Validity check
    return mt


def main(args):  # noqa: D103
    subsets = (
        args.subsets
    )  # TODO: Determine if splitting subset freq from whole callset agg
    hl.init(
        log=f"/generate_frequency_data{'.' + '_'.join(subsets) if subsets else ''}.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    vds = get_gnomad_v4_vds(
        split=True, test=args.test
    )  # TODO: Update to release = True once outlier filtering is complete
    meta_ht = meta.ht()

    if args.test:
        logger.info("Filtering to DRD2 in test VDS for testing purposes...")
        test_interval = [
            hl.parse_locus_interval("chr11:113409605-113475691")
        ]  # NOTE: test on gene DRD2 for now, will update to chr22 when original PR goes in
        vds = hl.vds.filter_intervals(vds, test_interval)
        variants, samples = vds.variant_data.count()
        logger.info(
            "Test VDS has %s variants in DRD2 in %s samples...", variants, samples
        )

    if args.calculate_gatk_af_diff:
        logger.info(
            "Densifying VDS and determining frequency differences across GATK"
            " versions..."
        )
        mt = hl.vds.to_dense_mt(vds)
        mt = annotate_gatk_version(mt)

        # Make this a module that can be recalled in the main for all chrs
        mt = mt.annotate_cols(
            sex_karyotype=meta_ht[mt.col_key].sex_imputation.sex_karyotype,
            pop=meta_ht[mt.col_key].population_inference.pop,
        )

        logger.info("Computing adj and sex adjusted genotypes...")
        mt = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(  # NOTE: Not needed just yet but we will
                mt.locus, mt.GT, mt.sex_karyotype
            ),
        )
        mt = annotate_adj(mt)
        logger.info("Annotating frequency...")
        mt = annotate_freq(
            mt,
            sex_expr=mt.sex_karyotype,
            pop_expr=mt.pop,
            additional_strata_expr={
                "gatk_version": mt.gatk_version
            },  # TODO: Flip to version and update gnomad_methods to expand this
        )
        # TODO: Add annotation for "if you applied homalt hot fix for this
        # variant, how many genotypes would change"

        mt = annotate_homalt_gt_change(mt)

    logger.info("Checkpointing MatrixTable...")
    mt = mt.checkpoint("gs://gnomad-tmp-4day/freq/test_freq_aggs.mt", overwrite=True)
    mt = mt.annotate_rows(high_ab_hets=hl.agg.count_where(mt.het_high_ab))
    mt = mt.checkpoint("gs://gnomad-tmp-4day/freq/test_freq_aggs2.mt", overwrite=True)

    # Validity checking
    mt.rows().show()
    mt.cols().show()
    print(mt.aggregate_cols(hl.agg.counter(mt.gatk_fix)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the MT.", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--subsets", help="Subsets to run frequency calculation on.", choices=SUBSETS
    )
    parser.add_argument(
        "--calculate-gatk-af-diff",
        help="Determine AF differences between pre and post GATK version 4.1.4.1.",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
