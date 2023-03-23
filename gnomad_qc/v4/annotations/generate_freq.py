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


def annotate_gatk_version(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate MAtrixTable's samples with HaplotypeCaller's GATK version.

    :param mt: Hail MatrixTable to annotate.
    :return mt:
    """
    logger.info("Annotating MT with GATK version if it exists....")
    meta_ht = meta.ht()
    gatk_ht = hl.import_table(
        "gs://gnomad-tmp-4day/batch/mwilson/v4_sample_htc_versions_322.txt"
    ).key_by(
        "s"
    )  # TODO: Discuss if this should be a resource or annotation on metadata, misssing GATK version info for UKB samples
    fixed_versions = hl.set(["4.1.4.1", "4.1.8.0"])
    mt = mt.annotate_cols(
        gatk_version=hl.or_missing(
            hl.is_defined(gatk_ht[mt.col_key].gatk_version),
            gatk_ht[mt.col_key].gatk_version,
        )
    )
    mt = mt.annotate_cols(
        gatk_fix=hl.case()
        .when(fixed_versions.contains(mt.gatk_version), "has-gatk-fix")
        .when(~fixed_versions.contains(mt.gatk_version), "no-gatk-fix")
        .when(hl.is_defined(meta_ht[mt.col_key].project_meta.ukb_meta.ukb_batch), "ukb")
        .or_missing()
    )
    return mt
    # Ok cool so the annotate_freq function wants the entries to be annotated with "adj" -- this is the high quality GT I mentioned
    # And we don't have that yet, right?


def main(args):  # noqa: D103
    subsets = (
        args.subsets
    )  # TODO: Determine if splitting subset freq from whole callset agg
    include_non_release = (
        args.include_non_release
    )  # TODO: Discuss if non-release samples will be included in calculating AF for homalt fix

    hl.init(
        log=f"/generate_frequency_data{'.' + '_'.join(subsets) if subsets else ''}.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    vds = get_gnomad_v4_vds(
        split=True, release_only=not include_non_release, test=args.test
    )  # TODO: Do we need to filter to test here? Is this anymore performant?
    meta_ht = meta.ht()

    if args.test:
        logger.info("Filtering to DRD2 in test VDS for testing purposes...")
        test_interval = [
            hl.parse_locus_interval("chr11:113409605-113475691")
        ]  # NOTE: test on gene DRD2
        vds = hl.vds.filter_intervals(vds, test_interval)
        variants, samples = vds.variant_data.count()
        logger.info(
            "Test VDS has %s variants in DRD2 in %s samples...", variants, samples
        )

    # if include_non_release:
    #     all_sample_count = vds.variant_data.count()[1]
    #     hq_samples = sample_meta.filter(sample_meta.high_quality)
    #     vds = hl.vds.filter_samples(vds, hq_samples)
    #     filtered_sample_count = vds.variant_data.count()[1]
    #     logger.info(
    #         "VDS has been filtered from %s samples to %s high quality samples...",
    #         all_sample_count,
    #         filtered_sample_count,
    #     )

    if args.calculate_gatk_af_diff:
        logger.info(
            "Densifying VDS and determining frequency differences across GATK"
            " versions..."
        )
        mt = hl.vds.to_dense_mt(vds)
        mt = annotate_gatk_version(mt)
        mt = mt.annotate_cols(
            sex_karyotype=meta_ht[mt.col_key].sex_imputation.sex_karyotype,
            pop=meta_ht[mt.col_key].population_inference.pop,
        )

        logger.info("Computing adj and sex adjusted genotypes...")
        mt = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(  # TODO: Not needed just yet but we will
                mt.locus, mt.GT, mt.sex_karyotype
            ),
        )
        mt = annotate_adj(mt)
        logger.info("Annotating frequency...")
        mt = annotate_freq(
            mt,
            sex_expr=mt.sex_karyotype,
            pop_expr=mt.pop,
            additional_strata_expr={"gatk_fix": mt.gatk_fix},
        )

    logger.info("Checkpointing MatrixTable...")
    mt = mt.checkpoint("gs://gnomad-tmp-4day/freq/test_freq_aggs.mt", overwrite=True)

    # Validity checking
    mt.rows().show()
    mt.cols().show()
    mt.aggregate_cols(
        hl.agg.counter(mt.gatk_fix)
    )  # Qin: I don't see the counts, maybe we have to use .aggregate instead of aggregate_cols?


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
        "--include-non-release",
        help="Includes un-releasable samples in the frequency calculations.",
        action="store_true",
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
