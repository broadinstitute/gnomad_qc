import argparse
import logging

import hail as hl

from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad.utils.annotations import annotate_adj
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_logging_path,
    get_gnomad_v4_vds,
)
from gnomad_qc.v4.resources.sample_qc import (
    gnomad_v4_pre_ld_prune_qc_sites,
    gnomad_v4_qc_sites,
    qc,
)
from gnomad_qc.slack_creds import slack_token


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_qc_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(
        log="/generate_qc_mt.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default
    hl._set_flags(use_new_shuffle="1")

    overwrite = args.overwrite
    ld_r2 = args.ld_r2
    test = args.test

    try:
        if args.create_intermediate_qc_mt:
            vds = get_gnomad_v4_vds(split=True, remove_hard_filtered_samples=False)
            if test:
                logger.info("Filtering to the first five partitions")
                vds = hl.vds.VariantDataset(
                    vds.reference_data._filter_partitions(range(20)),
                    vds.variant_data._filter_partitions(range(20)),
                )

            logger.info("Filtering variants to predetermined QC variants...")
            vds = hl.vds.filter_variants(vds, gnomad_v4_pre_ld_prune_qc_sites.ht())

            logger.info("Densifying data...")
            mt = hl.vds.to_dense_mt(vds)

            logger.info("Writing temp MT...")
            hl._set_flags(no_whole_stage_codegen="1")
            mt = mt.select_entries("GT", "GQ", "DP", "AD")

            mt_cp_path = get_checkpoint_path(
                f"dense_pre_ld_prune_qc_sites{'.test' if test else ''}", mt=True
            )
            mt.write(mt_cp_path, overwrite=True)
            mt = hl.read_matrix_table(mt_cp_path, _n_partitions=args.n_partitions)
            logger.info(
                "Number of predetermined QC variants found in the VDS: %d...",
                mt.count_rows(),
            )

            logger.info(
                "Annotating with adj and filtering to pre LD-pruned QC sites..."
            )
            mt = annotate_adj(mt)
            mt = get_qc_mt(
                mt,
                bi_allelic_only=args.bi_allelic_only,
                min_af=args.min_af,
                min_callrate=args.min_callrate,
                min_inbreeding_coeff_threshold=args.min_inbreeding_coeff_threshold,
                min_hardy_weinberg_threshold=None,
                apply_hard_filters=False,
                ld_r2=None,  # Done below to avoid an error
                filter_lcr=False,  # Already filtered from the initial set of QC variants
                filter_decoy=False,  # Doesn't exist for hg38
                filter_segdup=False,  # Already filtered from the initial set of QC variants
            )
            mt = mt.checkpoint(
                get_checkpoint_path(
                    f"dense_pre_ld_prune_qc_mt{'.test' if test else ''}", mt=True
                ),
                overwrite=True,
            )
            logger.info("Number of pre LD-pruned QC sites: %d...", mt.count_rows())

        if args.ld_prune:
            logger.warning(
                "The LD-prune step of this function requires non-preemptible workers only!"
            )
            mt = hl.read_matrix_table(
                get_checkpoint_path(
                    f"dense_pre_ld_prune_qc_mt{'.test' if test else ''}", mt=True
                )
            )
            logger.info("Number of pre LD-pruned QC sites: %d...", mt.count_rows())

            logger.info("Performing LD prune...")
            hl._set_flags(no_whole_stage_codegen=None)
            unfiltered_qc_mt = mt.unfilter_entries()
            ht = hl.ld_prune(unfiltered_qc_mt.GT, r2=ld_r2)
            ht = ht.annotate_globals(ld_r2=ld_r2)
            ht = ht.checkpoint(
                get_checkpoint_path("qc_sites.test")
                if test
                else gnomad_v4_qc_sites.path,
                overwrite=overwrite,
            )
            logger.info("Final number of QC sites: %d...", ht.count())

            logger.info("Filtering dense MT to final QC sites...")
            mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
            mt = mt.annotate_globals(qc_mt_params=mt.qc_mt_params.annotate(ld_r2=ld_r2))
            mt.write(
                get_checkpoint_path("dense_ld_prune_qc_mt.test", mt=True)
                if test
                else qc.path,
                overwrite=overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("generate_qc_mt"))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the VDS.", action="store_true"
    )
    parser.add_argument(
        "--create-intermediate-qc-mt",
        help="Create a checkpointed intermediate dense QC MT with all filters except the LD prune step.",
        action="store_true",
    )
    parser.add_argument(
        "--bi-allelic-only",
        help="Filter to variants that are bi-allelic.",
        action="store_true",
    )
    parser.add_argument(
        "--min-af",
        help="Minimum variant allele frequency to retain variant in QC MatrixTable.",
        default=0.001,
        type=float,
    )
    parser.add_argument(
        "--min-callrate",
        help="Minimum variant callrate to retain variant in qc matrix table.",
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--min-inbreeding-coeff-threshold",
        help="Minimum site inbreeding coefficient to keep.",
        default=-0.8,
        type=float,
    )
    parser.add_argument(
        "--ld-prune",
        help="Perform ld prune of the checkpointed intermediate dense QC MT.",
        action="store_true",
    )
    parser.add_argument(
        "--ld-r2",
        help="LD pruning cutoff",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--n-partitions",
        help="Desired number of partitions for output QC MatrixTable",
        default=1000,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
