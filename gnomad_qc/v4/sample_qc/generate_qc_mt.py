import argparse
import logging

import hail as hl

from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad.utils.annotations import get_adj_expr, get_lowqual_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.file_utils import file_exists
from gnomad.utils.sparse_mt import densify_sites, get_as_info_expr, get_site_info_expr
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    get_ukbb_data,
    last_END_positions_ht_path,
    logging_path,
)
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE
from ukbb_qc.resources.sample_qc import qc_ht_path, qc_mt_path, qc_sites_path
from ukbb_qc.slack_creds import slack_token
from ukbb_qc.utils.utils import get_qc_mt_sites


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_qc_data")
logger.setLevel(logging.INFO)


def main(args):
    hl.init(log="/create_qc_data.log", default_reference="GRCh38")

    data_source = "broad"
    freeze = args.freeze

    try:
        if args.compute_qc_mt:
            logger.info("Reading in raw MT...")
            mt = get_ukbb_data(
                data_source,
                freeze,
                key_by_locus_and_alleles=True,
                split=False,
                raw=True,
                repartition=args.repartition,
                n_partitions=args.raw_partitions,
            )
            mt = mt.transmute_entries(**mt.gvcf_info)
            mt = mt.select_entries(
                "LGT",
                "GQ",
                "DP",
                "LAD",
                "LA",
                "END",
                "QUALapprox",
                "VarDP",
                "ReadPosRankSum",
                "MQRankSum",
                "SB",
                "RAW_MQandDP",
            )

            logger.info("Reading in QC MT sites from tranche 2/freeze 5...")
            if not file_exists(qc_sites_path()):
                get_qc_mt_sites()
            qc_sites_ht = hl.read_table(qc_sites_path())
            logger.info(f"Number of QC sites: {qc_sites_ht.count()}")

            logger.info("Densifying sites...")
            last_END_ht = hl.read_table(last_END_positions_ht_path(freeze))
            mt = densify_sites(mt, qc_sites_ht, last_END_ht)

            logger.info("Checkpointing densified MT")
            mt = mt.checkpoint(
                get_checkpoint_path(
                    data_source, freeze, name="dense_qc_mt_v2_sites", mt=True
                ),
                overwrite=True,
            )

            logger.info("Repartitioning densified MT")
            mt = mt.naive_coalesce(args.n_partitions)
            mt = mt.checkpoint(
                get_checkpoint_path(
                    data_source,
                    freeze,
                    name="dense_qc_mt_v2_sites.repartitioned",
                    mt=True,
                ),
                overwrite=True,
            )

            # NOTE: Need MQ, QD, FS for hard filters
            logger.info("Adding info and low QUAL annotations and filtering to adj...")
            info_expr = get_site_info_expr(mt)
            info_expr = info_expr.annotate(**get_as_info_expr(mt))
            mt = mt.annotate_rows(info=info_expr)
            mt = mt.select_entries(
                GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA),
                adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
            )

            logger.info("Checkpointing MT...")
            mt = mt.checkpoint(
                get_checkpoint_path(
                    data_source,
                    freeze,
                    name=f"{data_source}.freeze_{freeze}.qc_sites.mt",
                    mt=True,
                ),
                overwrite=True,
            )

            # NOTE: adding decoy and segdup hg38 resources/filters are still pending
            qc_mt = get_qc_mt(
                mt,
                min_af=args.min_af,
                min_callrate=args.min_callrate,
                apply_hard_filters=True,
                ld_r2=None,
                filter_lcr=False,
                filter_decoy=False,
                filter_segdup=False,
            )
            # NOTE: not using default indel_phred_het_prior to make sure lowqual values match with standard pre-VQSR lowqual values
            # NOTE: Using AS_QUALapprox[1] because these are biallelic sites
            mt = mt.annotate_rows(
                AS_lowqual=get_lowqual_expr(
                    mt.alleles, mt.info.AS_QUALapprox[1], indel_phred_het_prior=40,
                )
            )
            qc_mt = qc_mt.checkpoint(
                qc_mt_path(data_source, freeze), overwrite=args.overwrite,
            )
            logger.info(
                f"Total number of bi-allelic, high-callrate, common SNPs for sample QC: {qc_mt.count_rows()}"
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(data_source, freeze))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--freeze", help="Data freeze to use", default=CURRENT_FREEZE, type=int,
    )
    parser.add_argument(
        "--compute_qc_mt",
        help="Compute matrix to be used in sample qc",
        action="store_true",
    )
    parser.add_argument(
        "--min_af",
        help="Minimum variant allele frequency to retain variant in QC MatrixTable.",
        default=0.001,
        type=float,
    )
    parser.add_argument(
        "--min_callrate",
        help="Minimum variant callrate to retain variant in qc matrix table.",
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--compute_sample_qc_ht",
        help="Compute sample qc on qc matrix table",
        action="store_true",
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output QC MatrixTable",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
