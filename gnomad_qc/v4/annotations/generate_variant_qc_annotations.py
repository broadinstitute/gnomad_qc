"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import (
    INFO_INT32_SUM_AGG_FIELDS,
    INFO_SUM_AGG_FIELDS,
    default_compute_info,
    get_as_info_expr,
    get_site_info_expr,
    split_info_annotation,
    split_lowqual_annotation,
)
from gnomad.utils.vcf import adjust_vcf_incompatible_types

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_info, info_vcf_path
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def split_info(info_ht: hl.Table) -> hl.Table:
    """
    Generate an info Table with split multi-allelic sites from the multi-allelic info Table.

    :param info_ht: Info Table with unsplit multi-allelics.
    :return: Info Table with split multi-allelics.
    """
    info_ht = hl.split_multi(info_ht)

    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation(info_ht.info, info_ht.a_index),
        ),
        AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )
    return info_ht


def main(args):
    """Generate all variant annotations needed for variant QC."""
    hl.init(
        default_reference="GRCh38",
        log="/variant_qc_annotations.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test_dataset = args.test_dataset
    test_n_partitions = args.test_n_partitions
    test = test_dataset or test_n_partitions
    overwrite = args.overwrite
    mt = get_gnomad_v4_vds(test=test_dataset, high_quality_only=True).variant_data

    if test_n_partitions:
        mt = mt._filter_partitions(range(test_n_partitions))

    if args.compute_info:
        # TODO: is there any reason to also compute info per platform?
        default_compute_info(mt, site_annotations=True).write(
            get_info(split=False, test=test).path, overwrite=overwrite
        )

    if args.split_info:
        split_info(get_info(split=False, test=test).ht()).write(
            get_info(test=test).path, overwrite=overwrite
        )

    if args.export_info_vcf:
        hl.export_vcf(
            adjust_vcf_incompatible_types(get_info(split=False, test=test).ht()),
            info_vcf_path(test=test),
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--test-dataset", help="Use the test dataset as input.", action="store_true"
    )
    parser.add_argument(
        "--test-n-partitions",
        help="Use only 2 partitions of the VDS as input for testing purposes.",
        nargs="?",
        const=2,
        type=int,
    )
    parser.add_argument("--compute-info", help="Computes info HT.", action="store_true")
    parser.add_argument("--split-info", help="Splits info HT.", action="store_true")
    parser.add_argument(
        "--export-info-vcf", help="Export info as VCF.", action="store_true"
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
