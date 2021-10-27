import argparse
import logging
from typing import Optional

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import split_info_annotation
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_vqsr_filters

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("import_vqsr")
logger.setLevel(logging.INFO)


def import_vqsr(
    vqsr_path: str,
    vqsr_type: str = "alleleSpecificTrans",
    num_partitions: int = 5000,
    overwrite: bool = False,
    import_header_path: Optional[str] = None,
) -> None:
    """
    Imports vqsr site vcf into a HT
    :param vqsr_path: Path to input vqsr site vcf. This can be specified as Hadoop glob patterns
    :param vqsr_type: One of `classic`, `alleleSpecific` (allele specific) or `alleleSpecificTrans`
        (allele specific with transmitted singletons)
    :param num_partitions: Number of partitions to use for the VQSR HT
    :param overwrite: Whether to overwrite imported VQSR HT
    :param import_header_path: Optional path to a header file to use for import
    :return: None
    """

    logger.info(f"Importing VQSR annotations for {vqsr_type} VQSR...")
    mt = hl.import_vcf(
        vqsr_path,
        force_bgz=True,
        reference_genome="GRCh38",
        header_file=import_header_path,
    ).repartition(num_partitions)

    ht = mt.rows()

    ht = ht.annotate(
        info=ht.info.annotate(
            AS_VQSLOD=ht.info.AS_VQSLOD.map(lambda x: hl.float(x)),
            AS_QUALapprox=ht.info.AS_QUALapprox.split("\|")[1:].map(
                lambda x: hl.int(x)
            ),
            AS_VarDP=ht.info.AS_VarDP.split("\|")[1:].map(lambda x: hl.int(x)),
            AS_SB_TABLE=ht.info.AS_SB_TABLE.split("\|").map(
                lambda x: x.split(",").map(lambda y: hl.int(y))
            ),
        )
    )

    ht = ht.checkpoint(
        get_vqsr_filters(f"vqsr_{vqsr_type}", split=False, finalized=False).path,
        overwrite=overwrite,
    )

    unsplit_count = ht.count()
    ht = hl.split_multi_hts(ht)

    ht = ht.annotate(
        info=ht.info.annotate(**split_info_annotation(ht.info, ht.a_index)),
    )

    ht = ht.checkpoint(
        get_vqsr_filters(f"vqsr_{vqsr_type}", split=True, finalized=False).path,
        overwrite=overwrite,
    )
    split_count = ht.count()
    logger.info(
        f"Found {unsplit_count} unsplit and {split_count} split variants with VQSR annotations"
    )


def main(args):
    hl.init(log="/load_data.log", default_reference="GRCh38")

    import_vqsr(
        args.vqsr_vcf_path,
        args.vqsr_type,
        args.n_partitions,
        args.overwrite,
        args.header_path,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--vqsr_vcf_path", help="Path to VQSR VCF. Can be specified as Hadoop glob patterns")
    parser.add_argument(
        "--vqsr_type",
        help="Type of VQSR corresponding to the VQSR VCF being loaded",
        default="alleleSpecificTrans",
        choices=["classic", "alleleSpecific", "alleleSpecificTrans"],
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output Table",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--header_path",
        help="Optional path to a header file to use for importing VQSR VCF",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
