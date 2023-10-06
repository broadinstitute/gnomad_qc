# noqa: D100

import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import split_info_annotation

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.variant_qc import get_variant_qc_result

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("import_variant_qc_vcf")
logger.setLevel(logging.INFO)


def import_variant_qc_vcf(
    vcf_path: str,
    model_id: str,
    num_partitions: int = 5000,
    overwrite: bool = False,
    import_header_path: Optional[str] = None,
    array_elements_required: bool = False,
) -> None:
    """
    Import variant QC result site VCF into a HT.

    :param vcf_path: Path to input variant QC result site vcf. This can be specified
        as Hadoop glob patterns.
    :param model_id: Model ID for the variant QC results. Must start with 'rf_',
        'vqsr_', or 'if_'.
    :param num_partitions: Number of partitions to use for the output HT.
    :param overwrite: Whether to overwrite data already present in the output HT.
    :param import_header_path: Optional path to a header file to use for import.
    :return: None
    """
    model_type = model_id.split("_")[0]
    if model_type not in ["rf", "vqsr", "if"]:
        raise ValueError(
            f"Model ID must start with 'rf_', 'vqsr_', or 'if_', but got {model_id}"
        )

    logger.info(f"Importing variant QC annotations for model: {model_id}...")
    logger.info(f"Array elements field as {array_elements_required}")
    mt = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        reference_genome="GRCh38",
        array_elements_required=array_elements_required,
        header_file=import_header_path,
    ).repartition(num_partitions)

    ht = mt.rows()

    unsplit_count = None
    if not args.is_split:
        if model_type == "vqsr":
            as_vqslod_expr = {"AS_VQSLOD": ht.info.AS_VQSLOD.map(lambda x: hl.float(x))}
        else:
            as_vqslod_expr = {}
        ht = ht.annotate(
            info=ht.info.annotate(
                **as_vqslod_expr,
                AS_QUALapprox=ht.info.AS_QUALapprox.split("\|")[1:].map(
                    lambda x: hl.int(x)
                ),
                AS_VarDP=ht.info.AS_VarDP.split("\|")[1:].map(lambda x: hl.int(x)),
                AS_SB_TABLE=ht.info.AS_SB_TABLE.split("\|").map(
                    lambda x: hl.or_missing(
                        x != "", x.split(",").map(lambda y: hl.int(y))
                    )
                ),
            )
        )

        ht = ht.checkpoint(
            get_variant_qc_result(model_id, split=False).path,
            overwrite=overwrite,
        )

        unsplit_count = ht.count()
        ht = hl.split_multi_hts(ht)

        ht = ht.annotate(
            info=ht.info.annotate(**split_info_annotation(ht.info, ht.a_index)),
        )

    ht = ht.checkpoint(
        get_variant_qc_result(model_id, split=True).path,
        overwrite=overwrite,
    )

    split_count = ht.count()
    if unsplit_count is None:
        logger.info(f"Found {split_count} split variants in the VCF")
    else:
        logger.info(
            f"Found {unsplit_count} unsplit and {split_count} split variants in the VCF"
        )


def main(args):  # noqa: D103
    hl.init(log="/load_variant_qc_vcf.log", default_reference="GRCh38")

    logger.info(f"passed array elements required as: {args.array_elements_required}")

    import_variant_qc_vcf(
        args.vcf_path,
        args.model_id,
        args.n_partitions,
        args.overwrite,
        args.header_path,
        args.array_elements_required,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Whether to overwrite data already present in the output Table.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )
    parser.add_argument(
        "--vcf-path",
        help="Path to variant QC result VCF. Can be specified as Hadoop glob patterns.",
        required=True,
    )
    parser.add_argument(
        "--model-id",
        help="Model ID for the variant QC result HT.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--compute-info-method",
        help=(
            "Compute info method used to generate the variant QC results. Options are"
            " 'AS', 'quasi' or 'set_long_AS_missing_info'."
        ),
        type=str,
        required=True,
        choices=["AS", "quasi", "set_long_AS_missing_info"],
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of desired partitions for output Table.",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--header-path",
        help=(
            "Optional path to a header file to use for importing the variant QC result"
            " VCF."
        ),
    )
    parser.add_argument(
        "--array-elements-required",
        action="store_true",
        help="Pass if you would like array elements required in import_vcf to be true.",
    )
    parser.add_argument(
        "--is-split", action="store_true", help="Whether the VCF is already split."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
