"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import get_lowqual_expr

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.resources.annotations import (
    aou_annotated_sites_only_vcf,
    aou_vcf_header,
    get_info_ht,
)
from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_annotations")
logger.setLevel(logging.INFO)


def create_info_ht(
    vcf_path: str,
    header_path: str,
    lowqual_indel_phred_het_prior: int = 40,
) -> hl.Table:
    """
    Import a VCF of AoU annotated sites, reformat annotations, and add AS_lowqual.

    :param vcf_path: Path to the annotated sites-only VCF.
    :param header_path: Path to the header file for the VCF.
    :param lowqual_indel_phred_het_prior: "Phred-scaled prior for a het genotype at a site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to be more consistent with the filtering used by Broad's Data Sciences Platform for VQSR
    :return: Hail Table with reformatted annotations.
    """
    ht = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        header_file=header_path,
        reference_genome="GRCh38",
    ).rows()

    logger.info("Reformatting annotations...")

    array_annotations = ["AS_FS", "AS_MQ", "AS_MQRankSum", "AS_ReadPosRankSum"]

    # AS_VarDP is in the format of "ref|alt" and VarDP is an int of the alt value from this, so can just set AS_VarDP to VarDP.
    # AS_SB_TABLE is an alternate formatting (string ref1, ref2 | alt 1, alt2)
    # of SB_TABLE,  so can just set AS_SB_TABLE to SB_TABLE.
    info_updates = {
        # Convert single-element array annotations to float64.
        **{ann: hl.float64(ht.info[ann][0]) for ann in array_annotations},
        # Extract singular element from AS_QD array and convert to int32.
        "AS_QD": hl.int32(ht.info.AS_QD[0]),
        "AS_QUALapprox": hl.int64(ht.info.QUALapprox),
        "AS_VarDP": ht.info.VarDP,
        "AS_SB_TABLE": ht.info.SB_TABLE,
    }

    ht = ht.transmute(info=ht.info.annotate(**info_updates))
    ht = ht.annotate(
        info=ht.info.drop("SB", "QUALapprox", "VarDP", "SB_TABLE", "AS_RAW_MQ")
    )
    ht = ht.drop("rsid")

    ht = ht.annotate(
        AS_lowqual=get_lowqual_expr(
            ht.alleles,
            ht.info.AS_QUALapprox,
            indel_phred_het_prior=lowqual_indel_phred_het_prior,
        )
    )

    return ht


def main(args):
    """Generate all variant annotations needed for variant QC."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/variant_qc_annotations.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )

    hl.default_reference("GRCh38")

    test = args.test
    overwrite = args.overwrite

    try:
        if args.create_info_ht:
            info_ht_path = get_info_ht(test=test).path
            check_resource_existence(
                output_step_resources={
                    "info_ht": [info_ht_path],
                },
                overwrite=overwrite,
            )
            ht = create_info_ht(
                vcf_path=aou_annotated_sites_only_vcf,
                header_path=aou_vcf_header,
                lowqual_indel_phred_het_prior=args.lowqual_indel_phred_het_prior,
            )
            ht.write(info_ht_path, overwrite=overwrite)
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(
            get_logging_path("generate_variant_qc_annotations", environment="rwb")
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Write to test path.",
        action="store_true",
    )
    parser.add_argument(
        "--create-info-ht",
        help="Create the info ht containing annotations needed for variant QC.",
        action="store_true",
    )
    parser.add_argument(
        "--lowqual-indel-phred-het-prior",
        help="Phred-scaled prior for a het genotype at a site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to be more consistent with the filtering used by Broad's Data Sciences Platform for VQSR.",
        default=40,
        type=int,
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
