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
        check_resource_existence(
            output_step_resources={
                "info_ht": [get_info_ht(test=test).path],
            },
            overwrite=overwrite,
        )

        logger.info("Importing VCF...")
        # vcf_path = "gs://fc-secure-b25d1307-7763-48b8-8045-fcae9caadfa1/echo_full_gnomad_annotated.sites-only.vcf.gz"
        # header_path = "gs://fc-secure-b25d1307-7763-48b8-8045-fcae9caadfa1/tmp/4_day/aou_annotation_sites_only_header.vcf"
        ht = hl.import_vcf(
            aou_annotated_sites_only_vcf,
            force_bgz=True,
            header_file=aou_vcf_header,
            reference_genome="GRCh38",
        ).rows()

        logger.info("Reformatting annotations...")
        array_annotations = ["AS_FS", "AS_MQ", "AS_MQRankSum", "AS_ReadPosRankSum"]

        # Build a dictionary of info updates.
        info_updates = {
            # Convert single-element array annotations to float64.
            **{ann: hl.float64(ht.info[ann][0]) for ann in array_annotations},
            # Extract singular element from  AS_QD array and convert to int32.
            "AS_QD": hl.int32(ht.info.AS_QD[0]),
            # Convert AS_QUALapprox to int64.
            "AS_QUALapprox": hl.int64(ht.info.QUALapprox),
            # Extract ALT value from AS_VarDP string and convert to int32.
            "AS_VarDP": hl.int32(ht.info.AS_VarDP.split("|")[1]),
        }

        # Apply info updates.
        ht = ht.transmute(info=ht.info.annotate(**info_updates))

        ht = ht.annotate(
            info=ht.info.annotate(
                AS_SB_TABLE=ht.info.AS_SB_TABLE.split("\||,").map(hl.parse_int)
            )
        )

        ht = ht.annotate(info=ht.info.drop("SB", "VarDP"))

        # Add AS_lowqual annotation.
        ht = ht.annotate(
            AS_lowqual=get_lowqual_expr(
                ht.alleles,
                ht.info.AS_QUALapprox,
                indel_phred_het_prior=args.lowqual_indel_phred_het_prior,
            )
        )

        ht.write(get_info_ht(test=test).path, overwrite=overwrite)
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
