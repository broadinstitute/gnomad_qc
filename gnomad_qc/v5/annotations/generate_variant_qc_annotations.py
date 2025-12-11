"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import annotate_allele_info, get_lowqual_expr
from gnomad.utils.sparse_mt import split_info_annotation

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.annotations.annotation_utils import get_adj_expr
from gnomad_qc.v5.resources.annotations import (
    aou_annotated_sites_only_vcf,
    aou_vcf_header,
    get_info_ht,
)
from gnomad_qc.v5.resources.basics import get_aou_vds, get_logging_path
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_annotations")
logger.setLevel(logging.INFO)


def generate_ac_info_ht(vds: hl.vds.VariantDataset) -> hl.Table:
    """Compute AC and AC_raw annotations for each allele count filter group.

    :param vds: VariantDataset to use for computing AC and AC_raw annotations.
    :return: Table with AC and AC_raw annotations split by high quality, release, and unrelated.
    """
    mt = vds.variant_data

    ac_filter_groups = {
        "high_quality": mt.meta.high_quality,
        "release": mt.meta.release,
        "unrelated": ~mt.meta.relatedness_inference.relatedness_filters.related,
    }

    mt = mt.annotate_cols(_ac_filter_groups=ac_filter_groups)
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    ac_info_expr = hl.struct()
    # First compute ACs for each non-ref allele, grouped by adj.
    grp_ac_expr = {
        f: hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai) & mt._ac_filter_groups[f],
                hl.agg.group_by(
                    get_adj_expr(mt.LGT, mt.GQ, mt.LAD),
                    hl.agg.sum(
                        mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[
                            mt.LA.index(ai)
                        ]
                    ),
                ),
            ),
            mt.alt_alleles_range_array,
        )
        for f in ac_filter_groups
    }

    # Then, for each non-ref allele, compute
    # 'AC' as the adj group
    # 'AC_raw' as the sum of adj and non-adj groups
    ac_info_expr = ac_info_expr.annotate(
        **{
            f"AC{'_' + f if f else f}_raw": grp.map(
                lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))
            )
            for f, grp in grp_ac_expr.items()
        },
        **{
            f"AC{'_' + f if f else f}": grp.map(lambda i: hl.int32(i.get(True, 0)))
            for f, grp in grp_ac_expr.items()
        },
    )

    ac_info_ht = mt.select_rows(AC_info=ac_info_expr).rows()

    # Split multi-allelic sites.
    ac_info_ht = annotate_allele_info(ac_info_ht)
    ac_info_ht = ac_info_ht.annotate(
        **{
            a: ac_info_ht[a].annotate(
                **split_info_annotation(ac_info_ht[a], ac_info_ht.a_index),
            )
            for a in ["AC_info"]
        },
    )
    ac_info_ht = ac_info_ht.drop("allele_info")

    return ac_info_ht


def create_info_ht(
    vcf_path: str,
    header_path: str,
    lowqual_indel_phred_het_prior: int = 40,
    vds: hl.vds.VariantDataset = None,
    test: bool = False,
) -> hl.Table:
    """
    Import a VCF of AoU annotated sites, reformat annotations, and add AS_lowqual.

    :param vcf_path: Path to the annotated sites-only VCF.
    :param header_path: Path to the header file for the VCF.
    :param lowqual_indel_phred_het_prior: Phred-scaled prior for a het genotype at a site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to be more consistent with the filtering used by Broad's Data Sciences Platform for VQSR.
    :param vds: VariantDataset to use for computing AC and AC_raw annotations.
    :param test: Whether to write run a test using just the first two partitions of the loaded VCF.
    :return: Hail Table with reformatted annotations.
    """
    ht = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        header_file=header_path,
        reference_genome="GRCh38",
    ).rows()

    if test:
        ht = ht._filter_partitions(range(2))

    logger.info("Reformatting annotations...")
    array_annotations = [
        "AS_FS",
        "AS_MQ",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_SOR",
    ]

    # AS_VarDP is in the format of "ref|alt" and VarDP is an int of the alt value from this, so can just set AS_VarDP to VarDP.
    # AS_SB_TABLE is an alternate formatting (string ref1, ref2 | alt 1, alt2)
    # of SB_TABLE,  so can just set AS_SB_TABLE to SB_TABLE.
    info_updates = {
        # Convert single-element array annotations to float64.
        **{ann: hl.float64(ht.info[ann][0]) for ann in array_annotations},
        # Extract singular element from AS_QD array and convert to float32.
        "AS_QD": hl.float32(ht.info.AS_QD[0]),
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

    logger.info("Adding AC info annotations to info ht...")
    ac_info_ht = generate_ac_info_ht(vds)
    ht = ht.annotate(AC_info=ac_info_ht[ht.key].AC_info)
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

            aou_vds = get_aou_vds(
                test=test,
                high_quality_only=True,
                annotate_meta=True,
            )

            ht = create_info_ht(
                vcf_path=aou_annotated_sites_only_vcf,
                header_path=aou_vcf_header,
                lowqual_indel_phred_het_prior=args.lowqual_indel_phred_het_prior,
                vds=aou_vds,
                test=test,
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
