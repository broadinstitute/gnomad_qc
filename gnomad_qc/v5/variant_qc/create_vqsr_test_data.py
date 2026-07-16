"""Create a schema-faithful VQSR test input from gnomAD v4 data.

The v5 VQSR pipeline (``vqsr.py``) consumes a sites-only info VCF
(:func:`get_info_vcf_path`) and, optionally, a true-positive/singleton VCF
(:func:`get_true_positive_vcf_path`). Those are produced by
``generate_variant_qc_annotations.py`` (branch ``kl/info_edits``) from AoU data,
which is not yet available.

This script stands in for that step so the VQSR pipeline can be exercised
end-to-end before the real v5 info VCF exists. It reads the public gnomAD v4
exomes info Table, restricts it to a small test region, and remaps its
allele-specific annotations onto the exact INFO schema that the v5 info VCF
export produces, then writes the info VCF and a transmitted-singleton VCF to the
``test=True`` paths that ``vqsr.py --test`` reads.

The v4 ``AS_info`` struct already carries every field the v5 export keeps that
VQSR reads, so the mapping is a straight rename plus the two type/shape
adjustments the v5 export applies:

    - ``AS_QD`` is cast to ``float32`` (v5 stores it as ``float32``).
    - ``AS_SB_TABLE`` is reshaped from a flat 4-element array to the nested
      ``[[ref_fwd, ref_rev], [alt_fwd, alt_rev]]`` form expected by
      ``adjust_vcf_incompatible_types``.
"""

import argparse
import logging

import hail as hl
from gnomad.utils.vcf import adjust_vcf_incompatible_types

from gnomad_qc.v4.resources.annotations import get_info as get_v4_info
from gnomad_qc.v5.resources.annotations import (
    get_info_vcf_path,
    get_true_positive_vcf_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_vqsr_test_data")
logger.setLevel(logging.INFO)

# Fields the v5 info VCF export carries that VQSR (and the downstream VCF import)
# read. Every one is present in the v4 exomes 'AS_info' struct.
FLOAT64_FEATURES = [
    "AS_MQRankSum",
    "AS_ReadPosRankSum",
    "AS_FS",
    "AS_SOR",
    "AS_MQ",
]


def build_info_vcf(region: str, out_path: str) -> None:
    """
    Build a sites-only info VCF matching the v5 info VCF INFO schema.

    :param region: Locus interval string (e.g. "chr22") to restrict the v4 info
        Table to.
    :param out_path: Path to write the info VCF (bgzipped, tabix indexed).
    """
    logger.info("Reading v4 exomes split info HT...")
    ht = hl.read_table(get_v4_info(split=True).path)

    logger.info("Restricting to test region %s...", region)
    ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval(region, reference_genome="GRCh38")]
    )

    # Remap the v4 'AS_info' struct onto the v5 info VCF INFO schema.
    ht = ht.select(
        info=hl.struct(
            AS_QD=hl.float32(ht.AS_info.AS_QD),
            AS_QUALapprox=hl.int64(ht.AS_info.AS_QUALapprox),
            AS_VarDP=hl.int32(ht.AS_info.AS_VarDP),
            AS_SB_TABLE=ht.AS_info.AS_SB_TABLE,
            **{f: ht.AS_info[f] for f in FLOAT64_FEATURES},
        )
    ).select_globals()

    # Reshape AS_SB_TABLE into the nested array form the v5 export uses so
    # adjust_vcf_incompatible_types emits it correctly.
    ht = ht.annotate(
        info=ht.info.annotate(
            AS_SB_TABLE=hl.array([ht.info.AS_SB_TABLE[0:2], ht.info.AS_SB_TABLE[2:4]])
        )
    )

    ht = adjust_vcf_incompatible_types(ht, pipe_delimited_annotations=[])

    logger.info("Exporting info VCF to %s", out_path)
    hl.export_vcf(ht, out_path, tabix=True)


def build_true_positive_vcf(region: str) -> None:
    """
    Build sites-only transmitted-singleton VCFs (raw and adj) for VQSR training.

    Uses v4 singletons (AC == 1) in the test region as a stand-in for the v5
    transmitted-singleton true-positive set. VQSR only needs the sites, so the
    exported VCFs are sites-only.

    :param region: Locus interval string (e.g. "chr22") to restrict to.
    """
    ht = hl.read_table(get_v4_info(split=True).path)
    ht = hl.filter_intervals(
        ht, [hl.parse_locus_interval(region, reference_genome="GRCh38")]
    )
    tp_ht = ht.filter(ht.AC_info.AC == 1).select().select_globals()
    tp_ht = tp_ht.checkpoint(hl.utils.new_temp_file("vqsr_test_tp", "ht"))
    logger.info("True-positive (singleton) sites: %d", tp_ht.count())

    for adj in (False, True):
        out_path = get_true_positive_vcf_path(
            test=True, true_positive_type="transmitted_singleton", adj=adj
        )
        logger.info("Exporting true-positive VCF (adj=%s) to %s", adj, out_path)
        hl.export_vcf(tp_ht, out_path, tabix=True)


def main(args):
    """Create the VQSR test info VCF and true-positive VCF from v4 data."""
    hl.init(
        backend="batch",
        tmp_dir="gs://gnomad-tmp-4day",
        gcs_requester_pays_configuration=args.gcp_billing_project,
        default_reference="GRCh38",
        regions=["us-central1"],
    )

    build_info_vcf(args.test_region, get_info_vcf_path(test=True))

    if not args.skip_true_positive_vcf:
        build_true_positive_vcf(args.test_region)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--test-region",
        type=str,
        default="chr22",
        help="Locus interval to restrict the v4 info Table to. Default is chr22.",
    )
    parser.add_argument(
        "--skip-true-positive-vcf",
        action="store_true",
        help="Skip building the transmitted-singleton true-positive VCF.",
    )
    parser.add_argument(
        "--gcp-billing-project",
        type=str,
        default="broad-mpg-gnomad",
        help="Google Cloud billing project for reading requester-pays buckets.",
    )
    return parser


if __name__ == "__main__":
    main(get_script_argument_parser().parse_args())
