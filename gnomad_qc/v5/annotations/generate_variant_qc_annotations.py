"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl
from gnomad.sample_qc.relatedness import filter_to_trios
from gnomad.variant_qc.pipeline import generate_sib_stats, generate_trio_stats

from gnomad_qc.v5.annotations.annotation_utils import annotate_adj
from gnomad_qc.v5.resources.basics import get_aou_vds, get_logging_path
from gnomad_qc.v5.resources.constants import GNOMAD_TMP_BUCKET
from gnomad_qc.v5.resources.sample_qc import pedigree, relatedness
from gnomad_qc.v5.sample_qc.identify_trios import filter_relatedness_ht

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("generate_variant_qc_annotations")
logger.setLevel(logging.INFO)


def run_generate_trio_stats(
    vds: hl.vds.VariantDataset,
    fam_ped: hl.Pedigree,
    fam_ht: hl.Table,
) -> hl.Table:
    """
    Generate trio transmission stats from a VariantDataset and pedigree info.

    :param vds: VariantDataset to generate trio stats from.
    :param fam_ped: Pedigree containing trio info.
    :param fam_ht: Table containing trio info.
    :return: Table containing trio stats.
    """
    # Filter the VDS to autosomes and trios.
    vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)
    vds = filter_to_trios(vds, fam_ht)

    vmt = vds.variant_data
    rmt = vds.reference_data

    # Filter the variant data to bi-allelic sites.
    vmt = vmt.filter_rows(hl.len(vmt.alleles) == 2)

    # Densify and annotate adj.
    vds = hl.vds.VariantDataset(reference_data=rmt, variant_data=vmt)
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.transmute_entries(GT=mt.LGT)
    mt = annotate_adj(mt)

    # Create trio matrix and generate trio stats.
    mt = hl.trio_matrix(mt, pedigree=fam_ped, complete_trios=True)
    return generate_trio_stats(mt)


def run_generate_sib_stats(
    mt: hl.MatrixTable,
    rel_ht: hl.Table,
    filter_ht: hl.Table,
) -> hl.Table:
    """
    Generate stats for the number of alternate alleles in common between sibling pairs.

    :param mt: MatrixTable to generate sibling stats from.
    :param rel_ht: Table containing relatedness info for pairs in `mt`.
    :param filter_ht: Table containing outlier filtering info for samples in `mt`.
    :return: Table containing sibling stats.
    """
    # Filter relatedness Table to non-filtered AoU genomes.
    rel_ht = filter_relatedness_ht(rel_ht, filter_ht)
    return generate_sib_stats(mt.transmute_entries(GT=mt.LGT), rel_ht)


def main(args):
    """Generate all variant annotations needed for variant QC."""
    hl.init(
        tmp_dir=f"gs://{GNOMAD_TMP_BUCKET}/tmp/4_day",
        log="generate_variant_qc_annotations.log",
    )
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test = args.test_n_partitions

    # NOTE: VDS will have 'aou_' prefix on sample IDs.
    vds = get_aou_vds(
        high_quality_only=True,
        filter_partitions=range(2) if test else None,
        annotate_meta=True,
    )

    try:
        if args.generate_trio_stats:
            pass
        if args.generate_sibling_stats:
            pass

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("generate_variant_qc_annotations.log"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Create script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--test-n-partitions",
        help="Use only n partitions of the VDS as input for testing purposes (default: 2).",
        nargs="?",
        const=2,
        type=int,
    )

    parser.add_argument(
        "--generate-trio-stats", help="Calculates trio stats", action="store_true"
    )
    parser.add_argument(
        "--generate-sibling-stats", help="Calculates sibling stats", action="store_true"
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
