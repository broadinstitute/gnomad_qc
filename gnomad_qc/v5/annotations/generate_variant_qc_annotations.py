"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import bi_allelic_expr
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.variant_qc.pipeline import generate_sib_stats, generate_trio_stats

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj
from gnomad_qc.v5.resources.annotations import get_sib_stats, get_trio_stats
from gnomad_qc.v5.resources.basics import (
    WORKSPACE_BUCKET,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import GNOMAD_TMP_BUCKET
from gnomad_qc.v5.resources.sample_qc import (
    dense_trios,
    finalized_outlier_filtering,
    pedigree,
    relatedness,
)
from gnomad_qc.v5.sample_qc.identify_trios import filter_relatedness_ht

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("generate_variant_qc_annotations")
logger.setLevel(logging.INFO)


def run_generate_trio_stats(
    mt: hl.MatrixTable,
    fam_ped: hl.Pedigree,
) -> hl.Table:
    """
    Generate trio transmission stats from a VariantDataset and pedigree info.

    :param mt: Dense trio MatrixTable.
    :param fam_ped: Pedigree containing trio info.
    :return: Table containing trio stats.
    """
    # Filter to autosomes and bi-allelic sites, also annotate adj.
    mt = filter_to_autosomes(mt)
    mt = mt.filter_rows(bi_allelic_expr(mt))
    mt = hl.experimental.sparse_split_multi(mt)
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
    if args.rwb:
        environment = "rwb"
        hl.init(
            log="/home/jupyter/workspaces/gnomadproduction/identify_trios.log",
            tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
        )
    else:
        environment = "batch"
        hl.init(
            tmp_dir=f"gs://{GNOMAD_TMP_BUCKET}/tmp/4_day",
            log="generate_variant_qc_annotations.log",
        )
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test_n_partitions = args.test_n_partitions
    test = test_n_partitions is not None

    # NOTE: VDS will have 'aou_' prefix on sample IDs.
    vds = get_aou_vds(
        high_quality_only=True,
        filter_partitions=range(test_n_partitions) if test_n_partitions else None,
        annotate_meta=True,
    )
    mt = vds.variant_data

    try:
        if args.generate_trio_stats:
            logger.info("Generating trio stats...")
            trio_stats_ht_path = get_trio_stats(test=test).path
            check_resource_existence(
                output_step_resources={"trio_stats_ht": trio_stats_ht_path}
            )

            ht = run_generate_trio_stats(
                dense_trios(test=test).mt(), pedigree(test=test).pedigree()
            )
            ht.write(trio_stats_ht_path, overwrite=overwrite)

        if args.generate_sibling_stats:
            logger.info("Generating sibling stats...")
            sib_stats_ht_path = get_sib_stats(test=test).path
            check_resource_existence(
                output_step_resources={"sib_stats_ht": sib_stats_ht_path}
            )
            ht = run_generate_sib_stats(
                mt, relatedness().ht(), finalized_outlier_filtering().ht()
            )
            ht.write(sib_stats_ht_path, overwrite=overwrite)

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            get_logging_path("generate_variant_qc_annotations", environment=environment)
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Create script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rwb",
        help="Run the script in RWB environment.",
        action="store_true",
    )
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
