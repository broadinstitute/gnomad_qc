"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl
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
from gnomad_qc.v5.resources.sample_qc import dense_trios, pedigree, relatedness

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
    # Add adj annotation and convert LGT to GT since only
    # autosomal bi-allelics are used to calculate trio stats.
    mt = annotate_adj(mt)
    mt = mt.transmute_entries(GT=mt.LGT)
    mt = hl.trio_matrix(mt, pedigree=fam_ped, complete_trios=True)
    return generate_trio_stats(mt)


def run_generate_sib_stats(
    mt: hl.MatrixTable,
    relatedness_ht: hl.Table,
) -> hl.Table:
    """
    Generate sibling stats from a VariantDataset and relatedness info.

    :param mt: Input MatrixTable.
    :param relatedness_ht: Table containing relatedness info.
    :return: Table containing sibling stats.
    """
    mt = annotate_adj(mt)
    mt = hl.experimental.sparse_split_multi(mt)
    return generate_sib_stats(mt, relatedness_ht)


def main(args):
    """Generate all variant annotations needed for variant QC."""
    if args.rwb:
        environment = "rwb"
        hl.init(
            log="/home/jupyter/workspaces/gnomadproduction/generate_variant_qc_annotations.log",
            tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
        )
    else:
        environment = "batch"
        hl.init(
            tmp_dir=f"gs://{GNOMAD_TMP_BUCKET}-4day",
            log="generate_variant_qc_annotations.log",
        )
        # TODO: Add machine configurations for Batch.
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

    try:
        if args.generate_trio_stats:
            logger.info("Generating trio stats...")
            trio_stats_ht_path = get_trio_stats(test=test, environment=environment).path
            check_resource_existence(
                output_step_resources={"trio_stats_ht": [trio_stats_ht_path]},
                overwrite=overwrite,
            )

            ht = run_generate_trio_stats(
                dense_trios(test=test).mt(), pedigree(test=test).pedigree()
            )
            ht.write(trio_stats_ht_path, overwrite=overwrite)

        if args.generate_sibling_stats:
            logger.info("Generating sibling stats...")
            sib_stats_ht_path = get_sib_stats(test=test, environment=environment).path
            check_resource_existence(
                output_step_resources={"sib_stats_ht": sib_stats_ht_path},
                overwrite=overwrite,
            )
            # Note: Checked sibling IDs; none of them have sample ID collisions.
            ht = run_generate_sib_stats(vds.variant_data, relatedness().ht())
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
