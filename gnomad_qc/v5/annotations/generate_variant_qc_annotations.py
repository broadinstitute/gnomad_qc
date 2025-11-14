"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl

from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.constants import GNOMAD_TMP_BUCKET

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("generate_variant_qc_annotations")
logger.setLevel(logging.INFO)


def main(args):
    """Generate all variant annotations needed for variant QC."""
    hl.init(
        tmp_dir=f"gs://{GNOMAD_TMP_BUCKET}/tmp/4_day",
        log="generate_variant_qc_annotations.log",
    )
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test = args.test_n_partitions

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
