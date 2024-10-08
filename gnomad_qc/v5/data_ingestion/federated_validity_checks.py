"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl
from gnomad.assessment.validity_checks import summarize_variants
from gnomad.resources.grch38.gnomad import public_release

from gnomad_qc.v5.resources.basics import get_logging_path

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("federated_validity_checks")
logger.setLevel(logging.INFO)


def main(args):
    """Perform validity checks for federated data."""
    hl.init(
        log="/federated_validity_checks.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")
    test_n_partitions = args.test_n_partitions

    try:
        # TODO: Add resources to intake federated data once obtained.
        ht = public_release(data_type="exomes").ht()

        if test_n_partitions:
            ht = ht._filter_partitions(range(test_n_partitions))

        expected_contigs = [
            i
            for i in hl.get_reference("GRCh38").contigs
            if i in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        ]
        logger.info("Summarizing variants and checking contigs.")
        summarize_variants(ht, expected_contigs=expected_contigs)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("federated_validity_checks"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the input for testing purposes. Defaults"
            "to 2 if passed without a value."
        ),
        nargs="?",
        const=2,
        type=int,
    )

    args = parser.parse_args()
    main(args)
