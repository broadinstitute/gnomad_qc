"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl

from gnomad.assessment.validity_checks import summarize_variants
from gnomad.resources.grch38.gnomad import public_release


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
    test = args.test

    # TODO: Add resources to intake federated data once obtained.
    if test:
        ht = public_release(data_type="exomes").ht()
        ht = ht._filter_partitions(range(2))

    expected_contigs = [
        i
        for i in hl.get_reference("GRCh38").contigs
        if i in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    ]
    logger.info("Summarizing variants and checking contigs.")
    summarize_variants(ht, expected_contigs=expected_contigs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test",
        help="Use the first two parititons of the gnomAD v4 exomes public release for testing.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
