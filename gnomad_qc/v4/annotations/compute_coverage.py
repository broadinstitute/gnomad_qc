# noqa: D100

import argparse
import logging
from typing import List, Union

import hail as hl
from gnomad.resources.grch38.gnomad import (
    CURRENT_EXOME_COVERAGE_RELEASE,
    coverage,
    coverage_tsv_path,
)
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.sparse_mt import compute_coverage_stats

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.basics import calling_intervals, get_gnomad_v4_vds

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("coverage")
logger.setLevel(logging.INFO)


def main(args):  # noqa: D103
    hl.init(
        log="/coverage.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    try:
        coverage_version = (
            args.coverage_version
            if args.coverage_version
            else CURRENT_EXOME_COVERAGE_RELEASE
        )
        test = args.test

        if args.compute_coverage_ht:
            check_resource_existence(
                output_step_resources={
                    "--compute-coverage-ht": [
                        coverage("exomes").versions[coverage_version].path
                    ],
                },
                overwrite=args.overwrite,
            )

            logger.info("Running compute coverage...")
            # Read in context Table.
            ref_ht = vep_context.versions["105"].ht()
            if test:
                ref_ht = ref_ht._filter_partitions(range(50))

            # Read in calling intervals.
            interval_ht = calling_intervals(
                interval_name=args.calling_interval_name,
                calling_interval_padding=args.calling_interval_padding,
            ).ht()

            # Read in VDS.
            vds = get_gnomad_v4_vds(
                release_only=True,
                test=test,
                filter_partitions=range(2) if test else None,
            )

            # Compute coverage stats.
            coverage_ht = compute_coverage_stats(vds, ref_ht, interval_ht)

            # Checkpoint Table
            coverage_ht = coverage_ht.checkpoint(
                "gs://gnomad-tmp/gnomad.genomes_v4.coverage.summary.ht",
                overwrite=args.overwrite,
            )

            # Write final result if not testing.
            if not test:
                coverage_ht = coverage_ht.naive_coalesce(5000)

                coverage_ht.write(
                    coverage("exomes").versions[coverage_version].path,
                    overwrite=args.overwrite,
                )

        if args.export_coverage_tsv:
            logger.info("Exporting coverage tsv...")
            ht = coverage("exomes").versions[coverage_version].ht()
            ht.export(coverage_tsv_path("exomes", coverage_version))

    finally:
        hl.copy_log(f"gs://gnomad-tmp-4day/coverage/compute_coverage.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to run a test using only the first 2 partitions of the VDS test"
            " dataset."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-coverage-ht", help="Compute the coverage HT.", action="store_true"
    )
    parser.add_argument(
        "--coverage-version",
        type=str,
        help=(
            "Specifies coverage version to read/write. If not set,"
            " gnomad.resources.grch38.gnomad.CURRENT_EXOME_COVERAGE_RELEASE is used."
        ),
    )
    parser.add_argument(
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb',"
            " 'broad', 'intersection', or 'union'."
        ),
        type=str,
        choices=["ukb", "broad", "intersection", "union"],
        default="union",
    )
    parser.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or"
            " 50 bp."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )
    parser.add_argument(
        "--export-coverage-tsv", help="Exports coverage TSV file.", action="store_true"
    )

    main(parser.parse_args())
