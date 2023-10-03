"""Script to compute coverage statistics on gnomAD v4 exomes."""
import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.sparse_mt import compute_coverage_stats

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.resources.basics import calling_intervals, get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import release_coverage, release_coverage_tsv_path

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("coverage")
logger.setLevel(logging.INFO)


def get_coverage_resources(
    test: bool,
    overwrite: bool,
    calling_interval_name: Optional[str] = None,
    calling_interval_padding: Optional[int] = None,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the coverage pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :param calling_interval_name: Name of calling intervals to use.
    :param calling_interval_padding: Padding to use for calling intervals.
    :return: PipelineResourceCollection containing resources for all steps of the
        coverage pipeline.
    """
    # Initialize coverage pipeline resource collection.
    coverage_pipeline = PipelineResourceCollection(
        pipeline_name="coverage",
        overwrite=overwrite,
    )
    # Create resource collection for each step of the coverage pipeline.
    if calling_interval_name and calling_interval_padding:
        coverage_input_resources = {
            "interval list": {
                "interval_ht": calling_intervals(
                    interval_name=calling_interval_name,
                    calling_interval_padding=calling_interval_padding,
                )
            },
        }
    else:
        coverage_input_resources = {}

    compute_coverage_ht = PipelineStepResourceCollection(
        "--compute-coverage-ht",
        input_resources=coverage_input_resources,
        output_resources={"coverage_ht": release_coverage(public=False, test=test)},
    )
    export_coverage_tsv = PipelineStepResourceCollection(
        "--export-coverage-tsv",
        output_resources={
            "coverage_tsv": release_coverage_tsv_path("exomes", test=test)
        },
        pipeline_input_steps=[compute_coverage_ht],
    )

    # Add all steps to the coverage pipeline resource collection.
    coverage_pipeline.add_steps(
        {
            "compute_coverage_ht": compute_coverage_ht,
            "export_coverage_tsv": export_coverage_tsv,
        }
    )

    return coverage_pipeline


def main(args):
    """Compute coverage statistics, including mean, median_approx, and coverage over certain DPs."""
    hl.init(
        log="/coverage.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )

    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    test = args.test
    overwrite = args.overwrite

    coverage_resources = get_coverage_resources(
        test=test,
        overwrite=overwrite,
        calling_interval_name=args.calling_interval_name,
        calling_interval_padding=args.calling_interval_padding,
    )

    try:
        if args.compute_coverage_ht:
            logger.info("Running compute coverage...")
            res = coverage_resources.compute_coverage_ht
            res.check_resource_existence()
            # Read in context Table.
            ref_ht = vep_context.versions["105"].ht()
            if test:
                ref_ht = ref_ht._filter_partitions(range(50))

            # Retain only 'locus' annotation from context Table
            ref_ht = ref_ht.key_by("locus").select().distinct()

            # Read in VDS.
            vds = get_gnomad_v4_vds(
                release_only=True,
                test=test,
                filter_partitions=range(2) if test else None,
                annotate_meta=True,
            )

            if args.stratify_by_ukb_and_platform:
                meta_expr = vds.variant_data.meta
                strata = {
                    "subset": hl.if_else(
                        meta_expr.project_meta.ukb_sample, "ukb", "non_ukb"
                    ),
                    "platform": meta_expr.platform_inference.qc_platform,
                }
            else:
                strata = None

            # Compute coverage stats.
            coverage_ht = compute_coverage_stats(
                vds,
                ref_ht,
                interval_ht=res.interval_ht.ht(),
                strata=strata,
            )

            # Checkpoint Table.
            coverage_ht = coverage_ht.checkpoint(
                hl.utils.new_temp_file("coverage", extension="ht")
            )

            # Naive coalesce and write out final Table.
            coverage_ht = coverage_ht.naive_coalesce(5000)
            coverage_ht.write(res.coverage_ht.path, overwrite=overwrite)

        if args.export_coverage_tsv:
            logger.info("Exporting coverage tsv...")
            res = coverage_resources.export_coverage_tsv
            res.check_resource_existence()
            res.coverage_ht.ht().export(res.coverage_tsv)

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
        "--stratify-by-ukb-and-platform",
        help="Whether to compute coverage stratified by UKB/non-UKB and platform.",
        action="store_true",
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
