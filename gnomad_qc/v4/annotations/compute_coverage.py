"""Script to compute coverage statistics on gnomAD v4 exomes."""
import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.annotations import (
    build_freq_stratification_list,
    generate_freq_group_membership_array,
)
from gnomad.utils.sparse_mt import (
    compute_allele_number_per_ref_site,
    compute_coverage_stats,
)

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.resources.annotations import get_downsampling
from gnomad_qc.v4.resources.basics import calling_intervals, get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import (
    release_all_sites_an,
    release_coverage,
    release_coverage_tsv_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("coverage")
logger.setLevel(logging.INFO)


def get_group_membership_ht(
    meta_ht: hl.Table,
    ds_ht: hl.Table,
    non_ukb_ds_ht: hl.Table,
) -> hl.Table:
    """
    Get group membership HT for all sites allele number stratification.

    :param meta_ht: Metadata HT.
    :param ds_ht: Full frequency downsampling HT.
    :param non_ukb_ds_ht: Non-UKB frequency downsampling HT.
    :return: Group membership HT.
    """
    # Filter to release samples.
    meta_ht = meta_ht.filter(meta_ht.release)

    # Filter to non-UKB samples.
    non_ukb_meta_ht = meta_ht.filter(~meta_ht.project_meta.ukb_sample)

    # Create group membership HT.
    ht = generate_freq_group_membership_array(
        meta_ht,
        build_freq_stratification_list(
            sex_expr=meta_ht.sex_imputation.sex_karyotype,
            pop_expr=meta_ht.population_inference.pop,
            downsampling_expr=ds_ht[meta_ht.key].downsampling,
        ),
        downsamplings=hl.eval(ds_ht.downsamplings),
        ds_pop_counts=hl.eval(ds_ht.ds_pop_counts),
    )

    # Add non-UKB group membership to HT.
    non_ukb_ht = generate_freq_group_membership_array(
        non_ukb_meta_ht,
        build_freq_stratification_list(
            sex_expr=non_ukb_meta_ht.sex_imputation.sex_karyotype,
            pop_expr=non_ukb_meta_ht.population_inference.pop,
            downsampling_expr=non_ukb_ds_ht[non_ukb_meta_ht.key].downsampling,
        ),
        downsamplings=hl.eval(non_ukb_ds_ht.downsamplings),
        ds_pop_counts=hl.eval(non_ukb_ds_ht.ds_pop_counts),
    )
    n_non_ukb = hl.eval(hl.len(non_ukb_ht.freq_meta))
    ht = ht.annotate(
        group_membership=ht.group_membership.extend(
            hl.coalesce(
                non_ukb_ht[ht.key].group_membership,
                hl.range(n_non_ukb).map(lambda x: False),
            )
        )
    )

    # Add non-UKB group metadata and sample counts to HT.
    non_ukb_globals = non_ukb_ht.index_globals()
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.extend(
            non_ukb_globals.freq_meta.map(
                lambda d: hl.dict(d.items().append(("subset", "non_ukb")))
            )
        ).map(
            lambda d: hl.dict(
                d.items().map(lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x))
            )
        ),
        freq_meta_sample_count=ht.freq_meta_sample_count.extend(
            non_ukb_globals.freq_meta_sample_count
        ),
    )

    return ht


def adjust_interval_padding(ht: hl.Table, padding: int) -> hl.Table:
    """
    Adjust interval padding in HT.

    .. warning::

        This function can lead to overlapping intervals, so it is not recommended for
        most applications. For example, it can be used to filter a variant list to all
        variants within the returned interval list, but would not work for getting an
        aggregate statistic for each interval if the desired output is independent
        statistics.

    :param ht: HT to adjust.
    :param padding: Padding to use.
    :return: HT with adjusted interval padding.
    """
    return ht.key_by(
        interval=hl.locus_interval(
            ht.interval.start.contig,
            ht.interval.start.position - padding,
            ht.interval.end.position + padding,
            reference_genome=ht.interval.start.dtype.reference_genome,
        )
    )


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
    if calling_interval_name is not None and calling_interval_padding is not None:
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
    compute_allele_number_ht = PipelineStepResourceCollection(
        "--compute-allele-number-ht",
        input_resources={
            **coverage_input_resources,
            **{
                "metadata Table": {"meta_ht": meta()},
                "downsampling Tables": {
                    "ds_ht": get_downsampling(),
                    "non_ukb_ds_ht": get_downsampling(subset="non_ukb"),
                },
            },
        },
        output_resources={
            "allele_number_ht": release_all_sites_an(public=False, test=test)
        },
    )
    export_coverage_files = PipelineStepResourceCollection(
        "--export-release-files",
        output_resources={
            "coverage_tsv": release_coverage_tsv_path("exomes", test=test),
            "release_ht": release_coverage(public=False, test=test, stratify=False),
        },
        pipeline_input_steps=[compute_coverage_ht],
    )

    # Add all steps to the coverage pipeline resource collection.
    coverage_pipeline.add_steps(
        {
            "compute_coverage_ht": compute_coverage_ht,
            "compute_allele_number_ht": compute_allele_number_ht,
            "export_coverage_files": export_coverage_files,
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
    n_partitions = args.n_partitions
    calling_interval_padding = args.calling_interval_padding
    adjust_padding = calling_interval_padding not in [0, 50]

    coverage_resources = get_coverage_resources(
        test=test,
        overwrite=overwrite,
        calling_interval_name=args.calling_interval_name,
        # We only have interval lists with 0 or 50bp padding, so if a different padding
        # is requested, use 0bp and adjust it based on the desired padding.
        calling_interval_padding=0 if adjust_padding else calling_interval_padding,
    )

    try:
        if args.compute_coverage_ht or args.compute_allele_number_ht:
            # Read in context Table.
            ref_ht = vep_context.versions["105"].ht()
            if test:
                ref_ht = ref_ht._filter_partitions(range(50))

            # Retain only 'locus' annotation from context Table.
            ref_ht = ref_ht.key_by("locus").select().distinct()

            # Read in VDS.
            vds = get_gnomad_v4_vds(
                release_only=True,
                test=test,
                filter_partitions=range(2) if test else None,
                annotate_meta=True,
            )

        if args.compute_coverage_ht:
            logger.info("Running compute coverage...")
            res = coverage_resources.compute_coverage_ht
            res.check_resource_existence()

            if args.stratify_by_ukb_and_platform:
                meta_expr = vds.variant_data.meta
                ukb_expr = meta_expr.project_meta.ukb_sample
                strata = [
                    {"ukb_strata": hl.if_else(ukb_expr, "ukb", "non_ukb")},
                    {
                        "subset": hl.or_missing(~ukb_expr, "non_ukb"),
                        "platform": meta_expr.platform_inference.qc_platform,
                    },
                ]
            else:
                strata = None

            # Compute coverage stats.
            coverage_ht = compute_coverage_stats(
                vds,
                ref_ht,
                interval_ht=(
                    adjust_interval_padding(
                        res.interval_ht.ht(), calling_interval_padding
                    )
                    if adjust_padding
                    else res.interval_ht.ht()
                ),
                strata_expr=strata,
            )

            # Checkpoint Table.
            if args.stratify_by_ukb_and_platform:
                coverage_ht = coverage_ht.annotate_globals(
                    coverage_stats_meta=coverage_ht.coverage_stats_meta.map(
                        lambda x: hl.dict(
                            x.items().map(
                                lambda m: hl.if_else(
                                    m[0] == "ukb_strata", ("subset", m[1]), m
                                )
                            )
                        )
                    )
                )
            coverage_ht = coverage_ht.checkpoint(
                hl.utils.new_temp_file("coverage", "ht")
            )

            # Naive coalesce and write out final Table.
            coverage_ht = coverage_ht.naive_coalesce(n_partitions)
            coverage_ht.write(res.coverage_ht.path, overwrite=overwrite)

        if args.compute_allele_number_ht:
            logger.info("Running compute all sites allele number HT...")
            res = coverage_resources.compute_allele_number_ht
            res.check_resource_existence()

            an_ht = compute_allele_number_per_ref_site(
                vds,
                ref_ht,
                interval_ht=(
                    adjust_interval_padding(
                        res.interval_ht.ht(), calling_interval_padding
                    )
                    if adjust_padding
                    else res.interval_ht.ht()
                ),
                group_membership_ht=get_group_membership_ht(
                    res.meta_ht.ht(),
                    res.ds_ht.ht(),
                    res.non_ukb_ds_ht.ht(),
                ),
            )

            an_ht = an_ht.checkpoint(hl.utils.new_temp_file("an", "ht"))

            # Naive coalesce and write out final Table.
            an_ht = an_ht.naive_coalesce(n_partitions)
            an_ht.write(res.allele_number_ht.path, overwrite=overwrite)

        if args.export_release_files:
            logger.info("Exporting coverage tsv...")
            res = coverage_resources.export_coverage_files
            res.check_resource_existence()
            ht = res.coverage_ht.ht()
            if "coverage_stats" in ht.row:
                ht = ht.select(
                    **{k: ht.coverage_stats[0][k] for k in ht.coverage_stats[0]}
                )
            ht = ht.drop("coverage_stats_meta", "coverage_stats_meta_sample_count")
            ht = ht.checkpoint(res.release_ht.path, overwrite=overwrite)
            ht.export(res.coverage_tsv)

    finally:
        hl.copy_log(f"gs://gnomad-tmp-4day/coverage/compute_coverage.log")


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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
        "--n-partitions",
        help="Number of partitions to use for the output Table.",
        type=int,
        default=5000,
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
        "--compute-allele-number-ht",
        help="Compute the all sites allele number HT.",
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
            "Number of base pair padding to use on the calling intervals. One of 0, "
            " 50, or 150 bp."
        ),
        type=int,
        choices=[0, 50, 150],
        default=50,
    )
    parser.add_argument(
        "--export-release-files", help="Exports coverage TSV file.", action="store_true"
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
