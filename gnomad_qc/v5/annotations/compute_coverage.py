"""Script to compute coverage statistics on gnomAD v5 genomes."""

import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.annotations import (
    build_freq_stratification_list,
    generate_freq_group_membership_array,
)
from gnomad.utils.sparse_mt import compute_coverage_stats, compute_stats_per_ref_site

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.annotations.compute_coverage import (
    compute_an_and_qual_hists_per_ref_site,
)
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds
from gnomad_qc.v4.resources.meta import meta

# TODO: Switch from v4>v5 once v5 sample QC is complete
from gnomad_qc.v5.resources.basics import get_logging_path  # get_aou_vds
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET

# from gnomad_qc.resources.meta import meta
from gnomad_qc.v5.resources.release import (
    release_all_sites_an_tsv_path,
    release_coverage_path,
    release_coverage_tsv_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("coverage")
logger.setLevel(logging.INFO)


def get_genomes_group_membership_ht(
    meta_ht: hl.Table,
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    :param meta_ht: Metadata HT.
    :return: Group membership HT.
    """
    # Filter to release samples.
    meta_ht = meta_ht.filter(meta_ht.release)

    # Create group membership HT.
    ht = generate_freq_group_membership_array(
        meta_ht,
        build_freq_stratification_list(
            sex_expr=meta_ht.sex_imputation.sex_karyotype,
            pop_expr=meta_ht.gen_anc_inference.gen_anc,
        ),
    )

    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: hl.dict(
                d.items().map(lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x))
            )
        ),
    )

    return ht


def main(args):
    """Compute coverage statistics, including mean, median_approx, and coverage over certain DPs."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/compute_coverage.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")

    # TODO: Remove this?
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    test_2_partitions = args.test_2_partitions
    test_chr22_chrx_chry = args.test_chr22_chrx_chry
    test = test_2_partitions or test_chr22_chrx_chry
    overwrite = args.overwrite
    data_type = args.data_type
    n_partitions = args.n_partitions

    try:
        if args.compute_coverage_ht or args.compute_all_sites_an_and_qual_hist_ht:
            # Read in context Table.
            ref_ht = vep_context.versions["105"].ht()
            if test_chr22_chrx_chry:
                chrom = ["chr22", "chrX", "chrY"]
                ref_ht = hl.filter_intervals(
                    ref_ht, [hl.parse_locus_interval(c) for c in chrom]
                )
            elif test_2_partitions:
                ref_ht = ref_ht._filter_partitions(range(5))

            # Retain only 'locus' annotation from context Table.
            ref_ht = ref_ht.key_by("locus").select().distinct()

            # Read in VDS.
            vds = get_gnomad_v4_genomes_vds(
                release_only=True,
                test=test,
                filter_partitions=range(2) if test_2_partitions else None,
                annotate_meta=True,
                chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
            )
            # vds = get_aou_vds(
            #    release_only=True,
            #    test=test,
            #    filter_partitions=range(2) if test_2_partitions else None,
            #    annotate_meta=True,
            #    chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
            # )

        if args.compute_coverage_ht:
            cov_ht_path = release_coverage_path(
                public=False, test=test, include_meta=True, coverage_type="coverage"
            )
            check_resource_existence(output_step_resources={"coverage_ht": cov_ht_path})

            logger.info("Running compute coverage...")
            coverage_ht = compute_coverage_stats(
                vds,
                ref_ht,
            )

            # Checkpoint Table.
            coverage_ht = coverage_ht.checkpoint(
                hl.utils.new_temp_file("coverage", "ht")
            )

            # Naive coalesce and write out final Table.
            coverage_ht = coverage_ht.naive_coalesce(n_partitions)
            coverage_ht.write(cov_ht_path, overwrite=overwrite)

        if args.compute_all_sites_an_and_qual_hist_ht:
            an_ht_path = release_coverage_path(
                public=False,
                test=test,
                include_meta=True,
                coverage_type="allele_number",
            )
            check_resource_existence(
                output_step_resources={"all_sites_an_ht": an_ht_path}
            )

            logger.info(
                "Running compute all sites allele number and quality histograms HT..."
            )
            group_membership_ht = get_genomes_group_membership_ht(
                meta().ht(),
            )
            group_membership_ht = group_membership_ht.checkpoint(
                hl.utils.new_temp_file("group_membership", "ht")
            )

            an_ht = compute_an_and_qual_hists_per_ref_site(
                vds,
                ref_ht,
                group_membership_ht=group_membership_ht,
            )
            an_ht = an_ht.checkpoint(hl.utils.new_temp_file("an", "ht"))
            # Naive coalesce and write out the intermediate HT.
            an_ht = an_ht.naive_coalesce(n_partitions)
            an_ht.write(an_ht_path.path, overwrite=overwrite)

        if args.export_coverage_release_files:
            cov_ht_path = release_coverage_path(
                public=False, test=test, include_meta=True, coverage_type="coverage"
            )
            release_ht_path = release_coverage_path(
                public=False, test=test, include_meta=False, coverage_type="coverage"
            )
            cov_tsv_path = release_coverage_tsv_path(test=test)
            check_resource_existence(
                input_step_resources={"cov_ht": cov_ht_path},
                output_step_resources={
                    "cov_release_ht": release_ht_path,
                    "cov_tsv": cov_tsv_path,
                },
            )

            logger.info("Exporting coverage tsv...")
            ht = hl.read_table(cov_ht_path)
            if "coverage_stats" in ht.row:
                ht = ht.select(
                    **{k: ht.coverage_stats[0][k] for k in ht.coverage_stats[0]}
                )
            ht = ht.drop("coverage_stats_meta", "coverage_stats_meta_sample_count")
            ht = ht.checkpoint(release_ht_path, overwrite=overwrite)
            ht.export(cov_tsv_path)

        if args.export_all_sites_an_release_files:
            an_ht_path = release_coverage_path(
                public=False,
                test=test,
                include_meta=True,
                coverage_type="allele_number",
            )
            release_ht_path = release_coverage_path(
                public=False,
                test=test,
                include_meta=False,
                coverage_type="allele_number",
            )
            an_tsv_path = release_coverage_tsv_path(test=test)
            check_resource_existence(
                input_step_resources={"an_ht": an_ht_path},
                output_step_resources={
                    "an_release_ht": release_ht_path,
                    "an_tsv": an_tsv_path,
                },
            )

            # Select only the AN and write out final Table.
            logger.info("Exporting all sites AN HT...")
            ht = hl.read_table(an_ht_path)
            ht = ht.select("AN")

            ht = ht.checkpoint(release_ht_path, overwrite=overwrite)

            logger.info("Exporting all sites AN tsv...")
            # Only export the adj AN for all release samples.
            ht = ht.annotate(AN=ht.AN[0])
            ht.export(an_tsv_path)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compute_coverage", environment="rwb"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--test-2-partitions",
        help=(
            "Whether to run a test using only the first 2 partitions of the VDS test"
            " dataset."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-chr22-chrx-chry",
        help=(
            "Whether to run a test using only the chr22, chrX, and chrY chromosomes of"
            " the VDS test dataset."
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
        "--compute-all-sites-an-and-qual-hist-ht",
        help="Compute the all sites allele number and quality histogram HT.",
        action="store_true",
    )
    parser.add_argument(
        "--export-coverage-release-files",
        help="Exports coverage release HT and TSV file.",
        action="store_true",
    )
    parser.add_argument(
        "--export-all-sites-an-release-files",
        help="Export all sites AN release HT and TSV file.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
