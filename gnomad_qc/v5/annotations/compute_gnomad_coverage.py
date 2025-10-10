"""Script to compute coverage, allele number, and quality histograms on gnomAD v5 genomes (gnomAD v4 minus consent drop samples)."""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import (
    build_freq_stratification_list,
    generate_freq_group_membership_array,
)

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v5.resources.annotations import group_membership
from gnomad_qc.v5.resources.basics import get_gnomad_v5_genomes_vds, get_logging_path

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("gnomad_coverage_and_an")
logger.setLevel(logging.INFO)


def get_group_membership_ht(meta_ht: hl.Table) -> hl.Table:
    """
    Get group membership HT for gnomAD v5 genomes.

    This is different from the group membership HT for v4 due to the consent drop samples.

    :param meta_ht: Meta HT.
    :return: Group membership HT.
    """
    # Filter to v4 release samples and drop consent drop samples.
    # NOTE: Not using v5 project meta because this script will be run in Dataproc.
    ht = meta_ht.filter(
        meta_ht.release
        & (
            (meta_ht.project_meta.research_project_key != "RP-1061")
            & (meta_ht.project_meta.research_project_key != "RP-1411")
        )
    )
    return generate_freq_group_membership_array(
        ht,
        build_freq_stratification_list(
            sex_expr=ht.sex_imputation.sex_karyotype,
            gen_anc_expr=ht.population_inference.pop,
        ),
    )


def main(args):
    """Compute coverage, allele number, and quality histograms on gnomAD v5 genomes (gnomAD v4 minus consent drop samples)."""
    hl.init(
        log="compute_gnomad_coverage.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")

    test_2_partitions = args.test_2_partitions
    test_chr22_chrx_chry = args.test_chr22_chrx_chry
    test = test_2_partitions or test_chr22_chrx_chry
    overwrite = args.overwrite
    n_partitions = args.n_partitions

    try:
        if args.write_group_membership_ht:
            logger.info("Writing group membership HT...")
            meta_ht_path = meta(data_type="genomes").path
            group_membership_ht_path = group_membership(
                test=args.test, data_set="gnomad"
            ).path
            check_resource_existence(
                input_step_resources={"meta_ht": meta_ht_path},
                output_step_resources={"group_membership_ht": group_membership_ht_path},
                overwrite=overwrite,
            )

            ht = get_group_membership_ht(hl.read_table(meta_ht_path))
            ht.write(group_membership_ht_path, overwrite=overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compute_gnomad_coverage", environment="dataproc"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--data-set",
        help="Data set to compute coverage for.",
        type=str,
        choices=["aou", "gnomad"],
        default="aou",
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
    group_membership_args = parser.add_argument_group(
        "Get group membership HT.",
    )
    group_membership_args.add_argument(
        "--write-group-membership-ht",
        help="Write group membership HT.",
        action="store_true",
    )
    group_membership_args.add_argument(
        "--test",
        help="Write test group membership HT to test path.",
        action="store_true",
    )
    parser.add_argument(
        "--compute-coverage-an-ht",
        help="Compute the coverage and allele number HT.",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
