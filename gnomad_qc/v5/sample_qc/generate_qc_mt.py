"""
Script to create a joint gnomAD v4 + AoU v8 dense MatrixTable filtered to predetermined QC sites.

For more information about the predetermined QC sites, see:
https://gnomad.broadinstitute.org/news/2023-11-gnomad-v4-0/#high-quality-sites-definition,
https://github.com/broadinstitute/gnomad_qc/blob/main/gnomad_qc/v4/sample_qc/generate_qc_mt.py,
and https://github.com/Nealelab/ccdg_qc/blob/master/scripts/pca_variant_filter.py.

The joint MatrixTable created by this script is used for relatedness and genetic ancestry inference.
"""

import argparse
import logging
from typing import List

import hail as hl

from gnomad_qc.v4.resources import sample_qc as v4_sample_qc
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_checkpoint_path,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import AOU_WGS_BUCKET, WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import sample_id_collisions, samples_to_exclude
from gnomad_qc.v5.resources.sample_qc import get_joint_qc

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("generate_qc_mt")
logger.setLevel(logging.INFO)


def union_aou_mts(
    s_to_exclude: hl.expr.SetExpression,
    ht: hl.Table,
    n_partitions: int = 10000,
    entry_annotations: List[str] = ["GT"],
    overwrite: bool = False,
) -> None:
    """
    Filter AoU ACAF and exome MTs to QC MT sites, remove samples to exclude, and union MTs.

    .. note::
        Both AoU MTs have 145192 partitions.

    :param s_to_exclude: Set of sample IDs to exclude from the AoU MatrixTable.
    :param ht: Table containing the gnomAD QC sites.
    :param n_partitions: Number of desired partitions for the unioned MatrixTable.
        Default is 10000.
    :param entry_annotations: List of entry annotations to keep in the unioned MatrixTable.
        Default is ["GT"].
    :param overwrite: Whether to overwrite output data. Default is False.
    :return: None; function writes unioned MT to temporary path.
    """
    logger.info(
        "Loading AoU ACAF and exome MatrixTables and removing unnecessary annotations..."
    )
    acaf_mt = (
        hl.read_matrix_table(f"gs://{AOU_WGS_BUCKET}/acaf_threshold/splitMT/hail.mt")
        .select_rows()
        .select_entries(*entry_annotations)
    )
    exome_mt = (
        hl.read_matrix_table(f"gs://{AOU_WGS_BUCKET}/exome/splitMT/hail.mt")
        .select_rows()
        .select_entries(*entry_annotations)
    )

    def _filter_aou_mt(
        mt: hl.MatrixTable,
        s_to_exclude: hl.expr.SetExpression,
        ht: hl.Table,
    ) -> hl.MatrixTable:
        """
        Filter AoU MT to gnomAD QC MT sites, remove hard filtered or low quality samples, and filter to `adj`.

        .. note::
            This function uses AoU's `NO_HQ_GENOTYPES` filter to filter to `adj`.

        :param mt: AoU MatrixTable.
        :param s_to_exclude: Set of sample IDs to exclude from the AoU MatrixTable.
        :param ht: Table containing the gnomAD QC sites.
        :return: Filtered AoU MatrixTable.
        """
        # Filter to gnomAD QC sites and use `NO_HQ_GENOTYPES` to filter to `adj`.
        mt = mt.filter_rows(hl.is_defined(ht[mt.row_key]))
        mt = mt.filter_rows(~mt.filters.contains("NO_HQ_GENOTYPES"))

        # Remove samples to exclude and their variant sites.
        mt = mt.filter_cols(~s_to_exclude.contains(mt.s))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        return mt

    logger.info("Filtering AoU ACAF and exome MatrixTables and unioning (on rows)...")
    acaf_mt = _filter_aou_mt(acaf_mt, s_to_exclude, ht)
    exome_mt = _filter_aou_mt(exome_mt, s_to_exclude, ht)
    union_mt = acaf_mt.union_rows(exome_mt)

    logger.info("Decreasing partitions and checkpointing...")
    union_mt = union_mt.naive_coalesce(n_partitions)
    union_mt.write(
        get_checkpoint_path("union_aou_mts", mt=True, environment="rwb"),
        overwrite=overwrite,
    )


def generate_qc_mt(
    gnomad_mt: hl.MatrixTable,
    aou_mt: hl.MatrixTable,
    n_partitions: int = 5000,
) -> hl.MatrixTable:
    """
    Union gnomAD v4 QC MT and AoU v8 MTs.

    :param gnomad_mt: gnomAD v4 QC MatrixTable.
    :param aou_mt: Joint AoU (ACAF + exome) v8 MatrixTable.
    :param n_partitions: Number of desired partitions for the unioned MatrixTable.
        Default is 5000.
    :return: Joint gnomAD v4 (exomes + genomes) + AoU v8 (genomes) QC MT.
    """
    logger.info("Resolving sample ID collisions...")
    sample_collisions = sample_id_collisions.ht()
    gnomad_mt = add_project_prefix_to_sample_collisions(
        t=gnomad_mt,
        sample_collisions=sample_collisions,
        project="gnomad",
    )
    aou_mt = add_project_prefix_to_sample_collisions(
        t=aou_mt,
        sample_collisions=sample_collisions,
        project="aou",
    )

    joint_mt = gnomad_mt.union_cols(aou_mt)
    return joint_mt.naive_coalesce(n_partitions)


def main(args):
    """Create a joint gnomAD + AoU dense MT of a diverse set of variants for relatedness/genetic ancestry PCA."""
    hl.init(
        log="/generate_qc_mt.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test = args.test

    try:
        if args.union_aou_mts:
            logger.info("Generating joint AoU MatrixTable...")
            union_aou_mts(
                s_to_exclude=samples_to_exclude().he(),
                ht=get_joint_qc(test=test).mt().rows(),
                n_partitions=args.n_partitions,
                overwrite=overwrite,
            )

        if args.generate_qc_mt:
            logger.info("Loading gnomAD v4 QC MatrixTable...")
            # NOTE: Dropping extra column and entry annotations because `union_cols` requires
            # that the column keys/schemas and entry schemas match across the two MTs.
            # Also dropping row annotations because they are not needed for the joint
            # MT.
            gnomad_mt = (
                v4_sample_qc.get_joint_qc(test=test)
                .mt()
                .select_globals()
                .select_rows()
                .select_cols()
                .select_entries("GT")
            )

            logger.info("Loading AoU joint ACAF + exome MatrixTable...")
            aou_mt = hl.read_matrix_table(
                get_checkpoint_path("union_aou_mts", mt=True, environment="rwb")
            )

            logger.info("Generating joint gnomAD + AoU QC MatrixTable...")
            joint_mt = generate_qc_mt(
                gnomad_mt=gnomad_mt,
                aou_mt=aou_mt,
                n_partitions=args.n_partitions,
            )
            joint_mt.write(
                get_joint_qc(test=test).path,
                overwrite=overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("generate_qc_mt", environment="rwb"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Read in a test dataset rather than the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--union-aou-mts",
        help=(
            "Filter AoU ACAF and exome MTs to QC MT sites, remove samples to exclude, and union MTs."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--generate-qc-mt",
        help="Generate joint gnomAD v4 + AoU v8 dense MatrixTable.",
        action="store_true",
    )
    parser.add_argument(
        "--n-partitions",
        help="Desired number of partitions for data output by step.",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite data.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
