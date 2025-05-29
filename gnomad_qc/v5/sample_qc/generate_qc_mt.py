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

import hail as hl
from gnomad.utils.filtering import filter_to_autosomes

from gnomad_qc.v4.resources import sample_qc as v4_sample_qc
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    aou_acaf_mt,
    aou_exome_mt,
    get_logging_path,
    get_samples_to_exclude,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import sample_id_collisions
from gnomad_qc.v5.resources.sample_qc import (
    get_aou_mt_union,
    get_joint_qc,
    hard_filtered_samples,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("generate_qc_mt")
logger.setLevel(logging.INFO)


def union_aou_mts(
    acaf_mt: hl.MatrixTable,
    exome_mt: hl.MatrixTable,
    s_to_exclude: hl.expr.SetExpression,
    ht: hl.Table,
    test: bool = False,
) -> hl.MatrixTable:
    """
    Filter AoU ACAF and exome MTs to QC MT sites, remove samples to exclude, filter to `adj`, and union MTs.

    .. note::

        - Both v8 AoU MTs have 145192 partitions.
        - This function unions both MTs to reduce the number of missing sites in the joint MT; after unioning these two MTs,
            only 5080 sites from the v4 QC MT are missing from the joint AoU MT.
        - The v7 AoU exome MT was missing 35,621 sites from the v4 QC MT, and the v7 ACAF MT was missing 65,935 sites.
            After unioning the v7 MTs, the joint MT was only missing 8203 sites.

    :param acaf_mt: AoU ACAF MatrixTable.
    :param exome_mt: AoU exome MatrixTable.
    :param s_to_exclude: Set of sample IDs to exclude from the AoU MatrixTable.
    :param ht: Table containing the gnomAD QC sites.
    :param test: Whether to filter to the first 2 partitions for testing.
        Default is False.
    :return: MatrixTable containing the union of AoU ACAF and exome MTs.
    """

    def _filter_aou_mt(
        mt: hl.MatrixTable,
    ) -> hl.MatrixTable:
        """
        Filter AoU MT to autosomal gnomAD QC MT sites, remove hard filtered or low quality samples, and filter to `adj`.

        .. note::

            This function uses AoU's `NO_HQ_GENOTYPES` filter to filter to `adj`.

        :param mt: AoU MatrixTable.
        :return: Filtered AoU MatrixTable.
        """
        if test:
            logger.info("Filtering AoU MTs to the first 2 partitions for testing...")
            mt = mt._filter_partitions(range(2))

        # Filter to gnomAD QC sites and use `NO_HQ_GENOTYPES` to filter to `adj`.
        mt = filter_to_autosomes(mt)
        mt = mt.semi_join_rows(ht)
        mt = mt.filter_rows(~mt.filters.contains("NO_HQ_GENOTYPES"))
        mt = mt.select_rows()

        # Remove samples to exclude and their variant sites.
        mt = mt.filter_cols(~s_to_exclude.contains(mt.s))
        mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
        return mt

    logger.info("Filtering AoU ACAF and exome MatrixTables and unioning (on rows)...")
    # NOTE: Add checkpoints after each MT if run into errors in AoU workbench.
    acaf_mt = _filter_aou_mt(acaf_mt)
    exome_mt = _filter_aou_mt(exome_mt)
    return acaf_mt.union_rows(exome_mt)


def generate_qc_mt(
    gnomad_mt: hl.MatrixTable,
    aou_mt: hl.MatrixTable,
) -> hl.MatrixTable:
    """
    Union gnomAD v4 QC MT and AoU v8 MTs.

    :param gnomad_mt: gnomAD v4 QC MatrixTable.
    :param aou_mt: Joint AoU (ACAF + exome) v8 MatrixTable.
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

    return gnomad_mt.union_cols(aou_mt)


def main(args):
    """Create a joint gnomAD + AoU dense MT of a diverse set of variants for relatedness/genetic ancestry PCA."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/generate_qc_mt.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")

    n_partitions = args.n_partitions
    overwrite = args.overwrite
    test = args.test
    v4_qc_sites_ht = v4_sample_qc.get_joint_qc().mt().rows()

    try:
        if args.union_aou_mts:

            logger.info(
                "Loading AoU ACAF and exome MatrixTables and removing unnecessary entry annotations..."
            )
            entry_annotations = args.entry_annotations
            acaf_mt = aou_acaf_mt.mt().select_entries(*entry_annotations)
            exome_mt = aou_exome_mt.mt().select_entries(*entry_annotations)

            logger.info("Generating joint AoU MatrixTable...")
            mt = union_aou_mts(
                acaf_mt=acaf_mt,
                exome_mt=exome_mt,
                s_to_exclude=get_samples_to_exclude(
                    filter_samples=hard_filtered_samples.ht()
                ),
                ht=v4_qc_sites_ht,
                test=test,
            )

            logger.info("Decreasing partitions and checkpointing...")
            mt = mt.naive_coalesce(n_partitions)
            mt.write(get_aou_mt_union.path, overwrite=overwrite)

        if args.check_missing_sites:
            logger.info(
                "Checking how many sites from v4 QC MT are missing from joint AoU MT..."
            )
            aou_ht = get_aou_mt_union().mt().rows()
            ht = v4_qc_sites_ht.select_rows()
            missing_sites = ht.anti_join(aou_ht)
            logger.info(
                "Number of QC MT sites missing from joint AoU MT: %i",
                missing_sites.count(),
            )

        if args.generate_qc_mt:
            logger.info("Loading gnomAD v4 QC MatrixTable...")
            # NOTE: Dropping extra column and entry annotations because `union_cols` requires
            # that the column keys/schemas and entry schemas match across the two MTs.
            # Also dropping row annotations because they are not needed for the joint
            # MT.
            gnomad_mt = (
                v4_sample_qc.get_joint_qc()
                .mt()
                .select_globals()
                .select_rows()
                .select_cols()
                .select_entries("GT")
            )

            if test:
                logger.info(
                    "Filtering gnomAD MT to the first 2 partitions for testing..."
                )
                gnomad_mt = gnomad_mt._filter_partitions(range(2))

            logger.info("Loading AoU joint ACAF + exome MatrixTable...")
            aou_mt = hl.read_matrix_table(
                get_checkpoint_path("union_aou_mts", mt=True, environment="rwb")
            )

            logger.info("Generating joint gnomAD + AoU QC MatrixTable...")
            joint_mt = generate_qc_mt(
                gnomad_mt=gnomad_mt,
                aou_mt=aou_mt,
            )
            joint_mt = joint_mt.naive_coalesce(n_partitions)
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
        help="Filter to the first 2 partitions for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--entry-annotations",
        help=(
            "Entry annotations to keep in the unioned AoU MatrixTable. "
            "Default is ['GT']."
        ),
        type=str,
        default=["GT"],
        nargs="+",
    )
    parser.add_argument(
        "--union-aou-mts",
        help=(
            "Filter AoU ACAF and exome MTs to QC MT sites, remove samples to exclude, and union MTs."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--check-missing-sites",
        help=("Check how many sites from v4 QC MT are missing from joint AoU MT."),
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
