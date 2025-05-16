"""Script to determine samples that fail hard filtering thresholds."""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.utils.annotations import bi_allelic_expr

from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_aou_vds,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import sample_id_collisions
from gnomad_qc.v5.resources.sample_qc import get_sample_qc

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hard_filters")
logger.setLevel(logging.INFO)


def compute_aou_sample_qc(
    n_partitions: int = 500,
    test: bool = False,
) -> hl.Table:
    """
    Perform sample QC on AoU VDS.

    .. note::

        We are not including the `n_alt_alleles_strata` parameter in this function—as we did for v4 exomes—
        because the distribution of alternate alleles in whole genome sequencing data is not as skewed as in exomes.
        For example, in AoU v8 genomes, 77.06% of variants are bi-allelic, compared to 76.65% in v4 genomes and
        only 35.77% in v4 exomes.

    :param n_partitions: Number of partitions to use when writing the sample QC table.
    :param test: If true, test the function on a smaller subset of the data.
    :return: Table containing sample QC metrics
    """
    if test:
        logger.info("Filtering to first interval that contains n_alt_alleles > 100...")
        vds = get_aou_vds(filter_intervals=["chr1:10440-10626"], split=True)
    else:
        logger.info("Loading AoU VDS...")
        vds = get_aou_vds(
            autosomes_only=True,
            split=True,
        )

    logger.info(
        "Excluding telomeres and centromeres from VDS (redundant but acts as a safety check)..."
    )
    # AoU pipeline has already filtered out the centromeres and telomeres, but
    # this serves as an additional safeguard.
    vds = hl.vds.filter_intervals(
        vds, intervals=telomeres_and_centromeres.ht(), keep=False
    )

    logger.info("Excluding loci with more than 100 alternative alleles...")
    vmt = vds.variant_data

    if test:
        logger.info("Number of rows before filtering: %s", vmt.count_rows())
    vmt = vmt.filter_rows(vmt.n_unsplit_alleles < 101)
    if test:
        logger.info("Number of rows after filtering: %s", vmt.count_rows())

    logger.info("Computing sample QC metrics...")
    sample_qc_ht = compute_stratified_sample_qc(
        vmt,
        strata={
            "bi_allelic": bi_allelic_expr(vmt),
            "multi_allelic": ~bi_allelic_expr(vmt),
        },
        tmp_ht_prefix=get_sample_qc(test=test).path[:-3],
        gt_col="GT",
    )

    return sample_qc_ht.naive_coalesce(n_partitions)


def main(args):
    """Determine samples that fail hard filtering thresholds."""
    hl.init(tmp_dir=f"{WORKSPACE_BUCKET}/tmp/4_day")
    hl.default_reference("GRCh38")

    test = args.test
    overwrite = args.overwrite

    if args.aou_sample_qc:
        ht = compute_aou_sample_qc(
            n_partitions=args.sample_qc_n_partitions,
            test=test,
        )
        ht = add_project_prefix_to_sample_collisions(
            ht, project="aou", sample_collisions=sample_id_collisions.ht()
        )
        ht.write(get_sample_qc(test=test).path, overwrite=overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all Matrixtables/Tables. (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Use the AoU test dataset instead of the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--aou-sample-qc",
        help="Compute Hail's VDS sample QC metrics on AoU.",
        action="store_true",
    )
    parser.add_argument(
        "--sample-qc-n-partitions",
        help="Number of partitions to use when writing the sample QC table.",
        type=int,
        default=500,
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
