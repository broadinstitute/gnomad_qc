"""Script to determine samples that fail hard filtering thresholds."""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.utils.annotations import bi_allelic_expr

from gnomad_qc.v5.resources.basics import get_aou_vds, get_logging_path
from gnomad_qc.v5.resources.sample_qc import get_aou_sample_qc

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hard_filters")
logger.setLevel(logging.INFO)


def compute_aou_sample_qc(
    n_partitions: int = 500,
    test: bool = False,
) -> hl.Table:
    """
    Perform sample QC on the VDS.

    :param n_partitions: Number of partitions to use when writing the sample QC table.
    :param test: If true, test the function on a smaller subset of the data.
    :return: Table containing sample QC metrics
    :rtype: hl.Table
    """
    logger.info("Computing sample QC")

    # Filter the first 100 samples for testing or filter the first 10 partitions?
    vds = get_aou_vds(
        autosomes_only=True, split=True, n_partitions=10 if test else None
    )

    # Remove centromeres and telomeres
    vds = hl.vds.filter_intervals(
        vds, intervals=telomeres_and_centromeres.ht(), keep=False
    )

    sample_qc_ht = compute_stratified_sample_qc(
        vds,
        strata={
            "bi_allelic": bi_allelic_expr(vds.variant_data),
            "multi_allelic": ~bi_allelic_expr(vds.variant_data),
        },
        tmp_ht_prefix=get_aou_sample_qc().path[:-3],
        gt_col="GT",
    )

    return sample_qc_ht.repartition(n_partitions)


def main(args):
    """Determine samples that fail hard filtering thresholds."""
    hl.init(
        log="/hard_filters.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default.
    hl._set_flags(use_new_shuffle="1")

    calling_interval_name = args.calling_interval_name
    calling_interval_padding = args.calling_interval_padding
    test = args.test
    overwrite = args.overwrite

    try:
        if args.aou_sample_qc:
            compute_aou_sample_qc(
                n_partitions=args.sample_qc_n_partitions,
                test=test,
            ).write(get_aou_sample_qc(test=test).path, overwrite=overwrite)

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("hard_filters"))


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
        help="Use the v4 test dataset instead of the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--aou-sample-qc",
        help="Compute Hail's VDS sample QC metrics on AoU.",
        action="store_true",
    )


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
