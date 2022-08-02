import argparse
import logging
from typing import Union

import hail as hl

from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad.utils.annotations import annotate_adj, get_adj_expr
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_logging_path,
    get_gnomad_v4_vds,
)
from gnomad_qc.v4.resources.sample_qc import (
    get_dense_predetermined_qc,
    hard_filtered_samples_no_sex,
    # hard_filtered_samples,
    predetermined_qc_sites,
    qc,
)
from gnomad_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_qc_data")
logger.setLevel(logging.INFO)


def create_filtered_dense_mt(
    mtds: Union[hl.vds.VariantDataset, hl.MatrixTable],
    test: bool = False,
    split: bool = False,
) -> hl.MatrixTable:
    """
    Filter a sparse MatrixTable or VariantDataset to a set of predetermined QC sites and return a dense MatrixTable.

    :param mtds: Input MatrixTable or VariantDataset.
    :param test: Whether to filter `mtds` to the first 20 partitions.
    :param split: Whether `mtds` should have multi-allelics split before filtering variants.
    :return: Filtered and densified MatrixTable.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if test:
        logger.info("Filtering to the first 20 partitions")
        if is_vds:
            mtds = hl.vds.VariantDataset(
                mtds.reference_data._filter_partitions(range(20)),
                mtds.variant_data._filter_partitions(range(20)),
            )
        else:
            mtds = mtds._filter_partitions(range(20))

    if is_vds:
        vds = mtds
    else:
        logger.info("Converting MatrixTable to VariantDataset...")
        mtds = mtds.select_entries("END", "LA", "LGT", "GQ", "DP", "LAD")
        vds = hl.vds.VariantDataset.from_merged_representation(mtds)

    if split:
        logger.info("Performing split-multi...")
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info("Filtering variants to predetermined QC variants...")
    vds = hl.vds.filter_variants(vds, predetermined_qc_sites.ht())

    logger.info("Densifying data...")
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Writing dense MatrixTable...")
    mt = mt.select_entries("GT", "GQ", "DP", "AD")

    return mt


def generate_qc_mt(
    v3_mt: hl.MatrixTable,
    v4_mt: hl.MatrixTable,
    bi_allelic_only: bool = False,
    min_af: float = 0.0001,
    min_callrate: float = 0.99,
    min_inbreeding_coeff_threshold: float = -0.8,
    ld_r2: float = 0.1,
) -> hl.MatrixTable:
    """
    Generate combined gnomAD v3 and v4 QC MatrixTable for use in relatedness and ancestry inference.

    :param v3_mt: Dense gnomAD v3 MatrixTable filtered to predetermined sites.
    :param v4_mt: Dense gnomAD v4 MatrixTable filtered to predetermined sites.
    :param bi_allelic_only: Whether to filter to bi-allelic variants.
    :param min_af: Minimum variant allele frequency to retain variant in QC MatrixTable.
    :param min_callrate: Minimum variant callrate to retain variant in QC MatrixTable.
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to retain variant in QC MatrixTable.
    :param ld_r2: LD-pruning cutoff.
    :return: MatrixTable of sites that pass QC filters.
    """
    # TODO: Should we only keep v3 release samples or all samples?
    logger.info(
        "Number of (variants, samples) in the v3.1 MatrixTable: %s...", v3_mt.count()
    )
    logger.info(
        "Number of (variants, samples) in the v4 MatrixTable: %s...", v4_mt.count()
    )
    # Remove v4 hard filtered samples
    # v4_mt = v4_mt.anti_join_cols(hard_filtered_samples.ht())
    v4_mt = v4_mt.anti_join_cols(hard_filtered_samples_no_sex.ht())

    samples_in_both = v4_mt.cols().semi_join_cols(v3_mt.cols())
    n_samples_in_both = samples_in_both.count()
    if n_samples_in_both > 0:
        logger.info("The following samples are found in both MatrixTables:")
        samples_in_both.show(n_samples_in_both)
        raise ValueError("Overlapping sample names were found in the MatrixTables!")

    logger.info(
        "Performing a union of the v3.1 and v4 predetermined QC site MatrixTable columns..."
    )
    mt = v3_mt.union_cols(v4_mt)

    logger.info("Annotating with adj and filtering to pre LD-pruned QC sites...")
    mt = annotate_adj(mt)

    cp_path = get_checkpoint_path("combined_pre_ld_prune_qc_mt", mt=True)
    # TODO: As I have it now, v3 and v4 are used to determine AF and callrate, is that OK, or should they be considered
    #  independently and then shared passing variants merged?
    mt = get_qc_mt(
        mt,
        bi_allelic_only=bi_allelic_only,
        min_af=min_af,
        min_callrate=min_callrate,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=None,  # Already filtered from the initial set of QC variants
        apply_hard_filters=False,
        ld_r2=ld_r2,
        checkpoint_path=cp_path,
        filter_lcr=False,  # Already filtered from the initial set of QC variants
        filter_decoy=False,  # Doesn't exist for hg38
        filter_segdup=False,  # Already filtered from the initial set of QC variants
    )
    logger.info(
        "Number of pre LD-pruned QC sites in gnomAD v3 + v4: %d...",
        hl.read_matrix_table(cp_path).count_rows(),
    )

    return mt


def main(args):
    hl.init(
        log="/generate_qc_mt.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default
    hl._set_flags(use_new_shuffle="1")

    overwrite = args.overwrite
    ld_r2 = args.ld_r2
    test = args.test
    n_partitions = args.n_partitions

    try:
        if args.create_v3_filtered_dense_mt:
            # Note: This command removes hard filtered samples
            mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
            mt = create_filtered_dense_mt(mt, test, split=True)
            mt = mt.checkpoint(
                get_dense_predetermined_qc(version="3.1", test=test).path,
                overwrite=overwrite,
            )
            logger.info(
                "Number of predetermined QC variants found in the gnomAD v3 MatrixTable: %d...",
                mt.count_rows(),
            )

        if args.create_v4_filtered_dense_mt:
            # Note: This subset dense MatrixTable was created before the final hard filtering was determined
            # Hard filtering is performed in `generate_qc_mt` before applying variant filters
            vds = get_gnomad_v4_vds(split=True, remove_hard_filtered_samples=False)
            mt = create_filtered_dense_mt(vds, test)
            mt = mt.checkpoint(
                get_dense_predetermined_qc(test=test).path, overwrite=overwrite
            )
            logger.info(
                "Number of predetermined QC variants found in the gnomAD v4 VDS: %d...",
                mt.count_rows(),
            )

        if args.generate_qc_mt:
            v3_mt = hl.read_matrix_table(
                get_dense_predetermined_qc(version="v3.1", test=test).path,
                _n_partitions=n_partitions,
            )
            v4_mt = hl.read_matrix_table(
                get_dense_predetermined_qc(test=test).path,
                _n_partitions=n_partitions,
            )
            mt = generate_qc_mt(
                v3_mt,
                v4_mt,
                bi_allelic_only=args.bi_allelic_only,
                min_af=args.min_af,
                min_callrate=args.min_callrate,
                min_inbreeding_coeff_threshold=args.min_inbreeding_coeff_threshold,
                ld_r2=ld_r2,
            )
            mt.write(
                get_checkpoint_path("dense_ld_prune_qc_mt.test", mt=True)
                if test
                else qc.path,
                overwrite=overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("generate_qc_mt"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the VDS.", action="store_true"
    )
    parser.add_argument(
        "--create-v3-filtered-dense-mt",
        help=(
            "Create a dense MatrixTable from the raw gnomAD v3.1 sparse MatrixTable filtered to predetermined "
            "QC variants."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--create-v4-filtered-dense-mt",
        help="Create a dense MatrixTable from the raw gnomAD v4 VariantDataset filtered to predetermined QC variants.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-qc-mt",
        help="Create the final merged gnomAD v3 + v4 QC MatrixTable with all specified filters and LD-pruning.",
        action="store_true",
    )
    parser.add_argument(
        "--bi-allelic-only",
        help="Filter to variants that are bi-allelic.",
        action="store_true",
    )
    parser.add_argument(
        "--min-af",
        help="Minimum variant allele frequency to retain variant in QC MatrixTable.",
        default=0.0001,
        type=float,
    )
    parser.add_argument(
        "--min-callrate",
        help="Minimum variant callrate to retain variant in QC MatrixTable.",
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--min-inbreeding-coeff-threshold",
        help="Minimum site inbreeding coefficient to retain variant in QC MatrixTable.",
        default=-0.8,
        type=float,
    )
    parser.add_argument(
        "--ld-r2",
        help="LD-pruning cutoff.",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--n-partitions",
        help="Desired number of partitions for output QC MatrixTable.",
        default=1000,
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
