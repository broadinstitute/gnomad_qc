"""Script to perform logistic regression on the unified MatrixTable of v4 exomess and genomes including ancestry PCs."""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import annotate_adj

from gnomad_qc.v4.resources.basics import (
    get_gnomad_v4_genomes_vds,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v4.resources.sample_qc import ancestry_pca_scores

# Steps to perform logistic regression:
# 1. Load gnomAD v4 exomes and genomes VDS, and split_multiallelics
# 2. Union exomes and genomes matrix tables
# 3. Filter and densify the v3 genomes matrix table
# 4. Filter to the test lists

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_test_intervals(ht: hl.Table, test_partitions: bool = False) -> hl.Table:
    """
    Get the test lists for the exomes and genomes MatrixTables.

    :param ht: Joint Table of v4 exomes and genomes.
    :param test_partitions: Whether to filter to test partitions.
    :return: Table of test intervals
    """
    # Filter to chr22
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chr22")])

    if test_partitions:
        logger.info("Filtering to 2 partitions on chr22 for test...")
        ht = ht._filter_partitions(range(2))
    else:
        logger.info(
            "Filtering to sets of variants on chr22 that have questionable CMH"
            " p-values..."
        )
        # Filter to variants with AF >= 0.1 and p-value < 10e-4
        ht1 = ht.filter(
            (ht.joint.freq[0].AF >= 0.1)
            & (ht.freq_comparison_stats.stat_union.p_value < 10e-4)
            & (
                ht.freq_comparison_stats.stat_union.stat_test_name
                == "cochran_mantel_haenszel_test"
            )
        )
        # Filter to variants with AF < 0.1 and p-value < 10e-4
        ht2 = ht.filter(
            (ht.joint.freq[0].AF < 0.1)
            & (ht.joint.freq[0].AC > 0)
            & (ht.freq_comparison_stats.stat_union.p_value < 10e-4)
            & (
                ht.freq_comparison_stats.stat_union.stat_test_name
                == "cochran_mantel_haenszel_test"
            )
        )
        ht2 = ht2.naive_coalesce(50)
        ht2 = ht2._filter_partitions(range(5))
        # Filter to ~8000 variants with p-value >= 10e-4
        ht3 = ht.filter(
            (ht.freq_comparison_stats.stat_union.p_value >= 10e-4)
            & (
                ht.freq_comparison_stats.stat_union.stat_test_name
                == "cochran_mantel_haenszel_test"
            )
        )
        ht3 = ht3._filter_partitions(range(5))
        ht = ht1.union(ht2).union(ht3)
    # make this an interval table
    ht = ht.annotate(
        interval=hl.locus_interval(
            ht.locus.contig, ht.locus.position, ht.locus.position + 1
        )
    )
    ht = ht.key_by(ht.interval)
    ht = ht.select()

    return ht


def filter_and_densify_vds(
    vds: hl.vds.VariantDataset,
    filtered_intervals: hl.Table,
    data_type: str = "exomes",
) -> hl.MatrixTable:
    """
    Filter and densify the VDS.

    :param vds: VDS to filter and densify.
    :param filtered_intervals: Table of filtered variants.
    :param data_type: exomes or genomes.
    :return: Densified MatrixTable
    """
    logger.info(f"Filtering and densifying {data_type} VDS...")
    vds = hl.vds.filter_intervals(vds, filtered_intervals, split_reference_blocks=True)
    rmt = vds.reference_data
    vmt = vds.variant_data
    # TODO: check where did 'LA', 'LAD', 'LGT' go? Why did we use them in
    # generate_freq? Was the VDS new?
    vmt = vmt.select_entries("AD", "DP", "GT", "GQ")
    vmt = vmt.annotate_cols(
        is_genome=True if data_type == "genomes" else False,
    )
    vds = hl.vds.VariantDataset(rmt, vmt)
    vds = vds.checkpoint(
        hl.utils.new_temp_file(f"temp_{data_type}_vds_filtered", "vds")
    )
    mt = hl.vds.to_dense_mt(vds)
    mt = annotate_adj(mt)
    # TODO: What do we do about the high ab het?
    # TODO: For testing, not adjust by sex karyotype yet
    mt = mt.select_entries("GT", "adj").select_cols("is_genome").select_rows()
    return mt


def main(args):
    """Perform logistic regression on the unioned exomes and genomes MatrixTable."""
    hl.init(
        log="/logistic_regression.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )

    test = args.test
    test_partitions = args.test_partitions
    overwrite = args.overwrite
    try:
        ht = release_sites(data_type="joint").ht()
        # Maybe only filter to variants with a union p-value < 10e-4 later

        if test:
            ht = get_test_intervals(ht, test_partitions=test_partitions)
            ht = ht.checkpoint(hl.utils.new_temp_file("test_intervals", "ht"))

        exomes_vds = get_gnomad_v4_vds(
            release_only=True,
            split=True,
        )
        genomes_vds = get_gnomad_v4_genomes_vds(
            release_only=True,
            split=True,
        )
        if args.filter_densify_exomes:
            exomes_mt = filter_and_densify_vds(exomes_vds, ht, data_type="exomes")
            exomes_mt.checkpoint(
                "gs://gnomad-tmp-30day/qin/mt/filtered_exomes_dense.mt"
            )
        if args.filter_densify_genomes:
            genomes_mt = filter_and_densify_vds(genomes_vds, ht, data_type="genomes")
            genomes_mt.checkpoint(
                "gs://gnomad-tmp-30day/qin/mt/filtered_genomes_dense.mt"
            )
        if args.union_dense_mts:
            exomes_mt = hl.read_matrix_table(
                "gs://gnomad-tmp-30day/qin/mt/filtered_exomes_dense.mt"
            )
            genomes_mt = hl.read_matrix_table(
                "gs://gnomad-tmp-30day/qin/mt/filtered_genomes_dense.mt"
            )
            mt = hl.MatrixTable.union_rows(exomes_mt, genomes_mt)
            mt.checkpoint(
                "gs://gnomad-tmp-30day/qin/mt/unioned_exomes_genomes_dense.mt"
            )
        if args.logistic_regression:
            mt = hl.read_matrix_table(
                "gs://gnomad-tmp-30day/qin/mt/unioned_exomes_genomes_dense.mt"
            )
            # Annotate with ancestry PCs
            pc_scores = ancestry_pca_scores().ht()
            mt = mt.annotate_cols(pc=pc_scores[mt.col_key].scores)
            # Perform logistic regression
            # Checked that number of PCs used to assign ancestry on the joint exomes and genomes
            # was 20
            ht = hl.logistic_regression_rows(
                "firth",
                y=mt.is_genome,
                x=mt.GT.n_alt_alleles(),
                covariates=[1] + [mt.pc[i] for i in range(21)],
            )

            ht.write(
                f"/{'test' if test else ''}_logistic_regression_results.ht",
                overwrite=overwrite,
            )

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("freq_logistic_regression"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Import VCFs and write MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--test-partitions",
        help="Filter to test partitions",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--filter-densify-exomes",
        help="Filter and densify exomes VDS",
        action="store_true",
    )
    parser.add_argument(
        "--filter-densify-genomes",
        help="Filter and densify genomes VDS",
        action="store_true",
    )
    parser.add_argument(
        "--union-dense-mts",
        help="Union dense exomes and genomes MatrixTables",
        action="store_true",
    )
    parser.add_argument(
        "--logistic-regression",
        help="Perform logistic regression on unioned exomes and genomes MatrixTable",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
