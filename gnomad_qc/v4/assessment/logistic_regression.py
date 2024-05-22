"""Script to perform logistic regression on the unified MatrixTable of v4 exomess and genomes including ancestry PCs."""

import argparse
import logging

import hail as hl

from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
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


def get_test_intervals(ht: hl.Table) -> hl.Table:
    """
    Get the test lists for the exomes and genomes MatrixTables.

    :param joint_ht: Joint HT of v4 exomes and genomes.
    :return: Test Table
    """
    # Filter to chr22
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chr22")])

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


def densify_union_exomes_genomes(
    exomes_vds: hl.vds.VariantDataset,
    genomes_vds: hl.vds.VariantDataset,
) -> hl.MatrixTable:
    """
    Union exomes and genomes MatrixTables.

    :param exomes_vds: VDS of exomes
    :param genomes_vds: VDS of genomes
    :return: Unioned MatrixTable
    """
    logger.info("Densifying exomes...")
    exomes_mt = hl.vds.to_dense_mt(exomes_vds)
    exomes_mt = exomes_mt.annotate_cols(is_genome=False)
    exomes_mt = exomes_mt.select_entries("GT").select_rows().select_cols("is_genome")

    logger.info("Densifying genomes...")
    genomes_mt = hl.vds.to_dense_mt(genomes_vds)
    genomes_mt = genomes_mt.annotate_cols(is_genome=True)
    genomes_mt = genomes_mt.select_entries("GT").select_rows().select_cols("is_genome")

    logger.info("Unioning exomes and genomes...")
    mt = exomes_mt.union_cols(genomes_mt)
    return mt


def main(args):
    """Perform logistic regression on the unioned exomes and genomes MatrixTable."""
    hl.init(
        log="/logistic_regression.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )

    test = args.test
    overwrite = args.overwrite

    exomes_vds = get_gnomad_v4_vds(
        release_only=True,
        split=True,
    )
    genomes_vds = get_gnomad_v4_genomes_vds(
        release_only=True,
        split=True,
    )

    ht = release_sites(data_type="joint").ht()
    # Maybe only filter to variants with a union p-value < 10e-4 later

    if test:
        logger.info("Filtering to test intervals...")
        ht = get_test_intervals(ht)
        ht = ht.checkpoint(hl.utils.new_temp_file("test_intervals", "ht"))
        exomes_vds = hl.vds.filter_intervals(
            exomes_vds, ht, split_reference_blocks=True
        )
        genomes_vds = hl.vds.filter_intervals(
            genomes_vds, ht, split_reference_blocks=True
        )

    mt = densify_union_exomes_genomes(exomes_vds, genomes_vds)
    mt = mt.checkpoint(
        hl.utils.new_temp_file("union_exomes_genomes", "mt"),
    )

    pc_scores = ancestry_pca_scores().ht()
    mt = mt.annotate_cols(pc=pc_scores[mt.col_key].scores)

    # Perform logistic regression
    ht = hl.logistic_regression_rows(
        "firth",
        y=mt.is_genome,
        x=mt.GT.n_alt_alleles(),
        covariates=[1] + [mt.pc[i] for i in range(10)],
    )

    ht.write(
        f"/{'test' if test else ''}_logistic_regression_results.ht", overwrite=overwrite
    )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Import VCFs and write MatrixTable",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing MatrixTable",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
