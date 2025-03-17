"""Script to assign global ancestry labels to gnomAD v4 genomes samples based on HGDP/TGP labels or AoU labels."""

import argparse
import logging
from typing import List

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file
from sympy.plotting.intervalmath import interval

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds
from gnomad_qc.v4.resources.sample_qc import get_joint_qc
from gnomad_qc.v5.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_union_dense_mt,
    hgdp_tgp_unrelateds_without_outliers_mt,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("ancestry_assignment")
logger.setLevel(logging.INFO)

# Function 1: get HGDP/TGP PC loadings

# Function 2: merge variants from AoU, HGDP/TGP and gnomAD v4 joint_qc_mt

# Function 3: make a dense MT with the variants (including entry fields: GT, GQ, DP,
# AD), will be (488412, 149623) for gnomAD v4 genomes

# get_predetermined_qc(version="3.1", test=False).mt() have (259482, 149625),
# fewer samples because it removed the hard-filtered samples. The update meta for v4
# genomes will 149623 samples.
# TODO: Do we need to get annotations like AF, site_callrate,
#  and site_inbreeding_coefficient?


def write_pca_results(
    pop_pca_eigenvalues: List[float],
    pop_pca_scores_ht: hl.Table,
    pop_pca_loadings_ht: hl.Table,
    data_set: str = "hgdp_tgp",
    overwrite: hl.bool = False,
    test: hl.bool = False,
):
    """
    Write out the eigenvalue hail Table, scores hail Table, and loadings hail Table returned by run_pca().

    :param pop_pca_eigenvalues: List of eigenvalues returned by run_pca.
    :param pop_pca_scores_ht: Table of scores returned by run_pca.
    :param pop_pca_loadings_ht: Table of loadings returned by run_pca.
    :param data_set: Data set used in sample QC, e.g. "aou" or "hgdp_tgp".
    :param overwrite: Whether to overwrite an existing file.
    :param test: Whether the test QC MT was used in the PCA.
    :return: None
    """
    pop_pca_eigenvalues_ht = hl.Table.parallelize(
        hl.literal(
            [{"PC": i + 1, "eigenvalue": x} for i, x in enumerate(pop_pca_eigenvalues)],
            "array<struct{PC: int, eigenvalue: float}>",
        )
    )
    pop_pca_eigenvalues_ht.write(
        ancestry_pca_eigenvalues(test, data_set=data_set).path,
        overwrite=overwrite,
    )
    pop_pca_scores_ht.write(
        ancestry_pca_scores(test, data_set=data_set).path,
        overwrite=overwrite,
    )
    pop_pca_loadings_ht.write(
        ancestry_pca_loadings(test, data_set=data_set).path,
        overwrite=overwrite,
    )


def run_hgdp_tgp_pca(test: bool, overwrite: bool, n_pcs: int = 20):
    """
    Run PCA on HGDP/TGP unrelated samples and export the results.

    ..note::

       This is getting the same results as in this notebook:
       https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb2.ipynb

    :param test: Whether to run in test mode (restricting to chr20).
    :param overwrite: Whether to overwrite existing PCA results.
    :param n_pcs: Number of PCs to compute. Default is 20.
    """
    # Load the HGDP/TGP unrelated samples matrix table
    mt = hgdp_tgp_unrelateds_without_outliers_mt.mt()

    # If test mode is enabled, restrict to chromosome 20
    if test:
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])

    # Perform PCA
    pop_eigenvalues, pop_scores_ht, pop_loadings_ht = run_pca_with_relateds(
        mt, n_pcs=n_pcs, autosomes_only=False
    )

    # Write PCA results
    write_pca_results(
        pop_eigenvalues,
        pop_scores_ht,
        pop_loadings_ht,
        data_set="hgdp_tgp",
        overwrite=overwrite,
        test=test,
    )


def get_aou_pc_loadings_ht(tsv_path: str) -> hl.Table:
    """
    Read the loadings from a TSV file and return a hail Table.

    :param tsv_path: Path to the TSV file.
    :return: Hail Table with the loadings.
    """
    aou_ht = hl.import_table(
        tsv_path,
        force_bgz=True,
        types={"loadings": hl.tarray(hl.tfloat64), "af": hl.tfloat64},
    )

    aou_ht = aou_ht.annotate(
        locus=hl.parse_locus(aou_ht.locus),
        alleles=hl.array(
            aou_ht.alleles.replace("\[", "")
            .replace("\]", "")
            .replace("'", "")
            .split(delim=", ")
        ),
    )

    aou_ht = aou_ht.key_by("locus", "alleles")

    return aou_ht


def union_variant_tables(
    hgdp_tgp_ht: hl.Table, aou_ht: hl.Table, test: bool
) -> hl.Table:
    """
    Merge the variants from AoU, HGDP/TGP and gnomAD v4 joint_qc_mt.

    :param hgdp_tgp_ht: Table with HGDP/TGP PC loadings.
    :param aou_ht: Table with AoU PC loadings.
    :param test: If True, run a test on chr20 only. Default is False.
    :return: Table with the merged list of variants.
    """
    # Load gnomAD v4 variants
    v4_ht = get_joint_qc().mt().rows().select()

    # Select only the necessary fields (removes annotations)
    tables = [hgdp_tgp_ht.select(), aou_ht.select(), v4_ht]

    # Apply filtering for test mode
    if test:
        interval = hl.parse_locus_interval("chr20")
        tables = [hl.filter_intervals(ht, [interval]) for ht in tables]

    # Perform union and enforce distinct loci
    merged_ht = hl.Table.union(*tables).key_by().order_by("locus")
    merged_ht = merged_ht.key_by("locus", "alleles").distinct()

    return merged_ht


def make_interval_table(variant_ht: hl.Table) -> hl.Table:
    """
    Create an interval table from a variant Hail Table (HT) with `locus`.

    :param variant_ht: Hail Table containing `locus`.
    :return: Hail Table with intervals covering each variant locus.
    """
    interval_ht = variant_ht.select(
        interval=hl.locus_interval(
            variant_ht.locus.contig,
            variant_ht.locus.position,
            variant_ht.locus.position,
            includes_end=True,
        )
    )
    interval_ht = interval_ht.key_by("interval")
    return interval_ht


def generate_dense_mt(variants_ht: hl.Table, test: bool = False) -> hl.MatrixTable:
    """
    Make a dense MT with the variants (including entry fields: GT, GQ, DP, AD), will be (488412, 149623) for gnomAD v4 genomes.

    :param variants_ht: Table with variants from AoU pc loadings, HGDP/TGP PC
       loadings and gnomAD v4 joint_qc_mt
    :param test: If True, run a test on chr20 only. Default is False.
    :return: MatrixTable with the variants
    """
    if test:
        variants_ht = hl.filter_intervals(
            variants_ht, [hl.parse_locus_interval("chr20")]
        )
    intervals_ht = make_interval_table(variants_ht)

    vds = get_gnomad_v4_genomes_vds(
        split=True,
        annotate_meta=True,
        remove_hard_filtered_samples=True,
        filter_variant_ht=variants_ht,
        filter_intervals=intervals_ht,
        entries_to_keep=["GT", "GQ", "DP", "AD"],
    )

    mt = hl.vds.to_dense_mt(vds)
    return mt


def main(args):
    """Assign global ancestry labels to samples."""
    hl.init(
        log="/assign_ancestry.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    test = args.test
    overwrite = args.overwrite

    if args.run_hgdp_tgp_pca:
        logger.info("Running PCA on HGDP/TGP unrelated samples...")
        run_hgdp_tgp_pca(test=test, overwrite=overwrite)

    if args.get_aou_pca_loadings:
        logger.info("Importing AoU PCA loadings...")
        tsv_path = ancestry_pca_loadings(data_set="aou").path.replace(".ht", ".tsv.gz")
        aou_ht = get_aou_pc_loadings_ht(tsv_path).repartition(500)
        aou_ht.write(ancestry_pca_loadings(data_set="aou").path, overwrite=overwrite)

    if args.generate_dense_mt:
        logger.info("Merging variants from AoU, HGDP/TGP and gnomAD v4 joint_qc_mt...")
        hgdp_tgp_ht = ancestry_pca_loadings(data_set="hgdp_tgp").ht()
        aou_ht = ancestry_pca_loadings(data_set="aou").ht()

        variants_ht = union_variant_tables(hgdp_tgp_ht, aou_ht, test).checkpoint(
            new_temp_file("merged_variants", "ht")
        )
        mt = generate_dense_mt(variants_ht, test=test)
        mt.write(get_union_dense_mt(test).path, overwrite=overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )
    parser.add_argument(
        "--test", help="Run script on test dataset.", action="store_true"
    )
    parser.add_argument(
        "--overwrite", help="Overwrite output files.", action="store_true"
    )
    parser.add_argument(
        "--run-hgdp-tgp-pca",
        help="Run PCA on HGDP/TGP unrelated samples.",
        action="store_true",
    )
    parser.add_argument(
        "--get-aou-pca-loadings",
        help="Import downloaded AoU PCA loadings to a hail Table.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-dense-mt",
        help="Make a dense MT with the variants.",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
