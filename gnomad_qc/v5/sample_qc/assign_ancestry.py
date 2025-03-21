"""Script to assign global ancestry labels to gnomAD v4 genomes samples based on HGDP/TGP labels or AoU labels."""

import argparse
import logging
import pickle
from typing import Any, List, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import (
    assign_population_pcs,
    pc_project,
    run_pca_with_relateds,
)
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, meta
from gnomad_qc.v4.resources.sample_qc import get_joint_qc
from gnomad_qc.v4.sample_qc.assign_ancestry import compute_precision_recall
from gnomad_qc.v5.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_pop_ht,
    get_pop_pr_ht,
    get_union_dense_mt,
    hgdp_tgp_unrelateds_without_outliers_mt,
    pop_rf_path,
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

    # Load and filter metadata to HGDP or TGP samples
    meta_ht = meta(data_type="genomes").ht()
    meta_ht = meta_ht.filter((meta_ht.subsets.hgdp | meta_ht.subsets.tgp))
    meta_ht = meta_ht.key_by(meta_ht.project_meta.sample_id)

    mt = mt.annotate_cols(new_s=meta_ht[mt.s].s).key_cols_by()
    mt = mt.transmute_cols(s=mt.new_s).key_cols_by("s")

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


def assign_pops_by_hgdp_tgp(
    min_prob: float,
    pcs: List[int] = list(range(1, 21)),
    missing_label: str = "remaining",
    test: bool = False,
) -> Tuple[hl.Table, Any]:
    """
    Assign global ancestry labels to gnomAD v4 genomes samples based on HGDP/TGP labels.

    ..note::

       This function will first project the gnomAD v4 genomes samples onto the
       HGDP/TGP PCA loadings and then assign global ancestry labels using a random
       forest model trained on the HGDP/TGP samples.

    :param min_prob: Minimum probability for a sample to be assigned a population.
    :param pcs: List of PCs to use for RF predictions. Default is [1, 2, ..., 20].
    :param missing_label: Label to assign to samples with missing predictions. Default is "remaining".
    :param test: Whether to run in test mode (restricting to chr20).
    :return: Tuple of the Table with population assignments and the RF model.
    """
    pop_field = "pop"

    pca_loadings_ht = ancestry_pca_loadings(data_set="hgdp_tgp").ht()
    pca_scores_ht = ancestry_pca_scores(data_set="hgdp_tgp").ht()
    meta_ht = meta(data_type="genomes").ht()

    mt = get_union_dense_mt(test=test).mt()

    # Collect sample IDs from the HGDP/TGP dataset
    hgdp_tgp_samples = set(pca_scores_ht.s.collect())

    # Filter the matrix table to keep samples in release but not in the
    # HGDP/TGP pc scores
    mt = mt.filter_cols(~hl.literal(hgdp_tgp_samples).contains(mt.s) & mt.meta.release)

    # fitler to adj
    mt = filter_to_adj(mt)

    logger.info(
        "Proecting %d genome samples with HGDP/TGP PCA loadings", mt.count_cols()
    )
    projections_ht = pc_project(mt, pca_loadings_ht)

    pca_scores_ht = pca_scores_ht.annotate(
        training_pop=meta_ht[pca_scores_ht.s].project_meta.project_pop
    )
    # Remove the 'oth' population from the training set
    pca_scores_ht = pca_scores_ht.transmute(
        training_pop=hl.if_else(
            pca_scores_ht.training_pop == "oth",
            hl.missing(hl.tstr),
            pca_scores_ht.training_pop,
        )
    )
    projections_ht = projections_ht.annotate(training_pop=hl.missing(hl.tstr))
    pca_scores_ht = pca_scores_ht.union(projections_ht).checkpoint(
        new_temp_file("pca_scores", "ht")
    )
    # TODO: Do we need to write the pca_scores_ht?

    # Run the pop RF.
    pop_ht, pops_rf_model = assign_population_pcs(
        pca_scores_ht,
        pc_cols=pcs,
        known_col="training_pop",
        output_col=pop_field,
        min_prob=min_prob,
        missing_label=missing_label,
    )

    pop_ht = pop_ht.annotate_globals(
        min_prob=min_prob,
        pcs=pcs,
    )

    return pop_ht, pops_rf_model


def get_aou_pc_loadings_ht(tsv_path: str) -> hl.Table:
    """
    Read the loadings from a TSV file and return a hail Table.

    As of 2025-03-18, it seems impossible to download a HT folder from AoU
    environment, so we read the HT and export it as a TSV and then read it back.

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

    if args.assign_pops_by_hgdp_tgp:
        pop_pcs = args.pop_pcs
        pop_pcs = list(range(1, pop_pcs[0] + 1)) if len(pop_pcs) == 1 else pop_pcs
        logger.info("Using following PCs: %s", pop_pcs)
        pop_ht, pops_rf_model = assign_pops_by_hgdp_tgp(
            args.min_pop_prob,
            pcs=pop_pcs,
            test=test,
        )

        logger.info("Writing pop ht...")
        pop_ht.write(get_pop_ht(test=test).path, overwrite=overwrite)

        with hl.hadoop_open(
            pop_rf_path(test=test),
            "wb",
        ) as out:
            pickle.dump(pops_rf_model, out)

    if args.compute_precision_recall:
        ht = compute_precision_recall(
            get_pop_ht(test=test).ht(), num_pr_points=args.number_pr_points
        )
        ht.write(get_pop_pr_ht(test=test).path, overwrite=overwrite)


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
    parser.add_argument(
        "--assign-pops-by-hgdp-tgp",
        help="Assign global ancestry labels to samples using HGDP/TGP labels.",
        action="store_true",
    )
    parser.add_argument(
        "--min-pop-prob",
        help="Minimum probability for a sample to be assigned a population.",
        type=float,
        default=0.75,
    )
    parser.add_argument(
        "--pop-pcs",
        help=(
            "List of PCs to use for ancestry assignment. The values provided should be"
            " 1-based. If a single integer is passed, the script assumes this"
            " represents the total PCs to use e.g. --pop-pcs=6 will use PCs"
            " 1,2,3,4,5,and 6. Defaults to 20 PCs."
        ),
        default=[20],
        type=int,
        nargs="+",
    )
    parser.add_argument(
        "--number-pr-points",
        help=(
            "Number of min prob cutoffs to compute PR metrics for. e.g. 100 will "
            "compute PR metrics for min prob of 0 to 1 in increments of 0.01. Default "
            "is 100."
        ),
        default=100,
        type=int,
    )
    parser.add_argument(
        "--compute-precision-recall",
        help=(
            "Compute precision and recall for the RF model using evaluation samples. "
            "This is computed for all evaluation samples as well as per population."
        ),
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
