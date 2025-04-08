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
from gnomad_qc.v3.resources.release import hgdp_tgp_subset
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, meta
from gnomad_qc.v4.resources.sample_qc import (
    hgdp_tgp_pop_outliers,  # this list is the same as gs://gcp-public-data--gnomad/release/3.1/secondary_analyses/hgdp_1kg_v2/pca/pca_outliers.txt
)
from gnomad_qc.v4.resources.sample_qc import get_joint_qc
from gnomad_qc.v4.sample_qc.assign_ancestry import compute_precision_recall
from gnomad_qc.v5.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_dense_mt_for_ancestry_inference,
    get_pop_ht,
    get_pop_pr_ht,
    hgdp_tgp_unrelateds_without_outliers_mt,
    pop_rf_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("ancestry_assignment")
logger.setLevel(logging.INFO)


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


def run_hgdp_tgp_pca(
    test: bool,
    overwrite: bool,
    include_relateds: bool = False,
    data_set: str = "hgdp_tgp",
    n_pcs: int = 20,
):
    """
    Run PCA on HGDP/TGP samples using high-quality variants identified by AoU.

    By default, uses unrelated samples. Set `include_relateds=True` to include relateds
    (minus outliers, hard-filtered, and 'oth' ancestry).

    .. note::
       The list of unrelated samples without outliers is based on:
       https://nbviewer.org/github/atgu/hgdp_tgp/blob/master/tutorials/nb2.ipynb

    :param test: If True, restrict to chr20 for testing purposes.
    :param overwrite: If True, overwrite existing PCA output.
    :param include_relateds: If True, include related samples.
    :param data_set: Data set used to run PCA, e.g. "hgdp_tgp" or 'hgdp_tgp_hq'.
    :param n_pcs: Number of principal components to compute. Default is 20.
    """
    mt = hgdp_tgp_subset(dense=True).mt()

    # Some HGDP/TGP sample IDs had a prefix 'v3.1::' in gnomAD releases but removed
    # in the subset release, we need to use the ones with the prefix to match gnomAD
    # and AoU genomes
    meta_ht = meta(data_type="genomes").ht()
    meta_ht = meta_ht.filter((meta_ht.subsets.hgdp | meta_ht.subsets.tgp))
    meta_ht = meta_ht.key_by(meta_ht.project_meta.sample_id)
    mt = mt.annotate_cols(new_s=meta_ht[mt.s].s).key_cols_by()
    mt = mt.transmute_cols(s=mt.new_s).key_cols_by("s")

    mt = mt.annotate_cols(
        hard_filtered=meta_ht[mt.s].sample_filters.hard_filtered,
        project_pop=meta_ht[mt.s].project_meta.project_pop,
    )

    if include_relateds:
        # Get the list of HGDP/TGP samples that are outliers
        hgdp_tgp_outliers = hgdp_tgp_pop_outliers.ht().s.collect()
        mt = mt.filter_cols(
            ~hl.literal(hgdp_tgp_outliers).contains(mt.s)
            & ~mt.hard_filtered
            & (mt.project_pop != "oth")
        )
        logger.info(
            "Filtering to %d high-quality samples (relateds included) that are not in 'oth' ancestry",
            mt.count_cols(),
        )
    else:
        # Get the list of HGDP/TGP unrelated samples without outliers
        unrelated_mt = hgdp_tgp_unrelateds_without_outliers_mt.mt()
        unrelated_samples_ht = unrelated_mt.cols()

        # Filter to unrelated samples
        mt = mt.filter_cols(hl.is_defined(unrelated_samples_ht[mt.col_key]))
        logger.info("Filtering to %d unrelated samples", mt.count_cols())

    # Get the list of high quality variants that AoU identified and used in their PCA
    aou_loadings_ht = ancestry_pca_loadings(data_set="aou").ht()

    # Filter to AoU/HGDP/TGP high quality variants
    mt = mt.semi_join_rows(aou_loadings_ht.select())

    # If test mode is enabled, restrict to chromosome 20
    if test:
        logger.info("Filtering to chr20 for testing purposes...")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])

    mt = mt.checkpoint(new_temp_file("hgdp_tgp_pca", "mt"))

    pop_eigenvalues, pop_scores_ht, pop_loadings_ht = run_pca_with_relateds(
        mt, n_pcs=n_pcs, autosomes_only=False
    )

    # Write PCA results
    write_pca_results(
        pop_eigenvalues,
        pop_scores_ht,
        pop_loadings_ht,
        data_set=data_set,
        overwrite=overwrite,
        test=test,
    )


def project_to_hgdp_tgp_pcs(
    mt: hl.MatrixTable, is_gnomad: bool = True, test: bool = False
) -> hl.Table:
    """
    Project gnomAD v4 genomes or AoU genomes onto HGDP/TGP PCA loadings.

    :param mt: MatrixTable containing the samples to project.
    :param is_gnomad: Whether the input MT is gnomAD v4 genomes. Default is True.
    :param test: If True, run a test on chr20 only. Default is False.
    :return: Table with the projections.
    """
    pca_loadings_ht = ancestry_pca_loadings(data_set="hgdp_tgp").ht()
    pca_scores_ht = ancestry_pca_scores(data_set="hgdp_tgp").ht()

    # Collect sample IDs from the HGDP/TGP dataset
    hgdp_tgp_samples = set(pca_scores_ht.s.collect())

    if test:
        logger.info("Filtering to chr20 for testing purposes...")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])

    # Filter out samples present in the HGDP/TGP dataset
    mt = mt.filter_cols(~hl.literal(hgdp_tgp_samples).contains(mt.s))
    if is_gnomad:
        # mt = mt.filter_cols(mt.meta.release)
        mt = filter_to_adj(mt)

    logger.info(
        "Projecting %d genome samples with HGDP/TGP PCA loadings", mt.count_cols()
    )
    projections_ht = pc_project(mt, pca_loadings_ht)
    pca_scores_ht = pca_scores_ht.union(projections_ht)

    return pca_scores_ht


def assign_pops_by_hgdp_tgp(
    pca_scores_ht: hl.Table,
    min_prob: float,
    pcs: List[int] = list(range(1, 21)),
    missing_label: str = "remaining",
) -> Tuple[hl.Table, Any]:
    """
    Assign global ancestry labels to gnomAD v4 genomes samples based on HGDP/TGP labels.

    ..note::

       This function will first project the gnomAD v4 genomes samples onto the
       HGDP/TGP PCA loadings and then assign global ancestry labels using a random
       forest model trained on the HGDP/TGP samples.

    :param mt: MatrixTable containing the samples to assign populations to.
    :param pca_scores_ht: Table with the PCA scores for the samples.
    :param min_prob: Minimum probability for a sample to be assigned a population.
    :param pcs: List of PCs to use for RF predictions. Default is [1, 2, ..., 20].
    :param missing_label: Label to assign to samples with missing predictions. Default is "remaining".
    :param adj_filters: Whether to filter to adj variants. Default is True.
    :param data_set: Data set used in sample QC, e.g. "hgdp_tgp".
    :return: Tuple of the Table with population assignments and the RF model.
    """
    pop_field = "pop"

    # Get the list of HGDP/TGP samples that are outliers
    hgdp_tgp_outliers = hgdp_tgp_pop_outliers.ht().s.collect()

    meta_ht = meta(data_type="genomes").ht()
    meta_ht = meta_ht.annotate(
        training_pop=hl.if_else(
            (meta_ht.subsets.hgdp | meta_ht.subsets.tgp)
            & ~hl.literal(hgdp_tgp_outliers).contains(meta_ht.s)
            & ~meta_ht.sample_filters.hard_filtered
            & (meta_ht.project_meta.project_pop != "oth"),
            meta_ht.project_meta.project_pop,
            hl.missing(hl.tstr),
        )
    )

    pca_scores_ht = pca_scores_ht.annotate(
        training_pop=meta_ht[pca_scores_ht.key].training_pop
    ).checkpoint(new_temp_file("pca_scores", "ht"))

    count_training_labels = pca_scores_ht.aggregate(
        hl.agg.count_where(~hl.is_missing(pca_scores_ht.training_pop))
    )
    logger.info(
        "%s HGDP/TGP samples are used to train the Random Forest Model",
        count_training_labels,
    )

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

    variants_ht = variants_ht.repartition(115376, shuffle=True).checkpoint(
        new_temp_file("variants", "ht")
    )

    vds = get_gnomad_v4_genomes_vds(
        split=True,
        annotate_meta=True,
        remove_hard_filtered_samples=True,
        filter_variant_ht=variants_ht,
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
    data_set = args.data_set

    if args.run_hgdp_tgp_pca:
        logger.info("Running PCA on HGDP/TGP unrelated samples...")
        run_hgdp_tgp_pca(
            test=test,
            overwrite=overwrite,
            include_relateds=args.include_relateds,
            data_set=data_set,
            n_pcs=args.n_pcs,
        )

    if args.get_aou_pca_loadings:
        logger.info("Importing AoU PCA loadings...")
        tsv_path = ancestry_pca_loadings(data_set="aou").path.replace(".ht", ".tsv.gz")
        aou_ht = get_aou_pc_loadings_ht(tsv_path).repartition(500)
        aou_ht.write(ancestry_pca_loadings(data_set="aou").path, overwrite=overwrite)

    if args.generate_dense_mt:
        logger.info(
            "Densifying gnomAD v4 genomes VDS using AoU/HGDP/TGP high-quality sites..."
        )

        aou_ht = ancestry_pca_loadings(data_set="aou").ht()
        mt = generate_dense_mt(aou_ht, test=test)

        mt.naive_coalesce(1000).write(
            get_dense_mt_for_ancestry_inference(test).path, overwrite=overwrite
        )

    if args.assign_pops_by_hgdp_tgp:
        pop_pcs = args.pop_pcs
        pop_pcs = list(range(1, pop_pcs[0] + 1)) if len(pop_pcs) == 1 else pop_pcs
        logger.info("Using following PCs: %s", pop_pcs)
        mt = get_dense_mt_for_ancestry_inference().mt()
        pca_scores_ht = project_to_hgdp_tgp_pcs(mt).checkpoint(
            ancestry_pca_scores(test=test, data_set=data_set).path,
            overwrite=overwrite,
        )
        pop_ht, pops_rf_model = assign_pops_by_hgdp_tgp(
            pca_scores_ht=pca_scores_ht,
            min_prob=args.min_pop_prob,
            pcs=pop_pcs,
        )

        logger.info("Writing pop ht...")
        pop_ht.write(get_pop_ht(test=test, data_set=data_set).path, overwrite=overwrite)

        with hl.hadoop_open(
            pop_rf_path(test=test, data_set=data_set),
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
        "--include-relateds",
        help="Include related samples in the PCA.",
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
        "--data-set",
        help="Data set used, e.g. 'hgdp_tgp' or 'aou'.",
        default="hgdp_tgp",
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
