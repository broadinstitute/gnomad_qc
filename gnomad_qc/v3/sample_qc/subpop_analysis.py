import argparse
import hail as hl
import logging

from gnomad.sample_qc.ancestry import run_pca_with_relateds

from gnomad_qc.v3.resources.basics import get_checkpoint_path
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    filtered_subpop_qc_mt,
    pca_related_samples_to_drop,
    subpop_qc,
)

from gnomad.resources.grch38.reference_data import lcr_intervals
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import densify_sites
from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad_qc.v3.resources import release_sites
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.annotations import get_info, last_END_position


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subpop_analysis")
logger.setLevel(logging.INFO)


def compute_subpop_qc_mt(
    mt: hl.MatrixTable, min_popmax_af: float = 0.001,
) -> hl.MatrixTable:
    """
    Generate the subpop QC MT to be used for all subpop analyses.
    
    Filter MT to bi-allelic SNVs at sites with popmax allele frequency greater than `min_popmax_af`, densify, and write MT.

    :param mt: Raw MatrixTable to use for the subpop analysis
    :param min_popmax_af: Minimum population max variant allele frequency to retain variant for the subpop QC MatrixTable
    :return: MatrixTable filtered to variants for subpop analysis
    """
    release_ht = hl.read_table(release_sites().path)

    # Filter to biallelic SNVs not in low-confidence regions and with a popmax above min_popmax_af
    qc_sites = release_ht.filter(
        (release_ht.popmax.AF > min_popmax_af)
        & ~release_ht.was_split
        & hl.is_snp(release_ht.alleles[0], release_ht.alleles[1])
        & hl.is_missing(lcr_intervals.ht()[release_ht.locus])
    ).key_by("locus")

    # Densify the MT and include only sites defined in the qc_sites HT
    mt = mt.select_entries(
        "END", GT=mt.LGT, adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
    )
    mt = densify_sites(mt, qc_sites, hl.read_table(last_END_position.path))

    return mt


def filter_subpop_qc(
    mt: hl.MatrixTable,
    pop: str,
    min_af: float = 0.001,
    min_inbreeding_coeff_threshold: float = -0.25,
    min_hardy_weinberg_threshold: float = 1e-8,
    ld_r2: float = 0.1,
) -> hl.MatrixTable:
    """
    Generate the QC MT per specified population.
    
    .. note::

        Hard filtered samples are removed before running `get_qc_mt`

    :param mt: The QC MT output by the 'compute_subpop_qc_mt' function
    :param pop: Population to which the QC MT should be filtered
    :param min_af: Minimum population variant allele frequency to retain variant in QC MT
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to retain variant in QC MT
    :param min_hardy_weinberg_threshold: Minimum site HW test p-value to keep
    :param ld_r2: Minimum r2 to keep when LD-pruning (set to `None` for no LD pruning)
    :return: Filtered QC MT with sample metadata for the specified population to use for subpop analysis
    """
    meta_ht = meta.ht()

    # Add info to the MT
    info_ht = get_info(split=False).ht()
    info_ht = info_ht.annotate(
        info=info_ht.info.select(
            # No need for AS_annotations since it's bi-allelic sites only
            **{x: info_ht.info[x] for x in info_ht.info if not x.startswith("AS_")}
        )
    )
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    # Add sample metadata to the QC MT
    mt = mt.annotate_cols(**meta_ht[mt.col_key])

    # Filter the  QC MT to the specified pop
    pop_mt = mt.filter_cols(mt.population_inference.pop == pop)

    # Remove hard filtered samples
    pop_mt = pop_mt.filter_cols(~pop_mt.sample_filters.hard_filtered)
    pop_mt = pop_mt.filter_rows(hl.agg.any(pop_mt.GT.is_non_ref()))

    # Generate a QC MT for the given pop
    pop_qc_mt = get_qc_mt(
        pop_mt,
        min_af=min_af,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=min_hardy_weinberg_threshold,
        ld_r2=ld_r2,
    )

    return pop_qc_mt


def annotate_subpop_meta(ht: hl.Table) -> hl.Table:
    meta_ht = meta.ht()
    ht = ht.annotate(**meta_ht[ht.key])
    ht = ht.annotate(
        subpop_description=ht.project_meta.subpop_description,
        v2_pop=ht.project_meta.v2_pop,
        v2_subpop=ht.project_meta.v2_subpop,
        pop=ht.population_inference.pop,
        project_id=ht.project_meta.project_id,
        project_pop=ht.project_meta.project_pop,
    )

    return ht


def main(args):  # noqa: D103
    pop = args.pop
    include_unreleasable_samples = args.include_unreleasable_samples
    high_quality = args.high_quality
    # Read in the raw mt
    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)

    # Filter to test partitions if specified
    if args.test:
        logger.info("Filtering to the first two partitions of the MT")
        mt = mt._filter_partitions(range(2))

    # Write out the densified MT
    if args.make_full_subpop_qc_mt:
        logger.info("Generating densified MT to use for all subpop analyses...")
        mt = compute_subpop_qc_mt(mt, args.min_popmax_af)
        mt = mt.write(
            get_checkpoint_path("test_make_full_subpop_qc", mt=True)
            if args.test
            else subpop_qc.path,
            overwrite=args.overwrite,
        )

    if args.run_subpop_pca:
        logger.info("Filtering subpop QC MT...")
        mt = subpop_qc.mt()
        mt = filter_subpop_qc(
            mt,
            pop,
            args.min_af,
            args.min_inbreeding_coeff_threshold,
            args.min_hardy_weinberg_threshold,
            args.ld_r2,
        )

        if args.checkpoint_filtered_subpop_qc:
            mt = mt.checkpoint(
                get_checkpoint_path("test_checkpoint_filtered_subpop_qc", mt=True)
                if args.test
                else filtered_subpop_qc_mt(pop),
                overwrite=args.overwrite,
            )

        if not include_unreleasable_samples:
            mt = mt.filter_cols(~mt.project_meta.releasable | mt.project_meta.exclude)
        if high_quality:
            mt = mt.filter_cols(mt.high_quality)
        if args.outlier_ht_path is not None:
            outliers = hl.read_table(args.outlier_ht_path)
            mt = mt.filter_cols(hl.is_missing(outliers[mt.col_key]))

        logger.info("Generating PCs for subpops...")
        relateds = pca_related_samples_to_drop.ht()
        pop_pca_evals, pop_pca_scores, pop_pca_loadings = run_pca_with_relateds(
            mt, relateds
        )

        pop_pca_evals_ht = hl.Table.parallelize(
            hl.literal(
                [{"PC": i + 1, "eigenvalue": x} for i, x in enumerate(pop_pca_evals)],
                "array<struct{PC: int, eigenvalue: float}>",
            )
        )
        pop_pca_evals_ht.write(
            get_checkpoint_path("test_pop_pca_evals_ht")
            if args.test
            else ancestry_pca_eigenvalues(
                include_unreleasable_samples, high_quality, pop
            ).path,
            overwrite=args.overwrite,
        )
        pop_pca_scores_ht.write(
            get_checkpoint_path("test_pop_pca_scores_ht")
            if args.test
            else ancestry_pca_scores(
                include_unreleasable_samples, high_quality, pop
            ).path,
            overwrite=args.overwrite,
        )
        pop_pca_loadings_ht.write(
            get_checkpoint_path("test_pop_pca_loadings_ht")
            if args.test
            else ancestry_pca_loadings(
                include_unreleasable_samples, high_quality, pop
            ).path,
            overwrite=args.overwrite,
        )

    # TODO: Need to use annotate_subpop_meta before subpop assignment
    if args.assign_subpops:
        raise NotImplementedError(
            "Sub-population assignment is not currently implemented."
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script generates a QC MT and PCA scores to use for subpop analyses"
    )
    parser.add_argument(
        "--slack-token", help="Slack token that allows integration with slack",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to",
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Runs a test of the code on only two partitions of the raw gnomAD v3 MT",
        action="store_true",
    )
    parser.add_argument(
        "--make-full-subpop-qc-mt",
        help="Runs function to create dense MT to use as QC MT for all subpop analyses. This uses --min-popmax-af to determine variants that need to be retained",
        action="store_true",
    )
    parser.add_argument(
        "--min-popmax-af",
        help="Minimum population max variant allele frequency to retain a variant for the subpop QC MT",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--run-subpop-pca",
        help="Runs function to generate PCA data for a certain population specified by --pop argument",
        action="store_true",
    )
    parser.add_argument(
        "--pop",
        help="Population to which the subpop QC MT should be filtered when generating the PCA data",
        type=str,
    )
    parser.add_argument(
        "--min-af",
        help="Minimum population variant allele frequency to retain variant in QC MT when generating the PCA data",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--min-inbreeding-coeff-threshold",
        help="Minimum site inbreeding coefficient to keep a variant when generating the PCA data",
        type=float,
        default=-0.25,
    )
    parser.add_argument(
        "--min-hardy-weinberg-threshold",
        help="Minimum site HW test p-value to keep when generating the PCA data",
        type=float,
        default=1e-8,
    )
    parser.add_argument(
        "--ld-r2", help="Minimum r2 to keep when LD-pruning", type=float, default=0.1,
    )
    parser.add_argument(
        "--outlier-ht-path",
        help="Path to Table keyed by column containing outliers to remove when generating the PCA data",
        type=str,
    )
    parser.add_argument(
        "--include-unreleasable-samples",
        help="Includes unreleasable samples for computing PCA",
        action="store_true",
    )
    parser.add_argument(
        "--high-quality",
        help="Filter to only high-quality samples when generating the PCA data",
        action="store_true",
    )
    parser.add_argument(
        "--run-subpop-pca",
        help="Runs function to generate PCA data for a certain population specified by --pop argument",
        action="store_true",
    )
    parser.add_argument(
        "--checkpoint-filtered-subpop-qc",
        help="Checkpoint the filtered subpop QC MT",
        action="store_true",
    )

    args = parser.parse_args()

    # Both a slack token and slack channel must be supplied to receive notifications on slack
    if args.slack_channel and args.slack_token:
        with slack_notifications(args.slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
