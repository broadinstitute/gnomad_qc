# noqa: D100

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import lcr_intervals
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.ancestry import run_pca_with_relateds
from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.file_utils import check_file_exists_raise_error
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v3.resources import release_sites
from gnomad_qc.v3.resources.annotations import get_info
from gnomad_qc.v3.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v3_mt,
    get_logging_path,
)
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    assigned_subpops,
    filtered_subpop_qc_mt,
    pca_related_samples_to_drop,
    subpop_outliers,
    subpop_qc,
)
from gnomad_qc.v3.sample_qc.sample_qc import assign_pops

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subpop_analysis")
logger.setLevel(logging.INFO)

CURATED_SUBPOPS = {
    # NOTE: 'Dai' is a known subpop label within eas but is removed as only 8 samples have defined subpop labels in this group, and it is similar to "Chinese Dai" which has more samples with defined subpop labels (92 samples) # noqa
    # NOTE: 'Han' is a known subpop label within eas but is removed as it is already encompassed by more distinct subpops, "Han Chinese" and "Southern Han Chinese" ("Han" overlaps both "Han Chinese" and "Southern Han Chinese" in PCA plots) # noqa
    # NOTE: 'Utah Residents (European Ancestry)' is a known subpop label within nfe but is removed as it is not a descriptive/accurate label and clusters near (0,0) on many PCs, which is potentially a result of missing data # noqa
    "eas": [
        "Cambodian",
        "Chinese Dai",
        "Daur",
        "Han Chinese",
        "Hezhen",
        "Japanese",
        "Kinh",
        "Lahu",
        "Miaozu",
        "Mongola",
        "Naxi",
        "Oroqen",
        "She",
        "Southern Han Chinese",
        "Tu",
        "Tujia",
        "Uygur",
        "Xibo",
        "Yakut",
        "Yizu",
    ],
    "sas": [
        "Balochi",
        "Bengali",
        "Brahui",
        "Burusho",
        "Gujarati",
        "Hazara",
        "Indian Telugu",
        "Kalash",
        "Makrani",
        "Pakistani",
        "Pathan",
        "Punjabi",
        "Sindhi",
        "Sri Lankan Tamil",
    ],
    "mid": ["Bedouin", "Druze", "Mozabite", "Palestinian"],
    "amr": [
        "Colombian",
        "Costa Rican",
        "Hawaiian",
        "Indigenous American",
        "Karitiana",
        "Maya",
        "Mexican-American",
        "Peruvian",
        "Pima",
        "Puerto Rican",
        "Rapa Nui from Easter Island",
        "Surui",
    ],
    "afr": [
        "African Caribbean",
        "African-American",
        "Bantu Kenya",
        "Bantu S Africa",
        "Biaka Pygmy",
        "Continental African",
        "Esan",
        "Gambian",
        "Luhya",
        "Mandenka",
        "Mbuti Pygmy",
        "Mende",
        "San",
        "Yoruba",
    ],
    "nfe": [
        "Adygei",
        "Basque",
        "British",
        "Dutch",
        "Estonian",
        "French",
        "German",
        "Iberian",
        "Italian",
        "Norwegian",
        "Orcadian",
        "Russian",
        "Sardinian",
        "Swedish",
        "Toscani",
    ],
}
"""
Manually curated list of which subpopulations to include within each global population. Manual curation was needed as some known subpop labels may be inaccurate, especially those that are cohort-level annotations from collaborator metadata.
For example, curation is needed to remove "Costa Rican" and "Indigenous American" subpops from the list of known subpops for the East Asian population.
"""


def compute_subpop_qc_mt(
    mt: hl.MatrixTable,
    min_popmax_af: float = 0.001,
) -> hl.MatrixTable:
    """
    Generate the subpop QC MT to be used for all subpop analyses.

    Filter MT to bi-allelic SNVs at sites with popmax allele frequency greater than `min_popmax_af`, densify, and write MT.

    :param mt: Raw MatrixTable to use for the subpop analysis
    :param min_popmax_af: Minimum population max variant allele frequency to retain variant for the subpop QC MatrixTable
    :return: MatrixTable filtered to variants for subpop analysis
    """
    release_ht = hl.read_table(release_sites().path)

    # Filter to biallelic SNVs not in low-confidence regions and with a popmax
    # above min_popmax_af
    qc_sites = release_ht.filter(
        (release_ht.popmax.AF > min_popmax_af)
        & ~release_ht.was_split
        & hl.is_snp(release_ht.alleles[0], release_ht.alleles[1])
        & hl.is_missing(lcr_intervals.ht()[release_ht.locus])
    )

    logger.info(
        "Checkpointing the QC sites HT of biallelic SNVs that are not in"
        " low-confidence regions and have a popmax above the specified minimum allele"
        f" frequency of {min_popmax_af}."
    )
    qc_sites = qc_sites.checkpoint(
        get_checkpoint_path("qc_sites"),
        overwrite=args.overwrite,
        _read_if_exists=not args.overwrite,
    )

    logger.info(
        "Converting MT to VDS..."
    )  # Convert to VDS to more efficiently densify the data
    mt = mt.select_entries(
        "END", "LA", "LGT", adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
    )

    vds = hl.vds.VariantDataset.from_merged_representation(mt)

    logger.info("Filtering to QC sites...")
    vds = hl.vds.filter_variants(vds, qc_sites)

    logger.info("Densifying data...")
    mt = hl.vds.to_dense_mt(vds)

    return mt


def filter_subpop_qc(
    mt: hl.MatrixTable,
    pop: str,
    min_af: float = 0.001,
    min_inbreeding_coeff_threshold: float = -0.25,
    min_hardy_weinberg_threshold: float = 1e-8,
    ld_r2: float = 0.1,
    n_partitions: int = None,
    block_size: int = None,
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
    :param n_partitions: Number of partitions to repartition the MT to before LD pruning
    :param block_size: If given, set the block size to this value when LD pruning
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
    mt = mt.transmute_entries(GT=mt.LGT)

    # Add sample metadata to the QC MT
    mt = mt.annotate_cols(**meta_ht[mt.col_key])

    # Filter the  QC MT to the specified pop
    pop_mt = mt.filter_cols(mt.population_inference.pop == pop)

    # Remove hard filtered samples
    pop_mt = pop_mt.filter_cols(~pop_mt.sample_filters.hard_filtered)
    pop_mt = pop_mt.filter_rows(hl.agg.any(pop_mt.GT.is_non_ref()))
    pop_mt = pop_mt.checkpoint(
        get_checkpoint_path("pop_mt", mt=True),
        overwrite=args.overwrite,
        _read_if_exists=not args.overwrite,
    )

    # Generate a QC MT for the given pop
    pop_qc_mt = get_qc_mt(
        pop_mt,
        min_af=min_af,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=min_hardy_weinberg_threshold,
        ld_r2=ld_r2,
        filter_decoy=False,
        checkpoint_path=get_checkpoint_path("intermediate_qc_mt", mt=True),
        n_partitions=n_partitions,
        block_size=block_size,
    )

    return pop_qc_mt


def drop_small_subpops(
    ht: hl.Table,
    min_additional_subpop_samples: int = None,
    unassigned_label: str = "remaining",
) -> hl.Table:
    """
    Drop small subpops from final subpop inference results.

    .. note::

        Samples within a dropped subpop will be reassigned to `unassigned_label`. Small subpops are those below `min_additional_subpop_samples` and criteria is based on the source of the known labels:
            - If the known labels came exclusively from HGDP and/or TGP(1000 Genomes Project), require newly assigned samples (n_newly_assigned) >= min_additional_subpop_samples
            - If the known labels came exclusively from a collaborator metadata source, require the sum of confirmed and newly assigned samples (n_assigned) >= min_additional_subpop_samples
            - If known labels came from both HGDP/TGP and collaborator metadata sources, require the sum of confirmed and newly assigned samples (n_assigned) minus confirmed HGDP/TGP samples (n_confirmed_hgdp_tgp) >= min_additional_subpop_samples

    :param ht: Table containing subpop inference results (under annotation name 'subpop').
    :param min_additional_subpop_samples: Minimum additional samples required to include subpop in final inference results. Default is 100.
    :param unassigned_label: Label to use for samples for which inferred subpop label will be dropped. Default is 'remaining'.
    :return: Table with final inference results in which samples belonging to small subpops have been reassigned to `unassigned_label`.
    """
    # For each training_pop, count the number of samples with known labels and the number of hgdp_or_tgp samples correctly assigned to their known label
    known_label_counts_ht = ht.group_by(ht.training_pop).aggregate(
        n_known_labels=hl.agg.count(),
        n_confirmed_hgdp_tgp=hl.agg.count_where(
            hl.is_defined(ht.training_pop)
            & hl.is_defined(ht.subpop)
            & (ht.subpop == ht.training_pop)
            & ht.hgdp_or_tgp
        ),
    )

    # For each subpop, count the number of samples assigned to that subpop overall, as well as newly assigned samples (did not have a known label)
    assigned_counts_ht = ht.group_by(ht.subpop).aggregate(
        n_assigned=hl.agg.count(),
        n_newly_assigned=hl.agg.count_where(
            ~ht.training_sample & ~ht.evaluation_sample & hl.is_defined(ht.subpop)
        ),
    )

    # Annotate whether known labels came from an hgdp_or_tgp sample, a collaborator metadata annotation, or a "mix" of both
    label_source_ht = ht.group_by(ht.training_pop).aggregate(
        known_label_source=(
            hl.case()
            .when(hl.agg.all(ht.hgdp_or_tgp), "hgdp_tgp_only")
            .when(hl.agg.all(~ht.hgdp_or_tgp), "collab_source")
            .default("mix")
        )
    )

    # Join the Table of known label counts with the Table of assigned subpop counts and label source
    counts_subpops_ht = known_label_counts_ht.join(
        assigned_counts_ht, how="outer"
    ).join(label_source_ht, how="outer")

    # Annotate which subpops to keep in the final inference results
    counts_subpops_ht = counts_subpops_ht.annotate(
        keep_subpop=(
            hl.case()
            .when(
                (counts_subpops_ht.known_label_source == "hgdp_tgp_only")
                & (counts_subpops_ht.n_newly_assigned >= min_additional_subpop_samples),
                True,
            )
            .when(
                (counts_subpops_ht.known_label_source == "collab_source")
                & (counts_subpops_ht.n_assigned >= min_additional_subpop_samples),
                True,
            )
            .when(
                (counts_subpops_ht.known_label_source == "mix")
                & (
                    counts_subpops_ht.n_assigned
                    - counts_subpops_ht.n_confirmed_hgdp_tgp
                    >= min_additional_subpop_samples
                ),
                True,
            )
            .default(False)
        )
    )
    counts_subpops_ht = counts_subpops_ht.filter(
        hl.is_defined(counts_subpops_ht.training_pop)
        & (counts_subpops_ht.training_pop != unassigned_label)
    )

    # Create a dictionary with 'training_pop' as key and whether or not to keep samples inferred as belonging to the subpop as values
    subpop_decisions = hl.dict(
        hl.tuple(
            [
                counts_subpops_ht.training_pop,
                hl.coalesce(counts_subpops_ht.keep_subpop, False),
            ]
        ).collect()
    )

    # Keep inferred labels only if the number of additional subpops exceeds or equals 'min_additional_subpop_samples'
    ht = ht.annotate(
        subpop=(
            hl.if_else(subpop_decisions.get(ht.subpop), ht.subpop, unassigned_label)
        )
    )

    ht = ht.annotate_globals(
        min_additional_subpop_samples=min_additional_subpop_samples
    )

    return ht


def main(args):  # noqa: D103
    pop = args.pop
    include_unreleasable_samples = args.include_unreleasable_samples
    high_quality = args.high_quality
    min_additional_subpop_samples = args.min_additional_subpop_samples
    unassigned_label = args.unassigned_label

    try:
        # Read in the raw mt
        mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)

        # Filter to test partitions if specified
        if args.test:
            logger.info("Filtering MT to chromosome 20")
            mt = mt.filter_rows(mt.locus.contig == "chr20")

        if args.remove_outliers:
            check_file_exists_raise_error(
                subpop_outliers(pop).path,
                error_if_not_exists=True,
                error_if_not_exists_msg=(
                    "The --remove-outliers option was used, but a Table of outlier"
                    f" samples does not exist for population {pop} at"
                    f" {subpop_outliers(pop).path}. Outliers should be manually"
                    " determined after visualizing the output of --run_subpop_pca."
                ),
            )
            outliers_ht = subpop_outliers(pop).ht()
        else:
            outliers_ht = None

        # Write out the densified MT
        if args.make_full_subpop_qc_mt:
            logger.info("Generating densified MT to use for all subpop analyses...")
            mt = compute_subpop_qc_mt(mt, args.min_popmax_af)
            mt.write(
                get_checkpoint_path("test_make_full_subpop_qc", mt=True)
                if args.test
                else subpop_qc.path,
                overwrite=args.overwrite,
            )

        if args.run_filter_subpop_qc:
            logger.info("Filtering subpop QC MT...")
            if args.test:
                mt = hl.read_matrix_table(
                    get_checkpoint_path("test_make_full_subpop_qc", mt=True)
                )
            else:
                mt = subpop_qc.mt()
            mt = filter_subpop_qc(
                mt,
                pop,
                args.min_af,
                args.min_inbreeding_coeff_threshold,
                args.min_hardy_weinberg_threshold,
                args.ld_r2,
                args.n_partitions,
                args.block_size,
            )

            mt.write(
                get_checkpoint_path("test_checkpoint_filtered_subpop_qc", mt=True)
                if args.test
                else filtered_subpop_qc_mt(pop),
                overwrite=args.overwrite,
            )

        if args.run_subpop_pca:
            # Read in the QC MT for a specified subpop and filter samples based on
            # user parameters
            mt = hl.read_matrix_table(
                get_checkpoint_path("test_checkpoint_filtered_subpop_qc", mt=True)
                if args.test
                else filtered_subpop_qc_mt(pop)
            )
            if not include_unreleasable_samples:
                mt = mt.filter_cols(
                    mt.project_meta.releasable
                    & ~mt.project_meta.exclude  # See https://github.com/broadinstitute/gnomad_meta/tree/master/v3.1#gnomad-project-metadata-annotation-definitions for further explanation of 'exclude'
                )
            if high_quality:
                mt = mt.filter_cols(mt.high_quality)

            logger.info("Generating PCs for subpops...")
            relateds_ht = pca_related_samples_to_drop.ht()
            pop_pca_evals, pop_pca_scores, pop_pca_loadings = run_pca_with_relateds(
                qc_mt=mt,
                related_samples_to_drop=relateds_ht,
                additional_samples_to_drop=outliers_ht,
                n_pcs=args.n_pcs,
            )

            pop_pca_evals_ht = hl.Table.parallelize(
                hl.literal(
                    [
                        {"PC": i + 1, "eigenvalue": x}
                        for i, x in enumerate(pop_pca_evals)
                    ],
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
            pop_pca_scores.write(
                get_checkpoint_path("test_pop_pca_scores_ht")
                if args.test
                else ancestry_pca_scores(
                    include_unreleasable_samples, high_quality, pop
                ).path,
                overwrite=args.overwrite,
            )
            pop_pca_loadings.write(
                get_checkpoint_path("test_pop_pca_loadings_ht")
                if args.test
                else ancestry_pca_loadings(
                    include_unreleasable_samples, high_quality, pop
                ).path,
                overwrite=args.overwrite,
            )

        if args.assign_subpops:
            logger.info("Assigning subpops...")
            joint_pca_ht, joint_pca_fit = assign_pops(
                min_prob=args.min_prob,
                include_unreleasable_samples=False,
                max_number_mislabeled_training_samples=args.max_number_mislabeled_training_samples,
                max_proportion_mislabeled_training_samples=args.max_proportion_mislabeled_training_samples,
                pcs=args.pcs,
                withhold_prop=args.withhold_prop,
                pop=pop,
                curated_subpops=CURATED_SUBPOPS[pop],
                additional_samples_to_drop=outliers_ht,
                high_quality=high_quality,
                missing_label=unassigned_label,
            )

            # Add metadata to subpop inference results
            meta_ht = meta.ht()

            meta_ht = meta_ht.select(
                meta_ht.project_meta.subpop_description,
                meta_ht.project_meta.research_project,
                meta_ht.project_meta.project_id,
                meta_ht.project_meta.title,
                meta_ht.project_meta.broad_external,
                meta_ht.project_meta.releasable,
                hgdp_or_tgp=meta_ht.subsets.hgdp | meta_ht.subsets.tgp,
            )

            joint_pca_ht = joint_pca_ht.annotate(**meta_ht[joint_pca_ht.key])

            if min_additional_subpop_samples:
                logger.info(
                    "Dropping small subpops (subpops with < %d samples added) ...",
                    min_additional_subpop_samples,
                )
                joint_pca_ht = drop_small_subpops(
                    joint_pca_ht, min_additional_subpop_samples, unassigned_label
                )

            # Annotate final subpop inference results, always keeping known labels from HGDP/1KG
            joint_pca_ht = joint_pca_ht.annotate(
                subpop=(
                    hl.case()
                    .when(joint_pca_ht.hgdp_or_tgp, joint_pca_ht.subpop_description)
                    .when(hl.is_defined(joint_pca_ht.subpop), joint_pca_ht.subpop)
                    .default(unassigned_label)
                )
            )

            joint_pca_ht.write(assigned_subpops(pop).path, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("subpop"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "This script generates a QC MT and PCA scores to use for subpop analyses."
        )
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="Runs a test of the code on only chr20 of the raw gnomAD v3 MT.",
        action="store_true",
    )
    parser.add_argument(
        "--make-full-subpop-qc-mt",
        help=(
            "Runs function to create dense MT to use as QC MT for all subpop analyses."
            " This uses --min-popmax-af to determine variants that need to be retained."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-popmax-af",
        help=(
            "Minimum population max variant allele frequency to retain a variant for"
            " the subpop QC MT."
        ),
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--run-filter-subpop-qc",
        help=(
            "Runs function to filter the QC MT to a certain subpop specified by --pop"
            " argument."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--pop",
        help=(
            "Population to which the subpop QC MT should be filtered when generating"
            " the PCA data."
        ),
        type=str,
    )
    parser.add_argument(
        "--run-subpop-pca",
        help=(
            "Runs function to generate PCA data for a certain population specified by"
            " --pop argument and based on variants in the subpop QC MT created using"
            " the --run-filter-subpop-qc argument."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--n-pcs",
        help="Number of PCs to compute.",
        type=int,
        default=20,
    )
    parser.add_argument(
        "--min-af",
        help=(
            "Minimum population variant allele frequency to retain variant in QC MT"
            " when generating the PCA data."
        ),
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--min-inbreeding-coeff-threshold",
        help=(
            "Minimum site inbreeding coefficient to keep a variant when generating the"
            " PCA data."
        ),
        type=float,
        default=-0.25,
    )
    parser.add_argument(
        "--min-hardy-weinberg-threshold",
        help="Minimum site HW test p-value to keep when generating the PCA data.",
        type=float,
        default=1e-8,
    )
    parser.add_argument(
        "--ld-r2",
        help="Minimum r2 to keep when LD-pruning.",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--n-partitions",
        help=(
            "Optional number of partitions to repartition the MT to before running LD"
            " pruning. Repartitioning to fewer partitions is useful after filtering out"
            " many variants to avoid errors regarding 'maximal_independent_set may run"
            " out of memory' while LD pruning."
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--block-size",
        help="Optional block size to use for LD pruning.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--remove-outliers",
        help=(
            "Whether to remove outliers when generating the PCA data. Outliers are"
            " manually determined after visualizing the PC plots."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--include-unreleasable-samples",
        help="Includes unreleasable samples for computing PCA.",
        action="store_true",
    )
    parser.add_argument(
        "--high-quality",
        help="Filter to only high-quality samples when generating the PCA data.",
        action="store_true",
    )
    parser.add_argument(
        "--assign-subpops",
        help="Runs function to assign subpops.",
        action="store_true",
    )
    parser.add_argument(
        "--min-prob",
        help="Minimum RF probability for subpop assignment.",
        type=float,
        default=0.90,
    )
    parser.add_argument(
        "--withhold-prop",
        help=(
            "Proportion of training pop samples to withhold from training, all samples"
            " will be kept if this flag is not used."
        ),
        type=float,
        default=None,
    )
    parser.add_argument(
        "--pcs",
        help="List of PCs to use in the subpop RF assignment.",
        type=int,
        nargs="+",
        default=list(range(1, 21)),
    )
    mislabel_parser = parser.add_mutually_exclusive_group(required=True)
    mislabel_parser.add_argument(
        "--max-number-mislabeled-training-samples",
        help=(
            "Maximum number of training samples that can be mislabelled. Can't be used"
            " if `max-proportion-mislabeled-training-samples` is already set."
        ),
        type=int,
        default=None,
    )
    mislabel_parser.add_argument(
        "--max-proportion-mislabeled-training-samples",
        help=(
            "Maximum proportion of training samples that can be mislabelled. Can't be"
            " used if `max-number-mislabeled-training-samples` is already set."
        ),
        type=float,
        default=None,
    )
    parser.add_argument(
        "--min-additional-subpop-samples",
        help=(
            "Minimum additional samples required to include subpop in final inference"
            " results. Default is None."
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--unassigned-label",
        help=(
            "Label for samples for that are not classified into a particular subpop."
            " Default is 'remaining'."
        ),
        type=str,
        default="remaining",
    )

    args = parser.parse_args()

    if args.slack_channel:
        from gnomad_qc.slack_creds import slack_token

        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
