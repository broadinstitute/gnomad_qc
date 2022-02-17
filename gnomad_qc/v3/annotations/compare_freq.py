"""
Compare frequencies for two gnomAD versions.

Generate a Hail Table containing frequencies for two gnomAD versions (specified using `args.versions_to_compare`),
and one of the following tests comparing the two frequencies:
    - Hail's contingency table test -- chi-squared or Fisher’s exact test of independence depending on min cell count
    - A logistic regression on the rows of a unified MatrixTable: version ~ genotype + [ancestry_pcs]. This second test
    is currently only implemented for the comparison of v3 genomes to v2 exomes and is therefore this:
    genome_vs_exome_status ~ genotype + [ancestry_pcs].

The Table is written out to the location of the more recent of the two gnomAD versions being compared.
"""
import argparse
import logging
from typing import List

import hail as hl

from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.resources.resource_utils import ResourceNotAvailable
from gnomad.resources.grch37.gnomad import (
    CURRENT_EXOME_RELEASE,
    CURRENT_GENOME_RELEASE,
    EXOME_POPS,
    GENOME_POPS,
    liftover,
    public_pca_loadings,
)
from gnomad.resources.grch37.gnomad import public_release as v2_public_release
from gnomad.resources.grch38.gnomad import POPS
from gnomad.resources.grch38.gnomad import public_release as v3_public_release
from gnomad.resources.grch38.gnomad import CURRENT_GENOME_RELEASE as V3_CURRENT_RELEASE
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.liftover import default_lift_data
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.basics import get_gnomad_data, get_gnomad_meta
from gnomad_qc.v3.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v3_mt,
    get_logging_path,
)
from gnomad_qc.v3.resources.annotations import get_freq_comparison
from gnomad_qc.v3.resources.sample_qc import (
    v2_v3_pc_project_pca_scores,
    v2_v3_pc_relate_pca_scores,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("compare_freq")
logger.setLevel(logging.INFO)

POPS_MAP = {
    "v2_exomes": {pop.lower() for pop in EXOME_POPS},
    "v2_genomes": {pop.lower() for pop in GENOME_POPS},
    "v3_genomes": {pop.lower() for pop in POPS},
}
POP_FORMAT = {
    "v2_exomes": "gnomad_{}",
    "v2_genomes": "gnomad_{}",
    "v3_genomes": "{}-adj",
}
CURRENT_RELEASE_MAP = {
    "v2_exomes": CURRENT_EXOME_RELEASE,
    "v2_genomes": CURRENT_GENOME_RELEASE,
    "v3_genomes": V3_CURRENT_RELEASE,
}


def extract_freq_info(version1: str, version2: str, pops: List[str]) -> hl.Table:
    """
    Merge frequency info for two of the gnomAD versions.

    This only keeps variants that PASS filtering in both callsets.

    :param version1: First gnomAD version to extract frequency information from (v3_genomes, v2_exomes, v2_genomes)
    :param version2: Second gnomAD version to extract frequency information from (v3_genomes, v2_exomes, v2_genomes)
    :param pops: List of populations in the order of the frequency lists, excluding the first 2 entries (adj, raw)
    :return: Table with combined and filtered frequency information from two gnomAD versions
    """
    versions = [version1, version2]
    release_resource_map = {"v2": v2_public_release, "v3": v3_public_release}

    if "v3_genomes" in versions:
        release_resource_map["v2"] = liftover

    logger.info("Reading release HTs and extracting frequency info...")
    hts = []
    for version in versions:
        v, d = version.split("_")
        ht = release_resource_map[v](d).ht()

        logger.info("Filtering %s to only variants that are PASS...", version)
        ht = ht.filter(hl.len(ht.filters) == 0)
        freq_index_dict = hl.eval(ht.freq_index_dict)

        # Keep full gnomAD callset adj [0], raw [1], and ancestry-specific adj frequencies
        freq_idx = [0, 1] + [
            freq_index_dict[POP_FORMAT[version].format(pop)] for pop in pops
        ]

        logger.info(
            "Keeping only frequencies that are needed (adj, raw, adj by pop)..."
        )
        ht = ht.select(**{f"freq_{version}": [ht.freq[i] for i in freq_idx]})
        ht = ht.select_globals(
            **{f"freq_meta_{version}": [ht.freq_meta[i] for i in freq_idx]}
        )

        hts.append(ht)

    logger.info("Performing an inner join on frequency HTs...")
    ht = hts[0].join(hts[1], how="inner")

    return ht


def perform_contingency_table_test(
    ht: hl.Table,
    version1: str,
    version2: str,
    pops: List[str],
    min_cell_count: int = 1000,
) -> hl.Table:
    """
    Perform Hail's `contingency_table_test` on the alleles counts between two gnomAD versions in `ht`.

    This is done on the 2x2 matrix of reference and alternate allele counts. The chi-squared test is used for any case
    where all cells of the 2x2 matrix are greater than `min_cell_count`. Otherwise Fisher’s exact test is used.

    Table must include two struct annotations with AN and AC counts for versions of gnomAD labeled `freq_version`.

    :param ht: Table containing frequency information
    :param version1: First gnomAD version to extract frequency information from (v3_genomes, v2_exomes, v2_genomes)
    :param version2: Second gnomAD version to extract frequency information from (v3_genomes, v2_exomes, v2_genomes)
    :param pops: List of populations in the order of the frequency lists, excluding the first 2 entries (adj, raw)
    :param min_cell_count: Minimum count in every cell to use the chi-squared test
    :return: Table with contingency table test results added
    """
    versions = [version1, version2]

    # First two entries of the frequency lists are excluded because they are adj and raw, we only need it by pop
    ht = ht.annotate(
        n_ref=hl.struct(
            **{
                version: ht[f"freq_{version}"][2:].AN - ht[f"freq_{version}"][2:].AC
                for version in versions
            }
        ),
        n_alt=hl.struct(
            **{version: ht[f"freq_{version}"][2:].AC for version in versions}
        ),
    )

    logger.info(
        "Computing chi squared and fisher exact tests on per population frequencies..."
    )
    ht = ht.transmute(
        contingency_table_test=hl.struct(
            **{
                pop: hl.contingency_table_test(
                    ht.n_alt[version1][i],
                    ht.n_ref[version1][i],
                    ht.n_alt[version2][i],
                    ht.n_ref[version2][i],
                    min_cell_count,
                )
                for i, pop in enumerate(pops)
            }
        ),
    )

    return ht


def filter_and_densify_v3_mt(filter_ht: hl.Table, test: bool = False):
    """
    Load, filter (to variants in `filter_ht`), and densify the gnomAD v3.1 MatrixTable.

    :param filter_ht: Table including variants the v3.1 MatrixTable should be filtered to
    :param test: Whether to filter the v3.1 MatrixTable to the first 5 partitions for testing
    :return: Dense gnomAD v3.1 MatrixTable containing only variants in `filter_ht`
    """
    # Loading unsplit v3 because this is needed for the conversion to VDS before the necessary densify
    logger.info("Loading v3.1 gnomAD MatrixTable...")
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=False,
        samples_meta=True,
        release_only=True,
    )

    if test:
        mt = mt._filter_partitions(range(5))

    logger.info("Converting MatrixTable to VDS...")
    mt = mt.select_entries(
        "END", "LA", "LGT", adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
    )
    vds = hl.vds.VariantDataset.from_merged_representation(mt)

    logger.info("Performing split-multi and filtering variants...")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)
    vds = hl.vds.filter_variants(vds, filter_ht)

    logger.info("Densifying data...")
    mt = hl.vds.to_dense_mt(vds)

    return mt


def project_on_exome_pop_pcs(test: bool = False) -> hl.Table:
    """
    Performs `pc_project` on gnomAD v3.1 using gnomAD v2 population PCA loadings.

    :param test: Whether to filter the v3.1 MatrixTable to the first 5 partitions for testing
    :return: Table with PC project scores of gnomAD v3.1 projected onto v2 PCs combined with the v2 PCs
    """
    # Load gnomAD metadata and population pc loadings
    v2_exome_meta_ht = get_gnomad_meta("exomes", full_meta=True)
    v2_exome_meta_ht = v2_exome_meta_ht.key_by("s")
    v2_exome_meta_ht = v2_exome_meta_ht.select(
        scores=[v2_exome_meta_ht[f"PC{i + 1}"] for i in range(0, 10)],
        version="v2_exomes",
    )
    v2_exome_loadings_ht = public_pca_loadings().ht()

    v2_exome_loadings_ht = default_lift_data(v2_exome_loadings_ht)
    v2_exome_loadings_ht = v2_exome_loadings_ht.filter(
        ~v2_exome_loadings_ht.ref_allele_mismatch
    )

    logger.info(f"Checkpointing the liftover of gnomAD v2 PCA loading sites.")
    v2_exome_loadings_ht = v2_exome_loadings_ht.checkpoint(
        get_checkpoint_path("v2_exome_loadings_liftover"), overwrite=True,
    )

    # Need to densify before running pc project
    logger.info("Loading v3.1 gnomAD MatrixTable...")
    mt = filter_and_densify_v3_mt(v2_exome_loadings_ht, test)

    logger.info(
        "Performing PC projection of gnomAD v3 samples into gnomAD v2 PC space..."
    )
    ht = hl.experimental.pc_project(
        mt.GT, v2_exome_loadings_ht.loadings, v2_exome_loadings_ht.pca_af,
    )
    ht = ht.annotate(version="v3_genomes")

    logger.info("Merging gnomAD v2 PCs with gnomAD v3 projected PCs...")
    joint_scores_ht = v2_exome_meta_ht.union(ht)

    return joint_scores_ht


def perform_logistic_regression(
    ht: hl.Table,
    version1: str,
    version2: str,
    test: bool,
    use_pc_project: bool,
    use_v3_qc_pc_scores: bool,
    read_checkpoint_if_exists: bool = False,
) -> hl.Table:
    """
    Perform Hail's `logistic_regression_rows` on the unified MatrixTable of two gnomAD versions including ancestry PCs.

    .. warning::

        This is currently implemented to only work with gnomAD v2 exomes and v3 genomes.

    The test performed is: genome_vs_exome_status ~ genotype + [ancestry_pcs]

    :param ht: Table containing frequency information
    :param version1: First gnomAD version to extract frequency information from (v3_genomes, v2_exomes)
    :param version2: Second gnomAD version to extract frequency information from (v3_genomes, v2_exomes)
    :param test: Whether to filter both MatrixTables to the first 5 partitions for testing
    :param use_pc_project: Whether to use PC project ancestry PCs
    :param use_v3_qc_pc_scores: Whether to use precomputed v2/v3 ancestry PCs (included relateds, doesn't include v3.1 samples and uses 10% of v3 ancestry PCA variants)
    :param read_checkpoint_if_exists: Whether to read the checkpointed v2/v3 unified MatrixTable if it exists
    :return: Table with logistic regression results
    """
    if version1 != "v3_genomes" and version2 != "v2_exomes":
        raise NotImplementedError(
            "Logistic regression with ancestry pca inclusion is currently only implemented for the comparison between "
            "gnomAD v3 genomes and v2 exomes."
        )

    if use_pc_project and not use_v3_qc_pc_scores:
        pc_ht = v2_v3_pc_project_pca_scores.ht()
    elif use_v3_qc_pc_scores:
        pc_ht = v2_v3_pc_relate_pca_scores.versions["3"].ht()
    else:
        ValueError(
            "One and only one of use_pc_project or use_v3_qc_pc_scores must be True!"
        )

    logger.info("Loading v3 gnomAD and v2 gnomAD exomes MatrixTables...")
    v2_mt = get_gnomad_data(data_type="exomes", split=True, release_samples=True)

    if test:
        v2_mt = v2_mt._filter_partitions(range(5))

    v3_mt = filter_and_densify_v3_mt(ht, test)

    logger.info("Adding is_genome annotation to the v3 MatrixTable...")
    v3_mt = v3_mt.annotate_cols(is_genome=True)
    v3_mt = v3_mt.select_entries("GT")

    logger.info(
        "Loading v2 exomes liftover HT for annotation onto the gnomAD v2 exome MatrixTable..."
    )
    v2_liftover_ht = liftover("exomes").ht()
    v2_liftover_ht = v2_liftover_ht.key_by("original_locus", "original_alleles")
    v2_liftover = v2_liftover_ht[v2_mt.row_key]

    logger.info(
        "Rekeying the gnomAD v2 exome MatrixTable with liftover locus and alleles..."
    )
    v2_mt = v2_mt.key_rows_by(locus=v2_liftover.locus, alleles=v2_liftover.alleles)

    logger.info(
        "Filtering gnomAD v2 MatrixTable to the variants found in the joined frequency HT and adding "
        "is_genome annotation..."
    )
    v2_mt = v2_mt.filter_rows(hl.is_defined(ht[v2_mt.row_key]))
    v2_mt = v2_mt.annotate_cols(is_genome=False)
    v2_mt = v2_mt.select_entries("GT")

    logger.info(
        "Performing a union of the gnomAD v2 exomes MatrixTable and gnomAD v3 MatrixTable columns..."
    )
    # Selecting only needed info in order to do the union on MatrixTables
    v3_mt = v3_mt.select_rows().select_cols("is_genome")
    v2_mt = v2_mt.select_rows().select_cols("is_genome")
    mt = v2_mt.union_cols(v3_mt)
    mt = mt.checkpoint(
        get_checkpoint_path("v2_exomes_v3_joint", mt=True),
        _read_if_exists=read_checkpoint_if_exists,
        overwrite=not read_checkpoint_if_exists,
    )

    logger.info(
        "Running a logistic regression of the rows: genome_vs_exome_status ~ genotype + [ancestry_pcs]..."
    )
    mt = mt.annotate_cols(pc=pc_ht[mt.col_key].scores)
    ht = hl.logistic_regression_rows(
        "firth",
        y=mt.is_genome,
        x=mt.GT.n_alt_alleles(),
        covariates=[1] + [mt.pc[i] for i in range(10)],
    )
    ht = ht.annotate_globals(
        pc_type="pc_project" if use_pc_project else "v3_qc_pc_scores"
    )

    return ht


def main(args):
    hl.init(log="/compare_freq.log")

    try:
        # Reorder so that v3_genomes is version1. Forces the output location to be in the more recent release location
        version1, version2 = [
            v
            for v in ["v3_genomes", "v2_genomes", "v2_exomes"]
            if v in args.versions_to_compare
        ]
        pops = list(POPS_MAP[version1] & POPS_MAP[version2])
        ht = extract_freq_info(version1, version2, pops)

        if args.test:
            ht = ht._filter_partitions(range(5))

        if args.contingency_table_test:
            ht = perform_contingency_table_test(
                ht, version1, version2, pops, min_cell_count=args.min_cell_count
            )

            ht.write(
                get_checkpoint_path(f"{version1}_{version2}.compare_freq.test")
                if args.test
                else get_freq_comparison(
                    CURRENT_RELEASE_MAP[version1],
                    version1.split("_")[1],
                    CURRENT_RELEASE_MAP[version2],
                    version2.split("_")[1],
                ).path,
                overwrite=args.overwrite,
            )

        if args.run_pc_project_v3_on_v2:
            ht = project_on_exome_pop_pcs(args.test)
            ht.write(
                get_checkpoint_path(f"v3_pc_project_on_v2.test")
                if args.test
                else v2_v3_pc_project_pca_scores.versions["3.1"].path,
                overwrite=args.overwrite,
            )

        if args.logistic_regression:
            logic_ht = perform_logistic_regression(
                ht,
                version1,
                version2,
                args.test,
                args.use_pc_project,
                args.use_v3_qc_pc_scores,
                args.read_checkpoint_if_exists,
            )

            ht = ht.annotate(logistic_regression=logic_ht[ht.key])

            ht.write(
                get_checkpoint_path(
                    f"{version1}_{version2}.compare_freq.logistic_regression.test"
                )
                if args.test
                else get_freq_comparison(
                    CURRENT_RELEASE_MAP[version1],
                    version1.split("_")[1],
                    CURRENT_RELEASE_MAP[version2],
                    version2.split("_")[1],
                    logistic_regression=True,
                ).path,
                overwrite=args.overwrite,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compare_freq"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite output files.", action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Filter Tables to only the first 5 partitions for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--versions-to-compare",
        help="Two gnomAD versions to compare allele frequencies.",
        nargs=2,
        choices=["v3_genomes", "v2_exomes", "v2_genomes"],
        default=["v3_genomes", "v2_exomes"],
    )
    parser.add_argument(
        "--contingency-table-test",
        help=(
            "Perform chi-squared or Fisher’s exact test of independence on the allele frequencies based "
            "on `min_cell_count`."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--min-cell-count",
        help="Minimum count in every cell to use the chi-squared test.",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--run-pc-project-v3-on-v2",
        help=(
            "Run the projection of v3.1 samples onto v2 PCs using v2 loadings. This is needed for '--logistic-regression' "
            "if '--use-pc-project' is wanted."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--logistic-regression",
        help=(
            "Perform a logistic regression of the rows: genome_vs_exome_status ~ genotype + [ancestry_pcs]. Only "
            "implemented for the comparison of v3 genomes to v2 exomes."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--use-pc-project",
        help="Use PC project scores as the ancestry PCs in the logistic regression '--logistic-regression'.",
        action="store_true",
    )
    parser.add_argument(
        "--use-v3-qc-pc-scores",
        help=(
            "Use the joint v3/v2_exome PC scores from an old run of PC relate (includes relateds and only used 10% of "
            "QC variants)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--read-checkpoint-if-exists",
        help="Whether to read the checkpointed v2/v3 unified MatrixTable if it exists. Used in '--logistic-regression'.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
