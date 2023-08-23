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
from gnomad.resources.grch37.gnomad import (
    CURRENT_EXOME_RELEASE,
    CURRENT_GENOME_RELEASE,
    EXOME_POPS,
    GENOME_POPS,
    liftover,
)
from gnomad.resources.grch37.gnomad import public_release as v2_public_release
from gnomad.resources.grch38.gnomad import CURRENT_GENOME_RELEASE as V3_CURRENT_RELEASE
from gnomad.resources.grch38.gnomad import POPS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.resources.grch38.gnomad import public_release as v3_public_release
from gnomad.utils.annotations import faf_expr, merge_freq_arrays, pop_max_expr
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict
from gnomad.utils.slack import slack_notifications
from statsmodels.stats.contingency_tables import StratifiedTable

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq_comparison
from gnomad_qc.v3.resources.basics import get_checkpoint_path, get_logging_path
from gnomad_qc.v4.resources.release import release_sites as v4_release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("compare_freq")
logger.setLevel(logging.INFO)

POPS_MAP = {
    "v2_exomes": {pop.lower() for pop in EXOME_POPS},
    "v2_genomes": {pop.lower() for pop in GENOME_POPS},
    "v3_genomes": {pop.lower() for pop in POPS["v3"]},
    "v4_exomes": {pop.lower() for pop in POPS["v4"]},
    "v4_genomes": {pop.lower() for pop in POPS["v4"]},
}
CURRENT_RELEASE_MAP = {
    "v2_exomes": CURRENT_EXOME_RELEASE,
    "v2_genomes": CURRENT_GENOME_RELEASE,
    "v3_genomes": V3_CURRENT_RELEASE,
    "v4_exomes": "4.0",  # TODO: Update once v4 release is created.
    "v4_genomes": "4.0",  # TODO: Update once v4 release is created.
}
POP_FORMAT = {
    "v2_exomes": "gnomad_{}",
    "v2_genomes": "gnomad_{}",
    "v3_genomes": "{}-adj",
}


def extract_freq_info(version1: str, version2: str, pops: List[str]) -> hl.Table:
    """
    Merge frequency info for two of the gnomAD versions.

    This only keeps variants that PASS filtering in both callsets.

    :param version1: First gnomAD version to extract frequency information from.
    :param version2: Second gnomAD version to extract frequency information from.
    :param pops: List of populations to include in the reduced frequency arrays for
        both gnomAD versions. The reduced frequency arrays will be in the order of
        `pops` after the first 2 entries (adj, raw).
    :return: Table with combined and filtered frequency information from two gnomAD
        versions.
    """
    versions = [version1, version2]
    versions_set = set([v.split("_")[0] for v in versions])
    release_resource_map = {
        "v2": v2_public_release,
        "v3": v3_public_release,
        "v4": v4_release_sites,  # v4_public_release,
    }

    if "v2" in versions_set and len(versions_set) > 1:
        release_resource_map["v2"] = liftover

    faf_pops = [pop for pop in pops if pop not in POPS_TO_REMOVE_FOR_POPMAX]

    logger.info("Reading release HTs and extracting frequency info...")
    hts = []
    for version in versions:
        v, d = version.split("_")
        ht = release_resource_map[v](d).ht()

        logger.info("Filtering %s to only variants that are PASS...", version)
        ht = ht.filter(hl.len(ht.filters) == 0)
        freq_meta = hl.eval(ht.freq_meta)

        if "grpmax" in ht.row:
            grpmax_label = "grpmax"
        else:
            grpmax_label = "popmax"
        grpmax_expr = ht[grpmax_label]

        if f"{grpmax_label}_index_dict" in ht.globals:
            grpmax_expr = grpmax_expr[ht[f"{grpmax_label}_index_dict"]["gnomad"]]

        if "faf_meta" in ht.globals:
            faf_meta = hl.eval(ht.faf_meta)
            faf_idx = [0, 1]
            faf_idx.extend(
                [
                    i
                    for i, m in enumerate(faf_meta)
                    if m.get("pop") in faf_pops and len(m) == 2
                ]
            )
            faf_meta = [faf_meta[i] for i in faf_idx]
        else:
            faf_index_dict = hl.eval(ht.faf_index_dict)
            pop_formatted = [POP_FORMAT[version].format(pop) for pop in faf_pops]
            faf_idx = [0, 1] + [faf_index_dict[pop] for pop in pop_formatted]
            faf_meta = [{"group": "adj"}, {"group": "raw"}]
            faf_meta.extend([{"group": "adj", "pop": pop} for pop in faf_pops])

        # Keep full gnomAD callset adj [0], raw [1], and ancestry-specific adj
        # frequencies.
        freq_idx = [0, 1]
        freq_idx.extend(
            [i for i, m in enumerate(freq_meta) if m.get("pop") in pops and len(m) == 2]
        )

        logger.info(
            "Keeping only frequencies that are needed (adj, raw, adj by pop)..."
        )
        ht = ht.select(
            **{
                f"freq_{version}": [ht.freq[i] for i in freq_idx],
                f"faf_{version}": [ht.faf[i] for i in faf_idx],
                f"grpmax_{version}": grpmax_expr,
            }
        )
        ht = ht.select_globals(
            **{
                f"freq_meta_{version}": [ht.freq_meta[i] for i in freq_idx],
                f"faf_meta_{version}": faf_meta,
            }
        )

        hts.append(ht)

    logger.info("Performing an inner join on frequency HTs...")
    ht = hts[0].join(hts[1], how="inner")

    joint_freq_expr, joint_freq_meta_expr = merge_freq_arrays(
        farrays=[ht[f"freq_{v}"] for v in versions],
        fmeta=[ht[f"freq_meta_{v}"] for v in versions],
    )
    joint_faf_expr, joint_faf_meta_expr = faf_expr(
        joint_freq_expr,
        hl.literal(joint_freq_meta_expr),
        ht.locus,
        POPS_TO_REMOVE_FOR_POPMAX,
    )
    joint_grpmax_expr = pop_max_expr(
        joint_freq_expr, hl.literal(joint_freq_meta_expr), POPS_TO_REMOVE_FOR_POPMAX
    )
    joint_faf_meta_by_pop = hl.literal(
        {
            m.get("pop"): i
            for i, m in enumerate(joint_faf_meta_expr)
            if m.get("pop") is not None
        }
    )
    joint_grpmax_expr = joint_grpmax_expr.annotate(
        faf95=joint_faf_expr[joint_faf_meta_by_pop.get(joint_grpmax_expr.pop)].faf95
    )
    ht = ht.annotate(
        joint_freq=joint_freq_expr,
        joint_faf=joint_faf_expr,
        joint_grpmax=joint_grpmax_expr,
    )
    ht = ht.annotate_globals(
        joint_freq_meta=joint_freq_meta_expr,
        joint_freq_index_dict=make_freq_index_dict(
            joint_freq_meta_expr, label_delimiter="-"
        ),
        joint_faf_meta=joint_faf_meta_expr,
        joint_faf_index_dict=make_faf_index_dict(
            joint_faf_meta_expr, label_delimiter="-"
        ),
    )

    return ht


def perform_contingency_table_test(
    freq1_expr: hl.expr.ArrayExpression,
    freq2_expr: hl.expr.ArrayExpression,
    min_cell_count: int = 1000,
) -> hl.expr.ArrayExpression:
    """
    Perform Hail's `contingency_table_test` on the alleles counts between two gnomAD versions in `ht`.

    This is done on the 2x2 matrix of reference and alternate allele counts. The
    chi-squared test is used for any case where all cells of the 2x2 matrix are greater
    than `min_cell_count`. Otherwise, Fisher’s exact test is used.

    Table must include two struct annotations with AN and AC counts for versions of
    gnomAD labeled `freq_version`.

    :param freq1_expr: Expression for the first version of gnomAD frequencies.
    :param freq2_expr: Expression for the second version of gnomAD frequencies.
    :param min_cell_count: Minimum count in every cell to use the chi-squared test
    :return: Table with contingency table test results added
    """
    logger.info("Computing chi squared and fisher exact tests on frequencies...")
    return hl.map(
        lambda x, y: hl.contingency_table_test(
            x.AC, x.AN - x.AC, y.AC, y.AN - y.AC, min_cell_count
        ),
        freq1_expr,
        freq2_expr,
    )


def perform_cmh_test(
    ht: hl.Table,
    freq1_expr: hl.expr.ArrayExpression,
    freq2_expr: hl.expr.ArrayExpression,
    freq1_meta_expr: hl.expr.ArrayExpression,
    freq2_meta_expr: hl.expr.ArrayExpression,
    pops: List[str],
) -> hl.Table:
    """
    Take the set of N populations with sufficient number of individuals in genomes
    (say > 5K or 10K) and exomes (say >10K or 20K), run a CMH test across the N 2x2
    tables (Exome/Genome allele1 and allele 2 counts) which would provide a p-value and
    an odds ratio. Then, since p-value alone might still flag common sites as
    significant based on ancestry mismatch within groups (though computing across all
    will likely alleviate this somewhat) I would guess we would want to flag those with
    ORs >2 (specifically OR<0.5 or OR>2) as inconsistent This is fast (mantelhaen.test
    in R even).
    """

    def _get_freq_by_pop(freq_expr, meta_expr):
        return {
            m.get("pop"): freq_expr[i]
            for i, m in enumerate(hl.eval(meta_expr))
            if m.get("pop") in pops
        }

    freq1_by_pop = _get_freq_by_pop(freq1_expr, freq1_meta_expr)
    freq2_by_pop = _get_freq_by_pop(freq2_expr, freq2_meta_expr)

    df = ht.select(
        an=[
            hl.bind(
                lambda x, y: [[x.AC, x.AN - x.AC], [y.AC, y.AN - y.AC]],
                freq1_by_pop[pop],
                freq2_by_pop[pop],
            )
            for pop in pops
        ]
    ).to_pandas()
    df["cmh"] = df.apply(lambda x: StratifiedTable(x.an).test_null_odds(), axis=1)
    df["pvalue"] = df.apply(lambda x: x.cmh.pvalue, axis=1)
    df["statistic"] = df.apply(lambda x: x.cmh.statistic, axis=1)
    df = df.drop(["cmh", "an"], axis=1)
    cmh_ht = hl.Table.from_pandas(df)

    cmh_ht = cmh_ht.key_by("locus", "alleles")
    cmh_ht = cmh_ht.annotate(
        cmh=hl.struct(odds_ratio=cmh_ht.statistic, p_value=cmh_ht.pvalue)
    )

    return cmh_ht[ht.key].cmh


def main(args):  # noqa: D103
    hl.init(log="/compare_freq.log")

    try:
        # Reorder so that v3_genomes is version1. Forces the output location to be
        # in the more recent release location
        version1, version2 = [
            v
            for v in [
                "v4_genomes",
                "v4_exomes",
                "v3_genomes",
                "v2_genomes",
                "v2_exomes",
            ]
            if v in args.versions_to_compare
        ]
        pops = list(POPS_MAP[version1] & POPS_MAP[version2])
        faf_pops = [pop for pop in pops if pop not in POPS_TO_REMOVE_FOR_POPMAX]
        ht = extract_freq_info(version1, version2, pops)

        if args.test:
            ht = ht._filter_partitions(range(20))
        ht = ht.checkpoint("gs://gnomad-tmp/julia/combined_faf/freq.ht", overwrite=True)

        if args.contingency_table_test:
            ht = ht.annotate(
                contingency_table_test=perform_contingency_table_test(
                    ht[f"freq_{version1}"],
                    ht[f"freq_{version2}"],
                    min_cell_count=args.min_cell_count,
                ),
                cochran_mantel_haenszel_test=perform_cmh_test(
                    ht,
                    ht[f"freq_{version1}"],
                    ht[f"freq_{version2}"],
                    ht[f"freq_meta_{version1}"],
                    ht[f"freq_meta_{version2}"],
                    pops=faf_pops,
                ),
            )
            ht = ht.annotate(
                **{
                    ann: ht[ann].rename({"pop": "grp"})
                    for ann in ht.row_value
                    if ann.startswith("grpmax")
                }
            )
            ht = ht.checkpoint(
                (
                    get_checkpoint_path(f"{version1}_{version2}.compare_freq.test")
                    if args.test
                    else get_freq_comparison(
                        CURRENT_RELEASE_MAP[version1],
                        version1.split("_")[1],
                        CURRENT_RELEASE_MAP[version2],
                        version2.split("_")[1],
                    ).path
                ),
                overwrite=args.overwrite,
            )
            ht.describe()

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compare_freq"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
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
        choices=["v3_genomes", "v2_exomes", "v2_genomes", "v4_genomes", "v4_exomes"],
        default=["v4_genomes", "v4_exomes"],
    )
    parser.add_argument(
        "--contingency-table-test",
        help=(
            "Perform chi-squared or Fisher’s exact test of independence on the allele"
            " frequencies based on `min_cell_count`."
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
            "Run the projection of v3.1 samples onto v2 PCs using v2 loadings. This is"
            " needed for '--logistic-regression' if '--use-pc-project' is wanted."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--logistic-regression",
        help=(
            "Perform a logistic regression of the rows: genome_vs_exome_status ~"
            " genotype + [ancestry_pcs]. Only implemented for the comparison of v3"
            " genomes to v2 exomes."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--use-pc-project",
        help=(
            "Use PC project scores as the ancestry PCs in the logistic regression"
            " '--logistic-regression'."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--use-v3-qc-pc-scores",
        help=(
            "Use the joint v3/v2_exome PC scores from an old run of PC relate (includes"
            " relateds and only used 10% of QC variants)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--read-checkpoint-if-exists",
        help=(
            "Whether to read the checkpointed v2/v3 unified MatrixTable if it exists."
            " Used in '--logistic-regression'."
        ),
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
