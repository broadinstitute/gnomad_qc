"""
Create a joint gnomAD v4 exome and genome frequency and FAF.

Generate a Hail Table containing frequencies for exomes and genomes in gnomAD v4, a
joint frequency, a joint FAF, and the following tests comparing the two frequencies:
    - Hail's contingency table test -- chi-squared or Fisher’s exact test of
      independence depending on min cell count.
    - Cochran–Mantel–Haenszel test -- stratified test of independence for 2x2xK
      contingency tables.

"""
import argparse
import logging
from typing import Dict, List, Set, Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import POPS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.utils.annotations import faf_expr, merge_freq_arrays, pop_max_expr
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications
from statsmodels.stats.contingency_tables import StratifiedTable

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token

# TODO: change to freq_ht when the freq_ht is finalized
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v4.create_release.create_release_sites_ht_genomes import (
    replace_oth_with_remaining,
)
from gnomad_qc.v4.resources.annotations import (
    get_combined_frequency,
    get_freq,
    get_freq_comparison,
)
from gnomad_qc.v4.resources.basics import get_logging_path
from gnomad_qc.v4.resources.release import get_combined_faf_release

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("compute_combined_faf")
logger.setLevel(logging.INFO)

faf_pops_to_exclude = POPS_TO_REMOVE_FOR_POPMAX


# TODO: not sure if we should put this function in gnomad.utils.annotations
def pop_max_for_faf_expr(
    faf: hl.expr.ArrayExpression,
    faf_meta: hl.expr.ArrayExpression,
) -> hl.expr.StructExpression:
    """

    Create an expression containing the information about the population that has the highest FAF in `faf_meta`.

    This resulting struct contains the following fields:

        - faf95_max: float64
        - faf95_max_pop: str
        - faf99_max: float64
        - faf99_max_pop: str

    :param faf: ArrayExpression of Structs with fields ['faf95', 'faf99']
    :param faf_meta: ArrayExpression of meta dictionaries corresponding to faf (as returned by faf_expr)

    :return: Popmax struct for faf
    """
    popmax_faf_indices = hl.range(0, hl.len(faf_meta)).filter(
        lambda i: (hl.set(faf_meta[i].keys()) == {"group", "pop"})
        & (faf_meta[i]["group"] == "adj")
    )
    popmax_faf_indices = hl.eval(popmax_faf_indices)

    faf95 = hl.rbind(
        hl.sorted(
            hl.array(
                [
                    hl.struct(faf=faf[i].faf95, population=faf_meta[i]["pop"])
                    for i in popmax_faf_indices
                ]
            ),
            key=lambda f: (-f.faf, f.population),
        ),
        lambda fafs: hl.if_else(
            (hl.len(fafs) > 0) & (fafs[0].faf > 0),
            hl.struct(faf95_max=fafs[0].faf, faf95_max_pop=fafs[0].population),
            hl.struct(
                faf95_max=hl.missing(hl.tfloat), faf95_max_pop=hl.missing(hl.tstr)
            ),
        ),
    )

    faf99 = hl.rbind(
        hl.sorted(
            hl.array(
                [
                    hl.struct(faf=faf[i].faf99, population=faf_meta[i]["pop"])
                    for i in popmax_faf_indices
                ]
            ),
            key=lambda f: (-f.faf, f.population),
        ),
        lambda fafs: hl.if_else(
            (hl.len(fafs) > 0) & (fafs[0].faf > 0),
            hl.struct(faf99_max=fafs[0].faf, faf99_max_pop=fafs[0].population),
            hl.struct(
                faf99_max=hl.missing(hl.tfloat), faf99_max_pop=hl.missing(hl.tstr)
            ),
        ),
    )
    faf_popmax = faf95.annotate(**faf99)
    return faf_popmax


def extract_freq_info(
    ht: hl.Table,
    pops: List[str],
    faf_pops: List[str],
    prefix: str,
) -> hl.Table:
    """
    Extract frequencies and FAF for populations in `pops` and `faf_pops` respectively.

    .. note::

        The reduced frequency arrays include adj and raw in addition to the `pops`.

    The following annotations are filtered and renamed:
        - freq: {prefix}_freq
        - faf: {prefix}_faf
        - grpmax: {prefix}_grpmax

    The following global annotations are filtered and renamed:
        - freq_meta: {prefix}_freq_meta
        - faf_meta: {prefix}_faf_meta

    This only keeps variants that PASS filtering.

    :param ht: Table with frequency and FAF information.
    :param pops: List of populations to include in the reduced frequency arrays.
    :param faf_pops: List of populations to include in the reduced FAF arrays.
    :param prefix: Prefix to add to each of the filtered annotations.
    :return: Table with filtered frequency and FAF information.
    """
    # logger.info("Filtering Table to only variants that are PASS...")
    # ht = ht.filter(hl.len(ht.filters) == 0)
    # TODO: we won't have this in the final freq_ht

    def _get_pop_meta_indices(
        meta: hl.ArrayExpression, pop_list: List[str]
    ) -> Tuple[List[int], List[dict]]:
        """
        Keep full gnomAD callset adj [0], raw [1], and ancestry-specific adj frequencies.

        :param meta: Array of frequency metadata.
        :param pop_list: List of populations to keep.
        :return: Indices of populations to keep and their metadata.
        """
        meta = hl.eval(meta)
        # for adj and raw, but raw is not in faf_meta for genomes release HT,
        # so we will remove the duplicate idx with list(set(idx))
        idx = [0, 1]
        idx.extend(
            [i for i, m in enumerate(meta) if m.get("pop") in pop_list and len(m) == 2]
        )
        meta = [meta[i] for i in list(set(idx))]

        return idx, meta

    logger.info("Keeping only frequencies that are needed (adj, raw, adj by pop)...")
    # if faf and faf_meta doesn't exist, then we calculate it with freq
    # TODO: remove this once v4 exomes freq_ht is finalized
    if "faf_meta" not in ht.globals.keys():
        faf, faf_meta = faf_expr(
            ht.freq, ht.freq_meta, ht.locus, pops_to_exclude=faf_pops_to_exclude
        )
        faf_meta_by_pop = {
            m.get("pop"): i for i, m in enumerate(faf_meta) if m.get("pop")
        }
        faf_meta_by_pop = hl.literal(faf_meta_by_pop)

        # Compute group max (popmax)
        grpmax = pop_max_expr(
            ht.freq, ht.freq_meta, pops_to_exclude=faf_pops_to_exclude
        )
        grpmax = grpmax.annotate(faf95=faf[faf_meta_by_pop.get(grpmax.pop)].faf95)

        ht = ht.annotate(faf=faf, grpmax=grpmax)
        ht = ht.annotate_globals(faf_meta=faf_meta)

    # Rename popmax to grpmax because it's called 'popmax' in v3 release HT.
    # TODO: remove this once v4 genomes freq_ht is finalized
    if "popmax" in ht.row.keys():
        ht = ht.transmute(grpmax=ht.popmax)

    freq_idx, freq_meta = _get_pop_meta_indices(ht.freq_meta, pops)
    faf_idx, faf_meta = _get_pop_meta_indices(ht.faf_meta, faf_pops)

    # Compute FAF max (fafmax)
    fafmax = pop_max_for_faf_expr(ht.faf, ht.faf_meta)
    ht = ht.annotate(fafmax=fafmax)

    # Rename filtered annotations with supplied prefix.
    ht = ht.select(
        **{
            f"{prefix}_freq": [ht.freq[i] for i in freq_idx],
            f"{prefix}_faf": [ht.faf[i] for i in faf_idx],
            f"{prefix}_fafmax": ht.fafmax,
            f"{prefix}_grpmax": ht.grpmax,
        }
    )
    ht = ht.select_globals(
        **{
            f"{prefix}_freq_meta": freq_meta,
            f"{prefix}_faf_meta": faf_meta,
        }
    )

    return ht


def get_joint_freq_and_faf(
    genomes_ht: hl.Table,
    exomes_ht: hl.Table,
    faf_pops_to_exclude: Set[str] = POPS_TO_REMOVE_FOR_POPMAX,
) -> hl.Table:
    """
    Get joint genomes and exomes frequency and FAF information.

    :param genomes_ht: Table with genomes frequency and FAF information.
    :param exomes_ht: Table with exomes frequency and FAF information.
    :param faf_pops_to_exclude: Set of populations to exclude from the FAF calculation.
    :return: Table with joint genomes and exomes frequency and FAF information.
    """
    logger.info("Performing an inner join on frequency HTs...")
    ht = genomes_ht.join(exomes_ht, how="inner")

    # Merge exomes and genomes frequencies.
    freq, freq_meta = merge_freq_arrays(
        farrays=[ht.genomes_freq, ht.exomes_freq],
        fmeta=[ht.genomes_freq_meta, ht.exomes_freq_meta],
    )
    freq_meta = hl.literal(freq_meta)

    # Compute FAF on the merged exomes + genomes frequencies.
    faf, faf_meta = faf_expr(
        freq, freq_meta, ht.locus, pops_to_exclude=faf_pops_to_exclude
    )
    faf_meta_by_pop = {m.get("pop"): i for i, m in enumerate(faf_meta) if m.get("pop")}
    faf_meta_by_pop = hl.literal(faf_meta_by_pop)

    # Compute group max (popmax) on the merged exomes + genomes frequencies.
    grpmax = pop_max_expr(freq, freq_meta, pops_to_exclude=faf_pops_to_exclude)
    grpmax = grpmax.annotate(faf95=faf[faf_meta_by_pop.get(grpmax.pop)].faf95)

    # Annotate Table with all joint exomes + genomes computations.
    ht = ht.annotate(joint_freq=freq, joint_faf=faf, joint_grpmax=grpmax)

    ht = ht.annotate_globals(
        joint_freq_meta=freq_meta,
        joint_freq_index_dict=make_freq_index_dict_from_meta(
            freq_meta, label_delimiter="-"
        ),
        joint_faf_meta=faf_meta,
        joint_faf_index_dict=make_faf_index_dict(faf_meta, label_delimiter="-"),
    )

    # Compute FAF max (fafmax) on the merged exomes + genomes frequencies.
    fafmax = pop_max_for_faf_expr(ht.joint_faf, ht.joint_faf_meta)
    ht = ht.annotate(joint_fafmax=fafmax)

    return ht


def perform_contingency_table_test(
    freq1_expr: hl.expr.ArrayExpression,
    freq2_expr: hl.expr.ArrayExpression,
    min_cell_count: int = 1000,
) -> hl.expr.ArrayExpression:
    """
    Perform Hail's `contingency_table_test` on the alleles counts between two frequency expressions.

    This is done on the 2x2 matrix of reference and alternate allele counts. The
    chi-squared test is used for any case where all cells of the 2x2 matrix are greater
    than `min_cell_count`. Otherwise, Fisher’s exact test is used.

    `freq1_expr` and `freq2_expr` should be ArrayExpressions of structs with 'AN' and
    'AC' annotations.

    :param freq1_expr: First ArrayExpression of frequencies to combine.
    :param freq2_expr: Second ArrayExpression of frequencies to combine.
    :param min_cell_count: Minimum count in every cell to use the chi-squared test.
    :return: ArrayExpression for contingency table test results.
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
    Perform the Cochran–Mantel–Haenszel test on the alleles counts between two frequency expressions using population as the stratification.

    This is done by creating a list of 2x2 matrices of freq1/freq2 reference and
    alternate allele counts for each population in pops.

    `freq1_expr` and `freq2_expr` should be ArrayExpressions of structs with 'AN' and
    'AC' annotations.

    :param ht: Table with joint exomes and genomes frequency and FAF information.
    :param freq1_expr: First ArrayExpression of frequencies to combine.
    :param freq2_expr: Second ArrayExpression of frequencies to combine.
    :param freq1_meta_expr: Frequency metadata for `freq1_expr`.
    :param freq2_meta_expr: Frequency metadata for `freq2_expr`.
    :param pops: List of populations to  include in the CMH test.
    :return: ArrayExpression for Cochran–Mantel–Haenszel test results.
    """

    def _get_freq_by_pop(
        freq_expr: hl.expr.ArrayExpression, meta_expr: hl.expr.ArrayExpression
    ) -> Dict[str, hl.expr.StructExpression]:
        """
        Get a dictionary of frequency StructExpressions by population.

        :param freq_expr: ArrayExpression of frequencies to combine.
        :param meta_expr: Frequency metadata for `freq_expr`.
        :return: Dictionary of frequency StructExpressions by population.
        """
        return {
            m.get("pop"): freq_expr[i]
            for i, m in enumerate(hl.eval(meta_expr))
            if m.get("pop") in pops
        }

    freq1_by_pop = _get_freq_by_pop(freq1_expr, freq1_meta_expr)
    freq2_by_pop = _get_freq_by_pop(freq2_expr, freq2_meta_expr)

    # Create list on 2x2 matrices of reference and alternate allele counts for each
    # population in pops and export to pandas for CMH test.
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

    # Convert CMH result pandas DataFrame to a Table and restructure the annotation.
    cmh_ht = hl.Table.from_pandas(df)
    cmh_ht = cmh_ht.key_by("locus", "alleles")
    cmh_ht = cmh_ht.annotate(
        cmh=hl.struct(odds_ratio=cmh_ht.statistic, p_value=cmh_ht.pvalue)
    )

    return cmh_ht[ht.key].cmh


def get_combine_faf_resources(
    overwrite: bool, test: bool, public: bool = False
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the combined FAF resource creation pipeline.

    :param overwrite: Whether to overwrite existing resources.
    :param test: Whether to use test resources.
    :param public: Whether to use the public finalized combined FAF resource.
    :return: PipelineResourceCollection containing resources for all steps of the
        combined FAF resource creation pipeline.
    """
    # Initialize pipeline resource collection.
    combine_faf_pipeline = PipelineResourceCollection(
        pipeline_name="combine_faf",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the pipeline.
    combined_frequency = PipelineStepResourceCollection(
        "--create-combined-frequency-table",
        output_resources={"freq_ht": get_combined_frequency(test=test)},
        input_resources={
            # TODO: rename when release scripts are finalized.
            "generate_freq.py": {
                # use the unfinalized freq_ht for now
                "exomes_ht": get_freq(finalized=False)
            },
            "create_release_sites_ht_genomes.py": {
                # TODO: "genomes_freq_ht": get_freq(test=test, data_type='genomes',
                # finalized=True)},
                "genomes_ht": release_sites()
            },
        },
    )
    contingency_table_test = PipelineStepResourceCollection(
        "--perform-contingency-table-test",
        output_resources={
            "contingency_table_ht": get_freq_comparison(
                "contingency_table_test", test=test
            )
        },
        pipeline_input_steps=[combined_frequency],
    )
    cmh_test = PipelineStepResourceCollection(
        "--perform-cochran-mantel-haenszel-test",
        output_resources={"cmh_ht": get_freq_comparison("cmh_test", test=test)},
        pipeline_input_steps=[combined_frequency],
    )
    finalize_faf = PipelineStepResourceCollection(
        "--finalize-combined-faf-release",
        output_resources={
            "final_combined_faf_ht": get_combined_faf_release(public=public),
        },
        pipeline_input_steps=[combined_frequency, contingency_table_test, cmh_test],
    )

    # Add all steps to the combined FAF pipeline resource collection.
    combine_faf_pipeline.add_steps(
        {
            "combined_frequency": combined_frequency,
            "contingency_table_test": contingency_table_test,
            "cmh_test": cmh_test,
            "finalize_faf": finalize_faf,
        }
    )

    return combine_faf_pipeline


def filter_gene_to_test(ht: hl.Table) -> hl.Table:
    """
    Filter to PCSK9 1:55039447-55064852 for testing.

    :param ht: Table with frequency and FAF information.
    :return: Table with frequency and FAF information of the filtered interval of a gene
    """
    return hl.filter_intervals(
        ht,
        [hl.parse_locus_interval("chr1:55039447-55064852", reference_genome="GRCh38")],
    )


def main(args):
    """Create combined FAF resource."""
    hl.init(log="/compute_combined_faf.log")
    test = args.test
    overwrite = args.overwrite
    pops = POPS["v3"] + POPS["v4"]
    faf_pops = [pop for pop in pops if pop not in POPS_TO_REMOVE_FOR_POPMAX]
    combine_faf_resources = get_combine_faf_resources(overwrite, test, args.public)

    try:
        if args.create_combined_frequency_table:
            res = combine_faf_resources.combined_frequency
            res.check_resource_existence()
            exomes_ht = res.exomes_ht.ht()
            genomes_ht = res.genomes_ht.ht()

            # replace oth with remaining in freq_meta and freq_index_dict for genomes_ht
            genomes_ht = replace_oth_with_remaining(genomes_ht)
            if test:
                # exomes_ht = res.exomes_freq_ht.ht()._filter_partitions(range(20))
                # genomes_ht = res.genomes_freq_ht.ht()._filter_partitions(range(20))
                # filter to PCSK9 1:55039447-55064852 for testing
                exomes_ht = filter_gene_to_test(exomes_ht)
                genomes_ht = filter_gene_to_test(genomes_ht)

            exomes_ht = extract_freq_info(exomes_ht, pops, faf_pops, "exomes")
            genomes_ht = extract_freq_info(genomes_ht, pops, faf_pops, "genomes")
            ht = get_joint_freq_and_faf(genomes_ht, exomes_ht)
            ht.write(res.freq_ht.path, overwrite=overwrite)

        if args.perform_contingency_table_test:
            res = combine_faf_resources.contingency_table_test
            res.check_resource_existence()
            ht = res.freq_ht.ht()
            ht = ht.select(
                # TODO: this step is very slow, need to optimize
                contingency_table_test=perform_contingency_table_test(
                    ht.genomes_freq,
                    ht.exomes_freq,
                    min_cell_count=args.min_cell_count,
                )
            )
            ht.write(res.contingency_table_ht.path, overwrite=overwrite)

        if args.perform_cochran_mantel_haenszel_test:
            res = combine_faf_resources.cmh_test
            res.check_resource_existence()
            ht = res.freq_ht.ht()
            ht = ht.select(
                cochran_mantel_haenszel_test=perform_cmh_test(
                    ht,
                    ht.genomes_freq,
                    ht.exomes_freq,
                    ht.exomes_freq_meta,
                    ht.exomes_freq_meta,
                    pops=faf_pops,
                ),
            )
            ht.write(res.cmh_ht.path, overwrite=overwrite)

        if args.finalize_combined_faf_release:
            raise NotImplementedError(
                "Finalizing combined FAF release is not yet implemented."
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compute_combined_faf"))


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
        help="Filter Tables to only the first 20 partitions for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--create-combined-frequency-table",
        help=(
            "Create a Table with frequency information for exomes, genomes, and the "
            "joint exome + genome frequencies. Included frequencies are adj, raw, and "
            "adj for all populations found in both the exomes and genomes. The table "
            "also includes FAF computed on the joint frequencies."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--perform-contingency-table-test",
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
        "--perform-cochran-mantel-haenszel-test",
        help=(
            "Perform the Cochran–Mantel–Haenszel test, a stratified test of "
            "independence for 2x2xK contingency tables, on the allele"
            " frequencies where K is the number of populations with FAF computed."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--finalize-combined-faf-release",
        help="Finalize the combined FAF Table for release.",
        action="store_true",
    )
    parser.add_argument(
        "--public",
        help="Whether to write the finalized Table to a public release path.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
