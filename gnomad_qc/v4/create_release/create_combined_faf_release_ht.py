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
from typing import Dict, List, Set

import hail as hl
from gnomad.resources.grch38.gnomad import POPS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.utils.annotations import (
    faf_expr,
    gen_anc_faf_max_expr,
    merge_freq_arrays,
    pop_max_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications
from statsmodels.stats.contingency_tables import StratifiedTable

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
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


def extract_freq_info(
    ht: hl.Table,
    prefix: str,
) -> hl.Table:
    """
    Extract frequencies and FAF for adj, raw (only for frequencies), adj by pop, adj by sex, and adj by pop/sex.

    The following annotations are filtered and renamed:
        - freq: {prefix}_freq
        - faf: {prefix}_faf
        - grpmax: {prefix}_grpmax
        - fafmax: {prefix}_fafmax

    The following global annotations are filtered and renamed:
        - freq_meta: {prefix}_freq_meta
        - faf_meta: {prefix}_faf_meta

    This only keeps variants that pass variant QC filtering.

    :param ht: Table with frequency and FAF information.
    :param prefix: Prefix to add to each of the filtered annotations.
    :return: Table with filtered frequency and FAF information.
    """
    logger.info(
        "Keeping only frequencies for adj, raw, adj by pop, adj by sex, and adj by "
        "pop/sex..."
    )
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            "freq": ht.freq,
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        ["group", "gen_anc", "sex"],
        combine_operator="or",
        exact_match=True,
    )
    faf_meta, faf = filter_arrays_by_meta(
        ht.faf_meta,
        {"faf": ht.faf},
        ["group", "gen_anc"],
        combine_operator="or",
        exact_match=True,
    )
    logger.info("Filtering to only adj frequencies for FAF...")
    faf_meta, faf = filter_arrays_by_meta(
        hl.literal(faf_meta),
        {"faf": faf["faf"]},
        {"group": ["adj"]},
        combine_operator="or",
    )

    # Select grpmax and fafmax
    grpmax_expr = ht.grpmax
    fafmax_expr = ht.gen_anc_faf_max
    if prefix == "exomes":
        # Note: The `grpmax` and `fafmax` structs in the exomes freq HT have two nested structs:
        # `gnomad` and `non_ukb`. This section selects only the `gnomad` values (values across full
        # v4 exomes release)
        grpmax_expr = grpmax_expr.gnomad
        fafmax_expr = fafmax_expr.gnomad

    # Rename filtered annotations with supplied prefix.
    ht = ht.select(
        **{
            f"{prefix}_freq": array_exprs["freq"],
            f"{prefix}_faf": faf["faf"],
            f"{prefix}_grpmax": grpmax_expr,
            f"{prefix}_fafmax": fafmax_expr,
        }
    )
    ht = ht.select_globals(
        **{
            f"{prefix}_freq_meta": freq_meta,
            f"{prefix}_freq_meta_sample_count": array_exprs["freq_meta_sample_count"],
            f"{prefix}_faf_meta": faf_meta,
        }
    )

    ht = ht.checkpoint(hl.utils.new_temp_file("extract_freq_info", "ht"))

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
    logger.info("Performing an outer join on frequency HTs...")
    ht = genomes_ht.join(exomes_ht, how="outer")

    # Merge exomes and genomes frequencies.
    freq, freq_meta, count_arrays_dict = merge_freq_arrays(
        farrays=[ht.genomes_freq, ht.exomes_freq],
        fmeta=[
            ht.index_globals().genomes_freq_meta,
            ht.index_globals().exomes_freq_meta,
        ],
        count_arrays={
            "counts": [
                ht.index_globals().genomes_freq_meta_sample_count,
                ht.index_globals().exomes_freq_meta_sample_count,
            ],
        },
    )
    freq_meta = hl.literal(freq_meta)

    # Compute FAF on the merged exomes + genomes frequencies.
    faf, faf_meta = faf_expr(
        freq,
        freq_meta,
        ht.locus,
        pops_to_exclude=faf_pops_to_exclude,
        pop_label="gen_anc",
    )
    faf_meta_by_pop = {
        m.get("gen_anc"): i for i, m in enumerate(faf_meta) if m.get("gen_anc")
    }
    faf_meta_by_pop = hl.literal(faf_meta_by_pop)
    # Compute group max (popmax) on the merged exomes + genomes frequencies.
    grpmax = pop_max_expr(
        freq, freq_meta, pops_to_exclude=faf_pops_to_exclude, pop_label="gen_anc"
    )
    grpmax = grpmax.annotate(faf95=faf[faf_meta_by_pop.get(grpmax.gen_anc)].faf95)

    # Annotate Table with all joint exomes + genomes computations.
    ht = ht.annotate(
        joint_freq=freq,
        joint_faf=faf,
        joint_fafmax=gen_anc_faf_max_expr(
            faf, hl.literal(faf_meta), pop_label="gen_anc"
        ),
        joint_grpmax=grpmax,
    )
    ht = ht.checkpoint(hl.utils.new_temp_file("combine_faf", "ht"))

    ht = ht.annotate_globals(
        joint_freq_meta=freq_meta,
        joint_freq_index_dict=make_freq_index_dict_from_meta(freq_meta),
        joint_faf_meta=faf_meta,
        joint_faf_index_dict=make_freq_index_dict_from_meta(hl.literal(faf_meta)),
        joint_freq_meta_sample_count=count_arrays_dict["counts"],
    )

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
    alternate allele counts for each population in pops. The stats used in
    `perform_contingency_table_test` can only be used on 2x2 matrices, so we perform
    that per population to get one statistic per population. The CMH test allows for
    multiple 2x2 matrices for a specific stratification, giving a single statistic
    across all populations.

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
        Get a dictionary of frequency StructExpressions by genetic ancestries.

        :param freq_expr: ArrayExpression of frequencies to combine.
        :param meta_expr: Frequency metadata for `freq_expr`.
        :return: Dictionary of frequency StructExpressions by population.
        """
        return {
            m.get("gen_anc"): freq_expr[i]
            for i, m in enumerate(hl.eval(meta_expr))
            if m.get("gen_anc") in pops
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
    overwrite: bool = False, test: bool = False
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the combined FAF resource creation pipeline.

    :param overwrite: Whether to overwrite existing resources. Default is False.
    :param test: Whether to use test resources. Default is False.
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
        output_resources={"comb_freq_ht": get_combined_frequency(test=test)},
        input_resources={
            "generate_freq.py": {"exomes_ht": get_freq(test=test)},
            "generate_freq_genomes.py": {
                "genomes_ht": get_freq(test=test, data_type="genomes")
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
        output_resources={"final_combined_faf_ht": get_combined_faf_release(test=test)},
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


def main(args):
    """Create combined FAF resource."""
    hl.init(
        log="/compute_combined_faf.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test_gene = args.test_gene
    overwrite = args.overwrite
    pops = list(set(POPS["v3"] + POPS["v4"]))
    faf_pops = [pop for pop in pops if pop not in POPS_TO_REMOVE_FOR_POPMAX]
    combine_faf_resources = get_combine_faf_resources(overwrite, test_gene)

    try:
        if args.create_combined_frequency_table:
            res = combine_faf_resources.combined_frequency
            res.check_resource_existence()
            exomes_ht = res.exomes_ht.ht()
            genomes_ht = res.genomes_ht.ht()

            if test_gene:
                # filter to PCSK9 1:55039447-55064852 for testing.
                exomes_ht = filter_gene_to_test(exomes_ht)
                genomes_ht = filter_gene_to_test(genomes_ht)

            # TODO: Need to resolve the type difference.
            genomes_ht = genomes_ht.annotate(
                freq=genomes_ht.freq.map(
                    lambda x: x.annotate(homozygote_count=hl.int32(x.homozygote_count))
                )
            )
            exomes_ht = exomes_ht.annotate(
                freq=exomes_ht.freq.map(
                    lambda x: x.annotate(homozygote_count=hl.int32(x.homozygote_count))
                )
            )
            exomes_ht = extract_freq_info(exomes_ht, "exomes")
            genomes_ht = extract_freq_info(genomes_ht, "genomes")

            ht = get_joint_freq_and_faf(genomes_ht, exomes_ht)
            ht = ht.annotate_globals(
                versions=hl.struct(
                    exomes=res.exomes_ht.default_version,
                    genomes=res.genomes_ht.default_version,
                )
            )
            ht.describe()
            ht.write(res.comb_freq_ht.path, overwrite=overwrite)

        if args.perform_contingency_table_test:
            res = combine_faf_resources.contingency_table_test
            res.check_resource_existence()
            ht = res.comb_freq_ht.ht()
            ht = ht.select(
                contingency_table_test=perform_contingency_table_test(
                    ht.genomes_freq,
                    ht.exomes_freq,
                    min_cell_count=args.min_cell_count,
                )
            )
            ht.describe()
            ht.write(res.contingency_table_ht.path, overwrite=overwrite)

        if args.perform_cochran_mantel_haenszel_test:
            res = combine_faf_resources.cmh_test
            res.check_resource_existence()
            ht = res.comb_freq_ht.ht()
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
            ht.describe()
            ht.write(res.cmh_ht.path, overwrite=overwrite)

        if args.finalize_combined_faf_release:
            res = combine_faf_resources.finalize_faf
            res.check_resource_existence()

            ht = res.comb_freq_ht.ht()
            ht = ht.annotate(
                contingency_table_test=res.contingency_table_ht.ht()[
                    ht.key
                ].contingency_table_test,
                cochran_mantel_haenszel_test=res.cmh_ht.ht()[
                    ht.key
                ].cochran_mantel_haenszel_test,
                joint_metric_data_type=hl.case()
                .when(
                    (hl.is_defined(ht.genomes_grpmax.AC))
                    & hl.is_defined(ht.exomes_grpmax.AC),
                    "both",
                )
                .when(hl.is_defined(ht.genomes_grpmax.AC), "genomes")
                .when(hl.is_defined(ht.exomes_grpmax.AC), "exomes")
                .default(hl.missing(hl.tstr)),
                joint_fafmax=ht.joint_fafmax.annotate(
                    joint_fafmax_data_type=hl.case()
                    .when(
                        (hl.is_defined(ht.genomes_fafmax.faf95_max))
                        & hl.is_defined(ht.exomes_fafmax.faf95_max),
                        "both",
                    )
                    .when(hl.is_defined(ht.genomes_fafmax.faf95_max), "genomes")
                    .when(hl.is_defined(ht.exomes_fafmax.faf95_max), "exomes")
                    .default(hl.missing(hl.tstr)),
                ),
            )
            ht.describe()
            ht.write(res.final_combined_faf_ht.path, overwrite=overwrite)

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
        "--test-gene",
        help="Filter Tables to only the PCSK9 gene for testing.",
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
        default=100,
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

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
