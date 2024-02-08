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
import pandas as pd
from gnomad.resources.grch38.gnomad import POPS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.utils.annotations import (
    faf_expr,
    gen_anc_faf_max_expr,
    merge_freq_arrays,
    pop_max_expr,
    set_female_y_metrics_to_na_expr,
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
from gnomad_qc.v4.resources.variant_qc import final_filter

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
    ht: hl.Table, prefix: str, apply_release_filters: bool = False
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
    :param apply_release_filters: Whether to apply the final release filters to the
        Table. Default is False.
    :return: Table with filtered frequency and FAF information.
    """
    if apply_release_filters:
        # Filter out chrM, AS_lowqual sites (these sites are dropped in the
        # final_filters HT so will not have information in `filters`) and AC_raw == 0.
        ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
        ht = ht.filter(hl.is_defined(ht.filters) & (ht.freq[1].AC > 0))

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
        # Note: The `grpmax` and `fafmax` structs in the exomes freq HT have two nested
        # structs: `gnomad` and `non_ukb`. This section selects only the `gnomad`
        # values (values across full v4 exomes release)
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
    faf_pops_to_exclude: Set[str] = POPS_TO_REMOVE_FOR_POPMAX["v4"],
) -> hl.Table:
    """
    Get joint genomes and exomes frequency and FAF information.

    :param genomes_ht: Table with genomes frequency and FAF information.
    :param exomes_ht: Table with exomes frequency and FAF information.
    :param faf_pops_to_exclude: Set of genetic ancestry groups to exclude from the FAF
        calculation.
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

    ht = ht.annotate(joint_freq=freq)
    ht = ht.annotate_globals(
        joint_freq_meta=freq_meta,
        joint_freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta)),
    )
    # NOTE: This checkpoint prevents a Class Too Large error in hail 0.2.122
    ht = ht.checkpoint(hl.utils.new_temp_file("combine_faf", "ht"))

    logger.info("Setting Y metrics to NA for XX groups...")
    freq = set_female_y_metrics_to_na_expr(
        ht,
        freq_expr=ht.joint_freq,
        freq_meta_expr=ht.joint_freq_meta,
        freq_index_dict_expr=ht.joint_freq_index_dict,
    )
    ht = ht.annotate(joint_freq=freq)

    # Compute FAF on the merged exomes + genomes frequencies.
    faf, faf_meta = faf_expr(
        ht.joint_freq,
        ht.joint_freq_meta,
        ht.locus,
        pops_to_exclude=faf_pops_to_exclude,
        pop_label="gen_anc",
    )
    faf_meta_by_pop = {
        m.get("gen_anc"): i
        for i, m in enumerate(faf_meta)
        if m.get("gen_anc") and len(m) == 2
    }
    faf_meta_by_pop = hl.literal(faf_meta_by_pop)
    # Compute group max (popmax) on the merged exomes + genomes frequencies.
    grpmax = pop_max_expr(
        ht.joint_freq,
        ht.joint_freq_meta,
        pops_to_exclude=faf_pops_to_exclude,
        pop_label="gen_anc",
    )
    grpmax = grpmax.annotate(faf95=faf[faf_meta_by_pop.get(grpmax.gen_anc)].faf95)

    # Annotate Table with all joint exomes + genomes computations.
    ht = ht.annotate(
        joint_faf=faf,
        joint_fafmax=gen_anc_faf_max_expr(
            faf, hl.literal(faf_meta), pop_label="gen_anc"
        ),
        joint_grpmax=grpmax,
    )

    ht = ht.annotate_globals(
        joint_faf_meta=faf_meta,
        joint_faf_index_dict=make_freq_index_dict_from_meta(hl.literal(faf_meta)),
        joint_freq_meta_sample_count=count_arrays_dict["counts"],
    )

    return ht


def perform_contingency_table_test(
    freq1_expr: hl.expr.ArrayExpression,
    freq2_expr: hl.expr.ArrayExpression,
    freq1_meta_expr: hl.expr.ArrayExpression,
    freq2_meta_expr: hl.expr.ArrayExpression,
    joint_meta_expr: hl.expr.ArrayExpression,
    min_cell_count: int = 5,
) -> hl.expr.ArrayExpression:
    """
    Perform Hail's `contingency_table_test` on the alleles counts between two frequency expressions.

    This is done on the 2x2 matrix of reference and alternate allele counts. The
    chi-squared test is used for any case where all cells of the 2x2 matrix are greater
    than `min_cell_count`. Otherwise, Fisher’s exact test is used.

    `freq1_expr` and `freq2_expr` should be ArrayExpressions of structs with 'AN' and
    'AC' annotations.

    .. note::

        The order of the output array expression will be the same as `joint_meta_expr`
        and any frequency group with missing or zero AC in both `freq1_expr` and
        `freq2_expr` (based on `freq1_meta_expr` and `freq2_meta_expr`) will be set to
        missing. Any frequency group in `freq1_meta_expr` or `freq2_meta_expr` that is
        not in `joint_meta_expr` will be excluded from tests.

    :param freq1_expr: First ArrayExpression of frequencies to combine.
    :param freq2_expr: Second ArrayExpression of frequencies to combine.
    :param freq1_meta_expr: Frequency metadata for `freq1_expr`.
    :param freq2_meta_expr: Frequency metadata for `freq2_expr`.
    :param joint_meta_expr: Joint frequency metadata, only used for ordering the output
        array expression.
    :param min_cell_count: Minimum count in every cell to use the chi-squared test.
        Default is 5.
    :return: ArrayExpression for contingency table test results.
    """
    # Using the joint_meta_expr to get the indexes of the two frequency expressions.
    freq_meta_idx = joint_meta_expr.map(
        lambda x: (freq1_meta_expr.index(x), freq2_meta_expr.index(x))
    )

    logger.info("Computing chi squared and fisher exact tests on frequencies...")
    # If both frequency structs are defined and at least one AC is greater than 0,
    # compute the chi-squared test otherwise return missing.
    return freq_meta_idx.map(
        lambda x: hl.or_missing(
            hl.is_defined(x[0]) & hl.is_defined(x[1]),
            hl.bind(
                lambda f1, f2: hl.or_missing(
                    (f1.AC > 0) | (f2.AC > 0),
                    hl.contingency_table_test(
                        f1.AC, f1.AN - f1.AC, f2.AC, f2.AN - f2.AC, min_cell_count
                    ),
                ),
                freq1_expr[x[0]],
                freq2_expr[x[1]],
            ),
        )
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
    Perform the Cochran–Mantel–Haenszel test on the alleles counts between two frequency expressions using genetic ancestry group as the stratification.

    This is done by creating a list of 2x2 matrices of freq1/freq2 reference and
    alternate allele counts for each genetic ancestry group in pops. The stats used in
    `perform_contingency_table_test` can only be used on 2x2 matrices, so we perform
    that per genetic ancestry group to get one statistic per genetic ancestry group.
    The CMH test allows for multiple 2x2 matrices for a specific stratification, giving
    a single statistic across all genetic ancestry groups.

    `freq1_expr` and `freq2_expr` should be ArrayExpressions of structs with 'AN' and
    'AC' annotations.

    .. note::

        Any genetic ancestry group with zero AC in both `freq1_expr` and `freq2_expr`
        will be excluded from the test.

    :param ht: Table with joint exomes and genomes frequency and FAF information.
    :param freq1_expr: First ArrayExpression of frequencies to combine.
    :param freq2_expr: Second ArrayExpression of frequencies to combine.
    :param freq1_meta_expr: Frequency metadata for `freq1_expr`.
    :param freq2_meta_expr: Frequency metadata for `freq2_expr`.
    :param pops: List of genetic ancestry groups to include in the CMH test.
    :return: ArrayExpression for Cochran–Mantel–Haenszel test results.
    """

    def _get_freq_by_pop(
        freq_expr: hl.expr.ArrayExpression, meta_expr: hl.expr.ArrayExpression
    ) -> Dict[str, hl.expr.StructExpression]:
        """
        Get a dictionary of frequency StructExpressions by genetic ancestries.

        :param freq_expr: ArrayExpression of frequencies to combine.
        :param meta_expr: Frequency metadata for `freq_expr`.
        :return: Dictionary of frequency StructExpressions by genetic ancestry group.
        """
        return {
            m.get("gen_anc"): freq_expr[i]
            for i, m in enumerate(hl.eval(meta_expr))
            if m.get("gen_anc") in pops
        }

    freq1_by_pop = _get_freq_by_pop(freq1_expr, freq1_meta_expr)
    freq2_by_pop = _get_freq_by_pop(freq2_expr, freq2_meta_expr)

    # Create list of 2x2 matrices of reference and alternate allele counts for each
    # genetic ancestry group in pops and export to pandas for CMH test.
    _ht = ht.select(
        an=[
            hl.bind(
                lambda x, y: hl.or_missing(
                    (x.AC > 0) | (y.AC > 0), [[x.AC, x.AN - x.AC], [y.AC, y.AN - y.AC]]
                ),
                freq1_by_pop[pop],
                freq2_by_pop[pop],
            )
            for pop in pops
        ]
    )

    # Remove any missing values from the list of 2x2 matrices and filter rows where
    # there are no 2x2 matrices.
    _ht = _ht.annotate(an=_ht.an.filter(lambda x: hl.is_defined(x)))
    _ht = _ht.filter(hl.len(_ht.an) > 0)

    # Add a temporary index to the Table to easily map values back to variants after
    # converting to a pandas DataFrame and back to a Hail Table.
    _ht = _ht.add_index(name="tmp_idx")
    _ht = _ht.annotate(tmp_idx=hl.int32(_ht.tmp_idx))
    tmp_path = hl.utils.new_temp_file("cmh_test", "parquet")

    # Export one compressed Parquet file per partition.
    df = _ht.to_spark().write.option("compression", "zstd")
    df.mode("overwrite").parquet(tmp_path)
    df = pd.read_parquet(tmp_path)

    # Perform CMH test on the list of 2x2 matrices (list of lists of 2 lists with 2
    # elements each).
    df["cmh"] = df.apply(
        lambda x: StratifiedTable(
            [list([list(a) for a in p]) for p in x.an]
        ).test_null_odds(),
        axis=1,
    )
    df["pvalue"] = df.apply(lambda x: x.cmh.pvalue, axis=1)
    df["statistic"] = df.apply(lambda x: x.cmh.statistic, axis=1)
    df = df[["tmp_idx", "pvalue", "statistic"]]

    # Convert CMH result pandas DataFrame to a Table and restructure the annotation.
    cmh_ht = hl.Table.from_pandas(df)
    cmh_ht = cmh_ht.key_by("tmp_idx")
    cmh_ht = _ht.select(
        cmh=hl.struct(
            odds_ratio=cmh_ht[_ht.tmp_idx].statistic, p_value=cmh_ht[_ht.tmp_idx].pvalue
        )
    )

    return cmh_ht[ht.key].cmh


def get_combine_faf_resources(
    overwrite: bool = False,
    test: bool = False,
    apply_release_filters: bool = False,
    include_contingency_table_test: bool = False,
    include_cochran_mantel_haenszel_test: bool = False,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the combined FAF resource creation pipeline.

    :param overwrite: Whether to overwrite existing resources. Default is False.
    :param test: Whether to use test resources. Default is False.
    :param apply_release_filters: Whether to get the resources for the filtered Tables.
    :param include_contingency_table_test: Whether to include the results from
        '--perform-contingency-table-test' on the combined FAF release Table.
    :param include_cochran_mantel_haenszel_test: Whether to include the results from
        '--perform-cochran-mantel-haenszel-test' on the combined FAF release Table.
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
        output_resources={
            "comb_freq_ht": get_combined_frequency(
                test=test, filtered=apply_release_filters
            )
        },
        input_resources={
            "generate_freq.py": {"exomes_ht": get_freq()},
            "generate_freq_genomes.py": {"genomes_ht": get_freq(data_type="genomes")},
            "final_filter.py": {"exomes_filter_ht": final_filter(data_type="exomes")},
            "final_filter_genomes.py": {
                "genomes_filter_ht": final_filter(data_type="genomes")
            },
        },
    )
    contingency_table_test = PipelineStepResourceCollection(
        "--perform-contingency-table-test",
        output_resources={
            "contingency_table_ht": get_freq_comparison(
                "contingency_table_test", test=test, filtered=apply_release_filters
            )
        },
        pipeline_input_steps=[combined_frequency],
    )
    cmh_test = PipelineStepResourceCollection(
        "--perform-cochran-mantel-haenszel-test",
        output_resources={
            "cmh_ht": get_freq_comparison(
                "cmh_test", test=test, filtered=apply_release_filters
            )
        },
        pipeline_input_steps=[combined_frequency],
    )
    finalize_faf_input_steps = [combined_frequency]
    if include_contingency_table_test:
        finalize_faf_input_steps.append(contingency_table_test)
    if include_cochran_mantel_haenszel_test:
        finalize_faf_input_steps.append(cmh_test)
    finalize_faf = PipelineStepResourceCollection(
        "--finalize-combined-faf-release",
        output_resources={
            "final_combined_faf_ht": get_combined_faf_release(
                test=test, filtered=apply_release_filters
            )
        },
        pipeline_input_steps=finalize_faf_input_steps,
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
    hl._set_flags(use_ssa_logs="1")
    test_gene = args.test_gene
    overwrite = args.overwrite
    apply_release_filters = args.apply_release_filters
    pops = list(set(POPS["v3"] + POPS["v4"]))
    faf_pops = [pop for pop in pops if pop not in POPS_TO_REMOVE_FOR_POPMAX["v4"]]
    combine_faf_resources = get_combine_faf_resources(
        overwrite,
        test_gene,
        apply_release_filters,
        args.include_contingency_table_test,
        args.include_cochran_mantel_haenszel_test,
    )

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
            if apply_release_filters:
                genomes_ht = genomes_ht.annotate(
                    filters=res.genomes_filter_ht.ht()[genomes_ht.key].filters
                )
                exomes_ht = exomes_ht.annotate(
                    filters=res.exomes_filter_ht.ht()[exomes_ht.key].filters
                )

            exomes_ht = extract_freq_info(exomes_ht, "exomes", apply_release_filters)
            genomes_ht = extract_freq_info(genomes_ht, "genomes", apply_release_filters)

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
            ht = ht.filter(
                hl.is_defined(ht.genomes_freq) & hl.is_defined(ht.exomes_freq)
            )
            ht = ht.select(
                contingency_table_test=perform_contingency_table_test(
                    ht.genomes_freq,
                    ht.exomes_freq,
                    ht.genomes_freq_meta,
                    ht.exomes_freq_meta,
                    ht.joint_freq_meta,
                    min_cell_count=args.min_cell_count,
                )
            )
            ht.write(res.contingency_table_ht.path, overwrite=overwrite)

        if args.perform_cochran_mantel_haenszel_test:
            res = combine_faf_resources.cmh_test
            res.check_resource_existence()
            ht = res.comb_freq_ht.ht()
            ht = ht.filter(
                hl.is_defined(ht.genomes_freq) & hl.is_defined(ht.exomes_freq)
            )
            ht = ht.select(
                cochran_mantel_haenszel_test=perform_cmh_test(
                    ht,
                    ht.genomes_freq,
                    ht.exomes_freq,
                    ht.genomes_freq_meta,
                    ht.exomes_freq_meta,
                    pops=faf_pops,
                ),
            )
            ht.write(res.cmh_ht.path, overwrite=overwrite)

        if args.finalize_combined_faf_release:
            res = combine_faf_resources.finalize_faf
            res.check_resource_existence()

            ht = res.comb_freq_ht.ht()
            stats_expr = {}
            if args.include_contingency_table_test:
                stats_expr["contingency_table_test"] = res.contingency_table_ht.ht()[
                    ht.key
                ].contingency_table_test
            if args.include_cochran_mantel_haenszel_test:
                stats_expr["cochran_mantel_haenszel_test"] = res.cmh_ht.ht()[
                    ht.key
                ].cochran_mantel_haenszel_test
            ht = ht.annotate(
                **stats_expr,
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
            ht.naive_coalesce(args.n_partitions).write(
                res.final_combined_faf_ht.path, overwrite=overwrite
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compute_combined_faf"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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
            "adj for all genetic ancestry groups found in both the exomes and genomes. "
            "The table also includes FAF computed on the joint frequencies."
        ),
        action="store_true",
    )
    parser.add_argument(
        # TODO: consider flipping this to be --skip-apply-release-filters
        "--apply-release-filters",
        help="Whether to apply the final release filters to the Table.",
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
        default=5,
    )
    parser.add_argument(
        "--perform-cochran-mantel-haenszel-test",
        help=(
            "Perform the Cochran–Mantel–Haenszel test, a stratified test of "
            "independence for 2x2xK contingency tables, on the allele "
            "frequencies where K is the number of genetic ancestry groups with FAF "
            "computed."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--finalize-combined-faf-release",
        help="Finalize the combined FAF Table for release.",
        action="store_true",
    )
    finalize_combined_faf = parser.add_argument_group(
        "Create finalized combined FAF release Table.",
        "Arguments for finalizing the combined FAF release Table.",
    )
    finalize_combined_faf.add_argument(
        "--include-contingency-table-test",
        help=(
            "Whether to include the results from '--perform-contingency-table-test' on "
            "the combined FAF release Table"
        ),
        action="store_true",
    )
    finalize_combined_faf.add_argument(
        "--include-cochran-mantel-haenszel-test",
        help=(
            "Whether to include the results from"
            " '--perform-cochran-mantel-haenszel-test' on the combined FAF release"
            " Table"
        ),
        action="store_true",
    )
    finalize_combined_faf.add_argument(
        "--n-partitions",
        help=(
            "Number of partitions to repartition the finalized combined FAF release "
            "Table to."
        ),
        type=int,
        default=10000,
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
