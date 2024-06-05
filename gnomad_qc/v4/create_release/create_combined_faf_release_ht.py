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
from typing import Callable, Dict, List, Optional, Set

import hail as hl
from gnomad.resources.grch38.gnomad import POPS, POPS_TO_REMOVE_FOR_POPMAX
from gnomad.resources.resource_utils import TableResource
from gnomad.utils.annotations import (
    faf_expr,
    gen_anc_faf_max_expr,
    merge_freq_arrays,
    merge_histograms,
    pop_max_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.annotations.compute_coverage import adjust_interval_padding
from gnomad_qc.v4.resources.annotations import (
    get_all_sites_an_and_qual_hists,
    get_combined_frequency,
    get_freq,
    get_freq_comparison,
)
from gnomad_qc.v4.resources.basics import (
    calling_intervals,
    get_checkpoint_path,
    get_logging_path,
)
from gnomad_qc.v4.resources.constants import (
    CURRENT_FREQ_VERSION,
    DATA_TYPES,
    RELEASE_DATA_TYPES,
)
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v4.resources.sample_qc import interval_qc_pass
from gnomad_qc.v4.resources.variant_qc import final_filter

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("compute_combined_faf")
logger.setLevel(logging.INFO)

CHR_LIST = [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]
"""List of chromosomes in the combined FAF release."""


def filter_gene_to_test(ht: hl.Table, pcsk9: bool, zfy: bool) -> hl.Table:
    """
    Filter to PCSK9 1:55039447-55064852 and/or ZFY Y:2935281-2982506 for testing.

    :param ht: Table with frequency and FAF information.
    :param pcsk9: Whether to filter to PCSK9 1:55039447-55064852.
    :param zfy: Whether to filter to ZFY Y:2935281-2982506.
    :return: Table with frequency and FAF information of the filtered interval of a gene
    """
    filter_loci = []
    if pcsk9:
        logger.info("Filtering to a subset of variants on PCSK9 on chr1...")
        filter_loci.append("chr1:55039447-55064852")
    if zfy:
        logger.info("Filtering to a subset of variants on ZFY on chrY...")
        filter_loci.append("chrY:2935281-2982506")

    return hl.filter_intervals(
        ht,
        [hl.parse_locus_interval(l, reference_genome="GRCh38") for l in filter_loci],
    )


def extract_freq_info(
    ht: hl.Table, prefix: str, apply_release_filters: bool = True
) -> hl.Table:
    """
    Extract frequencies and FAF for adj, raw (only for frequencies), adj by pop, adj by sex, and adj by pop/sex.

    The following annotations are renamed and where applicable, filtered:
        - freq: {prefix}_freq
        - faf: {prefix}_faf
        - grpmax: {prefix}_grpmax
        - fafmax: {prefix}_fafmax
        - qual_hists: {prefix}_qual_hists
        - raw_qual_hists: {prefix}_raw_qual_hists
        - age_hists: {prefix}_age_hists

    The following global annotations are filtered and renamed:
        - freq_meta: {prefix}_freq_meta
        - freq_index_dict: {prefix}_freq_index_dict
        - faf_meta: {prefix}_faf_meta
        - faf_index_dict: {prefix}_faf_index_dict
        - age_distribution: {prefix}_age_distribution

    If `apply_release_filters` is True, a {prefix}_filters annotation is added to the Table and the following variants are filtered:
        - chrM
        - AS_lowqual sites (these sites are dropped in the
          final_filters HT so will not have information in `filters`,
          hl.is_defined(ht.filters) is used)
        - AC_raw == 0

    :param ht: Table with frequency and FAF information.
    :param prefix: Prefix to add to each of the filtered annotations.
    :param apply_release_filters: Whether to apply the final release filters to the
        Table. Default is True.
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
    freq_index_dict = make_freq_index_dict_from_meta(hl.literal(freq_meta))
    logger.info(
        "Keeping only FAF for adj, adj by pop, adj by sex, and adj by pop/sex..."
    )
    # Exomes FAF was calculated for the whole gnomad and the ukb subset, we only
    # want to include the whole dataset in the joint release.
    faf_meta, faf = filter_arrays_by_meta(
        ht.faf_meta,
        {"faf": ht.faf},
        ["group", "gen_anc", "sex"],
        combine_operator="or",
        exact_match=True,
    )
    faf_index_dict = make_freq_index_dict_from_meta(hl.literal(faf_meta))

    # Select grpmax and fafmax
    grpmax_expr = ht.grpmax
    fafmax_expr = ht.gen_anc_faf_max
    if prefix == "exomes":
        # Note: The `grpmax` and `fafmax` structs in the exomes freq HT have two nested
        # structs: `gnomad` and `non_ukb`. This section selects only the `gnomad`
        # values (values across full v4 exomes release).
        grpmax_expr = grpmax_expr.gnomad
        fafmax_expr = fafmax_expr.gnomad

    # Drop 'faf95' annotation in the grpmax_expr if it exists. We no longer add this to
    # releases.
    if "faf95" in grpmax_expr:
        grpmax_expr = grpmax_expr.drop("faf95")

    # Rename filtered annotations with supplied prefix.
    ht = ht.select(
        **{f"{prefix}_filters": ht.filters} if apply_release_filters else {},
        **{
            f"{prefix}_freq": array_exprs["freq"],
            f"{prefix}_faf": faf["faf"],
            f"{prefix}_grpmax": grpmax_expr,
            f"{prefix}_fafmax": fafmax_expr,
            f"{prefix}_qual_hists": ht.histograms.qual_hists,
            f"{prefix}_raw_qual_hists": ht.histograms.raw_qual_hists,
            f"{prefix}_age_hists": ht.histograms.age_hists,
        },
    )
    ht = ht.select_globals(
        **{
            f"{prefix}_freq_meta": freq_meta,
            f"{prefix}_freq_index_dict": freq_index_dict,
            f"{prefix}_freq_meta_sample_count": array_exprs["freq_meta_sample_count"],
            f"{prefix}_faf_meta": faf_meta,
            f"{prefix}_faf_index_dict": faf_index_dict,
            f"{prefix}_age_distribution": ht.age_distribution,
        }
    )

    ht = ht.checkpoint(hl.utils.new_temp_file("extract_freq_info", "ht"))

    return ht


def add_all_sites_an_and_qual_hists(
    ht: hl.Table,
    exomes_all_sites_ht: hl.Table,
    genomes_all_sites_ht: hl.Table,
) -> hl.Table:
    """
    Add all sites AN and qual hists to the Table.

    :param ht: Table with frequency and FAF information.
    :param exomes_all_sites_ht: Table with all sites AN and qual hists for exomes.
    :param genomes_all_sites_ht: Table with all sites AN and qual hists for genomes.
    :return: Table with all sites AN and qual hists added.
    """

    def _get_freq_from_an_expr(
        freq_expr: hl.expr.ArrayExpression,
        all_sites_an_expr: hl.expr.ArrayExpression,
        freq_meta_expr: hl.expr.ArrayExpression,
        freq_index_dict_expr: hl.expr.DictExpression,
        an_meta_expr: hl.expr.ArrayExpression,
    ) -> hl.expr.ArrayExpression:
        """
        Get freq expression with all sites AN added if `freq_expr` is missing.

        :param freq_expr: ArrayExpression of frequencies.
        :param all_sites_an_expr: ArrayExpression of all sites AN.
        :param freq_meta_expr: Frequency metadata.
        :param freq_index_dict_expr: Frequency index dictionary.
        :param an_meta_expr: AN metadata.
        :return: ArrayExpression of frequencies with all sites AN added if `freq_expr`
            is missing.
        """
        meta_map = hl.dict(hl.enumerate(an_meta_expr).map(lambda x: (x[1], x[0])))
        all_sites_an_expr = freq_meta_expr.map(
            lambda x: hl.bind(
                lambda an: hl.struct(
                    AC=hl.or_missing(hl.is_defined(an), 0),
                    AF=hl.or_missing(an > 0, hl.float64(0.0)),
                    AN=hl.int32(an),
                    homozygote_count=hl.or_missing(hl.is_defined(an), 0),
                ),
                all_sites_an_expr[meta_map[x]],
            )
        )
        logger.info("Setting XX samples call stats to missing on chrY...")
        all_sites_an_expr = set_female_y_metrics_to_na_expr(
            ht,
            freq_expr=all_sites_an_expr,
            freq_meta_expr=freq_meta_expr,
            freq_index_dict_expr=freq_index_dict_expr,
        )

        return hl.coalesce(freq_expr, all_sites_an_expr)

    logger.info("Adding all sites AN and qual hists information where missing...")
    all_sites_by_data_type = {
        "exomes": exomes_all_sites_ht[ht.locus],
        "genomes": genomes_all_sites_ht[ht.locus],
    }
    all_sites_meta_by_data_type = {
        "exomes": exomes_all_sites_ht.index_globals().strata_meta,
        "genomes": genomes_all_sites_ht.index_globals().strata_meta,
    }
    ht = ht.annotate(
        **{
            f"{data_type}_freq": _get_freq_from_an_expr(
                ht[f"{data_type}_freq"],
                all_sites.AN,
                ht.index_globals()[f"{data_type}_freq_meta"],
                ht.index_globals()[f"{data_type}_freq_index_dict"],
                all_sites_meta_by_data_type[data_type],
            )
            for data_type, all_sites in all_sites_by_data_type.items()
        },
        **{
            f"{data_type}_{adj}qual_hists": hl.struct(
                **{
                    k: hl.coalesce(
                        ht[f"{data_type}_{adj}qual_hists"][k],
                        all_sites.qual_hists[f"{adj}qual_hists"].get(
                            k,
                            hl.missing(
                                hl.tstruct(
                                    bin_edges=hl.tarray(hl.tfloat64),
                                    bin_freq=hl.tarray(hl.tint64),
                                    n_smaller=hl.tint64,
                                    n_larger=hl.tint64,
                                )
                            ),
                        ),
                    )
                    for k in ht[f"{data_type}_{adj}qual_hists"]
                }
            )
            for data_type, all_sites in all_sites_by_data_type.items()
            for adj in ["", "raw_"]
        },
    )

    return ht


def get_joint_freq_and_faf(
    genomes_ht: hl.Table,
    exomes_ht: hl.Table,
    genomes_all_sites_ht: hl.Table,
    exomes_all_sites_ht: hl.Table,
    faf_pops_to_exclude: Set[str] = POPS_TO_REMOVE_FOR_POPMAX["v4"],
) -> hl.Table:
    """
    Get joint genomes and exomes frequency and FAF information.

    :param genomes_ht: Table with genomes frequency and FAF information.
    :param exomes_ht: Table with exomes frequency and FAF information.
    :param genomes_all_sites_ht: Table with all sites AN and qual hists for genomes.
    :param exomes_all_sites_ht: Table with all sites AN and qual hists for exomes.
    :param faf_pops_to_exclude: Set of genetic ancestry groups to exclude from the FAF
        calculation.
    :return: Table with joint genomes and exomes frequency and FAF information.
    """
    logger.info("Performing an outer join on frequency HTs...")
    ht = genomes_ht.join(exomes_ht, how="outer")
    ht = add_all_sites_an_and_qual_hists(ht, exomes_all_sites_ht, genomes_all_sites_ht)

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

    ht = ht.annotate(
        joint_freq=freq,
        **{
            f"joint_{hist_struct}": hl.struct(
                **{
                    h: merge_histograms(
                        [ht[f"{d}_{hist_struct}"][h] for d in DATA_TYPES]
                    )
                    for h in ht[f"genomes_{hist_struct}"]
                }
            )
            for hist_struct in ["qual_hists", "raw_qual_hists", "age_hists"]
        },
    )
    ht = ht.annotate_globals(
        joint_freq_meta=freq_meta,
        joint_freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta)),
    )
    # NOTE: This checkpoint prevents a Class Too Large error in hail 0.2.122.
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

    # Compute group max (popmax) on the merged exomes + genomes frequencies.
    grpmax = pop_max_expr(
        ht.joint_freq,
        ht.joint_freq_meta,
        pops_to_exclude=faf_pops_to_exclude,
        pop_label="gen_anc",
    )

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
        joint_age_distribution=merge_histograms(
            [ht.genomes_age_distribution, ht.exomes_age_distribution]
        ),
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
                    (f1.AC > 0) | (f2.AC > 0) & (f1.AN > 0) & (f2.AN > 0),
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
            if m.get("gen_anc") in pops and len(m) == 2
        }

    freq1_by_pop = _get_freq_by_pop(freq1_expr, freq1_meta_expr)
    freq2_by_pop = _get_freq_by_pop(freq2_expr, freq2_meta_expr)

    # Create list of the info for the 2x2 matrices of reference and alternate allele
    # counts for each genetic ancestry group in pops and remove any missing pops
    # from the list.
    cmh_input_expr = hl.array(
        [
            hl.bind(
                lambda x, y: hl.or_missing(
                    ((x.AC > 0) | (y.AC > 0)) & (x.AN > 0) & (y.AN > 0),
                    ([x.AC, x.AN - x.AC, y.AC, y.AN - y.AC], pop),
                ),
                freq1_by_pop[pop],
                freq2_by_pop[pop],
            )
            for pop in pops
        ]
    ).filter(hl.is_defined)

    return hl.or_missing(
        hl.len(cmh_input_expr) > 0,
        hl.cochran_mantel_haenszel_test(
            *[cmh_input_expr.map(lambda x: x[0][i]) for i in range(4)]
        )
        .annotate(gen_ancs=cmh_input_expr.map(lambda x: x[1]))
        .rename({"test_statistic": "chisq"}),
    )


def create_final_combined_faf_release(
    ht,
    contingency_table_ht: hl.Table,
    cmh_ht: hl.Table,
) -> hl.Table:
    """
    Create the final combined FAF release Table.

    :param ht: Table with joint exomes and genomes frequency and FAF information.
    :param contingency_table_ht: Table with contingency table test results to include
        on the final Table.
    :param cmh_ht: Table with Cochran–Mantel–Haenszel test results to include on the
        final Table.
    :return: Table with final combined FAF release information.
    """
    logger.info("Creating final combined FAF release Table...")

    def _get_struct_by_data_type(
        ht: hl.Table,
        annotation_expr: Optional[hl.expr.StructExpression] = None,
        condition_func: Optional[Callable[[str], hl.expr.BooleanExpression]] = None,
        postfix: str = "",
    ) -> hl.expr.StructExpression:
        """
        Group all annotations from the same data type into a struct.

        :param ht: Table with annotations to loop through.
        :param annotation_expr: StructExpression with annotations to loop through. If
            None, use ht.row.
        :param condition_func: Function that returns a BooleanExpression for a
            condition to check in addition to data type.
        :param postfix: String to add to the end of the new annotation following data
            type.
        :return: StructExpression of joint data type and the corresponding data type
            struct.
        """
        if condition_func is None:
            condition_func = lambda x: True
        if annotation_expr is None:
            annotation_expr = ht.row

        return hl.struct(
            **{
                f"{data_type}{postfix}": hl.struct(
                    **{
                        x.replace(f"{data_type}_", ""): ht[x]
                        for x in annotation_expr
                        if x.startswith(data_type) and condition_func(x)
                    }
                )
                for data_type in RELEASE_DATA_TYPES
            }
        )

    # Group all the histograms into a single struct per data type.
    ht = ht.transmute(
        **_get_struct_by_data_type(
            ht, condition_func=lambda x: x.endswith("_hists"), postfix="_histograms"
        )
    )

    # Add capture and calling intervals as region flags.
    intervals = [(c, calling_intervals(c, 0).ht()) for c in ["broad", "ukb"]]
    intervals = [(f"{c}_capture", i) for c, i in intervals] + [
        (f"{c}_calling", adjust_interval_padding(i, 150)) for c, i in intervals
    ]

    region_expr = {
        "fail_interval_qc": (
            ~interval_qc_pass(all_platforms=True).ht()[ht.locus].pass_interval_qc
        ),
        **{f"outside_{c}_region": hl.is_missing(i[ht.locus]) for c, i in intervals},
        **{
            f"not_called_in_{d}": hl.is_missing(ht[f"{d}_freq"]) | hl.is_missing(
                ht[f"{d}_freq"][1].AN
            )
            for d in DATA_TYPES
        },
    }

    # Get the stats expressions for the final Table.
    def _prep_stats_annotation(ht: hl.Table, stat_ht: hl.Table, stat_field: str):
        """
        Prepare stats expressions for the final Table.

        :param ht: Table to add stats to.
        :param stat_ht: Table with stats to add to the final Table.
        :param stat_field: Field to add to the final Table.
        :return: Dictionary of stats expressions for the final Table.
        """
        stat_expr = stat_ht[ht.key][stat_field]
        is_struct = False
        if isinstance(stat_expr, hl.expr.StructExpression):
            is_struct = True
            stat_expr = hl.array([stat_expr])

        stat_expr = stat_expr.map(
            lambda x: x.select(
                **{
                    a: hl.or_missing(~hl.is_nan(x[a]), x[a])
                    for a in x
                    if a != "gen_ancs"
                }
            )
        )

        if is_struct:
            stat_expr = stat_expr[0]

        return stat_expr

    # Prepare stats annotations for the final Table.
    ct_expr = _prep_stats_annotation(ht, contingency_table_ht, "contingency_table_test")
    cmh_expr = _prep_stats_annotation(ht, cmh_ht, "cochran_mantel_haenszel_test")
    anc_expr = cmh_ht[ht.key]["cochran_mantel_haenszel_test"].gen_ancs
    ct_union_expr = ct_expr[ht.joint_freq_index_dict[anc_expr[0] + "_adj"]]

    # Create a union stats annotation that uses the contingency table test if there is
    # only one genetic ancestry group, otherwise use the CMH test.
    union_expr = hl.if_else(
        hl.len(anc_expr) == 1,
        ct_union_expr.select("p_value").annotate(
            stat_test_name="contingency_table_test", gen_ancs=anc_expr
        ),
        cmh_expr.select("p_value").annotate(
            stat_test_name="cochran_mantel_haenszel_test", gen_ancs=anc_expr
        ),
    )
    stats_expr = {
        "contingency_table_test": ct_expr,
        "cochran_mantel_haenszel_test": cmh_expr,
        "stat_union": union_expr,
    }

    # Select the final columns for the Table, grouping all annotations from the same
    # data type into a struct.
    ht = ht.select(
        region_flags=hl.struct(**region_expr),
        **_get_struct_by_data_type(ht),
        freq_comparison_stats=hl.struct(**stats_expr),
    )

    # Group all global annotations from the same data type into a struct.
    ht = ht.transmute_globals(
        **_get_struct_by_data_type(ht, annotation_expr=ht.globals, postfix="_globals")
    )

    return ht.drop("versions")


def get_combine_faf_resources(
    overwrite: bool = False,
    test: bool = False,
    filtered: bool = True,
    stats_chr: str = None,
    stats_combine_all_chr: bool = False,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the combined FAF resource creation pipeline.

    :param overwrite: Whether to overwrite existing resources. Default is False.
    :param test: Whether to use test resources. Default is False.
    :param filtered: Whether to get the resources for the filtered Tables. Default is
        True.
    :param stats_chr: Chromosome to get temp stats resource for. Default is None, which
        will return the resources for the stats on the full exome/genome.
    :param stats_combine_all_chr: Whether to also get the stats resources for all
        chromosomes to be combined. Default is False.
    :return: PipelineResourceCollection containing resources for all steps of the
        combined FAF resource creation pipeline.
    """
    if stats_chr is not None and stats_combine_all_chr:
        raise ValueError(
            "Cannot specify both `stats_chr` and `stats_combine_all_chr`. Please "
            "specify only one."
        )

    # Initialize pipeline resource collection.
    combine_faf_pipeline = PipelineResourceCollection(
        pipeline_name="combine_faf",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the pipeline.
    combined_frequency = PipelineStepResourceCollection(
        "--create-combined-frequency-table",
        output_resources={
            "comb_freq_ht": get_combined_frequency(test=test, filtered=filtered)
        },
        input_resources={
            "generate_freq.py": {"exomes_ht": get_freq()},
            "generate_freq_genomes.py": {"genomes_ht": get_freq(data_type="genomes")},
            "final_filter.py": {"exomes_filter_ht": final_filter(data_type="exomes")},
            "final_filter_genomes.py": {
                "genomes_filter_ht": final_filter(data_type="genomes")
            },
            "all sites AN and qual hists HTs": {
                "exomes_all_sites_ht": get_all_sites_an_and_qual_hists(
                    data_type="exomes"
                ),
                "genomes_all_sites_ht": get_all_sites_an_and_qual_hists(
                    data_type="genomes"
                ),
            },
        },
    )

    cmh_resources = {}
    contingency_test_resources = {}
    if stats_chr is not None:
        cmh_resources["output_resources"] = {
            "cmh_ht": TableResource(get_checkpoint_path(f"cmh_{stats_chr}"))
        }
        contingency_test_resources["output_resources"] = {
            "contingency_table_ht": TableResource(
                get_checkpoint_path(f"contingency_table_{stats_chr}")
            )
        }
    else:
        cmh_resources["output_resources"] = {
            "cmh_ht": get_freq_comparison("cmh_test", test=test, filtered=filtered)
        }
        contingency_test_resources["output_resources"] = {
            "contingency_table_ht": get_freq_comparison(
                "contingency_table_test", test=test, filtered=filtered
            )
        }

    if stats_combine_all_chr:
        cmh_resources["input_resources"] = {
            "--perform-cochran-mantel-haenszel-test --stats-chr (any chr)": {
                f"{c}_cmh_ht": TableResource(get_checkpoint_path(f"cmh_{c}"))
                for c in CHR_LIST
            }
        }
        contingency_test_resources["input_resources"] = {
            "--perform-contingency-table-test --stats-chr (all chr)": {
                f"{c}_contingency_table_ht": TableResource(
                    get_checkpoint_path(f"contingency_table_{c}")
                )
                for c in CHR_LIST
            }
        }
    else:
        cmh_resources["pipeline_input_steps"] = [combined_frequency]
        contingency_test_resources["pipeline_input_steps"] = [combined_frequency]

    contingency_table_test = PipelineStepResourceCollection(
        "--perform-contingency-table-test", **contingency_test_resources
    )
    cmh_test = PipelineStepResourceCollection(
        "--perform-cochran-mantel-haenszel-test", **cmh_resources
    )
    finalize_faf = PipelineStepResourceCollection(
        "--finalize-combined-faf-release",
        output_resources={
            "final_combined_faf_ht": release_sites(data_type="joint", test=test)
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


def main(args):
    """Create combined FAF resource."""
    hl.init(
        log="/compute_combined_faf.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl._set_flags(use_ssa_logs="1")
    test = args.test_gene or args.test_y_gene
    overwrite = args.overwrite
    apply_release_filters = not args.skip_apply_release_filters
    pops = list(set(POPS["v3"]["genomes"] + POPS["v4"]["exomes"]))
    faf_pops = [pop for pop in pops if pop not in POPS_TO_REMOVE_FOR_POPMAX["v4"]]
    stats_chr = args.stats_chr
    stats_combine_all_chr = args.stats_combine_all_chr
    combine_faf_resources = get_combine_faf_resources(
        overwrite,
        test,
        apply_release_filters,
        stats_chr,
        stats_combine_all_chr,
    )

    try:
        if args.create_combined_frequency_table:
            res = combine_faf_resources.combined_frequency
            res.check_resource_existence()
            exomes_ht = res.exomes_ht.ht()
            genomes_ht = res.genomes_ht.ht()

            if test:
                exomes_ht = filter_gene_to_test(
                    exomes_ht, pcsk9=args.test_gene, zfy=args.test_y_gene
                )
                genomes_ht = filter_gene_to_test(
                    genomes_ht, pcsk9=args.test_gene, zfy=args.test_y_gene
                )

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

            ht = get_joint_freq_and_faf(
                genomes_ht,
                exomes_ht,
                res.genomes_all_sites_ht.ht(),
                res.exomes_all_sites_ht.ht(),
            )
            ht = ht.annotate_globals(
                versions=hl.struct(
                    exomes=CURRENT_FREQ_VERSION["exomes"],
                    genomes=CURRENT_FREQ_VERSION["genomes"],
                )
            )
            ht.describe()
            ht.write(res.comb_freq_ht.path, overwrite=overwrite)

        if args.perform_contingency_table_test:
            res = combine_faf_resources.contingency_table_test
            res.check_resource_existence()
            if stats_combine_all_chr:
                ht = hl.Table.union(
                    *[
                        getattr(res, f"{c}_contingency_table_ht").ht()
                        for c in [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]
                    ]
                )
            else:
                ht = res.comb_freq_ht.ht()
                ht = ht.filter(
                    hl.is_defined(ht.genomes_freq) & hl.is_defined(ht.exomes_freq)
                )
                if stats_chr is not None:
                    ht = hl.filter_intervals(ht, [hl.parse_locus_interval(stats_chr)])

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
            if stats_combine_all_chr:
                ht = hl.Table.union(
                    *[
                        getattr(res, f"{c}_cmh_ht").ht()
                        for c in [f"chr{c}" for c in range(1, 23)] + ["chrX", "chrY"]
                    ]
                )
            else:
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
                    )
                )
            ht.naive_coalesce(args.n_partitions).write(
                res.cmh_ht.path, overwrite=overwrite
            )

        if args.finalize_combined_faf_release:
            res = combine_faf_resources.finalize_faf
            res.check_resource_existence()

            ht = create_final_combined_faf_release(
                res.comb_freq_ht.ht(), res.contingency_table_ht.ht(), res.cmh_ht.ht()
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
        "--test-y-gene",
        help="Test on a subset of variants in ZFY on chrY.",
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
        "--skip-apply-release-filters",
        help="Whether to skip applying the final release filters to the Table.",
        action="store_true",
    )
    parser.add_argument(
        "--stats-chr",
        help="Chromosome to compute stats on.",
        type=str,
        choices=CHR_LIST,
    )
    parser.add_argument(
        "--stats-combine-all-chr",
        help=(
            "Whether to combined all chromosome stats. The stats calculations must "
            "have been performed for each chromosome using the `--stats-chr` argument."
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
