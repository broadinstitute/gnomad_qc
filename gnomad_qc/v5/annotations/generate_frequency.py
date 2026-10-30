"""
Generate gnomAD v5 genome frequencies.

gnomAD v5 adds no new gnomAD samples, so ``--process-gnomad`` reuses the v4 counts:
it computes AC, homozygote counts and age histograms for the samples that withdrew
consent after v4 and subtracts them, and keeps the v4 quality histograms as they are.
All of Us (AoU) via ``--process-aou`` aggregates AC, homozygote counts and histograms
from the AoU VDS variant data. ``--merge-datasets`` combines the two tables, takes AN
for every stratum from the all-sites AN release ``compute_coverage.py`` wrote, and
adds AF, FAF, grpmax and the inbreeding coefficient.

Genome-wide, ``--use-batch-fanout`` fans the AoU step out in Query-on-Batch using the
interval chunk layout ``compute_coverage.py`` already wrote beside its AN table, and
``--merge-freq-chunks`` tree-merges the chunk tables.

Examples::

    # AoU: fan out, then merge.
    generate_frequency.py --process-aou --use-batch-fanout --environment batch
    generate_frequency.py --process-aou --merge-freq-chunks --environment batch
    # AoU as a single job on a test region.
    generate_frequency.py --process-aou --test-region chr22:10510002-11110002
    # gnomAD, then the merge.
    generate_frequency.py --process-gnomad --environment batch
    generate_frequency.py --merge-datasets --environment batch
"""

import argparse
import logging
from collections.abc import Callable
from functools import partial
from typing import Any

import hail as hl
from gnomad.resources.grch38.gnomad import GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX
from gnomad.sample_qc.sex import adjust_sex_ploidy, adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    _read_reduction_globals,
    age_hists_expr,
    annotate_adj,
    bi_allelic_site_inbreeding_expr,
    expand_strata_array_from_leaves,
    faf_expr,
    gen_anc_faf_max_expr,
    grpmax_expr,
    merge_array_expressions,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.resources.release import release_sites as v3_release_sites
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
from gnomad_qc.v5.annotations.batch_fanout import (
    RelayJobSpec,
    apply_path_suffix,
    chunk_intervals_hash,
    chunk_intervals_path,
    chunk_path,
    chunk_read_intervals,
    combine_suffix,
    fold_chunk_layout,
    group_path,
    list_present_chunk_indices,
    load_chunk_layout,
    orchestrate_tree_merge,
    parse_region_interval,
    qob_sizing_flags,
    relay_context,
    resolve_commit,
    run_chunk_fanout,
    select_chunks_for_contig,
    spans_sex_chromosome,
    split_intervals,
    submit_relay_batch,
    test_region_hash,
    union_and_write_hts,
)
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_ac_homalt_and_hists,
    get_freq,
    group_membership,
)
from gnomad_qc.v5.resources.basics import (
    _BATCH_RESOURCE_PARAMS,
    _file_exists_for_env,
    _get_batch_resource_kwargs,
    _init_hail,
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.meta import meta
from gnomad_qc.v5.resources.release import release_coverage_path

# force=True: the gnomad imports already called basicConfig, so without it this call
# is a no-op and logger.info is dropped.
logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("v5_frequency")

# Hail 0.2.137 runs on the Batch workers' Java 21.
DEFAULT_BATCH_IMAGE = (
    "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:0.2.137"
)
AGE_HIST_BINS = (30, 80, 10)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _release_meta(project: str, environment: str) -> hl.Table:
    """
    Return sex karyotype and age for one project's release samples, at 10 partitions.

    The v5 meta table has ~330 partitions, so every join or aggregate against it would
    fan ~330 tiny QoB tasks.

    :param project: ``"aou"`` or ``"gnomad"``.
    :param environment: Compute environment.
    :return: Table keyed by ``s`` with ``sex_karyotype`` and ``age``.
    """
    ht = meta(data_type="genomes", environment=environment).ht()
    ht = ht.filter((ht.project_meta.project == project) & ht.release)
    return ht.select("sex_karyotype", age=ht.project_meta.age).naive_coalesce(
        10
    )  # TODO: Confirm that the naive_coalesce to 10 instead of checkpointing then naive_coalesce is ok. Seems like its small enough and only done once?


def _age_distribution(ht: hl.Table) -> hl.Struct:
    """
    Return the age histogram global for the release-meta table.

    :param ht: Table from :func:`_release_meta`.
    :return: Age histogram struct.
    """
    return ht.aggregate(hl.agg.hist(ht.age, *AGE_HIST_BINS))


def _group_membership_ht(data_set: str, test: bool, environment: str) -> hl.Table:
    """
    Read the group-membership HT ``compute_coverage.py`` built the AN with.

    The merge matches AC to AN by stratum, so the strata and samples must be the
    ones the AN was built with; a missing HT is an error rather than something to
    rebuild here. The gnomAD HT holds the consent-withdrawal samples, and its
    coverage-only ``v3_release`` stratum is dropped.

    :param data_set: ``"aou"`` or ``"gnomad"``.
    :param test: Read the test-scoped path.
    :param environment: Compute environment.
    :return: HT with a ``group_membership`` array and ``freq_meta`` and
        ``freq_meta_sample_count`` globals. AoU's is cell-reduced, with ``_full``
        globals that map cells back to the full strata.
    """
    path = apply_path_suffix(
        group_membership(test=test, data_set=data_set, environment=environment).path,
        "cells",
    )
    if not _file_exists_for_env(path, environment):
        raise FileNotFoundError(
            f"{data_set} group-membership HT not found at {path}; run"
            " compute_coverage.py --write-group-membership-ht first."
        )
    ht = hl.read_table(path)
    if data_set == "gnomad":
        freq_meta, arrays = filter_arrays_by_meta(
            ht.freq_meta,
            {
                "group_membership": ht.group_membership,
                "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
            },
            items_to_filter=["v3_release"],
            keep=False,
        )
        ht = ht.annotate(group_membership=arrays["group_membership"]).annotate_globals(
            freq_meta=freq_meta,
            freq_meta_sample_count=arrays["freq_meta_sample_count"],
        )
    return ht


def _carrier_hists(c: hl.expr.StructExpression) -> hl.expr.StructExpression:
    """
    Quality and age histogram aggregators over the carriers of one split row.

    v5 keeps only the adj quality histograms; passing ``adj_expr`` would also build
    the raw ones.

    :param c: Carrier struct with ``GT``, ``GQ``, ``AD``, ``adj`` and ``age``.
    :return: Struct of ``qual_hists`` and ``age_hists`` aggregators.
    """
    return hl.struct(
        qual_hists=hl.agg.filter(
            c.adj,
            qual_hist_expr(
                gt_expr=c.GT,
                gq_expr=c.GQ,
                dp_expr=hl.sum(c.AD),
                ab_expr=c.AD[1] / hl.sum(c.AD),
            ),
        ),
        age_hists=age_hists_expr(c.adj, c.GT, c.age),
    )


def _merge_hist_struct(
    hist1: hl.expr.StructExpression,
    hist2: hl.expr.StructExpression,
    operation: str = "sum",
    set_negatives_to_zero: bool = False,
) -> hl.expr.StructExpression:
    """
    Sum or subtract every histogram field of two histogram structs.

    A histogram missing from ``hist2`` leaves ``hist1``'s. One missing from ``hist1``
    takes ``hist2``'s for ``"sum"`` and stays missing for ``"diff"``, since there is
    nothing to subtract from.

    :param hist1: First struct of histograms.
    :param hist2: Second struct with the same fields.
    :param operation: ``"sum"`` or ``"diff"``.
    :param set_negatives_to_zero: Set a negative ``"diff"`` bin, ``n_smaller`` or
        ``n_larger`` to 0 instead of erroring. Default is False.
    :return: Struct of merged histograms.
    """

    def _merge(h1, h2):
        keep = hl.is_defined(h1)
        if operation == "sum":
            keep |= hl.is_defined(h2)
        return hl.or_missing(
            keep,
            merge_histograms(
                [h1, h2],
                operation=operation,
                set_negatives_to_zero=set_negatives_to_zero,
            ),
        )

    return hl.struct(
        **{field: _merge(hist1[field], hist2[field]) for field in hist1.dtype.fields}
    )


def _localize_carriers(mt: hl.MatrixTable, group_membership_ht: hl.Table) -> hl.Table:
    """
    Localize each row to its carriers, each listing the strata it counts toward.

    A carrier is a defined entry. It counts toward every stratum its sample is in,
    except the adj strata when the entry fails adj, so one explode and group_by in
    :func:`_aggregate_carriers` covers every stratum.

    :param mt: Variant MT with an ``adj`` entry field and an ``age`` column field.
    :param group_membership_ht: Membership HT whose ``group_membership`` array and
        ``freq_meta`` global define the strata.
    :return: Table of ``mt``'s rows with a ``carriers`` array of the defined entries,
        each with ``age`` and ``counted_in`` (stratum indices) added.
    """
    freq_meta = [dict(m) for m in hl.eval(group_membership_ht.freq_meta)]
    adj_groups = hl.literal([m.get("group", "NA") == "adj" for m in freq_meta])
    gm = group_membership_ht.select(
        strata=hl.enumerate(group_membership_ht.group_membership)
        .filter(lambda t: t[1])
        .map(lambda t: t[0])
    )
    mt = mt.annotate_cols(strata=gm[mt.col_key].strata)
    lt = mt.localize_entries("entries", "cols")
    carriers = (
        hl.enumerate(lt.entries)
        .filter(lambda t: hl.is_defined(t[1]))
        .map(
            lambda t: t[1].annotate(
                age=lt.cols[t[0]].age,
                counted_in=lt.cols[t[0]].strata.filter(
                    lambda g: ~adj_groups[g] | hl.coalesce(t[1].adj, False)
                ),
            )
        )
    )
    return lt.select(carriers=carriers).drop("cols")


def _aggregate_carriers(
    lt: hl.Table,
    group_membership_ht: hl.Table,
    hists: Callable[[hl.expr.StructExpression], hl.expr.StructExpression],
) -> hl.Table:
    """
    Aggregate per-stratum AC and homozygote counts and histograms over carriers.

    Only diploid hom-alt calls are homozygotes, as in ``hl.agg.call_stats`` and so
    in v4, so hemizygous XY calls on chrX non-PAR and chrY are not.

    :param lt: Table from :func:`_localize_carriers`, keyed by the split
        ``locus, alleles``, with a ``GT`` on each carrier.
    :param group_membership_ht: The membership HT ``lt``'s ``counted_in`` indexes.
    :param hists: Returns the histogram aggregators over one carrier.
    :return: Table with ``hist_fields`` and ``strata`` (array of ``AC`` /
        ``homozygote_count`` structs, one per stratum).
    """
    n_groups = hl.eval(hl.len(group_membership_ht.freq_meta))
    empty = hl.struct(AC=hl.int64(0), homozygote_count=hl.int64(0))
    lt = lt.annotate(
        hist_fields=lt.carriers.aggregate(hists),
        _by_stratum=lt.carriers.aggregate(
            lambda c: hl.agg.explode(
                lambda g: hl.agg.group_by(
                    g,
                    hl.struct(
                        AC=hl.agg.sum(c.GT.n_alt_alleles()),
                        homozygote_count=hl.agg.count_where(
                            c.GT.is_diploid() & c.GT.is_hom_var()
                        ),
                    ),
                ),
                c.counted_in,
            )
        ),
    )
    return lt.select(
        "hist_fields",
        strata=hl.range(n_groups).map(lambda g: lt._by_stratum.get(g, empty)),
    )


def _int32_counts(strata: hl.expr.ArrayExpression) -> hl.expr.ArrayExpression:
    """
    Cast per-stratum counts to int32, the v4 release's type.

    :param strata: Array of ``AC`` / ``homozygote_count`` structs.
    :return: The array with int32 fields.
    """
    return strata.map(
        lambda s: hl.struct(
            AC=hl.int32(s.AC), homozygote_count=hl.int32(s.homozygote_count)
        )
    )


def _merge_ac_homalt(
    arrays: list[hl.expr.ArrayExpression],
    metas: list[hl.expr.ArrayExpression],
    sample_counts: list[hl.expr.ArrayExpression],
    operation: str = "sum",
    set_negatives_to_zero: bool = False,
) -> tuple[hl.expr.ArrayExpression, list[dict[str, str]], hl.expr.ArrayExpression]:
    """
    Sum or subtract per-stratum AC and homozygote counts, matching strata by meta.

    A missing array or stratum counts as 0. The output keeps the first array's
    strata order; ``"sum"`` appends strata found only in later arrays, ``"diff"``
    keeps only the first array's strata and errors on a negative result.

    :param arrays: ``freq`` arrays of ``AC`` / ``homozygote_count`` structs.
    :param metas: ``freq_meta`` of each array.
    :param sample_counts: ``freq_meta_sample_count`` of each array.
    :param operation: ``"sum"`` or ``"diff"``.
    :param set_negatives_to_zero: Set a negative ``"diff"`` AC or homozygote count to
        0 instead of erroring; a negative sample count still errors. Default is False.
    :return: ``(freq, freq_meta, freq_meta_sample_count)``.
    """
    freq, freq_meta = merge_array_expressions(
        arrays,
        metas,
        operation=operation,
        set_negatives_to_zero=set_negatives_to_zero,
        struct_fields=["AC", "homozygote_count"],
    )
    # A negative sample count is a strata or sample mismatch, so it always errors.
    counts, _ = merge_array_expressions(sample_counts, metas, operation=operation)
    return freq, freq_meta, counts


# ---------------------------------------------------------------------------
# AoU
# ---------------------------------------------------------------------------


def _prepare_aou_vds(
    aou_vds: hl.vds.VariantDataset,
    group_membership_ht: hl.Table,
    environment: str,
    skip_sex_ploidy: bool,
) -> hl.MatrixTable:
    """
    Prepare the AoU variant data for the frequency aggregation.

    Rows stay unsplit, so adj is computed on ``LGT``/``LAD``;
    :func:`_sparse_split_strata_and_hists` compacts to carriers before it splits.

    :param aou_vds: AoU VDS filtered to release samples.
    :param group_membership_ht: The cells HT from :func:`_group_membership_ht`.
    :param environment: Compute environment.
    :param skip_sex_ploidy: Skip the sex-ploidy adjustment; only valid when the
        read scope has no chrX/chrY loci.
    :return: Unsplit variant MT with ``LGT``, ``GQ``, ``LAD``, ``LA``, ``adj``
        entries, ``sex_karyotype`` and ``age`` on the columns, and the strata globals.
    """
    vmt = aou_vds.variant_data
    release = _release_meta("aou", environment)
    m = release[vmt.col_key]
    vmt = vmt.select_cols(sex_karyotype=m.sex_karyotype, age=m.age)
    if skip_sex_ploidy:
        lgt = vmt.LGT
    else:
        lgt = adjusted_sex_ploidy_expr(vmt.locus, vmt.LGT, vmt.sex_karyotype)
    vmt = vmt.select_entries(LGT=lgt, GQ=vmt.GQ, LAD=vmt.LAD, LA=vmt.LA)
    vmt = annotate_adj_no_dp(vmt)
    gg = group_membership_ht.index_globals()
    return vmt.select_globals(
        freq_meta=gg.freq_meta_full,
        freq_meta_sample_count=gg.freq_meta_sample_count_full,
        age_distribution=_age_distribution(release),
        downsamplings=gg.downsamplings,
    )


def _sparse_split_strata_and_hists(
    mt: hl.MatrixTable, group_membership_ht: hl.Table
) -> hl.Table:
    """
    Split multi-allelics and aggregate per-stratum counts and histograms over carriers.

    The entries array has a slot per sample but most variants have few carriers, so
    splitting first would process mostly missing entries once per split row.
    Compacting to carriers first makes every later step proportional to carriers.
    The split matches ``sparse_split_multi`` with ``filter_changed_loci=True``, except
    that bi-allelic rows are min-repped too.

    :param mt: Prepared, unsplit variant MT from :func:`_prepare_aou_vds`.
    :param group_membership_ht: Membership HT whose ``group_membership`` array and
        ``freq_meta`` global define the strata (cells).
    :return: Table keyed by the split ``locus, alleles`` with ``hist_fields`` and
        ``strata`` (array of ``AC`` / ``homozygote_count`` structs, one per stratum)
        and ``mt``'s globals.
    """
    lt = _localize_carriers(mt, group_membership_ht)

    def _split_struct(i):
        return hl.bind(
            lambda mr: hl.or_missing(
                mr.locus == lt.locus,
                hl.struct(locus=lt.locus, alleles=mr.alleles, a_index=i),
            ),
            hl.min_rep(lt.locus, [lt.alleles[0], lt.alleles[i]]),
        )

    splits = hl.if_else(
        hl.len(lt.alleles) < 2,
        [hl.struct(locus=lt.locus, alleles=lt.alleles, a_index=1)],
        hl._sort_by(
            hl.range(1, hl.len(lt.alleles)).map(_split_struct).filter(hl.is_defined),
            lambda l, r: hl._compare(l.alleles, r.alleles) < 0,
        ),
    )

    def _split_carrier(e, a_index):
        lai = hl.fold(
            lambda acc, k: hl.if_else(e.LA[k] == a_index, k, acc),
            hl.missing(hl.tint32),
            hl.range(hl.len(e.LA)),
        )

        def _with_lai(lai):
            non_ref_ad = hl.or_else(e.LAD[lai], 0)
            return e.annotate(
                GT=hl.if_else(
                    e.LGT.is_non_ref(),
                    hl.downcode(e.LGT, hl.or_else(lai, hl.len(e.LA))),
                    e.LGT,
                ),
                AD=hl.or_missing(
                    hl.is_defined(e.LAD), [hl.sum(e.LAD) - non_ref_ad, non_ref_ad]
                ),
            ).drop("LGT", "LAD", "LA")

        return hl.bind(_with_lai, lai)

    lt = lt.annotate(_split=splits).explode("_split")
    # Re-key as sparse_split_multi does so no sort is triggered.
    lt = lt._key_by_assert_sorted("locus")
    lt = lt.transmute(
        alleles=lt._split.alleles,
        carriers=lt.carriers.map(lambda e: _split_carrier(e, lt._split.a_index)),
    )
    lt = lt._key_by_assert_sorted("locus", "alleles")
    return _aggregate_carriers(lt, group_membership_ht, _carrier_hists)


def _region_read_intervals(
    test_region: list[str] | None, read_subintervals: int | None
) -> list[hl.utils.Interval] | None:
    """
    Turn ``--test-region`` strings into read intervals, tiled for parallelism.

    :param test_region: Region strings, or None.
    :param read_subintervals: Slices per region (see :func:`split_intervals`).
    :return: Read intervals, or None when no region was given.
    """
    if not test_region:
        return None
    intervals = [parse_region_interval(r) for r in test_region]
    return split_intervals(intervals, read_subintervals or 1)


def process_aou_dataset(
    environment: str,
    test: bool = False,
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
    read_intervals: list[hl.utils.Interval] | None = None,
) -> hl.Table:
    """
    Compute AoU per-stratum AC and homozygote counts and histograms over carriers.

    One pipeline for both the single job and each fan-out chunk; only the read
    scope differs. AN and AF are added by the merge from the all-sites AN release,
    because AN needs the reference blocks while these counts only need the carriers.

    :param environment: Compute environment.
    :param test: Use test paths.
    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading VDS partitions.
    :param chrom: Single contig to scope the reads to.
    :param read_intervals: Read-time intervals (a region or a layout chunk).
    :return: Table with ``freq`` and ``histograms``.
    """
    group_membership_ht = _group_membership_ht("aou", test, environment)
    aou_vds = get_aou_vds(
        annotate_meta=False,
        release_only=True,
        test=test_vds,
        filter_partitions=(
            list(range(test_partitions))
            if (test_partitions and not read_intervals)
            else None
        ),
        chrom=chrom,
        read_intervals=read_intervals,
        log_sample_counts=False,
        environment=environment,
    )
    vmt = _prepare_aou_vds(
        aou_vds,
        group_membership_ht,
        environment,
        skip_sex_ploidy=not spans_sex_chromosome(read_intervals, chrom),
    )
    ht = _sparse_split_strata_and_hists(vmt, group_membership_ht)

    # Cells back to the full strata, which the merge matches to the AN release.
    r = _read_reduction_globals(group_membership_ht.globals)
    ht = ht.select(
        freq=_int32_counts(
            expand_strata_array_from_leaves(
                ht.strata, r["leaf_indices"], r["decomposition"], r["n_full"]
            )
        ),
        histograms=hl.struct(
            qual_hists=ht.hist_fields.qual_hists,
            age_hists=ht.hist_fields.age_hists,
        ),
    )
    return select_final_dataset_fields(ht, dataset="aou")


# ---------------------------------------------------------------------------
# gnomAD
# ---------------------------------------------------------------------------


def _prepare_consent_mt(
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
) -> hl.MatrixTable:
    """
    Load and split the variant data of the samples that withdrew consent after v4.

    The counts subtracted from v4 come only from carriers, so the reference data is
    never read.

    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading partitions.
    :param chrom: Single contig to scope the load to.
    :return: Split variant MT of the withdrawn samples.
    """
    vds = get_gnomad_v5_genomes_vds(
        release_only=True,
        test=test_vds,
        consent_drop_only=True,
        annotate_meta=True,
        filter_partitions=(
            list(range(test_partitions))
            if (test_partitions is not None and chrom is None)
            else None
        ),
        chrom=chrom,
    )
    logger.info(
        "VDS filtered to %s consent-withdrawal samples.", vds.variant_data.count_cols()
    )
    vmt = vds.variant_data
    # Fewer, fuller partitions are cheaper for so few samples, and the QoB driver,
    # which collects one result per partition, died at 1,480 on chr22 (batch 8447491).
    vmt = vmt.naive_coalesce(max(1, vmt.n_partitions() // 15))
    meta = vmt.meta
    vmt = vmt.select_cols(
        sex_karyotype=meta.sex_imputation.sex_karyotype,
        # v3.1 fell back to the bin midpoint when the exact age was missing.
        age=hl.coalesce(meta.project_meta.age, meta.project_meta.age_alt),
    )
    # LPL stays so the split recalculates GQ from the downcoded PL, as it did for
    # v3.1; without it the site-level GQ is kept and adj changes at multi-allelic
    # sites.
    vmt = vmt.select_entries(
        "LA", "LAD", "LPL", "DP", "GQ", "LGT", _het_non_ref=vmt.LGT.is_het_non_ref()
    )
    vmt = vmt.checkpoint(new_temp_file("consent_samples_variant_data", "mt"))
    return hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)


def _consent_freq_ht(mt: hl.MatrixTable, test: bool, environment: str) -> hl.Table:
    """
    AC, homozygote counts and age histograms of the consent-withdrawal samples.

    All three come only from non-reference genotypes, which are all entries of the
    variant data, so nothing is densified. The strata and samples are the ones
    ``compute_coverage.py`` subtracted from the v4 AN.

    :param mt: Split variant MT from :func:`_prepare_consent_mt`.
    :param test: Use test paths.
    :param environment: Compute environment.
    :return: Table with ``freq`` (``AC`` / ``homozygote_count`` per stratum),
        ``age_hists`` and the strata globals.
    """
    # v3.1's order, which v4 inherited: adj from the diploid call, then sex ploidy,
    # then the hom-alt fix using the v3.0 release AF.
    # v4's counts used the pre-0.9.1 het-ref AB reading, so the subtraction must too.
    mt = annotate_adj(mt, legacy_het_ab=True)
    mt = adjust_sex_ploidy(mt, mt.sex_karyotype, xy_str="XY", xx_str="XX")
    v30 = v3_release_sites(public=True).versions["3.0"].ht()[mt.row_key]
    # Alleles no v3.0 release sample carried have no row in the public v3.0 release;
    # here they come only from the 21 withdrawn samples outside that release. The
    # internal v3.0 table v3.1's hom-alt fix read gave them AF 0, so use 0. It is
    # wrong only where no v3.0 release sample had an adj call (AF was missing); the
    # subtraction sets the resulting negatives to 0 and counts them.
    mt = mt.annotate_rows(_v30_af=hl.if_else(hl.is_defined(v30), v30.freq[0].AF, 0.0))
    mt = mt.annotate_entries(
        GT=hom_alt_depletion_fix(
            mt.GT,
            # Where the v3.0 AF is missing, v4 has v3.1.1's fix, which turned every het
            # with AB above 0.9 or missing into a missing call, het non-ref or not.
            het_non_ref_expr=mt._het_non_ref & hl.is_defined(mt._v30_af),
            af_expr=mt._v30_af,
            ab_expr=mt.AD[1] / mt.DP,
            use_v3_1_correction=True,
        ),
        # v4's age histograms are v3.1.1's, which the v3.1.2 het non-ref patch did not
        # recompute, so they were built with the fix before its het non-ref check.
        _age_GT=hom_alt_depletion_fix(
            mt.GT,
            het_non_ref_expr=hl.bool(False),
            af_expr=mt._v30_af,
            ab_expr=mt.AD[1] / mt.DP,
            use_v3_1_correction=True,
        ),
    )
    gm_ht = _group_membership_ht("gnomad", test, environment)
    ht = _aggregate_carriers(
        _localize_carriers(mt, gm_ht),
        gm_ht,
        lambda c: hl.struct(age_hists=age_hists_expr(c.adj, c._age_GT, c.age)),
    )
    gg = gm_ht.index_globals()
    ht = ht.select(
        freq=_int32_counts(ht.strata), age_hists=ht.hist_fields.age_hists
    ).select_globals(
        freq_meta=gg.freq_meta, freq_meta_sample_count=gg.freq_meta_sample_count
    )
    return ht.checkpoint(new_temp_file("consent_freq_and_hists", "ht"))


def _subtract_consent(v4_ht: hl.Table, consent_ht: hl.Table) -> hl.Table:
    """
    Subtract the consent-withdrawal counts and age histograms from the v4 table.

    Sites the withdrawn samples do not carry subtract nothing, so every v4 row
    passes through with its values intact. A count or age-histogram bin can only go
    negative where a withdrawn genotype v4 never counted is subtracted; it is set to
    0, as v4 did for the samples it subtracted. ``n_negative_set_to_zero`` counts the
    strata whose AC or homozygote count was set to 0, not the age-histogram bins.

    :param v4_ht: v4 release table with ``freq`` of ``AC`` / ``homozygote_count``
        structs.
    :param consent_ht: Table from :func:`_consent_freq_ht`.
    :return: v4 table with ``freq``, ``histograms.age_hists`` and the strata
        globals updated, and ``n_negative_set_to_zero`` added.
    """
    # The diff drops withdrawn strata v4 lacks without an error, so check first.
    v4_meta = hl.eval(v4_ht.freq_meta)
    missing = [m for m in hl.eval(consent_ht.freq_meta) if m not in v4_meta]
    if missing:
        raise ValueError(f"Withdrawn-sample strata missing from v4: {missing}")
    consent = consent_ht[v4_ht.key]
    ht = v4_ht.annotate(consent_freq=consent.freq, consent_age_hists=consent.age_hists)
    ht = ht.annotate_globals(
        consent_freq_meta=consent_ht.index_globals().freq_meta,
        consent_freq_meta_sample_count=consent_ht.index_globals().freq_meta_sample_count,
    )
    g = ht.index_globals()
    freq, freq_meta, counts = _merge_ac_homalt(
        [ht.freq, ht.consent_freq],
        [g.freq_meta, g.consent_freq_meta],
        [g.freq_meta_sample_count, g.consent_freq_meta_sample_count],
        operation="diff",
        set_negatives_to_zero=True,
    )
    # Strata where the withdrawn counts exceed v4's, so the diff above was set to 0.
    excess, _ = merge_array_expressions(
        [ht.consent_freq, ht.freq],
        [g.consent_freq_meta, g.freq_meta],
        operation="diff",
        set_negatives_to_zero=True,
        struct_fields=["AC", "homozygote_count"],
    )
    ht = ht.annotate(
        freq=freq,
        n_negative_set_to_zero=hl.if_else(
            hl.is_defined(ht.consent_freq),
            hl.len(excess.filter(lambda x: (x.AC > 0) | (x.homozygote_count > 0))),
            0,
        ),
        histograms=ht.histograms.annotate(
            age_hists=_merge_hist_struct(
                ht.histograms.age_hists,
                ht.consent_age_hists,
                operation="diff",
                set_negatives_to_zero=True,
            )
        ),
    ).drop("consent_freq", "consent_age_hists")
    return ht.annotate_globals(freq_meta=freq_meta, freq_meta_sample_count=counts)


def process_gnomad_dataset(
    environment: str,
    test: bool = False,
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
) -> hl.Table:
    """
    Build the gnomAD v5 AC, homozygote-count and histogram table from the v4 release.

    :param environment: Compute environment.
    :param test: Use test paths and restrict v4 to the consent sites.
    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading partitions.
    :param chrom: Single contig; the v4 table is scoped to it so a per-contig output
        holds only that contig.
    :return: Table with ``freq`` (``AC`` / ``homozygote_count`` per stratum) and
        ``histograms``.
    """
    v4_ht = release_sites(data_type="genomes").ht()
    # The HGDP/1KG subset strata have no all-sites AN to merge with AoU, so v5 drops
    # them, as the v4.1 joint release did.
    freq_meta, arrays = filter_arrays_by_meta(
        v4_ht.freq_meta,
        {
            "freq": v4_ht.freq,
            "freq_meta_sample_count": v4_ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=["subset"],
        keep=False,
    )
    v4_ht = v4_ht.select(
        freq=arrays["freq"].map(lambda f: f.select("AC", "homozygote_count")),
        # v4 copied these from the v3.1.2 release, so variants new in v4 have none.
        # Only the age histograms are adjusted for the withdrawn samples, below.
        histograms=v4_ht.histograms.drop("raw_qual_hists"),
    ).select_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=arrays["freq_meta_sample_count"],
    )
    if chrom:
        v4_ht = hl.filter_intervals(v4_ht, [hl.parse_locus_interval(chrom)])
    mt = _prepare_consent_mt(test_vds, test_partitions, chrom)
    consent_ht = _consent_freq_ht(mt, test, environment)
    if test:
        v4_ht = v4_ht.filter(hl.is_defined(consent_ht[v4_ht.key]))
        v4_ht = v4_ht.naive_coalesce(100).checkpoint(new_temp_file("v4_test", "ht"))
    ht = _subtract_consent(v4_ht, consent_ht)
    # Recomputed from the v5 meta, which already excludes the withdrawn samples,
    # rather than subtracted from v4's global.
    ht = ht.annotate_globals(
        age_distribution=_age_distribution(_release_meta("gnomad", environment))
    )
    return select_final_dataset_fields(ht, dataset="gnomad")


# ---------------------------------------------------------------------------
# Merge, FAF and final fields
# ---------------------------------------------------------------------------


def select_final_dataset_fields(ht: hl.Table, dataset: str = "gnomad") -> hl.Table:
    """
    Keep only the fields the next step reads.

    :param ht: Table to trim.
    :param dataset: ``"gnomad"``, ``"aou"`` or ``"merged"``.
    :return: Trimmed table.
    """
    if dataset not in ["gnomad", "aou", "merged"]:
        raise ValueError(f"Invalid dataset: {dataset}")
    if dataset in ["gnomad", "aou"]:
        final_globals = ["freq_meta", "freq_meta_sample_count", "age_distribution"]
        if dataset == "aou":
            # Only AoU has downsamplings; the v4 genomes release carries none.
            final_globals.append("downsamplings")
        # --process-gnomad logs it after the write; the merge drops it.
        fields = ["n_negative_set_to_zero"] if dataset == "gnomad" else []
        return ht.select("freq", "histograms", *fields).select_globals(*final_globals)
    ht = ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(ht.freq_meta),
        faf_index_dict=make_freq_index_dict_from_meta(ht.faf_meta),
    )
    return ht.select(
        "freq", "faf", "grpmax", "fafmax", "inbreeding_coeff", "histograms"
    ).select_globals(
        "freq_meta",
        "freq_index_dict",
        "freq_meta_sample_count",
        "faf_meta",
        "faf_index_dict",
        "age_distribution",
        "aou_downsamplings",
    )


def merge_gnomad_and_aou_frequencies(
    gnomad_freq_ht: hl.Table, aou_freq_ht: hl.Table, an_ht: hl.Table
) -> hl.Table:
    """
    Combine the gnomAD and AoU AC and homozygote counts into one ``freq`` array.

    The tables are outer-joined so a variant seen in only one dataset is kept. AC and
    homozygote counts of the strata the two tables share are summed, and the AoU
    strata are kept again as an ``aou`` subset so AoU can be reported on its own. AoU
    downsamplings only mean something within AoU, so they appear only in the subset.
    AN for every stratum comes from the all-sites AN release, so a variant absent from
    one dataset still counts that dataset's samples in AN. Histograms and the age
    distribution are summed too; a histogram one dataset lacks takes the other's.
    Variants with raw AC 0 in both datasets are dropped.

    :param gnomad_freq_ht: gnomAD AC, homozygote-count and histogram table.
    :param aou_freq_ht: AoU AC, homozygote-count and histogram table.
    :param an_ht: v5 all-sites AN release table, keyed by locus, built over the same
        strata and samples.
    :return: Merged table with ``freq``, ``histograms`` and the strata globals.
    """
    aou = aou_freq_ht.select(
        aou_freq=aou_freq_ht.freq, aou_histograms=aou_freq_ht.histograms
    ).select_globals()
    ht = gnomad_freq_ht.join(aou, how="outer")
    ag = aou_freq_ht.index_globals()
    ht = ht.annotate_globals(
        aou_freq_meta=ag.freq_meta,
        aou_freq_meta_sample_count=ag.freq_meta_sample_count,
        aou_age_distribution=ag.age_distribution,
        aou_downsamplings=ag.downsamplings,
    )
    g = ht.index_globals()
    merged_freq, merged_meta, counts = _merge_ac_homalt(
        [ht.freq, ht.aou_freq],
        [g.freq_meta, g.aou_freq_meta],
        [g.freq_meta_sample_count, g.aou_freq_meta_sample_count],
    )
    merged_meta = [dict(m) for m in merged_meta]
    global_idx = [i for i, d in enumerate(merged_meta) if "downsampling" not in d]
    aou_meta_list = hl.eval(aou_freq_ht.freq_meta)
    freq_meta = [merged_meta[i] for i in global_idx] + [
        {**d, "subset": "aou"} for d in aou_meta_list
    ]
    sample_counts = hl.eval(
        hl.array([counts[i] for i in global_idx]).extend(g.aou_freq_meta_sample_count)
    )
    an_meta, an_counts = hl.eval((an_ht.strata_meta, an_ht.strata_sample_count))
    an_meta = [dict(m) for m in an_meta]
    no_an = [m for m in freq_meta if m not in an_meta]
    if no_an:
        raise ValueError(f"Strata missing from the all-sites AN table: {no_an}")
    an_idx = [an_meta.index(m) for m in freq_meta]
    if sample_counts != [an_counts[i] for i in an_idx]:
        raise ValueError(
            "Freq and all-sites AN strata sample counts differ, so AC and AN were"
            " computed over different samples."
        )

    # A dataset with no row at a variant has no carriers there.
    no_aou_row = hl.range(len(aou_meta_list)).map(
        lambda _: hl.struct(AC=0, homozygote_count=0)
    )
    ac_hom = hl.array([merged_freq[i] for i in global_idx]).extend(
        hl.or_else(ht.aou_freq, no_aou_row)
    )
    an = an_ht[ht.locus].AN
    ht = ht.annotate(
        freq=hl.zip(ac_hom, hl.literal(an_idx)).map(
            lambda x: hl.bind(
                lambda an_i: hl.struct(
                    AC=x[0].AC,
                    AF=hl.or_missing(an_i > 0, x[0].AC / an_i),
                    AN=an_i,
                    homozygote_count=x[0].homozygote_count,
                ),
                hl.int32(an[x[1]]),
            )
        ),
        histograms=hl.struct(
            qual_hists=_merge_hist_struct(
                ht.histograms.qual_hists, ht.aou_histograms.qual_hists
            ),
            age_hists=_merge_hist_struct(
                ht.histograms.age_hists, ht.aou_histograms.age_hists
            ),
        ),
    )
    # Variants left with no raw carrier in either dataset (dead alleles, or carried
    # only by the withdrawn samples) are dropped, as v4's release filtered raw AC > 0.
    ht = ht.filter(ht.freq[freq_meta.index({"group": "raw"})].AC > 0)
    return ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=sample_counts,
        age_distribution=merge_histograms(
            [g.age_distribution, g.aou_age_distribution], operation="sum"
        ),
    )


def calculate_faf_and_grpmax_annotations(ht: hl.Table) -> hl.Table:
    """
    Add FAF, grpmax and fafmax overall and for the AoU subset, and the inbreeding coefficient.

    ``faf_expr`` and ``grpmax_expr`` select strata by exact key set, so the
    ``subset`` key is removed from the AoU meta for the call and put back on the
    returned ``faf_meta``.

    :param ht: Merged frequency table.
    :return: Table with ``faf``, ``grpmax``, ``fafmax``, ``inbreeding_coeff`` and
        the ``faf_meta`` global.
    """
    aou_meta, aou_arrays = filter_arrays_by_meta(
        ht.freq_meta,
        {"freq": ht.freq},
        items_to_filter={"subset": ["aou"]},
        keep=True,
        combine_operator="or",
    )
    aou_meta = aou_meta.map(
        lambda d: hl.dict(d.items().filter(lambda x: x[0] != "subset"))
    )
    freq_metas = {
        "gnomad": (ht.freq, ht.index_globals().freq_meta),
        "aou": (aou_arrays["freq"], aou_meta),
    }
    faf_exprs, faf_meta_exprs, grpmax_exprs, fafmax_exprs = [], [], {}, {}
    for dataset, (freq, freq_meta) in freq_metas.items():
        faf, faf_meta = faf_expr(
            freq, freq_meta, ht.locus, GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v5"]
        )
        grpmax_exprs[dataset] = grpmax_expr(
            freq, freq_meta, GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v5"]
        )
        fafmax_exprs[dataset] = gen_anc_faf_max_expr(faf, faf_meta)
        if dataset == "aou":
            faf_meta = [{**x, "subset": "aou"} for x in faf_meta]
        faf_exprs.append(faf)
        faf_meta_exprs.append(faf_meta)
    ht = ht.annotate(
        faf=hl.flatten(faf_exprs),
        grpmax=hl.struct(**grpmax_exprs),
        fafmax=hl.struct(**fafmax_exprs),
        inbreeding_coeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1]),
    )
    return ht.annotate_globals(faf_meta=hl.flatten(faf_meta_exprs))


# ---------------------------------------------------------------------------
# Hail Batch relay fan-out for AoU
#
# The relay machinery is batch_fanout.py, shared with compute_coverage.py; only
# freq policy lives here.
# ---------------------------------------------------------------------------


def _load_chunk_layout(
    an_environment: str, test: bool, n_chunks: int | None
) -> tuple[dict[str, Any], str]:
    """
    Load the layout coverage wrote beside the AN table, folded by ``--n-chunks``.

    Every orchestrator and worker must pass the same ``n_chunks`` so they agree on
    the chunk indices and the layout hash.

    :param an_environment: Environment (bucket) the AN table lives in.
    :param test: Read the test-scoped path.
    :param n_chunks: ``--n-chunks`` fold target.
    :return: ``(layout, layout hash)``.
    :raises FileNotFoundError: if the layout has not been written.
    """
    path = chunk_intervals_path(
        coverage_and_an_path(test=test, environment=an_environment).path
    )
    data = fold_chunk_layout(load_chunk_layout(path), n_chunks)
    logger.info("Using chunk layout %s (%d freq chunks).", path, len(data["chunks"]))
    return data, chunk_intervals_hash(data)


def _chunk_scope(args: argparse.Namespace) -> tuple[list[int], int, str]:
    """
    Return the chunks this run covers: eligible indices, total count and layout hash.

    Shared by the fan-out and the merge so they can never disagree on which chunks
    exist.

    :param args: Parsed CLI args.
    :return: ``(eligible indices, total chunks, layout hash)``.
    """
    if args.test_region:
        return [0], 1, test_region_hash(args.test_region)
    data, intervals_hash = _load_chunk_layout(
        args.an_environment, args.test, args.n_chunks
    )
    contigs = [c["contig"] for c in data["chunks"]]
    return select_chunks_for_contig(contigs, args.chrom), len(contigs), intervals_hash


def _relay_flags(args: argparse.Namespace, chunk: bool) -> str:
    """
    Return the flags every relay gets; chunk relays also get the read flags.

    :param args: Parsed CLI args.
    :param chunk: Include the chunk-only flags.
    :return: Space-joined flag string.
    """
    flags = [
        "--environment batch",
        f"--an-environment {args.an_environment}",
        f"--billing-project {args.billing_project}",
        f"--tmp-dir-days {args.tmp_dir_days}",
    ]
    flags += qob_sizing_flags(args)
    valued = [("app_name", "--app-name"), ("n_chunks", "--n-chunks")]
    if chunk:
        valued += [("read_subintervals", "--read-subintervals")]
    for attr, flag in valued:
        value = getattr(args, attr)
        if value is not None and value != "":
            flags.append(f"{flag} {value}")
    if chunk:
        if args.test_region:
            flags.append(f"--test-region {' '.join(args.test_region)}")
        if args.test_vds:
            flags.append("--test-vds")
        if args.test:
            flags.append("--test")
    return " ".join(flags)


def _relay_context(
    args: argparse.Namespace,
) -> tuple[str, str, str, Callable[[str, list[RelayJobSpec], str], int | None]]:
    """
    Return the commit, setup command, script and submitter every relay shares.

    The Hail version comes from the image, so no pin is passed.

    :param args: Parsed CLI args.
    :return: ``(commit, setup_cmd, script, submit)``. ``submit`` takes a batch name,
        job specs and a log label, and honors ``--batch-dry-run``.
    """
    commit = resolve_commit()
    setup_cmd, backend_kwargs, script = relay_context(
        commit,
        args.methods_branch,
        args.billing_project,
        args.batch_remote_tmpdir,
        "generate_frequency.py",
    )
    submit = partial(
        submit_relay_batch, args.batch_image, backend_kwargs, dry_run=args.batch_dry_run
    )
    return commit, setup_cmd, script, submit


def _orchestrate_freq_fanout(args: argparse.Namespace, ht_path: str) -> None:
    """
    Fan the AoU freq run out as one relay job per chunk.

    :param args: Parsed CLI args.
    :param ht_path: Final AoU freq HT path the chunks are namespaced under.
    :return: None.
    """
    eligible, n_total, intervals_hash = _chunk_scope(args)
    commit, setup_cmd, script, submit = _relay_context(args)
    common = _relay_flags(args, chunk=True)
    scope = "region" if args.test_region else f"{len(eligible)}of{n_total}c"

    def make_job(i: int) -> RelayJobSpec:
        return RelayJobSpec(
            name=f"freq_chunk_{i:06d}",
            cpu=args.chunk_cpu,
            memory=args.chunk_memory,
            storage=args.chunk_storage,
            attempts=1,
            command=(
                f"{setup_cmd}{script} --run-chunk --chunk-index {i}"
                f" --chunk-output {chunk_path(ht_path, i, intervals_hash)} {common}"
            ),
        )

    run_chunk_fanout(
        ht_path,
        intervals_hash,
        eligible,
        make_job=make_job,
        submit=submit,
        batch_base=combine_suffix(
            f"v5_freq_aou_chunk_{scope}_{intervals_hash[:8]}",
            args.aou_freq_ht_suffix,
            args.chrom,
        ),
        wave_size=args.wave_size,
        commit=commit,
        app_name=args.app_name,
        overwrite=args.overwrite,
        dry_run=args.batch_dry_run,
    )


def _run_freq_chunk(args: argparse.Namespace) -> None:
    """
    Compute and write one chunk's AoU freq HT, inside a relay under a nested QoB driver.

    The read intervals come from the layout JSON or the region strings, so the
    worker never opens the VDS to find them.

    :param args: Parsed CLI args (``chunk_index``, ``chunk_output``).
    :return: None.
    """
    _initialize_hail(args)
    if args.test_region:
        intervals_hash = test_region_hash(args.test_region)
        sub_intervals = _region_read_intervals(args.test_region, args.read_subintervals)
    else:
        data, intervals_hash = _load_chunk_layout(
            args.an_environment, args.test, args.n_chunks
        )
        sub_intervals = chunk_read_intervals(data, args.chunk_index)
    logger.info("Chunk %d: %d read intervals.", args.chunk_index, len(sub_intervals))
    ht = process_aou_dataset(
        environment=args.environment,
        test=args.test,
        test_vds=args.test_vds,
        read_intervals=sub_intervals,
    )
    ht = ht.annotate_globals(freq_chunk_intervals_hash=intervals_hash)
    ht.write(args.chunk_output, overwrite=True)
    logger.info("Wrote chunk %d to %s", args.chunk_index, args.chunk_output)


def _orchestrate_freq_merge(args: argparse.Namespace, ht_path: str) -> None:
    """
    Union the chunk HTs into the final AoU freq HT with a tree of relay merge jobs.

    :param args: Parsed CLI args.
    :param ht_path: Final AoU freq HT path.
    :return: None.
    """
    eligible, _n_total, intervals_hash = _chunk_scope(args)
    present = list_present_chunk_indices(ht_path, intervals_hash)
    missing = [i for i in eligible if i not in present]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} of {len(eligible)} chunks are missing (first:"
            f" {missing[:5]}); rerun --use-batch-fanout first."
        )
    _commit, setup_cmd, script, submit = _relay_context(args)
    gs = args.merge_group_size
    orchestrate_tree_merge(
        inputs=[chunk_path(ht_path, i, intervals_hash) for i in eligible],
        final_output_path=ht_path,
        group_path_fn=lambda level, idx: group_path(
            ht_path, level, idx, f"gs{gs}", intervals_hash
        ),
        merge_group_size=gs,
        overwrite=args.overwrite,
        final_coalesce_to=args.n_partitions,
        batch_prefix=combine_suffix(
            "v5_freq_merge", args.aou_freq_ht_suffix, args.chrom
        ),
        job_prefix="freq_merge",
        setup_cmd=setup_cmd,
        script=script,
        common_flags=_relay_flags(args, chunk=False),
        merge_cpu=args.merge_cpu,
        merge_memory=args.merge_memory,
        merge_storage=args.merge_storage,
        final_merge_storage=args.final_merge_storage,
        submit=submit,
    )


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------


def _initialize_hail(args: argparse.Namespace) -> None:
    """
    Start Hail the same way for every role that runs Hail in-process.

    A script-specific tmp directory keeps freq and coverage temp files apart.

    :param args: Parsed CLI args.
    :return: None.
    """
    tmp = qc_temp_prefix(environment=args.environment, days=args.tmp_dir_days)
    _init_hail(
        "v5_frequency_generation",
        args.environment,
        batch_billing_project=args.billing_project,
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{tmp}frequency_generation",
        **_get_batch_resource_kwargs(args),
    )


def main(args):
    """Generate v5 frequency data."""
    args.an_environment = args.an_environment or (
        args.environment if args.test else "dataproc"
    )
    if args.environment == "batch" and args.app_name is None:
        steps = [
            s
            for s, on in (
                ("gnomad", args.process_gnomad),
                (
                    "aou",
                    args.process_aou or args.use_batch_fanout or args.merge_freq_chunks,
                ),
                ("merge", args.merge_datasets),
            )
            if on
        ]
        args.app_name = f"v5_freq_{'_'.join(steps)}" if steps else None

    if args.run_chunk:
        _run_freq_chunk(args)
        return
    if args.run_merge:
        _initialize_hail(args)
        union_and_write_hts(
            args.merge_inputs, args.merge_output, args.merge_coalesce_to
        )
        return

    aou_suffix = combine_suffix(args.aou_freq_ht_suffix, args.chrom)
    aou_out_suffix = combine_suffix(aou_suffix, args.freq_output_suffix)
    gnomad_out_suffix = combine_suffix(args.chrom, args.freq_output_suffix)
    freq_kwargs = {
        "test": args.test,
        "data_type": "genomes",
        "environment": args.environment,
    }

    # Orchestrators submit relays and never start Hail here.
    if args.use_batch_fanout or args.merge_freq_chunks:
        aou_freq = get_ac_homalt_and_hists(
            data_set="aou", suffix=aou_out_suffix, **freq_kwargs
        )
        if args.use_batch_fanout:
            # Not gated on the final HT: the fan-out never writes it, and a rerun
            # after a merge must stay possible.
            _orchestrate_freq_fanout(args, aou_freq.path)
            return
        check_resource_existence(
            output_step_resources={"merge-freq-chunks": [aou_freq]},
            overwrite=args.overwrite,
        )
        _orchestrate_freq_merge(args, aou_freq.path)
        return

    _initialize_hail(args)
    try:
        if args.assemble_chrom_freq:
            data_set = args.assemble_data_set
            base_suffix = args.aou_freq_ht_suffix if data_set == "aou" else None
            canonical = get_ac_homalt_and_hists(
                data_set=data_set, suffix=base_suffix, **freq_kwargs
            )
            per_contig = [
                get_ac_homalt_and_hists(
                    data_set=data_set,
                    suffix=combine_suffix(base_suffix, c),
                    **freq_kwargs,
                )
                for c in args.contigs
            ]
            check_resource_existence(
                input_step_resources={"per-contig-freq": per_contig},
                output_step_resources={"assemble-chrom-freq": [canonical]},
                overwrite=args.overwrite,
            )
            union_and_write_hts(
                [p.path for p in per_contig], canonical.path, args.n_partitions
            )

        if args.process_gnomad:
            gnomad_freq = get_ac_homalt_and_hists(
                data_set="gnomad", suffix=gnomad_out_suffix, **freq_kwargs
            )
            check_resource_existence(
                output_step_resources={"process-gnomad": [gnomad_freq]},
                overwrite=args.overwrite,
            )
            ht = process_gnomad_dataset(
                environment=args.environment,
                test=args.test,
                test_vds=args.test_vds,
                test_partitions=args.test_partitions,
                chrom=args.chrom,
            )
            ht.write(gnomad_freq.path, overwrite=args.overwrite)
            ht = gnomad_freq.ht()
            n_neg = ht.aggregate(hl.agg.sum(ht.n_negative_set_to_zero))
            if n_neg > args.max_negative_set_to_zero:
                raise ValueError(
                    f"{n_neg} variant strata had withdrawn counts above v4's and were"
                    f" set to 0, more than --max-negative-set-to-zero"
                    f" ({args.max_negative_set_to_zero}); check the strata and"
                    f" samples before using {gnomad_freq.path}."
                )
            if n_neg:
                logger.warning(
                    "%d variant strata had withdrawn counts above v4's and were set"
                    " to 0.",
                    n_neg,
                )

        if args.process_aou:
            aou_freq = get_ac_homalt_and_hists(
                data_set="aou", suffix=aou_out_suffix, **freq_kwargs
            )
            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq]},
                overwrite=args.overwrite,
            )
            ht = process_aou_dataset(
                environment=args.environment,
                test=args.test,
                test_vds=args.test_vds,
                test_partitions=args.test_partitions,
                chrom=args.chrom,
                read_intervals=_region_read_intervals(
                    args.test_region, args.read_subintervals
                ),
            )
            ht.write(aou_freq.path, overwrite=args.overwrite)

        if args.merge_datasets:
            merged = get_freq(
                data_set="merged", suffix=args.freq_output_suffix, **freq_kwargs
            )
            check_resource_existence(
                output_step_resources={"merge-datasets": [merged]},
                overwrite=args.overwrite,
            )
            gnomad_freq = get_ac_homalt_and_hists(data_set="gnomad", **freq_kwargs)
            aou_freq = get_ac_homalt_and_hists(
                data_set="aou", suffix=args.aou_freq_ht_suffix, **freq_kwargs
            )
            an_path = release_coverage_path(
                test=args.test,
                coverage_type="allele_number",
                environment=args.an_environment,
            )
            check_resource_existence(
                input_step_resources={
                    "process-gnomad": [gnomad_freq],
                    "process-aou": [aou_freq],
                    "export-an-release-files": [an_path],
                }
            )
            ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq.ht(), aou_freq.ht(), hl.read_table(an_path)
            )
            ht = ht.checkpoint(new_temp_file("merged_freq", "ht"))
            n_no_an = ht.aggregate(hl.agg.count_where(hl.is_missing(ht.freq[0].AN)))
            if n_no_an:
                logger.warning(
                    "%d merged variants have no all-sites AN at their locus; their AN"
                    " and AF are missing.",
                    n_no_an,
                )
            ht = calculate_faf_and_grpmax_annotations(ht)
            ht = select_final_dataset_fields(ht, dataset="merged")
            ht.write(merged.path, overwrite=args.overwrite)
    finally:
        # On the batch backend the driver's log lives in Batch, not locally.
        if args.environment != "batch":
            hl.copy_log(
                get_logging_path(
                    "v5_frequency_run",
                    environment=args.environment,
                    tmp_dir_days=args.tmp_dir_days,
                )
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate frequency data for gnomAD v5."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing Hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--aou-freq-ht-suffix",
        default=None,
        help=(
            "Tag for the AoU freq HT so an experimental run does not overwrite the"
            " default; applied to the write and to the read in --merge-datasets."
        ),
    )
    parser.add_argument(
        "--freq-output-suffix",
        default=None,
        help=(
            "Tag the HT this run writes without changing any input it reads, so a"
            " trial output never feeds --merge-datasets."
        ),
    )

    steps = parser.add_argument_group("processing steps")
    steps.add_argument(
        "--process-gnomad",
        action="store_true",
        help="gnomAD v5 frequencies: v4 minus the samples that withdrew consent.",
    )
    steps.add_argument(
        "--max-negative-set-to-zero",
        type=int,
        default=1000,
        help=(
            "Fail --process-gnomad after its write if more variant strata than this"
            " went negative and were set to 0. A few are expected where v4 never"
            " counted a withdrawn genotype; many point to a strata or sample"
            " mismatch. Default 1000."
        ),
    )
    steps.add_argument(
        "--process-aou",
        action="store_true",
        help=(
            "AoU as one job: AC, homozygote counts and histograms from the VDS variant"
            " data. Use --use-batch-fanout for the genome."
        ),
    )
    steps.add_argument(
        "--use-batch-fanout",
        action="store_true",
        help=(
            "Run the AoU frequencies as one Batch job per coverage chunk so a single"
            " crash does not lose the run; combine the chunk HTs with"
            " --merge-freq-chunks. Rerunning resumes from the chunks that are missing."
        ),
    )
    steps.add_argument(
        "--merge-freq-chunks",
        action="store_true",
        help="Union the fan-out's chunk HTs into the final AoU freq HT.",
    )
    steps.add_argument(
        "--merge-datasets",
        action="store_true",
        help=(
            "Combine the gnomAD and AoU tables and add FAF, grpmax and inbreeding. AN"
            " comes from the compute_coverage all-sites AN release."
        ),
    )
    steps.add_argument(
        "--chrom",
        default=None,
        help=(
            "Run one contig at a time so a failure costs only that contig; each contig"
            " writes its own HT, joined later with --assemble-chrom-freq."
        ),
    )
    steps.add_argument(
        "--assemble-chrom-freq",
        action="store_true",
        help="Union the per-contig HTs of --assemble-data-set into the canonical HT.",
    )
    steps.add_argument(
        "--assemble-data-set",
        choices=["aou", "gnomad"],
        default="aou",
        help="Data set whose per-contig HTs --assemble-chrom-freq unions. Default aou.",
    )
    steps.add_argument(
        "--contigs",
        nargs="+",
        default=None,
        help="Contigs to union with --assemble-chrom-freq.",
    )

    test_group = parser.add_argument_group("testing options")
    test_group.add_argument(
        "--test-vds",
        action="store_true",
        help="Read the test VDS instead of the full one and turn on test paths.",
    )
    test_group.add_argument(
        "--test-partitions",
        type=int,
        default=None,
        help="Run on the first N VDS partitions and turn on test paths.",
    )
    test_group.add_argument(
        "--test-region",
        nargs="+",
        default=None,
        help=(
            "Restrict the AoU run to these intervals (chr:start-end, half-open like"
            " compute_coverage --test-region). Turns on test paths."
        ),
    )
    test_group.add_argument(
        "--read-subintervals",
        type=int,
        default=48,
        help=(
            "Split each --test-region into this many read intervals so the aggregation"
            " runs in parallel instead of in one partition. Default 48; 1 disables."
        ),
    )
    test_group.add_argument("--test", action="store_true", help=argparse.SUPPRESS)

    env_group = parser.add_argument_group("environment configuration")
    env_group.add_argument(
        "--environment",
        choices=["rwb", "batch"],
        default="batch",
        help="Environment to run in.",
    )
    env_group.add_argument(
        "--an-environment",
        choices=["rwb", "batch", "dataproc"],
        default=None,
        help=(
            "Bucket the AN tables and the chunk layout are read from. Default:"
            " dataproc, where prod compute_coverage writes them, or --environment"
            " with --test."
        ),
    )
    env_group.add_argument(
        "--n-chunks",
        type=int,
        default=None,
        help=(
            "Fold neighbouring coverage chunks into about this many freq jobs to cut"
            " per-job overhead; results are unchanged. Default: one job per coverage"
            " chunk."
        ),
    )
    env_group.add_argument(
        "--tmp-dir-days",
        type=int,
        default=4,
        help="Temp directory retention in days. Default 4.",
    )
    env_group.add_argument(
        "--billing-project",
        default="gnomad-production",
        help="Hail Batch billing project for the QoB driver and every relay job.",
    )

    batch_group = parser.add_argument_group(
        "batch configuration", "QoB sizing for this process (--environment batch only)."
    )
    batch_group.add_argument("--app-name", default=None, help="Hail Batch name.")
    batch_group.add_argument(
        "--driver-cores", type=int, default=None, help="QoB driver cores."
    )
    batch_group.add_argument(
        "--driver-memory", default=None, help="QoB driver memory class, e.g. highmem."
    )
    batch_group.add_argument(
        "--worker-cores",
        default=None,
        help="Cores per QoB worker; Hail Batch accepts 1, 2, 4 or 8 for JVM jobs.",
    )
    batch_group.add_argument(
        "--worker-memory", default=None, help="QoB worker memory class, e.g. highmem."
    )

    fanout = parser.add_argument_group(
        "batch fan-out configuration",
        "Relay settings for --use-batch-fanout and --merge-freq-chunks.",
    )
    fanout.add_argument(
        "--wave-size",
        type=int,
        default=1000,
        help="Chunks per sequential wave, bounding concurrent relays. <=0 means one wave. Default 1000.",
    )
    fanout.add_argument(
        "--chunk-driver-cores",
        type=int,
        default=1,
        help="Cores for each chunk's QoB driver; it only builds the query, so 1 is enough.",
    )
    fanout.add_argument(
        "--chunk-driver-memory",
        default="highmem",
        help="Memory class for each chunk's QoB driver; highmem is the smallest that compiles the query.",
    )
    fanout.add_argument(
        "--chunk-worker-cores",
        default=None,
        help="Cores per chunk QoB worker; accepts the same values as --worker-cores.",
    )
    fanout.add_argument(
        "--chunk-worker-memory", default=None, help="Memory class per chunk QoB worker."
    )
    fanout.add_argument(
        "--chunk-cpu",
        type=float,
        default=0.5,
        help="CPU per chunk relay; it only waits on its QoB batch. Default 0.5.",
    )
    fanout.add_argument(
        "--chunk-memory",
        default="standard",
        help="Memory class per chunk relay. Default standard.",
    )
    fanout.add_argument(
        "--chunk-storage", default="25Gi", help="Storage per chunk relay. Default 25Gi."
    )
    fanout.add_argument(
        "--merge-group-size",
        type=int,
        default=500,
        help="Chunk HTs unioned per merge job. Default 500.",
    )
    fanout.add_argument(
        "--n-partitions",
        type=int,
        default=None,
        help="Partitions of the final merged or assembled HT. Default: no coalesce.",
    )
    fanout.add_argument(
        "--merge-cpu", type=int, default=4, help="CPU per merge job. Default 4."
    )
    fanout.add_argument(
        "--merge-memory",
        default="standard",
        help="Memory class per merge job. Default standard.",
    )
    fanout.add_argument(
        "--merge-storage",
        default="50Gi",
        help="Storage per intermediate merge job. Default 50Gi.",
    )
    fanout.add_argument(
        "--final-merge-storage",
        default="100Gi",
        help="Storage for the final merge job. Default 100Gi.",
    )
    fanout.add_argument(
        "--batch-image",
        default=DEFAULT_BATCH_IMAGE,
        help="Relay image; its Hail version is the one the run uses.",
    )
    fanout.add_argument(
        "--batch-remote-tmpdir",
        default=None,
        help="gs:// scratch for the Batch ServiceBackend. Default: Hail's configured one.",
    )
    fanout.add_argument(
        "--methods-branch",
        default="main",
        help="gnomad_methods branch or commit the relays pull. Default main.",
    )
    fanout.add_argument(
        "--batch-dry-run",
        action="store_true",
        help="Validate the first wave's DAG without running it.",
    )

    worker = parser.add_argument_group(
        "batch worker subcommands",
        "Internal; passed by the orchestrator to its relays.",
    )
    worker.add_argument("--run-chunk", action="store_true", help=argparse.SUPPRESS)
    worker.add_argument("--chunk-index", type=int, default=None, help=argparse.SUPPRESS)
    worker.add_argument("--chunk-output", default=None, help=argparse.SUPPRESS)
    worker.add_argument("--run-merge", action="store_true", help=argparse.SUPPRESS)
    worker.add_argument(
        "--merge-inputs", nargs="+", default=None, help=argparse.SUPPRESS
    )
    worker.add_argument("--merge-output", default=None, help=argparse.SUPPRESS)
    worker.add_argument(
        "--merge-coalesce-to", type=int, default=None, help=argparse.SUPPRESS
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    args.test = (
        args.test
        or args.test_vds
        or args.test_partitions is not None
        or args.test_region is not None
    )
    provided = [a for a in _BATCH_RESOURCE_PARAMS if getattr(args, a, None) is not None]
    if provided and args.environment != "batch":
        parser.error(
            "Batch arguments ("
            + ", ".join("--" + a.replace("_", "-") for a in provided)
            + ") require --environment batch."
        )
    if args.test_region is not None and args.test_partitions is not None:
        parser.error("--test-region and --test-partitions are mutually exclusive.")
    if args.test_region and args.chrom:
        parser.error("--test-region already scopes the contig; do not pass --chrom.")
    if args.assemble_chrom_freq:
        if not args.contigs:
            parser.error("--assemble-chrom-freq requires --contigs.")
        if args.chrom:
            parser.error(
                "--assemble-chrom-freq writes the contig-unscoped HT; do not pass --chrom."
            )
    if args.run_chunk and (args.chunk_index is None or args.chunk_output is None):
        parser.error("--run-chunk requires --chunk-index and --chunk-output.")
    if args.run_merge and not (args.merge_inputs and args.merge_output):
        parser.error("--run-merge requires --merge-inputs and --merge-output.")

    main(args)
