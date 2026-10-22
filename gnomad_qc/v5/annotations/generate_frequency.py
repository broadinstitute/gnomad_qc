"""
Generate gnomAD v5 genome frequencies.

gnomAD v5 adds no new gnomAD samples, so ``--process-gnomad`` reuses the v4
frequencies: it computes frequencies and age histograms for the samples that
withdrew consent after v4 and subtracts them. All of Us (AoU) has nothing to reuse,
so ``--process-aou`` aggregates AC, homozygote counts and histograms from the AoU
VDS variant data and joins the allele number that ``compute_coverage.py`` already
computed at every site. ``--merge-datasets`` combines the two tables and adds FAF,
grpmax and the inbreeding coefficient, which only make sense on the combined cohort.

The AoU run is too large for one Query-on-Batch job, so it fans out over the chunk
layout ``compute_coverage.py`` wrote beside its AN table (``--use-batch-fanout``) and
tree-merges the chunk tables afterwards (``--merge-freq-chunks``). Reusing coverage's
layout means every freq chunk reads exactly the AN partitions written for it.

Examples::

    # AoU: fan out, then merge.
    generate_frequency.py --use-batch-fanout --environment batch
    generate_frequency.py --merge-freq-chunks --environment batch
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
    compute_freq_by_strata,
    expand_strata_array_from_leaves,
    faf_expr,
    gen_anc_faf_max_expr,
    grpmax_expr,
    merge_freq_arrays,
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

# force=True so this handler wins over the ones Hail installs at import time;
# without it every logger.info call is dropped when running against the batch backend.
logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("v5_frequency")

# Hail 0.2.137 is the lowest version that runs on the Batch workers' Java 21.
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

    The meta table has ~330 partitions, so every join or aggregate against it fans
    ~330 tiny QoB tasks; a lazy select-and-coalesce of the two fields freq needs
    makes that 10, with nothing written.

    :param project: ``"aou"`` or ``"gnomad"``.
    :param environment: Compute environment.
    :return: Table keyed by ``s`` with ``sex_karyotype`` and ``age``.
    """
    ht = meta(data_type="genomes", environment=environment).ht()
    ht = ht.filter((ht.project_meta.project == project) & ht.release)
    return ht.select("sex_karyotype", age=ht.project_meta.age).naive_coalesce(10)


def _age_distribution(ht: hl.Table) -> hl.Struct:
    """
    Return the age histogram global over a release-meta table.

    Computed from the meta table rather than from the variant MT's columns so it
    never forces an extra pass over the variant data.

    :param ht: Table from :func:`_release_meta`.
    :return: Age histogram struct.
    """
    return ht.aggregate(hl.agg.hist(ht.age, *AGE_HIST_BINS))


def _aou_group_membership_ht(test: bool, environment: str) -> hl.Table:
    """
    Read the cell-reduced AoU group-membership HT ``compute_coverage.py`` built the AN with.

    Freq's AC must be stratified over exactly the sample sets the AN was, so the join
    lines up index for index; reading any other membership table would be wrong, so
    a missing HT is an error rather than a fallback.

    :param test: Read the test-scoped path.
    :param environment: Compute environment.
    :return: The cells HT (one boolean per cell; ``_full`` globals map cells back to
        the full strata).
    :raises FileNotFoundError: if the cells HT has not been written.
    """
    path = apply_path_suffix(
        group_membership(test=test, data_set="aou", environment=environment).path,
        "cells",
    )
    if not _file_exists_for_env(path, environment):
        raise FileNotFoundError(
            f"AoU group-membership cells HT not found at {path}; run"
            " compute_coverage.py --write-group-membership-ht first."
        )
    return hl.read_table(path)


def _carrier_hists(c: hl.expr.StructExpression) -> hl.expr.StructExpression:
    """
    Quality and age histogram aggregators over the carriers of one split row.

    Only the adj quality histograms are kept in v5, so they are filtered at the
    source instead of also building the raw ones.

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
) -> hl.expr.StructExpression:
    """
    Merge every histogram field of two histogram structs.

    :param hist1: First struct of histograms.
    :param hist2: Second struct with the same fields.
    :param operation: ``"sum"`` or ``"diff"``.
    :return: Struct of merged histograms.
    """
    return hl.struct(
        **{
            field: merge_histograms([hist1[field], hist2[field]], operation=operation)
            for field in hist1.dtype.fields
        }
    )


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

    Rows are left unsplit on purpose: adj is computed on the local ``LGT``/``LAD``
    fields, and :func:`_sparse_split_strata_and_hists` splits each row only after
    compacting its entries to the carriers, which is what keeps the per-row cost
    proportional to carriers rather than samples.

    :param aou_vds: AoU VDS filtered to release samples.
    :param group_membership_ht: The cells HT from :func:`_aou_group_membership_ht`.
    :param environment: Compute environment.
    :param skip_sex_ploidy: Skip the sex-ploidy adjustment; only valid when the
        read scope has no chrX/chrY loci.
    :return: Unsplit variant MT with ``LGT``, ``GQ``, ``LAD``, ``LA``, ``adj``
        entries, ``age`` on the columns, and the strata globals.
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

    The entries array has one slot per sample (~365k) but the median variant has
    three carriers, so ``split_multi`` followed by ``agg_by_strata`` spends almost all
    of its time on missing entries, once per split row and per stratum. Compacting
    to carriers first and splitting afterwards makes every later step proportional
    to carriers. Output is identical to ``hl.vds.split_multi(filter_changed_loci=True)``
    plus ``agg_by_strata``: the split reproduces ``sparse_split_multi``'s allele
    min-rep, row ordering and ``LGT``/``LAD`` downcoding, applied to bi-allelic rows
    as well because dead-allele removal can leave a bi-allelic row non-minimal.

    :param mt: Prepared, unsplit variant MT from :func:`_prepare_aou_vds`.
    :param group_membership_ht: Membership HT whose ``group_membership`` array and
        ``freq_meta`` global define the strata (cells).
    :return: Table keyed by the split ``locus, alleles`` with ``hist_fields`` and
        ``strata`` (array of ``AC`` / ``homozygote_count`` structs, one per stratum)
        and ``mt``'s globals.
    """
    freq_meta = [dict(m) for m in hl.eval(group_membership_ht.freq_meta)]
    n_groups = len(freq_meta)
    adj_groups = hl.literal([m.get("group", "NA") == "adj" for m in freq_meta])
    gm = group_membership_ht.select(
        strata=hl.enumerate(group_membership_ht.group_membership)
        .filter(lambda t: t[1])
        .map(lambda t: t[0])
    )
    mt = mt.annotate_cols(strata=gm[mt.col_key].strata)
    lt = mt.localize_entries("entries", "cols")
    # Keep only carriers. Each carries the strata it counts toward (adj strata only
    # when the call is adj) so one explode+group_by does every stratum at once.
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
    lt = lt.select(carriers=carriers)

    # Same split structs as sparse_split_multi; an alt whose min_rep moves the locus
    # is dropped (filter_changed_loci).
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

    # Same entry transform as sparse_split_multi for the fields present (LGT, LAD,
    # LA, GQ).
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
    # Re-key exactly as sparse_split_multi does so no sort is triggered.
    lt = lt._key_by_assert_sorted("locus")
    lt = lt.transmute(
        alleles=lt._split.alleles,
        carriers=lt.carriers.map(lambda e: _split_carrier(e, lt._split.a_index)),
    )
    lt = lt._key_by_assert_sorted("locus", "alleles")

    empty = hl.struct(AC=hl.int64(0), homozygote_count=hl.int64(0))
    lt = lt.annotate(
        hist_fields=lt.carriers.aggregate(_carrier_hists),
        _by_stratum=lt.carriers.aggregate(
            lambda c: hl.agg.explode(
                lambda g: hl.agg.group_by(
                    g,
                    hl.struct(
                        AC=hl.agg.sum(c.GT.n_alt_alleles()),
                        homozygote_count=hl.agg.count_where(c.GT.is_hom_var()),
                    ),
                ),
                c.counted_in,
            )
        ),
    )
    lt = lt.annotate(
        strata=hl.range(n_groups).map(lambda g: lt._by_stratum.get(g, empty))
    )
    return lt.select("hist_fields", "strata").drop("cols")


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
    an_environment: str,
    test: bool = False,
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
    read_intervals: list[hl.utils.Interval] | None = None,
    all_sites_an_suffix: str | None = None,
) -> hl.Table:
    """
    Compute AoU frequencies, with AN joined from the all-sites-AN table.

    One pipeline for both the single job and each fan-out chunk; only the read
    scope differs. AN comes from ``compute_coverage.py`` because it needs the
    reference blocks, while freq only needs the carriers.

    :param environment: Compute environment.
    :param an_environment: Environment (bucket) the AN table is read from.
    :param test: Use test paths.
    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading VDS partitions.
    :param chrom: Single contig to scope the reads to.
    :param read_intervals: Read-time intervals (a region or a layout chunk); these
        also scope the AN read, so a region test matches an AN table made over the
        same region.
    :param all_sites_an_suffix: Suffix of the AN table written with
        ``compute_coverage.py --cov-and-an-output-suffix``.
    :return: Table with ``freq`` and ``histograms``.
    """
    group_membership_ht = _aou_group_membership_ht(test, environment)
    # Read-time interval pruning rather than a post-read filter: the latter leaves a
    # small region fanned across thousands of empty partitions.
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

    # Cells back to the full strata, which is what the AN array is indexed by.
    r = _read_reduction_globals(group_membership_ht.globals)
    ht = ht.annotate(
        strata=expand_strata_array_from_leaves(
            ht.strata, r["leaf_indices"], r["decomposition"], r["n_full"]
        )
    )
    # Checkpoint before the AN join: one fused query from split to write compiles
    # to more JVM bytecode than a standard driver can hold.
    ht = ht.checkpoint(new_temp_file("aou_freq_agg", "ht"))

    an_path = apply_path_suffix(
        coverage_and_an_path(test=test, environment=an_environment).path,
        all_sites_an_suffix,
    )
    logger.info("Joining AN from %s...", an_path)
    an_ht = hl.read_table(an_path)
    if read_intervals:
        an_ht = hl.filter_intervals(an_ht, read_intervals)
    elif chrom:
        an_ht = hl.filter_intervals(an_ht, [hl.parse_locus_interval(chrom)])
    ht = ht.annotate(_an=an_ht[ht.locus].AN)
    # A variant with no AN at its locus gets a missing freq array below, which
    # nothing downstream would flag. The AN table should cover every AoU variant
    # locus (checked 2026-10-02), so any count here means the two runs disagree.
    n_missing_an = ht.aggregate(hl.agg.count_where(hl.is_missing(ht._an)))
    if n_missing_an:
        logger.warning(
            "%d variant rows have no AN at their locus (outside the all-sites AN"
            " table); their freq array will be missing.",
            n_missing_an,
        )
    ht = ht.select(
        freq=hl.map(
            lambda s, an: hl.struct(
                AC=hl.int32(s.AC),
                AF=hl.if_else(an > 0, s.AC / an, hl.missing(hl.tfloat64)),
                AN=hl.int32(an),
                homozygote_count=hl.int32(s.homozygote_count),
            ),
            ht.strata,
            ht._an,
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


def _prepare_consent_vds(
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
) -> hl.vds.VariantDataset:
    """
    Load and split the VDS of the samples that withdrew consent after v4.

    Only the load, the column and entry narrowing and the split happen here. adj,
    sex ploidy and the hom-alt fix are applied on the dense matrix in
    :func:`_consent_freq_ht`, in v3.1's order, so reference blocks get the same
    treatment they had when v4's counts were made.

    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading partitions.
    :param chrom: Single contig to scope the load to.
    :return: Split VDS of the withdrawn samples.
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
    # A few hundred samples do not need the full VDS's partitioning (1,480 on chr22).
    # Fewer, fuller partitions make the extraction cheaper, and the QoB driver, which
    # collects one result per partition, died at 1,480 (batch 8447491).
    n_parts = max(1, vds.variant_data.n_partitions() // 15)
    vds = hl.vds.VariantDataset(
        vds.reference_data.naive_coalesce(n_parts),
        vds.variant_data.naive_coalesce(n_parts),
    )
    vmt = vds.variant_data
    meta = vmt.meta
    vmt = vmt.select_cols(
        gen_anc=meta.population_inference.pop,
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
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = vds.checkpoint(new_temp_file("consent_samples_vds", "vds"))
    return hl.vds.split_multi(vds, filter_changed_loci=True)


def _consent_freq_ht(
    vds: hl.vds.VariantDataset, test: bool, environment: str
) -> hl.Table:
    """
    Frequencies and age histograms of the consent-withdrawal samples.

    Densifying is affordable because the VDS holds only the withdrawn samples. The
    strata come from the gnomAD membership table, so the subtracted arrays line up
    with v4's.

    :param vds: Split VDS from :func:`_prepare_consent_vds`.
    :param test: Use test paths.
    :param environment: Compute environment.
    :return: Table with ``freq`` and ``age_hists``.
    """
    mt = hl.vds.to_dense_mt(vds)
    # v3.1's order, which v4 inherited: adj from the diploid call (reference blocks
    # included, from the block's GQ and DP), then sex ploidy, then the hom-alt fix
    # using the v3.0 release AF. The hom-alt fix only changes hets, so reference
    # blocks pass through it untouched.
    mt = annotate_adj(mt)
    mt = adjust_sex_ploidy(mt, mt.sex_karyotype, xy_str="XY", xx_str="XX")
    v30_ht = v3_release_sites(public=True).versions["3.0"].ht()
    mt = mt.annotate_rows(_v30_af=v30_ht[mt.row_key].freq[0].AF)
    mt = mt.annotate_entries(
        GT=hom_alt_depletion_fix(
            mt.GT,
            het_non_ref_expr=mt._het_non_ref,
            af_expr=mt._v30_af,
            ab_expr=mt.AD[1] / mt.DP,
            use_v3_1_correction=True,
        )
    )
    gm_ht = group_membership(test=test, data_set="gnomad", environment=environment).ht()
    mt = mt.annotate_cols(group_membership=gm_ht[mt.col_key].group_membership)
    mt = mt.annotate_globals(
        freq_meta=gm_ht.index_globals().freq_meta,
        freq_meta_sample_count=gm_ht.index_globals().freq_meta_sample_count,
    )
    # A row field, so freq and age hists come from one pass over the dense MT.
    mt = mt.annotate_rows(age_hists=age_hists_expr(mt.adj, mt.GT, mt.age))
    ht = compute_freq_by_strata(mt, select_fields=["age_hists"])
    return ht.checkpoint(new_temp_file("consent_freq_and_hists", "ht"))


def _subtract_consent(v4_ht: hl.Table, consent_ht: hl.Table) -> hl.Table:
    """
    Subtract the consent-withdrawal counts and age histograms from the v4 table.

    Sites the withdrawn samples do not carry subtract nothing, so every v4 row
    passes through with its values intact.

    :param v4_ht: v4 release table.
    :param consent_ht: Table from :func:`_consent_freq_ht`.
    :return: v4 table with ``freq``, ``histograms.age_hists`` and the strata
        globals updated.
    """
    consent = consent_ht[v4_ht.key]
    ht = v4_ht.annotate(consent_freq=consent.freq, consent_age_hists=consent.age_hists)
    ht = ht.annotate_globals(
        consent_freq_meta=consent_ht.index_globals().freq_meta,
        consent_freq_meta_sample_count=consent_ht.index_globals().freq_meta_sample_count,
    )
    g = ht.index_globals()
    freq, freq_meta, counts = merge_freq_arrays(
        [ht.freq, ht.consent_freq],
        [g.freq_meta, g.consent_freq_meta],
        operation="diff",
        count_arrays={
            "freq_meta_sample_count": [
                g.freq_meta_sample_count,
                g.consent_freq_meta_sample_count,
            ]
        },
    )
    ht = ht.annotate(
        freq=freq,
        histograms=ht.histograms.annotate(
            age_hists=_merge_hist_struct(
                ht.histograms.age_hists, ht.consent_age_hists, operation="diff"
            )
        ),
    ).drop("consent_freq", "consent_age_hists")
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=counts["freq_meta_sample_count"],
    )
    return ht.checkpoint(new_temp_file("gnomad_freq_minus_consent", "ht"))


def process_gnomad_dataset(
    environment: str,
    test: bool = False,
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
) -> hl.Table:
    """
    Build the gnomAD v5 frequency table from the v4 release.

    v5 adds no gnomAD samples, so instead of recomputing from the VDS the withdrawn
    samples' frequencies and age histograms are subtracted from v4. FAF, grpmax and
    inbreeding are left to the merge step because they depend on the combined strata.

    :param environment: Compute environment.
    :param test: Use test paths and restrict v4 to the consent sites.
    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading partitions.
    :param chrom: Single contig; the v4 table is scoped to it so a per-contig output
        holds only that contig.
    :return: gnomAD v5 frequency table.
    """
    v4_ht = release_sites(data_type="genomes").ht()
    if chrom:
        v4_ht = hl.filter_intervals(v4_ht, [hl.parse_locus_interval(chrom)])
    vds = _prepare_consent_vds(test_vds, test_partitions, chrom)
    consent_ht = _consent_freq_ht(vds, test, environment)
    if test:
        v4_ht = v4_ht.filter(hl.is_defined(consent_ht[v4_ht.key]))
        v4_ht = v4_ht.naive_coalesce(100).checkpoint(new_temp_file("v4_test", "ht"))
    ht = _subtract_consent(v4_ht, consent_ht)
    # raw_qual_hists was approved for removal from v5; the AoU side never builds it.
    ht = ht.annotate(histograms=ht.histograms.drop("raw_qual_hists"))
    # Recomputed from the v5 meta rather than adjusted from v4's global: v5 both
    # drops the withdrawn samples and carries corrected ages.
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

    gnomAD and AoU tables carry just ``freq`` and ``histograms`` because FAF, grpmax
    and the inbreeding coefficient are only meaningful on the merged cohort; the
    merged table also gets the index dicts the release and VCF export look fields
    up by.

    :param ht: Table to trim.
    :param dataset: ``"gnomad"``, ``"aou"`` or ``"merged"``.
    :return: Trimmed table.
    """
    if dataset not in ["gnomad", "aou", "merged"]:
        raise ValueError(f"Invalid dataset: {dataset}")
    if dataset in ["gnomad", "aou"]:
        final_globals = ["freq_meta", "freq_meta_sample_count", "age_distribution"]
        if dataset == "aou":
            # Only AoU has downsamplings; gnomAD's were dropped with the v4 subset.
            final_globals.append("downsamplings")
        return ht.select("freq", "histograms").select_globals(*final_globals)
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
    gnomad_freq_ht: hl.Table, aou_freq_ht: hl.Table
) -> hl.Table:
    """
    Combine the gnomAD and AoU frequency tables into one ``freq`` array.

    Matching strata are summed into the global (all-sample) strata, and the AoU
    strata are kept again as an ``aou`` subset so the release can report AoU alone.
    AoU downsamplings only mean something within AoU, so they appear in the subset
    but not in the global strata. Histograms and the age distribution are summed
    the same way.

    :param gnomad_freq_ht: gnomAD v5 frequency table.
    :param aou_freq_ht: AoU frequency table.
    :return: Merged table with ``freq``, ``histograms`` and the strata globals.
    """
    # One index of the AoU table for both fields, so the compiled query has one join.
    aou = aou_freq_ht[gnomad_freq_ht.key]
    ht = gnomad_freq_ht.annotate(aou_freq=aou.freq, aou_histograms=aou.histograms)
    ag = aou_freq_ht.index_globals()
    ht = ht.annotate_globals(
        aou_freq_meta=ag.freq_meta,
        aou_freq_meta_sample_count=ag.freq_meta_sample_count,
        aou_age_distribution=ag.age_distribution,
        aou_downsamplings=ag.downsamplings,
    )
    g = ht.index_globals()
    merged_freq, merged_meta, counts = merge_freq_arrays(
        [ht.freq, ht.aou_freq],
        [g.freq_meta, g.aou_freq_meta],
        operation="sum",
        count_arrays={
            "counts": [g.freq_meta_sample_count, g.aou_freq_meta_sample_count]
        },
    )
    merged_meta_list = hl.eval(merged_meta)
    global_idx = [i for i, d in enumerate(merged_meta_list) if "downsampling" not in d]
    aou_subset_meta = g.aou_freq_meta.map(
        lambda d: hl.dict(d.items().append(("subset", "aou")))
    )
    ht = ht.annotate(
        freq=hl.array([merged_freq[i] for i in global_idx]).extend(ht.aou_freq),
        histograms=hl.struct(
            qual_hists=_merge_hist_struct(
                ht.histograms.qual_hists, ht.aou_histograms.qual_hists
            ),
            age_hists=_merge_hist_struct(
                ht.histograms.age_hists, ht.aou_histograms.age_hists
            ),
        ),
    )
    return ht.annotate_globals(
        freq_meta=hl.literal([merged_meta_list[i] for i in global_idx]).extend(
            aou_subset_meta
        ),
        freq_meta_sample_count=hl.array(
            [counts["counts"][i] for i in global_idx]
        ).extend(g.aou_freq_meta_sample_count),
        age_distribution=merge_histograms(
            [g.age_distribution, g.aou_age_distribution], operation="sum"
        ),
    )


def calculate_faf_and_grpmax_annotations(ht: hl.Table) -> hl.Table:
    """
    Add FAF, grpmax, fafmax and the inbreeding coefficient, overall and for the AoU subset.

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
# The relay machinery is batch_fanout.py, shared with compute_coverage.py. What stays
# here is freq policy: where the layout lives and how far it is folded, the flags a
# freq relay gets, the relay command, and the chunk worker itself.
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
    exist. A ``--test-region`` run is one chunk with no layout JSON.

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
        valued += [
            ("read_subintervals", "--read-subintervals"),
            ("all_sites_an_suffix", "--all-sites-an-suffix"),
        ]
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

    Resolved once per orchestrator run so every relay it submits checks out the same
    commit. The Hail version comes from the image, so no pin is passed.

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

    :func:`run_chunk_fanout` skips the chunks with a ``_SUCCESS``, dispatches the
    rest in waves of ``--wave-size`` and records the misses, so rerunning this step
    resumes a partial run. This function only names the batch and builds the freq
    relay command.

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

    The chunk's read intervals come from the coverage layout (looked up by index with
    ``hailtop.fs``, so no QoB job or VDS open is needed to find them) or, for a
    ``--test-region`` run, from the region itself. The layout hash is stamped into a
    global for provenance.

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
        an_environment=args.an_environment,
        test=args.test,
        test_vds=args.test_vds,
        read_intervals=sub_intervals,
        all_sites_an_suffix=args.all_sites_an_suffix,
    )
    ht = ht.annotate_globals(freq_chunk_intervals_hash=intervals_hash)
    ht.write(args.chunk_output, overwrite=True)
    logger.info("Wrote chunk %d to %s", args.chunk_index, args.chunk_output)


def _orchestrate_freq_merge(args: argparse.Namespace, ht_path: str) -> None:
    """
    Union the chunk HTs into the final AoU freq HT with a tree of relay merge jobs.

    Chunks are enumerated with the same helper as the fan-out, and a missing chunk
    stops the merge rather than silently dropping loci. :func:`orchestrate_tree_merge`
    builds the tree: intermediate groups that already exist are skipped unless
    ``--overwrite``, and the final job always writes.

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
    args.an_environment = args.an_environment or args.environment
    if args.environment == "batch" and args.app_name is None:
        # Name the Batch after the steps so runs are findable in the Batch UI.
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

    # Relay workers.
    if args.run_chunk:
        _run_freq_chunk(args)
        return
    if args.run_merge:
        _initialize_hail(args)
        union_and_write_hts(
            args.merge_inputs, args.merge_output, args.merge_coalesce_to
        )
        return

    # Suffixes keep per-contig and tagged outputs on their own paths;
    # --merge-datasets reads the untagged, assembled AoU HT.
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
        aou_freq = get_freq(data_set="aou", suffix=aou_out_suffix, **freq_kwargs)
        if args.use_batch_fanout:
            # Not gated on the final HT: the fan-out writes chunks, and a rerun
            # after a merge must still be allowed.
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
            canonical = get_freq(data_set=data_set, suffix=base_suffix, **freq_kwargs)
            per_contig = [
                get_freq(
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
            gnomad_freq = get_freq(
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

        if args.process_aou:
            aou_freq = get_freq(data_set="aou", suffix=aou_out_suffix, **freq_kwargs)
            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq]},
                overwrite=args.overwrite,
            )
            ht = process_aou_dataset(
                environment=args.environment,
                an_environment=args.an_environment,
                test=args.test,
                test_vds=args.test_vds,
                test_partitions=args.test_partitions,
                chrom=args.chrom,
                read_intervals=_region_read_intervals(
                    args.test_region, args.read_subintervals
                ),
                all_sites_an_suffix=args.all_sites_an_suffix,
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
            gnomad_freq = get_freq(data_set="gnomad", **freq_kwargs)
            aou_freq = get_freq(
                data_set="aou", suffix=args.aou_freq_ht_suffix, **freq_kwargs
            )
            check_resource_existence(
                input_step_resources={
                    "process-gnomad": [gnomad_freq],
                    "process-aou": [aou_freq],
                }
            )
            ht = merge_gnomad_and_aou_frequencies(gnomad_freq.ht(), aou_freq.ht())
            ht = ht.checkpoint(new_temp_file("merged_freq", "ht"))
            ht = calculate_faf_and_grpmax_annotations(ht)
            ht = select_final_dataset_fields(ht, dataset="merged")
            ht.write(merged.path, overwrite=args.overwrite)
    finally:
        # Batch keeps its own logs; there is no local file to copy.
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
        "--all-sites-an-suffix",
        default=None,
        help=(
            "Read the AN table compute_coverage wrote with --cov-and-an-output-suffix"
            " <suffix> instead of the default."
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
        "--process-aou",
        action="store_true",
        help=(
            "AoU frequencies as one job: AC, homozygote counts and histograms from the"
            " VDS, AN joined from the coverage table. Use --use-batch-fanout for the"
            " genome."
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
        help="Combine the gnomAD and AoU tables and add FAF, grpmax and inbreeding.",
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
            " compute_coverage --test-region) so its output can be checked against an"
            " AN table made over the same region. Turns on test paths."
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
            "Bucket the AN table and the chunk layout beside it are read from; pass"
            " dataproc when compute_coverage ran with --results-environment dataproc."
            " Default: the value of --environment."
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
        "Relay sizing for --use-batch-fanout and --merge-freq-chunks.",
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
        help="Cores per chunk QoB worker; Hail Batch accepts 1, 2, 4 or 8 for JVM jobs.",
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

    # The testing scopes all switch the run to test paths.
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
