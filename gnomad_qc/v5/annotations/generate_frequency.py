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
import hashlib
import json
import logging
import os
import re
import subprocess
from collections.abc import Sequence
from datetime import datetime, timezone
from typing import Any, NamedTuple

import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs
from gnomad.resources.grch38.gnomad import GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    _read_reduction_globals,
    age_hists_expr,
    bi_allelic_site_inbreeding_expr,
    compute_freq_by_strata,
    expand_strata_array_from_leaves,
    faf_expr,
    gen_anc_faf_max_expr,
    get_adj_expr,
    grpmax_expr,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
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

# Relay jobs and their nested QoB jobs stay in the one region that can reach the
# AoU VDS.
BATCH_REGIONS = ["us-central1"]
BATCH_REMOTE_TMPDIR = "gs://fc-11093c2b-590e-424a-91ac-0cc040d562fc/batch-tmp"
# Hail 0.2.137 is the lowest version that runs on the Batch workers' Java 21.
DEFAULT_BATCH_IMAGE = (
    "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:0.2.137"
)
RELAY_SCRIPT = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/generate_frequency.py"
AGE_HIST_BINS = (30, 80, 10)


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------


def _apply_path_suffix(path: str, suffix: str | None) -> str:
    """
    Insert ``_<suffix>`` before the ``.ht`` extension; unchanged when ``suffix`` is empty.

    Same convention as ``compute_coverage.py``, so a suffixed AN table can be read.

    :param path: HT path ending in ``.ht``.
    :param suffix: Suffix without the leading underscore.
    :return: Suffixed path.
    """
    if not suffix:
        return path
    return path.rstrip("/").removesuffix(".ht") + f"_{suffix}.ht"


def _combine_suffix(*parts: str | None) -> str | None:
    """
    Join the non-empty suffix parts with underscores so each scope writes its own path.

    :param parts: Suffix fragments; empty ones are skipped.
    :return: Joined suffix, or None when every part is empty.
    """
    kept = [p for p in parts if p]
    return "_".join(kept) if kept else None


def _parse_region_interval(
    s: str, reference_genome: str = "GRCh38"
) -> hl.utils.Interval:
    """
    Parse ``contig:start-end`` into a half-open locus interval.

    Half-open, like ``compute_coverage.py``'s regions, so adjacent regions never share
    a locus and a region test lines up with the AN table made over the same region.

    :param s: Interval string; commas in positions are allowed.
    :param reference_genome: Reference genome name.
    :return: ``[start, end)`` locus interval.
    """
    contig, span = s.split(":")
    start_pos, end_pos = (int(p.replace(",", "")) for p in span.split("-"))
    return hl.Interval(
        hl.Locus(contig, start_pos, reference_genome=reference_genome),
        hl.Locus(contig, end_pos, reference_genome=reference_genome),
        includes_start=True,
        includes_end=False,
    )


def _split_intervals(
    intervals: list[hl.utils.Interval], n: int
) -> list[hl.utils.Interval]:
    """
    Tile each half-open interval into ``n`` equal position slices.

    ``read_vds`` partitions by the intervals it is given, so a region read as one
    interval lands in one partition and its whole aggregation runs on one core.
    Slices tile the input exactly, so no locus is dropped or counted twice.

    :param intervals: Half-open locus intervals.
    :param n: Slices per interval; ``<= 1`` returns the input unchanged.
    :return: The slices, in order.
    """
    if n <= 1:
        return intervals
    out = []
    for iv in intervals:
        contig, rg = iv.start.contig, iv.start.reference_genome
        start, end = iv.start.position, iv.end.position
        if end - start <= n:
            out.append(iv)
            continue
        step = (end - start) // n
        bounds = [start + i * step for i in range(n)] + [end]
        out.extend(
            hl.Interval(
                hl.Locus(contig, bounds[i], reference_genome=rg),
                hl.Locus(contig, bounds[i + 1], reference_genome=rg),
                includes_start=True,
                includes_end=False,
            )
            for i in range(n)
        )
    return out


def _spans_sex_chromosome(
    intervals: list[hl.utils.Interval] | None, chrom: str | None = None
) -> bool:
    """
    Return whether a read scope can include chrX or chrY.

    The sex-ploidy adjustment is the identity on autosomes, so a scope that cannot
    touch chrX or chrY skips it and its self-join. An unscoped read counts as spanning.

    :param intervals: Read intervals, or None.
    :param chrom: Single contig, or None.
    :return: True unless the scope is known to exclude both sex contigs.
    """
    if intervals:
        rg = intervals[0].start.reference_genome
        sex = set(rg.x_contigs) | set(rg.y_contigs)
        order = {c: k for k, c in enumerate(rg.contigs)}
        contigs: set[str] = set()
        for iv in intervals:
            contigs.update(
                rg.contigs[order[iv.start.contig] : order[iv.end.contig] + 1]
            )
        return bool(contigs & sex)
    if chrom:
        rg = hl.get_reference("GRCh38")
        return chrom in set(rg.x_contigs) | set(rg.y_contigs)
    return True


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
    path = _apply_path_suffix(
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
    :param read_subintervals: Slices per region (see :func:`_split_intervals`).
    :return: Read intervals, or None when no region was given.
    """
    if not test_region:
        return None
    intervals = [_parse_region_interval(r) for r in test_region]
    return _split_intervals(intervals, read_subintervals or 1)


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
        skip_sex_ploidy=not _spans_sex_chromosome(read_intervals, chrom),
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

    an_path = _apply_path_suffix(
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
    v4_ht: hl.Table,
    test_vds: bool = False,
    test_partitions: int | None = None,
    chrom: str | None = None,
) -> hl.vds.VariantDataset:
    """
    Build the VDS of samples that withdrew consent after v4, prepared as v4 was.

    adj, sex ploidy and the hom-alt depletion fix are applied in the same order as
    the v4 genomes release so the subtracted counts match what v4 added.

    :param v4_ht: v4 release table, for the AF the hom-alt fix needs.
    :param test_vds: Read the test VDS.
    :param test_partitions: Keep only this many leading partitions.
    :param chrom: Single contig to scope the load to.
    :return: Prepared, split VDS.
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
    vmt = vmt.select_cols(
        gen_anc=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        age=vmt.meta.project_meta.age,
    )
    vmt = vmt.select_entries(
        "LA", "LAD", "DP", "GQ", "LGT", _het_non_ref=vmt.LGT.is_het_non_ref()
    )
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = vds.checkpoint(new_temp_file("consent_samples_vds", "vds"))
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    vmt = vds.variant_data
    vmt = vmt.annotate_rows(v4_af=v4_ht[vmt.row_key].freq[0].AF)
    # adj on the diploid call, then ploidy, then the hom-alt fix: the v4 genomes order.
    vmt = vmt.annotate_entries(adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD))
    vmt = vmt.select_entries(
        "AD",
        "DP",
        "GQ",
        "_het_non_ref",
        "adj",
        GT=adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype),
    )
    vmt = vmt.annotate_entries(
        GT=hom_alt_depletion_fix(
            vmt.GT,
            het_non_ref_expr=vmt._het_non_ref,
            af_expr=vmt.v4_af,
            ab_expr=vmt.AD[1] / vmt.DP,
            use_v3_1_correction=True,
        )
    )
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    return vds.checkpoint(new_temp_file("consent_samples_vds_prepared", "vds"))


def _consent_freq_ht(
    vds: hl.vds.VariantDataset, test: bool, environment: str
) -> hl.Table:
    """
    Frequencies and age histograms of the consent-withdrawal samples.

    Densifying is affordable because the VDS holds only the withdrawn samples. The
    strata come from the gnomAD membership table coverage was built with, so the
    subtracted arrays line up with v4 and with the AN table.

    :param vds: Prepared VDS from :func:`_prepare_consent_vds`.
    :param test: Use test paths.
    :param environment: Compute environment.
    :return: Table with ``freq`` and ``age_hists``.
    """
    mt = hl.vds.to_dense_mt(vds)
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
    vds = _prepare_consent_vds(v4_ht, test_vds, test_partitions, chrom)
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
# One genome-wide QoB query is one driver, and a driver crash loses everything, so the
# AoU run is cut into chunks. Each chunk runs from a small non-spot Batch job (a
# "relay") that starts its own QoB driver; a crash costs one chunk, a rerun skips the
# chunks that already have a _SUCCESS marker, and the orchestrator never starts Hail.
# Chunk outputs are namespaced by a hash of the chunk layout because coverage samples
# its layout without a seed: a chunk left by an older layout must never be merged.
# ---------------------------------------------------------------------------


class _RelayJobSpec(NamedTuple):
    """One relay job for :func:`_submit_relay_batch`."""

    name: str
    cpu: float
    memory: str
    storage: str
    command: str


def _resolve_commit() -> str:
    """
    Return the gnomad_qc commit the relays check out.

    ``GNOMAD_QC_COMMIT`` wins so an orchestrator launched from a tarball checkout
    (no ``.git``) can still pin the relays.

    :return: Full commit hash.
    """
    return os.getenv("GNOMAD_QC_COMMIT") or (
        subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    )


def _build_setup_command(
    commit: str, methods_branch: str, gcp_billing_project: str = "broad-mpg-gnomad"
) -> str:
    """
    Return the shell prefix every relay runs before the script.

    Both repos change faster than the image, so the relay pulls them at the pinned
    commit and branch on start-up. The Hail config and the ``quota_project_id``
    patch exist because a bare container has no gcloud config to hand the
    requester-pays project to Hail's Java GCS client. The Hail version comes from the
    image.

    :param commit: gnomad_qc commit to check out.
    :param methods_branch: gnomad_methods branch or commit to check out.
    :param gcp_billing_project: Requester-pays project.
    :return: Shell command string ending in a newline.
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = (
        "https://github.com/broadinstitute/gnomad_methods/archive/"
        f"{methods_branch}.tar.gz"
    )
    config_body = (
        "[batch]\n"
        "billing_project = gnomad-production\n"
        f"remote_tmpdir = {BATCH_REMOTE_TMPDIR}\n"
        "[gcs_requester_pays]\n"
        f"project = {gcp_billing_project}\n"
    )
    return (
        "set -euxo pipefail\n"
        "mkdir -p ~/.config/hail ~/.hail\n"
        "cat > ~/.config/hail/config.ini <<'HAILCFG'\n"
        f"{config_body}"
        "HAILCFG\n"
        "cp ~/.config/hail/config.ini ~/.hail/config.ini\n"
        f"python3 -c \"import json, os; p='/gsa-key/key.json';"
        f" d=json.load(open(p)); d['quota_project_id']='{gcp_billing_project}';"
        f" json.dump(d, open(p+'.new','w')); os.replace(p+'.new', p)\"\n"
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_branch.replace('/', '-')} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def _relay_context(args: argparse.Namespace) -> tuple[str, str, dict]:
    """
    Return the commit, setup command and backend kwargs every relay submission shares.

    Built once per orchestrator so the fan-out and the merge cannot drift onto
    different commits or images.

    :param args: Parsed CLI args.
    :return: ``(commit, setup_cmd, backend_kwargs)``.
    """
    commit = _resolve_commit()
    setup_cmd = _build_setup_command(commit, args.methods_branch)
    backend_kwargs = {"billing_project": args.billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir
    return commit, setup_cmd, backend_kwargs


def _submit_relay_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    batch_name: str,
    job_specs: Sequence[_RelayJobSpec],
    log_label: str,
) -> int | None:
    """
    Submit one Hail Batch of relay jobs and wait for it.

    Relays are non-spot because a preempted relay orphans the QoB batch it waits
    on, and they run once because a Batch retry cannot cancel that orphan and would
    race it. Nothing is submitted when ``job_specs`` is empty.

    :param args: Parsed CLI args (``batch_image``, ``batch_dry_run``).
    :param backend_kwargs: kwargs for ``hb.ServiceBackend``.
    :param batch_name: Hail Batch name.
    :param job_specs: Jobs to submit.
    :param log_label: Noun for log messages ("chunk", "merge").
    :return: Batch id, or None when nothing ran.
    """
    if not job_specs:
        logger.info(
            "No pending %s jobs for %s; nothing submitted.", log_label, batch_name
        )
        return None
    backend = hb.ServiceBackend(**backend_kwargs)
    try:
        batch = hb.Batch(name=batch_name, backend=backend)
        for spec in job_specs:
            j = batch.new_job(name=spec.name)
            j.image(args.batch_image)
            j.cpu(spec.cpu)
            j.memory(spec.memory)
            j.storage(spec.storage)
            j.regions(BATCH_REGIONS)
            j.spot(False)
            j.n_max_attempts(1)
            j.command(spec.command)
        logger.info(
            "Submitting Hail Batch '%s': %d %s jobs (dry_run=%s)",
            batch_name,
            len(job_specs),
            log_label,
            args.batch_dry_run,
        )
        submitted = batch.run(dry_run=args.batch_dry_run)
        return getattr(submitted, "id", None)
    finally:
        backend.close()


def _chunk_intervals_hash(data: dict[str, Any]) -> str:
    """
    Return a 16-hex-char content hash of a chunk layout.

    :param data: Parsed chunk-intervals JSON.
    :return: First 16 hex chars of the SHA-256 of its canonical serialization.
    """
    payload = {k: v for k, v in data.items() if k != "intervals_hash"}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


def _test_region_hash(test_region: Sequence[str]) -> str:
    """
    Return the layout hash of a ``--test-region`` run, which has no layout JSON.

    Hashing the region keeps two region tests under one output path from seeing
    each other's chunk as already present.

    :param test_region: Region strings as given on the command line.
    :return: ``test_region_`` plus a 16-hex-char hash.
    """
    return f"test_region_{_chunk_intervals_hash({'test_region': list(test_region)})}"


def _interval_from_list(t: Sequence, reference_genome: str) -> hl.utils.Interval:
    """
    Rebuild a locus interval from the chunk-intervals JSON form.

    :param t: ``[start_contig, start_pos, end_contig, end_pos, includes_start,
        includes_end]``.
    :param reference_genome: Reference genome name.
    :return: Locus interval.
    """
    sc, sp, ec, ep, incs, ince = t
    return hl.Interval(
        hl.Locus(sc, sp, reference_genome=reference_genome),
        hl.Locus(ec, ep, reference_genome=reference_genome),
        includes_start=incs,
        includes_end=ince,
    )


def _group_chunk_layout(data: dict[str, Any], n_chunks: int | None) -> dict[str, Any]:
    """
    Fold consecutive coverage chunks into about ``n_chunks`` freq chunks.

    Coverage sizes its chunks for its own heavy compute; freq's work per chunk is so
    small that the fixed relay plus QoB driver cost per job dominates, so fewer,
    larger jobs are cheaper. A folded chunk reads the union of its parts, never
    across a contig, so results and the AN join are unchanged. The folded layout
    hashes differently, so its chunks live in their own directory.

    :param data: Parsed chunk-intervals JSON.
    :param n_chunks: Target chunk count; None, or at least the current count, keeps
        the layout as is.
    :return: Layout with the folded ``chunks``.
    """
    n_coverage = len(data["chunks"])
    if n_chunks is None or n_chunks >= n_coverage:
        return data
    if n_chunks < 1:
        raise ValueError(f"--n-chunks must be >= 1, got {n_chunks}")
    per_job = -(-n_coverage // n_chunks)
    grouped: list[dict[str, Any]] = []
    for c in data["chunks"]:
        last = grouped[-1] if grouped else None
        if last and last["contig"] == c["contig"] and last["_n"] < per_job:
            last["intervals"].extend(c["intervals"])
            last["_n"] += 1
        else:
            grouped.append(
                {"contig": c["contig"], "intervals": list(c["intervals"]), "_n": 1}
            )
    for g in grouped:
        del g["_n"]
    return {**data, "chunks": grouped}


def _load_chunk_layout(
    an_environment: str, test: bool, n_chunks: int | None
) -> tuple[dict[str, Any], str]:
    """
    Load the chunk layout ``compute_coverage.py`` wrote beside the AN table.

    Read with ``hailtop.fs`` so the orchestrator and every worker find their chunks
    without a QoB job. All of them must pass the same ``n_chunks`` so they agree on
    the chunk indices and the layout hash.

    :param an_environment: Environment (bucket) the AN table lives in.
    :param test: Read the test-scoped path.
    :param n_chunks: ``--n-chunks`` fold target.
    :return: ``(layout, layout hash)``.
    :raises FileNotFoundError: if the layout has not been written.
    """
    base = coverage_and_an_path(test=test, environment=an_environment).path
    path = base.rstrip("/").removesuffix(".ht") + "_chunk_intervals.json"
    if not file_exists(path):
        raise FileNotFoundError(
            f"Chunk layout not found at {path}; run compute_coverage.py"
            " --write-chunk-intervals first (and pass --an-environment if the AN"
            " table lives in another bucket)."
        )
    with hfs.open(path) as f:
        data = _group_chunk_layout(json.load(f), n_chunks)
    logger.info("Using chunk layout %s (%d freq chunks).", path, len(data["chunks"]))
    return data, _chunk_intervals_hash(data)


def _chunk_path(ht_path: str, idx: int, intervals_hash: str) -> str:
    """
    Return a chunk's HT path, beside the final HT and namespaced by the layout hash.

    :param ht_path: Final AoU freq HT path.
    :param idx: Chunk index.
    :param intervals_hash: Layout hash.
    :return: ``<ht>_chunks/<hash>/<idx:08d>.chunk.ht``.
    """
    base = ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_chunks/{intervals_hash}/{idx:08d}.chunk.ht"


def _group_path(
    ht_path: str, level: int, group_idx: int, merge_group_size: int, intervals_hash: str
) -> str:
    """
    Return a merge-tree intermediate HT path.

    Level, tree shape and layout hash are in the directory so a rerun with a
    different tree or layout writes fresh instead of reusing stale groups.

    :param ht_path: Final AoU freq HT path.
    :param level: Merge-tree level, 1-indexed.
    :param group_idx: Group index within the level.
    :param merge_group_size: Inputs per merge job.
    :param intervals_hash: Layout hash.
    :return: ``<ht>_merge_groups_gs<size>/<hash>/L<level>_<group>.ht``.
    """
    base = ht_path.rstrip("/").removesuffix(".ht")
    return (
        f"{base}_merge_groups_gs{merge_group_size}/{intervals_hash}/"
        f"L{level:02d}_{group_idx:08d}.ht"
    )


def _list_present_chunk_indices(ht_path: str, intervals_hash: str) -> set[int]:
    """
    Return the chunk indices that have a ``_SUCCESS`` marker.

    One listing instead of one GCS stat per chunk, and keyed on ``_SUCCESS`` so a
    half-written chunk is rerun. A missing directory lists as empty.

    :param ht_path: Final AoU freq HT path.
    :param intervals_hash: Layout hash.
    :return: Completed chunk indices.
    """
    base = ht_path.rstrip("/").removesuffix(".ht")
    present: set[int] = set()
    for entry in hfs.ls(f"{base}_chunks/{intervals_hash}/*/_SUCCESS"):
        m = re.search(r"/(\d+)\.chunk\.ht/_SUCCESS$", entry.path)
        if m:
            present.add(int(m.group(1)))
    return present


def _write_failed_chunks_manifest(
    ht_path: str,
    intervals_hash: str,
    failed: Sequence[int],
    n_dispatched: int,
    commit: str,
    app_name: str | None,
    waves: Sequence[dict[str, Any]],
) -> str | None:
    """
    Record which dispatched chunks did not land, or clear the record when all did.

    Rerunning the fan-out already resumes from the missing chunks; the manifest is
    the durable record of what failed and when, since log scrollback is not.

    :param ht_path: Final AoU freq HT path.
    :param intervals_hash: Layout hash.
    :param failed: Dispatched chunk indices with no ``_SUCCESS``.
    :param n_dispatched: Chunks dispatched by this run.
    :param commit: gnomad_qc commit the relays ran.
    :param app_name: ``--app-name`` the relays used.
    :param waves: Per-wave records.
    :return: Manifest path when one was written, else None.
    """
    base = ht_path.rstrip("/").removesuffix(".ht")
    path = f"{base}_chunks/{intervals_hash}/_failed_chunks.json"
    if not failed:
        if file_exists(path):
            hfs.remove(path)
        return None
    payload = {
        "written_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "commit": commit,
        "app_name": app_name,
        "intervals_hash": intervals_hash,
        "n_dispatched": n_dispatched,
        "n_failed": len(failed),
        "failed_chunk_indices": sorted(failed),
        "waves": list(waves),
    }
    with hfs.open(path, "w") as f:
        f.write(json.dumps(payload, indent=2) + "\n")
    return path


def _chunk_scope(args: argparse.Namespace) -> tuple[list[int], int, str]:
    """
    Return the chunks this run covers: eligible indices, total count and layout hash.

    Shared by the fan-out and the merge so they can never disagree on which chunks
    exist. A ``--test-region`` run is one chunk with no layout JSON.

    :param args: Parsed CLI args.
    :return: ``(eligible indices, total chunks, layout hash)``.
    """
    if args.test_region:
        return [0], 1, _test_region_hash(args.test_region)
    data, intervals_hash = _load_chunk_layout(
        args.an_environment, args.test, args.n_chunks
    )
    contigs = [c["contig"] for c in data["chunks"]]
    if not args.chrom:
        return list(range(len(contigs))), len(contigs), intervals_hash
    eligible = [i for i, c in enumerate(contigs) if c == args.chrom]
    if not eligible:
        raise ValueError(
            f"No chunks on --chrom {args.chrom}; the layout covers {sorted(set(contigs))}."
        )
    return eligible, len(contigs), intervals_hash


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
    valued = [
        ("app_name", "--app-name"),
        ("n_chunks", "--n-chunks"),
        ("chunk_driver_cores", "--driver-cores"),
        ("chunk_driver_memory", "--driver-memory"),
        ("chunk_worker_cores", "--worker-cores"),
        ("chunk_worker_memory", "--worker-memory"),
    ]
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


def _orchestrate_freq_fanout(args: argparse.Namespace, ht_path: str) -> None:
    """
    Submit one relay job per pending chunk, in sequential waves of ``--wave-size``.

    Waves bound how many relays and nested QoB drivers run at once, and keep each
    ``batch.run()`` sequential, which Hail Batch's progress display requires. Chunks
    with a ``_SUCCESS`` are skipped, so rerunning this step resumes a partial run. A
    batch does not raise on job failure, so every wave is re-listed afterwards and
    the misses are recorded in a manifest.

    :param args: Parsed CLI args.
    :param ht_path: Final AoU freq HT path the chunks are namespaced under.
    :return: None.
    """
    eligible, n_total, intervals_hash = _chunk_scope(args)
    present = (
        set()
        if args.overwrite
        else _list_present_chunk_indices(ht_path, intervals_hash)
    )
    pending = [i for i in eligible if i not in present]
    logger.info(
        "Fan-out: %d chunks, %d eligible, %d pending, %d already complete.",
        n_total,
        len(eligible),
        len(pending),
        len(eligible) - len(pending),
    )
    if not pending:
        return

    commit, setup_cmd, backend_kwargs = _relay_context(args)
    common = _relay_flags(args, chunk=True)
    scope = "region" if args.test_region else f"{len(pending)}of{n_total}c"
    batch_base = _combine_suffix(
        f"v5_freq_aou_chunk_{scope}_{intervals_hash[:8]}",
        args.aou_freq_ht_suffix,
        args.chrom,
    )
    wave_size = args.wave_size
    if wave_size <= 0 or wave_size >= len(pending):
        waves = [pending]
    else:
        waves = [pending[i : i + wave_size] for i in range(0, len(pending), wave_size)]
    failed: list[int] = []
    records: list[dict[str, Any]] = []
    for wi, wave in enumerate(waves, start=1):
        label = f"_w{wi:03d}of{len(waves):03d}" if len(waves) > 1 else ""
        specs = [
            _RelayJobSpec(
                name=f"freq_chunk_{i:06d}",
                cpu=args.chunk_cpu,
                memory=args.chunk_memory,
                storage=args.chunk_storage,
                command=(
                    f"{setup_cmd}{RELAY_SCRIPT} --run-chunk --chunk-index {i}"
                    f" --chunk-output {_chunk_path(ht_path, i, intervals_hash)} {common}"
                ),
            )
            for i in wave
        ]
        logger.info(
            "Wave %d/%d: %d chunks (%d..%d).",
            wi,
            len(waves),
            len(wave),
            wave[0],
            wave[-1],
        )
        batch_id = _submit_relay_batch(
            args, backend_kwargs, batch_base + label, specs, "chunk"
        )
        if args.batch_dry_run:
            logger.info("Dry run: wave DAG validated; stopping.")
            return
        present = _list_present_chunk_indices(ht_path, intervals_hash)
        wave_failed = [i for i in wave if i not in present]
        failed.extend(wave_failed)
        records.append({"wave": wi, "batch_id": batch_id, "failed": wave_failed})
        if wave_failed:
            logger.warning(
                "Wave %d/%d: %d/%d chunk(s) MISSING after run: %s%s",
                wi,
                len(waves),
                len(wave_failed),
                len(wave),
                wave_failed[:25],
                " ..." if len(wave_failed) > 25 else "",
            )
        else:
            logger.info(
                "Wave %d/%d complete; all %d chunks present.", wi, len(waves), len(wave)
            )
    manifest = _write_failed_chunks_manifest(
        ht_path, intervals_hash, failed, len(pending), commit, args.app_name, records
    )
    if failed:
        logger.warning(
            "Fan-out finished with %d/%d chunk(s) missing; see %s. Rerun"
            " --use-batch-fanout to retry them.",
            len(failed),
            len(pending),
            manifest,
        )
    else:
        logger.info(
            "Fan-out finished: all %d dispatched chunk(s) present.", len(pending)
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
        intervals_hash = _test_region_hash(args.test_region)
        sub_intervals = _region_read_intervals(args.test_region, args.read_subintervals)
    else:
        data, intervals_hash = _load_chunk_layout(
            args.an_environment, args.test, args.n_chunks
        )
        chunks = data["chunks"]
        if not 0 <= args.chunk_index < len(chunks):
            raise ValueError(
                f"chunk index {args.chunk_index} is out of range for a layout of"
                f" {len(chunks)} chunks; the fan-out and the layout are out of sync."
            )
        rg = data["reference_genome"]
        sub_intervals = [
            _interval_from_list(t, rg) for t in chunks[args.chunk_index]["intervals"]
        ]
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


def _union_and_write(
    input_paths: Sequence[str], output_path: str, coalesce_to: int | None = None
) -> None:
    """
    Union HTs that share a schema and globals, and write the result.

    Globals come from the first input, which is right because every input was built
    from the same strata tables.

    :param input_paths: HT paths to union.
    :param output_path: Destination path, overwritten.
    :param coalesce_to: ``naive_coalesce`` target before writing, or None.
    :return: None.
    """
    logger.info("Merging %d HTs -> %s", len(input_paths), output_path)
    merged = hl.Table.union(*[hl.read_table(p) for p in input_paths])
    if coalesce_to is not None:
        merged = merged.naive_coalesce(coalesce_to)
    merged.write(output_path, overwrite=True)


def _orchestrate_freq_merge(args: argparse.Namespace, ht_path: str) -> None:
    """
    Union the chunk HTs into the final AoU freq HT with a tree of relay merge jobs.

    A tree keeps every job's input list and memory bounded at thousands of chunks.
    Chunks are enumerated with the same helper as the fan-out, and a missing chunk
    stops the merge rather than silently dropping loci. Intermediate groups that
    already exist are skipped unless ``--overwrite``.

    :param args: Parsed CLI args.
    :param ht_path: Final AoU freq HT path.
    :return: None.
    """
    eligible, _n_total, intervals_hash = _chunk_scope(args)
    present = _list_present_chunk_indices(ht_path, intervals_hash)
    missing = [i for i in eligible if i not in present]
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} of {len(eligible)} chunks are missing (first:"
            f" {missing[:5]}); rerun --use-batch-fanout first."
        )
    _commit, setup_cmd, backend_kwargs = _relay_context(args)
    common = _relay_flags(args, chunk=False)
    batch_base = _combine_suffix("v5_freq_merge", args.aou_freq_ht_suffix, args.chrom)
    gs = args.merge_group_size

    def merge_job(name, output, inputs, coalesce_to, storage):
        coalesce = (
            f" --merge-coalesce-to {coalesce_to}" if coalesce_to is not None else ""
        )
        return _RelayJobSpec(
            name=name,
            cpu=args.merge_cpu,
            memory=args.merge_memory,
            storage=storage,
            command=(
                f"{setup_cmd}{RELAY_SCRIPT} --run-merge --merge-output {output}{coalesce}"
                f" --merge-inputs {' '.join(inputs)} {common}"
            ),
        )

    inputs = [_chunk_path(ht_path, i, intervals_hash) for i in eligible]
    level = 1
    while len(inputs) > gs:
        groups = [inputs[i : i + gs] for i in range(0, len(inputs), gs)]
        outputs = [
            _group_path(ht_path, level, g, gs, intervals_hash)
            for g in range(len(groups))
        ]
        pending = [
            g
            for g in range(len(groups))
            if args.overwrite or not file_exists(outputs[g])
        ]
        logger.info(
            "Merge level %d: %d inputs -> %d groups, %d pending.",
            level,
            len(inputs),
            len(groups),
            len(pending),
        )
        specs = [
            merge_job(
                f"freq_merge_L{level:02d}_{g:06d}",
                outputs[g],
                groups[g],
                len(groups[g]),
                args.merge_storage,
            )
            for g in pending
        ]
        _submit_relay_batch(
            args, backend_kwargs, f"{batch_base}_L{level:02d}", specs, "merge"
        )
        inputs = outputs
        level += 1
    logger.info("Final merge: %d inputs -> %s", len(inputs), ht_path)
    final = merge_job(
        "freq_merge_final", ht_path, inputs, args.n_partitions, args.final_merge_storage
    )
    _submit_relay_batch(
        args, backend_kwargs, f"{batch_base}_final", [final], "final-merge"
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
        billing_project=args.billing_project,
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
        _union_and_write(args.merge_inputs, args.merge_output, args.merge_coalesce_to)
        return

    # Suffixes keep per-contig and tagged outputs on their own paths;
    # --merge-datasets reads the untagged, assembled AoU HT.
    aou_suffix = _combine_suffix(args.aou_freq_ht_suffix, args.chrom)
    aou_out_suffix = _combine_suffix(aou_suffix, args.freq_output_suffix)
    gnomad_out_suffix = _combine_suffix(args.chrom, args.freq_output_suffix)
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
                    suffix=_combine_suffix(base_suffix, c),
                    **freq_kwargs,
                )
                for c in args.contigs
            ]
            check_resource_existence(
                input_step_resources={"per-contig-freq": per_contig},
                output_step_resources={"assemble-chrom-freq": [canonical]},
                overwrite=args.overwrite,
            )
            _union_and_write(
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
