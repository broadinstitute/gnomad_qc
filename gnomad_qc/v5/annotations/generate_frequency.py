"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies and histograms for:
1. gnomAD dataset - updating v4 frequencies by subtracting consent withdrawal samples
2. AoU dataset - AC / homozygote counts / histograms aggregated from the VDS variant
   data and joined to the all-sites allele number written by compute_coverage.py

Processing Workflow:
--------------------
gnomAD (--process-gnomad):
1. Load v4 frequency table (contains frequencies and age histograms)
2. Prepare consent withdrawal VDS (split multiallelics, annotate metadata)
3. Calculate frequencies and age histograms for consent samples
4. Subtract from v4 frequencies to get updated gnomAD v5 frequencies

AoU (--process-aou):
1. Fan out over the chunk layout compute_coverage.py wrote beside the all-sites-AN
   HT (--use-batch-fanout): one relay per chunk aggregates AC / homozygote_count /
   histograms per stratum from the variant data and joins the chunk's AN.
2. Tree-merge the per-chunk HTs (--merge-freq-chunks).
   Without --use-batch-fanout the same compute runs as a single job (test scopes).

Merged dataset (--merge-datasets):
1. Merge frequency data and histograms from both gnomAD and AoU datasets.
2. Calculate FAF, grpmax, and other post-processing annotations on merged dataset.

Usage Examples:
---------------
# Process AoU (relay fan-out over the coverage chunk layout, then merge).
python generate_frequency.py --process-aou --use-batch-fanout --environment batch
python generate_frequency.py --process-aou --merge-freq-chunks --environment batch

# Process AoU as a single job over a test region.
python generate_frequency.py --process-aou --test-region chr22:10510002-11110002

# Process gnomAD consent withdrawals
python generate_frequency.py --process-gnomad --environment batch

# Run gnomAD in test mode
python generate_frequency.py --process-gnomad --test --test-partitions 2

# Merge both datasets
python generate_frequency.py --merge-datasets --environment batch --app-name "merged_freq" --driver-cores 8 --worker-memory highmem
"""

import argparse
import copy
import hashlib
import json
import logging
import re
import subprocess
from typing import Any, NamedTuple

import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs
from gnomad.resources.grch38.gnomad import GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
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
from gnomad.utils.vcf import SORT_ORDER
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_aou_freq_chunk_path,
    get_freq,
    group_membership,
)
from gnomad_qc.v5.resources.basics import (
    _file_exists_for_env,
    _get_batch_resource_kwargs,
    _init_hail,
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.meta import meta

# Use force=True so that our root handler wins over any handler that
# Hail / hailtop / absl / other deps installed during their imports above.
# Without force=True, basicConfig is a no-op when the root logger already has
# a handler, and all `logger.info` calls get silently dropped when running
# locally against the batch backend (Mode 1 QoB).
logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("v5_frequency")

# Hail Batch regions for relay chunk/merge jobs (mirrors compute_coverage).
BATCH_REGIONS = ["us-central1"]
# Chunk/merge relay image: Hail 0.2.137 + gnomad_methods deps baked in.
# 0.2.137 is the floor for QoB: the Batch worker JVM is Java 21, which the
# 0.2.128 JAR cannot run on.
DEFAULT_BATCH_IMAGE = (
    "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:0.2.137"
)
logger.setLevel(logging.INFO)


def _combine_freq_suffix(suffix: str | None, chrom: str | None) -> str | None:
    """
    Fold an optional ``--chrom`` contig into the freq HT suffix.

    A ``--chrom`` run appends the contig to the suffix so every per-contig output
    (AoU freq HT, chunk/group HTs) is stored at its own path -- a per-contig run can
    go all the way through fan-out -> merge without colliding with or clobbering
    another contig's outputs, and a late failure only loses that contig. Assemble the
    per-contig HTs with ``--assemble-chrom-freq``. Folded in one place (not by mutating
    ``args``) so the orchestrator and the relay workers it spawns -- which each
    re-resolve from the same ``suffix`` + ``chrom`` -- always agree, with no risk of
    double-folding.

    :param suffix: Optional base suffix (e.g. ``--aou-freq-ht-suffix``).
    :param chrom: Optional single contig (``--chrom``).
    :return: ``"{suffix}_{chrom}"`` when both set, else whichever is set, else None.
    """
    if suffix and chrom:
        return f"{suffix}_{chrom}"
    return suffix or chrom


def _apply_path_suffix(path: str, suffix: str | None) -> str:
    """
    Insert ``_<suffix>`` before the ``.ht`` extension, or return unchanged if no suffix.

    Mirrors ``compute_coverage.py``'s ``_apply_path_suffix`` (underscore convention), so
    a suffixed all-sites-AN HT written by ``compute_coverage.py``
    (``--cov-and-an-output-suffix``, e.g. ``coverage_and_an_chunkhash_test.ht``) can be
    targeted by the frequency all-sites-AN read.

    :param path: HT path ending in ``.ht``.
    :param suffix: Optional suffix string (no leading underscore). If falsy, ``path`` is
        returned unchanged.
    :return: Suffix-applied path.
    """
    if not suffix:
        return path
    return path.rstrip("/").removesuffix(".ht") + f"_{suffix}.ht"


def _parse_region_interval(
    s: str, reference_genome: str = "GRCh38"
) -> hl.utils.Interval:
    """
    Parse a ``contig:start-end`` string into a half-open Python locus interval.

    Mirrors ``compute_coverage.py``'s ``_parse_region_interval`` so a ``--test-region``
    freq run scopes to exactly the same loci an all-sites-AN test HT was generated over.
    Returns a concrete ``hl.utils.Interval`` (not an ``hl.parse_locus_interval``
    expression) with half-open ``[start, end)`` bounds, so adjacent regions (e.g. two
    consecutive 50 kb intervals) stay disjoint -- no locus is double-counted at the
    shared boundary.

    :param s: Interval string, e.g. ``chr1:55058666-55108666`` (commas allowed).
    :param reference_genome: Reference-genome name. Default "GRCh38".
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


def _split_intervals_for_read(
    intervals: list[hl.utils.Interval], n: int
) -> list[hl.utils.Interval]:
    """
    Split each half-open locus interval into ``n`` roughly equal position sub-intervals.

    ``hl.vds.read_vds(intervals=...)`` partitions the read BY the intervals it is given,
    so passing N tiling sub-intervals for a region makes the variant read land in N
    partitions instead of 1 -- giving the downstream per-variant aggregation N-way row
    parallelism with no shuffle (each sub-interval reads its own slice). Safe for the
    variant-only all-sites-AN aggregation (point loci; no reference-block straddle to
    undercount). Sub-intervals stay half-open and tile the original exactly, so no locus
    is dropped or double-counted. Position-based (not variant-density-balanced), so
    per-partition variant counts may be uneven -- fine for parallelism, and simple.

    :param intervals: Half-open locus intervals to split.
    :param n: Number of sub-intervals per input interval (``<= 1`` returns unchanged).
    :return: Flattened list of ``n`` sub-intervals per input interval.
    """
    if n <= 1:
        return intervals
    out = []
    for iv in intervals:
        contig = iv.start.contig
        rg = iv.start.reference_genome
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


def mt_hist_fields(mt: hl.MatrixTable) -> hl.StructExpression:
    """
    Annotate allele balance quality metrics histograms and age histograms onto MatrixTable.

    :param mt: Input MatrixTable.
    :return: Struct with allele balance, quality metrics histograms, and age histograms.
    """
    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    qual_hists = qual_hist_expr(
        gt_expr=mt.GT,
        gq_expr=mt.GQ,
        dp_expr=hl.sum(mt.AD),
        adj_expr=mt.adj,
        ab_expr=(mt.AD[1] / hl.sum(mt.AD)),
        split_adj_and_raw=True,
    )
    return hl.struct(
        qual_hists=qual_hists,
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
    )


def _aou_release_meta_small(environment: str = "batch") -> hl.Table:
    """
    Return sex karyotype and age for the AoU release samples as a 10-partition Table.

    The sample meta HT is ~330 partitions for ~365k rows, so joining or aggregating it
    directly fans ~330 tiny QoB tasks per use. Selecting the two fields freq needs and
    coalescing (lazily; nothing is written) cuts that to 10 tasks per use.

    :param environment: Environment to use. Default is "batch".
    :return: Table keyed by ``s`` with ``sex_karyotype`` and ``age``.
    """
    meta_ht = meta(data_type="genomes", environment=environment).ht()
    meta_ht = meta_ht.filter((meta_ht.project_meta.project == "aou") & meta_ht.release)
    return meta_ht.select("sex_karyotype", age=meta_ht.project_meta.age).naive_coalesce(
        10
    )


def _aou_age_distribution(meta_small: hl.Table):
    """
    Return the AoU release age-distribution histogram.

    Chunk-independent, so each fan-out chunk recomputes the same value; over the
    10-partition table from :func:`_aou_release_meta_small` that is one small aggregate.

    :param meta_small: Table from :func:`_aou_release_meta_small`.
    :return: Age histogram struct over the AoU release samples.
    """
    return meta_small.aggregate(hl.agg.hist(meta_small.age, 30, 80, 10))


def _load_release_aou_vds(
    environment: str,
    test_vds: bool = False,
    filter_partitions: list[int] | None = None,
    read_intervals: list[hl.utils.Interval] | None = None,
    chrom: str | None = None,
) -> hl.vds.VariantDataset:
    """
    Load the AoU VDS restricted to the release samples, with no per-row bookkeeping.

    ``get_aou_vds(release_only=True)`` removes the hard-filtered samples with
    ``remove_dead_alleles=True`` and applies the release filter through
    ``hl.vds.filter_samples``; each of those walks every row's full sample array (a
    dead-allele count plus an LA/LAD rewrite, and an "any entries left" count). The
    release samples are already free of hard-filtered samples, and the frequency
    calculation drops rows with raw AC 0 over the release samples anyway, so neither
    pass changes the output. Here the columns are filtered directly to the release
    list (the permanent JSON ``write_aou_vds_sample_jsons`` writes) and rows are left
    alone; the empty ones fall out at the raw-AC filter.

    :param environment: Compute environment.
    :param test_vds: Whether to load the test VDS.
    :param filter_partitions: Optional partition indices to read.
    :param read_intervals: Optional locus intervals to prune the read to.
    :param chrom: Optional single contig.
    :return: VDS whose variant data holds only release samples.
    """
    from gnomad_qc.v5.resources.meta import load_aou_sample_artifact_json

    release = load_aou_sample_artifact_json("release_samples.json", environment)
    if release is None:
        raise FileNotFoundError(
            "release_samples.json is missing; run write_aou_vds_sample_jsons first."
        )
    vds = get_aou_vds(
        release_only=False,
        remove_hard_filtered_samples=False,
        remove_dead_alleles=False,
        # Prefix colliding sample IDs so they match the post-prefix release list.
        add_project_prefix=True,
        annotate_meta=False,
        log_sample_counts=False,
        test=test_vds,
        filter_partitions=filter_partitions,
        read_intervals=read_intervals,
        chrom=chrom,
        environment=environment,
    )
    keep = hl.Table.parallelize([hl.struct(s=s) for s in release], key="s")
    vmt = vds.variant_data
    vmt = vmt.filter_cols(hl.is_defined(keep[vmt.col_key]))
    return hl.vds.VariantDataset(vds.reference_data, vmt)


def _aou_group_membership_ht(
    test: bool = False, environment: str = "batch"
) -> tuple[hl.Table, bool]:
    """
    Read the AoU group-membership HT ``compute_coverage.py`` built the all-sites AN with.

    ``compute_coverage.py --write-group-membership-ht`` writes the cell-reduced HT
    (``generate_freq_group_membership_array(reduce_to_cells=True)``) at the
    ``_cells``-suffixed path beside the ``group_membership`` resource: one boolean per
    distinct membership pattern, so each sample is in exactly one adj cell and one raw
    cell, and globals ``freq_meta_full`` / ``freq_meta_sample_count_full`` /
    ``freq_leaf_indices`` / ``freq_group_decomposition`` map the cells back to the
    full strata list. Reading that HT (not a copy) keeps freq's AC on the same sample
    sets as the AN it is joined to. It is already written at 10 partitions.

    :param test: Whether to read the test-scoped path.
    :param environment: Compute environment.
    :return: ``(ht, reduced)`` -- ``reduced`` is True when the cells HT was read;
        False when only the full (one boolean per stratum) HT exists.
    """
    gm_path = group_membership(test=test, data_set="aou", environment=environment).path
    cells_path = _apply_path_suffix(gm_path, "cells")
    if _file_exists_for_env(cells_path, environment):
        logger.info("Using cell-reduced AoU group membership HT: %s", cells_path)
        return hl.read_table(cells_path), True
    logger.warning(
        "Cell-reduced AoU group membership HT not found at %s; using the full HT at %s.",
        cells_path,
        gm_path,
    )
    return hl.read_table(gm_path), False


def _spans_sex_chromosome(
    intervals: list[hl.utils.Interval] | None, chrom: str | None = None
) -> bool:
    """
    Return whether the read scope touches chrX or chrY.

    ``adjusted_sex_ploidy_expr`` only changes genotypes on those contigs (its first
    case returns the genotype unchanged on autosomes), so a chunk, region, or
    ``--chrom`` scope that touches neither can skip it. Skipping saves the self-join
    of the variant rows and columns the expression builds
    (``annotate_and_index_source_mt_for_sex_ploidy``). An unscoped (whole-genome)
    read is treated as spanning.

    :param intervals: Locus intervals of the read scope, or None.
    :param chrom: Single contig of the read scope, or None.
    :return: True if the scope may include chrX or chrY.
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


def _prepare_aou_vds(
    aou_vds: hl.vds.VariantDataset,
    test: bool = False,
    environment: str = "batch",
    age_distribution=None,
    skip_sex_ploidy: bool = False,
) -> hl.MatrixTable:
    """
    Prepare the AoU variant data for the all-sites-AN frequency calculation.

    Joins sex karyotype and age onto the columns from a small meta view
    (:func:`_aou_release_meta_small`), adjusts sex ploidy, annotates adj, splits
    multi-allelics, and sets the strata globals from the group membership HT
    ``compute_coverage.py`` built the all-sites AN with. Order matters: ploidy ->
    adj -> split, so adj is computed on the local (``LGT``/``LAD``) fields.

    :param aou_vds: AoU VariantDataset (release samples).
    :param test: Whether running in test mode.
    :param environment: Environment being used. Default is "batch". Must be one of "rwb"
        or "batch".
    :param age_distribution: Precomputed AoU age-distribution histogram to set as the
        global. When None, computed here via ``_aou_age_distribution``.
    :param skip_sex_ploidy: Skip the sex-ploidy adjustment. Only valid when the read
        scope has no chrX/chrY loci (see :func:`_spans_sex_chromosome`), where the
        adjustment is the identity. Default False.
    :return: Prepared, split AoU variant MatrixTable.
    """
    aou_vmt = aou_vds.variant_data
    # The group membership HT compute_coverage built the all-sites AN with (its
    # cell-reduced form when present); freq's strata must match the AN's. Only its
    # globals are used here -- the aggregation step joins the membership itself.
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht, gm_reduced = _aou_group_membership_ht(
        test=test, environment=environment
    )

    logger.info("Selecting cols for frequency stratification...")
    # Only sex_karyotype (sex-ploidy adjustment) and age (age_hists) are needed on the
    # columns; join them from the small coalesced meta table, not the 330-partition one.
    meta_small = _aou_release_meta_small(environment)
    meta_indexed = meta_small[aou_vmt.col_key]
    aou_vmt = aou_vmt.select_cols(
        sex_karyotype=meta_indexed.sex_karyotype, age=meta_indexed.age
    )
    if skip_sex_ploidy:
        logger.info("Read scope has no chrX/chrY loci; skipping sex-ploidy adjustment.")
        lgt = aou_vmt.LGT
    else:
        lgt = adjusted_sex_ploidy_expr(
            aou_vmt.locus, aou_vmt.LGT, aou_vmt.sex_karyotype
        )
    aou_vmt = aou_vmt.select_entries(
        LGT=lgt, GQ=aou_vmt.GQ, LAD=aou_vmt.LAD, LA=aou_vmt.LA
    )
    # AoU adj uses the shared helper (also used by coverage and variant QC):
    # the usual gnomAD cutoffs with DP approximated as sum(LAD).
    aou_vmt = annotate_adj_no_dp(aou_vmt)
    aou_vds = hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt)
    aou_vds = hl.vds.split_multi(aou_vds, filter_changed_loci=True)
    aou_vmt = aou_vds.variant_data

    logger.info("Annotating globals...")
    # Age-distribution global comes from the sample metadata table (one scan of the
    # sample-keyed meta table), NOT aggregate_cols on the prepared variant MT -- the
    # latter would force a full variant-MT pass per downstream eager action. Set as a
    # literal global, so it survives the aggregation step.
    if age_distribution is None:
        age_distribution = _aou_age_distribution(meta_small)
    # The cells HT keeps the full strata list under the `_full` globals.
    gg = group_membership_ht.index_globals()
    return aou_vmt.select_globals(
        freq_meta=gg.freq_meta_full if gm_reduced else gg.freq_meta,
        freq_meta_sample_count=(
            gg.freq_meta_sample_count_full if gm_reduced else gg.freq_meta_sample_count
        ),
        age_distribution=age_distribution,
        downsamplings=gg.downsamplings,
    )


def _sparse_strata_and_hists(
    mt: hl.MatrixTable, group_membership_ht: hl.Table
) -> hl.Table:
    """
    Aggregate AC / homozygote_count per stratum and the histograms over defined entries only.

    Output-identical to ``agg_by_strata`` over the same MatrixTable plus
    ``mt_hist_fields``, but each row visits only its defined entries. VDS variant data
    stores an entry only for samples with a non-reference call, so a row's entries
    array (one slot per sample, ~365k) is ~99.9% missing: the median variant has 3
    carriers. Aggregating over the full array walks every slot once per aggregation
    pass; here the array is compacted once and every aggregator runs over the
    carriers.

    Semantics reproduced exactly:
      - a call is counted in each stratum its sample belongs to; for adj strata only
        when ``adj`` is True (a missing ``adj`` counts as False, as ``hl.agg.filter``
        does);
      - ``AC`` = sum of alt alleles, ``homozygote_count`` = number of hom-var calls
        (missing genotypes contribute nothing, as with ``hl.agg.sum`` /
        ``count_where``);
      - the histograms come from the same ``qual_hist_expr`` / ``age_hists_expr``
        expressions as :func:`mt_hist_fields`, so binning is unchanged.

    :param mt: Prepared, split variant MatrixTable with entry fields ``GT``, ``GQ``,
        ``AD``, ``adj`` and column field ``age``.
    :param group_membership_ht: Group membership HT (full or cells); its
        ``group_membership`` array and ``freq_meta`` global define the strata.
    :return: Table keyed like ``mt`` with row fields ``hist_fields``, ``AC``,
        ``homozygote_count`` (arrays over the HT's strata) and ``mt``'s globals.
    """
    freq_meta = [dict(m) for m in hl.eval(group_membership_ht.freq_meta)]
    n_groups = len(freq_meta)
    adj_groups = hl.literal([m.get("group", "NA") == "adj" for m in freq_meta])
    # Per sample: the indices of the strata it belongs to.
    gm = group_membership_ht.select(
        strata=hl.enumerate(group_membership_ht.group_membership)
        .filter(lambda t: t[1])
        .map(lambda t: t[0])
    )
    mt = mt.annotate_cols(strata=gm[mt.col_key].strata)
    lt = mt.localize_entries("entries", "cols")
    # Compact to the defined entries, carrying each carrier's age and the strata its
    # call is counted in (raw strata always, adj strata only when the call is adj).
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
    empty = hl.struct(AC=hl.int64(0), homozygote_count=hl.int64(0))
    lt = lt.annotate(
        hist_fields=carriers.aggregate(lambda c: mt_hist_fields(c)),
        _by_stratum=carriers.aggregate(
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
        AC=hl.range(n_groups).map(lambda g: lt._by_stratum.get(g, empty).AC),
        homozygote_count=hl.range(n_groups).map(
            lambda g: lt._by_stratum.get(g, empty).homozygote_count
        ),
    )
    return lt.select("hist_fields", "AC", "homozygote_count").select_globals(
        *[g for g in lt.globals.dtype if g != "cols"]
    )


def _calculate_aou_frequencies_and_hists_using_all_sites_ans(
    aou_variant_mt: hl.MatrixTable,
    test: bool = False,
    environment: str = "batch",
    chrom: str | None = None,
    region_intervals: list[hl.utils.Interval] | None = None,
    all_sites_an_suffix: str | None = None,
    an_environment: str | None = None,
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using all sites ANs.

    :param aou_variant_mt: Prepared variant MatrixTable.
    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :param an_environment: Environment the all-sites-AN HT is read from (selects its
        bucket; see ``--an-environment``). Default None (same as ``environment``).
    :param chrom: Optional single contig; when set, the all-sites-AN HT is
        interval-filtered to it so only that contig's AN is read (the variant MT is
        already contig-scoped by the VDS load, so this is a read-pruning optimization).
    :param region_intervals: TESTING ONLY: optional list of half-open ``hl.Interval``
        objects (from ``--test-region``); when set, the all-sites-AN HT is filtered to
        exactly these intervals (the same ones used to scope the variant MT read), so
        the AN join matches an all-sites-AN test HT generated for the same region.
        Takes precedence over ``chrom`` for the AN-HT filter.
    :param all_sites_an_suffix: Optional suffix (no leading underscore) inserted before
        the ``.ht`` extension of the all-sites-AN HT path, to read a suffixed AN HT
        written by ``compute_coverage.py --cov-and-an-output-suffix`` (e.g.
        ``chunkhash_test`` -> ``...coverage_and_an_chunkhash_test.ht``). Default None
        (the base ``coverage_and_an.ht``).
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    all_sites_an_path = _apply_path_suffix(
        coverage_and_an_path(test=test, environment=an_environment or environment).path,
        all_sites_an_suffix,
    )
    logger.info("Reading all sites AN HT from %s...", all_sites_an_path)
    all_sites_an_ht = hl.read_table(all_sites_an_path)
    if region_intervals:
        all_sites_an_ht = hl.filter_intervals(all_sites_an_ht, region_intervals)
    elif chrom:
        all_sites_an_ht = hl.filter_intervals(
            all_sites_an_ht, [hl.parse_locus_interval(chrom)]
        )
    logger.info("Annotating frequencies with all sites ANs...")
    # Read the SAME group membership the all-sites-AN HT was built with, so the AC /
    # homozygote_count strata line up with the AN array on the join below. That is the
    # cell-reduced HT when compute_coverage wrote one (each sample in one adj cell and
    # one raw cell); the cells are summed back to the full strata below using the
    # decomposition globals the HT carries. The full HT is ~330 partitions, so it is
    # coalesced lazily to 10 before the column join.
    group_membership_ht, reduced = _aou_group_membership_ht(
        test=test, environment=environment
    )
    if not reduced:
        group_membership_ht = group_membership_ht.naive_coalesce(10)

    # Per-stratum AC / homozygote_count and the histograms, aggregated over each row's
    # defined entries only (see _sparse_strata_and_hists). Only the adj qual_hists is
    # kept in v5 (raw_qual_hists was approved for removal); dropping it here, before the
    # checkpoint, lets Hail prune its aggregators.
    aou_variant_freq_ht = _sparse_strata_and_hists(aou_variant_mt, group_membership_ht)
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        hist_fields=aou_variant_freq_ht.hist_fields.annotate(
            qual_hists=aou_variant_freq_ht.hist_fields.qual_hists.drop("raw_qual_hists")
        )
    )

    # With the cells HT, AC / homozygote_count come back over cells only; expand each
    # back to the full strata by summation (both are summable; a stratum's cells are
    # disjoint) using the decomposition the HT carries, and restore the full freq_meta
    # so the arrays line up with the full-strata all-sites AN joined below.
    if reduced:
        gg = group_membership_ht.index_globals()
        # Batch the four global reads into ONE hl.eval round-trip -- each hl.eval is a
        # separate QoB driver execution.
        ev = hl.eval(
            hl.struct(
                leaf_indices=gg.freq_leaf_indices,
                decomposition=gg.freq_group_decomposition,
                freq_meta_full=gg.freq_meta_full,
                freq_meta_sample_count_full=gg.freq_meta_sample_count_full,
            )
        )
        leaf_indices = ev.leaf_indices
        decomposition = {i: d for i, d in enumerate(ev.decomposition) if d}
        freq_meta_full = ev.freq_meta_full
        freq_meta_sample_count_full = ev.freq_meta_sample_count_full
        n_full = len(freq_meta_full)
        aou_variant_freq_ht = aou_variant_freq_ht.annotate(
            AC=expand_strata_array_from_leaves(
                aou_variant_freq_ht.AC, leaf_indices, decomposition, n_full
            ),
            homozygote_count=expand_strata_array_from_leaves(
                aou_variant_freq_ht.homozygote_count,
                leaf_indices,
                decomposition,
                n_full,
            ),
        )
        aou_variant_freq_ht = aou_variant_freq_ht.annotate_globals(
            freq_meta=freq_meta_full,
            freq_meta_sample_count=freq_meta_sample_count_full,
        )

    # Keep only rows with at least one alt allele among the release samples. This is
    # the row-set rule that makes the loader's per-row passes unnecessary (see
    # _load_release_aou_vds): alleles carried only by non-release samples, which used
    # to appear as AC 0 rows, are dropped.
    freq_meta_now = [dict(m) for m in hl.eval(aou_variant_freq_ht.freq_meta)]
    raw_idx = freq_meta_now.index({"group": "raw"})
    aou_variant_freq_ht = aou_variant_freq_ht.filter(
        aou_variant_freq_ht.AC[raw_idx] > 0
    )

    # Checkpoint the per-variant aggregated freq HT (AC/homozygote_count arrays +
    # hist_fields) before the AN join, freq-struct build, and write. Running the whole
    # pipeline (split_multi -> sparse hists + strata aggregation -> AN left-join -> write)
    # as ONE fused query compiles to a ~3,500-node IR and a huge volume of generated JVM
    # bytecode on the DRIVER, OOMing a standard (~8GB) driver at write time even though
    # the data itself is tiny (~130MB RegionPool). Materializing this compact
    # (rows-only, no sample columns) aggregated HT splits the query in two -- the heavy
    # hist+strata aggregation on one side, the AN join + finalize + write on the other
    # -- so each compiles far less code and a standard driver suffices.
    aou_variant_freq_ht = aou_variant_freq_ht.checkpoint(
        new_temp_file("aou_allsites_freq_agg", "ht")
    )

    # Load AN values from all sites ANs table (calculated by another script but used
    # same group membership HT so same strata order).
    logger.info("Annotating AN values from all sites ANs...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        all_sites_an=all_sites_an_ht[aou_variant_freq_ht.locus].AN
    )

    logger.info("Building complete frequency struct with imported AN values...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        freq=hl.map(
            lambda AC, hom_alt, AN: hl.struct(
                AC=hl.int32(AC),
                AF=hl.if_else(AN > 0, AC / AN, hl.missing(hl.tfloat64)),
                AN=hl.int32(AN),
                homozygote_count=hl.int32(hom_alt),
            ),
            aou_variant_freq_ht.AC,
            aou_variant_freq_ht.homozygote_count,
            aou_variant_freq_ht.all_sites_an,
        ),
    ).drop("all_sites_an")

    # Nest histograms to match gnomAD structure. Only the adj-filtered ``qual_hists``
    # is kept; ``raw_qual_hists`` was already dropped before the checkpoint above (so its
    # aggregators are pruned, ~halving the per-variant quality-hist cost).
    aou_variant_freq_ht = aou_variant_freq_ht.select(
        freq=aou_variant_freq_ht.freq,
        histograms=hl.struct(
            qual_hists=aou_variant_freq_ht.hist_fields.qual_hists.qual_hists,
            age_hists=aou_variant_freq_ht.hist_fields.age_hists,
        ),
    )

    return aou_variant_freq_ht


def process_aou_dataset(
    test_vds: bool = False,
    test_partitions: int = None,
    environment: str = "batch",
    chrom: str | None = None,
    test_region: list[str] | None = None,
    all_sites_an_suffix: str | None = None,
    read_subintervals: int | None = None,
    an_environment: str | None = None,
) -> hl.Table:
    """
    Process the All of Us dataset for frequency calculations and age histograms.

    Single-job all-sites-AN path: load the AoU VDS, prepare it, aggregate AC /
    homozygote_count and the histograms per stratum, and join the all-sites AN
    ``compute_coverage.py`` wrote. The production run uses the relay fan-out
    (``--use-batch-fanout``) over the same compute per chunk.

    :param test_vds: Whether to run in test mode on test VDS.
    :param test_partitions: Number of partitions to use in test mode. Default is None.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :param chrom: Optional single contig to scope the AoU VDS load (and the
        all-sites-AN read) to. Default is None (whole genome).
    :param test_region: TESTING ONLY: optional list of ``contig:start-end`` strings
        (half-open) to scope the AoU VDS load (and the all-sites-AN read) to, so the
        output can be matched against an all-sites-AN test HT for the same region.
        Auto-enables test mode. Default None.
    :param all_sites_an_suffix: Optional suffix for the all-sites-AN HT path (to read a
        suffixed AN HT from ``compute_coverage.py --cov-and-an-output-suffix``).
        Default None.
    :param an_environment: Environment the all-sites-AN HT is read from (see
        ``--an-environment``). Default None (same as ``environment``).
    :return: Table with freq and age_hists annotations for AoU dataset.
    """
    # --test-region is a testing-only scope, so it enables test mode (test paths).
    test = test_vds or test_partitions is not None or test_region is not None
    # Parse the region strings once into half-open intervals reused for both the VDS
    # load and the all-sites-AN HT filter (same semantics as the AN test HT).
    region_intervals = (
        [_parse_region_interval(r) for r in test_region] if test_region else None
    )
    # Optionally split the region into sub-intervals so the variant read lands in
    # ``read_subintervals`` partitions (N-way row parallelism for the aggregation)
    # instead of the single partition a whole-region read prunes to.
    if region_intervals and read_subintervals:
        region_intervals = _split_intervals_for_read(
            region_intervals, read_subintervals
        )
        logger.info(
            "Split --test-region into %d read sub-interval(s) for parallelism.",
            len(region_intervals),
        )

    # --test-region scopes via read-time interval pruning (read_intervals), NOT a
    # post-read filter_intervals: the latter leaves the variant_data read fanned across
    # thousands of (empty) partitions for a small region. Columns are filtered to the
    # release samples with no per-row passes (see _load_release_aou_vds).
    aou_vds = _load_release_aou_vds(
        environment,
        test_vds=test_vds,
        filter_partitions=(
            list(range(test_partitions))
            if (test_partitions and not region_intervals)
            else None
        ),
        read_intervals=region_intervals,
        chrom=chrom,
    )
    aou_vmt = _prepare_aou_vds(
        aou_vds,
        test=test,
        environment=environment,
        skip_sex_ploidy=not _spans_sex_chromosome(region_intervals, chrom),
    )

    logger.info("Calculating AoU frequencies and age histograms...")
    aou_freq_ht = _calculate_aou_frequencies_and_hists_using_all_sites_ans(
        aou_vmt,
        test=test,
        environment=environment,
        chrom=chrom,
        region_intervals=region_intervals,
        all_sites_an_suffix=all_sites_an_suffix,
        an_environment=an_environment,
    )
    return select_final_dataset_fields(aou_freq_ht, dataset="aou")


# ---------------------------------------------------------------------------
# Hail Batch fan-out for --process-aou
# ---------------------------------------------------------------------------


def _list_present_freq_chunk_indices(sample_chunk_path: str) -> set[int]:
    """
    Return the set of chunk indices with a completed (``_SUCCESS``) output.

    One ``hailtop.fs.ls`` glob of ``<freq_chunks_dir>/*.freq.chunk_*.ht/_SUCCESS``
    instead of a per-chunk ``file_exists`` probe (tens of thousands of serial GCS
    stats at prod scale). Keying on ``_SUCCESS`` (not the directory) means a
    partially-written chunk is correctly treated as absent; a missing directory
    returns an empty set. Only ``chunk_`` outputs are matched, so per-group merge
    HTs in the same directory are ignored.

    :param sample_chunk_path: Any per-chunk HT path (e.g. from
        ``get_aou_freq_chunk_path(0, kind="chunk", ...)``); its parent directory
        is the freq_chunks directory that is globbed.
    :return: Set of completed chunk indices.
    """
    chunk_dir = sample_chunk_path.rstrip("/").rsplit("/", 1)[0]
    present: set[int] = set()
    try:
        entries = hfs.ls(f"{chunk_dir}/*.freq.chunk_*.ht/_SUCCESS")
    except (FileNotFoundError, OSError):
        # Directory does not exist yet (no chunks written) -> nothing present.
        return present
    for entry in entries:
        m = re.search(r"\.freq\.chunk_(\d+)\.ht/_SUCCESS$", entry.path)
        if m:
            present.add(int(m.group(1)))
    return present


def _build_setup_command(
    commit: str,
    gcp_billing_project: str = "broad-mpg-gnomad",
    methods_branch: str = "main",
) -> str:
    """
    Build shell commands to download gnomad_qc and gnomad_methods and configure Hail.

    Both repos are actively developed, so we pull them at runtime rather
    than relying on what's baked into the Docker image. The image provides
    ``hail`` and system dependencies (g++, curl).

    Mirrors the hardened setup used by ``compute_coverage.py`` (adopted for
    parity):

    - Writes the Hail ``config.ini`` to BOTH the canonical XDG path
      (``~/.config/hail/config.ini``, what ``hailtop.config.get_user_config_path``
      returns) AND the legacy ``~/.hail/config.ini``, so ``hl.init`` finds the
      Batch billing project, ``remote_tmpdir``, and GCS requester-pays project on
      any Hail version.
    - Patches ``/gsa-key/key.json`` with a ``quota_project_id`` field so
      requester-pays reads (the AoU VDS) have a billing-project fallback for the
      container's Java GCS client (works on a laptop via gcloud config, but not in
      a bare container without this).
    - Relies on the batch image (``--batch-image``, default the
      ``v5_freq_batch:0.2.137`` image) to supply the validated Hail version --
      0.2.137 baked in -- instead of a per-job runtime pip reinstall. (0.2.137
      predates the 0.2.138 requester-pays-propagation regression for VDS metadata
      reads; pass a different ``--batch-image`` to change the Hail version.)

    :param commit: Git commit hash to pin gnomad_qc to.
    :param gcp_billing_project: GCP project for requester-pays reads; patched into
        the GSA key as ``quota_project_id`` and written as the requester-pays
        project. Default ``"broad-mpg-gnomad"``.
    :param methods_branch: Branch/commit of gnomad_methods to pull.
        Default is ``"main"``.
    :return: Shell command string.
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = f"https://github.com/broadinstitute/gnomad_methods/archive/{methods_branch}.tar.gz"
    methods_dir_suffix = methods_branch.replace("/", "-")
    config_body = (
        "[batch]\n"
        "billing_project = gnomad-production\n"
        "remote_tmpdir = gs://fc-11093c2b-590e-424a-91ac-0cc040d562fc/batch-tmp\n"
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
        # TODO: Remove this GSA-key patch once Hail's requester-pays propagation
        # to the container's Java GCS client is fixed (likely 0.2.139+).
        f"python3 -c \"import json, os; p='/gsa-key/key.json';"
        f" d=json.load(open(p)); d['quota_project_id']='{gcp_billing_project}';"
        f" json.dump(d, open(p+'.new','w')); os.replace(p+'.new', p)\"\n"
        # Hail version comes from the batch image (default v5_freq_batch:0.2.137),
        # not a runtime reinstall -- keeps container startup fast and the version
        # reproducible. Pass a different --batch-image to change it.
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_dir_suffix} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def _prepare_consent_vds(
    v4_ht: hl.Table,
    test_vds: bool = False,
    test_partitions: int = 2,
    chrom: str | None = None,
) -> hl.vds.VariantDataset:
    """
    Load and prepare VDS for consent withdrawal sample processing.

    :param v4_ht: v4 release table for AF annotation.
    :param test_vds: Whether running in test mode.
    :param test_partitions: Number of partitions to use in test mode. Default is 2.
    :param chrom: Optional single contig to scope the consent VDS load to. Default
        is None (whole genome).
    :return: Prepared VDS with consent samples, split multiallelics, and annotations.
    """
    logger.info("Loading and preparing VDS for consent withdrawal samples...")

    vds = get_gnomad_v5_genomes_vds(
        release_only=True,
        test=test_vds,
        consent_drop_only=True,
        annotate_meta=True,
        # Partition-slice only for a test run with no --chrom; a --chrom run scopes via
        # the contig filter instead, and a full prod run reads all partitions.
        filter_partitions=(
            list(range(test_partitions))
            if (test_partitions is not None and chrom is None)
            else None
        ),
        chrom=chrom,
    )

    logger.info(
        "VDS has been filtered to %s consent withdrawal samples...",
        vds.variant_data.count_cols(),
    )

    vmt = vds.variant_data
    vmt = vmt.select_cols(
        gen_anc=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        age=vmt.meta.project_meta.age,
    )

    logger.info("Selecting entries and annotating non_ref hets pre-split...")
    vmt = vmt.select_entries(
        "LA", "LAD", "DP", "GQ", "LGT", _het_non_ref=vmt.LGT.is_het_non_ref()
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = vds.checkpoint(new_temp_file("consent_samples_vds", "vds"))

    logger.info("Splitting multiallelics in gnomAD sample withdrawal VDS...")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    # Annotate with v4 frequencies for hom alt depletion fix
    vmt = vds.variant_data
    vmt = vmt.annotate_rows(v4_af=v4_ht[vmt.row_key].freq[0].AF)

    # This follows the v3/v4 genomes workflow for adj and sex adjusted genotypes which
    # were added before the hom alt depletion fix.
    # The correct order is to do the hom alt fix before adjusting sex ploidy and before
    # determining the adj annotation because haploid GTs have different adj filtering
    # criteria, but the option to adjust ploidy after adj is included for consistency
    # with v3.1, where we added the adj annotation before adjusting for sex ploidy.
    logger.info("Computing sex adjusted genotypes and quality annotations...")
    vmt = vmt.annotate_entries(
        adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD),
    )
    vmt = vmt.select_entries(
        "AD",
        "DP",
        "GQ",
        "_het_non_ref",
        "adj",
        GT=adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype),
    )

    # We set use_v3_1_correction to True to mimic the v4 genomes approach.
    logger.info("Applying v4 genomes hom alt depletion fix...")
    vmt = vmt.annotate_entries(
        GT=hom_alt_depletion_fix(
            vmt.GT,
            het_non_ref_expr=vmt._het_non_ref,
            af_expr=vmt.v4_af,
            ab_expr=vmt.AD[1] / vmt.DP,
            use_v3_1_correction=True,
        )
    )
    logger.info("Annotating age distribution...")
    vmt = vmt.annotate_globals(
        age_distribution=vmt.aggregate_cols(hl.agg.hist(vmt.age, 30, 80, 10))
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    return vds.checkpoint(new_temp_file("consent_samples_vds_prepared", "vds"))


def _calculate_consent_frequencies_and_age_histograms(
    vds: hl.vds.VariantDataset,
    test_run: bool = False,
    environment: str = "batch",
) -> hl.Table:
    """
    Calculate frequencies and age histograms for consent withdrawal samples.

    :param vds: Prepared VDS with consent samples.
    :param test_run: Whether running in test mode.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Table with freq and age_hists annotations for consent samples.
    """
    logger.info("Densifying VDS for frequency calculations...")
    mt = hl.vds.to_dense_mt(vds)
    # Group membership table is already filtered to consent drop samples and is in GCS.
    group_membership_ht = group_membership(
        test=test_run, data_set="gnomad", environment=environment
    ).ht()

    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership,
    )
    mt = mt.annotate_globals(
        freq_meta=group_membership_ht.index_globals().freq_meta,
        freq_meta_sample_count=group_membership_ht.index_globals().freq_meta_sample_count,
    )

    logger.info(
        "Calculating frequencies and age histograms using compute_freq_by_strata..."
    )

    # Annotate hists_fields on MatrixTable rows before calling compute_freq_by_strata so
    # to keep the age_hists annotation on the frequency table.
    logger.info("Annotating hists_fields on MatrixTable rows...")
    mt = mt.annotate_rows(
        hists_fields=hl.struct(
            age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        )
    )

    logger.info(
        "Computing frequencies for consent samples using compute_freq_by_strata..."
    )
    consent_freq_ht = compute_freq_by_strata(
        mt,
        select_fields=["hists_fields"],
    )

    consent_freq_ht = consent_freq_ht.transmute(
        age_hists=consent_freq_ht.hists_fields.age_hists,
    )

    return consent_freq_ht.checkpoint(new_temp_file("consent_freq_and_hists", "ht"))


def _subtract_consent_frequencies_and_age_histograms(
    v4_ht: hl.Table,
    consent_freq_ht: hl.Table,
) -> hl.Table:
    """
    Subtract consent withdrawal frequencies and age histograms from v4 frequency table.

    :param v4_ht: v4 release table (contains both freq and histograms.age_hists).
    :param consent_freq_ht: Consent withdrawal table with freq and age_hists annotations.
    :return: Updated frequency table with consent frequencies and age histograms subtracted.
    """
    logger.info(
        "Subtracting consent withdrawal frequencies and age histograms from v4 release table..."
    )

    joined_freq_ht = v4_ht.annotate(
        consent_freq=consent_freq_ht[v4_ht.key].freq,
        consent_age_hists=consent_freq_ht[v4_ht.key].age_hists,
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=consent_freq_ht.index_globals().freq_meta,
        consent_freq_meta_sample_count=consent_freq_ht.index_globals().freq_meta_sample_count,
    )

    logger.info("Subtracting consent frequencies...")
    updated_freq_expr, updated_freq_meta, updated_sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.consent_freq],
        [
            joined_freq_ht.index_globals().freq_meta,
            joined_freq_ht.index_globals().consent_freq_meta,
        ],
        operation="diff",
        count_arrays={
            "freq_meta_sample_count": [
                joined_freq_ht.index_globals().freq_meta_sample_count,
                joined_freq_ht.index_globals().consent_freq_meta_sample_count,
            ],
        },
    )
    # Update the frequency table with freq changes.
    joined_freq_ht = joined_freq_ht.annotate(freq=updated_freq_expr)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        freq_meta=updated_freq_meta,
        freq_meta_sample_count=updated_sample_counts["freq_meta_sample_count"],
    )

    logger.info("Subtracting consent age histograms...")
    updated_age_hist_het = merge_histograms(
        [
            joined_freq_ht.histograms.age_hists.age_hist_het,
            joined_freq_ht.consent_age_hists.age_hist_het,
        ],
        operation="diff",
    )
    updated_age_hist_hom = merge_histograms(
        [
            joined_freq_ht.histograms.age_hists.age_hist_hom,
            joined_freq_ht.consent_age_hists.age_hist_hom,
        ],
        operation="diff",
    )

    # Update the frequency table with age hist changes.
    joined_freq_ht = joined_freq_ht.annotate(
        histograms=joined_freq_ht.histograms.annotate(
            age_hists=joined_freq_ht.histograms.age_hists.annotate(
                age_hist_het=updated_age_hist_het,
                age_hist_hom=updated_age_hist_hom,
            )
        ),
    )

    return joined_freq_ht.checkpoint(new_temp_file("merged_freq_and_hists", "ht"))


def select_final_dataset_fields(ht: hl.Table, dataset: str = "gnomad") -> hl.Table:
    """
    Create final dataset freq Table with only desired annotations.

    :param ht: Hail Table containing all annotations.
    :param dataset: Dataset to select final fields, either "gnomad", "aou" or "merged".
    :return: Hail Table with final annotations.
    """
    if dataset not in ["gnomad", "aou", "merged"]:
        raise ValueError(f"Invalid dataset: {dataset}")

    if dataset in ["gnomad", "aou"]:
        final_globals = ["freq_meta", "freq_meta_sample_count", "age_distribution"]
        final_fields = ["freq", "histograms"]

        if dataset == "aou":
            # AoU has one extra 'downsamplings' global field that is not present in
            # gnomAD.
            final_globals.append("downsamplings")

        # Convert all int64 annotations in the freq struct to int32s for merging type
        # compatibility.
        ht = ht.annotate(
            freq=ht.freq.map(
                lambda x: x.annotate(
                    **{k: hl.int32(v) for k, v in x.items() if v.dtype == hl.tint64}
                )
            )
        )
    else:
        sort_order = copy.deepcopy(SORT_ORDER)

        ht = ht.annotate_globals(
            freq_index_dict=make_freq_index_dict_from_meta(
                ht.freq_meta, sort_order=sort_order
            ),
            faf_index_dict=make_freq_index_dict_from_meta(
                ht.faf_meta, sort_order=sort_order
            ),
        )
        final_globals = [
            "freq_meta",
            "freq_index_dict",
            "freq_meta_sample_count",
            "faf_meta",
            "faf_index_dict",
            "age_distribution",
            "aou_downsamplings",
        ]
        final_fields = [
            "freq",
            "faf",
            "grpmax",
            "fafmax",
            "inbreeding_coeff",
            "histograms",
        ]

    return ht.select(*final_fields).select_globals(*final_globals)


def _fix_v4_global_age_distribution(
    freq_ht: hl.Table, environment: str = "batch"
) -> hl.Table:
    """
    Fix the age distribution global annotation in the frequency table.

    :param freq_ht: Frequency table to annotate with the age distribution.
    :param environment: Environment to use. Default is "batch".
    :return: Frequency table with the age distribution global annotation fixed.
    """
    # Use v5 meta as age is already fixed in the v5 project metadata as are the consent
    # withdrawal samples' releasable field.
    meta_ht = meta(environment=environment).ht()
    meta_ht = meta_ht.filter(
        (meta_ht.release) & (meta_ht.project_meta.project == "gnomad")
    )
    meta_ht = meta_ht.annotate_globals(
        age_distribution=meta_ht.aggregate(
            hl.agg.hist(meta_ht.project_meta.age, 30, 80, 10)
        )
    )
    freq_ht = freq_ht.annotate_globals(
        age_distribution=meta_ht.index_globals().age_distribution
    )

    return freq_ht


def process_gnomad_dataset(
    test_vds: bool = False,
    test_partitions: int = 2,
    environment: str = "batch",
    chrom: str | None = None,
) -> hl.Table:
    """
    Process gnomAD dataset to update v4 frequency HT by removing consent withdrawal samples.

    This function performs frequency adjustment by:
    1. Loading v4 frequency HT (contains both frequencies and age histograms)
    2. Loading consent withdrawal VDS
    3. Filtering to sites present in BOTH consent VDS AND v4 frequency table
    4. Calculating frequencies and age histograms for consent withdrawal samples
    5. Subtracting both frequencies and age histograms from v4 frequency HT
    6. Only overwriting fields that were actually updated in the final output

    :param test_vds: Whether to run on test vds. Default is False.
    :param test_partitions: Number of partitions to filter to in test mode. Default is 2.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :param chrom: Optional single contig to scope the consent VDS and the v4 release
        HT to, so the output freq HT covers only that contig. Default is None.
    :return: Updated frequency HT with updated frequencies and age histograms for gnomAD dataset.
    """
    test_run = test_vds or test_partitions is not None
    v4_ht = release_sites(data_type="genomes").ht()
    # Scope the v4 release HT to --chrom so the per-contig output contains ONLY that
    # contig (else the whole-genome v4 freq would pass through unsubtracted for
    # off-contig sites and --assemble-chrom-freq would double-count).
    if chrom:
        v4_ht = hl.filter_intervals(v4_ht, [hl.parse_locus_interval(chrom)])

    vds = _prepare_consent_vds(
        v4_ht,
        test_vds=test_vds,
        test_partitions=test_partitions,
        chrom=chrom,
    )

    logger.info("Calculating frequencies and age histograms for consent samples...")
    consent_freq_ht = _calculate_consent_frequencies_and_age_histograms(
        vds, test_run=test_run, environment=environment
    )

    if test_run:
        v4_ht = v4_ht.filter(hl.is_defined(consent_freq_ht[v4_ht.key]))
        v4_ht = v4_ht.naive_coalesce(100).checkpoint(
            new_temp_file("v4_ht_filtered_test", "ht")
        )

    logger.info("Subtracting consent frequencies and age histograms from v4...")
    updated_freq_ht = _subtract_consent_frequencies_and_age_histograms(
        v4_ht, consent_freq_ht
    )

    logger.info("Merging updated frequency fields...")
    freq_ht = _merge_updated_frequency_fields(v4_ht, updated_freq_ht)

    logger.info("Reannotating gnomAD's age distribution global annotation...")
    freq_ht = _fix_v4_global_age_distribution(freq_ht, environment=environment)

    # The "aou" subset is built at the AoU+gnomAD merge.

    # Select only the fields that were updated as FAF/grpmax/inbreeding_coeff annotations
    # will be calculated on the final merged dataset.
    logger.info("Selecting gnomAD freq HT final fields...")
    final_freq_ht = select_final_dataset_fields(freq_ht, dataset="gnomad")

    return final_freq_ht


def _merge_updated_frequency_fields(
    v4_release_ht: hl.Table, updated_freq_ht: hl.Table
) -> hl.Table:
    """
    Merge frequency tables, only overwriting fields that were actually updated.

    For sites that exist in updated_freq_ht, use the updated values.
    For sites that don't exist in updated_freq_ht, keep original values.

    Note: FAF/grpmax/inbreeding_coeff annotations are not calculated during consent
    withdrawal processing and will be calculated later on the final merged dataset.

    :param v4_release_ht: Original v4 release table.
    :param updated_freq_ht: Updated frequency table with consent withdrawals subtracted.
    :return: Final frequency table with selective field updates.
    """
    logger.info("Merging frequency tables with selective field updates...")

    # Bring in updated values with a single lookup.
    updated_row = updated_freq_ht[v4_release_ht.key]

    # Update freq and age_hists in a single annotate to avoid source mismatch:
    # - freq: use updated if present, otherwise keep original
    # - histograms.age_hists: update only age_hists, preserving qual_hists
    # Drop v4's raw_qual_hists: it was approved for removal from v5 (matching the AoU
    # compute and compute_coverage.py), so the merged output carries only adj
    # qual_hists.
    final_freq_ht = v4_release_ht.annotate(
        freq=hl.coalesce(updated_row.freq, v4_release_ht.freq),
        histograms=v4_release_ht.histograms.annotate(
            age_hists=hl.coalesce(
                updated_row.histograms.age_hists,
                v4_release_ht.histograms.age_hists,
            )
        ).drop("raw_qual_hists"),
    )

    # Update globals from updated table.
    updated_globals = {}
    for global_field in ["freq_meta", "freq_meta_sample_count"]:
        if global_field in updated_freq_ht.globals:
            updated_globals[global_field] = updated_freq_ht.index_globals()[
                global_field
            ]

    if updated_globals:
        final_freq_ht = final_freq_ht.annotate_globals(**updated_globals)

    return final_freq_ht


def merge_gnomad_and_aou_frequencies(
    gnomad_freq_ht: hl.Table,
    aou_freq_ht: hl.Table,
) -> hl.Table:
    """
    Merge frequency data and histograms from gnomAD and All of Us datasets.

    :param gnomad_freq_ht: Frequency Table for gnomAD.
    :param aou_freq_ht: Frequency Table for AoU.
    :return: Merged frequency Table with combined frequencies and histograms.
    """
    # Index the AoU freq HT ONCE and pull both freq and histograms off the same indexed
    # struct so the compiled IR emits a single join (one scan of the genome-wide AoU HT)
    # instead of two. aou_histograms doesn't depend on any field added between here and
    # the old second annotate, and annotate preserves the key, so this is
    # output-identical.
    aou_indexed = aou_freq_ht[gnomad_freq_ht.key]
    joined_freq_ht = gnomad_freq_ht.annotate(
        aou_freq=aou_indexed.freq,
        aou_histograms=aou_indexed.histograms,
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        aou_freq_meta=aou_freq_ht.index_globals().freq_meta,
        aou_freq_meta_sample_count=aou_freq_ht.index_globals().freq_meta_sample_count,
        aou_age_distribution=aou_freq_ht.index_globals().age_distribution,
        aou_downsamplings=aou_freq_ht.index_globals().downsamplings,
    )

    merged_freq, merged_meta, sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.aou_freq],
        [
            joined_freq_ht.index_globals().freq_meta,
            joined_freq_ht.index_globals().aou_freq_meta,
        ],
        operation="sum",
        count_arrays={
            "counts": [
                joined_freq_ht.index_globals().freq_meta_sample_count,
                joined_freq_ht.index_globals().aou_freq_meta_sample_count,
            ],
        },
    )

    # AoU-only downsampling is dropped from global; it's exposed only in the
    # "aou" subset.
    merged_meta_list = hl.eval(merged_meta)
    global_idx = [i for i, d in enumerate(merged_meta_list) if "downsampling" not in d]
    global_meta = [merged_meta_list[i] for i in global_idx]
    global_freq = hl.array([merged_freq[i] for i in global_idx])
    global_counts = hl.array([sample_counts["counts"][i] for i in global_idx])

    # AoU is defined at every site, so the parallel arrays extend without fill.
    aou_subset_meta = joined_freq_ht.index_globals().aou_freq_meta.map(
        lambda d: hl.dict(d.items().append(("subset", "aou")))
    )
    freq = global_freq.extend(joined_freq_ht.aou_freq)
    freq_meta = hl.literal(global_meta).extend(aou_subset_meta)
    freq_meta_sample_count = global_counts.extend(
        joined_freq_ht.index_globals().aou_freq_meta_sample_count
    )

    joined_freq_ht = joined_freq_ht.annotate(freq=freq).annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=freq_meta_sample_count,
        freq_index_dict=make_freq_index_dict_from_meta(freq_meta),
    )

    logger.info("Merging quality histograms and age histograms from both datasets...")
    # aou_histograms was pulled off the single aou_freq_ht index at the top of this
    # function (with aou_freq), so no second join here.

    def _merge_hist_struct(hist1, hist2, operation="sum"):
        """Merge all fields of two histogram structs."""
        return hl.struct(
            **{
                field: merge_histograms(
                    [hist1[field], hist2[field]], operation=operation
                )
                for field in hist1.dtype.fields
            }
        )

    # raw_qual_hists was dropped from v5 (both the AoU compute and the gnomAD v4-reuse
    # side), so only adj qual_hists and age_hists are merged.
    merged_histograms = hl.struct(
        qual_hists=_merge_hist_struct(
            joined_freq_ht.histograms.qual_hists,
            joined_freq_ht.aou_histograms.qual_hists,
        ),
        age_hists=_merge_hist_struct(
            joined_freq_ht.histograms.age_hists,
            joined_freq_ht.aou_histograms.age_hists,
        ),
    )

    joined_freq_ht = joined_freq_ht.annotate(histograms=merged_histograms)

    # Merge age_distribution global histograms (single histogram per dataset, not
    # per-strata like freq_meta_sample_count).
    merged_age_distribution = merge_histograms(
        [
            joined_freq_ht.index_globals().age_distribution,
            joined_freq_ht.index_globals().aou_age_distribution,
        ],
        operation="sum",
    )
    joined_freq_ht = joined_freq_ht.annotate_globals(
        age_distribution=merged_age_distribution
    )

    return joined_freq_ht


def calculate_faf_and_grpmax_annotations(
    ht: hl.Table,
) -> hl.Table:
    """
    Calculate FAF, grpmax, gen_anc_faf_max, and inbreeding coefficient annotations.

    Computes filtering allele frequencies and grpmax for both the full dataset
    (gnomad) and the AoU-only "aou" subset.

    :param ht: Merged frequency table for AoU and gnomAD data containing 'freq' and
        'freq_meta' annotations.
    :return: Table with 'faf', 'grpmax', 'gen_anc_faf_max', and 'inbreeding_coeff'.
    """
    logger.info(
        "Filtering frequencies to just 'aou' subset entries for 'faf' calculations..."
    )
    # Filter to the aou subset and remove the "subset" key from freq_meta so faf_expr
    # pulls the correct indices.
    aou_freq_meta, aou_array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {"freq": ht.freq},
        items_to_filter={"subset": ["aou"]},
        keep=True,
        combine_operator="or",
    )
    aou_freq_meta = aou_freq_meta.map(
        lambda d: hl.dict(d.items().filter(lambda x: x[0] != "subset"))
    )
    # Use the filtered freq expression (already a per-row expr on ht) and the filtered
    # meta directly, instead of wrapping into a table and self-joining it back -- the
    # self-join compiled to a second scan of the merged checkpoint.
    freq_metas = {
        "gnomad": (ht.freq, ht.index_globals().freq_meta),
        "aou": (
            aou_array_exprs["freq"],
            aou_freq_meta,
        ),
    }

    faf_exprs = []
    faf_meta_exprs = []
    grpmax_exprs = {}
    gen_anc_faf_max_exprs = {}

    for dataset, (freq, meta) in freq_metas.items():
        faf, faf_meta = faf_expr(
            freq, meta, ht.locus, GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v5"]
        )
        grpmax = grpmax_expr(freq, meta, GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v5"])
        gen_anc_faf_max = gen_anc_faf_max_expr(faf, faf_meta)

        # Add subset back to aou faf meta.
        if dataset == "aou":
            faf_meta = [{**x, "subset": "aou"} for x in faf_meta]

        faf_exprs.append(faf)
        faf_meta_exprs.append(faf_meta)
        grpmax_exprs[dataset] = grpmax
        gen_anc_faf_max_exprs[dataset] = gen_anc_faf_max

    logger.info(
        "Annotating 'faf', 'grpmax', 'gen_anc_faf_max', and 'inbreeding_coeff'..."
    )
    ht = ht.annotate(
        faf=hl.flatten(faf_exprs),
        grpmax=hl.struct(**grpmax_exprs),
        fafmax=hl.struct(**gen_anc_faf_max_exprs),
        inbreeding_coeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1]),
    )
    faf_meta_exprs = hl.flatten(faf_meta_exprs)
    ht = ht.annotate_globals(
        faf_meta=faf_meta_exprs,
        faf_index_dict=make_freq_index_dict_from_meta(faf_meta_exprs),
    )

    ht = ht.checkpoint(new_temp_file("freq_with_faf", "ht"))

    return ht


def _get_default_app_name(args) -> str | None:
    """
    Derive a Hail Batch app name from the processing-step args.

    Combines the requested processing steps (``--process-gnomad``,
    ``--process-aou``, ``--merge-datasets``) into a single name so the run is
    easy to identify in the Hail Batch UI. Returns ``None`` when no processing
    step is selected.

    :param args: Parsed command-line arguments.
    :return: Default app name, or ``None`` if no steps were requested.
    """
    steps = []
    if args.process_gnomad:
        steps.append("gnomad")
    if args.process_aou:
        steps.append("aou")
    if args.merge_datasets:
        steps.append("merge")

    if not steps:
        return None

    return f"v5_freq_{'_'.join(steps)}"


def _initialize_hail(args) -> None:
    """
    Initialize Hail with appropriate configuration for the environment.

    :param args: Parsed command-line arguments.
    """
    # Auto-derive an app name from the requested processing steps if the user
    # did not pass one explicitly. Only relevant for the batch backend.
    if args.environment == "batch" and args.app_name is None:
        args.app_name = _get_default_app_name(args)

    _init_hail(
        "v5_frequency_generation",
        args.environment,
        billing_project=getattr(args, "gcp_billing_project", None),
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{qc_temp_prefix(environment=args.environment, days=args.tmp_dir_days)}frequency_generation",
        **_get_batch_resource_kwargs(args),
    )


# ===========================================================================
# 1. Chunk layout (read from compute_coverage's chunk-intervals JSON) + path/hash helpers
# ===========================================================================


def _coverage_chunk_intervals_path(an_environment: str, test: bool = False) -> str:
    """
    Return the path of the chunk-intervals JSON ``compute_coverage.py`` writes beside the AN HT.

    ``compute_coverage.py --write-chunk-intervals`` stores its chunk layout as
    ``<coverage_and_an>_chunk_intervals.json`` next to the canonical (unsuffixed)
    all-sites-AN HT, in the bucket its ``--results-environment`` selects. The freq
    fan-out reads that same file, so every freq chunk covers exactly one coverage
    chunk and each chunk's AN read lands on the AN HT partitions written for it. Freq
    has no chunk precompute of its own.

    :param an_environment: Environment the all-sites-AN HT is read from
        (``--an-environment``).
    :param test: If True, return the test-scoped path.
    :return: GCS path to the coverage chunk-intervals JSON.
    """
    base = coverage_and_an_path(test=test, environment=an_environment).path
    return base.rstrip("/").removesuffix(".ht") + "_chunk_intervals.json"


def _normalize_chunk_intervals(data: dict[str, Any]) -> dict[str, Any]:
    """
    Return the chunk-intervals data with freq's key names.

    ``compute_coverage.py`` names the per-chunk interval list ``intervals`` and the
    per-chunk count ``read_subintervals_per_chunk``; freq uses ``sub_intervals`` and
    ``read_subintervals``. The interval lists themselves have the same shape
    (``[start_contig, start_pos, end_contig, end_pos, includes_start, includes_end]``).
    ``ref_block_max_length`` is dropped: the variant-only freq read has no reference
    blocks to widen for. Data already in freq's layout is returned unchanged.

    :param data: Parsed chunk-intervals JSON in either layout.
    :return: Dict with ``read_subintervals``, ``reference_genome``, and ``chunks``
        (each ``{"contig", "sub_intervals"}``).
    """
    if "sub_intervals" in data["chunks"][0]:
        return data
    return {
        "read_subintervals": data["read_subintervals_per_chunk"],
        "reference_genome": data["reference_genome"],
        "chunks": [
            {"contig": c["contig"], "sub_intervals": c["intervals"]}
            for c in data["chunks"]
        ],
    }


def _load_freq_chunk_intervals(
    an_environment: str, test: bool = False
) -> tuple[dict[str, Any], str, str]:
    """
    Load the chunk layout the fan-out, merge, and workers share.

    Read driver-side with ``hailtop.fs`` (no QoB job) from the coverage chunk-intervals
    JSON beside the all-sites-AN HT (see :func:`_coverage_chunk_intervals_path`).

    :param an_environment: Environment the all-sites-AN HT is read from.
    :param test: If True, read the test-scoped path.
    :return: ``(data in freq's layout, path it was read from, layout hash)``.
    :raises FileNotFoundError: if the JSON is absent.
    """
    path = _coverage_chunk_intervals_path(an_environment, test)
    if not file_exists(path):
        raise FileNotFoundError(
            f"chunk-intervals JSON not found at {path}. The freq fan-out uses the"
            " layout compute_coverage.py --write-chunk-intervals writes beside the"
            " all-sites-AN HT; run that first (and pass --an-environment if the AN HT"
            " lives in another bucket)."
        )
    with hfs.open(path) as f:
        data = _normalize_chunk_intervals(json.load(f))
    logger.info("Using chunk layout from %s (%d chunks).", path, len(data["chunks"]))
    return data, path, _freq_chunk_intervals_hash(data)


def _freq_chunk_intervals_hash(data: dict[str, Any]) -> str:
    """
    Return a stable short content hash of the chunk-intervals JSON's boundaries.

    ``compute_coverage.py`` derives its chunk boundaries by sampling the key
    distribution WITHOUT a fixed seed, so identical inputs produce different cut points
    on each of its ``--write-chunk-intervals`` runs. Because chunk outputs are keyed only
    by index and the ``_SUCCESS`` skip-check / merge are existence-only, a chunk left at
    index N by a PRIOR layout (different boundaries) would otherwise pass as "present"
    and be merged with this layout's chunks -> overlapping/duplicate loci. Namespacing
    every chunk output by this hash (folded into the freq chunk-path suffix) keeps
    layouts in separate directories, so a stale chunk is recomputed rather than silently
    merged.

    Computed over the normalized content (see :func:`_normalize_chunk_intervals`), so
    re-loading the same layout always yields the same value.

    :param data: Chunk-intervals data in freq's layout.
    :return: First 16 hex chars of the SHA-256 of the canonical serialization.
    """
    payload = {k: v for k, v in data.items() if k != "intervals_hash"}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


def _combine_suffix(*parts: str | None) -> str | None:
    """
    Join the non-empty suffix parts with underscores (freq uses ``.``/``_`` suffixes).

    Used to fold the run's ``--aou-freq-ht-suffix`` and the chunk-intervals content hash
    into the single ``suffix=`` that ``get_aou_freq_chunk_path`` accepts, so every
    per-chunk / per-group HT is namespaced by its layout.

    :param parts: Suffix fragments (any may be None/empty).
    :return: Underscore-joined suffix, or None if all parts are empty.
    """
    kept = [p for p in parts if p]
    return "_".join(kept) if kept else None


def _freq_chunk_path(
    idx: int,
    intervals_hash: str,
    test: bool = False,
    environment: str = "batch",
    suffix: str | None = None,
) -> str:
    """
    Return the per-chunk freq HT path, namespaced by the chunk-intervals layout hash.

    Reuses freq's ``get_aou_freq_chunk_path(kind="chunk")`` and folds ``intervals_hash``
    (plus any run ``suffix``) into its ``suffix=`` so chunks from a different
    (non-reproducible) layout land in a distinct ``freq_chunks.<...>/`` directory and can
    never be mixed at merge time (compute_coverage's ``_chunk_path`` role).

    :param idx: Chunk index (zero-based).
    :param intervals_hash: Content hash of the chunk-intervals layout that produced it.
    :param test: Whether to use the test-scoped path.
    :param environment: Compute environment.
    :param suffix: Optional run suffix (e.g. ``--aou-freq-ht-suffix`` folded with
        ``--chrom``).
    :return: GCS path for this chunk's partial freq HT.
    """
    return get_aou_freq_chunk_path(
        idx,
        kind="chunk",
        test=test,
        environment=environment,
        suffix=_combine_suffix(suffix, intervals_hash),
    )


def _freq_group_path(
    level: int,
    group_idx: int,
    intervals_hash: str,
    merge_group_size: int,
    test: bool = False,
    environment: str = "batch",
    suffix: str | None = None,
) -> str:
    """
    Return a per-group merged HT path (tree-reduce intermediate level output).

    Level-tagged so a recursive merge tree doesn't overwrite earlier-level outputs, and
    tree-shape-tagged (``gs<gs>``) + layout-hash-tagged so a rerun with a different
    group size or a regenerated layout writes to a fresh directory instead of reusing
    stale group HTs (compute_coverage's ``_group_path`` role, expressed through freq's
    ``kind="group"`` resolver + ``suffix=``). The chunk shape itself is part of the
    layout hash.

    :param level: Merge-tree level (1-indexed); level N output feeds level N+1.
    :param group_idx: Group index within this level (zero-based).
    :param intervals_hash: Content hash of the chunk-intervals layout.
    :param merge_group_size: Chunk HTs per group-merge job (tree fan-in).
    :param test: Whether to use the test-scoped path.
    :param environment: Compute environment.
    :param suffix: Optional run suffix.
    :return: Per-group HT path.
    """
    tree = f"gs{merge_group_size}_L{level:02d}"
    return get_aou_freq_chunk_path(
        group_idx,
        kind="group",
        test=test,
        environment=environment,
        suffix=_combine_suffix(suffix, intervals_hash, tree),
    )


def _interval_from_list(t: list, reference_genome: str) -> hl.utils.Interval:
    """
    Reconstruct a locus interval from its chunk-intervals JSON serialization.

    :param t: ``[start_contig, start_pos, end_contig, end_pos, includes_start,
        includes_end]`` (the layout ``compute_coverage.py`` writes).
    :param reference_genome: Reference-genome name (e.g. "GRCh38").
    :return: Locus interval.
    """
    sc, sp, ec, ep, incs, ince = t
    return hl.Interval(
        hl.Locus(sc, sp, reference_genome=reference_genome),
        hl.Locus(ec, ep, reference_genome=reference_genome),
        includes_start=incs,
        includes_end=ince,
    )


# ===========================================================================
# 2. Relay skeleton (shared submit machinery + eligibility)
# ===========================================================================


class _RelayJobSpec(NamedTuple):
    """One relay job's per-job config for :func:`_submit_relay_batch`."""

    name: str
    cpu: float
    memory: str
    storage: str
    attempts: int
    command: str


def _submit_relay_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    batch_name: str,
    job_specs: list[_RelayJobSpec],
    log_label: str,
) -> None:
    """
    Build and submit one Hail Batch of relay jobs sharing the same config.

    Shared skeleton for the chunk and merge submitters: each relay job is a NON-SPOT
    coordinator container (preemption mid-wait would orphan its inner QoB job) pinned to
    ``BATCH_REGIONS`` and run with a single attempt (no retry -- an OOM re-run only
    re-crashes, and non-spot means no preemption to recover from; failed chunks are
    collected into a manifest instead). Parallelism comes from Hail Batch's own
    scheduler running the N jobs concurrently. No-ops (skips ``batch.run()``)
    when ``job_specs`` is empty.

    :param args: Parsed CLI args (reads ``batch_image``, ``batch_dry_run``).
    :param backend_kwargs: kwargs for the ``hb.ServiceBackend(...)`` constructor.
    :param batch_name: Hail Batch name.
    :param job_specs: Per-job config (name + sizing + retry count + command).
    :param log_label: Noun for log messages ("chunk" / "merge").
    :return: None.
    """
    if not job_specs:
        logger.info(
            "  no pending %s jobs for %s; skipping batch.run()", log_label, batch_name
        )
        return

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
            # Relay is a coordinator waiting on its inner QoB job; preemption mid-wait
            # orphans that inner job, so relays are non-spot.
            j.spot(False)
            j.n_max_attempts(spec.attempts)
            j.command(spec.command)

        logger.info(
            "Submitting Hail Batch '%s': %d %s jobs (dry_run=%s)",
            batch_name,
            len(job_specs),
            log_label,
            args.batch_dry_run,
        )
        batch.run(dry_run=args.batch_dry_run)
    finally:
        backend.close()


def _build_freq_relay_common_flags(args: argparse.Namespace, *, chunk: bool) -> str:
    """
    Build the CLI flag string shared by per-chunk / per-merge relay invocations.

    Translates the orchestrator's ``args`` into a CLI string appended to each relay's
    ``--run-chunk`` / ``--run-merge`` command. The shared base (environment / AN
    environment / billing / tmp-dir / app-name / nested-QoB driver sizing) is common to
    both. With ``chunk=True`` the read/compute flags a chunk relay also needs
    (``--read-subintervals``, ``--all-sites-an-suffix``, ``--aou-freq-ht-suffix``,
    ``--chrom``, ``--test-vds`` / ``--test-partitions``) are appended; the merge relay
    (which only unions HTs) omits them. Per-job flags (``--chunk-*`` / ``--merge-*``)
    are added by the submit helpers. Flag order is irrelevant to argparse.

    :param args: Parsed CLI args.
    :param chunk: If True, include the chunk-only read/compute flags; if False, build
        the (smaller) merge relay flag set.
    :return: Space-joined ``--flag value`` string.
    """
    flags = [
        "--environment batch",
        f"--an-environment {args.an_environment}",
        f"--billing-project {args.billing_project}",
        f"--tmp-dir-days {args.tmp_dir_days}",
    ]
    if args.app_name:
        flags.append(f"--app-name {args.app_name}")
    # Nested-QoB driver/worker sizing (decoupled from the orchestrator's own
    # --driver-*/--worker-*).
    if args.chunk_driver_cores is not None:
        flags.append(f"--driver-cores {args.chunk_driver_cores}")
    if args.chunk_driver_memory:
        flags.append(f"--driver-memory {args.chunk_driver_memory}")
    if args.chunk_worker_memory:
        flags.append(f"--worker-memory {args.chunk_worker_memory}")
    if args.chunk_worker_cores:
        flags.append(f"--worker-cores {args.chunk_worker_cores}")
    if chunk:
        flags.append(f"--read-subintervals {args.read_subintervals}")
        if args.test_region:
            # The region IS this chunk's read intervals (no layout JSON); the worker
            # resolves them from --test-region and namespaces outputs as "test_region".
            flags.append(f"--test-region {' '.join(args.test_region)}")
        if args.all_sites_an_suffix:
            flags.append(f"--all-sites-an-suffix {args.all_sites_an_suffix}")
        if args.aou_freq_ht_suffix:
            flags.append(f"--aou-freq-ht-suffix {args.aou_freq_ht_suffix}")
        if args.chrom:
            flags.append(f"--chrom {args.chrom}")
        if args.test_vds:
            flags.append("--test-vds")
        if args.test_partitions is not None:
            flags.append(f"--test-partitions {args.test_partitions}")
    return " ".join(flags)


def _eligible_freq_chunk_indices(
    args: argparse.Namespace,
) -> tuple[list[str | None], list[int], str]:
    """
    Enumerate fan-out chunks and the subset selected by ``--chrom``.

    Single source of truth shared by the chunk orchestrator and the merge so they always
    agree on which chunks exist. Chunks come from the compute_coverage chunk layout
    beside the all-sites-AN HT (:func:`_load_freq_chunk_intervals`): it is read
    driver-side (no VDS access -- which is VPC-SC-blocked outside the perimeter) to get
    each chunk's contig, then only chunks whose contig matches ``--chrom`` are kept (all
    chunks when ``--chrom`` unset).

    :param args: Parsed CLI args (reads ``an_environment``, ``test``, ``chrom``).
    :return: ``(chunk_contigs, eligible, intervals_hash)`` -- the per-chunk contig indexed
        by chunk index, the list of eligible chunk indices after the ``--chrom`` filter,
        and the content hash of the chunk-intervals layout (namespaces chunk outputs).
    """
    if args.test_region:
        # A --test-region run is a single chunk whose intervals ARE the region (no JSON);
        # "test_region" namespaces its outputs, matching compute_coverage. Do not combine
        # with --chrom (the region already scopes the contig).
        chunk_contigs: list[str | None] = [None]
        intervals_hash = "test_region"
    else:
        data, _path, intervals_hash = _load_freq_chunk_intervals(
            args.an_environment, args.test
        )
        chunk_contigs = [c["contig"] for c in data["chunks"]]
    n_chunks = len(chunk_contigs)
    if args.chrom:
        eligible = [i for i in range(n_chunks) if chunk_contigs[i] == args.chrom]
        if not eligible:
            raise ValueError(
                f"No chunks match --chrom {args.chrom}; contigs in the precompute:"
                f" {sorted(set(chunk_contigs))}."
            )
    else:
        eligible = list(range(n_chunks))
    return chunk_contigs, eligible, intervals_hash


def _submit_freq_chunk_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    chunk_indices: list[int],
    intervals_hash: str,
    setup_cmd: str,
    common_flags_str: str,
    script: str,
    suffix: str | None,
    wave_label: str | None = None,
    n_chunks: int | None = None,
) -> None:
    """
    Build and submit one Hail Batch containing all pending chunk jobs.

    Each chunk job is a relay container that runs ``--run-chunk`` (which does
    ``hl.init(backend="batch")`` and spawns its own QoB driver). Parallelism comes from
    Hail Batch's own scheduler. Existence checks happen in the orchestrator before this
    is called (``chunk_indices`` is the already-filtered pending set).

    :param args: Parsed CLI args (reads ``chunk_cpu/memory/storage``, ``batch_image``,
        ``batch_dry_run``).
    :param backend_kwargs: kwargs for the per-call ``hb.ServiceBackend(...)``.
    :param chunk_indices: Pending chunk indices to submit.
    :param intervals_hash: Content hash of the chunk-intervals layout, used to namespace
        each chunk's output path.
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags from
        ``_build_freq_relay_common_flags(args, chunk=True)``.
    :param script: Path to the script inside the relay container.
    :param suffix: Run suffix folded into the chunk output path.
    :param wave_label: Optional suffix appended to the batch name per wave.
    :param n_chunks: Total chunks in the layout (batch name only).
    :return: None.
    """
    scope = "region" if args.test_region else f"{len(chunk_indices)}of{n_chunks}c"
    batch_name = f"v5_freq_aou_chunk_{scope}_{intervals_hash[:8]}"
    if suffix:
        batch_name += f"_{suffix}"
    if wave_label:
        batch_name += f"_{wave_label}"

    job_specs = []
    for idx in chunk_indices:
        path = _freq_chunk_path(
            idx,
            intervals_hash,
            test=args.test,
            environment=args.environment,
            suffix=suffix,
        )
        # Chunk identity is the chunk INDEX: the worker looks itself up in the per-chunk
        # JSON by --chunk-start (= idx); --chunk-stop is idx+1. The VDS read is
        # interval-based (read_intervals from the chunk's sub-intervals).
        command = (
            f"{setup_cmd}{script} --run-chunk"
            f" --chunk-start {idx} --chunk-stop {idx + 1}"
            f" --chunk-output {path}"
            f" {common_flags_str}"
        )
        job_specs.append(
            _RelayJobSpec(
                name=f"freq_chunk_{idx:06d}",
                cpu=args.chunk_cpu,
                memory=args.chunk_memory,
                storage=args.chunk_storage,
                attempts=1,
                # No retry: an OOM re-run just re-crashes; failures go to the manifest.
                command=command,
            )
        )
    _submit_relay_batch(args, backend_kwargs, batch_name, job_specs, "chunk")


def _freq_failed_chunks_path(sample_chunk_path: str) -> str:
    """
    Return the path to the fan-out's failed-chunk manifest.

    Co-located in the freq_chunks directory (``_failed_chunks.json``) so it is
    namespaced by the same intervals-hash + run suffix embedded in the chunk paths;
    a ``--chrom``/``--test``/suffix-scoped run therefore gets its own manifest and
    never collides with another.

    :param sample_chunk_path: Any per-chunk HT path (e.g. ``_freq_chunk_path(0, ...)``);
        its parent directory is the freq_chunks directory the manifest lives in.
    :return: GCS path to the failed-chunk manifest JSON.
    """
    chunk_dir = sample_chunk_path.rstrip("/").rsplit("/", 1)[0]
    return f"{chunk_dir}/_failed_chunks.json"


def _write_failed_chunks_manifest(
    args: argparse.Namespace,
    suffix: str | None,
    intervals_hash: str,
    failed_indices: list[int],
) -> str | None:
    """
    Write (or clear) the manifest of chunks that finished without a ``_SUCCESS``.

    Relays run with a single attempt (an OOM re-run just re-crashes; see
    ``_submit_relay_batch``), so failures are collected and recorded rather than retried
    in place. The manifest is self-describing -- each failed chunk's index, contig, and
    read sub-intervals -- so the user can eyeball what failed and, if needed, bump
    ``--chunk-driver-memory`` and rerun just those via ``--rerun-failed``. On a run with
    no failures any stale manifest is removed so it never lies about the latest state.

    :param args: Parsed CLI args (reads ``environment``, ``test``, ``chrom``).
    :param suffix: Run suffix folded into chunk output paths.
    :param intervals_hash: Content hash namespacing the chunk layout.
    :param failed_indices: Chunk indices missing ``_SUCCESS`` after the fan-out.
    :return: The manifest path if one was written, else None.
    """
    sample = _freq_chunk_path(
        0, intervals_hash, test=args.test, environment=args.environment, suffix=suffix
    )
    manifest_path = _freq_failed_chunks_path(sample)
    if not failed_indices:
        if file_exists(manifest_path):
            hfs.remove(manifest_path)
        return None
    if intervals_hash == "test_region":
        # A --test-region run has no layout JSON; the region itself is the chunk.
        records = [
            {"index": i, "contig": None, "sub_intervals": args.test_region}
            for i in sorted(failed_indices)
        ]
    else:
        chunks = _load_freq_chunk_intervals(args.an_environment, args.test)[0]["chunks"]
        records = [
            {
                "index": i,
                "contig": chunks[i]["contig"],
                "sub_intervals": chunks[i]["sub_intervals"],
            }
            for i in sorted(failed_indices)
        ]
    manifest = {
        "intervals_hash": intervals_hash,
        "suffix": suffix,
        "chrom": args.chrom,
        "n_failed": len(records),
        "failed_indices": [r["index"] for r in records],
        "chunks": records,
    }
    with hfs.open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)
    return manifest_path


def _read_failed_chunks_manifest(path: str) -> set[int]:
    """
    Read a failed-chunk manifest into a set of chunk indices for ``--rerun-failed``.

    Accepts either a manifest written by ``_write_failed_chunks_manifest`` (a JSON object
    with ``failed_indices``) or a hand-curated plain-text file of one integer index per
    line (``#`` comments allowed), so a user can trim the auto-written list to a subset.

    :param path: Path to the manifest (JSON or newline-delimited ints).
    :return: Set of chunk indices to reprocess.
    """
    if not file_exists(path):
        raise FileNotFoundError(f"--rerun-failed manifest not found: {path}")
    with hfs.open(path) as f:
        raw = f.read()
    if isinstance(raw, bytes):
        raw = raw.decode("utf-8")
    try:
        return set(json.loads(raw)["failed_indices"])
    except (json.JSONDecodeError, KeyError, TypeError):
        indices: set[int] = set()
        for line in raw.splitlines():
            line = line.split("#", 1)[0].strip()
            if line:
                indices.add(int(line))
        return indices


def _orchestrate_freq_batch(
    args: argparse.Namespace,
    setup_cmd: str,
    backend_kwargs: dict,
    script: str,
    suffix: str | None,
) -> None:
    """
    Fan the all-sites-AN freq compute out as relay chunk jobs (one job per chunk).

    Pending chunks are split into sequential **waves** of ``--wave-size`` chunks. Each
    wave is its own Hail Batch (one relay job per chunk; each relay spawns its own QoB
    driver via ``hl.init(backend="batch")``) and runs to completion before the next is
    submitted -- which (a) bounds the number of concurrently-running relays and their
    nested QoB drivers to ``--wave-size``, and (b) avoids the process-global rich.Live
    crash from concurrent ``batch.run()`` calls. Within a wave, Hail Batch's own
    scheduler runs the relays in parallel.

    Idempotency: chunks whose ``_SUCCESS`` already exists are skipped (ONE directory
    listing), so a failed/partial wave is resumed by simply rerunning this step. Merge is
    a separate step (``--merge-freq-chunks``).

    :param args: Parsed CLI args (reads ``wave_size``, ``overwrite``, chunk sizing, ...).
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param backend_kwargs: kwargs for the per-call ``hb.ServiceBackend(...)``.
    :param script: Path to the script inside the relay container.
    :param suffix: Run suffix folded into chunk output paths.
    :return: None.
    """
    chunk_contigs, eligible, intervals_hash = _eligible_freq_chunk_indices(args)
    n_chunks = len(chunk_contigs)

    if args.overwrite:
        pending_indices = list(eligible)
    else:
        # ONE directory listing instead of per-chunk serial existence probes. Reuse
        # freq's _list_present_freq_chunk_indices (globs *.freq.chunk_*.ht/_SUCCESS).
        present = _list_present_freq_chunk_indices(
            _freq_chunk_path(
                0,
                intervals_hash,
                test=args.test,
                environment=args.environment,
                suffix=suffix,
            )
        )
        pending_indices = [idx for idx in eligible if idx not in present]
    if args.rerun_failed:
        # Reprocess only the chunks listed in a failed-chunk manifest (intersected with
        # what is still missing _SUCCESS). Lets the user bump --chunk-driver-memory and
        # rerun just the OOM chunks from a prior run.
        rerun_set = _read_failed_chunks_manifest(args.rerun_failed)
        before = len(pending_indices)
        pending_indices = [idx for idx in pending_indices if idx in rerun_set]
        logger.info(
            "--rerun-failed %s: restricted %d pending -> %d listed chunk(s).",
            args.rerun_failed,
            before,
            len(pending_indices),
        )
    logger.info(
        "Freq fan-out: %d chunks total, %d eligible%s, %d pending, %d skipped"
        " (overwrite=%s)",
        n_chunks,
        len(eligible),
        f" (--chrom {args.chrom})" if args.chrom else "",
        len(pending_indices),
        len(eligible) - len(pending_indices),
        args.overwrite,
    )
    if not pending_indices:
        logger.info("All chunks already complete; nothing to submit.")
        return

    common_flags_str = _build_freq_relay_common_flags(args, chunk=True)

    wave_size = args.wave_size
    if wave_size <= 0 or wave_size >= len(pending_indices):
        waves = [pending_indices]
    else:
        waves = [
            pending_indices[i : i + wave_size]
            for i in range(0, len(pending_indices), wave_size)
        ]
    n_waves = len(waves)
    logger.info(
        "Dispatching %d pending chunks in %d sequential wave(s) of up to %d each.",
        len(pending_indices),
        n_waves,
        wave_size if wave_size > 0 else len(pending_indices),
    )

    all_failed: list[int] = []
    for wi, wave_indices in enumerate(waves, start=1):
        wave_label = f"w{wi:03d}of{n_waves:03d}" if n_waves > 1 else None
        logger.info(
            "Wave %d/%d: submitting %d chunks (indices %d..%d).",
            wi,
            n_waves,
            len(wave_indices),
            wave_indices[0],
            wave_indices[-1],
        )
        _submit_freq_chunk_batch(
            args=args,
            backend_kwargs=backend_kwargs,
            chunk_indices=wave_indices,
            intervals_hash=intervals_hash,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
            suffix=suffix,
            wave_label=wave_label,
            n_chunks=n_chunks,
        )
        # A dry run submits nothing, so an output check would spuriously flag every
        # chunk as failed (and write a manifest); skip it.
        if args.batch_dry_run:
            continue
        # batch.run() does not raise on per-job failure; re-check this wave's outputs
        # (one listing) and surface any that did not land (rather than only discovering
        # them at --merge-freq-chunks time).
        present = _list_present_freq_chunk_indices(
            _freq_chunk_path(
                0,
                intervals_hash,
                test=args.test,
                environment=args.environment,
                suffix=suffix,
            )
        )
        failed = [idx for idx in wave_indices if idx not in present]
        all_failed.extend(failed)
        if failed:
            logger.warning(
                "Wave %d/%d complete but %d/%d chunk(s) MISSING after run (single"
                " attempt, no retry); missing indices: %s%s",
                wi,
                n_waves,
                len(failed),
                len(wave_indices),
                failed[:25],
                " ..." if len(failed) > 25 else "",
            )
        else:
            logger.info(
                "Wave %d/%d complete; all %d chunks present.",
                wi,
                n_waves,
                len(wave_indices),
            )

    if args.batch_dry_run:
        logger.info("Dry run: submitted no jobs; skipping the failed-chunk manifest.")
        return
    manifest_path = _write_failed_chunks_manifest(
        args, suffix, intervals_hash, all_failed
    )
    if all_failed:
        logger.warning(
            "Fan-out finished: %d/%d dispatched chunk(s) missing _SUCCESS. Wrote"
            " failed-chunk manifest %s. Rerun with --rerun-failed %s (optionally"
            " --chunk-driver-memory highmem) to reprocess just these.",
            len(all_failed),
            len(pending_indices),
            manifest_path,
            manifest_path,
        )
    else:
        logger.info(
            "Fan-out finished: all %d dispatched chunk(s) present.",
            len(pending_indices),
        )


def _orchestrate_freq_list_failed(
    args: argparse.Namespace,
    suffix: str | None,
) -> None:
    """
    Scan the fan-out outputs and (re)write the failed-chunk manifest, no jobs submitted.

    A standalone, Hail-free counterpart to the manifest the fan-out writes itself:
    globs ``_SUCCESS`` once (``_list_present_freq_chunk_indices``) against the eligible
    set (``_eligible_freq_chunk_indices``) and records every chunk still missing. Use it
    to recover the failed list when the orchestrator process itself died mid-run (and so
    never wrote the manifest), or just to re-check status before ``--rerun-failed`` /
    ``--merge-freq-chunks``.

    :param args: Parsed CLI args (reads ``environment``, ``test``, ``chrom``).
    :param suffix: Run suffix folded into chunk output paths.
    :return: None.
    """
    _chunk_contigs, eligible, intervals_hash = _eligible_freq_chunk_indices(args)
    present = _list_present_freq_chunk_indices(
        _freq_chunk_path(
            0,
            intervals_hash,
            test=args.test,
            environment=args.environment,
            suffix=suffix,
        )
    )
    failed = [idx for idx in eligible if idx not in present]
    manifest_path = _write_failed_chunks_manifest(args, suffix, intervals_hash, failed)
    if failed:
        logger.warning(
            "%d/%d eligible chunk(s) missing _SUCCESS; wrote manifest %s.",
            len(failed),
            len(eligible),
            manifest_path,
        )
    else:
        logger.info(
            "All %d eligible chunk(s) present; no failed-chunk manifest written.",
            len(eligible),
        )


def _orchestrate_freq_fanout(
    args: argparse.Namespace,
    suffix: str | None,
) -> None:
    """
    Top-level orchestrator: layout check -> chunk fan-out.

    Never initializes Hail (its checks/listing use ``hailtop.fs``): builds the setup
    command + backend kwargs and delegates to ``_orchestrate_freq_batch``. The final
    merged HT is produced by the separate ``--merge-freq-chunks`` step.

    :param args: Parsed CLI args (reads ``environment``, ``test``, ``methods_branch``,
        ``billing_project``, ``batch_remote_tmpdir``, ``batch_image``, plus everything
        ``_orchestrate_freq_batch`` reads).
    :param suffix: Run suffix (``--aou-freq-ht-suffix`` folded with ``--chrom``) used to
        namespace chunk outputs; must match the merge step.
    :return: None.
    """
    # Fail fast (before submitting jobs) if the coverage chunk layout is missing -- the
    # fan-out and merge enumerate chunks from it. A --test-region run has no JSON (the
    # region is the single chunk), so skip the check there.
    if not args.test_region:
        _load_freq_chunk_intervals(args.an_environment, args.test)

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    # Reuse freq's _build_setup_command (targets v5_freq_batch:0.2.137, no
    # hail reinstall).
    setup_cmd = _build_setup_command(commit, methods_branch=args.methods_branch)

    backend_kwargs = {"billing_project": args.billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir

    script = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/generate_frequency.py"
    _orchestrate_freq_batch(args, setup_cmd, backend_kwargs, script, suffix)


# ===========================================================================
# 3. Chunk worker (all-sites-AN under a NESTED QoB driver)
# ===========================================================================


def _run_aou_freq_chunk_all_sites_ans(args: argparse.Namespace) -> None:
    """
    Compute the all-sites-AN AoU freq HT for a SINGLE chunk and write it.

    Runs under a nested QoB driver inside a relay container the orchestrator submitted.

    Steps:
      1. Init Hail with the BATCH (QoB) backend -- ``hl.init(backend="batch")`` via
         ``_initialize_hail`` -- NOT local Spark.
      2. Look this chunk up by index (``--chunk-start``) in the compute_coverage chunk
         layout beside the all-sites-AN HT, read driver-side (``hailtop.fs``, no QoB
         job) so we skip a full VDS re-open and the query-on-batch cold-start a
         Hail-Table read would incur. A chunk never spans a contig boundary.
      3. Read the AoU VDS via ``read_intervals=sub_intervals`` (one partition per
         sub-interval, no shuffle), with metadata + release filter.
      4. ``_prepare_aou_vds`` (ploidy adjust -> adj -> split; strata globals).
      5. ``_calculate_aou_frequencies_and_hists_using_all_sites_ans(..., region_intervals=
         sub_intervals, ...)`` -- passing ``region_intervals`` so the all-sites-AN HT read
         is pruned to the chunk. (This helper already annotates ``mt_hist_fields`` and
         drops ``raw_qual_hists`` internally, exactly as the single-job path does.)
      6. ``select_final_dataset_fields(ht, dataset="aou")``; stamp the layout hash into a
         global for provenance; write to ``--chunk-output``.

    :param args: Parsed CLI args (reads ``chunk_start``, ``chunk_output``, ``environment``,
        ``an_environment``, ``test_vds``, ``test_partitions``, ``chrom``,
        ``all_sites_an_suffix``, ``aou_freq_ht_suffix``).
    :return: None.
    """
    # (1) Nested QoB init. _initialize_hail derives the batch app name and calls
    # _init_hail(backend="batch").
    _initialize_hail(args)

    environment = args.environment
    test = (
        args.test_vds
        or args.test_partitions is not None
        or args.test_region is not None
    )
    start = args.chunk_start

    # (2) Resolve this chunk's read sub-intervals. An explicit --test-region IS the single
    # chunk -- its intervals are the region, "test_region" namespaces its outputs, and no
    # layout JSON is read (matching compute_coverage's --test-region layout). This is
    # what makes a region test land where a region-scoped AN table actually has data.
    # Otherwise look the chunk up by index in the coverage chunk layout, read
    # driver-side (hailtop.fs, no QoB job) so we skip a full VDS re-open.
    if args.test_region:
        intervals_hash = "test_region"
        # Split into read_subintervals so the read lands in that many partitions,
        # matching the fanned-out chunks -- a single-partition read collects the whole
        # region's aggregation on the driver and OOMs it.
        sub_intervals = _split_intervals_for_read(
            [_parse_region_interval(r) for r in args.test_region],
            args.read_subintervals,
        )
        logger.info(
            "Split --test-region into %d read sub-intervals for chunk %d (no precompute"
            " JSON).",
            len(sub_intervals),
            start,
        )
    else:
        data, intervals_path, intervals_hash = _load_freq_chunk_intervals(
            args.an_environment, test
        )
        chunk_meta = data["chunks"]
        if not 0 <= start < len(chunk_meta):
            raise ValueError(
                f"chunk index {start} is out of range [0, {len(chunk_meta)}) in"
                f" {intervals_path}; the fan-out and the coverage chunk layout are out"
                " of sync."
            )
        entry = chunk_meta[start]
        rg = data["reference_genome"]
        sub_intervals = [_interval_from_list(t, rg) for t in entry["sub_intervals"]]
        logger.info(
            "Read %d sub-intervals for chunk %d (contig %s) from the precompute.",
            len(sub_intervals),
            start,
            entry["contig"],
        )

    # Auto-derive --chunk-output if omitted (manual single-chunk run); the orchestrator
    # always passes it explicitly. Namespaced by the layout hash so a stale chunk is
    # never skipped-as-present or merged in.
    output_path = args.chunk_output
    if output_path is None:
        output_path = _freq_chunk_path(
            start,
            intervals_hash,
            test=test,
            environment=environment,
            suffix=_combine_suffix(args.aou_freq_ht_suffix, args.chrom),
        )
        logger.info("Auto-derived --chunk-output: %s", output_path)

    # (3) Read the AoU VDS pruned to this chunk's sub-intervals (no straddle-widening --
    # all-sites-AN is variant-only), columns filtered to the release samples with no
    # per-row passes (see _load_release_aou_vds); sex_karyotype / age come from the
    # small meta table joined in _prepare_aou_vds.
    vds = _load_release_aou_vds(
        environment,
        test_vds=args.test_vds,
        read_intervals=sub_intervals,
        chrom=args.chrom,
    )

    # (4) Prepare (ploidy adjust -> adj -> split, LGT-preserving). A chunk never
    # spans a contig, so an autosomal chunk skips the (identity) ploidy adjustment.
    vmt = _prepare_aou_vds(
        vds,
        test=test,
        environment=environment,
        skip_sex_ploidy=not _spans_sex_chromosome(sub_intervals),
    )

    # (5) All-sites-AN compute, with the AN HT read pruned to this chunk's intervals.
    freq_ht = _calculate_aou_frequencies_and_hists_using_all_sites_ans(
        vmt,
        test=test,
        environment=environment,
        chrom=args.chrom,
        region_intervals=sub_intervals,
        all_sites_an_suffix=args.all_sites_an_suffix,
        an_environment=args.an_environment,
    )

    # (6) Final field select, provenance stamp, write.
    freq_ht = select_final_dataset_fields(freq_ht, dataset="aou")
    freq_ht = freq_ht.annotate_globals(freq_chunk_intervals_hash=intervals_hash)
    freq_ht.write(output_path, overwrite=True)
    logger.info("Wrote freq chunk %d to %s", start, output_path)


# ===========================================================================
# 4. Merge (tree-reduce union of chunk HTs)
# ===========================================================================


def _run_freq_merge(
    input_paths: list[str],
    output_path: str,
    coalesce_to: int | None = None,
) -> None:
    """
    Union a list of partial AoU freq HTs and write the result.

    Globals are identical across all inputs (built from the same group_membership HT), so
    the union inherits them from the first HT. Used for both group-level and final merges
    in the ``--merge-freq-chunks`` pipeline.

    :param input_paths: GCS paths of the partial HTs to union.
    :param output_path: GCS path to write the merged HT.
    :param coalesce_to: If set, ``naive_coalesce`` to this many partitions before writing
        (one output partition per input for group merges; ``--n-partitions`` for the
        final merge). Default None (natural sum-of-input partition count).
    :return: None.
    """
    logger.info(
        "Merging %d HTs -> %s (coalesce_to=%s)",
        len(input_paths),
        output_path,
        coalesce_to,
    )
    hts = [hl.read_table(p) for p in input_paths]
    merged = hl.Table.union(*hts) if len(hts) > 1 else hts[0]
    if coalesce_to is not None:
        merged = merged.naive_coalesce(coalesce_to)
    merged.write(output_path, overwrite=True)
    logger.info("Wrote merged HT to %s", output_path)


def _submit_freq_merge_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    group_indices: list[int],
    groups: list[list[str]],
    group_output_paths: list[str],
    setup_cmd: str,
    common_flags_str: str,
    script: str,
    level: int,
    suffix: str | None,
) -> None:
    """
    Build and submit one Hail Batch containing all pending group-merge jobs.

    Each job runs ``--run-merge`` over its assigned inputs and writes an intermediate
    per-group HT (via ``_freq_group_path``). Only intermediate levels go through this
    function; the final union is submitted separately by ``_orchestrate_freq_merge`` and
    is the sole job that writes the canonical final freq HT. Per-group coalesce target is
    the number of inputs in that group (one output partition per input).

    :param args: Parsed CLI args (reads ``merge_cpu/memory/storage``, ``batch_image``,
        ``batch_dry_run``).
    :param backend_kwargs: kwargs for the per-call ``hb.ServiceBackend(...)``.
    :param group_indices: Pending group indices to submit.
    :param groups: For each group index, the list of input HT paths to union.
    :param group_output_paths: For each group index, the output HT path.
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags from
        ``_build_freq_relay_common_flags(args, chunk=False)``.
    :param script: Path to the script inside the relay container.
    :param level: Merge-tree level (1-indexed); used in the batch and job names.
    :param suffix: Run suffix (batch-name only).
    :return: None.
    """
    batch_name = f"v5_freq_merge_L{level:02d}"
    if suffix:
        batch_name += f"_{suffix}"

    job_specs = []
    for group_idx in group_indices:
        group_inputs = groups[group_idx]
        command = (
            f"{setup_cmd}{script} --run-merge"
            f" --merge-output {group_output_paths[group_idx]}"
            f" --merge-coalesce-to {len(group_inputs)}"
            f" --merge-inputs {' '.join(group_inputs)}"
            f" {common_flags_str}"
        )
        job_specs.append(
            _RelayJobSpec(
                name=f"freq_merge_L{level:02d}_{group_idx:06d}",
                cpu=args.merge_cpu,
                memory=args.merge_memory,
                storage=args.merge_storage,
                attempts=1,
                # No retry: an OOM re-run just re-crashes; failures go to the manifest.
                command=command,
            )
        )
    _submit_relay_batch(args, backend_kwargs, batch_name, job_specs, "merge")


def _orchestrate_freq_merge(
    args: argparse.Namespace,
    final_output_path: str,
    suffix: str | None,
) -> None:
    """
    Recursive tree-reduce merge of per-chunk HTs into ``final_output_path``.

    Counterpart to ``_orchestrate_freq_fanout``: runs AFTER the fan-out has produced all
    per-chunk HTs. Submits Hail Batch jobs (QoB-from-container per job) and exits without
    initializing Hail in this process.

    Discovery uses the same ``_eligible_freq_chunk_indices`` the fan-out uses (so the two
    agree), filtered by ``--chrom``. Every expected chunk must have a ``_SUCCESS`` marker;
    missing chunks fail loudly so the user re-runs the fan-out before merging.

    Recursive tree: from N chunks, each level groups inputs into windows of
    ``--merge-group-size`` and emits one ``--run-merge`` job per group; level-k>1 inputs
    are level-(k-1) outputs. Iteration stops when <= one group remains; that group is the
    final-merge job that writes ``final_output_path``. Safe to re-run: group HTs whose
    ``_SUCCESS`` exists are skipped; the final HT is skipped if it exists (unless
    ``--overwrite``).

    :param args: Parsed CLI args (reads ``merge_group_size``, ``overwrite``,
        ``n_partitions``, merge sizing, billing/tmpdir/image, ...).
    :param final_output_path: GCS path for the final merged AoU freq HT.
    :param suffix: Run suffix used to namespace group/chunk HTs; must match the fan-out.
    :return: None.
    """
    chunk_contigs, eligible, intervals_hash = _eligible_freq_chunk_indices(args)
    n_chunks = len(chunk_contigs)
    logger.info(
        "Verifying %d expected chunk HTs exist (of %d total)...",
        len(eligible),
        n_chunks,
    )
    present = _list_present_freq_chunk_indices(
        _freq_chunk_path(
            0,
            intervals_hash,
            test=args.test,
            environment=args.environment,
            suffix=suffix,
        )
    )
    missing = [i for i in eligible if i not in present]
    if missing:
        raise FileNotFoundError(
            f"--merge-freq-chunks: {len(missing)} of {len(eligible)} expected chunks"
            f" missing (first few idx: {missing[:5]}). Run the fan-out to (re)compute"
            " missing chunks first."
        )
    logger.info("All %d expected chunks present.", len(eligible))

    gs = args.merge_group_size

    # Precompute the level shape so we can log the full plan upfront.
    shape = [len(eligible)]
    while shape[-1] > gs:
        shape.append((shape[-1] + gs - 1) // gs)
    logger.info(
        "Merge tree (group_size=%d): %s -> 1 final HT (%d intermediate level(s))",
        gs,
        " -> ".join(str(n) for n in shape),
        len(shape) - 1,
    )

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    setup_cmd = _build_setup_command(commit, methods_branch=args.methods_branch)

    backend_kwargs = {"billing_project": args.billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir

    script = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/generate_frequency.py"
    common_flags_str = _build_freq_relay_common_flags(args, chunk=False)

    # Intermediate levels: each emits a level-tagged group HT per output. Iterates while
    # #inputs > gs; stops when one final merge can union everything remaining.
    inputs = [
        _freq_chunk_path(
            i,
            intervals_hash,
            test=args.test,
            environment=args.environment,
            suffix=suffix,
        )
        for i in eligible
    ]
    level = 1
    while len(inputs) > gs:
        n_in = len(inputs)
        n_out = (n_in + gs - 1) // gs
        groups = [inputs[i : i + gs] for i in range(0, n_in, gs)]
        out_paths = [
            _freq_group_path(
                level,
                idx,
                intervals_hash,
                args.merge_group_size,
                test=args.test,
                environment=args.environment,
                suffix=suffix,
            )
            for idx in range(n_out)
        ]

        if args.overwrite:
            pending = list(range(n_out))
        else:
            pending = []
            for idx in range(n_out):
                if file_exists(out_paths[idx]):
                    logger.info(
                        "Skipping already-complete L%d group %d at %s",
                        level,
                        idx,
                        out_paths[idx],
                    )
                else:
                    pending.append(idx)
        logger.info(
            "Level %d dispatch: %d groups total, %d pending, %d skipped (%d -> %d)",
            level,
            n_out,
            len(pending),
            n_out - len(pending),
            n_in,
            n_out,
        )
        # _submit_freq_merge_batch no-ops when there are no pending merges.
        _submit_freq_merge_batch(
            args=args,
            backend_kwargs=backend_kwargs,
            group_indices=pending,
            groups=groups,
            group_output_paths=out_paths,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
            level=level,
            suffix=suffix,
        )
        inputs = out_paths
        level += 1

    # Final merge: a single job that unions the remaining (<= gs) inputs and writes the
    # canonical output.
    if not args.overwrite and file_exists(final_output_path):
        logger.info(
            "Final merge HT exists at %s; skipping (pass --overwrite to rewrite).",
            final_output_path,
        )
        return

    final_batch_name = "v5_freq_merge_final"
    if suffix:
        final_batch_name += f"_{suffix}"
    coalesce_flag = (
        f" --merge-coalesce-to {args.n_partitions}"
        if args.n_partitions is not None
        else ""
    )
    logger.info("Final merge: %d inputs -> %s", len(inputs), final_output_path)
    final_spec = _RelayJobSpec(
        name="freq_merge_final",
        cpu=args.merge_cpu,
        memory=args.merge_memory,
        storage=args.final_merge_storage,
        attempts=1,
        # No retry: an OOM re-run just re-crashes; failures go to the manifest.
        command=(
            f"{setup_cmd}{script} --run-merge"
            f" --merge-output {final_output_path}"
            f"{coalesce_flag}"
            f" --merge-inputs {' '.join(inputs)}"
            f" {common_flags_str}"
        ),
    )
    _submit_relay_batch(
        args, backend_kwargs, final_batch_name, [final_spec], "final-merge"
    )


def main(args):
    """Generate v5 frequency data."""
    # The all-sites-AN HT (and the coverage chunk layout beside it) may live in another
    # bucket than the rest of the run; default to the run's own environment.
    args.an_environment = args.an_environment or args.environment

    # --- Batch worker subcommands (run inside Hail Batch containers) ---
    if args.run_chunk:
        args.test = (
            args.test_vds
            or args.test_partitions is not None
            or args.test_region is not None
        )
        # Relay chunk worker: nested QoB, reads its own sub-intervals from the
        # coverage chunk layout and initializes Hail itself.
        _run_aou_freq_chunk_all_sites_ans(args)
        return

    if args.run_merge:
        # Relay tree-merge worker (union + optional coalesce).
        _initialize_hail(args)
        _run_freq_merge(
            input_paths=args.merge_inputs,
            output_path=args.merge_output,
            coalesce_to=args.merge_coalesce_to,
        )
        return

    # --- Normal orchestrator flow ---
    environment = args.environment
    test_vds = args.test_vds
    test_partitions = args.test_partitions
    overwrite = args.overwrite
    aou_freq_suffix = args.aou_freq_ht_suffix
    chrom = args.chrom
    # A --chrom run writes to a contig-scoped output path so per-contig runs don't
    # collide; assemble them later with --assemble-chrom-freq. The merge-datasets step
    # reads the assembled (contig-unscoped) AoU HT, so it keeps ``aou_freq_suffix``.
    aou_freq_suffix_chrom = _combine_freq_suffix(aou_freq_suffix, chrom)
    # Output-only tag appended to whichever freq HT this run writes (process-aou /
    # process-gnomad / merge-datasets). Not applied to any input read, so a tagged test
    # output never feeds --merge-datasets.
    freq_output_suffix = args.freq_output_suffix
    aou_out_suffix = _combine_freq_suffix(aou_freq_suffix_chrom, freq_output_suffix)
    gnomad_out_suffix = _combine_freq_suffix(
        _combine_freq_suffix(None, chrom), freq_output_suffix
    )
    tmp_dir_days = args.tmp_dir_days
    # --test-region is a testing-only scope, so it auto-enables test mode (test paths).
    test_run = test_vds or test_partitions is not None or args.test_region is not None

    _initialize_hail(args)

    try:
        logger.info("Running generate_frequency.py...")

        if args.assemble_chrom_freq:
            # Union the per-contig freq HTs (each written by a separate
            # --process-{aou,gnomad} --chrom <contig> run) into the canonical
            # (contig-unscoped) freq HT. Globals are identical across contigs (same
            # group_membership HT), so the union inherits them from the first input.
            data_set = args.assemble_data_set
            base_suffix = aou_freq_suffix if data_set == "aou" else None
            canonical_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set=data_set,
                environment=environment,
                suffix=base_suffix,
            )
            per_contig_paths = [
                get_freq(
                    test=test_run,
                    data_type="genomes",
                    data_set=data_set,
                    environment=environment,
                    suffix=_combine_freq_suffix(base_suffix, c),
                )
                for c in args.contigs
            ]
            check_resource_existence(
                input_step_resources={"per-contig-freq": per_contig_paths},
                output_step_resources={"assemble-chrom-freq": [canonical_freq]},
                overwrite=overwrite,
            )
            logger.info(
                "Assembling %d per-contig %s freq HTs into %s...",
                len(per_contig_paths),
                data_set,
                canonical_freq.path,
            )
            hts = [p.ht() for p in per_contig_paths]
            assembled_ht = hl.Table.union(*hts) if len(hts) > 1 else hts[0]
            if args.n_partitions is not None:
                assembled_ht = assembled_ht.naive_coalesce(args.n_partitions)
            assembled_ht.write(canonical_freq.path, overwrite=overwrite)

        if args.process_gnomad:
            logger.info("Processing gnomAD dataset...")

            gnomad_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set="gnomad",
                environment=environment,
                suffix=gnomad_out_suffix,
            )

            check_resource_existence(
                output_step_resources={"process-gnomad": [gnomad_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = process_gnomad_dataset(
                test_vds=test_vds,
                test_partitions=test_partitions,
                environment=environment,
                chrom=chrom,
            )

            logger.info(
                "Writing gnomAD frequency HT (with embedded age histograms) to %s...",
                gnomad_freq.path,
            )
            gnomad_freq_ht.write(gnomad_freq.path, overwrite=overwrite)

        if args.process_aou:
            logger.info("Processing All of Us dataset...")
            aou_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set="aou",
                environment=environment,
                suffix=aou_out_suffix,
            )

            if args.list_failed_chunks:
                # Scan fan-out outputs and (re)write the failed-chunk manifest;
                # no jobs submitted, writes no HT -- so no output-existence gate.
                args.test = test_run
                _orchestrate_freq_list_failed(args, aou_freq_suffix_chrom)
                return

            if args.use_batch_fanout and not args.merge_freq_chunks:
                # Crash-resilient relay fan-out: one non-spot coordinator per chunk
                # runs the all-sites-AN compute via nested QoB, so a single chunk crash
                # can't wreck the whole run. It writes per-chunk HTs (NOT the final
                # freq HT), so it is not gated on the final output existing -- that
                # would block resuming a run whose merge already wrote a final. Assemble
                # with --merge-freq-chunks after this completes.
                args.test = test_run
                _orchestrate_freq_fanout(args, aou_freq_suffix_chrom)
                return

            # From here on every path writes the final freq HT, so gate on it.
            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq]},
                overwrite=overwrite,
            )

            if args.merge_freq_chunks:
                # Tree-reduce the per-chunk HTs produced by the relay fan-out
                # into the final freq HT (run AFTER the fan-out finishes).
                args.test = test_run
                _orchestrate_freq_merge(args, aou_freq.path, aou_freq_suffix_chrom)
                return

            # Single-job all-sites-AN path (test regions / small scopes).
            aou_freq_ht = process_aou_dataset(
                test_vds=test_vds,
                test_partitions=test_partitions,
                environment=environment,
                chrom=chrom,
                test_region=args.test_region,
                all_sites_an_suffix=args.all_sites_an_suffix,
                read_subintervals=args.read_subintervals,
                an_environment=args.an_environment,
            )
            logger.info("Writing AoU frequency HT to %s...", aou_freq.path)
            aou_freq_ht.write(aou_freq.path, overwrite=overwrite)

        if args.merge_datasets:
            logger.info(
                "Merging frequency data and age histograms from both datasets..."
            )

            merged_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set="merged",
                environment=environment,
                suffix=freq_output_suffix,
            )

            check_resource_existence(
                output_step_resources={"merge-datasets": [merged_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = get_freq(
                data_type="genomes",
                test=test_run,
                data_set="gnomad",
                environment=environment,
            )
            aou_freq_ht = get_freq(
                data_type="genomes",
                test=test_run,
                data_set="aou",
                environment=environment,
                suffix=aou_freq_suffix,
            )

            check_resource_existence(
                input_step_resources={
                    "process-gnomad": [gnomad_freq_ht],
                    "process-aou": [aou_freq_ht],
                }
            )
            merged_freq_ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq_ht.ht(),
                aou_freq_ht.ht(),
            )
            merged_freq_ht = merged_freq_ht.checkpoint(
                new_temp_file("merged_freq", "ht")
            )

            logger.info(
                "Calculating FAF, grpmax, and other annotations on merged dataset..."
            )
            merged_freq_ht = calculate_faf_and_grpmax_annotations(merged_freq_ht)

            merged_freq_ht = select_final_dataset_fields(
                merged_freq_ht, dataset="merged"
            )

            logger.info("Writing merged frequency HT to %s...", merged_freq.path)
            merged_freq_ht.write(merged_freq.path, overwrite=overwrite)

    finally:
        # Skip in batch mode: the JVM runs remotely so there is no local log
        # file to copy, and Hail Batch already retains the full driver/worker
        # logs (accessible via the Batch UI or `hailctl batch log`).
        if environment != "batch":
            hl.copy_log(
                get_logging_path(
                    "v5_frequency_run",
                    environment=environment,
                    tmp_dir_days=tmp_dir_days,
                )
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate frequency data for gnomAD v5."
    )

    # General arguments.
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--aou-freq-ht-suffix",
        type=str,
        default=None,
        help=(
            "Optional filename suffix inserted before the '.ht' extension on "
            "the AoU frequency HT, e.g. 'split_vds' -> "
            "'...frequencies.split_vds.ht'. Useful for tagging test / "
            "experimental runs so they don't collide with the default output "
            "path. Applied to both --process-aou and --merge-datasets reads "
            "of the AoU freq HT so the merge picks up the same file. Default "
            "is None (no suffix)."
        ),
    )
    parser.add_argument(
        "--all-sites-an-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix (no leading underscore) inserted before the '.ht' "
            "extension of the all-sites-AN HT read on the --use-all-sites-ans path, "
            "to target a suffixed AN HT written by 'compute_coverage.py "
            "--cov-and-an-output-suffix' (e.g. 'chunkhash_test' -> "
            "'...coverage_and_an_chunkhash_test.ht'). Uses the same underscore "
            "convention as compute_coverage. Default is None (base coverage_and_an.ht)."
        ),
    )
    parser.add_argument(
        "--freq-output-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix appended to the WRITTEN freq HT for whichever step runs "
            "(--process-aou / --process-gnomad / --merge-datasets), e.g. "
            "'...frequencies.chunkhash_test.ht'. Output-only: unlike "
            "--aou-freq-ht-suffix it does NOT change any input read, so a tagged test "
            "output is a dead-end artifact and never feeds --merge-datasets. Combines "
            "with --aou-freq-ht-suffix and --chrom on the output path. Default is None."
        ),
    )

    # Test/debug arguments.
    test_group = parser.add_argument_group("testing options")
    test_group.add_argument(
        "--test-vds",
        help="Use the test VDS for given project.",
        action="store_true",
    )
    test_group.add_argument(
        "--test-partitions",
        type=int,
        default=None,
        help=(
            "Filter the VDS to this many partitions for testing. When supplied,"
            " enables test mode. The all-sites-AN fan-out ignores the value itself"
            " (its chunks come from the compute_coverage chunk layout, which is already"
            " test-scoped) and uses the flag only to select test paths. Default is None"
            " (process all partitions)."
        ),
    )
    test_group.add_argument(
        "--chrom",
        default=None,
        help=(
            "Filter input data to a single contig (e.g. --chrom chr22) and append it"
            " to the output freq HT suffix, so a per-contig run goes all the way"
            " through processing (fan-out -> merge) without clobbering another"
            " contig -- a late failure only loses that contig, reducing the chance of"
            " a catastrophic whole-genome loss. Applies to both the AoU and gnomAD"
            " paths. Assemble the per-contig AoU HTs with --assemble-chrom-freq."
            " Independent of --test-*."
        ),
    )
    test_group.add_argument(
        "--test-region",
        nargs="+",
        default=None,
        help=(
            "TESTING ONLY: scope the AoU read to these explicit locus interval(s), e.g."
            " --test-region chr1:55058666-55108666 chr1:55108666-55158666. Parsed as"
            " half-open [start, end) intervals (same as compute_coverage.py's"
            " --test-region), so the frequency output can be matched against an"
            " all-sites-AN test HT generated over the same region in QoB testing; the"
            " AoU VDS and the all-sites-AN HT are both filtered to these intervals."
            " Auto-enables test mode (reads the test all-sites-AN HT and writes test"
            " output paths) while still using the real AoU VDS scoped to the region."
            " The all-sites-AN QoB test region was chr1:55058666-55158666. Mutually"
            " exclusive with --test-partitions."
        ),
    )
    test_group.add_argument(
        "--read-subintervals",
        type=int,
        default=48,
        help=(
            "PERF: split each --test-region into this many read sub-intervals so the"
            " variant read lands in this many partitions (N-way row parallelism for the"
            " all-sites-AN aggregation) instead of the single partition a whole-region read"
            " prunes to. A single partition serializes the ~245k-sample aggregation into"
            " one task -- slow, and it OOMs a standard driver -- so this is on by default."
            " Position-based tiling; only affects the --use-all-sites-ans path when a"
            " --test-region is set (whole-genome / --chrom runs already parallelize via the"
            " VDS's native partitions). ~48 was the measured cost/runtime knee for a"
            " ~70k-variant region and the exact count barely matters past it. Pass 1 to"
            " disable. Default 48."
        ),
    )

    # Processing step arguments.
    processing_group = parser.add_argument_group("processing steps")
    processing_group.add_argument(
        "--process-gnomad",
        help="Process gnomAD dataset for frequency calculations.",
        action="store_true",
    )
    processing_group.add_argument(
        "--process-aou",
        help=(
            "Process the AoU dataset for frequency calculations: AC, homozygote counts,"
            " and histograms are aggregated from the AoU VDS variant data and joined to"
            " the all-sites AN written by compute_coverage.py (no densify). Runs as a"
            " single job, or as a Hail Batch relay fan-out with --use-batch-fanout."
        ),
        action="store_true",
    )
    processing_group.add_argument(
        "--merge-datasets",
        help="Merge frequency data from both gnomAD and AoU datasets.",
        action="store_true",
    )
    processing_group.add_argument(
        "--assemble-chrom-freq",
        help=(
            "Union the per-contig freq HTs (each written by a separate"
            " --process-{aou,gnomad} --chrom <contig> run) into the canonical freq HT"
            " for --assemble-data-set. Provide the contigs to union with --contigs."
            " Runs in-process (do not pass --chrom)."
        ),
        action="store_true",
    )
    processing_group.add_argument(
        "--assemble-data-set",
        choices=["aou", "gnomad"],
        default="aou",
        help=(
            "Data set whose per-contig freq HTs --assemble-chrom-freq unions. Default"
            " 'aou'."
        ),
    )
    processing_group.add_argument(
        "--contigs",
        nargs="+",
        default=None,
        help=(
            "Contigs to union for --assemble-chrom-freq (e.g. --contigs chr1 chr2 ...)."
            " Each must have a completed per-contig freq HT. Required with"
            " --assemble-chrom-freq."
        ),
    )

    # Batch worker subcommands (invoked by Hail Batch jobs, not by users).
    worker_group = parser.add_argument_group(
        "batch worker subcommands",
        "Internal subcommands run by Hail Batch chunk/merge jobs.",
    )
    worker_group.add_argument(
        "--run-chunk",
        help=argparse.SUPPRESS,
        action="store_true",
    )
    worker_group.add_argument(
        "--chunk-start",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--chunk-stop",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--chunk-output",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--run-merge",
        help=argparse.SUPPRESS,
        action="store_true",
    )
    worker_group.add_argument(
        "--merge-inputs",
        nargs="+",
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--merge-output",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )

    # Environment configuration.
    env_group = parser.add_argument_group("environment configuration")
    env_group.add_argument(
        "--environment",
        help="Environment to run in.",
        choices=["rwb", "batch"],
        default="batch",
    )
    env_group.add_argument(
        "--an-environment",
        choices=["rwb", "batch", "dataproc"],
        default=None,
        help=(
            "Environment (bucket) the all-sites-AN HT and the compute_coverage chunk"
            " layout beside it are read from. Everything else (VDS, metadata, group"
            " membership, chunk outputs) stays on --environment. Pass 'dataproc' when"
            " compute_coverage ran with --results-environment dataproc. Default is the"
            " value of --environment."
        ),
    )
    env_group.add_argument(
        "--tmp-dir-days",
        type=int,
        default=4,
        help="Number of days for temp directory retention. Default is 4.",
    )
    env_group.add_argument(
        "--billing-project",
        type=str,
        default="gnomad-production",
        help="Google Cloud billing project for reading requester pays buckets.",
    )

    # Batch-specific configuration.
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
    )
    batch_group.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_group.add_argument(
        "--driver-cores",
        type=int,
        default=None,
        help="Number of cores for driver node.",
    )
    batch_group.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_group.add_argument(
        "--worker-cores",
        type=str,
        default=None,
        help=(
            "Cores per QoB worker job. Hail Batch requires 1, 2, 4, or 8 for JVM"
            " jobs (fractional cores are rejected with a 400)."
        ),
    )
    batch_group.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
    )

    # Batch fan-out configuration (--process-aou --use-batch-fanout).
    fanout_group = parser.add_argument_group(
        "batch fan-out configuration",
        "Parameters for the --process-aou Hail Batch relay fan-out and its merge.",
    )
    fanout_group.add_argument(
        "--wave-size",
        type=int,
        default=1000,
        help=(
            "All-sites-AN relay fan-out: chunks per sequential wave (bounds concurrent"
            " relays). <=0 submits a single batch. Default 1000."
        ),
    )
    fanout_group.add_argument(
        "--chunk-driver-cores",
        type=int,
        default=1,
        help=(
            "Relay: nested-QoB driver cores per chunk. The driver only builds the IR"
            " and waits on the workers, so 1 core suffices. Default 1."
        ),
    )
    fanout_group.add_argument(
        "--chunk-driver-memory",
        type=str,
        default="highmem",
        help=(
            "Relay: nested-QoB driver memory per chunk. 1-core 'highmem' (6.5 GB) is"
            " the measured minimum: 1-core 'standard' (4 GB) OOMs while lowering the"
            " fused compute IR. Default 'highmem'."
        ),
    )
    fanout_group.add_argument(
        "--chunk-worker-memory",
        type=str,
        default=None,
        help=(
            "Relay: nested-QoB worker memory class per chunk ('standard'/'highmem')."
            " Chunk workers peak under 200 MB, but Hail Batch rejects 'lowmem' for JVM"
            " jobs. Default None (Hail's default, 'standard')."
        ),
    )
    fanout_group.add_argument(
        "--chunk-worker-cores",
        type=str,
        default=None,
        help=(
            "Relay: cores per nested-QoB worker job. Hail Batch requires 1, 2, 4, or 8"
            " for JVM jobs (a fractional value is rejected with a 400; measured"
            " 2026-09-30). The chunk tasks are CPU-bound, so more cores per worker"
            " do not lower cost. Default None (Hail's default, 1)."
        ),
    )
    fanout_group.add_argument(
        "--list-failed-chunks",
        action="store_true",
        help=(
            "All-sites-AN relay fan-out: scan the fan-out outputs and (re)write the"
            " failed-chunk manifest (chunks missing _SUCCESS). No jobs submitted; use to"
            " recover the failed list if the orchestrator process died mid-run."
        ),
    )
    fanout_group.add_argument(
        "--rerun-failed",
        type=str,
        default=None,
        help=(
            "All-sites-AN relay fan-out: path to a failed-chunk manifest (from a prior"
            " fan-out or --list-failed-chunks). Restricts the fan-out to just those"
            " chunks -- pair with --chunk-driver-memory highmem to reprocess OOM chunks."
        ),
    )
    fanout_group.add_argument(
        "--merge-freq-chunks",
        action="store_true",
        help=(
            "All-sites-AN relay fan-out: tree-reduce the per-chunk HTs into the final"
            " freq HT (run AFTER the fan-out finishes)."
        ),
    )
    fanout_group.add_argument(
        "--merge-coalesce-to",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    fanout_group.add_argument(
        "--n-partitions",
        type=int,
        default=None,
        help=(
            "Number of partitions the final merged freq HT is coalesced to"
            " (--merge-freq-chunks final merge and --assemble-chrom-freq). The fan-out"
            " takes its chunks from the compute_coverage chunk layout, not from this."
            " Default None (no coalesce)."
        ),
    )
    fanout_group.add_argument(
        "--use-batch-fanout",
        help=(
            "Run --process-aou as a Hail Batch relay fan-out (one non-spot relay per"
            " chunk of the compute_coverage layout, each driving a nested QoB job)"
            " instead of a single job. Writes per-chunk HTs; assemble them with"
            " --merge-freq-chunks."
        ),
        action="store_true",
    )
    fanout_group.add_argument(
        "--merge-group-size",
        type=int,
        default=500,
        help="Number of chunk HTs unioned per group-merge job. Default 500.",
    )
    fanout_group.add_argument(
        "--chunk-cpu",
        type=float,
        default=0.5,
        help=(
            "CPU request per chunk relay job; fractional values are allowed -- the"
            " relay is a QoB client that builds the query and waits on the nested"
            " batch, and 0.5 core is the measured working size. Default 0.5."
        ),
    )
    fanout_group.add_argument(
        "--chunk-memory",
        type=str,
        default="standard",
        help="Memory preset per chunk relay job. Default 'standard'.",
    )
    fanout_group.add_argument(
        "--chunk-storage",
        type=str,
        default="25Gi",
        help="Extra /io storage per chunk job. Default '25Gi'.",
    )
    fanout_group.add_argument(
        "--merge-cpu",
        type=int,
        default=4,
        help="CPU request per merge job. Default 4.",
    )
    fanout_group.add_argument(
        "--merge-memory",
        type=str,
        default="standard",
        help="Memory preset for merge jobs. Default 'standard'.",
    )
    fanout_group.add_argument(
        "--merge-storage",
        type=str,
        default="50Gi",
        help="Extra storage per group-merge job. Default '50Gi'.",
    )
    fanout_group.add_argument(
        "--final-merge-storage",
        type=str,
        default="100Gi",
        help="Extra storage for the final merge job. Default '100Gi'.",
    )
    fanout_group.add_argument(
        "--batch-image",
        type=str,
        default=DEFAULT_BATCH_IMAGE,
        help=(
            "Docker image for chunk and merge jobs. Defaults to the"
            " v5_freq_batch:0.2.137 image (Hail 0.2.137 + gnomad_methods deps baked"
            " in). Pass another image to change the Hail version -- the version comes"
            " from the image, not a runtime reinstall."
        ),
    )
    fanout_group.add_argument(
        "--batch-remote-tmpdir",
        type=str,
        default=None,
        help=(
            "gs:// path used by the Hail Batch ServiceBackend for its"
            " scratch. Defaults to the standard batch tmp bucket."
        ),
    )
    fanout_group.add_argument(
        "--methods-branch",
        type=str,
        default="main",
        help=(
            "Branch or commit of gnomad_methods to pull at runtime." " Default 'main'."
        ),
    )
    fanout_group.add_argument(
        "--batch-dry-run",
        help="Validate the fan-out DAG without actually running jobs.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    batch_args = [
        "app_name",
        "driver_cores",
        "driver_memory",
        "worker_cores",
        "worker_memory",
    ]
    provided_batch_args = [arg for arg in batch_args if getattr(args, arg) is not None]
    if provided_batch_args and args.environment != "batch":
        parser.error(
            f"Batch configuration arguments ({', '.join('--' + a.replace('_', '-') for a in provided_batch_args)}) "
            f"require --environment=batch"
        )

    # --test-region and --test-partitions are two different test-scoping mechanisms
    # (explicit intervals vs a partition slice); combining them is contradictory.
    if args.test_region is not None and args.test_partitions is not None:
        parser.error("--test-region and --test-partitions are mutually exclusive.")

    # --assemble-chrom-freq unions the per-contig HTs into the canonical
    # (contig-unscoped) path, so it needs the contig list and must not itself be
    # contig-scoped.
    if args.assemble_chrom_freq:
        if not args.contigs:
            parser.error("--assemble-chrom-freq requires --contigs.")
        if args.chrom:
            parser.error(
                "--assemble-chrom-freq unions every per-contig HT into the canonical"
                " (contig-unscoped) path; do not pass --chrom."
            )

    main(args)
