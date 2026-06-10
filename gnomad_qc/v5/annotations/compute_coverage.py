r"""
Compute coverage, allele number, and quality histograms for gnomAD v5 genomes.

v5 genomes = AoU v8 (new) + gnomAD v4 genomes minus the consent-drop set.
This script computes per-reference-site coverage, AN, and (AoU only) qual
histograms for each project track and joins them into the v5 release HT/TSV.

Execution roles:
----------------
One entry-point (this script), three process roles dispatched by ``main``
on the CLI flags. They are mutually exclusive (the ``__main__`` block
validates this):

1. ORCHESTRATOR (``--use-batch-fanout`` or ``--merge-cov-chunks``):
   runs wherever you launch it, does not call ``_init_hail`` (but the
   ``_SUCCESS``/output skip-checks via ``_file_exists_for_env`` do
   initialize the Hail backend), builds and submits a Hail Batch of
   relay jobs, then returns. Used for the prod-scale AoU compute (the
   dataset is too big for one job).
2. WORKER (``--run-chunk`` or ``--run-merge``): each relay container the
   orchestrator submitted re-invokes THIS script with ``--run-chunk`` /
   ``--run-merge``. Initializes Hail, does exactly one chunk/merge unit,
   and returns BEFORE the try/finally. The chunk worker is what actually
   calls ``compute_all_release_stats_per_ref_site``.
3. IN-PROCESS PIPELINE (none of the above flags): initializes Hail and
   runs the try/finally step chain. Reached only when not orchestrating
   or working — i.e. the strict single-job compute (gnomAD consent-drop,
   or dev/test AoU) and the AoU release assembly (the ``--export-*`` /
   ``--merge-qual-hists`` steps, run as a separate invocation after
   fan-out + ``--merge-cov-chunks``).

So a full prod AoU release is three *separate* invocations: (1) fan-out
compute, (2) merge-cov-chunks, (3) the in-process assembly invocation.
The group-membership / downsampling HTs are written by their own
in-process invocation first; fan-out then reads them from disk.

Processing Workflow:
--------------------
Per-project setup (--project-name {aou,gnomad})::

    1. (AoU only) Write the downsampling HT (--write-aou-downsampling-ht).
    2. Write the group membership HT (--write-group-membership-ht): strata =
       sex_karyotype x genetic_ancestry x (downsampling for AoU). Use
       --reduce-min-aggs to write the leaf-only variant.
    2b. Write the preprocessed vep_context sites HT (--write-vep-context-sites):
        deduped/telomere-stripped/locus-keyed sites read co-partitioned by every
        chunk (so the dedup + strip runs once, not per chunk).
    3. Compute the per-ref-site coverage/AN/qual-hists HT — the dense step.
       Strict single-job via --compute-all-cov-release-stats-ht, or the
       prod-scale fan-out (--use-batch-fanout then --merge-cov-chunks). See
       "Execution roles" above for how those dispatch.
    3b. (full runs) Validate the merged HT covers every vep_context site
        (--validate-cov-and-an) — anti-join, fail on any dropped site.

Per-project merge (--project-name gnomad)::

    4. Build the gnomAD v5 coverage HT (--merge-gnomad-coverage) and AN HT
       (--merge-gnomad-an) by subtracting the consent-drop samples from the
       v4 release HTs.

Release assembly (--project-name aou)::

    5. Join AoU + gnomAD v5 outputs and export release files:
       --export-coverage-release-files (coverage HT + TSV)
       --export-an-release-files (AN HT + TSV)
       --merge-qual-hists (qual hists HT; gnomAD v4 hists are reused as-is —
       they were not recomputed for v5).

Usage Examples::

    # Per-chunk fan-out compute + merge on Hail Batch (the prod path):
    python compute_coverage.py --project-name aou --environment batch \\
        --use-batch-fanout --partitions-per-chunk 3
    python compute_coverage.py --project-name aou --environment batch \\
        --merge-cov-chunks --n-partitions 10000

    # Strict whole-VDS compute on rwb:
    python compute_coverage.py --project-name aou --environment rwb \\
        --compute-all-cov-release-stats-ht

    # gnomAD v4 consent-drop subtraction + AoU/gnomAD merge:
    python compute_coverage.py --project-name gnomad --environment dataproc \\
        --merge-gnomad-coverage --merge-gnomad-an
    python compute_coverage.py --project-name aou --environment batch \\
        --export-coverage-release-files --export-an-release-files \\
        --merge-qual-hists
"""

import argparse
import logging
import subprocess
from functools import reduce
from typing import List, NamedTuple, Optional, Sequence, Tuple

import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.gnomad import CURRENT_GENOME_AN_RELEASE as v4_AN_RELEASE
from gnomad.resources.grch38.gnomad import (
    CURRENT_GENOME_COVERAGE_RELEASE as v4_COVERAGE_RELEASE,
)
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, public_release
from gnomad.resources.grch38.reference_data import (
    telomeres_and_centromeres,
    vep_context,
)
from gnomad.utils.annotations import (
    annotate_downsamplings,
    build_freq_stratification_list,
    generate_freq_group_membership_array,
    merge_array_expressions,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.file_utils import repartition_for_join
from gnomad.utils.sparse_mt import (
    compute_stats_per_ref_site,
    get_allele_number_agg_func,
    get_coverage_agg_func,
)
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_aou_downsampling,
    group_membership,
    qual_hists,
)
from gnomad_qc.v5.resources.basics import (
    _get_batch_resource_kwargs,
    _init_hail,
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.meta import meta
from gnomad_qc.v5.resources.release import (
    release_all_sites_an_tsv_path,
    release_coverage_path,
    release_coverage_tsv_path,
)

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("v5_coverage_and_an")
logger.setLevel(logging.INFO)

COVERAGE_OVER_X_BINS = (1, 5, 10, 15, 20, 25, 30, 50, 100)

# Overall groups that non-reducible aggregations (coverage_stats →
# {"group": "adj"}, qual_hists → {"group": "raw"}) are pinned to via
# `entry_agg_group_membership`. Under `--reduce-min-aggs` these must be
# retained as leaves in the reduced group_membership HT so they stay
# directly computed (a cheap index lookup) instead of reconstructed
# per-row from their leaf-children, which is what made the reduced AoU
# run more expensive than the non-reduced one.
PINNED_LEAF_GROUPS = [{"group": "adj"}, {"group": "raw"}]

# All Hail Batch jobs (chunk relays + merges) are pinned to the region the
# input data lives in: the AoU VDS, vep_context, and outputs are all
# us-central1. Co-locating the workers with the data avoids inter-region GCS
# egress, and since the relays are non-spot (j.spot(False)) there's no
# spot-capacity reason to spread across regions.
BATCH_REGIONS = ["us-central1"]


def get_downsampling_ht(ht: hl.Table) -> hl.Table:
    """
    Get Table with downsampling groups for all samples.

    v5 downsampling is only applied to the AoU dataset.
    Desired groups:
    - 10,000
    - 100,000
    - Genetic ancestry group sizes for AFR, AMR, NFE

    Only AFR, AMR, and NFE per-group downsamplings are generated: the gen_anc
    expression is restricted to those groups before ``annotate_downsamplings``,
    which adds each present group's sample count to the downsamplings list and a
    per-group index. The global 10k/100k still cover all samples (they use a
    gen-anc-independent global index).

    :param ht: Input Table.
    :return: Table with downsampling groups.
    """
    logger.info(
        "Determining downsampling groups for AoU...",
    )
    downsamplings = DOWNSAMPLINGS["v5"]
    # Restrict per-group downsamplings to the desired genetic ancestry groups so
    # annotate_downsamplings doesn't generate them for every other group too.
    gen_anc = ht.genetic_ancestry_inference.gen_anc
    gen_anc = hl.or_missing(
        hl.literal({"afr", "amr", "nfe"}).contains(gen_anc), gen_anc
    )
    ht = annotate_downsamplings(ht, downsamplings, gen_anc)
    return ht


def get_group_membership_ht(
    meta_ht: hl.Table,
    project: str,
    ds_ht: Optional[hl.Table] = None,
    reduce_min_aggs: bool = False,
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    :param meta_ht: Meta HT.
    :param project: Project name. Must be "aou" or "gnomad". If "gnomad", function will filter meta HT to only consent drop samples.
    :param ds_ht: Optional downsampling HT. Only used for AoU.
    :param reduce_min_aggs: If True, build the HT with `reduce_to_minimal_groups=True`. Downstream `compute_stats_per_ref_site` calls must then list every non-reducible annotation under `entry_agg_group_membership` and pass `reducible_aggs={"AN"}` (or similar) so AN gets expanded back to full strata.
    :return: Group membership HT.
    """
    if project == "aou":
        ht = generate_freq_group_membership_array(
            meta_ht,
            build_freq_stratification_list(
                sex_expr=meta_ht.sex_karyotype,
                gen_anc_expr=meta_ht.genetic_ancestry_inference.gen_anc,
                downsampling_expr=ds_ht[meta_ht.key].downsampling,
            ),
            reduce_to_minimal_groups=reduce_min_aggs,
            force_leaf_groups=PINNED_LEAF_GROUPS,
            downsamplings=hl.eval(ds_ht.downsamplings),
            ds_gen_anc_counts=hl.eval(ds_ht.ds_gen_anc_counts),
        )

        ht = ht.annotate_globals(
            freq_meta=ht.freq_meta.map(
                lambda d: hl.dict(
                    d.items().map(
                        lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x)
                    )
                )
            ),
        )

    elif project == "gnomad":
        # Filter to v4 consent drop samples.
        # NOTE: Not using v5 project meta here because this part will be run in
        # Dataproc.
        ht = meta_ht.filter(
            meta_ht.release
            & (
                (meta_ht.project_meta.research_project_key == "RP-1061")
                | (meta_ht.project_meta.research_project_key == "RP-1411")
            )
        )
        ht = generate_freq_group_membership_array(
            ht,
            build_freq_stratification_list(
                sex_expr=ht.sex_imputation.sex_karyotype,
                gen_anc_expr=ht.population_inference.pop,
            ),
            reduce_to_minimal_groups=reduce_min_aggs,
            force_leaf_groups=PINNED_LEAF_GROUPS,
        )
    else:
        raise ValueError(f"project must be 'aou' or 'gnomad', got {project!r}")

    return ht


def _file_exists_for_env(path: str, environment: str) -> bool:
    """
    Check if a path exists, tolerant of permission errors in batch mode.

    On the batch backend, anonymous file probes against requester-pays
    buckets can raise permission errors before getting to "exists / does
    not exist." Treat those *specifically* as "exists" so the chunk is
    skipped rather than re-run; the next stage will surface a real error
    if the file is actually broken. Any other exception (network,
    asyncio, etc.) is re-raised so we fail loud instead of silently
    skipping work.

    :param path: GCS path to probe.
    :param environment: Compute environment. Only "batch" activates
        the permission-error fallback; other environments propagate.
    :return: True if the file exists (or is assumed to under batch
        permission errors), False otherwise.
    """
    from gnomad.utils.file_utils import file_exists

    try:
        return file_exists(path)
    except PermissionError as e:
        # This PermissionError was isolated to Hail 0.2.138's anonymous probe of
        # requester-pays paths. Treating it as "exists" is a workaround, but it is
        # dangerous: a genuine creds/bucket misconfig would mark a never-written
        # chunk as present and let the merge run over missing data. Logged at ERROR
        # so a broad perms problem (many "assumed exists") is visible, not silent.
        if environment == "batch":
            logger.error(
                "file_exists on %s raised PermissionError (%s); ASSUMING EXISTS"
                " (Hail 0.2.138 workaround). If many of these appear, a real"
                " permissions problem may be silently skipping unwritten chunks.",
                path,
                e,
            )
            return True
        raise


def _chunk_path(cov_and_an_ht_path: str, idx: int) -> str:
    """
    Return a sibling per-chunk HT path under ``<cov_and_an_path>_chunks/``.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param idx: Chunk index (zero-based).
    :return: ``<cov_and_an_path>_chunks/<idx:08d>.chunk.ht``.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_chunks/{idx:08d}.chunk.ht"


def _group_path(cov_and_an_ht_path: str, level: int, group_idx: int) -> str:
    """
    Return a per-group merged HT path under ``<cov_and_an_path>_merge_groups/``.

    Level-tagged so a recursive merge tree (multiple levels of grouping
    when chunk count exceeds ``merge_group_size`` once-over) doesn't
    overwrite earlier-level outputs.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param level: Merge-tree level (1-indexed); level N output is the
        input to level N+1.
    :param group_idx: Group index within this level (zero-based).
    :return: ``<cov_and_an_path>_merge_groups/L<level:02d>_<idx:08d>.ht``.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_merge_groups/L{level:02d}_{group_idx:08d}.ht"


def _apply_path_suffix(path: str, suffix: Optional[str]) -> str:
    """
    Insert ``_<suffix>`` before the ``.ht`` extension, or return unchanged if no suffix.

    :param path: HT path ending in ``.ht``.
    :param suffix: Optional suffix string (no leading underscore). If
        falsy, ``path`` is returned unchanged.
    :return: Suffix-applied path.
    """
    if not suffix:
        return path
    return path.rstrip("/").removesuffix(".ht") + f"_{suffix}.ht"


def _resolve_cov_and_an_ht_path(
    project: str,
    environment: str,
    test: bool,
    suffix: Optional[str],
) -> str:
    """
    Return the cov_and_an HT path, applying ``suffix`` when set.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param test: Whether to route to the test path.
    :param suffix: Optional suffix appended before the ``.ht`` extension
        (for A/B comparison runs).
    :return: Fully resolved cov_and_an HT path.
    """
    path = coverage_and_an_path(
        test=test, data_set=project, environment=environment
    ).path
    return _apply_path_suffix(path, suffix)


def _resolve_group_membership_ht_path(
    project: str,
    environment: str,
    test: bool,
    reduce_min_aggs: bool,
) -> str:
    """
    Return the group_membership HT path, applying ``_reduce.ht`` under ``reduce_min_aggs``.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param test: Whether to route to the test path.
    :param reduce_min_aggs: If True, append ``_reduce`` to the path —
        used to keep the leaf-only-reduction HT separate from the
        full-strata HT.
    :return: Fully resolved group_membership HT path.
    """
    path = group_membership(test=test, data_set=project, environment=environment).path
    return _apply_path_suffix(path, "reduce" if reduce_min_aggs else None)


def _log_name_for_run(
    run_chunk: bool,
    run_merge: bool,
    chunk_start: int,
    chunk_stop: int,
    merge_output: Optional[str],
) -> str:
    """
    Build a per-worker log name so concurrent chunk/merge workers don't clobber the monolithic log.

    :param run_chunk: ``--run-chunk`` flag.
    :param run_merge: ``--run-merge`` flag.
    :param chunk_start: Chunk start partition index (used when
        ``run_chunk`` is True).
    :param chunk_stop: Chunk stop partition index (used when
        ``run_chunk`` is True).
    :param merge_output: Merge output HT path (used to derive a slug
        when ``run_merge`` is True).
    :return: Log name; suitable for use as a log-file basename.
    """
    if run_chunk:
        return f"v5_cov_chunk_{chunk_start:08d}_{chunk_stop:08d}"
    if run_merge:
        merge_slug = "merge"
        if merge_output:
            merge_slug = merge_output.rstrip("/").split("/")[-1].removesuffix(".ht")
        return f"v5_cov_{merge_slug}"
    return "v5_coverage_and_an_generation"


def _derive_chunk_locus_intervals(
    vds_filtered: hl.vds.VariantDataset,
    n_subdivisions: int = 1,
    reference_genome: str = "GRCh38",
) -> List[hl.utils.Interval]:
    """
    Derive per-contig locus sub-intervals covering the filtered VDS reference_data.

    The vep_context HT (used as ``ref_ht``) and the AoU VDS have independent
    partition layouts, so chunking each by its own partition index would not
    yield the same locus range. Instead, we derive the locus extent of the
    VDS chunk and use it to align ``ref_ht`` to the same range.

    When ``n_subdivisions > 1``, each contig's ``lo..hi`` extent is split into
    that many equal-position sub-intervals. The returned objects are concrete
    ``hl.Interval`` instances (NOT ``IntervalExpression`` s — ``read_vds`` /
    ``read_matrix_table``'s ``_intervals=`` reader path requires Python-native
    intervals so the IR can serialize them as partition boundaries). Pass to
    ``hl.vds.read_vds(intervals=...)`` and ``hl.read_table(_intervals=...)``
    to co-partition both sides of the densify join at read time, avoiding the
    shuffle that a post-read ``repartition`` would cost. With
    ``n_subdivisions=1`` (default), returns one big interval per contig.

    :param vds_filtered: VDS already filtered to the chunk's partitions
        (the locus extent is derived from its ``reference_data``).
    :param n_subdivisions: Number of equal-position sub-intervals per
        contig. Default 1 (one interval per contig).
    :param reference_genome: Reference genome name for the constructed
        ``hl.Locus`` objects. Default ``"GRCh38"``.
    :return: List of ``hl.Interval`` objects covering the VDS chunk.
    """
    rd = vds_filtered.reference_data
    bounds = rd.aggregate_rows(
        hl.agg.group_by(
            rd.locus.contig,
            hl.struct(
                lo=hl.agg.min(rd.locus.position),
                hi=hl.agg.max(rd.locus.position),
            ),
        )
    )
    n = max(n_subdivisions, 1)
    sub_intervals: List[hl.utils.Interval] = []
    for contig, b in bounds.items():
        lo, hi = b.lo, b.hi + 1
        total = max(hi - lo, 1)
        step = max(total // n, 1)
        for i in range(n):
            sub_lo = lo + i * step
            sub_hi = lo + (i + 1) * step if i < n - 1 else hi
            if sub_hi <= sub_lo:
                continue
            sub_intervals.append(
                hl.Interval(
                    hl.Locus(contig, sub_lo, reference_genome=reference_genome),
                    hl.Locus(contig, sub_hi, reference_genome=reference_genome),
                    includes_start=True,
                    includes_end=False,
                )
            )
    return sub_intervals


def _derive_ref_partition_intervals(
    n_partitions: int, chrom: Optional[List[str]] = None
) -> List[hl.utils.Interval]:
    """
    Derive balanced locus intervals from the vep_context sites table.

    Thin wrapper over ``gnomad.utils.file_utils.repartition_for_join`` (balanced
    ``_calculate_new_partitions`` intervals, per Hail-team guidance "read all the
    tables with the same partitions") used by the strict single-job path to
    co-partition the VDS and vep_context reads. Restricting to ``chrom`` first
    keeps the intervals aligned with a ``--chrom``-filtered run.

    :param n_partitions: Number of balanced partitions to derive when positive; a
        falsy value (0) uses ``repartition_for_join``'s default
        1.1x-of-native-partitions multiplier (no number to guess for the full run).
    :param chrom: Optional list of contigs to restrict vep_context to before
        deriving partitions, so the intervals match a ``--chrom``-filtered run.
    :return: List of bare-``Locus`` ``hl.Interval`` objects.
    """
    ref_ht = vep_context.versions["105"].ht().key_by("locus")
    if chrom:
        ref_ht = hl.filter_intervals(
            ref_ht, [hl.parse_locus_interval(c) for c in chrom]
        )
    return repartition_for_join(
        ref_ht,
        locus_intervals=True,
        **({"n_partitions": n_partitions} if n_partitions else {}),
    )


def compute_all_release_stats_per_ref_site(
    vds: hl.vds.VariantDataset,
    ref_ht: hl.Table,
    sex_karyotype_field: str,
    project: str,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
    interval_ht: Optional[hl.Table] = None,
    group_membership_ht: Optional[hl.Table] = None,
    reduce_min_aggs: bool = False,
) -> hl.Table:
    """
    Compute coverage, allele number, and quality histograms per reference site.

    .. note::

        Running this function prior to calculating frequencies removes the need for an additional
        densify for frequency calculations.

    :param vds: Input VDS. Reference data must carry ``END``/``GQ``/``DP``;
        a ``LEN`` field is added if missing.
    :param ref_ht: Locus-only sites Table (typically derived from
        ``vep_context`` and stripped of telomeres/centromeres) defining
        the reference positions at which to aggregate.
    :param sex_karyotype_field: Dotted path on the variant_data column
        struct to the sample's sex karyotype (e.g.
        ``"meta.sex_karyotype"``). Used by ``compute_stats_per_ref_site``
        to set per-sample ploidy on sex chromosomes.
    :param project: "aou" or "gnomad". When "aou", additionally fans out
        per-site qual histograms (``qual_hists``); "gnomad" computes only
        coverage and AN.
    :param coverage_over_x_bins: DP thresholds for the ``over_X``
        cumulative-sample-count fields written into the output HT.
        Default is :data:`COVERAGE_OVER_X_BINS`.
    :param interval_ht: Optional interval Table forwarded to
        ``compute_stats_per_ref_site`` for partition pruning. Unused on
        the v5 path.
    :param group_membership_ht: Group-membership Table built by
        ``get_group_membership_ht``. Defines the per-stratum sample sets
        AN is fanned out across.
    :param reduce_min_aggs: If True, pass ``reducible_aggs={"AN"}`` to
        ``compute_stats_per_ref_site`` so AN goes through the
        leaf-reduction path. Requires ``group_membership_ht`` to have
        been built with ``reduce_to_minimal_groups=True``.
    :return: HT keyed by locus with per-stratum ``AN``, flat
        ``mean``/``over_X``/``median_approx``/``total_DP`` coverage
        fields (from the global adj-filtered group), and
        ``qual_hists`` when ``project == "aou"``.
    """

    def _get_hists(qual_expr) -> hl.expr.Expression:
        """
        Build adj-only GQ + DP histograms from a [GQ, DP, adj] expression triple.

        Selecting ``.qual_hists`` (the adj-filtered struct) drops the unfiltered
        ``raw_qual_hists``: those aggregations are never referenced, so Hail prunes
        them. Raw qual hists are not loaded into the browser and were approved for
        removal from v5 (saves the extra per-site histogram aggregation).
        """
        return hl.struct(
            qual_hists=qual_hist_expr(
                gq_expr=qual_expr[0],
                dp_expr=qual_expr[1],
                adj_expr=qual_expr[2] == 1,
                split_adj_and_raw=True,
            ).qual_hists
        )

    # Set up coverage bins.
    cov_bins = sorted(coverage_over_x_bins)
    rev_cov_bins = list(reversed(cov_bins))
    max_cov_bin = cov_bins[-1]
    cov_bins = hl.array(cov_bins)

    entry_agg_funcs = {
        "AN": get_allele_number_agg_func("LGT"),
        "coverage_stats": get_coverage_agg_func(dp_field="DP", max_cov_bin=max_cov_bin),
    }
    # Restrict `coverage_stats` to the global adj-filtered group only —
    # downstream code uses `coverage_stats[0]` exclusively (the global adj
    # group; transmuted into flat mean/sum/over_X fields), so per-strata
    # coverage aggregation is wasted work. We use `{"group": "adj"}` (and
    # not `"raw"`) because `compute_stats_per_ref_site` does NOT rewrite
    # group labels to "raw" when a pre-built `group_membership_ht` is
    # passed (the rewrite only happens on the `strata_expr` code path), so
    # `freq_meta[0]` is `{"group": "adj"}` and we need to match it.
    # AN is intentionally omitted so it still fans out across all strata,
    # which is what downstream consumers need.
    entry_agg_group_membership = {
        "coverage_stats": [{"group": "adj"}],
    }
    # Only compute qual hists for AoU.
    if project == "aou":
        entry_agg_funcs["qual_hists"] = (lambda t: [t.GQ, t.DP, t.adj], _get_hists)

        # Below we use just the raw group for qual hist computations because
        # qual_hist_expr does its own adj filtering when adj is passed as an
        # argument; we keep only the adj-filtered histograms (see _get_hists --
        # raw qual hists were approved for removal from v5).
        entry_agg_group_membership["qual_hists"] = [{"group": "raw"}]

    logger.info(
        "Computing coverage, allele number, and optionally qual hists per reference site..."
    )

    vmt = vds.variant_data
    sex_expr = reduce(lambda x, field: x[field], sex_karyotype_field.split("."), vmt)
    vmt = vmt.annotate_cols(sex_karyotype=sex_expr)
    rmt = vds.reference_data
    # Computing LEN is only needed for VDS written before Hail 0.2.134, which
    # added LEN to the reference data; future VDS versions ship it already, so
    # this branch can be dropped once the inputs are all >= 0.2.134.
    # https://hail.is/docs/0.2/change_log.html#version-0-2-134
    if "LEN" not in rmt.entry:
        rmt = rmt.annotate_entries(LEN=rmt.END - rmt.locus.position + 1)
    vds = hl.vds.VariantDataset(rmt, vmt)

    ht = compute_stats_per_ref_site(
        vds,
        ref_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["GQ", "DP"],
        reducible_aggs={"AN"} if reduce_min_aggs else None,
        entry_agg_group_membership=entry_agg_group_membership,
        sex_karyotype_field="sex_karyotype",
    )

    # This expression aggregates the DP counter in reverse order of the cov_bins and
    # computes the cumulative sum over them. It needs to be in reverse order because we
    # want the sum over samples covered by > X.
    def _cov_stats(
        cov_stat: hl.expr.StructExpression,
    ) -> hl.expr.StructExpression:
        """Convert per-DP coverage_counter into cumulative ``over_X`` sample counts."""
        # The coverage was already floored to the max_coverage_bin, so no more
        # aggregation is needed for the max bin.
        count_expr = cov_stat.coverage_counter
        max_bin_expr = hl.int32(count_expr.get(max_cov_bin, 0))

        # For each of the other bins, coverage is summed between the boundaries.
        bin_expr = hl.range(hl.len(cov_bins) - 1, 0, step=-1)
        bin_expr = bin_expr.map(
            lambda i: hl.sum(
                hl.range(cov_bins[i - 1], cov_bins[i]).map(
                    lambda j: hl.int32(count_expr.get(j, 0))
                )
            )
        )
        bin_expr = hl.cumulative_sum(hl.array([max_bin_expr]).extend(bin_expr))
        # NOTE: Keeping these as sample counts rather than fractions to join with
        # gnomAD v4 genomes.
        bin_expr = {f"over_{x}": bin_expr[i] for i, x in enumerate(rev_cov_bins)}
        return cov_stat.annotate(**bin_expr).drop("coverage_counter")

    # Keep coverage stats from global adj grouping (index 0) only.
    ht = ht.annotate_globals(
        coverage_stats_meta_sample_count=ht.strata_sample_count[0],
    )
    cov_stats_expr = _cov_stats(ht.coverage_stats[0])

    ht = ht.transmute(**cov_stats_expr)

    if project == "aou":
        # `qual_hists` as returned by `compute_stats_per_ref_site` is an array of length 1 so we drop the array here.
        ht = ht.annotate(qual_hists=ht.qual_hists[0])
    return ht


def _rename_cov_annotations(
    ht: hl.Table,
    project: str,
    sample_count: int,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
) -> hl.Table:
    """
    Rename coverage annotations prior to merging Tables.

    Function transforms mean back into sum using `sample_count` argument.

    :param ht: Input HT.
    :param project: Project name.
    :param sample_count: Number of samples in HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :return: Renamed HT.
    """
    # Transform mean back into sum.
    ht = ht.transmute(sum=ht.mean * sample_count)

    # Rename annotations to include project.
    row_fields = list(ht.row_value)
    rename_dict = {f: f"{f}_{project}" for f in row_fields}
    ht = ht.rename(rename_dict)
    if project == "gnomad" or project == "gnomad_release":
        # Revert v4 genomes fraction over X bins to sample count over X bins.
        ht = ht.transmute(
            **{
                f"over_{x}_{project}": (
                    hl.float64(ht[f"over_{x}_{project}"]) * sample_count
                )
                for x in coverage_over_x_bins
            }
        )
    return ht


def _merge_coverage_fields(
    ht: hl.Table,
    project_1: str,
    project_2: str,
    sample_count: int,
    operation: str,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
) -> hl.expr.DictExpression:
    """
    Merge coverage fields from two Tables.

    .. note::
        - Function does not merge `median_approx` fields.

    :param ht: Input HT. Must have annotations from both projects.
    :param project_1: First project name.
    :param project_2: Second project name.
    :param sample_count: Total sample count.
    :param operation: Operation to perform on the coverage fields. Must be "sum" or "diff".
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :return: Merged fields.
    """
    if operation == "diff":
        # Make sure over_x fields are float64s.
        merged_fields = {
            "sum_gnomad": ht[f"sum_{project_1}"] - ht[f"sum_{project_2}"],
            "total_DP_gnomad": (
                ht[f"total_DP_{project_1}"] - ht[f"total_DP_{project_2}"]
            ),
        }
        merged_fields.update(
            {
                f"over_{x}_gnomad": (
                    ht[f"over_{x}_{project_1}"] - ht[f"over_{x}_{project_2}"]
                )
                for x in coverage_over_x_bins
            }
        )
    else:
        merged_fields = {
            "mean": (ht[f"sum_{project_1}"] + ht[f"sum_{project_2}"]) / sample_count,
        }
        merged_fields.update(
            {
                f"over_{x}": (
                    (ht[f"over_{x}_{project_1}"] + ht[f"over_{x}_{project_2}"])
                    / sample_count
                )
                for x in coverage_over_x_bins
            }
        )
    return merged_fields


def merge_gnomad_coverage_hts(
    gnomad_ht: hl.Table,
    gnomad_release_ht: hl.Table,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
    v4_count: int = 76215,
    consent_drop_count: int = 866,
) -> hl.Table:
    """
    Subtract consent drop samples from gnomAD v4 genomes release HT to create gnomAD v5 genomes coverage HT.

    :param gnomad_ht: gnomAD coverage HT (contains coverage for consent drop samples only).
    :param gnomad_release_ht: gnomAD v4 genomes coverage release HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :param v4_count: Number of release gnomAD v4 genome samples. Default is 76215.
    :param consent_drop_count: Number of consent drop gnomAD v4 genome samples. Default is 866.
    :return: gnomAD v5 genomes coverage HT.
    """
    logger.info(
        "Subtracting gnomAD v4 consent drop samples from gnomAD v4 genomes release HT..."
    )
    gnomad_ht = _rename_cov_annotations(
        gnomad_ht, "gnomad", consent_drop_count, coverage_over_x_bins
    )
    gnomad_release_ht = _rename_cov_annotations(
        gnomad_release_ht, "gnomad_release", v4_count, coverage_over_x_bins
    )
    gnomad_v5_count = v4_count - consent_drop_count
    logger.info("Total number of gnomAD v5 release genomes: %s", gnomad_v5_count)

    gnomad_ht = gnomad_ht.join(gnomad_release_ht, "right")
    merged_fields = _merge_coverage_fields(
        ht=gnomad_ht,
        project_1="gnomad_release",
        project_2="gnomad",
        sample_count=gnomad_v5_count,
        operation="diff",
    )
    gnomad_ht = gnomad_ht.transmute(**merged_fields)

    # Keep median_approx from v4 release.
    gnomad_ht = gnomad_ht.transmute(
        median_approx_gnomad=gnomad_ht.median_approx_gnomad_release
    )

    # Drop unnecessary globals and add back v5 gnomAD genomes count.
    gnomad_ht = gnomad_ht.select_globals()
    gnomad_ht = gnomad_ht.annotate_globals(
        coverage_stats_meta_sample_count=gnomad_v5_count,
    )
    return gnomad_ht


def join_aou_and_gnomad_coverage_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
    gnomad_v5_count: int = 76215 - 866,
) -> hl.Table:
    """
    Join AoU and gnomAD coverage HTs for release.

    :param aou_ht: AoU coverage HT.
    :param gnomad_ht: gnomAD v5 genomes coverage HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :param gnomad_v5_count: Number of release gnomAD v5 genome samples. Default is 76215 - 866.
    :return: Joined HT.
    """
    aou_count = hl.eval(aou_ht.coverage_stats_meta_sample_count)
    logger.info("Total number of AoU v8 release samples: %s", aou_count)

    logger.info("Merging AoU and gnomAD v5 coverage HTs...")
    aou_ht = _rename_cov_annotations(aou_ht, "aou", aou_count, coverage_over_x_bins)
    v5_count = aou_count + gnomad_v5_count
    logger.info("Total number of AoU + gnomAD v5 release genomes: %s", v5_count)
    ht = aou_ht.join(gnomad_ht, "left")
    merged_fields = _merge_coverage_fields(
        ht=ht,
        project_1="aou",
        project_2="gnomad",
        sample_count=v5_count,
        operation="sum",
    )
    ht = ht.transmute(**merged_fields)
    return ht.select_globals()


def _rename_fields(
    ht: hl.Table, field_name: str, project: str, rename_globals: bool
) -> hl.Table:
    """
    Rename fields by adding project name prior to merging Tables.

    Used for AN and qual hists merging but not coverage because
    coverage does not have globals and also requires extra transformations
    to transform v4 mean back into sum.

    :param ht: Input HT.
    :param project: Project name.
    :param rename_globals: Whether to rename globals.
    :return: Renamed HT.
    """
    if rename_globals:
        rename_globals = {
            f"strata_meta_{project}": ht.strata_meta,
            f"strata_sample_count_{project}": ht.strata_sample_count,
        }
        ht = ht.transmute_globals(**rename_globals)
    rename_dict = {f"{field_name}": f"{field_name}_{project}"}
    return ht.rename(rename_dict)


def _merge_an_fields(
    ht: hl.Table, project_1: str, project_2: str, operation: str
) -> Tuple[
    hl.expr.ArrayExpression,
    hl.expr.ArrayExpression,
    hl.expr.DictExpression,
]:
    """
    Merge AN fields from two projects.

    :param ht: Input HT. Must have annotations from both projects.
    :param project_1: First project name.
    :param project_2: Second project name.
    :param operation: Operation to perform on the AN fields. Must be "sum" or "diff".
    :return: Joined AN, strata meta, and count arrays.
    """
    joint_an, joint_strata_meta, count_arrays_dict = merge_array_expressions(
        arrays=[ht[f"AN_{project_1}"], ht[f"AN_{project_2}"]],
        meta=[
            ht.index_globals()[f"strata_meta_{project_1}"],
            ht.index_globals()[f"strata_meta_{project_2}"],
        ],
        count_arrays={
            "counts": [
                ht.index_globals()[f"strata_sample_count_{project_1}"],
                ht.index_globals()[f"strata_sample_count_{project_2}"],
            ],
        },
        operation=operation,
    )
    return joint_an, joint_strata_meta, count_arrays_dict


def merge_gnomad_an_hts(
    gnomad_ht: hl.Table,
    gnomad_release_ht: hl.Table,
) -> hl.Table:
    """
    Subtract consent drop samples from gnomAD v4 genomes release HT to create gnomAD v5 genomes AN HT.

    :param gnomad_ht: gnomAD AN HT (contains AN for consent drop samples only).
    :param gnomad_release_ht: gnomAD v4 genomes release AN HT.
    :return: gnomAD v5 genomes AN HT.
    """
    logger.info(
        "Subtracting gnomAD v4 consent drop samples from gnomAD v4 genomes release HT..."
    )
    gnomad_ht = _rename_fields(gnomad_ht, "AN", "gnomad", rename_globals=True)
    gnomad_release_ht = _rename_fields(
        gnomad_release_ht, "AN", "gnomad_release", rename_globals=True
    )

    gnomad_ht = gnomad_ht.join(gnomad_release_ht, "right")
    joint_an, joint_strata_meta, count_arrays_dict = _merge_an_fields(
        ht=gnomad_ht,
        project_1="gnomad_release",
        project_2="gnomad",
        operation="diff",
    )
    gnomad_ht = gnomad_ht.annotate(AN_gnomad=joint_an)
    gnomad_ht = gnomad_ht.annotate_globals(
        strata_meta_gnomad=joint_strata_meta,
        strata_sample_count_gnomad=count_arrays_dict["counts"],
    )
    gnomad_ht = gnomad_ht.drop(
        "strata_meta_gnomad_release",
        "strata_sample_count_gnomad_release",
        "coverage_stats_meta_sample_count",
    )
    gnomad_ht = gnomad_ht.select("AN_gnomad")
    return gnomad_ht


def join_aou_and_gnomad_an_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
) -> hl.Table:
    """
    Join AoU and gnomAD AN HTs for release.

    :param aou_ht: AoU AN HT.
    :param gnomad_ht: gnomAD v5 genomes AN HT.
    :return: Joined HT.
    """
    aou_ht = _rename_fields(aou_ht, "AN", "aou", rename_globals=True)

    logger.info("Merging AoU and gnomAD v5 AN HTs...")
    ht = aou_ht.join(gnomad_ht, "left")
    ht = ht.checkpoint(new_temp_file("aou_and_gnomad_join", "ht"))

    joint_an, joint_strata_meta, count_arrays_dict = _merge_an_fields(
        ht=ht,
        project_1="aou",
        project_2="gnomad",
        operation="sum",
    )
    ht = ht.annotate(AN=joint_an)
    ht = ht.annotate_globals(
        strata_meta=joint_strata_meta,
        strata_sample_count=count_arrays_dict["counts"],
    )
    return ht


def join_aou_and_gnomad_qual_hists_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
) -> hl.Table:
    """
    Join AoU and gnomAD qual hists HTs for release.

    .. note::
        We did not compute qual hists for the gnomAD v4 genomes release
        (https://github.com/broadinstitute/gnomad_qc/blob/e65bdbb5768113c0129199a875d845da245690e2/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1139).
        This means we will also not recompute hists on the gnomAD v4 genomes for v5,
        which also means we will not subtract values from the samples to drop for consent reasons.

    :param aou_ht: AoU qual hists HT.
    :param gnomad_ht: gnomAD qual hists HT.
    :return: Joined HT.
    """
    aou_ht = _rename_fields(aou_ht, "qual_hists", "aou", rename_globals=False)
    gnomad_ht = _rename_fields(gnomad_ht, "qual_hists", "gnomad", rename_globals=False)
    ht = aou_ht.join(gnomad_ht, "left")
    qual_hists = [
        "gq_hist_all",
        "dp_hist_all",
    ]
    # Only adj qual_hists are kept in v5 (raw_qual_hists were dropped at compute
    # time; see _get_hists). The reused gnomAD v4 hists still carry raw_qual_hists,
    # but we don't merge or emit it.
    hist_structs = {
        "qual_hists": qual_hists,
    }
    hists_expr = {
        hist_struct: hl.struct(
            **{
                h: merge_histograms(
                    [
                        ht.qual_hists_aou[hist_struct][h],
                        ht.qual_hists_gnomad[hist_struct][h],
                    ],
                    operation="sum",
                )
                for h in hists
            }
        )
        for hist_struct, hists in hist_structs.items()
    }
    ht = ht.annotate(**hists_expr)
    return ht


def _filter_to_locus_bounds(target_ht: hl.Table, source_ht: hl.Table) -> hl.Table:
    """
    Filter `target_ht` to per-contig locus bounds derived from `source_ht`.

    Used in test runs where `target_ht` (e.g., a public release HT) and
    `source_ht` (e.g., a freshly-computed test HT) do not share partition
    layouts, so their first-N-partition slices do not overlap. Filtering to
    the per-contig min..max locus range of `source_ht` guarantees test
    join/merge overlap at low cost (one small aggregation + partition-pruned
    interval filter).

    :param target_ht: Table to filter.
    :param source_ht: Table whose locus ranges define the filter intervals.
    :return: `target_ht` filtered to the per-contig locus ranges of
        `source_ht`.
    """
    bounds = source_ht.aggregate(
        hl.agg.group_by(
            source_ht.locus.contig,
            hl.tuple(
                [
                    hl.agg.min(source_ht.locus.position),
                    hl.agg.max(source_ht.locus.position),
                ]
            ),
        )
    )
    intervals = [
        hl.parse_locus_interval(f"{contig}:{lo}-{hi + 1}", reference_genome="GRCh38")
        for contig, (lo, hi) in bounds.items()
    ]
    return hl.filter_intervals(target_ht, intervals)


def _load_project_vds(
    project: str,
    environment: str,
    partition_range: Optional[List[int]] = None,
    sub_intervals: Optional[List[hl.utils.Interval]] = None,
    chrom: Optional[List[str]] = None,
    test: bool = False,
    test_sample_subset: bool = False,
) -> Tuple[hl.vds.VariantDataset, str]:
    """
    Load the per-project VDS with consistent test/subsample handling.

    Centralizes the project-conditional VDS load that's needed by both the
    chunk worker (``_run_coverage_chunk``) and the strict path (``main``):

    - filter_partitions vs read_intervals routing (both AoU and gnomAD):
      ``sub_intervals`` takes precedence; ``partition_range`` is the
      fallback. When ``sub_intervals`` is provided the VDS is read with
      ``read_intervals`` so it's co-partitioned with the vep_context read
      on the same intervals (a shuffle-free densify join).
    - AoU-specific DP synthesis from LAD + ``annotate_adj_no_dp`` (the AoU
      v8 VDS lacks DP, which ``compute_stats_per_ref_site`` requires).
    - ``test_sample_subset`` application: AoU only, reads ``meta`` from
      disk inside the helper so it's self-contained.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param partition_range: VDS partition indices (e.g. ``list(range(3))``)
        or ``None`` for the full VDS.
    :param sub_intervals: Locus sub-intervals for read-time partitioning
        (AoU only; ignored for gnomAD).
    :param chrom: Optional list of contigs to filter to.
    :param test: Whether this is a test run (gates ``test_sample_subset``).
    :param test_sample_subset: If True (AoU only, and ``test``), subsample
        ~0.1% of the AoU sample set for cheap tiny-cohort runs.
    :return: Tuple of ``(vds, sex_karyotype_field)``.
    """
    if project == "aou":
        sex_karyotype_field = "meta.sex_karyotype"
        vds = get_aou_vds(
            release_only=True,
            filter_partitions=None if sub_intervals else partition_range,
            read_intervals=sub_intervals,
            annotate_meta=True,
            chrom=chrom,
            environment=environment,
        )
        vmt = vds.variant_data
        vmt = annotate_adj_no_dp(vmt)
        vmt = vmt.annotate_entries(DP=hl.sum(vmt.LAD))
        vds = hl.vds.VariantDataset(vds.reference_data, vmt)

        if test and test_sample_subset:
            meta_ht = hl.read_table(
                meta(data_type="genomes", environment=environment).path
            )
            meta_ht = meta_ht.filter(
                (meta_ht.project_meta.project == project) & (meta_ht.release)
            ).select()
            meta_ht = meta_ht.sample(0.001, seed=42)
            vds = hl.vds.filter_samples(vds, meta_ht)
    else:
        sex_karyotype_field = "meta.sex_imputation.sex_karyotype"
        vds = get_gnomad_v5_genomes_vds(
            release_only=True,
            consent_drop_only=True,
            filter_partitions=None if sub_intervals else partition_range,
            read_intervals=sub_intervals,
            annotate_meta=True,
            chrom=chrom,
        )
    return vds, sex_karyotype_field


def _probe_vds(
    project: str,
    environment: str,
    partition_range: Optional[List[int]],
    chrom: Optional[List[str]],
) -> hl.vds.VariantDataset:
    """
    Cheap reference_data-bounds probe-load of the per-project VDS for a partition range.

    Loads via ``filter_partitions`` only (no sample filtering / DP synthesis); the
    caller derives a locus extent from ``reference_data`` (e.g. via
    ``repartition_for_join``). Used by both the chunk worker and the strict path's
    ``--test-n-partitions`` co-partitioning branch.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param partition_range: VDS partition indices to probe (e.g.
        ``list(range(2))``) or ``None`` for the full VDS.
    :param chrom: Optional list of contigs to filter to.
    :return: Probe VDS (reference/variant data only).
    """
    if project == "aou":
        return get_aou_vds(
            filter_partitions=partition_range,
            chrom=chrom,
            environment=environment,
            remove_hard_filtered_samples=False,
            log_sample_counts=False,
        )
    return get_gnomad_v5_genomes_vds(
        filter_partitions=partition_range,
        chrom=chrom,
    )


def _vep_context_sites_path(environment: str) -> str:
    """
    Return the path to the preprocessed vep_context sites HT (30-day storage).

    The content (deduped, telomere/centromere-stripped, locus-keyed vep_context)
    is the same for test and prod, so a single per-environment path is used.

    :param environment: Compute environment.
    :return: GCS path to the preprocessed sites HT.
    """
    return f"{qc_temp_prefix(environment=environment, days=30)}vep_context_sites.ht"


def _build_vep_context_sites_ht() -> hl.Table:
    """
    Build the locus-keyed, deduped, telomere/centromere-stripped vep_context sites HT.

    This is the per-ref-site set the coverage/AN compute aggregates at. Previously
    this dedup + strip ran inside every chunk (``.distinct()`` shuffle + a driver
    ``collect()`` of the telomere intervals, ~48k times); it is now done once via
    ``--write-vep-context-sites`` and read co-partitioned per chunk.

    :return: Locus-keyed sites HT.
    """
    ref_ht = vep_context.versions["105"].ht().key_by("locus").select().distinct()
    ref_ht = hl.filter_intervals(
        ref_ht,
        telomeres_and_centromeres.ht().interval.collect(),
        keep=False,
    )
    return ref_ht


def _build_chunk_ref_ht(
    vds_filtered: hl.vds.VariantDataset,
    partition_count: int,
    chrom: Optional[List[str]],
    sites_path: str,
    sub_intervals: Optional[List[hl.utils.Interval]] = None,
) -> hl.Table:
    """
    Build the per-chunk ``ref_ht`` from the preprocessed vep_context sites HT.

    Reads ``sites_path`` (already locus-keyed, deduped, and telomere/centromere-
    stripped by ``--write-vep-context-sites``). When ``sub_intervals`` is provided,
    reads with ``_intervals=sub_intervals`` so the ref_ht is co-partitioned with the
    VDS read on the same intervals (shuffle-free densify join). Otherwise (the
    strict path / single-interval chunks) reads the whole sites HT and filters to
    ``chrom`` (when set) or the VDS chunk's per-contig locus extent.

    :param vds_filtered: VDS already filtered to the chunk's partitions
        (the locus extent is derived from its ``reference_data``).
    :param partition_count: When > 0 and ``sub_intervals``/``chrom`` are None,
        triggers the locus-extent filter against the VDS chunk. 0 means "whole".
    :param chrom: Optional list of contigs to filter to (takes precedence over
        the ``partition_count`` locus-extent filter).
    :param sites_path: Path to the preprocessed vep_context sites HT.
    :param sub_intervals: Optional read-time intervals; take precedence over both
        ``chrom`` and the locus-extent filter.
    :return: Reference HT keyed by locus.
    """
    if sub_intervals is not None:
        logger.info(
            "Reading vep_context sites with %d sub-intervals (one partition each).",
            len(sub_intervals),
        )
        return hl.read_table(sites_path, _intervals=sub_intervals)

    ref_ht = hl.read_table(sites_path)
    if chrom:
        ref_ht = hl.filter_intervals(
            ref_ht, [hl.parse_locus_interval(c) for c in chrom]
        )
    elif partition_count > 0:
        chunk_intervals = _derive_chunk_locus_intervals(vds_filtered)
        ref_ht = hl.filter_intervals(ref_ht, chunk_intervals)
    return ref_ht


def _run_coverage_chunk(args: argparse.Namespace) -> None:
    """
    Compute one chunk of the coverage/AN HT and write it to ``args.chunk_output``.

    Invoked by the per-chunk Hail Batch relay container (see
    ``_orchestrate_coverage_batch``). Hail must already be initialized
    (``main`` does this via ``_init_hail`` before dispatching here).

    Steps:
      1. Probe-load the VDS via
         ``filter_partitions=range(chunk_start, chunk_stop)`` and derive
         ``--read-subintervals-per-chunk`` balanced locus sub-intervals over the
         chunk via ``repartition_for_join`` (``_calculate_new_partitions``, not an
         ``aggregate_rows`` scan).
      2. Re-read the VDS via ``hl.vds.read_vds(intervals=sub_intervals)`` — one
         VDS partition per sub-interval (no shuffle).
      3. Read ``vep_context`` with the same ``_intervals=sub_intervals`` —
         co-partitioned with the VDS, so the densify join is shuffle-free.
      4. Run ``compute_all_release_stats_per_ref_site``.
      5. Write to ``args.chunk_output``.

    When ``--read-subintervals-per-chunk=1`` (legacy behavior), skips the
    probe/re-read and reads the VDS once via ``filter_partitions``, then
    builds ``ref_ht`` via the locus-extent filter path.

    :param args: Parsed CLI args; reads ``project_name``, ``environment``,
        ``chunk_start``, ``chunk_stop``, ``test``, ``chrom``,
        ``read_subintervals_per_chunk``, ``chunk_output``,
        ``cov_and_an_output_suffix``, ``reduce_min_aggs``,
        ``test_sample_subset``. Hail must already be initialized.
    :return: None.
    """
    project = args.project_name
    environment = args.environment
    # __main__ validation guarantees --chunk-stop > --chunk-start, so
    # partition_range is always a non-empty list.
    start, stop = args.chunk_start, args.chunk_stop
    partition_range = list(range(start, stop))
    n = stop - start

    test = args.test
    chrom = args.chrom
    n_sub = max(args.read_subintervals_per_chunk, 1)

    if args.chunk_output is None:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
        )
        args.chunk_output = _chunk_path(cov_and_an_ht_path, start)
        logger.info("Auto-derived --chunk-output: %s", args.chunk_output)

    group_membership_ht_path = _resolve_group_membership_ht_path(
        project,
        environment,
        test=test,
        reduce_min_aggs=args.reduce_min_aggs,
    )

    sub_intervals: Optional[List[hl.utils.Interval]] = None
    if n_sub > 1 and chrom:
        logger.warning(
            "--chrom is set, so sub-interval co-partitioning (read-subintervals-"
            "per-chunk=%d) is disabled; this chunk uses the slower filter_partitions"
            " + locus-extent path (test/prod path divergence).",
            n_sub,
        )
    if n_sub > 1 and not chrom:
        # Probe: cheap reference_data-bounds load via filter_partitions.
        vds_probe = _probe_vds(project, environment, partition_range, chrom)
        sub_intervals = repartition_for_join(
            vds_probe.reference_data.rows(), n_partitions=n_sub, locus_intervals=True
        )
        logger.info(
            "Derived %d balanced sub-intervals from chunk [%d, %d) for read-time"
            " partitioning.",
            len(sub_intervals),
            start,
            stop,
        )

    vds, sex_karyotype_field = _load_project_vds(
        project=project,
        environment=environment,
        partition_range=partition_range,
        sub_intervals=sub_intervals,
        chrom=chrom,
        test=test,
        test_sample_subset=args.test_sample_subset,
    )

    ref_ht = _build_chunk_ref_ht(
        vds_filtered=vds,
        partition_count=n,
        chrom=chrom,
        sites_path=_vep_context_sites_path(environment),
        sub_intervals=sub_intervals,
    )

    cov_and_an_ht = compute_all_release_stats_per_ref_site(
        vds,
        ref_ht,
        sex_karyotype_field=sex_karyotype_field,
        project=project,
        group_membership_ht=hl.read_table(group_membership_ht_path),
        reduce_min_aggs=args.reduce_min_aggs,
    )
    cov_and_an_ht.write(args.chunk_output, overwrite=True)
    logger.info("Wrote chunk [%d, %d) to %s", start, stop, args.chunk_output)


def _run_coverage_merge(
    input_paths: List[str],
    output_path: str,
    coalesce_to: Optional[int] = None,
) -> None:
    """
    Union per-chunk coverage HTs and write the merged HT to ``output_path``.

    Globals (``strata_meta``, etc.) are identical across all chunks (computed
    from the same group_membership HT), so the union inherits them from the
    first input. Used for both group-level and final merges in the
    ``--merge-cov-chunks`` pipeline.

    :param input_paths: HT paths to union.
    :param output_path: Destination HT path.
    :param coalesce_to: If set, ``naive_coalesce`` to this many partitions
        before writing. For group-level merges, set to roughly the per-chunk
        partition count to avoid the partition-count blowup
        (sum-of-input-partitions). For the final merge, set to the desired
        final partition count (e.g. ``--n-partitions``).
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


def _build_setup_command(
    commit: str,
    gcp_billing_project: str,
    methods_branch: str = "main",
) -> str:
    """
    Build shell commands to download gnomad_qc and gnomad_methods.

    Both repos are actively developed, so chunk containers pull them at
    runtime rather than relying on what's baked into the Docker image. The
    image provides ``hail`` and system dependencies. Mirrors the equivalent
    helper in ``generate_frequency.py`` from earlier this branch.

    Also writes the Hail ``config.ini`` (Batch billing project,
    ``remote_tmpdir``, and GCS requester-pays project) so the relay's
    ``hl.init(backend="batch")`` finds them; see the inline comment below
    for why both the XDG and legacy paths are written.

    Also patches ``/gsa-key/key.json`` with a ``quota_project_id`` field
    so requester-pays reads (e.g., the AoU VDS) succeed from the QoB
    driver pod the relay spawns. Without this field Hail's QoB doesn't
    fully propagate ``gcs_requester_pays_configuration`` through to the
    driver pod's Java GCS client, and the read 400s with no fallback
    (works fine from a laptop because gcloud config supplies the quota
    project there).

    :param commit: Git commit hash to pin gnomad_qc to.
    :param gcp_billing_project: GCP project for requester-pays reads;
        patched into the GSA key as ``quota_project_id``.
    :param methods_branch: Branch/commit of gnomad_methods to pull.
    :return: Shell command string (terminated with newline).
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = (
        "https://github.com/broadinstitute/gnomad_methods/archive/"
        f"{methods_branch}.tar.gz"
    )
    methods_dir_suffix = methods_branch.replace("/", "-")
    # Hail config so `hl.init(backend="batch")` finds the Batch billing
    # project and remote_tmpdir, plus the GCS requester-pays project.
    # The canonical path is the XDG-style `~/.config/hail/config.ini`
    # (what `hailtop.config.get_user_config_path()` returns); also write
    # `~/.hail/config.ini` as a legacy fallback for older Hail versions.
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
        # TODO: Remove this GSA-key patch once Hail's
        # gcs_requester_pays_configuration propagation to the QoB driver
        # pod is fixed (likely 0.2.139+). Patches quota_project_id into
        # /gsa-key/key.json so the QoB driver pod's Java GCS client has
        # a billing-project fallback for requester-pays reads.
        f"python3 -c \"import json, os; p='/gsa-key/key.json';"
        f" d=json.load(open(p)); d['quota_project_id']='{gcp_billing_project}';"
        f" json.dump(d, open(p+'.new','w')); os.replace(p+'.new', p)\"\n"
        # TODO: Remove this Hail pin once 0.2.139+ fixes the
        # requester-pays propagation regression that 0.2.138 introduced
        # for `load_references_from_dataset` (AoU VDS metadata reads
        # 400'd with "no user project" until we pinned back to 0.2.137).
        # The relay's Hail Python version determines the JAR the QoB
        # driver downloads, so pinning here pins the entire pipeline.
        "/opt/venv/bin/pip install --quiet --upgrade --force-reinstall"
        " --no-deps hail==0.2.137\n"
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_dir_suffix} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def _build_relay_common_flags(args: argparse.Namespace, *, chunk: bool) -> str:
    """
    Build the CLI flag string shared by per-chunk / per-merge relay invocations.

    Translates the orchestrator's ``args.Namespace`` into a CLI string appended to
    each relay's ``--run-chunk`` / ``--run-merge`` command. The shared base
    (project / environment / billing / tmp-dir / n-partitions / experimental /
    app-name / output-suffix / QoB driver+worker sizing) is common to both. With
    ``chunk=True`` the read/compute flags a chunk relay also needs
    (``--read-subintervals-per-chunk``, ``--reduce-min-aggs``,
    ``--test-sample-subset``, ``--chrom``, ``--test``) are appended; the merge
    relay (which only unions HTs) omits them. Per-job flags (``--chunk-*`` /
    ``--merge-*``) are added by the submit helpers. Flag order is irrelevant to
    argparse.

    :param args: Parsed CLI args.
    :param chunk: If True, include the chunk-only read/compute flags; if False,
        build the (smaller) merge relay flag set.
    :return: Space-joined ``--flag value`` string.
    """
    flags = [
        f"--project-name {args.project_name}",
        "--environment batch",
        f"--gcp-billing-project {args.gcp_billing_project}",
        f"--tmp-dir-days {args.tmp_dir_days}",
    ]
    if args.n_partitions is not None:
        flags.append(f"--n-partitions {args.n_partitions}")
    if args.experimental:
        # With --experimental the inner QoB driver attaches to this outer batch
        # (HAIL_BATCH_ID); without it, each relay's QoB creates its own Hail Batch.
        flags.append("--experimental")
    if args.app_name:
        flags.append(f"--app-name {args.app_name}")
    if args.cov_and_an_output_suffix:
        flags.append(f"--cov-and-an-output-suffix {args.cov_and_an_output_suffix}")
    # QoB driver/worker sizing; merge workers reuse the chunk sizing (same
    # QoB-from-container shape).
    if args.chunk_driver_cores is not None:
        flags.append(f"--driver-cores {args.chunk_driver_cores}")
    if args.chunk_driver_memory:
        flags.append(f"--driver-memory {args.chunk_driver_memory}")
    if args.chunk_worker_cores is not None:
        flags.append(f"--worker-cores {args.chunk_worker_cores}")
    if args.chunk_worker_memory:
        flags.append(f"--worker-memory {args.chunk_worker_memory}")
    if chunk:
        flags.append(
            f"--read-subintervals-per-chunk {args.read_subintervals_per_chunk}"
        )
        if args.reduce_min_aggs:
            flags.append("--reduce-min-aggs")
        if args.test_sample_subset:
            flags.append("--test-sample-subset")
        if args.chrom:
            flags.append(f"--chrom {' '.join(args.chrom)}")
        if args.test:
            flags.append("--test")
    return " ".join(flags)


class _RelayJobSpec(NamedTuple):
    """One relay job's per-job config for :func:`_submit_relay_batch`."""

    name: str
    cpu: float
    memory: str
    storage: str
    command: str


def _submit_relay_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    batch_name: str,
    job_specs: List[_RelayJobSpec],
    log_label: str,
) -> None:
    """
    Build and submit one Hail Batch of relay jobs sharing the same config.

    Shared skeleton for ``_submit_chunk_batch`` and ``_submit_merge_batch``: each
    relay job is a non-spot coordinator container (preemption mid-wait would orphan
    its inner QoB job) pinned to ``BATCH_REGIONS``, with ``--chunk-attempts``
    retries. Parallelism comes from Hail Batch's own scheduler running the N jobs
    concurrently. No-ops (skips ``batch.run()``) when ``job_specs`` is empty.

    :param args: Parsed CLI args (reads ``batch_image``, ``chunk_attempts``,
        ``batch_dry_run``).
    :param backend_kwargs: kwargs for the ``hb.ServiceBackend(...)`` constructor.
    :param batch_name: Hail Batch name.
    :param job_specs: Per-job config (name + QoB sizing + command).
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
            # Relay is a coordinator waiting on its inner QoB job; preemption
            # mid-wait orphans that inner job, so relays are non-spot.
            j.spot(False)
            j.n_max_attempts(args.chunk_attempts)
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


def _submit_chunk_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    chunk_indices: List[int],
    cov_and_an_ht_path: str,
    setup_cmd: str,
    common_flags_str: str,
    script: str,
    wave_label: Optional[str] = None,
) -> None:
    """
    Build and submit one Hail Batch containing all pending chunk jobs.

    Each chunk job is a relay container that runs ``--run-chunk`` (which
    spawns its own QoB driver). Parallelism comes from Hail Batch's
    own job scheduler — N parallel jobs in one batch run concurrently
    against Hail Batch's worker pool.

    Existence checks happen in the orchestrator's main thread before this
    function is called (see ``_orchestrate_coverage_batch``);
    ``chunk_indices`` is the already-filtered pending set.

    :param args: Parsed CLI args (reads ``project_name``,
        ``total_partitions``, ``partitions_per_chunk``,
        ``read_subintervals_per_chunk``, ``cov_and_an_output_suffix``,
        ``batch_image``, ``chunk_cpu/memory/storage``,
        ``chunk_attempts``, ``batch_dry_run``).
    :param backend_kwargs: kwargs for the per-call
        ``hb.ServiceBackend(...)`` constructor (``billing_project``,
        ``remote_tmpdir``).
    :param chunk_indices: Pending chunk indices to submit.
    :param cov_and_an_ht_path: Resolved canonical output HT path; each
        chunk writes to a sibling ``_chunks/<idx>.chunk.ht``.
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags from
        ``_build_relay_common_flags(args, chunk=True)``.
    :param script: Path to the script inside the relay container.
    :param wave_label: Optional suffix appended to the Hail Batch name to
        distinguish per-wave batches (e.g. ``"w03of49"``). None for a
        single (non-waved) batch.
    :return: None.
    """
    project = args.project_name
    total = args.total_partitions
    pp = args.partitions_per_chunk
    batch_name = (
        f"v5_cov_{project}_{total}p_{pp}ppc_sub{args.read_subintervals_per_chunk}"
    )
    if args.cov_and_an_output_suffix:
        batch_name += f"_{args.cov_and_an_output_suffix}"
    if wave_label:
        batch_name += f"_{wave_label}"

    job_specs = []
    for idx in chunk_indices:
        path = _chunk_path(cov_and_an_ht_path, idx)
        start = idx * pp
        stop = min(start + pp, total)
        command = (
            f"{setup_cmd}{script} --run-chunk"
            f" --chunk-start {start} --chunk-stop {stop}"
            f" --chunk-output {path}"
            f" {common_flags_str}"
        )
        job_specs.append(
            _RelayJobSpec(
                name=f"cov_chunk_{idx:06d}_{start}_{stop}",
                cpu=args.chunk_cpu,
                memory=args.chunk_memory,
                storage=args.chunk_storage,
                command=command,
            )
        )
    _submit_relay_batch(args, backend_kwargs, batch_name, job_specs, "chunk")


def _orchestrate_coverage_batch(
    args: argparse.Namespace, cov_and_an_ht_path: str
) -> None:
    """
    Fan per-partition coverage/AN compute out as relay chunk jobs in Hail Batch.

    Pending chunks are split into sequential **waves** of
    ``--wave-size`` chunks. Each wave is its own Hail Batch (one relay
    job per chunk; each relay spawns its own QoB driver via
    ``hl.init(backend="batch")``) and is run to completion before the
    next wave is submitted. Waves are sequential — ``batch.run()`` blocks
    until the wave finishes — which (a) bounds the number of
    concurrently-running chunk relays and their nested QoB drivers to
    ``--wave-size`` so a ~tens-of-thousands-chunk prod run doesn't
    overwhelm the shared Hail Batch service, and (b) is the
    only safe option: Hail Batch's progress display is a process-global
    rich.Live, so *concurrent* ``batch.run()`` calls crash with ``Only
    one live display may be active at once`` (sequential calls are
    fine). Within a wave, Hail Batch's own scheduler runs the relay jobs
    in parallel up to worker-pool capacity. Set ``--wave-size <= 0`` (or
    >= the pending count) for a single batch (legacy behavior).

    Note on the binding constraint: these relays and their nested QoB
    drivers/workers run on the **Hail Batch service's** autoscaled worker
    pool (provisioned in the Hail team's GCP project), not on
    gnomad-production Compute/Dataproc. So ``--wave-size`` is bounded by
    the shared service's capacity and by any per-billing-project
    concurrent-core cap the Hail team configures for ``gnomad-production``
    — *not* by gnomad-production's own GCP CPU / in-use-IP quota. Confirm
    that cap (and acceptable concurrency on the shared deployment) with
    the Hail team before the prod run.

    By default each chunk's inner QoB creates its own separate Hail Batch;
    pass ``--experimental`` to attach the inner QoB to this outer batch
    via ``HAIL_BATCH_ID``.

    Idempotency: chunks whose ``_SUCCESS`` already exists are skipped, so
    a failed/partial wave is resumed by simply rerunning this step (only
    not-yet-complete chunks are re-submitted). ``batch.run()`` does not
    raise on per-job failure, so a wave with some failed chunks does not
    abort the remaining waves; rerun to retry the failures.

    The merge step is intentionally separate; run ``--merge-cov-chunks``
    after this finishes.

    :param args: Parsed CLI args (reads ``project_name``,
        ``total_partitions``, ``partitions_per_chunk``, ``wave_size``,
        ``overwrite``, ``methods_branch``, ``gcp_billing_project``,
        ``batch_billing_project``, ``batch_remote_tmpdir``, plus
        everything ``_build_relay_common_flags`` and
        ``_submit_chunk_batch`` read).
    :param cov_and_an_ht_path: Canonical output cov_and_an HT path
        (resolved by ``main``). Each chunk writes to
        ``<cov_and_an_path>_chunks/<idx>.chunk.ht``; the final HT at
        ``cov_and_an_ht_path`` is produced by ``--merge-cov-chunks``.
    :return: None.
    """
    project = args.project_name
    total = args.total_partitions
    pp = args.partitions_per_chunk
    if total <= 0 or pp <= 0:
        raise ValueError(
            "--use-batch-fanout requires --total-partitions and"
            " --partitions-per-chunk to be > 0."
        )

    # Fail fast (before submitting thousands of jobs) if the preprocessed
    # vep_context sites HT the chunks read isn't there yet.
    sites_path = _vep_context_sites_path(args.environment)
    if not _file_exists_for_env(sites_path, "batch"):
        raise FileNotFoundError(
            f"vep_context sites HT not found at {sites_path}. Run"
            " --write-vep-context-sites first (it is a prerequisite for the"
            " chunk compute)."
        )

    n_chunks = (total + pp - 1) // pp

    # Pre-compute pending chunk indices so the summary log is accurate
    # and so we can short-circuit when nothing is pending.
    if args.overwrite:
        pending_indices = list(range(n_chunks))
    else:
        pending_indices = []
        for idx in range(n_chunks):
            path = _chunk_path(cov_and_an_ht_path, idx)
            if _file_exists_for_env(path, "batch"):
                logger.info("Skipping already-complete chunk %d at %s", idx, path)
            else:
                pending_indices.append(idx)
    logger.info(
        "Coverage fan-out: %d chunks total, %d pending, %d skipped"
        " (overwrite=%s, project=%s)",
        n_chunks,
        len(pending_indices),
        n_chunks - len(pending_indices),
        args.overwrite,
        project,
    )

    if not pending_indices:
        logger.info("All chunks already complete; nothing to submit.")
        return

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    setup_cmd = _build_setup_command(
        commit,
        gcp_billing_project=args.gcp_billing_project,
        methods_branch=args.methods_branch,
    )

    backend_kwargs = {"billing_project": args.batch_billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir

    common_flags_str = _build_relay_common_flags(args, chunk=True)
    script = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/compute_coverage.py"

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
        "Dispatching %d pending chunks in %d sequential wave(s) of up to"
        " %d chunks each.",
        len(pending_indices),
        n_waves,
        wave_size if wave_size > 0 else len(pending_indices),
    )

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
        _submit_chunk_batch(
            args=args,
            backend_kwargs=backend_kwargs,
            chunk_indices=wave_indices,
            cov_and_an_ht_path=cov_and_an_ht_path,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
            wave_label=wave_label,
        )
        # batch.run() does not raise on per-job failure, so re-probe this wave's
        # chunk outputs and surface any that did not land (rather than only
        # discovering them at --merge-cov-chunks time).
        failed = [
            idx
            for idx in wave_indices
            if not _file_exists_for_env(_chunk_path(cov_and_an_ht_path, idx), "batch")
        ]
        if failed:
            logger.warning(
                "Wave %d/%d complete but %d/%d chunk(s) MISSING after run"
                " (rerun --use-batch-fanout to retry); missing indices: %s%s",
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


def _submit_merge_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    group_indices: List[int],
    groups: List[List[str]],
    group_output_paths: List[str],
    setup_cmd: str,
    common_flags_str: str,
    script: str,
    level: int,
) -> None:
    """
    Build and submit one Hail Batch containing all pending group-merge jobs.

    Each job runs ``--run-merge`` over its assigned inputs and writes to
    ``<cov_and_an_path>_merge_groups/L<level>_<group_idx>.ht``.
    Parallelism is Hail Batch's own job scheduler running the N jobs
    concurrently against the worker pool.

    Per-group coalesce target is the number of inputs in that group
    (one output partition per chunk-equivalent of input work), so the
    last (possibly smaller) group still produces a valid coalesce.

    :param args: Parsed CLI args (reads ``project_name``,
        ``cov_and_an_output_suffix``, ``batch_image``,
        ``merge_cpu``, ``merge_memory``, ``merge_storage``,
        ``chunk_attempts``, ``batch_dry_run``).
    :param backend_kwargs: kwargs for the per-call
        ``hb.ServiceBackend(...)`` constructor.
    :param group_indices: Pending group indices to submit.
    :param groups: For each group index, the list of input HT paths to
        union.
    :param group_output_paths: For each group index, the output HT path
        this level's merge writes to.
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags from
        ``_build_relay_common_flags(args, chunk=False)``.
    :param script: Path to the script inside the relay container.
    :param level: Merge-tree level (1-indexed); used in the batch and
        job names.
    :return: None.
    """
    project = args.project_name
    batch_name = f"v5_cov_merge_L{level:02d}_{project}"
    if args.cov_and_an_output_suffix:
        batch_name += f"_{args.cov_and_an_output_suffix}"

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
                name=f"cov_merge_L{level:02d}_{group_idx:06d}",
                cpu=args.merge_cpu,
                memory=args.merge_memory,
                storage=args.merge_storage,
                command=command,
            )
        )
    _submit_relay_batch(args, backend_kwargs, batch_name, job_specs, "merge")


def _orchestrate_coverage_merge(
    args: argparse.Namespace, cov_and_an_ht_path: str
) -> None:
    """
    Two-level tree-reduce merge of per-chunk HTs into ``cov_and_an_ht_path``.

    Counterpart to ``_orchestrate_coverage_batch``: this orchestrator runs
    *after* ``--use-batch-fanout`` has produced all per-chunk HTs. It
    submits Hail Batch jobs (QoB-from-container per job) and exits without
    initializing Hail in this process.

    Discovery: chunk paths are reconstructed from ``--total-partitions`` /
    ``--partitions-per-chunk``. Every expected chunk must have a
    ``_SUCCESS`` marker; missing chunks fail loudly so the user re-runs
    the chunk orchestrator before merging.

    Recursive tree: starting from N chunks, each level groups its
    inputs into windows of ``--merge-group-size`` and emits one
    ``--run-merge`` job per group. Inputs of level 1 are the chunk HTs;
    inputs of level k>1 are the outputs of level k-1. Iteration stops
    when at most one group remains; that group becomes the final-merge
    job that writes to ``cov_and_an_ht_path``.

    Per-job coalesce target is the number of inputs in the group (one
    partition per input). The final merge additionally accepts
    ``--n-partitions`` to coalesce the final HT to a specific count;
    unset = keep natural sum-of-input partition count.

    Each level is one Hail Batch containing N parallel jobs (Hail
    Batch's scheduler handles intra-level parallelism). Levels run
    sequentially; the orchestrator blocks on each level via
    ``batch.run(wait=True)``.

    Idempotency: at each level, group HTs whose ``_SUCCESS`` exists are
    skipped; the final HT is skipped if it exists (unless ``--overwrite``).

    :param args: Parsed CLI args (reads ``project_name``,
        ``total_partitions``, ``partitions_per_chunk``,
        ``merge_group_size``, ``overwrite``, ``n_partitions``,
        ``methods_branch``, ``gcp_billing_project``,
        ``batch_billing_project``, ``batch_remote_tmpdir``,
        ``batch_image``, ``merge_cpu``, ``merge_memory``,
        ``final_merge_storage``, ``chunk_attempts``,
        ``cov_and_an_output_suffix``, ``batch_dry_run``).
    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :return: None.
    """
    project = args.project_name
    total = args.total_partitions
    pp = args.partitions_per_chunk
    if total <= 0 or pp <= 0:
        raise ValueError(
            "--merge-cov-chunks requires --total-partitions and"
            " --partitions-per-chunk > 0."
        )

    n_chunks = (total + pp - 1) // pp
    chunk_paths = [_chunk_path(cov_and_an_ht_path, i) for i in range(n_chunks)]

    logger.info("Verifying %d expected chunk HTs exist...", n_chunks)
    missing = [p for p in chunk_paths if not _file_exists_for_env(p, "batch")]
    if missing:
        raise FileNotFoundError(
            f"--merge-cov-chunks: {len(missing)} of {n_chunks} chunks"
            f" missing (first few: {missing[:5]}). Run --use-batch-fanout"
            " to (re)compute missing chunks first."
        )
    logger.info("All %d chunks present.", n_chunks)

    gs = args.merge_group_size

    # Precompute the level shape so we can log the full plan upfront.
    shape = [n_chunks]
    while shape[-1] > gs:
        shape.append((shape[-1] + gs - 1) // gs)
    # shape[-1] is the input count for the final merge (≤ gs). The
    # final merge writes 1 HT, which is cov_and_an_ht_path.
    logger.info(
        "Merge tree (group_size=%d): %s -> 1 final HT (%d intermediate level(s))",
        gs,
        " -> ".join(str(n) for n in shape),
        len(shape) - 1,
    )

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    setup_cmd = _build_setup_command(
        commit,
        gcp_billing_project=args.gcp_billing_project,
        methods_branch=args.methods_branch,
    )

    backend_kwargs = {"billing_project": args.batch_billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir

    script = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/compute_coverage.py"
    common_flags_str = _build_relay_common_flags(args, chunk=False)

    # Intermediate levels: each emits a level-tagged group HT per output.
    # Iterates while #inputs > gs; stops when one final merge can union
    # everything remaining.
    inputs = chunk_paths
    level = 1
    while len(inputs) > gs:
        n_in = len(inputs)
        n_out = (n_in + gs - 1) // gs
        groups = [inputs[i : i + gs] for i in range(0, n_in, gs)]
        out_paths = [
            _group_path(cov_and_an_ht_path, level, idx) for idx in range(n_out)
        ]

        if args.overwrite:
            pending = list(range(n_out))
        else:
            pending = []
            for idx in range(n_out):
                if _file_exists_for_env(out_paths[idx], "batch"):
                    logger.info(
                        "Skipping already-complete L%d group %d at %s",
                        level,
                        idx,
                        out_paths[idx],
                    )
                else:
                    pending.append(idx)
        logger.info(
            "Level %d dispatch: %d groups total, %d pending, %d skipped" " (%d -> %d)",
            level,
            n_out,
            len(pending),
            n_out - len(pending),
            n_in,
            n_out,
        )
        # _submit_merge_batch no-ops when there are no pending merges.
        _submit_merge_batch(
            args=args,
            backend_kwargs=backend_kwargs,
            group_indices=pending,
            groups=groups,
            group_output_paths=out_paths,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
            level=level,
        )
        inputs = out_paths
        level += 1

    # Final merge: a single job that unions the remaining (<= gs) inputs
    # and writes the canonical output.
    if not args.overwrite and _file_exists_for_env(cov_and_an_ht_path, "batch"):
        logger.info(
            "Final merge HT exists at %s; skipping (pass --overwrite to rewrite).",
            cov_and_an_ht_path,
        )
        return

    final_batch_name = f"v5_cov_merge_final_{project}"
    if args.cov_and_an_output_suffix:
        final_batch_name += f"_{args.cov_and_an_output_suffix}"
    coalesce_flag = (
        f" --merge-coalesce-to {args.n_partitions}"
        if args.n_partitions is not None
        else ""
    )
    logger.info("Final merge: %d inputs -> %s", len(inputs), cov_and_an_ht_path)
    final_spec = _RelayJobSpec(
        name="cov_merge_final",
        cpu=args.merge_cpu,
        memory=args.merge_memory,
        storage=args.final_merge_storage,
        command=(
            f"{setup_cmd}{script} --run-merge"
            f" --merge-output {cov_and_an_ht_path}"
            f"{coalesce_flag}"
            f" --merge-inputs {' '.join(inputs)}"
            f" {common_flags_str}"
        ),
    )
    _submit_relay_batch(
        args, backend_kwargs, final_batch_name, [final_spec], "final-merge"
    )


def main(args):
    """
    Compute all sites coverage, AN, and quality histograms for v5 genomes.

    Flat dispatcher over three mutually-exclusive process ROLES (see the
    "Execution roles" section in the module docstring for the full
    picture):

      - ROLE 1 ORCHESTRATOR (``--use-batch-fanout`` /
        ``--merge-cov-chunks``): submit a Hail Batch, NO Hail init,
        return.
      - ROLE 2 WORKER (``--run-chunk`` / ``--run-merge``): init Hail,
        run one chunk/merge unit, return — *before* the try/finally.
      - ROLE 3 IN-PROCESS PIPELINE (no role flag): init Hail, run the
        try/finally step chain (SETUP -> COMPUTE -> ASSEMBLE).
    """
    project = args.project_name
    environment = args.environment
    test = args.test or args.test_n_partitions is not None
    chrom = args.chrom
    overwrite = args.overwrite

    # ===================================================================
    # ROLE 1: ORCHESTRATOR — submit a Hail Batch of relay jobs and return
    # WITHOUT initializing Hail (each per-chunk relay spawns its own QoB
    # driver). Returns here; never reaches _init_hail or the try/finally.
    #   --use-batch-fanout : scatter the per-chunk compute.
    #   --merge-cov-chunks : tree-reduce the chunk HTs (run after fanout).
    # ===================================================================
    if args.use_batch_fanout:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
        )
        _orchestrate_coverage_batch(args, cov_and_an_ht_path)
        return

    if args.merge_cov_chunks:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
        )
        _orchestrate_coverage_merge(args, cov_and_an_ht_path)
        return

    log_name = _log_name_for_run(
        run_chunk=args.run_chunk,
        run_merge=args.run_merge,
        chunk_start=args.chunk_start,
        chunk_stop=args.chunk_stop,
        merge_output=args.merge_output,
    )

    # QoB / dataproc init — reached by ROLE 2 (worker) and ROLE 3
    # (in-process pipeline) only; ROLE 1 returned above. ``--experimental``
    # attaches the QoB driver to an existing Hail Batch (batch_id
    # auto-resolved from HAIL_BATCH_ID inside _init_hail); otherwise each
    # run creates its own Hail Batch.
    _init_hail(
        log_name,
        environment,
        billing_project=getattr(args, "gcp_billing_project", None),
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{qc_temp_prefix(environment=environment, days=args.tmp_dir_days)}coverage_and_an_generation",
        experimental=args.experimental,
        **_get_batch_resource_kwargs(args),
    )

    # ===================================================================
    # ROLE 2: WORKER — one relay container the ROLE 1 orchestrator
    # submitted, re-invoking this script with --run-chunk / --run-merge.
    # Run exactly one unit, then return (BEFORE the try/finally below).
    # _run_coverage_chunk is what calls compute_all_release_stats_per_ref_site.
    # ===================================================================
    if args.run_chunk:
        _run_coverage_chunk(args)
        return
    if args.run_merge:
        _run_coverage_merge(
            args.merge_inputs, args.merge_output, coalesce_to=args.merge_coalesce_to
        )
        return

    # ===================================================================
    # ROLE 3: IN-PROCESS PIPELINE — reached only when NO orchestrator/
    # worker flag is set: the strict single-job compute (gnomAD
    # consent-drop, or dev/test AoU) and the AoU release assembly
    # (--export-* / --merge-qual-hists, run as a separate invocation
    # after fan-out + --merge-cov-chunks). Phases below run in order,
    # each gated by its CLI flag: SETUP -> COMPUTE -> ASSEMBLE.
    # ===================================================================
    try:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
        )
        group_membership_ht_path = _resolve_group_membership_ht_path(
            project,
            environment,
            test=test,
            reduce_min_aggs=args.reduce_min_aggs,
        )

        # AoU downsampling + sample-meta are AoU-only resources (the
        # helpers only allow rwb/batch). Skip them for gnomad runs so the
        # dataproc comparison can resolve resources.
        downsampling_ht_path = None
        meta_ht = None
        if project == "aou":
            downsampling_ht_path = get_aou_downsampling(
                test=test, environment=environment
            ).path
            meta_ht = meta(data_type="genomes", environment=environment).ht()

        # --- SETUP (prerequisites): write the input HTs the COMPUTE
        # phase (and the fan-out workers) read from disk. ---
        if args.write_vep_context_sites:
            logger.info("Writing preprocessed vep_context sites HT...")
            sites_path = _vep_context_sites_path(environment)
            check_resource_existence(
                output_step_resources={"vep_context_sites": [sites_path]},
                overwrite=overwrite,
            )
            _build_vep_context_sites_ht().write(sites_path, overwrite=overwrite)

        if args.write_aou_downsampling_ht:
            logger.info("Writing AoU downsampling HT...")
            check_resource_existence(
                output_step_resources={"downsampling_ht": [downsampling_ht_path]},
                overwrite=overwrite,
            )
            ds_meta_ht = meta_ht.filter(
                (meta_ht.project_meta.project == project) & (meta_ht.release)
            )
            ds_ht = get_downsampling_ht(ds_meta_ht)
            ds_ht.write(downsampling_ht_path, overwrite=overwrite)

        if args.write_group_membership_ht:
            check_resource_existence(
                output_step_resources={
                    "group_membership_ht": [group_membership_ht_path]
                },
                overwrite=overwrite,
            )
            if project == "gnomad":
                logger.info("Writing gnomAD group membership HT...")
                group_membership_ht = get_group_membership_ht(
                    v4_meta(data_type="genomes").ht(),
                    project=project,
                    reduce_min_aggs=args.reduce_min_aggs,
                )
            else:
                logger.info("Writing AoU group membership HT...")
                aou_meta_ht = meta_ht.filter(
                    (meta_ht.project_meta.project == project) & (meta_ht.release)
                )
                group_membership_ht = get_group_membership_ht(
                    meta_ht=aou_meta_ht,
                    project=project,
                    ds_ht=hl.read_table(downsampling_ht_path),
                    reduce_min_aggs=args.reduce_min_aggs,
                )
            group_membership_ht.write(group_membership_ht_path, overwrite=overwrite)

        # --- VALIDATE: every intended vep_context site got AN/coverage after
        # the fan-out + merge. Run only on a FULL run (a partial/--test run
        # covers only its partitions' sites, so the anti-join would flag the
        # rest as missing). ---
        if args.validate_cov_and_an:
            logger.info("Validating cov_and_an HT covers all vep_context sites...")
            sites_ht = hl.read_table(_vep_context_sites_path(environment))
            merged_ht = hl.read_table(cov_and_an_ht_path)
            n_missing = sites_ht.anti_join(merged_ht).count()
            if n_missing:
                raise ValueError(
                    f"cov_and_an HT at {cov_and_an_ht_path} is MISSING {n_missing}"
                    " vep_context site(s) after fan-out + merge — sites were"
                    " dropped. Re-run --use-batch-fanout / --merge-cov-chunks."
                )
            logger.info(
                "Validation passed: all %d vep_context sites present.",
                sites_ht.count(),
            )

        # --- COMPUTE (strict, single job): whole-VDS per-ref-site
        # coverage/AN/qual-hists. Prod AoU does this via ROLE 1 fan-out
        # instead and never reaches here; this is gnomAD and dev/test. ---
        if args.compute_all_cov_release_stats_ht:
            logger.info(
                "Computing coverage, all sites allele number, and optionally"
                " quality histograms HT for %s...",
                project,
            )
            check_resource_existence(
                output_step_resources={"coverage_and_an_ht": [cov_and_an_ht_path]},
                overwrite=overwrite,
            )
            # --test-n-partitions: read only the first N partitions of the VDS
            # for a cheap test (the strict path otherwise reads the whole VDS).
            # partition_count scopes the ref_ht to those partitions' extent when
            # not co-partitioning.
            test_n = args.test_n_partitions
            strict_partition_range = list(range(test_n)) if test_n else None

            # Optionally co-partition the VDS and vep_context reads so the densify
            # join is shuffle-free (Hail-team guidance: read all tables with the
            # same partitions). The full run uses repartition_for_join's 1.1x
            # native-partition multiplier (bare flag); a value gives an explicit
            # count. With --test-n-partitions the intervals are derived from those
            # partitions' extent; otherwise from the (chrom-restricted) vep_context.
            strict_sub_intervals = None
            if args.partitions_for_rep_on_read is not None:
                # 0 (bare flag) -> repartition_for_join's 1.1x multiplier;
                # a positive int -> that explicit partition count.
                n = args.partitions_for_rep_on_read
                logger.info(
                    "Co-partitioning the VDS and vep_context reads (%s).",
                    (
                        f"{n} balanced partitions"
                        if n
                        else "repartition_for_join 1.1x native-partition multiplier"
                    ),
                )
                if test_n:
                    probe = _probe_vds(
                        project, environment, strict_partition_range, chrom
                    )
                    strict_sub_intervals = repartition_for_join(
                        probe.reference_data.rows(),
                        locus_intervals=True,
                        **({"n_partitions": n} if n else {}),
                    )
                else:
                    strict_sub_intervals = _derive_ref_partition_intervals(
                        n, chrom=chrom
                    )
            vds, sex_karyotype_field = _load_project_vds(
                project=project,
                environment=environment,
                partition_range=strict_partition_range,
                sub_intervals=strict_sub_intervals,
                chrom=chrom,
                test=test,
                test_sample_subset=args.test_sample_subset,
            )
            ref_ht = _build_chunk_ref_ht(
                vds_filtered=vds,
                partition_count=test_n or 0,
                chrom=chrom,
                sites_path=_vep_context_sites_path(environment),
                sub_intervals=strict_sub_intervals,
            )

            cov_and_an_ht = compute_all_release_stats_per_ref_site(
                vds,
                ref_ht,
                sex_karyotype_field=sex_karyotype_field,
                project=project,
                group_membership_ht=hl.read_table(group_membership_ht_path),
                reduce_min_aggs=args.reduce_min_aggs,
            )
            cov_and_an_ht = cov_and_an_ht.checkpoint(
                new_temp_file(f"{project}_cov_and_an", "ht")
            )
            if args.n_partitions is not None:
                cov_and_an_ht = cov_and_an_ht.naive_coalesce(args.n_partitions)
            cov_and_an_ht.write(cov_and_an_ht_path, overwrite=overwrite)

        # --- ASSEMBLE (gnomAD): subtract the consent-drop cohort from
        # the v4 release coverage/AN HTs to get gnomAD v5. ---
        if args.merge_gnomad_coverage:
            logger.info("Building gnomAD v5 coverage HT (subtracting consent-drop)...")
            merged_gnomad_coverage_ht_path = (
                f"{qc_temp_prefix(environment=environment)}"
                "gnomad_v5_genomes_coverage.ht"
            )
            check_resource_existence(
                output_step_resources={
                    "merged_gnomad_coverage_ht": [merged_gnomad_coverage_ht_path],
                },
                overwrite=overwrite,
            )
            gnomad_ht = hl.read_table(cov_and_an_ht_path).drop("AN")
            gnomad_release_ht = hl.read_table(
                release_coverage_path(release_version=v4_COVERAGE_RELEASE, public=True)
            )
            if test:
                gnomad_release_ht = _filter_to_locus_bounds(
                    gnomad_release_ht, gnomad_ht
                )
            ht = merge_gnomad_coverage_hts(gnomad_ht, gnomad_release_ht)
            ht.write(merged_gnomad_coverage_ht_path, overwrite=overwrite)

        if args.merge_gnomad_an:
            logger.info("Building gnomAD v5 AN HT (subtracting consent-drop)...")
            merged_gnomad_an_ht_path = (
                f"{qc_temp_prefix(environment=environment)}gnomad_v5_genomes_an.ht"
            )
            check_resource_existence(
                output_step_resources={
                    "merged_gnomad_an_ht": [merged_gnomad_an_ht_path]
                },
                overwrite=overwrite,
            )
            gnomad_ht = hl.read_table(cov_and_an_ht_path).select("AN")
            gnomad_release_ht = hl.read_table(
                release_coverage_path(
                    release_version=v4_AN_RELEASE,
                    public=True,
                    coverage_type="allele_number",
                )
            )
            if test:
                gnomad_release_ht = _filter_to_locus_bounds(
                    gnomad_release_ht, gnomad_ht
                )
            ht = merge_gnomad_an_hts(gnomad_ht, gnomad_release_ht)
            ht.write(merged_gnomad_an_ht_path, overwrite=overwrite)

        # --- ASSEMBLE (AoU): join AoU + gnomAD v5 and export the release
        # HT/TSVs + merge qual hists. Run as a separate invocation after
        # ROLE 1 fan-out + --merge-cov-chunks have produced the AoU
        # cov_and_an HT. ---
        if args.export_coverage_release_files:
            logger.info("Exporting coverage release HT and TSV...")
            cov_ht_path = release_coverage_path(
                public=False,
                test=test,
                coverage_type="coverage",
                environment=environment,
            )
            cov_tsv_path = release_coverage_tsv_path(test=test, environment=environment)
            gnomad_coverage_ht_path = (
                f"{qc_temp_prefix(environment=environment)}"
                "gnomad_v5_genomes_coverage.ht"
            )
            check_resource_existence(
                input_step_resources={"gnomad_coverage_ht": [gnomad_coverage_ht_path]},
                output_step_resources={
                    "cov_release_ht": [cov_ht_path],
                    "cov_tsv": [cov_tsv_path],
                },
                overwrite=overwrite,
            )
            aou_ht = hl.read_table(cov_and_an_ht_path).drop("AN", "qual_hists")
            gnomad_ht = hl.read_table(gnomad_coverage_ht_path)
            ht = join_aou_and_gnomad_coverage_ht(aou_ht, gnomad_ht)
            ht = ht.checkpoint(new_temp_file("aou_and_gnomad_cov_join", "ht"))
            if args.n_partitions is not None:
                ht = ht.naive_coalesce(args.n_partitions)
            ht = ht.checkpoint(cov_ht_path, overwrite=overwrite)
            ht.export(cov_tsv_path)

        if args.export_an_release_files:
            logger.info("Exporting AN release HT and TSV...")
            an_ht_path = release_coverage_path(
                public=False,
                test=test,
                coverage_type="allele_number",
                environment=environment,
            )
            an_tsv_path = release_all_sites_an_tsv_path(
                test=test, environment=environment
            )
            gnomad_an_ht_path = (
                f"{qc_temp_prefix(environment=environment)}gnomad_v5_genomes_an.ht"
            )
            check_resource_existence(
                input_step_resources={"gnomad_an_ht": [gnomad_an_ht_path]},
                output_step_resources={
                    "an_release_ht": [an_ht_path],
                    "an_release_tsv": [an_tsv_path],
                },
                overwrite=overwrite,
            )
            aou_ht = hl.read_table(cov_and_an_ht_path).select("AN")
            gnomad_ht = hl.read_table(gnomad_an_ht_path)
            ht = join_aou_and_gnomad_an_ht(aou_ht, gnomad_ht)
            ht = ht.checkpoint(new_temp_file("aou_and_gnomad_an_join", "ht"))
            ht = ht.select("AN")
            ht = ht.select_globals(
                strata_meta=ht.strata_meta,
                strata_sample_count=ht.strata_sample_count,
            )
            if args.n_partitions is not None:
                ht = ht.naive_coalesce(args.n_partitions)
            ht = ht.checkpoint(an_ht_path, overwrite=overwrite)
            ht = ht.transmute(AN=ht.AN[0])
            ht.export(an_tsv_path)

        if args.merge_qual_hists:
            logger.info("Merging AoU + gnomAD v4 qual hists...")
            qual_hists_path = qual_hists(test=test, environment=environment).path
            qual_hists_path = _apply_path_suffix(
                qual_hists_path, args.qual_hists_output_suffix
            )
            check_resource_existence(
                output_step_resources={"qual_hists_ht": [qual_hists_path]},
                overwrite=overwrite,
            )
            aou_ht = (
                hl.read_table(cov_and_an_ht_path).select("qual_hists").select_globals()
            )
            gnomad_ht = (
                public_release(data_type="genomes")
                .ht()
                .select("histograms")
                .select_globals()
            )
            if test:
                gnomad_ht = _filter_to_locus_bounds(gnomad_ht, aou_ht)
            # Drop age hists (handled in the frequency script) and rekey
            # by locus. `distinct()` deduplicates the multi-row-per-locus
            # rows that result from rekeying a (locus, alleles)-keyed HT
            # to locus-only; the `_all` hists are locus-level and
            # identical across split rows, so distinct keeps a representative.
            gnomad_ht = gnomad_ht.select(
                qual_hists=gnomad_ht.histograms.drop("age_hists")
            )
            gnomad_ht = gnomad_ht.key_by("locus").distinct()
            ht = join_aou_and_gnomad_qual_hists_ht(aou_ht, gnomad_ht)
            ht.write(qual_hists_path, overwrite=overwrite)
    finally:
        if environment != "batch":
            logger.info("Copying hail log to logging bucket...")
            hl.copy_log(
                get_logging_path(
                    "v5_coverage_and_an_generation",
                    environment=environment,
                    tmp_dir_days=args.tmp_dir_days,
                )
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-name",
        help="Project name.",
        default="aou",
        type=str,
        choices=["aou", "gnomad"],
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )

    # Environment configuration.
    env_group = parser.add_argument_group("environment configuration")
    env_group.add_argument(
        "--environment",
        help="Compute environment.",
        choices=["rwb", "batch", "dataproc"],
        default="rwb",
    )
    env_group.add_argument(
        "--tmp-dir-days",
        type=int,
        default=4,
        help="Number of days for temp directory retention. Default is 4.",
    )
    env_group.add_argument(
        "--gcp-billing-project",
        type=str,
        default="broad-mpg-gnomad",
        help="Google Cloud billing project for reading requester pays buckets.",
    )

    # Batch-specific configuration.
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
    )
    batch_group.add_argument(
        "--experimental",
        action="store_true",
        help=(
            "Route the QoB init through `hl.experimental.init` instead of"
            " `hl.init` and attach the QoB driver to an existing Hail Batch."
            " Requires HAIL_BATCH_ID to be set in the env (Hail Batch"
            " injects this inside batch jobs); raises if not. Without this"
            " flag, each invocation creates its own Hail Batch."
        ),
    )
    batch_group.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_group.add_argument(
        "--driver-cores",
        type=str,
        default=None,
        help=(
            "Number of cores for the driver node. Pass a power of two"
            " between 0.25 and 16 (as a string, e.g. '2' or '0.5')."
        ),
    )
    batch_group.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_group.add_argument(
        "--jvm-heap-size",
        type=str,
        default=None,
        help=(
            "Max JVM heap size (-Xmx) for the in-process QoB driver under"
            " --experimental, e.g. '5g' or '2500m'. Plumbed to"
            " hl.experimental.init(jvm_heap_size=...). Set to ~50-70%% of"
            " container memory; the rest is for native off-heap"
            " (RegionPool), Python, and OS overhead. Ignored without"
            " --experimental."
        ),
    )
    batch_group.add_argument(
        "--worker-cores",
        type=str,
        default=None,
        help=(
            "Number of cores per worker node. Pass a power of two"
            " between 0.25 and 16 (as a string, e.g. '1' or '0.5')."
        ),
    )
    batch_group.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
    )
    parser.add_argument(
        "--n-partitions",
        help=(
            "Number of partitions to use for the output Table. If unset"
            " (default), no final naive_coalesce is applied — the output"
            " keeps the IR's natural partition count (strict path) or"
            " the merge tree's per-level coalesce shape (--merge-cov-chunks)."
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--cov-and-an-output-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix appended to the cov_and_an HT path (before the"
            " .ht extension). Use to write the output to a sibling location"
            " for A/B comparison. Applies to both writes (step 3) and reads"
            " (downstream steps), so pass the same suffix consistently."
        ),
    )
    parser.add_argument(
        "--qual-hists-output-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix appended to the qual_hists HT path (before the"
            " .ht extension) — analogous to --cov-and-an-output-suffix but"
            " for the --merge-qual-hists output. Use for A/B comparison."
        ),
    )
    parser.add_argument(
        "--reduce-min-aggs",
        action="store_true",
        help=(
            "Use the leaf-only `reduce_to_minimal_groups` optimization: build"
            " the group_membership_ht with `reduce_to_minimal_groups=True`"
            " and pass `reducible_aggs={'AN'}` to `compute_stats_per_ref_site`."
            " Must be set consistently for both --write-group-membership-ht"
            " and --compute-all-cov-release-stats-ht runs (the gmh's"
            " `freq_reduced` global determines whether reduction takes effect"
            " at compute time)."
        ),
    )
    parser.add_argument(
        "--test-sample-subset",
        action="store_true",
        help=(
            "When `--test` is set on AoU, additionally subsample the AoU"
            " sample set to ~0.1%% via `meta_ht.sample(0.001)`. Default"
            " off: a `--test` run uses all samples but the partition /"
            " chrom subset, so AN values are stable and comparable"
            " across runs. Useful for a cheap tiny-cohort sanity check"
            " rather than a real-scale slice."
        ),
    )

    test_group = parser.add_argument_group(
        "Test mode",
        "Route inputs/outputs to test paths and (optionally) filter data.",
    )
    test_group.add_argument(
        "--test",
        action="store_true",
        help=(
            "Route reads/writes to test paths (group_membership,"
            " cov_and_an, release HTs, etc.). Independent of"
            " ``--chrom``: combine for a chrom-filtered test run."
        ),
    )
    parser.add_argument(
        "--chrom",
        nargs="+",
        default=None,
        help=(
            "Filter input data to these contigs (e.g."
            " ``--chrom chr22 chrX chrY``). Default: no chrom filter."
            " Independent of --test."
        ),
    )

    group_membership_args = parser.add_argument_group(
        "Get gnomAD genomes group membership HT.",
    )
    group_membership_args.add_argument(
        "--write-group-membership-ht",
        help="Write group membership HT.",
        action="store_true",
    )

    parser.add_argument(
        "--write-aou-downsampling-ht",
        help="Write v5 downsampling HT.",
        action="store_true",
    )
    parser.add_argument(
        "--write-vep-context-sites",
        help=(
            "Preprocess vep_context once into a deduped, telomere/centromere-"
            "stripped, locus-keyed sites HT (30-day storage) that every chunk reads"
            " co-partitioned. Prerequisite for the compute; run before"
            " --use-batch-fanout / --compute-all-cov-release-stats-ht."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--validate-cov-and-an",
        help=(
            "Anti-join the vep_context sites HT against the merged cov_and_an HT and"
            " fail if any site is missing (i.e. dropped during fan-out + merge)."
            " Run after --merge-cov-chunks on a FULL run only (a partial/--test run"
            " covers only its partitions' sites)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-all-cov-release-stats-ht",
        help="Compute the all sites coverage, allele number, and quality histogram HT.",
        action="store_true",
    )
    parser.add_argument(
        "--partitions-for-rep-on-read",
        type=int,
        nargs="?",
        const=0,
        default=None,
        help=(
            "Strict single-job compute (--compute-all-cov-release-stats-ht) only:"
            " co-partition the VDS and vep_context reads (via repartition_for_join)"
            " so the densify join is shuffle-free. Pass the flag with NO value to"
            " use repartition_for_join's default 1.1x-of-native-partitions"
            " multiplier (recommended for the full run -- no number to guess), or"
            " pass an explicit integer partition count (handy for small tests)."
            " Omit the flag entirely to read each with its native partitioning."
            " Honors --chrom. Used for the single gnomAD dataproc job, which is not"
            " chunked/fanned out."
        ),
    )
    parser.add_argument(
        "--test-n-partitions",
        type=int,
        default=None,
        help=(
            "Strict single-job compute (--compute-all-cov-release-stats-ht) only:"
            " read just the first N partitions of the VDS for a cheap test (the"
            " strict path otherwise reads the whole VDS; replaces the dropped"
            " --test-2-partitions). The ref_ht is scoped to those partitions'"
            " locus extent. Combine with --partitions-for-rep-on-read to co-partition just"
            " those N partitions. Default None reads the whole VDS."
        ),
    )

    coverage_args = parser.add_argument_group(
        "Compute coverage release stats HT.",
    )
    coverage_args.add_argument(
        "--merge-gnomad-coverage",
        help="Subtract consent drop samples from v4 release HT to create gnomAD v5 genomes coverage HT.",
        action="store_true",
    )
    coverage_args.add_argument(
        "--export-coverage-release-files",
        help="Join and export AoU + gnomAD v4 coverage release HT and TSV file.",
        action="store_true",
    )

    an_args = parser.add_argument_group(
        "Compute AN release stats HT.",
    )
    an_args.add_argument(
        "--merge-gnomad-an",
        help="Subtract consent drop samples from v4 release HT to create gnomAD v5 genomes AN HT.",
        action="store_true",
    )
    an_args.add_argument(
        "--export-an-release-files",
        help="Exports joint AoU + gnomAD v4 AN release HT and TSV file.",
        action="store_true",
    )

    parser.add_argument(
        "--merge-qual-hists",
        help="Merge variant quality histograms from AoU v8 and gnomAD v4 genomes.",
        action="store_true",
    )

    fanout_group = parser.add_argument_group(
        "Hail Batch fan-out for `--compute-all-cov-release-stats-ht`.",
        "Submit one Hail Batch job per partition chunk; each chunk is a tiny"
        " relay container that spawns its own QoB driver to do the densify."
        " Idempotent: chunks whose `_SUCCESS` exists are skipped on rerun.",
    )
    fanout_group.add_argument(
        "--use-batch-fanout",
        action="store_true",
        help=(
            "Build and submit a Hail Batch pipeline that fans out the"
            " coverage/AN compute by partition chunks (QoB-from-container"
            " per chunk). Mutually exclusive with"
            " --compute-all-cov-release-stats-ht, --merge-cov-chunks,"
            " --run-chunk, and --run-merge. After all chunks complete,"
            " run --merge-cov-chunks to union them into the final HT."
        ),
    )
    fanout_group.add_argument(
        "--merge-cov-chunks",
        action="store_true",
        help=(
            "Build and submit the two-level tree-reduce merge of per-chunk"
            " HTs produced by --use-batch-fanout. Level 1 (group merges,"
            " sized by --merge-group-size) writes to"
            " <cov_and_an_path>_merge_groups/<idx>.ht; level 2 unions all"
            " group HTs into <cov_and_an_path>. Chunk discovery uses"
            " --total-partitions / --partitions-per-chunk (same as"
            " --use-batch-fanout); missing chunks fail loudly. Requires"
            " --environment=batch."
        ),
    )
    fanout_group.add_argument(
        "--partitions-per-chunk",
        type=int,
        default=3,
        help=(
            "Number of VDS partitions per chunk job. Default 3 (the"
            " driver-comfort load used for the prod AoU runs). IMPORTANT:"
            " this must be identical across the compute"
            " (--use-batch-fanout) and merge (--merge-cov-chunks) reruns."
            " Chunk discovery derives the expected chunk count from"
            " --total-partitions / --partitions-per-chunk, so passing a"
            " different value on the merge silently unions the wrong set of"
            " chunk HTs (no error) and skews the expected-chunk check."
        ),
    )
    fanout_group.add_argument(
        "--read-subintervals-per-chunk",
        type=int,
        default=50,
        help=(
            "Per-chunk: subdivide the chunk's locus extent into this many"
            " equal-position sub-intervals, then re-read both the VDS"
            " (`hl.vds.read_vds(intervals=...)`) and `vep_context`"
            " (`hl.read_table(_intervals=...)`) with those intervals. One"
            " partition per sub-interval on both sides, co-partitioned,"
            " no shuffle. Default 50. Set to 1 to skip the probe/re-read"
            " and use plain `filter_partitions` (legacy behavior)."
        ),
    )
    fanout_group.add_argument(
        "--total-partitions",
        type=int,
        default=145192,
        help=(
            "Total VDS partition count to scatter across. Default 145192"
            " (prod AoU VDS partition count). Used by --use-batch-fanout."
        ),
    )
    fanout_group.add_argument(
        "--wave-size",
        type=int,
        default=1000,
        help=(
            "Max chunks (= relay jobs) submitted per Hail Batch under"
            " --use-batch-fanout. Pending chunks are split into"
            " sequential waves of this size; each wave's Hail Batch runs"
            " to completion before the next is submitted, so the number"
            " of concurrently-running chunk relays (and their nested QoB"
            " drivers) is bounded by --wave-size. Prevents a single"
            " ~tens-of-thousands-job batch from overwhelming the shared"
            " Hail Batch service. Note: these jobs run on the Hail Batch"
            " service's worker pool (the Hail team's GCP project), not"
            " gnomad-production Compute/Dataproc, so the binding limit is"
            " the shared service's capacity / any per-billing-project core"
            " cap set by the Hail team — not gnomad-production's own GCP"
            " quota. Waves are restartable: a"
            " rerun skips chunks whose _SUCCESS already exists. Set <= 0"
            " or >= total chunk count for a single batch (legacy"
            " behavior). Default 1000."
        ),
    )
    fanout_group.add_argument(
        "--merge-group-size",
        type=int,
        default=500,
        help="Chunk HTs per group-merge job. Default 500.",
    )
    fanout_group.add_argument(
        "--methods-branch",
        type=str,
        default="main",
        help="Branch/commit of gnomad_methods to clone in chunk containers.",
    )
    fanout_group.add_argument(
        "--batch-image",
        type=str,
        default=(
            "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:latest"
        ),
        help="Docker image for the chunk + merge BashJob containers.",
    )
    fanout_group.add_argument(
        "--batch-billing-project",
        type=str,
        default="gnomad-production",
        help="Hail Batch billing project for the orchestrator pipeline.",
    )
    fanout_group.add_argument(
        "--batch-remote-tmpdir",
        type=str,
        default=None,
        help=(
            "GCS scratch path for hb.ServiceBackend. Defaults to the value"
            " encoded in `_build_setup_command` config."
        ),
    )
    fanout_group.add_argument(
        "--chunk-cpu",
        type=float,
        default=0.5,
        help="CPU per chunk relay container. Default 0.5 (relay only).",
    )
    fanout_group.add_argument(
        "--chunk-memory",
        type=str,
        default="standard",
        choices=["lowmem", "standard", "highmem"],
        help="Hail Batch memory preset per chunk relay container.",
    )
    fanout_group.add_argument(
        "--chunk-storage",
        type=str,
        default="5Gi",
        help="/io storage per chunk relay container. Default 5Gi.",
    )
    fanout_group.add_argument(
        "--chunk-driver-cores",
        type=str,
        default=None,
        help=(
            "Cores for the QoB driver pod each chunk/merge relay spawns."
            " Power of two between 0.25 and 16 (as a string, e.g. '2' or"
            " '0.5'). Forwarded to the relay's --driver-cores; decoupled"
            " from the orchestrator's own --driver-cores."
        ),
    )
    fanout_group.add_argument(
        "--chunk-driver-memory",
        type=str,
        default=None,
        help=(
            "Memory profile for the QoB driver pod each chunk/merge"
            " relay spawns, e.g. 'highmem'. Forwarded to the relay's"
            " --driver-memory; decoupled from the orchestrator's own"
            " --driver-memory."
        ),
    )
    fanout_group.add_argument(
        "--chunk-worker-cores",
        type=str,
        default=None,
        help=(
            "Cores per QoB worker pod the chunk/merge driver dispatches."
            " Power of two between 0.25 and 16 (as a string, e.g. '1' or"
            " '0.5'). Forwarded to the relay's --worker-cores; decoupled"
            " from the orchestrator's own --worker-cores."
        ),
    )
    fanout_group.add_argument(
        "--chunk-worker-memory",
        type=str,
        default=None,
        help=(
            "Memory profile per QoB worker pod the chunk/merge driver"
            " dispatches, e.g. 'lowmem'. Forwarded to the relay's"
            " --worker-memory; decoupled from the orchestrator's own"
            " --worker-memory."
        ),
    )
    fanout_group.add_argument(
        "--chunk-attempts",
        type=int,
        default=5,
        help="Max retry attempts (n_max_attempts) per chunk job.",
    )
    fanout_group.add_argument(
        "--merge-cpu",
        type=int,
        default=4,
        help="CPU per merge container.",
    )
    fanout_group.add_argument(
        "--merge-memory",
        type=str,
        default="standard",
        help="Hail Batch memory preset per merge container.",
    )
    fanout_group.add_argument(
        "--merge-storage",
        type=str,
        default="50Gi",
        help="Extra /io storage per group-merge container.",
    )
    fanout_group.add_argument(
        "--final-merge-storage",
        type=str,
        default="100Gi",
        help="Extra /io storage for the final merge container.",
    )
    fanout_group.add_argument(
        "--batch-dry-run",
        action="store_true",
        help="If set, validate the Hail Batch DAG without submitting.",
    )

    worker_group = parser.add_argument_group(
        "Per-chunk worker entry-points for --use-batch-fanout (set by the"
        " orchestrator; not for direct use unless debugging a single chunk).",
    )
    worker_group.add_argument(
        "--run-chunk",
        action="store_true",
        help=(
            "Worker entry-point: compute one chunk and write to"
            " --chunk-output. Requires --chunk-start, --chunk-stop,"
            " --chunk-output."
        ),
    )
    worker_group.add_argument(
        "--chunk-start",
        type=int,
        default=0,
        help="First VDS partition index for --run-chunk (inclusive).",
    )
    worker_group.add_argument(
        "--chunk-stop",
        type=int,
        default=0,
        help="Last VDS partition index for --run-chunk (exclusive).",
    )
    worker_group.add_argument(
        "--chunk-output",
        type=str,
        default=None,
        help=(
            "GCS output path for --run-chunk's per-chunk HT. Optional: when"
            " unset, auto-derived from the resolved cov_and_an HT path"
            " (plus `--cov-and-an-output-suffix` if set) and the chunk's"
            " `--chunk-start` index (`_chunks/<chunk_start:08d>.chunk.ht`)."
            " The orchestrator (`--use-batch-fanout`) always passes this"
            " explicitly; this default fires only for standalone"
            " `--run-chunk` invocations (smoketests)."
        ),
    )
    worker_group.add_argument(
        "--run-merge",
        action="store_true",
        help=(
            "Worker entry-point: union the HTs at --merge-inputs and write"
            " to --merge-output."
        ),
    )
    worker_group.add_argument(
        "--merge-inputs",
        type=str,
        nargs="+",
        default=None,
        help="GCS paths of per-chunk (or per-group) HTs to union.",
    )
    worker_group.add_argument(
        "--merge-output",
        type=str,
        default=None,
        help="GCS output path for --run-merge's merged HT.",
    )
    worker_group.add_argument(
        "--merge-coalesce-to",
        type=int,
        default=None,
        help=(
            "If set, the --run-merge worker naive_coalesces the unioned"
            " HT to this many partitions before writing. Set by the merge"
            " orchestrator to bound per-level partition count (group"
            " merges coalesce to len(inputs); the final merge coalesces"
            " to --n-partitions when that is set)."
        ),
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    # Every argument that only makes sense with --environment=batch lives in
    # this one list so the validation is complete and easy to scan. "Provided"
    # is detected by comparing against the parser default, which works
    # uniformly for None-default args, value-default args (--chunk-cpu /
    # -memory / -storage), and store_true flags (--experimental and the
    # --use-batch-fanout / --merge-cov-chunks mode flags) — so no arg has to
    # be left out just because it has a non-None default.
    batch_only_args = [
        "app_name",
        "driver_cores",
        "driver_memory",
        "jvm_heap_size",
        "chunk_driver_cores",
        "chunk_driver_memory",
        "worker_cores",
        "worker_memory",
        "chunk_worker_cores",
        "chunk_worker_memory",
        "chunk_cpu",
        "chunk_memory",
        "chunk_storage",
        "experimental",
        "use_batch_fanout",
        "merge_cov_chunks",
    ]
    provided_batch_args = [
        a for a in batch_only_args if getattr(args, a) != parser.get_default(a)
    ]
    if provided_batch_args and args.environment != "batch":
        parser.error(
            "Batch-only arguments ("
            + ", ".join("--" + a.replace("_", "-") for a in provided_batch_args)
            + ") require --environment=batch"
        )

    # --jvm-heap-size only applies to the in-process JVM under --experimental.
    if args.jvm_heap_size is not None and not args.experimental:
        parser.error(
            "--jvm-heap-size requires --experimental (it controls the"
            " in-process JVM heap for hl.experimental.init)."
        )

    if args.project_name == "aou" and args.environment not in ("rwb", "batch"):
        parser.error(
            "--project-name=aou requires --environment to be 'rwb' or 'batch'."
        )
    if args.project_name == "gnomad" and args.environment not in ("dataproc", "batch"):
        parser.error(
            "--project-name=gnomad requires --environment to be 'dataproc' or 'batch'."
        )

    if args.write_aou_downsampling_ht and args.project_name != "aou":
        raise ValueError(
            "--write-aou-downsampling-ht requires --project-name to be 'aou'."
        )

    if args.merge_gnomad_coverage and args.project_name != "gnomad":
        raise ValueError(
            "--merge-gnomad-coverage requires --project-name to be 'gnomad'."
        )

    if args.merge_gnomad_an and args.project_name != "gnomad":
        raise ValueError("--merge-gnomad-an requires --project-name to be 'gnomad'.")

    if args.export_coverage_release_files and args.project_name != "aou":
        raise ValueError(
            "--export-coverage-release-files requires --project-name to be 'aou'."
        )

    if args.export_an_release_files and args.project_name != "aou":
        raise ValueError(
            "--export-an-release-files requires --project-name to be 'aou'."
        )

    if args.merge_qual_hists and args.project_name != "aou":
        raise ValueError("--merge-qual-hists requires --project-name to be 'aou'.")

    # --use-batch-fanout / --merge-cov-chunks / --run-chunk / --run-merge
    # are mutually exclusive entry-points. The orchestrators set
    # --run-chunk / --run-merge inside their relay containers; users
    # should not pass the worker flags from the command line alongside
    # the orchestrator flags.
    fanout_modes = sum(
        bool(x)
        for x in (
            args.use_batch_fanout,
            args.merge_cov_chunks,
            args.run_chunk,
            args.run_merge,
        )
    )
    if fanout_modes > 1:
        parser.error(
            "--use-batch-fanout, --merge-cov-chunks, --run-chunk, and"
            " --run-merge are mutually exclusive."
        )
    if args.use_batch_fanout and args.compute_all_cov_release_stats_ht:
        parser.error(
            "--use-batch-fanout is an alternative to"
            " --compute-all-cov-release-stats-ht; do not pass both."
        )
    if args.merge_cov_chunks and args.compute_all_cov_release_stats_ht:
        parser.error(
            "--merge-cov-chunks is a separate step from"
            " --compute-all-cov-release-stats-ht; do not pass both."
        )
    if args.run_chunk:
        if args.chunk_stop <= args.chunk_start:
            parser.error("--run-chunk requires --chunk-stop > --chunk-start.")
        # ``--chunk-output`` is optional: when not set it's auto-derived in
        # ``_run_coverage_chunk`` from the resolved cov_and_an HT path +
        # `--cov-and-an-output-suffix` + chunk_start. Useful for one-off
        # `--run-chunk` invocations (smoketests) — the orchestrator path
        # always passes `--chunk-output` explicitly.
    if args.run_merge:
        if args.merge_output is None or not args.merge_inputs:
            parser.error("--run-merge requires --merge-output and --merge-inputs.")

    main(args)
