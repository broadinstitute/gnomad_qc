r"""
Compute coverage, allele number, and quality histograms for gnomAD v5 genomes.

v5 genomes = AoU v8 (new) + gnomAD v4 genomes minus the consent-drop set.
This script computes per-reference-site coverage, AN, and (AoU only) qual
histograms for each project and joins them to create the v5 coverage and AN HT/TSV.

Execution roles:
----------------
This script dispatches three process roles by ``main`` on the CLI flags.
They are mutually exclusive -- dispatched by early return in ``main`` (the
``__main__`` block also rejects conflicting entry-point flags):

1. ORCHESTRATOR (``--use-batch-fanout`` or ``--merge-cov-chunks``):
   runs wherever you launch it and never initializes Hail (its
   prerequisite / ``_SUCCESS`` skip-checks use ``file_exists``, which is
   hailtop.fs -- no Hail backend), builds and submits a Hail Batch of
   relay jobs, then returns. Used for the prod-scale AoU compute (the
   dataset is too big for one job).
2. WORKER (``--run-chunk`` or ``--run-merge``): each relay container the
   orchestrator submitted re-invokes this script with ``--run-chunk`` /
   ``--run-merge``. Initializes Hail, does exactly one chunk/merge unit,
   and returns before the try/finally. The chunk worker is what actually
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
    2b. (AoU only) Write the sample artifacts (--write-sample-artifacts):
        collisions.json + release_samples.json. Without them every chunk worker
        scans the ~330-partition sample-collisions HT and collects the
        ~330-partition meta HT inside get_aou_vds. Optional (the chunk falls back
        to those scans), but at fan-out scale it is paid per chunk.
    3. Write the preprocessed vep_context sites HT (--write-vep-context-sites):
        deduped/telomere+centromere+chrM-stripped/locus-keyed sites read co-partitioned
        by every chunk (so the dedup + strip runs once, not per chunk). This table is
        the definition of "every site the compute must cover" -- step 5 validates the
        merged HT against it exactly.
    4. Compute the per-ref-site coverage/AN/qual-hists HT — the dense step.
       Strict single-job via --compute-all-cov-release-stats-ht, or the
       prod-scale fan-out: --write-chunk-intervals first (one VDS open
       precomputes every chunk's read sub-intervals; required -- the fan-out
       and merge read this JSON to enumerate chunks), then --use-batch-fanout,
       then --merge-cov-chunks. See "Execution roles" above for how those dispatch.
       Chunk boundaries are sampled without a fixed seed, so each
       --write-chunk-intervals run yields a different layout; chunk outputs are
       therefore namespaced by a content hash of that JSON
       (``_chunks/<hash>/<idx>.chunk.ht``). Re-running the fan-out reuses the
       existing JSON (same hash) and retries only the missing chunks; regenerating
       the JSON produces a new hash, so chunks from the old layout are recomputed
       under the new one rather than silently merged into overlapping loci.
    5. Validate the merged HT covers every vep_context site
        (--validate-cov-and-an) — anti-join, fail on any dropped site.

Per-project merge (--project-name gnomad)::

    6. Build the gnomAD v5 coverage HT (--merge-gnomad-coverage) and AN HT
       (--merge-gnomad-an) by subtracting the consent-drop samples from the
       v4 release HTs.

Release assembly (--project-name aou)::

    7. Join AoU + gnomAD v5 outputs and export release files:
       --export-coverage-release-files (coverage HT + TSV; gnomAD-only --
         AoU computes no coverage stats, its ref blocks carry no DP)
       --export-an-release-files (AN HT + TSV)
       --merge-qual-hists (qual hists HT; gnomAD v4 hists are reused as-is —
       they were not recomputed for v5).

Usage Examples::

    # Per-chunk fan-out compute + merge (on batch). Run setup steps 1-3 +
    # --write-chunk-intervals first; the fan-out reads those outputs:
    python compute_coverage.py --project-name aou --environment batch \\
        --use-batch-fanout --partitions-per-chunk 3
    python compute_coverage.py --project-name aou --environment batch \\
        --merge-cov-chunks --n-partitions 10000

    # gnomAD v4 consent-drop subtraction (on dataproc):
    hailctl dataproc submit cluster-name compute_coverage.py \\
        --project-name gnomad --environment dataproc \\
        --merge-gnomad-coverage --merge-gnomad-an

    # AoU/gnomAD release assembly (on batch):
    python compute_coverage.py --project-name aou --environment batch \\
        --export-coverage-release-files --export-an-release-files \\
        --merge-qual-hists
"""

import argparse
import hashlib
import json
import logging
import os
import re
import shlex
import subprocess
import sys
from collections.abc import Sequence
from datetime import datetime, timezone
from functools import reduce
from itertools import groupby
from typing import Any, NamedTuple

import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs
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
    COVERAGE_OVER_X_BINS,
    annotate_downsamplings,
    build_freq_stratification_list,
    generate_freq_group_membership_array,
    merge_array_expressions,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.file_utils import file_exists, repartition_for_join
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

# Overall groups that non-reducible aggregations (qual_hists →
# {"group": "raw"}; coverage_stats → {"group": "adj"}, gnomAD runs only)
# are pinned to via ``entry_agg_group_membership``. Under ``--reduce-min-aggs`` these must be
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

# gnomAD sample counts
GNOMAD_SAMPLE_COUNT = 71702
GNOMAD_CONSENT_DROP_SAMPLE_COUNT = 849

# Contigs excluded from the all-sites compute. The mitochondria are called by a
# separate pipeline, and the v4 genomes all-sites AN release contains no chrM, so the
# gnomAD side of the release join has nothing to merge there. ``vep_context`` v105 DOES
# include chrM, so both the sites HT (:func:`_build_vep_context_sites_ht`) and the
# per-chunk interval precompute (:func:`_build_chunk_intervals`) must drop these:
# leaving chrM in the sites HT would demand coverage the fan-out never produces and
# ``--validate-cov-and-an`` would fail at the end of the run, while leaving it in the
# chunk set would spawn chunks whose sites read returns nothing.
EXCLUDED_CONTIGS = ["chrM"]

# Partition count for the written group_membership HT. It is one row per sample
# (~365k for AoU) and exists only to build a broadcast per-strata index, so it needs
# far fewer partitions than the ~330 it inherits from the meta HT it is built from.
GROUP_MEMBERSHIP_N_PARTITIONS = 10


def get_downsampling_ht(ht: hl.Table) -> hl.Table:
    """
    Get Table with downsampling groups for all samples.

    v5 downsampling is only applied to the AoU dataset.
    Desired groups:
    - 10,000
    - 100,000
    - Genetic ancestry group sizes for AFR, AMR, NFE

    Only AFR, AMR, and NFE per-group downsamplings are generated: the full gen_anc
    expression is passed to ``annotate_downsamplings`` with
    ``gen_ancs_to_downsample=["afr", "amr", "nfe"]``, which adds those groups'
    sample counts to the downsamplings list. Every group still gets a per-group
    index; the global 10k/100k cover all samples (gen-anc-independent global index).

    :param ht: Input Table.
    :return: Table with downsampling groups.
    """
    logger.info(
        "Determining downsampling groups for AoU...",
    )
    downsamplings = DOWNSAMPLINGS["v5"]
    # Restrict the per-group downsamplings to the desired genetic ancestry groups
    # via gen_ancs_to_downsample so annotate_downsamplings doesn't generate them for
    # every other group too. The global 10k/100k still cover all samples, and every
    # group still gets a per-group index (just no per-group downsampling).
    ht = annotate_downsamplings(
        ht,
        downsamplings,
        ht.genetic_ancestry_inference.gen_anc,
        gen_ancs_to_downsample=["afr", "amr", "nfe"],
    )
    return ht


def get_group_membership_ht(
    meta_ht: hl.Table,
    project: str,
    ds_ht: hl.Table | None = None,
    reduce_min_aggs: bool = False,
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    :param meta_ht: Meta HT.
    :param project: Project name. Must be "aou" or "gnomad". If "gnomad", function will filter meta HT to only consent drop samples.
    :param ds_ht: Optional downsampling HT. Only used for AoU.
    :param reduce_min_aggs: If True, build the HT with ``reduce_to_minimal_groups=True``. Downstream ``compute_stats_per_ref_site`` calls must then list every non-reducible annotation under ``entry_agg_group_membership`` and pass ``reducible_aggs={"AN"}`` (or similar) so AN gets expanded back to full strata.
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

    else:
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

    return ht


def _chunk_intervals_hash(data: dict[str, Any]) -> str:
    """
    Return a stable short content hash of the chunk-intervals JSON's boundaries.

    The chunk boundaries come from ``repartition_for_join`` ->
    ``ht._calculate_new_partitions``, which samples the key distribution WITHOUT a
    fixed seed, so identical ``(flags, project, VDS)`` inputs produce different cut
    points on each ``--write-chunk-intervals`` run. Because chunk outputs are keyed
    only by index and the ``_SUCCESS`` skip-check / merge are existence-only, a chunk
    left at index N by a PRIOR run (different boundaries) would otherwise pass as
    "present" and be merged with this run's chunks -> overlapping/duplicate loci,
    caught only by the final row-count validation. Namespacing every chunk output by
    this hash (see :func:`_chunk_path`) keeps layouts from different runs in separate
    directories, so the skip-check and merge -- which list only the current hash's
    directory -- never mix them; a stale chunk is recomputed rather than silently
    merged.

    Computed over the boundary-determining content (everything except any embedded
    ``intervals_hash``), so re-loading the same JSON always yields the same value.

    :param data: Parsed chunk-intervals JSON (from ``--write-chunk-intervals``).
    :return: First 16 hex chars of the SHA-256 of the canonical serialization.
    """
    payload = {k: v for k, v in data.items() if k != "intervals_hash"}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


def _chunk_path(cov_and_an_ht_path: str, idx: int, intervals_hash: str) -> str:
    """
    Return a sibling per-chunk HT path under ``<cov_and_an_path>_chunks/<hash>/``.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param idx: Chunk index (zero-based).
    :param intervals_hash: Content hash of the chunk-intervals layout that produced
        this chunk (see :func:`_chunk_intervals_hash`); namespaces the output so
        chunks from a different layout can never be mixed at merge time.
    :return: ``<cov_and_an_path>_chunks/<hash>/<idx:08d>.chunk.ht``.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_chunks/{intervals_hash}/{idx:08d}.chunk.ht"


def _list_present_chunk_indices(
    cov_and_an_ht_path: str, intervals_hash: str
) -> set[int]:
    """
    Return the set of chunk indices with a completed (``_SUCCESS``) output.

    One ``hailtop.fs.ls`` glob of ``<…>_chunks/<hash>/*/_SUCCESS`` instead of a
    per-chunk ``file_exists`` probe (tens of thousands of serial GCS stats at prod
    scale). Scoping the glob to ``intervals_hash`` means only chunks from the CURRENT
    layout are seen -- a stale chunk from a prior run (different boundaries, different
    hash) is invisible, so it is recomputed rather than merged (see
    :func:`_chunk_intervals_hash`). Keying on ``_SUCCESS`` (not the dir) means a
    partially-written chunk is correctly treated as absent; a missing directory
    returns an empty set.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param intervals_hash: Content hash of the current chunk-intervals layout.
    :return: Set of completed chunk indices for this layout.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    present: set[int] = set()
    for entry in hfs.ls(f"{base}_chunks/{intervals_hash}/*/_SUCCESS"):
        m = re.search(r"/(\d+)\.chunk\.ht/_SUCCESS$", entry.path)
        if m:
            present.add(int(m.group(1)))
    return present


def _failed_chunks_path(cov_and_an_ht_path: str, intervals_hash: str) -> str:
    """
    Return the path of the failed-chunk manifest for this layout.

    Sibling of the chunk outputs and namespaced by ``intervals_hash`` for the same
    reason they are: a manifest of indices only means something against the layout
    that produced those indices.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param intervals_hash: Content hash of the current chunk-intervals layout.
    :return: ``<cov_and_an_path>_chunks/<hash>/_failed_chunks.json``.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_chunks/{intervals_hash}/_failed_chunks.json"


def _write_failed_chunks_manifest(
    cov_and_an_ht_path: str,
    intervals_hash: str,
    failed: Sequence[int],
    n_dispatched: int,
    run_id: str,
    commit: str,
    app_name: str | None,
    waves: Sequence[dict[str, Any]],
) -> str:
    """
    Record the chunk indices that did not land, for a later targeted rerun.

    Relay jobs do not retry themselves (see ``--chunk-attempts``), so the set of
    chunks still missing when the orchestrator exits is the rerun list. The
    orchestrator's own log holds it too, but a log scrollback is not a durable
    record -- this writes it next to the chunks so the rerun list survives the
    orchestrator process.

    Rerunning ``--use-batch-fanout`` unchanged picks these up automatically (it
    dispatches whatever is not already present), so the manifest is informational:
    it tells you how much is left and which indices to look at, and is rewritten
    (not appended to) on each wave so it always reflects the latest attempt.

    Because it is one fixed path per layout, successive runs overwrite each other
    -- so every manifest is stamped with the orchestrator ``run_id``, the wall-clock
    write time, the code commit, the ``--app-name``, and the Hail Batch id of each
    wave. Without those, an overwritten manifest is indistinguishable from a stale
    one left by an earlier attempt at the same region.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param intervals_hash: Content hash of the current chunk-intervals layout.
    :param failed: Chunk indices dispatched in this run that have no ``_SUCCESS``.
    :param n_dispatched: Number of chunks dispatched so far in this run.
    :param run_id: Identifier for this orchestrator run (see ``_new_run_id``).
    :param commit: gnomad_qc commit the relays ran.
    :param app_name: ``--app-name`` passed to the relays (names their QoB batches).
    :param waves: Per-wave records: ``wave``, ``batch_id``, ``n_dispatched``,
        ``failed_chunk_indices``.
    :return: Path the manifest was written to.
    """
    path = _failed_chunks_path(cov_and_an_ht_path, intervals_hash)
    payload = {
        "run_id": run_id,
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


def _new_run_id() -> str:
    """
    Return a UTC-timestamp identifier for one orchestrator run.

    Distinguishes the artifacts of successive attempts at the same chunk layout
    (which share one manifest path); seconds resolution is enough because a run
    submits at least one Hail Batch and so can never repeat within a second.

    :return: ``run-YYYYmmddTHHMMSSZ``.
    """
    return "run-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def _group_path(
    cov_and_an_ht_path: str,
    level: int,
    group_idx: int,
    partitions_per_chunk: int,
    merge_group_size: int,
) -> str:
    """
    Return a per-group merged HT path under ``<cov_and_an_path>_merge_groups_*/``.

    Level-tagged so a recursive merge tree (multiple levels of grouping when chunk
    count exceeds ``merge_group_size`` once-over) doesn't overwrite earlier-level
    outputs. The merge-groups directory also encodes the tree-shape params
    (``partitions_per_chunk`` + ``merge_group_size``): a rerun with a different
    value writes to a fresh directory rather than silently reusing / mixing stale
    group HTs from an incompatible tree shape.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param level: Merge-tree level (1-indexed); level N output is the
        input to level N+1.
    :param group_idx: Group index within this level (zero-based).
    :param partitions_per_chunk: Partitions per chunk (tree-base shape).
    :param merge_group_size: Chunk HTs per group-merge job (tree fan-in).
    :return: Per-group HT path.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    tree = f"pp{partitions_per_chunk}_gs{merge_group_size}"
    return f"{base}_merge_groups_{tree}/L{level:02d}_{group_idx:08d}.ht"


def _apply_path_suffix(path: str, suffix: str | None) -> str:
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


def _results_environment(
    environment: str, test: bool, override: str | None = None
) -> str:
    """
    Return the bucket-selecting environment for coverage release artifacts.

    Release artifacts go to the gnomAD (``dataproc``) bucket regardless of
    compute environment, so an AoU run in ``batch`` writes to the gnomAD bucket.
    A test run keeps results in the compute environment's tmp bucket, unless
    ``override`` (``--results-environment``) forces a specific bucket -- e.g. to
    smoke-test the batch-compute -> dataproc-bucket write before a prod run.

    :param environment: Compute environment.
    :param test: Whether this is a test run.
    :param override: Explicit results environment (``--results-environment``);
        takes precedence over the test/prod default when set.
    :return: ``override`` if set; else ``"dataproc"`` for prod, ``environment`` for test.
    """
    if override is not None:
        return override
    return environment if test else "dataproc"


def _resolve_cov_and_an_ht_path(
    project: str,
    environment: str,
    test: bool,
    suffix: str | None,
    chrom: str | None = None,
) -> str:
    """
    Return the cov_and_an HT path, applying ``suffix`` and ``chrom`` when set.

    A ``--chrom`` run appends the contig to the suffix so every derived path
    (per-chunk, merge-group, final HT, and the validate read) is contig-scoped and a
    per-contig run can go all the way through fan-out -> merge -> validate without
    colliding with or clobbering another contig's outputs. Folded here (not by
    mutating ``args``) so the orchestrator and the relay workers it spawns -- which
    each re-resolve from the same ``suffix`` + ``chrom`` -- always agree, with no
    risk of double-folding. Assemble the per-contig HTs with
    ``--assemble-chrom-coverage``.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param test: Whether to route to the test path.
    :param suffix: Optional suffix appended before the ``.ht`` extension
        (for A/B comparison runs).
    :param chrom: Optional single contig (``--chrom``); appended to ``suffix`` so the
        path is contig-scoped.
    :return: Fully resolved cov_and_an HT path.
    """
    path = coverage_and_an_path(
        test=test, data_set=project, environment=environment
    ).path
    full_suffix = f"{suffix}_{chrom}" if suffix and chrom else (suffix or chrom)
    return _apply_path_suffix(path, full_suffix)


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


def _gnomad_v5_merged_path(environment: str, coverage_type: str, test: bool) -> str:
    """
    Return the gnomAD v5 merged (consent-drop-subtracted) coverage/AN HT path.

    Written by ``--merge-gnomad-coverage`` / ``--merge-gnomad-an`` and read back
    by the ``--export-*`` join, so the write and read sides must resolve to the
    same path. ``test`` routes to a separate ``_test`` path so a test run never
    clobbers (or reads) the prod intermediate.

    :param environment: Compute environment.
    :param coverage_type: "coverage" or "allele_number".
    :param test: If True, append ``_test`` to keep test isolated from prod.
    :return: GCS path to the gnomAD v5 merged HT.
    """
    name = "coverage" if coverage_type == "coverage" else "an"
    return (
        f"{qc_temp_prefix(environment=environment)}"
        f"gnomad_v5_genomes_{name}{'_test' if test else ''}.ht"
    )


def _log_name_for_run(
    run_chunk: bool,
    run_merge: bool,
    chunk_start: int,
    chunk_stop: int,
    merge_output: str | None,
) -> str:
    """
    Build a per-worker log name so concurrent chunk/merge workers don't clobber the monolithic log.

    :param run_chunk: ``--run-chunk`` flag.
    :param run_merge: ``--run-merge`` flag.
    :param chunk_start: Chunk start index (used when ``run_chunk`` is True).
    :param chunk_stop: Chunk stop index, exclusive (used when ``run_chunk`` is True).
    :param merge_output: Merge output HT path (used to derive a label when ``run_merge`` is True).
    :return: Log name; suitable for use as a log-file basename.
    """
    if run_chunk:
        return f"v5_cov_chunk_{chunk_start:08d}_{chunk_stop:08d}"
    if run_merge:
        merge_label = "merge"
        if merge_output:
            merge_label = merge_output.rstrip("/").split("/")[-1].removesuffix(".ht")
        return f"v5_cov_{merge_label}"
    return "v5_coverage_and_an_generation"


def _derive_chunk_locus_intervals(
    vds_filtered: hl.vds.VariantDataset,
    n_subdivisions: int = 1,
    reference_genome: str = "GRCh38",
) -> list[hl.utils.Interval]:
    """
    Derive per-contig locus sub-intervals covering the filtered VDS reference_data.

    The vep_context HT (used as ``ref_ht``) and the AoU VDS have independent
    partition layouts, so chunking each by its own partition index would not
    yield the same loci. Instead, we derive the loci within the VDS chunk and
    use them to align ``ref_ht`` to those same loci.

    Returns one concrete ``hl.Interval`` per contig (with ``n_subdivisions=1``,
    the default and the only value any caller passes). Used to scope sites to the
    loci within the probed VDS partitions: the strict single-job compute applies
    it via ``hl.filter_intervals`` to scope ``ref_ht`` (``_build_chunk_ref_ht``),
    and ``--write-vep-context-sites --test`` feeds it to
    ``_build_vep_context_sites_ht``, which reads vep_context with ``_intervals=``.
    ``n_subdivisions > 1`` would split each contig's ``lo..hi`` range into that
    many equal-position sub-intervals, but no caller does so.

    :param vds_filtered: VDS already filtered to the chunk's partitions
        (the loci within the chunk are derived from its ``reference_data``).
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
    sub_intervals: list[hl.utils.Interval] = []
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
    n_partitions: int, chrom: str | None = None
) -> list[hl.utils.Interval]:
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
    :param chrom: Optional single contig to restrict vep_context to before
        deriving partitions, so the intervals match a ``--chrom``-filtered run.
    :return: List of bare-``Locus`` ``hl.Interval`` objects.
    """
    ref_ht = vep_context.versions["105"].ht().key_by("locus")
    if chrom:
        ref_ht = hl.filter_intervals(ref_ht, [hl.parse_locus_interval(chrom)])
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
    interval_ht: hl.Table | None = None,
    group_membership_ht: hl.Table | None = None,
    reduce_min_aggs: bool = False,
    experimental_densify: bool = False,
) -> hl.Table:
    """
    Compute coverage, allele number, and quality histograms per reference site.

    .. note::

        Running this function prior to calculating frequencies removes the need for an additional
        densify for frequency calculations.

    :param vds: Input VDS. Reference data must carry ``END`` (and ``GQ`` for
        adj/qual hists); a ``LEN`` field is added if missing. DP is read from
        variant data only -- missing DP on reference blocks is treated as 0.
    :param ref_ht: Locus-only sites Table (typically derived from
        ``vep_context`` and stripped of telomeres/centromeres/chrM) defining
        the reference positions at which to aggregate.
    :param sex_karyotype_field: Dotted path on the variant_data column
        struct to the sample's sex karyotype (e.g.
        ``"meta.sex_karyotype"``). Used by ``compute_stats_per_ref_site``
        to set per-sample ploidy on sex chromosomes.
    :param project: "aou" or "gnomad". When "aou", additionally fans out
        per-site qual histograms (``qual_hists``) and computes no coverage
        stats (see :return:); "gnomad" computes only coverage and AN.
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
    :param experimental_densify: If True, merge the VDS into one sparse MT
        (``to_merged_sparse_mt``) and unify the genotype into ``LGT`` so
        ``compute_stats_per_ref_site`` takes the  ``hl.experimental.densify`` path
        instead of ``to_dense_mt``. Default is False.
    :return: HT keyed by locus with per-stratum ``AN``; for
        ``project == "gnomad"``, flat
        ``mean``/``over_X``/``median_approx``/``total_DP`` coverage fields
        (from the global adj-filtered group); for ``project == "aou"``,
        ``qual_hists`` and no coverage fields. AoU reference blocks carry
        GQ/END/GT/LEN and no DP, so DP is defined only where a sample
        carries a variant and coverage would aggregate to ~0 genome-wide
        (measured on chr1:55.06-55.16Mb: mean 0.048x, ``over_1`` 432 of
        365,318 samples), so it is never computed for AoU.
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

    entry_agg_funcs = {"AN": get_allele_number_agg_func("LGT")}
    # Coverage stats are gnomAD-only: AoU coverage would measure ~0 (no DP on
    # reference blocks -- see the docstring), so AoU skips the 101-bin
    # ``hl.agg.counter`` plus ``hl.agg.approx_median`` sketch per reference site.
    if project == "gnomad":
        entry_agg_funcs["coverage_stats"] = get_coverage_agg_func(
            dp_field="DP", max_cov_bin=max_cov_bin
        )
    # Restrict ``coverage_stats`` to the global adj-filtered group only —
    # downstream code uses ``coverage_stats[0]`` exclusively (the global adj
    # group; transmuted into flat mean/sum/over_X fields), so per-strata
    # coverage aggregation is wasted work. We use ``{"group": "adj"}`` (and
    # not ``"raw"``) because ``compute_stats_per_ref_site`` does NOT rewrite
    # group labels to "raw" when a pre-built ``group_membership_ht`` is
    # passed (the rewrite only happens on the ``strata_expr`` code path), so
    # ``freq_meta[0]`` is ``{"group": "adj"}`` and we need to match it.
    # AN is intentionally omitted so it still fans out across all strata,
    # which is what downstream consumers need.
    # NOTE: with no AoU coverage stats, nothing on the AoU path needs
    # ``{"group": "adj"}`` to be a directly-computed leaf anymore, but that is
    # fixed in the group_membership HT (``PINNED_LEAF_GROUPS``) at write time, so
    # unpinning it is a separate change requiring a group_membership regen -- an
    # additional saving not measured here.
    entry_agg_group_membership = {}
    if project == "gnomad":
        entry_agg_group_membership["coverage_stats"] = [{"group": "adj"}]
    # Only compute qual hists for AoU.
    if project == "aou":
        entry_agg_funcs["qual_hists"] = (lambda t: [t.GQ, t.DP, t.adj], _get_hists)

        # Below we use just the raw group for qual hist computations because qual hists
        # has its own built-in adj filtering when adj is passed as an argument and will
        # produce both adj and raw histograms. We keep only the adj-filtered histograms
        # (see _get_hists -- raw qual hists were approved for removal from v5).
        entry_agg_group_membership["qual_hists"] = [{"group": "raw"}]

    logger.info(
        "Computing coverage, allele number, and optionally qual hists per reference site..."
    )

    vmt = vds.variant_data
    sex_expr = reduce(lambda x, field: x[field], sex_karyotype_field.split("."), vmt)
    vmt = vmt.annotate_cols(sex_karyotype=sex_expr)
    rmt = vds.reference_data
    # Computing LEN is only needed for VDS written before Hail 0.2.134, which
    # added LEN to the reference data; VDS written using Hail 0.2.134 or greater have it
    # already, so this branch can be dropped once the inputs are all >= 0.2.134.
    # https://hail.is/docs/0.2/change_log.html#version-0-2-134
    if "LEN" not in rmt.entry:
        rmt = rmt.annotate_entries(LEN=rmt.END - rmt.locus.position + 1)
    vds = hl.vds.VariantDataset(rmt, vmt)

    # Merge the VDS into one sparse MT so compute_stats_per_ref_site takes the lean
    # ``hl.experimental.densify`` scan (the is_vds=False branch) instead of
    # ``to_dense_mt``'s variant+reference outer-join + per-cell ``coalesce_join``. The
    # genotype is unified per cell (ref blocks carry GT, variant sites LGT); AN sums
    # ploidy, which is local/global- and multiallelic-invariant, and coalesce mirrors
    # what ``to_dense_mt`` does internally, so outputs match. Passing a missing ref
    # allele makes ``to_merged_sparse_mt`` skip its ``has_sequence``/``sequence_context``
    # branch, so no reference-genome FASTA sequence is needed.
    mtds = vds
    if experimental_densify:
        mtds = hl.vds.to_merged_sparse_mt(
            vds, ref_allele_function=lambda locus: hl.missing("str")
        )
        # Unify the genotype into LGT: ref blocks carry GT, variant sites carry LGT.
        # compute_stats_per_ref_site resolves gt_field = GT if present else LGT, and the
        # AN aggregator reads LGT.ploidy, so we expose LGT alone (drop GT) to match the
        # VDS/variant_data schema (where only LGT exists). Ploidy is identical for GT and
        # LGT, so AN is unchanged and ref sites are not undercounted.
        mtds = mtds.annotate_entries(LGT=hl.coalesce(mtds.LGT, mtds.GT))
        mtds = mtds.select_entries(
            *[f for f in ("LGT", "GQ", "DP", "adj", "END") if f in mtds.entry]
        )

    ht = compute_stats_per_ref_site(
        mtds,
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

    # Keep coverage stats from global adj grouping (index 0) only. The sample-count
    # global is annotated either way -- it is just ``strata_sample_count[0]`` (the
    # full cohort) and downstream readers expect it.
    ht = ht.annotate_globals(
        coverage_stats_meta_sample_count=ht.strata_sample_count[0],
    )
    if project == "gnomad":
        cov_stats_expr = _cov_stats(ht.coverage_stats[0])
        ht = ht.transmute(**cov_stats_expr)

    if project == "aou":
        # ``qual_hists`` as returned by ``compute_stats_per_ref_site`` is an array of length 1 so we drop the array here.
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

    Function transforms mean back into sum using ``sample_count`` argument.

    :param ht: Input HT.
    :param project: Project name.
    :param sample_count: Number of samples in HT.
    :param coverage_over_x_bins: Boundaries for the samples-over-X fields.
        Default is :data:`COVERAGE_OVER_X_BINS`.
    :return: Renamed HT.
    """
    # Transform mean back into sum.
    ht = ht.transmute(sum=ht.mean * sample_count)

    # Rename annotations to include project.
    row_fields = list(ht.row_value)
    rename_dict = {f: f"{f}_{project}" for f in row_fields}
    ht = ht.rename(rename_dict)
    if project == "gnomad_release":
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
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
) -> hl.expr.DictExpression:
    """
    Subtract ``project_2``'s coverage fields from ``project_1``'s.

    Only the consent-drop subtraction remains as a caller: AoU computes no
    coverage stats, so there is no AoU + gnomAD coverage sum any more.

    .. note::
        - Function does not merge ``median_approx`` fields.

    :param ht: Input HT. Must have annotations from both projects.
    :param project_1: Project subtracted FROM (the minuend).
    :param project_2: Project being subtracted (the subtrahend).
    :param coverage_over_x_bins: Boundaries for the samples-over-X fields.
        Default is :data:`COVERAGE_OVER_X_BINS`.
    :return: Merged fields.
    """
    # Make sure over_x fields are float64s.
    merged_fields = {
        "sum_gnomad": ht[f"sum_{project_1}"] - ht[f"sum_{project_2}"],
        "total_DP_gnomad": (ht[f"total_DP_{project_1}"] - ht[f"total_DP_{project_2}"]),
    }
    merged_fields.update(
        {
            f"over_{x}_gnomad": (
                ht[f"over_{x}_{project_1}"] - ht[f"over_{x}_{project_2}"]
            )
            for x in coverage_over_x_bins
        }
    )
    return merged_fields


def merge_gnomad_coverage_hts(
    gnomad_ht: hl.Table,
    gnomad_release_ht: hl.Table,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
    gnomad_sample_count: int = GNOMAD_SAMPLE_COUNT,
    consent_drop_count: int = GNOMAD_CONSENT_DROP_SAMPLE_COUNT,
) -> hl.Table:
    """
    Subtract consent drop samples from gnomAD v4 genomes release HT to create gnomAD v5 genomes coverage HT.

    :param gnomad_ht: gnomAD coverage HT (contains coverage for consent drop samples only).
    :param gnomad_release_ht: gnomAD v4 genomes coverage release HT.
    :param coverage_over_x_bins: Boundaries for the samples-over-X fields.
        Default is :data:`COVERAGE_OVER_X_BINS`.
    :param gnomad_sample_count: Number of release gnomAD v4 genome samples. Default is `GNOMAD_SAMPLE_COUNT`.
    :param consent_drop_count: Number of consent drop gnomAD v4 genome samples. Default is `GNOMAD_CONSENT_DROP_SAMPLE_COUNT`.
    :return: gnomAD v5 genomes coverage HT.
    """
    logger.info(
        "Subtracting gnomAD v4 consent drop samples from gnomAD v4 genomes release HT..."
    )
    gnomad_ht = _rename_cov_annotations(
        gnomad_ht, "gnomad", consent_drop_count, coverage_over_x_bins
    )
    gnomad_release_ht = _rename_cov_annotations(
        gnomad_release_ht, "gnomad_release", gnomad_sample_count, coverage_over_x_bins
    )
    gnomad_v5_count = gnomad_sample_count - consent_drop_count
    logger.info("Total number of gnomAD v5 release genomes: %s", gnomad_v5_count)

    gnomad_ht = gnomad_ht.join(gnomad_release_ht, "right")
    merged_fields = _merge_coverage_fields(
        ht=gnomad_ht,
        project_1="gnomad_release",
        project_2="gnomad",
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


def _rename_fields(
    ht: hl.Table, field_name: str, project: str, rename_globals: bool
) -> hl.Table:
    """
    Rename fields by adding project name prior to merging Tables.

    Used for AN and qual hists, not coverage: coverage goes through
    ``_rename_cov_annotations`` instead (it needs the v4 mean->sum transform) and
    doesn't need the ``strata_meta``/``strata_sample_count`` global renaming this
    function does.

    :param ht: Input HT.
    :param field_name: Row field to suffix with ``project`` (e.g. ``"AN"``).
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
) -> tuple[
    hl.expr.ArrayExpression,
    list[dict[str, str]],
    dict[str, hl.expr.ArrayExpression],
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

    The combined (AoU + gnomAD) per-stratum AN for the strata both projects share is
    written first as the unlabeled "global" entries (the AoU-only downsampling strata
    are NOT global). The AoU-only AN is then appended as an ``"aou"`` subset (each
    appended ``strata_meta`` entry carries ``{"subset": "aou"}``), including the
    AoU-only downsampling strata -- which thus need no ``"aou-downsampling"`` rename.

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
    # Global = the AoU + gnomAD combined AN for the strata both projects share (the
    # standard adj/raw x gen_anc x sex strata). The AoU-only downsampling strata are
    # NOT global -- they are exposed only in the "aou" subset below -- so drop them
    # here (joint_strata_meta is a Python list; keep the matching AN / count entries).
    global_idx = [i for i, d in enumerate(joint_strata_meta) if "downsampling" not in d]
    global_meta = [joint_strata_meta[i] for i in global_idx]
    ht = ht.annotate(AN=hl.array([joint_an[i] for i in global_idx]))
    ht = ht.annotate_globals(
        strata_meta=global_meta,
        strata_sample_count=hl.array(
            [count_arrays_dict["counts"][i] for i in global_idx]
        ),
    )

    # Append the AoU-only AN as the "aou" subset. The global entries set above stay
    # unlabeled; each appended entry gets {"subset": "aou"} added to its strata_meta,
    # including the AoU-only downsampling strata (so they carry the subset tag and
    # need no "aou-downsampling" key rename). The AoU side (row AN_aou and the
    # strata_meta_aou / strata_sample_count_aou globals) survives the join, and AN_aou
    # is always defined (AoU is the left side of the join), so we extend the three
    # parallel arrays in lockstep with no fill.
    subset_meta = ht.index_globals().strata_meta_aou.map(
        lambda d: hl.dict(d.items().append(("subset", "aou")))
    )
    ht = ht.annotate(AN=ht.AN.extend(ht.AN_aou))
    ht = ht.annotate_globals(
        strata_meta=ht.strata_meta.extend(subset_meta),
        strata_sample_count=ht.strata_sample_count.extend(
            ht.index_globals().strata_sample_count_aou
        ),
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
    # time; see ``_get_hists``). The reused gnomAD v4 hists still carry raw_qual_hists,
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
    Filter ``target_ht`` to per-contig locus bounds derived from ``source_ht``.

    Used in test runs where ``target_ht`` (e.g., a public release HT) and
    ``source_ht`` (e.g., a freshly-computed test HT) do not share partition
    layouts, so their first-N-partition slices do not overlap. Filtering to
    the per-contig min..max locus range of ``source_ht`` guarantees test
    join/merge overlap at low cost (one small aggregation + partition-pruned
    interval filter).

    :param target_ht: Table to filter.
    :param source_ht: Table whose locus ranges define the filter intervals.
    :return: ``target_ht`` filtered to the per-contig locus ranges of
        ``source_ht``.
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


def _aou_sample_artifact_path(name: str, environment: str, test: bool = False) -> str:
    """
    Return the path to a chunk-invariant AoU sample-level artifact (30-day storage).

    These are small files written once per run so each chunk worker reads them
    instead of rescanning the ~330-partition sample tables: ``collisions.json``
    (sample IDs colliding with gnomAD) and ``release_samples.json`` (AoU release
    sample IDs).

    NOTE: this path scheme is shared with ``generate_frequency.py``, which writes
    these files under ``--write-sample-artifacts`` (and also writes ``meta_small.ht``
    and ``gm_coalesced[_reduce].ht``, which are frequency-specific). The two are
    deliberately byte-identical schemes so one writer serves both pipelines; when
    those branches merge, hoist this helper into ``basics.py`` and have both import
    it rather than keeping two copies in sync by hand.

    :param name: Artifact file name (with extension).
    :param environment: Compute environment.
    :param test: If True, return the test-scoped path.
    :return: GCS path to the artifact.
    """
    scope = "test" if test else "full"
    return (
        f"{qc_temp_prefix(environment=environment, days=30)}"
        f"aou_sample_artifacts_{scope}/{name}"
    )


def _read_sample_artifact_json(
    name: str, environment: str, test: bool = False
) -> Any | None:
    """
    Return the parsed JSON sample artifact, or None when it has not been precomputed.

    Absence is not an error: every caller falls back to computing the value itself,
    so a run without precomputed artifacts is correct, only more expensive.

    :param name: Artifact file name (with extension).
    :param environment: Compute environment.
    :param test: If True, read the test-scoped path.
    :return: Parsed JSON value, or None if the file is absent.
    """
    path = _aou_sample_artifact_path(name, environment, test)
    if not file_exists(path):
        return None
    with hfs.open(path) as f:
        return json.load(f)


def _load_project_vds(
    project: str,
    environment: str,
    partition_range: list[int] | None = None,
    sub_intervals: list[hl.utils.Interval] | None = None,
    filter_intervals: list[hl.utils.Interval] | None = None,
    chrom: str | None = None,
    test: bool = False,
    test_sample_subset: bool = False,
) -> tuple[hl.vds.VariantDataset, str]:
    """
    Load the per-project VDS with consistent test/subsample handling.

    Centralizes the project-conditional VDS load that's needed by both the
    chunk worker (``_run_coverage_chunk``) and the strict path (``main``):

    - filter_partitions vs read_intervals routing (both AoU and gnomAD):
      ``sub_intervals`` takes precedence; ``partition_range`` is the
      fallback. When ``sub_intervals`` is provided the VDS is read with
      ``read_intervals`` so it's co-partitioned with the vep_context read
      (callers may widen the VDS read's leading edge for straddling
      reference blocks, so the endpoints aren't necessarily identical).
    - AoU-specific DP synthesis from LAD + ``annotate_adj_no_dp`` (the AoU
      v8 VDS lacks DP, which ``compute_stats_per_ref_site`` requires).
    - ``test_sample_subset`` application: AoU only, reads ``meta`` from
      disk inside the helper so it's self-contained.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param partition_range: VDS partition indices (e.g. ``list(range(3))``)
        or ``None`` for the full VDS.
    :param sub_intervals: Locus sub-intervals for read-time partitioning
        (applies to both AoU and gnomAD).
    :param filter_intervals: Locus intervals passed to ``get_*_vds`` as
        ``filter_intervals`` (``split_reference_blocks=True``): subsets the VDS to
        these intervals and SPLITS reference blocks straddling the edges, so the
        in-region piece keeps its coverage with a bounded read (no backing the read
        up by ``ref_block_max_length``). Used by ``--test-region``.
    :param chrom: Optional single contig to filter to.
    :param test: Whether this is a test run (gates ``test_sample_subset``).
    :param test_sample_subset: If True (AoU only, and ``test``), subsample
        ~0.1% of the AoU sample set for cheap tiny-cohort runs.
    :return: Tuple of ``(vds, sex_karyotype_field)``.
    """
    if project == "aou":
        sex_karyotype_field = "meta.sex_karyotype"
        # Each get_aou_vds call otherwise pays a full aggregation of the
        # ~329-partition sample-collisions HT (to fetch 27 IDs) plus a collect of the
        # ~330-partition meta HT for release_only (365k IDs) -- per chunk process, so
        # the fan-out pays them thousands of times over. Read the precomputed JSON
        # instead when it exists; fall back to the scans when it doesn't.
        collisions = _read_sample_artifact_json("collisions.json", environment, test)
        release_samples = _read_sample_artifact_json(
            "release_samples.json", environment, test
        )
        vds = get_aou_vds(
            release_only=release_samples is None,
            filter_samples=release_samples,
            # The release list holds POST-prefix IDs, and prefixing only runs when one
            # of release_only/annotate_meta/add_project_prefix is set -- annotate_meta
            # covers it here, but set it explicitly so the filter can't silently miss
            # the prefixed samples if annotate_meta ever goes away.
            add_project_prefix=release_samples is not None,
            sample_collisions=set(collisions) if collisions is not None else None,
            filter_partitions=(
                None if (sub_intervals or filter_intervals) else partition_range
            ),
            read_intervals=sub_intervals,
            filter_intervals=filter_intervals,
            annotate_meta=True,
            chrom=chrom,
            environment=environment,
        )
        vmt = vds.variant_data
        vmt = annotate_adj_no_dp(vmt)
        vmt = vmt.annotate_entries(DP=hl.sum(vmt.LAD))
        rmt = vds.reference_data
        # Reference blocks are hom-ref with no AD/LAD; gq_only computes adj from GQ
        # alone (the AB test is vacuous for non-het calls, so this equals get_adj_expr
        # for ref blocks). Without adj on ref data, ref-block samples weren't marked
        # adj, undercounting adj AN/coverage at reference sites.
        rmt = annotate_adj_no_dp(rmt, gq_only=True)
        vds = hl.vds.VariantDataset(rmt, vmt)

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
            filter_partitions=(
                None if (sub_intervals or filter_intervals) else partition_range
            ),
            read_intervals=sub_intervals,
            filter_intervals=filter_intervals,
            annotate_meta=True,
            chrom=chrom,
        )
    return vds, sex_karyotype_field


def _probe_vds(
    project: str,
    environment: str,
    partition_range: list[int] | None,
    chrom: str | None,
    filter_intervals: list[hl.utils.Interval] | None = None,
) -> hl.vds.VariantDataset:
    """
    Cheap reference_data-bounds probe-load of the per-project VDS.

    Loads via ``filter_partitions`` (or ``filter_intervals``) only (no
    release/subset filtering or DP synthesis, though AoU's base sample exclusion
    still applies); the caller derives the loci from ``reference_data``
    (e.g. via ``repartition_for_join``). Used by the chunk worker (both the normal
    partition-range path and the ``--test-region`` fan-out) and the strict path's
    ``--test-n-partitions`` co-partitioning branch.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param partition_range: VDS partition indices to probe (e.g.
        ``list(range(2))``) or ``None`` for the full VDS.
    :param chrom: Optional single contig to filter to.
    :param filter_intervals: Optional locus intervals to bound the probe to (via
        ``hl.vds.filter_intervals``, splitting straddling reference blocks). Used by
        the ``--test-region`` fan-out to probe just the region's reference-data loci.
    :return: Probe VDS (reference/variant data only).
    """
    if project == "aou":
        return get_aou_vds(
            filter_partitions=partition_range,
            chrom=chrom,
            environment=environment,
            remove_hard_filtered_samples=False,
            log_sample_counts=False,
            filter_intervals=filter_intervals,
        )
    return get_gnomad_v5_genomes_vds(
        filter_partitions=partition_range,
        chrom=chrom,
        filter_intervals=filter_intervals,
    )


def _vep_context_sites_path(test: bool = False) -> str:
    """
    Return the path to the preprocessed vep_context sites HT (30-day storage).

    Prod is genome-wide; a ``--test`` write is scoped to the loci within the test
    VDS partitions and written to a separate ``_test`` path so it never clobbers
    the genome-wide prod table.

    :param test: If True, return the test-scoped sites path.
    :return: GCS path to the preprocessed sites HT.
    """
    # Pin to dataproc so it is accessible everywhere
    return f"{qc_temp_prefix(environment='dataproc', days=30)}{'vep_context_sites_test.ht' if test else 'vep_context_sites.ht'}"


def _build_vep_context_sites_ht(
    intervals: list[hl.utils.Interval] | None = None,
) -> hl.Table:
    """
    Build the locus-keyed, deduped, telomere/centromere/chrM-stripped vep_context sites HT.

    This is the per-ref-site set the coverage/AN compute aggregates at. Previously
    this dedup + strip ran inside every chunk (``.distinct()`` shuffle + a driver
    ``collect()`` of the telomere intervals, ~48k times); it is now done once via
    ``--write-vep-context-sites`` and read co-partitioned per chunk.

    :param intervals: Optional locus intervals to scope the table to (read-time
        pruned); used for a cheap ``--test`` write scoped to the loci within the
        test partitions. None preprocesses the whole genome.
    :return: Locus-keyed sites HT.
    """
    if intervals is not None:
        ref_ht = vep_context.versions["105"].ht(read_args={"_intervals": intervals})
    else:
        ref_ht = vep_context.versions["105"].ht()
    ref_ht = ref_ht.key_by("locus").select().distinct()
    # Drop telomeres/centromeres and EXCLUDED_CONTIGS (chrM) in one pass. Dropped here
    # rather than downstream so this table stays the single definition of "every site
    # the compute must cover" and --validate-cov-and-an remains an exact-equality
    # check. Concrete intervals (hl.eval) so the list stays homogeneous with the
    # collected telomere/centromere intervals.
    drop_intervals = telomeres_and_centromeres.ht().interval.collect()
    drop_intervals += [
        hl.eval(hl.parse_locus_interval(c, reference_genome="GRCh38"))
        for c in EXCLUDED_CONTIGS
    ]
    ref_ht = hl.filter_intervals(ref_ht, drop_intervals, keep=False)
    return ref_ht


def _chunk_intervals_path(environment: str, test: bool = False) -> str:
    """
    Return the path to the precomputed per-chunk read sub-intervals JSON (30-day storage).

    Written once by ``--write-chunk-intervals``: a single VDS open derives the balanced
    read sub-intervals for *every* chunk, replacing the per-chunk probe (a redundant
    ~400s VDS open per chunk). Stored as a plain JSON file (not a Hail Table) so the
    orchestrator/merge (to enumerate and --chrom-filter chunks) and the worker (to look
    up its own chunk) read it driver-side with ``hailtop.fs`` -- no QoB job, so no
    query-on-batch cold-start. The JSON holds ``chunks`` (a list indexed by chunk index;
    each entry ``{"contig", "intervals"}``, never spanning a contig boundary), the VDS
    ``ref_block_max_length``, the reference-genome name (for the leading-edge boundary
    fix, see :func:``_expand_leading_edges``), and ``read_subintervals_per_chunk``.

    :param environment: Compute environment.
    :param test: If True, return the test-scoped intervals path.
    :return: GCS path to the precomputed chunk-intervals JSON.
    """
    name = "chunk_intervals_test.json" if test else "chunk_intervals.json"
    return f"{qc_temp_prefix(environment=environment, days=30)}{name}"


def _interval_to_list(iv: hl.utils.Interval) -> list[str | int | bool]:
    """
    Serialize a locus interval to a JSON-friendly list.

    :param iv: Locus interval (Python ``hl.Interval`` with ``hl.Locus`` endpoints).
    :return: ``[start_contig, start_pos, end_contig, end_pos, includes_start,
        includes_end]``.
    """
    return [
        iv.start.contig,
        iv.start.position,
        iv.end.contig,
        iv.end.position,
        iv.includes_start,
        iv.includes_end,
    ]


def _interval_from_list(
    t: list[str | int | bool], reference_genome: str
) -> hl.utils.Interval:
    """
    Reconstruct a locus interval from its :func:``_interval_to_list`` serialization.

    :param t: ``[start_contig, start_pos, end_contig, end_pos, includes_start,
        includes_end]``.
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


def _parse_region_interval(
    s: str, reference_genome: str = "GRCh38"
) -> hl.utils.Interval:
    """
    Parse a ``contig:start-end`` string into a half-open Python locus interval.

    Returns a concrete ``hl.utils.Interval`` (not a Hail ``IntervalExpression``,
    which ``hl.parse_locus_interval`` would give) so it can be passed to
    ``read_args={"_intervals": ...}`` for read-time partition pruning -- that path
    requires Python-bool ``includes_start`` / ``includes_end``. Half-open ``[start,
    end)`` so adjacent intervals stay disjoint (no double-counted loci). Used by
    ``--test-region``.

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


def _split_intervals_at_contigs(
    intervals: list[hl.utils.Interval], reference_genome: str
) -> list[hl.utils.Interval]:
    """
    Split any locus interval that straddles a contig boundary into one per contig.

    ``repartition_for_join`` returns a contiguous partition of the keyspace, so the
    interval covering a contig transition is ``[chrA:x, chrB:y)`` -- it spans the tail
    of chrA and the head of chrB. Contig-keyed chunking requires every interval (hence
    every chunk) to belong to exactly one contig, so such intervals are cut at contig
    lengths; single-contig intervals pass through unchanged. The pieces tile the
    original interval exactly (no gap or overlap), so coverage is preserved.

    :param intervals: Locus intervals (e.g. from ``repartition_for_join``), sorted.
    :param reference_genome: Reference-genome name (e.g. "GRCh38").
    :return: Intervals, each with ``start.contig == end.contig``.
    """
    rg = hl.get_reference(reference_genome)
    contigs = rg.contigs
    lengths = rg.lengths
    out: list[hl.utils.Interval] = []
    for iv in intervals:
        if iv.start.contig == iv.end.contig:
            out.append(iv)
            continue
        si = contigs.index(iv.start.contig)
        ei = contigs.index(iv.end.contig)
        for k in range(si, ei + 1):
            contig = contigs[k]
            first, last = k == si, k == ei
            lo = iv.start.position if first else 1
            hi = iv.end.position if last else lengths[contig]
            inc_s = iv.includes_start if first else True
            inc_e = iv.includes_end if last else True
            # Drop an empty tail piece (e.g. end at chrB:1 exclusive => nothing on
            # chrB).
            if hi < lo or (hi == lo and not (inc_s and inc_e)):
                continue
            out.append(
                hl.Interval(
                    hl.Locus(contig, lo, reference_genome=reference_genome),
                    hl.Locus(contig, hi, reference_genome=reference_genome),
                    includes_start=inc_s,
                    includes_end=inc_e,
                )
            )
    return out


def _build_chunk_intervals(
    project: str,
    environment: str,
    total_partitions: int,
    partitions_per_chunk: int,
    n_sub: int,
    chrom: str | None = None,
) -> dict[str, Any]:
    """
    Precompute every chunk's balanced read sub-intervals in ONE VDS open, by contig.

    Replaces the per-chunk probe (``_probe_vds`` + ``repartition_for_join`` on each
    chunk). Opens the VDS once over the leading ``total_partitions``, derives
    ``n_chunks * n_sub`` balanced sub-intervals over all reference-data loci (the
    pre-split count; the contig split below can add a few), then splits any that
    straddle a contig boundary (:func:`_split_intervals_at_contigs`), groups them by
    contig, and slices each contig's run into chunks of ``n_sub`` consecutive
    sub-intervals. So no chunk crosses a contig boundary, and ``--chrom`` can select a
    single contig's chunks; chunks remain disjoint and their union covers all loci.

    Returns a JSON-serializable dict (intervals as plain lists), read driver-side by the
    orchestrator (to enumerate/filter chunks by contig) and the worker (to look up its
    own chunk by index) without a QoB job.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param total_partitions: Number of leading VDS partitions to derive intervals over
        (the probe extent); with ``partitions_per_chunk`` and ``n_sub`` it sets the
        sub-interval granularity (``n_chunks * n_sub`` balanced sub-intervals). The
        caller passes ``args.test_n_partitions or args.total_partitions``, so for a
        test ``--test-n-partitions`` overrides this to the first N VDS partitions
        (combine with ``chrom`` to slice a single contig).
    :param partitions_per_chunk: Sets ``n_chunks = ceil(total_partitions / this)``, the
        sub-interval-count divisor (no longer a chunk key -- chunks are contig-grouped).
    :param n_sub: Sub-intervals per chunk (``--read-subintervals-per-chunk``).
    :param chrom: Optional single contig to restrict the precompute to.
    :return: Dict with ``read_subintervals_per_chunk``, ``ref_block_max_length``,
        ``reference_genome``, and ``chunks`` (a list indexed by chunk index; each entry
        ``{"contig": str, "intervals": [serialized sub-intervals]}``).
    """
    n_chunks = (total_partitions + partitions_per_chunk - 1) // partitions_per_chunk
    vds = _probe_vds(project, environment, list(range(total_partitions)), chrom)
    rbml_field = hl.vds.VariantDataset.ref_block_max_length_field
    if rbml_field not in vds.reference_data.globals:
        raise ValueError(
            f"VDS reference_data lacks the '{rbml_field}' global, so blocks"
            " straddling chunk boundaries cannot be bounded; run"
            " hl.vds.truncate_reference_blocks / store_ref_block_max_length"
            " on the VDS first."
        )
    max_ref_block_len = hl.eval(vds.reference_data.index_globals()[rbml_field])
    all_subs = repartition_for_join(
        vds.reference_data.rows(),
        n_partitions=n_chunks * n_sub,
        locus_intervals=True,
    )
    if not all_subs:
        raise ValueError(
            "repartition_for_join returned no sub-intervals; the VDS reference_data"
            " is empty for the requested partitions/contigs."
        )
    rg = all_subs[0].start.reference_genome
    contig_subs = _split_intervals_at_contigs(all_subs, rg.name)
    # Keep the chunk set on the same contigs as the sites HT (see EXCLUDED_CONTIGS).
    # Safe to drop after the contig split, when every interval belongs to exactly one
    # contig. If the VDS has no reference data on these contigs this is a no-op.
    dropped = [iv for iv in contig_subs if iv.start.contig in EXCLUDED_CONTIGS]
    if dropped:
        logger.info(
            "Dropping %d sub-interval(s) on excluded contig(s) %s.",
            len(dropped),
            sorted({iv.start.contig for iv in dropped}),
        )
        contig_subs = [
            iv for iv in contig_subs if iv.start.contig not in EXCLUDED_CONTIGS
        ]
    chunks: list[dict[str, Any]] = []
    for contig, group in groupby(contig_subs, key=lambda iv: iv.start.contig):
        contig_ivs = list(group)
        for j in range(0, len(contig_ivs), n_sub):
            chunks.append(
                {
                    "contig": contig,
                    "intervals": [
                        _interval_to_list(iv) for iv in contig_ivs[j : j + n_sub]
                    ],
                }
            )
    # Defensive: no chunk may mix contigs (--chrom filtering and per-contig disjointness
    # both depend on it). _interval_to_list = [s_contig, s_pos, e_contig, e_pos, ...].
    assert all(
        iv[0] == iv[2] == c["contig"] for c in chunks for iv in c["intervals"]
    ), "a chunk interval crosses a contig boundary"
    logger.info(
        "Built chunk-intervals: %d chunks across %d contigs (%d sub-intervals;"
        " ref_block_max_length=%d).",
        len(chunks),
        len({c["contig"] for c in chunks}),
        len(contig_subs),
        max_ref_block_len,
    )
    return {
        "read_subintervals_per_chunk": n_sub,
        "ref_block_max_length": max_ref_block_len,
        "reference_genome": rg.name,
        "chunks": chunks,
    }


def _build_chunk_ref_ht(
    vds_filtered: hl.vds.VariantDataset,
    partition_count: int,
    chrom: str | None,
    sites_path: str,
    sub_intervals: list[hl.utils.Interval] | None = None,
) -> hl.Table:
    """
    Build the per-chunk ``ref_ht`` from the preprocessed vep_context sites HT.

    Reads ``sites_path`` (already locus-keyed, deduped, and telomere/centromere/chrM-
    stripped by ``--write-vep-context-sites``). When ``sub_intervals`` is provided,
    reads with ``_intervals=sub_intervals`` so the ref_ht is co-partitioned with the
    VDS read (the VDS read's leading edge may be widened for straddling reference
    blocks; the sites stay un-widened so they're disjoint per chunk). Otherwise (no
    read-time sub-intervals -- the partition-range / whole-VDS reads) reads the whole
    sites HT and filters to ``chrom`` (when set) or the per-contig loci within the
    VDS chunk.

    :param vds_filtered: VDS already filtered to the chunk's partitions
        (the loci within the chunk are derived from its ``reference_data``).
    :param partition_count: When > 0 and ``sub_intervals``/``chrom`` are None,
        filters ``ref_ht`` to the loci within the VDS chunk. 0 means "whole".
    :param chrom: Optional single contig to filter to (takes precedence over
        the ``partition_count``-driven filter to the chunk's loci).
    :param sites_path: Path to the preprocessed vep_context sites HT.
    :param sub_intervals: Optional read-time intervals; take precedence over both
        ``chrom`` and the filter to the chunk's loci.
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
        ref_ht = hl.filter_intervals(ref_ht, [hl.parse_locus_interval(chrom)])
    elif partition_count > 0:
        chunk_intervals = _derive_chunk_locus_intervals(vds_filtered)
        ref_ht = hl.filter_intervals(ref_ht, chunk_intervals)
    return ref_ht


def _expand_leading_edges(
    intervals: list[hl.utils.Interval], max_ref_block_len: int
) -> list[hl.utils.Interval]:
    """
    Back up each contig's first interval start by ``max_ref_block_len - 1``.

    Each chunk is a separate ``read_vds(intervals=...)`` that gates reference
    blocks by their START locus. A reference block straddling in from the
    PREVIOUS chunk (start < this chunk's first locus, END inside it) is therefore
    dropped, silently undercounting coverage/AN at the chunk's leading sites for
    every sample whose block straddles. Widening ONLY the VDS read's leading edge
    per contig pulls those straddlers back in; the cumulative densify scan then
    carries them forward to this chunk's sites. The sites (``ref_ht``) read is
    left unchanged, so the chunk still emits only its own (disjoint) sites — no
    site is computed by two chunks, so there is no double-counting. ``ref_block_
    max_length`` is a hard cap, so a block live at the leading edge starts within
    ``max_ref_block_len - 1`` bases before it; backing up exactly that catches
    every straddler and no earlier (irrelevant) block.

    :param intervals: Sorted locus intervals (one VDS read-partition each).
    :param max_ref_block_len: The VDS ``ref_block_max_length`` global value.
    :return: Intervals with each contig's leading start backed up (clamped at 1).
    """
    seen_contigs: set[str] = set()
    expanded: list[hl.utils.Interval] = []
    for interval in intervals:
        contig = interval.start.contig
        if contig in seen_contigs:
            expanded.append(interval)
            continue
        seen_contigs.add(contig)
        new_pos = max(interval.start.position - (max_ref_block_len - 1), 1)
        expanded.append(
            hl.Interval(
                hl.Locus(
                    contig, new_pos, reference_genome=interval.start.reference_genome
                ),
                interval.end,
                includes_start=interval.includes_start,
                includes_end=interval.includes_end,
            )
        )
    return expanded


def _run_coverage_chunk(args: argparse.Namespace) -> None:
    """
    Compute one chunk of the coverage/AN HT and write it to ``args.chunk_output``.

    Invoked by the per-chunk Hail Batch relay container (see
    ``_orchestrate_coverage_batch``). Hail must already be initialized
    (``main`` does this via ``_init_hail`` before dispatching here).

    Steps:
      1. Obtain this chunk's balanced locus sub-intervals: for a --test-region run, the
         passed regions; otherwise this chunk's entry (by index, ``--chunk-start``) in
         the per-chunk precompute JSON (``--write-chunk-intervals``, required), read
         driver-side. A chunk never spans a contig boundary.
      2. Re-read the VDS via ``hl.vds.read_vds(intervals=...)`` — one partition per
         sub-interval (no shuffle) — with each contig's leading edge backed up by
         ``ref_block_max_length - 1`` (``_expand_leading_edges``) so reference
         blocks straddling in from the previous chunk are not gated out.
      3. Read the preprocessed vep_context sites HT (``_build_chunk_ref_ht``) with the
         same (un-widened) ``_intervals=sub_intervals`` — co-partitioned with the VDS
         so the densify join is shuffle-free, and the chunk emits only its own disjoint
         sites.
      4. Run ``compute_all_release_stats_per_ref_site``.
      5. Write to ``args.chunk_output``.

    :param args: Parsed CLI args; reads ``project_name``, ``environment``,
        ``chunk_start``, ``chunk_stop``, ``test``, ``test_region``, ``chrom``,
        ``read_subintervals_per_chunk``, ``chunk_output``, ``results_environment``,
        ``cov_and_an_output_suffix``, ``reduce_min_aggs``, ``test_sample_subset``,
        ``experimental_densify``. Hail must already be initialized.
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
    results_environment = _results_environment(
        environment, test, args.results_environment
    )

    group_membership_ht_path = _resolve_group_membership_ht_path(
        project,
        environment,
        test=test,
        reduce_min_aggs=args.reduce_min_aggs,
    )

    sub_intervals: list[hl.utils.Interval] | None = None
    # The VDS read uses the same intervals EXCEPT each contig's leading edge is
    # backed up by ref_block_max_length-1 to capture reference blocks straddling
    # in from the previous chunk (see _expand_leading_edges). The sites/ref_ht
    # read keeps the un-widened sub_intervals so sites stay disjoint per chunk.
    vds_read_intervals: list[hl.utils.Interval] | None = None
    # For --test-region the VDS is read via filter_intervals (which splits
    # straddling reference blocks) instead of read_intervals; this holds those.
    vds_filter_intervals: list[hl.utils.Interval] | None = None
    # Content hash of the chunk-intervals layout this chunk belongs to; namespaces
    # the output path so chunks from a different (non-reproducible) layout are never
    # mixed at merge time. "test_region" for a --test-region run (no JSON).
    intervals_hash: str
    if args.test_region:
        # Explicit region (--test-region): scope the chunk to these intervals, then
        # balance them into n_sub sub-intervals using the same probe-and-
        # repartition_for_join machinery the prod precompute (_build_chunk_intervals)
        # uses, so the densify + per-site aggregation is spread across many small
        # partitions. Reading the whole region as one partition per interval pins a
        # single worker at its memory limit and gets OOM-killed on AoU's large sample
        # count (the dense per-site agg over ~365k samples does not fit one worker).
        intervals_hash = "test_region"
        region_intervals = [_parse_region_interval(r) for r in args.test_region]
        if n_sub > 1:
            # Probe the VDS bounded to the region (filter_intervals splits straddling
            # reference blocks, so the probe's reference-data loci span exactly the
            # region), balance those loci into n_sub sub-intervals, and read them via
            # read_intervals with each contig's leading edge widened -- the same
            # machinery (and the same no-data-loss guarantee) the prod fan-out uses
            # and validates, just scoped to the region instead of a partition range.
            vds_probe = _probe_vds(
                project, environment, None, chrom, filter_intervals=region_intervals
            )
            sub_intervals = repartition_for_join(
                vds_probe.reference_data.rows(),
                n_partitions=n_sub,
                locus_intervals=True,
            )
            rbml_field = hl.vds.VariantDataset.ref_block_max_length_field
            if rbml_field not in vds_probe.reference_data.globals:
                raise ValueError(
                    f"VDS reference_data lacks the '{rbml_field}' global, so blocks"
                    " straddling sub-interval boundaries cannot be bounded; run"
                    " hl.vds.truncate_reference_blocks / store_ref_block_max_length"
                    " on the VDS first."
                )
            max_ref_block_len = hl.eval(
                vds_probe.reference_data.index_globals()[rbml_field]
            )
            vds_read_intervals = _expand_leading_edges(sub_intervals, max_ref_block_len)
            logger.info(
                "Region fan-out: balanced %d --test-region interval(s) into %d read"
                " sub-intervals (read via read_intervals + leading-edge widening).",
                len(region_intervals),
                len(sub_intervals),
            )
        else:
            # n_sub == 1: no fan-out -- read the region as one partition per interval
            # via filter_intervals (splits straddling reference blocks; bounded).
            sub_intervals = region_intervals
            vds_filter_intervals = region_intervals
    else:
        # Look this chunk up by index in the per-chunk precompute JSON
        # (--write-chunk-intervals, required for the fan-out). Read it driver-side
        # (hailtop.fs, NOT a QoB job) so we skip both the ~400s VDS re-open AND the
        # query-on-batch cold-start a Hail-Table read would incur as this worker's
        # first action. --chrom is handled by the orchestrator (it only fans out the
        # selected contigs' chunks), so the worker just processes chunks[start].
        intervals_path = _chunk_intervals_path(environment, test)
        if not file_exists(intervals_path):
            raise FileNotFoundError(
                f"chunk-intervals JSON not found at {intervals_path};"
                " --write-chunk-intervals is required before the fan-out."
            )
        with hfs.open(intervals_path) as f:
            data = json.load(f)
        intervals_hash = _chunk_intervals_hash(data)
        chunk_meta = data["chunks"]
        if not 0 <= start < len(chunk_meta):
            raise ValueError(
                f"chunk index {start} is out of range [0, {len(chunk_meta)}) in"
                f" {intervals_path}; the fan-out and precompute are out of sync --"
                " regenerate --write-chunk-intervals."
            )
        entry = chunk_meta[start]
        rg = data["reference_genome"]
        sub_intervals = [_interval_from_list(t, rg) for t in entry["intervals"]]
        max_ref_block_len = data["ref_block_max_length"]
        logger.info(
            "Read %d sub-intervals for chunk %d (contig %s) from the precompute.",
            len(sub_intervals),
            start,
            entry["contig"],
        )
        # Back up the VDS read's leading edge per contig so reference blocks straddling
        # in from the previous chunk are not gated out (else coverage/AN is
        # undercounted at the chunk's first <=max_len sites). The sites/ref_ht read
        # keeps the un-widened sub_intervals so sites stay disjoint per chunk.
        vds_read_intervals = _expand_leading_edges(sub_intervals, max_ref_block_len)

    # Deferred until the layout hash is known: a chunk output must be namespaced by
    # its layout so a stale chunk from a prior run is never skipped-as-present or
    # merged in (the fan-out passes --chunk-output explicitly; this covers a manual
    # single-chunk run where it is omitted).
    if args.chunk_output is None:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            results_environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
            chrom=args.chrom,
        )
        args.chunk_output = _chunk_path(cov_and_an_ht_path, start, intervals_hash)
        logger.info("Auto-derived --chunk-output: %s", args.chunk_output)

    vds, sex_karyotype_field = _load_project_vds(
        project=project,
        environment=environment,
        partition_range=partition_range,
        sub_intervals=vds_read_intervals,
        filter_intervals=vds_filter_intervals,
        chrom=chrom,
        test=test,
        test_sample_subset=args.test_sample_subset,
    )

    ref_ht = _build_chunk_ref_ht(
        vds_filtered=vds,
        partition_count=n,
        chrom=chrom,
        sites_path=_vep_context_sites_path(test),
        sub_intervals=sub_intervals,
    )

    cov_and_an_ht = compute_all_release_stats_per_ref_site(
        vds,
        ref_ht,
        sex_karyotype_field=sex_karyotype_field,
        project=project,
        group_membership_ht=hl.read_table(group_membership_ht_path),
        reduce_min_aggs=args.reduce_min_aggs,
        experimental_densify=args.experimental_densify,
    )
    # Stamp the layout hash into globals for provenance: the union at merge time
    # inherits the first input's globals, so the final release HT records which
    # chunk-intervals layout produced it.
    cov_and_an_ht = cov_and_an_ht.annotate_globals(chunk_intervals_hash=intervals_hash)
    cov_and_an_ht.write(args.chunk_output, overwrite=True)
    logger.info("Wrote chunk [%d, %d) to %s", start, stop, args.chunk_output)


def _run_coverage_merge(
    input_paths: list[str],
    output_path: str,
    coalesce_to: int | None = None,
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
        before writing. For group-level merges, set to the number of input HTs
        in the group (one output partition per input) to avoid the
        partition-count blowup (sum-of-input-partitions). For the final merge,
        set to the desired
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


def _resolve_commit() -> str:
    """
    Return the gnomad_qc commit that relay containers should check out.

    Prefers the ``GNOMAD_QC_COMMIT`` env var, which ``--submit-orchestrator``
    sets in the orchestrator job: the in-job checkout comes from a GitHub
    archive tarball, which has no ``.git`` directory for ``git rev-parse`` to
    read. Falls back to ``git rev-parse HEAD`` when running from a real clone
    (a laptop or worktree).

    :return: Full commit hash.
    """
    commit = os.getenv("GNOMAD_QC_COMMIT")
    if commit:
        return commit
    return subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()


def _build_setup_command(
    commit: str,
    gcp_billing_project: str,
    methods_branch: str = "main",
) -> str:
    """
    Build commands to download gnomad_qc and gnomad_methods and set up the pipeline's hail configuration.

    Both repos are actively developed, so chunk containers pull them at
    runtime rather than relying on what's baked into the Docker image. The
    image provides ``hail`` and system dependencies.

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
    # Hail config so ``hl.init(backend="batch")`` finds the Batch billing
    # project and remote_tmpdir, plus the GCS requester-pays project.
    # The canonical path is the XDG-style ``~/.config/hail/config.ini``
    # (what ``hailtop.config.get_user_config_path()`` returns); also write
    # ``~/.hail/config.ini`` as a legacy fallback for older Hail versions.
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
        # Pin the pipeline's Hail version (the relay's Hail Python version
        # determines the JAR the QoB driver downloads, so pinning here pins the
        # entire pipeline). 0.2.139 is past the 0.2.138 requester-pays propagation
        # regression that forced an earlier 0.2.137 pin, and is what the freq
        # fan-out runs (via the v5_freq_batch:0.2.139 image); the GSA-key patch
        # above is kept until confirmed unnecessary on 0.2.139.
        #
        # WARNING: this version is a floor for every Table this pipeline READS, not
        # just a preference. Hail's table format is not backward compatible -- an
        # older reader hard-fails with "incompatible file format" on a table written
        # by a newer Hail. So any HT an input step writes (the vep_context sites HT,
        # the chunk-intervals JSON is exempt, group_membership, downsampling) must be
        # written by an environment whose Hail version is <= the version pinned here.
        # Lowering this pin without rewriting those inputs breaks the fan-out; JSON
        # artifacts are immune, which is one reason the sample artifacts are JSON.
        "/opt/venv/bin/pip install --quiet --upgrade --force-reinstall"
        " --no-deps hail==0.2.139\n"
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
    app-name / output-suffix / results-environment / QoB driver+worker sizing) is
    common to both. With
    ``chunk=True`` the read/compute flags a chunk relay also needs
    (``--read-subintervals-per-chunk``, ``--reduce-min-aggs``,
    ``--experimental-densify``, ``--test-sample-subset``, ``--chrom``,
    ``--test-region``, ``--test``) are
    appended; the merge relay (which only unions HTs) omits them. Per-job flags (``--chunk-*`` /
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
    if args.results_environment:
        # Forward the results-bucket override so a worker's fallback path resolution
        # matches the orchestrator's (the orchestrator already passes explicit
        # --chunk-output/--merge-output resolved with it).
        flags.append(f"--results-environment {args.results_environment}")
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
        if args.experimental_densify:
            flags.append("--experimental-densify")
        if args.test_sample_subset:
            flags.append("--test-sample-subset")
        if args.chrom:
            flags.append(f"--chrom {args.chrom}")
        if args.test_region:
            flags.append(f"--test-region {' '.join(args.test_region)}")
        # args.test is the canonical test flag (normalized at parse time to include
        # --test-n-partitions / --test-region). The worker reads its chunk from the
        # JSON by index, so it needs --test (for test-scoped paths), not the partition
        # count -- forward --test so the worker resolves the same test paths.
        if args.test:
            flags.append("--test")
    return " ".join(flags)


class _RelayJobSpec(NamedTuple):
    """One relay job's per-job config for :func:``_submit_relay_batch``."""

    name: str
    cpu: float
    memory: str
    storage: str
    chunk_attempts: int
    command: str


def _submit_relay_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    batch_name: str,
    job_specs: list[_RelayJobSpec],
    log_label: str,
) -> int | None:
    """
    Build and submit one Hail Batch of relay jobs sharing the same config.

    Shared skeleton for ``_submit_chunk_batch`` and ``_submit_merge_batch``: each
    relay job is a non-spot coordinator container (preemption mid-wait would orphan
    its inner QoB job) pinned to ``BATCH_REGIONS``, with ``--chunk-attempts``
    attempts (1 by default -- no auto-retry, see that flag's help). Parallelism
    comes from Hail Batch's own scheduler running the N jobs
    concurrently. No-ops (skips ``batch.run()``) when ``job_specs`` is empty.

    :param args: Parsed CLI args (reads ``batch_image``, ``batch_dry_run``).
    :param backend_kwargs: kwargs for the ``hb.ServiceBackend(...)`` constructor.
    :param batch_name: Hail Batch name.
    :param job_specs: Per-job config (name + QoB sizing + retry count + command).
    :param log_label: Noun for log messages ("chunk" / "merge").
    :return: Hail Batch id of the submitted batch, or None if nothing was submitted
        (empty ``job_specs``) or the batch was a ``--batch-dry-run``.
    """
    if not job_specs:
        logger.info(
            "  no pending %s jobs for %s; skipping batch.run()", log_label, batch_name
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
            # Relay is a coordinator waiting on its inner QoB job; preemption
            # mid-wait orphans that inner job, so relays are non-spot.
            j.spot(False)
            # Defaults to 1 attempt. Non-spot does not make an attempt unlosable
            # (a worker can still go 'not_responding'), and a retried relay cannot
            # cancel the orphaned inner QoB batch -- it just submits a second one
            # racing the first. See --chunk-attempts.
            j.n_max_attempts(spec.chunk_attempts)
            j.command(spec.command)

        logger.info(
            "Submitting Hail Batch '%s': %d %s jobs (dry_run=%s)",
            batch_name,
            len(job_specs),
            log_label,
            args.batch_dry_run,
        )
        submitted = batch.run(dry_run=args.batch_dry_run)
        # None under --batch-dry-run (nothing was submitted, so there is no id).
        return getattr(submitted, "id", None)
    finally:
        backend.close()


def _submit_orchestrator_batch(args: argparse.Namespace) -> None:
    """
    Submit THIS orchestrator invocation as one small non-spot Hail Batch job.

    Driving a full fan-out locally means keeping one process alive for the
    whole run (waves block in ``batch.run()``), which ties the campaign to a
    laptop staying awake. This instead re-invokes the same command line (minus
    ``--submit-orchestrator``) inside a Hail Batch job and returns immediately
    (``batch.run(wait=False)``). The orchestrator never initializes Hail — it
    only uses hailtop to submit relay batches, which works from inside a Batch
    job the same way each relay's nested QoB submission does — so the job needs
    no JVM sizing and fits the smallest container.

    The commit for relay checkouts is resolved here, where ``.git`` exists, and
    forwarded via the ``GNOMAD_QC_COMMIT`` env var (see ``_resolve_commit``).

    If the orchestrator job dies, resubmit this same command: completed chunks
    are skipped via their ``_SUCCESS`` markers and the ``_failed_chunks.json``
    manifest records what the dead run had dispatched.

    :param args: Parsed CLI args (reads the batch/billing config; everything
        else is forwarded verbatim to the wrapped run).
    :return: None.
    """
    commit = _resolve_commit()
    setup_cmd = _build_setup_command(
        commit,
        gcp_billing_project=args.gcp_billing_project,
        methods_branch=args.methods_branch,
    )
    forwarded = [a for a in sys.argv[1:] if a != "--submit-orchestrator"]
    command = (
        f"export GNOMAD_QC_COMMIT={commit}\n"
        "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/compute_coverage.py "
        + " ".join(shlex.quote(a) for a in forwarded)
    )

    backend_kwargs = {"billing_project": args.batch_billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir
    backend = hb.ServiceBackend(**backend_kwargs)
    try:
        batch = hb.Batch(name=f"{args.app_name}_orchestrator", backend=backend)
        j = batch.new_job(name="orchestrator")
        j.image(args.batch_image)
        j.cpu(0.5)
        j.memory("standard")
        j.storage("10Gi")
        j.regions(BATCH_REGIONS)
        # Long-lived, mostly idle coordinator: the worst profile to preempt
        # (losing it mid-campaign stops wave dispatch), and non-spot costs
        # cents at this size.
        j.spot(False)
        # A second attempt would dispatch waves concurrently with any still
        # unfinished from the first; recovery is manual resubmission against
        # the _SUCCESS/manifest state instead.
        j.n_max_attempts(1)
        j.command(setup_cmd + command)
        logger.info(
            "Submitting orchestrator job '%s_orchestrator' (commit %s," " dry_run=%s)",
            args.app_name,
            commit[:12],
            args.batch_dry_run,
        )
        submitted = batch.run(wait=False, dry_run=args.batch_dry_run)
        batch_id = getattr(submitted, "id", None)
        logger.info(
            "Orchestrator running as Hail Batch %s. Monitor via the Batch UI;"
            " chunk state lives in the _SUCCESS markers and"
            " _failed_chunks.json. Resubmit this same command to resume if the"
            " job dies.",
            batch_id,
        )
    finally:
        backend.close()


def _submit_chunk_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    chunk_indices: list[int],
    cov_and_an_ht_path: str,
    intervals_hash: str,
    setup_cmd: str,
    common_flags_str: str,
    script: str,
    wave_label: str | None = None,
) -> int | None:
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
        chunk writes to a sibling ``_chunks/<intervals_hash>/<idx>.chunk.ht``.
    :param intervals_hash: Content hash of the chunk-intervals layout, used to
        namespace each chunk's output path (see :func:`_chunk_path`).
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags from
        ``_build_relay_common_flags(args, chunk=True)``.
    :param script: Path to the script inside the relay container.
    :param wave_label: Optional suffix appended to the Hail Batch name to
        distinguish per-wave batches (e.g. ``"w003of049"``). None for a
        single (non-waved) batch.
    :return: Hail Batch id of the submitted wave (None under
        ``--batch-dry-run``), for the failed-chunk manifest.
    """
    project = args.project_name
    total = args.total_partitions
    ppc = args.partitions_per_chunk
    batch_name = (
        f"v5_cov_{project}_{total}p_{ppc}ppc_sub{args.read_subintervals_per_chunk}"
    )
    if args.cov_and_an_output_suffix:
        batch_name += f"_{args.cov_and_an_output_suffix}"
    if wave_label:
        batch_name += f"_{wave_label}"

    job_specs = []
    for idx in chunk_indices:
        path = _chunk_path(cov_and_an_ht_path, idx, intervals_hash)
        # Chunk identity is the chunk INDEX: the worker looks itself up in the
        # per-chunk JSON by --chunk-start (= idx), and --chunk-stop is idx+1. The
        # VDS read is interval-based (read_intervals from those sub-intervals), so the
        # old native-partition range is no longer used.
        command = (
            f"{setup_cmd}{script} --run-chunk"
            f" --chunk-start {idx} --chunk-stop {idx + 1}"
            f" --chunk-output {path}"
            f" {common_flags_str}"
        )
        job_specs.append(
            _RelayJobSpec(
                name=f"cov_chunk_{idx:06d}",
                cpu=args.chunk_cpu,
                memory=args.chunk_memory,
                storage=args.chunk_storage,
                chunk_attempts=args.chunk_attempts,
                command=command,
            )
        )
    return _submit_relay_batch(args, backend_kwargs, batch_name, job_specs, "chunk")


def _eligible_chunk_indices(
    args: argparse.Namespace,
) -> tuple[list[str | None], list[int], str]:
    """
    Enumerate fan-out chunks and the subset selected by ``--chrom``.

    Single source of truth shared by the chunk orchestrator and the merge so they
    always agree on which chunks exist. A ``--test-region`` run is one chunk (the worker
    derives its sub-intervals from the regions). Otherwise chunks come from the
    per-chunk precompute JSON (``--write-chunk-intervals``, required): it is read
    driver-side (no VDS access -- which is VPC-SC-blocked outside the perimeter) to get
    each chunk's contig, then only chunks whose contig is in ``--chrom`` are kept (all
    chunks when ``--chrom`` is unset).

    :param args: Parsed CLI args (reads ``test_region``, ``environment``, ``test``,
        ``chrom``).
    :return: ``(chunk_contigs, eligible, intervals_hash)`` -- the per-chunk contig
        (``None`` for the single ``--test-region`` chunk) indexed by chunk index, the
        list of eligible chunk indices after the ``--chrom`` filter, and the content
        hash of the chunk-intervals layout (``"test_region"`` for a ``--test-region``
        run, which has no JSON) used to namespace chunk outputs.
    """
    if args.test_region:
        chunk_contigs: list[str | None] = [None]
        intervals_hash = "test_region"
    else:
        intervals_path = _chunk_intervals_path(args.environment, args.test)
        if not file_exists(intervals_path):
            raise FileNotFoundError(
                f"chunk-intervals JSON not found at {intervals_path}. Run"
                " --write-chunk-intervals first (required: the fan-out and merge"
                " enumerate chunks by contig from it)."
            )
        with hfs.open(intervals_path) as f:
            data = json.load(f)
        chunk_contigs = [c["contig"] for c in data["chunks"]]
        intervals_hash = _chunk_intervals_hash(data)
    n_chunks = len(chunk_contigs)
    if args.chrom:
        eligible = [i for i in range(n_chunks) if chunk_contigs[i] == args.chrom]
        if not eligible:
            raise ValueError(
                f"No chunks match --chrom {args.chrom}; contigs in the precompute:"
                f" {sorted(c for c in set(chunk_contigs) if c is not None)}."
            )
    else:
        eligible = list(range(n_chunks))
    return chunk_contigs, eligible, intervals_hash


def _orchestrate_coverage_batch(
    args: argparse.Namespace, cov_and_an_ht_path: str
) -> None:
    """
    Fan coverage/AN compute out as relay chunk jobs in Hail Batch (one job per chunk).

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
    one live display may be active at once``. Within a wave, Hail Batch's
    own scheduler runs the relay jobs in parallel up to worker-pool capacity.
    Set ``--wave-size <= 0`` (or >= the pending count) for a single batch
    (legacy behavior).

    Note on the binding constraint: these relays and their nested QoB
    drivers/workers run on the **Hail Batch service's** autoscaled worker
    pool (provisioned in the Hail team's GCP project). So ``--wave-size``
    is bounded by the shared service's capacity and by any per-billing-project
    concurrent-core cap the Hail team configures for ``gnomad-production``
    — *not* by gnomad-production's own GCP CPU / in-use-IP quota.

    By default each chunk's inner QoB creates its own separate Hail Batch;
    pass ``--experimental`` to attach the inner QoB to this outer batch
    via ``HAIL_BATCH_ID``.

    Idempotency: chunks whose ``_SUCCESS`` already exists are skipped, so
    a failed/partial wave is resumed by simply rerunning this step (only
    not-yet-complete chunks are re-submitted). ``batch.run()`` does not
    raise on per-job failure, so a wave with some failed chunks does not
    abort the remaining waves. After the last wave, chunks still missing
    are automatically re-dispatched for up to ``--fanout-retry-passes``
    passes; whatever remains after that is recorded in
    ``_failed_chunks.json`` for a manual rerun.

    The merge step is intentionally separate; run ``--merge-cov-chunks``
    after this finishes.

    :param args: Parsed CLI args (reads ``project_name``, ``wave_size``,
        ``overwrite``, ``methods_branch``, ``gcp_billing_project``,
        ``batch_billing_project``, ``batch_remote_tmpdir``, plus everything
        ``_build_relay_common_flags`` and ``_submit_chunk_batch`` read --
        the chunk set comes from the chunk-intervals JSON, not
        ``total_partitions``/``partitions_per_chunk``).
    :param cov_and_an_ht_path: Canonical output cov_and_an HT path
        (resolved by ``main``). Each chunk writes to
        ``<cov_and_an_path>_chunks/<intervals_hash>/<idx>.chunk.ht``; the final HT
        at ``cov_and_an_ht_path`` is produced by ``--merge-cov-chunks``.
    :return: None.
    """
    project = args.project_name

    # Fail fast (before submitting thousands of jobs) if the preprocessed
    # vep_context sites HT the chunks read isn't there yet.
    sites_path = _vep_context_sites_path(args.test)
    if not file_exists(sites_path):
        raise FileNotFoundError(
            f"vep_context sites HT not found at {sites_path}. Run"
            " --write-vep-context-sites first (it is a prerequisite for the"
            " chunk compute)."
        )

    # Enumerate chunks (and the --chrom subset) via the shared helper, so the
    # orchestrator and the merge always agree on which chunks exist. intervals_hash
    # namespaces every chunk output by its layout so a stale chunk from a prior
    # --write-chunk-intervals run can never be skipped-as-present or merged in.
    chunk_contigs, eligible, intervals_hash = _eligible_chunk_indices(args)
    n_chunks = len(chunk_contigs)

    # Pre-compute pending chunk indices so the summary log is accurate
    # and so we can short-circuit when nothing is pending.
    if args.overwrite:
        pending_indices = list(eligible)
    else:
        # One directory listing instead of per-chunk serial existence probes.
        present = _list_present_chunk_indices(cov_and_an_ht_path, intervals_hash)
        pending_indices = [idx for idx in eligible if idx not in present]
    logger.info(
        "Coverage fan-out: %d chunks total, %d eligible%s, %d pending, %d skipped"
        " (overwrite=%s, project=%s)",
        n_chunks,
        len(eligible),
        f" (--chrom {args.chrom})" if args.chrom else "",
        len(pending_indices),
        len(eligible) - len(pending_indices),
        args.overwrite,
        project,
    )

    if not pending_indices:
        logger.info("All chunks already complete; nothing to submit.")
        return

    commit = _resolve_commit()
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

    run_id = _new_run_id()
    logger.info(
        "Orchestrator run id: %s (stamped into the failed-chunk manifest)", run_id
    )

    all_failed: list[int] = []
    wave_records: list[dict[str, Any]] = []
    n_dispatched = 0
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
        wave_batch_id = _submit_chunk_batch(
            args=args,
            backend_kwargs=backend_kwargs,
            chunk_indices=wave_indices,
            cov_and_an_ht_path=cov_and_an_ht_path,
            intervals_hash=intervals_hash,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
            wave_label=wave_label,
        )
        # batch.run() does not raise on per-job failure, so re-check this wave's
        # chunk outputs (one listing) and surface any that did not land (rather
        # than only discovering them at --merge-cov-chunks time).
        present = _list_present_chunk_indices(cov_and_an_ht_path, intervals_hash)
        failed = [idx for idx in wave_indices if idx not in present]
        n_dispatched += len(wave_indices)
        all_failed.extend(failed)
        wave_records.append(
            {
                "wave": wi,
                "batch_id": wave_batch_id,
                "n_dispatched": len(wave_indices),
                "failed_chunk_indices": failed,
            }
        )
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
        # Rewrite after every wave (not only at the end) so the rerun list is
        # already on disk if the orchestrator itself dies mid-run.
        manifest = _write_failed_chunks_manifest(
            cov_and_an_ht_path=cov_and_an_ht_path,
            intervals_hash=intervals_hash,
            failed=all_failed,
            n_dispatched=n_dispatched,
            run_id=run_id,
            commit=commit,
            app_name=args.app_name,
            waves=wave_records,
        )
        if all_failed:
            logger.warning(
                "%d of %d dispatched chunk(s) still missing; rerun list written to"
                " %s (run_id=%s)",
                len(all_failed),
                n_dispatched,
                manifest,
                run_id,
            )

    # Automatic retry passes for still-missing chunks. Safe where Batch-level
    # attempts (--chunk-attempts) are not: waves are sequential and
    # batch.run() blocks, so by the time a pass runs, every relay from the
    # previous dispatch has exited — a re-dispatch cannot race a live relay's
    # inner QoB batch for the same chunk output. Each pass is one batch (a
    # missing set large enough to need waves indicates a systematic failure
    # that retries won't fix). Chunks missing after the last pass stay in
    # _failed_chunks.json for follow-up.
    for pass_no in range(1, max(args.fanout_retry_passes, 0) + 1):
        if not all_failed:
            break
        retry_indices = sorted(set(all_failed))
        logger.info(
            "Retry pass %d/%d: re-dispatching %d missing chunk(s): %s%s",
            pass_no,
            args.fanout_retry_passes,
            len(retry_indices),
            retry_indices[:25],
            " ..." if len(retry_indices) > 25 else "",
        )
        retry_batch_id = _submit_chunk_batch(
            args=args,
            backend_kwargs=backend_kwargs,
            chunk_indices=retry_indices,
            cov_and_an_ht_path=cov_and_an_ht_path,
            intervals_hash=intervals_hash,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
            wave_label=f"retry{pass_no}",
        )
        present = _list_present_chunk_indices(cov_and_an_ht_path, intervals_hash)
        all_failed = [idx for idx in retry_indices if idx not in present]
        n_dispatched += len(retry_indices)
        wave_records.append(
            {
                "wave": f"retry{pass_no}",
                "batch_id": retry_batch_id,
                "n_dispatched": len(retry_indices),
                "failed_chunk_indices": all_failed,
            }
        )
        manifest = _write_failed_chunks_manifest(
            cov_and_an_ht_path=cov_and_an_ht_path,
            intervals_hash=intervals_hash,
            failed=all_failed,
            n_dispatched=n_dispatched,
            run_id=run_id,
            commit=commit,
            app_name=args.app_name,
            waves=wave_records,
        )
        if all_failed:
            logger.warning(
                "Retry pass %d/%d complete; %d chunk(s) STILL missing (see %s).",
                pass_no,
                args.fanout_retry_passes,
                len(all_failed),
                manifest,
            )
        else:
            logger.info(
                "Retry pass %d/%d complete; all chunks present.",
                pass_no,
                args.fanout_retry_passes,
            )


def _submit_merge_batch(
    args: argparse.Namespace,
    backend_kwargs: dict,
    group_indices: list[int],
    groups: list[list[str]],
    group_output_paths: list[str],
    setup_cmd: str,
    common_flags_str: str,
    script: str,
    level: int,
) -> None:
    """
    Build and submit one Hail Batch containing all pending group-merge jobs.

    Each job runs ``--run-merge`` over its assigned inputs and writes an
    intermediate per-group HT at
    ``<cov_and_an_path>_merge_groups_pp<ppc>_gs<group_size>/L<level>_<group_idx>.ht``
    (via ``_group_path``). Only intermediate levels go through this function; the
    final union is submitted separately by ``_orchestrate_coverage_merge`` and is
    the sole job that writes the canonical ``cov_and_an_ht_path``.
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
                chunk_attempts=args.chunk_attempts,
                command=command,
            )
        )
    _submit_relay_batch(args, backend_kwargs, batch_name, job_specs, "merge")


def _orchestrate_coverage_merge(
    args: argparse.Namespace, cov_and_an_ht_path: str
) -> None:
    """
    Recursive tree-reduce merge of per-chunk HTs into ``cov_and_an_ht_path``.

    Counterpart to ``_orchestrate_coverage_batch``: this orchestrator runs
    *after* ``--use-batch-fanout`` has produced all per-chunk HTs. It
    submits Hail Batch jobs (QoB-from-container per job) and exits without
    initializing Hail in this process.

    Discovery: chunks are enumerated from the per-chunk precompute JSON (via the
    same ``_eligible_chunk_indices`` the fan-out uses, so the two agree), filtered by
    ``--chrom`` when set. Every expected chunk must have a ``_SUCCESS`` marker; missing
    chunks fail loudly so the user re-runs the chunk orchestrator before merging.

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

    Safe to re-run: at each level, group HTs whose ``_SUCCESS`` exists are
    skipped; the final HT is skipped if it exists (unless ``--overwrite`` is passed).

    :param args: Parsed CLI args (reads ``project_name``,
        ``partitions_per_chunk``,
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

    # Enumerate chunks (and the --chrom subset) with the same helper the orchestrator
    # uses, so the merge expects exactly the chunks the fan-out produced. A --chrom run
    # merges only that contig's chunks into cov_and_an_ht_path.
    chunk_contigs, eligible, intervals_hash = _eligible_chunk_indices(args)
    n_chunks = len(chunk_contigs)
    logger.info(
        "Verifying %d expected chunk HTs exist (of %d total)...",
        len(eligible),
        n_chunks,
    )
    present = _list_present_chunk_indices(cov_and_an_ht_path, intervals_hash)
    missing = [i for i in eligible if i not in present]
    if missing:
        raise FileNotFoundError(
            f"--merge-cov-chunks: {len(missing)} of {len(eligible)} expected chunks"
            f" missing (first few idx: {missing[:5]}). Run --use-batch-fanout"
            " to (re)compute missing chunks first."
        )
    logger.info("All %d expected chunks present.", len(eligible))

    gs = args.merge_group_size

    # Precompute the level shape so we can log the full plan upfront.
    shape = [len(eligible)]
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

    commit = _resolve_commit()
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
    inputs = [_chunk_path(cov_and_an_ht_path, i, intervals_hash) for i in eligible]
    level = 1
    while len(inputs) > gs:
        n_in = len(inputs)
        n_out = (n_in + gs - 1) // gs
        groups = [inputs[i : i + gs] for i in range(0, n_in, gs)]
        out_paths = [
            _group_path(
                cov_and_an_ht_path,
                level,
                idx,
                args.partitions_per_chunk,
                args.merge_group_size,
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
    if not args.overwrite and file_exists(cov_and_an_ht_path):
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
        chunk_attempts=args.chunk_attempts,
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
        ``--merge-cov-chunks``): submit a Hail Batch and return.
      - ROLE 2 WORKER (``--run-chunk`` / ``--run-merge``): init Hail,
        run one chunk/merge unit, return — before the try/finally.
      - ROLE 3 IN-PROCESS PIPELINE (no role flag): init Hail, run the
        try/finally step chain (SETUP -> COMPUTE -> ASSEMBLE).
    """
    project = args.project_name
    environment = args.environment
    # args.test is normalized at parse time to fold in --test-n-partitions /
    # --test-region.
    test = args.test
    chrom = args.chrom
    overwrite = args.overwrite
    reduce_min_aggs = args.reduce_min_aggs
    results_environment = _results_environment(
        environment, test, args.results_environment
    )

    # ===================================================================
    # ROLE 1: ORCHESTRATOR — submit a Hail Batch of relay jobs, then return
    # straight out of main() (each per-chunk relay spawns its own QoB driver).
    # The configured _init_hail and the try/finally below are NEVER reached on
    # an orchestrator run; they belong to the separate ROLE 2 (worker) / ROLE 3
    # (in-process) invocations. The existence checks (file_exists) and chunk
    # listing both use hailtop.fs, so no Hail backend is started on an
    # orchestrator run.
    #   --use-batch-fanout : scatter the per-chunk compute.
    #   --merge-cov-chunks : tree-reduce the chunk HTs (run after fanout).
    #   --submit-orchestrator : run either of the above inside a Batch job
    #       instead of this process (submits and returns immediately).
    # ===================================================================
    if args.submit_orchestrator:
        _submit_orchestrator_batch(args)
        return

    if args.use_batch_fanout:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            results_environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
            chrom=args.chrom,
        )
        _orchestrate_coverage_batch(args, cov_and_an_ht_path)
        return

    if args.merge_cov_chunks:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            results_environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
            chrom=args.chrom,
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

    # Configured QoB / dataproc init — reached by ROLE 2 (worker) and ROLE 3
    # (in-process pipeline) only; ROLE 1 returned above without starting any Hail
    # backend (its existence checks / GCS listing use hailtop.fs).
    # ``--experimental`` attaches the QoB driver to an existing Hail Batch
    # (batch_id auto-resolved from HAIL_BATCH_ID inside _init_hail); otherwise
    # each run creates its own Hail Batch.
    #
    # Hail 0.2.139 workaround: the --experimental-densify scan crashes at compile
    # time ("scala.NotImplementedError ... SimpleSStream.settableTupleTypes")
    # whenever the scanned table has more partitions than Hail's branching factor
    # (default 50), because 0.2.139's tree-scan lowering zips a stream of streams
    # (LowerTableIR "table_scan_down_pass"; 0.2.137/0.2.138 are unaffected).
    # Raising the branching factor above the chunk's partition count keeps the
    # scan combine single-level, which never enters the broken lowering path. The
    # chunk reads at most --read-subintervals-per-chunk partitions (the interval
    # balancing can return fewer); the 2x margin only changes how many serialized
    # scan states the driver merges at once.
    branching_factor = None
    if args.run_chunk and args.experimental_densify:
        branching_factor = max(50, 2 * args.read_subintervals_per_chunk)

    _init_hail(
        log_name,
        environment,
        billing_project=getattr(args, "gcp_billing_project", None),
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{qc_temp_prefix(environment=environment, days=args.tmp_dir_days)}coverage_and_an_generation",
        experimental=args.experimental,
        branching_factor=branching_factor,
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
    # (--export-* / --merge-qual-hists, run as a separate job after fan-out
    #  + --merge-cov-chunks).
    # ===================================================================
    try:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            results_environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
            chrom=args.chrom,
        )

        # Union the per-contig HTs (each ..._<contig>.ht, written by a --chrom
        # merge) into the canonical HT. Contigs come from the precompute JSON;
        # each per-contig path is the canonical path with the contig appended.
        # Validate afterward with --validate-cov-and-an (no --chrom).
        if args.assemble_chrom_coverage:
            with hfs.open(_chunk_intervals_path(environment, test)) as f:
                contigs = sorted({c["contig"] for c in json.load(f)["chunks"]})
            inputs = [_apply_path_suffix(cov_and_an_ht_path, c) for c in contigs]
            _run_coverage_merge(
                inputs, cov_and_an_ht_path, coalesce_to=args.n_partitions
            )
            return

        group_membership_ht_path = _resolve_group_membership_ht_path(
            project,
            environment,
            test=test,
            reduce_min_aggs=reduce_min_aggs,
        )

        downsampling_ht_path = (
            get_aou_downsampling(test=test, environment=environment).path
            if project == "aou"
            else None
        )
        meta_ht = (
            meta(data_type="genomes", environment=environment).ht()
            if project == "aou"
            else None
        )

        # Materialize the chunk-invariant AoU sample artifacts once per run so each
        # chunk worker reads two small JSON files instead of rescanning the
        # ~330-partition sample tables (see _aou_sample_artifact_path). Cheap enough
        # to be worth running before any fan-out; generate_frequency reads the same
        # files, so whichever pipeline runs first serves both.
        if args.write_sample_artifacts:
            if project != "aou":
                raise ValueError(
                    "--write-sample-artifacts only applies to --project-name aou"
                    " (gnomAD has no sample-collision or AoU release artifact)."
                )
            # Imported here to match basics.get_aou_vds, which imports it lazily to
            # avoid a circular import.
            from gnomad_qc.v5.resources.meta import get_sample_id_collisions

            sc_ht = get_sample_id_collisions(environment=environment).ht()
            collisions = sorted(sc_ht.aggregate(hl.agg.collect_as_set(sc_ht.s)))
            path = _aou_sample_artifact_path("collisions.json", environment, test)
            with hfs.open(path, "w") as f:
                json.dump(collisions, f)
            logger.info("Wrote %d collision sample IDs: %s", len(collisions), path)

            # Must match get_aou_vds's release_only filter exactly, since it replaces
            # it: post-prefix AoU release sample IDs.
            release_samples = meta_ht.filter(
                (meta_ht.project_meta.project == "aou") & meta_ht.release
            ).s.collect()
            path = _aou_sample_artifact_path("release_samples.json", environment, test)
            with hfs.open(path, "w") as f:
                json.dump(release_samples, f)
            logger.info("Wrote %d release sample IDs: %s", len(release_samples), path)

        # Process the context HT once so the compute (fan-out workers especially) don't
        # have to re-compute it each time.
        if args.write_vep_context_sites:
            logger.info("Writing preprocessed vep_context sites HT...")
            sites_path = _vep_context_sites_path(test=test)
            check_resource_existence(
                output_step_resources={"vep_context_sites": [sites_path]},
                overwrite=overwrite,
            )
            site_intervals = None
            if test and args.test_region:
                # Explicit region (--test-region): scope the sites to these intervals
                # directly. vep_context is reference data, so no VDS probe is needed;
                # this preprocesses a dense, non-telomeric test region (the first VDS
                # partitions are the chr1 telomere).
                site_intervals = [_parse_region_interval(r) for r in args.test_region]
                logger.info(
                    "Test sites scoped to %d explicit --test-region interval(s): %s",
                    len(site_intervals),
                    args.test_region,
                )
            elif test:
                # Scope the test sites to the loci the test compute actually reads, so we
                # don't preprocess the whole genome for a tiny test. For a fan-out test,
                # use the chunk-intervals JSON's sub-intervals -- the same loci the
                # workers read (range(N) native partitions is the genome start = chr1, so
                # it cannot compose with --chrom). Falls back to the first N native VDS
                # partitions for the strict single-job path (no JSON). Validation confirms
                # the test sites cover the same loci as the compute.
                intervals_path = _chunk_intervals_path(environment, test)
                if file_exists(intervals_path):
                    _, eligible, _ = _eligible_chunk_indices(args)
                    with hfs.open(intervals_path) as f:
                        cm = json.load(f)
                    rg_name = cm["reference_genome"]
                    site_intervals = [
                        _interval_from_list(t, rg_name)
                        for i in eligible
                        for t in cm["chunks"][i]["intervals"]
                    ]
                    logger.info(
                        "Test sites scoped to %d sub-intervals (%d chunks) from the"
                        " chunk-intervals JSON.",
                        len(site_intervals),
                        len(eligible),
                    )
                elif chrom:
                    raise ValueError(
                        "--write-vep-context-sites --test --chrom requires the"
                        " chunk-intervals JSON (run --write-chunk-intervals first):"
                        " range(N) native partitions is the genome start (chr1) and"
                        " cannot scope to --chrom."
                    )
                else:
                    n_test_parts = args.test_n_partitions
                    if not n_test_parts or n_test_parts > 1000:
                        raise ValueError(
                            "--write-vep-context-sites --test (strict single-job"
                            " path) requires --test-n-partitions in 1..1000 to scope"
                            " the test sites; pass it, or run --write-chunk-intervals"
                            " first (fan-out test)."
                        )
                    probe = _probe_vds(
                        project, environment, list(range(n_test_parts)), chrom
                    )
                    site_intervals = _derive_chunk_locus_intervals(probe)
                    logger.info(
                        "Test sites scoped to the loci within the first %d partitions.",
                        n_test_parts,
                    )
            _build_vep_context_sites_ht(site_intervals).write(
                sites_path, overwrite=overwrite
            )

        if args.write_chunk_intervals:
            logger.info("Precomputing per-chunk read sub-intervals...")
            intervals_path = _chunk_intervals_path(environment, test=test)
            check_resource_existence(
                output_step_resources={"chunk_intervals": [intervals_path]},
                overwrite=overwrite,
            )
            chunk_intervals = _build_chunk_intervals(
                project=project,
                environment=environment,
                # --test-n-partitions scopes the probe to the first N VDS partitions
                # (overriding the full-VDS --total-partitions default) for a cheap test.
                total_partitions=args.test_n_partitions or args.total_partitions,
                partitions_per_chunk=args.partitions_per_chunk,
                n_sub=max(args.read_subintervals_per_chunk, 1),
                chrom=chrom,
            )
            with hfs.open(intervals_path, "w") as f:
                json.dump(chunk_intervals, f)
            logger.info("Wrote chunk-intervals JSON: %s", intervals_path)

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
                    reduce_min_aggs=reduce_min_aggs,
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
                    reduce_min_aggs=reduce_min_aggs,
                )
            # Coalesce on write. This HT is one row per sample (~365k for AoU) but
            # inherits the meta HT's ~330 partitions, and every consumer either joins
            # or key-indexes it -- which fans that consumer's read/collect stages to
            # ~330 tasks, per chunk, across the whole fan-out. Consumers cannot fix it
            # themselves: a lazy naive_coalesce does not propagate through a
            # key-indexed join, so the coalesced copy has to be materialized, and doing
            # it here means every consumer benefits instead of each re-checkpointing
            # its own copy.
            group_membership_ht.naive_coalesce(GROUP_MEMBERSHIP_N_PARTITIONS).write(
                group_membership_ht_path, overwrite=overwrite
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
            test_n = args.test_n_partitions
            strict_partition_range = list(range(test_n)) if test_n else None

            # Optionally co-partition the VDS and vep_context reads for more
            # efficient joins. --test-region (handled in the branch just below)
            # overrides this and reads via filter_intervals. Otherwise, the bare
            # --partitions-for-rep-on-read flag uses repartition_for_join's 1.1x
            # native-partition multiplier (the natural choice for a full run, where
            # there's no count to guess); a value gives an explicit partition count.
            # With --test-n-partitions the intervals come from the loci within those
            # partitions; otherwise from the (chrom-restricted) vep_context.
            strict_sub_intervals = None
            strict_filter_intervals = None
            if args.test_region:
                # Explicit region: read the VDS via filter_intervals (splits
                # straddling reference blocks; see the chunk worker) and the sites at
                # these intervals; skip the partition-based scoping entirely.
                strict_partition_range = None
                strict_sub_intervals = [
                    _parse_region_interval(r) for r in args.test_region
                ]
                strict_filter_intervals = strict_sub_intervals
            elif args.partitions_for_rep_on_read is not None:
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

            # The VDS read uses the same intervals as the sites via read_intervals,
            # EXCEPT for --test-region, which reads via filter_intervals (splitting
            # straddling reference blocks; see the chunk worker) -- so read_intervals
            # is dropped there.
            strict_vds_read_intervals = (
                None if args.test_region else strict_sub_intervals
            )
            vds, sex_karyotype_field = _load_project_vds(
                project=project,
                environment=environment,
                partition_range=strict_partition_range,
                sub_intervals=strict_vds_read_intervals,
                filter_intervals=strict_filter_intervals,
                chrom=chrom,
                test=test,
                test_sample_subset=args.test_sample_subset,
            )
            ref_ht = _build_chunk_ref_ht(
                vds_filtered=vds,
                partition_count=test_n or 0,
                chrom=chrom,
                sites_path=_vep_context_sites_path(test),
                sub_intervals=strict_sub_intervals,
            )

            cov_and_an_ht = compute_all_release_stats_per_ref_site(
                vds,
                ref_ht,
                sex_karyotype_field=sex_karyotype_field,
                project=project,
                group_membership_ht=hl.read_table(group_membership_ht_path),
                reduce_min_aggs=reduce_min_aggs,
            )
            if args.n_partitions is not None:
                # Checkpoint to temp only to stabilize before the coalesce; when
                # not coalescing, write straight to the destination (skip a full
                # temp write+read of the whole-VDS result).
                cov_and_an_ht = cov_and_an_ht.checkpoint(
                    new_temp_file(f"{project}_cov_and_an", "ht")
                ).naive_coalesce(args.n_partitions)
            cov_and_an_ht.write(cov_and_an_ht_path, overwrite=overwrite)

        # Validate that every intended vep_context site got AN/coverage. Runs
        # after cov_and_an is produced (by the compute just above, or by a
        # prior fan-out + merge) and before the gnomAD subtraction reads it.
        if args.validate_cov_and_an:
            logger.info("Validating cov_and_an HT covers all vep_context sites...")
            sites_ht = hl.read_table(_vep_context_sites_path(test))
            # A scoped run produces cov_and_an for only part of the genome, so validate
            # against just what ran (else the rest reads as "missing"). A fan-out
            # --test-n-partitions run covers only the first N sub-intervals (a partial
            # contig), so scope to exactly the chunks that ran, from the JSON; otherwise
            # --chrom scopes to whole contigs.
            intervals_path = _chunk_intervals_path(environment, test)
            if args.test_n_partitions and file_exists(intervals_path):
                _, eligible, _ = _eligible_chunk_indices(args)
                with hfs.open(intervals_path) as f:
                    cm = json.load(f)
                rg = cm["reference_genome"]
                sites_ht = hl.filter_intervals(
                    sites_ht,
                    [
                        _interval_from_list(t, rg)
                        for i in eligible
                        for t in cm["chunks"][i]["intervals"]
                    ],
                )
            elif args.chrom:
                sites_ht = hl.filter_intervals(
                    sites_ht, [hl.parse_locus_interval(args.chrom)]
                )
            merged_ht = hl.read_table(cov_and_an_ht_path)
            missing_ht = sites_ht.anti_join(merged_ht).checkpoint(
                new_temp_file("cov_an_missing_sites", "ht")
            )
            n_missing = missing_ht.count()
            n_sites = sites_ht.count()
            n_merged = merged_ht.count()
            if n_missing:
                # The compute emits an AN=0 row for every vep_context site it
                # processes (AN=0 where there's no VDS reference data), so within the
                # scoped set a site is never legitimately absent. Report a sample of
                # the missing loci to aid debugging:
                # a clustered count points at a dropped chunk, a scattered one at
                # the compute failing to emit AN=0 rows.
                sample = [str(x) for x in missing_ht.head(15).locus.collect()]
                raise ValueError(
                    f"cov_and_an HT at {cov_and_an_ht_path} is MISSING"
                    f" {n_missing}/{n_sites} vep_context site(s) after fan-out +"
                    " merge -- every site must be present (AN=0 if uncovered), so"
                    " this is a real loss, not a boundary artifact."
                    f" Sample missing loci: {sample}."
                    " Clustered -> a dropped chunk (re-run the fan-out/merge);"
                    " scattered -> the coverage compute failed to emit AN=0 rows."
                )
            # Check for duplicate sites
            if n_merged != n_sites:
                raise ValueError(
                    f"cov_and_an HT at {cov_and_an_ht_path} has {n_merged} rows"
                    f" but covers all {n_sites} vep_context sites — the"
                    f" {n_merged - n_sites} extra row(s) indicate duplicate loci"
                    " (e.g. overlapping chunk sub-intervals). Inspect before"
                    " releasing."
                )
            logger.info(
                "Validation passed: cov_and_an HT == vep_context sites exactly"
                " (%d rows, no missing, no extras).",
                n_sites,
            )

        # --- ASSEMBLE (gnomAD): subtract the consent-drop cohort from
        # the v4 release coverage/AN HTs to get gnomAD v5. ---
        if args.merge_gnomad_coverage:
            logger.info("Building gnomAD v5 coverage HT (subtracting consent-drop)...")
            merged_gnomad_coverage_ht_path = _gnomad_v5_merged_path(
                results_environment, "coverage", test
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
            merged_gnomad_an_ht_path = _gnomad_v5_merged_path(
                results_environment, "allele_number", test
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

        # --- ASSEMBLE (AoU): join AoU + gnomAD v5 and export the coverage and AN
        # HT/TSVs + merge qual hists. Run as a separate invocation after
        # ROLE 1 fan-out + --merge-cov-chunks have produced the AoU
        # cov_and_an HT. ---
        if args.export_coverage_release_files:
            logger.info("Exporting coverage release HT and TSV...")
            cov_ht_path = release_coverage_path(
                public=False,
                test=test,
                coverage_type="coverage",
                environment=results_environment,
            )
            cov_tsv_path = release_coverage_tsv_path(
                test=test, environment=results_environment
            )
            gnomad_coverage_ht_path = _gnomad_v5_merged_path(
                results_environment, "coverage", test
            )
            check_resource_existence(
                input_step_resources={"gnomad_coverage_ht": [gnomad_coverage_ht_path]},
                output_step_resources={
                    "cov_release_ht": [cov_ht_path],
                    "cov_tsv": [cov_tsv_path],
                },
                overwrite=overwrite,
            )
            # Coverage is gnomAD-only (AoU computes no coverage stats -- its
            # reference blocks carry no DP), so the release HT is the gnomAD
            # consent-drop result published as-is, with no AoU join.
            ht = hl.read_table(gnomad_coverage_ht_path)
            if args.n_partitions is not None:
                # Temp checkpoint only to stabilize before the coalesce;
                # otherwise the destination checkpoint below materializes it once.
                ht = ht.checkpoint(new_temp_file("gnomad_cov_release", "ht"))
                ht = ht.naive_coalesce(args.n_partitions)
            ht = ht.checkpoint(cov_ht_path, overwrite=overwrite)
            ht.export(cov_tsv_path)

        if args.export_an_release_files:
            logger.info("Exporting AN release HT and TSV...")
            an_ht_path = release_coverage_path(
                public=False,
                test=test,
                coverage_type="allele_number",
                environment=results_environment,
            )
            an_tsv_path = release_all_sites_an_tsv_path(
                test=test, environment=results_environment
            )
            gnomad_an_ht_path = _gnomad_v5_merged_path(
                results_environment, "allele_number", test
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
            ht = ht.select("AN")
            ht = ht.select_globals(
                strata_meta=ht.strata_meta,
                strata_sample_count=ht.strata_sample_count,
            )
            if args.n_partitions is not None:
                # join_aou_and_gnomad_an_ht already checkpoints; the extra temp
                # checkpoint is only to stabilize before the coalesce.
                ht = ht.checkpoint(new_temp_file("aou_and_gnomad_an_join", "ht"))
                ht = ht.naive_coalesce(args.n_partitions)
            ht = ht.checkpoint(an_ht_path, overwrite=overwrite)
            ht = ht.transmute(AN=ht.AN[0])
            ht.export(an_tsv_path)

        if args.merge_qual_hists:
            logger.info("Merging AoU + gnomAD v4 qual hists...")
            qual_hists_path = qual_hists(
                test=test, environment=results_environment
            ).path
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
            # by locus. ``distinct()`` deduplicates the multi-row-per-locus
            # rows that result from rekeying a (locus, alleles)-keyed HT
            # to locus-only; the ``_all`` hists are locus-level and
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
        "--results-environment",
        help=(
            "Override the bucket the coverage release artifacts (cov_and_an,"
            " gnomAD v5 merged tables, release files, qual-hists) write to,"
            " independent of the compute --environment. Default (unset): prod"
            " writes to 'dataproc', a --test run stays in the compute environment."
            " Set to 'dataproc' to smoke-test the batch-compute -> dataproc-bucket"
            " write path on a --test run before a prod run."
        ),
        choices=["rwb", "batch", "dataproc"],
        default=None,
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
            "Route the QoB init through ``hl.experimental.init`` instead of"
            " ``hl.init`` and attach the QoB driver to an existing Hail Batch."
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
        choices=["standard", "highmem"],
        help=(
            "Memory type for driver node. 'lowmem' is not offered: Hail Batch rejects"
            " it for JVM jobs ('400: jvm jobs cannot be on lowmem machines')."
        ),
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
            " (RegionPool), Python, and OS overhead. Rejected without"
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
        choices=["standard", "highmem"],
        help=(
            "Memory type for worker nodes. Default None leaves QoB's own default"
            " ('standard'), which is already generous -- chunk workers measured a"
            " 272 MiB peak, 7%% of a 1-core standard pod. 'lowmem' is not offered:"
            " Hail Batch rejects it for JVM jobs ('400: jvm jobs cannot be on lowmem"
            " machines')."
        ),
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
            " for A/B comparison. Applies to both writes (the compute / fan-out"
            " / merge) and reads (downstream steps), so pass the same suffix"
            " consistently."
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
        "--experimental-densify",
        action="store_true",
        help=(
            "compute coverage/AN via ``hl.experimental.densify`` on a merged sparse MT"
            " instead of ``hl.vds.to_dense_mt``, skipping the variant+reference per-cell"
            " coalesce. Output is equivalent; lets us  A/B the densify cost. Chunk /"
            " --use-batch-fanout path only."
        ),
    )
    parser.add_argument(
        "--reduce-min-aggs",
        action="store_true",
        help=(
            "Use the leaf-only ``reduce_to_minimal_groups`` optimization: build"
            " the group_membership_ht with ``reduce_to_minimal_groups=True``"
            " and pass ``reducible_aggs={'AN'}`` to ``compute_stats_per_ref_site``."
            " Must be set consistently for both --write-group-membership-ht"
            " and --compute-all-cov-release-stats-ht runs (the gmh's"
            " ``freq_reduced`` global determines whether reduction takes effect"
            " at compute time)."
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
    test_group.add_argument(
        "--test-sample-subset",
        action="store_true",
        help=(
            "When ``--test`` is set on AoU, additionally subsample the AoU"
            " sample set to ~0.1%% via ``meta_ht.sample(0.001)``. Default"
            " off: a ``--test`` run uses all samples but the partition /"
            " chrom subset, so AN values are stable and comparable"
            " across runs. Useful for a cheap tiny-cohort sanity check"
            " rather than a real-scale slice."
        ),
    )
    test_group.add_argument(
        "--test-n-partitions",
        type=int,
        default=None,
        help=(
            "Cheap test over the first N VDS partitions. With"
            " --write-chunk-intervals, scopes the probe to those partitions (overriding"
            " --total-partitions) -> ceil(N / --partitions-per-chunk) chunks, so the"
            " fan-out / merge / validate run just those (and --write-vep-context-sites"
            " --test scopes to the same loci); combine with --chrom to slice a single"
            " contig. With the strict single-job --compute-all-cov-release-stats-ht,"
            " reads the first N native VDS partitions (the genome start, so this mode"
            " does not compose with --chrom) and scopes the ref_ht to their loci; combine"
            " with --partitions-for-rep-on-read to co-partition just those. Default None:"
            " full run."
        ),
    )
    test_group.add_argument(
        "--chrom",
        default=None,
        help=(
            "Filter input data to a single contig (e.g. ``--chrom chr22``). The"
            " contig is appended to the output path suffix, so a --chrom run goes"
            " all the way through fan-out -> merge -> validate without clobbering"
            " another contig; assemble the per-contig HTs with"
            " --assemble-chrom-coverage. Default: no chrom filter. Independent of"
            " --test."
        ),
    )
    test_group.add_argument(
        "--test-region",
        nargs="+",
        default=None,
        help=(
            "For --write-vep-context-sites, the strict"
            " --compute-all-cov-release-stats-ht compute, and the fan-out chunk"
            " worker: scope the test to these explicit locus intervals instead of"
            " the first N VDS partitions (e.g."
            " ``--test-region chr1:55058666-55108666 chr1:55108666-55158666``)."
            " The vep_context sites HT and the VDS read are scoped to the region;"
            " straddling reference blocks are handled per path: the strict compute"
            " (and the chunk worker when --read-subintervals-per-chunk=1) read via"
            " filter_intervals with split_reference_blocks; the chunk worker fan-out"
            " (--read-subintervals-per-chunk>1) balances the region into that many"
            " read sub-intervals and widens the per-contig leading edge instead"
            " (sites stay un-widened). Use to test a dense, non-telomeric region --"
            " the first VDS partitions are the chr1 telomere."
            " Auto-enables test mode; mutually exclusive with --test-n-partitions."
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
        "--write-sample-artifacts",
        help=(
            "AoU only. Write the chunk-invariant sample artifacts (collisions.json,"
            " release_samples.json) that chunk workers read instead of each rescanning"
            " the ~330-partition sample-collisions and meta tables. Run once per run"
            " before the fan-out; without it every chunk pays those scans. Shares its"
            " paths with generate_frequency's --write-sample-artifacts."
        ),
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
            "Preprocess vep_context once into a deduped, telomere/centromere/chrM-"
            "stripped, locus-keyed sites HT (30-day storage) that every chunk reads"
            " co-partitioned. Prerequisite for the compute; run before"
            " --use-batch-fanout / --compute-all-cov-release-stats-ht. With --test,"
            " writes a cheap, separate _test table scoped to the loci the test compute"
            " reads -- from the chunk-intervals JSON, or (strict single-job path) the"
            " first --test-n-partitions VDS partitions -- instead of the whole genome."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--write-chunk-intervals",
        help=(
            "Precompute every chunk's balanced read sub-intervals in ONE VDS open and"
            " store them as a JSON file (30-day storage): a list of chunks, each tagged"
            " with its contig and never spanning a contig boundary. REQUIRED before the"
            " fan-out, which enumerates and --chrom-filters chunks from this file"
            " driver-side (instead of re-probing the VDS, ~400s per chunk). Uses"
            " --total-partitions / --partitions-per-chunk /"
            " --read-subintervals-per-chunk for sub-interval granularity; pass --chrom"
            " to restrict to some contigs. With --test, writes a separate _test file and"
            " --test-n-partitions overrides --total-partitions as the probe extent."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--validate-cov-and-an",
        help=(
            "Anti-join the vep_context sites HT against the merged cov_and_an HT and"
            " fail if any site is missing (i.e. dropped during fan-out + merge)."
            " Run after --merge-cov-chunks. Works for --test too (the _test sites"
            " table is scoped to the loci within the same partitions as the test compute)."
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
        help=(
            "Export the gnomAD v5 coverage HT and TSV. Coverage is gnomAD-only:"
            " AoU computes no coverage stats (its reference blocks carry no DP)."
        ),
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
        help="Exports the joint AoU + gnomAD v5 AN release HT and TSV file.",
        action="store_true",
    )

    parser.add_argument(
        "--merge-qual-hists",
        help="Merge variant quality histograms from AoU v8 and gnomAD v4 genomes.",
        action="store_true",
    )

    fanout_group = parser.add_argument_group(
        "Hail Batch fan-out for ``--compute-all-cov-release-stats-ht``.",
        "Submit one Hail Batch job per partition chunk; each chunk is a tiny"
        " relay container that spawns its own QoB driver to do the densify."
        " Idempotent: chunks whose ``_SUCCESS`` exists are skipped on rerun.",
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
        "--submit-orchestrator",
        action="store_true",
        help=(
            "Run the orchestrator itself inside Hail Batch instead of"
            " locally: submit one small non-spot job that re-invokes this"
            " command (minus this flag) and return immediately, so a full"
            " fan-out doesn't require keeping a laptop process alive for"
            " the whole run. Requires --use-batch-fanout or"
            " --merge-cov-chunks. Progress lives in the Batch UI and"
            " _failed_chunks.json; if the orchestrator job dies, resubmit"
            " this same command to resume (completed chunks are skipped)."
        ),
    )
    fanout_group.add_argument(
        "--fanout-retry-passes",
        type=int,
        default=1,
        help=(
            "After the last wave, automatically re-dispatch still-missing"
            " chunks up to this many extra passes (each pass is one Hail"
            " Batch of just the missing chunks). Unlike --chunk-attempts"
            " (Batch-level job retry), a pass only starts after every"
            " previously dispatched relay has exited, so it cannot race an"
            " orphaned inner QoB batch for the same chunk output. Chunks"
            " still missing after the final pass are recorded in"
            " _failed_chunks.json. Set 0 to disable. Default 1."
        ),
    )
    fanout_group.add_argument(
        "--merge-cov-chunks",
        action="store_true",
        help=(
            "Build and submit the recursive tree-reduce merge of per-chunk"
            " HTs produced by --use-batch-fanout. Intermediate levels (group"
            " merges sized by --merge-group-size) write to"
            " <cov_and_an_path>_merge_groups_pp<ppc>_gs<group_size>/"
            "L<level>_<idx>.ht; the merge recurses until one final union writes"
            " <cov_and_an_path>. Chunk discovery reads the chunk-intervals JSON"
            " (same as --use-batch-fanout); missing chunks fail loudly. Requires"
            " --environment=batch."
        ),
    )
    fanout_group.add_argument(
        "--assemble-chrom-coverage",
        action="store_true",
        help=(
            "Union the per-contig cov_and_an HTs (each ..._<contig>.ht, produced"
            " by running --merge-cov-chunks --chrom <contig> per contig) into the"
            " canonical HT. Contigs are enumerated from the chunk-intervals JSON."
            " Runs in-process (no --chrom). Validate the result with"
            " --validate-cov-and-an (no --chrom)."
        ),
    )
    fanout_group.add_argument(
        "--partitions-per-chunk",
        type=int,
        default=3,
        help=(
            "Number of VDS partitions per chunk job. Default 3, a shape from"
            " early testing; each chunk pays fixed overhead (relay + QoB driver +"
            " scan combine), so larger values are substantially cheaper at"
            " production scale (2026-08 calibration: ~$0.35/chunk fixed +"
            " ~$0.068/partition). IMPORTANT:"
            " keep identical across --write-chunk-intervals / --use-batch-fanout"
            " / --merge-cov-chunks. The chunk set is fixed in the JSON at"
            " --write-chunk-intervals time (n_chunks = ceil(extent /"
            " --partitions-per-chunk)) and read back by fan-out and merge; a"
            " mismatch misaligns the merge tree's group-HT paths."
        ),
    )
    fanout_group.add_argument(
        "--read-subintervals-per-chunk",
        type=int,
        default=50,
        help=(
            "Per-chunk: subdivide the loci within a chunk into this many"
            " data-balanced sub-intervals (computed once at"
            " --write-chunk-intervals time and stored in the JSON), then re-read"
            " both the VDS (``hl.vds.read_vds(intervals=...)``) and"
            " ``vep_context`` (``hl.read_table(_intervals=...)``) with those"
            " intervals. One partition per sub-interval on both sides,"
            " co-partitioned, no shuffle. Default 50. The VDS is always read by"
            " intervals, never by partition index."
        ),
    )
    fanout_group.add_argument(
        "--total-partitions",
        type=int,
        default=145192,
        help=(
            "Total VDS partition count to scatter across (the full-VDS prod"
            " extent). Default 145192 (prod AoU VDS partition count). Consumed by"
            " --write-chunk-intervals to size the precompute (the fan-out then"
            " reads that JSON); for a test, --test-n-partitions overrides this to"
            " scope the probe to the first N partitions."
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
            "GCS scratch path for the orchestrator's hb.ServiceBackend. Default"
            " None: ServiceBackend uses the Hail config on the host running the"
            " orchestrator."
        ),
    )
    fanout_group.add_argument(
        "--chunk-cpu",
        type=float,
        default=0.5,
        help=(
            "CPU per chunk relay container. Default 0.5: the relay only builds IR and"
            " waits on its inner QoB batch, so 1 core measured as pure waste"
            " ($0.052 -> $0.028 per chunk on the freq fan-out)."
        ),
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
        default="2",
        help=(
            "Cores for the QoB driver pod each chunk/merge relay spawns."
            " Power of two between 0.25 and 16 (as a string, e.g. '2' or"
            " '0.5'). Forwarded to the relay's --driver-cores; decoupled"
            " from the orchestrator's own --driver-cores. Default '2': with the"
            " highmem memory default this is a 13 GB driver, the smallest"
            " configuration measured to survive a coverage chunk (2-core standard,"
            " 8 GB, died; 1-core highmem, 6.5 GB, is smaller than what died and"
            " untested on coverage). The driver is billed for the whole worker"
            " barrier while idle, so cores beyond the memory they bring are"
            " near-pure waste."
        ),
    )
    fanout_group.add_argument(
        "--chunk-driver-memory",
        type=str,
        default="highmem",
        choices=["standard", "highmem"],
        help=(
            "Memory profile for the QoB driver pod each chunk/merge relay spawns."
            " Forwarded to the relay's --driver-memory; decoupled from the"
            " orchestrator's own --driver-memory. Default 'highmem': the driver dies"
            " during IR lowering of the big fused execute on 'standard' (measured on a"
            " PCSK9 coverage chunk -- 2-core standard, 8 GB, died with"
            " 'unexpected end of stream in jvm'; 2-core highmem completed). This is a"
            " JVM-total limit, not a RegionPool one, so it does not show up as a Hail"
            " out-of-memory error. 'lowmem' is rejected by Hail Batch for JVM jobs."
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
        choices=["standard", "highmem"],
        help=(
            "Memory profile per QoB worker pod the chunk/merge driver"
            " dispatches. Forwarded to the relay's --worker-memory; decoupled from"
            " the orchestrator's own --worker-memory. Default None leaves QoB's own"
            " default ('standard'), which is already generous -- chunk workers"
            " measured a 272 MiB peak, 7%% of a 1-core standard pod. 'lowmem' is NOT"
            " offered: Hail Batch rejects it for JVM jobs ('400: jvm jobs cannot be"
            " on lowmem machines')."
        ),
    )
    fanout_group.add_argument(
        "--chunk-attempts",
        type=int,
        default=1,
        help=(
            "Max attempts (n_max_attempts) per relay job (chunk and merge)."
            " Default 1 = no auto-retry, which is deliberate: a relay is only a"
            " coordinator for its inner QoB batch, and Batch retries an attempt it"
            " considers lost (e.g. reason 'not_responding') without being able to"
            " cancel that inner batch -- so a retry submits a SECOND QoB batch"
            " racing the first for the same chunk output. Failed chunks are"
            " recorded in _failed_chunks.json next to the chunk outputs; rerun"
            " --use-batch-fanout to pick them up."
        ),
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
        help=(
            "Chunk index for --run-chunk: the chunk's entry in the per-chunk"
            " precompute JSON (a --test-region run uses 0). Set by the orchestrator."
        ),
    )
    worker_group.add_argument(
        "--chunk-stop",
        type=int,
        default=0,
        help=(
            "Chunk index end (exclusive); the orchestrator sets it to"
            " --chunk-start + 1."
        ),
    )
    worker_group.add_argument(
        "--chunk-output",
        type=str,
        default=None,
        help=(
            "GCS output path for --run-chunk's per-chunk HT. Optional: when"
            " unset, auto-derived from the resolved cov_and_an HT path"
            " (plus ``--cov-and-an-output-suffix`` if set) and the chunk's"
            " ``--chunk-start`` index"
            " (``_chunks/<intervals_hash>/<chunk_start:08d>.chunk.ht``)."
            " The orchestrator (``--use-batch-fanout``) always passes this"
            " explicitly; this default fires only for standalone"
            " ``--run-chunk`` invocations (smoketests)."
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

    # --test-n-partitions and --test-region are themselves test runs, so fold them
    # into the canonical args.test once here. Everything downstream (main, the
    # orchestrators, and the helpers that resolve test-scoped paths) reads args.test,
    # so this keeps the test/prod path choice consistent across all of them.
    args.test = (
        args.test or args.test_n_partitions is not None or args.test_region is not None
    )

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
        "submit_orchestrator",
        "fanout_retry_passes",
        # Fan-out / merge orchestration args (only consumed by the batch
        # orchestrator + relays); passing them in a non-batch env is a mistake.
        "total_partitions",
        "partitions_per_chunk",
        "read_subintervals_per_chunk",
        "wave_size",
        "merge_group_size",
        "methods_branch",
        "batch_image",
        "chunk_attempts",
        "merge_cpu",
        "merge_memory",
        "merge_storage",
        "final_merge_storage",
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

    # --submit-orchestrator only wraps the two orchestrator roles.
    if args.submit_orchestrator and not (
        args.use_batch_fanout or args.merge_cov_chunks
    ):
        parser.error(
            "--submit-orchestrator requires --use-batch-fanout or"
            " --merge-cov-chunks (it wraps an orchestrator run in a Batch job)."
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

    if args.write_sample_artifacts and args.project_name != "aou":
        parser.error("--write-sample-artifacts requires --project-name to be 'aou'.")

    # NOTE: 'lowmem' is excluded from the choices of every driver/worker memory flag
    # rather than checked here -- Hail Batch rejects it for JVM jobs outright ("400:
    # jvm jobs cannot be on lowmem machines") and every QoB pod this script starts is
    # a JVM job, so argparse refuses it before anything is submitted.

    # Single-action entry points: each dispatches exactly one unit in main() and
    # returns, so pairing one with another entry point or with an in-process step
    # flag would silently skip the other. --run-chunk / --run-merge are set by the
    # orchestrators inside relay containers; users shouldn't pass them from the CLI.
    entry_points = {
        "--use-batch-fanout": args.use_batch_fanout,
        "--merge-cov-chunks": args.merge_cov_chunks,
        "--assemble-chrom-coverage": args.assemble_chrom_coverage,
        "--run-chunk": args.run_chunk,
        "--run-merge": args.run_merge,
    }
    # In-process steps that run sequentially in main()'s ROLE 3 block; these may be
    # combined with each other (e.g. compute + validate), but not with an entry point.
    step_flags = {
        "--write-aou-downsampling-ht": args.write_aou_downsampling_ht,
        "--write-group-membership-ht": args.write_group_membership_ht,
        "--write-sample-artifacts": args.write_sample_artifacts,
        "--write-vep-context-sites": args.write_vep_context_sites,
        "--write-chunk-intervals": args.write_chunk_intervals,
        "--compute-all-cov-release-stats-ht": args.compute_all_cov_release_stats_ht,
        "--validate-cov-and-an": args.validate_cov_and_an,
        "--merge-gnomad-coverage": args.merge_gnomad_coverage,
        "--merge-gnomad-an": args.merge_gnomad_an,
        "--export-coverage-release-files": args.export_coverage_release_files,
        "--export-an-release-files": args.export_an_release_files,
        "--merge-qual-hists": args.merge_qual_hists,
    }
    active_entry = [f for f, v in entry_points.items() if v]
    if len(active_entry) > 1:
        parser.error(f"{', '.join(active_entry)} are mutually exclusive entry points.")
    if active_entry:
        active_steps = [f for f, v in step_flags.items() if v]
        if active_steps:
            parser.error(
                f"{active_entry[0]} runs on its own and returns before the"
                f" in-process steps; {', '.join(active_steps)} would be silently"
                " skipped. Run them in a separate invocation."
            )
    if args.assemble_chrom_coverage and args.chrom:
        parser.error(
            "--assemble-chrom-coverage unions every per-contig HT into the"
            " canonical (contig-unscoped) path; do not pass --chrom."
        )
    if args.test_region and args.test_n_partitions is not None:
        parser.error(
            "--test-region and --test-n-partitions are mutually exclusive:"
            " --test-region scopes to explicit intervals, --test-n-partitions"
            " to the first N VDS partitions."
        )
    if args.run_chunk:
        if args.chunk_stop <= args.chunk_start:
            parser.error("--run-chunk requires --chunk-stop > --chunk-start.")
        # ``--chunk-output`` is optional: when not set it's auto-derived in
        # ``_run_coverage_chunk`` from the resolved cov_and_an HT path +
        # ``--cov-and-an-output-suffix`` + chunk_start. Useful for one-off
        # ``--run-chunk`` invocations (smoketests) — the orchestrator path
        # always passes ``--chunk-output`` explicitly.
    if args.run_merge:
        if args.merge_output is None or not args.merge_inputs:
            parser.error("--run-merge requires --merge-output and --merge-inputs.")

    main(args)
