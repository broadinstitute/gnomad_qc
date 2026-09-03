r"""
Compute coverage, allele number, and quality histograms for gnomAD v5 genomes.

v5 genomes = AoU v8 (new) + gnomAD v4 genomes minus the consent-drop set. This
script computes per-reference-site coverage (gnomAD only), AN, and qual
histograms (AoU only) per project and joins them into the v5 release HT/TSVs.

Execution roles (mutually exclusive; dispatched by early return in ``main``):

1. ORCHESTRATOR (``--use-batch-fanout`` / ``--merge-cov-chunks``): never
   initializes Hail; submits a Hail Batch of relay jobs and returns. Used for
   the prod-scale AoU compute.
2. WORKER (``--run-chunk`` / ``--run-merge``): one relay container; initializes
   Hail and runs exactly one chunk/merge unit. The chunk worker is what calls
   ``compute_all_release_stats_per_ref_site``.
3. IN-PROCESS PIPELINE (neither): the strict single-job compute (gnomAD
   consent-drop, dev/test AoU), the setup writers, and the release assembly.

A full prod AoU release is three separate invocations: fan-out compute,
merge-cov-chunks, then the in-process assembly. The setup HTs (downsampling,
group membership, sites, chunk intervals) are written by their own
in-process invocations first; the fan-out reads them from disk. The AoU sample
artifact JSONs (release samples, ID collisions) come from
``create_sample_qc_metadata_ht.py --write-vds-sample-artifact-jsons`` and are
read inside ``get_aou_vds``.

Workflow::

    1. (AoU) --write-aou-downsampling-ht
    2. --write-group-membership-ht: one aggregation cell per distinct
       membership pattern; every AN group is reconstructed by summing cells.
    3. --write-vep-context-sites: deduped, telomere/centromere/chrM-stripped
       sites HT -- the definition of "every site the compute must cover".
    4. Compute: --compute-all-cov-release-stats-ht (strict single job), or
       --write-chunk-intervals, then --use-batch-fanout, then
       --merge-cov-chunks. Chunk outputs are namespaced by a content hash of
       the intervals JSON: rerunning the fan-out retries only missing chunks,
       and a regenerated JSON (new layout, new hash) never mixes with old
       chunks.
    5. --validate-cov-and-an: the merged HT must cover every site exactly.
    6. (gnomad) --merge-gnomad-coverage / --merge-gnomad-an: subtract the
       consent-drop cohort from the v4 release HTs.
    7. (aou) Release: --export-coverage-release-files (gnomAD-only; AoU
       computes no coverage stats), --export-an-release-files,
       --merge-qual-hists (gnomAD v4 hists reused as-is).

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
    get_adj_expr,
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

# All Batch jobs are pinned to the region the input data lives in (AoU VDS,
# vep_context, and outputs are us-central1) to avoid inter-region GCS egress.
BATCH_REGIONS = ["us-central1"]

# gnomAD sample counts
GNOMAD_SAMPLE_COUNT = 71702
GNOMAD_CONSENT_DROP_SAMPLE_COUNT = 849

# chrM is called by a separate pipeline and absent from the v4 all-sites AN
# release, but present in vep_context v105 -- both the sites HT and the chunk
# intervals must drop it, or --validate-cov-and-an fails on sites the fan-out
# never computes.
EXCLUDED_CONTIGS = ["chrM"]

# The group_membership HT is one row per sample (~365k for AoU) and only feeds
# a broadcast index; far fewer partitions than the ~330 it inherits from meta.
GROUP_MEMBERSHIP_N_PARTITIONS = 10


def get_downsampling_ht(ht: hl.Table) -> hl.Table:
    """
    Get Table with downsampling groups for all samples.

    v5 downsamplings (AoU only): global 10,000 and 100,000, plus per-group
    sizes for afr/amr/nfe only (``gen_ancs_to_downsample``). Every group still
    gets a per-group index.

    :param ht: Input Table.
    :return: Table with downsampling groups.
    """
    logger.info(
        "Determining downsampling groups for AoU...",
    )
    downsamplings = DOWNSAMPLINGS["v5"]
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
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    The HT is reduced to cells (``reduce_to_cells=True``): one aggregation per
    distinct membership pattern instead of one per group, with every group --
    downsamplings included -- reconstructed by summing cells. Compute calls
    must pass ``reducible_aggs={"AN"}`` and pin every non-summable annotation
    (qual_hists, coverage_stats) to a single group via
    ``entry_agg_group_membership``; a pinned group that is not itself a cell
    (gnomAD coverage_stats -> adj) is aggregated over the union of its cells.

    :param meta_ht: Meta HT.
    :param project: "aou" or "gnomad". For "gnomad", filters ``meta_ht`` to the consent-drop samples.
    :param ds_ht: Optional downsampling HT (AoU only).
    :return: Group membership HT.
    """
    reduce_kwargs = dict(reduce_to_cells=True)
    if project == "aou":
        ht = generate_freq_group_membership_array(
            meta_ht,
            build_freq_stratification_list(
                sex_expr=meta_ht.sex_karyotype,
                gen_anc_expr=meta_ht.genetic_ancestry_inference.gen_anc,
                downsampling_expr=ds_ht[meta_ht.key].downsampling,
            ),
            **reduce_kwargs,
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
        # Filter to v4 consent-drop samples (v4 meta, not v5 project meta: this
        # step runs on Dataproc).
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
            **reduce_kwargs,
        )

    # Coalesce: this is a small per-sample lookup (~365k rows) that otherwise
    # inherits ~330 partitions from the meta HT and fans every downstream
    # read/collect into ~330 tasks. The coalesced copy must be materialized
    # (the write below): a lazy naive_coalesce does not propagate through a
    # key-indexed join, so consumers cannot fix it themselves.
    return ht.naive_coalesce(GROUP_MEMBERSHIP_N_PARTITIONS)


def _chunk_intervals_hash(data: dict[str, Any]) -> str:
    """
    Return a stable 16-hex-char content hash of the chunk-intervals JSON.

    ``repartition_for_join`` samples partition boundaries without a fixed seed,
    so each ``--write-chunk-intervals`` run yields a different layout. Chunk
    outputs are keyed only by index, so without namespacing by layout a stale
    chunk from a prior run would pass the existence skip-check and merge into
    overlapping loci. Computed over everything except any embedded
    ``intervals_hash``.

    :param data: Parsed chunk-intervals JSON.
    :return: First 16 hex chars of the SHA-256 of the canonical serialization.
    """
    payload = {k: v for k, v in data.items() if k != "intervals_hash"}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


def _test_region_hash(args: argparse.Namespace) -> str:
    """
    Return the layout hash for a ``--test-region`` run.

    Such runs have no chunk-intervals JSON, so the hash covers the region
    strings instead. Two test runs on different regions under the same output
    path then never see each other's chunk as already present. The
    sub-interval count is deliberately not included: the merge and validate
    steps do not pass it, and it does not change the chunk's contents.

    :param args: Parsed CLI args with ``test_region`` set.
    :return: ``test_region_`` plus a 16-hex-char hash.
    """
    return (
        f"test_region_{_chunk_intervals_hash({'test_region': list(args.test_region)})}"
    )


def _chunk_path(cov_and_an_ht_path: str, idx: int, intervals_hash: str) -> str:
    """
    Return a sibling per-chunk HT path under ``<cov_and_an_path>_chunks/<hash>/``.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param idx: Chunk index (zero-based).
    :param intervals_hash: Layout hash namespacing the output (see :func:`_chunk_intervals_hash`).
    :return: ``<cov_and_an_path>_chunks/<hash>/<idx:08d>.chunk.ht``.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_chunks/{intervals_hash}/{idx:08d}.chunk.ht"


def _list_present_chunk_indices(
    cov_and_an_ht_path: str, intervals_hash: str
) -> set[int]:
    """
    Return the set of chunk indices with a completed (``_SUCCESS``) output.

    One ``ls`` glob instead of per-chunk serial existence probes. Scoped to
    ``intervals_hash`` so only the current layout's chunks count; keying on
    ``_SUCCESS`` treats a partially-written chunk as absent.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param intervals_hash: Current chunk-intervals layout hash.
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

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param intervals_hash: Current chunk-intervals layout hash.
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

    Written next to the chunks (log scrollback is not durable) and rewritten
    after every wave. Rerunning ``--use-batch-fanout`` picks missing chunks up
    automatically; the manifest is informational. Successive runs overwrite the
    one path per layout, so each write is stamped with run id, time, commit,
    app name, and per-wave batch ids to distinguish attempts.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param intervals_hash: Current chunk-intervals layout hash.
    :param failed: Dispatched chunk indices with no ``_SUCCESS``.
    :param n_dispatched: Number of chunks dispatched so far in this run.
    :param run_id: Identifier for this orchestrator run.
    :param commit: gnomad_qc commit the relays ran.
    :param app_name: ``--app-name`` passed to the relays.
    :param waves: Per-wave records.
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

    Level-tagged so recursive merge levels don't overwrite each other; the
    directory encodes the tree-shape params so a rerun with different values
    writes fresh rather than mixing stale group HTs.

    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :param level: Merge-tree level (1-indexed).
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

    Prod results go to the gnomAD (``dataproc``) bucket regardless of compute
    environment; a test run stays in the compute environment's bucket unless
    ``override`` forces one.

    :param environment: Compute environment.
    :param test: Whether this is a test run.
    :param override: Explicit results environment (``--results-environment``).
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

    ``--chrom`` is appended to the suffix so every derived path is contig-scoped
    and per-contig runs never collide. Folded here (not by mutating ``args``)
    so the orchestrator and its relay workers always resolve identically.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param test: Whether to route to the test path.
    :param suffix: Optional suffix appended before the ``.ht`` extension.
    :param chrom: Optional single contig, appended to the suffix.
    :return: Fully resolved cov_and_an HT path.
    """
    path = coverage_and_an_path(
        test=test, data_set=project, environment=environment
    ).path
    full_suffix = f"{suffix}_{chrom}" if suffix and chrom else (suffix or chrom)
    return _apply_path_suffix(path, full_suffix)


def _group_membership_ht_path(project: str, environment: str, test: bool) -> str:
    """
    Return the path of the cell-reduced group_membership HT this script writes and reads.

    It is kept apart (``_cells`` suffix) from the full-shape HT at the
    ``group_membership`` resource path, which ``generate_frequency.py`` reads
    and which does not expand cells.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param test: Whether to route to the test path.
    :return: Fully resolved group_membership HT path.
    """
    path = group_membership(test=test, data_set=project, environment=environment).path
    return _apply_path_suffix(path, "cells")


def _gnomad_v5_merged_path(environment: str, coverage_type: str, test: bool) -> str:
    """
    Return the gnomAD v5 merged (consent-drop-subtracted) coverage/AN HT path.

    Written by ``--merge-gnomad-coverage`` / ``--merge-gnomad-an`` and read
    back by the export steps; ``test`` routes to a separate path.

    :param environment: Compute environment.
    :param coverage_type: "coverage" or "allele_number".
    :param test: If True, append ``_test``.
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
    Build a per-worker log name so concurrent workers don't clobber one log.

    :param run_chunk: ``--run-chunk`` flag.
    :param run_merge: ``--run-merge`` flag.
    :param chunk_start: Chunk start index.
    :param chunk_stop: Chunk stop index (exclusive).
    :param merge_output: Merge output HT path (label source for ``run_merge``).
    :return: Log-file basename.
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
    Derive per-contig locus intervals covering the filtered VDS reference_data.

    The vep_context HT and the AoU VDS have independent partition layouts, so
    the loci are derived from the VDS chunk and used to align ``ref_ht`` reads
    to the same span. Every caller uses the default ``n_subdivisions=1`` (one
    interval per contig).

    :param vds_filtered: VDS already filtered to the chunk's partitions.
    :param n_subdivisions: Equal-position sub-intervals per contig. Default 1.
    :param reference_genome: Reference genome name. Default ``"GRCh38"``.
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

    Thin wrapper over ``repartition_for_join``, used by the strict single-job
    path to co-partition the VDS and vep_context reads.

    :param n_partitions: Balanced partition count when positive; 0 uses
        ``repartition_for_join``'s default 1.1x-of-native multiplier.
    :param chrom: Optional single contig to restrict to first.
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


def _spans_sex_chromosome(intervals: Sequence[hl.utils.Interval]) -> bool:
    """
    Return whether any locus interval touches chrX or chrY.

    The sex-karyotype ploidy adjustment in ``compute_stats_per_ref_site`` only
    changes genotypes on those contigs, so a chunk that touches neither can
    skip it.

    :param intervals: Locus intervals (Python ``hl.utils.Interval`` objects).
    :return: True if any interval covers any part of an X or Y contig.
    """
    rg = intervals[0].start.reference_genome
    order = {c: k for k, c in enumerate(rg.contigs)}
    # An interval may span several contigs (partition bounds are not
    # contig-aligned), so take every contig from its start to its end.
    contigs: set[str] = set()
    for i in intervals:
        contigs.update(rg.contigs[order[i.start.contig] : order[i.end.contig] + 1])
    return bool(contigs & (set(rg.x_contigs) | set(rg.y_contigs)))


def compute_all_release_stats_per_ref_site(
    vds: hl.vds.VariantDataset,
    ref_ht: hl.Table,
    sex_karyotype_field: str | None,
    project: str,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
    interval_ht: hl.Table | None = None,
    group_membership_ht: hl.Table | None = None,
) -> hl.Table:
    """
    Compute coverage, allele number, and quality histograms per reference site.

    .. note::

        Running this prior to calculating frequencies removes the need for an
        additional densify there.

    :param vds: Input VDS. Reference data must carry ``END`` (and ``GQ``); a
        ``LEN`` field is added if missing. DP is read from variant data only.
    :param ref_ht: Locus-only sites Table defining the positions to aggregate at.
    :param sex_karyotype_field: Dotted path on the variant_data column struct to
        the sample's sex karyotype (sets per-sample ploidy on sex chromosomes),
        or None to skip the adjustment. It is the identity on autosomes, so
        pass None when ``ref_ht`` has no chrX/chrY sites (see
        :func:`_spans_sex_chromosome`): that also skips the per-sample
        karyotype lookup and its column checkpoint.
    :param project: "aou" or "gnomad". AoU adds ``qual_hists`` and computes no
        coverage stats; gnomAD computes coverage and AN only.
    :param coverage_over_x_bins: DP thresholds for the ``over_X`` fields.
        Default is :data:`COVERAGE_OVER_X_BINS`.
    :param interval_ht: Optional interval Table for partition pruning. Unused
        on the v5 path.
    :param group_membership_ht: Group-membership Table defining the per-stratum
        sample sets AN is fanned out across. Must be the cell-reduced HT from
        :func:`get_group_membership_ht`: AN is passed as the only
        ``reducible_aggs`` entry, so it is aggregated once per cell and every
        group is rebuilt by summing cells.
    :return: HT keyed by locus with per-stratum ``AN``; flat coverage fields
        (global adj group) for gnomAD; ``qual_hists`` (adj GQ histogram only)
        for AoU. AoU reference blocks carry no DP, so its coverage would
        measure ~0 genome-wide (chr1:55.06-55.16Mb: mean 0.048x, ``over_1``
        432 of 365,318) and is never computed, and for the same reason its DP
        histogram at reference sites is empty and is not computed either.
    """

    def _get_hists(qual_expr) -> hl.expr.Expression:
        """
        Build the adj-only GQ histogram from a [GQ, adj] pair.

        Selecting ``.qual_hists`` drops ``raw_qual_hists`` (approved for
        removal from v5), so Hail prunes those aggregations. No DP histogram:
        AoU reference blocks carry no DP, so it would be empty at reference
        sites (also approved for removal).
        """
        return hl.struct(
            qual_hists=qual_hist_expr(
                gq_expr=qual_expr[0],
                adj_expr=qual_expr[1] == 1,
                split_adj_and_raw=True,
            ).qual_hists
        )

    # Set up coverage bins.
    cov_bins = sorted(coverage_over_x_bins)
    rev_cov_bins = list(reversed(cov_bins))
    max_cov_bin = cov_bins[-1]
    cov_bins = hl.array(cov_bins)

    entry_agg_funcs = {"AN": get_allele_number_agg_func("LGT")}
    # Coverage stats are gnomAD-only: AoU ref blocks carry no DP (see docstring).
    if project == "gnomad":
        entry_agg_funcs["coverage_stats"] = get_coverage_agg_func(
            dp_field="DP", max_cov_bin=max_cov_bin
        )
    # Pin coverage_stats to the global adj group only: downstream uses
    # coverage_stats[0] exclusively, so per-strata coverage is wasted work. The
    # label must be "adj" (not "raw") -- with a pre-built group_membership_ht,
    # compute_stats_per_ref_site does not rewrite labels, so freq_meta[0] is
    # {"group": "adj"}. AN is omitted so it fans out across all strata.
    entry_agg_group_membership = {}
    if project == "gnomad":
        entry_agg_group_membership["coverage_stats"] = [{"group": "adj"}]
    if project == "aou":
        entry_agg_funcs["qual_hists"] = (lambda t: [t.GQ, t.adj], _get_hists)

        # qual_hists does its own adj filtering (adj is an argument), so it runs
        # on the raw group; _get_hists keeps only the adj-filtered result.
        entry_agg_group_membership["qual_hists"] = [{"group": "raw"}]

    logger.info(
        "Computing coverage, allele number, and optionally qual hists per reference site..."
    )

    vmt = vds.variant_data
    if sex_karyotype_field is not None:
        sex_expr = reduce(
            lambda x, field: x[field], sex_karyotype_field.split("."), vmt
        )
        vmt = vmt.annotate_cols(sex_karyotype=sex_expr)
    rmt = vds.reference_data
    # Hail < 0.2.134 VDSes lack LEN on reference data; drop this branch once
    # all inputs are >= 0.2.134.
    if "LEN" not in rmt.entry:
        rmt = rmt.annotate_entries(LEN=rmt.END - rmt.locus.position + 1)
    vds = hl.vds.VariantDataset(rmt, vmt)

    # Merge the VDS into one sparse MT so compute_stats_per_ref_site takes the
    # lean hl.experimental.densify scan instead of to_dense_mt's outer-join +
    # per-cell coalesce; outputs match. A missing ref allele skips
    # to_merged_sparse_mt's sequence-context branch (no FASTA needed).
    mtds = hl.vds.to_merged_sparse_mt(
        vds, ref_allele_function=lambda locus: hl.missing("str")
    )
    # Unify the genotype into LGT (ref blocks carry GT, variant sites LGT)
    # and drop GT so gt_field resolves to LGT. GT/LGT ploidy is identical,
    # so AN is unchanged. AoU arrives with adj already annotated; gnomAD keeps
    # LAD so compute_stats_per_ref_site annotates adj itself, after the
    # sex-ploidy adjustment (haploid calls get the haploid DP cutoff).
    mtds = mtds.annotate_entries(LGT=hl.coalesce(mtds.LGT, mtds.GT))
    mtds = mtds.select_entries(
        *[f for f in ("LGT", "GQ", "DP", "LAD", "adj", "END") if f in mtds.entry]
    )
    # AoU has nothing DP-based left to compute (adj is already annotated,
    # no coverage stats, no DP histogram), so drop DP from the densify. gnomAD
    # keeps it for coverage_stats and adj.
    if project == "aou":
        mtds = mtds.drop("DP")

    ht = compute_stats_per_ref_site(
        mtds,
        ref_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["GQ", "DP"] if project == "gnomad" else ["GQ"],
        reducible_aggs={"AN"},
        entry_agg_group_membership=entry_agg_group_membership,
        sex_karyotype_field=(None if sex_karyotype_field is None else "sex_karyotype"),
    )

    def _cov_stats(
        cov_stat: hl.expr.StructExpression,
    ) -> hl.expr.StructExpression:
        """Convert per-DP coverage_counter into cumulative ``over_X`` sample counts."""
        # Coverage was floored to the max bin upstream.
        count_expr = cov_stat.coverage_counter
        max_bin_expr = hl.int32(count_expr.get(max_cov_bin, 0))

        bin_expr = hl.range(hl.len(cov_bins) - 1, 0, step=-1)
        bin_expr = bin_expr.map(
            lambda i: hl.sum(
                hl.range(cov_bins[i - 1], cov_bins[i]).map(
                    lambda j: hl.int32(count_expr.get(j, 0))
                )
            )
        )
        bin_expr = hl.cumulative_sum(hl.array([max_bin_expr]).extend(bin_expr))
        # Sample counts, not fractions, to join with gnomAD v4 genomes.
        bin_expr = {f"over_{x}": bin_expr[i] for i, x in enumerate(rev_cov_bins)}
        return cov_stat.annotate(**bin_expr).drop("coverage_counter")

    # The sample-count global is annotated for both projects; downstream
    # readers expect it.
    ht = ht.annotate_globals(
        coverage_stats_meta_sample_count=ht.strata_sample_count[0],
    )
    if project == "gnomad":
        cov_stats_expr = _cov_stats(ht.coverage_stats[0])
        ht = ht.transmute(**cov_stats_expr)

    if project == "aou":
        # qual_hists comes back as a length-1 array; unwrap it.
        ht = ht.annotate(qual_hists=ht.qual_hists[0])
    return ht


def _rename_cov_annotations(
    ht: hl.Table,
    project: str,
    sample_count: int,
    coverage_over_x_bins: Sequence[int] = COVERAGE_OVER_X_BINS,
) -> hl.Table:
    """
    Rename coverage annotations (suffix ``_<project>``) prior to merging.

    Transforms ``mean`` back into ``sum`` using ``sample_count``; for the v4
    release HT also converts fraction over-X bins back to sample counts.

    :param ht: Input HT.
    :param project: Project name.
    :param sample_count: Number of samples in ``ht``.
    :param coverage_over_x_bins: Boundaries for the over-X fields. Default is :data:`COVERAGE_OVER_X_BINS`.
    :return: Renamed HT.
    """
    ht = ht.transmute(sum=ht.mean * sample_count)

    row_fields = list(ht.row_value)
    rename_dict = {f: f"{f}_{project}" for f in row_fields}
    ht = ht.rename(rename_dict)
    if project == "gnomad_release":
        # Fraction over-X bins back to sample counts.
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

    Only the consent-drop subtraction remains as a caller (AoU computes no
    coverage stats, so there is no coverage sum). ``median_approx`` is not
    merged.

    :param ht: Input HT with both projects' annotations.
    :param project_1: Minuend project.
    :param project_2: Subtrahend project.
    :param coverage_over_x_bins: Boundaries for the over-X fields. Default is :data:`COVERAGE_OVER_X_BINS`.
    :return: Merged fields.
    """
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
    Subtract consent-drop samples from the v4 coverage release to create the gnomAD v5 coverage HT.

    :param gnomad_ht: Coverage HT of the consent-drop samples only.
    :param gnomad_release_ht: gnomAD v4 genomes coverage release HT.
    :param coverage_over_x_bins: Boundaries for the over-X fields. Default is :data:`COVERAGE_OVER_X_BINS`.
    :param gnomad_sample_count: v4 release genome count. Default `GNOMAD_SAMPLE_COUNT`.
    :param consent_drop_count: Consent-drop genome count. Default `GNOMAD_CONSENT_DROP_SAMPLE_COUNT`.
    :return: gnomAD v5 genomes coverage HT in the release schema (``mean``,
        ``median_approx``, ``total_DP``, fraction ``over_X``).
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

    # Back to the v4 release schema: mean and fraction over-X bins from the
    # subtracted sums; median_approx is kept from the v4 release as-is. This
    # HT is exported for release unchanged (AoU computes no coverage stats).
    gnomad_ht = gnomad_ht.select(
        mean=gnomad_ht.sum_gnomad / gnomad_v5_count,
        median_approx=gnomad_ht.median_approx_gnomad_release,
        total_DP=gnomad_ht.total_DP_gnomad,
        **{
            f"over_{x}": gnomad_ht[f"over_{x}_gnomad"] / gnomad_v5_count
            for x in coverage_over_x_bins
        },
    )

    gnomad_ht = gnomad_ht.select_globals()
    gnomad_ht = gnomad_ht.annotate_globals(
        coverage_stats_meta_sample_count=gnomad_v5_count,
    )
    return gnomad_ht


def _rename_fields(
    ht: hl.Table, field_name: str, project: str, rename_globals: bool
) -> hl.Table:
    """
    Rename ``field_name`` (and optionally the strata globals) with a ``_<project>`` suffix.

    Used for AN and qual hists; coverage goes through ``_rename_cov_annotations``.

    :param ht: Input HT.
    :param field_name: Row field to suffix (e.g. ``"AN"``).
    :param project: Project name.
    :param rename_globals: Whether to rename ``strata_meta``/``strata_sample_count`` too.
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

    The shared strata's combined AN is written first as the unlabeled "global"
    entries; the AoU-only AN (including the downsampling strata) is appended as
    an ``{"subset": "aou"}``-tagged subset.

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
    # Keep only the shared strata as the unlabeled globals; the AoU-only
    # downsampling strata appear solely in the "aou" subset below.
    global_idx = [i for i, d in enumerate(joint_strata_meta) if "downsampling" not in d]
    global_meta = [joint_strata_meta[i] for i in global_idx]
    ht = ht.annotate(AN=hl.array([joint_an[i] for i in global_idx]))
    ht = ht.annotate_globals(
        strata_meta=global_meta,
        strata_sample_count=hl.array(
            [count_arrays_dict["counts"][i] for i in global_idx]
        ),
    )

    # Append the full AoU AN as the "aou" subset. AN_aou is always defined
    # (AoU is the left side of the join), so the three parallel arrays extend
    # in lockstep with no fill.
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
        Qual hists were not computed for the gnomAD v4 genomes release, so v5
        reuses the v4 hists as-is (no consent-drop subtraction).

    :param aou_ht: AoU qual hists HT.
    :param gnomad_ht: gnomAD qual hists HT.
    :return: Joined HT.
    """
    aou_ht = _rename_fields(aou_ht, "qual_hists", "aou", rename_globals=False)
    gnomad_ht = _rename_fields(gnomad_ht, "qual_hists", "gnomad", rename_globals=False)
    ht = aou_ht.join(gnomad_ht, "left")
    # Only the adj GQ histogram survives in v5: the v4 raw_qual_hists are not
    # merged, and dp_hist_all is not computed for AoU (no DP on its reference
    # blocks), so the v4 one is dropped rather than released unmerged.
    qual_hists = [
        "gq_hist_all",
    ]
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

    Used in test runs where the two tables' first-N-partition slices would not
    overlap; filtering to the source's per-contig min..max guarantees overlap
    cheaply.

    :param target_ht: Table to filter.
    :param source_ht: Table whose locus ranges define the filter intervals.
    :return: Filtered ``target_ht``.
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
    partition_range: list[int] | None = None,
    sub_intervals: list[hl.utils.Interval] | None = None,
    filter_intervals: list[hl.utils.Interval] | None = None,
    chrom: str | None = None,
    test: bool = False,
    test_sample_subset: bool = False,
) -> tuple[hl.vds.VariantDataset, str]:
    """
    Load the per-project VDS with consistent test/subsample handling.

    Shared by the chunk worker and the strict path: routes ``sub_intervals``
    (read-time co-partitioning) over ``partition_range``, synthesizes AoU DP
    from LAD (the AoU v8 VDS lacks DP), annotates AoU adj (standard cutoffs at
    variant sites; all ref-site genotypes pass), and applies the optional AoU
    test subsample.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param partition_range: VDS partition indices, or None for the full VDS.
    :param sub_intervals: Locus intervals for read-time partitioning.
    :param filter_intervals: Intervals passed to ``get_*_vds`` as
        ``filter_intervals`` (splits straddling reference blocks). Used by
        ``--test-region``.
    :param chrom: Optional single contig to filter to.
    :param test: Whether this is a test run.
    :param test_sample_subset: If True (AoU and ``test``), subsample ~0.1% of
        the AoU samples.
    :return: Tuple of ``(vds, sex_karyotype_field)``.
    """
    if project == "aou":
        sex_karyotype_field = "meta.sex_karyotype"
        # get_aou_vds reads the precomputed sample artifact JSONs (release
        # samples, ID collisions; see write_aou_vds_sample_jsons) and falls
        # back to the sample-table scans only when they are absent.
        # Release samples are already free of hard-filtered samples, and AN,
        # GQ and DP never read alleles, so skip the hard-filter Table scan, the
        # dead-allele LA/LAD rewrite and the count_cols() logging -- each is
        # extra per-chunk work with no effect on the output.
        vds = get_aou_vds(
            release_only=True,
            remove_hard_filtered_samples=False,
            remove_dead_alleles=False,
            log_sample_counts=False,
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
        # Variant sites: the usual gnomAD adj cutoffs (GQ >= 20, DP >= 10 with
        # DP approximated as sum(LAD), AB >= 0.2 for het calls).
        vmt = vmt.annotate_entries(DP=hl.sum(vmt.LAD))
        vmt = vmt.annotate_entries(adj=get_adj_expr(vmt.LGT, vmt.GQ, vmt.DP, vmt.LAD))
        rmt = vds.reference_data
        # Ref sites: all genotypes pass. AoU ref blocks are GQ-banded to
        # {20, 30, 40} and GQ0 is never written, so there is no principled GQ
        # cutoff to apply.
        rmt = rmt.annotate_entries(adj=True)
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

    Loads via ``filter_partitions`` / ``filter_intervals`` only -- no release
    or subset filtering, no DP synthesis (AoU's base sample exclusion still
    applies); the caller derives loci from ``reference_data``.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param partition_range: VDS partition indices to probe, or None for the full VDS.
    :param chrom: Optional single contig to filter to.
    :param filter_intervals: Optional locus intervals to bound the probe to
        (splits straddling reference blocks). Used by ``--test-region``.
    :return: Probe VDS.
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


def _vep_context_sites_path(
    test: bool = False, test_region: list[str] | None = None
) -> str:
    """
    Return the path to the preprocessed vep_context sites HT (30-day storage).

    Prod is genome-wide; a ``--test`` write is scoped and goes to a separate
    ``_test`` path. A ``--test-region`` write is further namespaced by the
    region(s) so that test chains for different regions never read each
    other's sites HT (a stale sites HT silently yields empty chunks).

    :param test: If True, return the test-scoped sites path.
    :param test_region: The ``--test-region`` strings, if any (test only).
    :return: GCS path to the preprocessed sites HT.
    """
    name = "vep_context_sites"
    if test:
        name += "_test"
        if test_region:
            name += "_" + "__".join(
                re.sub(r"[^A-Za-z0-9]+", "_", r) for r in test_region
            )
    # Pin to dataproc so it is accessible everywhere
    return f"{qc_temp_prefix(environment='dataproc', days=30)}{name}.ht"


def _build_vep_context_sites_ht(
    intervals: list[hl.utils.Interval] | None = None,
) -> hl.Table:
    """
    Build the locus-keyed, deduped, telomere/centromere/chrM-stripped sites HT.

    The per-ref-site set the compute aggregates at; previously the dedup +
    strip ran inside every chunk, now once here.

    :param intervals: Optional locus intervals to scope the table to
        (read-time pruned); None preprocesses the whole genome.
    :return: Locus-keyed sites HT.
    """
    if intervals is not None:
        ref_ht = vep_context.versions["105"].ht(read_args={"_intervals": intervals})
    else:
        ref_ht = vep_context.versions["105"].ht()
    ref_ht = ref_ht.key_by("locus").select().distinct()
    # Drop telomeres/centromeres and EXCLUDED_CONTIGS here so this table stays
    # the single definition of "every site the compute must cover" and
    # --validate-cov-and-an remains an exact-equality check.
    drop_intervals = telomeres_and_centromeres.ht().interval.collect()
    drop_intervals += [
        hl.eval(hl.parse_locus_interval(c, reference_genome="GRCh38"))
        for c in EXCLUDED_CONTIGS
    ]
    ref_ht = hl.filter_intervals(ref_ht, drop_intervals, keep=False)
    return ref_ht


def _chunk_intervals_path(environment: str, test: bool = False) -> str:
    """
    Return the path of the precomputed per-chunk read sub-intervals JSON (30-day storage).

    Written once by ``--write-chunk-intervals``: one VDS open derives every
    chunk's balanced sub-intervals. Plain JSON (not a Hail Table) so the
    orchestrator, merge, and workers read it driver-side with no QoB job.
    Holds ``chunks`` (``{"contig", "intervals"}`` per chunk index, never
    spanning a contig boundary), ``ref_block_max_length``,
    ``reference_genome``, and ``read_subintervals_per_chunk``.

    :param environment: Compute environment.
    :param test: If True, return the test-scoped intervals path.
    :return: GCS path to the chunk-intervals JSON.
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

    A concrete ``hl.utils.Interval`` (not an ``IntervalExpression``) so it can
    go to ``read_args={"_intervals": ...}``; half-open so adjacent intervals
    stay disjoint. Used by ``--test-region``.

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

    ``repartition_for_join`` returns a contiguous partition of the keyspace, so
    an interval can span the tail of one contig and the head of the next;
    contig-keyed chunking needs every interval on exactly one contig. The
    pieces tile the original exactly.

    :param intervals: Sorted locus intervals.
    :param reference_genome: Reference-genome name.
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
            # Drop an empty tail piece (e.g. end at chrB:1 exclusive).
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


def _snap_edges_to_bounds(
    intervals: list[hl.utils.Interval],
    bounds: dict[str, tuple[int, int, bool]],
) -> list[hl.utils.Interval]:
    """
    Extend each contig's first interval back to a bound start and its last forward to a bound end.

    ``repartition_for_join`` bounds span the probe's first to last reference-block
    key, so sites before the first block (or after the last) in the requested span
    would never be read and would be dropped instead of emitted with AN=0.

    :param intervals: Sorted locus intervals, each on one contig.
    :param bounds: Per-contig ``(start position, end position, includes_end)``;
        contigs absent from ``bounds`` are left as is.
    :return: Intervals with per-contig edges snapped to ``bounds``.
    """
    out = list(intervals)
    by_contig: dict[str, list[int]] = {}
    for i, iv in enumerate(out):
        by_contig.setdefault(iv.start.contig, []).append(i)
    for contig, idxs in by_contig.items():
        if contig not in bounds:
            continue
        lo, hi, inc_hi = bounds[contig]
        rg = out[idxs[0]].start.reference_genome
        first = out[idxs[0]]
        if lo < first.start.position:
            out[idxs[0]] = hl.Interval(
                hl.Locus(contig, lo, reference_genome=rg),
                first.end,
                includes_start=True,
                includes_end=first.includes_end,
            )
        last = out[idxs[-1]]
        if hi > last.end.position or (
            hi == last.end.position and inc_hi and not last.includes_end
        ):
            out[idxs[-1]] = hl.Interval(
                last.start,
                hl.Locus(contig, hi, reference_genome=rg),
                includes_start=last.includes_start,
                includes_end=inc_hi,
            )
    return out


def _region_bounds(
    regions: list[hl.utils.Interval],
) -> dict[str, tuple[int, int, bool]]:
    """
    Per-contig ``(min start, max end, includes_end)`` of the ``--test-region`` intervals.

    :param regions: Parsed ``--test-region`` intervals.
    :return: Bounds for ``_snap_edges_to_bounds``.
    """
    bounds: dict[str, tuple[int, int, bool]] = {}
    for r in regions:
        c = r.start.contig
        lo, hi, inc = bounds.get(c, (r.start.position, r.end.position, r.includes_end))
        if r.end.position > hi or (r.end.position == hi and r.includes_end):
            hi, inc = r.end.position, r.includes_end
        bounds[c] = (min(lo, r.start.position), hi, inc)
    return bounds


def _build_chunk_intervals(
    project: str,
    environment: str,
    total_partitions: int,
    partitions_per_chunk: int,
    n_sub: int,
    chrom: str | None = None,
) -> dict[str, Any]:
    """
    Precompute every chunk's balanced read sub-intervals in one VDS open, by contig.

    Replaces the per-chunk probe. Derives ``n_chunks * n_sub`` balanced
    sub-intervals over all reference-data loci, splits any that straddle a
    contig boundary, groups by contig, and slices each contig's run into
    chunks of ``n_sub``. No chunk crosses a contig boundary, chunks are
    disjoint, and their union covers every contig end to end.

    :param project: "aou" or "gnomad".
    :param environment: Compute environment.
    :param total_partitions: Leading VDS partitions to derive intervals over
        (``--test-n-partitions`` overrides this for tests).
    :param partitions_per_chunk: Sets ``n_chunks = ceil(total_partitions / this)``.
    :param n_sub: Sub-intervals per chunk (``--read-subintervals-per-chunk``).
    :param chrom: Optional single contig to restrict the precompute to.
    :return: JSON-serializable dict: ``read_subintervals_per_chunk``,
        ``ref_block_max_length``, ``reference_genome``, ``chunks``.
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
    # Keep the chunk set on the same contigs as the sites HT (EXCLUDED_CONTIGS);
    # safe to drop only after the contig split.
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
    # The bounds span the first to last reference-block key; sites outside
    # that span on each contig must still be read (AN=0), not dropped.
    contig_subs = _snap_edges_to_bounds(
        contig_subs,
        {c: (1, rg.lengths[c], True) for c in {iv.start.contig for iv in contig_subs}},
    )
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
    # --chrom filtering and per-contig disjointness both require single-contig chunks.
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

    With ``sub_intervals``, reads co-partitioned with the VDS read (the VDS
    edge may be widened for straddling reference blocks; the sites stay
    un-widened so chunks emit disjoint sites). Otherwise reads the whole sites
    HT and filters to ``chrom`` or the chunk's loci.

    :param vds_filtered: VDS already filtered to the chunk's partitions.
    :param partition_count: When > 0 and no ``sub_intervals``/``chrom``, filter
        ``ref_ht`` to the chunk's loci. 0 means whole.
    :param chrom: Optional single contig (takes precedence over the
        ``partition_count`` filter).
    :param sites_path: Path to the preprocessed sites HT.
    :param sub_intervals: Optional read-time intervals; take precedence over both.
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

    Each chunk's ``read_vds(intervals=...)`` gates reference blocks by START
    locus, so a block straddling in from the previous chunk would be dropped,
    silently undercounting coverage/AN at the chunk's leading sites. Widening
    only the VDS read's leading edge pulls the straddlers back in; the sites
    read stays unchanged, so no site is computed twice.
    ``ref_block_max_length`` caps how far back a live block can start, so this
    catches every straddler and nothing earlier.

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

    Invoked by the per-chunk relay container; Hail must already be initialized.
    Steps: look up this chunk's balanced sub-intervals (``--test-region``
    regions, or by index in the precompute JSON); read the VDS with those
    intervals, each contig's leading edge widened (``_expand_leading_edges``);
    read the sites HT with the same un-widened intervals (co-partitioned, so
    the densify join is shuffle-free and sites stay disjoint); run
    ``compute_all_release_stats_per_ref_site``; write.

    :param args: Parsed CLI args.
    :return: None.
    """
    project = args.project_name
    environment = args.environment
    start, stop = args.chunk_start, args.chunk_stop
    partition_range = list(range(start, stop))
    n = stop - start

    test = args.test
    chrom = args.chrom
    n_sub = max(args.read_subintervals_per_chunk, 1)
    results_environment = _results_environment(
        environment, test, args.results_environment
    )

    group_membership_ht_path = _group_membership_ht_path(project, environment, test)

    sub_intervals: list[hl.utils.Interval] | None = None
    # sub_intervals with leading edges widened (see _expand_leading_edges).
    vds_read_intervals: list[hl.utils.Interval] | None = None
    # --test-region reads the VDS via filter_intervals instead; holds those.
    vds_filter_intervals: list[hl.utils.Interval] | None = None
    # Layout hash namespacing the output (see _test_region_hash for --test-region).
    intervals_hash: str
    if args.test_region:
        # Balance the region into n_sub sub-intervals with the same machinery as
        # the prod precompute: one partition per region OOMs a single worker at
        # AoU's ~365k samples.
        intervals_hash = _test_region_hash(args)
        region_intervals = [_parse_region_interval(r) for r in args.test_region]
        if n_sub > 1:
            vds_probe = _probe_vds(
                project, environment, None, chrom, filter_intervals=region_intervals
            )
            sub_intervals = repartition_for_join(
                vds_probe.reference_data.rows(),
                n_partitions=n_sub,
                locus_intervals=True,
            )
            # The bounds span the probe's first to last reference-block key;
            # region sites outside them must still be read (AN=0), not dropped.
            sub_intervals = _snap_edges_to_bounds(
                sub_intervals, _region_bounds(region_intervals)
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
            # n_sub == 1: read the region as one partition per interval.
            sub_intervals = region_intervals
            vds_filter_intervals = region_intervals
    else:
        # Look this chunk up by index in the precompute JSON, driver-side (no
        # QoB job, no ~400s VDS re-open). --chrom is handled by the
        # orchestrator, so the worker just processes chunks[start].
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
        vds_read_intervals = _expand_leading_edges(sub_intervals, max_ref_block_len)

    # Auto-derive for manual single-chunk runs (the fan-out always passes
    # --chunk-output); deferred to here because the layout hash namespaces it.
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
        sites_path=_vep_context_sites_path(test, args.test_region),
        sub_intervals=sub_intervals,
    )

    # No chunk crosses a contig (_build_chunk_intervals; --test-region is
    # per-contig), so an autosomal chunk skips the sex-ploidy adjustment.
    if not _spans_sex_chromosome(sub_intervals):
        logger.info("Autosomal chunk: skipping the sex-karyotype ploidy adjustment.")
        sex_karyotype_field = None
    cov_and_an_ht = compute_all_release_stats_per_ref_site(
        vds,
        ref_ht,
        sex_karyotype_field=sex_karyotype_field,
        project=project,
        group_membership_ht=hl.read_table(group_membership_ht_path),
    )
    # Provenance: the merge inherits globals from its first input, so the final
    # HT records the layout that produced it.
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

    Globals are identical across chunks, so the union inherits them from the
    first input.

    :param input_paths: HT paths to union.
    :param output_path: Destination HT path.
    :param coalesce_to: If set, ``naive_coalesce`` to this many partitions
        before writing (group merges use len(inputs); the final merge uses
        ``--n-partitions``).
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

    Prefers the ``GNOMAD_QC_COMMIT`` env var (set by ``--submit-orchestrator``:
    the in-job checkout is a GitHub tarball with no ``.git``); else
    ``git rev-parse HEAD``.

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
    Build the relay shell prefix: repo checkouts + Hail config + version pin.

    Both repos are pulled at runtime (the image only provides hail + system
    deps). Also writes the Hail config.ini (Batch billing project,
    remote_tmpdir, requester-pays project) and patches ``/gsa-key/key.json``
    with a ``quota_project_id`` so requester-pays reads succeed from the QoB
    driver pod (Hail's own propagation doesn't reach the driver's Java GCS
    client; works from a laptop only because gcloud supplies it there).

    :param commit: gnomad_qc commit to pin.
    :param gcp_billing_project: Requester-pays project; patched into the GSA key.
    :param methods_branch: Branch/commit of gnomad_methods to pull.
    :return: Shell command string (terminated with newline).
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = (
        "https://github.com/broadinstitute/gnomad_methods/archive/"
        f"{methods_branch}.tar.gz"
    )
    methods_dir_suffix = methods_branch.replace("/", "-")
    # Hail config for hl.init(backend="batch"); write both the XDG path and
    # the legacy ~/.hail path.
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
        # TODO: drop this GSA-key patch once Hail propagates
        # gcs_requester_pays_configuration to the QoB driver pod.
        f"python3 -c \"import json, os; p='/gsa-key/key.json';"
        f" d=json.load(open(p)); d['quota_project_id']='{gcp_billing_project}';"
        f" json.dump(d, open(p+'.new','w')); os.replace(p+'.new', p)\"\n"
        # Pin the pipeline's Hail version (the relay's Python version sets the
        # QoB JAR, so this pins everything). WARNING: also a floor for every HT
        # this pipeline reads -- Hail's table format is not backward
        # compatible, so lowering the pin without rewriting the input HTs
        # breaks the fan-out (JSON artifacts are immune).
        "/opt/venv/bin/pip install --quiet --upgrade --force-reinstall"
        " --no-deps hail==0.2.128\n"
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_dir_suffix} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def _build_relay_common_flags(args: argparse.Namespace, *, chunk: bool) -> str:
    """
    Build the CLI flag string shared by per-chunk / per-merge relay invocations.

    ``chunk=True`` appends the read/compute flags a chunk relay needs; the
    merge relay (which only unions HTs) omits them. Per-job flags are added by
    the submit helpers.

    :param args: Parsed CLI args.
    :param chunk: Include the chunk-only read/compute flags.
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
        # --experimental: the inner QoB driver attaches to this outer batch.
        flags.append("--experimental")
    if args.app_name:
        flags.append(f"--app-name {args.app_name}")
    if args.cov_and_an_output_suffix:
        flags.append(f"--cov-and-an-output-suffix {args.cov_and_an_output_suffix}")
    if args.results_environment:
        # Keep worker fallback path resolution consistent with the orchestrator.
        flags.append(f"--results-environment {args.results_environment}")
    # QoB driver/worker sizing; merge workers reuse the chunk sizing.
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
        if args.test_sample_subset:
            flags.append("--test-sample-subset")
        if args.chrom:
            flags.append(f"--chrom {args.chrom}")
        if args.test_region:
            flags.append(f"--test-region {' '.join(args.test_region)}")
        # Forward the normalized --test so the worker resolves the same
        # test-scoped paths.
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

    Shared by chunk and merge submits. Each relay is a non-spot coordinator
    (preemption mid-wait would orphan its inner QoB batch), pinned to
    ``BATCH_REGIONS``, with ``--chunk-attempts`` attempts. No-ops on empty
    ``job_specs``.

    :param args: Parsed CLI args.
    :param backend_kwargs: kwargs for ``hb.ServiceBackend``.
    :param batch_name: Hail Batch name.
    :param job_specs: Per-job config.
    :param log_label: Noun for log messages ("chunk" / "merge").
    :return: Batch id, or None (empty specs or ``--batch-dry-run``).
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
            # Non-spot: preemption mid-wait would orphan the inner QoB batch.
            j.spot(False)
            # Default 1 attempt: a Batch-level retry cannot cancel the orphaned
            # inner QoB batch and would race it. See --chunk-attempts.
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
        return getattr(submitted, "id", None)
    finally:
        backend.close()


def _submit_orchestrator_batch(args: argparse.Namespace) -> None:
    """
    Submit THIS orchestrator invocation as one small non-spot Hail Batch job.

    Re-invokes the same command line (minus ``--submit-orchestrator``) inside a
    Batch job and returns immediately (``batch.run(wait=False)``), so a full
    fan-out doesn't tie up a laptop. The orchestrator never initializes Hail,
    so the job needs no JVM sizing. The commit is resolved here (where
    ``.git`` exists) and forwarded via ``GNOMAD_QC_COMMIT``. If the job dies,
    resubmit: completed chunks are skipped via ``_SUCCESS``.

    :param args: Parsed CLI args; forwarded verbatim to the wrapped run.
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
        # Non-spot: losing the coordinator mid-campaign stops wave dispatch.
        j.spot(False)
        # One attempt: a retry would dispatch waves concurrently with the
        # first's; recovery is manual resubmission.
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

    Each chunk job is a relay container running ``--run-chunk``.
    ``chunk_indices`` is the already-filtered pending set.

    :param args: Parsed CLI args.
    :param backend_kwargs: kwargs for ``hb.ServiceBackend``.
    :param chunk_indices: Pending chunk indices to submit.
    :param cov_and_an_ht_path: Resolved canonical output HT path.
    :param intervals_hash: Layout hash namespacing each chunk's output.
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags.
    :param script: Script path inside the relay container.
    :param wave_label: Optional batch-name suffix (e.g. ``"w003of049"``).
    :return: Batch id of the submitted wave (None under ``--batch-dry-run``).
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
        # Chunk identity is the index: the worker looks itself up in the JSON
        # by --chunk-start (= idx); --chunk-stop is idx+1.
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

    Shared by the chunk orchestrator and the merge so they always agree. A
    ``--test-region`` run is one chunk; otherwise chunks come from the
    precompute JSON, read driver-side.

    :param args: Parsed CLI args.
    :return: ``(chunk_contigs, eligible, intervals_hash)`` -- per-chunk contig
        (None for the single ``--test-region`` chunk), eligible indices after
        the ``--chrom`` filter, and the layout hash (:func:`_test_region_hash`
        when there is no JSON).
    """
    if args.test_region:
        chunk_contigs: list[str | None] = [None]
        intervals_hash = _test_region_hash(args)
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
    Fan coverage/AN compute out as relay chunk jobs in Hail Batch (one per chunk).

    Pending chunks are split into sequential waves of ``--wave-size``; each
    wave is one Hail Batch run to completion before the next. Sequential waves
    bound concurrency (relays + their nested QoB drivers) and avoid Hail
    Batch's process-global progress display, which crashes on concurrent
    ``batch.run()`` calls. The binding capacity is the shared Hail Batch
    service's worker pool, not gnomad-production's own GCP quota.

    Idempotent: chunks with a ``_SUCCESS`` are skipped, so rerunning resumes.
    ``batch.run()`` does not raise on per-job failure; after the last wave,
    missing chunks are re-dispatched for up to ``--fanout-retry-passes``
    passes, and whatever remains is recorded in ``_failed_chunks.json``. Run
    ``--merge-cov-chunks`` after this finishes.

    :param args: Parsed CLI args.
    :param cov_and_an_ht_path: Canonical output HT path; each chunk writes a
        sibling ``_chunks/<hash>/<idx>.chunk.ht``.
    :return: None.
    """
    project = args.project_name

    # Fail fast if the sites HT the chunks read isn't there yet.
    sites_path = _vep_context_sites_path(args.test, args.test_region)
    if not file_exists(sites_path):
        raise FileNotFoundError(
            f"vep_context sites HT not found at {sites_path}. Run"
            " --write-vep-context-sites first (it is a prerequisite for the"
            " chunk compute)."
        )

    chunk_contigs, eligible, intervals_hash = _eligible_chunk_indices(args)
    n_chunks = len(chunk_contigs)

    if args.overwrite:
        pending_indices = list(eligible)
    else:
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
        # batch.run() does not raise on per-job failure; re-check the outputs.
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
        # Rewrite after every wave so the list survives an orchestrator death.
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

    # Automatic retry passes. Safe where Batch-level attempts are not: by the
    # time a pass runs, every relay from the previous dispatch has exited, so a
    # re-dispatch cannot race a live relay's inner QoB batch. Chunks missing
    # after the last pass stay in _failed_chunks.json.
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

    Each job runs ``--run-merge`` and writes a per-group HT (``_group_path``).
    Only intermediate levels go through here; the final union is submitted
    separately. Per-group coalesce target is the group's input count.

    :param args: Parsed CLI args.
    :param backend_kwargs: kwargs for ``hb.ServiceBackend``.
    :param group_indices: Pending group indices to submit.
    :param groups: Input HT paths per group index.
    :param group_output_paths: Output HT path per group index.
    :param setup_cmd: Shell prefix from ``_build_setup_command``.
    :param common_flags_str: Shared CLI flags.
    :param script: Script path inside the relay container.
    :param level: Merge-tree level (1-indexed).
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

    Runs after the fan-out; submits Batch jobs and never initializes Hail.
    Chunks are enumerated with the same helper as the fan-out; missing chunks
    fail loudly. Each level groups its inputs into windows of
    ``--merge-group-size`` and emits one ``--run-merge`` job per group, until
    one final union writes the canonical path. Levels are sequential; safe to
    re-run (existing group HTs and final HT are skipped without
    ``--overwrite``).

    :param args: Parsed CLI args.
    :param cov_and_an_ht_path: Canonical output cov_and_an HT path.
    :return: None.
    """
    project = args.project_name

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

    # Level shape, to log the full plan upfront.
    shape = [len(eligible)]
    while shape[-1] > gs:
        shape.append((shape[-1] + gs - 1) // gs)
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

    # Intermediate levels: iterate while #inputs > gs.
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

    # Final merge: one job unions the remaining (<= gs) inputs.
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

    Dispatches the three mutually-exclusive roles (see the module docstring):
    orchestrator (submit and return), worker (one chunk/merge unit), or the
    in-process pipeline (setup -> compute -> assemble).
    """
    project = args.project_name
    environment = args.environment
    # args.test is normalized at parse time to fold in --test-n-partitions /
    # --test-region.
    test = args.test
    chrom = args.chrom
    overwrite = args.overwrite
    results_environment = _results_environment(
        environment, test, args.results_environment
    )

    # ===================================================================
    # ROLE 1: ORCHESTRATOR — submit a Hail Batch of relay jobs and return.
    # Never initializes Hail (existence checks use hailtop.fs).
    #   --use-batch-fanout : scatter the per-chunk compute.
    #   --merge-cov-chunks : tree-reduce the chunk HTs.
    #   --submit-orchestrator : run either of the above inside a Batch job.
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

    # QoB / dataproc init — ROLE 2 (worker) and ROLE 3 (in-process) only.
    #
    # Hail 0.2.139 workaround: the hl.experimental.densify scan crashes at
    # compile time (NotImplementedError, SimpleSStream.settableTupleTypes) when
    # the scanned table has more partitions than the branching factor (default
    # 50); 0.2.137/0.2.138 are unaffected. A branching factor above the chunk's
    # partition count (chunks read at most --read-subintervals-per-chunk) keeps
    # the scan combine single-level, avoiding the broken lowering. Chunk
    # workers only: the single-job path runs on dataproc (0.2.137). A
    # --test-region chunk with one sub-interval reads the region's native VDS
    # partitions (filter_intervals), so it gets a fixed high floor instead.
    branching_factor = None
    if args.run_chunk:
        branching_factor = max(
            50, 2 * args.read_subintervals_per_chunk, 1024 if args.test_region else 0
        )

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
    # ROLE 2: WORKER — run exactly one chunk/merge unit, then return.
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
    # ROLE 3: IN-PROCESS PIPELINE — the strict single-job compute (gnomAD
    # consent-drop, dev/test AoU), setup writers, and release assembly.
    # ===================================================================
    try:
        cov_and_an_ht_path = _resolve_cov_and_an_ht_path(
            project,
            results_environment,
            test=test,
            suffix=args.cov_and_an_output_suffix,
            chrom=args.chrom,
        )

        # Union the per-contig HTs (written by --chrom merges) into the
        # canonical HT; validate afterward with --validate-cov-and-an.
        if args.assemble_chrom_coverage:
            with hfs.open(_chunk_intervals_path(environment, test)) as f:
                contigs = sorted({c["contig"] for c in json.load(f)["chunks"]})
            inputs = [_apply_path_suffix(cov_and_an_ht_path, c) for c in contigs]
            _run_coverage_merge(
                inputs, cov_and_an_ht_path, coalesce_to=args.n_partitions
            )
            return

        group_membership_ht_path = _group_membership_ht_path(project, environment, test)

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

        if args.write_vep_context_sites:
            logger.info("Writing preprocessed vep_context sites HT...")
            sites_path = _vep_context_sites_path(test, args.test_region)
            check_resource_existence(
                output_step_resources={"vep_context_sites": [sites_path]},
                overwrite=overwrite,
            )
            site_intervals = None
            if test and args.test_region:
                # Scope the sites directly to the regions (no VDS probe needed).
                site_intervals = [_parse_region_interval(r) for r in args.test_region]
                logger.info(
                    "Test sites scoped to %d explicit --test-region interval(s): %s",
                    len(site_intervals),
                    args.test_region,
                )
            elif test:
                # Scope the test sites to the loci the test compute reads: the
                # chunk-intervals JSON's sub-intervals when present, else the
                # first N native VDS partitions (strict single-job path).
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
                )
            # Already coalesced to GROUP_MEMBERSHIP_N_PARTITIONS in
            # get_group_membership_ht; writing materializes it for every consumer.
            group_membership_ht.write(group_membership_ht_path, overwrite=overwrite)

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

            # Optionally co-partition the VDS and vep_context reads;
            # --test-region overrides and reads via filter_intervals.
            strict_sub_intervals = None
            strict_filter_intervals = None
            if args.test_region:
                # Explicit region: filter_intervals read; skip partition scoping.
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

            # --test-region reads via filter_intervals; read_intervals is dropped.
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
                sites_path=_vep_context_sites_path(test, args.test_region),
                sub_intervals=strict_sub_intervals,
            )
            # None (whole VDS) keeps the adjustment; a per-contig or --test-region
            # span on autosomes skips it.
            if strict_sub_intervals and not _spans_sex_chromosome(strict_sub_intervals):
                logger.info(
                    "Autosomal span: skipping the sex-karyotype ploidy adjustment."
                )
                sex_karyotype_field = None

            cov_and_an_ht = compute_all_release_stats_per_ref_site(
                vds,
                ref_ht,
                sex_karyotype_field=sex_karyotype_field,
                project=project,
                group_membership_ht=hl.read_table(group_membership_ht_path),
            )
            if args.n_partitions is not None:
                # Temp checkpoint only to stabilize before the coalesce.
                cov_and_an_ht = cov_and_an_ht.checkpoint(
                    new_temp_file(f"{project}_cov_and_an", "ht")
                ).naive_coalesce(args.n_partitions)
            cov_and_an_ht.write(cov_and_an_ht_path, overwrite=overwrite)

        if args.validate_cov_and_an:
            logger.info("Validating cov_and_an HT covers all vep_context sites...")
            sites_ht = hl.read_table(_vep_context_sites_path(test, args.test_region))
            # Scope validation to what actually ran: the JSON's chunks for a
            # --test-n-partitions run, whole contigs for --chrom.
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
                # The compute emits AN=0 rows for uncovered sites, so a missing
                # site is a real loss; the error text explains the diagnosis.
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

        # --- ASSEMBLE (gnomAD): v4 release minus consent-drop -> gnomAD v5. ---
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
            # The subtraction rebuilds DP sums as mean * n, so n must be the
            # cohort the HT was actually computed on, not the script constant.
            consent_drop_count = hl.eval(gnomad_ht.coverage_stats_meta_sample_count)
            if consent_drop_count != GNOMAD_CONSENT_DROP_SAMPLE_COUNT:
                logger.warning(
                    "Consent-drop coverage HT has %d samples; script constant is %d.",
                    consent_drop_count,
                    GNOMAD_CONSENT_DROP_SAMPLE_COUNT,
                )
            ht = merge_gnomad_coverage_hts(
                gnomad_ht, gnomad_release_ht, consent_drop_count=consent_drop_count
            )
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

        # --- ASSEMBLE (AoU): release exports + qual-hists merge. ---
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
            # Coverage is gnomAD-only (AoU computes no coverage stats), so the
            # release HT is the gnomAD consent-drop result as-is.
            ht = hl.read_table(gnomad_coverage_ht_path)
            if args.n_partitions is not None:
                # Temp checkpoint only to stabilize before the coalesce.
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
                # Temp checkpoint only to stabilize before the coalesce.
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
            # Drop age hists (the frequency script's job) and rekey by locus;
            # distinct() dedups split rows (the _all hists are locus-level).
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
            "Override the bucket the release artifacts write to, independent of"
            " the compute --environment. Default: prod writes to 'dataproc', a"
            " --test run stays in the compute environment."
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
            "Route QoB init through hl.experimental.init and attach the driver"
            " to an existing Hail Batch (HAIL_BATCH_ID; raises if unset)."
            " Without it, each invocation creates its own Batch."
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
            "Max JVM heap (-Xmx) for the in-process QoB driver under"
            " --experimental, e.g. '5g'. Set to ~50-70%% of container memory"
            " (the rest is off-heap/Python/OS). Rejected without --experimental."
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
            "Memory type for worker nodes. Default None = QoB's 'standard',"
            " already generous (272 MiB peak measured). 'lowmem' is rejected"
            " for JVM jobs."
        ),
    )
    parser.add_argument(
        "--n-partitions",
        help=(
            "Partitions for the output Table. Default None: no final" " naive_coalesce."
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--cov-and-an-output-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix before the .ht extension (A/B comparison)."
            " Applies to writes and downstream reads; pass consistently."
        ),
    )
    parser.add_argument(
        "--qual-hists-output-suffix",
        type=str,
        default=None,
        help=("Like --cov-and-an-output-suffix, for the --merge-qual-hists output."),
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
            "With --test on AoU, subsample to ~0.1%% of samples. Default off:"
            " a --test run keeps all samples so AN values stay comparable."
        ),
    )
    test_group.add_argument(
        "--test-n-partitions",
        type=int,
        default=None,
        help=(
            "Cheap test over the first N VDS partitions. With"
            " --write-chunk-intervals, scopes the probe (-> ceil(N/ppc) chunks;"
            " combine with --chrom). With the strict compute, reads the first N"
            " native partitions (genome start; does not compose with --chrom)."
            " Default None: full run."
        ),
    )
    test_group.add_argument(
        "--chrom",
        default=None,
        help=(
            "Filter to a single contig (e.g. chr22). Appended to the output"
            " suffix so per-contig runs never collide; assemble with"
            " --assemble-chrom-coverage."
        ),
    )
    test_group.add_argument(
        "--test-region",
        nargs="+",
        default=None,
        help=(
            "Scope a test to explicit locus intervals (e.g."
            " chr1:55058666-55108666) instead of the first N VDS partitions;"
            " auto-enables test mode. Sites and VDS reads are region-scoped;"
            " straddling ref blocks are split (strict path / n_sub=1) or"
            " handled by leading-edge widening (fan-out). Mutually exclusive"
            " with --test-n-partitions."
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
            "Preprocess vep_context once into the deduped, stripped, locus-keyed"
            " sites HT every chunk reads. Prerequisite for the compute. With"
            " --test, writes a scoped _test table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--write-chunk-intervals",
        help=(
            "Precompute every chunk's balanced read sub-intervals in one VDS"
            " open (JSON, 30-day storage). REQUIRED before the fan-out. Uses"
            " --total-partitions / --partitions-per-chunk /"
            " --read-subintervals-per-chunk; --test writes a _test file."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--validate-cov-and-an",
        help=(
            "Anti-join the sites HT against the merged cov_and_an HT; fail if"
            " any site is missing. Run after --merge-cov-chunks."
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
            "Strict single-job compute only: co-partition the VDS and"
            " vep_context reads (shuffle-free join). Flag alone = the"
            " 1.1x-native default (recommended for full runs); a value = an"
            " explicit count. Honors --chrom."
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
            "Fan the coverage/AN compute out as one Batch job per chunk."
            " Mutually exclusive with the other entry points; run"
            " --merge-cov-chunks afterward."
        ),
    )
    fanout_group.add_argument(
        "--submit-orchestrator",
        action="store_true",
        help=(
            "Run the orchestrator inside Hail Batch: submit one small non-spot"
            " job that re-invokes this command (minus this flag) and return."
            " Requires --use-batch-fanout or --merge-cov-chunks; resubmit the"
            " same command to resume if the job dies."
        ),
    )
    fanout_group.add_argument(
        "--fanout-retry-passes",
        type=int,
        default=1,
        help=(
            "Extra automatic re-dispatch passes for still-missing chunks after"
            " the last wave. Unlike --chunk-attempts, a pass starts only after"
            " every prior relay has exited, so it cannot race an orphaned inner"
            " batch. 0 disables. Default 1."
        ),
    )
    fanout_group.add_argument(
        "--merge-cov-chunks",
        action="store_true",
        help=(
            "Submit the recursive tree-reduce merge of the fan-out's chunk HTs"
            " (group size --merge-group-size) down to the canonical HT."
            " Missing chunks fail loudly."
        ),
    )
    fanout_group.add_argument(
        "--assemble-chrom-coverage",
        action="store_true",
        help=(
            "Union the per-contig HTs from --merge-cov-chunks --chrom runs into"
            " the canonical HT (contigs from the chunk-intervals JSON). Runs"
            " in-process."
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
            " ~$0.068/partition). IMPORTANT: keep identical across"
            " --write-chunk-intervals / --use-batch-fanout / --merge-cov-chunks"
            " (the JSON fixes the chunk set)."
        ),
    )
    fanout_group.add_argument(
        "--read-subintervals-per-chunk",
        type=int,
        default=50,
        help=(
            "Balanced sub-intervals per chunk (computed at"
            " --write-chunk-intervals time); the VDS and vep_context are both"
            " read with these intervals -- co-partitioned, no shuffle."
            " Default 50."
        ),
    )
    fanout_group.add_argument(
        "--total-partitions",
        type=int,
        default=145192,
        help=(
            "Total VDS partitions to scatter (full prod extent). Default 145192"
            " (prod AoU VDS). Consumed by --write-chunk-intervals;"
            " --test-n-partitions overrides for tests."
        ),
    )
    fanout_group.add_argument(
        "--wave-size",
        type=int,
        default=1000,
        help=(
            "Max chunks per Hail Batch under --use-batch-fanout; waves run"
            " sequentially, bounding concurrent relays and their QoB drivers."
            " Capacity is the shared Hail Batch service pool, not"
            " gnomad-production GCP quota. <= 0 for a single batch."
            " Default 1000."
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
            "GCS scratch for the orchestrator's ServiceBackend. Default None:"
            " host Hail config."
        ),
    )
    fanout_group.add_argument(
        "--chunk-cpu",
        type=float,
        default=0.5,
        help=(
            "CPU per chunk relay container. Default 0.5: the relay only builds"
            " IR and waits on its inner QoB batch."
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
            "Memory profile for the QoB driver pod. Default 'highmem': on"
            " 'standard' the driver died lowering the big fused execute"
            " ('unexpected end of stream in jvm' -- a JVM-total limit, not a"
            " RegionPool OOM). 'lowmem' is rejected for JVM jobs."
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
            "Memory per QoB worker pod. Default None = QoB's 'standard',"
            " already generous (272 MiB peak measured). 'lowmem' is rejected"
            " for JVM jobs."
        ),
    )
    fanout_group.add_argument(
        "--chunk-attempts",
        type=int,
        default=1,
        help=(
            "Max attempts per relay job. Default 1 = no auto-retry: Batch"
            " cannot cancel a lost attempt's inner QoB batch, so a retry would"
            " race it. Failures land in _failed_chunks.json; rerun the fan-out."
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
            "Output path for --run-chunk's HT. When unset, auto-derived"
            " (the orchestrator always passes it explicitly)."
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
            "naive_coalesce the unioned HT to this many partitions before"
            " writing (set by the merge orchestrator)."
        ),
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    # Fold --test-n-partitions / --test-region into the canonical args.test.
    args.test = (
        args.test or args.test_n_partitions is not None or args.test_region is not None
    )

    # Every batch-only argument lives in this one list; "provided" is detected
    # by comparing against the parser default.
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
        # Fan-out / merge orchestration args.
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

    if args.submit_orchestrator and not (
        args.use_batch_fanout or args.merge_cov_chunks
    ):
        parser.error(
            "--submit-orchestrator requires --use-batch-fanout or"
            " --merge-cov-chunks (it wraps an orchestrator run in a Batch job)."
        )

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

    # 'lowmem' is excluded via choices on every memory flag: Hail Batch rejects
    # it for JVM jobs, and every QoB pod this script starts is a JVM job.

    # Single-action entry points: each dispatches one unit and returns; pairing
    # with another entry point or an in-process step would silently skip work.
    entry_points = {
        "--use-batch-fanout": args.use_batch_fanout,
        "--merge-cov-chunks": args.merge_cov_chunks,
        "--assemble-chrom-coverage": args.assemble_chrom_coverage,
        "--run-chunk": args.run_chunk,
        "--run-merge": args.run_merge,
    }
    step_flags = {
        "--write-aou-downsampling-ht": args.write_aou_downsampling_ht,
        "--write-group-membership-ht": args.write_group_membership_ht,
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
    if args.run_merge:
        if args.merge_output is None or not args.merge_inputs:
            parser.error("--run-merge requires --merge-output and --merge-inputs.")

    main(args)
