"""
Hail Batch relay fan-out shared by ``compute_coverage.py`` and ``generate_frequency.py``.

A Query-on-Batch run has one driver, and a driver crash loses everything it was
computing. Both scripts therefore cut the genome into chunks and run each chunk from
its own small, non-spot Batch job (a "relay") that starts a QoB driver inside its
container. A crash then costs one chunk, and a rerun skips the chunks that already
have a ``_SUCCESS`` marker. The orchestrator that submits the relays never starts
Hail.

Chunk outputs are stored under a hash of the chunk layout. The layout is sampled
without a fixed seed, so each regeneration moves the cut points, and a chunk left
behind by an older layout must not be taken for a finished chunk of the new one.
"""

import hashlib
import json
import logging
import os
import re
import subprocess
from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, NamedTuple, Optional, Sequence, Set, Tuple

import hail as hl
import hailtop.batch as hb
import hailtop.fs as hfs
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v5.resources.constants import BATCH_REGIONS, BATCH_TMP_BUCKET

logger = logging.getLogger("v5_batch_fanout")

BATCH_REMOTE_TMPDIR = f"gs://{BATCH_TMP_BUCKET}"


class RelayJobSpec(NamedTuple):
    """One relay job for :func:`submit_relay_batch`."""

    name: str
    cpu: float
    memory: str
    storage: str
    attempts: int
    command: str


# ---------------------------------------------------------------------------
# Relay containers
# ---------------------------------------------------------------------------


def resolve_commit() -> str:
    """
    Return the gnomad_qc commit the relay containers should check out.

    The ``GNOMAD_QC_COMMIT`` environment variable takes precedence so an orchestrator
    that itself runs inside a Batch job, where the checkout is a tarball with no
    ``.git`` directory, can still pin the relays.

    :return: Full commit hash.
    """
    return os.getenv("GNOMAD_QC_COMMIT") or (
        subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    )


def build_setup_command(
    commit: str,
    gcp_billing_project: str = "broad-mpg-gnomad",
    methods_branch: str = "main",
    hail_version: Optional[str] = None,
) -> str:
    """
    Return the shell prefix every relay runs before the script.

    Both repos change faster than the image, so the relay downloads them at the
    pinned commit and branch when it starts. The Hail config file sets the Batch
    billing project and scratch directory for the nested QoB driver, and together
    with the ``quota_project_id`` patch gives Hail's GCS client a requester-pays
    project, which a bare container has no gcloud config to supply.

    :param commit: gnomad_qc commit to check out.
    :param gcp_billing_project: Requester-pays project.
    :param methods_branch: gnomad_methods branch or commit to check out.
    :param hail_version: If set, reinstall this Hail version in the container;
        otherwise the image's Hail is used. The relay's Python Hail sets the QoB
        JAR, so this pins everything. WARNING: the pin is also a floor for every
        HT the pipeline reads. Hail's table format is not backward compatible,
        so lowering the pin without rewriting the input HTs breaks the fan-out.
    :return: Shell command string ending in a newline.
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = (
        "https://github.com/broadinstitute/gnomad_methods/archive/"
        f"{methods_branch}.tar.gz"
    )
    methods_dir_suffix = methods_branch.replace("/", "-")
    config_body = (
        "[batch]\n"
        "billing_project = gnomad-production\n"
        f"remote_tmpdir = {BATCH_REMOTE_TMPDIR}\n"
        "[gcs_requester_pays]\n"
        f"project = {gcp_billing_project}\n"
    )
    pin = (
        "/opt/venv/bin/pip install --quiet --upgrade --force-reinstall --no-deps"
        f" hail=={hail_version}\n"
        if hail_version
        else ""
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
        f"{pin}"
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_dir_suffix} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def relay_context(
    commit: str,
    methods_branch: str,
    batch_billing_project: str,
    batch_remote_tmpdir: Optional[str],
    script_name: str,
    gcp_billing_project: str = "broad-mpg-gnomad",
    hail_version: Optional[str] = None,
) -> Tuple[str, Dict[str, str], str]:
    """
    Return what every relay submission needs: setup command, backend kwargs, script.

    Built once per orchestrator so the fan-out and the merge check out the same
    commits.

    :param commit: gnomad_qc commit the relays check out.
    :param methods_branch: gnomad_methods branch or commit the relays check out.
    :param batch_billing_project: Hail Batch billing project for the relay jobs.
    :param batch_remote_tmpdir: ``hb.ServiceBackend`` scratch, or None for Hail's
        configured one.
    :param script_name: Script file name under ``gnomad_qc/v5/annotations``.
    :param gcp_billing_project: Requester-pays project for the relays.
    :param hail_version: Optional Hail pin for :func:`build_setup_command`.
    :return: ``(setup_cmd, backend_kwargs, script_command)``.
    """
    setup_cmd = build_setup_command(
        commit, gcp_billing_project, methods_branch, hail_version
    )
    backend_kwargs = {"billing_project": batch_billing_project}
    if batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = batch_remote_tmpdir
    script = f"python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/{script_name}"
    return setup_cmd, backend_kwargs, script


def submit_relay_batch(
    batch_image: str,
    backend_kwargs: Dict[str, str],
    batch_name: str,
    job_specs: Sequence[RelayJobSpec],
    log_label: str,
    dry_run: bool = False,
) -> Optional[int]:
    """
    Submit one Hail Batch of relay jobs and wait for it.

    Relays are non-spot because a preempted relay leaves the QoB batch it was
    waiting on running with nobody to collect it. For the same reason callers
    normally ask for one attempt: a Batch retry cannot cancel that orphaned batch
    and would run alongside it.

    :param batch_image: Docker image for the relay containers.
    :param backend_kwargs: kwargs for ``hb.ServiceBackend``.
    :param batch_name: Hail Batch name.
    :param job_specs: Jobs to submit; nothing is submitted when empty.
    :param log_label: Noun for log messages ("chunk", "merge").
    :param dry_run: Validate the batch without running it.
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
            j.image(batch_image)
            j.cpu(spec.cpu)
            j.memory(spec.memory)
            j.storage(spec.storage)
            j.regions(BATCH_REGIONS)
            j.spot(False)
            j.n_max_attempts(spec.attempts)
            j.command(spec.command)
        logger.info(
            "Submitting Hail Batch '%s': %d %s jobs (dry_run=%s)",
            batch_name,
            len(job_specs),
            log_label,
            dry_run,
        )
        submitted = batch.run(dry_run=dry_run)
        return getattr(submitted, "id", None)
    finally:
        backend.close()


# ---------------------------------------------------------------------------
# Chunk layout and intervals
# ---------------------------------------------------------------------------


def chunk_intervals_hash(data: Dict[str, Any]) -> str:
    """
    Return a 16-hex-char content hash of a chunk layout.

    Every regeneration of the layout moves the cut points, and keying outputs by
    this hash keeps each layout's chunks apart.

    :param data: Parsed chunk-intervals JSON.
    :return: First 16 hex chars of the SHA-256 of its canonical serialization.
    """
    payload = {k: v for k, v in data.items() if k != "intervals_hash"}
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:16]


def test_region_hash(test_region: Sequence[str]) -> str:
    """
    Return the layout hash of a ``--test-region`` run, which has no layout JSON.

    Hashing the region strings keeps two region tests under the same output path
    from seeing each other's chunk as already present.

    :param test_region: Region strings as passed on the command line.
    :return: ``test_region_`` plus a 16-hex-char hash.
    """
    return f"test_region_{chunk_intervals_hash({'test_region': list(test_region)})}"


def interval_to_list(iv: hl.utils.Interval) -> list:
    """
    Serialize a locus interval for the chunk-intervals JSON.

    :param iv: Locus interval.
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


def interval_from_list(t: Sequence, reference_genome: str) -> hl.utils.Interval:
    """
    Rebuild a locus interval from its chunk-intervals JSON form.

    :param t: Output of :func:`interval_to_list`.
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


def parse_region_interval(
    interval_s: str, reference_genome: str = "GRCh38"
) -> hl.utils.Interval:
    """
    Parse ``contig:start-end`` into a half-open locus interval.

    Half-open so adjacent regions never share a locus.

    :param interval_s: Interval string; commas in positions are allowed.
    :param reference_genome: Reference genome name.
    :return: ``[start, end)`` locus interval.
    """
    contig, span = interval_s.split(":")
    start_pos, end_pos = (int(p.replace(",", "")) for p in span.split("-"))
    return hl.Interval(
        hl.Locus(contig, start_pos, reference_genome=reference_genome),
        hl.Locus(contig, end_pos, reference_genome=reference_genome),
        includes_start=True,
        includes_end=False,
    )


def spans_sex_chromosome(
    intervals: Optional[Sequence[hl.utils.Interval]], chrom: Optional[str] = None
) -> bool:
    """
    Return whether a read scope can include chrX or chrY.

    The sex-ploidy adjustment changes nothing on autosomes, so a scope that cannot
    touch chrX or chrY may skip it. A read with no scope counts as spanning.

    :param intervals: Read intervals, or None.
    :param chrom: Single contig, or None.
    :return: True unless the scope is known to exclude both sex contigs.
    """
    if intervals:
        rg = intervals[0].start.reference_genome
        sex = set(rg.x_contigs) | set(rg.y_contigs)
        order = {c: k for k, c in enumerate(rg.contigs)}
        contigs: Set[str] = set()
        for iv in intervals:
            contigs.update(
                rg.contigs[order[iv.start.contig] : order[iv.end.contig] + 1]
            )
        return bool(contigs & sex)
    if chrom:
        rg = hl.get_reference("GRCh38")
        return chrom in set(rg.x_contigs) | set(rg.y_contigs)
    return True


def select_chunks_for_contig(
    chunk_contigs: Sequence[Optional[str]], chrom: Optional[str]
) -> List[int]:
    """
    Return the chunk indices on ``chrom``, or every index when ``chrom`` is None.

    :param chunk_contigs: Contig of each chunk, by index.
    :param chrom: Contig to keep, or None.
    :return: Eligible chunk indices.
    """
    if not chrom:
        return list(range(len(chunk_contigs)))
    eligible = [i for i, c in enumerate(chunk_contigs) if c == chrom]
    if not eligible:
        raise ValueError(
            f"No chunks on --chrom {chrom}; the layout covers"
            f" {sorted(c for c in set(chunk_contigs) if c)}."
        )
    return eligible


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------


def apply_path_suffix(path: str, suffix: Optional[str]) -> str:
    """
    Insert ``_<suffix>`` before the ``.ht`` extension; unchanged when ``suffix`` is empty.

    :param path: HT path ending in ``.ht``.
    :param suffix: Suffix without the leading underscore.
    :return: Suffixed path.
    """
    if not suffix:
        return path
    return path.rstrip("/").removesuffix(".ht") + f"_{suffix}.ht"


def combine_suffix(*parts: Optional[str]) -> Optional[str]:
    """
    Join the non-empty parts with underscores, so each scope writes its own path.

    :param parts: Suffix fragments; empty ones are skipped.
    :return: Joined suffix, or None when every part is empty.
    """
    kept = [p for p in parts if p]
    return "_".join(kept) if kept else None


def _base(ht_path: str) -> str:
    return ht_path.rstrip("/").removesuffix(".ht")


def chunk_path(ht_path: str, idx: int, intervals_hash: str) -> str:
    """
    Return the per-chunk HT path beside the final HT.

    :param ht_path: Final HT path the chunks will be merged into.
    :param idx: Chunk index.
    :param intervals_hash: Layout hash namespacing the chunk directory.
    :return: ``<ht>_chunks/<hash>/<idx:08d>.chunk.ht``.
    """
    return f"{_base(ht_path)}_chunks/{intervals_hash}/{idx:08d}.chunk.ht"


def group_path(
    ht_path: str, level: int, group_idx: int, tree_tag: str, intervals_hash: str
) -> str:
    """
    Return a merge-tree intermediate HT path.

    Level, tree shape and layout hash are all in the directory so a rerun with a
    different tree or layout writes fresh instead of reusing stale groups.

    :param ht_path: Final HT path.
    :param level: Merge-tree level, 1-indexed.
    :param group_idx: Group index within the level.
    :param tree_tag: Tree-shape tag, e.g. ``gs500``.
    :param intervals_hash: Layout hash.
    :return: ``<ht>_merge_groups_<tag>/<hash>/L<level>_<group>.ht``.
    """
    return (
        f"{_base(ht_path)}_merge_groups_{tree_tag}/{intervals_hash}/"
        f"L{level:02d}_{group_idx:08d}.ht"
    )


def failed_chunks_path(ht_path: str, intervals_hash: str) -> str:
    """
    Return the failed-chunk manifest path for one layout.

    :param ht_path: Final HT path.
    :param intervals_hash: Layout hash.
    :return: ``<ht>_chunks/<hash>/_failed_chunks.json``.
    """
    return f"{_base(ht_path)}_chunks/{intervals_hash}/_failed_chunks.json"


def list_present_chunk_indices(ht_path: str, intervals_hash: str) -> Set[int]:
    """
    Return the chunk indices that have a ``_SUCCESS`` marker.

    One listing is cheaper than one GCS stat per chunk, and looking for ``_SUCCESS``
    rather than the chunk directory means a half-written chunk is rerun. A missing
    directory lists as empty.

    :param ht_path: Final HT path.
    :param intervals_hash: Layout hash.
    :return: Completed chunk indices.
    """
    present: Set[int] = set()
    for entry in hfs.ls(f"{_base(ht_path)}_chunks/{intervals_hash}/*/_SUCCESS"):
        m = re.search(r"/(\d+)\.chunk\.ht/_SUCCESS$", entry.path)
        if m:
            present.add(int(m.group(1)))
    return present


def new_run_id() -> str:
    """
    Return a UTC timestamp id for one orchestrator run.

    :return: ``run-YYYYmmddTHHMMSSZ``.
    """
    return "run-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def write_failed_chunks_manifest(
    ht_path: str,
    intervals_hash: str,
    failed: Sequence[int],
    n_dispatched: int,
    run_id: str,
    commit: str,
    app_name: Optional[str],
    waves: Sequence[Dict[str, Any]],
) -> Optional[str]:
    """
    Record which dispatched chunks did not land, or clear the record when all did.

    A rerun of the fan-out resumes from the missing chunks on its own. The manifest
    is a durable record of what failed and when, which the logs are not.

    :param ht_path: Final HT path.
    :param intervals_hash: Layout hash.
    :param failed: Dispatched chunk indices with no ``_SUCCESS``.
    :param n_dispatched: Chunks dispatched by this run.
    :param run_id: This orchestrator run's id.
    :param commit: gnomad_qc commit the relays ran.
    :param app_name: ``--app-name`` the relays used.
    :param waves: Per-wave records: wave number, batch id, chunks dispatched and
        chunks that failed.
    :return: Manifest path when one was written, else None.
    """
    path = failed_chunks_path(ht_path, intervals_hash)
    if not failed:
        if file_exists(path):
            hfs.remove(path)
        return None
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


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def dispatch_in_waves(
    pending: Sequence[int],
    wave_size: int,
    submit_wave: Callable[[List[int], Optional[str]], Optional[int]],
    present_indices: Callable[[], Set[int]],
    dry_run: bool = False,
    on_wave_complete: Optional[
        Callable[[List[int], List[Dict[str, Any]]], None]
    ] = None,
) -> Tuple[List[int], List[Dict[str, Any]]]:
    """
    Submit pending chunks in sequential waves and report the ones that did not land.

    Waves cap how many relays, and so how many nested QoB drivers, run at once. A
    batch run does not raise when a job fails, so the chunk directory is listed
    again after each wave to find the chunks that did not land.

    :param pending: Chunk indices to dispatch.
    :param wave_size: Chunks per wave; ``<= 0`` means one wave.
    :param submit_wave: Submits one wave and returns its batch id.
    :param present_indices: Lists the chunk indices with ``_SUCCESS``.
    :param dry_run: Validate the first wave's batch without running it, then stop.
    :param on_wave_complete: Called after each wave with the failed indices and
        wave records so far, e.g. to rewrite the failed-chunk manifest so it
        survives an orchestrator death. Not called on a dry run.
    :return: ``(failed indices, wave records)``.
    """
    if not pending:
        return [], []
    if wave_size <= 0 or wave_size >= len(pending):
        waves = [list(pending)]
    else:
        waves = [
            list(pending[i : i + wave_size]) for i in range(0, len(pending), wave_size)
        ]
    n_waves = len(waves)
    logger.info(
        "Dispatching %d pending chunks in %d sequential wave(s).", len(pending), n_waves
    )
    failed: List[int] = []
    records: List[Dict[str, Any]] = []
    for wi, wave in enumerate(waves, start=1):
        label = f"w{wi:03d}of{n_waves:03d}" if n_waves > 1 else None
        logger.info(
            "Wave %d/%d: submitting %d chunks (indices %d..%d).",
            wi,
            n_waves,
            len(wave),
            wave[0],
            wave[-1],
        )
        batch_id = submit_wave(wave, label)
        if dry_run:
            logger.info("Dry run: wave DAG validated; stopping.")
            return [], []
        present = present_indices()
        wave_failed = [i for i in wave if i not in present]
        failed.extend(wave_failed)
        records.append(
            {
                "wave": wi,
                "batch_id": batch_id,
                "n_dispatched": len(wave),
                "failed_chunk_indices": wave_failed,
            }
        )
        if wave_failed:
            logger.warning(
                "Wave %d/%d complete but %d/%d chunk(s) MISSING after run: %s%s",
                wi,
                n_waves,
                len(wave_failed),
                len(wave),
                wave_failed[:25],
                " ..." if len(wave_failed) > 25 else "",
            )
        else:
            logger.info(
                "Wave %d/%d complete; all %d chunks present.", wi, n_waves, len(wave)
            )
        if on_wave_complete is not None:
            on_wave_complete(failed, records)
    return failed, records


def orchestrate_tree_merge(
    inputs: Sequence[str],
    final_output_path: str,
    group_path_fn: Callable[[int, int], str],
    merge_group_size: int,
    overwrite: bool,
    final_coalesce_to: Optional[int],
    batch_prefix: str,
    job_prefix: str,
    setup_cmd: str,
    script: str,
    common_flags: str,
    merge_cpu: float,
    merge_memory: str,
    merge_storage: str,
    final_merge_storage: str,
    submit: Callable[[str, List[RelayJobSpec], str], Optional[int]],
    attempts: int = 1,
) -> None:
    """
    Union ``inputs`` into ``final_output_path`` with a tree of relay merge jobs.

    With thousands of chunks, one job that unions everything would need a huge
    input list and too much memory, so the merge is a tree. Each level groups its
    inputs ``merge_group_size`` at a time and submits one ``--run-merge`` job per
    group; a group whose HT already exists is skipped unless ``overwrite`` is set.
    The last level is a single job that writes the final HT, and it always runs and
    overwrites, so a caller that wants to keep an existing final HT checks for it
    before calling.

    :param inputs: Chunk HT paths, in order.
    :param final_output_path: Final HT path.
    :param group_path_fn: ``(level, group_idx) -> intermediate HT path``.
    :param merge_group_size: Inputs per merge job.
    :param overwrite: Recompute intermediate groups that already exist.
    :param final_coalesce_to: ``--merge-coalesce-to`` for the final job, or None.
    :param batch_prefix: Batch name prefix; level and ``final`` are appended.
    :param job_prefix: Job name prefix.
    :param setup_cmd: Relay shell prefix from :func:`build_setup_command`.
    :param script: Script command inside the relay container.
    :param common_flags: Flags appended to every ``--run-merge`` command.
    :param merge_cpu: CPU per merge job.
    :param merge_memory: Memory class per merge job.
    :param merge_storage: Storage per intermediate merge job.
    :param final_merge_storage: Storage for the final merge job.
    :param submit: ``(batch_name, job_specs, log_label) -> batch id``.
    :param attempts: Batch attempts per merge job. Default 1.
    :return: None.
    """
    gs = merge_group_size
    shape = [len(inputs)]
    while shape[-1] > gs:
        shape.append((shape[-1] + gs - 1) // gs)
    logger.info(
        "Merge tree (group_size=%d): %s -> 1 final HT (%d intermediate level(s))",
        gs,
        " -> ".join(str(n) for n in shape),
        len(shape) - 1,
    )

    def merge_job(
        name: str, output: str, group_inputs: Sequence[str], coalesce, storage
    ):
        coalesce_flag = (
            f" --merge-coalesce-to {coalesce}" if coalesce is not None else ""
        )
        return RelayJobSpec(
            name=name,
            cpu=merge_cpu,
            memory=merge_memory,
            storage=storage,
            attempts=attempts,
            command=(
                f"{setup_cmd}{script} --run-merge --merge-output {output}{coalesce_flag}"
                f" --merge-inputs {' '.join(group_inputs)} {common_flags}"
            ),
        )

    inputs = list(inputs)
    level = 1
    while len(inputs) > gs:
        groups = [inputs[i : i + gs] for i in range(0, len(inputs), gs)]
        out_paths = [group_path_fn(level, idx) for idx in range(len(groups))]
        pending = [
            idx
            for idx in range(len(groups))
            if overwrite or not file_exists(out_paths[idx])
        ]
        logger.info(
            "Level %d: %d groups, %d pending, %d already complete (%d -> %d inputs)",
            level,
            len(groups),
            len(pending),
            len(groups) - len(pending),
            len(inputs),
            len(groups),
        )
        specs = [
            merge_job(
                f"{job_prefix}_L{level:02d}_{idx:06d}",
                out_paths[idx],
                groups[idx],
                len(groups[idx]),
                merge_storage,
            )
            for idx in pending
        ]
        submit(f"{batch_prefix}_L{level:02d}", specs, "merge")
        inputs = out_paths
        level += 1

    logger.info("Final merge: %d inputs -> %s", len(inputs), final_output_path)
    final = merge_job(
        f"{job_prefix}_final",
        final_output_path,
        inputs,
        final_coalesce_to,
        final_merge_storage,
    )
    submit(f"{batch_prefix}_final", [final], "final-merge")


def union_and_write_hts(
    input_paths: Sequence[str], output_path: str, coalesce_to: Optional[int] = None
) -> None:
    """
    Union HTs that share a schema and globals, and write the result.

    Needs a running Hail context. ``Table.union`` keeps the globals of the first
    input, which loses nothing because every chunk was built from the same inputs
    and carries the same globals.

    :param input_paths: HT paths to union.
    :param output_path: Destination path, overwritten.
    :param coalesce_to: ``naive_coalesce`` target before writing, or None.
    :return: None.
    """
    logger.info(
        "Merging %d HTs -> %s (coalesce_to=%s)",
        len(input_paths),
        output_path,
        coalesce_to,
    )
    merged = hl.Table.union(*[hl.read_table(p) for p in input_paths])
    if coalesce_to is not None:
        merged = merged.naive_coalesce(coalesce_to)
    merged.write(output_path, overwrite=True)
    logger.info("Wrote merged HT to %s", output_path)
