"""
Compute genome-wide coverage from a stratified subset of AoU CRAMs.

The AoU v8 VDS carries no DP field (ref blocks have only banded GQ; variant
rows carry LAD, whose sum is informative-read depth only), so release coverage
is measured directly from CRAMs on a 2.5-5k stratified sample subset. Depth
matches the historical definition (samtools depth: -a, base quality >= 10,
MAPQ >= 20, duplicates excluded) over the primary contigs minus
telomeres/centromeres. Results are RAW (uncapped): `mean` is exact raw;
`median` is raw, right-censored at --hist-ceiling (default 1000; censored
sites are visible via the `n_above_ceiling` column and reported by
--validate). `mean_capped` additionally applies the historical per-sample
100x cap, computed exactly from the in-job histogram — a per-sample cap
cannot be applied after per-site averaging, so it cannot be deferred to the
HT stage without shipping full histograms.

All CRAM- and person_id-touching steps run as Hail Batch jobs (workstations
cannot reach the AoU bucket); only per-site aggregates leave the jobs. Regions
are processed independently (one job per ~10 Mb region reads that region from
every CRAM via the .crai), so no per-sample intermediates are ever stored and
a preemption only re-reads one region.

Steps:
  --run-unit-tests      Local, no cloud access: checks the histogram summary
                        math against brute force, stratum allocation, and
                        region construction. Run before review / any submission.
  --select-samples      Stratified draw (release wave x sequencing center,
                        proportional with a per-stratum floor) from the CRAM
                        manifest + genomic_metrics; writes the person_id-keyed
                        selection to the batch bucket and prints an ID-free
                        stratum/coverage-distribution summary to the job log.
  --compute-coverage    One batch job per region: per-site histogram over all
                        selected samples -> mean, median, over_X fractions;
                        writes region TSVs (aggregates only).
  --aggregate-regions   Import region TSVs into a locus-keyed Hail Table.
  --validate            Sanity-check the coverage HT (row count vs footprint,
                        constant n, over_X monotonicity/ranges) and print
                        summary quantiles for review.

Test run (cheap end-to-end, ~50 samples x 2 chr22 regions, tmp paths):
  python cram_coverage.py --test --select-samples
  python cram_coverage.py --test --compute-coverage
  python cram_coverage.py --test --aggregate-regions
  python cram_coverage.py --test --validate
"""

import argparse
import inspect
import logging
import urllib.request
from typing import Dict, List, Tuple

import hailtop.batch as hb

from gnomad_qc.v5.resources.annotations import (
    cram_coverage_path,
    cram_coverage_region_root,
    cram_coverage_selection_path,
)
from gnomad_qc.v5.resources.constants import AOU_GENOMIC_METRICS_PATH

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("cram_coverage")
logger.setLevel(logging.INFO)

BATCH_IMAGE = "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:0.2.137"
RP_PROJECT = "broad-mpg-gnomad"
AOU_CRAM_MANIFEST = "gs://fc-aou-datasets-controlled/v8/wgs/cram/manifest.csv"
REF_FASTA = (
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta"
)
TELOMERES_CENTROMERES_BED = (
    "https://storage.googleapis.com/gcp-public-data--gnomad/resources/grch38/"
    "telomeres_and_centromeres/hg38.telomeresAndMergedCentromeres.bed"
)
OVER_THRESHOLDS = [1, 5, 10, 15, 20, 25, 30, 50, 100]
# Historical per-sample cap, applied only to the `mean_capped` column.
DP_CAP = 100
# Default histogram bin ceiling (--hist-ceiling): raw depths above it share the
# top bin, right-censoring `median` there. Memory is (ceiling+1) x region_len
# x 2 bytes. 1000 is ~25x the genome-wide median depth; only pathological
# collapsed-repeat sites can censor, and they are flagged via n_above_ceiling.
HIST_CEILING = 1000

GRCH38_CONTIG_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}


def allocate_strata(sizes, n_samples: int, floor_n: int):
    """
    Allocate samples across strata proportionally with a per-stratum floor.

    Floors are capped at the stratum size; the total is adjusted to exactly
    `n_samples` by adding to the strata with the most headroom / removing from
    the largest above-floor allocations.

    :param sizes: pandas Series of stratum sizes.
    :param n_samples: Total samples to allocate.
    :param floor_n: Minimum per stratum (capped at stratum size).
    :return: pandas Series of per-stratum allocations.
    """
    import numpy as np

    assert n_samples <= sizes.sum(), "n_samples exceeds cohort size"
    floors = np.minimum(floor_n, sizes)
    alloc = (sizes / sizes.sum() * n_samples).round().astype(int)
    alloc = np.maximum(alloc, floors)
    alloc = np.minimum(alloc, sizes)
    while alloc.sum() != n_samples:
        if alloc.sum() < n_samples:
            headroom = sizes - alloc
            target = headroom.idxmax()
            assert headroom[target] > 0, "no headroom left"
            alloc[target] += 1
        else:
            above = alloc[alloc > floors]
            target = above.idxmax() if len(above) else alloc.idxmax()
            alloc[target] -= 1
    return alloc


def hist_to_summary(hist, n_samples: int, thresholds, cap: int):
    """
    Convert a per-site raw-depth histogram to coverage summary arrays.

    :param hist: (ceiling+1, n_sites) array; hist[d, i] = samples with raw
        depth d (clipped to the ceiling) at site i. The ceiling must be >= cap,
        which makes `mean_capped` exact.
    :param n_samples: Expected sample count at every site.
    :param thresholds: over_X thresholds (each <= cap, so cap-invariant).
    :param cap: Per-sample depth cap for `mean_capped`.
    :return: (mean_capped, median, {threshold: fraction >= threshold},
        n_above_ceiling) arrays; `median` is raw, right-censored at the ceiling.
    """
    import numpy as np

    assert hist.shape[0] - 1 >= cap, "histogram ceiling must be >= cap"
    n = hist.sum(axis=0)
    assert (n == n_samples).all(), "histogram counts != sample count"
    capped_depths = np.minimum(np.arange(hist.shape[0], dtype=np.float64), cap)
    mean_capped = (capped_depths[:, None] * hist).sum(axis=0) / n_samples
    cum = hist.cumsum(axis=0)
    # Lower median: smallest depth d with cumulative count >= n/2.
    median = (cum < (n_samples / 2)).sum(axis=0)
    over = {t: (n_samples - cum[t - 1]) / n_samples for t in thresholds}
    return mean_capped, median, over, hist[-1]


SELECT_PY = (
    r"""
import sys

import numpy as np
import pandas as pd
import hailtop.fs as hfs

"""
    + inspect.getsource(allocate_strata)
    + r"""

out_path, n_samples, floor_n, seed = (
    sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
)

with hfs.open(sys.argv[5], "r") as f:
    m = pd.read_csv(f, dtype=str)
with hfs.open(sys.argv[6], "r") as f:
    met = pd.read_csv(f, sep="\t", dtype=str)

m["wave"] = m.cram_uri.str.extract(r"/cram/([^/]+)/")
met = met[["research_id", "site_id", "mean_coverage"]]
df = m.merge(met, left_on="person_id", right_on="research_id", how="left")
df["batch"] = df.wave.fillna("unknown") + "|" + df.site_id.fillna("unknown")
df["mean_coverage"] = pd.to_numeric(df.mean_coverage, errors="coerce")

sizes = df.batch.value_counts()
alloc = allocate_strata(sizes, n_samples, floor_n)
sel = pd.concat(
    [
        df[df.batch == b].sample(n, random_state=seed)
        for b, n in alloc.items()
        if n > 0
    ]
)
assert len(sel) == n_samples, len(sel)
assert sel.cram_uri.notna().all() and sel.cram_index_uri.notna().all()
sel = sel.sort_values(["batch", "person_id"])

print("=== per-stratum allocation (cohort n, selected n)")
for b in sizes.index:
    print(f"{b}\t{sizes[b]}\t{alloc.get(b, 0)}")
print("=== mean_coverage deciles: cohort vs selected")
q = [i / 10 for i in range(1, 10)]
print(df.mean_coverage.quantile(q).round(2).to_string())
print(sel.mean_coverage.quantile(q).round(2).to_string())

sel[["person_id", "cram_uri", "cram_index_uri", "batch"]].to_csv(
    out_path, sep="\t", header=False, index=False
)
"""
)

REF_PY = r"""
import os

import hailtop.fs as hfs

ref = os.environ["REF_GS"]
for src, dst in [(ref, "ref.fa"), (ref + ".fai", "ref.fa.fai")]:
    hfs.copy(src, "file://" + os.path.abspath(dst))
print("ref.fa bytes:", os.path.getsize("ref.fa"), flush=True)
"""

# One region, all samples: workers stream capped per-site depth vectors; the
# main process accumulates a single (DP_CAP+1 x region_len) histogram, then
# writes per-site aggregates. (d, arange) index pairs are unique per sample,
# so plain fancy-indexed += is a correct single-pass increment.
DEPTH_PY = (
    r"""
import multiprocessing as mp
import os
import sys

import numpy as np
import pandas as pd

"""
    + inspect.getsource(hist_to_summary)
    + r"""

sel_path, contig, start, end, cap, out_tsv = (
    sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]),
    int(sys.argv[5]), sys.argv[6],
)
thresholds = [int(t) for t in sys.argv[7].split(",")]
ceiling = int(sys.argv[8])
region = f"{contig}:{start}-{end}"
length = end - start + 1

rows = [l.split("\t") for l in open(sel_path).read().splitlines() if l.strip()]
print(f"{region}: {len(rows)} samples", flush=True)


def sample_depth(row):
    import google.auth.transport.requests
    import pysam
    from google.oauth2 import service_account

    person, cram, crai, batch = row
    os.environ["GCS_REQUESTER_PAYS_PROJECT"] = os.environ["RP_PROJECT"]
    creds = service_account.Credentials.from_service_account_file(
        "/gsa-key/key.json",
        scopes=["https://www.googleapis.com/auth/devstorage.read_only"],
    )
    creds.refresh(google.auth.transport.requests.Request())
    os.environ["GCS_OAUTH_TOKEN"] = creds.token

    bam = f"slice_{person}.bam"
    pysam.view("-b", "-T", "ref.fa", "-X", cram, crai, region, "-o", bam,
               catch_stdout=False)
    pysam.index(bam)
    dtsv = f"depth_{person}.tsv"
    pysam.depth("-a", "-r", region, "-q", "10", "-Q", "20", "-o", dtsv, bam,
                catch_stdout=False)
    d = pd.read_csv(dtsv, sep="\t", header=None, usecols=[1, 2],
                    names=["pos", "DP"], dtype=np.int32)
    os.remove(bam); os.remove(bam + ".bai"); os.remove(dtsv)
    vec = np.zeros(length, dtype=np.int32)
    vec[d.pos.to_numpy() - start] = d.DP.to_numpy()
    return vec


hist = np.zeros((ceiling + 1, length), dtype=np.uint16)
raw_sum = np.zeros(length, dtype=np.int64)
idx = np.arange(length)
with mp.Pool(processes=int(os.environ.get("N_PROCS", "8"))) as pool:
    for i, vec in enumerate(pool.imap_unordered(sample_depth, rows, chunksize=1)):
        raw_sum += vec
        hist[np.minimum(vec, ceiling), idx] += 1
        if (i + 1) % 250 == 0:
            print(f"{region}: {i + 1}/{len(rows)} samples done", flush=True)

mean_capped, median, over, n_above = hist_to_summary(hist, len(rows), thresholds, cap)
out = pd.DataFrame({"contig": contig, "pos": idx + start,
                    "mean": (raw_sum / len(rows)).round(4),
                    "mean_capped": mean_capped.round(4), "median": median})
for t in thresholds:
    out[f"over_{t}"] = over[t].round(6)
out["n"] = len(rows)
out["n_above_ceiling"] = n_above
out.to_csv("region.tsv", sep="\t", index=False)

import pysam
pysam.tabix_compress("region.tsv", "region.tsv.bgz", force=True)
import hailtop.fs as hfs
hfs.copy("file://" + os.path.abspath("region.tsv.bgz"), out_tsv)
print(f"{region}: wrote {out_tsv}", flush=True)
"""
)


def build_regions(region_mb: int, contigs: List[str]) -> List[Tuple[str, int, int]]:
    """
    Build mask-complement regions chopped to `region_mb`.

    :param region_mb: Maximum region size in Mb.
    :param contigs: Contigs to include.
    :return: List of (contig, start, end) 1-based inclusive regions.
    """
    mask = load_mask()
    step = region_mb * 1_000_000
    regions = []
    for contig in contigs:
        pos = 1
        breaks = sorted(mask.get(contig, [])) + [
            (GRCH38_CONTIG_LENGTHS[contig] + 1, GRCH38_CONTIG_LENGTHS[contig] + 1)
        ]
        for m_start, m_end in breaks:
            keep_start, keep_end = pos, min(m_start - 1, GRCH38_CONTIG_LENGTHS[contig])
            while keep_start <= keep_end:
                regions.append(
                    (contig, keep_start, min(keep_start + step - 1, keep_end))
                )
                keep_start += step
            pos = max(pos, m_end + 1)
    return regions


def load_mask() -> Dict[str, List[Tuple[int, int]]]:
    """
    Load telomere/centromere intervals (1-based inclusive) from the public BED.

    :return: Dict of contig -> sorted (start, end) intervals.
    """
    mask: Dict[str, List[Tuple[int, int]]] = {}
    with urllib.request.urlopen(TELOMERES_CENTROMERES_BED) as f:
        for line in f.read().decode().splitlines():
            if not line.strip():
                continue
            c, s, e = line.split()[:3]
            mask.setdefault(c, []).append((int(s) + 1, int(e)))
    return mask


def run_unit_tests() -> None:
    """Run local unit tests (no cloud access); raises on failure."""
    import numpy as np
    import pandas as pd

    # hist_to_summary vs brute force on random raw depth data (odd n for an
    # exact lower median == np.median). Depths exceed both cap and ceiling to
    # exercise capping and censoring.
    rng = np.random.default_rng(0)
    n_samples, n_sites, cap, ceiling = 37, 500, 20, 60
    depth = rng.integers(0, ceiling + 30, size=(n_samples, n_sites))
    hist = np.zeros((ceiling + 1, n_sites), dtype=np.uint16)
    raw_sum = np.zeros(n_sites, dtype=np.int64)
    idx = np.arange(n_sites)
    for vec in depth:  # same accumulation as DEPTH_PY
        raw_sum += vec
        hist[np.minimum(vec, ceiling), idx] += 1
    thresholds = [1, 5, 10, cap]
    mean_capped, median, over, n_above = hist_to_summary(
        hist, n_samples, thresholds, cap
    )
    assert np.allclose(raw_sum / n_samples, depth.mean(axis=0)), "raw mean mismatch"
    assert np.allclose(
        mean_capped, np.minimum(depth, cap).mean(axis=0)
    ), "capped mean mismatch"
    assert np.array_equal(
        median, np.minimum(np.median(depth, axis=0), ceiling)
    ), "median mismatch"
    for t in thresholds:
        assert np.allclose(over[t], (depth >= t).mean(axis=0)), f"over_{t} mismatch"
    assert np.array_equal(n_above, (depth >= ceiling).sum(axis=0)), "n_above mismatch"

    # allocate_strata: exact total, floors respected, sizes never exceeded.
    for n, floor_n in [(5000, 250), (2500, 250), (100, 5), (7, 2)]:
        sizes = pd.Series({"a": 200_000, "b": 50_000, "c": 900, "d": 12})
        if n > sizes.sum():
            continue
        alloc = allocate_strata(sizes, n, floor_n)
        assert alloc.sum() == n, (n, alloc.sum())
        assert (alloc <= sizes).all(), "allocation exceeds stratum size"
        feasible = np.minimum(floor_n, sizes)
        if feasible.sum() <= n:
            assert (alloc >= feasible).all(), "floor violated"

    # build_regions: within bounds, <= step, disjoint from the mask, and the
    # footprint equals contig length minus in-bounds mask bases.
    contig = "chr22"
    mask = load_mask()
    regions = [r for r in build_regions(10, [contig])]
    length = GRCH38_CONTIG_LENGTHS[contig]
    masked = sum(
        min(e, length) - max(s, 1) + 1
        for s, e in mask[contig]
        if s <= length and e >= 1
    )
    total = sum(e - s + 1 for _, s, e in regions)
    assert total == length - masked, (total, length - masked)
    for _, s, e in regions:
        assert 1 <= s <= e <= length
        assert e - s + 1 <= 10_000_000
        for m_s, m_e in mask[contig]:
            assert e < m_s or s > m_e, "region overlaps mask"

    logger.info("All unit tests passed.")


def get_batch(args) -> hb.Batch:
    """
    Create a Hail Batch with the standard v5 backend.

    :param args: Argparse arguments (billing project, remote tmpdir).
    :return: Batch object.
    """
    backend = hb.ServiceBackend(
        billing_project=args.billing_project, remote_tmpdir=args.remote_tmpdir
    )
    return hb.Batch(name="cram-coverage", backend=backend)


def main(args):
    """Compute CRAM-based coverage on a stratified AoU sample subset."""
    if args.run_unit_tests:
        run_unit_tests()
        return

    test = args.test
    n_samples = args.n_samples or (50 if test else 5000)
    stratum_floor = args.stratum_floor if not test else min(args.stratum_floor, 5)
    selection_path = cram_coverage_selection_path(test=test)
    region_root = cram_coverage_region_root(test=test)
    contigs = ["chr22"] if test else list(GRCH38_CONTIG_LENGTHS)

    if args.select_samples:
        b = get_batch(args)
        j = b.new_bash_job(name="select-samples")
        j.image(BATCH_IMAGE)
        j.cpu(1)
        j.memory("standard")
        j.regions(["us-central1"])
        j.command("set -euo pipefail")
        j.command(f"hailctl config set gcs_requester_pays/project {RP_PROJECT}")
        j.command(f"cat > select.py <<'PYEOF'{SELECT_PY}PYEOF")
        j.command(
            f'python3 select.py "{j.selected}" {n_samples} '
            f'{stratum_floor} {args.seed} "{AOU_CRAM_MANIFEST}" '
            f'"{AOU_GENOMIC_METRICS_PATH}"'
        )
        b.write_output(j.selected, selection_path)
        b.run(wait=False)
        logger.info("Submitted sample selection; selection -> %s", selection_path)

    if args.compute_coverage:
        assert args.hist_ceiling >= DP_CAP, "--hist-ceiling must be >= 100"
        regions = build_regions(args.region_mb, contigs)
        if test:
            regions = regions[: args.test_n_regions]
        logger.info(
            "%d regions, %.0f Mb footprint",
            len(regions),
            sum(e - s + 1 for _, s, e in regions) / 1e6,
        )
        b = get_batch(args)
        selected = b.read_input(selection_path)
        for contig, start, end in regions:
            j = b.new_bash_job(name=f"depth-{contig}-{start}")
            j.image(BATCH_IMAGE)
            j.cpu(args.cpu_per_job)
            j.memory("standard")
            j.storage("30Gi")
            j.regions(["us-central1"])
            if args.non_preemptible and hasattr(j, "spot"):
                j.spot(False)
            j.command("set -euo pipefail")
            j.command("pip install --quiet pysam 2>&1 | tail -1")
            j.command(
                f'export REF_GS="{REF_FASTA}" RP_PROJECT="{RP_PROJECT}" '
                f'N_PROCS="{args.cpu_per_job}"'
            )
            j.command(f"cat > fetch_ref.py <<'PYEOF'{REF_PY}PYEOF")
            j.command("python3 fetch_ref.py")
            j.command(f"cat > depth.py <<'PYEOF'{DEPTH_PY}PYEOF")
            out_tsv = f"{region_root}/region_{contig}_{start}_{end}.tsv.bgz"
            j.command(
                f'python3 depth.py "{selected}" {contig} {start} {end} '
                f'{DP_CAP} "{out_tsv}" {",".join(map(str, OVER_THRESHOLDS))} '
                f"{args.hist_ceiling}"
            )
        b.run(wait=False)
        logger.info(
            "Submitted %d region jobs; outputs -> %s", len(regions), region_root
        )

    if args.aggregate_regions or args.validate:
        import hail as hl

        coverage_ht_path = cram_coverage_path(
            test=test, environment=args.environment
        ).path

    if args.aggregate_regions:
        from gnomad_qc.resource_utils import check_resource_existence

        check_resource_existence(
            output_step_resources={"--aggregate-regions": [coverage_ht_path]},
            overwrite=args.overwrite,
        )
        ht = hl.import_table(
            f"{region_root}/region_*.tsv.bgz",
            force_bgz=True,
            types={
                "pos": hl.tint32,
                "mean": hl.tfloat64,
                "mean_capped": hl.tfloat64,
                "median": hl.tint32,
                "n": hl.tint32,
                "n_above_ceiling": hl.tint32,
                **{f"over_{t}": hl.tfloat64 for t in OVER_THRESHOLDS},
            },
        )
        ht = ht.transmute(locus=hl.locus(ht.contig, ht.pos, reference_genome="GRCh38"))
        ht = ht.key_by("locus")
        ht.write(coverage_ht_path, overwrite=args.overwrite)

    if args.validate:
        ht = hl.read_table(coverage_ht_path)
        regions = build_regions(args.region_mb, contigs)
        if test:
            regions = regions[: args.test_n_regions]
        expected = sum(e - s + 1 for _, s, e in regions)
        over_cols = [f"over_{t}" for t in OVER_THRESHOLDS]
        monotone_expr = hl.literal(True)
        for hi, lo in zip(over_cols, over_cols[1:]):
            monotone_expr = monotone_expr & (ht[hi] >= ht[lo])
        in_range_expr = (
            (ht.mean >= 0)
            & (ht.mean_capped >= 0)
            & (ht.mean_capped <= DP_CAP)
            & (ht.mean_capped <= ht.mean)
            & (ht.median >= 0)
            & (ht.median <= args.hist_ceiling)
            & (ht.n_above_ceiling >= 0)
            & (ht.n_above_ceiling <= ht.n)
        )
        for c in over_cols:
            in_range_expr = in_range_expr & (ht[c] >= 0) & (ht[c] <= 1)
        stats = ht.aggregate(
            hl.struct(
                n_rows=hl.agg.count(),
                n_values=hl.agg.collect_as_set(ht.n),
                monotone=hl.agg.all(monotone_expr),
                in_range=hl.agg.all(in_range_expr),
                mean_q=hl.agg.approx_quantiles(ht.mean, [0.1, 0.5, 0.9]),
                median_q=hl.agg.approx_quantiles(ht.median, [0.1, 0.5, 0.9]),
                over_20_mean=hl.agg.mean(ht.over_20),
                n_censored_median=hl.agg.count_where(ht.median == args.hist_ceiling),
                n_sites_above_ceiling=hl.agg.count_where(ht.n_above_ceiling > 0),
            )
        )
        assert stats.n_rows == expected, f"rows {stats.n_rows} != footprint {expected}"
        assert len(stats.n_values) == 1, f"non-constant n: {stats.n_values}"
        assert stats.monotone, "over_X not monotone decreasing"
        assert stats.in_range, "column values out of range"
        logger.info(
            "Validation passed: %d sites, n=%s, mean deciles %s, median deciles %s, "
            "mean over_20 %.4f, censored medians %d, sites with any sample above "
            "ceiling %d",
            stats.n_rows,
            stats.n_values,
            stats.mean_q,
            stats.median_q,
            stats.over_20_mean,
            stats.n_censored_median,
            stats.n_sites_above_ceiling,
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run-unit-tests",
        help="Run local unit tests (no cloud access) and exit.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Use tmp paths, chr22 regions only, and 50 samples by default.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-regions", type=int, default=2, help="Regions to run with --test."
    )
    parser.add_argument("--overwrite", help="Overwrite data.", action="store_true")
    parser.add_argument("--billing-project", default="gnomad-production")
    parser.add_argument(
        "--remote-tmpdir",
        default="gs://fc-11093c2b-590e-424a-91ac-0cc040d562fc/batch-tmp",
    )
    parser.add_argument(
        "--environment",
        default="batch",
        choices=["rwb", "batch", "dataproc"],
        help="Environment for the aggregated coverage HT path.",
    )

    select_args = parser.add_argument_group("Sample selection")
    select_args.add_argument(
        "--select-samples",
        help="Submit the stratified sample-selection batch job.",
        action="store_true",
    )
    select_args.add_argument(
        "--n-samples",
        type=int,
        help="Number of samples to select. Default 5000 (50 with --test).",
    )
    select_args.add_argument(
        "--stratum-floor",
        type=int,
        default=250,
        help="Minimum samples per wave x center stratum (capped at stratum size).",
    )
    select_args.add_argument("--seed", type=int, default=42)

    depth_args = parser.add_argument_group("Coverage computation")
    depth_args.add_argument(
        "--compute-coverage",
        help="Submit one batch job per region computing per-site aggregates.",
        action="store_true",
    )
    depth_args.add_argument(
        "--region-mb", type=int, default=5, help="Region size in Mb."
    )
    depth_args.add_argument(
        "--hist-ceiling",
        type=int,
        default=HIST_CEILING,
        help="Raw-depth histogram bin ceiling; must be >= 100. Median is "
        "right-censored here (see n_above_ceiling). Job histogram memory is "
        "(ceiling+1) x region length x 2 bytes.",
    )
    depth_args.add_argument(
        "--cpu-per-job", type=int, default=16, help="Cores per region job."
    )
    depth_args.add_argument(
        "--non-preemptible",
        help="Run region jobs non-preemptible (preemption re-reads one region).",
        action="store_true",
    )

    agg_args = parser.add_argument_group("Aggregation and validation")
    agg_args.add_argument(
        "--aggregate-regions",
        help="Import region TSVs into the locus-keyed coverage HT.",
        action="store_true",
    )
    agg_args.add_argument(
        "--validate",
        help="Sanity-check the coverage HT and print summary quantiles.",
        action="store_true",
    )
    return parser


if __name__ == "__main__":
    main(get_script_argument_parser().parse_args())
