"""
Compute per-pixel (grid-bin) allele frequencies for gnomAD v4.

Input: a bin assignment TSV with columns idx, bin_id, eval_fold, bin_size.
  - idx maps to gnomAD sample IDs via the gridmax_sample_key resource HT.
  - bin_id is an A_20 lattice cell (21 '_'-separated integers summing to 0).
  - eval_fold: 0 = training fold (used for AF and privacy filtering), 1 = held-out.
  - Samples absent from the TSV have bin_id=missing and are excluded from all counts.

Output: adj call stats per (bin, fold), keyed by bin_id in the `freq` (training)
and `freq_eval` (held-out) dicts. --aggregate-hierarchy adds coarser dyadic levels
(freq_x{factor} / freq_eval_x{factor}); --finalize suppresses groups with fewer
than --min-bin-size samples at every level, filtering each dict by its own fold's
sample count (training-fold count for freq/freq_x{factor}, eval-fold count for
freq_eval/freq_eval_x{factor}).

Pipeline steps (run all by default; use flags to run individual steps):
  --load-bin-assignments  Load TSV, map idx→s, add coarsened bin_id_x{f} columns, write training_samples.ht
  --export-bin-summary    Export raw per-(bin, fold) pop×data_type counts to bin_summary.tsv
  --prep-vds              Filter VDS to exome samples with bin assignments (both folds), split multiallelics
  --compute-freq          Densify VDS and compute per-pixel adj frequencies for both eval folds
  --correct-high-ab-hets  Apply high-AB-het GATK correction to both fold freq dicts
  --aggregate-hierarchy   Sum fine-bin call stats into coarser lattice levels; add sample-count globals
  --finalize              Write privacy-filtered final HT (suppress small groups at all levels)
"""

import argparse
import logging
import re
from typing import Optional

import hail as hl

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.annotations.generate_freq import (
    densify_and_prep_vds_for_freq,
    get_vds_for_freq,
)
from gnomad_qc.v4.resources.annotations import (
    get_gridmax_bin_summary,
    get_gridmax_freq,
    get_gridmax_split_vds,
    get_gridmax_training_samples,
    gridmax_sample_key,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gridmax_freq")
logger.setLevel(logging.INFO)

DEFAULT_MIN_BIN_SIZE = 50
DEFAULT_AB_CUTOFF = 0.9
DEFAULT_AF_THRESHOLD = 0.01


def get_gridmax_pipeline_resources(
    overwrite: bool = False,
    test: bool = False,
    chrom: Optional[str] = None,
    test_gene: bool = False,
) -> PipelineResourceCollection:
    """
    Get the gridmax pipeline resource collection.

    :param overwrite: Whether to overwrite existing outputs.
    :param test: Whether to use tmp paths for testing.
    :param chrom: Optional chromosome suffix for output paths.
    :param test_gene: Whether to use the test-gene (DRD2) suffix in output paths.
    :return: PipelineResourceCollection for the gridmax pipeline.
    """
    gridmax_pipeline = PipelineResourceCollection(
        pipeline_name="gridmax_freq",
        overwrite=overwrite,
    )
    load_bin_assignments_step = PipelineStepResourceCollection(
        "--load-bin-assignments",
        input_resources={
            "compute_gridmax_freq.py (sample key HT)": {
                "gridmax_sample_key": gridmax_sample_key
            }
        },
        output_resources={"ht": get_gridmax_training_samples()},
    )
    export_bin_summary_step = PipelineStepResourceCollection(
        "--export-bin-summary",
        pipeline_input_steps=[load_bin_assignments_step],
        output_resources={"bin_summary": get_gridmax_bin_summary()},
    )
    prep_vds_step = PipelineStepResourceCollection(
        "--prep-vds",
        pipeline_input_steps=[load_bin_assignments_step],
        output_resources={
            "split_vds": get_gridmax_split_vds(
                test=test, chrom=chrom, test_gene=test_gene
            )
        },
    )
    compute_freq_step = PipelineStepResourceCollection(
        "--compute-freq",
        pipeline_input_steps=[prep_vds_step, load_bin_assignments_step],
        output_resources={
            "freq_ht": get_gridmax_freq(
                step="freq", test=test, chrom=chrom, test_gene=test_gene
            )
        },
    )
    correct_high_ab_hets_step = PipelineStepResourceCollection(
        "--correct-high-ab-hets",
        pipeline_input_steps=[compute_freq_step],
        output_resources={
            "corrected_freq_ht": get_gridmax_freq(
                step="corrected", test=test, chrom=chrom, test_gene=test_gene
            )
        },
    )
    aggregate_hierarchy_step = PipelineStepResourceCollection(
        "--aggregate-hierarchy",
        pipeline_input_steps=[correct_high_ab_hets_step, load_bin_assignments_step],
        output_resources={
            "hierarchy_freq_ht": get_gridmax_freq(
                step="hierarchy", test=test, chrom=chrom, test_gene=test_gene
            )
        },
    )
    finalize_step = PipelineStepResourceCollection(
        "--finalize",
        pipeline_input_steps=[aggregate_hierarchy_step],
        output_resources={
            "final_freq_ht": get_gridmax_freq(
                step="final", test=test, chrom=chrom, test_gene=test_gene
            )
        },
    )
    gridmax_pipeline.add_steps(
        {
            "load_bin_assignments": load_bin_assignments_step,
            "export_bin_summary": export_bin_summary_step,
            "prep_vds": prep_vds_step,
            "compute_freq": compute_freq_step,
            "correct_high_ab_hets": correct_high_ab_hets_step,
            "aggregate_hierarchy": aggregate_hierarchy_step,
            "finalize": finalize_step,
        }
    )
    return gridmax_pipeline


def load_bin_assignments(bin_tsv_path: str) -> hl.Table:
    """
    Load bin assignments and return the sample key HT with bin columns added.

    Starts from the full gridmax sample key HT (all samples keyed by idx) and
    annotates each sample with bin_id, bin_size, eval_fold, and one bin_id_x{f}
    column per dyadic coarsening level (parent IDs at each hierarchy level). The
    set of levels is derived from the data: powers of 2 from /2 up to and
    including the factor that collapses every bin to the single root cell.

    :param bin_tsv_path: GCS path to bin assignments TSV.
    :return: Table keyed by 's' with fields bin_id, bin_id_x{f}, bin_size, eval_fold, pop, data_type.
    """
    logger.info("Reading sample key from %s...", gridmax_sample_key.path)
    sample_key_ht = gridmax_sample_key.ht()

    logger.info("Loading bin assignments from %s...", bin_tsv_path)
    bin_ht = hl.import_table(
        bin_tsv_path,
        delimiter="\t",
        types={"idx": hl.tint64, "eval_fold": hl.tint32, "bin_size": hl.tint32},
        key="idx",
    )

    max_abs = bin_ht.aggregate(
        hl.agg.max(hl.max(bin_ht.bin_id.split("_").map(lambda s: hl.abs(hl.int32(s)))))
    )
    factors = _dyadic_factors(max_abs)
    logger.info(
        "Max |coordinate| = %d; hierarchy coarsening factors = %s", max_abs, factors
    )

    bin_row = bin_ht[sample_key_ht.key]
    sample_key_ht = sample_key_ht.annotate(
        bin_id=bin_row.bin_id,
        bin_size=bin_row.bin_size,
        eval_fold=bin_row.eval_fold,
        **{f"bin_id_x{f}": _coarsen_bin_id_hl(bin_row.bin_id, f) for f in factors},
    )
    return sample_key_ht.key_by("s").select(
        "bin_id",
        *[f"bin_id_x{f}" for f in factors],
        "bin_size",
        "eval_fold",
        "pop",
        "data_type",
    )


def export_bin_summary(ht: hl.Table, output_path: str) -> None:
    """
    Export raw per-bin sample counts broken down by pop, data_type, and eval_fold.

    Coarsened bin_id_x{factor} columns present in ht are included in the output,
    enabling the plotter to aggregate at any resolution. No privacy filter is
    applied here; the plotter (and --finalize) suppress small groups using
    training-fold counts.

    :param ht: Training HT keyed by 's' with bin_id, bin_id_x{factor}, pop, data_type, eval_fold.
    :param output_path: Full GCS path for the output TSV.
    """
    coarse_cols = sorted(
        [c for c in ht.row if re.match(r"bin_id_x\d+$", c)],
        key=lambda c: int(c[len("bin_id_x") :]),
    )
    ht = ht.filter(hl.is_defined(ht.bin_id))
    ht = ht.group_by(
        "bin_id",
        *coarse_cols,
        "pop",
        "data_type",
        "eval_fold",
    ).aggregate(n=hl.agg.count())
    logger.info("Exporting bin summary to %s...", output_path)
    ht.export(output_path)


def prep_vds_for_gridmax(
    ht: hl.Table,
    chrom: Optional[str] = None,
    test_gene: bool = False,
) -> hl.vds.VariantDataset:
    """
    Prep gnomAD v4 VDS for gridmax frequency computation.

    Calls get_vds_for_freq for standard VDS prep (multi-allelic split, entry
    annotations), then filters cols to exome samples with bin assignments (both
    eval folds) and annotates bin_id and eval_fold.

    :param ht: Table keyed by 's' with fields bin_id and eval_fold.
    :param chrom: Optional chromosome to restrict to (e.g. 'chr20').
    :param test_gene: If True, restrict to DRD2 for quick testing.
    :return: Prepped VDS with bin_id and eval_fold col annotations, both folds included.
    """
    vds = get_vds_for_freq(test_gene=test_gene, chrom=chrom)

    counts = ht.aggregate(
        hl.struct(
            n_exome=hl.agg.filter(ht.data_type == "exome", hl.agg.count()),
            total=hl.agg.count(),
        )
    )
    logger.info(
        "VDS is exomes-only; restricting to %d exome samples "
        "(%d genome samples in joint key will not contribute to AC/AN).",
        counts.n_exome,
        counts.total - counts.n_exome,
    )
    ht = ht.filter(hl.is_defined(ht.bin_id) & (ht.data_type == "exome"))

    logger.info(
        "Filtering VDS to exome samples with bin assignments (both eval folds)..."
    )
    vmt = vds.variant_data
    rmt = vds.reference_data
    vmt = vmt.annotate_cols(bin_id=ht[vmt.s].bin_id, eval_fold=ht[vmt.s].eval_fold)
    vmt = vmt.filter_cols(hl.is_defined(ht[vmt.s]))
    vmt = vmt.select_cols("sex_karyotype", "fixed_homalt_model", "bin_id", "eval_fold")
    rmt = rmt.filter_cols(hl.is_defined(ht[rmt.s]))
    rmt = rmt.filter_rows(hl.agg.count() > 0)

    return hl.vds.VariantDataset(rmt, vmt)


def compute_pixel_freq(
    vds: hl.vds.VariantDataset,
    ab_cutoff: float = DEFAULT_AB_CUTOFF,
) -> hl.Table:
    """
    Densify VDS and compute per-bin adj frequencies using hl.agg.group_by.

    Uses a single group_by aggregation to avoid JVM ClassTooLargeException. All
    bins (including small ones) are included so that they can contribute to coarser
    hierarchy levels in the --aggregate-hierarchy step.

    Call stats are stored as alt-only scalars (AC, AF, homozygote_count for the
    alt allele, plus AN), matching the release freq schema. The reference count is
    redundant given AN (AC_ref = AN - AC, AF_ref = 1 - AF) and is not stored.

    :param vds: VariantDataset with bin_id and eval_fold col annotations (both folds).
    :param ab_cutoff: Allele balance cutoff for high-AB-het annotation.
    :return: Table with freq, freq_eval, high_ab_hets, high_ab_hets_eval dicts keyed by bin_id.
    """
    logger.info("Densifying VDS for frequency computation...")
    mt = densify_and_prep_vds_for_freq(vds, ab_cutoff=ab_cutoff)
    logger.info("Aggregating per-bin frequencies with group_by (both eval folds)...")
    high_ab_het_flag = mt.GT.is_het_ref() & ~mt.fixed_homalt_model & mt._high_ab_het_ref
    ht = mt.select_rows(
        freq=hl.agg.filter(
            mt.adj & (mt.eval_fold == 0),
            hl.agg.group_by(mt.bin_id, hl.agg.call_stats(mt.GT, mt.alleles)),
        ),
        freq_eval=hl.agg.filter(
            mt.adj & (mt.eval_fold == 1),
            hl.agg.group_by(mt.bin_id, hl.agg.call_stats(mt.GT, mt.alleles)),
        ),
        high_ab_hets=hl.agg.filter(
            mt.adj & (mt.eval_fold == 0),
            hl.agg.group_by(mt.bin_id, hl.agg.sum(hl.int(high_ab_het_flag))),
        ),
        high_ab_hets_eval=hl.agg.filter(
            mt.adj & (mt.eval_fold == 1),
            hl.agg.group_by(mt.bin_id, hl.agg.sum(hl.int(high_ab_het_flag))),
        ),
    ).rows()
    # hl.agg.call_stats returns [ref, alt] arrays; keep alt-only scalars.
    return ht.annotate(
        freq=_flatten_call_stats(ht.freq),
        freq_eval=_flatten_call_stats(ht.freq_eval),
    )


def _flatten_call_stats(freq_expr: hl.DictExpression) -> hl.DictExpression:
    """Reduce a bin_id -> array-form call_stats dict to alt-only scalar structs."""
    return hl.dict(
        freq_expr.items().map(
            lambda kv: hl.tuple(
                [
                    kv[0],
                    hl.struct(
                        AC=kv[1].AC[1],
                        AF=kv[1].AF[1],
                        AN=kv[1].AN,
                        homozygote_count=kv[1].homozygote_count[1],
                    ),
                ]
            )
        )
    )


def _apply_high_ab_het_correction(stats, n, af_threshold):
    # Reclassifying n high-AB het calls as hom-alt adds one alt allele per sample
    # (AC += n) with AN unchanged; the alt hom count rises by n.
    new_ac = stats.AC + hl.int32(n)
    new_af = hl.if_else(
        stats.AN > 0, hl.float64(new_ac) / stats.AN, hl.missing(hl.tfloat64)
    )
    corrected = hl.struct(
        AC=new_ac,
        AF=new_af,
        AN=stats.AN,
        homozygote_count=stats.homozygote_count + hl.int32(n),
    )
    return hl.if_else(stats.AF > af_threshold, corrected, stats)


def correct_high_ab_hets(
    ht: hl.Table, af_threshold: float = DEFAULT_AF_THRESHOLD
) -> hl.Table:
    """
    Correct call statistics for high-AB-het GATK artifact.

    For each bin, if adj AF exceeds af_threshold, adjusts AC and homozygote_count
    upward by the per-bin high-AB-het count and recomputes AF.

    :param ht: Frequency Table with freq and high_ab_hets dicts keyed by bin_id.
    :param af_threshold: Per-bin AF above which to apply correction.
    :return: Table with corrected freq dict; high_ab_hets dropped.
    """

    def _correct_dict(freq_expr, hets_expr):
        return hl.dict(
            freq_expr.items().map(
                lambda kv: hl.tuple(
                    [
                        kv[0],
                        hl.bind(
                            lambda stats, n: _apply_high_ab_het_correction(
                                stats, n, af_threshold
                            ),
                            kv[1],
                            hets_expr[kv[0]],
                        ),
                    ]
                )
            )
        )

    ht = ht.annotate(
        freq=_correct_dict(ht.freq, ht.high_ab_hets),
        freq_eval=_correct_dict(ht.freq_eval, ht.high_ab_hets_eval),
    )
    return ht.drop("high_ab_hets", "high_ab_hets_eval").annotate_globals(
        af_threshold_for_correction=af_threshold
    )


def _dyadic_factors(max_abs: int) -> list:
    """
    Return dyadic coarsening factors [2, 4, 8, ...] up to the root.

    Includes the first power of 2 that exceeds max_abs, which collapses every
    coordinate to 0 (the single root cell). Coarsening uses truncation toward
    zero, so this is the point at which all bins merge.

    :param max_abs: Largest absolute coordinate value across all bin IDs.
    :return: List of powers of 2 defining the hierarchy levels.
    """
    factors = []
    f = 2
    while f // 2 <= max_abs:
        factors.append(f)
        f *= 2
    return factors


def _coarsen_bin_id_hl(
    bin_id_expr: hl.StringExpression, factor: int
) -> hl.StringExpression:
    """
    Coarsen a bin ID by dividing each coordinate by factor, truncating toward zero.

    Truncation toward zero (not floor division) ensures negative coordinates
    collapse to 0 at large factors, so the hierarchy reaches a single root cell.
    """
    return hl.if_else(
        hl.is_defined(bin_id_expr),
        hl.delimit(
            bin_id_expr.split("_").map(
                lambda s: hl.str(hl.int32(hl.int32(s) / factor))
            ),
            "_",
        ),
        hl.missing(hl.tstr),
    )


def _sum_call_stats(a, b):
    return hl.struct(
        AC=a.AC + b.AC,
        AF=hl.missing(hl.tfloat64),
        AN=a.AN + b.AN,
        homozygote_count=a.homozygote_count + b.homozygote_count,
    )


def _recompute_af(stats):
    return stats.annotate(
        AF=hl.if_else(
            stats.AN > 0, hl.float64(stats.AC) / stats.AN, hl.missing(hl.tfloat64)
        )
    )


def _aggregate_to_coarser_level(
    freq_expr: hl.DictExpression,
    children_map: hl.DictExpression,
) -> hl.DictExpression:
    zero = hl.struct(
        AC=hl.int32(0),
        AF=hl.missing(hl.tfloat64),
        AN=hl.int32(0),
        homozygote_count=hl.int32(0),
    )
    return hl.dict(
        children_map.items().map(
            lambda kv: hl.tuple(
                [
                    kv[0],
                    _recompute_af(
                        kv[1].fold(
                            lambda acc, child_id: _sum_call_stats(
                                acc, hl.coalesce(freq_expr.get(child_id), zero)
                            ),
                            zero,
                        )
                    ),
                ]
            )
        )
    )


def aggregate_hierarchy(
    freq_ht: hl.Table,
    training_ht: hl.Table,
) -> hl.Table:
    """
    Aggregate per-bin call stats to coarser hierarchy levels.

    Each A_20 bin ID is a 21-dimensional integer coordinate vector (coordinates
    sum to 0). Coarser bins are the bin_id_x{factor} columns precomputed by
    --load-bin-assignments; children at the finer level are summed into the
    parent named by those columns. The parent assignments are read directly from
    the training HT (never recomputed here), so this always matches the bin
    summary and the training HT exactly.

    Levels are discovered from the bin_id_x{factor} columns present in training_ht.

    Adds freq_x{factor} / freq_eval_x{factor} dict fields and per-fold, exome-only
    sample-count globals (sample_counts / sample_counts_eval and their _x{factor}
    variants) for use by --finalize.

    :param freq_ht: Corrected frequency Table with freq dict keyed by bin_id.
    :param training_ht: Training samples Table keyed by 's' with bin_id and bin_id_x{factor} fields.
    :return: freq_ht with additional coarser-level freq dicts and sample count globals.
    """
    factors = sorted(
        int(m.group(1))
        for c in training_ht.row
        for m in [re.match(r"bin_id_x(\d+)$", c)]
        if m
    )
    logger.info("Hierarchy factors discovered from training HT: %s", factors)
    coarse_cols = [f"bin_id_x{f}" for f in factors]

    # Group by bin_id and its precomputed parent IDs (constant within a bin_id).
    # The bin_id_x{f} column values ARE the parent assignments, so the hierarchy
    # here matches --load-bin-assignments and the bin summary exactly.
    #
    # Count exome samples only: the freq/freq_eval dicts are computed from exome
    # samples (prep_vds_for_gridmax filters to data_type == "exome"), so the
    # privacy sample counts must reflect exactly the samples in those arrays.
    filt_ht = training_ht.filter(
        hl.is_defined(training_ht.bin_id) & (training_ht.data_type == "exome")
    )
    counts_rows = (
        filt_ht.group_by("bin_id", *coarse_cols)
        .aggregate(
            n_train=hl.agg.filter(filt_ht.eval_fold == 0, hl.agg.count()),
            n_eval=hl.agg.filter(filt_ht.eval_fold == 1, hl.agg.count()),
        )
        .collect()
    )
    l0_train = {row.bin_id: row.n_train for row in counts_rows}
    l0_eval = {row.bin_id: row.n_eval for row in counts_rows}

    level_fields = {}
    global_fields = {"sample_counts": l0_train, "sample_counts_eval": l0_eval}

    for factor in factors:
        col = f"bin_id_x{factor}"
        children_map = {}
        coarse_train = {}
        coarse_eval = {}
        for row in counts_rows:
            parent = row[col]
            children_map.setdefault(parent, []).append(row.bin_id)
            coarse_train[parent] = coarse_train.get(parent, 0) + row.n_train
            coarse_eval[parent] = coarse_eval.get(parent, 0) + row.n_eval

        logger.info("Factor /%d: %d coarser groups", factor, len(children_map))
        children_map_lit = hl.literal(
            children_map, hl.tdict(hl.tstr, hl.tarray(hl.tstr))
        )
        level_fields[f"freq_x{factor}"] = _aggregate_to_coarser_level(
            freq_ht.freq, children_map_lit
        )
        level_fields[f"freq_eval_x{factor}"] = _aggregate_to_coarser_level(
            freq_ht.freq_eval, children_map_lit
        )
        global_fields[f"sample_counts_x{factor}"] = coarse_train
        global_fields[f"sample_counts_eval_x{factor}"] = coarse_eval

    freq_ht = freq_ht.annotate(**level_fields)
    return freq_ht.annotate_globals(
        **{k: hl.literal(v) for k, v in global_fields.items()}
    )


def finalize_freq_ht(freq_ht: hl.Table, min_bin_size: int) -> hl.Table:
    """
    Write a privacy-filtered freq HT, suppressing groups below min_bin_size at all levels.

    Sample counts at each level are taken from globals added by --aggregate-hierarchy.

    :param freq_ht: Hierarchy freq HT with freq, freq_x2, ... and sample_counts* globals.
    :param min_bin_size: Groups with fewer than this many samples are suppressed.
    :return: freq_ht with small groups filtered from all freq dicts and from the
        sample_counts* globals, so suppressed group sizes are not exposed.
    """
    globals_eval = hl.eval(freq_ht.index_globals())
    factors = sorted(
        int(m.group(1))
        for k in globals_eval
        for m in [re.match(r"sample_counts_x(\d+)$", k)]
        if m
    )

    def large_set(counts_dict):
        return hl.literal(
            {k for k, v in counts_dict.items() if v >= min_bin_size},
            hl.tset(hl.tstr),
        )

    def scrub(counts_dict):
        return hl.literal(
            {k: v for k, v in counts_dict.items() if v >= min_bin_size},
            hl.tdict(hl.tstr, hl.tint64),
        )

    large_train = large_set(globals_eval.sample_counts)
    large_eval = large_set(globals_eval.sample_counts_eval)
    result = freq_ht.annotate(
        freq=hl.dict(
            freq_ht.freq.items().filter(lambda kv: large_train.contains(kv[0]))
        ),
        freq_eval=hl.dict(
            freq_ht.freq_eval.items().filter(lambda kv: large_eval.contains(kv[0]))
        ),
    )

    # Scrub the sample-count globals to match the row-level suppression: keep counts
    # only for released (large) groups, so the exact sizes of suppressed small groups
    # do not leak via the "privacy-filtered" HT's globals.
    scrubbed_globals = {
        "sample_counts": scrub(globals_eval.sample_counts),
        "sample_counts_eval": scrub(globals_eval.sample_counts_eval),
    }
    for factor in factors:
        large = large_set(getattr(globals_eval, f"sample_counts_x{factor}"))
        large_ev = large_set(getattr(globals_eval, f"sample_counts_eval_x{factor}"))
        result = result.annotate(
            **{
                f"freq_x{factor}": hl.dict(
                    getattr(result, f"freq_x{factor}")
                    .items()
                    .filter(lambda kv: large.contains(kv[0]))
                ),
                f"freq_eval_x{factor}": hl.dict(
                    getattr(result, f"freq_eval_x{factor}")
                    .items()
                    .filter(lambda kv: large_ev.contains(kv[0]))
                ),
            }
        )
        scrubbed_globals[f"sample_counts_x{factor}"] = scrub(
            getattr(globals_eval, f"sample_counts_x{factor}")
        )
        scrubbed_globals[f"sample_counts_eval_x{factor}"] = scrub(
            getattr(globals_eval, f"sample_counts_eval_x{factor}")
        )
        logger.info(
            "factor /%d: %d train / %d eval groups pass min_bin_size=%d",
            factor,
            hl.eval(hl.len(large)),
            hl.eval(hl.len(large_ev)),
            min_bin_size,
        )
    return result.annotate_globals(min_bin_size=min_bin_size, **scrubbed_globals)


def main(args):
    """Run the gridmax frequency pipeline."""
    hl.init(log="/tmp/gridmax.log", tmp_dir=args.tmp_dir)

    overwrite = args.overwrite
    chrom = args.chrom
    test_gene = args.test_gene
    test = test_gene

    resources = get_gridmax_pipeline_resources(
        overwrite=overwrite, test=test, chrom=chrom, test_gene=test_gene
    )

    run_all = not any(
        [
            args.load_bin_assignments,
            args.export_bin_summary,
            args.prep_vds,
            args.compute_freq,
            args.correct_high_ab_hets,
            args.aggregate_hierarchy,
            args.finalize,
        ]
    )

    if run_all or args.load_bin_assignments:
        if not args.bin_assignments:
            raise ValueError("--bin-assignments is required for --load-bin-assignments")
        res = resources.load_bin_assignments
        res.check_resource_existence()
        logger.info("Step 1: Load bin assignments and map to gnomAD sample IDs...")
        ht = load_bin_assignments(args.bin_assignments)
        ht.checkpoint(res.ht.path, overwrite=overwrite)
        logger.info("Samples in key HT: %d", res.ht.ht().count())

    if run_all or args.export_bin_summary:
        res = resources.export_bin_summary
        res.check_resource_existence()
        export_bin_summary(
            res.ht.ht(),
            res.bin_summary,
        )

    if run_all or args.prep_vds:
        res = resources.prep_vds
        res.check_resource_existence()
        logger.info("Step 2: Prep and split VDS...")
        vds = prep_vds_for_gridmax(
            res.ht.ht(),
            chrom=chrom,
            test_gene=test_gene,
        )
        vds.checkpoint(res.split_vds.path, overwrite=overwrite)

    if run_all or args.compute_freq:
        res = resources.compute_freq
        res.check_resource_existence()
        logger.info("Step 3: Densify and compute per-pixel frequencies...")
        freq_ht = compute_pixel_freq(
            res.split_vds.vds(),
            ab_cutoff=args.ab_cutoff,
        )
        freq_ht.checkpoint(res.freq_ht.path, overwrite=overwrite)

    if run_all or args.correct_high_ab_hets:
        res = resources.correct_high_ab_hets
        res.check_resource_existence()
        logger.info("Step 4: Apply high-AB-het correction...")
        freq_ht = correct_high_ab_hets(res.freq_ht.ht(), af_threshold=args.af_threshold)
        freq_ht.checkpoint(res.corrected_freq_ht.path, overwrite=overwrite)

    if run_all or args.aggregate_hierarchy:
        res = resources.aggregate_hierarchy
        res.check_resource_existence()
        logger.info("Step 5: Aggregate call stats across hierarchy levels...")
        freq_ht = aggregate_hierarchy(res.corrected_freq_ht.ht(), res.ht.ht())
        freq_ht.checkpoint(res.hierarchy_freq_ht.path, overwrite=overwrite)

    if run_all or args.finalize:
        res = resources.finalize
        res.check_resource_existence()
        logger.info("Step 6: Write privacy-filtered final HT...")
        freq_ht = finalize_freq_ht(res.hierarchy_freq_ht.ht(), args.min_bin_size)
        freq_ht.checkpoint(res.final_freq_ht.path, overwrite=overwrite)
        freq_ht.describe()
        logger.info("Done. Final output: %s", res.final_freq_ht.path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    step_group = parser.add_argument_group(
        "pipeline steps",
        "Select one or more steps to run. If none are specified, all steps run in order.",
    )
    step_group.add_argument(
        "--load-bin-assignments",
        action="store_true",
        help="Step 1: Load bin TSV, map idx→s, checkpoint training_samples.ht.",
    )
    step_group.add_argument(
        "--export-bin-summary",
        action="store_true",
        help=(
            "Export raw per-bin counts by pop, data_type, and eval_fold (with "
            "coarsened bin_id_x{factor} columns) to bin_summary.tsv. No privacy "
            "filter applied; downstream code suppresses small groups."
        ),
    )
    step_group.add_argument(
        "--prep-vds",
        action="store_true",
        help=(
            "Step 2: Filter VDS to exome samples with bin assignments (both eval "
            "folds), split multiallelics, checkpoint split_vds.vds."
        ),
    )
    step_group.add_argument(
        "--compute-freq",
        action="store_true",
        help=(
            "Step 3: Densify VDS and compute per-pixel adj frequencies for both "
            "eval folds (freq, freq_eval), checkpoint freq.ht."
        ),
    )
    step_group.add_argument(
        "--correct-high-ab-hets",
        action="store_true",
        help="Step 4: Apply high-AB-het correction, write freq.corrected.ht.",
    )
    step_group.add_argument(
        "--aggregate-hierarchy",
        action="store_true",
        help=(
            "Step 5: Aggregate per-bin call stats to coarser A_20 lattice levels "
            "(dyadic factors up to the root, set at --load-bin-assignments). Adds "
            "freq_x{factor} fields and sample_counts globals. Writes freq.hierarchy.ht."
        ),
    )
    step_group.add_argument(
        "--finalize",
        action="store_true",
        help=(
            "Step 6: Suppress groups below --min-bin-size at all hierarchy levels "
            "and write freq.final.ht."
        ),
    )
    load_group = parser.add_argument_group("--load-bin-assignments (step 1) parameters")
    load_group.add_argument(
        "--bin-assignments",
        help=(
            "GCS path to bin assignments TSV "
            "(columns: idx, bin_id, eval_fold, bin_size). "
            "Required for --load-bin-assignments."
        ),
    )

    compute_group = parser.add_argument_group("--compute-freq (step 3) parameters")
    compute_group.add_argument(
        "--ab-cutoff",
        type=float,
        default=DEFAULT_AB_CUTOFF,
        help=f"Allele balance cutoff for high-AB-het annotation (default: {DEFAULT_AB_CUTOFF}).",
    )

    correct_group = parser.add_argument_group(
        "--correct-high-ab-hets (step 4) parameters"
    )
    correct_group.add_argument(
        "--af-threshold",
        type=float,
        default=DEFAULT_AF_THRESHOLD,
        help=(
            f"AF threshold above which to apply high-AB-het correction "
            f"(default: {DEFAULT_AF_THRESHOLD})."
        ),
    )

    finalize_group = parser.add_argument_group("--finalize (step 6) parameters")
    finalize_group.add_argument(
        "--min-bin-size",
        type=int,
        default=DEFAULT_MIN_BIN_SIZE,
        help=(
            f"Privacy threshold for group size (default: {DEFAULT_MIN_BIN_SIZE}). "
            "In --finalize, groups with fewer than this many samples are removed at "
            "every hierarchy level, filtering each dict by its own fold's sample "
            "count (training-fold count for freq/freq_x{factor}, eval-fold count for "
            "freq_eval/freq_eval_x{factor}). freq.hierarchy.ht is unaffected."
        ),
    )

    scope_group = parser.add_argument_group(
        "scope",
        "Restrict the run and set the output-path suffix; apply to --prep-vds onward.",
    )
    scope_group.add_argument(
        "--chrom",
        help="Restrict computation to this chromosome (e.g. chr20).",
    )
    scope_group.add_argument(
        "--test-gene",
        action="store_true",
        help=(
            "Restrict --prep-vds onward to DRD2 (chr11:113409605-113475691) for quick "
            "testing, writing those outputs to gnomad-tmp. --load-bin-assignments and "
            "--export-bin-summary are per-sample and gene-independent, so they always "
            "write to the production paths and are not isolated by this flag."
        ),
    )

    general_group = parser.add_argument_group("general")
    general_group.add_argument(
        "--tmp-dir",
        default="gs://gnomad-tmp-4day",
        help="Hail temporary directory (default: gs://gnomad-tmp-4day).",
    )
    general_group.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs.",
    )

    args = parser.parse_args()
    main(args)
