"""
Compute per-pixel (grid-bin) allele frequencies for gnomAD v4.

Input: a bin assignment TSV with columns idx, bin_id, eval_fold, bin_size.
  - idx maps to gnomAD sample IDs via the gridmax_sample_key resource HT.
  - All samples in the key HT are used; bin columns are null for samples absent
    from the TSV, which causes them to contribute to adj-overall only.

Pipeline steps (run all by default; use flags to run individual steps):
  --load-bin-assignments  Load TSV, map idx→s, write training_samples.ht
  --export-bin-summary    Export pop×data_type counts per bin to bin_summary.tsv
  --prep-vds              Filter VDS to exome samples with bin assignments, split multiallelics
  --compute-freq          Densify VDS and compute per-pixel frequencies
  --correct-high-ab-hets  Apply high-AB-het GATK correction
  --finalize              Write privacy-filtered final HT (suppress small bins)

Small bins (bin_size < --min-bin-size) are excluded from per-bin strata during
--compute-freq; those samples contribute to adj-overall only. --finalize writes
freq.final.ht with freq_index_dict recomputed from the filtered freq_meta.
"""

import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.utils.annotations import (
    build_freq_stratification_list,
    compute_freq_by_strata,
    generate_freq_group_membership_array,
)
from gnomad.utils.release import make_freq_index_dict_from_meta

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


def high_ab_het(mt: hl.MatrixTable) -> hl.Int32Expression:
    """
    Return 1 if a call is a high allele balance heterozygous call, else 0.

    compute_freq_by_strata calls entry_agg_funcs as f(mt), so this takes the
    full MT rather than separate (entry, col) structs as in generate_freq.py.

    :param mt: Input MatrixTable with GT, adj, fixed_homalt_model, _high_ab_het_ref.
    :return: 1 if high-AB het call, else 0.
    """
    return hl.int(
        mt.GT.is_het_ref() & mt.adj & ~mt.fixed_homalt_model & mt._high_ab_het_ref
    )


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
        output_resources={"ht": get_gridmax_training_samples(test=test)},
    )
    export_bin_summary_step = PipelineStepResourceCollection(
        "--export-bin-summary",
        pipeline_input_steps=[load_bin_assignments_step],
        output_resources={"bin_summary": get_gridmax_bin_summary(test=test)},
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
    finalize_step = PipelineStepResourceCollection(
        "--finalize",
        pipeline_input_steps=[correct_high_ab_hets_step, load_bin_assignments_step],
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
            "finalize": finalize_step,
        }
    )
    return gridmax_pipeline


def load_bin_assignments(bin_tsv_path: str) -> hl.Table:
    """
    Load bin assignments and return the sample key HT with bin columns added.

    Starts from the full gridmax sample key HT (all samples keyed by idx) and
    annotates each sample with bin_id, bin_size, and eval_fold from the TSV.
    Samples not present in the TSV receive missing values for those columns and
    will contribute to the adj-overall group only in frequency computation.

    :param bin_tsv_path: GCS path to bin assignments TSV.
    :return: Table keyed by 's' with fields bin_id, bin_size, eval_fold, pop, data_type.
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

    bin_row = bin_ht[sample_key_ht.key]
    sample_key_ht = sample_key_ht.annotate(
        bin_id=bin_row.bin_id,
        bin_size=bin_row.bin_size,
        eval_fold=bin_row.eval_fold,
    )
    return sample_key_ht.key_by("s").select(
        "bin_id", "bin_size", "eval_fold", "pop", "data_type"
    )


def export_bin_summary(ht: hl.Table, output_path: str, min_bin_size: int) -> None:
    """
    Export per-bin sample counts broken down by pop and data_type.

    Produces a TSV with columns: bin_id, pop, data_type, eval_fold, n.
    Bins with bin_size < min_bin_size are labeled 'suppressed' for display;
    the internal training HT retains their true bin_id.
    Intended for local plotting via plot_bin_summary.py.

    :param ht: Table keyed by 's' with fields bin_id, bin_size, pop, data_type.
    :param output_path: Full GCS path for the output TSV.
    :param min_bin_size: Bins below this size are labeled 'suppressed' in the TSV.
    """
    ht = ht.filter(hl.is_defined(ht.bin_id))
    ht = ht.annotate(
        bin_id=hl.if_else(ht.bin_size < min_bin_size, "suppressed", ht.bin_id)
    )
    ht = ht.group_by("bin_id", "pop", "data_type", "eval_fold").aggregate(
        n=hl.agg.count()
    )
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
    annotations), then filters cols to training-fold samples and annotates bin_id.

    :param ht: Table keyed by 's' with field bin_id (possibly missing).
    :param chrom: Optional chromosome to restrict to (e.g. 'chr20').
    :param test_gene: If True, restrict to DRD2 for quick testing.
    :return: Prepped VDS filtered to training samples with bin_id col annotation.
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
    ht = ht.filter(ht.data_type == "exome")

    logger.info(
        "Filtering VDS to samples with bin assignments and annotating bin_id..."
    )
    vmt = vds.variant_data
    rmt = vds.reference_data
    vmt = vmt.annotate_cols(bin_id=ht[vmt.s].bin_id)
    vmt = vmt.filter_cols(hl.is_defined(ht[vmt.s]))
    vmt = vmt.select_cols("sex_karyotype", "fixed_homalt_model", "bin_id")
    rmt = rmt.filter_cols(hl.is_defined(ht[rmt.s]))
    rmt = rmt.filter_rows(hl.agg.count() > 0)

    return hl.vds.VariantDataset(rmt, vmt)


def compute_pixel_freq(
    vds: hl.vds.VariantDataset,
    ht: hl.Table,
    min_bin_size: int,
    ab_cutoff: float = DEFAULT_AB_CUTOFF,
) -> hl.Table:
    """
    Densify VDS and compute per-pixel frequencies using bin_id as the pop stratification.

    Samples in bins with fewer than min_bin_size training samples have their
    bin_id set to missing and contribute to the adj-overall count only.
    This keeps the number of strata within Hail's JVM class-size limits.

    :param vds: VariantDataset filtered to training samples with bin_id col annotation.
    :param ht: Table keyed by 's' with field bin_size.
    :param min_bin_size: Bins below this size are excluded from per-bin freq computation.
    :param ab_cutoff: Allele balance cutoff for high-AB-het annotation.
    :return: Frequency Table with adj + per-bin freq array and freq_meta global.
    """
    logger.info("Densifying VDS for frequency computation...")
    mt = densify_and_prep_vds_for_freq(vds, ab_cutoff=ab_cutoff)
    # TODO: bin_size is from the joint (exome+genome) TSV but AC/AN is exome-only;
    #   suppression thresholds should use exome-only counts. Revisit after chr20 test.
    mt = mt.annotate_cols(
        bin_id=hl.or_missing(ht[mt.s].bin_size >= min_bin_size, mt.bin_id)
    )
    meta_ht = mt.cols()
    logger.info("Building frequency stratification by bin_id...")
    strata_expr = build_freq_stratification_list(gen_anc_expr=meta_ht.bin_id)
    group_membership_ht = generate_freq_group_membership_array(meta_ht, strata_expr)

    group_membership_globals = group_membership_ht.index_globals().annotate(
        non_ukb_downsamplings=hl.missing(hl.tarray(hl.tint)),
        non_ukb_ds_pop_counts=hl.missing(hl.tdict(hl.tstr, hl.tint)),
    )

    logger.info("Aggregating per-pixel frequencies...")
    freq_ht = compute_freq_by_strata(
        mt.annotate_cols(
            group_membership=group_membership_ht[mt.col_key].group_membership
        ),
        entry_agg_funcs={"high_ab_hets_by_group": (high_ab_het, hl.agg.sum)},
    )
    freq_ht = freq_ht.annotate_globals(**group_membership_globals)
    return freq_ht


def correct_high_ab_hets(
    ht: hl.Table, af_threshold: float = DEFAULT_AF_THRESHOLD
) -> hl.Table:
    """
    Correct call statistics for high-AB-het GATK artifact.

    At sites where the raw adj AF exceeds af_threshold, adjusts AC and
    homozygote_count upward by the per-group high-AB-het count.

    Note: the existing correct_for_high_ab_hets in generate_freq.py also corrects
    histograms, which we don't compute here, so this is a freq-only variant.

    :param ht: Frequency Table with freq and high_ab_hets_by_group arrays.
    :param af_threshold: AF above which to apply correction.
    :return: Table with corrected freq array and freq_index_dict global.
    """
    corrected = hl.map(
        lambda f, g: hl.struct(
            AC=hl.int32(f.AC + g),
            AF=hl.if_else(f.AN > 0, (f.AC + g) / f.AN, hl.missing(hl.tfloat64)),
            AN=f.AN,
            homozygote_count=f.homozygote_count + g,
        ),
        ht.freq,
        ht.high_ab_hets_by_group,
    )
    # Raw group (index 1) must remain unadjusted; restore it from the original freq.
    corrected = hl.zip(hl.range(hl.len(corrected)), corrected).map(
        lambda pair: hl.if_else(pair[0] == 1, ht.freq[pair[0]], pair[1])
    )
    ht = ht.annotate(
        freq=hl.if_else(
            ht.freq[0].AF > af_threshold, corrected, ht.freq, missing_false=True
        )
    )
    ht = ht.drop("high_ab_hets_by_group")
    ht = ht.annotate_globals(
        freq_index_dict=make_freq_index_dict_from_meta(ht.freq_meta),
        af_threshold_for_correction=af_threshold,
    )
    return ht


def finalize_freq_ht(freq_ht: hl.Table, ht: hl.Table, min_bin_size: int) -> hl.Table:
    """
    Write a privacy-filtered freq HT by removing small-bin entries.

    Collects bin_ids where bin_size < min_bin_size from the training HT, then
    removes the corresponding entries from freq, freq_meta, and
    freq_meta_sample_count, keeping all three arrays index-aligned. The overall
    adj entry (no gen_anc) and all large-bin entries are unchanged.
    freq_index_dict is recomputed from the filtered freq_meta.

    :param freq_ht: Corrected frequency Table with freq_meta global.
    :param ht: Table keyed by 's' with fields bin_id, bin_size.
    :param min_bin_size: Bins below this size are removed from the output.
    :return: freq_ht with small-bin entries dropped from freq, freq_meta, and
        freq_meta_sample_count, and freq_index_dict recomputed.
    """
    small_bins = ht.aggregate(
        hl.agg.filter(
            ht.bin_size < min_bin_size,
            hl.agg.collect_as_set(ht.bin_id),
        )
    )
    logger.info(
        "Suppressing %d bins with bin_size < %d in final output...",
        len(small_bins),
        min_bin_size,
    )
    small_bins_lit = hl.literal(small_bins, hl.tset(hl.tstr))
    ht_globals = freq_ht.index_globals()
    freq_meta = ht_globals.freq_meta
    suppress_mask = freq_meta.map(
        lambda m: hl.coalesce(small_bins_lit.contains(m.get("gen_anc")), False)
    )
    freq_filtered = (
        hl.zip(freq_ht.freq, suppress_mask)
        .filter(lambda pair: ~pair[1])
        .map(lambda pair: pair[0])
    )
    freq_meta_filtered = (
        hl.zip(freq_meta, suppress_mask)
        .filter(lambda pair: ~pair[1])
        .map(lambda pair: pair[0])
    )
    freq_meta_sample_count_filtered = (
        hl.zip(ht_globals.freq_meta_sample_count, suppress_mask)
        .filter(lambda pair: ~pair[1])
        .map(lambda pair: pair[0])
    )
    return freq_ht.annotate(freq=freq_filtered).annotate_globals(
        freq_meta=freq_meta_filtered,
        freq_meta_sample_count=freq_meta_sample_count_filtered,
        freq_index_dict=make_freq_index_dict_from_meta(freq_meta_filtered),
    )


def main(args):
    """Run the gridmax frequency pipeline."""
    hl.init(tmp_dir=args.tmp_dir)

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
            args.finalize,
        ]
    )

    if run_all or args.load_bin_assignments:
        res = resources.load_bin_assignments
        res.check_resource_existence()
        logger.info("Step 1: Load bin assignments and map to gnomAD sample IDs...")
        ht = load_bin_assignments(args.bin_assignments)
        ht.checkpoint(res.ht.path, overwrite=overwrite)
        logger.info("Samples in key HT: %d", res.ht.ht().count())

    if run_all or args.export_bin_summary:
        res = resources.export_bin_summary
        res.check_resource_existence()
        export_bin_summary(res.ht.ht(), res.bin_summary, args.min_bin_size)

    if run_all or args.prep_vds:
        res = resources.prep_vds
        res.check_resource_existence()
        logger.info("Step 2: Prep and split VDS...")
        vds = prep_vds_for_gridmax(res.ht.ht(), chrom=chrom, test_gene=test_gene)
        vds.checkpoint(res.split_vds.path, overwrite=overwrite)

    if run_all or args.compute_freq:
        res = resources.compute_freq
        res.check_resource_existence()
        logger.info("Step 3: Densify and compute per-pixel frequencies...")
        freq_ht = compute_pixel_freq(
            res.split_vds.vds(),
            res.ht.ht(),
            min_bin_size=args.min_bin_size,
            ab_cutoff=args.ab_cutoff,
        )
        freq_ht.checkpoint(res.freq_ht.path, overwrite=overwrite)

    if run_all or args.correct_high_ab_hets:
        res = resources.correct_high_ab_hets
        res.check_resource_existence()
        logger.info("Step 4: Apply high-AB-het correction...")
        freq_ht = correct_high_ab_hets(res.freq_ht.ht(), af_threshold=args.af_threshold)
        freq_ht.checkpoint(res.corrected_freq_ht.path, overwrite=overwrite)

    if run_all or args.finalize:
        res = resources.finalize
        res.check_resource_existence()
        logger.info("Step 5: Write privacy-filtered final HT...")
        freq_ht = finalize_freq_ht(
            res.corrected_freq_ht.ht(), res.ht.ht(), args.min_bin_size
        )
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
        help="Export per-bin pop×data_type counts to bin_summary.tsv alongside training_samples.ht.",
    )
    step_group.add_argument(
        "--prep-vds",
        action="store_true",
        help="Step 2: Filter VDS to training samples, split multiallelics, checkpoint split_vds.vds.",
    )
    step_group.add_argument(
        "--compute-freq",
        action="store_true",
        help="Step 3: Densify VDS and compute per-pixel frequencies, checkpoint freq.ht.",
    )
    step_group.add_argument(
        "--correct-high-ab-hets",
        action="store_true",
        help="Step 4: Apply high-AB-het correction, write freq.corrected.ht.",
    )
    step_group.add_argument(
        "--finalize",
        action="store_true",
        help=(
            "Step 5: Suppress small-bin freq entries (bin_size < --min-bin-size) "
            "and write freq.final.ht. The internal corrected HT is unchanged."
        ),
    )
    parser.add_argument(
        "--bin-assignments",
        required=True,
        help=(
            "GCS path to bin assignments TSV "
            "(columns: idx, bin_id, eval_fold, bin_size)."
        ),
    )
    parser.add_argument(
        "--chrom",
        help="Restrict computation to this chromosome (e.g. chr20).",
    )
    parser.add_argument(
        "--test-gene",
        action="store_true",
        help="Restrict to DRD2 (chr11:113409605-113475691) for quick testing. Outputs go to gnomad-tmp.",
    )
    parser.add_argument(
        "--min-bin-size",
        type=int,
        default=DEFAULT_MIN_BIN_SIZE,
        help=(
            f"Privacy threshold for bin size (default: {DEFAULT_MIN_BIN_SIZE}). "
            "Bins below this size are labeled 'suppressed' in bin_summary.tsv and "
            "are removed from freq, freq_meta, and freq_meta_sample_count in "
            "freq.final.ht (arrays shortened, freq_index_dict recomputed). "
            "The internal training HT and freq.corrected.ht are unaffected."
        ),
    )
    parser.add_argument(
        "--ab-cutoff",
        type=float,
        default=DEFAULT_AB_CUTOFF,
        help=f"Allele balance cutoff for high-AB-het annotation (default: {DEFAULT_AB_CUTOFF}).",
    )
    parser.add_argument(
        "--af-threshold",
        type=float,
        default=DEFAULT_AF_THRESHOLD,
        help=(
            f"AF threshold above which to apply high-AB-het correction "
            f"(default: {DEFAULT_AF_THRESHOLD})."
        ),
    )
    parser.add_argument(
        "--tmp-dir",
        default="gs://gnomad-tmp-4day",
        help="Hail temporary directory (default: gs://gnomad-tmp-4day).",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing outputs.",
    )
    args = parser.parse_args()
    main(args)
