"""
Compare VEP 105 vs VEP 115 runs on the gnomAD context HT.

Main Checks:
- schema: Schema comparison between VEP versions
- counts: Variant count comparison per chromosome
- centromere: Chr18 centromere patch region check
- missingness: VEP field missingness rates
- consequences: Most severe consequence distribution
- consistency: Annotation consistency spot check
- loftee: LOFTEE annotation comparison
- sift_polyphen: SIFT/PolyPhen prediction comparison
- mane: MANE transcript coverage (VEP 115)
- biotype: Transcript biotype distribution
- colocated: Colocated variants (dbSNP, ClinVar)
- regulatory: Regulatory/motif consequences

Debug Checks:
- debug_tx_counts: Transcript count distribution analysis
- debug_lof_filter: LoF filter string difference analysis
- debug_mane: MANE field value inspection
- debug_regulatory_motif: Regulatory/motif feature inspection
- debug_sift_polyphen: SIFT/PolyPhen category transitions

Usage:
    # Run all checks
    python temp_vep_check.py

    # Run all checks with debug
    python temp_vep_check.py --include-debug

    # Run only debug checks (faster, for investigating specific issues)
    python temp_vep_check.py --debug-only

    # Run specific check groups
    python temp_vep_check.py --checks schema counts loftee

    # List available checks
    python temp_vep_check.py --list-checks
"""

import argparse
import logging
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.assessment.validity_checks import count_vep_annotated_variants_per_interval
from gnomad.resources.grch38.reference_data import ensembl_interval, vep_context

from gnomad_qc.v5.resources.basics import qc_temp_prefix

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def create_global_sample(
    ht105: hl.Table, ht115: hl.Table, sample_frac: float = 0.01, seed: int = 42
) -> hl.Table:
    """
    Create a globally sampled and joined table for reuse across multiple checks.

    Samples from ht105, joins with ht115, and checkpoints the result for efficient
    reuse. This avoids repeated expensive sampling operations.

    :param ht105: VEP 105 table.
    :param ht115: VEP 115 table.
    :param sample_frac: Fraction of variants to sample (default 1%).
    :param seed: Random seed for reproducible sampling.
    :return: Checkpointed joined table with vep105 and vep115 annotations.
    """
    logger.info(f"Creating global sample ({sample_frac:.2%}) and checkpointing...")

    # Sample from one table first (cheap operation).
    sampled_ht = ht105.sample(sample_frac, seed=seed)

    # Use semi-join pattern: filter ht105 to sampled keys, then join ht115.
    ht = sampled_ht.select(vep105=sampled_ht.vep)
    ht = ht.annotate(vep115=ht115[ht.key].vep)
    # ht = ht.filter(hl.is_defined(ht.vep115))

    # Checkpoint for efficient reuse across multiple checks.
    checkpoint_path = qc_temp_prefix() + "vep_comparison_sample.ht"
    ht = ht.checkpoint(checkpoint_path, overwrite=True)

    logger.info(f"Global sample created and checkpointed: {checkpoint_path}")
    return ht


def compare_schemas(ht105: hl.Table, ht115: hl.Table) -> Dict[str, List[str]]:
    """
    Compare VEP schemas between two tables.

    :param ht105: VEP 105 annotated table.
    :param ht115: VEP 115 annotated table.
    :return: Dict with 'only_in_105', 'only_in_115', 'in_both' field lists.
    """
    logger.info("Comparing VEP schemas...")

    def get_vep_fields(ht: hl.Table, path: str = "vep") -> set:
        """Recursively get all field paths in the VEP struct, handling nested arrays."""
        fields = set()
        vep_type = ht[path].dtype

        def traverse(dtype, prefix):
            if isinstance(dtype, hl.tstruct):
                # Recurse into struct fields.
                for name, subtype in dtype.items():
                    traverse(subtype, f"{prefix}.{name}")
            elif isinstance(dtype, hl.tarray):
                # For arrays, recurse into element type with [] notation.
                element_type = dtype.element_type
                if isinstance(element_type, hl.tstruct):
                    # Array of structs - recurse into struct fields.
                    for name, subtype in element_type.items():
                        traverse(subtype, f"{prefix}[].{name}")
                elif isinstance(element_type, hl.tarray):
                    # Nested array (e.g., Array[Array[...]]) - recurse with nested [].
                    traverse(element_type, f"{prefix}[]")
                else:
                    # Array of primitives - just add the array itself.
                    fields.add(f"{prefix}[]")
            else:
                # Primitive type - add the field.
                fields.add(prefix)

        traverse(vep_type, path)
        return fields

    fields_105 = get_vep_fields(ht105)
    fields_115 = get_vep_fields(ht115)

    result = {
        "only_in_105": sorted(fields_105 - fields_115),
        "only_in_115": sorted(fields_115 - fields_105),
        "in_both": sorted(fields_105 & fields_115),
    }

    logger.info(f"Fields only in VEP 105: {result['only_in_105']}")
    logger.info(f"Fields only in VEP 115: {result['only_in_115']}")
    logger.info(f"Shared fields: {len(result['in_both'])}")

    return result


def count_variants_per_contig(ht: hl.Table) -> Dict[str, int]:
    """
    Count variants per chromosome.

    :param ht: Input table with locus key.
    :return: Dict of contig -> count.
    """
    logger.info("Counting variants per contig...")
    counts = ht.aggregate(hl.agg.counter(ht.locus.contig))
    return dict(sorted(counts.items(), key=lambda x: x[0]))


def compare_variant_counts(
    ht105: hl.Table,
    ht115: hl.Table,
) -> Tuple[Dict[str, int], Dict[str, int], Dict[str, int]]:
    """
    Compare variant counts per chromosome between two tables.

    :param ht105: VEP 105 table.
    :param ht115: VEP 115 table.
    :return: Tuple of (counts_105, counts_115, differences).
    """
    logger.info("Comparing variant counts per contig...")

    counts_105 = count_variants_per_contig(ht105)
    counts_115 = count_variants_per_contig(ht115)

    all_contigs = set(counts_105.keys()) | set(counts_115.keys())
    differences = {}

    for contig in sorted(all_contigs):
        c105 = counts_105.get(contig, 0)
        c115 = counts_115.get(contig, 0)
        if c105 != c115:
            differences[contig] = c115 - c105
            logger.warning(
                f"Contig {contig}: VEP105={c105:,}, VEP115={c115:,}, diff={c115-c105:,}"
            )

    logger.info("Variant counts per contig:")
    logger.info("Contig        VEP105         VEP115         Difference")
    logger.info("-" * 55)
    for contig in sorted(all_contigs):
        c105 = counts_105.get(contig, 0)
        c115 = counts_115.get(contig, 0)
        diff = c115 - c105
        logger.info(f"{contig:<12} {c105:>12,} {c115:>12,} {diff:>+12,}")

    if not differences:
        logger.info("All contigs have matching variant counts!")

    return counts_105, counts_115, differences


def check_vep_missingness_from_sample(
    sampled_ht: hl.Table,
) -> Tuple[Dict[str, float], Dict[str, float]]:
    """
    Check missingness rates for key VEP fields using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Tuple of (missingness_105, missingness_115) dicts.
    """
    logger.info("Checking VEP missingness from global sample...")

    def get_missingness(vep_expr, version_label):
        """Calculate missingness for a VEP expression."""
        return sampled_ht.aggregate(
            hl.struct(
                vep_missing=hl.agg.fraction(hl.is_missing(vep_expr)),
                transcript_consequences_missing=hl.agg.fraction(
                    hl.is_missing(vep_expr.transcript_consequences)
                ),
                transcript_consequences_empty=hl.agg.fraction(
                    hl.is_defined(vep_expr.transcript_consequences)
                    & (hl.len(vep_expr.transcript_consequences) == 0)
                ),
                most_severe_consequence_missing=hl.agg.fraction(
                    hl.is_missing(vep_expr.most_severe_consequence)
                ),
                # Check for variants with no canonical transcript.
                no_canonical=hl.agg.fraction(
                    hl.is_missing(vep_expr.transcript_consequences)
                    | ~hl.any(
                        vep_expr.transcript_consequences.map(
                            lambda x: hl.coalesce(x.canonical, 0) == 1
                        )
                    )
                ),
                # Check for lof annotation presence.
                has_any_lof=hl.agg.fraction(
                    hl.is_defined(vep_expr.transcript_consequences)
                    & hl.any(
                        vep_expr.transcript_consequences.map(
                            lambda x: hl.is_defined(x.lof)
                        )
                    )
                ),
            )
        )

    missingness_105 = get_missingness(sampled_ht.vep105, "VEP105")
    missingness_115 = get_missingness(sampled_ht.vep115, "VEP115")

    result_105 = {k: v for k, v in missingness_105.items()}
    result_115 = {k: v for k, v in missingness_115.items()}

    for field in result_105.keys():
        logger.info(
            f"  VEP105 {field}: {result_105[field]:.4%}, "
            f"VEP115 {field}: {result_115[field]:.4%}"
        )

    return result_105, result_115


def compare_consequence_distributions_from_sample(
    sampled_ht: hl.Table,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Compare most_severe_consequence distributions using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Tuple of consequence counters for each version.
    """
    logger.info("Comparing most_severe_consequence distributions from global sample...")

    dist_105 = sampled_ht.aggregate(
        hl.agg.counter(sampled_ht.vep105.most_severe_consequence)
    )
    dist_115 = sampled_ht.aggregate(
        hl.agg.counter(sampled_ht.vep115.most_severe_consequence)
    )

    # Log differences.
    all_consequences = set(dist_105.keys()) | set(dist_115.keys())
    for csq in sorted(all_consequences):
        c105 = dist_105.get(csq, 0)
        c115 = dist_115.get(csq, 0)
        if c105 != c115:
            pct_diff = ((c115 - c105) / max(c105, 1)) * 100
            logger.info(f"  {csq}: VEP105={c105:,}, VEP115={c115:,} ({pct_diff:+.1f}%)")

    return dist_105, dist_115


def compare_biotype_distribution_from_sample(
    sampled_ht: hl.Table,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    Compare transcript biotype distributions using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Tuple of biotype counters for each version.
    """
    logger.info("Comparing transcript biotype distributions from global sample...")

    biotypes_105 = sampled_ht.aggregate(
        hl.agg.explode(
            lambda x: hl.agg.counter(x.biotype),
            sampled_ht.vep105.transcript_consequences,
        )
    )
    biotypes_115 = sampled_ht.aggregate(
        hl.agg.explode(
            lambda x: hl.agg.counter(x.biotype),
            sampled_ht.vep115.transcript_consequences,
        )
    )

    all_biotypes = set(biotypes_105.keys()) | set(biotypes_115.keys())
    logger.info(f"  Found {len(all_biotypes)} unique biotypes")

    # Log significant differences.
    for bt in sorted(all_biotypes):
        c105 = biotypes_105.get(bt, 0)
        c115 = biotypes_115.get(bt, 0)
        if c105 > 0 or c115 > 0:
            pct_diff = (
                ((c115 - c105) / max(c105, 1)) * 100 if c105 > 0 else float("inf")
            )
            if abs(pct_diff) > 10 or c105 == 0 or c115 == 0:
                logger.info(
                    f"    {bt}: VEP105={c105:,}, VEP115={c115:,} ({pct_diff:+.1f}%)"
                )

    return biotypes_105, biotypes_115


def check_chr18_centromere(ht115: hl.Table) -> Dict[str, int]:
    """
    Specifically check chr18 centromere region that was patched.

    The centromere of chr18 is approximately at 17.2-19.0 Mb.

    :param ht115: VEP 115 table (the patched one).
    :return: Dict with counts for centromere region.
    """
    logger.info("Checking chr18 centromere region (patch area)...")

    # Centromere region - approximate bounds.
    centromere_interval = hl.parse_locus_interval(
        "chr18:15000000-21000000", reference_genome="GRCh38"
    )

    ht_centro = ht115.filter(centromere_interval.contains(ht115.locus))

    stats = ht_centro.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            n_vep_missing=hl.agg.count_where(hl.is_missing(ht_centro.vep)),
            n_transcript_missing=hl.agg.count_where(
                hl.is_defined(ht_centro.vep)
                & hl.is_missing(ht_centro.vep.transcript_consequences)
            ),
            n_transcript_empty=hl.agg.count_where(
                hl.is_defined(ht_centro.vep)
                & hl.is_defined(ht_centro.vep.transcript_consequences)
                & (hl.len(ht_centro.vep.transcript_consequences) == 0)
            ),
        )
    )

    logger.info(f"  chr18 centromere variants: {stats.n_variants:,}")
    logger.info(f"  VEP missing: {stats.n_vep_missing:,}")
    logger.info(f"  transcript_consequences missing: {stats.n_transcript_missing:,}")
    logger.info(f"  transcript_consequences empty: {stats.n_transcript_empty:,}")

    return dict(stats)


def spot_check_shared_annotations(
    ht105: hl.Table, ht115: hl.Table, n_samples: int = 10_000
) -> Dict[str, float]:
    """
    Spot check that shared annotation fields have consistent values.

    Randomly samples variants and compares core fields that shouldn't change.

    :param ht105: VEP 105 table.
    :param ht115: VEP 115 table.
    :param n_samples: Number of variants to sample.
    :return: Dict of field -> agreement rate.
    """
    logger.info(f"Spot checking ~{n_samples:,} variants for annotation consistency...")

    # Join tables on key.
    ht = ht105.select(vep105=ht105.vep)
    ht = ht.annotate(vep115=ht115[ht.key].vep)

    # Filter to variants present in both with non-empty consequences.
    ht = ht.filter(
        hl.is_defined(ht.vep115)
        & (hl.len(ht.vep105.transcript_consequences) > 0)
        & (hl.len(ht.vep115.transcript_consequences) > 0)
    )

    # Random sample for better coverage across genome.
    ht = ht.sample(0.0001, seed=42).head(n_samples)

    # Helper to get canonical transcript gene_id (first one if multiple).
    def get_canonical_gene_id(vep):
        canonical = vep.transcript_consequences.filter(
            lambda x: hl.coalesce(x.canonical, 0) == 1
        )
        return hl.or_missing(hl.len(canonical) > 0, canonical[0].gene_id)

    ht = ht.annotate(
        canonical_gene_105=get_canonical_gene_id(ht.vep105),
        canonical_gene_115=get_canonical_gene_id(ht.vep115),
    )

    # Compare most_severe_consequence.
    agreement = ht.aggregate(
        hl.struct(
            most_severe_csq_match=hl.agg.fraction(
                ht.vep105.most_severe_consequence == ht.vep115.most_severe_consequence
            ),
            # Compare canonical transcript gene_id.
            canonical_gene_match=hl.agg.fraction(
                ht.canonical_gene_105 == ht.canonical_gene_115
            ),
            # Compare number of transcript consequences.
            n_transcripts_match=hl.agg.fraction(
                hl.len(ht.vep105.transcript_consequences)
                == hl.len(ht.vep115.transcript_consequences)
            ),
        )
    )

    for field, rate in agreement.items():
        logger.info(f"  {field}: {rate:.2%}")

    return dict(agreement)


def run_gene_coverage_check(
    ht105: hl.Table, ht115: hl.Table, interval_ht: hl.Table
) -> Tuple[hl.Table, hl.Table]:
    """
    Run gene-level coverage check using gnomad utility.

    :param ht105: VEP 105 table.
    :param ht115: VEP 115 table.
    :param interval_ht: Ensembl gene intervals.
    :return: Tuple of annotated interval tables for each version.
    """
    logger.info("Running gene-level VEP coverage check...")

    interval_105 = count_vep_annotated_variants_per_interval(ht105, interval_ht)
    interval_115 = count_vep_annotated_variants_per_interval(ht115, interval_ht)

    return interval_105, interval_115


def compare_loftee_annotations_from_sample(sampled_ht: hl.Table) -> Dict[str, Any]:
    """
    Compare LOFTEE (LoF) annotations between versions using a pre-sampled table.

    LOFTEE is critical for variant interpretation - should be consistent.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict with comparison stats.
    """
    logger.info("Comparing LOFTEE annotations from global sample...")

    ht = sampled_ht

    # Get canonical transcript LoF for each.
    def get_canonical_lof(vep):
        canonical = vep.transcript_consequences.filter(
            lambda x: hl.coalesce(x.canonical, 0) == 1
        )
        return hl.or_missing(
            hl.len(canonical) > 0,
            hl.struct(
                lof=canonical[0].lof,
                lof_filter=canonical[0].lof_filter,
                lof_flags=canonical[0].lof_flags,
            ),
        )

    ht = ht.annotate(
        lof105=get_canonical_lof(ht.vep105),
        lof115=get_canonical_lof(ht.vep115),
    )

    # Filter to variants with LoF in at least one version.
    ht = ht.filter(
        (hl.is_defined(ht.lof105) & hl.is_defined(ht.lof105.lof))
        | (hl.is_defined(ht.lof115) & hl.is_defined(ht.lof115.lof))
    )

    stats = ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            lof_match=hl.agg.fraction(ht.lof105.lof == ht.lof115.lof),
            lof_filter_match=hl.agg.fraction(
                ht.lof105.lof_filter == ht.lof115.lof_filter
            ),
            # Count LoF categories.
            lof105_hc=hl.agg.count_where(ht.lof105.lof == "HC"),
            lof115_hc=hl.agg.count_where(ht.lof115.lof == "HC"),
            lof105_lc=hl.agg.count_where(ht.lof105.lof == "LC"),
            lof115_lc=hl.agg.count_where(ht.lof115.lof == "LC"),
            # Discordant cases.
            hc_to_lc=hl.agg.count_where(
                (ht.lof105.lof == "HC") & (ht.lof115.lof == "LC")
            ),
            lc_to_hc=hl.agg.count_where(
                (ht.lof105.lof == "LC") & (ht.lof115.lof == "HC")
            ),
            gained_lof=hl.agg.count_where(
                hl.is_missing(ht.lof105.lof) & hl.is_defined(ht.lof115.lof)
            ),
            lost_lof=hl.agg.count_where(
                hl.is_defined(ht.lof105.lof) & hl.is_missing(ht.lof115.lof)
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} variants with LoF annotations")
    logger.info(f"  LoF category match rate: {stats.lof_match:.2%}")
    logger.info(f"  LoF filter match rate: {stats.lof_filter_match:.2%}")
    logger.info(f"  VEP105 HC: {stats.lof105_hc:,}, VEP115 HC: {stats.lof115_hc:,}")
    logger.info(
        f"  HC->LC changes: {stats.hc_to_lc:,}, LC->HC changes: {stats.lc_to_hc:,}"
    )
    logger.info(f"  Gained LoF: {stats.gained_lof:,}, Lost LoF: {stats.lost_lof:,}")

    return dict(stats)


def compare_sift_polyphen_from_sample(sampled_ht: hl.Table) -> Dict[str, Any]:
    """
    Compare SIFT and PolyPhen predictions between versions using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict with agreement rates.
    """
    logger.info("Comparing SIFT/PolyPhen predictions from global sample...")

    ht = sampled_ht

    # Filter to missense variants.
    ht = ht.filter(
        (ht.vep105.most_severe_consequence == "missense_variant")
        | (ht.vep115.most_severe_consequence == "missense_variant")
    )

    def get_canonical_predictions(vep):
        canonical = vep.transcript_consequences.filter(
            lambda x: hl.coalesce(x.canonical, 0) == 1
        )
        return hl.or_missing(
            hl.len(canonical) > 0,
            hl.struct(
                sift_prediction=canonical[0].sift_prediction,
                sift_score=canonical[0].sift_score,
                polyphen_prediction=canonical[0].polyphen_prediction,
                polyphen_score=canonical[0].polyphen_score,
            ),
        )

    ht = ht.annotate(
        pred105=get_canonical_predictions(ht.vep105),
        pred115=get_canonical_predictions(ht.vep115),
    )

    ht = ht.filter(hl.is_defined(ht.pred105) & hl.is_defined(ht.pred115))

    stats = ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            sift_pred_match=hl.agg.fraction(
                ht.pred105.sift_prediction == ht.pred115.sift_prediction
            ),
            polyphen_pred_match=hl.agg.fraction(
                ht.pred105.polyphen_prediction == ht.pred115.polyphen_prediction
            ),
            # Score correlation (for non-missing).
            sift_score_corr=hl.agg.filter(
                hl.is_defined(ht.pred105.sift_score)
                & hl.is_defined(ht.pred115.sift_score),
                hl.agg.corr(ht.pred105.sift_score, ht.pred115.sift_score),
            ),
            polyphen_score_corr=hl.agg.filter(
                hl.is_defined(ht.pred105.polyphen_score)
                & hl.is_defined(ht.pred115.polyphen_score),
                hl.agg.corr(ht.pred105.polyphen_score, ht.pred115.polyphen_score),
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} missense variants")
    logger.info(f"  SIFT prediction match: {stats.sift_pred_match:.2%}")
    logger.info(f"  PolyPhen prediction match: {stats.polyphen_pred_match:.2%}")
    logger.info(f"  SIFT score correlation: {stats.sift_score_corr:.4f}")
    logger.info(f"  PolyPhen score correlation: {stats.polyphen_score_corr:.4f}")

    return dict(stats)


def compare_mane_annotations_from_sample(sampled_ht: hl.Table) -> Dict[str, Any]:
    """
    Check MANE transcript annotation coverage in VEP 115 using a pre-sampled table.

    MANE (Matched Annotation from NCBI and EBI) is important for clinical interpretation.
    VEP 115 has improved MANE support with the new 'mane' array field.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict with MANE coverage stats.
    """
    logger.info("Checking MANE transcript coverage (VEP 115) from global sample...")

    # Filter to variants with protein_coding transcripts in VEP 115 (where
    # MANE would be).
    ht = sampled_ht.filter(
        (
            hl.any(
                sampled_ht.vep115.transcript_consequences.map(
                    lambda x: x.biotype == "protein_coding"
                )
            )
        )
        | (
            hl.any(
                sampled_ht.vep105.transcript_consequences.map(
                    lambda x: x.biotype == "protein_coding"
                )
            )
        )
    )

    # MANE annotations are a VEP 115 feature, so only VEP 115 stats are meaningful
    # VEP 105 won't have these fields, so we just show VEP 115 stats
    stats = ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            # VEP 115 stats only (MANE is a VEP 115 feature)
            vep115_mane_select=hl.agg.fraction(
                hl.any(
                    ht.vep115.transcript_consequences.map(
                        lambda x: hl.is_defined(x.mane_select)
                    )
                )
            ),
            vep115_mane_plus_clinical=hl.agg.fraction(
                hl.any(
                    ht.vep115.transcript_consequences.map(
                        lambda x: hl.is_defined(x.mane_plus_clinical)
                    )
                )
            ),
            vep115_mane_array=hl.agg.fraction(
                hl.any(
                    ht.vep115.transcript_consequences.map(
                        lambda x: hl.is_defined(x.mane) & (hl.len(x.mane) > 0)
                    )
                )
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} protein-coding variants")
    logger.info(f"  VEP115 - MANE Select: {stats.vep115_mane_select:.2%}")
    logger.info(f"  VEP115 - MANE Plus Clinical: {stats.vep115_mane_plus_clinical:.2%}")
    logger.info(f"  VEP115 - MANE array: {stats.vep115_mane_array:.2%}")
    logger.info(
        "  Note: MANE annotations are a VEP 115 feature (not available in VEP 105)"
    )

    return dict(stats)


def check_colocated_variants(
    ht105: hl.Table, ht115: hl.Table, n_samples: int = 50_000
) -> Dict[str, Any]:
    """
    Compare colocated_variants (dbSNP rsIDs, ClinVar) between versions.

    :param ht105: VEP 105 table.
    :param ht115: VEP 115 table.
    :param n_samples: Number of variants to sample.
    :return: Dict with comparison stats.
    """
    logger.info("Comparing colocated_variants annotations...")

    ht = ht105.select(vep105=ht105.vep)
    ht = ht.annotate(vep115=ht115[ht.key].vep)
    ht = ht.filter(hl.is_defined(ht.vep115))

    # Random sample for better coverage.
    ht = ht.sample(0.0001, seed=42).head(n_samples)

    stats = ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            # Has colocated variants.
            has_colocated_105=hl.agg.fraction(
                hl.is_defined(ht.vep105.colocated_variants)
                & (hl.len(ht.vep105.colocated_variants) > 0)
            ),
            has_colocated_115=hl.agg.fraction(
                hl.is_defined(ht.vep115.colocated_variants)
                & (hl.len(ht.vep115.colocated_variants) > 0)
            ),
            # Has rsID.
            has_rsid_105=hl.agg.fraction(
                hl.any(
                    ht.vep105.colocated_variants.map(
                        lambda x: hl.is_defined(x.id) & x.id.startswith("rs")
                    )
                )
            ),
            has_rsid_115=hl.agg.fraction(
                hl.any(
                    ht.vep115.colocated_variants.map(
                        lambda x: hl.is_defined(x.id) & x.id.startswith("rs")
                    )
                )
            ),
            # Has ClinVar significance.
            has_clin_sig_105=hl.agg.fraction(
                hl.any(
                    ht.vep105.colocated_variants.map(
                        lambda x: hl.is_defined(x.clin_sig) & (hl.len(x.clin_sig) > 0)
                    )
                )
            ),
            has_clin_sig_115=hl.agg.fraction(
                hl.any(
                    ht.vep115.colocated_variants.map(
                        lambda x: hl.is_defined(x.clin_sig) & (hl.len(x.clin_sig) > 0)
                    )
                )
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} variants")
    logger.info(
        f"  Has colocated variants: VEP105={stats.has_colocated_105:.2%}, "
        f"VEP115={stats.has_colocated_115:.2%}"
    )
    logger.info(
        f"  Has rsID: VEP105={stats.has_rsid_105:.2%}, VEP115={stats.has_rsid_115:.2%}"
    )
    logger.info(
        f"  Has ClinVar: VEP105={stats.has_clin_sig_105:.2%}, "
        f"VEP115={stats.has_clin_sig_115:.2%}"
    )

    return dict(stats)


def spot_check_shared_annotations_from_sample(sampled_ht: hl.Table) -> Dict[str, float]:
    """
    Spot check that shared annotation fields have consistent values using pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict of field -> agreement rate.
    """
    logger.info(
        f"Spot checking ~{sampled_ht.count():,} variants for annotation consistency..."
    )

    # Helper to get canonical transcript gene_id (first one if multiple).
    def get_canonical_gene_id(vep):
        canonical = vep.transcript_consequences.filter(
            lambda x: hl.coalesce(x.canonical, 0) == 1
        )
        return hl.or_missing(hl.len(canonical) > 0, canonical[0].gene_id)

    ht = sampled_ht.annotate(
        canonical_gene_105=get_canonical_gene_id(sampled_ht.vep105),
        canonical_gene_115=get_canonical_gene_id(sampled_ht.vep115),
    )

    # Compare most_severe_consequence.
    agreement = ht.aggregate(
        hl.struct(
            most_severe_csq_match=hl.agg.fraction(
                ht.vep105.most_severe_consequence == ht.vep115.most_severe_consequence
            ),
            # Compare canonical transcript gene_id.
            canonical_gene_match=hl.agg.fraction(
                ht.canonical_gene_105 == ht.canonical_gene_115
            ),
            # Compare number of transcript consequences.
            n_transcripts_match=hl.agg.fraction(
                hl.len(ht.vep105.transcript_consequences)
                == hl.len(ht.vep115.transcript_consequences)
            ),
        )
    )

    for field, rate in agreement.items():
        logger.info(f"  {field}: {rate:.2%}")

    return dict(agreement)


def check_colocated_variants_from_sample(sampled_ht: hl.Table) -> Dict[str, Any]:
    """
    Compare colocated_variants using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict with comparison stats.
    """
    logger.info("Comparing colocated_variants annotations from global sample...")

    stats = sampled_ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            # Has colocated variants.
            has_colocated_105=hl.agg.fraction(
                hl.is_defined(sampled_ht.vep105.colocated_variants)
                & (hl.len(sampled_ht.vep105.colocated_variants) > 0)
            ),
            has_colocated_115=hl.agg.fraction(
                hl.is_defined(sampled_ht.vep115.colocated_variants)
                & (hl.len(sampled_ht.vep115.colocated_variants) > 0)
            ),
            # Has rsID.
            has_rsid_105=hl.agg.fraction(
                hl.any(
                    sampled_ht.vep105.colocated_variants.map(
                        lambda x: hl.is_defined(x.id) & x.id.startswith("rs")
                    )
                )
            ),
            has_rsid_115=hl.agg.fraction(
                hl.any(
                    sampled_ht.vep115.colocated_variants.map(
                        lambda x: hl.is_defined(x.id) & x.id.startswith("rs")
                    )
                )
            ),
            # Has ClinVar significance.
            has_clin_sig_105=hl.agg.fraction(
                hl.any(
                    sampled_ht.vep105.colocated_variants.map(
                        lambda x: hl.is_defined(x.clin_sig) & (hl.len(x.clin_sig) > 0)
                    )
                )
            ),
            has_clin_sig_115=hl.agg.fraction(
                hl.any(
                    sampled_ht.vep115.colocated_variants.map(
                        lambda x: hl.is_defined(x.clin_sig) & (hl.len(x.clin_sig) > 0)
                    )
                )
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} variants")
    logger.info(
        f"  Has colocated variants: VEP105={stats.has_colocated_105:.2%}, "
        f"VEP115={stats.has_colocated_115:.2%}"
    )
    logger.info(
        f"  Has rsID: VEP105={stats.has_rsid_105:.2%}, VEP115={stats.has_rsid_115:.2%}"
    )
    logger.info(
        f"  Has ClinVar: VEP105={stats.has_clin_sig_105:.2%}, "
        f"VEP115={stats.has_clin_sig_115:.2%}"
    )

    return dict(stats)


def compare_regulatory_consequences_from_sample(sampled_ht: hl.Table) -> Dict[str, Any]:
    """
    Compare regulatory consequences using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict with comparison stats.
    """
    logger.info("Comparing regulatory feature consequences from global sample...")

    stats = sampled_ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            has_regulatory_105=hl.agg.fraction(
                hl.is_defined(sampled_ht.vep105.regulatory_feature_consequences)
                & (hl.len(sampled_ht.vep105.regulatory_feature_consequences) > 0)
            ),
            has_regulatory_115=hl.agg.fraction(
                hl.is_defined(sampled_ht.vep115.regulatory_feature_consequences)
                & (hl.len(sampled_ht.vep115.regulatory_feature_consequences) > 0)
            ),
            has_motif_105=hl.agg.fraction(
                hl.is_defined(sampled_ht.vep105.motif_feature_consequences)
                & (hl.len(sampled_ht.vep105.motif_feature_consequences) > 0)
            ),
            has_motif_115=hl.agg.fraction(
                hl.is_defined(sampled_ht.vep115.motif_feature_consequences)
                & (hl.len(sampled_ht.vep115.motif_feature_consequences) > 0)
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} variants")
    logger.info(
        f"  Has regulatory: VEP105={stats.has_regulatory_105:.2%}, "
        f"VEP115={stats.has_regulatory_115:.2%}"
    )
    logger.info(
        f"  Has motif: VEP105={stats.has_motif_105:.2%}, "
        f"VEP115={stats.has_motif_115:.2%}"
    )

    return dict(stats)


def compare_transcript_counts_distribution_from_sample(
    sampled_ht: hl.Table,
) -> Dict[str, Any]:
    """
    Compare transcript count distributions using a pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    :return: Dict with distribution stats.
    """
    logger.info("Comparing transcript count distributions from global sample...")

    ht = sampled_ht.annotate(
        n_tx_105=hl.len(sampled_ht.vep105.transcript_consequences),
        n_tx_115=hl.len(sampled_ht.vep115.transcript_consequences),
    )

    stats = ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            mean_tx_105=hl.agg.mean(ht.n_tx_105),
            mean_tx_115=hl.agg.mean(ht.n_tx_115),
            median_tx_105=hl.agg.approx_quantiles(ht.n_tx_105, 0.5),
            median_tx_115=hl.agg.approx_quantiles(ht.n_tx_115, 0.5),
            max_tx_105=hl.agg.max(ht.n_tx_105),
            max_tx_115=hl.agg.max(ht.n_tx_115),
            # Correlation between counts.
            tx_count_corr=hl.agg.corr(hl.float64(ht.n_tx_105), hl.float64(ht.n_tx_115)),
            # What fraction have MORE transcripts in 115?
            more_in_115=hl.agg.fraction(ht.n_tx_115 > ht.n_tx_105),
            fewer_in_115=hl.agg.fraction(ht.n_tx_115 < ht.n_tx_105),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} variants")
    logger.info(
        f"  Mean transcripts: VEP105={stats.mean_tx_105:.1f}, "
        f"VEP115={stats.mean_tx_115:.1f}"
    )
    logger.info(
        f"  Median transcripts: VEP105={stats.median_tx_105:.0f}, "
        f"VEP115={stats.median_tx_115:.0f}"
    )
    logger.info(
        f"  Max transcripts: VEP105={stats.max_tx_105:,}, "
        f"VEP115={stats.max_tx_115:,}"
    )
    logger.info(f"  Transcript count correlation: {stats.tx_count_corr:.4f}")
    logger.info(
        f"  More transcripts in 115: {stats.more_in_115:.2%}, "
        f"Fewer: {stats.fewer_in_115:.2%}"
    )

    return dict(stats)


def check_specific_genes(
    ht105: hl.Table, ht115: hl.Table, genes: List[str] = None
) -> Dict[str, Dict]:
    """
    Spot check specific clinically important genes.

    :param ht105: VEP 105 table.
    :param ht115: VEP 115 table.
    :param genes: List of gene symbols to check.
    :return: Dict with per-gene comparison stats.
    """
    if genes is None:
        # Clinically important genes to spot check.
        genes = ["BRCA1", "BRCA2", "TP53", "CFTR", "TTN", "APOB", "LDLR"]

    logger.info(f"Spot checking specific genes: {genes}")

    results = {}
    for gene in genes:
        # Filter to variants annotated with this gene.
        ht = ht105.filter(
            hl.any(
                ht105.vep.transcript_consequences.map(lambda x: x.gene_symbol == gene)
            )
        )
        ht = ht.select(vep105=ht.vep)
        ht = ht.annotate(vep115=ht115[ht.key].vep)
        ht = ht.filter(hl.is_defined(ht.vep115))

        stats = ht.aggregate(
            hl.struct(
                n_variants=hl.agg.count(),
                csq_match=hl.agg.fraction(
                    ht.vep105.most_severe_consequence
                    == ht.vep115.most_severe_consequence
                ),
            )
        )

        results[gene] = dict(stats)
        logger.info(
            f"  {gene}: {stats.n_variants:,} variants, "
            f"consequence match: {stats.csq_match:.2%}"
        )

    return results


def debug_lof_filter_differences_from_sample(sampled_ht: hl.Table) -> None:
    """
    Debug why lof_filter match rate is so low despite lof category matching using pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    """
    logger.info("Debugging LoF filter differences from global sample...")

    ht = sampled_ht

    def get_canonical_lof_info(vep):
        canonical = vep.transcript_consequences.filter(
            lambda x: hl.coalesce(x.canonical, 0) == 1
        )
        return hl.or_missing(
            hl.len(canonical) > 0,
            hl.struct(
                lof=canonical[0].lof,
                lof_filter=canonical[0].lof_filter,
            ),
        )

    ht = ht.annotate(
        lof105=get_canonical_lof_info(ht.vep105),
        lof115=get_canonical_lof_info(ht.vep115),
    )

    # Filter to variants with LoF in at least one version.
    ht = ht.filter(
        (hl.is_defined(ht.lof105) & hl.is_defined(ht.lof105.lof))
        | (hl.is_defined(ht.lof115) & hl.is_defined(ht.lof115.lof))
    )

    # Overall LoF category transitions.
    lof_transitions = ht.aggregate(
        hl.agg.counter(hl.struct(from_=ht.lof105.lof, to=ht.lof115.lof))
    )

    logger.info("  LoF category transitions:")
    for combo, count in sorted(lof_transitions.items(), key=lambda x: -x[1])[:10]:
        logger.info(f"    '{combo.from_}' -> '{combo.to}': {count:,}")

    # Filter string transitions (for variants with LoF in both).
    ht_both = ht.filter(hl.is_defined(ht.lof105.lof) & hl.is_defined(ht.lof115.lof))

    filter_combos = ht_both.aggregate(
        hl.agg.counter(
            hl.struct(
                filter105=ht_both.lof105.lof_filter, filter115=ht_both.lof115.lof_filter
            )
        )
    )

    logger.info("  LoF filter string transitions (when both have LoF):")
    for combo, count in sorted(filter_combos.items(), key=lambda x: -x[1])[:15]:
        match = "✓" if combo.filter105 == combo.filter115 else "✗"
        logger.info(
            f"    {match} '{combo.filter105}' -> '{combo.filter115}': {count:,}"
        )


def debug_mane_annotations(ht115: hl.Table, n_samples: int = 1000) -> None:
    """
    Debug why MANE annotations show 0%.

    :param ht115: VEP 115 table.
    :param n_samples: Number of variants to sample.
    """
    logger.info("Debugging MANE annotations in VEP 115...")

    # Filter to variants with protein_coding transcripts (where MANE would be).
    ht = ht115.filter(
        hl.any(
            ht115.vep.transcript_consequences.map(
                lambda x: x.biotype == "protein_coding"
            )
        )
    )

    # Sample randomly, not from beginning.
    ht = ht.sample(0.001, seed=42).head(n_samples)

    # Check what MANE fields actually contain.
    stats = ht.aggregate(
        hl.struct(
            n_variants=hl.agg.count(),
            # Check canonical protein_coding transcripts specifically.
            has_any_mane_select=hl.agg.fraction(
                hl.any(
                    ht.vep.transcript_consequences.map(
                        lambda x: hl.is_defined(x.mane_select)
                        & (x.biotype == "protein_coding")
                    )
                )
            ),
            has_any_mane_plus=hl.agg.fraction(
                hl.any(
                    ht.vep.transcript_consequences.map(
                        lambda x: hl.is_defined(x.mane_plus_clinical)
                        & (x.biotype == "protein_coding")
                    )
                )
            ),
            # Sample actual MANE values.
            sample_mane_select=hl.agg.filter(
                hl.any(
                    ht.vep.transcript_consequences.map(
                        lambda x: hl.is_defined(x.mane_select)
                    )
                ),
                hl.agg.take(
                    ht.vep.transcript_consequences.filter(
                        lambda x: hl.is_defined(x.mane_select)
                    )[0].mane_select,
                    3,
                ),
            ),
        )
    )

    logger.info(f"  Sampled {stats.n_variants:,} protein-coding variants")
    logger.info(f"  Has any mane_select: {stats.has_any_mane_select:.2%}")
    logger.info(f"  Has any mane_plus_clinical: {stats.has_any_mane_plus:.2%}")
    logger.info(f"  Sample mane_select values: {stats.sample_mane_select}")

    # Also check the VEP config to see if MANE was enabled.
    if "vep_config" in ht115.globals.dtype.fields:
        vep_config = hl.eval(ht115.vep_config)
        logger.info(
            f"  VEP config contains 'mane': {'mane' in str(vep_config).lower()}"
        )
        logger.info(
            f"  VEP config contains 'regulatory': {'regulatory' in str(vep_config).lower()}"
        )
        logger.info(
            f"  VEP config contains 'motif': {'motif' in str(vep_config).lower()}"
        )


def debug_regulatory_motif(ht115: hl.Table, n_samples: int = 1000) -> None:
    """
    Debug regulatory and motif feature annotations in VEP 115.

    Aggregates per-contig in a single pass for efficiency, logging progress by chromosome.

    :param ht115: VEP 115 table.
    :param n_samples: Number of variants to sample for examples.
    """
    logger.info("Debugging regulatory and motif features in VEP 115...")
    logger.info("Aggregating per-contig (single pass through table)...")

    # Define boolean expressions for has_motif and has_regulatory.
    has_motif = hl.is_defined(ht115.vep.motif_feature_consequences) & (
        hl.len(ht115.vep.motif_feature_consequences) > 0
    )
    has_regulatory = hl.is_defined(ht115.vep.regulatory_feature_consequences) & (
        hl.len(ht115.vep.regulatory_feature_consequences) > 0
    )

    # Single-pass aggregation grouped by contig.
    per_contig_stats = ht115.aggregate(
        hl.agg.group_by(
            ht115.locus.contig,
            hl.struct(
                n_variants=hl.agg.count(),
                n_with_motif=hl.agg.count_where(has_motif),
                n_with_regulatory=hl.agg.count_where(has_regulatory),
            ),
        )
    )

    # Log per-contig results and accumulate totals.
    logger.info("=" * 60)
    logger.info("PER-CONTIG MOTIF/REGULATORY RESULTS")
    logger.info("=" * 60)
    logger.info(f"{'Contig':<10} {'Variants':>15} {'Motif':>12} {'Regulatory':>12}")
    logger.info("-" * 60)

    total_variants = 0
    total_motif = 0
    total_regulatory = 0

    # Sort contigs for readable output.
    def contig_sort_key(c):
        if c.startswith("chr"):
            c = c[3:]
        if c.isdigit():
            return (0, int(c))
        return (1, c)

    for contig in sorted(per_contig_stats.keys(), key=contig_sort_key):
        stats = per_contig_stats[contig]
        total_variants += stats.n_variants
        total_motif += stats.n_with_motif
        total_regulatory += stats.n_with_regulatory

        # Only log if there's something interesting (motif or regulatory found).
        if stats.n_with_motif > 0 or stats.n_with_regulatory > 0:
            logger.info(
                f"{contig:<10} {stats.n_variants:>15,} {stats.n_with_motif:>12,} "
                f"{stats.n_with_regulatory:>12,}"
            )
        else:
            # Log abbreviated for contigs with no hits.
            logger.info(f"{contig:<10} {stats.n_variants:>15,} {'0':>12} {'0':>12}")

    logger.info("-" * 60)
    logger.info(
        f"{'TOTAL':<10} {total_variants:>15,} {total_motif:>12,} {total_regulatory:>12,}"
    )
    logger.info("=" * 60)

    # Summary.
    logger.info("SUMMARY")
    logger.info("=" * 60)
    logger.info(f"  Total variants: {total_variants:,}")
    logger.info(f"  ANY variant has motif features: {total_motif > 0}")
    logger.info(f"  ANY variant has regulatory features: {total_regulatory > 0}")

    if total_motif > 0:
        pct_motif = total_motif / total_variants * 100
        logger.info(f"  Motif count: {total_motif:,} ({pct_motif:.4f}%)")
    else:
        logger.warning("  MOTIF FIELD IS EMPTY - no variants have motif features!")

    if total_regulatory > 0:
        pct_reg = total_regulatory / total_variants * 100
        logger.info(f"  Regulatory count: {total_regulatory:,} ({pct_reg:.4f}%)")
    else:
        logger.warning(
            "  REGULATORY FIELD IS EMPTY - no variants have regulatory features!"
        )

    # Grab sample values if any exist (quick filter + take).
    if total_motif > 0:
        logger.info("  Fetching sample motif values...")
        sample_motif = ht115.filter(has_motif).head(3).collect()
        for i, row in enumerate(sample_motif):
            logger.info(
                f"    Motif example {i + 1}: {row.vep.motif_feature_consequences}"
            )

    if total_regulatory > 0:
        logger.info("  Fetching sample regulatory values...")
        sample_reg = ht115.filter(has_regulatory).head(3).collect()
        for i, row in enumerate(sample_reg):
            logger.info(
                f"    Regulatory example {i + 1}: {row.vep.regulatory_feature_consequences}"
            )

    # Check VEP config for regulatory/motif plugins.
    if "vep_config" in ht115.globals.dtype.fields:
        vep_config = hl.eval(ht115.vep_config)
        logger.info(
            f"  VEP config contains 'regulatory': {'regulatory' in str(vep_config).lower()}"
        )
        logger.info(
            f"  VEP config contains 'motif': {'motif' in str(vep_config).lower()}"
        )


def check_sift_polyphen_category_changes_from_sample(sampled_ht: hl.Table) -> None:
    """
    Debug SIFT/PolyPhen category changes despite score correlation using pre-sampled table.

    :param sampled_ht: Pre-sampled and joined table with vep105 and vep115 annotations.
    """
    logger.info("Analyzing SIFT/PolyPhen category transitions from global sample...")

    ht = sampled_ht.filter(
        (ht.vep105.most_severe_consequence == "missense_variant")
        | (ht.vep115.most_severe_consequence == "missense_variant")
    )

    def get_canonical_predictions(vep):
        canonical = vep.transcript_consequences.filter(
            lambda x: hl.coalesce(x.canonical, 0) == 1
        )
        return hl.or_missing(
            hl.len(canonical) > 0,
            hl.struct(
                sift=canonical[0].sift_prediction,
                polyphen=canonical[0].polyphen_prediction,
            ),
        )

    ht = ht.annotate(
        pred105=get_canonical_predictions(ht.vep105),
        pred115=get_canonical_predictions(ht.vep115),
    )
    ht = ht.filter(hl.is_defined(ht.pred105) & hl.is_defined(ht.pred115))

    # SIFT transitions.
    sift_transitions = ht.aggregate(
        hl.agg.counter(hl.struct(from_=ht.pred105.sift, to=ht.pred115.sift))
    )

    logger.info("  SIFT category transitions:")
    for combo, count in sorted(sift_transitions.items(), key=lambda x: -x[1])[:10]:
        if combo.from_ != combo.to:
            logger.info(f"    '{combo.from_}' -> '{combo.to}': {count:,}")

    # PolyPhen transitions.
    polyphen_transitions = ht.aggregate(
        hl.agg.counter(hl.struct(from_=ht.pred105.polyphen, to=ht.pred115.polyphen))
    )

    logger.info("  PolyPhen category transitions:")
    for combo, count in sorted(polyphen_transitions.items(), key=lambda x: -x[1])[:10]:
        if combo.from_ != combo.to:
            logger.info(f"    '{combo.from_}' -> '{combo.to}': {count:,}")


# Available check groups.
AVAILABLE_CHECKS = {
    "schema": "Schema comparison between VEP versions",
    "counts": "Variant count comparison per chromosome",
    "centromere": "Chr18 centromere patch region check",
    "missingness": "VEP field missingness rates",
    "consequences": "Most severe consequence distribution",
    "consistency": "Annotation consistency spot check",
    "loftee": "LOFTEE annotation comparison",
    "sift_polyphen": "SIFT/PolyPhen prediction comparison",
    "mane": "MANE transcript coverage (VEP 115)",
    "biotype": "Transcript biotype distribution",
    "colocated": "Colocated variants (dbSNP, ClinVar)",
    "regulatory": "Regulatory/motif consequences",
    "genes": "Specific gene spot checks (slow)",
    "gene_coverage": "Gene-level coverage (expensive)",
}

DEBUG_CHECKS = {
    "debug_tx_counts": "Transcript count distribution analysis",
    "debug_lof_filter": "LoF filter string difference analysis",
    "debug_mane": "MANE field value inspection",
    "debug_regulatory_motif": "Regulatory/motif feature inspection",
    "debug_sift_polyphen": "SIFT/PolyPhen category transitions",
}


def get_parser() -> argparse.ArgumentParser:
    """Create argument parser."""
    parser = argparse.ArgumentParser(
        description="Compare VEP 105 vs VEP 115 annotations on gnomAD context HT.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--vep105-path",
        help="Path to VEP 105 table. Default: gnomad public vep_context resource.",
        default=None,
    )
    parser.add_argument(
        "--vep115-path",
        help="Path to VEP 115 table.",
        default="gs://gnomad-tmp-30day/completed_vep_context.ht",
    )
    parser.add_argument(
        "--checks",
        nargs="+",
        choices=list(AVAILABLE_CHECKS.keys()) + list(DEBUG_CHECKS.keys()),
        help="Specific checks to run. Default: all main checks.",
    )
    parser.add_argument(
        "--debug-only",
        action="store_true",
        help="Run only debug checks (faster, for investigating specific issues).",
    )
    parser.add_argument(
        "--include-debug",
        action="store_true",
        help="Include debug checks when running main checks.",
    )
    parser.add_argument(
        "--list-checks",
        action="store_true",
        help="List all available checks and exit.",
    )
    parser.add_argument(
        "--tmp-dir",
        default="gs://gnomad-tmp-30day/vep_comparison",
        help="Temporary directory for Hail.",
    )
    parser.add_argument(
        "--no-global-sample",
        action="store_true",
        help=(
            "Skip creating the global sample. This speeds up startup but disables "
            "checks that require paired comparisons (missingness, consequences, "
            "consistency, loftee, sift_polyphen, biotype, colocated, regulatory, "
            "debug_lof_filter, debug_sift_polyphen)."
        ),
    )

    return parser


CHECKS_REQUIRING_GLOBAL_SAMPLE = {
    "missingness",
    "consequences",
    "consistency",
    "loftee",
    "sift_polyphen",
    "biotype",
    "colocated",
    "regulatory",
    "debug_lof_filter",
    "debug_sift_polyphen",
}


def run_check(
    name: str,
    ht105: hl.Table,
    ht115: hl.Table,
    global_sample: Optional[hl.Table],
    results: dict,
) -> None:
    """Run a single check by name."""
    # Skip checks that require global_sample if it's not available.
    if name in CHECKS_REQUIRING_GLOBAL_SAMPLE and global_sample is None:
        logger.warning(
            f"SKIPPING '{name}' check - requires global sample (use without "
            "--no-global-sample to enable)"
        )
        return

    logger.info("=" * 60)

    if name == "schema":
        logger.info("SCHEMA COMPARISON")
        logger.info("=" * 60)
        results["schema_diff"] = compare_schemas(ht105, ht115)

    elif name == "counts":
        logger.info("VARIANT COUNT COMPARISON")
        logger.info("=" * 60)
        counts_105, counts_115, count_diffs = compare_variant_counts(ht105, ht115)
        results["count_diffs"] = count_diffs

    elif name == "centromere":
        logger.info("CHR18 CENTROMERE CHECK (PATCH REGION)")
        logger.info("=" * 60)
        results["centro_stats"] = check_chr18_centromere(ht115)

    elif name == "missingness":
        logger.info("MISSINGNESS CHECK")
        logger.info("=" * 60)
        results["miss_105"], results["miss_115"] = check_vep_missingness_from_sample(
            global_sample
        )

    elif name == "consequences":
        logger.info("CONSEQUENCE DISTRIBUTION COMPARISON")
        logger.info("=" * 60)
        results["csq_105"], results["csq_115"] = (
            compare_consequence_distributions_from_sample(global_sample)
        )

    elif name == "consistency":
        logger.info("ANNOTATION CONSISTENCY SPOT CHECK")
        logger.info("=" * 60)
        results["agreement"] = spot_check_shared_annotations_from_sample(global_sample)

    elif name == "loftee":
        logger.info("LOFTEE ANNOTATION COMPARISON")
        logger.info("=" * 60)
        results["loftee_stats"] = compare_loftee_annotations_from_sample(global_sample)

    elif name == "sift_polyphen":
        logger.info("SIFT/POLYPHEN COMPARISON")
        logger.info("=" * 60)
        results["sift_polyphen_stats"] = compare_sift_polyphen_from_sample(
            global_sample
        )

    elif name == "mane":
        logger.info("MANE TRANSCRIPT COVERAGE (VEP 115)")
        logger.info("=" * 60)
        results["mane_stats"] = compare_mane_annotations_from_sample(global_sample)

    elif name == "biotype":
        logger.info("BIOTYPE DISTRIBUTION COMPARISON")
        logger.info("=" * 60)
        results["biotypes_105"], results["biotypes_115"] = (
            compare_biotype_distribution_from_sample(global_sample)
        )

    elif name == "colocated":
        logger.info("COLOCATED VARIANTS COMPARISON")
        logger.info("=" * 60)
        results["colocated_stats"] = check_colocated_variants_from_sample(global_sample)

    elif name == "regulatory":
        logger.info("REGULATORY CONSEQUENCES COMPARISON")
        logger.info("=" * 60)
        results["regulatory_stats"] = compare_regulatory_consequences_from_sample(
            global_sample
        )

    #    elif name == "genes":
    #        logger.info("SPECIFIC GENE SPOT CHECKS")
    #        logger.info("=" * 60)
    #        results["gene_stats"] = check_specific_genes(ht105, ht115)

    elif name == "gene_coverage":
        logger.info("GENE-LEVEL COVERAGE CHECK")
        logger.info("=" * 60)
        interval_ht = ensembl_interval.ht()
        results["interval_105"], results["interval_115"] = run_gene_coverage_check(
            ht105, ht115, interval_ht
        )

    # Debug checks.
    elif name == "debug_tx_counts":
        logger.info("TRANSCRIPT COUNT DISTRIBUTION")
        logger.info("=" * 60)
        results["tx_count_stats"] = compare_transcript_counts_distribution_from_sample(
            global_sample
        )

    elif name == "debug_lof_filter":
        logger.info("LOF FILTER DIFFERENCE DEBUG")
        logger.info("=" * 60)
        debug_lof_filter_differences_from_sample(global_sample)

    elif name == "debug_mane":
        logger.info("MANE ANNOTATION DEBUG")
        logger.info("=" * 60)
        debug_mane_annotations(ht115)

    elif name == "debug_regulatory_motif":
        logger.info("REGULATORY/MOTIF FEATURE DEBUG")
        logger.info("=" * 60)
        debug_regulatory_motif(ht115)

    elif name == "debug_sift_polyphen":
        logger.info("SIFT/POLYPHEN CATEGORY TRANSITIONS")
        logger.info("=" * 60)
        check_sift_polyphen_category_changes_from_sample(global_sample)


def print_summary(results: dict) -> None:
    """Print summary of results."""
    logger.info("=" * 60)
    logger.info("SUMMARY")
    logger.info("=" * 60)

    if "schema_diff" in results:
        logger.info(
            f"Schema fields only in VEP105: {len(results['schema_diff']['only_in_105'])}"
        )
        logger.info(
            f"Schema fields only in VEP115: {len(results['schema_diff']['only_in_115'])}"
        )

    if "count_diffs" in results:
        logger.info(f"Contigs with count differences: {len(results['count_diffs'])}")

    if "miss_105" in results:
        logger.info(
            f"VEP105 vep_missing rate: {results['miss_105'].get('vep_missing', 'N/A')}"
        )
    if "miss_115" in results:
        logger.info(
            f"VEP115 vep_missing rate: {results['miss_115'].get('vep_missing', 'N/A')}"
        )

    if "agreement" in results:
        logger.info(
            f"most_severe_consequence agreement: "
            f"{results['agreement'].get('most_severe_csq_match', 'N/A')}"
        )

    if "loftee_stats" in results:
        logger.info(
            f"LOFTEE category match: {results['loftee_stats'].get('lof_match', 'N/A')}"
        )

    if "sift_polyphen_stats" in results:
        logger.info(
            f"SIFT prediction match: "
            f"{results['sift_polyphen_stats'].get('sift_pred_match', 'N/A')}"
        )
        logger.info(
            f"PolyPhen prediction match: "
            f"{results['sift_polyphen_stats'].get('polyphen_pred_match', 'N/A')}"
        )


def main(args: argparse.Namespace) -> None:
    """Run VEP 105 vs 115 comparison."""
    if args.list_checks:
        print("\nAvailable main checks:")
        for name, desc in AVAILABLE_CHECKS.items():
            print(f"  {name:15} - {desc}")
        print("\nDebug checks (use --debug-only or --include-debug):")
        for name, desc in DEBUG_CHECKS.items():
            print(f"  {name:15} - {desc}")
        return

    hl.init(
        quiet=True,
        log="/vep_comparison.log",
        tmp_dir=args.tmp_dir,
    )
    hl.default_reference("GRCh38")

    logger.info("Loading VEP tables...")
    ht115 = hl.read_table(args.vep115_path)
    ht105 = hl.read_table(args.vep105_path) if args.vep105_path else vep_context.ht()

    # Create global sample once for all checks that need it (unless disabled).
    if args.no_global_sample:
        logger.warning(
            "Global sample creation DISABLED (--no-global-sample). "
            "Checks requiring paired comparisons will be skipped!"
        )
        global_sample = None
    else:
        logger.info("Creating global sample for efficient reuse...")
        global_sample = create_global_sample(ht105, ht115, sample_frac=0.01)

    results = {}

    # Determine which checks to run.
    if args.debug_only:
        checks_to_run = list(DEBUG_CHECKS.keys())
    elif args.checks:
        checks_to_run = args.checks
        if args.include_debug:
            checks_to_run.extend(DEBUG_CHECKS.keys())
    else:
        # Default: all main checks.
        checks_to_run = list(AVAILABLE_CHECKS.keys())
        # Exclude expensive optional checks by default.
        checks_to_run.remove("genes")
        checks_to_run.remove("gene_coverage")
        if args.include_debug:
            checks_to_run.extend(DEBUG_CHECKS.keys())

    logger.info(f"Running checks: {checks_to_run}")

    for check in checks_to_run:
        run_check(check, ht105, ht115, global_sample, results)

    print_summary(results)
    logger.info("Comparison complete!")


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    main(args)
