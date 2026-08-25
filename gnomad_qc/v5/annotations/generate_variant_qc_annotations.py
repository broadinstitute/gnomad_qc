"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import json
import logging
import math
from collections import OrderedDict
from typing import Dict, List

import hail as hl
import hailtop.fs as hfs
from gnomad.resources.grch38.gnomad import GROUPS
from gnomad.resources.grch38.reference_data import get_truth_ht
from gnomad.utils.annotations import (
    annotate_allele_info,
    get_lowqual_expr,
    pab_max_expr,
)
from gnomad.utils.file_utils import file_exists, repartition_for_join
from gnomad.utils.sparse_mt import split_info_annotation
from gnomad.utils.vcf import adjust_vcf_incompatible_types
from gnomad.variant_qc.pipeline import (
    INFO_FEATURES,
    NON_INFO_FEATURES,
    TRUTH_DATA,
    generate_sib_stats,
    generate_trio_stats,
)
from gnomad.variant_qc.random_forest import median_impute_features

from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp, get_adj_expr
from gnomad_qc.v5.resources.annotations import (
    get_ac_info_ht_checkpoint_path,
    get_vcf_ht_checkpoint_path,
    get_aou_annotated_sites_only_vcf,
    get_aou_vcf_header,
    get_info_ht,
    get_sib_stats,
    get_trio_stats,
    get_true_positive_vcf_path,
    get_variant_qc_annotations,
    info_vcf_path,
)
from gnomad_qc.v5.resources.basics import (
    _check_resource_existence,
    qc_temp_prefix,
    _get_batch_resource_kwargs,
    _init_hail,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.sample_qc import dense_trios, pedigree, relatedness

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("generate_variant_qc_annotations")
logger.setLevel(logging.INFO)


def _aou_sample_artifact_path(name: str, environment: str, test: bool = False) -> str:
    """
    Return the path to a chunk-invariant sample-level artifact (30-day storage).

    These are written once by ``--write-sample-artifacts`` so later runs read a
    small file instead of each rescanning the ~330-partition sample tables. This
    pipeline writes and reads two artifacts: ``collisions.json`` (colliding sample
    IDs) and ``high_quality_samples.json`` (post-prefix high-quality sample IDs).
    Readers fall back to computing each one when its file is absent.

    :param name: Artifact file name (with extension).
    :param environment: Compute environment.
    :param test: If True, return the test-scoped path.
    :return: GCS path to the artifact.
    """
    scope = "test" if test else "full"
    return f"{qc_temp_prefix(environment=environment, days=30)}aou_sample_artifacts_{scope}/{name}"


def _read_sample_artifact_json(name: str, environment: str, test: bool = False):
    """
    Return the parsed JSON sample artifact, or None when it has not been precomputed.

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


def _get_sample_artifact_vds_kwargs(environment: str, test: bool) -> dict:
    """
    Build `get_aou_vds` sample-filtering kwargs from precomputed sample artifacts.

    Reads the run-invariant artifacts written by ``--write-sample-artifacts`` and
    falls back to the original in-loader scans for any that are absent:
    ``collisions.json`` replaces the sample-collisions Table scan, and
    ``high_quality_samples.json`` replaces the meta-table filter behind
    ``high_quality_only=True``. The writer uses the exact filter this substitutes
    for (``meta_ht.high_quality``, collected post-prefix), so the two cannot
    drift; ``add_project_prefix=True`` is set explicitly because the JSON holds
    post-prefix IDs.

    :param environment: Compute environment.
    :param test: If True, read the test-scoped artifacts (the 10-sample test VDS).
    :return: Keyword arguments for `get_aou_vds`.
    """
    kwargs = {"high_quality_only": True}
    collisions = _read_sample_artifact_json("collisions.json", environment, test)
    if collisions is not None:
        logger.info(
            "Using precomputed collisions.json (%d sample IDs).", len(collisions)
        )
        kwargs["sample_collisions"] = set(collisions)
    hq_samples = _read_sample_artifact_json(
        "high_quality_samples.json", environment, test
    )
    if hq_samples is not None:
        logger.info(
            "Using precomputed high_quality_samples.json (%d sample IDs).",
            len(hq_samples),
        )
        kwargs.update(
            high_quality_only=False,
            filter_samples=hq_samples,
            add_project_prefix=True,
        )
    return kwargs


def _optional_global(value, dtype) -> hl.expr.Expression:
    """
    Return `value` as a Hail literal of `dtype`, or a typed missing when None.

    Used to build parameter-provenance global structs from optional CLI args,
    whose None values Hail cannot type-infer on its own.

    :param value: Python value or None.
    :param dtype: Hail type for the literal / missing value.
    :return: Hail expression of `dtype`.
    """
    return hl.missing(dtype) if value is None else hl.literal(value, dtype)


def _grp_ac_expr_original(mt: hl.MatrixTable, ac_filter_groups: Dict) -> Dict:
    """
    Compute per-alt AC dicts by scanning every global alt allele (original).

    Per-genotype cost scales with the number of global alts at the locus.

    :param mt: MatrixTable with `_ac_filter_groups` and `alt_alleles_range_array`
        already annotated.
    :param ac_filter_groups: Mapping of filter-group name to boolean col expr.
    :return: Mapping of filter-group name to an array (indexed by alt) of
        `dict<adj, AC>`.
    """
    return {
        f: hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai) & mt._ac_filter_groups[f],
                hl.agg.group_by(
                    get_adj_expr(mt.LGT, mt.GQ, mt.LAD),
                    hl.agg.sum(
                        mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[
                            mt.LA.index(ai)
                        ]
                    ),
                ),
            ),
            mt.alt_alleles_range_array,
        )
        for f in ac_filter_groups
    }


def _grp_ac_expr_local(mt: hl.MatrixTable, ac_filter_groups: Dict) -> Dict:
    """
    Compute per-alt AC dicts in local-allele space (the `use_local_allele_agg` AC path).

    Instead of scanning every global alt allele for every genotype (cost ~
    n_global_alts per genotype), each genotype scatters only its local alt
    alleles as (global_idx, dosage) pairs (cost ~ n_local_alleles per genotype).
    This is much faster at high-allele loci while producing the same
    ``{group: array<dict<adj, AC>>}`` shape as the original aggregation.

    :param mt: MatrixTable with `_ac_filter_groups` and `alt_alleles_range_array`
        already annotated.
    :param ac_filter_groups: Mapping of filter-group name to boolean col expr.
    :return: Mapping of filter-group name to an array (indexed by alt) of
        `dict<adj, AC>`.
    """
    adj = get_adj_expr(mt.LGT, mt.GQ, mt.LAD)
    # Dosage in LOCAL space: local_dosage[j] = copies (0/1/2) of local allele j.
    local_dosage = mt.LGT.one_hot_alleles(hl.len(mt.LA))
    # Per-genotype (global_idx, dosage) pairs for its local alleles, ref removed.
    pairs = (
        hl.range(hl.len(mt.LA))
        .map(lambda j: hl.struct(g=mt.LA[j], d=local_dosage[j]))
        .filter(lambda s: s.g != 0)
    )

    def per_group(f):
        # dict<global_idx, dict<adj, AC>>
        scattered = hl.agg.filter(
            mt._ac_filter_groups[f],
            hl.agg.explode(
                lambda s: hl.agg.group_by(s.g, hl.agg.group_by(adj, hl.agg.sum(s.d))),
                pairs,
            ),
        )
        # Re-expand to the array-indexed-by-alt shape the rest of the pipeline
        # expects: position (ai - 1) holds dict<adj, AC> for global alt ai.
        empty = hl.empty_dict(hl.tbool, hl.tint64)
        return mt.alt_alleles_range_array.map(lambda ai: scattered.get(ai, empty))

    return {f: per_group(f) for f in ac_filter_groups}


def _local_pab_max_expr(mt: hl.MatrixTable) -> hl.expr.ArrayExpression:
    """
    Compute `AS_pab_max` per alt allele in local-allele space.

    `pab_max_expr` expands every genotype's AD to global space (via
    ``local_to_global``, length n_alleles) and array-aggregates over all global
    alts, costing ~n_global_alts per genotype -- the dominant cost at high-allele
    loci that the local-allele AC aggregation alone leaves untouched. This
    computes the identical result by scattering only each genotype's local
    alleles.

    ``pab_max[a] = max over het genotypes of binom_test(AD_global[a], sum(AD), .5)``.
    Because ``binom_test(0, n, .5)`` is the minimum over the success count, every
    genotype's term for an allele it does not carry equals a shared per-genotype
    baseline, and a carrier's real term is always >= its own baseline term. Hence

        ``pab_max[a] == max(scattered_max[a], baseline)``

    where ``baseline`` is the max over het genotypes of ``binom_test(0, sum(AD), .5)``
    and ``scattered_max[a]`` is the max over het carriers of allele ``a``. This
    identity is exact, so the output matches `pab_max_expr` bit-for-bit.

    :param mt: MatrixTable with `alt_alleles_range_array` annotated and LA/LGT/LAD
        entries.
    :return: Array (indexed by alt) of maximum allele-balance binomial p-values.
    """
    is_het = mt.LGT.is_het()
    # Total depth (incl. ref); equals hl.sum of the global AD used by pab_max_expr.
    sum_ad = hl.sum(mt.LAD)
    # Shared per-genotype term for any allele the genotype does not carry.
    baseline = hl.agg.filter(
        is_het, hl.agg.max(hl.binom_test(0, sum_ad, 0.5, "two-sided"))
    )
    # (global_idx, pab) for each local ALT allele the genotype carries.
    pairs = (
        hl.range(hl.len(mt.LA))
        .map(
            lambda j: hl.struct(
                g=mt.LA[j],
                p=hl.binom_test(mt.LAD[j], sum_ad, 0.5, "two-sided"),
            )
        )
        .filter(lambda s: s.g != 0)
    )
    scattered_max = hl.agg.filter(
        is_het,
        hl.agg.explode(lambda s: hl.agg.group_by(s.g, hl.agg.max(s.p)), pairs),
    )
    return mt.alt_alleles_range_array.map(
        lambda ai: hl.if_else(
            scattered_max.contains(ai),
            hl.max(scattered_max[ai], baseline),
            baseline,
        )
    )


def generate_ac_info_ht(
    vds: hl.vds.VariantDataset,
    max_alleles: int = None,
    min_alleles: int = None,
    use_local_allele_agg: bool = False,
) -> hl.Table:
    """
    Compute AC and AC_raw annotations for each allele count filter group.

    Function also adds `AS_pab_max` and `allele_info` annotations.

    With `use_local_allele_agg`, both the AC aggregation and `AS_pab_max` are
    computed in local-allele space: per-genotype cost scales with the number of
    local alleles instead of global alts, which is much faster at high-allele
    loci and produces output identical to the original. The local-allele path
    adds per-genotype hash-map overhead that can be a net loss at low-allele
    loci, so the original dense aggregation remains the default.

    :param vds: VariantDataset to use for computing AC and AC_raw annotations.
    :param max_alleles: Optional maximum number of alleles (including the reference) allowed at a locus. Variants at loci with more than this many alleles are filtered out. Default is None (no filtering).
    :param min_alleles: Optional minimum number of alleles (including the reference) allowed at a locus. Variants at loci with fewer than this many alleles are filtered out. Default is None (no filtering).
    :param use_local_allele_agg: If True, compute the AC aggregation AND `AS_pab_max` in local-allele space. This removes the ~n_global_alts-per-genotype costs in the AC aggregation and `pab_max_expr`, the dominant costs at high-allele loci. Recommended for high-allele strata (e.g. min_alleles >= 10). Default is False (original aggregation).
    :return: Table with AC and AC_raw annotations split by high quality, release, and unrelated.
    """
    mt = vds.variant_data

    if max_alleles is not None:
        logger.info("Filtering out variants with more than %d alleles...", max_alleles)
        mt = mt.filter_rows(hl.len(mt.alleles) <= max_alleles)

    if min_alleles is not None:
        logger.info("Filtering out variants with fewer than %d alleles...", min_alleles)
        mt = mt.filter_rows(hl.len(mt.alleles) >= min_alleles)

    ac_filter_groups = {
        "high_quality": mt.meta.high_quality,
        "release": mt.meta.release,
        "unrelated": ~mt.meta.relatedness_inference.relatedness_filters.related,
    }

    mt = mt.annotate_cols(_ac_filter_groups=ac_filter_groups)
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    ac_info_expr = hl.struct()
    # First compute ACs for each non-ref allele, grouped by adj.
    if use_local_allele_agg:
        logger.info("Using local-allele-space AC aggregation...")
        grp_ac_expr = _grp_ac_expr_local(mt, ac_filter_groups)
    else:
        grp_ac_expr = _grp_ac_expr_original(mt, ac_filter_groups)

    # Then, for each non-ref allele, compute
    # 'AC' as the adj group
    # 'AC_raw' as the sum of adj and non-adj groups
    ac_info_expr = ac_info_expr.annotate(
        **{
            f"AC{'_' + f if f else f}_raw": grp.map(
                lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))
            )
            for f, grp in grp_ac_expr.items()
        },
        **{
            f"AC{'_' + f if f else f}": grp.map(lambda i: hl.int32(i.get(True, 0)))
            for f, grp in grp_ac_expr.items()
        },
    )

    if use_local_allele_agg:
        logger.info("Using local-allele-space AS_pab_max...")
        pab_expr = _local_pab_max_expr(mt)
    else:
        pab_expr = pab_max_expr(mt.LGT, mt.LAD, mt.LA, hl.len(mt.alleles))
    ac_info_expr = ac_info_expr.annotate(AS_pab_max=pab_expr)

    ac_info_ht = mt.select_rows(AC_info=ac_info_expr).rows()

    # Split multi-allelic sites.
    ac_info_ht = annotate_allele_info(ac_info_ht)
    ac_info_ht = ac_info_ht.annotate(
        **{
            a: (
                ac_info_ht[a].annotate(
                    **split_info_annotation(ac_info_ht[a], ac_info_ht.a_index),
                )
            )
            for a in ["AC_info"]
        },
    )

    return ac_info_ht


def validate_local_allele_agg(
    vds: hl.vds.VariantDataset,
    n_partitions: int = 5,
    max_alleles: int = None,
    min_alleles: int = None,
    use_local_allele_agg: bool = False,
    pab_atol: float = 1e-9,
) -> bool:
    """
    Diff the local-allele-space aggregation against the original on a subset.

    Builds the pre-split AC / AC_raw arrays and `AS_pab_max` with both the
    original and the local-allele strategy and asserts they match at every
    locus.

    Comparisons handle missingness and NaN explicitly (rather than letting
    `hl.agg.all` silently skip them), verify the key sets are identical in both
    directions (so row-set drift is caught, not swallowed), and fail on empty
    input (so a zero-locus run is never a vacuous PASS). AC arrays (integers) are
    compared exactly; `AS_pab_max` (float) is compared to within `pab_atol`.

    :param vds: VariantDataset to validate on.
    :param n_partitions: Number of leading `variant_data` partitions to sample.
        Pass None to validate every loaded partition.
    :param max_alleles: Optional max alleles filter applied before aggregation.
    :param min_alleles: Optional min alleles filter applied before aggregation.
    :param use_local_allele_agg: If True, validate the local-allele AC aggregation
        and `AS_pab_max` against the original aggregation.
    :param pab_atol: Absolute tolerance for the `AS_pab_max` float comparison.
        Default is 1e-9.
    :return: True iff the key sets match and every compared field matches at every
        sampled locus (and the input was non-empty).
    """
    if not use_local_allele_agg:
        logger.warning(
            "use_local_allele_agg is False; nothing to validate."
        )
        return True

    def _ac_arrays(candidate: bool) -> hl.Table:
        mt = vds.variant_data
        if n_partitions is not None:
            mt = mt._filter_partitions(range(n_partitions))
        if max_alleles is not None:
            mt = mt.filter_rows(hl.len(mt.alleles) <= max_alleles)
        if min_alleles is not None:
            mt = mt.filter_rows(hl.len(mt.alleles) >= min_alleles)

        ac_filter_groups = {
            "high_quality": mt.meta.high_quality,
            "release": mt.meta.release,
            "unrelated": ~mt.meta.relatedness_inference.relatedness_filters.related,
        }
        mt = mt.annotate_cols(_ac_filter_groups=ac_filter_groups)
        mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

        if candidate:
            grp_ac_expr = _grp_ac_expr_local(mt, ac_filter_groups)
        else:
            grp_ac_expr = _grp_ac_expr_original(mt, ac_filter_groups)
        ac_info_expr = hl.struct().annotate(
            **{
                f"AC{'_' + f if f else f}_raw": grp.map(
                    lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))
                )
                for f, grp in grp_ac_expr.items()
            },
            **{
                f"AC{'_' + f if f else f}": grp.map(lambda i: hl.int32(i.get(True, 0)))
                for f, grp in grp_ac_expr.items()
            },
        )
        ac_info_expr = ac_info_expr.annotate(
            AS_pab_max=(
                _local_pab_max_expr(mt)
                if candidate
                else pab_max_expr(mt.LGT, mt.LAD, mt.LA, hl.len(mt.alleles))
            )
        )
        return mt.select_rows(AC_info=ac_info_expr).rows()

    # Checkpoint both sides so the counts, anti-joins, and field diff below don't
    # each re-run the (expensive) aggregation.
    orig = _ac_arrays(candidate=False).checkpoint(
        hl.utils.new_temp_file("validate_orig", "ht")
    )
    cand = _ac_arrays(candidate=True).checkpoint(
        hl.utils.new_temp_file("validate_cand", "ht")
    )

    # (2) Key-set parity in both directions; (3) empty-input guard.
    n_orig = orig.count()
    n_cand = cand.count()
    if n_orig == 0 and n_cand == 0:
        logger.error(
            "Validation input is empty (0 loci) after partition/allele filters; "
            "nothing was validated. Failing to avoid a vacuous PASS."
        )
        return False
    n_orig_only = orig.anti_join(cand).count()
    n_cand_only = cand.anti_join(orig).count()
    keys_match = n_orig_only == 0 and n_cand_only == 0

    # (4) Compare every field the local-allele path changes: the AC arrays and
    # AS_pab_max.
    compare_fields = list(orig.AC_info)

    # (1) Element-wise array match that treats both-missing / both-NaN as equal
    # (so divergence to NA is caught, not silently skipped) and never yields a
    # missing predicate that `hl.agg.all` would skip. Unequal lengths => mismatch.
    def _elem_int(x, y):
        return hl.if_else(
            hl.is_missing(x) | hl.is_missing(y),
            hl.is_missing(x) & hl.is_missing(y),
            x == y,
        )

    def _elem_float(x, y):
        return (
            hl.case()
            .when(hl.is_missing(x) | hl.is_missing(y), hl.is_missing(x) & hl.is_missing(y))
            .when(hl.is_nan(x) | hl.is_nan(y), hl.is_nan(x) & hl.is_nan(y))
            .default(hl.abs(x - y) <= pab_atol)
        )

    def _arr_match(a, b, is_float):
        elem = _elem_float if is_float else _elem_int
        return hl.if_else(
            hl.is_missing(a) | hl.is_missing(b),
            hl.is_missing(a) & hl.is_missing(b),
            (hl.len(a) == hl.len(b))
            & hl.all(hl.zip(a, b).map(lambda p: elem(p[0], p[1]))),
        )

    joined = orig.annotate(_cand=cand[orig.key].AC_info)
    res = joined.aggregate(
        hl.struct(
            n=hl.agg.count(),
            **{
                fld: hl.agg.all(
                    _arr_match(
                        joined.AC_info[fld],
                        joined._cand[fld],
                        is_float=(fld == "AS_pab_max"),
                    )
                )
                for fld in compare_fields
            },
        )
    )

    label = "local-allele aggregation"
    part_desc = "all loaded" if n_partitions is None else f"{n_partitions}"
    fields_ok = all(bool(res[fld]) for fld in compare_fields)
    all_ok = keys_match and res.n > 0 and fields_ok

    logger.info(
        "Validated %d loci across %s partitions (orig=%d, cand=%d).",
        res.n,
        part_desc,
        n_orig,
        n_cand,
    )
    if not keys_match:
        logger.warning(
            "Key-set mismatch: %d loci only in original, %d only in candidate.",
            n_orig_only,
            n_cand_only,
        )
    for fld in compare_fields:
        logger.info("  %-20s match=%s", fld, bool(res[fld]))
    if all_ok:
        logger.info(
            "PASS: %s matches original on %d loci (%s partitions).",
            label,
            res.n,
            part_desc,
        )
    else:
        mismatched = [f for f in compare_fields if not res[f]]
        logger.warning(
            "FAIL: %s — key_sets_match=%s, mismatched fields=%s",
            label,
            keys_match,
            mismatched,
        )
    return all_ok


def create_sites_vcf_ht(
    vcf_path: str,
    header_path: str,
    lowqual_indel_phred_het_prior: int = 40,
    test: bool = False,
) -> hl.Table:
    """
    Import the AoU annotated sites-only VCF and reformat it into a HT.

    Converts the array-typed AS_* annotations to scalars, drops fields redundant
    with their AS_* versions, and adds `AS_lowqual`.

    :param vcf_path: Path to the annotated sites-only VCF.
    :param header_path: Path to the header file for the VCF.
    :param lowqual_indel_phred_het_prior: Phred-scaled prior for a het genotype at a site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to be more consistent with the filtering used by Broad's Data Sciences Platform for VQSR.
    :param test: Whether to use just the first two partitions of the loaded VCF.
    :return: Reformatted sites VCF Table.
    """
    ht = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        header_file=header_path,
        reference_genome="GRCh38",
    ).rows()

    if test:
        ht = ht._filter_partitions(range(2))

    # NOTE: This VCF is already split to biallelics, so `hl.len(ht.alleles)` is
    # always 2. Allele-count filtering is therefore applied only on the unsplit
    # VDS inside `generate_ac_info_ht`; the final rows come from that side.
    logger.info("Reformatting annotations...")
    array_annotations = [
        "AS_FS",
        "AS_MQ",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_SOR",
    ]

    # AS_VarDP is in the format of "ref|alt" and VarDP is an int of the alt value from this, so can just set AS_VarDP to VarDP.
    # AS_SB_TABLE is an alternate formatting (string ref1, ref2 | alt 1, alt2)
    # of SB_TABLE,  so can just set AS_SB_TABLE to SB_TABLE.
    info_updates = {
        # Convert single-element array annotations to float64.
        **{ann: hl.float64(ht.info[ann][0]) for ann in array_annotations},
        # Extract singular element from AS_QD array and convert to float32.
        "AS_QD": hl.float32(ht.info.AS_QD[0]),
        "AS_QUALapprox": hl.int64(ht.info.QUALapprox),
        "AS_VarDP": ht.info.VarDP,
        "AS_SB_TABLE": ht.info.SB_TABLE,
    }

    ht = ht.transmute(info=ht.info.annotate(**info_updates))
    ht = ht.annotate(
        info=ht.info.drop("SB", "QUALapprox", "VarDP", "SB_TABLE", "AS_RAW_MQ")
    )
    ht = ht.drop("rsid")

    ht = ht.annotate(
        AS_lowqual=get_lowqual_expr(
            ht.alleles,
            ht.info.AS_QUALapprox,
            indel_phred_het_prior=lowqual_indel_phred_het_prior,
        )
    )

    return ht


def join_vcf_and_ac_info_hts(
    ac_info_ht: hl.Table,
    vcf_ht: hl.Table,
    ac_info_ht_path: str = None,
    vcf_ht_path: str = None,
    use_repartition_for_join: bool = False,
) -> hl.Table:
    """
    Annotate the AC info HT with the reformatted sites VCF HT's annotations.

    The join is driven by the AC info side's keys (the AoU VCF has variants not
    present in the samples we consider high quality, and those rows are dropped);
    `AS_pab_max` is moved from `AC_info` into `info` to match the final schema.

    :param ac_info_ht: AC info HT (split, e.g. a `generate_ac_info_ht` checkpoint or a `union_ac_info_hts` output).
    :param vcf_ht: Reformatted sites VCF HT (from `create_sites_vcf_ht`).
    :param ac_info_ht_path: Readable path for `ac_info_ht`. Required when `use_repartition_for_join` is True (both sides are re-read with shared `_intervals`).
    :param vcf_ht_path: Readable path for `vcf_ht`. If None and `use_repartition_for_join` is True, the VCF HT is first written to a temp path.
    :param use_repartition_for_join: If True, co-partition the VCF and AC info HTs onto identical interval boundaries (via `repartition_for_join`) so the join is partition-aligned and requires no shuffle. Best for full-genome runs. Default is False.
    :return: Info HT.
    """
    if use_repartition_for_join:
        if ac_info_ht_path is None:
            raise ValueError(
                "use_repartition_for_join=True requires ac_info_ht_path (a readable "
                "path to re-read the AC info HT with shared _intervals)."
            )
        logger.info("Co-partitioning VCF and AC info HTs for a shuffle-free join...")
        if vcf_ht_path is None:
            vcf_ht_path = hl.utils.new_temp_file("info_vcf_ht", "ht")
            vcf_ht.write(vcf_ht_path, overwrite=True)
        intervals = repartition_for_join(vcf_ht_path)
        vcf_ht = hl.read_table(vcf_ht_path, _intervals=intervals)
        ac_info_ht = hl.read_table(ac_info_ht_path, _intervals=intervals)
    info_ht = ac_info_ht.annotate(**vcf_ht[ac_info_ht.key])
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(AS_pab_max=info_ht.AC_info.AS_pab_max),
        AC_info=info_ht.AC_info.drop("AS_pab_max"),
    )
    return info_ht


def _parse_allele_stratum(spec: str) -> tuple:
    """
    Parse a MIN:MAX allele-stratum spec into (min_alleles, max_alleles).

    Bounds are inclusive and either side may be omitted for an open bound,
    matching the --min-alleles/--max-alleles of the generate runs: ":9" ->
    (None, 9), "10:100" -> (10, 100), "101:" -> (101, None).

    :param spec: Stratum spec string.
    :return: Tuple of (min_alleles, max_alleles), each possibly None.
    """
    if ":" not in spec:
        raise ValueError(
            f"Invalid allele stratum spec {spec!r}; expected MIN:MAX with either "
            "side optional, e.g. ':9', '10:100', '101:'."
        )
    lo, hi = spec.split(":", 1)
    min_alleles = int(lo) if lo else None
    max_alleles = int(hi) if hi else None
    if min_alleles is None and max_alleles is None:
        raise ValueError(f"Allele stratum spec {spec!r} has neither bound.")
    return min_alleles, max_alleles


def union_ac_info_hts(
    ac_info_ht_paths: List[str],
    n_partitions: int = None,
    expected_strata: List[Dict] = None,
) -> hl.Table:
    """
    Union per-allele-count-stratum AC info HTs into a single AC info HT.

    The per-stratum AC info HTs come from separate `--generate-ac-info-ht` runs over
    disjoint allele-count ranges (e.g. `--max-alleles 9`, `--min-alleles 10
    --max-alleles 100`, `--min-alleles 101`), so the result is a row union, not a
    key join. All inputs must share identical row and key schemas. When every input
    carries the `ac_info_ht_parameters` global (written by --generate-ac-info-ht),
    the union replaces it with an `ac_info_ht_strata` array holding each input's
    parameters in input order; otherwise globals are taken from the first input.

    After the union, asserts that the total row count equals the sum of the input row
    counts and that all keys are distinct, so accidental locus overlap between strata
    (e.g. overlapping min/max allele ranges -- both bounds are inclusive) fails loudly
    instead of silently duplicating rows.

    :param ac_info_ht_paths: Paths to the per-stratum AC info HTs to union. Must contain at least two paths.
    :param n_partitions: Optional number of partitions for the unioned Table. Default is None (keep the union's partitioning).
    :param expected_strata: Optional list (aligned with `ac_info_ht_paths`) of dicts of expected `ac_info_ht_parameters` global fields (e.g. contig, min_alleles, max_alleles, test). Each input's recorded parameters must match its expectation -- a path name is only a convention, but the globals are written by the run itself, so this catches a table sitting at the right path with the wrong contents. Inputs without the provenance global fail verification. Default is None (no verification).
    :return: Unioned AC info HT (checkpointed to a temp path).
    """
    if len(ac_info_ht_paths) < 2:
        raise ValueError("union_ac_info_hts requires at least two AC info HT paths.")

    hts = [hl.read_table(p) for p in ac_info_ht_paths]

    # Evaluate each input's provenance global once; used for verification here
    # and for the ac_info_ht_strata global at the end.
    input_params = [
        hl.eval(h.globals.ac_info_ht_parameters)
        if "ac_info_ht_parameters" in list(h.globals)
        else None
        for h in hts
    ]

    if expected_strata is not None:
        if len(expected_strata) != len(ac_info_ht_paths):
            raise ValueError(
                "expected_strata must align 1:1 with ac_info_ht_paths "
                f"({len(expected_strata)} vs {len(ac_info_ht_paths)})."
            )
        mismatches = []
        for path, params, expected in zip(
            ac_info_ht_paths, input_params, expected_strata
        ):
            if params is None:
                mismatches.append(
                    f"{path}: no ac_info_ht_parameters global to verify against"
                )
                continue
            for field, want in expected.items():
                got = params[field]
                if got != want:
                    mismatches.append(f"{path}: {field}={got!r}, expected {want!r}")
        if mismatches:
            raise ValueError(
                "Union input parameter verification failed (each table's recorded "
                "ac_info_ht_parameters must match the stratum its path claims):\n  "
                + "\n  ".join(mismatches)
            )
        logger.info(
            "Verified recorded parameters of all %d union inputs against their "
            "expected strata.",
            len(hts),
        )

    first = hts[0]
    for path, ht in zip(ac_info_ht_paths[1:], hts[1:]):
        if ht.row.dtype != first.row.dtype or ht.key.dtype != first.key.dtype:
            raise ValueError(
                f"Schema mismatch between {ac_info_ht_paths[0]} and {path}:\n"
                f"  {ac_info_ht_paths[0]}: key={first.key.dtype}, row={first.row.dtype}\n"
                f"  {path}: key={ht.key.dtype}, row={ht.row.dtype}"
            )
        if ht.globals.dtype != first.globals.dtype:
            logger.warning(
                "Global schemas differ between %s and %s; globals from the first "
                "input are kept in the union.",
                ac_info_ht_paths[0],
                path,
            )

    part_counts = [ht.count() for ht in hts]
    for path, n in zip(ac_info_ht_paths, part_counts):
        logger.info("Input AC info HT %s: %d rows", path, n)

    ht = first.union(*hts[1:])
    tmp_path = hl.utils.new_temp_file("union_ac_info_ht", "ht")
    ht = ht.checkpoint(tmp_path)

    n_total = ht.count()
    n_expected = sum(part_counts)
    if n_total != n_expected:
        raise ValueError(
            f"Union row count ({n_total}) does not equal the sum of the input row "
            f"counts ({n_expected})."
        )
    n_distinct = ht.distinct().count()
    if n_distinct != n_total:
        raise ValueError(
            f"Union contains {n_total - n_distinct} duplicate keys; the input "
            "strata overlap. Check the --min-alleles/--max-alleles ranges used for "
            "each input run."
        )
    logger.info(
        "Unioned AC info HT contains %d rows across %d inputs.", n_total, len(hts)
    )

    if n_partitions is not None:
        ht = hl.read_table(tmp_path, _n_partitions=n_partitions)

    # Table.union keeps only the first input's globals, which would silently
    # claim one stratum's parameters for the whole table. Replace the single
    # `ac_info_ht_parameters` struct with the full per-stratum list.
    if all(params is not None for params in input_params):
        params_dtype = hts[0].globals.ac_info_ht_parameters.dtype
        ht = ht.drop("ac_info_ht_parameters")
        ht = ht.annotate_globals(
            ac_info_ht_strata=hl.literal(
                input_params, dtype=hl.tarray(params_dtype)
            )
        )
    else:
        logger.warning(
            "Not all union inputs carry the ac_info_ht_parameters global "
            "(pre-provenance checkpoints?); the union keeps the first input's "
            "globals unchanged."
        )

    return ht


def run_generate_trio_stats(
    mt: hl.MatrixTable,
    fam_ped: hl.Pedigree,
) -> hl.Table:
    """
    Generate trio transmission stats from a VariantDataset and pedigree info.

    :param mt: Dense trio MatrixTable.
    :param fam_ped: Pedigree containing trio info.
    :return: Table containing trio stats.
    """
    # Add adj annotation and convert LGT to GT since only
    # autosomal bi-allelics are used to calculate trio stats.
    mt = annotate_adj_no_dp(mt)
    mt = mt.transmute_entries(GT=mt.LGT)
    mt = hl.trio_matrix(mt, pedigree=fam_ped, complete_trios=True)
    return generate_trio_stats(mt)


def run_generate_sib_stats(
    mt: hl.MatrixTable,
    relatedness_ht: hl.Table,
) -> hl.Table:
    """
    Generate sibling stats from a VariantDataset and relatedness info.

    :param mt: Input MatrixTable.
    :param relatedness_ht: Table containing relatedness info.
    :return: Table containing sibling stats.
    """
    mt = annotate_adj_no_dp(mt)
    mt = hl.experimental.sparse_split_multi(mt)
    return generate_sib_stats(mt, relatedness_ht)


def create_variant_qc_annotation_ht(
    info_ht: hl.Table,
    trio_stats_ht: hl.Table,
    sib_stats_ht: hl.Table,
    impute_features: bool = True,
    n_partitions: int = 5000,
) -> hl.Table:
    """
    Create a Table with all necessary annotations for variant QC.

    Annotations that are included:

        Features for RF:
            - variant_type
            - allele_type
            - n_alt_alleles
            - has_star
            - AS_QD
            - AS_pab_max
            - AS_MQRankSum
            - AS_SOR
            - AS_ReadPosRankSum

        Training sites (bool):
            - transmitted_singleton
            - sibling_singleton
            - fail_hard_filters - (ht.AS_QD < 0.5) | (ht.AS_FS > 60) | (ht.AS_MQ < 30)

    :param info_ht: Info Table with split multi-allelics.
    :param trio_stats_ht: Table with trio statistics.
    :param sib_stats_ht: Table with sibling statistics.
    :param impute_features: Whether to impute features using feature medians (this is
        done by variant type).
    :param n_partitions: Number of partitions to use for final annotated Table.
    :return: Hail Table with all annotations needed for variant QC.
    """
    truth_data_ht = get_truth_ht()

    ht = info_ht.transmute(**info_ht.AC_info, **info_ht.allele_info)
    ht = ht.annotate(**ht.info.select(*INFO_FEATURES))

    if impute_features:
        impute_ht = ht.select("variant_type", **ht.info)
        impute_ht = median_impute_features(
            impute_ht, {"variant_type": impute_ht.variant_type}
        ).checkpoint(hl.utils.new_temp_file("median_impute", "ht"))

        impute_result = impute_ht[ht.key]
        ht = ht.annotate(
            info=impute_result.drop("feature_imputed", "variant_type"),
            feature_imputed=impute_result.feature_imputed,
        )
        ht = ht.annotate_globals(feature_medians=hl.eval(impute_ht.feature_medians))

    logger.info("Annotating Table with trio and sibling stats and reference truth data")
    trio_stats_ht = trio_stats_ht.select(
        *[f"{a}_{group}" for a in ["n_transmitted", "ac_children"] for group in GROUPS]
    )
    ht = ht.annotate(
        **trio_stats_ht[ht.key],
        **sib_stats_ht[ht.key],
        **truth_data_ht[ht.key],
    )
    tp_map = {
        "transmitted_singleton": "n_transmitted",
        "sibling_singleton": "n_sib_shared_variants",
    }

    # Filter to only variants found in high quality samples and are not lowqual.
    ht = ht.filter((ht.AC_high_quality_raw > 0) & ~ht.AS_lowqual)

    select_dict = {tp: hl.or_else(ht[tp], False) for tp in TRUTH_DATA}
    select_dict.update(
        {
            f"{tp}_{group}": hl.or_else(
                (ht[f"{n}_{group}"] == 1)
                & (ht[f"AC_high_quality{'' if group == 'adj' else '_raw'}"] == 2),
                False,
            )
            for tp, n in tp_map.items()
            for group in GROUPS
        }
    )
    select_dict.update(
        {
            # NOTE: Previous versions used QD < 2, but we decided this was too
            # stringent based on the distribution of this metric in AoU.
            "fail_hard_filters": (ht.AS_QD < 0.5) | (ht.AS_FS > 60) | (ht.AS_MQ < 30),
            "singleton": ht.AC_release_raw == 1,
            "ac_raw": ht.AC_high_quality_raw,
            "ac": ht.AC_release,
            "ac_unrelated_raw": ht.AC_unrelated_raw,
        }
    )

    if impute_features:
        select_dict["feature_imputed"] = ht.feature_imputed

    ht = ht.select(
        "a_index",
        "was_split",
        *NON_INFO_FEATURES,
        "info",
        **select_dict,
    )

    temp_path = hl.utils.new_temp_file("variant_qc_annotations", "ht")
    ht.write(temp_path)
    ht = hl.read_table(temp_path, _n_partitions=n_partitions)
    ht.describe()

    summary = ht.group_by(
        *TRUTH_DATA, *[f"{tp}_{group}" for tp in tp_map for group in GROUPS]
    ).aggregate(n=hl.agg.count())
    logger.info("Summary of truth data annotations:")
    summary.show(-1)

    return ht


def get_tp_ht_for_vcf_export(
    ht: hl.Table,
    transmitted_singletons: bool = False,
    sibling_singletons: bool = False,
) -> Dict[str, hl.Table]:
    """
    Get Tables with raw and adj true positive variants to export as a VCF for use in VQSR.

    :param ht: Input Table with transmitted singleton and sibling singleton information.
    :param transmitted_singletons: Whether to include transmitted singletons in the
        true positive variants.
    :param sibling_singletons: Whether to include sibling singletons in the true
        positive variants.
    :return: Dictionary of 'raw' and 'adj' true positive variant sites Tables.
    """
    if not transmitted_singletons and not sibling_singletons:
        raise ValueError(
            "At least one of transmitted_singletons or sibling_singletons must be set "
            "to True"
        )
    tp_hts = {}
    for group in GROUPS:
        filter_expr = False
        if transmitted_singletons:
            filter_expr = ht[f"transmitted_singleton_{group}"]
        if sibling_singletons:
            filter_expr = filter_expr | ht[f"sibling_singleton_{group}"]

        filtered_ht = ht.filter(filter_expr).select().select_globals()
        filtered_ht = filtered_ht.checkpoint(
            hl.utils.new_temp_file("true_positive_variants", "ht"),
        )
        logger.info(
            "True positive %s Table for VCF export contains %d variants",
            group,
            filtered_ht.count(),
        )
        tp_hts[group] = filtered_ht

    return tp_hts


def _derive_chunk_locus_intervals(
    vds_filtered: hl.vds.VariantDataset,
    n_subdivisions: int = 1,
    total_subdivisions: int = None,
    reference_genome: str = "GRCh38",
) -> List[hl.utils.Interval]:
    """
    Derive per-contig locus sub-intervals covering the filtered VDS variant_data.

    Bounds are measured on the chunk's VARIANT rows (min/max variant locus per
    contig), not on reference blocks: variant partitions are contiguous locus
    ranges, so adjacent chunks' variant-derived spans provably leave no variant
    in a seam between them, and the measured table is the same one
    `generate_ac_info_ht` later aggregates.

    :param vds_filtered: VDS already filtered to the chunk's partitions.
    :param n_subdivisions: Number of equal-position sub-intervals per contig. Ignored when `total_subdivisions` is set.
    :param total_subdivisions: Optional total number of sub-intervals to spread across all contigs in the chunk, allocated proportionally to each contig's position span (each contig gets at least one). Takes precedence over `n_subdivisions`.
    :param reference_genome: Reference genome for constructed `hl.Locus` objects.
    :return: List of `hl.Interval` objects covering the VDS chunk.
    """
    vd = vds_filtered.variant_data
    bounds = vd.aggregate_rows(
        hl.agg.group_by(
            vd.locus.contig,
            hl.struct(
                lo=hl.agg.min(vd.locus.position),
                hi=hl.agg.max(vd.locus.position),
            ),
        )
    )
    spans = {contig: max(b.hi + 1 - b.lo, 1) for contig, b in bounds.items()}
    total_span = sum(spans.values())
    # Iterate contigs in genomic order: the aggregation dict comes back in
    # lexicographic key order (chr1, chr10, ..., chr2, ...), which read_vds's
    # partitioner rejects for multi-contig chunks.
    contig_order = {
        c: i for i, c in enumerate(hl.get_reference(reference_genome).contigs)
    }
    sub_intervals: List[hl.utils.Interval] = []
    for contig in sorted(bounds, key=lambda c: contig_order[c]):
        b = bounds[contig]
        if total_subdivisions is not None:
            n = max(1, round(total_subdivisions * spans[contig] / total_span))
        else:
            n = max(n_subdivisions, 1)
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


def scout_target_loci(
    vds: hl.vds.VariantDataset,
    checkpoint_path: str,
    min_alleles: int = None,
    max_alleles: int = None,
    contig: str = None,
) -> hl.Table:
    """
    Find loci whose allele count falls within a target range (cheap rows-only pass).

    This is the first pass of a two-pass "scout then re-read" strategy: it scans only
    the row metadata of the VDS variant data (entries are pruned, so no genotype
    aggregation runs) to identify the loci of interest, and writes them to a
    checkpointed Table that stays distributed (no driver collect of every locus).

    :param vds: VariantDataset to scout. Only `variant_data` row keys are read.
    :param checkpoint_path: Path to checkpoint the target loci Table.
    :param min_alleles: Optional minimum number of alleles (including the reference)
        at a locus. Loci with fewer are excluded. Default is None (no minimum).
    :param max_alleles: Optional maximum number of alleles (including the reference)
        at a locus. Loci with more are excluded. Default is None (no maximum).
    :param contig: Optional contig (e.g. "chr1") to restrict the scout to.
        Default is None (whole input).
    :return: Table keyed by `locus, alleles` containing only the target loci.
    """
    ht = vds.variant_data.rows().select()

    if contig is not None:
        ht = hl.filter_intervals(
            ht, [hl.parse_locus_interval(contig, reference_genome="GRCh38")]
        )

    if min_alleles is not None:
        ht = ht.filter(hl.len(ht.alleles) >= min_alleles)
    if max_alleles is not None:
        ht = ht.filter(hl.len(ht.alleles) <= max_alleles)

    return ht.checkpoint(checkpoint_path, overwrite=True)


def group_scout_loci_into_intervals(
    target_ht: hl.Table,
    n_partitions: int = None,
    rows_per_partition: int = None,
    reference_genome: str = "GRCh38",
) -> List[hl.utils.Interval]:
    """
    Group scouted target loci into a bounded set of locus intervals.

    Rather than one interval (and therefore one partition) per target locus, the
    sorted target loci are grouped into chunks so that each resulting interval -- and
    thus each partition produced when the VDS is re-read with these intervals --
    holds many target loci. This keeps the partition count in a healthy range
    regardless of how many target loci there are.

    Intervals never span contigs.

    :param target_ht: Table of target loci (e.g. from `scout_target_loci`).
    :param n_partitions: Desired (approximate) number of intervals/partitions. If
        set, takes precedence over `rows_per_partition`.
    :param rows_per_partition: Target number of loci per interval/partition. Used
        only when `n_partitions` is None. Defaults to 2000 if neither is set.
    :param reference_genome: Reference genome used to order contigs so the returned
        intervals are in genomic order (required by `read_vds`'s partitioner).
    :return: List of `hl.Interval` objects, in genomic order, each covering a
        contiguous chunk of target loci within a single contig.
    """
    loci = target_ht.key_by("locus").select().distinct()
    collected = loci.locus.collect()
    n_loci = len(collected)
    if n_loci == 0:
        return []

    if n_partitions is None:
        rows_per_partition = rows_per_partition or 2000
        n_partitions = max(1, math.ceil(n_loci / rows_per_partition))
    n_partitions = max(1, min(n_partitions, n_loci))
    chunk = max(1, math.ceil(n_loci / n_partitions))

    # Group loci by contig; chunk within contig so no interval spans a contig
    # boundary. `collect()` order is not guaranteed, so contigs and positions are
    # explicitly sorted into genomic order for the read_vds partitioner.
    contig_order = {
        c: i for i, c in enumerate(hl.get_reference(reference_genome).contigs)
    }
    by_contig: "OrderedDict[str, list]" = OrderedDict()
    for loc in collected:
        by_contig.setdefault(loc.contig, []).append(loc)

    intervals: List[hl.utils.Interval] = []
    for contig in sorted(by_contig, key=lambda c: contig_order[c]):
        locs = sorted(by_contig[contig], key=lambda loc: loc.position)
        for i in range(0, len(locs), chunk):
            group = locs[i : i + chunk]
            intervals.append(
                hl.Interval(
                    group[0],
                    group[-1],
                    includes_start=True,
                    includes_end=True,
                )
            )
    return intervals


def compute_scout_intervals(args) -> List[hl.utils.Interval]:
    """
    Run the scout pass and derive partition-scaled re-read intervals.

    :param args: Parsed CLI args.
    :return: List of locus intervals covering the scouted target loci.
    """
    environment = args.environment

    scout_input_partitions = args.test_n_partitions
    vds = get_aou_vds(
        **_get_sample_artifact_vds_kwargs(environment, args.test),
        filter_partitions=(
            range(scout_input_partitions) if scout_input_partitions else None
        ),
        annotate_meta=False,
        test=args.test,
        environment=environment,
    )
    if scout_input_partitions:
        logger.info(
            "Scouting only the first %d input partitions of the VDS",
            scout_input_partitions,
        )

    target_ht = scout_target_loci(
        vds,
        checkpoint_path=hl.utils.new_temp_file("scout_target_loci", "ht"),
        min_alleles=args.min_alleles,
        max_alleles=args.max_alleles,
        contig=args.contig,
    )
    n_targets = target_ht.count()
    logger.info("Scout found %d target loci", n_targets)

    intervals = group_scout_loci_into_intervals(
        target_ht,
        n_partitions=args.scout_n_partitions,
        rows_per_partition=args.scout_rows_per_partition,
    )
    logger.info(
        "Scout derived %d intervals for re-read (avg ~%.1f target loci/interval)",
        len(intervals),
        (n_targets / len(intervals)) if intervals else 0.0,
    )
    return intervals


def compute_chunks(args):
    """
    Derive per-contig locus sub-intervals for a chunk of the VDS.

    :param args: Parsed CLI args.
    :return: List of locus intervals covering the chunk.
    """
    environment = args.environment
    start, stop = args.chunk_start, args.chunk_stop
    partition_range = list(range(start, stop))

    vds_probe = get_aou_vds(
        filter_partitions=partition_range,
        environment=environment,
        remove_hard_filtered_samples=False,
        log_sample_counts=False,
    )

    if args.read_subintervals_scale is not None:
        n_parts = vds_probe.variant_data.n_partitions()
        total_sub = max(1, math.ceil(args.read_subintervals_scale * n_parts))
        logger.info(
            "Scaling explode: %d chunk partitions x %s = %d total sub-intervals",
            n_parts,
            args.read_subintervals_scale,
            total_sub,
        )
        sub_intervals = _derive_chunk_locus_intervals(
            vds_probe, total_subdivisions=total_sub
        )
    else:
        n_sub = max(args.read_subintervals_per_chunk or 1, 1)
        sub_intervals = _derive_chunk_locus_intervals(vds_probe, n_subdivisions=n_sub)

    logger.info(
        "Chunk partitions %d-%d: derived %d sub-intervals", start, stop, len(sub_intervals)
    )

    return sub_intervals


def compute_contig_intervals(args) -> List[hl.utils.Interval]:
    """
    Derive equal-width read intervals covering a single contig.

    The contig's span is known from the reference genome, so unlike chunk mode no
    bounds aggregation over the VDS reference data is needed. With
    --read-subintervals-scale, the sub-interval count is the scale times the
    number of VDS partitions overlapping the contig (a metadata-only lookup);
    otherwise --read-subintervals-per-chunk gives the absolute count.

    :param args: Parsed CLI args.
    :return: List of locus intervals covering the contig, in genomic order.
    """
    contig = args.contig
    reference_genome = "GRCh38"
    rg = hl.get_reference(reference_genome)
    contig_len = rg.contig_length(contig)

    if args.read_subintervals_scale is not None:
        vds_probe = get_aou_vds(
            environment=args.environment,
            remove_hard_filtered_samples=False,
            log_sample_counts=False,
        )
        n_parts = hl.vds.filter_intervals(
            vds_probe,
            [hl.parse_locus_interval(contig, reference_genome=reference_genome)],
        ).variant_data.n_partitions()
        n_sub = max(1, math.ceil(args.read_subintervals_scale * n_parts))
        logger.info(
            "Scaling explode: %d %s partitions x %s = %d total sub-intervals",
            n_parts,
            contig,
            args.read_subintervals_scale,
            n_sub,
        )
    else:
        n_sub = max(args.read_subintervals_per_chunk or 1, 1)
        if n_sub == 1:
            logger.warning(
                "Deriving a single read interval for all of %s (one interval = one "
                "partition on re-read); pass --read-subintervals-per-chunk or "
                "--read-subintervals-scale to subdivide it.",
                contig,
            )

    n_sub = min(n_sub, contig_len)
    step = max(contig_len // n_sub, 1)
    intervals: List[hl.utils.Interval] = []
    for i in range(n_sub):
        sub_lo = 1 + i * step
        last = i == n_sub - 1
        sub_hi = contig_len if last else 1 + (i + 1) * step
        intervals.append(
            hl.Interval(
                hl.Locus(contig, sub_lo, reference_genome=reference_genome),
                hl.Locus(contig, sub_hi, reference_genome=reference_genome),
                includes_start=True,
                includes_end=last,
            )
        )
    return intervals


def _validate_args(args, test: bool) -> None:
    """
    Raise ValueError for invalid CLI flag combinations.

    Only pure flag checks live here; validations coupled to derived values
    (e.g. the union input count) stay with the code that derives them.

    :param args: Parsed CLI args.
    :param test: Whether this is a test run (--test or --test-n-partitions).
    """
    if args.test_n_partitions is not None and (
        args.explode_partitions
        or args.chunk_start is not None
        or args.chunk_stop is not None
    ):
        raise ValueError(
            "--test-n-partitions cannot be combined with chunk processing "
            "(--explode-partitions/--chunk-start/--chunk-stop): the chunk "
            "intervals override the partition filter, so the run would process "
            "the full chunk while silently writing to test output paths. Use one "
            "or the other."
        )

    if (
        args.read_subintervals_scale is not None
        and args.read_subintervals_per_chunk is not None
    ):
        raise ValueError(
            "--read-subintervals-scale and --read-subintervals-per-chunk are "
            "mutually exclusive: use the scale to derive the sub-interval count "
            "from the chunk's partition count, or the absolute per-contig count, "
            "not both."
        )

    if args.contig is not None:
        if args.contig not in hl.get_reference("GRCh38").contigs:
            raise ValueError(
                f"--contig {args.contig!r} is not a GRCh38 contig."
            )
        if (
            args.explode_partitions
            or args.chunk_start is not None
            or args.chunk_stop is not None
        ):
            raise ValueError(
                "--contig and chunk processing (--explode-partitions/"
                "--chunk-start/--chunk-stop) are mutually exclusive: use one "
                "read-restriction mode."
            )
        if args.test_n_partitions is not None:
            raise ValueError(
                "--contig cannot be combined with --test-n-partitions; use a "
                "small contig (e.g. chr21) as the test unit instead."
            )

    if args.scout_alleles and args.min_alleles is None and args.max_alleles is None:
        raise ValueError(
            "--scout-alleles requires at least one of --min-alleles or --max-alleles"
        )

    if args.export_true_positive_vcfs and not (
        args.transmitted_singletons or args.sibling_singletons
    ):
        raise ValueError(
            "--export-true-positive-vcfs requires at least one of"
            " --transmitted-singletons or --sibling-singletons"
        )

    if args.union_ac_info_hts and args.generate_ac_info_ht:
        raise ValueError(
            "--union-ac-info-hts and --generate-ac-info-ht are mutually exclusive "
            "in a single invocation: both write to the AC info HT checkpoint "
            "path. Generate the per-stratum AC info HTs first, then union them "
            "in a separate run."
        )
    if args.use_stable_info_paths and test:
        raise ValueError(
            "--use-stable-info-paths is for production runs; it cannot be "
            "combined with --test/--test-n-partitions (test AC info HTs belong "
            "in the temp bucket)."
        )


def _derive_read_intervals(args, need_intervals: bool):
    """
    Derive the VDS read intervals for the selected read-restriction mode.

    Dispatches to chunk (--explode-partitions), scout (--scout-alleles), or
    contig (--contig) derivation; returns None when no mode is selected or when
    no step in this invocation reads the VDS.

    :param args: Parsed CLI args.
    :param need_intervals: Whether any selected step reads the VDS.
    :return: List of locus intervals, or None.
    """
    sub_intervals = None
    if not need_intervals and (
        args.explode_partitions or args.scout_alleles or args.contig
    ):
        logger.info(
            "Skipping interval derivation: the selected steps read neither the VDS "
            "nor the VCF by interval."
        )
    elif args.explode_partitions:
        logger.info("Explode partitions...")
        sub_intervals = compute_chunks(args)
        logger.info("Derived %d sub-intervals from chunk", len(sub_intervals))
    elif args.scout_alleles:
        logger.info("Scouting target loci by allele count...")
        sub_intervals = compute_scout_intervals(args)
        logger.info("Derived %d scout intervals", len(sub_intervals))
    elif args.contig:
        logger.info("Deriving read intervals for contig %s...", args.contig)
        sub_intervals = compute_contig_intervals(args)
        logger.info("Derived %d contig sub-intervals", len(sub_intervals))

    if sub_intervals is not None and len(sub_intervals) == 0:
        # No intervals is a legitimate outcome (e.g. a scout finding no loci in
        # the allele range, or an empty chunk). Proceed with the empty list --
        # get_aou_vds honors it as a zero-row read -- so the run writes an
        # EMPTY AC info HT with full provenance globals, positively recording
        # that this stratum was attempted. The union's derived-input existence
        # check then passes and the empty member contributes zero rows.
        logger.warning(
            "Interval derivation returned no intervals; writing an empty AC "
            "info HT for this run."
        )
    return sub_intervals


def _resolve_union_inputs(args, test: bool, environment: str):
    """
    Resolve the union step's input paths and expected per-stratum parameters.

    Inputs are either listed explicitly (--union-input-ac-info-ht-paths) or
    derived as the cross-product of --union-contigs and --union-allele-strata
    using the same path derivation as --generate-ac-info-ht. Derived inputs come
    with an expectations list that union_ac_info_hts verifies against each
    table's recorded ac_info_ht_parameters global.

    :param args: Parsed CLI args.
    :param test: Whether this is a test run (--test or --test-n-partitions).
    :param environment: Compute environment.
    :return: Tuple of (input paths or None, expected strata or None).
    """
    union_input_ac_info_ht_paths = args.union_input_ac_info_ht_paths
    union_expected_strata = None
    if args.union_contigs or args.union_allele_strata:
        if args.union_input_ac_info_ht_paths is not None:
            raise ValueError(
                "--union-input-ac-info-ht-paths is mutually exclusive with "
                "--union-contigs/--union-allele-strata: list the inputs "
                "explicitly or derive them, not both."
            )
        if not (args.union_contigs and args.union_allele_strata):
            raise ValueError(
                "--union-contigs and --union-allele-strata must be used together."
            )
        rg_contigs = hl.get_reference("GRCh38").contigs
        contigs = args.union_contigs
        if contigs == ["autosomes"]:
            contigs = [f"chr{i}" for i in range(1, 23)]
        bad_contigs = [c for c in contigs if c not in rg_contigs]
        if bad_contigs:
            raise ValueError(
                f"--union-contigs contains non-GRCh38 contigs: {bad_contigs}"
            )
        strata = [_parse_allele_stratum(s) for s in args.union_allele_strata]
        union_input_ac_info_ht_paths = [
            get_ac_info_ht_checkpoint_path(
                test=test,
                environment=environment,
                stable=args.use_stable_info_paths,
                contig=contig,
                min_alleles=stratum_min,
                max_alleles=stratum_max,
            )
            for contig in contigs
            for stratum_min, stratum_max in strata
        ]
        # Verified against each input's recorded ac_info_ht_parameters global in
        # union_ac_info_hts: the path locates the table, the globals prove it is
        # the stratum the path claims.
        union_expected_strata = [
            {
                "contig": contig,
                "min_alleles": stratum_min,
                "max_alleles": stratum_max,
                "test": test,
            }
            for contig in contigs
            for stratum_min, stratum_max in strata
        ]
        logger.info(
            "Derived %d union input paths (%d contigs x %d allele strata):",
            len(union_input_ac_info_ht_paths),
            len(contigs),
            len(strata),
        )
        for path in union_input_ac_info_ht_paths:
            logger.info("  %s", path)
    if args.union_ac_info_hts and (
        union_input_ac_info_ht_paths is None
        or len(union_input_ac_info_ht_paths) < 2
    ):
        raise ValueError(
            "--union-ac-info-hts requires input paths: pass "
            "--union-input-ac-info-ht-paths (at least two), or derive them with "
            "--union-contigs and --union-allele-strata."
        )
    return union_input_ac_info_ht_paths, union_expected_strata


def main(args):
    """Generate all variant annotations needed for variant QC."""
    environment = args.environment
    _init_hail(
        "generate_variant_qc_annotations",
        environment,
        **_get_batch_resource_kwargs(args),
    )

    overwrite = args.overwrite
    test_n_partitions = args.test_n_partitions
    test = args.test or test_n_partitions is not None

    _validate_args(args, test)

    # Interval derivation is only needed to read the VDS.
    need_vds = args.generate_ac_info_ht or args.generate_sibling_stats
    sub_intervals = _derive_read_intervals(args, need_intervals=need_vds)

    info_ht_path = args.info_ht_path_override or get_info_ht(
        test=test, environment=environment
    ).path
    trio_stats_ht_path = get_trio_stats(test=test, environment=environment).path
    sib_stats_ht_path = get_sib_stats(test=test, environment=environment).path
    variant_qc_annotation_ht_path = get_variant_qc_annotations(
        test=test, environment=environment
    ).path
    # The union output and the final join's input use the same derived "union"
    # path, so the two steps line up without manual path plumbing; per-stratum
    # generate runs get paths derived from their partition/allele parameters.
    ac_info_ht_checkpoint_path = args.ac_info_ht_checkpoint_path_override or (
        get_ac_info_ht_checkpoint_path(
            test=test,
            environment=environment,
            stable=args.use_stable_info_paths,
            test_n_partitions=test_n_partitions,
            contig=args.contig,
            chunk_start=args.chunk_start,
            chunk_stop=args.chunk_stop,
            min_alleles=args.min_alleles,
            max_alleles=args.max_alleles,
            union=args.union_ac_info_hts or args.create_final_info_ht,
        )
    )
    if args.generate_ac_info_ht or args.union_ac_info_hts or args.create_final_info_ht:
        logger.info("AC info HT checkpoint path: %s", ac_info_ht_checkpoint_path)
    vcf_ht_checkpoint_path = args.vcf_ht_checkpoint_path_override or (
        get_vcf_ht_checkpoint_path(
            test=test,
            environment=environment,
            stable=args.use_stable_info_paths,
        )
    )
    if args.create_sites_vcf_ht or args.create_final_info_ht:
        logger.info("Sites VCF HT checkpoint path: %s", vcf_ht_checkpoint_path)

    union_input_ac_info_ht_paths, union_expected_strata = _resolve_union_inputs(
        args, test, environment
    )
    # Load VDS only if needed; the union, sites-VCF, and final-join steps never
    # read it.
    vds = None
    if need_vds:
        # NOTE: VDS will have 'aou_' prefix on sample IDs.
        vds = get_aou_vds(
            **_get_sample_artifact_vds_kwargs(environment, args.test),
            filter_partitions=(
                None
                if sub_intervals is not None
                else range(test_n_partitions) if test_n_partitions else None
            ),
            read_intervals=sub_intervals,
            annotate_meta=True,
            # NOTE: Using args.test here so that sibling stats test can be calculated from
            # a few partitions of the full (not test) VDS).
            test=args.test,
            environment=environment,
        )

    try:
        if args.write_sample_artifacts:
            logger.info("Writing run-invariant sample artifacts...")
            # Import here to avoid circular imports.
            from gnomad_qc.v5.resources.meta import get_sample_id_collisions, meta

            sc_ht = get_sample_id_collisions(environment=environment).ht()
            collisions = sorted(sc_ht.aggregate(hl.agg.collect_as_set(sc_ht.s)))
            path = _aou_sample_artifact_path("collisions.json", environment, args.test)
            with hfs.open(path, "w") as f:
                json.dump(collisions, f)
            logger.info("Wrote %d collision sample IDs: %s", len(collisions), path)

            # Must be the exact filter high_quality_only=True applies in
            # get_aou_vds; this JSON substitutes for that filter.
            meta_ht = meta(data_type="genomes", environment=environment).ht()
            hq_samples = meta_ht.filter(meta_ht.high_quality).s.collect()
            path = _aou_sample_artifact_path(
                "high_quality_samples.json", environment, args.test
            )
            with hfs.open(path, "w") as f:
                json.dump(hq_samples, f)
            logger.info(
                "Wrote %d high-quality sample IDs: %s", len(hq_samples), path
            )

        if args.validate_local_allele_agg and args.generate_ac_info_ht:
            logger.info(
                "Validating local-allele-space aggregation against original..."
            )
            if not validate_local_allele_agg(
                vds,
                n_partitions=None,
                max_alleles=args.max_alleles,
                min_alleles=args.min_alleles,
                use_local_allele_agg=args.use_local_allele_agg,
            ):
                raise ValueError(
                    "Aggregation validation failed; aborting before write."
                )
        if args.generate_ac_info_ht:
            logger.info("Generating AC info HT (no VCF join)...")
            _check_resource_existence(
                environment=environment,
                output_step_resources={
                    "ac_info_ht": [ac_info_ht_checkpoint_path],
                },
                overwrite=overwrite,
            )
            ht = generate_ac_info_ht(
                vds,
                max_alleles=args.max_alleles,
                min_alleles=args.min_alleles,
                use_local_allele_agg=args.use_local_allele_agg,
            )
            # Record the parameters that determined this stratum's contents so
            # the checkpoint is self-describing; the union step gathers these
            # into an `ac_info_ht_strata` array on the combined table.
            ht = ht.annotate_globals(
                ac_info_ht_parameters=hl.struct(
                    environment=environment,
                    test=test,
                    test_n_partitions=_optional_global(test_n_partitions, hl.tint32),
                    contig=_optional_global(args.contig, hl.tstr),
                    explode_partitions=args.explode_partitions,
                    chunk_start=_optional_global(args.chunk_start, hl.tint32),
                    chunk_stop=_optional_global(args.chunk_stop, hl.tint32),
                    read_subintervals_per_chunk=_optional_global(
                        args.read_subintervals_per_chunk, hl.tint32
                    ),
                    read_subintervals_scale=_optional_global(
                        args.read_subintervals_scale, hl.tfloat64
                    ),
                    scout_alleles=args.scout_alleles,
                    scout_n_partitions=_optional_global(
                        args.scout_n_partitions, hl.tint32
                    ),
                    scout_rows_per_partition=_optional_global(
                        args.scout_rows_per_partition, hl.tint32
                    ),
                    min_alleles=_optional_global(args.min_alleles, hl.tint32),
                    max_alleles=_optional_global(args.max_alleles, hl.tint32),
                    use_local_allele_agg=args.use_local_allele_agg,
                )
            )
            ht.write(ac_info_ht_checkpoint_path, overwrite=overwrite)

        if args.union_ac_info_hts:
            logger.info("Unioning per-stratum AC info HTs...")
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "union_input_ac_info_hts": union_input_ac_info_ht_paths,
                },
                output_step_resources={
                    "ac_info_ht": [ac_info_ht_checkpoint_path],
                },
                overwrite=overwrite,
            )
            ht = union_ac_info_hts(
                union_input_ac_info_ht_paths,
                n_partitions=args.union_n_partitions,
                expected_strata=union_expected_strata,
            )
            ht.write(ac_info_ht_checkpoint_path, overwrite=overwrite)

        if args.create_sites_vcf_ht:
            logger.info("Importing and reformatting the sites-only VCF into a HT...")
            _check_resource_existence(
                environment=environment,
                output_step_resources={
                    "vcf_ht": [vcf_ht_checkpoint_path],
                },
                overwrite=overwrite,
            )
            ht = create_sites_vcf_ht(
                vcf_path=get_aou_annotated_sites_only_vcf(environment=environment),
                header_path=get_aou_vcf_header(environment=environment),
                lowqual_indel_phred_het_prior=args.lowqual_indel_phred_het_prior,
                test=test,
            )
            ht.write(vcf_ht_checkpoint_path, overwrite=overwrite)

        if args.create_final_info_ht:
            logger.info(
                "Joining the AC info HT with the reformatted VCF HT to create the "
                "final info HT..."
            )
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "ac_info_ht": [ac_info_ht_checkpoint_path],
                    "vcf_ht": [vcf_ht_checkpoint_path],
                },
                output_step_resources={"info_ht": [info_ht_path]},
                overwrite=overwrite,
            )
            ac_info_ht = hl.read_table(ac_info_ht_checkpoint_path)
            logger.info(
                "Reusing reformatted VCF HT at %s...", vcf_ht_checkpoint_path
            )
            vcf_ht = hl.read_table(vcf_ht_checkpoint_path)
            ht = join_vcf_and_ac_info_hts(
                ac_info_ht,
                vcf_ht,
                ac_info_ht_path=ac_info_ht_checkpoint_path,
                vcf_ht_path=vcf_ht_checkpoint_path,
                use_repartition_for_join=args.repartition_for_join,
            )
            ht.write(info_ht_path, overwrite=overwrite)

        if args.export_info_vcf:
            logger.info("Exporting info ht as VCF...")
            out_info_vcf_path = info_vcf_path(test=test, environment=environment)
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "info_ht": [info_ht_path],
                },
                output_step_resources={
                    "info_vcf_path": [out_info_vcf_path],
                },
                overwrite=overwrite,
            )
            info_ht = hl.read_table(info_ht_path)
            # TODO: Check if AS_QUALapprox and AS_VarDP are needed for v5 (not used for v4) and if so need preceeded pipe.
            # Reformat AS_SB_TABLE to be a nested array of arrays for proper use
            # within the 'adjust_vcf_incompatible_types' function.
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    AS_SB_TABLE=hl.array(
                        [info_ht.info.AS_SB_TABLE[0:2], info_ht.info.AS_SB_TABLE[2:4]]
                    )
                )
            )
            info_ht = adjust_vcf_incompatible_types(
                info_ht, pipe_delimited_annotations=[]
            )
            hl.export_vcf(info_ht, out_info_vcf_path, tabix=True)

        if args.generate_trio_stats:
            logger.info("Generating trio stats...")
            _check_resource_existence(
                environment=environment,
                output_step_resources={"trio_stats_ht": [trio_stats_ht_path]},
                overwrite=overwrite,
            )

            ht = run_generate_trio_stats(
                dense_trios(test=test, environment=environment).mt(),
                pedigree(test=test, environment=environment).pedigree(),
            )
            ht.write(trio_stats_ht_path, overwrite=overwrite)

        if args.generate_sibling_stats:
            logger.info("Generating sibling stats...")
            _check_resource_existence(
                environment=environment,
                output_step_resources={"sib_stats_ht": [sib_stats_ht_path]},
                overwrite=overwrite,
            )
            # Note: Checked sibling IDs; none of them have sample ID collisions.
            ht = run_generate_sib_stats(
                vds.variant_data, relatedness(environment=environment).ht()
            )
            ht.write(sib_stats_ht_path, overwrite=overwrite)

        if args.create_variant_qc_annotation_ht:
            logger.info("Creating variant QC annotation HT...")
            _check_resource_existence(
                environment=environment,
                output_step_resources={
                    "variant_qc_annotation_ht": [variant_qc_annotation_ht_path]
                },
                overwrite=overwrite,
            )
            ht = create_variant_qc_annotation_ht(
                hl.read_table(info_ht_path),
                hl.read_table(trio_stats_ht_path),
                hl.read_table(sib_stats_ht_path),
                impute_features=args.impute_features,
                n_partitions=args.n_partitions,
            )
            ht.write(variant_qc_annotation_ht_path, overwrite=overwrite)

        if args.export_true_positive_vcfs:
            logger.info("Exporting true positive VCFs...")

            tp_parts = []
            if args.transmitted_singletons:
                tp_parts.append("transmitted_singleton")
            if args.sibling_singletons:
                tp_parts.append("sibling_singleton")
            tp_type = ".".join(tp_parts)

            vcf_path_kwargs = dict(
                test=test, environment=environment, true_positive_type=tp_type
            )
            raw_tp_vcf_path = get_true_positive_vcf_path(**vcf_path_kwargs, adj=False)
            adj_tp_vcf_path = get_true_positive_vcf_path(**vcf_path_kwargs, adj=True)

            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "variant_qc_annotation_ht": [variant_qc_annotation_ht_path],
                },
                output_step_resources={
                    "raw_true_positive_vcf_path": [raw_tp_vcf_path],
                    "adj_true_positive_vcf_path": [adj_tp_vcf_path],
                },
                overwrite=overwrite,
            )

            tp_hts = get_tp_ht_for_vcf_export(
                hl.read_table(variant_qc_annotation_ht_path),
                transmitted_singletons=args.transmitted_singletons,
                sibling_singletons=args.sibling_singletons,
            )
            hl.export_vcf(tp_hts["raw"], raw_tp_vcf_path, tabix=True)
            hl.export_vcf(tp_hts["adj"], adj_tp_vcf_path, tabix=True)

    finally:
        if environment == "rwb":
            logger.info("Copying log to logging bucket...")
            hl.copy_log(
                get_logging_path(
                    "generate_variant_qc_annotations", environment=environment
                )
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--environment",
        help="Environment where script will run.",
        choices=["rwb", "batch"],
        default="rwb",
    )
    parser.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    parser.add_argument(
        "--driver-cores",
        help="Number of cores. Applies to Batch environment only. Hail default is 1 if unspecified.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--driver-memory",
        help="Memory for driver node. Applies to Batch environment only. Hail default is 'standard' if unspecified.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--worker-cores",
        help="Number of cores. Applies to Batch environment only. Hail default is 1 if unspecified.",
        type=int,
        default=None,
    )
    parser.add_argument(
        "--worker-memory",
        help="Memory for worker nodes. Applies to Batch environment only. Hail default is 'standard' if unspecified.",
        type=str,
        default=None,
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Write to test path.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-partitions",
        help="Use only n partitions of the VDS as input for testing purposes (default: 2).",
        nargs="?",
        const=2,
        type=int,
    )
    parser.add_argument(
        "--write-sample-artifacts",
        help=(
            "Precompute run-invariant sample artifacts (collisions.json and "
            "high_quality_samples.json) to the shared 30-day artifact path, so "
            "every subsequent VDS load reads two small JSONs instead of "
            "rescanning the ~330-partition sample tables. Run once before a "
            "fan-out of --generate-ac-info-ht jobs; loads fall back to the "
            "original scans when the artifacts are absent. The path scheme "
            "matches the pending --write-sample-artifacts on the frequency "
            "pipeline branch, so either writer's collisions.json serves both."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--generate-trio-stats",
        help="Calculates trio stats.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-sibling-stats",
        help="Calculates sibling stats.",
        action="store_true",
    )
    parser.add_argument(
        "--vcf-ht-checkpoint-path-override",
        help=(
            "Optional override path for the reformatted sites VCF HT. By "
            "default the path is derived (sites_vcf[_test].ht in the 30-day temp "
            "bucket, or the durable annotations bucket with "
            "--use-stable-info-paths); written by --create-sites-vcf-ht and read "
            "by --create-final-info-ht."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--repartition-for-join",
        help=(
            "With --create-final-info-ht, co-partition the VCF and AC info HTs "
            "onto identical interval boundaries before the join so it requires no "
            "shuffle (removes the join shuffle stage). Best for full-genome runs."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--info-ht-path-override",
        help=(
            "Optional override path for the info HT output. "
            "If set, this path is used instead of the default resource path."
        ),
        type=str,
        default=None,
    )
    split_workflow_args = parser.add_argument_group(
        "Split info HT workflow",
        "Create the info HT in independent steps: per-stratum AC info HTs "
        "(--generate-ac-info-ht, one run per --min-alleles/--max-alleles range), "
        "union them (--union-ac-info-hts), reformat the sites-only VCF once "
        "(--create-sites-vcf-ht), and join to write the final info HT "
        "(--create-final-info-ht). The VCF is never interval-filtered, so split "
        "(min-repped) VCF loci can never be dropped at interval boundaries.",
    )
    split_workflow_args.add_argument(
        "--generate-ac-info-ht",
        help=(
            "Run only the AC info aggregation on the VDS (no VCF import or join). "
            "Combine with --min-alleles/--max-alleles, --use-local-allele-agg, and "
            "scout/chunk processing. The checkpoint path is derived automatically "
            "from the run's partition/allele parameters (see "
            "--ac-info-ht-checkpoint-path-override to override)."
        ),
        action="store_true",
    )
    split_workflow_args.add_argument(
        "--union-ac-info-hts",
        help=(
            "Union per-allele-count-stratum AC info HTs (from separate "
            "--generate-ac-info-ht runs) and write the result to "
            "--ac-info-ht-checkpoint-path-override. Fails if the input strata overlap "
            "(both allele bounds are inclusive, so e.g. max 9 pairs with min 10). "
            "Mutually exclusive with --generate-ac-info-ht."
        ),
        action="store_true",
    )
    split_workflow_args.add_argument(
        "--union-input-ac-info-ht-paths",
        help=(
            "Space-separated paths to the per-stratum AC info HTs to union with "
            "--union-ac-info-hts. At least two paths are required. Alternatively "
            "derive the list with --union-contigs/--union-allele-strata."
        ),
        nargs="+",
        type=str,
        default=None,
    )
    split_workflow_args.add_argument(
        "--union-contigs",
        help=(
            "With --union-ac-info-hts, derive the input paths instead of listing "
            "them: one per contig x allele stratum, using the same path "
            "derivation as --generate-ac-info-ht (so run flags like --test must "
            "match the generate runs). Pass contig names (e.g. chr1 chr2 chr3) "
            "or 'autosomes' for chr1-chr22. Requires "
            "--union-allele-strata; mutually exclusive with "
            "--union-input-ac-info-ht-paths. Derived inputs are existence-checked "
            "up front, so a forgotten stratum run fails before any compute."
        ),
        nargs="+",
        type=str,
        default=None,
    )
    split_workflow_args.add_argument(
        "--union-allele-strata",
        help=(
            "Allele strata as MIN:MAX specs (inclusive, either side optional) "
            "matching the --min-alleles/--max-alleles of the generate runs, e.g. "
            "':9' '10:100' '101:'. Used with --union-contigs."
        ),
        nargs="+",
        type=str,
        default=None,
    )
    split_workflow_args.add_argument(
        "--union-n-partitions",
        help=(
            "Optional number of partitions for the unioned AC info HT. Default is "
            "None (keep the union's partitioning)."
        ),
        type=int,
        default=None,
    )
    split_workflow_args.add_argument(
        "--create-sites-vcf-ht",
        help=(
            "Import and reformat the full AoU annotated sites-only VCF into a HT "
            "(adds AS_lowqual) and write it to the derived sites VCF HT path "
            "(see --vcf-ht-checkpoint-path-override to override). No VDS is loaded and no "
            "interval filtering is applied."
        ),
        action="store_true",
    )
    split_workflow_args.add_argument(
        "--create-final-info-ht",
        help=(
            "Join the AC info HT at --ac-info-ht-checkpoint-path-override (e.g. the "
            "--union-ac-info-hts output) with the reformatted VCF HT to write the "
            "final info HT to --info-ht-path-override (or the default resource path). "
            "Reads the reformatted VCF HT at the derived path (or "
            "--vcf-ht-checkpoint-path-override override), so run --create-sites-vcf-ht "
            "first. Combine with --repartition-for-join for a shuffle-free "
            "full-genome join."
        ),
        action="store_true",
    )
    split_workflow_args.add_argument(
        "--use-stable-info-paths",
        help=(
            "Derive AC info HT and sites VCF HT checkpoint paths under the "
            "durable annotations bucket instead of the 30-day temp bucket, so "
            "production artifacts are kept. Every step of a workflow "
            "(--generate-ac-info-ht, --union-ac-info-hts, --create-sites-vcf-ht, "
            "--create-final-info-ht) must agree on this flag for the derived "
            "paths to line up. Not allowed with --test/--test-n-partitions."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--ac-info-ht-checkpoint-path-override",
        help=(
            "Optional override path for the AC info HT checkpoint. By default the "
            "path is derived automatically from the run's parameters "
            "(test/--test-n-partitions, --chunk-start/--chunk-stop, "
            "--min-alleles/--max-alleles, union), landing in the 30-day temp "
            "bucket, or in the durable annotations bucket with "
            "--use-stable-info-paths; the --union-ac-info-hts output and "
            "the --create-final-info-ht input share the derived 'union' path."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--contig",
        help=(
            "Restrict the --generate-ac-info-ht VDS read to a single contig "
            "(e.g. chr1). Without --scout-alleles, the contig span is subdivided "
            "into read intervals using --read-subintervals-scale (times the "
            "contig's VDS partition count) or --read-subintervals-per-chunk "
            "(absolute count); no bounds aggregation is needed since the contig "
            "span is known. With --scout-alleles, the scout pass is restricted "
            "to the contig. Mutually exclusive with chunk processing and "
            "--test-n-partitions."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--max-alleles",
        help=(
            "Optional maximum number of alleles (including the reference) allowed at "
            "a locus when creating the info HT. Variants at loci with more than this "
            "many alleles are filtered out. Default is None (no filtering)."
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--min-alleles",
        help=(
            "Optional minimum number of alleles (including the reference) allowed at "
            "a locus when creating the info HT. Variants at loci with fewer than this "
            "many alleles are filtered out. Default is None (no filtering)."
        ),
        type=int,
        default=None,
    )
    parser.add_argument(
        "--lowqual-indel-phred-het-prior",
        help="Phred-scaled prior for a het genotype at a site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to be more consistent with the filtering used by Broad's Data Sciences Platform for VQSR.",
        default=40,
        type=int,
    )
    parser.add_argument(
        "--export-info-vcf", help="Export info ht as VCF.", action="store_true"
    )
    parser.add_argument(
        "--use-local-allele-agg",
        help=(
            "Compute the AC aggregation and AS_pab_max in local-allele space, so "
            "per-genotype cost scales with each genotype's local alleles instead "
            "of every global alt. Output is identical to the original "
            "aggregation. Recommended for high-allele strata (e.g. --min-alleles "
            "10 and above); at low-allele loci the original default is slightly "
            "cheaper. (Formerly --run-opt2.)"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--validate-local-allele-agg",
        help=(
            "Before writing, diff the --use-local-allele-agg aggregation (AC "
            "arrays and AS_pab_max) against the original on the loaded VDS and "
            "abort if they do not match. Use with a small input (e.g. "
            "--scout-alleles or --test-n-partitions)."
        ),
        action="store_true",
    )

    chunk_args = parser.add_argument_group("Chunk processing parameters")
    chunk_args.add_argument(
        "--explode-partitions",
        help="Derive locus sub-intervals for chunk processing.",
        action="store_true",
    )
    chunk_args.add_argument(
        "--chunk-start",
        help="Start partition index for chunk processing.",
        type=int,
        default=None,
    )
    chunk_args.add_argument(
        "--chunk-stop",
        help="Stop partition index for chunk processing (exclusive).",
        type=int,
        default=None,
    )
    chunk_args.add_argument(
        "--read-subintervals-per-chunk",
        help=(
            "Number of locus sub-intervals to subdivide each chunk into, per "
            "contig. Default is None (1 sub-interval per contig, legacy "
            "behavior). Mutually exclusive with --read-subintervals-scale."
        ),
        type=int,
        default=None,
    )
    chunk_args.add_argument(
        "--read-subintervals-scale",
        help=(
            "Derive the total sub-interval count automatically as this multiplier "
            "times the number of VDS partitions in the chunk (e.g. 15 with a "
            "3-partition chunk gives ~45 sub-intervals), allocated across contigs "
            "proportionally to their position span. Use instead of "
            "--read-subintervals-per-chunk when the chunk's partition count is not "
            "known up front (e.g. a --contig run). Mutually exclusive with "
            "--read-subintervals-per-chunk."
        ),
        type=float,
        default=None,
    )

    scout_args = parser.add_argument_group(
        "Scout target loci by allele count",
        "Two-pass mode: a cheap rows-only pass finds loci within the "
        "--min-alleles/--max-alleles range, then the VDS is re-read targeting only "
        "those loci with a scaled (not per-row) number of partitions.",
    )
    scout_args.add_argument(
        "--scout-alleles",
        help=(
            "Enable the two-pass scout mode when creating the info HT. Requires at "
            "least one of --min-alleles or --max-alleles."
        ),
        action="store_true",
    )
    scout_args.add_argument(
        "--scout-n-partitions",
        help=(
            "Approximate number of intervals/partitions to spread the scouted target "
            "loci across on re-read. Takes precedence over "
            "--scout-rows-per-partition."
        ),
        type=int,
        default=None,
    )
    scout_args.add_argument(
        "--scout-rows-per-partition",
        help=(
            "Target number of scouted loci per interval/partition on re-read. Used "
            "only when --scout-n-partitions is not set. Default is 2000."
        ),
        type=int,
        default=None,
    )

    variant_qc_annotation_args = parser.add_argument_group(
        "Variant QC annotation HT parameters."
    )
    variant_qc_annotation_args.add_argument(
        "--create-variant-qc-annotation-ht",
        help="Creates an annotated HT with features for variant QC.",
        action="store_true",
    )
    variant_qc_annotation_args.add_argument(
        "--impute-features",
        help="If set, imputation is performed for variant QC features.",
        action="store_true",
    )
    variant_qc_annotation_args.add_argument(
        "--n-partitions",
        help="Desired number of partitions for variant QC annotation HT.",
        type=int,
        default=5000,
    )
    tp_vcf_args = parser.add_argument_group(
        "Export true positive VCFs",
        "Arguments used to define true positive variant set.",
    )
    tp_vcf_args.add_argument(
        "--export-true-positive-vcfs",
        help=(
            "Exports true positive variants (--transmitted-singletons and/or"
            " --sibling-singletons) to VCF files."
        ),
        action="store_true",
    )
    tp_vcf_args.add_argument(
        "--transmitted-singletons",
        help=(
            "Include transmitted singletons in the exports of true positive variants to"
            " VCF files."
        ),
        action="store_true",
    )
    tp_vcf_args.add_argument(
        "--sibling-singletons",
        help=(
            "Include sibling singletons in the exports of true positive variants to VCF"
            " files."
        ),
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
