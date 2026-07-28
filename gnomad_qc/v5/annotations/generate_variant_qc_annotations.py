"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging
from typing import Dict, List

import hail as hl
from gnomad.resources.grch38.gnomad import GROUPS
from gnomad.resources.grch38.reference_data import get_truth_ht
from gnomad.utils.annotations import (
    annotate_allele_info,
    get_lowqual_expr,
    pab_max_expr,
)
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
    _get_batch_resource_kwargs,
    _init_hail,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.sample_qc import (
    DENSE_TRIO_TEST_CHROMS,
    dense_trios,
    get_dense_trio_chroms,
    pedigree,
    relatedness,
)

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("generate_variant_qc_annotations")
logger.setLevel(logging.INFO)


def generate_ac_info_ht(vds: hl.vds.VariantDataset) -> hl.Table:
    """
    Compute AC and AC_raw annotations for each allele count filter group.

    Function also adds `AS_pab_max` and `allele_info` annotations.

    :param vds: VariantDataset to use for computing AC and AC_raw annotations.
    :return: Table with AC and AC_raw annotations split by high quality, release, and unrelated.
    """
    mt = vds.variant_data

    ac_filter_groups = {
        "high_quality": mt.meta.high_quality,
        "release": mt.meta.release,
        "unrelated": ~mt.meta.relatedness_inference.relatedness_filters.related,
    }

    mt = mt.annotate_cols(_ac_filter_groups=ac_filter_groups)
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    ac_info_expr = hl.struct()
    # First compute ACs for each non-ref allele, grouped by adj.
    grp_ac_expr = {
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

    ac_info_expr = ac_info_expr.annotate(
        AS_pab_max=pab_max_expr(mt.LGT, mt.LAD, mt.LA, hl.len(mt.alleles))
    )

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


def create_info_ht(
    vcf_path: str,
    header_path: str,
    lowqual_indel_phred_het_prior: int = 40,
    vds: hl.vds.VariantDataset = None,
    test: bool = False,
) -> hl.Table:
    """
    Import a VCF of AoU annotated sites, reformat annotations, and add AS_lowqual.

    :param vcf_path: Path to the annotated sites-only VCF.
    :param header_path: Path to the header file for the VCF.
    :param lowqual_indel_phred_het_prior: Phred-scaled prior for a het genotype at a site with a low quality indel. Default is 40. We use 1/10k bases (phred=40) to be more consistent with the filtering used by Broad's Data Sciences Platform for VQSR.
    :param vds: VariantDataset to use for computing AC and AC_raw annotations.
    :param test: Whether to write run a test using just the first two partitions of the loaded VCF.
    :return: Hail Table with reformatted annotations.
    """
    ht = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        header_file=header_path,
        reference_genome="GRCh38",
    ).rows()

    if test:
        ht = ht._filter_partitions(range(2))

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

    logger.info("Adding AC info annotations to info ht...")
    ac_info_ht = generate_ac_info_ht(vds)
    ac_info_ht = ac_info_ht.checkpoint(hl.utils.new_temp_file("ac_info_ht", "ht"))
    # Annotate ac_info_ht with info annotations because info VCF
    # provided by AoU has variants not present in samples we consider high quality.
    info_ht = ac_info_ht.annotate(**ht[ac_info_ht.key])
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(AS_pab_max=info_ht.AC_info.AS_pab_max),
        AC_info=info_ht.AC_info.drop("AS_pab_max"),
    )
    return info_ht


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


# Per-chromosome trio stats inherit the dense trio MT's fine partitioning
# (~2k-12k partitions each for a few tens of MB), so each is coalesced down to
# this value -- both at generation and again as it's read into the union.
# Across the 22 autosomes that's ~22 * 25 = ~550 partitions combined, versus
# the ~138k we'd get from their native partitionings.
TRIO_STATS_N_PARTITIONS_PER_CHROM = 25

# adj/raw paired count fields produced by generate_trio_stats_expr.
TRIO_STATS_PAIRED_FIELDS = [
    "n_transmitted",
    "n_untransmitted",
    "n_de_novos",
    "ac_parents",
    "an_parents",
    "ac_children",
    "an_children",
]


def validate_trio_stats(
    ht: hl.Table,
    per_chrom_hts: Dict[str, hl.Table],
    expected_chroms: List[str],
) -> None:
    """
    Validate the combined trio stats HT for union integrity and value sanity.

    Raises on hard failures (row/contig loss, schema mismatch, invariant
    violations); logs biological plausibility metrics for eyeballing.

    :param ht: Combined (unioned) trio stats HT.
    :param per_chrom_hts: Per-chromosome trio stats HTs keyed by contig.
    :param expected_chroms: Contigs the union should cover (autosomes).
    """
    # Structural checks (cheap; catch a broken union first).
    ref_ht = per_chrom_hts[expected_chroms[0]]
    if ht.row.dtype != ref_ht.row.dtype or list(ht.key) != list(ref_ht.key):
        raise ValueError("Combined trio stats schema/key differs from per-chrom HT.")

    per_chrom_counts = {c: h.count() for c, h in per_chrom_hts.items()}

    # Single pass over the combined HT for everything else.
    agg = ht.aggregate(
        hl.struct(
            n_rows=hl.agg.count(),
            contigs=hl.agg.counter(ht.locus.contig),
            n_multiallelic=hl.agg.count_where(hl.len(ht.alleles) != 2),
            n_non_autosome=hl.agg.count_where(~ht.locus.in_autosome()),
            n_missing=hl.agg.count_where(
                hl.any(
                    [
                        hl.is_missing(ht[f"{f}_{s}"])
                        for f in TRIO_STATS_PAIRED_FIELDS
                        for s in ["raw", "adj"]
                    ]
                )
            ),
            n_negative=hl.agg.count_where(
                hl.any(
                    [
                        ht[f"{f}_{s}"] < 0
                        for f in TRIO_STATS_PAIRED_FIELDS
                        for s in ["raw", "adj"]
                    ]
                )
            ),
            # adj must be a subset of raw for every paired field.
            n_adj_gt_raw=hl.agg.count_where(
                hl.any(
                    [ht[f"{f}_adj"] > ht[f"{f}_raw"] for f in TRIO_STATS_PAIRED_FIELDS]
                )
            ),
            # AC <= AN for parents and children, both strata.
            n_ac_gt_an=hl.agg.count_where(
                hl.any(
                    [
                        ht[f"ac_{grp}_{s}"] > ht[f"an_{grp}_{s}"]
                        for grp in ["parents", "children"]
                        for s in ["raw", "adj"]
                    ]
                )
            ),
            sum_t_adj=hl.agg.sum(ht.n_transmitted_adj),
            sum_u_adj=hl.agg.sum(ht.n_untransmitted_adj),
            sum_dnm_raw=hl.agg.sum(ht.n_de_novos_raw),
            sum_dnm_adj=hl.agg.sum(ht.n_de_novos_adj),
        )
    )

    observed_contigs = set(agg.contigs)
    problems = []

    # Row-count conservation (union concatenates, so exact).
    total_per_chrom = sum(per_chrom_counts.values())
    if agg.n_rows != total_per_chrom:
        problems.append(
            f"row count {agg.n_rows} != sum of per-chrom counts {total_per_chrom}"
        )
    # Per-contig conservation + contig set.
    if observed_contigs != set(expected_chroms):
        problems.append(
            f"contigs {sorted(observed_contigs)} != expected {sorted(expected_chroms)}"
        )
    for c in expected_chroms:
        if agg.contigs.get(c, 0) != per_chrom_counts.get(c):
            problems.append(
                f"{c}: combined {agg.contigs.get(c, 0)} != per-chrom "
                f"{per_chrom_counts.get(c)}"
            )
    # Content invariants (all must be zero).
    for label, n in [
        ("multiallelic rows", agg.n_multiallelic),
        ("non-autosomal rows", agg.n_non_autosome),
        ("rows with missing counts", agg.n_missing),
        ("rows with negative counts", agg.n_negative),
        ("rows where adj > raw", agg.n_adj_gt_raw),
        ("rows where AC > AN", agg.n_ac_gt_an),
    ]:
        if n != 0:
            problems.append(f"{n} {label}")

    # Biological plausibility (logged, not fatal).
    denom = agg.sum_t_adj + agg.sum_u_adj
    ratio = agg.sum_t_adj / denom if denom > 0 else float("nan")
    logger.info(
        "Combined trio stats: %d rows across %d contigs.",
        agg.n_rows,
        len(observed_contigs),
    )
    logger.info("Transmission ratio (adj) = %.4f (expect ~0.5).", ratio)
    logger.info("De novos: raw=%d, adj=%d.", agg.sum_dnm_raw, agg.sum_dnm_adj)

    if problems:
        raise ValueError(
            "Trio stats validation FAILED:\n  - " + "\n  - ".join(problems)
        )
    logger.info("Trio stats validation PASSED.")


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
    test_chrom = args.test_chrom
    test = args.test or test_n_partitions is not None or test_chrom is not None

    info_ht_path = get_info_ht(test=test, environment=environment).path
    trio_stats_ht_path = get_trio_stats(test=test, environment=environment).path
    sib_stats_ht_path = get_sib_stats(test=test, environment=environment).path
    variant_qc_annotation_ht_path = get_variant_qc_annotations(
        test=test, environment=environment
    ).path

    if args.export_true_positive_vcfs and not (
        args.transmitted_singletons or args.sibling_singletons
    ):
        raise ValueError(
            "--export-true-positive-vcfs requires at least one of"
            " --transmitted-singletons or --sibling-singletons"
        )

    # Load VDS only if needed for create_info_ht or generate_sibling_stats.
    if args.create_info_ht or args.generate_sibling_stats:
        # NOTE: VDS will have 'aou_' prefix on sample IDs.
        vds = get_aou_vds(
            high_quality_only=True,
            filter_partitions=range(test_n_partitions) if test_n_partitions else None,
            annotate_meta=True,
            # NOTE: Using args.test here so that sibling stats test can be calculated from
            # a few partitions of the full (not test) VDS).
            test=args.test,
            environment=environment,
        )

    try:
        if args.create_info_ht:
            _check_resource_existence(
                environment=environment,
                output_step_resources={
                    "info_ht": [info_ht_path],
                },
                overwrite=overwrite,
            )

            ht = create_info_ht(
                vcf_path=get_aou_annotated_sites_only_vcf(environment=environment),
                header_path=get_aou_vcf_header(environment=environment),
                lowqual_indel_phred_het_prior=args.lowqual_indel_phred_het_prior,
                vds=vds,
                test=test,
            )
            ht.write(info_ht_path, overwrite=overwrite)
        if args.export_info_vcf:
            logger.info("Exporting info ht as VCF...")
            out_info_vcf_path = info_vcf_path(test=test, environment=environment)
            # --test-chrom cuts the full info HT to the requested contigs. The test info
            # HT is partition-filtered, so its contigs won't match a test run's -L args.
            in_info_ht_path = (
                get_info_ht(test=False, environment=environment).path
                if test_chrom
                else info_ht_path
            )
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "info_ht": [in_info_ht_path],
                },
                output_step_resources={
                    "info_vcf_path": [out_info_vcf_path],
                },
                overwrite=overwrite,
            )
            info_ht = hl.read_table(in_info_ht_path)
            if test_chrom:
                info_ht = hl.filter_intervals(
                    info_ht,
                    [
                        hl.parse_locus_interval(c, reference_genome="GRCh38")
                        for c in test_chrom
                    ],
                )
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
            # GATK's isolation forest infers allele-specific mode from Number=A in the
            # header (the explicit --use-allele-specific-annotations flag was replaced
            # by this header check in GATK 4.5.0.0); Hail exports arrays as Number=. so
            # declare the AS features Number=A. VQSR (VariantRecalibrator) is flag-driven
            # and works with either.
            vcf_metadata = {
                "info": {
                    f: {"Number": "A", "Type": "Float", "Description": ""}
                    for f in INFO_FEATURES
                }
            }
            hl.export_vcf(info_ht, out_info_vcf_path, tabix=True, metadata=vcf_metadata)

        if args.generate_trio_stats:
            chrom = args.chrom
            if chrom is None:
                raise ValueError(
                    "--chrom is required for --generate-trio-stats; trio stats are "
                    "computed one chromosome at a time (combine with --union-trio-stats)."
                )
            if test and chrom not in DENSE_TRIO_TEST_CHROMS:
                raise ValueError(
                    f"--test computes trio stats only for {DENSE_TRIO_TEST_CHROMS}; got "
                    f"--chrom {chrom}."
                )
            logger.info("Generating trio stats for %s...", chrom)
            per_chrom_trio_stats_path = get_trio_stats(
                test=test, environment=environment, chrom=chrom
            ).path
            _check_resource_existence(
                environment=environment,
                output_step_resources={"trio_stats_ht": [per_chrom_trio_stats_path]},
                overwrite=overwrite,
            )
            ht = run_generate_trio_stats(
                dense_trios(test=test, chrom=chrom, environment=environment).mt(),
                pedigree(test=test, environment=environment).pedigree(),
            )
            # The dense trio MT is finely partitioned, so the per-chrom trio stats
            # inherit ~thousands of tiny partitions. Checkpoint the full result,
            # then naive_coalesce off the materialized partitioning before writing.
            ht = ht.checkpoint(hl.utils.new_temp_file(f"trio_stats_{chrom}", "ht"))
            ht = ht.naive_coalesce(TRIO_STATS_N_PARTITIONS_PER_CHROM)
            ht.write(per_chrom_trio_stats_path, overwrite=overwrite)

        if args.union_trio_stats:
            logger.info("Unioning per-chromosome trio stats...")
            chroms = get_dense_trio_chroms(test)
            per_chrom_trio_stats_paths = [
                get_trio_stats(test=test, environment=environment, chrom=c).path
                for c in chroms
            ]
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "trio_stats_per_chrom": per_chrom_trio_stats_paths
                },
                output_step_resources={"trio_stats_ht": [trio_stats_ht_path]},
                overwrite=overwrite,
            )
            # The per-chrom HTs inherit the dense trio MTs' fine partitioning
            # (~2k-12k partitions each), so naive_coalesce each one as it is read;
            # union() then concatenate to get len(chroms) *
            # TRIO_STATS_N_PARTITIONS_PER_CHROM partitions instead of the sum of every
            # chrom's native count (~138k).
            hts = [
                get_trio_stats(test=test, environment=environment, chrom=c)
                .ht()
                .naive_coalesce(TRIO_STATS_N_PARTITIONS_PER_CHROM)
                for c in chroms
            ]
            ht = hts[0] if len(hts) == 1 else hts[0].union(*hts[1:])
            ht.write(trio_stats_ht_path, overwrite=overwrite)

        if args.validate_trio_stats:
            logger.info("Validating combined trio stats...")
            chroms = get_dense_trio_chroms(test)
            per_chrom = {
                c: get_trio_stats(test=test, environment=environment, chrom=c)
                for c in chroms
            }
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "trio_stats_ht": [trio_stats_ht_path],
                    "trio_stats_per_chrom": [r.path for r in per_chrom.values()],
                },
            )
            validate_trio_stats(
                get_trio_stats(test=test, environment=environment).ht(),
                {c: r.ht() for c, r in per_chrom.items()},
                chroms,
            )

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
        "--generate-trio-stats",
        help=(
            "Calculate trio stats for a single chromosome (--chrom). Run once per "
            "chromosome, then combine with --union-trio-stats."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--union-trio-stats",
        help="Union the per-chromosome trio stats HTs into the combined trio stats HT.",
        action="store_true",
    )
    parser.add_argument(
        "--validate-trio-stats",
        help="Validate the combined trio stats HT (union integrity + value sanity).",
        action="store_true",
    )
    parser.add_argument(
        "--chrom",
        help=(
            "Single chromosome (e.g. 'chr20'). Required for --generate-trio-stats; trio "
            "stats are computed one chromosome at a time."
        ),
        type=str,
        default=None,
    )
    parser.add_argument(
        "--test-chrom",
        help=(
            "Contig(s) (e.g. 'chr22') to cut the full info HT to for --export-info-vcf. "
            "Use to match the contigs a downstream test run scores. Implies --test."
        ),
        type=str,
        nargs="+",
        default=None,
    )
    parser.add_argument(
        "--generate-sibling-stats",
        help="Calculates sibling stats.",
        action="store_true",
    )
    parser.add_argument(
        "--create-info-ht",
        help="Create the info ht containing annotations needed for variant QC.",
        action="store_true",
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
