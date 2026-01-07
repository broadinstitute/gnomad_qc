"""Script to generate annotations for variant QC on gnomAD v5."""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import GROUPS
from gnomad.resources.grch38.reference_data import get_truth_ht
from gnomad.utils.annotations import (
    annotate_allele_info,
    get_lowqual_expr,
    pab_max_expr,
)
from gnomad.utils.sparse_mt import split_info_annotation
from gnomad.variant_qc.pipeline import generate_sib_stats, generate_trio_stats
from gnomad.variant_qc.random_forest import median_impute_features

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.annotations.generate_variant_qc_annotations import (
    INFO_FEATURES,
    NON_INFO_FEATURES,
    TRUTH_DATA,
)
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp, get_adj_expr
from gnomad_qc.v5.resources.annotations import (
    aou_annotated_sites_only_vcf,
    aou_vcf_header,
    get_info_ht,
    get_sib_stats,
    get_trio_stats,
    get_variant_qc_annotations,
)
from gnomad_qc.v5.resources.basics import get_aou_vds, get_logging_path
from gnomad_qc.v5.resources.constants import GNOMAD_TMP_BUCKET, WORKSPACE_BUCKET
from gnomad_qc.v5.resources.sample_qc import dense_trios, pedigree, relatedness

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
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
    ht = ht.annotate(**ac_info_ht[ht.key])
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
            - fail_hard_filters - (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30)

    :param info_ht: Info Table with split multi-allelics.
    :param trio_stats_ht: Table with trio statistics.
    :param sib_stats_ht: Table with sibling statistics.
    :param impute_features: Whether to impute features using feature medians (this is
        done by variant type).
    :param n_partitions: Number of partitions to use for final annotated Table.
    :return: Hail Table with all annotations needed for variant QC.
    """
    truth_data_ht = get_truth_ht()

    ht = info_ht.transmute(
        **info_ht.AC_info, **info_ht.allele_info, **info_ht.site_info
    )
    ht = ht.annotate(info=ht.info.select(*INFO_FEATURES))

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
                & (ht[f"AC_high_quality{'' if group == 'adj' else f'_raw'}"] == 2),
                False,
            )
            for tp, n in tp_map.items()
            for group in GROUPS
        }
    )
    select_dict.update(
        {
            "fail_hard_filters": (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30),
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


def main(args):
    """Generate all variant annotations needed for variant QC."""
    if args.rwb:
        environment = "rwb"
        hl.init(
            log="/home/jupyter/workspaces/gnomadproduction/generate_variant_qc_annotations.log",
            tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
        )
    else:
        environment = "batch"
        hl.init(
            tmp_dir=f"gs://{GNOMAD_TMP_BUCKET}-4day",
            log="generate_variant_qc_annotations.log",
        )
        # TODO: Add machine configurations for Batch.
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test_n_partitions = args.test_n_partitions
    test = args.test or test_n_partitions is not None

    info_ht_path = get_info_ht(test=test, environment=environment).path
    trio_stats_ht_path = get_trio_stats(test=test, environment=environment).path
    sib_stats_ht_path = get_sib_stats(test=test, environment=environment).path
    variant_qc_annotation_ht_path = get_variant_qc_annotations(
        test=test, environment=environment
    ).path

    # NOTE: VDS will have 'aou_' prefix on sample IDs.
    vds = get_aou_vds(
        high_quality_only=True,
        filter_partitions=range(test_n_partitions) if test_n_partitions else None,
        annotate_meta=True,
        # NOTE: Using args.test here so that sibling stats test can be calculated from
        # a few partitions of the full (not test) VDS).
        test=args.test,
    )

    try:
        if args.create_info_ht:
            check_resource_existence(
                output_step_resources={
                    "info_ht": [info_ht_path],
                },
                overwrite=overwrite,
            )

            ht = create_info_ht(
                vcf_path=aou_annotated_sites_only_vcf,
                header_path=aou_vcf_header,
                lowqual_indel_phred_het_prior=args.lowqual_indel_phred_het_prior,
                vds=vds,
                test=test,
            )
            ht.write(info_ht_path, overwrite=overwrite)

        if args.generate_trio_stats:
            logger.info("Generating trio stats...")
            check_resource_existence(
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
            check_resource_existence(
                output_step_resources={"sib_stats_ht": [sib_stats_ht_path]},
                overwrite=overwrite,
            )
            # Note: Checked sibling IDs; none of them have sample ID collisions.
            ht = run_generate_sib_stats(vds.variant_data, relatedness().ht())
            ht.write(sib_stats_ht_path, overwrite=overwrite)

        if args.create_variant_qc_annotation_ht:
            logger.info("Creating variant QC annotation HT...")
            check_resource_existence(
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

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(
            get_logging_path("generate_variant_qc_annotations", environment=environment)
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rwb",
        help="Run the script in RWB environment.",
        action="store_true",
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
        help="Calculates trio stats.",
        action="store_true",
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

    variant_qc_annotation_args = parser.add_argument_group(
        "Variant QC annotation HT parameters"
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
        help="Desired number of partitions for variant QC annotation HT .",
        type=int,
        default=5000,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
