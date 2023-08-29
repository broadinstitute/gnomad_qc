"""Script to generate annotations for variant QC on gnomAD v4."""
import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.assessment.validity_checks import count_vep_annotated_variants_per_interval
from gnomad.resources.grch38.reference_data import ensembl_interval
from gnomad.resources.resource_utils import TableResource
from gnomad.sample_qc.relatedness import filter_mt_to_trios
from gnomad.utils.annotations import annotate_allele_info
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import (
    default_compute_info,
    get_as_info_expr,
    split_info_annotation,
    split_lowqual_annotation,
)
from gnomad.utils.vcf import adjust_vcf_incompatible_types
from gnomad.utils.vep import vep_or_lookup_vep
from gnomad.variant_qc.pipeline import generate_sib_stats, generate_trio_stats

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import (
    get_info,
    get_sib_stats,
    get_trio_stats,
    get_vep,
    info_vcf_path,
    validate_vep_path,
)
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.sample_qc import pedigree, relatedness

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def extract_as_pls(
    lpl_expr: hl.expr.ArrayExpression,
    allele_idx: hl.expr.Int32Expression,
) -> hl.expr.ArrayExpression:
    """
    Extract PLs for a specific allele from an LPL array expression.

    PL/LPL represents the normalized Phred-scaled likelihoods of the possible
    genotypes from all considered alleles (or local alleles).

    If three alleles are considered, LPL genotype indexes are:
    [0/0, 0/1, 1/1, 0/2, 1/2, 2/2].

    If we want to extract the PLs for each alternate allele, we need to extract:
        - allele 1: [0/0, 0/1, 1/1]
        - allele 2: [0/0, 0/2, 2/2]

    Example:
        - LPL: [138, 98, 154, 26, 0, 14]
        - Extract allele 1 PLs: [0/0, 0/1, 1/1] -> [138, 98, 154]
        - Extract allele 2 PLs: [0/0, 0/2, 2/2] -> [138, 26, 14]

    :param lpl_expr: LPL ArrayExpression.
    :param allele_idx: The index of the alternate allele to extract PLs for.
    :return: ArrayExpression of PLs for the specified allele.
    """
    calls_to_keep = hl.array(
        [hl.call(0, 0), hl.call(0, allele_idx), hl.call(allele_idx, allele_idx)]
    )
    return calls_to_keep.map(lambda c: lpl_expr[c.unphased_diploid_gt_index()])


def recompute_as_qualapprox_from_lpl(mt: hl.MatrixTable) -> hl.expr.ArrayExpression:
    """
    Recompute AS_QUALapprox from LPL.

    QUALapprox is the (Phred-scaled) probability that all reads at the site are hom-ref,
    so QUALapprox is PL[0]. To get the QUALapprox for just one allele, pull out the
    PLs for just that allele, then normalize by subtracting the smallest element from
    all the entries (so the best genotype is 0) and then use the normalized PL[0]
    value for that allele's QUALapprox.

    .. note::

        - The first element of AS_QUALapprox is always None.
        - If the allele is a star allele, we set QUALapprox for that allele to 0.
        - If GQ == 0 and PL[0] for the allele == 1, we set QUALapprox for the allele
          to 0.

    Example:
        Starting Values:
            - alleles: [‘G’, ‘*’, ‘A’, ‘C’, ‘GCTT’, ‘GT’, ‘T’]
            - LGT: 1/2
            - LA: [0, 1, 6]
            - LPL: [138, 98, 154, 26, 0, 14]
            - QUALapprox: 138

        Use `extract_as_pls` to get PLs for each allele:
            - allele 1: [138, 98, 154]
            - allele 2: [138, 26, 14]

        Normalize PLs by subtracting the smallest element from all the PLs:
            - allele 1: [138-98, 98-98, 154-98] -> [40, 0, 56]
            - allele 2: [138-14, 26-14, 14-14] -> [124, 12, 0]

        Use the first element of the allele specific PLs to generate AS_QUALapprox:
        [None, 40, 124]

        Set QUALapprox to 0 for the star allele: [None, 0, 124]

    :param mt: Input MatrixTable.
    :return: AS_QUALapprox ArrayExpression recomputed from LPL.
    """
    return hl.enumerate(mt.LA).map(
        lambda i: (
            hl.case()
            .when(mt.alleles[i[1]] == "*", 0)
            .when(
                i[0] > 0,
                hl.bind(
                    lambda pl_0: hl.if_else((mt.GQ == 0) & (pl_0 == 1), 0, pl_0),
                    hl.bind(lambda x: x[0] - hl.min(x), extract_as_pls(mt.LPL, i[0])),
                ),
            )
            .or_missing()
        )
    )


def correct_as_annotations(
    mt: hl.MatrixTable,
    set_to_missing: bool = False,
) -> hl.expr.StructExpression:
    """
    Correct allele specific annotations that are longer than the length of LA.

    For some entries in the MatrixTable, the following annotations are longer than LA,
    when they should be the same length as LA:
        - AS_SB_TABLE
        - AS_RAW_MQ
        - AS_RAW_ReadPosRankSum
        - AS_RAW_MQRankSum

    This function corrects these annotations by either dropping the alternate allele
    with the index corresponding to the min value of AS_RAW_MQ, or setting them to
    missing if `set_to_missing` is True.

    :param mt: Input MatrixTable.
    :param set_to_missing: Whether to set the annotations to missing instead of
        correcting them.
    :return: StructExpression with corrected allele specific annotations.
    """
    annotations_to_correct = [
        "AS_SB_TABLE",
        "AS_RAW_MQ",
        "AS_RAW_ReadPosRankSum",
        "AS_RAW_MQRankSum",
    ]
    annotations_to_correct = {a: mt.gvcf_info[a] for a in annotations_to_correct}

    # Identify index corresponding to min of AS_RAW_MQ, skipping the reference allele.
    as_raw_mq_no_ref = mt.gvcf_info.AS_RAW_MQ[1:]
    idx_remove = as_raw_mq_no_ref.index(hl.min(as_raw_mq_no_ref)) + 1

    corrected_annotations = {
        a: hl.if_else(
            hl.len(expr) > hl.len(mt.LA),
            hl.or_missing(
                not set_to_missing, expr[:idx_remove].extend(expr[idx_remove + 1 :])
            ),
            expr,
        )
        for a, expr in annotations_to_correct.items()
    }

    return hl.struct(**corrected_annotations)


def run_compute_info(
    mt: hl.MatrixTable,
    max_n_alleles: Optional[int] = None,
    min_n_alleles: Optional[int] = None,
) -> hl.Table:
    """
    Run compute info on a MatrixTable.

    ..note::

        Adds a fix for AS_QUALapprox by recomputing from LPL because some were found to
        have different lengths than LA.

    :param mt: Input MatrixTable.
    :param max_n_alleles: Maximum number of alleles for the site to be included in
        computations.
    :param min_n_alleles: Minimum number of alleles for the site to be included in
        computations.
    :return: Table with info annotations.
    """
    if max_n_alleles:
        mt = mt.filter_rows(hl.len(mt.alleles) <= max_n_alleles)
    if min_n_alleles:
        mt = mt.filter_rows(hl.len(mt.alleles) >= min_n_alleles)

    mt = mt.annotate_entries(
        gvcf_info=mt.gvcf_info.annotate(
            AS_QUALapprox=recompute_as_qualapprox_from_lpl(mt),
            **correct_as_annotations(mt),
        )
    )
    info_ht = default_compute_info(
        mt,
        site_annotations=True,
        as_annotations=True,
        ac_filter_groups={
            "high_quality": mt.meta.high_quality,
            "release": mt.meta.release,
            "unrelated": ~mt.meta.sample_filters.relatedness_filters.related,
        },
        n_partitions=None,
    )

    mt = mt.annotate_entries(gvcf_info=correct_as_annotations(mt, set_to_missing=True))
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))
    ht = mt.select_rows(
        **get_as_info_expr(
            mt,
            sum_agg_fields=["AS_RAW_MQ"],
            int32_sum_agg_fields=[],
            median_agg_fields=["AS_RAW_ReadPosRankSum", "AS_RAW_MQRankSum"],
            array_sum_agg_fields=["AS_SB_TABLE"],
            treat_fields_as_allele_specific=True,
        )
    ).rows()

    info_ht = info_ht.annotate(set_long_AS_missing_info=ht[info_ht.key])

    return info_ht


def split_info(info_ht: hl.Table) -> hl.Table:
    """
    Generate an info Table with split multi-allelic sites from the multi-allelic info Table.

    .. note::

        gnomad_methods' `annotate_allele_info` splits multi-allelic sites before the
        `info` annotation is split to ensure that all sites in the returned Table are
        annotated with allele info.

    :param info_ht: Info Table with unsplit multi-allelics.
    :return: Info Table with split multi-allelics.
    """
    info_ht = annotate_allele_info(info_ht)
    info_annotations_to_split = ["info", "quasi_info", "set_long_AS_missing_info"]
    info_ht = info_ht.annotate(
        **{
            a: info_ht[a].annotate(
                **split_info_annotation(info_ht[a], info_ht.a_index),
            )
            for a in info_annotations_to_split
        },
        AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )

    return info_ht


def run_generate_trio_stats(
    vds: hl.vds.VariantDataset,
    fam_ped: hl.Pedigree,
    fam_ht: hl.Table,
) -> hl.Table:
    """
    Generate trio transmission stats from a VariantDataset and pedigree info.

    :param vds: VariantDataset to generate trio stats from.
    :param fam_ped: Pedigree containing trio info.
    :param fam_ht: Table containing trio info.
    :return: Table containing trio stats.
    """
    # Filter the VDS to autosomes.
    vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)
    vmt = vds.variant_data
    rmt = vds.reference_data

    # Filter the variant data to bi-allelic sites.
    vmt = vmt.filter_rows(hl.len(vmt.alleles) == 2)

    # Filter the variant data and reference data to only the trios.
    vmt = filter_mt_to_trios(vmt, fam_ht)
    rmt = rmt.filter_cols(hl.is_defined(vmt.cols()[rmt.col_key]))

    mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(rmt, vmt))
    mt = mt.transmute_entries(GT=mt.LGT)
    mt = hl.trio_matrix(mt, pedigree=fam_ped, complete_trios=True)

    return generate_trio_stats(mt, bi_allelic_only=False)


def run_generate_sib_stats(
    mt: hl.MatrixTable,
    rel_ht: hl.Table,
) -> hl.Table:
    """
    Generate stats for the number of alternate alleles in common between sibling pairs.

    :param mt: MatrixTable to generate sibling stats from.
    :param rel_ht: Table containing relatedness info for pairs in `mt`.
    :return: Table containing sibling stats.
    """
    # Filter relatedness Table to only exomes-exome pairs.
    rel_ht = rel_ht.filter(
        (rel_ht.i.data_type == "exomes") & (rel_ht.j.data_type == "exomes")
    )
    return generate_sib_stats(mt.transmute_entries(GT=mt.LGT), rel_ht)


def get_variant_qc_annotation_resources(
    test: bool,
    overwrite: bool,
    over_n_alleles: Optional[bool] = None,
    combine_compute_info: bool = False,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC annotation pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :param over_n_alleles: Whether to use a temporary info TableResource for results.
        When True, use temporary info TableResource for only sites that have more
        than the passed arg --compute-info-split-n-alleles alleles. When False, use
        temporary info TableResource for only sites with fewer alleles. When None,
        the finalize info ht is used instead of a temporary location. Default is None.
    :param combine_compute_info: Whether the input for --compute-info should be the two
        temporary files (with and without the --compute-info-over-split-n-alleles flag)
        produced by running --compute-info with --compute-info-split-n-alleles.
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC annotation pipeline.
    """
    if over_n_alleles is None or combine_compute_info:
        info_ht = get_info(split=False, test=test)
    else:
        info_ht = TableResource(
            get_checkpoint_path(
                f"compute_info{'.test' if test else ''}.{'over_n_alleles' if over_n_alleles else 'under_n_alleles'}"
            )
        )
    compute_info_input_resources = {}
    if combine_compute_info:
        res_key = "--compute-info --compute-info-split-n-alleles"
        for n_alleles in ["under", "over"]:
            if n_alleles == "over":
                res_key += " --compute-info-over-split-n-alleles"
            compute_info_input_resources[res_key] = {
                f"{n_alleles}_info_ht": TableResource(
                    get_checkpoint_path(
                        f"compute_info{'.test' if test else ''}.{n_alleles}_n_alleles"
                    )
                )
            }

    # Initialize variant QC annotation pipeline resource collection.
    ann_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc_annotation",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the variant QC annotation pipeline.
    compute_info = PipelineStepResourceCollection(
        "--compute-info",
        input_resources=compute_info_input_resources,
        output_resources={"info_ht": info_ht},
    )
    split_info_ann = PipelineStepResourceCollection(
        "--split-info",
        output_resources={"split_info_ht": get_info(test=test)},
        pipeline_input_steps=[compute_info],
    )
    export_info_vcf = PipelineStepResourceCollection(
        "--export-info-vcf",
        output_resources={"info_vcf": info_vcf_path(test=test)},
        pipeline_input_steps=[compute_info],
    )
    run_vep = PipelineStepResourceCollection(
        "--run-vep",
        output_resources={"vep_ht": get_vep(data_type="exomes", test=test)},
    )
    validate_vep = PipelineStepResourceCollection(
        "--validate-vep",
        output_resources={"vep_count_ht": validate_vep_path(test=test)},
        pipeline_input_steps=[run_vep],
    )
    trio_stats = PipelineStepResourceCollection(
        "--generate-trio-stats",
        output_resources={"trio_stats_ht": get_trio_stats(test=test)},
        input_resources={"identify_trios.py --finalize-ped": {"final_ped": pedigree()}},
    )
    sib_stats = PipelineStepResourceCollection(
        "--generate-sib-stats",
        output_resources={"sib_stats_ht": get_sib_stats(test=test)},
        input_resources={
            "relatedness.py --finalize-relatedness-ht": {"rel_ht": relatedness()}
        },
    )

    # Add all steps to the variant QC annotation pipeline resource collection.
    ann_pipeline.add_steps(
        {
            "compute_info": compute_info,
            "split_info": split_info_ann,
            "export_info_vcf": export_info_vcf,
            "run_vep": run_vep,
            "validate_vep": validate_vep,
            "generate_trio_stats": trio_stats,
            "generate_sib_stats": sib_stats,
        }
    )

    return ann_pipeline


def main(args):
    """Generate all variant annotations needed for variant QC."""
    hl.init(
        default_reference="GRCh38",
        log="/variant_qc_annotations.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test_dataset = args.test_dataset
    test_n_partitions = args.test_n_partitions
    test = test_dataset or test_n_partitions
    over_split_n_alleles = args.compute_info_over_split_n_alleles
    split_n_alleles = args.compute_info_split_n_alleles
    combine_compute_info = args.combine_compute_info
    run_vep = args.run_vep
    overwrite = args.overwrite

    max_n_alleles = min_n_alleles = over_n_alleles = None
    if split_n_alleles is not None:
        if over_split_n_alleles:
            min_n_alleles = split_n_alleles
            over_n_alleles = True
        else:
            max_n_alleles = split_n_alleles - 1
            over_n_alleles = False

    vqc_resources = get_variant_qc_annotation_resources(
        test=test,
        overwrite=overwrite,
        over_n_alleles=over_n_alleles,
        combine_compute_info=combine_compute_info,
    )
    vds = get_gnomad_v4_vds(
        test=test_dataset,
        high_quality_only=True,
        # Keep control/truth samples because they are used in variant QC.
        keep_controls=True,
        annotate_meta=True,
    )
    mt = vds.variant_data

    if test_n_partitions:
        mt = mt._filter_partitions(range(test_n_partitions))

    try:
        if args.compute_info:
            # TODO: is there any reason to also compute info per platform?
            res = vqc_resources.compute_info
            res.check_resource_existence()
            if combine_compute_info:
                ht = mt.select_rows().rows()
                lt_ht = res.under_info_ht.ht()
                gt_eq_ht = res.over_info_ht.ht()
                ht = ht.annotate(**hl.coalesce(lt_ht[ht.key], gt_eq_ht[ht.key]))
            else:
                ht = run_compute_info(mt, max_n_alleles, min_n_alleles)

            if split_n_alleles is None or combine_compute_info:
                ht = ht.naive_coalesce(args.compute_info_n_partitions)

            ht.write(res.info_ht.path, overwrite=overwrite)

        if args.split_info:
            res = vqc_resources.split_info
            res.check_resource_existence()
            split_info(res.info_ht.ht()).write(
                res.split_info_ht.path, overwrite=overwrite
            )

        if args.export_info_vcf:
            res = vqc_resources.export_info_vcf
            res.check_resource_existence()
            hl.export_vcf(adjust_vcf_incompatible_types(res.info_ht.ht()), res.info_vcf)

        if run_vep:
            res = vqc_resources.run_vep
            res.check_resource_existence()
            ht = hl.split_multi(
                get_gnomad_v4_vds(test=test_dataset).variant_dataset.rows()
            )
            ht = vep_or_lookup_vep(ht, vep_version=args.vep_version)
            ht.write(res.vep_ht.path, overwrite=overwrite)

        if args.validate_vep:
            res = vqc_resources.validate_vep
            res.check_resource_existence()
            count_ht = count_vep_annotated_variants_per_interval(
                res.vep_ht.ht(), ensembl_interval.ht()
            )
            count_ht.write(res.vep_count_ht.path, overwrite=args.overwrite)

        if args.generate_trio_stats:
            res = vqc_resources.generate_trio_stats
            res.check_resource_existence()
            ht = run_generate_trio_stats(
                vds, res.final_ped.pedigree(), res.final_ped.ht()
            )
            ht.write(res.trio_stats_ht.path, overwrite=overwrite)

        if args.generate_sibling_stats:
            res = vqc_resources.generate_sib_stats
            res.check_resource_existence()
            ht = run_generate_sib_stats(mt, res.rel_ht.ht())
            ht.write(res.sib_stats_ht.path, overwrite=overwrite)

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("variant_qc_annotations.log"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--test-dataset", help="Use the test dataset as input.", action="store_true"
    )
    parser.add_argument(
        "--test-n-partitions",
        help="Use only 2 partitions of the VDS as input for testing purposes.",
        nargs="?",
        const=2,
        type=int,
    )
    compute_info_args = parser.add_argument_group(
        "Compute info HT.",
        "Arguments relevant to computing the info HT..",
    )
    compute_info_args.add_argument(
        "--compute-info", help="Compute info HT.", action="store_true"
    )
    compute_info_args.add_argument(
        "--compute-info-split-n-alleles",
        help=(
            "Number of alleles at a site to filter results on. If "
            "--compute-info-over-split-n-alleles is used, the results are filtered to "
            "sites with greater than or equal to the value supplied, otherwise the "
            "results are filtered to sites with less than the value supplied. By "
            "default no sites are filtered."
        ),
        default=None,
        type=int,
    )
    compute_info_args.add_argument(
        "--compute-info-over-split-n-alleles",
        help=(
            "Whether to filter to sites greater than or equal to the value supplied to"
            "--compute-info-split-n-alleles. By default, sites are filtered to sites "
            "less than that value, or None if --compute-info-split-n-alleles is not "
            "supplied."
        ),
        action="store_true",
    )
    compute_info_args.add_argument(
        "--combine-compute-info",
        help=(
            "Whether to combine the output from running --compute-info "
            "--compute-info-split-n-alleles with and without the "
            "--compute-info-over-split-n-alleles flag."
        ),
        action="store_true",
    )
    compute_info_args.add_argument(
        "--compute-info-n-partitions",
        help="Number of desired partitions for the info HT.",
        default=5000,
        type=int,
    )

    parser.add_argument("--split-info", help="Split info HT.", action="store_true")
    parser.add_argument(
        "--export-info-vcf", help="Export info as VCF.", action="store_true"
    )
    parser.add_argument(
        "--run-vep", help="Generates vep annotations.", action="store_true"
    )
    parser.add_argument(
        "--validate-vep",
        help=(
            "Validate that variants in protein-coding genes are correctly annotated by"
            " VEP."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--vep-version",
        help="Version of VEPed context Table to use in vep_or_lookup_vep.",
        default="105",
    )
    parser.add_argument(
        "--generate-trio-stats", help="Calculates trio stats", action="store_true"
    )
    parser.add_argument(
        "--generate-sibling-stats",
        help="Calculated sibling variant sharing stats",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
