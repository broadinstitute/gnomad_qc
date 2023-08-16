"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl

from gnomad.assessment.validity_checks import count_vep_annotated_variants_per_interval
from gnomad.resources.grch38.reference_data import ensembl_interval
from gnomad.sample_qc.relatedness import filter_mt_to_trios
from gnomad.utils.annotations import annotate_allele_info
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import (
    default_compute_info,
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
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
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
        - If the allele is a star allele, we set QUALapprox to 0.

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
                    lambda pl_0: hl.if_else((mt.GQ < 2) & (pl_0 == 1), 0, pl_0),
                    hl.bind(lambda x: x[0] - hl.min(x), extract_as_pls(mt.LPL, i[0])),
                ),
            )
            .or_missing()
        )
    )


def run_compute_info(mt: hl.MatrixTable, test: bool = False) -> hl.Table:
    """
    Run compute info on a MatrixTable.

    ..note::

        Adds a fix for AS_QUALapprox by recomputing from LPL because some were found to
        have different lengths than LA.

    :param mt: Input MatrixTable.
    :param test: Whether to use the test dataset. Default is False.
    :return: Table with info annotations.
    """
    mt = mt.annotate_entries(
        gvcf_info=mt.gvcf_info.annotate(
            AS_QUALapprox=recompute_as_qualapprox_from_lpl(mt)
        )
    )

    if test:
        unrelated_expr = ~mt.meta.rand_sampling_meta.related
    else:
        unrelated_expr = ~mt.meta.sample_filters.relatedness_filters.related

    return default_compute_info(
        mt,
        site_annotations=True,
        as_annotations=True,
        ac_filter_groups={"release": mt.meta.release, "unrelated": unrelated_expr},
    )


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
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation(info_ht.info, info_ht.a_index),
        ),
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


def export_tp_vcf(
    transmitted_singletons: bool = True,
    sibling_singletons: bool = True,
):
    """
    Export true positive variants to VCF for use in VQSR.

    :param transmitted_singletons: Should transmitted singletons be included
    :param sibling_singletons: Should sibling singletons be included
    :return: None
    """
    if not (transmitted_singletons | sibling_singletons):
        raise ValueError(
            "At least one of transmitted_singletons or sibling_singletons must be set "
            "to True"
        )

    trio_stats_ht = hl.read_table(
        var_annotations_ht_path("trio_stats", data_source, freeze)
    )
    sib_stats_ht = hl.read_table(
        var_annotations_ht_path("sib_stats", data_source, freeze)
    )

    for transmission_confidence in ["raw", "adj"]:
        filter_expr = False
        true_positive_type = ""
        if transmitted_singletons:
            filter_expr = filter_expr | (
                trio_stats_ht[qc_ac_ht.key][f"n_transmitted_{transmission_confidence}"]
                == 1
            ) & (qc_ac_ht.ac_qc_samples_raw == 2)
            true_positive_type = true_positive_type + "ts_"

        if sibling_singletons:
            filter_expr = filter_expr | (
                sib_stats_ht[qc_ac_ht.key][
                    f"n_sib_shared_variants_{transmission_confidence}"
                ]
                == 1
            ) & (qc_ac_ht.ac_qc_samples_raw == 2)
            true_positive_type = true_positive_type + "ss_"

        ht = qc_ac_ht.filter(filter_expr)
        mt = hl.MatrixTable.from_rows_table(ht)
        logger.info(
            f"Exporting {transmission_confidence} transmitted singleton VCF with"
            f" {mt.count()} variants..."
        )
        hl.export_vcf(
            mt,
            get_true_positive_vcf_path(
                true_positive_type=true_positive_type,
                adj=(transmission_confidence == "adj"),
                data_source=data_source,
                freeze=freeze,
            ),
            tabix=True,
        )


def get_variant_qc_annotation_resources(
    test: bool, overwrite: bool
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC annotation pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC annotation pipeline.
    """
    # Initialize variant QC annotation pipeline resource collection.
    ann_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc_annotation",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the variant QC annotation pipeline.
    compute_info = PipelineStepResourceCollection(
        "--compute-info",
        output_resources={"info_ht": get_info(split=False, test=test)},
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
    run_vep = args.run_vep
    overwrite = args.overwrite
    vqc_resources = get_variant_qc_annotation_resources(test=test, overwrite=overwrite)
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

    if args.compute_info:
        # TODO: is there any reason to also compute info per platform?
        res = vqc_resources.compute_info
        res.check_resource_existence()
        ht = run_compute_info(mt, test=test_dataset)
        ht.write(res.info_ht.path, overwrite=overwrite)

    if args.split_info:
        res = vqc_resources.split_info
        res.check_resource_existence()
        split_info(res.info_ht.ht()).write(res.split_info_ht.path, overwrite=overwrite)

    if args.export_info_vcf:
        res = vqc_resources.export_info_vcf
        res.check_resource_existence()
        hl.export_vcf(adjust_vcf_incompatible_types(res.info_ht.ht()), res.info_vcf)

    if run_vep:
        res = vqc_resources.run_vep
        res.check_resource_existence()
        ht = hl.split_multi(get_gnomad_v4_vds(test=test_dataset).variant_dataset.rows())
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
        ht = run_generate_trio_stats(vds, res.final_ped.pedigree(), res.final_ped.ht())
        ht.write(res.trio_stats_ht.path, overwrite=overwrite)

    if args.generate_sibling_stats:
        res = vqc_resources.generate_sib_stats
        res.check_resource_existence()
        ht = run_generate_sib_stats(mt, res.rel_ht.ht())
        ht.write(res.sib_stats_ht.path, overwrite=overwrite)

    if args.annotate_truth_data:
        logger.info("Joining truth data annotations...")
        ht = get_ukbb_data(data_source, freeze).rows().select()
        truth_ht = get_truth_ht()
        ht = ht.join(truth_ht, how="left")
        ht = ht.checkpoint(
            var_annotations_ht_path("truth_data", data_source, freeze),
            overwrite=overwrite,
        )
        ht.summarize()

    if args.export_true_positive_vcfs:
        logger.info("Exporting true positive variants to VCFs...")
        export_tp_vcf(
            data_source,
            freeze,
            transmitted_singletons=args.transmitted_singletons,
            sibling_singletons=args.sibling_singletons,
        )


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
    parser.add_argument("--compute-info", help="Compute info HT.", action="store_true")
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
    parser.add_argument(
        "--annotate_truth_data",
        help="Creates a HT of UKBB variants annotated with truth sites",
        action="store_true",
    )
    parser.add_argument(
        "--export_true_positive_vcfs",
        help="Exports true positive variants to VCF files.",
        action="store_true",
    )
    parser.add_argument(
        "--transmitted_singletons",
        help=(
            "Include transmitted singletons in the exports of true positive variants to"
            " VCF files."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--sibling_singletons",
        help=(
            "Include sibling singletons in the exports of true positive variants to VCF"
            " files."
        ),
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
