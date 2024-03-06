"""Script to validate and export only the false duplication Hail Table of 3 genes to a VCF."""
import argparse
import logging

import hail as hl
from gnomad.utils.annotations import (
    faf_expr,
    gen_anc_faf_max_expr,
    merge_freq_arrays,
    pop_max_expr,
)
from gnomad.utils.release import make_freq_index_dict_from_meta

from gnomad_qc.v2.annotations.generate_frequency_data import POPS_TO_REMOVE_FOR_POPMAX
from gnomad_qc.v2.resources.basics import get_gnomad_liftover_data_path
from gnomad_qc.v4.resources.release import get_false_dup_genes_path

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("false_dup_genes")
logger.setLevel(logging.INFO)

FALSE_DUP_GENES = ["KCNE1", "CBS", "CRYAA"]

from typing import Dict, List, Optional

from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    AS_VQSR_FIELDS,
    FAF_POPS,
    FORMAT_DICT,
    HISTS,
    IN_SILICO_ANNOTATIONS_INFO_DICT,
    INFO_DICT,
    REGION_FLAG_FIELDS,
    SEXES,
    SITE_FIELDS,
    VRS_FIELDS_DICT,
    add_as_info_dict,
    adjust_vcf_incompatible_types,
    build_vcf_export_reference,
    create_label_groups,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
    make_vcf_filter_dict,
    rekey_new_reference,
)
from gnomad.utils.vep import (
    # CURRENT_VEP_VERSION,
    VEP_CSQ_FIELDS,
    VEP_CSQ_HEADER,
    vep_struct_to_csq,
)

from gnomad_qc.v4.create_release.validate_and_export_vcf import process_vep_csq_header
from gnomad_qc.v4.resources.release import get_false_dup_genes_path

FALSE_DUP_INFO = [
    "pab_max",
    "info_MQRankSum",
    "info_SOR",
    "info_InbreedingCoeff",
    "info_ReadPosRankSum",
    "info_FS",
    "info_QD",
    "info_MQ",
    "info_DP",
    "qual",
]

FALSE_DUP_AS = [
    "DB",
    "DP",
    "DS",
    "FS",
    "InbreedingCoeff",
    "MQ",
    "MQRankSum",
    "QD",
    "ReadPosRankSum",
    "SOR",
]

FALSE_DUP_AS_VQSR = ["culprit", "VQSLOD", "NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE"]

FALSE_DUP_AGE_HISTS = [
    "age_hist_het_controls_exomes",
    "age_hist_het_gnomad_exomes",
    "age_hist_het_non_cancer_exomes",
    "age_hist_het_non_neuro_exomes",
    "age_hist_het_non_topmed_exomes",
    "age_hist_hom_controls_exomes",
    "age_hist_hom_gnomad_exomes",
    "age_hist_hom_non_cancer_exomes",
    "age_hist_hom_non_neuro_exomes",
    "age_hist_hom_non_topmed_exomes",
    "age_hist_het_controls_genomes",
    "age_hist_het_gnomad_genomes",
    "age_hist_het_non_neuro_genomes",
    "age_hist_het_non_topmed_genomes",
    "age_hist_hom_controls_genomes",
    "age_hist_hom_gnomad_genomes",
    "age_hist_hom_non_neuro_genomes",
    "age_hist_hom_non_topmed_genomes",
]


def v4_false_dup_unfurl_annotations(
    ht: hl.Table,
) -> hl.expr.StructExpression:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays.
    Adapted from work in /v4/create_release/validate_and_export_vcf.py

    The values of the returned dictionary are Hail Expressions describing how to access
    the corresponding values.

    :param ht: Table containing the nested variant annotation arrays to be unfurled.
    :return: StructExpression containing variant annotations and their corresponding
        expressions and updated entries and set of frequency entries to remove from the
        VCF.
    """
    expr_dict = {}
    data_types = ["exomes", "genomes"]

    # Unfurl exome and genomefreq index dict
    # Cycles through each key and index (e.g., k=afr_XX, i=31)
    for dt in data_types:
        logger.info(f"Unfurling freq data in {dt}...")
        freq_idx = hl.eval(ht[f"freq_index_dict_{dt}"])
        expr_dict.update(
            {
                f"{f if f != 'homozygote_count' else 'nhomalt'}_{k}_{dt}": ht[
                    f"v2_{dt}"
                ].freq[i][f]
                for k, i in freq_idx.items()
                for f in ht[f"v2_{dt}"].freq[0].keys()
            }
        )

    logger.info("Unfurling joint freq data...")
    joint_freq_idx = hl.eval(ht.joint_freq_index_dict)
    expr_dict.update(
        {
            f"{f if f != 'homozygote_count' else 'nhomalt'}_{k}_joint": (
                ht.v2_joint.joint_freq[i][f]
            )
            for k, i in joint_freq_idx.items()
            for f in ht.v2_joint.joint_freq[0].keys()
        }
    )

    for dt in data_types:
        logger.info(f"Adding grpmax data from {dt}...")
        grpmax_idx = hl.eval(ht[f"popmax_index_dict_{dt}"])
        grpmax_dict = {}
        grpmax_dict.update(
            {
                f"{f if f != 'homozygote_count' else 'nhomalt'}{k}_grpmax_{dt}": ht[
                    f"v2_{dt}"
                ].popmax[i][f]
                for k, i in grpmax_idx.items()
                for f in ht[f"v2_{dt}"].popmax[0].keys()
            }
        )

        expr_dict.update(grpmax_dict)

    # Create the fields for joint grpmax
    logger.info("Adding joint grpmax data...")
    joint_grpmax_idx = ht.v2_joint.joint_grpmax
    joint_grpmax_dict = {"grpmax_joint": joint_grpmax_idx.gen_anc}
    joint_grpmax_dict.update(
        {
            f"{f if f != 'homozygote_count' else 'nhomalt'}_grpmax_joint": (
                joint_grpmax_idx[f]
            )
            for f in [f for f in joint_grpmax_idx._fields if f != "gen_anc"]
        }
    )
    expr_dict.update(joint_grpmax_dict)

    for dt in data_types:
        logger.info(f"Unfurling faf data from {dt}...")
        faf_idx = hl.eval(ht[f"faf_index_dict_{dt}"])
        expr_dict.update(
            {
                f"{f}_{k}_{dt}": ht[f"v2_{dt}"].faf[i][f]
                for f in ht[f"v2_{dt}"].faf[0].keys()
                if "meta" not in f
                for k, i in faf_idx.items()
            }
        )

    logger.info("Unfurling joint faf data...")
    joint_faf_idx = hl.eval(ht.joint_faf_index_dict)
    expr_dict.update(
        {
            f"{f}_joint_{k}": ht.v2_joint.joint_faf[i][f]
            for f in ht.v2_joint.joint_faf[0].keys()
            for k, i in joint_faf_idx.items()
        }
    )

    # Note that v2->GRCh38 liftover files do not have fafmax data

    logger.info("Unfurling joint fafmax data...")
    joint_fafmax_idx = ht.v2_joint.joint_fafmax
    joint_fafmax_dict = {
        f"fafmax_{f if f != 'joint_fafmax_data_type' else 'data_type'}_joint": (
            joint_fafmax_idx[f]
        )
        for f in joint_fafmax_idx.keys()
    }
    expr_dict.update(joint_fafmax_dict)

    for dt in data_types:
        logger.info(f"Unfurling age hists data in {dt}...")
        hist_idx = ht[
            f"age_index_dict_{dt}"
        ]  # index of which hists are which #ht[f"v2_{dt}"]
        age_hists = ["age_hist_het", "age_hist_hom"]  # in ht.v2_{dt}
        age_hist_dict = {
            f"{hist}_{key_name}_{dt}": (
                hl.delimit(ht[f"v2_{dt}"][hist][index_val], delimiter="|")
                if "bin" in key_name
                else ht[f"v2_{dt}"][hist][index_val]
            )
            for hist in age_hists
            for key_name, index_val in hl.eval(hist_idx.items())
        }
        expr_dict.update(age_hist_dict)

    return hl.struct(**expr_dict)


def custom_make_info_expr(
    t: hl.Table,
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expression for variant annotations to be included in VCF INFO field.

    :param t: Table containing variant annotations to be reformatted for VCF export.
    :param hist_prefix: Prefix to use for histograms.
    :param data_type: Data type to make info expression for. One of "exomes" or
        "genomes". Default is "exomes".
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    """
    vcf_info_dict = {}

    for dt in ["exomes", "genomes"]:
        for field in FALSE_DUP_INFO:  # .remove('QUALapprox'):
            vcf_info_dict[f"{field.replace('info_','')}_{dt}"] = t[f"v2_{dt}"][
                f"{field}"
            ]
        for as_field in FALSE_DUP_AS:
            vcf_info_dict[f"as_{as_field}_{dt}"] = t[f"v2_{dt}"].allele_info[
                f"{as_field}"
            ]
        for vqsr_field in FALSE_DUP_AS_VQSR:
            vcf_info_dict[f"as_vqsr_{vqsr_field}_{dt}"] = t[f"v2_{dt}"].allele_info[
                f"{vqsr_field}"
            ]
        # Add region_flag and allele_info fields to info dict
        for at_field in ALLELE_TYPE_FIELDS:
            vcf_info_dict[f"{at_field}_{dt}"] = t[f"v2_{dt}"][f"{at_field}"]
        for r_field in REGION_FLAG_FIELDS:
            if "non_par" not in r_field:
                vcf_info_dict[f"{r_field}_{dt}"] = t[f"v2_{dt}"][f"{r_field}"]
            else:
                pass

        for hist in HISTS:
            hist_dict = {
                f"{hist}_bin_freq_{dt}": hl.delimit(
                    t[f"v2_{dt}"][hist].bin_freq, delimiter="|"
                ),
            }
            vcf_info_dict.update(hist_dict)

            if "dp" in hist:
                vcf_info_dict.update(
                    {f"{hist}_n_larger_{dt}": t[f"v2_{dt}"][hist].n_larger},
                )

        for age_hist in FALSE_DUP_AGE_HISTS:
            if dt in age_hist:
                age_hist_dict = {
                    f"{age_hist.replace(dt,'')}_bin_freq_{dt}": hl.delimit(
                        t.info[f"{age_hist}"].bin_freq, delimiter="|"
                    ),
                }

            vcf_info_dict.update(age_hist_dict)

        vcf_info_dict.update({f"vep_{dt}": t[f"v2_{dt}_vep"]})

    return vcf_info_dict


def _joint_filters(ht: hl.Table) -> hl.Table:
    """ "
    Annotate Hail Table with new 'filters' field, for status in exomes and genomes.

    :param ht: Hail Table with ht.v2_exomes.filters and h2.v2_genomes.filters
    :return: Hail Table with added ht.filters of type set<str>
    """
    ht = ht.annotate(
        exome_pass=(hl.len(ht.v2_exomes.filters) == 0)
        & (~hl.is_missing(ht.v2_exomes.filters)),
        genome_pass=(hl.len(ht.v2_genomes.filters) == 0)
        & (~hl.is_missing(ht.v2_genomes.filters)),
    )

    ht = ht.annotate(
        filters=hl.if_else(
            (ht.exome_pass) & (ht.genome_pass),
            {"PASS_EXOMES", "PASS_GENOMES"},
            hl.if_else(
                ht.exome_pass,
                {"PASS_EXOMES", "NO_PASS_GENOMES"},
                hl.if_else(
                    ht.genome_pass,
                    {"NO_PASS_EXOMES", "PASS_GENOMES"},
                    {"NO_PASS_EXOMES", "NO_PASS_GENOMES"},
                ),
            ),
        )
    )

    return ht


def prepare_false_dup_ht_for_validation(
    ht: hl.Table,
    vcf_info_reorder: Optional[List[str]] = None,
) -> hl.Table:
    """
    Prepare HT for validity checks and export.

    :param ht: Release Hail Table
    :param vcf_info_reorder: Order of VCF INFO fields
    :return: Hail Table prepared for validity checks and export
    """
    logger.info(
        "Unfurling nested gnomAD frequency annotations and add to INFO field..."
    )
    info_struct = v4_false_dup_unfurl_annotations(ht)

    logger.info("Constructing INFO field")
    # Remove SIFT and Polyphen from CSQ fields or they will be inserted with
    csq_fields = "|".join(
        [c for c in VEP_CSQ_FIELDS["101"].split("|") if c not in ["SIFT", "PolyPhen"]]
    )

    ht = ht.annotate(
        release_v2_exomes=ht.v2_exomes,
        release_v2_genomes=ht.v2_genomes,
        release_v2_joint=ht.v2_joint,
        info=info_struct,  # THIS IS THE BIG PART WOO
        v2_exomes_vep=vep_struct_to_csq(
            ht.v2_exomes.vep, csq_fields=csq_fields, has_polyphen_sift=False
        ),
        v2_genomes_vep=vep_struct_to_csq(
            ht.v2_genomes.vep, csq_fields=csq_fields, has_polyphen_sift=False
        ),
    )

    # Add variant annotations to INFO field
    ht = ht.annotate(info=ht.info.annotate(**custom_make_info_expr(ht)))

    # Add VEP-specific globals
    ht = ht.annotate_globals(
        vep_csq_header=process_vep_csq_header(VEP_CSQ_HEADER),
    )

    # Remove
    ht = ht.annotate(info=ht.info.drop(*FALSE_DUP_AGE_HISTS))

    ht = _joint_filters(ht)

    ht = ht.annotate(
        rsid=hl.if_else(
            ~hl.is_missing(ht.v2_exomes.rsid), ht.v2_exomes.rsid, ht.v2_genomes.rsid
        )
    )

    # Select relevant fields for VCF export
    ht = ht.select(
        "info", "filters", "rsid"
    )  # we'd want filters information, that is IN

    return ht


def main(args):
    ht = hl.read_table(get_false_dup_genes_path())
    ht = prepare_false_dup_ht_for_validation(ht)
    hl.export_vcf(ht, "gs://gnomad-tmp-4day/tester_with_filters.vcf.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Option to overwrite existing custom liftover table.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Option to store output at test paths.",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
