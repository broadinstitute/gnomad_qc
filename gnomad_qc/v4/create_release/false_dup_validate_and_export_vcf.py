"""Script to validate and export only the false duplication Hail Table of 3 genes to a VCF."""

import argparse
import logging
from copy import deepcopy

import hail as hl
from gnomad.resources.grch38.gnomad import HGDP_POPS, POPS, SUBSETS, TGP_POPS
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.annotations import (
    faf_expr,
    gen_anc_faf_max_expr,
    merge_freq_arrays,
    pop_max_expr,
)
from gnomad.utils.release import make_freq_index_dict_from_meta

from gnomad_qc.v2.annotations.generate_frequency_data import POPS_TO_REMOVE_FOR_POPMAX
from gnomad_qc.v2.resources.basics import get_gnomad_liftover_data_path
from gnomad_qc.v4.create_release.create_false_dup_liftover import FALSE_DUP_GENES
from gnomad_qc.v4.resources.release import get_false_dup_genes_path

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("false_dup_genes")
logger.setLevel(logging.INFO)

LIFTOVER_POPS = {
    "afr": "African/African-American",
    "amr": "Admixed American",
    "asj": "Ashkenazi Jewish",
    "eas": "East Asian",
    "eas_jpn": "Japanese",
    "eas_kor": "Korean",
    "eas_oea": "non-Korean, non-Japanese East Asian",
    "fin": "Finnish",
    "nfe": "Non-Finnish European",
    "nfe_bgr": "Non-Finnish European - Bulgarian",
    "nfe_est": "Non-Finnish European - Estonian",
    "nfe_nwe": "Non-Finnish European - North-Western European",
    "nfe_onf": "Non-Finnish European - Other Non-Finnish",
    "nfe_seu": "Non-Finnish European - Southern Europe",
    "nfe_swe": "Non-Finnish European - Swedish",
    "oth": "Other",
    "sas": "South Asian",
}

# TODO: is this different for exomes vs genomes ?
LIFTOVER_SUBSETS = {
    "v2_liftover": ["controls", "gnomad", "non_cancer", "non_neuro", "non_topmed"]
}

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
from gnomad.utils.vep import (  # CURRENT_VEP_VERSION,
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
                f"{f if f != 'homozygote_count' else 'nhomalt'}_{k}_{dt}": [
                    ht[f"v2_{dt}"].freq[i][f]
                ]
                for k, i in freq_idx.items()
                for f in ht[f"v2_{dt}"].freq[0].keys()
            }
        )

    logger.info("Unfurling joint freq data...")
    joint_freq_idx = hl.eval(ht.joint_freq_index_dict)
    expr_dict.update(
        {
            f"{f if f != 'homozygote_count' else 'nhomalt'}_{k}_joint": [
                (ht.v2_joint.joint_freq[i][f])
            ]
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
                f"{f if f != 'homozygote_count' else 'nhomalt'}{k}_grpmax_{dt}": [
                    ht[f"v2_{dt}"].popmax[i][f]
                ]
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
            f"{f if f != 'homozygote_count' else 'nhomalt'}_grpmax_joint": [
                (joint_grpmax_idx[f])
            ]
            for f in [f for f in joint_grpmax_idx._fields if f != "gen_anc"]
        }
    )
    expr_dict.update(joint_grpmax_dict)

    for dt in data_types:
        logger.info(f"Unfurling faf data from {dt}...")
        faf_idx = hl.eval(ht[f"faf_index_dict_{dt}"])
        expr_dict.update(
            {
                f"{f}_{k}_{dt}": [ht[f"v2_{dt}"].faf[i][f]]
                for f in ht[f"v2_{dt}"].faf[0].keys()
                if "meta" not in f
                for k, i in faf_idx.items()
            }
        )

    logger.info("Unfurling joint faf data...")
    joint_faf_idx = hl.eval(ht.joint_faf_index_dict)
    expr_dict.update(
        {
            f"{f}_joint_{k}": [ht.v2_joint.joint_faf[i][f]]
            for f in ht.v2_joint.joint_faf[0].keys()
            for k, i in joint_faf_idx.items()
        }
    )

    # Note that v2->GRCh38 liftover files do not have fafmax data

    logger.info("Unfurling joint fafmax data...")
    joint_fafmax_idx = ht.v2_joint.joint_fafmax
    joint_fafmax_dict = {
        f"fafmax_{f if f != 'joint_fafmax_data_type' else 'data_type'}_joint": [
            (joint_fafmax_idx[f])
        ]
        for f in joint_fafmax_idx.keys()
    }
    expr_dict.update(joint_fafmax_dict)

    for dt in data_types:
        logger.info(f"Unfurling age hists data in {dt}...")
        hist_idx = ht[
            f"age_index_dict_{dt}"
        ]  # index of which hists are which #ht[f"v2_{dt}"]
        age_hists = ["age_hist_het", "age_hist_hom"]
        for hist in age_hists:
            for key_name, index_val in hl.eval(hist_idx.items()):
                age_hist_dict = {
                    f"{hist}_{key_name}_bin_edges_{dt}": [
                        hl.delimit(
                            ht[f"v2_{dt}"][hist][index_val].bin_edges, delimiter="|"
                        )
                    ],
                    f"{hist}_{key_name}_bin_freq_{dt}": [
                        hl.delimit(
                            ht[f"v2_{dt}"][hist][index_val].bin_freq, delimiter="|"
                        )
                    ],
                    f"{hist}_{key_name}_n_smaller_{dt}": [
                        ht[f"v2_{dt}"][hist][index_val].n_smaller
                    ],
                    f"{hist}_{key_name}_n_larger_{dt}": [
                        ht[f"v2_{dt}"][hist][index_val].n_larger
                    ],
                }
                expr_dict.update(age_hist_dict)

    return hl.struct(**expr_dict)


def make_info_expr_false_dups(
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
            if field != "qual":
                vcf_info_dict[f"{field.replace('info_','')}_{dt}"] = [
                    t[f"v2_{dt}"][f"{field}"]
                ]
            else:
                vcf_info_dict[f"{field.replace('info_','')}_{dt}"] = t[f"v2_{dt}"][
                    f"{field}"
                ]
        for as_field in FALSE_DUP_AS:
            if as_field not in ["DS", "DB"]:
                vcf_info_dict[f"as_{as_field}_{dt}"] = [
                    t[f"v2_{dt}"].allele_info[f"{as_field}"]
                ]
            else:
                vcf_info_dict[f"as_{as_field}_{dt}"] = t[f"v2_{dt}"].allele_info[
                    f"{as_field}"
                ]
        # Not made into arrays. These are bools.
        for vqsr_field in FALSE_DUP_AS_VQSR:
            vcf_info_dict[f"as_vqsr_{vqsr_field}_{dt}"] = t[f"v2_{dt}"].allele_info[
                f"{vqsr_field}"
            ]
        # Add region_flag and allele_info fields to info dict
        # Not array. These are bools.
        for at_field in ALLELE_TYPE_FIELDS:
            vcf_info_dict[f"{at_field}_{dt}"] = t[f"v2_{dt}"][f"{at_field}"]
        # Also not array. These are bool.
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

        vcf_info_dict.update({f"vep_{dt}": t[f"v2_{dt}_vep"]})

        vcf_info_dict.update({f"filters_{dt}": t[f"v2_{dt}"].filters})

    return vcf_info_dict


def prepare_joint_filters(ht: hl.Table) -> hl.Table:
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
            {"PASS"},
            hl.if_else(
                ht.exome_pass,
                {"GENOMES_FILTERED"},
                hl.if_else(
                    ht.genome_pass,
                    {"EXOMES_FILTERED"},
                    {"EXOMES_FILTERED", "GENOMES_FILTERED"},
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
    ht = ht.annotate(info=ht.info.annotate(**make_info_expr_false_dups(ht)))

    # Add VEP-specific globals
    ht = ht.annotate_globals(
        vep_csq_header=process_vep_csq_header(VEP_CSQ_HEADER),
    )

    ht = prepare_joint_filters(ht)

    ht = ht.annotate(
        rsid=hl.if_else(
            ~hl.is_missing(ht.v2_exomes.rsid), ht.v2_exomes.rsid, ht.v2_genomes.rsid
        )
    )

    # Select relevant fields for VCF export
    ht = ht.select("info", "filters", "rsid")

    return ht


def prepare_vcf_liftover_header(unfurled_ht: hl.Table) -> Dict[str, Dict]:
    """
    Make custom dictionary to populate VCF Header with inteligible information.
    """
    v2_exomes_filter_dict = make_vcf_filter_dict(
        hl.eval(unfurled_ht.rf_exomes.rf_snv_cutoff.min_score),
        hl.eval(unfurled_ht.rf_exomes.rf_indel_cutoff.min_score),
        inbreeding_cutoff=-0.3,  # established
        variant_qc_filter="RF",
    )

    v2_genomes_filter_dict = make_vcf_filter_dict(
        hl.eval(unfurled_ht.rf_genomes.rf_snv_cutoff.min_score),
        hl.eval(unfurled_ht.rf_genomes.rf_indel_cutoff.min_score),
        inbreeding_cutoff=-0.3,  # established
        variant_qc_filter="RF",
    )

    ac0_text = v2_exomes_filter_dict["AC0"]["Description"]
    inbreeding_coeff_text = v2_exomes_filter_dict["InbreedingCoeff"]["Description"]
    rf_exomes_text = v2_exomes_filter_dict["RF"]["Description"]
    rf_genomes_text = v2_genomes_filter_dict["RF"]["Description"]

    custom_filter_dict = {
        "EXOMES_FILTERED": {
            "Description": (
                "Filtered out or not present in Exomes. Cause stored in"
                " the info annotation filters_exomes."
            )
        },
        "GENOMES_FILTERED": {
            "Description": (
                "Filtered out or not present in Genomes. Cause stored in"
                " the info annotation filters_genomes."
            )
        },
        "PASS": {"Description": "Passed all variant filters"},
        "AC0": {
            "Description": f"In INFO exomes_filters and genomes_filters: {ac0_text}"
        },
        "InbreedingCoeff": {
            "Description": f"In INFO exomes_filters and genomes_filters: {inbreeding_coeff_text}"
        },
        "RF": {
            "Description": f"In exomes for exomes_filters: {rf_exomes_text}. In genomes for genomes_filters: {rf_genomes_text}"
        },
    }

    vcf_info_dict = populate_info_dict(list(unfurled_ht.info))

    for dt in ["exomes", "genomes"]:
        vcf_info_dict[f"filters_{dt}"] = {
            "Description:": f"Reasons for variant failure in {dt}. See #FILTER for more",
            "Type": "String",
            "Number": ".",
        }

    return {
        "info": vcf_info_dict,
        "filter": custom_filter_dict,
    }


def populate_subset_info_dict(
    subset: str,
    description_text: str = "",
    data_type: str = "exomes",
    pops: Dict[str, str] = LIFTOVER_POPS,
    faf_pops: Dict[str, List[str]] = FAF_POPS,
    sexes: List[str] = ["male", "female"],
    label_delimiter: str = "_",
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` to populate INFO dictionary for the requested `subset`.

    Creates:
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample genetic
            ancestry group, sex both for adj and raw data.
        - INFO fields for filtering allele frequency (faf) annotations.

    :param subset: Sample subset in dataset. "" is used as a placeholder for the full
        dataset.
    :param description_text: Text describing the sample subset that should be added to
        the INFO description.
    :param data_type: One of "exomes" or "genomes". Default is "exomes".
    :param pops: Dict of sample global genetic ancestry names for the gnomAD data type.
        Default is POPS.
    :param faf_pops: Dict with gnomAD version (keys) and faf genentic ancestry group
        names (values). Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param label_delimiter: String to use as delimiter when making group label
        combinations. Default is '_'.
    :return: Dictionary containing Subset specific INFO header fields.
    """
    vcf_info_dict = {}
    # Remove unnecessary pop names from FAF_POPS dict depending on data type
    # and version of FAF_POPS.
    faf_pops_version = "v4" if data_type == "exomes" else "v3"
    faf_pops = {pop: POP_NAMES[pop] for pop in faf_pops[faf_pops_version]}

    # Add FAF fields to dict.
    faf_label_groups = create_label_groups(
        pops=faf_pops, sexes=sexes, all_groups=["adj"]
    )
    for label_group in faf_label_groups:
        faf_dict = {}
        faf_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                pop_names=faf_pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                faf=True,
                description_text=description_text,
                suffix=data_type,
            )
        )
        faf_dict_fix = {k.replace("_adj_", "_"): v for k, v in faf_dict.items()}
        vcf_info_dict.update(faf_dict_fix)
    # Add AC, AN, AF, nhomalt fields to dict.
    label_groups = create_label_groups(pops=pops, sexes=sexes)
    for label_group in label_groups:
        label_group_dict = {}
        label_group_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                pop_names=pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                description_text=description_text,
                callstats=True,
                suffix=data_type,
            )
        )
        vcf_info_dict.update(
            {k.replace("_adj_", "_"): v for k, v in label_group_dict.items()}
        )

    # Add grpmax.
    grpmax_dict = make_info_dict(
        suffix=f"{subset}_{data_type}",
        label_delimiter=label_delimiter,
        pop_names=pops,
        grpmax=True,
        description_text=description_text,
    )
    vcf_info_dict.update(
        {
            k.replace(
                f"grpmax_{subset}",
                f"{subset}_grpmax",
            ): v
            for k, v in grpmax_dict.items()
        }
    )

    return vcf_info_dict


def populate_info_dict(
    info_fields: List[str],
    bin_edges: Dict[str, str] = None,
    age_hist_distribution: str = None,
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT,
    subset_list: List[str] = LIFTOVER_SUBSETS,
    pops: Dict[str, str] = LIFTOVER_POPS,
    faf_pops: Dict[str, List[str]] = FAF_POPS,
    sexes: List[str] = SEXES,
    in_silico_dict: Dict[str, Dict[str, str]] = IN_SILICO_ANNOTATIONS_INFO_DICT,
    vrs_fields_dict: Dict[str, Dict[str, str]] = VRS_FIELDS_DICT,
    label_delimiter: str = "_",
    data_type: str = "exomes",
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` and `make_hist_dict` to populate INFO dictionary.

    Used during VCF export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for
          heterozygous and homozygous variant carriers).
        - INFO fields for grpmax AC, AN, AF, nhomalt, and grpmax genetic ancestry group.
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample genetic
          ancestry group, sex both for adj and raw data.
        - INFO fields for filtering allele frequency (faf) annotations.
        - INFO fields for variant histograms (hist_bin_freq for each histogram and
          hist_n_larger for DP histograms).

    :param info_fields: List of info fields to add to the info dict. Default is None.
    :param bin_edges: Dictionary of variant annotation histograms and their associated
        bin edges.
    :param age_hist_distribution: Pipe-delimited string of overall age histogram bin
        frequency.
    :param info_dict: INFO dict to be populated.
    :param subset_list: List of sample subsets in dataset. Default is SUBSETS.
    :param pops: Dict of sample global genetic ancestry names for the gnomAD data type.
    :param faf_pops: Dict with gnomAD version (keys) and faf genentic ancestry group
        names (values). Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param in_silico_dict: Dictionary of in silico predictor score descriptions.
    :param vrs_fields_dict: Dictionary with VRS annotations.
    :param label_delimiter: String to use as delimiter when making group label
        combinations.
    :param data_type: Data type to populate info dict for. One of "exomes" or
        "genomes". Default is "exomes".
    :return: Updated INFO dictionary for VCF export.
    """
    # Get existing info fields from predefined info_dict, e.g. `FS`,
    # `non_par`, `negative_train_site`...
    vcf_info_dict = info_dict.copy()  # this does NOT contain "_exomes" or "_genomes" in

    # Add allele-specific fields to info dict, including AS_VQSR_FIELDS
    # this ought to be fine as is, I think!
    vcf_info_dict.update(
        add_as_info_dict(info_dict=info_dict, as_fields=AS_FIELDS + AS_VQSR_FIELDS)
    )

    for dt in ["exomes", "genomes", "joint"]:
        # so the names don't exist with "_exomes" or "_genomes" inside of the vcf_info_dict
        # but they do inside of info_fields

        # for each row in the MT:
        if dt != "joint":
            for f in info_fields:

                # de-data type it
                f_sans_dt = f.replace(
                    f"_{dt}", ""
                )  # this is how it appears in the texty vcf_info_dict

                # if it's non-AS and we have a description for it, add it in!
                if f_sans_dt in list(vcf_info_dict.keys()):
                    appending_struct = vcf_info_dict[f_sans_dt].copy()
                    appending_struct["Description"] += f" in {dt}"
                    vcf_info_dict[f] = appending_struct

                # AS fix, needed to grab AS info
                if f_sans_dt[:3] == "as_":
                    f_sans_dt = f_sans_dt.replace(
                        "as_", "AS_", 1
                    )  # only replace first instance

                    if f_sans_dt in list(vcf_info_dict.keys()):
                        appending_struct = vcf_info_dict[f_sans_dt].copy()
                        appending_struct["Description"] += f" - in {dt}"
                        vcf_info_dict[f] = appending_struct

            for subset in subset_list["v2_liftover"]:
                subset_pops = deepcopy(pops)
                if (subset == "joint") | (data_type == "genomes"):
                    subset_pops.update({"ami": "Amish"})
                description_text = "" if subset == "" else f" in {subset} subset"

                subset_info = populate_subset_info_dict(
                    subset=subset,
                    description_text=description_text,
                    data_type=dt,
                    pops=subset_pops,
                    faf_pops=faf_pops,
                    sexes=["male", "female"],
                    label_delimiter=label_delimiter,
                )

                vcf_info_dict.update(
                    {k.replace("v2_liftover_", ""): v for k, v in subset_info.items()}
                )

    if age_hist_distribution:
        age_hist_distribution = "|".join(str(x) for x in age_hist_distribution)
    else:
        logger.info("Punting on age_hists...")

    return vcf_info_dict


def main(args):
    ht = hl.read_table(get_false_dup_genes_path(release_version="4.0"))
    ht = prepare_false_dup_ht_for_validation(ht)
    appended_header = prepare_vcf_liftover_header(ht)

    export_path = get_false_dup_genes_path(test=args.test).replace(".ht", ".vcf.bgz")
    logger.info(f"Writing to: {export_path}...")
    hl.export_vcf(ht, export_path, metadata=appended_header)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Option to store output at test paths.",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
