"""Script to validate and export gnomAD VCFs."""

import argparse
import logging
import pickle
from copy import deepcopy
from pprint import pprint
from typing import Dict, List, Optional, Set, Tuple

import hail as hl
from gnomad.assessment.validity_checks import (
    check_global_and_row_annot_lengths,
    pprint_global_anns,
    validate_release_t,
    vcf_field_check,
)
from gnomad.resources.grch38.gnomad import HGDP_POPS, POPS, SUBSETS, TGP_POPS
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.filtering import remove_fields_from_constant
from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    AS_VQSR_FIELDS,
    FAF_POPS,
    HISTS,
    IN_SILICO_ANNOTATIONS_INFO_DICT,
    INFO_DICT,
    JOINT_REGION_FLAG_FIELDS,
    JOINT_REGION_FLAGS_INFO_DICT,
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
    CURRENT_VEP_VERSION,
    VEP_CSQ_FIELDS,
    VEP_CSQ_HEADER,
    vep_struct_to_csq,
)

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.resources.basics import get_logging_path, qc_temp_prefix
from gnomad_qc.v4.resources.release import (
    append_to_vcf_header_path,
    release_header_path,
    release_sites,
    release_vcf_path,
    validated_release_ht,
)

# Add new site fields
NEW_SITE_FIELDS = [
    "monoallelic",
    "only_het",
    "transmitted_singleton",
]
SITE_FIELDS = deepcopy(SITE_FIELDS)
SITE_FIELDS.extend(NEW_SITE_FIELDS)
SITE_FIELDS = {
    "exomes": SITE_FIELDS + ["sibling_singleton"],
    "genomes": SITE_FIELDS,
}

# Add updatedVQSR fields
NEW_AS_VQSR_FIELDS = ["negative_train_site", "positive_train_site"]
AS_VQSR_FIELDSS = deepcopy(AS_VQSR_FIELDS)
AS_VQSR_FIELDS.extend(NEW_AS_VQSR_FIELDS)

# Update inbreeding coeff
MISSING_AS_FIELDS = ["InbreedingCoeff"]
AS_FIELDS = remove_fields_from_constant(AS_FIELDS, ["InbreedingCoeff"])
AS_FIELDS.append("inbreeding_coeff")

# Drop decoy, still doesn't exist on 38
REGION_FLAG_FIELDS = deepcopy(REGION_FLAG_FIELDS)
REGION_FLAG_FIELDS = remove_fields_from_constant(
    REGION_FLAG_FIELDS, ["decoy", "nonpar"]
)
REGION_FLAG_FIELDS = {
    "exomes": REGION_FLAG_FIELDS + [
        "fail_interval_qc",
        "outside_ukb_capture_region",
        "outside_broad_capture_region",
    ],
    "genomes": REGION_FLAG_FIELDS,
    "joint": JOINT_REGION_FLAG_FIELDS,
}

# Remove original alleles for containing non-releasable alleles
ALLELE_TYPE_FIELDS = deepcopy(ALLELE_TYPE_FIELDS)
ALLELE_TYPE_FIELDS = remove_fields_from_constant(
    ALLELE_TYPE_FIELDS, ["original_alleles"]
)
# Remove original alleles for containing non-releasable alleles
ALLELE_TYPE_FIELDS = {
    "exomes": ALLELE_TYPE_FIELDS,
    "genomes": remove_fields_from_constant(ALLELE_TYPE_FIELDS, ["has_star"]),
}

INSILICO_FIELDS = [
    "cadd",
    "revel_max",
    "spliceai_ds_max",
    "pangolin_largest_ds",
    "phylop",
    "sift_max",
    "polyphen_max",
]

GENOME_SUBSETS_TO_DROP = remove_fields_from_constant(
    deepcopy(SUBSETS["v3"]), ["hgdp", "tgp"]
)
SUBSETS = {
    "exomes": deepcopy(SUBSETS["v4"]),
    "genomes": remove_fields_from_constant(
        deepcopy(SUBSETS["v3"]), GENOME_SUBSETS_TO_DROP
    ),
    "joint": [""],
}

POPS = deepcopy(POPS["v4"])
POPS["joint"] = set(POPS["exomes"]) | set(POPS["genomes"])

# Remove unnecessary pop names from POP_NAMES dict
POPS = {d: {pop: POP_NAMES[pop] for pop in pops} for d, pops in POPS.items()}

SAMPLE_SUM_SETS_AND_POPS = {
    "exomes": {"non_ukb": POPS["exomes"]},
    "genomes": {"hgdp": HGDP_POPS, "tgp": TGP_POPS},
    "joint": None,
}

# Row annotaions and their associated global annotations for length comparison
LEN_COMP_GLOBAL_ROWS = {
    "freq": ["freq_meta", "freq_index_dict", "freq_meta_sample_count"],
    "faf": ["faf_meta", "faf_index_dict"],
}

# Joint row annotations and their associated global annotations for length comparison
LEN_COMP_JOINT_GLOBAL_ROWS = {
    "joint_freq": [
        "joint_freq_meta",
        "joint_freq_index_dict",
        "joint_freq_meta_sample_count",
    ],
    "joint_faf": ["joint_faf_meta", "joint_faf_index_dict"],
}


# VCF INFO fields to reorder.
VCF_INFO_REORDER = [
    "AC",
    "AN",
    "AF",
    "grpmax",
    "fafmax_faf95_max",
    "fafmax_faf95_max_gen_anc",
]


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_export_resources(
    overwrite: bool = False,
    data_type: str = "exomes",
    test: Optional[bool] = False,
    contig: Optional[str] = None,
) -> PipelineResourceCollection:
    """
    Get export resources.

    :param overwrite: Whether to overwrite existing files.
    :param data_type: Data type to get resources for. One of "exomes" or "genomes".
        Default is "exomes".
    :param test: Whether to use test resources.
    :param contig: Contig to get resources for. Default is None.
    :return: Export resources.
    """
    export_pipeline = PipelineResourceCollection(
        pipeline_name="validate_and_export_vcf",
        overwrite=overwrite,
    )
    validate_release_ht = PipelineStepResourceCollection(
        "--validate-release-ht",
        input_resources={
            "create_release_sites_ht.py": {
                "release_ht": release_sites(data_type=data_type, test=test)
            }
        },
        output_resources={
            "validated_ht": validated_release_ht(test=test, data_type=data_type)
        },
    )
    prepare_vcf_header = PipelineStepResourceCollection(
        "--prepare-vcf-header",
        pipeline_input_steps=[validate_release_ht],
        add_input_resources={
            "create_release_sites_ht.py": {
                "release_ht": release_sites(data_type=data_type)
            }
        },
        output_resources={
            "vcf_header_path": release_header_path(test=test, data_type=data_type)
        },
    )
    export_vcf = PipelineStepResourceCollection(
        "--export-vcf",
        pipeline_input_steps=[validate_release_ht, prepare_vcf_header],
        output_resources={
            "release_vcf": release_vcf_path(
                test=test,
                data_type=data_type,
                contig=contig,
            )
        },
    )

    export_pipeline.add_steps(
        {
            "validate_release_ht": validate_release_ht,
            "prepare_vcf_header": prepare_vcf_header,
            "export_vcf": export_vcf,
        }
    )
    return export_pipeline


def filter_to_test(ht: hl.Table, num_partitions: int = 2) -> hl.Table:
    """
    Filter Table to `num_partitions` partitions on chr20, chrX, and chrY for testing.

    :param ht: Input Table to filter.
    :param num_partitions: Number of partitions to grab from each chromosome.
    :return: Input Table filtered to `num_partitions` on chr20, chrX, and chrY.
    """
    logger.info(
        "Filtering to %d partitions on chr20, chrX, and chrY (for tests only)...",
        num_partitions,
    )
    ht_chr20 = hl.filter_intervals(ht, [hl.parse_locus_interval("chr20")])
    ht_chr20 = ht_chr20._filter_partitions(range(num_partitions))

    ht_chrx = hl.filter_intervals(ht, [hl.parse_locus_interval("chrX")])
    ht_chrx = ht_chrx._filter_partitions(range(num_partitions))

    ht_chry = hl.filter_intervals(ht, [hl.parse_locus_interval("chrY")])
    ht_chry = ht_chry._filter_partitions(range(num_partitions))
    ht = ht_chr20.union(ht_chrx, ht_chry)

    return ht


def select_type_from_joint_ht(ht: hl.Table, data_type: str) -> hl.Table:
    """
    Select all fields from the joint HT that are relevant to `data_type`.

    :param ht: Joint release HT.
    :param data_type: Data type to select in joint HT. One of "exomes", "genomes", or
        "joint".
    :return: Joint HT with fields relevant to `data_type`.
    """
    global_fields = [f"{data_type}_globals"]
    row_fields = [data_type, "region_flags"]
    if data_type == "joint":
        row_fields.append("freq_comparison_stats")
    ht = ht.select_globals(*global_fields)
    ht = ht.select(*row_fields)
    ht = ht.transmute_globals(**ht[f"{data_type}_globals"])
    ht = ht.transmute(**ht[data_type])
    return ht


def unfurl_nested_annotations(
    ht: hl.Table,
    entries_to_remove: Set[str] = None,
    data_type: str = "exomes",
    joint_included: bool = False,
    hist_prefix: str = "",
    freq_comparison_included: bool = False,
    for_joint_validation: bool = False,
) -> [hl.expr.StructExpression, Set[str], Dict[str, str]]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays.

    The values of the returned dictionary are Hail Expressions describing how to access
    the corresponding values.

    :param ht: Table containing the nested variant annotation arrays to be unfurled.
    :param entries_to_remove: Optional Set of frequency entries to remove for vcf_export.
    :param data_type: Data type to unfurl nested annotations for. One of "exomes",
        "genomes", or "joint". Default is "exomes".
    :param joint_included: Whether joint frequency data is included in the exome or
        genome HT. Default is False.
    :param hist_prefix: Prefix to use for histograms. Default is "".
    :param freq_comparison_included: Whether frequency comparison data is included in
        the HT. Default is False.
    :param for_joint_validation: Whether to prepare HT for joint validation. Default is
        False.
    :return: StructExpression containing variant annotations and their corresponding
        expressions and updated entries, set of frequency entries to remove from the
        VCF, and a dict of fields to rename when `for_joint_validation` is True.
    """
    expr_dict = {}
    rename_dict = {}

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=afr_XX, i=31)
    logger.info("Unfurling freq data...")
    freq_idx = hl.eval(ht.freq_index_dict)
    for k, i in freq_idx.items():
        for f in ht.freq[0].keys():
            field_name = f if f != "homozygote_count" else "nhomalt"
            expr_dict[f"{field_name}_{k}"] = ht.freq[i][f]
            if for_joint_validation:
                rename_dict[f"{field_name}_{k}"] = f"{field_name}_{data_type}_{k}"

    if joint_included:
        logger.info("Unfurling joint freq data...")
        joint_freq_idx = hl.eval(ht.joint_freq_index_dict)
        expr_dict.update(
            {
                f"{f if f != 'homozygote_count' else 'nhomalt'}_joint_{k}": (
                    ht.joint_freq[i][f]
                )
                for k, i in joint_freq_idx.items()
                for f in ht.joint_freq[0].keys()
            }
        )

    # This creates fields like grpmax, AC_grpmax_non_ukb...
    logger.info("Adding grpmax data...")
    grpmax_idx = ht.grpmax
    if data_type == "exomes" and not for_joint_validation:
        grpmax_dict = {}
        for s in grpmax_idx.keys():
            grpmax_dict.update(
                {f"grpmax{'_'+s if s != 'gnomad' else ''}": grpmax_idx[s].gen_anc}
            )
            grpmax_dict.update(
                {
                    f"{f if f != 'homozygote_count' else 'nhomalt'}_grpmax{'_'+s if s != 'gnomad' else ''}": grpmax_idx[
                        s
                    ][
                        f
                    ]
                    for f in grpmax_idx[s].keys()
                    if f != "gen_anc"
                }
            )
    else:
        grpmax_dict = {"grpmax": grpmax_idx.gen_anc}
        grpmax_rename = {
            f: f if f != "homozygote_count" else "nhomalt"
            for f in grpmax_idx.keys()
            if f != "gen_anc"
        }
        grpmax_dict.update(
            {f"{v}_grpmax": grpmax_idx[k] for k, v in grpmax_rename.items()}
        )
        if for_joint_validation:
            rename_dict["grpmax"] = f"grpmax_{data_type}"
            rename_dict.update(
                {
                    f"{v}_grpmax": f"{v}_grpmax_{data_type}"
                    for v in grpmax_rename.values()
                }
            )
    expr_dict.update(grpmax_dict)

    if joint_included:
        logger.info("Adding joint grpmax data...")
        joint_grpmax_idx = ht.joint_grpmax
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

    logger.info("Unfurling faf data...")
    faf_idx = hl.eval(ht.faf_index_dict)
    for k, i in faf_idx.items():
        for f in ht.faf[0].keys():
            expr_dict[f"{f}_{k}"] = ht.faf[i][f]
            if for_joint_validation:
                rename_dict[f"{f}_{k}"] = f"{f}_{data_type}_{k}"

    logger.info("Unfurling fafmax data...")
    fafmax_idx = ht.fafmax
    if data_type == "exomes" and not for_joint_validation:
        fafmax_dict = {
            f"fafmax_{f}{'_'+s if s != 'gnomad' else ''}": fafmax_idx[s][f]
            for s in fafmax_idx.keys()
            for f in fafmax_idx[s].keys()
        }
    else:
        fafmax_dict = {f"fafmax_{f}": fafmax_idx[f] for f in fafmax_idx.keys()}
        if for_joint_validation:
            rename_dict.update(
                {f"fafmax_{f}": f"fafmax_{f}_{data_type}" for f in fafmax_idx.keys()}
            )
    expr_dict.update(fafmax_dict)

    if joint_included:
        logger.info("Unfurling joint faf data...")
        joint_faf_idx = hl.eval(ht.joint_faf_index_dict)
        expr_dict.update(
            {
                f"{f}_joint_{k}": ht.joint_faf[i][f]
                for f in ht.joint_faf[0].keys()
                for k, i in joint_faf_idx.items()
            }
        )

        logger.info("Unfurling joint fafmax data...")
        joint_fafmax_idx = ht.joint_fafmax
        joint_fafmax_dict = {
            f"fafmax_{f if f != 'joint_fafmax_data_type' else 'data_type'}_joint": (
                joint_fafmax_idx[f]
            )
            for f in joint_fafmax_idx.keys()
        }
        expr_dict.update(joint_fafmax_dict)

    logger.info("Unfurling age hists...")
    hist_idx = ht.histograms.age_hists
    age_hists = ["age_hist_het", "age_hist_hom"]
    for hist in age_hists:
        for f in hist_idx[hist].keys():
            expr_dict[f"{hist}_{f}"] = (
                hl.delimit(hist_idx[hist][f], delimiter="|")
                if "bin" in f
                else hist_idx[hist][f]
            )
            if for_joint_validation:
                rename_dict[f"{hist}_{f}"] = f"{hist}_{f}_{data_type}"

    logger.info("Unfurling variant quality histograms...")
    # Add underscore to hist_prefix if it isn't empty
    if hist_prefix != "":
        hist_prefix += "_"

    # Histograms to export are:
    # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
    # We previously dropped:
    # _n_smaller for all hists
    # _bin_edges for all hists
    # _n_larger for all hists EXCEPT DP hists
    for hist in HISTS:
        hist_dict = {
            f"{hist}_bin_freq": hl.delimit(
                ht.histograms.qual_hists[hist].bin_freq, delimiter="|"
            ),
        }
        expr_dict.update(hist_dict)
        if for_joint_validation:
            rename_dict.update(
                {f"{hist}_bin_freq": f"{hist_prefix}{hist}_bin_freq_{data_type}"}
            )

        if "dp" in hist:
            expr_dict.update(
                {f"{hist}_n_larger": ht.histograms.qual_hists[hist].n_larger},
            )
            if for_joint_validation:
                rename_dict.update(
                    {f"{hist}_n_larger": f"{hist_prefix}{hist}_n_larger_{data_type}"}
                )

    if freq_comparison_included:
        logger.info("Unfurling contingency table test results...")
        contingency_idx = hl.eval(ht.freq_index_dict)
        for k, i in contingency_idx.items():
            for f in ht.freq_comparison_stats.contingency_table_test[0].keys():
                key = f"CTT_{f}_{k}"
                expr = ht.freq_comparison_stats.contingency_table_test[i][f]
                expr_dict[key] = expr

        logger.info("Unfurling Cochran-Mantel-Haenszel test results...")
        expr_dict["CMH_chisq"] = (
            ht.freq_comparison_stats.cochran_mantel_haenszel_test.chisq
        )
        expr_dict["CMH_p_value"] = (
            ht.freq_comparison_stats.cochran_mantel_haenszel_test.p_value
        )
        logger.info("Unfurling unionized stats...")
        expr_dict["stat_union_p_value"] = ht.freq_comparison_stats.stat_union.p_value
        expr_dict["stat_union_test_name"] = (
            ht.freq_comparison_stats.stat_union.stat_test_name
        )
        expr_dict["stat_union_gen_ancs"] = ht.freq_comparison_stats.stat_union.gen_ancs

    return hl.struct(**expr_dict), entries_to_remove, rename_dict


def make_info_expr(
    t: hl.Table,
    data_type: str = "exomes",
    for_joint_validation: bool = False,
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expression for variant annotations to be included in VCF INFO field.

    :param t: Table containing variant annotations to be reformatted for VCF export.
    :param data_type: Data type to make info expression for. One of "exomes", "genomes",
        or "joint". Default is "exomes".
    :param for_joint_validation: Whether to prepare HT for joint validation. Default is False.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    """
    vcf_info_dict = {}

    # Set data type to joint if for_joint_validation is True so the correct region flag
    # fields are used.
    if for_joint_validation:
        data_type = "joint"

    if "region_flags" in t.row:
        # Add region_flag to info dict
        for field in REGION_FLAG_FIELDS[data_type]:
            vcf_info_dict[field] = t["region_flags"][f"{field}"]

    if for_joint_validation:
        return vcf_info_dict

    # Add site-level annotations and AS annotations to vcf_info_dict
    for field in SITE_FIELDS[data_type] + AS_FIELDS:
        vcf_info_dict[field] = t["release_ht_info"][f"{field}"]

    for field in AS_VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr_results"][f"{field}"]

    # Add allele_info fields to info dict
    for field in ALLELE_TYPE_FIELDS[data_type]:
        vcf_info_dict[field] = t["allele_info"][f"{field}"]

    # Add in silico annotations to info dict
    insilico_idx = t.in_silico_predictors
    for field in INSILICO_FIELDS:
        if field == "cadd":
            vcf_info_dict[f"{field}_raw_score"] = insilico_idx[field]["raw_score"]
            vcf_info_dict[f"{field}_phred"] = insilico_idx[field]["phred"]
        else:
            vcf_info_dict[field] = insilico_idx[field]

    # Add VRS annotations to info dict
    for field in VRS_FIELDS_DICT:
        vcf_info_dict[field] = t["release_ht_info"]["vrs"][field]

    # Add vep annotations to info dict
    vcf_info_dict["vep"] = t["vep"]

    return vcf_info_dict


def prepare_ht_for_validation(
    ht: hl.Table,
    data_type: str = "exomes",
    freq_entries_to_remove: Optional[List[str]] = None,
    vcf_info_reorder: Optional[List[str]] = None,
    joint_included: bool = False,
    for_joint_validation: bool = True,
    freq_comparison_included: bool = False,
) -> hl.Table:
    """
    Prepare HT for validity checks and export.

    :param ht: Release Hail Table.
    :param data_type: Data type to prepare HT for. One of "exomes", "genomes", or
        "joint". Default is "exomes".
    :param freq_entries_to_remove: List of entries to remove from freq.
    :param vcf_info_reorder: Order of VCF INFO fields.
    :param joint_included: Whether joint frequency data is included in the HT. Default
        is False.
    :param for_joint_validation: Whether to prepare HT for joint validation. Default is
        True.
    :param freq_comparison_included: Whether frequency comparison data is included in
        the HT. Default is False.
    :return: Hail Table prepared for validity checks and export and a dictionary of
        fields to rename when `for_joint_validation` is True.
    """
    logger.info(
        "Unfurling nested gnomAD frequency annotations and add to INFO field..."
    )
    info_struct, freq_entries_to_remove, rename_dict = unfurl_nested_annotations(
        ht,
        entries_to_remove=freq_entries_to_remove,
        data_type=data_type,
        joint_included=joint_included,
        freq_comparison_included=freq_comparison_included,
        for_joint_validation=for_joint_validation,
    )

    logger.info("Constructing INFO field")
    # Remove SIFT and Polyphen from CSQ fields or they will be inserted with
    # missing values by vep_struct_to_csq. These fields are processed separately
    # as in silico annotations.
    csq_fields = "|".join(
        [
            c
            for c in VEP_CSQ_FIELDS[CURRENT_VEP_VERSION].split("|")
            if c != "PolyPhen" and c != "SIFT"
        ]
    )

    if for_joint_validation:
        ann_expr = {"info": info_struct}
        if "region_flags" in ht.row:
            ann_expr["region_flags"] = ht.region_flags
    else:
        ann_expr = {
            "region_flag": ht.region_flags,
            "release_ht_info": ht.info,
            "info": info_struct,
            "rsid": hl.str(";").join(ht.rsid),
            "vep": vep_struct_to_csq(
                ht.vep, csq_fields=csq_fields, has_polyphen_sift=False
            ),
        }

    ht = ht.annotate(**ann_expr)

    # Add variant annotations to INFO field
    # This adds the following:
    #   region flag for problematic regions
    #   annotations in ht.release_ht_info (site and allele-specific annotations),
    #   info struct (unfurled data obtained above),
    #   dbSNP rsIDs
    #   all VEP annotations
    ht = ht.annotate(
        info=ht.info.annotate(
            **make_info_expr(
                ht,
                data_type=data_type,
                for_joint_validation=for_joint_validation,
            )
        )
    )

    if for_joint_validation:
        ht = ht.annotate_globals(
            freq_entries_to_remove=(
                freq_entries_to_remove
                if freq_entries_to_remove
                else hl.empty_set(hl.tstr)
            ),
        )
    else:
        ht = ht.annotate_globals(
            vep_csq_header=process_vep_csq_header(VEP_CSQ_HEADER),
            freq_entries_to_remove=(
                freq_entries_to_remove
                if freq_entries_to_remove
                else hl.empty_set(hl.tstr)
            ),
        )

    # Select relevant fields for VCF export
    if for_joint_validation:
        if "filters" in ht.row:
            filters_expr = ht.filters
        else:
            filters_expr = hl.missing(hl.tset(hl.tstr))
        ht = ht.select("info", filters=filters_expr)
    else:
        ht = ht.select("info", "filters", "rsid")

    if vcf_info_reorder:
        logger.info("Rearranging fields to desired order...")
        ht = ht.annotate(
            info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder))
        )

    return ht, rename_dict


def populate_subset_info_dict(
    subset: str,
    description_text: str,
    data_type: str = "exomes",
    pops: Dict[str, str] = POPS["exomes"],
    faf_pops: Dict[str, List[str]] = FAF_POPS,
    sexes: List[str] = SEXES,
    label_delimiter: str = "_",
    freq_comparison_included: bool = False,
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
    :param data_type: One of "exomes", "genomes", or "joint". Default is "exomes".
    :param faf_pops: Dict with gnomAD version (keys) and faf genentic ancestry group
        names (values). Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param label_delimiter: String to use as delimiter when making group label
        combinations. Default is '_'.
    :param freq_comparison_included: Whether frequency comparison data is included in the HT.
    :return: Dictionary containing Subset specific INFO header fields.
    """
    vcf_info_dict = {}
    # Remove unnecessary pop names from FAF_POPS dict depending on data type
    # and version of FAF_POPS.
    faf_pops_version = "v3" if data_type == "genomes" or subset == "genomes" else "v4"
    faf_pops = {pop: POP_NAMES[pop] for pop in faf_pops[faf_pops_version]}

    # Add FAF fields to dict.
    faf_label_groups = create_label_groups(
        pops=faf_pops, sexes=sexes, all_groups=["adj"]
    )
    for label_group in faf_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                pop_names=faf_pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                faf=True,
                description_text=description_text,
            )
        )
    # Add AC, AN, AF, nhomalt fields to dict.
    label_groups = create_label_groups(pops=pops, sexes=sexes)
    for label_group in label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                pop_names=pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                description_text=description_text,
                callstats=True,
            )
        )

    # Add grpmax.
    vcf_info_dict.update(
        make_info_dict(
            suffix=subset,
            label_delimiter=label_delimiter,
            pop_names=pops,
            grpmax=True,
            description_text=description_text,
        )
    )

    # Add fafmax.
    vcf_info_dict.update(
        make_info_dict(
            suffix=subset,
            label_delimiter=label_delimiter,
            pop_names=pops,
            fafmax=True,
            description_text=description_text,
        )
    )
    if freq_comparison_included:
        ctt_label_groups = create_label_groups(pops=pops, sexes=sexes)
        for label_group in ctt_label_groups:
            vcf_info_dict.update(
                make_info_dict(
                    freq_ctt=True,
                    label_groups=label_group,
                )
            )
        vcf_info_dict.update(
            make_info_dict(
                freq_cmh=True,
            )
        )
        vcf_info_dict.update(
            make_info_dict(
                freq_stat_union=True,
            )
        )

    return vcf_info_dict


def populate_info_dict(
    info_fields: List[str],
    bin_edges: Dict[str, str],
    age_hist_distribution: str = None,
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT,
    subset_list: List[str] = SUBSETS["exomes"],
    pops: Dict[str, str] = POPS["exomes"],
    faf_pops: Dict[str, List[str]] = FAF_POPS,
    sexes: List[str] = SEXES,
    in_silico_dict: Dict[str, Dict[str, str]] = IN_SILICO_ANNOTATIONS_INFO_DICT,
    vrs_fields_dict: Dict[str, Dict[str, str]] = VRS_FIELDS_DICT,
    label_delimiter: str = "_",
    data_type: str = "exomes",
    freq_comparison_included: bool = False,
    extra_suffix: str = None,
    extra_description_text: str = None,
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
    :param subset_list: List of sample subsets in dataset. Default is SUBSETS["exomes"].
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
    :param freq_comparison_included: Whether frequency comparison data is included in the HT.
    :param extra_suffix: Suffix to add to INFO field.
    :param extra_description_text: Extra description text to add to INFO field.
    :return: Updated INFO dictionary for VCF export.
    """
    vcf_info_dict = {}
    if data_type == "joint" and "joint" in subset_list:
        vcf_info_dict.update(JOINT_REGION_FLAGS_INFO_DICT)
    else:
        # Get existing info fields from predefined info_dict, e.g. `FS`,
        # `non_par`, `negative_train_site`...
        vcf_info_dict.update(info_dict)
        vcf_info_dict = {f: vcf_info_dict[f] for f in info_fields if f in vcf_info_dict}
        # Add allele-specific fields to info dict, including AS_VQSR_FIELDS
        vcf_info_dict.update(
            add_as_info_dict(info_dict=info_dict, as_fields=AS_FIELDS + AS_VQSR_FIELDS)
        )

    for subset in subset_list:
        subset_pops = deepcopy(pops)
        if data_type == "joint":
            description_text = f" in {subset} dataset " if subset != "" else ""
        else:
            description_text = "" if subset == "" else f" in {subset} subset"

        vcf_info_dict.update(
            populate_subset_info_dict(
                subset=subset,
                description_text=description_text,
                data_type=data_type,
                pops=subset_pops,
                faf_pops=faf_pops,
                sexes=sexes,
                label_delimiter=label_delimiter,
                freq_comparison_included=freq_comparison_included,
            )
        )

    if age_hist_distribution:
        age_hist_distribution = "|".join(str(x) for x in age_hist_distribution)

    # Add age histogram data to info dict.
    vcf_info_dict.update(
        make_info_dict(
            suffix=extra_suffix,
            label_delimiter=label_delimiter,
            bin_edges=bin_edges,
            age_hist_distribution=age_hist_distribution,
            description_text=extra_description_text,
        )
    )

    # Add variant quality histograms to info dict.
    vcf_info_dict.update(
        make_hist_dict(
            bin_edges,
            adj=True,
            drop_n_smaller_larger=True,
            suffix=extra_suffix,
            description_text=extra_description_text,
        )
    )
    if data_type == "joint":
        vcf_info_dict["exomes_filters"] = {
            "Description": "Filters annotation for exomes data type."
        }
        vcf_info_dict["genomes_filters"] = {
            "Description": "Filters annotation for genomes data type."
        }
        return vcf_info_dict

    # Add in silico prediction annotations to info_dict.
    vcf_info_dict.update(in_silico_dict)

    # Add VRS annotations to info_dict.
    vcf_info_dict.update(vrs_fields_dict)

    return vcf_info_dict


def prepare_vcf_header_dict(
    ht: hl.Table,
    validated_ht: Optional[hl.Table],
    info_fields: List[str],
    bin_edges: Dict[str, str],
    age_hist_distribution: str,
    subset_list: List[str],
    pops: Dict[str, str],
    data_type: str = "exomes",
    joint_included: bool = False,
    freq_comparison_included: bool = False,
    extra_suffix: str = None,
    extra_description_text: str = None,
) -> Dict[str, Dict[str, str]]:
    """
    Prepare VCF header dictionary.

    :param ht: Input Table
    :param validated_ht: Validated HT with unfurled info fields.
    :param info_fields: List of info fields to add to the info dict.
    :param bin_edges: Dictionary of variant annotation histograms and their associated
        bin edges.
    :param age_hist_distribution: Pipe-delimited string of overall age histogram bin
        frequency.
    :param subset_list: List of sample subsets in dataset.
    :param pops: List of sample global genetic ancestry group names for gnomAD data type.
    :param data_type: Data type to prepare VCF header for. One of "exomes" or "genomes".
        Default is "exomes".
    :param joint_included: Whether joint frequency data is included in the HT. Default is False.
    :param freq_comparison_included: Whether frequency comparison data is included in the HT.
    :param extra_suffix: Suffix to add to INFO field.
    :param extra_description_text: Extra description text to add to INFO field.
    :return: Prepared VCF header dictionary.
    """
    if data_type != "joint":
        logger.info("Making FILTER dict for VCF...")
        filter_dict = make_vcf_filter_dict(
            hl.eval(ht.filtering_model.snv_cutoff.min_score),
            hl.eval(ht.filtering_model.indel_cutoff.min_score),
            inbreeding_cutoff=hl.eval(ht.inbreeding_coeff_cutoff),
            variant_qc_filter=hl.eval(ht.filtering_model.filter_name),
        )
        # subset = "" represents full dataset in VCF header construction, the
        # logic in gnomad_methods is built around this.
        subset_list.extend(["", "joint"] if joint_included else [""])

    logger.info("Making INFO dict for VCF...")
    vcf_info_dict = populate_info_dict(
        info_fields=info_fields,
        bin_edges=bin_edges,
        age_hist_distribution=age_hist_distribution,
        subset_list=subset_list,
        pops=pops,
        data_type=data_type,
        freq_comparison_included=freq_comparison_included,
        extra_suffix=extra_suffix,
        extra_description_text=extra_description_text,
    )

    if data_type != "joint":
        vcf_info_dict.update(
            {"vep": {"Description": hl.eval(validated_ht.vep_csq_header)}}
        )

    # Adjust keys to remove adj tags before exporting to VCF.
    new_vcf_info_dict = {i.replace("_adj", ""): j for i, j in vcf_info_dict.items()}

    if data_type == "joint":
        header_dict = new_vcf_info_dict
    else:
        header_dict = {
            "info": new_vcf_info_dict,
            "filter": filter_dict,
        }

    return header_dict


def get_downsamplings_fields(ht: hl.Table) -> List[str]:
    """
    Get downsampling specific annotations from info struct.

    .. note::

        This retrieves any annotation that contains any downsampling value in its name.

    :param ht: Input Table.
    :return: List of downsampling specific annotations to drop from info struct.
    """
    ds = hl.set(hl.flatten(ht.downsamplings.values()))
    ds = hl.eval(ds.map(lambda d: hl.str(d)))
    ds_fields = list(
        set([field for d in ds for field in list(ht.info) if str(d) in field])
    )

    return ds_fields


def format_validated_ht_for_export(
    ht: hl.Table,
    data_type: str = "exomes",
    vcf_info_reorder: Optional[List[str]] = VCF_INFO_REORDER,
    info_fields_to_drop: Optional[List[str]] = None,
) -> Tuple[hl.Table, List[str]]:
    """
    Format validated HT for export.

    Drop downsamplings frequency stats from info, rearrange info, and make sure fields
    are VCF compatible.

    :param ht: Validated HT.
    :param data_type: Data type to format validated HT for. One of "exomes" or "genomes".
        Default is "exomes".
    :param vcf_info_reorder: Order of VCF INFO fields. These will be placed in front of
        all other fields in the order specified.
    :param info_fields_to_drop: List of info fields to drop from the info struct.
    :return: Formatted HT and list rename row annotations.
    """
    if info_fields_to_drop is None:
        info_fields_to_drop = []

    if data_type == "exomes":
        logger.info("Getting downsampling annotations to drop from info struct...")
        ds_fields = get_downsamplings_fields(ht)
        info_fields_to_drop.extend(ds_fields)

    # Drop any info annotation with "hgdp" or "tgp" in the name for genomes.
    if data_type == "genomes":
        logger.info("Dropping hgdp and tgp annotations from info struct...")
        info_fields_to_drop.extend(
            [f for f in list(ht.info) if (("hgdp" in f) | ("tgp" in f))]
        )

    logger.info("Add age_histogram bin edges to info fields to drop...")
    info_fields_to_drop.extend(["age_hist_het_bin_edges", "age_hist_hom_bin_edges"])

    logger.info("Adding 'SB' to info fields to drop...")
    info_fields_to_drop.append("SB")

    logger.info(
        "Dropping the following fields from info struct: %s...",
        pprint(info_fields_to_drop),
    )
    ht = ht.annotate(info=ht.info.drop(*info_fields_to_drop))

    logger.info("Dropping _'adj' from info fields...")
    row_annots = list(ht.info)
    new_row_annots = [x.replace("_adj", "") for x in row_annots]
    info_annot_mapping = dict(
        zip(new_row_annots, [ht.info[f"{x}"] for x in row_annots])
    )
    ht = ht.transmute(info=hl.struct(**info_annot_mapping))

    logger.info("Adjusting VCF incompatible types...")
    # Reformat AS_SB_TABLE for use in adjust_vcf_incompatible_types
    ht = ht.annotate(
        info=ht.info.annotate(
            AS_SB_TABLE=hl.array([ht.info.AS_SB_TABLE[:2], ht.info.AS_SB_TABLE[2:]])
        )
    )
    # The Table is already split so there are no annotations that need to be
    # pipe delimited.
    ht = adjust_vcf_incompatible_types(ht, pipe_delimited_annotations=[])

    logger.info("Rearranging fields to desired order...")
    ht = ht.annotate(
        info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder))
    )
    return ht, new_row_annots


def process_vep_csq_header(vep_csq_header: str = VEP_CSQ_HEADER) -> str:
    """
    Process VEP CSQ header string, delimited by '|', to remove polyphen and sift annotations.

    :param vep_csq_header: VEP CSQ header.
    :return: Processed VEP CSQ header.
    """
    logger.info("Processing VEP CSQ header...")
    vep_csq_header = vep_csq_header.split("|")
    vep_csq_header = [f for f in vep_csq_header if f not in ["PolyPhen", "SIFT"]]
    vep_csq_header = "|".join(vep_csq_header)
    return vep_csq_header


def check_globals_for_retired_terms(ht: hl.Table) -> None:
    """
    Check list of dictionaries to see if the keys in the dictionaries contain either 'pop and 'oth'.

    :param ht: Input Table
    """
    logger.info("Checking globals for retired terms...")
    errors = []

    for field in ht.globals:
        if field.endswith("meta"):
            for d in hl.eval(ht[field]):
                if "pop" in d.keys():
                    errors.append(
                        f"Found retired term 'pop' in global {field} annotation: {d}"
                    )
                if "oth" in d.values():
                    errors.append(
                        f"Found retired term 'oth' in global {field} annotation: {d}"
                    )
        if "index_dict" in field:
            for k in hl.eval(ht[field]).keys():
                if "oth" in k:
                    errors.append(
                        f"Found retired term 'oth' in global {field} annotation: {k}"
                    )

    if len(errors) > 0:
        logger.info("Failed retired term check")
        pprint(errors)
    else:
        logger.info("Passed retired term check: No retired terms found in globals.")


def main(args):
    """Validate release Table and export VCFs."""
    hl.init(
        log="/validate_and_export_vcf.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with and this hits a
    # class too large error without this no whole stage codegen flag
    # after I added joint_fafmax -- just the straw that broke hails back
    hl._set_flags(use_ssa_logs="1", no_whole_stage_codegen="1")

    overwrite = args.overwrite
    test = args.test
    data_type = args.data_type
    contig = args.contig
    joint_included = args.joint_included
    for_joint = data_type == "joint"

    resources = get_export_resources(
        overwrite,
        data_type,
        test,
    )

    if contig and test:
        raise ValueError(
            "Test argument cannot be used with contig argument as test filters"
            " to chr20, X, and Y."
        )
    if data_type == "joint" and joint_included:
        raise ValueError(
            "Joint data type and joint included cannot both be true. joint_included"
            "is only used for exomes and genomes data types."
        )
    try:
        if args.validate_release_ht:
            logger.info("Running %s release HT validation...", data_type)
            res = resources.validate_release_ht
            res.check_resource_existence()
            ht = res.release_ht.ht()

            logger.info(
                "Checking globals for retired terms and checking their associated row"
                " annotation lengths..."
            )
            check_globals_for_retired_terms(ht)
            pprint_global_anns(ht)
            validate_hts = {}
            if data_type == "joint":
                iter_data_types = ["exomes", "genomes", "joint"]
            else:
                iter_data_types = [data_type]

            for dt in iter_data_types:
                if for_joint:
                    dt_ht = select_type_from_joint_ht(ht, dt)
                else:
                    dt_ht = ht

                len_comp_global_rows = (
                    LEN_COMP_GLOBAL_ROWS + LEN_COMP_JOINT_GLOBAL_ROWS
                    if joint_included
                    else LEN_COMP_GLOBAL_ROWS
                )
                check_global_and_row_annot_lengths(dt_ht, len_comp_global_rows)

                logger.info("Preparing %s HT for validity checks and export...", dt)
                dt_ht, rename_dict = prepare_ht_for_validation(
                    dt_ht,
                    data_type=dt,
                    joint_included=joint_included,
                    freq_comparison_included=(dt == "joint"),
                    for_joint_validation=for_joint,
                )
                dt_ht = dt_ht.checkpoint(hl.utils.new_temp_file(f"validate_{dt}", "ht"))
                site_gt_check_expr = None

                if data_type != "joint":
                    site_gt_check_expr = {
                        "monoallelic": dt_ht.info.monoallelic,
                        "only_het": dt_ht.info.only_het,
                    }

                validate_release_t(
                    dt_ht,
                    subsets=SUBSETS[data_type],
                    pops=POPS[dt],
                    site_gt_check_expr=site_gt_check_expr,
                    verbose=args.verbose,
                    delimiter="_",
                    sample_sum_sets_and_pops=SAMPLE_SUM_SETS_AND_POPS[data_type],
                    variant_filter_field="AS_VQSR",
                    problematic_regions=REGION_FLAG_FIELDS[data_type],
                    single_filter_count=True,
                    filters_check=False if dt == "joint" else True,
                )
                if for_joint:
                    ordered_rename_dict = {
                        key: rename_dict.get(key, key) for key in dt_ht.info.keys()
                    }
                    dt_ht = dt_ht.annotate(info=dt_ht.info.rename(ordered_rename_dict))
                    if dt != "joint":
                        dt_ht = dt_ht.select("info")

                    dt_ht = dt_ht.select_globals(
                        **{f"{dt}_{f}": dt_ht[f] for f in dt_ht.globals}
                    )

                validate_hts[dt] = dt_ht

            ht = validate_hts[data_type]
            if for_joint:
                in_joint_ht = set(ht.info.keys())
                for dt in ["exomes", "genomes"]:
                    info_expr = validate_hts[dt][ht.key].info
                    info_expr = info_expr.select(
                        *[f for f in info_expr if f not in in_joint_ht]
                    )
                    ht = ht.annotate(info=ht.info.annotate(**info_expr))
                    ht = ht.annotate_globals(**validate_hts[dt].index_globals())

            ht.describe()

            # Note: Checkpoint saves time in validity checks and final export by not
            # needing to run the VCF HT prep on each chromosome -- more needs to happen
            # before ready for export, but this is an intermediate write.
            logger.info("Writing prepared VCF HT for validity checks and export...")
            ht.write(res.validated_ht.path, overwrite=overwrite)

        if args.prepare_vcf_header:
            logger.info("Preparing VCF header dict...")
            res = resources.prepare_vcf_header
            res.check_resource_existence()
            ht = res.release_ht.ht()
            validated_ht = res.validated_ht.ht()

            if not for_joint:
                # v4 Genomes drops subsets from VCF
                subsets = SUBSETS["exomes"] if data_type == "exomes" else []
                header_dict = prepare_vcf_header_dict(
                    ht,
                    validated_ht=validated_ht,
                    info_fields=list(validated_ht.info),
                    bin_edges=make_hist_bin_edges_expr(
                        ht,
                        include_age_hists=True,
                    ),
                    age_hist_distribution=hl.eval(ht.age_distribution.bin_freq),
                    subset_list=subsets,
                    pops=POPS[data_type],
                    data_type=data_type,
                    joint_included=joint_included,
                )
            else:
                header_dict = {"info": {}}
                for dt in ["exomes", "genomes", "joint"]:
                    dt_ht = select_type_from_joint_ht(ht, dt)
                    temp_header_dict = prepare_vcf_header_dict(
                        dt_ht,
                        validated_ht=validated_ht,
                        info_fields=[f for f in validated_ht.info.keys() if dt in f],
                        bin_edges=make_hist_bin_edges_expr(
                            dt_ht, include_age_hists=True
                        ),
                        age_hist_distribution=hl.eval(dt_ht.age_distribution.bin_freq),
                        subset_list=[dt],
                        pops=POPS[dt],
                        data_type="joint",
                        joint_included=joint_included,
                        freq_comparison_included=(dt == "joint"),
                        extra_suffix=dt,
                        extra_description_text=f" in {dt} dataset",
                    )
                    header_dict["info"].update(temp_header_dict)

            logger.info("Writing VCF header dict...")
            with hl.hadoop_open(res.vcf_header_path, "wb") as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)

        # NOTE: The following step is not yet implemented for joint and should not
        # be reviewed for it
        if args.export_vcf:
            if data_type == "joint":
                raise ValueError("Joint data type is not yet supported for VCF export.")
            if contig and test:
                raise ValueError(
                    "Test argument cannot be used with contig argument as test filters"
                    " to chr20, X, and Y."
                )
            contig = f"chr{contig}" if contig else None
            logger.info(f"Exporting VCF{f' for {contig}' if contig else ''}...")
            res = resources.export_vcf
            res.check_resource_existence()
            ht = res.validated_ht.ht()
            logger.info("Reading release HT from %s...", res.validated_ht.path)
            with hl.hadoop_open(res.vcf_header_path, "rb") as f:
                header_dict = pickle.load(f)

            if test:
                logger.info("Filtering to test partitions on chr20, X, and Y...")
                ht = filter_to_test(ht)
            if contig:
                logger.info(f"Filtering to {contig}...")
                ht = hl.filter_intervals(
                    ht, [hl.parse_locus_interval(contig, reference_genome="GRCh38")]
                )

            ht, new_row_annots = format_validated_ht_for_export(ht, data_type=data_type)

            logger.info("Running check on VCF fields and info dict...")
            if not vcf_field_check(ht, header_dict, new_row_annots):
                raise ValueError("Did not pass VCF field check")

            output_path = (
                f"{qc_temp_prefix(data_type=data_type)}gnomad.{data_type}.test{'.chr'+contig if contig else ''}.vcf.bgz"
                if test
                else release_vcf_path(contig=contig, data_type=data_type)
            )
            # Match header order to info field order
            logger.info("Matching header order to info field order...")
            ordered_vcf_info_dict = {
                f: header_dict["info"][f]
                for f in list(ht.info)
                if f in header_dict["info"]
            }
            header_dict.update({"info": ordered_vcf_info_dict})

            logger.info("Exporting VCF...")
            export_reference = build_vcf_export_reference(
                "gnomAD_GRCh38", keep_chrM=False
            )
            hl.export_vcf(
                rekey_new_reference(ht, export_reference),
                output_path,
                metadata=header_dict,
                append_to_header=append_to_vcf_header_path(data_type=data_type),
                tabix=True,
            )

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("validity_checks_and_export"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validate-release-ht",
        help="Run release HT validation",
        action="store_true",
    )
    parser.add_argument(
        "--verbose",
        help="Log successes in addition to failures during validation",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data and start from raw inputs",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            "For validation, test will run on chr20, chrX, and chrY. For VCF export,"
            " test runs on PCSK9 region. Outputs to test bucket."
        ),
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--data-type",
        help="Data type to run validity checks on.",
        default="exomes",
        choices=["exomes", "genomes", "joint"],
    )
    parser.add_argument(
        "--prepare-vcf-header",
        help="Prepare VCF header dict.",
        action="store_true",
    )
    parser.add_argument(
        "--export-vcf",
        help="Export VCF.",
        action="store_true",
    )
    parser.add_argument(
        "--contig",
        help="Contig to export VCF for.",
        default=None,
    )
    parser.add_argument(
        "--joint-included",
        help="Whether joint frequency data is included in the release HT.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
