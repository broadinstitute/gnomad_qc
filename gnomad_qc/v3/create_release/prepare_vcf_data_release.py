import argparse
import logging
import pickle
from typing import Dict, List, Optional, Set, Union

import hail as hl

from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.resources.grch38.gnomad import (
    COHORTS_WITH_POP_STORED_AS_SUBPOP,
    HGDP_POPS,
    KG_POPS,
    KG_POP_NAMES,
    POPS,
    SEXES,
    SUBSETS,
)
from gnomad.resources.resource_utils import DataException
from gnomad.utils.file_utils import file_exists
from gnomad.utils.vep import VEP_CSQ_HEADER, vep_struct_to_csq
from gnomad.utils.vcf import (
    add_as_info_dict,
    adjust_vcf_incompatible_types,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    AS_VQSR_FIELDS,
    FAF_POPS,
    FORMAT_DICT,
    GROUPS,
    HISTS,
    INFO_DICT,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
    make_vcf_filter_dict,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
)
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF

from gnomad_qc.v3.create_release.sanity_checks import (
    sanity_check_release_ht,
    vcf_field_check,
)
from gnomad_qc.v3.resources.basics import get_checkpoint_path, qc_temp_prefix

# TODO: Uncomment when this resource goes in
# from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.release import (
    append_to_vcf_header_path,
    release_header_path,
    release_subset,
    release_vcf_path,
)
from gnomad_qc.v3.utils import (
    build_export_reference,
    rekey_new_reference,
    remove_fields_from_globals,
)

hl.stop()

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)

# Add monoallelic, QUALapprox, and AS_SB_TABLE to vcf_info_dict
VCF_INFO_DICT = INFO_DICT
VCF_INFO_DICT["monoallelic"] = {
    "Description": "All samples are all homozygous alternate for the variant"
}
VCF_INFO_DICT["QUALapprox"] = {
    "Number": "1",
    "Description": "Sum of PL[0] values; used to approximate the QUAL score",
}
VCF_INFO_DICT["AS_SB_TABLE"] = {
    "Number": ".",
    "Description": "Allele-specific forward/reverse read counts for strand bias tests",
}

# Add new site fields
NEW_SITE_FIELDS = [
    "monoallelic",
    "QUALapprox",
    "transmitted_singleton",
]
SITE_FIELDS.extend(NEW_SITE_FIELDS)
AS_FIELDS.append("AS_SB_TABLE")

# Remove original alleles for containing non-releasable alleles
MISSING_ALLELE_TYPE_FIELDS = ["original_alleles", "has_star"]
remove_fields_from_globals(ALLELE_TYPE_FIELDS, MISSING_ALLELE_TYPE_FIELDS)

# Remove BaseQRankSum and SOR from site fields (doesn't exist in v3.1)
MISSING_SITES_FIELDS = ["BaseQRankSum", "SOR"]
remove_fields_from_globals(SITE_FIELDS, MISSING_SITES_FIELDS)

# Remove AS_BaseQRankSum and AS_SOR from AS fields
MISSING_AS_FIELDS = ["AS_BaseQRankSum", "AS_VarDP"]
remove_fields_from_globals(AS_FIELDS, MISSING_AS_FIELDS)

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST_FOR_VCF = SUBSETS.copy()
SUBSET_LIST_FOR_VCF.append("")
remove_fields_from_globals(SUBSET_LIST_FOR_VCF, COHORTS_WITH_POP_STORED_AS_SUBPOP)

# Remove decoy from region field flag
MISSING_REGION_FIELDS = ["decoy"]
remove_fields_from_globals(REGION_FLAG_FIELDS, MISSING_REGION_FIELDS)

# All missing fields to remove from vcf info dict
MISSING_INFO_FIELDS = (
    MISSING_ALLELE_TYPE_FIELDS
    + MISSING_AS_FIELDS
    + MISSING_REGION_FIELDS
    + MISSING_SITES_FIELDS
    + RF_FIELDS
)

# Remove unnecessary pop names from POPS dict
POPS = {pop: POP_NAMES[pop] for pop in POPS}

# Remove unnecessary pop names from FAF_POPS dict
FAF_POPS = {pop: POP_NAMES[pop] for pop in FAF_POPS}

# Get HGDP + TGP(KG) subset pop names
HGDP_KG_KEEP_POPS = KG_POPS + HGDP_POPS
HGDP_KG_POPS = {}
for pop in HGDP_KG_KEEP_POPS:
    if pop in KG_POP_NAMES:
        HGDP_KG_POPS[pop] = KG_POP_NAMES[pop]
    else:
        HGDP_KG_POPS[pop] = pop.capitalize()


# Histograms to exclude from the VCF export
DROP_HISTS = (
    [x + "_n_smaller" for x in HISTS if "dp_hist_all" not in x]
    + [x + "_n_larger" for x in HISTS if "dp_" not in x]
    + [x + "_raw_n_smaller" for x in HISTS]
    + [x + "_raw_bin_edges" for x in HISTS]
    + [x + "_raw_n_larger" for x in HISTS]
    + [x + "_raw_bin_freq" for x in HISTS]
    + [x + "_bin_edges" for x in HISTS]
)

# Dictionary with in silico score descriptions to include in the VCF INFO header
IN_SILICO_ANNOTATIONS_INFO_DICT = {
    "cadd_raw_score": {
        "Number": "1",
        "Description": "Raw CADD scores are interpretable as the extent to which the annotation profile for a given variant suggests that the variant is likely to be 'observed' (negative values) vs 'simulated' (positive values). Larger values are more deleterious.",
    },
    "cadd_phred": {
        "Number": "1",
        "Description": "Cadd Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on the rank of each variant relative to all possible 8.6 billion substitutions in the human reference genome. Larger values are more deleterious",
    },
    "revel_score": {
        "Number": "1",
        "Description": "dbNSFP's Revel score from 0 to 1. Variants with higher scores are predicted to be more likely to be deleterious.",
    },
    "splice_ai_max_ds": {
        "Number": "1",
        "Description": "Illumina's SpliceAI max delta score; interpreted as the probability of the variant being splice-altering.",
    },
    "splice_ai_consequence": {
        "Description": "The consequence term associated with the max delta score in 'splice_ai_max_ds'.",
    },
    "primate_ai_score": {
        "Number": "1",
        "Description": "PrimateAI's deleteriousness score from 0 (less deleterious) to 1 (more deleterious).",
    },
}

# Used for HGDP + TGP subset MT VCF output only
FORMAT_DICT.update(
    {
        "RGQ": {
            "Number": "1",
            "Type": "Integer",
            "Description": "Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)",
        }
    }
)

# VCF INFO field reordering
VCF_INFO_REORDER = ["AC-adj", "AN-adj", "AF-adj", "popmax", "faf95_popmax"]
HGDP_KG_VCF_INFO_REORDER = [
    "AC-adj",
    "AN-adj",
    "AF-adj",
    "AC-raw",
    "AN-raw",
    "AF-raw",
    "gnomad-AC-adj",
    "gnomad-AN-adj",
    "gnomad-AF-adj",
    "gnomad-popmax",
    "gnomad-faf95_popmax",
    "gnomad-AC-raw",
    "gnomad-AN-raw",
    "gnomad-AF-raw",
]


# TODO: USE RESOURCES once it is in, looks like it will be called release_sites so just remove this
def release_ht_path():
    return "gs://gnomad/release/3.1.1/ht/genomes/gnomad.genomes.v3.1.1.sites.ht"


def populate_subset_info_dict(
    subset: str,
    description_text: str,
    groups: List[str] = GROUPS,
    pops: List[str] = POPS,
    faf_pops: List[str] = FAF_POPS,
    sexes: List[str] = SEXES,
    label_delimiter: str = "_",
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` to populate INFO dictionary with specific sexes, population names, and filtering allele
    frequency (faf) pops for the requested subset.

    :param subset: Sample subset in dataset.
    :param description_text: Text describing the sample subset that should be added to the INFO description.
    :param groups: List of sample groups [adj, raw]. Default is GROUPS.
    :param pops: List of sample global population names for gnomAD genomes. Default is POPS.
    :param faf_pops: List of faf population names. Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :return: Dictionary containing Subset specific INFO header fields.
    """
    def _create_label_groups(
        pops: Union[Dict[str, str], List[str]],
        sexes: List[str],
        group: List[str] = ["adj"],
    ) -> List[Dict[str, List[str]]]:
        """
        Generates list of label group dictionaries needed to populate info dictionary.

        Label dictionaries are passed as input to `make_info_dict`.

        :param pops: Dict or list of population names.
        :param sexes: List of sample sexes.
        :param group: List of data types (adj, raw). Default is ["adj"].
        :return: List of label group dictionaries.
        """
        return [
            dict(group=groups),  # this is to capture raw fields
            dict(group=group, sex=sexes),
            dict(group=group, pop=pops),
            dict(group=group, pop=pops, sex=sexes),
        ]

    vcf_info_dict = {}
    faf_label_groups = _create_label_groups(pops=faf_pops, sexes=sexes)
    for label_group in faf_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                pop_names=pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                faf=True,
                description_text=description_text,
            )
        )

    label_groups = _create_label_groups(pops=pops, sexes=sexes)
    for label_group in label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                prefix_before_metric=False,
                pop_names=pops,
                label_groups=label_group,
                label_delimiter=label_delimiter,
                description_text=description_text,
            )
        )

    # Add popmax to info dict
    vcf_info_dict.update(
        make_info_dict(
            prefix=subset,
            prefix_before_metric=False,
            label_delimiter=label_delimiter,
            pop_names=pops,
            popmax=True,
            description_text=description_text,
        )
    )

    return vcf_info_dict


def populate_info_dict(
    bin_edges: Dict[str, str],
    age_hist_data: str = None,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    subset_list: List[str] = SUBSETS,
    groups: List[str] = GROUPS,
    pops: List[str] = POPS,
    gnomad_pops: List[str] = POPS,
    faf_pops: Dict[str, str] = FAF_POPS,
    sexes: List[str] = SEXES,
    in_silico_dict: Dict[str, Dict[str, str]] = IN_SILICO_ANNOTATIONS_INFO_DICT,
    label_delimiter="_",
) -> Dict[str, Dict[str, str]]:
    """
    Call `make_info_dict` and `make_hist_dict` to populate INFO dictionary with specific sexes, population names, and filtering allele frequency (faf) pops.

    Used during VCF export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations
        - INFO fields for variant histograms (hist_bin_freq, hist_n_smaller, hist_n_larger for each histogram)

    :param bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param info_dict: INFO dict to be populated.
    :param subset_list: List of sample subsets in dataset. Default is SUBSETS.
    :param groups: List of sample groups [adj, raw]. Default is GROUPS.
    :param pops: List of sample global population names for gnomAD genomes. Default is POPS.
    :param gnomad_pops: List of sample global population names for gnomAD genomes. Default is POPS.
    :param faf_pops: List of faf population names. Default is FAF_POPS.
    :param sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :param in_silico_dict: Dictionary of in silico predictor score descriptions.
    :return: Updated INFO dictionary for VCF export.
    """
    vcf_info_dict = info_dict.copy()

    # Remove MISSING_INFO_FIELDS from info dict
    for field in MISSING_INFO_FIELDS:
        vcf_info_dict.pop(field, None)

    # Add allele-specific fields to info dict, including AS_VQSR_FIELDS
    vcf_info_dict.update(
        add_as_info_dict(info_dict=info_dict, as_fields=AS_FIELDS + AS_VQSR_FIELDS)
    )

    for subset in subset_list:
        if "gnomad" in subset:
            description_text = " in gnomAD"
            pops = gnomad_pops
        else:
            description_text = "" if subset == "" else f" in {subset} subset"
            pops = pops

        vcf_info_dict.update(
            populate_subset_info_dict(
                subset=subset,
                description_text=description_text,
                groups=groups,
                pops=pops,
                faf_pops=faf_pops,
                sexes=sexes,
                label_delimiter=label_delimiter,
            )
        )

    if age_hist_data:
        age_hist_data = "|".join(str(x) for x in age_hist_data)

    vcf_info_dict.update(
        make_info_dict(
            prefix="",
            prefix_before_metric=False,
            label_delimiter=label_delimiter,
            bin_edges=bin_edges,
            popmax=True,
            age_hist_data=age_hist_data,
        )
    )

    # Add variant quality histograms to info dict
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=True))

    # Add Analyst annotations to info_dict
    vcf_info_dict.update(in_silico_dict)

    return vcf_info_dict


def make_info_expr(t: Union[hl.MatrixTable, hl.Table]) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expression for variant annotations to be included in VCF INFO field.

    :param Table/MatrixTable t: Table/MatrixTable containing variant annotations to be reformatted for VCF export.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    :rtype: Dict[str, hl.expr.Expression]
    """
    vcf_info_dict = {}
    # Add site-level annotations to vcf_info_dict
    for field in SITE_FIELDS:
        vcf_info_dict[field] = t["release_ht_info"][f"{field}"]

    # Add AS annotations to info dict
    for field in AS_FIELDS:
        vcf_info_dict[field] = t["release_ht_info"][f"{field}"]
    for field in VQSR_FIELDS + AS_VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr"][f"{field}"]

    # Add region_flag and allele_info fields to info dict
    for field in ALLELE_TYPE_FIELDS:
        vcf_info_dict[field] = t["allele_info"][f"{field}"]
    for field in REGION_FLAG_FIELDS:
        vcf_info_dict[field] = t["region_flag"][f"{field}"]

    # Add histograms to info dict
    for hist in HISTS:
        for prefix in ["qual_hists", "raw_qual_hists"]:
            hist_name = hist
            if "raw" in prefix:
                hist_name = f"{hist}_raw"

            hist_dict = {
                f"{hist_name}_bin_freq": hl.delimit(
                    t[prefix][hist].bin_freq, delimiter="|"
                ),
                f"{hist_name}_bin_edges": hl.delimit(
                    t[prefix][hist].bin_edges, delimiter="|"
                ),
                f"{hist_name}_n_smaller": t[prefix][hist].n_smaller,
                f"{hist_name}_n_larger": t[prefix][hist].n_larger,
            }
            vcf_info_dict.update(hist_dict)

    # Add in silico annotations to info dict
    vcf_info_dict["cadd_raw_score"] = t["cadd"]["raw_score"]
    vcf_info_dict["cadd_phred"] = t["cadd"]["phred"]

    vcf_info_dict["revel_score"] = t["revel"]["revel_score"]

    vcf_info_dict["splice_ai_max_ds"] = t["splice_ai"]["splice_ai_score"]
    vcf_info_dict["splice_ai_consequence"] = t["splice_ai"]["splice_consequence"]

    vcf_info_dict["primate_ai_score"] = t["primate_ai"]["primate_ai_score"]

    return vcf_info_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table],
    is_subset: bool = False,
    add_gnomad_release: bool = False,
    entries_to_remove: Set[str] = None,
) -> [hl.struct, List]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.

    :param t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param add_gnomad_release: Should the gnomAD release frequencies be unfurled.
    :param is_subset: Is this for the release of a subset.
    :param entries_to_remove: Frequency entries to remove for vcf_export.
    :return: Dictionary containing variant annotations and their corresponding values.
    """
    expr_dict = dict()

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict
    if is_subset & add_gnomad_release:
        gnomad_prefix = f"gnomad"
        popmax = f"{gnomad_prefix}_popmax"
        faf = f"{gnomad_prefix}_faf"
        freq = f"{gnomad_prefix}_freq"
        faf_idx = hl.eval(t.globals[f"{gnomad_prefix}_faf_index_dict"])
        freq_idx = hl.eval(t.globals[f"{gnomad_prefix}_freq_index_dict"])
    elif is_subset:
        freq = "cohort_freq"
        freq_idx = hl.eval(t.globals[f"cohort_freq_index_dict"])
    else:
        popmax = "popmax"
        freq = "freq"
        freq_idx = hl.eval(t.freq_index_dict)
        faf = "faf"
        faf_idx = hl.eval(t.faf_index_dict)

    # freqs to remove for vcf_export
    freq_entries_to_remove_vcf = []
    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    logger.info("Unfurling freq data...")
    for k, i in freq_idx.items():
        prefix = ""
        if add_gnomad_release:
            prefix = f"{gnomad_prefix}-"

        # Set combination to key
        # e.g., set entry of 'afr_adj' to combo
        combo = k
        combo_dict = {
            f"{prefix}AC-{combo}": t[freq][i].AC,
            f"{prefix}AN-{combo}": t[freq][i].AN,
            f"{prefix}AF-{combo}": t[freq][i].AF,
            f"{prefix}nhomalt-{combo}": t[freq][i].homozygote_count,
        }
        expr_dict.update(combo_dict)

        if k.split("-")[0] in entries_to_remove:
            freq_entries_to_remove_vcf.extend(combo_dict.keys())

    # Add popmax
    if is_subset & add_gnomad_release:
        prefix = f"{gnomad_prefix}-"
    else:
        prefix = ""

    combo_dict = {
        f"{prefix}popmax": t[popmax].pop,
        f"{prefix}AC_popmax": t[popmax].AC,
        f"{prefix}AN_popmax": t[popmax].AN,
        f"{prefix}AF_popmax": t[popmax].AF,
        f"{prefix}nhomalt_popmax": t[popmax].homozygote_count,
        f"{prefix}faf95_popmax": t[popmax].faf95,
    }
    expr_dict.update(combo_dict)

    ## Unfurl FAF index dict
    logger.info("Unfurling faf data...")
    for (
        k,
        i,
    ) in faf_idx.items():  # NOTE: faf annotations are all done on adj-only groupings
        combo_dict = {
            f"{prefix}faf95-{k}": t[faf][i].faf95,
            f"{prefix}faf99-{k}": t[faf][i].faf99,
        }
        expr_dict.update(combo_dict)

    # Unfurl ages
    age_hist_dict = {
        "age_hist_het_bin_freq": hl.delimit(t.age_hist_het.bin_freq, delimiter="|"),
        "age_hist_het_bin_edges": hl.delimit(t.age_hist_het.bin_edges, delimiter="|"),
        "age_hist_het_n_smaller": t.age_hist_het.n_smaller,
        "age_hist_het_n_larger": t.age_hist_het.n_larger,
        "age_hist_hom_bin_freq": hl.delimit(t.age_hist_hom.bin_freq, delimiter="|"),
        "age_hist_hom_bin_edges": hl.delimit(t.age_hist_hom.bin_edges, delimiter="|"),
        "age_hist_hom_n_smaller": t.age_hist_hom.n_smaller,
        "age_hist_hom_n_larger": t.age_hist_hom.n_larger,
    }
    expr_dict.update(age_hist_dict)

    return hl.struct(**expr_dict), freq_entries_to_remove_vcf


def filter_to_test(
    t: Union[hl.Table, hl.MatrixTable], num_partitions: int = 1
) -> Union[hl.Table, hl.MatrixTable]:
    """
    Filter Table/MatrixTable to `num_partitions` partitions on chr20, chrX, and chrY for testing.

    :param t: Input Table/MatrixTable to filter.
    :param num_partitions: Number of partitions to grab from each chromosome.
    :return: Input Table/MatrixTable filtered to `num_partitions` on chr20, chrX, and chrY.
    """
    logger.info(
        f"Filtering to {num_partitions} partitions on chr20, chrX, and chrY (for tests only)..."
    )
    t_chr20 = hl.filter_intervals(t, [hl.parse_locus_interval("chr20")])
    t_chr20 = t_chr20._filter_partitions(range(num_partitions))

    t_chrx = hl.filter_intervals(t, [hl.parse_locus_interval("chrX")])
    t_chrx = t_chrx._filter_partitions(range(num_partitions))

    t_chry = hl.filter_intervals(t, [hl.parse_locus_interval("chrY")])
    t_chry = t_chry._filter_partitions(range(num_partitions))

    t = t_chr20.union(t_chrx, t_chry)

    return t


def prepare_vcf_ht(
    ht: hl.Table,
    is_subset: bool,
    add_gnomad_release: bool,
    freq_entries_to_remove: List[str],
    field_reorder: Optional[List[str]] = None,
) -> hl.Table:
    """
    Prepare the Table used for sanity checks and VCF export

    :param ht: Table containing the
    nested variant annotation arrays to be unfurled.
    :param is_subset: Is this for the release of a subset.
    :param add_gnomad_release: Should the gnomAD release frequencies be unfurled.
    :param freq_entries_to_remove: Frequency entries to remove for vcf_export.
    :param field_reorder: Optional list of INFO fields to reorder, the rest of the fields are added after this list.
    :return: Prepared HT for sanity checks and VCF export
    """
    logger.info("Starting preparation of VCF HT...")
    logger.info("Adding non-PAR annotation...")
    region_flag_expr = ht.region_flag.annotate(
        nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar())
    )

    logger.info(
        "Unfurling nested gnomAD frequency annotations and add to INFO field..."
    )
    info_struct, freq_entries_to_remove = unfurl_nested_annotations(
        ht,
        entries_to_remove=freq_entries_to_remove,
        is_subset=is_subset,
        add_gnomad_release=add_gnomad_release,
    )

    # NOTE: Merging rsid set into a semi-colon delimited string
    # dbsnp might have multiple identifiers for one variant
    # thus, rsid is a set annotation, starting with version b154 for dbsnp resource:
    # https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/resources/grch38/reference_data.py#L136
    # `export_vcf` expects this field to be a string, and vcf specs
    # say this field may be delimited by a semi-colon:
    # https://samtools.github.io/hts-specs/VCFv4.2.pdf
    logger.info("Reformatting rsid...")
    rsid_expr = hl.str(";").join(ht.rsid)

    logger.info("Reformatting VEP annotation...")
    vep_expr = vep_struct_to_csq(ht.vep)

    # Add variant annotations to INFO field
    # This adds annotations from:
    #   RF struct, VQSR struct, allele_info struct,
    #   info struct (site and allele-specific annotations),
    #   region_flag struct, and
    #   raw_qual_hists/qual_hists structs.
    logger.info("Constructing INFO field")
    ht = ht.annotate(
        region_flag=region_flag_expr,
        release_ht_info=ht.info,
        info=info_struct,
        rsid=rsid_expr,
        vep=vep_expr,
    )
    ht = ht.annotate(info=ht.info.annotate(**make_info_expr(ht), vep=ht.vep))
    ht = ht.annotate_globals(
        vep_csq_header=VEP_CSQ_HEADER, freq_entries_to_remove=freq_entries_to_remove,
    )

    # Select relevant fields for VCF export
    ht = ht.select("info", "filters", "rsid")

    logger.info("Rearranging fields to desired order...")
    ht = ht.annotate(info=ht.info.select(*field_reorder, *ht.info.drop(*field_reorder)))

    return ht


def prepare_vcf_header_dict(
    t: Union[hl.Table, hl.MatrixTable],
    bin_edges: Dict[str, str],
    age_hist_data: str,
    subset_list: List[str],
    pops: List[str],
) -> Dict[str, Dict[str, str]]:
    """
    Prepare VCF header dictionary.

    :param t: Input MatrixTable/Table
    :param bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param subset_list: List of sample subsets in dataset.
    :param pops: List of sample global population names for gnomAD genomes.
    :return: Prepared VCF header dictionary
    """
    logger.info("Making FILTER dict for VCF...")
    filter_dict = make_vcf_filter_dict(
        hl.eval(t.filtering_model.snv_cutoff.min_score),
        hl.eval(t.filtering_model.indel_cutoff.min_score),
        inbreeding_cutoff=INBREEDING_COEFF_HARD_CUTOFF,
        variant_qc_filter="AS_VQSR",
    )

    logger.info("Making INFO dict for VCF...")
    vcf_info_dict = populate_info_dict(
        bin_edges=bin_edges,
        age_hist_data=age_hist_data,
        subset_list=subset_list,
        pops=pops,
    )
    vcf_info_dict.update({"vep": {"Description": hl.eval(t.vep_csq_header)}})

    # Adjust keys to remove adj tags before exporting to VCF
    new_vcf_info_dict = {}
    for i, j in vcf_info_dict.items():
        i = i.replace("_adj", "")
        i = i.replace("-", "_")  # VCF 4.3 specs do not allow hyphens in info fields
        new_vcf_info_dict[i] = j

    header_dict = {
        "info": new_vcf_info_dict,
        "filter": filter_dict,
    }

    return header_dict


def main(args):

    hl.init(
        log="/vcf_release.log",
        default_reference="GRCh38",
        tmp_dir="hdfs:///vcf_write.tmp/",
    )
    chromosome = args.export_chromosome
    export_reference = build_export_reference()

    if chromosome and args.test:
        raise ValueError("chromosome argument doesn't work with the test flag.")

    # Setup of parameters and Table/MatrixTable based on hgdp_1kg_subset flag
    if args.hgdp_1kg_subset:
        pops = HGDP_KG_POPS
        subsets = ["", "gnomad"]
        is_subset = True
        add_gnomad_release = True
        temp_ht_path = get_checkpoint_path(
            f"vcf_prep_{'test' if args.test else ''}_hgdp_1kg"
        )

        logger.info("Reading in HGDP + TGP subset release MT, keeping only rows...")
        ht = release_subset(subset="hgdp_1kg", dense=True).mt().rows()
        age_hist_data = None
        field_reorder = [
            "AC-adj",
            "AN-adj",
            "AF-adj",
            "AC-raw",
            "AN-raw",
            "AF-raw",
            "gnomad-AC-adj",
            "gnomad-AN-adj",
            "gnomad-AF-adj",
            "gnomad-popmax",
            "gnomad-faf95_popmax",
            "gnomad-AC-raw",
            "gnomad-AN-raw",
            "gnomad-AF-raw",
        ]
    else:
        pops = POPS
        subsets = SUBSET_LIST_FOR_VCF
        is_subset = False
        add_gnomad_release = False
        temp_ht_path = get_checkpoint_path(f"vcf_prep_{'test' if args.test else ''}")

        logger.info("Reading in release HT...")
        ht = hl.read_table(release_ht_path())  # TODO: Change to release_sites().ht()
        logger.info("Getting age hist data...")
        age_hist_data = hl.eval(ht.age_distribution)
        field_reorder = ["AC-adj", "AN-adj", "AF-adj", "popmax", "faf95_popmax"]

    # Downsampling and subset entries to remove from VCF's freq export
    # Note: Need to extract the non-standard downsamplings from the freq_meta struct to the FREQ_ENTRIES_TO_REMOVE
    freq_entries_to_remove = {
        str(x["downsampling"]) for x in hl.eval(ht.freq_meta) if "downsampling" in x
    }
    freq_entries_to_remove.update(set(COHORTS_WITH_POP_STORED_AS_SUBPOP))

    try:
        if args.test:
            ht = filter_to_test(ht)

        if args.prepare_vcf_ht:
            ht = prepare_vcf_ht(
                ht, is_subset, add_gnomad_release, freq_entries_to_remove, field_reorder
            )

            # Note: Checkpoint saves time for the final export by not needing to run the VCF HT prep on each chromosome
            logger.info("Checkpointing prepared VCF HT for sanity checks and export...")
            ht = ht.checkpoint(temp_ht_path, overwrite=True)
            ht.describe()

        if args.prepare_vcf_header_dict or args.sanity_check or args.export_vcf:
            if not file_exists(temp_ht_path):
                raise DataException(
                    "The intermediate HT output doesn't exist, 'prepare_vcf_ht' needs to be run to create this file"
                )
            prepared_vcf_ht = hl.read_table(temp_ht_path)

        if args.prepare_vcf_header_dict:
            logger.info("Making histogram bin edges...")
            bin_edges = make_hist_bin_edges_expr(ht, prefix="")

            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(
                release_header_path(
                    subset="hgdp_1kg" if args.hgdp_1kg_subset else None
                ),
                "wb",
            ) as p:
                pickle.dump(
                    prepare_vcf_header_dict(
                        prepared_vcf_ht, bin_edges, age_hist_data, subsets, pops
                    ),
                    p,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )

        if args.sanity_check:
            sanity_check_release_ht(
                prepared_vcf_ht,
                SUBSETS,
                missingness_threshold=0.5,
                verbose=args.verbose,
            )

        if args.export_vcf:
            with hl.hadoop_open(
                release_header_path(
                    subset="hgdp_1kg" if args.hgdp_1kg_subset else None
                ),
                "rb",
            ) as f:
                header_dict = pickle.load(f)
            logger.info(
                "Dropping histograms and frequency entries that are not needed in VCF..."
            )
            prepared_vcf_ht = prepared_vcf_ht.annotate(
                info=prepared_vcf_ht.info.drop(
                    *DROP_HISTS, *hl.eval(prepared_vcf_ht.freq_entries_to_remove)
                )
            )

            # Reformat names to remove "adj" pre-export
            # e.g, renaming "AC-adj" to "AC"
            # All unlabeled frequency information is assumed to be adj
            logger.info("Dropping 'adj' from info annotations...")
            row_annots = list(prepared_vcf_ht.info)
            new_row_annots = []
            for x in row_annots:
                x = x.replace("-adj", "")
                x = x.replace(
                    "-", "_"
                )  # VCF 4.3 specs do not allow hyphens in info fields
                new_row_annots.append(x)

            info_annot_mapping = dict(
                zip(new_row_annots, [prepared_vcf_ht.info[f"{x}"] for x in row_annots])
            )
            prepared_vcf_ht = prepared_vcf_ht.transmute(
                info=hl.struct(**info_annot_mapping)
            )

            logger.info("Adjusting VCF incompatiable types...")
            # Reformat AS_SB_TABLE for use in ht_to_vcf_mt
            prepared_vcf_ht = prepared_vcf_ht.annotate(
                info=prepared_vcf_ht.info.annotate(
                    AS_SB_TABLE=hl.array(
                        [
                            prepared_vcf_ht.info.AS_SB_TABLE[:2],
                            prepared_vcf_ht.info.AS_SB_TABLE[2:],
                        ]
                    )
                )
            )
            prepared_vcf_ht = adjust_vcf_incompatible_types(prepared_vcf_ht)

            logger.info("Rearranging fields to desired order...")
            if args.hgdp_1kg_subset:
                t = (
                    release_subset(subset="hgdp_1kg", dense=True)
                    .mt()
                    .select_rows()
                    .annotate_rows(**prepared_vcf_ht)
                )
            else:
                t = prepared_vcf_ht

            if args.test:
                t = filter_to_test(t)

            if chromosome:
                t = hl.filter_intervals(t, [hl.parse_locus_interval(chromosome)])

            logger.info("Running check on VCF fields and info dict...")
            if not vcf_field_check(prepared_vcf_ht, header_dict, new_row_annots):
                logger.error("Did not pass VCF field check")

            logger.info(f"Export chromosome {chromosome}....")
            if args.test:
                output_path = f"{qc_temp_prefix}/gnomad.genomes_vcf_test{'hgdp_1kg_subset' if args.hgdp_1kg_subset else None}.vcf.bgz"
            else:
                output_path = release_vcf_path(
                    contig=chromosome,
                    subset="hgdp_1kg_subset" if args.hgdp_1kg_subset else None,
                )

            hl.export_vcf(
                rekey_new_reference(t, export_reference),
                output_path,
                append_to_header=append_to_vcf_header_path(),
                metadata=header_dict,
                tabix=True,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(f"{qc_temp_prefix}/logs/vcf_export.log")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--hgdp_1kg_subset",
        help="Use HGDP + TGP subset Matrix table instead of the gnomAD release HT",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Create release files using only 5 partitions on chr20, chrX, and chrY for testing purposes",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_vcf_ht",
        help="Use release Table or MatrixTable to create vcf HT",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_vcf_header_dict",
        help="Prepare the VCF header dictionary.",
        action="store_true",
    )
    parser.add_argument(
        "--sanity_check", help="Run sanity checks function", action="store_true"
    )
    parser.add_argument("--export_vcf", help="Export VCF", action="store_true")
    parser.add_argument("--export_chromosome", help="Which chromosome to export as VCF")
    parser.add_argument(
        "--verbose",
        help="Run sanity checks function with verbose output",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)
