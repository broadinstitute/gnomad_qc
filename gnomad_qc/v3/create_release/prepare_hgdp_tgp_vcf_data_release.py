import argparse
import json
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import vep_struct_to_csq, VEP_CSQ_HEADER
from gnomad.utils.vcf import (
    add_as_info_dict,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    RF_FIELDS,
    ENTRIES,
    FAF_POPS,
    FORMAT_DICT,
    GROUPS,
    HISTS,
    ht_to_vcf_mt,
    INFO_DICT,
    make_hist_dict,
    make_info_dict,
    make_label_combos,
    make_vcf_filter_dict,
    REGION_FLAG_FIELDS,
    SITE_FIELDS,
    SPARSE_ENTRIES,
    set_female_y_metrics_to_na,
    VQSR_FIELDS,
    remove_fields_from_globals,
    make_combo_header_text,
)

from gnomad_qc.v3.create_release.sanity_checks_hgdp_tgp import (
    sanity_check_release_mt,
    vcf_field_check,
)

from gnomad_qc.v3.resources.annotations import get_info

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)

# Add capture region and sibling singletons to vcf_info_dict
VCF_INFO_DICT = INFO_DICT

# Remove original alleles for containing non-releasable alleles
MISSING_ALLELE_TYPE_FIELDS = ["original_alleles", "has_star"]
remove_fields_from_globals(ALLELE_TYPE_FIELDS, MISSING_ALLELE_TYPE_FIELDS)

# Remove decoy from region field flag
MISSING_REGION_FIELDS = ["decoy"]
remove_fields_from_globals(REGION_FLAG_FIELDS, MISSING_REGION_FIELDS)

# Remove BaseQRankSum and SOR from site fields (doesn't exist in v3.1)
MISSING_SITES_FIELDS = ["BaseQRankSum", "SOR"]
remove_fields_from_globals(SITE_FIELDS, MISSING_SITES_FIELDS)

# Remove AS_BaseQRankSum and AS_SOR from AS fields
MISSING_AS_FIELDS = ["AS_BaseQRankSum", "AS_VarDP"]
remove_fields_from_globals(AS_FIELDS, MISSING_AS_FIELDS)

# All missing fields to remove from vcf info dict
MISSING_INFO_FIELDS = (
    MISSING_ALLELE_TYPE_FIELDS
    + MISSING_AS_FIELDS
    + MISSING_REGION_FIELDS
    + MISSING_SITES_FIELDS
    + RF_FIELDS
)

SEXES = ["XX", "XY"]

FORMAT_DICT.update({'RGQ':{"Number":"1","Type":"Integer","Description":"Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)"}})

TGP_POPS= {'CHB': 'Han Chinese','JPT': 'Japanese','CHS': 'Southern Han Chinese','CDX': 'Chinese Dai','KHV': 'Kinh','CEU': 'Utah Residents (European Ancestry)','TSI': 'Toscani','FIN': 'Finnish','GBR': 'British','IBS': 'Iberian','YRI': 'Yoruba','LWK': 'Luhya','GWD': 'Gambian','MSL': 'Mende','ESN': 'Esan','ASW': 'African-American','ACB': 'African Caribbean','MXL': 'Mexican-American','PUR': 'Puerto Rican','CLM': 'Colombian','PEL': 'Peruvian','GIH': 'Gujarati','PJL': 'Punjabi',
'BEB': 'Bengali','STU': 'Sri Lankan Tamil','ITU': 'Indian Telugu'}


# TODO: WHAT TO DO WITH POPS?????? Can we just use what is in subpop, or should we be using the longer descriptions?

# Select populations to keep from the list of population names in POP_NAMES
# This removes pop names we don't want to use
# (e.g., "uniform", "consanguineous") to reduce clutter
KEEP_POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]
GNOMAD_POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}
GNOMAD_POPS["mid"] = "Middle Eastern"


KG_POPS = ['esn', 'pur', 'pjl', 'clm', 'jpt', 'chb', 'stu', 'itu', 'tsi', 'mxl', 'ceu', 'msl', 'yri', 'beb', 'fin', 'khv', 'cdx', 'lwk', 'acb', 'asw', 'ibs', 'gbr', 'pel', 'gih', 'chs', 'gwd']
HGDP_POPS = ['japanese', 'papuan', 'adygei', 'orcadian', 'biakapygmy', 'yakut', 'han', 'uygur', 'miaozu', 'mongola', 'balochi', 'bedouin', 'russian', 'daur', 'pima', 'hezhen', 'sindhi', 'yizu', 'oroqen', 'san', 'tuscan', 'tu', 'palestinian', 'tujia', 'druze', 'pathan', 'basque', 'makrani', 'italian', 'naxi', 'karitiana', 'sardinian', 'mbutipygmy', 'mozabite', 'yoruba', 'lahu', 'dai', 'cambodian', 'melanesian', 'french', 'brahui', 'hazara', 'bantusafrica', 'surui', 'mandenka', 'kalash', 'xibo', 'colombian', 'bantukenya', 'she', 'burusho', 'maya']
SUBSET_KEEP_POPS = KG_POPS + HGDP_POPS

SUBSET_POPS = {}
for pop in SUBSET_KEEP_POPS:
    if pop.upper() in TGP_POPS:
        SUBSET_POPS[pop] = TGP_POPS[pop.upper()]
    else:
        SUBSET_POPS[pop] = pop.capitalize()

EXPORT_HISTS = [
    "gq_hist_alt_bin_freq",
    "gq_hist_all_bin_freq",
    "dp_hist_alt_bin_freq",
    "dp_hist_alt_n_larger",
    "dp_hist_all_bin_freq",
    "dp_hist_all_n_larger",
    "ab_hist_alt_bin_freq",
]

SORT_ORDER = [
    "popmax",
    "pop",
    "subpop",
    "sex",
    "group",
]
"""
Order to sort subgroupings during VCF export.
Ensures that INFO labels in VCF are in desired order (e.g., raw_AC_afr_female).
"""

INBREEDING_CUTOFF = -0.3

ANALYST_ANOTATIONS_INFO_DICT = {
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


def make_hist_bin_edges_expr(
    ht: hl.Table, hists: List[str] = HISTS, prefix: str = ""
) -> Dict[str, str]:
    """
    Create dictionaries containing variant histogram annotations and their associated bin edges, formatted into a string
    separated by pipe delimiters.

    :param ht: Table containing histogram variant annotations.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :param prefix: Prefix text for age histogram bin edges.  Default is empty string.
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values.
    """
    # Add underscore to prefix if it isn't empty
    if prefix != "":
        prefix += "_"

    edges_dict = {}

    for hist in hists:

        # Parse hists calculated on both raw and adj-filtered data
        for hist_type in [f"{prefix}raw_qual_hists", f"{prefix}qual_hists"]:

            hist_name = hist
            if "raw" in hist_type:
                hist_name = f"{hist}_raw"

            edges_dict[hist_name] = "|".join(
                map(
                    lambda x: f"{x:.2f}" if "ab" in hist else str(int(x)),
                    ht.head(1)[hist_type][hist].collect()[0].bin_edges,
                )
            )
    return edges_dict

SUBSET_LIST = ["", "gnomad"]

def populate_info_dict( # TODO: Need to make sure this handles pops correctly and that it adds in both cohort and gnomad info
    bin_edges: Dict[str, str],
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    groups: List[str] = GROUPS,
    subset_list = SUBSET_LIST,
    subset_pops: Dict[str, str] = SUBSET_POPS,
    gnomad_pops: Dict[str, str] = GNOMAD_POPS,
    faf_pops: List[str] = FAF_POPS,
    sexes: List[str] = SEXES,
    analyst_dict: Dict[str, Dict[str, str]] = ANALYST_ANOTATIONS_INFO_DICT,
) -> Dict[str, Dict[str, str]]:
    """
    Calls `make_info_dict` and `make_hist_dict` to populate INFO dictionary with specific sexes, population names, and filtering allele frequency (faf) pops.
    Used during VCF export.
    Creates:
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations
        - INFO fields for variant histograms (hist_bin_freq, hist_n_smaller, hist_n_larger for each histogram)

    :param Dict[str, str] bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param Dict[str, Dict[str, str]] info_dict: INFO dict to be populated.
    :param List[str] groups: List of sample groups [adj, raw]. Default is GROUPS.
    :param Dict[str, str] pops: List of sample global population names for gnomAD genomes. Default is POPS.
    :param List[str] faf_pops: List of faf population names. Default is FAF_POPS.
    :param List[str] sexes: gnomAD sample sexes used in VCF export. Default is SEXES.
    :rtype: Dict[str, Dict[str, str]]
    """
    vcf_info_dict = info_dict.copy()

    # Remove MISSING_REGION_FIELDS from info dict
    for field in MISSING_INFO_FIELDS:
        vcf_info_dict.pop(field, None)

    # Add allele-specific fields to info dict, including AS_VQSLOD and AS_culprit
    # NOTE: need to think about how to resolve AS VQSR fields to avoid having to make temp_AS_fields variable in the future
    temp_AS_fields = AS_FIELDS.copy()
    temp_AS_fields.extend(["AS_culprit", "AS_VQSLOD"])
    vcf_info_dict.update(
        add_as_info_dict(info_dict=info_dict, as_fields=temp_AS_fields)
    )

    def _create_label_groups(
        pops: Dict[str, str], sexes: List[str], group: List[str] = ["adj"],
    ) -> List[Dict[str, List[str]]]:
        """
        Generates list of label group dictionaries needed to populate info dictionary.
        Label dictionaries are passed as input to `make_info_dict`.
        :param Dict[str, str] pops: List of population names.
        :param List[str] sexes: List of sample sexes.
        :param List[str] group: List of data types (adj, raw). Default is ["adj"].
        :return: List of label group dictionaries.
        :rtype: List[Dict[str, List[str]]]
        """
        return [
            dict(group=groups),  # this is to capture raw fields
            dict(group=group, sex=sexes),
            dict(group=group, pop=pops),
            dict(group=group, pop=pops, sex=sexes),
        ]

    for subset in subset_list:

        if "gnomad" in subset:
            description_text = " in gnomAD"

            # faf labels are the same for exomes/genomes
            faf_label_groups = _create_label_groups(pops=faf_pops, sexes=sexes)
            for label_group in faf_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_pops,
                        label_groups=label_group,
                        faf=True,
                        description_text=description_text,
                    )
                )
            gnomad_label_groups = _create_label_groups(
                pops=gnomad_pops, sexes=sexes
            )
            for label_group in gnomad_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset,
                        pop_names=gnomad_pops,
                        label_groups=label_group,
                        description_text=description_text,
                    )
                )

            # Add popmax to info dict
            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    pop_names=gnomad_pops,
                    popmax=True,
                    description_text=description_text,
                )
            )

        else:
            subset_label_groups = _create_label_groups(pops=subset_pops, sexes=sexes)
            for label_group in subset_label_groups:
                vcf_info_dict.update(
                    make_info_dict(
                        prefix=subset, pop_names=subset_pops, label_groups=label_group,
                    )
                )

    # Add variant quality histograms to info dict
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=True))

    # Add Analyst annotations to info_dict
    vcf_info_dict.update(analyst_dict)

    return vcf_info_dict


def make_info_expr(t: Union[hl.MatrixTable, hl.Table]) -> Dict[str, hl.expr.Expression]:
    """
    Makes Hail expression for variant annotations to be included in VCF INFO field.
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
    for field in VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr"][f"{field}"]
    # Add region_flag and allele_info fields to info dict
    for field in ALLELE_TYPE_FIELDS:
        vcf_info_dict[field] = t["allele_info"][f"{field}"]
    for field in REGION_FLAG_FIELDS:
        vcf_info_dict[field] = t["region_flag"][f"{field}"]

    # Add histograms to info dict
    for hist in HISTS:
        for prefix in ["gnomad_qual_hists", "gnomad_raw_qual_hists"]:
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

    # Add analyst annotations to info dict
    vcf_info_dict["cadd_raw_score"] = t["cadd"]["raw_score"]
    vcf_info_dict["cadd_phred"] = t["cadd"]["phred"]

    vcf_info_dict["revel_score"] = t["revel"]["revel_score"]

    vcf_info_dict["splice_ai_max_ds"] = t["splice_ai"]["max_ds"]
    vcf_info_dict["splice_ai_consequence"] = t["splice_ai"]["splice_consequence"]

    vcf_info_dict["primate_ai_score"] = t["primate_ai"]["primate_ai_score"]

    return vcf_info_dict


def make_info_dict(
        prefix: str = "",
        pop_names: Dict[str, str] = POP_NAMES,
        label_groups: Dict[str, str] = None,
        bin_edges: Dict[str, str] = None,
        faf: bool = False,
        popmax: bool = False,
        description_text: str = "",
        age_hist_data: str = None,
        sort_order: List[str] = SORT_ORDER,
) -> Dict[str, Dict[str, str]]:
    """
    Generate dictionary of Number and Description attributes of VCF INFO fields.

    Used to populate the INFO fields of the VCF header during export.

    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex, and subpopulation, both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations

    :param prefix: Prefix string for data, e.g. "gnomAD". Default is empty string.
    :param pop_names: Dict with global population names (keys) and population descriptions (values). Default is POP_NAMES.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations.
    :param popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: Dictionary keyed by VCF INFO annotations, where values are dictionaries of Number and Description attributes.
    """
    if prefix != "":
        prefix = f"{prefix}-"

    info_dict = dict()

    if age_hist_data:
        age_hist_dict = {
            f"{prefix}age_hist_het_bin_freq": {
                "Number": "A",
                "Description": f"Histogram of ages of heterozygous individuals{description_text}; bin edges are: {bin_edges['het']}; total number of individuals of any genotype bin: {age_hist_data}",
            },
            f"{prefix}age_hist_het_n_smaller": {
                "Number": "A",
                "Description": f"Count of age values falling below lowest histogram bin edge for heterozygous individuals{description_text}",
            },
            f"{prefix}age_hist_het_n_larger": {
                "Number": "A",
                "Description": f"Count of age values falling above highest histogram bin edge for heterozygous individuals{description_text}",
            },
            f"{prefix}age_hist_hom_bin_freq": {
                "Number": "A",
                "Description": f"Histogram of ages of homozygous alternate individuals{description_text}; bin edges are: {bin_edges['hom']}; total number of individuals of any genotype bin: {age_hist_data}",
            },
            f"{prefix}age_hist_hom_n_smaller": {
                "Number": "A",
                "Description": f"Count of age values falling below lowest histogram bin edge for homozygous alternate individuals{description_text}",
            },
            f"{prefix}age_hist_hom_n_larger": {
                "Number": "A",
                "Description": f"Count of age values falling above highest histogram bin edge for homozygous alternate individuals{description_text}",
            },
        }
        info_dict.update(age_hist_dict)

    if popmax:
        popmax_dict = {
            f"{prefix}popmax": {
                "Number": "A",
                "Description": f"Population with maximum AF{description_text}",
            },
            f"{prefix}AC_popmax": {
                "Number": "A",
                "Description": f"Allele count in the population with the maximum AF{description_text}",
            },
            f"{prefix}AN_popmax": {
                "Number": "A",
                "Description": f"Total number of alleles in the population with the maximum AF{description_text}",
            },
            f"{prefix}AF_popmax": {
                "Number": "A",
                "Description": f"Maximum allele frequency across populations{description_text}",
            },
            f"{prefix}nhomalt_popmax": {
                "Number": "A",
                "Description": f"Count of homozygous individuals in the population with the maximum allele frequency{description_text}",
            },
            f"{prefix}faf95_popmax": {
                "Number": "A",
                "Description": f"Filtering allele frequency (using Poisson 95% CI) for the population with the maximum allele frequency{description_text}",
            },
        }
        info_dict.update(popmax_dict)

    else:
        group_types = sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
        combos = make_label_combos(label_groups)

        for combo in combos:
            loop_description_text = description_text
            combo_fields = combo.split("-")
            group_dict = dict(zip(group_types, combo_fields))

            for_combo = make_combo_header_text("for", group_dict, prefix, pop_names)
            in_combo = make_combo_header_text("in", group_dict, prefix, pop_names)

            if not faf:
                combo_dict = {
                    f"{prefix}AC-{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele count{for_combo}",
                    },
                    f"{prefix}AN-{combo}": {
                        "Number": "1",
                        "Description": f"Total number of alleles{in_combo}",
                    },
                    f"{prefix}AF-{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele frequency{in_combo}",
                    },
                    f"{prefix}nhomalt-{combo}": {
                        "Number": "A",
                        "Description": f"Count of homozygous individuals{in_combo}",
                    },
                }
            else:
                if ("XX" in combo_fields) | ("XY" in combo_fields):
                    loop_description_text = loop_description_text + " in non-PAR regions of sex chromosomes only"
                combo_dict = {
                    f"{prefix}faf95-{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 95% CI){for_combo}{loop_description_text}",
                    },
                    f"{prefix}faf99-{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 99% CI){for_combo}{loop_description_text}",
                    },
                }
            info_dict.update(combo_dict)
    return info_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table], gnomad_release: bool
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.
    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param List[str] pops: List of global populations in frequency array.
    :return: Dictionary containing variant annotations and their corresponding values.
    :rtype: Dict[str, hl.expr.Expression]
    """
    expr_dict = dict()

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict for UKBB
    if gnomad_release:
        gnomad_prefix = f"gnomad"
        popmax = f"{gnomad_prefix}_popmax"
        faf = f"{gnomad_prefix}_faf"
        freq = f"{gnomad_prefix}_freq"
        faf_idx = hl.eval(t.globals[f"{gnomad_prefix}_faf_index_dict"])
        freq_idx = hl.eval(t.globals[f"{gnomad_prefix}_freq_index_dict"])
    else:
        freq = "cohort_freq"
        freq_idx = hl.eval(t.globals[f"cohort_freq_index_dict"])

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    logger.info("Unfurling freq data...")
    for k, i in freq_idx.items():
        prefix = ""
        if gnomad_release:
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

    # Add popmax
    if gnomad_release:
        prefix = f"{gnomad_prefix}-"

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

    return expr_dict


def build_export_reference() -> hl.ReferenceGenome:
    """
    Creates export reference based on GRCh38. Eliminates all non-standard contigs
    :return: Reference for VCF export containing chr1-22,X,Y, and M
    :rtype: hl.ReferenceGenome
    """
    ref = hl.get_reference("GRCh38")
    my_contigs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
    export_reference = hl.ReferenceGenome(
        name="export_reference",
        contigs=my_contigs,
        lengths={my_contig: ref.lengths[my_contig] for my_contig in my_contigs},
        x_contigs=ref.x_contigs,
        y_contigs=ref.y_contigs,
        par=[(interval.start.contig,
            interval.start.position,
            interval.end.position)
            for interval in ref.par],
        mt_contigs=ref.mt_contigs
    )
    return export_reference


def main(args):

    hl.init(log="/vcf_release.log", default_reference="GRCh38", tmp_dir='hdfs:///vcf_write.tmp/')
    chromosome = args.export_chromosome
    if args.prepare_vcf_ht:
        logger.info("Starting VCF process...")
        logger.info("Getting raw MT and dropping all unnecessary entries...")
        mt = hl.read_matrix_table("gs://gnomad/release/3.1/mt/genomes/hgdp_1kg.gnomad.v3.1.mt")  # TODO: Need to add toe resources

        analyst_ht = hl.read_table(
                'gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1/gnomad_genomes_v3.1.analyst_annotations.ht')  # TODO: update path to use resources
        mt = mt.annotate_rows(**analyst_ht[mt.row_key])
        mt = mt.transmute_cols(sex_karyotype=mt.sex_imputation.sex_karyotype)

        logger.info("Removing chrM...")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chrM")], keep=False)

        if chromosome:
            if chromosome:
                mt = hl.filter_intervals(
                    mt, [hl.parse_locus_interval(chromosome, reference_genome="GRCh38")]
                )

        if args.test:
            logger.info("Filtering to chr20 and chrX (for tests only)...")
            # Using chr20 to test a small autosome and chrX to test a sex chromosome
            # Some annotations (like FAF) are 100% missing on autosomes
            mt = hl.filter_intervals(
                mt,
                [hl.parse_locus_interval("chr20"), hl.parse_locus_interval("chrX")],
            )

        logger.info("Making histogram bin edges...")
        # NOTE: using release HT here because histograms aren't necessarily defined
        # in the first row of the subset MT (we may have filtered that row)
        bin_edges = make_hist_bin_edges_expr(hl.read_table("gs://gnomad/release/3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht"), prefix="")  # TODO: update path to use resources

        logger.info("Making INFO dict for VCF...")
        vcf_info_dict = populate_info_dict(
            bin_edges=bin_edges,
        )

        # Add non-PAR annotation
        mt = mt.annotate_rows(
            region_flag=mt.region_flag.annotate(
                nonpar=(mt.locus.in_x_nonpar() | mt.locus.in_y_nonpar())
            )
        )

        # Unfurl nested subset frequency annotations and add to INFO field
        mt = mt.annotate_rows(release_ht_info=mt.info)
        mt = mt.annotate_rows(
            info=hl.struct(
                **unfurl_nested_annotations(
                    mt,
                    gnomad_release=False,
                )
            )
        )
        # Unfurl nested gnomAD genome frequency annotations and add to info field
        mt = mt.annotate_rows(
            info=mt.info.annotate(
                **unfurl_nested_annotations(
                    mt,
                    gnomad_release=True,
                )
            )
        )

        logger.info("Constructing INFO field")
        # Add variant annotations to INFO field
        # This adds annotations from:
        #   RF struct, VQSR struct, allele_info struct,
        #   info struct (site and allele-specific annotations),
        #   region_flag struct, and
        #   raw_qual_hists/qual_hists structs.

        mt = mt.annotate_rows(info=mt.info.annotate(**make_info_expr(mt)))

        # Reformat vep annotation
        mt = mt.annotate_globals(vep_csq_header=VEP_CSQ_HEADER)
        mt = mt.annotate_rows(vep=vep_struct_to_csq(mt.vep))
        mt = mt.annotate_rows(info=mt.info.annotate(vep=mt.vep))

        # Select relevant fields for VCF export
        mt = mt.select_rows("info", "filters", "rsid")
        vcf_info_dict.update({"vep": {"Description": hl.eval(mt.vep_csq_header)}})

        # Make filter dict
        filter_dict = make_vcf_filter_dict(
            hl.eval(mt.filtering_model.snv_cutoff.min_score),
            hl.eval(mt.filtering_model.indel_cutoff.min_score),
            inbreeding_cutoff=INBREEDING_CUTOFF,
        )

        new_vcf_info_dict = {  # Adjust keys to remove adj tags before exporting to VCF
            i.replace("-adj", ""): j for i, j in vcf_info_dict.items()
        }
        header_dict = {
            "info": new_vcf_info_dict,
            "filter": filter_dict,
            "format": FORMAT_DICT,
        }

        logger.info("Saving header dict to pickle...")
        with hl.hadoop_open(
            f"gs://gnomad-tmp/gnomad_v3.1_vcf_header{chromosome}", "wb"  # TODO: update path to use resources
        ) as p:
            pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)

    if args.sanity_check:
        info_ht = get_info(split=True).ht()
        mt2 = mt.filter_rows(
            ~info_ht[mt.row_key].AS_lowqual
            & ~hl.is_defined(telomeres_and_centromeres.ht()[mt.locus])
            & (hl.len(mt.alleles) > 1)
        )
        sanity_check_release_mt(
            mt2, ["", "gnomad"], missingness_threshold=0.5, verbose=args.verbose
        )

    if args.export_vcf:
        logger.warning(
            "VCF export will densify! Make sure you have an autoscaling cluster."
        )
        logger.info("Reading header dict from pickle...")
        with hl.hadoop_open(f"gs://gnomad-tmp/gnomad_v3.1_vcf_header{chromosome}", "rb") as p:  # TODO: update path to use resources
            header_dict = pickle.load(p)


        logger.info("Dropping histograms that are not needed in VCF...")
        # Drop unnecessary histograms
        drop_hists = (
            [x + "_n_smaller" for x in HISTS if "dp_hist_all" not in x]
            + [x + "_n_larger" for x in HISTS if "dp_" not in x]
            + [x + "_raw_n_smaller" for x in HISTS]
            + [x + "_raw_bin_edges" for x in HISTS]
            + [x + "_raw_n_larger" for x in HISTS]
            + [x + "_raw_bin_freq" for x in HISTS]
            + [x + "_bin_edges" for x in HISTS]
        )
        mt = mt.annotate_rows(info=mt.info.drop(*drop_hists))

        # Reformat names to remove "adj" pre-export
        # e.g, renaming "AC-adj" to "AC"
        # All unlabeled frequency information is assumed to be adj
        logger.info("Dropping 'adj' from info annotations...")
        row_annots = list(mt.info)
        new_row_annots = []
        for x in row_annots:
            x = x.replace("-adj", "")
            new_row_annots.append(x)

        info_annot_mapping = dict(
            zip(new_row_annots, [mt.info[f"{x}"] for x in row_annots])
        )

        mt = mt.transmute_rows(info=hl.struct(**info_annot_mapping))

        logger.info("Running check on VCF fields and info dict...")
        if not vcf_field_check(mt, header_dict, new_row_annots, list(mt.entry)):
            logger.error("Did not pass VCF field check.")

        # Convert int64 fields to int32 (int64 isn't supported by VCF)
        for f, ft in mt.info.dtype.items():
            if ft == hl.dtype("int64"):
                logger.warning(
                    f"Coercing field info.{f} from int64 to int32 for VCF output. Value will be capped at int32 max value."
                )
                mt = mt.annotate_rows(
                    info=mt.info.annotate(
                        **{f: hl.int32(hl.min(2 ** 31 - 1, mt.info[f]))}
                    )
                )
            elif ft == hl.dtype("array<int64>"):
                logger.warning(
                    f"Coercing field info.{f} from array<int64> to array<int32> for VCF output. Array values will be capped at int32 max value."
                )
                mt = mt.annotate_rows(
                    info=mt.info.annotate(
                        **{
                            f: mt.info[f].map(
                                lambda x: hl.int32(hl.min(2 ** 31 - 1, x))
                            )
                        }
                    )
                )


        logger.info("Rearranging fields to desired order...")
        mt = mt.annotate_rows(
            info=mt.info.select(
                "AC",
                "AN",
                "AF",
                "AC-raw",
                "AN-raw",
                "AF-raw",
                "gnomad-AC",
                "gnomad-AN",
                "gnomad-AF",
                "gnomad-popmax",
                "gnomad-faf95_popmax",
                "gnomad-AC-raw",
                "gnomad-AN-raw",
                "gnomad-AF-raw",
                *mt.info.drop(
                    "AC",
                    "AN",
                    "AF",
                    "AC-raw",
                    "AN-raw",
                    "AF-raw",
                    "gnomad-AC",
                    "gnomad-AN",
                    "gnomad-AF",
                    "gnomad-popmax",
                    "gnomad-faf95_popmax",
                    "gnomad-AC-raw",
                    "gnomad-AN-raw",
                    "gnomad-AF-raw",
                ),
            )
        )
        logger.info(f"Export chromosome {chromosome}....")

        logger.info("Densifying and exporting VCF...")
        mt = hl.experimental.densify(mt)

        logger.info("Removing low QUAL variants alleles...")
        info_ht = get_info(split=True).ht()
        mt = mt.filter_rows(
            ~info_ht[mt.row_key].AS_lowqual
            & ~hl.is_defined(telomeres_and_centromeres.ht()[mt.locus])
            & (hl.len(mt.alleles) > 1)
        )

        logger.info("Adjusting sex ploidy...")
        mt = adjust_sex_ploidy(
            mt, mt.sex_karyotype, male_str="XY", female_str="XX"
        )
        mt = mt.select_cols()


        ref = hl.get_reference("GRCh38")
        my_contigs = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY', 'chrM']
        assert (all(my_contig in ref.contigs for my_contig in my_contigs))
        assert (all(x_contig in my_contigs for x_contig in ref.x_contigs))
        assert (all(y_contig in my_contigs for y_contig in ref.y_contigs))
        assert (all(mt_contig in my_contigs for mt_contig in ref.mt_contigs))
        my_reference = hl.ReferenceGenome(
            name="gnomAD_GRCh38",
            contigs=my_contigs,
            lengths={my_contig: ref.lengths[my_contig]
                        for my_contig in my_contigs},
            x_contigs=ref.x_contigs,
            y_contigs=ref.y_contigs,
            par=[(interval.start.contig,
                  interval.start.position,
                  interval.end.position)
                 for interval in ref.par],
            mt_contigs=ref.mt_contigs
        )
        assert (
                mt.aggregate_rows(
                    hl.agg.count_where(
                        ~hl.set(my_reference.contigs).contains(mt.locus.contig)
                    )
                ) == 0
        )
        mt = mt.rename({"locus": "locus_original"})
        mt = mt.annotate_rows(
            locus=hl.locus(
                mt.locus_original.contig,
                mt.locus_original.position,
                reference_genome=my_reference
            )
        )
        mt = mt.key_rows_by("locus", "alleles").drop("locus_original")

        hl.export_vcf(
            mt,
            f"gs://gnomad/release/3.1/vcf/genomes/gnomad.genomes.v3.1.hgdp_1kg_subset.{chromosome}.vcf.bgz",  # TODO: update path to use resources
            append_to_header="gs://gnomad/release/3.1/vcf/genomes/extra_fields_for_header_HGDP_TGP.tsv",  # TODO: update path to use resources
            metadata=header_dict,
            tabix=True,
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test",
        help="Create release files using only chr20 and chrX for testing purposes",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_vcf_ht", help="Use release ht to create vcf ht", action="store_true"
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