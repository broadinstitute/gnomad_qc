import argparse
import copy
import itertools
import logging
import pickle
from typing import Dict, List, Union

import hail as hl

from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.resources.grch38.gnomad import (
    COHORTS_WITH_POP_STORED_AS_SUBPOP,
    DOWNSAMPLINGS,
    POPS,
    REGION_FLAG_FIELDS,
    SEXES,
    SUBSETS,
)
from gnomad.utils.vep import VEP_CSQ_HEADER, vep_struct_to_csq
from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    FAF_POPS,
    GROUPS,
    HISTS,
    INFO_DICT,
    INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
    make_combo_header_text,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_label_combos,
    SORT_ORDER,
    SITE_FIELDS,
    VQSR_FIELDS,
)
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF

from gnomad_qc.v3.create_release.sanity_checks import (
    sanity_check_release_ht,
    vcf_field_check,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)

# Add capture region and sibling singletons to vcf_info_dict
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


def remove_fields_from_globals(global_field: List[str], fields_to_remove: List[str]):
    """
    Removes fields from the pre-defined global field variables.

    :param global_field: Global list of fields
    :param fields_to_remove: List of fields to remove from global (they must be in the global list)
    """
    for field in fields_to_remove:
        if field in global_field:
            global_field.remove(field)
        else:
            logger.info(f"'{field}'' missing from {global_field}")


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
SUBSET_LIST_FOR_VCF = SUBSETS
SUBSET_LIST_FOR_VCF.append("")
remove_fields_from_globals(SUBSET_LIST_FOR_VCF, COHORTS_WITH_POP_STORED_AS_SUBPOP)

# Remove unnecessary pop names from pops dict
POPS = {pop: POP_NAMES[pop] for pop in POPS}

# downsampling and subset entries to remove from VCF's freq export
FREQ_ENTRIES_TO_REMOVE = DOWNSAMPLINGS + COHORTS_WITH_POP_STORED_AS_SUBPOP

ANALYST_ANNOTATIONS_INFO_DICT = {
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


def ht_to_vcf_mt(
    info_ht: hl.Table,
    create_mt: hl.bool = False,
    pipe_delimited_annotations: List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.MatrixTable:
    """
    Creates a MT ready for vcf export from a HT. In particular, the following conversions are done:
    - All int64 are coerced to int32
    - Fields specified by `pipe_delimited_annotations` will be converted from arrays to pipe-delimited strings

    .. note::

        The MT returned has no cols.

    :param info_ht: Input HT
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct)
    :return: MatrixTable ready for VCF export
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    info_ht = info_ht.key_by("locus", "alleles")

    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in info_ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            logger.warning(
                f"Coercing field info.{f} from int64 to int32 for VCF output. Value will be capped at int32 max value."
            )
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{f: hl.int32(hl.min(2 ** 31 - 1, info_ht.info[f]))}
                )
            )
        elif ft == hl.dtype("array<int64>"):
            logger.warning(
                f"Coercing field info.{f} from array<int64> to array<int32> for VCF output. Array values will be capped at int32 max value."
            )
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{
                        f: info_ht.info[f].map(
                            lambda x: hl.int32(hl.min(2 ** 31 - 1, x))
                        )
                    }
                )
            )

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, but it is overwritten below
    for f in pipe_delimited_annotations:
        if f in info_ht.info and f != "AS_SB_TABLE":
            info_expr[f] = "|" + get_pipe_expr(info_ht.info[f])

    # Flatten SB if it is an array of arrays
    if "SB" in info_ht.info and not isinstance(
        info_ht.info.SB, hl.expr.ArrayNumericExpression
    ):
        info_expr["SB"] = info_ht.info.SB[0].extend(info_ht.info.SB[1])

    if "AS_SB_TABLE" in info_ht.info:
        info_expr["AS_SB_TABLE"] = get_pipe_expr(
            hl.array([info_ht.info.AS_SB_TABLE[:2], info_ht.info.AS_SB_TABLE[2:]]).map(lambda x: hl.delimit(x, ","))
        )

    # Annotate with new expression
    info_ht = info_ht.annotate(info=info_ht.info.annotate(**info_expr))
    info_t = info_ht

    if create_mt:
        # Add 's' empty string field required to cast HT to MT
        info_t = info_t.annotate(s=hl.null(hl.tstr))

        # Create an MT with no cols so that we can export to VCF
        info_t = info_t.to_matrix_table_row_major(columns=["s"], entry_field_name="s")
        info_t = info_t.filter_cols(False)

    return info_t


def make_info_dict(
    prefix: str = "",
    pop_names: Dict[str, str] = POP_NAMES,
    label_groups: Dict[str, str] = None,
    label_delimiter: str = "_",
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
    :param label_delimiter: String to use as delimiter when making group label combinations.
    :param bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation.
    :param faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations.
    :param popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations.
    :param description_text: Optional text to append to the end of descriptions. Needs to start with a space if specified.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: Dictionary keyed by VCF INFO annotations, where values are dictionaries of Number and Description attributes.
    """
    if prefix != "":
        prefix = f"{prefix}_" #TODO: DOES THIS NEED TO BE {label_delimiter}?

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
        combos = make_label_combos(label_groups, label_delimiter=label_delimiter)

        for combo in combos:
            loop_description_text = description_text
            combo_fields = combo.split(label_delimiter)
            group_dict = dict(zip(group_types, combo_fields))

            for_combo = make_combo_header_text("for", group_dict, prefix, pop_names)
            in_combo = make_combo_header_text("in", group_dict, prefix, pop_names)

            #TODO: IN COMPARISON WITH COMMON CODE THERE IS DIFFERENT ORDERING, WHICH DO WE WANT MOVING FORWARD?
            if not faf:
                combo_dict = {
                    f"AC-{prefix}{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele count{for_combo}{loop_description_text}",
                    },
                    f"AN-{prefix}{combo}": {
                        "Number": "1",
                        "Description": f"Total number of alleles{in_combo}{loop_description_text}",
                    },
                    f"AF-{prefix}{combo}": {
                        "Number": "A",
                        "Description": f"Alternate allele frequency{in_combo}{loop_description_text}",
                    },
                    f"nhomalt-{prefix}{combo}": {
                        "Number": "A",
                        "Description": f"Count of homozygous individuals{in_combo}{loop_description_text}",
                    },
                }
            else:
                # TODO: THE IF STATEMENT BELOW IS NOT IN THE COMMON CODE, WHICH DO WE WANT?
                if ("XX" in combo_fields) | ("XY" in combo_fields):
                    loop_description_text = (
                        loop_description_text
                        + " in non-PAR regions of sex chromosomes only"
                    )
                combo_dict = {
                    f"faf95-{prefix}{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 95% CI){for_combo}{loop_description_text}",
                    },
                    f"faf99-{prefix}{combo}": {
                        "Number": "A",
                        "Description": f"Filtering allele frequency (using Poisson 99% CI){for_combo}{loop_description_text}",
                    },
                }
            info_dict.update(combo_dict)
    return info_dict


def make_vcf_filter_dict(
    snp_cutoff: float, indel_cutoff: float, inbreeding_cutoff: float
) -> Dict[str, str]:
    """
    Generates dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations.

    Generates descriptions for:
        - AC0 filter
        - InbreedingCoeff filter
        - RF filter #TODO: GENERALIZE WE USED AS_VQSR FOR V#
        - PASS (passed all variant filters)

    :param snp_cutoff: Minimum SNP cutoff score from random forest model.
    :param indel_cutoff: Minimum indel cutoff score from random forest model.
    :param inbreeding_cutoff: Inbreeding coefficient hard cutoff.
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes.
    """
    filter_dict = {
        "AC0": {
            "Description": "Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)"
        },
        "InbreedingCoeff": {
            "Description": f"InbreedingCoeff < {inbreeding_cutoff}"
        },
        #TODO: CHANGE COMMON CODE TO TAKE DIFFERENT FILTERING METHOD
        "RF": {
            "Description": f"Failed random forest filtering thresholds of {snp_cutoff} for SNPs and {indel_cutoff} for indels (probabilities of being a true positive variant)"
        },
        "PASS": {"Description": "Passed all variant filters"},
    }
    return filter_dict


#TODO: THIS LOOKS EXACTLY LIKE COMMON CODE JUST USE THAT AS IS
def add_as_info_dict(
    info_dict: Dict[str, Dict[str, str]] = INFO_DICT, as_fields: List[str] = AS_FIELDS
) -> Dict[str, Dict[str, str]]:
    """
    Updates info dictionary with allele-specific terms and their descriptions.

    Used in VCF export.

    :param info_dict: Dictionary containing site-level annotations and their descriptions. Default is INFO_DICT.
    :param as_fields: List containing allele-specific fields to be added to info_dict. Default is AS_FIELDS.
    :return: Dictionary with allele specific annotations, their descriptions, and their VCF number field.
    """
    as_dict = {}
    for field in as_fields:

        try:
            # Strip AS_ from field name
            site_field = field[3:]

            # Get site description from info dictionary and make first letter lower case
            first_letter = info_dict[site_field]["Description"][0].lower()
            rest_of_description = info_dict[site_field]["Description"][1:]

            as_dict[field] = {}
            as_dict[field]["Number"] = "A"
            as_dict[field][
                "Description"
            ] = f"Allele-specific {first_letter}{rest_of_description}"

        except KeyError:
            logger.warning(f"{field} is not present in input info dictionary!")

    return as_dict


#TODO:Is this something that should go in gnomad_methods? I am not sure of the standards around globals,
# I always had this notion that you don’t modify globals, but this is an odd case because it isn’t modifying the
# global in the file that contains the global and basically making a global for the current file so maybe it is OK
def remove_fields_from_globals(global_field: List[str], fields_to_remove: List[str]):
    """
    Removes fields from the pre-defined global field variables.

    :param global_field: Global list of fields
    :param fields_to_remove: List of fields to remove from global (they must be in the global list)
    """
    for field in fields_to_remove:
        if field in global_field:
            global_field.remove(field)
        else:
            logger.info(f"'{field}'' missing from {global_field}")


#TODO: SHOULD THIS BE IN A MORE COMMON LOCATION?
def build_export_reference() -> hl.ReferenceGenome:
    """
    Creates export reference based on GRCh38. Eliminates all non-standard contigs
    :return: Reference for VCF export containing chr1-22,X,Y, and M
    :rtype: hl.ReferenceGenome
    """
    ref = hl.get_reference("GRCh38")
    my_contigs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]
    export_reference = hl.ReferenceGenome(
        name="gnomAD_GRCh38",
        contigs=my_contigs,
        lengths={my_contig: ref.lengths[my_contig] for my_contig in my_contigs},
        x_contigs=ref.x_contigs,
        y_contigs=ref.y_contigs,
        par=[
            (interval.start.contig, interval.start.position, interval.end.position)
            for interval in ref.par
        ],
        mt_contigs=ref.mt_contigs,
    )
    return export_reference

#TODO:USE RESOURCES
def release_ht_path():
    return "gs://gnomad/release/3.1/ht/genomes/gnomad.genomes.v3.1.sites.reference_fixed.ht"


def populate_info_dict(
    bin_edges: Dict[str, str],
    age_hist_data: str,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    subset_list: List[str] = SUBSETS,
    groups: List[str] = GROUPS,
    pops: Dict[str, str] = POPS,
    faf_pops: List[str] = FAF_POPS,
    sexes: List[str] = SEXES,
    analyst_dict: Dict[str, Dict[str, str]] = ANALYST_ANNOTATIONS_INFO_DICT,
) -> Dict[str, Dict[str, str]]:
    """
    Calls `make_info_dict` and `make_hist_dict` to populate INFO dictionary with specific sexes, population names, and filtering allele frequency (faf) pops.
    Used during VCF export.
    Creates:
        - INFO fields for age histograms (bin freq, n_smaller, and n_larger for heterozygous and homozygous variant carriers)
        - INFO fields for popmax AC, AN, AF, nhomalt, and popmax population
        - INFO fields for AC, AN, AF, nhomalt for each combination of sample population, sex both for adj and raw data
        - INFO fields for filtering allele frequency (faf) annotations
        - INFO fields for variant histograms (hist_bin_freq, hist_n_smaller, hist_n_larger for each histogram)
    :param Dict[str, List[str]] subpops: Dictionary of global population names (keys)
        and all hybrid population cluster names associated with that global pop (values).
    :param Dict[str, str] bin_edges: Dictionary of variant annotation histograms and their associated bin edges.
    :param str age_hist_data: Pipe-delimited string of age histograms, from `get_age_distributions`.
    :param Dict[str, Dict[str, str]] info_dict: INFO dict to be populated.
    :param List[str] subset_list: List of sample subsets in dataset. Default is SUBSETS.
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
            dict(group=groups),  # this is to capture high level raw fields
            dict(group=group, sex=sexes),
            dict(group=group, pop=pops),
            dict(group=group, pop=pops, sex=sexes),
        ]

    for subset in subset_list:
        description_text = "" if subset == "" else f" in {subset} subset"

        label_groups = _create_label_groups(pops=pops, sexes=sexes)
        for label_group in label_groups:
            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset,
                    pop_names=pops,
                    label_groups=label_group,
                    description_text=description_text,
                )
            )

    faf_label_groups = _create_label_groups(pops=faf_pops, sexes=sexes)
    for label_group in faf_label_groups:
        vcf_info_dict.update(
            make_info_dict(
                prefix="", pop_names=pops, label_groups=label_group, faf=True,
            )
        )

    vcf_info_dict.update(
        make_info_dict(
            prefix="",
            bin_edges=bin_edges,
            popmax=True,
            age_hist_data="|".join(str(x) for x in age_hist_data),
        )
    )

    # Add variant quality histograms to info dict
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=True, label_delimiter="-"))

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

    # Add analyst annotations to info dict
    vcf_info_dict["cadd_raw_score"] = t["cadd"]["raw_score"]
    vcf_info_dict["cadd_phred"] = t["cadd"]["phred"]

    vcf_info_dict["revel_score"] = t["revel"]["revel_score"]

    vcf_info_dict["splice_ai_max_ds"] = t["splice_ai"]["splice_ai_score"]
    vcf_info_dict["splice_ai_consequence"] = t["splice_ai"]["splice_consequence"]

    vcf_info_dict["primate_ai_score"] = t["primate_ai"]["primate_ai_score"]

    return vcf_info_dict


def unfurl_nested_annotations(
        t: Union[hl.MatrixTable, hl.Table], pops: List[str], entries_to_remove: List[str] = None
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

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict
    popmax = "popmax"

    freq = "freq"
    freq_idx = hl.eval(t.freq_index_dict)

    if entries_to_remove:
        # freqs to remove for vcf_export
        new_freq_idx = freq_idx
        for entry in entries_to_remove:
            new_freq_idx = {
                key: val
                for key, val in new_freq_idx.items()
                if not key.startswith(entry)
            }
        freq_idx = new_freq_idx

    faf = "faf"
    faf_idx = hl.eval(t.faf_index_dict)

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    logger.info("Unfurling freq data...")
    for k, i in freq_idx.items():
        # Set combination to key
        # e.g., set entry of 'afr_adj' to combo
        combo = k
        combo_dict = {
            f"AC-{combo}": t[freq][i].AC,
            f"AN-{combo}": t[freq][i].AN,
            f"AF-{combo}": t[freq][i].AF,
            f"nhomalt-{combo}": t[freq][i].homozygote_count,
        }
        expr_dict.update(combo_dict)

    # Add popmax
    combo_dict = {
        "popmax": t[popmax].pop,
        "AC_popmax": t[popmax].AC,
        "AN_popmax": t[popmax].AN,
        "AF_popmax": t[popmax].AF,
        "nhomalt_popmax": t[popmax].homozygote_count,
        "faf95_popmax": t[popmax].faf95,
    }
    expr_dict.update(combo_dict)

    ## Unfurl FAF index dict
    logger.info("Unfurling faf data...")
    for (
        k,
        i,
    ) in faf_idx.items():  # NOTE: faf annotations are all done on adj-only groupings
        combo_dict = {
            f"faf95-{k}": t[faf][i].faf95,
            f"faf99-{k}": t[faf][i].faf99,
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

    return expr_dict


def main(args):

    hl.init(
        log="/vcf_release.log",
        default_reference="GRCh38",
        tmp_dir="hdfs:///vcf_write.tmp/",
    )
    try:

        if args.prepare_vcf_ht:
            chromosome = args.export_chromosome
            logger.info("Starting VCF process...")
            logger.info("Reading in release HT...")
            ht = hl.read_table(release_ht_path())
            export_reference = build_export_reference()

            #TODO: Confirm that this is needed
            ht = ht.rename({"locus": "locus_original"})
            ht = ht.annotate(
                locus=hl.locus(
                    ht.locus_original.contig,
                    ht.locus_original.position,
                    reference_genome=export_reference
                )
            )
            ht = ht.key_by("locus", "alleles").drop("locus_original")

            if chromosome:
                ht = hl.filter_intervals(
                    ht,
                    [
                        hl.parse_locus_interval(
                            chromosome, reference_genome=export_reference
                        )
                    ],
                )
                #TODO: MAKE this smaller test an option!!!
                # ht = ht._filter_partitions(range(1))

            if args.test:
                logger.info("Filtering to chr20, chrX, and chrY (for tests only)...")
                # Using chr20 to test a small autosome and chrX/chrY to test sex chromosomes
                ht = hl.filter_intervals(
                    ht,
                    [
                        hl.parse_locus_interval(
                            "chr20", reference_genome=export_reference
                        ),
                        hl.parse_locus_interval(
                            "chrX", reference_genome=export_reference
                        ),
                        hl.parse_locus_interval(
                            "chrY", reference_genome=export_reference
                        ),
                    ],
                )

            logger.info("Making histogram bin edges...")
            bin_edges = make_hist_bin_edges_expr(ht, prefix="")

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(ht.age_distribution)

            logger.info("Making INFO dict for VCF...")
            vcf_info_dict = populate_info_dict(
                bin_edges=bin_edges,
                age_hist_data=age_hist_data,
                subset_list=SUBSET_LIST_FOR_VCF,
            )

            # Add non-PAR annotation
            ht = ht.annotate(
                region_flag=ht.region_flag.annotate(
                    nonpar=(ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar())
                )
            )

            # Unfurl nested gnomAD frequency annotations and add to info field
            ht = ht.annotate(release_ht_info=ht.info)
            ht = ht.annotate(
                info=hl.struct(
                    **unfurl_nested_annotations(
                        ht, pops=POPS, entries_to_remove=FREQ_ENTRIES_TO_REMOVE
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

            ht = ht.annotate(info=ht.info.annotate(**make_info_expr(ht)))

            # Reformat vep annotation
            ht = ht.annotate_globals(vep_csq_header=VEP_CSQ_HEADER)
            ht = ht.annotate(vep=vep_struct_to_csq(ht.vep))
            ht = ht.annotate(info=ht.info.annotate(vep=ht.vep))

            # NOTE: Merging rsid set into a semi-colon delimited string
            # dbsnp might have multiple identifiers for one variant
            # thus, rsid is a set annotation, starting with version b154 for dbsnp resource:
            # https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/resources/grch38/reference_data.py#L136
            # `export_vcf` expects this field to be a string, and vcf specs
            # say this field may be delimited by a semi-colon:
            # https://samtools.github.io/hts-specs/VCFv4.2.pdf
            logger.info("Reformatting rsid...")
            ht = ht.annotate(rsid=hl.str(";").join(ht.rsid))

            # Select relevant fields for VCF export
            ht = ht.select("info", "filters", "rsid")
            vcf_info_dict.update({"vep": {"Description": hl.eval(ht.vep_csq_header)}})

            #TODO: ADD TO RESOURCES, but this checkpoint really helps for export!
            ht = ht.checkpoint(f"gs://gnomad-tmp/gnomad_v3.1_vcfs/vcf_ht_checkpoint_chr_all.ht", overwrite=True)

            # Make filter dict
            filter_dict = make_vcf_filter_dict(
                hl.eval(ht.filtering_model.snv_cutoff.min_score),
                hl.eval(ht.filtering_model.indel_cutoff.min_score),
                inbreeding_cutoff=INBREEDING_COEFF_HARD_CUTOFF,
            )

            # Adjust keys to remove adj tags before exporting to VCF
            new_vcf_info_dict = {}
            for i, j in vcf_info_dict.items():
                i = i.replace("-adj", "")
                i = i.replace(
                    "-", "_"
                )  # VCF 4.3 specs do not allow hyphens in info fields
                new_vcf_info_dict[i] = j

            header_dict = {
                "info": new_vcf_info_dict,
                "filter": filter_dict,
            }

            # TODO: CHANGE TO RESOURCE LOCATION
            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(
                    "gs://gnomad-mwilson/v3.1.1/release/vcf_header", "wb"
            ) as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)

        if args.sanity_check:
            #TODO: MIGHT NEED TO ADD TO sanity_check_release_ht, I think I had trouble with this somethimes, maybe needing to use or not use the reference_genome? will need to test
            sanity_check_release_ht(
                ht, SUBSETS, missingness_threshold=0.5, verbose=args.verbose, reference_genome=export_reference
            )

        if args.export_vcf:
            chromosome = args.export_chromosome
            #TODO: AGAIN ADD TO RESOURCES FOR BOTH BELOW
            ht = hl.read_table(f"gs://gnomad-tmp/gnomad_v3.1_vcfs/vcf_ht_checkpoint_chr_all.ht")
            with hl.hadoop_open(
                    "gs://gnomad-mwilson/v3.1.1/release/vcf_header",
                    "rb",
            ) as f:
                header_dict = pickle.load(f)

            if chromosome:
                export_reference = build_export_reference()
                ht = hl.filter_intervals(
                    ht,
                    [
                        hl.parse_locus_interval(
                            chromosome, reference_genome=export_reference
                        )
                    ],
                )

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
                + ["age_hist_het_bin_edges", "age_hist_hom_bin_edges"]
            )
            ht = ht.annotate(info=ht.info.drop(*drop_hists))

            # Reformat names to remove "adj" pre-export
            # e.g, renaming "AC-adj" to "AC"
            # All unlabeled frequency information is assumed to be adj
            logger.info("Dropping 'adj' from info annotations...")
            row_annots = list(ht.info)
            new_row_annots = []
            for x in row_annots:
                x = x.replace("-adj", "")
                x = x.replace(
                    "-", "_"
                )  # VCF 4.3 specs do not allow hyphens in info fields
                new_row_annots.append(x)

            info_annot_mapping = dict(
                zip(new_row_annots, [ht.info[f"{x}"] for x in row_annots])
            )

            ht = ht.transmute(info=hl.struct(**info_annot_mapping))

            logger.info("Running check on VCF fields and info dict...")
            if not vcf_field_check(ht, header_dict, new_row_annots):
                logger.error("Did not pass VCF field check")

            logger.info("Adjusting VCF incompatiable types...")
            ht = ht_to_vcf_mt(ht)

            logger.info("Rearranging fields to desired order...")
            ht = ht.annotate(
                info=ht.info.select(
                    "AC",
                    "AN",
                    "AF",
                    "popmax",
                    "faf95_popmax",
                    *ht.info.drop("AC", "AN", "AF", "popmax", "faf95_popmax"),
                )
            )
            ht.describe()
            logger.info(f"Export chromosome {chromosome}....")
            #TODO: ADD BOTH PATHS BELOW TO RESOURCES!!
            hl.export_vcf(
                ht,
                f"gs://gnomad/release/3.1.1/vcf/genomes/gnomad.genomes.v3.1.1.sites.{chromosome}.vcf.bgz",
                append_to_header="gs://gnomad/release/3.1.1/vcf/genomes/extra_fields_for_header.tsv",
                metadata=header_dict,
                tabix=True,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        #TODO: ADD TO RESOURCES!
        hl.copy_log("gs://gnomad-mwilson/logs/v3.1.1/vcf_export.log")


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
