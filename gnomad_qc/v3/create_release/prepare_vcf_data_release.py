import argparse
import json
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import vep_struct_to_csq
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
    INFO_DICT,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
    make_label_combos,
    make_vcf_filter_dict,
    REGION_FLAG_FIELDS,
    SITE_FIELDS,
    set_female_y_metrics_to_na,
    VQSR_FIELDS,
    remove_fields_from_globals,
)

from gnomad_qc.v3.create_release.sanity_checks import sanity_check_release_ht

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

VEP_CSQ_FIELDS = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|LoF_flags|LoF_info"
"""
Constant that defines the order of VEP annotations used in VCF export.
"""

VEP_CSQ_HEADER = f"Consequence annotations from Ensembl VEP. Format: {VEP_CSQ_FIELDS}"
"""
Constant that contains description for VEP used in VCF export.
"""

SEXES = ["XX", "XY"]

# Make subset list (used in properly filling out VCF header descriptions and naming VCF info fields)
SUBSET_LIST = [
    "controls_and_biobanks",
    "non_cancer",
    "non_neuro",
    "non_topmed",
    "non_v2",
    "hgdp",
    "tgp",
]

SUBSET_LIST_FOR_VCF = [
    "controls_and_biobanks",
    "non_cancer",
    "non_neuro",
    "non_topmed",
    "non_v2",
    "",
]

# Select populations to keep from the list of population names in POP_NAMES
# This removes pop names we don't want to use
# (e.g., "uniform", "consanguineous") to reduce clutter
KEEP_POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]

# Remove unnecessary pop names from pops dict
POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}
POPS["mid"] = "Middle Eastern"

# downsampling and subset entries to remove from VCF's freq export
FREQ_ENTRIES_TO_REMOVE = [
    "10",
    "20",
    "50",
    "100",
    "158",
    "200",
    "456",
    "500",
    "1000",
    "1047",
    "1736",
    "2000",
    "2419",
    "2604",
    "5000",
    "5316",
    "7647",
    "10000",
    "15000",
    "20000",
    "20744",
    "25000",
    "30000",
    "34029",
    "40000",
    "60000",
    "70000",
    "75000",
    "hgdp",
    "tgp",
]

EXPORT_HISTS = [
    "gq_hist_alt_bin_freq",
    "gq_hist_all_bin_freq",
    "dp_hist_alt_bin_freq",
    "dp_hist_alt_n_larger",
    "dp_hist_all_bin_freq",
    "dp_hist_all_n_larger",
    "ab_hist_alt_bin_freq",
]


INBREEDING_CUTOFF = -0.3


def release_ht_path():
    return "gs://gnomad-mwilson/untitled-folder/release_test.ht"  # "gs://gnomad/release/v3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht"


def populate_info_dict(
    bin_edges: Dict[str, str],
    age_hist_data: str,
    info_dict: Dict[str, Dict[str, str]] = VCF_INFO_DICT,
    subset_list: List[str] = SUBSET_LIST,
    groups: List[str] = GROUPS,
    pops: Dict[str, str] = POPS,
    faf_pops: List[str] = FAF_POPS,
    sexes: List[str] = SEXES,
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
    :param List[str] subset_list: List of sample subsets in dataset. Default is SUBSET_LIST.
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
                prefix="",
                pop_names=pops,
                label_groups=label_group,
                faf=True,
                description_text=description_text,
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
    vcf_info_dict.update(
        make_hist_dict(bin_edges, adj=True)
    )  # , dict_hists=EXPORT_HISTS)) <-this fails because bins go with the prefix of the hist, not specific bin
    return vcf_info_dict


def make_info_expr(t: Union[hl.MatrixTable, hl.Table]) -> Dict[str, hl.expr.Expression]:
    """
    Makes Hail expression for variant annotations to be included in VCF INFO field.
    :param Table/MatrixTable t: Table/MatrixTable containing variant annotations to be reformatted for VCF export.
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    :rtype: Dict[str, hl.expr.Expression]
    """
    # Start info dict with region_flag and allele_info fields
    vcf_info_dict = {}
    for field in ALLELE_TYPE_FIELDS:
        vcf_info_dict[field] = t["allele_info"][f"{field}"]
    for field in REGION_FLAG_FIELDS:
        vcf_info_dict[field] = t["region_flag"][f"{field}"]

    # Add site-level annotations to vcf_info_dict
    for field in SITE_FIELDS:
        vcf_info_dict[field] = t["info"][f"{field}"]
    for field in VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr"][f"{field}"]

    # Add AS annotations to info dict
    for field in AS_FIELDS:
        vcf_info_dict[field] = t["info"][f"{field}"]

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

    vcf_info_dict["splice_ai_max_ds"] = t["splice_ai"]["max_ds"]
    vcf_info_dict["splice_ai_consequence"] = t["splice_ai"]["splice_consequence"]

    vcf_info_dict["primate_ai_score"] = t["primate_ai"]["primate_ai_score"]

    return vcf_info_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table], pops: List[str], entries_to_remove: List[str]
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

    combo_dict = {
        "popmax": t[popmax].pop,
        "AC_popmax": t[popmax].AC,
        "AN_popmax": t[popmax].AN,
        "AF_popmax": t[popmax].AF,
        "nhomalt_popmax": t[popmax].homozygote_count,
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

    hl.init(log="/vcf_release.log", default_reference="GRCh38")
    try:

        if args.prepare_vcf_ht:
            logger.info("Starting VCF process...")
            logger.info("Reading in release HT...")
            ht = hl.read_table(release_ht_path())

            logger.info("Removing chrM...")
            ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)

            if args.export_chromosome:
                ht = hl.filter_intervals(
                    ht, [hl.parse_locus_interval(args.export_chromosome)]
                )

            if args.test:
                logger.info("Filtering to chr20 and chrX (for tests only)...")
                # Using chr20 to test a small autosome and chrX to test a sex chromosome
                # Some annotations (like FAF) are 100% missing on autosomes
                ht = hl.filter_intervals(
                    ht,
                    [hl.parse_locus_interval("chr20"), hl.parse_locus_interval("chrX")],
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

            logger.info("Constructing INFO field")
            # Add variant annotations to INFO field
            # This adds annotations from:
            #   RF struct, VQSR struct, allele_info struct,
            #   info struct (site and allele-specific annotations),
            #   region_flag struct, and
            #   raw_qual_hists/qual_hists structs.

            ht = ht.annotate(info=hl.struct(**make_info_expr(ht)))

            # Unfurl nested gnomAD frequency annotations and add to info field
            ht = ht.annotate(
                info=ht.info.annotate(
                    **unfurl_nested_annotations(
                        ht, pops=POPS, entries_to_remove=FREQ_ENTRIES_TO_REMOVE
                    )
                )
            )
            # ht = ht.annotate(**set_female_y_metrics_to_na(ht))

            # Reformat vep annotation
            ht = ht.annotate_globals(vep_csq_header=VEP_CSQ_HEADER)
            ht = ht.annotate(vep=vep_struct_to_csq(ht.vep))
            ht = ht.annotate(info=ht.info.annotate(vep=ht.vep))

            # Select relevant fields for VCF export
            ht = ht.select("info", "filters", "rsid")
            vcf_info_dict.update({"vep": {"Description": hl.eval(ht.vep_csq_header)}})

            # Make filter dict
            filter_dict = make_vcf_filter_dict(
                hl.eval(ht.filtering_model.snv_cutoff.min_score),
                hl.eval(ht.filtering_model.indel_cutoff.min_score),
                inbreeding_cutoff=INBREEDING_CUTOFF,
            )

            new_vcf_info_dict = {  # Adjust keys to remove adj tags before exporting to VCF
                i.replace("-adj", ""): j for i, j in vcf_info_dict.items()
            }
            header_dict = {
                "info": new_vcf_info_dict,
                "filter": filter_dict,
            }

            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(
                "gs://gnomad-mwilson/v3.1/release/vcf_header", "wb"
            ) as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)

        if args.sanity_check:
            sanity_check_release_ht(
                ht, SUBSET_LIST, missingness_threshold=0.5, verbose=args.verbose
            )

        if args.export_vcf:
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
            row_annots = list(ht.info)
            new_row_annots = []
            for x in row_annots:
                x = x.replace("-adj", "")
                new_row_annots.append(x)

            info_annot_mapping = dict(
                zip(new_row_annots, [ht.info[f"{x}"] for x in row_annots])
            )

            ht = ht.transmute(info=hl.struct(**info_annot_mapping))

            ht.describe()
            print(header_dict)
            hl.export_vcf(
                ht,
                "gs://gnomad-mwilson/untitled-folder/release_test_w_append.vcf.bgz",
                append_to_header="gs://gnomad-mwilson/untitled-folder/vcf_header_non_info.tsv",
                metadata=header_dict,
                tabix=True,
            )

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log("gs://gnomad-mwilson/logs/v3.1/sanity_logs/10202020.log")


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
