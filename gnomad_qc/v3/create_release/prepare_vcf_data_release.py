import argparse
import json
import logging
import pickle
import sys
from typing import Dict, List, Union

import hail as hl

from gnomad_methods.gnomad.sample_qc.ancestry import POP_NAMES
from gnomad_methods.gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad_methods.gnomad.utils.reference_genome import get_reference_genome
from gnomad_methods.gnomad.utils.slack import slack_notifications
from gnomad_methods.gnomad.utils.vep import vep_struct_to_csq
from gnomad_methods.gnomad.utils.vcf import (
    add_as_info_dict,
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
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
    RF_FIELDS,
    SITE_FIELDS,
    set_female_y_metrics_to_na,
    SEXES,
    VQSR_FIELDS,
)

from gnomad_qc.gnomad_qc.v3.resources.release import release_ht_path

from ukbb_qc.assessment.sanity_checks import (
    sanity_check_release_mt,
    vcf_field_check,
)
from ukbb_qc.resources.basics import (
    get_checkpoint_path,
    get_ukbb_data,
    logging_path,
    release_header_path,
    release_mt_path,
    release_vcf_path,
)

from ukbb_qc.resources.variant_qc import info_ht_path
from ukbb_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)

# Add capture region and sibling singletons to vcf_info_dict
VCF_INFO_DICT = INFO_DICT

# Remove decoy from region flag field
REGION_FLAG_FIELDS.remove("decoy")

# Missing region field flag to remove from vcf info dict
MISSING_REGION_FIELDS = "decoy"

# Remove BaseQRankSum from site fields (doesn't exist in UKBB 300K)
SITE_FIELDS.remove("BaseQRankSum")

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

# Select populations to keep from the list of population names in POP_NAMES
# This removes pop names we don't want to use
# (e.g., "uniform", "consanguineous") to reduce clutter
KEEP_POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"]

# Remove unnecessary pop names from pops dict
POPS = {pop: POP_NAMES[pop] for pop in KEEP_POPS}
POPS["mid"] = "Middle Eastern"


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
    :param List[str] gnomad_sexes: gnomAD sample sexes used in VCF export. Default is SEXES. 

    :rtype: Dict[str, Dict[str, str]]
    """
    vcf_info_dict = info_dict

    # Remove MISSING_REGION_FIELDS from info dict
    for field in MISSING_REGION_FIELDS:
        vcf_info_dict.pop(field, None)

    # Add allele-specific fields to info dict, including AS_VQSLOD and AS_culprit
    # NOTE: need to think about how to resolve AS VQSR fields to avoid having to make temp_AS_fields variable in the future
    temp_AS_fields = AS_FIELDS.copy()
    temp_AS_fields.extend(["AS_culprit", "AS_VQSLOD"])
    vcf_info_dict.update(
        add_as_info_dict(info_dict=vcf_info_dict, as_fields=temp_AS_fields)
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
        label_groups = _create_label_groups(pops=pops, sexes=sexes)
        for label_group in label_groups:
            vcf_info_dict.update(
                make_info_dict(prefix=subset, pop_names=pops, label_groups=label_group,)
            )

        faf_label_groups = _create_label_groups(pops=faf_pops, sexes=sexes)
        for label_group in faf_label_groups:
            vcf_info_dict.update(
                make_info_dict(
                    prefix=subset, pop_names=pops, label_groups=label_group, faf=True,
                )
            )

        vcf_info_dict.update(
            make_info_dict(
                prefix=subset,
                bin_edges=bin_edges,
                popmax=True,
                age_hist_data="|".join(str(x) for x in age_hist_data),
            )
        )

    # Add variant quality histograms to info dict
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=True))
    vcf_info_dict.update(make_hist_dict(bin_edges, adj=False))
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
    for field in RF_FIELDS:
        vcf_info_dict[field] = t["rf"][f"{field}"]
    for field in VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr"][f"{field}"]

    # Add AS annotations to info dict
    for field in AS_FIELDS:
        vcf_info_dict[field] = t["info"][f"{field}"]

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
    return vcf_info_dict


def index_globals(
    freq_meta: List[Dict[str, str]], label_groups: Dict
) -> Dict[str, int]:
    """
    Create a dictionary keyed by the specified label groupings with values describing the corresponding index of each grouping entry in the freq_meta array annotation.
    e.g., if freq_meta is [{'group': 'adj'}, {'group': 'raw'}, {'group': 'adj', 'pop': 'nfe'}], then this function will return 
    {'adj': 0, 'raw': 1, 'nfe_adj': 2}.
    :param List[Dict[str, str]] freq_meta: Ordered list containing dictionaries describing all the grouping combinations contained in the
        frequency row annotation.
    param Dict[str, List[str]] label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :return: Dictionary keyed by specified label grouping combinations, with values describing the corresponding index
        of each grouping entry in the frequency array.
    :rtype: Dict[str, int]
    """
    combos = make_label_combos(label_groups)
    index_dict = {}
    for combo in combos:
        combo_fields = combo.split("_")
        for i, v in enumerate(freq_meta):
            if set(v.values()) == set(combo_fields):
                # NOTE: this is UKBB specific
                # This is to avoid updating the dictionary when the subpop is the same as the pop
                # E.g., this is to avoid updating the index for {'group': 'adj', 'pop': 'nfe'} with the
                # index for {'group': 'adj', 'pop': 'nfe', 'subpop': 'nfe'},
                if "subpop" in v and v["subpop"] == v["pop"]:
                    continue
                index_dict.update({f"{combo}": i})
    return index_dict


def make_freq_meta_index_dict(
    freq_meta: List[str],
    pops: List[str],
    groups: List[str] = GROUPS,
    sexes: List[str] = SEXES,
) -> Dict[str, int]:
    """
    Makes a dictionary of the entries in the frequency array annotation, where keys are the grouping combinations and the values
    are the 0-based integer indices.
    :param List[str] freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation.
    :param bool gnomad: Whether to index a gnomAD sample freq_meta list.
    :param List[str] pops: List of global populations in frequency array. Used for both gnomAD and UKBB. 
        Can handle populations to unique to gnomAD/UKBB or a union of all population names.
    :param List[str] groups: Group names used to generate labels for high quality genotypes and all raw genotypes. Default is GROUPS.
    :param List[str] sexes: gnomAD sample sexes used in VCF export. Default is SEXES. 
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict[str, int]
    """
    index_dict = index_globals(freq_meta, dict(group=groups))
    index_dict.update(index_globals(freq_meta, dict(group=groups, pop=pops)))
    index_dict.update(index_globals(freq_meta, dict(group=groups, sex=sexes)))
    index_dict.update(index_globals(freq_meta, dict(group=groups, pop=pops, sex=sexes)))

    return index_dict


def make_index_dict(
    t: Union[hl.MatrixTable, hl.Table], freq_meta_str: str, pops: List[str],
) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array.
    :param Table ht: Table or MatrixTable containing freq_meta global annotation to be indexed
    :param str freq_meta: freq_meta global annotation to be indexed (freq_meta, gnomad_exomes_freq_meta, or gnomad_genomes_freq_meta)
    :param List[str] pops: List of global populations in frequency array. Used for both gnomAD and UKBB. 
        Can handle populations to unique to gnomAD/UKBB or a union of all population names.
    :param List[str] subpops: List of subpops in frequency array.
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    """
    freq_meta = hl.eval(t.globals[freq_meta_str])
    # check if indexing gnomAD data
    if "gnomad" in freq_meta_str:
        index_dict = make_freq_meta_index_dict(freq_meta, pops=pops)
    else:
        index_dict = make_freq_meta_index_dict(freq_meta, pops=pops)
    return index_dict


def unfurl_nested_annotations(
    t: Union[hl.MatrixTable, hl.Table], pops: List[str],
) -> Dict[str, hl.expr.Expression]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values.
    :param Table/MatrixTable t: Table/MatrixTable containing the nested variant annotation arrays to be unfurled.
    :param bool gnomad: Whether the annotations are from gnomAD.
    :param bool genome: Whether the annotations are from genome data (relevant only to gnomAD data).
    :param List[str] pops: List of global populations in frequency array. 
    :param List[str] subpops: List of all UKBB subpops (possible hybrid population cluster names). 
    :return: Dictionary containing variant annotations and their corresponding values.
    :rtype: Dict[str, hl.expr.Expression]
    """
    expr_dict = dict()

    # Set variables to locate necessary fields, compute freq index dicts, and compute faf index dict for UKBB
    faf = "faf"
    freq = "freq"
    faf_idx = make_index_dict(t=t, freq_meta_str="faf_meta", pops=pops)
    popmax = "popmax"
    freq_idx = make_index_dict(t=t, freq_meta_str="freq_meta", pops=pops)

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    for k, i in freq_idx.items():
        # Set combination to key
        # e.g., set entry of 'afr_adj' to combo
        combo = k
        combo_dict = {
            f"AC_{combo}": t[freq][i].AC,
            f"AN_{combo}": t[freq][i].AN,
            f"AF_{combo}": t[freq][i].AF,
            f"nhomalt_{combo}": t[freq][i].homozygote_count,
        }
        expr_dict.update(combo_dict)

    ## Unfurl FAF index dict
    for (
        k,
        i,
    ) in faf_idx.items():  # NOTE: faf annotations are all done on adj-only groupings
        entry = k.split("_")
        # Set combo to equal entry
        combo_fields = entry
        combo = k

        # NOTE: need to compute UKBB separately because UKBB no longer has faf meta bundled into faf
        combo_dict = {
            f"faf95_{combo}": hl.or_missing(
                hl.set(hl.eval(t.faf_meta[i].values())) == set(combo_fields),
                t[faf][i].faf95,
            ),
            f"faf99_{combo}": hl.or_missing(
                hl.set(hl.eval(t.faf_meta[i].values())) == set(combo_fields),
                t[faf][i].faf99,
            ),
        }
        expr_dict.update(combo_dict)

    prefix = ""

    combo_dict = {
        f"{prefix}popmax": t[popmax].pop,
        f"{prefix}AC_popmax": t[popmax].AC,
        f"{prefix}AN_popmax": t[popmax].AN,
        f"{prefix}AF_popmax": t[popmax].AF,
        f"{prefix}nhomalt_popmax": t[popmax].homozygote_count,
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

        if args.prepare_vcf_mt:
            logger.info("Starting VCF process...")
            logger.info("Reading in release HT...")
            ht = hl.read_table(release_ht_path())

            logger.info("Removing chrM...")
            ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)

            if args.test:
                logger.info("Filtering to chr20 and chrX (for tests only)...")
                # Using filter intervals to keep all the work done by get_ukbb_data
                # (removing sample with withdrawn consent/their ref blocks/variants,
                # also keeping meta col annotations)
                # Using chr20 to test a small autosome and chrX to test a sex chromosome
                # Some annotations (like FAF) are 100% missing on autosomes
                mt = hl.filter_intervals(
                    ht,
                    [hl.parse_locus_interval("chr20"), hl.parse_locus_interval("chrX")],
                )

            logger.info("Making histogram bin edges...")
            bin_edges = make_hist_bin_edges_expr(ht, prefix="")

            logger.info("Getting age hist data...")
            age_hist_data = hl.eval(ht.age_distribution)

            logger.info("Making INFO dict for VCF...")
            vcf_info_dict = populate_info_dict(
                bin_edges=bin_edges, age_hist_data=age_hist_data,
            )

            # Adjust keys to remove adj tags before exporting to VCF
            new_vcf_info_dict = {
                i.replace("_adj", ""): j for i, j in vcf_info_dict.items()
            }

            # Add non-PAR annotation
            ht = ht.annotate_rows(
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


            # Unfurl nested gnomAD  frequency annotations and add to info field
            ht = ht.annotate(
                info=ht.info.annotate(
                    **unfurl_nested_annotations(
                        ht, pops=POPS
                    )
                )
            )
            ht = mt.annotate_rows(**set_female_y_metrics_to_na(ht)) #Don't do this

            # Reformat vep annotation
            mt = mt.annotate_rows(vep=vep_struct_to_csq(mt.vep))
            mt = mt.annotate_rows(info=mt.info.annotate(vep=mt.vep))

            # Select relevant fields for VCF export
            mt = mt.select_rows("info", "filters", "rsid", "qual")
            new_vcf_info_dict.update(
                {"vep": {"Description": hl.eval(mt.vep_csq_header)}}
            )

            # Make filter dict and add field for MonoAllelic filter
            filter_dict = make_vcf_filter_dict(
                hl.eval(mt.rf_globals.rf_snv_cutoff.min_score),
                hl.eval(mt.rf_globals.rf_indel_cutoff.min_score),
                hl.eval(mt.rf_globals.inbreeding_cutoff),
            )
            filter_dict["MonoAllelic"] = {
                "Description": "Samples are all homozygous reference or all homozygous alternate for the variant"
            }
            header_dict = {
                "info": new_vcf_info_dict,
                "filter": filter_dict,
                "format": FORMAT_DICT,
            }

            logger.info("Saving header dict to pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "wb") as p:
                pickle.dump(header_dict, p, protocol=pickle.HIGHEST_PROTOCOL)
            mt.write(release_mt_path(*tranche_data), args.overwrite)

        if args.sanity_check:
            mt = hl.read_matrix_table(release_mt_path(*tranche_data))
            # NOTE: removing lowqual and star alleles here to avoid having additional failed missingness checks
            info_ht = hl.read_table(info_ht_path(data_source, freeze))
            mt = mt.filter_rows(
                (~info_ht[mt.row_key].AS_lowqual)
                & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
            )
            sanity_check_release_mt(
                mt, SUBSET_LIST, missingness_threshold=0.5, verbose=args.verbose
            )

        if args.prepare_release_vcf:

            logger.warning(
                "VCF export will densify! Make sure you have an autoscaling cluster."
            )
            if not args.per_chromosome and not args.parallelize:
                logger.error("Need to choose how to export the release VCF. Exiting...")
                sys.exit(1)

            mt = hl.read_matrix_table(release_mt_path(*tranche_data))
            logger.info("Reading header dict from pickle...")
            with hl.hadoop_open(release_header_path(*tranche_data), "rb") as p:
                header_dict = pickle.load(p)

            # Drop unnecessary histograms
            # TODO: figure out which hists we want to export and only create those for 500K
            drop_hists = (
                [x + "_n_smaller" for x in HISTS]
                + [x + "_bin_edges" for x in HISTS]
                + [x + "_n_larger" for x in HISTS if "dp_" not in x]
            )
            drop_hists.extend(
                [x + "_raw_n_smaller" for x in HISTS]
                + [x + "_raw_bin_edges" for x in HISTS]
                + [x + "_raw_n_larger" for x in HISTS if "dp_" not in x]
                + ["age_hist_hom_bin_edges", "age_hist_het_bin_edges"]
            )
            mt = mt.annotate_rows(info=mt.info.drop(*drop_hists))

            # Reformat names to remove "adj" pre-export
            # e.g, renaming "AC_adj" to "AC"
            # All unlabeled frequency information is assumed to be adj
            row_annots = list(mt.row.info)
            new_row_annots = []
            for x in row_annots:
                x = x.replace("_adj", "")
                new_row_annots.append(x)

            info_annot_mapping = dict(
                zip(new_row_annots, [mt.info[f"{x}"] for x in row_annots])
            )

            # Confirm all VCF fields and descriptions are present
            if not vcf_field_check(mt, header_dict, new_row_annots, list(mt.entry)):
                logger.error("Did not pass VCF field check.")
                return

            mt = mt.transmute_rows(info=hl.struct(**info_annot_mapping))

            # Rearrange INFO field in desired ordering
            mt = mt.annotate_rows(
                info=mt.info.select(
                    "AC",
                    "AN",
                    "AF",
                    "rf_tp_probability",
                    *mt.info.drop("AC", "AN", "AF", "rf_tp_probability"),
                )
            )

            # Export VCFs by chromosome
            if args.per_chromosome:
                ht = mt.rows().checkpoint(
                    get_checkpoint_path(*tranche_data, name="flat_vcf_ready", mt=False),
                    overwrite=args.overwrite,
                )

                rg = get_reference_genome(mt.locus)
                contigs = rg.contigs[:24]  # autosomes + X/Y
                logger.info(f"Contigs: {contigs}")

                for contig in contigs:
                    # Checked with Hail team about the fastest way to filter to a contig
                    # This method shouldn't be any faster than `filter_intervals`: the same amount of data is read in both cases
                    # `_calculate_new_partitions` might give us more parallelism downstream
                    # Decided to stick with `_calculate_new_partitions` method because it felt much faster on the 300K tranche
                    mt = hl.read_matrix_table(release_mt_path(*tranche_data))
                    mt = hl.filter_intervals(mt, [hl.parse_locus_interval(contig)])
                    intervals = mt._calculate_new_partitions(10000)
                    mt = hl.read_matrix_table(
                        release_mt_path(*tranche_data), _intervals=intervals
                    )
                    mt = mt.annotate_rows(**ht[mt.row_key])

                    logger.info("Densifying and exporting VCF...")
                    mt = hl.experimental.densify(mt)

                    logger.info("Removing low QUAL variants and * alleles...")
                    info_ht = hl.read_table(info_ht_path(data_source, freeze))
                    mt = mt.filter_rows(
                        (~info_ht[mt.row_key].AS_lowqual)
                        & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
                    )

                    logger.info("Adjusting sex ploidy...")
                    mt = adjust_sex_ploidy(
                        mt, mt.sex_karyotype, male_str="XY", female_str="XX"
                    )

                    hl.export_vcf(
                        mt.select_cols(),
                        release_vcf_path(*tranche_data, contig=contig),
                        metadata=header_dict,
                    )

            # Export sharded VCF
            if args.parallelize:

                logger.info("Densifying...")
                mt = hl.experimental.densify(mt)

                logger.info("Removing low QUAL variants and * alleles...")
                info_ht = hl.read_table(info_ht_path(data_source, freeze))
                mt = mt.filter_rows(
                    (~info_ht[mt.row_key].AS_lowqual)
                    & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
                )

                logger.info("Adjusting sex ploidy...")
                mt = adjust_sex_ploidy(
                    mt, mt.sex_karyotype, male_str="XY", female_str="XX"
                )
                mt = mt.select_cols()

                if args.n_shards:
                    mt = mt.naive_coalesce(args.n_shards)

                hl.export_vcf(
                    mt,
                    release_vcf_path(*tranche_data),
                    parallel="header_per_shard",
                    metadata=header_dict,
                )
    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(logging_path(*tranche_data))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--key_by_locus_and_alleles",
        help="Whether to key raw MT by locus and alleles.",
        action="store_true",
    )
    parser.add_argument(
        "--repartition",
        help="Repartition raw MT on read. Needs to be true for tranche 3/freeze 6/300K.",
        action="store_true",
    )
    parser.add_argument(
        "--raw_partitions",
        help="Number of desired partitions for the raw MT. Necessary only for tranche 3/freeze 6/300K. Used only if --repartition is also specified",
        default=30000,
        type=int,
    )
    parser.add_argument(
        "--test",
        help="Create release files using only chr20 and chrX for testing purposes",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_vcf_mt", help="Use release mt to create vcf mt", action="store_true"
    )
    parser.add_argument(
        "--sanity_check", help="Run sanity checks function", action="store_true"
    )
    parser.add_argument(
        "--verbose",
        help="Run sanity checks function with verbose output",
        action="store_true",
    )
    parser.add_argument(
        "--prepare_release_vcf", help="Prepare release VCF", action="store_true"
    )
    export_opts = parser.add_mutually_exclusive_group()
    export_opts.add_argument(
        "--per_chromosome",
        help="Prepare release VCFs per chromosome",
        action="store_true",
    )
    export_opts.add_argument(
        "--parallelize",
        help="Parallelize VCF export by exporting sharded VCF",
        action="store_true",
    )
    parser.add_argument(
        "--n_shards",
        help="Desired number of shards for output VCF (if --parallelize is set)",
        type=int,
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
