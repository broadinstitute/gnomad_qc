import logging
from typing import Dict, List, Optional

import copy
import itertools
import hail as hl

from gnomad.assessment.sanity_checks import (
    generic_field_check,
    make_filters_sanity_check_expr,
)

from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    HISTS,
    INFO_DICT,
    REGION_FLAG_FIELDS,
    SORT_ORDER,
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
    "Number": "1",
    "Description": "Allele-specific forward/reverse read counts for strand bias tests",
}

# Add new site fields
NEW_SITE_FIELDS = [
    "monoallelic",
    "QUALapprox",
    "transmitted_singleton",
]
SITE_FIELDS.extend(NEW_SITE_FIELDS)
AS_FIELDS.extend(["AS_SB_TABLE"])

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

def make_label_combos(
    label_groups: Dict[str, List[str]], sort_order: List[str] = SORT_ORDER,
) -> List[str]:
    """
    Make combinations of all possible labels for a supplied dictionary of label groups.

    For example, if label_groups is `{"sex": ["male", "female"], "pop": ["afr", "nfe", "amr"]}`,
    this function will return `["afr_male", "afr_female", "nfe_male", "nfe_female", "amr_male", "amr_female']`

    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: list of all possible combinations of values for the supplied label groupings.
    """
    copy_label_groups = copy.deepcopy(label_groups)
    if len(copy_label_groups) == 1:
        return [item for sublist in copy_label_groups.values() for item in sublist]
    anchor_group = sorted(copy_label_groups.keys(), key=lambda x: sort_order.index(x))[
        0
    ]
    anchor_val = copy_label_groups.pop(anchor_group)
    combos = []
    for x, y in itertools.product(anchor_val, make_label_combos(copy_label_groups)):
        combos.append("{0}-{1}".format(x, y))
    return combos



POPS = ["afr", "ami", "amr", "asj", "eas", "fin", "mid", "nfe", "oth", "sas"]
SEXES = ["XX", "XY"]
TGP_POPS = [
    "cdx",
    "gih",
    "msl",
    "ibs",
    "beb",
    "itu",
    "ceu",
    "tsi",
    "fin",
    "pjl",
    "gwd",
    "pel",
    "pur",
    "acb",
    "chs",
    "mxl",
    "asw",
    "clm",
    "yri",
    "lwk",
    "gbr",
    "chb",
    "esn",
    "jpt",
    "khv",
    "stu",
]
HGDP_POPS = [
    "adygei",
    "burusho",
    "palestinian",
    "han",
    "bantukenya",
    "french",
    "bantusafrica",
    "basque",
    "tu",
    "dai",
    "makrani",
    "pima",
    "tujia",
    "mozabite",
    "mandenka",
    "surui",
    "orcadian",
    "bedouin",
    "she",
    "brahui",
    "naxi",
    "daur",
    "lahu",
    "miaozu",
    "xibo",
    "russian",
    "balochi",
    "maya",
    "sardinian",
    "yoruba",
    "colombian",
    "karitiana",
    "yizu",
    "pathan",
    "hazara",
    "cambodian",
    "kalash",
    "yakut",
    "italian",
    "hezhen",
    "mongola",
    "tuscan",
    "uygur",
    "sindhi",
    "druze",
    "japanese",
    "oroqen",
]
LABEL_GROUP_SORT_ORDER = [
    "subset",
    "downsampling",
    "popmax",
    "pop",
    "subpop",
    "sex",
    "group",
]
DOWNSAMPLINGS = {
    "10": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas", "mid"],
    "20": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas", "mid"],
    "50": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas", "mid"],
    "100": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas", "mid"],
    "158": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas", "mid"],
    "200": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"],
    "456": ["afr", "ami", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"],
    "500": ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"],
    "1000": ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"],
    "1047": ["afr", "amr", "asj", "eas", "fin", "nfe", "oth", "sas"],
    "1736": ["afr", "amr", "asj", "eas", "fin", "nfe", "sas"],
    "2000": ["afr", "amr", "eas", "fin", "nfe", "sas"],
    "2419": ["afr", "amr", "eas", "fin", "nfe", "sas"],
    "2604": ["afr", "amr", "eas", "fin", "nfe"],
    "5000": ["afr", "amr", "fin", "nfe"],
    "7647": ["afr", "amr", "nfe"],
    "10000": ["afr", "nfe"],
    "15000": ["afr", "nfe"],
    "20000": ["afr", "nfe"],
    "20744": ["afr", "nfe"],
    "25000": ["nfe"],
    "30000": ["nfe"],
    "34029": ["nfe"],
}
logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("sanity_checks")
logger.setLevel(logging.INFO)


# Resources to check raw HT upon loading
def summarize_ht(ht: hl.Table, monoallelic_check: hl.bool) -> hl.Struct:
    """
    Gets a summary of variants in a MatrixTable.
    Prints number of variants to stdout, and checks that each chromosome has variant calls
    :param MatrixTable ht: Raw MatrixTable to be checked
    :return: Struct of MatrixTable variant summary
    :rtype: Struct
    """

    var_summary = hl.summarize_variants(ht, show=False)
    logger.info(
        f"Dataset has {var_summary.n_variants} variants in {var_summary.contigs}"
    )

    # check that all contigs have variant calls
    for contig in var_summary.contigs:
        if var_summary.contigs[contig] == 0:
            logger.warning(f"{contig} has no variants called")

    if monoallelic_check:
        logger.info(
            f"There are {ht.filter(ht.info.monoallelic).count()} monoallelic sites in the dataset."
        )

    return var_summary


# Resources to check HT upon VCF export
def filters_sanity_check(ht: hl.Table) -> None:
    """
    Summarizes variants filtered under various conditions in input Table.
    Summarizes counts for:
        - Total number of variants
        - Fraction of variants removed due to:
            - Any filter
            - Inbreeding coefficient filter in combination with any other filter
            - AC0 filter in combination with any other filter
            - Random forest filtering in combination with any other filter
            - VQSR filtering in combination with any other filter
            - Only inbreeding coefficient filter
            - Only AC0 filter
            - Only RF filtering
            - Only VQSR filtering
    :param hl.Table ht: Input Table.
    :return: None
    :rtype: None
    """
    ht_explode = ht.explode(ht.filters)
    logger.info(
        f"hl.agg.counter filters: {ht_explode.aggregate(hl.agg.counter(ht_explode.filters))}"
    )
    # NOTE: in_problematic_region check will need to be updated if we get hg38 decoy
    ht = ht.annotate(
        is_filtered=hl.len(ht.filters) > 0,
        in_problematic_region=hl.any(
            lambda x: x, [ht.info.lcr, ht.info.segdup, ht.info.nonpar]
        ),
    )

    def _filter_agg_order(
        ht: hl.Table,
        group_exprs: Dict[str, hl.expr.Expression],
        n_rows: int = None,
        n_cols: int = None,
        extra_filter_checks: Optional[Dict[str, hl.expr.Expression]] = None,
    ) -> None:
        """
        Performs sanity checks to measure percentages of variants filtered under different conditions.
        :param hl.Table ht: Input Table.
        :param hl.expr.Expression group_exprs: Dictionary of expressions to group the table by.
        :param int n_rows: Number of rows to show.
        :param int n_cols: Number of columns to show.
        :return: None
        """
        # NOTE: make_filters_sanity_check_expr returns a dict with %ages of variants filtered
        ht.group_by(**group_exprs).aggregate(
            **make_filters_sanity_check_expr(ht, extra_filter_checks)
        ).order_by(hl.desc("n")).show(n_rows, n_cols)

    logger.info(
        "Checking distributions of filtered variants amongst variant filters..."
    )

    new_filters_dict = {
        "frac_vqsr": hl.agg.fraction(ht.filters.contains("AS_VQSR")),
        "frac_vqsr_only": hl.agg.fraction(
            ht.filters.contains("AS_VQSR") & (ht.filters.length() == 1)
        ),
    }
    _filter_agg_order(
        ht, {"is_filtered": ht.is_filtered}, extra_filter_checks=new_filters_dict
    )

    logger.info("Checking distributions of variant type amongst variant filters...")
    _filter_agg_order(
        ht, {"allele_type": ht.info.allele_type}, extra_filter_checks=new_filters_dict
    )

    logger.info(
        "Checking distributions of variant type and region type amongst variant filters..."
    )
    _filter_agg_order(
        ht,
        {
            "allele_type": ht.info.allele_type,
            "in_problematic_region": ht.in_problematic_region,
        },
        50,
        140,
        extra_filter_checks=new_filters_dict,
    )

    logger.info(
        "Checking distributions of variant type, region type, and number of alt alleles amongst variant filters..."
    )
    _filter_agg_order(
        ht,
        {
            "allele_type": ht.info.allele_type,
            "in_problematic_region": ht.in_problematic_region,
            "n_alt_alleles": ht.info.n_alt_alleles,
        },
        50,
        140,
        extra_filter_checks=new_filters_dict,
    )


def histograms_sanity_check(
    ht: hl.Table, verbose: bool, hists: List[str] = HISTS
) -> None:
    """
    Checks the number of variants that have nonzero values in their n_smaller and n_larger bins of quality histograms (both raw and adj).
    :param hl.Table ht: Input Table.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param List[str] hists: List of variant annotation histograms.
    :return: None
    :rtype: None
    """
    for hist in hists:
        for suffix in ["", "raw"]:
            if suffix == "raw":
                logger.info("Checking raw qual hists...")
                hist = f"{hist}_{suffix}"
            else:
                logger.info("Checking adj qual hists...")

            # Check subfield == 0
            generic_field_check(
                ht,
                cond_expr=(ht.info[f"{hist}_n_smaller"] != 0),
                check_description=f"{hist}_n_smaller == 0",
                display_fields=[f"info.{hist}_n_smaller"],
                verbose=verbose,
            )
            if hist not in [
                "dp_hist_alt",
                "dp_hist_all",
            ]:  # NOTE: DP hists can have nonzero values in n_larger bin
                generic_field_check(
                    ht,
                    cond_expr=(ht.info[f"{hist}_n_larger"] != 0),
                    check_description=f"{hist}_n_larger == 0",
                    display_fields=[f"info.{hist}_n_larger"],
                    verbose=verbose,
                )


def raw_and_adj_sanity_checks(ht: hl.Table, subsets: List[str], verbose: bool):
    """
    Performs sanity checks on raw and adj data in input Table.
    Checks that:
        - Raw AC, AN, AF are not 0
        - Adj AN is not 0 and AC and AF are not negative
        - Raw values for AC, AN, nhomalt in each sample subset are greater than or equal to their corresponding adj values
    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """
    for subfield in ["AC", "AF"]:
        # Check raw AC, AF > 0

        generic_field_check(
            ht,
            cond_expr=(ht.info[f"{subfield}-raw"] <= 0),
            check_description=f"{subfield}-raw > 0",
            display_fields=[f"info.{subfield}-raw"],
            verbose=verbose,
        )
        # Check adj AC, AF >=0
        generic_field_check(
            ht,
            cond_expr=(ht.info[f"{subfield}-adj"] < 0),
            check_description=f"{subfield}-adj >= 0",
            display_fields=[f"info.{subfield}_-dj", "filters"],
            verbose=verbose,
        )

    # Check raw AN > 0
    generic_field_check(
        ht,
        cond_expr=(ht.info["AN-raw"] <= 0),
        check_description="AN-raw > 0",
        display_fields=["info.AN-raw"],
        verbose=verbose,
    )

    # Check adj AN >= 0
    generic_field_check(
        ht,
        cond_expr=(ht.info["AN-adj"] < 0),
        check_description="AN-adj >= 0",
        display_fields=["info.AN-adj"],
        verbose=verbose,
    )
    # Check overall gnomad's raw subfields >= adj
    for subfield in ["AC", "AN", "nhomalt"]:
        generic_field_check(
            ht,
            cond_expr=(ht.info[f"{subfield}-raw"] < ht.info[f"{subfield}-adj"]),
            check_description=f"{subfield}-raw >= {subfield}-adj",
            display_fields=[f"info.{subfield}-raw", f"info.{subfield}-adj",],
            verbose=verbose,
        )

    for subset in subsets:
        for subfield in ["AC", "AN", "nhomalt"]:
            # Check AC_raw >= AC adj
            generic_field_check(
                ht,
                cond_expr=(
                    ht.info[f"{subfield}-{subset}-raw"]
                    < ht.info[f"{subfield}-{subset}-adj"]
                ),
                check_description=f"{subfield}-{subset}-raw >= {subfield}-{subset}-adj",
                display_fields=[
                    f"info.{subfield}-{subset}-raw",
                    f"info.{subfield}-{subset}-adj",
                ],
                verbose=verbose,
            )


def frequency_sanity_checks(ht: hl.Table, subsets: List[str], verbose: bool) -> None:
    """
    Performs sanity checks on frequency data in input Table.
    Checks:
        - Number of sites where gnomAD callset frequency is equal to a gnomAD subset frequency (both raw and adj)
    
    Also performs small spot checks:
        - Counts total number of sites where the gnomAD allele count annotation is defined (both raw and adj)
        
    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """
    for subset in subsets:
        for subfield in ["AC", "AN", "nhomalt"]:
            logger.info("adj checks")
            generic_field_check(
                ht,
                cond_expr=(
                    ht.info[f"{subfield}-adj"] == ht.info[f"{subfield}-{subset}-adj"]
                ),
                check_description=f"{subfield}-adj != {subfield}-{subset}-adj",
                display_fields=[
                    f"info.{subfield}-adj",
                    f"info.{subfield}-{subset}-adj",
                ],
                verbose=verbose,
            )

    freq_counts = ht.aggregate(
        hl.struct(
            total_defined_AC=hl.agg.count_where(hl.is_defined(ht.info["AC-adj"])),
            total_defined_AC_raw=hl.agg.count_where(hl.is_defined(ht.info["AC-raw"])),
        )
    )
    logger.info(f"Frequency spot check counts: {freq_counts}")


def sample_sum_check(
    ht: hl.Table,
    prefix: str,
    label_groups: Dict[str, List[str]],
    verbose: bool,
    subpop: bool = None,
    sort_order: List[str] = LABEL_GROUP_SORT_ORDER,
) -> None:
    """
    Compute afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    display results from checking the sum of the specified annotations in the terminal.

    :param ht: Table containing annotations to be summed.
    :param prefix: String indicating sample subset.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :return: None
    """
    if prefix != "":
        prefix += "-"

    label_combos = make_label_combos(label_groups)
    combo_AC = [ht.info[f"AC-{prefix}{x}"] for x in label_combos]
    combo_AN = [ht.info[f"AN-{prefix}{x}"] for x in label_combos]
    combo_nhomalt = [ht.info[f"nhomalt-{prefix}{x}"] for x in label_combos]

    group = label_groups.pop("group")[0]
    alt_groups = "_".join(
        sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
    )
    annot_dict = {
        f"sum_AC-{group}-{alt_groups}": hl.sum(combo_AC),
        f"sum_AN-{group}-{alt_groups}": hl.sum(combo_AN),
        f"sum_nhomalt-{group}-{alt_groups}": hl.sum(combo_nhomalt),
    }

    ht = ht.annotate(**annot_dict)

    for subfield in ["AC", "AN", "nhomalt"]:
        generic_field_check(
            ht,
            (
                ht.info[f"{subfield}-{prefix}{group}"]
                != ht[f"sum_{subfield}-{group}-{alt_groups}"]
            ),
            f"{subfield}-{prefix}{group} = sum({subfield}-{group}-{alt_groups})",
            [
                f"info.{subfield}-{prefix}{group}",
                f"sum_{subfield}-{group}-{alt_groups}",
            ],
            verbose,
        )


def sample_sum_sanity_checks(
    ht: hl.Table,
    subsets: List[str],
    info_metrics: List[str],
    verbose: bool,
    pops: List[str] = POPS,
    sexes: List[str] = SEXES,
) -> None:
    """
    Performs sanity checks on sample sums in input Table.
    Computes afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    displays results from checking the sum of the specified annotations in the terminal.
    Also checks that annotations for all expected sample populations are present.
    :param hl.Table ht: Input Table.
    :param List[str] subsets: List of sample subsets.
    :param List[str] info_metrics: List of metrics in info struct of input Table.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param pop_names: Dict with global population names (keys) and population descriptions (values).
    :return: None
    :rtype: None
    """
    # Add "" for sum checks on entire callset
    subsets.append("")
    # Perform sample sum checks per subset
    for subset in subsets:
        pop_names = pops
        if subset == "hgdp":
            pop_names = HGDP_POPS
        if subset == "tgp":
            pop_names = TGP_POPS

        sample_sum_check(ht, subset, dict(group=["adj"], pop=pop_names), verbose)
        sample_sum_check(ht, subset, dict(group=["adj"], sex=sexes), verbose)
        sample_sum_check(
            ht,
            subset,
            dict(group=["adj"], pop=list(set(pop_names)), sex=sexes),
            verbose,
        )


def sex_chr_sanity_checks(
    ht: hl.Table, info_metrics: List[str], contigs: List[str], verbose: bool, reference_genome
) -> None:
    """
    Performs sanity checks for annotations on the sex chromosomes.
    Checks:
        - That metrics for chrY variants in female samples are NA and not 0
        - That nhomalt counts are equal to female nhomalt counts for all non-PAR chrX variants
    :param hl.Table ht: Input Table.
    :param List[str] info_metrics: List of metrics in info struct of input Table.
    :param List[str] contigs: List of contigs present in input Table.
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """

    female_metrics = [x for x in info_metrics if "-female" in x or "-XX" in x]

    if "chrY" in contigs:
        logger.info("Check values of female metrics for Y variants are NA:")
        ht_y = hl.filter_intervals(ht, [hl.parse_locus_interval("chrY", reference_genome=reference_genome)])
        metrics_values = {}
        for metric in female_metrics:
            metrics_values[metric] = hl.agg.any(hl.is_defined(ht_y.info[metric]))
        output = ht_y.aggregate(hl.struct(**metrics_values))
        for metric, value in dict(output).items():
            if value:
                values_found = ht_y.aggregate(
                    hl.agg.filter(
                        hl.is_defined(ht_y.info[metric]),
                        hl.agg.take(ht_y.info[metric], 1),
                    )
                )
                logger.info(
                    f"FAILED {metric} = {None} check for Y variants. Values found: {values_found}"
                )
            else:
                logger.info(f"PASSED {metric} = {None} check for Y variants")

    ht_x = hl.filter_intervals(ht, [hl.parse_locus_interval("chrX", reference_genome=reference_genome)])
    ht_xnonpar = ht_x.filter(ht_x.locus.in_x_nonpar())
    n = ht_xnonpar.count()
    logger.info(f"Found {n} X nonpar sites")

    logger.info("Check (nhomalt == nhomalt_female) for X nonpar variants:")
    female_metrics = [x for x in female_metrics if "nhomalt" in x]
    for metric in female_metrics:
        standard_field = metric.replace("-female", "").replace("-XX", "")
        generic_field_check(
            ht_xnonpar,
            (ht_xnonpar.info[f"{metric}"] != ht_xnonpar.info[f"{standard_field}"]),
            f"{metric} == {standard_field}",
            [f"info.{metric}", f"info.{standard_field}"],
            verbose,
        )


def missingness_sanity_checks(
    ht: hl.Table,
    info_metrics: List[str],
    non_info_metrics: List[str],
    n_sites: int,
    missingness_threshold: float,
) -> None:
    """
    Checks amount of missingness in all row annotations.
    Prints metric to terminal if more than missingness_threshold% of annotations for that metric are missing.
    :param hl.Table ht: Input Table.
    :param List[str] info_metrics: List of metrics in info struct of input Table.
    :param List[str] non_info_metrics: List of row annotations minus info struct from input Table.
    :param int n_sites: Number of sites in input Table.
    :param float missingness_threshold: Upper cutoff for allowed amount of missingness.
    :return: None
    :rtype: None
    """
    import math

    ht = ht.checkpoint("gs://gnomad-tmp/v3.1.1_sanity_check_missingness.ht", overwrite=True)

    logger.info(
        f"Missingness threshold (upper cutoff for what is allowed for missingness checks): {missingness_threshold}"
    )
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht[x])) / n_sites

    num_metrics = len(metrics_frac_missing)
    logger.info(f"There are {num_metrics} metrics to calculate missingness on...")

    if num_metrics > 1000:
        num_split = math.ceil(num_metrics / 1000.0)
        num_factor = num_metrics / num_split
        for i in range(int(num_split)):
            metrics_frac_missing_split = list(metrics_frac_missing.items())[int(i * num_factor):int((i + 1) * num_factor)]
            output = ht.aggregate(hl.struct(**dict(metrics_frac_missing_split)))

            n_fail = 0
            for metric, value in dict(output).items():
                message = f"missingness check for {metric}: {100 * value}% missing"
                if value > missingness_threshold:
                    logger.info(f"FAILED {message}")
                    n_fail += 1
                else:
                    logger.info(f"Passed {message}")

    logger.info(f"{n_fail} missing metrics checks failed")


def vcf_field_check(
    ht: hl.Table,
    header_dict: Dict[str, Dict[str, Dict[str, str]]],
    row_annotations: List[str],
    hists: List[str] = HISTS,
) -> bool:
    """
    Checks that all VCF fields and descriptions are present in input Table and VCF header dictionary.

    :param hl.MatrixTable ht: Input Table to be exported to VCF.
    :param Dict[str, Dict[str, Dict[str, str]]] header_dict: VCF header dictionary.
    :param List[str] row_annotations: List of row annotations in MatrixTable.
    :param List[str] hists: List of variant histogram annotations. Default is HISTS.
    :return: Bool with whether all expected fields and descriptions are present.
    :rtype: bool
    """
    # Confirm all VCF fields/descriptions are present before exporting
    hist_fields = []
    for hist in hists:
        hist_fields.extend(
            [
                f"{hist}_bin_freq",
                f"{hist}_n_smaller",
                f"{hist}_n_larger",
                f"{hist}_raw_bin_freq",
                f"{hist}_raw_n_smaller",
                f"{hist}_raw_n_larger",
            ]
        )

    missing_fields = []
    missing_descriptions = []
    for item in ["info", "filter"]:
        if item == "info":
            annots = row_annotations
        else:
            annot_ht = ht.explode(ht.filters)
            annots = list(annot_ht.aggregate(hl.agg.collect_as_set(annot_ht.filters)))

        temp_missing_fields = []
        temp_missing_descriptions = []
        for field in annots:
            try:
                description = header_dict[item][field]
                if len(description) == 0:
                    logger.warning(
                        f"{field} in HT info field has empty description in VCF header!"
                    )
                    temp_missing_descriptions.append(field)
            except KeyError:
                logger.warning(
                    f"{field} in HT info field does not exist in VCF header!"
                )
                # NOTE: some hists are not exported, so ignoring here
                # END entry is also not exported (removed during densify)
                if (field not in hist_fields) and (field != "END"):
                    temp_missing_fields.append(field)

        missing_fields.extend(temp_missing_fields)
        missing_descriptions.extend(temp_missing_descriptions)

    if len(missing_fields) != 0 or len(missing_descriptions) != 0:
        logger.error(
            "Some fields are either missing or missing descriptions in the VCF header! Please reconcile."
        )
        logger.error(f"Missing fields: {missing_fields}")
        logger.error(f"Missing descriptions: {missing_descriptions}")
        return False

    logger.info("Passed VCF fields check!")
    return True


def sanity_check_release_ht(
    ht: hl.Table,
    subsets: List[str],
    missingness_threshold: float = 0.5,
    verbose: bool = True,
    reference_genome="GRCh38",
) -> None:
    """
    Perform a battery of sanity checks on a specified group of subsets in a MatrixTable containing variant annotations.
    Includes:
    - Summaries of % filter status for different partitions of variants
    - Histogram outlier bin checks
    - Checks on AC, AN, and AF annotations
    - Checks that subgroup annotation values add up to the supergroup annotation values
    - Checks on sex-chromosome annotations; and summaries of % missingness in variant annotations
    :param MatrixTable ht: MatrixTable containing variant annotations to check.
    :param List[str] subsets: List of subsets to be checked.
    :param float missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.5
    :param bool verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated.
    :return: None (terminal display of results from the battery of sanity checks).
    :rtype: None
    """

    # Perform basic checks -- number of variants, number of contigs, number of samples
    """logger.info("BASIC SUMMARY OF INPUT TABLE:")
    summarize_ht(ht, monoallelic_check=True)

    logger.info("VARIANT FILTER SUMMARIES:")
    filters_sanity_check(ht)

    logger.info("HISTOGRAM CHECKS:")
    histograms_sanity_check(ht, verbose=verbose)

    logger.info("RAW AND ADJ CHECKS:")
    raw_and_adj_sanity_checks(ht, subsets, verbose)

    logger.info("FREQUENCY CHECKS:")
    frequency_sanity_checks(ht, subsets, verbose)
    """
    # Pull row annotations from HT
    info_metrics = list(ht.row.info)
    non_info_metrics = list(ht.row)
    non_info_metrics.remove("info")

    """logger.info("SAMPLE SUM CHECKS:")
    sample_sum_sanity_checks(ht, subsets, info_metrics, verbose)

    logger.info("SEX CHROMOSOME ANNOTATION CHECKS:")
    contigs = ht.aggregate(hl.agg.collect_as_set(ht.locus.contig))
    sex_chr_sanity_checks(ht, info_metrics, contigs, verbose, reference_genome=reference_genome)
    """
    logger.info("MISSINGNESS CHECKS:")
    n_sites = ht.count()
    missingness_sanity_checks(
        ht, info_metrics, non_info_metrics, n_sites, missingness_threshold
    )
