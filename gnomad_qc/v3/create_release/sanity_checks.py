import logging
from typing import Dict, List, Optional, Union

import hail as hl
from gnomad.assessment.sanity_checks import (
    generic_field_check,
    make_filters_sanity_check_expr,
)
from gnomad.resources.grch38.gnomad import HGDP_POPS, KG_POPS, POPS, SEXES, SUBSETS
from gnomad.utils.vcf import HISTS, SORT_ORDER, make_label_combos

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


def summarize_t(
    t: Union[hl.MatrixTable, hl.Table], monoallelic_check: hl.bool = False
) -> hl.Struct:
    """
    Get summary of variants in a MatrixTable or Table.

    Print the number of variants to stdout and check that each chromosome has variant calls.

    :param t: Input MatrixTable or Table to be checked.
    :param monoallelic_check: Log how many monoallelic sites are in the Table; requires a monoallelic annotation within an info struct.
    :rtype: Struct
    """
    if isinstance(t, hl.MatrixTable):
        logger.info(f"Dataset has {t.count_cols()} samples.")
        t = t.rows()

    var_summary = hl.summarize_variants(t, show=False)
    logger.info(
        f"Dataset has {var_summary.n_variants} variants distributed across the following contigs: {var_summary.contigs}"
    )

    # check that all contigs have variant calls
    for contig in var_summary.contigs:
        if var_summary.contigs[contig] == 0:
            logger.warning(f"{contig} has no variants called")

    if monoallelic_check:
        mono_sites = t.filter(t.info.monoallelic).count()
        logger.info(f"There are {mono_sites} monoallelic sites in the dataset.")

    return var_summary


# Resources to check HT/MT upon VCF export
def filters_sanity_check(t: Union[hl.MatrixTable, hl.Table]) -> None:
    """
    Summarize variants filtered under various conditions in input MatrixTable or Table.

    Summarize counts for:
        - Total number of variants
        - Fraction of variants removed due to:
            - Any filter
            - Inbreeding coefficient filter in combination with any other filter
            - AC0 filter in combination with any other filter
            - VQSR filtering in combination with any other filter
            - Only inbreeding coefficient filter
            - Only AC0 filter
            - Only RF filtering
            - Only VQSR filtering

    :param t: Input MatrixTable or Table to be checked.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    filters = t.aggregate(hl.agg.counter(t.filters))
    logger.info(f"hl.agg.counter filters: {filters}")

    filtered_expr = hl.len(t.filters) > 0
    problematic_region_expr = hl.any(
        lambda x: x,
        [
            t.info.lcr,
            t.info.segdup,
            t.info.nonpar,
        ],  # NOTE: in_problematic_region check will need to be updated if we get hg38 decoy
    )

    t = t.annotate(
        is_filtered=filtered_expr, in_problematic_region=problematic_region_expr
    )

    def _filter_agg_order(
        t: Union[hl.MatrixTable, hl.Table],
        group_exprs: Dict[str, hl.expr.Expression],
        n_rows: int = None,
        n_cols: int = None,
        extra_filter_checks: Optional[Dict[str, hl.expr.Expression]] = None,
    ) -> None:
        """
        Perform sanity checks to measure percentages of variants filtered under different conditions.

        :param t: Input MatrixTable or Table.
        :param group_exprs: Dictionary of expressions to group the Table by.
        :param n_rows: Number of rows to show.
        :param n_cols: Number of columns to show.
        :return: None
        """
        t = t.rows() if isinstance(t, hl.MatrixTable) else t
        # NOTE: make_filters_sanity_check_expr returns a dict with %ages of variants filtered
        t.group_by(**group_exprs).aggregate(
            **make_filters_sanity_check_expr(t, extra_filter_checks)
        ).order_by(hl.desc("n")).show(n_rows, n_cols)

    logger.info(
        "Checking distributions of filtered variants amongst variant filters..."
    )

    new_filters_dict = {
        "frac_vqsr": hl.agg.fraction(t.filters.contains("AS_VQSR")),
        "frac_vqsr_only": hl.agg.fraction(
            t.filters.contains("AS_VQSR") & (t.filters.length() == 1)
        ),
    }
    _filter_agg_order(
        t, {"is_filtered": t.is_filtered}, extra_filter_checks=new_filters_dict
    )

    logger.info("Checking distributions of variant type amongst variant filters...")
    _filter_agg_order(
        t, {"allele_type": t.info.allele_type}, extra_filter_checks=new_filters_dict
    )

    logger.info(
        "Checking distributions of variant type and region type amongst variant filters..."
    )
    _filter_agg_order(
        t,
        {
            "allele_type": t.info.allele_type,
            "in_problematic_region": t.in_problematic_region,
        },
        50,
        140,
        extra_filter_checks=new_filters_dict,
    )

    logger.info(
        "Checking distributions of variant type, region type, and number of alt alleles amongst variant filters..."
    )
    _filter_agg_order(
        t,
        {
            "allele_type": t.info.allele_type,
            "in_problematic_region": t.in_problematic_region,
            "n_alt_alleles": t.info.n_alt_alleles,
        },
        50,
        140,
        extra_filter_checks=new_filters_dict,
    )


def histograms_sanity_check(
    t: Union[hl.MatrixTable, hl.Table], verbose: bool, hists: List[str] = HISTS
) -> None:
    """
    Check that variants have nonzero values in their n_smaller and, with the exeption of DP hist, n_larger bins of quality histograms (both raw and adj).

    All n_smaller and n_larger annotations must be within an info struct annotation. 

    :param t: Input MatrixTable or Table.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param hists: List of variant annotation histograms.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    for hist in hists:
        for suffix in ["", "raw"]:
            if suffix == "raw":
                logger.info("Checking raw qual hists...")
                hist = f"{hist}_{suffix}"
            else:
                logger.info("Checking adj qual hists...")

            # Check subfield == 0
            generic_field_check(
                t,
                cond_expr=(t.info[f"{hist}_n_smaller"] != 0),
                check_description=f"{hist}_n_smaller == 0",
                display_fields=[f"info.{hist}_n_smaller"],
                verbose=verbose,
            )
            if hist not in [
                "dp_hist_alt",
                "dp_hist_all",
            ]:  # NOTE: DP hists can have nonzero values in n_larger bin
                generic_field_check(
                    t,
                    cond_expr=(t.info[f"{hist}_n_larger"] != 0),
                    check_description=f"{hist}_n_larger == 0",
                    display_fields=[f"info.{hist}_n_larger"],
                    verbose=verbose,
                )


def raw_and_adj_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    verbose: bool,
    delimiter: str = "-",
) -> None:
    """
    Perform sanity checks on raw and adj data in input Table/MatrixTable.

    Check that:
        - Raw AC, AN, AF are not 0
        - Adj AN is not 0 and AC and AF are not negative
        - Raw values for AC, AN, nhomalt in each sample subset are greater than or equal to their corresponding adj values
    Raw and adj call stat annotations must be in an info struct annotation on the Table/MatrixTable. 

    :param t: Input MatrixTable or Table to check.
    :param subsets: List of sample subsets.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param delimiter: String to use as delimiter when making group label combinations.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    for subfield in ["AC", "AF"]:
        field = f"{subfield}{delimiter}"
        # Check raw AC, AF > 0
        generic_field_check(
            t,
            cond_expr=(t.info[f"{field}raw"] <= 0),
            check_description=f"{field}raw > 0",
            display_fields=[f"info.{field}raw"],
            verbose=verbose,
        )
        # Check adj AC, AF >=0
        generic_field_check(
            t,
            cond_expr=(t.info[f"{field}adj"] < 0),
            check_description=f"{field}adj >= 0",
            display_fields=[f"info.{field}adj", "filters"],
            verbose=verbose,
        )

    # Check raw AN > 0
    an_raw_field = f"AN{delimiter}raw"
    generic_field_check(
        t,
        cond_expr=(t.info[an_raw_field] <= 0),
        check_description=f"{an_raw_field} > 0",
        display_fields=[f"info.{an_raw_field}"],
        verbose=verbose,
    )

    an_adj_field = f"AN{delimiter}adj"
    # Check adj AN >= 0
    generic_field_check(
        t,
        cond_expr=(t.info[an_adj_field] < 0),
        check_description=f"{an_adj_field} >= 0",
        display_fields=[f"info.{an_adj_field}"],
        verbose=verbose,
    )

    # Check overall raw subfields >= adj
    for subfield in ["AC", "AN", "nhomalt"]:
        field = f"{subfield}{delimiter}"
        generic_field_check(
            t,
            cond_expr=(t.info[f"{field}raw"] < t.info[f"{field}adj"]),
            check_description=f"{field}raw >= {field}adj",
            display_fields=[f"info.{field}raw", f"info.{field}adj",],
            verbose=verbose,
        )

    for subset in subsets:
        for subfield in ["AC", "AN", "nhomalt"]:
            # Check AC_raw >= AC adj
            field = f"{subfield}{delimiter}{subset}{delimiter}"
            generic_field_check(
                t,
                cond_expr=(t.info[f"{field}raw"] < t.info[f"{field}adj"]),
                check_description=f"{field}raw >= {field}adj",
                display_fields=[f"info.{field}raw", f"info.{field}adj",],
                verbose=verbose,
            )


def frequency_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    verbose: bool,
    show_percent_sites: bool = True,
    delimiter: str = "-",
) -> None:
    """
    Perform sanity checks on frequency data in input Table.

    Check:
        - Number of sites where callset frequency is equal to a subset frequency (raw and adj)
        - Total number of sites where the allele count annotation is defined (raw and adj)
        
    :param t: Input MatrixTable or Table.
    :param subsets: List of sample subsets.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param show_percent_sites: If true, show the percentage and count of overall sites that fail; if False, only show the number of sites that fail.
    :param delimiter: String to use as delimiter when making group label combinations.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    for subset in subsets:
        if subset:
            subset += delimiter
            for subfield in ["AC", "AN", "nhomalt"]:
                for group in ["adj", "raw"]:
                    logger.info(f"Comparing subset {group} frequencies to entire callset")
                    subfield_label = f"{subfield}{delimiter}{group}"
                    subfield_subset_label = f"{subfield}{delimiter}{subset}{group}"

                    generic_field_check(
                        t,
                        cond_expr=(t.info[subfield_label] == t.info[subfield_subset_label]),
                        check_description=f"{subfield_label} != {subfield_subset_label}",
                        display_fields=[
                            f"info.{subfield_label}",
                            f"info.{subfield_subset_label}",
                        ],
                        verbose=verbose,
                        show_percent_sites=show_percent_sites,
                    )
            

    freq_counts = t.aggregate(
        hl.struct(
            total_defined_AC=hl.agg.count_where(hl.is_defined(t.info["AC-adj"])),
            total_defined_AC_raw=hl.agg.count_where(hl.is_defined(t.info["AC-raw"])),
        )
    )
    logger.info(f"Frequency spot check counts: {freq_counts}")


def sample_sum_check(
    t: Union[hl.MatrixTable, hl.Table],
    subset: str,
    label_groups: Dict[str, List[str]],
    verbose: bool,
    subpop: bool = None,
    sort_order: List[str] = SORT_ORDER,
    delimiter: str = "-",
) -> None:
    """
    Compute the sum of call stats annotations for a specified group of annotations, compare to the annotated version, and display the result in stdout.

    :param t: Input MatrixTable or Table containing annotations to be summed.
    :param subset: String indicating sample subset.
    :param label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"]).
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked.
    :param sort_order: List containing order to sort label group combinations. Default is SORT_ORDER.
    :param delimiter: String to use as delimiter when making group label combinations.
    :return: None
    """

    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    if subset != "":
        subset += delimiter

    label_combos = make_label_combos(label_groups, label_delimiter=delimiter)
    combo_AC = [t.info[f"AC-{subset}{x}"] for x in label_combos]
    combo_AN = [t.info[f"AN-{subset}{x}"] for x in label_combos]
    combo_nhomalt = [t.info[f"nhomalt-{subset}{x}"] for x in label_combos]

    group = label_groups.pop("group")[0]
    alt_groups = delimiter.join(
        sorted(label_groups.keys(), key=lambda x: sort_order.index(x))
    )
    group_expr = f"{subset}{group}{delimiter}{alt_groups}"
    annot_dict = {
        f"sum_AC-{group_expr}": hl.sum(combo_AC),
        f"sum_AN-{group_expr}": hl.sum(combo_AN),
        f"sum_nhomalt-{group_expr}": hl.sum(combo_nhomalt),
    }

    t = t.annotate(**annot_dict)

    for subfield in ["AC", "AN", "nhomalt"]:

        group_expr = f"{subfield}{delimiter}{subset}{group}"
        alt_groups_expr = f"{subfield}{delimiter}{subset}{group}{delimiter}{alt_groups}"

        generic_field_check(
            t,
            (t.info[f"{group_expr}"] != t[f"sum_{alt_groups_expr}"]),
            f"{group_expr} = sum({alt_groups_expr})",
            [f"info.{group_expr}", f"sum_{alt_groups_expr}",],
            verbose,
        )


def sample_sum_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    info_metrics: List[str],
    verbose: bool,
    pops: List[str] = POPS,
    sexes: List[str] = SEXES,
) -> None:
    """
    Compute the sum of annotations for a specified group of annotations, compare to the annotated version and display the results in stdout.
    Also check that annotations for all expected sample populations are present.

    :param t: Input MatrixTable or Table.
    :param subsets: List of sample subsets.
    :param info_metrics: List of metrics in info struct of input Table.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :param pop_names: Dict with global population names (keys) and population descriptions (values).
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t
    
    # Add an empty string for sum checks on entire callset
    subsets.append("")

    # Perform sample sum checks per subset
    for subset in subsets:
        pop_names = pops
        if subset == "hgdp":
            pop_names = HGDP_POPS
        if subset == "tgp":
            pop_names = KG_POPS

        sample_sum_check(t, subset, dict(group=["adj"], pop=pop_names), verbose)
        sample_sum_check(t, subset, dict(group=["adj"], sex=sexes), verbose)
        sample_sum_check(
            t,
            subset,
            dict(group=["adj"], pop=list(set(pop_names)), sex=sexes),
            verbose,
        )


def sex_chr_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    info_metrics: List[str],
    contigs: List[str],
    verbose: bool,
) -> None:
    """
    Performs sanity checks for annotations on the sex chromosomes.
    Check:
        - That metrics for chrY variants in female samples are NA and not 0
        - That nhomalt counts are equal to female nhomalt counts for all non-PAR chrX variants

    :param t: Input MatrixTable or Table.
    :param info_metrics: List of metrics in info struct of input Table.
    :param contigs: List of contigs present in input Table.
    :param verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    female_metrics = [x for x in info_metrics if "-female" in x or "-XX" in x]

    if "chrY" in contigs:
        logger.info("Check values of female metrics for Y variants are NA:")
        t_y = hl.filter_intervals(t, [hl.parse_locus_interval("chrY")])
        metrics_values = {}
        for metric in female_metrics:
            metrics_values[metric] = hl.agg.any(hl.is_defined(t_y.info[metric]))
        output = t_y.aggregate(hl.struct(**metrics_values))
        for metric, value in dict(output).items():
            if value:
                values_found = t_y.aggregate(
                    hl.agg.filter(
                        hl.is_defined(t_y.info[metric]),
                        hl.agg.take(t_y.info[metric], 1),
                    )
                )
                logger.info(
                    f"FAILED {metric} = {None} check for Y variants. Values found: {values_found}"
                )
            else:
                logger.info(f"PASSED {metric} = {None} check for Y variants")

    t_x = hl.filter_intervals(t, [hl.parse_locus_interval("chrX")])
    t_xnonpar = t_x.filter(t_x.locus.in_x_nonpar())
    n = t_xnonpar.count()
    logger.info(f"Found {n} X nonpar sites")

    logger.info("Check (nhomalt == nhomalt_female) for X nonpar variants:")
    female_metrics = [x for x in female_metrics if "nhomalt" in x]
    for metric in female_metrics:
        standard_field = metric.replace("-female", "").replace("-XX", "")
        generic_field_check(
            t_xnonpar,
            (t_xnonpar.info[f"{metric}"] != t_xnonpar.info[f"{standard_field}"]),
            f"{metric} == {standard_field}",
            [f"info.{metric}", f"info.{standard_field}"],
            verbose,
        )


def missingness_sanity_checks(
    t: Union[hl.MatrixTable, hl.Table],
    info_metrics: List[str],
    non_info_metrics: List[str],
    n_sites: int,
    missingness_threshold: float,
) -> None:
    """
    Check amount of missingness in all row annotations.
    Print metric to sdout if the metric annotations missingness exceeds the missingness_threshold.

    :param t: Input MatrixTable or Table.
    :param info_metrics: List of metrics in info struct of input Table.
    :param non_info_metrics: List of row annotations minus info struct from input Table.
    :param n_sites: Number of sites in input Table.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness.
    :return: None
    :rtype: None
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

    logger.info(
        f"Missingness threshold (upper cutoff for what is allowed for missingness checks): {missingness_threshold}"
    )
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(t.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(t[x])) / n_sites
    output = t.aggregate(hl.struct(**metrics_frac_missing))

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
    t: Union[hl.MatrixTable, hl.Table],
    header_dict: Dict[str, Dict[str, Dict[str, str]]],
    row_annotations: List[str],
    hists: List[str] = HISTS,
) -> bool:
    """
    Check that all VCF fields and descriptions are present in input Table and VCF header dictionary.

    :param t: Input MatrixTable or Tableto be exported to VCF.
    :param header_dict: VCF header dictionary.
    :param row_annotations: List of row annotations in MatrixTable.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :return: Bool with whether all expected fields and descriptions are present.
    :rtype: bool
    """
    t = t.rows() if isinstance(t, hl.MatrixTable) else t

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
            annot_t = t.explode(t.filters)
            annots = list(annot_t.aggregate(hl.agg.collect_as_set(annot_t.filters)))

        temp_missing_fields = []
        temp_missing_descriptions = []
        for field in annots:
            try:
                description = header_dict[item][field]
                if len(description) == 0:
                    logger.warning(
                        f"{field} in T info field has empty description in VCF header!"
                    )
                    temp_missing_descriptions.append(field)
            except KeyError:
                logger.warning(f"{field} in T info field does not exist in VCF header!")
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


def sanity_check_release_t(
    t: Union[hl.MatrixTable, hl.Table],
    subsets: List[str],
    missingness_threshold: float = 0.5,
    verbose: bool = True,
) -> None:
    """
    Perform a battery of sanity checks on a specified group of subsets in a MatrixTable containing variant annotations.
    Includes:
    - Summaries of % filter status for different partitions of variants
    - Histogram outlier bin checks
    - Checks on AC, AN, and AF annotations
    - Checks that subgroup annotation values add up to the supergroup annotation values
    - Checks on sex-chromosome annotations; and summaries of % missingness in variant annotations

    :param t: Input MatrixTable or Table containing variant annotations to check.
    :param subsets: List of subsets to be checked.
    :param missingness_threshold: Upper cutoff for allowed amount of missingness. Default is 0.5
    :param verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated.
    :return: None (terminal display of results from the battery of sanity checks).
    :rtype: None
    """

    # Perform basic checks -- number of variants, number of contigs, number of samples
    logger.info("BASIC SUMMARY OF INPUT TABLE:")
    summarize_t(t, monoallelic_check=True)

    logger.info("VARIANT FILTER SUMMARIES:")
    filters_sanity_check(t)

    logger.info("HISTOGRAM CHECKS:")
    histograms_sanity_check(t, verbose=verbose)

    logger.info("RAW AND ADJ CHECKS:")
    raw_and_adj_sanity_checks(t, subsets, verbose)

    logger.info("FREQUENCY CHECKS:")
    frequency_sanity_checks(t, subsets, verbose)

    # Pull row annotations from HT
    info_metrics = list(t.row.info)
    non_info_metrics = list(t.row)
    non_info_metrics.remove("info")

    logger.info("SAMPLE SUM CHECKS:")
    sample_sum_sanity_checks(t, subsets, info_metrics, verbose)

    logger.info("SEX CHROMOSOME ANNOTATION CHECKS:")
    contigs = t.aggregate(hl.agg.collect_as_set(t.locus.contig))
    sex_chr_sanity_checks(t, info_metrics, contigs, verbose)

    logger.info("MISSINGNESS CHECKS:")
    n_sites = t.count()
    missingness_sanity_checks(
        t, info_metrics, non_info_metrics, n_sites, missingness_threshold
    )
