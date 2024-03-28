"""Script to export joint release Table to VCF."""
import argparse
import logging
from copy import deepcopy
from pprint import pprint
from typing import Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.vcf import (
    HISTS,
    SEXES,
    adjust_vcf_incompatible_types,
    create_label_groups,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
)

from gnomad_qc.v4.create_release.validate_and_export_vcf import (
    FAF_POPS,
    POPS,
    populate_subset_info_dict,
)
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("export_joint_vcf")
logger.setLevel(logging.INFO)

expr_dict = {}

REGION_FLAG_FIELDS = [
    "fail_interval_qc",
    "outside_broad_capture_region",
    "outside_ukb_capture_region",
    "outside_broad_calling_region",
    "outside_ukb_calling_region",
    "not_called_in_exomes",
    "not_called_in_genomes",
]

# VCF INFO fields to reorder.
VCF_INFO_REORDER = [
    "AC",
    "AN",
    "AF",
    "grpmax",
    "fafmax_faf95_max",
    "fafmax_faf95_max_gen_anc",
]


def unfurl_nested_structs(ht: hl.Table) -> hl.expr.StructExpression:
    """
    Unfurl nested annotations in a Table to a flat expression.

    :param ht: Table with nested annotations
    :return: StructExpression containing variant annotations and their corresponding
        expressions and updated entries
    """
    expr_dict = {}

    for data_type in ["exomes", "genomes", "joint"]:
        for annotation in ["freq", "faf"]:
            # Generating index dictionary
            annotation_meta = getattr(ht[data_type + "_globals"], annotation + "_meta")
            annotation_index_dict = make_freq_index_dict_from_meta(annotation_meta)
            annotation_idx = hl.eval(annotation_index_dict)

            logger.info(f"Unfurling {annotation} data...")
            for k, i in annotation_idx.items():
                for f in ht[data_type][annotation][0].keys():
                    key = (
                        f"{f if f != 'homozygote_count' else 'nhomalt'}_{data_type}_{annotation}_{k}"
                    )
                    expr = ht[data_type][annotation][i][f]
                    expr_dict[key] = expr

        logger.info(f"Unfurling grpmax data...")
        grpmax_idx = ht[data_type].grpmax
        grpmax_dict = {"grpmax_" + data_type: grpmax_idx.gen_anc}
        for f in [field for field in grpmax_idx._fields if field != "gen_anc"]:
            key = f"{f if f != 'homozygote_count' else 'nhomalt'}_grpmax_{data_type}"
            grpmax_dict[key] = grpmax_idx[f]
        expr_dict.update(grpmax_dict)

        logger.info(f"Unfurling fafmax data...")
        fafmax_idx = ht[data_type].fafmax
        fafmax_dict = {
            f"fafmax_{f}_{data_type}": fafmax_idx[f] for f in fafmax_idx.keys()
        }
        expr_dict.update(fafmax_dict)

        logger.info(f"Unfurling age hists...")
        age_hist_idx = ht[data_type].histograms.age_hists
        AGE_HISTS = ["age_hist_het", "age_hist_hom"]
        for hist in AGE_HISTS:
            for f in age_hist_idx[hist].keys():
                key = f"{hist}_{f}_{data_type}"
                expr = (
                    hl.delimit(age_hist_idx[hist][f], delimiter="|")
                    if "bin" in f
                    else age_hist_idx[hist][f]
                )
                expr_dict[key] = expr

        # Histograms to export are:
        # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
        # We previously dropped:
        # _n_smaller for all hists
        # _bin_edges for all hists
        # _n_larger for all hists EXCEPT DP hists
        # TODO: Are we still dropping them from joint VCF?
        for hist in HISTS:
            hist_dict = {
                f"{hist}_bin_freq_{data_type}": hl.delimit(
                    ht[data_type].histograms.qual_hists[hist].bin_freq, delimiter="|"
                ),
            }
            expr_dict.update(hist_dict)

            if "dp" in hist:
                expr_dict.update(
                    {
                        f"{hist}_n_larger_{data_type}": (
                            ht[data_type].histograms.qual_hists[hist].n_larger
                        )
                    },
                )
    logger.info(
        "Unfurling contingency table test results of exomes and genomes freq"
        " comparison..."
    )
    contingency_idx = hl.eval(ht.joint_globals.freq_index_dict)
    for k, i in contingency_idx.items():
        for f in ht.freq_comparison_stats.contingency_table_test[0].keys():
            key = f"CTT_{f}_{k}"
            expr = ht.freq_comparison_stats.contingency_table_test[i][f]
            expr_dict[key] = expr

    return hl.struct(**expr_dict)


def prepare_ht_for_export(
    ht: hl.Table,
    vcf_info_reorder: Optional[List[str]] = VCF_INFO_REORDER,
    info_fields_to_drop: Optional[List[str]] = None,
) -> Tuple[hl.Table, List[str]]:
    """
    Format validated HT for export.

    Drop downsamplings frequency stats from info, rearrange info, and make sure fields
    are VCF compatible.

    :param ht: Inpuy joint release HT.
    :param vcf_info_reorder: Order of VCF INFO fields. These will be placed in front of
        all other fields in the order specified.
    :param info_fields_to_drop: List of info fields to drop from the info struct.
    :return: Formatted HT and list rename row annotations.
    """
    ht = ht.annotate(info=unfurl_nested_structs(ht))
    ht = ht.annotate(info=ht.info.annotate(**{f: ht[f] for f in REGION_FLAG_FIELDS}))
    ht = ht.annotate(
        info=ht.info.annotate(
            CMH_p_value=ht.freq_comparison_stats.cochran_mantel_haenszel_test.p_value,
            CMH_chisq=ht.freq_comparison_stats.cochran_mantel_haenszel_test.chisq,
        ),
    )

    if vcf_info_reorder:
        logger.info("Rearranging fields to desired order...")
        ht = ht.annotate(
            info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder))
        )

    # Select relevant fields for VCF export
    ht = ht.select("info")

    if info_fields_to_drop is None:
        info_fields_to_drop = []

    # TODO: Are we dropping these fields from joint VCF?
    logger.info("Add age_histogram bin edges to info fields to drop...")
    info_fields_to_drop.extend(["age_hist_het_bin_edges", "age_hist_hom_bin_edges"])

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
    # The Table is already split so there are no annotations that need to be
    # pipe delimited.
    ht = adjust_vcf_incompatible_types(ht, pipe_delimited_annotations=[])

    return ht, new_row_annots


def populate_info_dict(
    info_fields: List[str],
    bin_edges: Dict[str, str],
    age_hist_distribution: str = None,
    subset_list: List[str] = ["exomes", "genomes", "joint"],
    pops: Dict[str, str] = POPS,
    faf_pops: Dict[str, List[str]] = FAF_POPS,
    sexes: List[str] = SEXES,
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
    :param label_delimiter: String to use as delimiter when making group label
        combinations.
    :param data_type: Data type to populate info dict for. One of "exomes" or
        "genomes". Default is "exomes".
    :return: Updated INFO dictionary for VCF export.
    """
    vcf_info_dict = {}
    vcf_info_dict = {f: vcf_info_dict[f] for f in info_fields if f in vcf_info_dict}

    for subset in subset_list:
        subset_pops = deepcopy(pops)
        if (subset == "joint") | (data_type == "genomes"):
            subset_pops.update({"ami": "Amish"})
        description_text = "" if subset == "" else f" in {subset} subset"

        vcf_info_dict.update(
            populate_subset_info_dict(
                subset=subset,
                description_text=description_text,
                data_type=subset,  # Because we are using joint release
                pops=subset_pops,
                faf_pops=faf_pops,
                sexes=sexes,
                label_delimiter=label_delimiter,
            )
        )

    if age_hist_distribution:
        age_hist_distribution = "|".join(str(x) for x in age_hist_distribution)

    # Add age histogram data to info dict.
    vcf_info_dict.update(
        make_info_dict(
            label_delimiter=label_delimiter,
            bin_edges=bin_edges,
            age_hist_distribution=age_hist_distribution,
        )
    )

    # Add variant quality histograms to info dict.
    vcf_info_dict.update(
        make_hist_dict(bin_edges, adj=True, drop_n_smaller_larger=True)
    )

    return vcf_info_dict


def make_hist_bin_edges_expr(
    ht: hl.Table,
    hists: List[str] = HISTS,
    row_name: str = "",
    prefix: str = "",
    label_delimiter: str = "_",
    include_age_hists: bool = True,
) -> Dict[str, str]:
    """
    Create dictionaries containing variant histogram annotations and their associated bin edges, formatted into a string separated by pipe delimiters.

    :param ht: Table containing histogram variant annotations.
    :param hists: List of variant histogram annotations. Default is HISTS.
    :param row_name: Name of row annotation containing histogram data. In the joint releast HT, it could be "exomes", "genomes", or "joint".
    :param prefix: Prefix text for age histogram bin edges.  Default is empty string.
    :param label_delimiter: String used as delimiter between prefix and histogram annotation.
    :param include_age_hists: Include age histogram annotations.
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values.
    """
    # Add underscore to prefix if it isn't empty
    if prefix != "":
        prefix += label_delimiter

    edges_dict = {}
    if include_age_hists:
        edges_dict.update(
            {
                f"{prefix}{call_type}": "|".join(
                    map(
                        lambda x: f"{x:.1f}",
                        ht.head(1)[row_name]
                        .histograms.age_hists[f"age_hist_{call_type}"]
                        .collect()[0]
                        .bin_edges,
                    )
                )
                for call_type in ["het", "hom"]
            }
        )

    for hist in hists:
        # Parse hists calculated on both raw and adj-filtered data
        for hist_type in [f"{prefix}raw_qual_hists", f"{prefix}qual_hists"]:
            hist_name = hist
            if "raw" in hist_type:
                hist_name = f"{prefix}{hist}_raw"

            edges_dict[hist_name] = "|".join(
                map(
                    lambda x: f"{x:.2f}" if "ab" in hist else str(int(x)),
                    ht.head(1)[row_name]
                    .histograms[hist_type][hist]
                    .collect()[0]
                    .bin_edges,
                )
            )

    return edges_dict


def prepare_vcf_header_dict(
    ht: hl.Table,
    bin_edges: Dict[str, str],
    age_hist_distribution: str,
    subset_list: List[str],
    pops: Dict[str, str],
    data_type: str = "exomes",
    joint_included: bool = False,
) -> Dict[str, Dict[str, str]]:
    """
    Prepare VCF header dictionary.

    :param ht: Input Table
    :param validated_ht: Validated HT with unfurled info fields.
    :param bin_edges: Dictionary of variant annotation histograms and their associated
        bin edges.
    :param age_hist_distribution: Pipe-delimited string of overall age histogram bin
        frequency.
    :param subset_list: List of sample subsets in dataset.
    :param pops: List of sample global genetic ancestry group names for gnomAD data type.
    :param format_dict: Dictionary describing MatrixTable entries. Used in header for
        VCF export.
    :param data_type: Data type to prepare VCF header for. One of "exomes" or "genomes".
        Default is "exomes".
    :param joint_included: Whether joint frequency data is included in the HT.
    :return: Prepared VCF header dictionary.
    """
    # subset = "" represents full dataset in VCF header construction, the
    # logic in gnomad_methods is built around this.
    subset_list.extend(["", "joint"] if joint_included else [""])
    logger.info("Making INFO dict for VCF...")
    vcf_info_dict = populate_info_dict(
        info_fields=list(ht.info),
        bin_edges=bin_edges,
        age_hist_distribution=age_hist_distribution,
        subset_list=subset_list,
        pops=pops,
        data_type=data_type,
    )

    # Adjust keys to remove adj tags before exporting to VCF.
    new_vcf_info_dict = {i.replace("_adj", ""): j for i, j in vcf_info_dict.items()}

    header_dict = {
        "info": new_vcf_info_dict,
    }

    return header_dict


def main(args):
    """Export joint Table to VCF."""
    hl.init(
        log="/export_joint_vcf.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )

    overwrite = args.overwrite
    test = args.test
    contig = args.contig

    ht = release_sites(data_type="joint").ht()

    if contig and test:
        raise ValueError(
            "Test argument cannot be used with contig argument as test filters"
            " to chr20, X, and Y."
        )

    if test:
        logger.info("Filter to PCSK9: 1:55039447-55064852 for testing...")
        ht = hl.filter_intervals(
            ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )
    if contig:
        logger.info(f"Filtering to {contig}...")
        ht = hl.filter_intervals(
            ht, [hl.parse_locus_interval(contig, reference_genome="GRCh38")]
        )
