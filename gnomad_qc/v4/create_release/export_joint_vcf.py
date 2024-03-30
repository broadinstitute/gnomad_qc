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
    REGION_FLAGS_JOINT_INFO_DICT,
    SEXES,
    adjust_vcf_incompatible_types,
    create_label_groups,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
)

from gnomad_qc.v4.create_release.validate_and_export_vcf import FAF_POPS, POPS
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


def prepare_info_fields(ht: hl.Table) -> hl.expr.StructExpression:
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

            logger.info(f"Unfurling {annotation} data from {data_type}...")
            for k, i in annotation_idx.items():
                for f in ht[data_type][annotation][0].keys():
                    key = (
                        f"{f if f != 'homozygote_count' else 'nhomalt'}_{data_type}_{k}"
                    )
                    expr = ht[data_type][annotation][i][f]
                    expr_dict[key] = expr

        logger.info(f"Unfurling grpmax data from {data_type}...")
        grpmax_idx = ht[data_type].grpmax
        grpmax_dict = {"grpmax_" + data_type: grpmax_idx.gen_anc}
        for f in [field for field in grpmax_idx._fields if field != "gen_anc"]:
            key = f"{f if f != 'homozygote_count' else 'nhomalt'}_grpmax_{data_type}"
            grpmax_dict[key] = grpmax_idx[f]
        expr_dict.update(grpmax_dict)

        logger.info(f"Unfurling fafmax data from {data_type}...")
        fafmax_idx = ht[data_type].fafmax
        fafmax_dict = {
            f"fafmax_{f}_{data_type}": fafmax_idx[f] for f in fafmax_idx.keys()
        }
        expr_dict.update(fafmax_dict)

        logger.info(f"Unfurling age hists from {data_type}...")
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

        logger.info(f"Unfurling variant quality histograms from {data_type}...")
        # Histograms to export are:
        # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
        # We previously dropped:
        # _n_smaller for all hists
        # _bin_edges for all hists
        # _n_larger for all hists EXCEPT DP hists
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

    logger.info("Unfurling region flags...")
    for f in REGION_FLAG_FIELDS:
        expr_dict[f] = ht.region_flags[f]

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

    logger.info(
        "Unfurling Cochran-Mantel-Haenszel test results of exomes and genomes freq"
        " comparison..."
    )
    expr_dict["CMH_chisq"] = ht.freq_comparison_stats.cochran_mantel_haenszel_test.chisq
    expr_dict["CMH_p_value"] = (
        ht.freq_comparison_stats.cochran_mantel_haenszel_test.p_value
    )
    return hl.struct(**expr_dict)


def prepare_vcf_header_dict(
    ht: hl.Table,
    pops: Dict[str, str] = POPS,
    sexes: List[str] = SEXES,
) -> Dict[str, Dict[str, str]]:
    """
    Prepare VCF header dictionary.

    :param ht: Input Table
    :return: Prepared VCF header dictionary
    """
    vcf_info_dict = {}
    vcf_info_dict.update(REGION_FLAGS_JOINT_INFO_DICT)
    for data_type in ["exomes", "genomes", "joint"]:
        description_text = f" in {data_type} subset"  # TODO: do we call it subset?
        logger.info("Preparing freq VCF header for %s subset...", data_type)
        if data_type == "exomes":
            pops_to_use = pops.copy()
        else:
            pops_to_use = pops.copy()
            pops_to_use.update({"ami": "Amish"})
        freq_label_groups = create_label_groups(pops=pops_to_use, sexes=sexes)
        for label_group in freq_label_groups:
            vcf_info_dict.update(
                make_info_dict(
                    prefix=data_type,
                    prefix_before_metric=False,
                    pop_names=pops_to_use,
                    label_groups=label_group,
                    description_text=description_text,
                    callstats=True,
                )
            )
        logger.info("Preparing faf VCF header for %s subset...", data_type)
        if data_type == "exomes" or data_type == "joint":
            faf_pops_version = "v4"
        else:
            faf_pops_version = "v3"
        faf_pops = {pop: POP_NAMES[pop] for pop in FAF_POPS[faf_pops_version]}
        if data_type == "exomes" or data_type == "genomes":
            # There are no sex-specific FAFs in the joint dataset.
            faf_label_groups = create_label_groups(
                pops=faf_pops, sexes=[], all_groups=["adj"]
            )
        else:
            faf_label_groups = create_label_groups(
                pops=faf_pops, sexes=SEXES, all_groups=["adj"]
            )
        for label_group in faf_label_groups:
            vcf_info_dict.update(
                make_info_dict(
                    prefix=data_type,
                    prefix_before_metric=False,
                    pop_names=faf_pops,
                    label_groups=label_group,
                    faf=True,
                    description_text=description_text,
                )
            )
        logger.info(
            "Preparing grpmax and fafmax VCF header for %s subset...", data_type
        )
        vcf_info_dict.update(
            make_info_dict(
                suffix=data_type,
                pop_names=pops,
                grpmax=True,
                description_text=description_text,
            )
        )
        vcf_info_dict.update(
            make_info_dict(
                suffix=data_type,
                fafmax=True,
                description_text=description_text,
            )
        )

        logger.info("Preparing age histograms VCF header for %s subset...", data_type)
        vcf_info_dict.update(
            make_info_dict(
                suffix=data_type,
                bin_edges=make_hist_bin_edges_expr(ht, ann_with_hists=data_type),
                age_hist_distribution=hl.eval(
                    f"ht.{data_type}_globals.age_distribution)"
                ),
                description_text=description_text,
            ),
        )

        logger.info(
            "Preparing variant quality histograms VCF header for %s subset...",
            data_type,
        )
        vcf_info_dict.update(
            make_hist_dict(
                make_hist_bin_edges_expr(ht, ann_with_hists=data_type),
                adj=True,
                drop_n_smaller_larger=True,
                suffix=data_type,
                description_text=description_text,
            )
        )

    pops_ctt = pops.copy()
    pops_ctt.update({"ami": "Amish"})
    ctt_label_groups = create_label_groups(pops=pops_ctt, sexes=sexes)
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

    # Adjust keys to remove adj tags before exporting to VCF.
    new_vcf_info_dict = {i.replace("_adj", ""): j for i, j in vcf_info_dict.items()}

    vcf_header_dict = {
        "info": new_vcf_info_dict,
    }

    return vcf_header_dict


def prepare_ht_for_export(
    ht: hl.Table,
    vcf_info_reorder: Optional[List[str]] = VCF_INFO_REORDER,
    info_fields_to_drop: Optional[List[str]] = None,
) -> Tuple[hl.Table, List[str]]:
    """
    Format validated HT for export.

    Drop downsamplings frequency stats from info, rearrange info, and make sure fields
    are VCF compatible.

    :param ht: Input joint release HT.
    :param vcf_info_reorder: Order of VCF INFO fields. These will be placed in front of
        all other fields in the order specified.
    :param info_fields_to_drop: List of info fields to drop from the info struct.
    :return: Formatted HT and list rename row annotations.
    """
    if vcf_info_reorder:
        logger.info("Rearranging fields to desired order...")
        ht = ht.annotate(
            info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder))
        )

    # Select relevant fields for VCF export
    ht = ht.select("info")

    if info_fields_to_drop is None:
        info_fields_to_drop = []

    logger.info("Add age_histogram bin edges to info fields to drop...")
    for data_type in ["exomes", "genomes", "joint"]:
        info_fields_to_drop.extend(
            [
                f"age_hist_het_bin_edges_{data_type}",
                f"age_hist_hom_bin_edges_{data_type}",
            ]
        )

    logger.info(
        "Dropping the following fields from info struct: %s...",
        pprint(info_fields_to_drop),
    )
    ht = ht.annotate(info=ht.info.drop(*info_fields_to_drop))

    logger.info("Dropping '_adj' from info fields...")
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
