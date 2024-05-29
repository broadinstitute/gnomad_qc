"""Find grpmax stats for gnomAD v4 and v2."""

import argparse
import logging
from pprint import pprint
from typing import Dict

import hail as hl
from gnomad.utils.vep import (
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_LOW_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    CSQ_NON_CODING,
)
from tabulate import tabulate

from gnomad_qc.v2.resources.basics import get_gnomad_public_data
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("grpmax_comps")
logger.setLevel(logging.INFO)


DIVERSE_GRPS = hl.literal({"afr", "amr", "eas", "mid", "sas"})
EUR_GRPS = {"all_eur": ["nfe", "fin", "asj"], "nfe_only": ["nfe"]}
FILTER_VALUES_TO_DROP = hl.array(["AC0", "InbreedingCoeff"])
AF_THRESHOLDS = [0.0001, 0.001]
NS_CONSEQ_TERMS = hl.literal(CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT)


def get_eur_freq(ht: hl.Table, eur_grps: list, version: str = "v4"):
    """
    Calculate European genetic ancestry super group AF.

    :param ht: gnomAD release Table
    :param eur_grps: Set of genetic ancestry groups to be included in AF
    :param version: gnomAD version
    :return: Table with European AF annotation
    """
    if version == "v4":
        eur_indices = hl.literal(eur_grps).map(lambda g: ht.freq_index_dict[g + "_adj"])
    else:
        eur_indices = hl.literal(eur_grps).map(
            lambda g: ht.freq_index_dict["gnomad_" + g]
        )
    ht = ht.annotate(
        eur_AC=hl.sum(
            eur_indices.map(
                lambda f: hl.if_else(hl.is_defined(ht.freq[f].AC), ht.freq[f].AC, 0)
            )
        ),
        eur_AN=hl.sum(
            eur_indices.map(
                lambda f: hl.if_else(hl.is_defined(ht.freq[f].AN), ht.freq[f].AN, 0)
            )
        ),
    )
    ht = ht.annotate(eur_AF=hl.if_else(ht.eur_AN > 0, ht.eur_AC / ht.eur_AN, 0.0))
    return ht.drop("eur_AC", "eur_AN")


def filter_to_threshold(
    ht: hl.Table,
    af_threshold: float = 0.001,
    version: str = "v4",
):
    """
    Filter to variants where eur AF is < threshold while grpmax > threshold.

    :param ht: gnomAD release Table
    :param af_threshold: AF dividing threshold
    :param version: gnomAD version
    :return: Table filtered to variants meething threshold specifications
    """
    if version == "v4":
        grpmax_expr = ht.grpmax.gnomad
    else:
        grpmax_expr = ht.popmax[0]
    ht = ht.filter((grpmax_expr.AF > af_threshold) & (ht.eur_AF < af_threshold))
    return ht


def version_stats(
    ht: hl.Table, version: str = "v4", can_only: bool = False, drop_hard_filtered_variants: bool = False,
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Calculate grpmax stats for a given gnomAD version.

    :param ht: gnomAD release Table
    :param version: gnomAD version
    :param can_only: Only consider MANE Select and canonical transcripts.
    :param drop_hard_filtered_variants: Remove only hard filtered variants but keep all other variants regardless of variant QC status
    :return: Dictionary of grpmax stats
    """
    logger.info("Calculating grpmax stats for %s", version)
    # Group results into nfe and all eur grps
    results_by_eur_grping = {}

    logger.info(
        "Total number of variants in %s before filtering: %s", version, ht.count()
    )
    t_variants = ht.count()
    if can_only:
        logger.info("Filtering to only MANE Select and canonical transcripts")
        ht = ht.explode(ht.vep.transcript_consequences)
        # All MANE select transcripts in v4 are also the canonical transcript
        ht = ht.filter(hl.is_defined(ht.vep.transcript_consequences.canonical))
        logger.info(
            "Filtering on mane_select and canonical transcript consequence term to "
            "keep only non-synonymous variants..."
        )
        ht = ht.filter(
            hl.any(
                ht.vep.transcript_consequences.consequence_terms.map(
                    lambda x: ~hl.literal(
                        CSQ_NON_CODING + CSQ_CODING_LOW_IMPACT
                    ).contains(x)
                )
            )
        )
        # Remove duplicate sites after exploding, most_severe_consequence is the same for all
        # transcript consequences so we can use distinct to remove duplicates
        ht = ht.distinct()
        ht = ht.checkpoint(
            f"gs://gnomad-tmp-4day/grpmax_comps_{version}_canonical_non_syn.ht",
            overwrite=True,
        )
        ns_variants = ht.count()
        logger.info(
            "Total number of variants in %s after canonical/mane and consequence "
            "filtering: %s (%.2f%% of total variants)",
            version,
            ns_variants,
            (ns_variants / t_variants) * 100,
        )

    else:
        # Filter to only non-synonymous terms CSQ_CODING_HIGH_IMPACT +
        # CSQ_CODING_MEDIUM_IMPACT
        logger.info(
            "Filtering on most_severe_consequence to keep only non-synonymous"
            " variants..."
        )
        ht = ht.filter(
            ~hl.literal(CSQ_NON_CODING + CSQ_CODING_LOW_IMPACT).contains(
                ht.vep.most_severe_consequence
            )
        )
        ht = ht.checkpoint(
            f"gs://gnomad-tmp-4day/grpmax_comps_{version}_non_syn.ht", overwrite=True
        )
        ns_variants = ht.count()
        logger.info(
            "Total number of variants in %s after consequence filtering: %s (%.2f%% of "
            "total variants)",
            version,
            ns_variants,
            (ns_variants / t_variants) * 100,
        )
    if drop_hard_filtered_variants:
        # Filter out AC0 and InbreedingCoeff variants
        ht = ht.filter(~ht.filters.any(lambda x: FILTER_VALUES_TO_DROP.contains(x)))
        log_message = "ONLY hard filters (RF/VQSR are retained)"
    else:
        ht = ht.filter(hl.len(ht.filters) == 0)
        log_message = "all filters"

    p_ns_variants = ht.count()
    logger.info(
        "Total number of variants passing %s: %s (%.2f%% of"
        " non-synonymous variants, %.2f%% of total variants)",
        log_message,
        p_ns_variants,
        (p_ns_variants / ns_variants) * 100,
        (p_ns_variants / t_variants) * 100,
    )
    # Iterate through just nfe group and all eur grp for calculating eur AF
    for grp_id, grps in EUR_GRPS.items():
        # Get european frequency by calculating the cumulative AF in the passed
        # europe genetic ancestry list
        p_ht = get_eur_freq(ht, eur_grps=grps, version=version)

        grpmax_by_thresholds = {}
        for threshold in AF_THRESHOLDS:
            # Filter to keep only variants that split at the AF threshold
            t_ht = filter_to_threshold(p_ht, threshold, version=version)
            t_ht = t_ht.checkpoint(
                f"gs://gnomad-tmp-4day/grpmax_comps_{version}_{grp_id}_{threshold}.ht",
                overwrite=True,
            )
            logger.info(
                "Number of variants after filtering to AF threshold for %s at %s "
                "threshold in %s: %s (%.2f%% of passing non-synonymous variants, %.2f%%"
                " of total variants)",
                grp_id,
                threshold,
                version,
                t_ht.count(),
                t_ht.count() / p_ns_variants * 100,
                t_ht.count() / t_variants * 100,
            )
            if version == "v4":
                grpmax_ga_expr = t_ht.grpmax.gnomad.gen_anc
            else:
                grpmax_ga_expr = t_ht.popmax[0].pop

            t_ht = t_ht.filter(DIVERSE_GRPS.contains(grpmax_ga_expr))
            if version == "v4":
                grpmax_ga_expr = t_ht.grpmax.gnomad.gen_anc
            else:
                grpmax_ga_expr = t_ht.popmax[0].pop
            # For each diverse genetic ancestry group, aggregate the number of
            # variants where that group is grpmax
            grpmax_by_thresholds[threshold] = t_ht.aggregate(
                hl.agg.counter(grpmax_ga_expr)
            )
            # Add a total count for all groups
            grpmax_by_thresholds[threshold]["all"] = t_ht.count()

        results_by_eur_grping[grp_id] = grpmax_by_thresholds
    logger.info(f"Results for {version}")
    pprint(results_by_eur_grping, sort_dicts=False)
    return results_by_eur_grping


def create_table(
    version_dict: Dict[str, Dict[str, Dict[str, int]]],
    data_subset: str,
):
    """
    Create tables of grpmax stats.

    Tables have a column for each version as well as the difference for each grpmax group.

    :param version_dict: Dictionary of grpmax stats
    :param data_subset: Data subset
    """
    headers = ["grpmax_grp", "v2", "v4", "difference(v4-v2)"]
    for threshold in AF_THRESHOLDS:
        table = []
        for grpmax_group in hl.eval(hl.array(DIVERSE_GRPS)) + ["all"]:
            v2_val = version_dict["v2"][data_subset][threshold].get(grpmax_group, 0)
            v4_val = version_dict["v4"][data_subset][threshold].get(grpmax_group, 0)
            diff = v4_val - v2_val
            if diff < 0:
                diff = f"\033[91m{diff}\033[0m"
            else:
                diff = f"\033[92m{diff}\033[0m"
            table.append([grpmax_group, v2_val, v4_val, diff])
        logger.info(
            "\nNonsynonymous variant count by grpmax genetic ancestry group where the\n"
            f"grpmax AF is above {threshold} and the {data_subset} AF is below it...\n"
            f"{tabulate(table,headers=headers,tablefmt='fancy_grid')}"
        )


def main(args):
    """Find grpmax stats for gnomAD v4 and v2."""
    version_dict = {}
    for version in ["v4", "v2"]:
        if version == "v4":
            ht = release_sites().ht()
        else:
            ht = get_gnomad_public_data("exomes")
        version_dict[version] = version_stats(
            ht, version=version, can_only=args.canonical_only
        )

    # Create tables for "all_eur" and "nfe_only" data
    create_table(version_dict, data_subset="all_eur")
    create_table(version_dict, data_subset="nfe_only")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--canonical-only",
        action="store_true",
        help="Only consider MANE Select and canonical transcripts",
    )
    parser.add_argument(
        "--drop-hard-filtered-variants",
        action="store_true",
        help="Drop variants that have been hard filtered",
    )
    args = parser.parse_args()
    main(args)
