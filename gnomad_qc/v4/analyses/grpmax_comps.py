"""Find grpmax stats for gnomAD v4 and v2."""

import logging
from pprint import pprint
from typing import Dict

import hail as hl
from gnomad.utils.vep import CSQ_CODING_HIGH_IMPACT, CSQ_CODING_MEDIUM_IMPACT
from tabulate import tabulate

from gnomad_qc.v2.resources.basics import get_gnomad_public_data
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("grpmax_comps")
logger.setLevel(logging.INFO)


DIVERSE_GRPS = hl.literal({"afr", "amr", "eas", "mid", "sas"})
EUR_GRPS = {"all_eur": ["nfe", "fin", "asj"], "nfe_only": ["nfe"]}

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
        eur_AF=hl.sum(
            eur_indices.map(
                lambda f: hl.if_else(hl.is_defined(ht.freq[f].AC), ht.freq[f].AC, 0)
            )
        )
        / hl.sum(
            eur_indices.map(
                lambda f: hl.if_else(hl.is_defined(ht.freq[f].AN), ht.freq[f].AN, 0)
            )
        )
    )
    ht = ht.filter(hl.is_defined(ht.eur_AF))

    return ht


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
    ht: hl.Table, version: str = "v4"
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Calculate grpmax stats for a given gnomAD version.

    :param ht: gnomAD release Table
    :param version: gnomAD version
    :return: Dictionary of grpmax stats
    """
    logger.info(f"Calculating grpmax stats for {version}")
    # Group results into nfe and all eur grps
    results_by_eur_grping = {}
    # Filter to only non-synonymous terms CSQ_CODING_HIGH_IMPACT +
    # CSQ_CODING_MEDIUM_IMPACT
    ht = ht.filter(NS_CONSEQ_TERMS.contains(ht.vep.most_severe_consequence))
    # Filter to PASS variants only
    ht = ht.filter(ht.filters.length() == 0)
    # Iterate through just nfe group and all eur grp for calculating eur AF
    for grp_id, grps in EUR_GRPS.items():
        # Get european frequency by calculating the cumulative AF in the passed
        # europe genetic ancestry list
        p_ht = get_eur_freq(ht, eur_grps=grps, version=version)

        grpmax_by_thresholds = {}
        for threshold in AF_THRESHOLDS:
            # Filter to keep only variants that split at the AF threshold
            t_ht = filter_to_threshold(p_ht, threshold, version=version)
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
        for grpmax_group in hl.eval(DIVERSE_GRPS):
            v2_val = version_dict["v2"][data_subset][threshold].get(grpmax_group, 0)
            v4_val = version_dict["v4"][data_subset][threshold].get(grpmax_group, 0)
            diff = v4_val - v2_val
            if diff < 0:
                diff = f"\033[91m{diff}\033[0m"
            else:
                diff = f"\033[92m{diff}\033[0m"
            table.append([grpmax_group, v2_val, v4_val, diff])
        logger.info(
            "Nonsynonymous variant count by grpmax genetic ancestry group where the "
            f"grpmax AF is above {threshold} and the {data_subset} AF is below it..."
        )
        logger.info(
            tabulate(
                table,
                headers=headers,
                tablefmt="fancy_grid",
            )
        )


def main():
    """Find grpmax stats for gnomAD v4 and v2."""
    version_dict = {}
    for version in ["v4", "v2"]:
        if version == "v4":
            ht = release_sites().ht()
        else:
            ht = get_gnomad_public_data("exomes")
        version_dict[version] = version_stats(ht, version=version)

    # Create tables for "all_eur" and "nfe_only" data
    create_table(version_dict, data_subset="all_eur")
    create_table(version_dict, data_subset="nfe_only")


if __name__ == "__main__":
    main()
