"""Find grpmax stats for gnomAD v4 and v2."""

import argparse
import logging
from copy import deepcopy
from pprint import pprint
from typing import Dict

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.vep import (
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    filter_vep_to_canonical_transcripts,
    process_consequences,
)
from tabulate import tabulate

from gnomad_qc.v2.resources.basics import (
    get_gnomad_liftover_data_path,
    get_gnomad_public_data,
)
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("grpmax_comps")
logger.setLevel(logging.INFO)


DIVERSE_GRPS = ["afr", "amr", "eas", "mid", "sas"]
EUR_GRPS = {"all_eur": ["nfe", "fin", "asj"], "nfe_only": ["nfe"]}
FILTER_VALUES_TO_DROP = hl.array(["AC0", "InbreedingCoeff"])
AF_THRESHOLDS = [0.0001, 0.001, 0.01]
NS_CONSEQ_TERMS = CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT


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
    ht: hl.Table, af_threshold: float = 0.001, version: str = "v4", eur_filter=False
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
    filter_expr = grpmax_expr.AF > af_threshold
    if eur_filter:
        filter_expr &= ht.eur_AF < af_threshold
    ht = ht.filter(filter_expr)
    return ht


def version_stats(
    ht: hl.Table,
    t_variants: int,
    version: str = "v4",
    grpmax_counts: bool = False,
    eur_filter: bool = False,
) -> Dict[str, Dict[str, Dict[str, int]]]:
    """
    Calculate grpmax stats for a given gnomAD version.

    :param ht: gnomAD release Table
    :param t_variants: Total number of variants in the Table
    :param version: gnomAD version
    :param grpmax_counts: Calculate using only grpmax counts for diverse genetic ancestry groups
    :return: Dictionary of grpmax stats
    """
    results_dict = {}

    if grpmax_counts:
        if version == "v4":
            ht = ht.annotate(
                grpmax_ga=ht.grpmax.gnomad.gen_anc,
                grpmax_af=ht.grpmax.gnomad.AF,
            )
        else:
            ht = ht.annotate(
                grpmax_ga=ht.popmax[0].pop,
                grpmax_af=ht.popmax[0].AF,
            )

        if eur_filter:
            logger.info(
                "Calculating grpmax counts for diverse genetic ancestry groups to "
                "deduplicate variant counts..."
            )
            # Iterate through just nfe group and all eur grp for calculating eur AF
            for grp_id, grps in EUR_GRPS.items():
                # Get european frequency by calculating the cumulative AF in the passed
                # europe genetic ancestry list
                p_ht = get_eur_freq(ht, eur_grps=grps, version=version)

                counts_by_thresholds = {}
                for threshold in AF_THRESHOLDS:
                    t_ht = filter_to_threshold(
                        p_ht, threshold, version=version, eur_filter=eur_filter
                    )
                    t_ht = t_ht.checkpoint(
                        f"gs://gnomad-tmp-4day/grpmax_comps_{version}_{grp_id}_{threshold}.ht",
                        overwrite=True,
                    )
                    logger.info(
                        "Number of variants after filtering grpmax AF threshold for %s"
                        " at %s threshold in %s: %s ( %.2f%% of total variants)",
                        grp_id,
                        threshold,
                        version,
                        t_ht.count(),
                        t_ht.count() / t_variants * 100,
                    )

                    t_ht = t_ht.filter(
                        hl.literal(DIVERSE_GRPS).contains(t_ht.grpmax_ga)
                    )

                    # For each diverse genetic ancestry group, aggregate the number of
                    # variants where that group is grpmax
                    counts_by_thresholds[threshold] = t_ht.aggregate(
                        hl.agg.counter(t_ht.grpmax_ga)
                    )
                    # Add a total count for all groups
                    counts_by_thresholds[threshold]["all"] = t_ht.count()
                results_dict[grp_id] = counts_by_thresholds
        else:
            logger.info(
                "Aggregating grpmax variant counts for genetic ancestry groups..."
            )
            grps = deepcopy(DIVERSE_GRPS) + ["nfe"]
            for threshold in AF_THRESHOLDS:
                t_ht = filter_to_threshold(ht, threshold, version=version)
                results_dict[threshold] = t_ht.aggregate(
                    hl.dict(
                        {grp: hl.agg.count_where(t_ht.grpmax_ga == grp) for grp in grps}
                    )
                )
                results_dict[threshold]["all"] = t_ht.count()

    else:
        logger.info(
            "Aggregating variant counts for genetic ancestry groups, variants "
            "could be counted multiple times if they are above the AF threshold for "
            "multiple groups..."
        )
        grps = deepcopy(DIVERSE_GRPS) + ["nfe"]
        grp_index_dict = {
            g: ht.freq_index_dict.get(g + "_adj" if version == "v4" else "gnomad_" + g)
            for g in grps
        }
        for threshold in AF_THRESHOLDS:
            results_dict[threshold] = ht.aggregate(
                hl.dict(
                    {
                        grp: hl.agg.count_where(ht.freq[idx].AF > threshold)
                        for grp, idx in grp_index_dict.items()
                    }
                )
            )

    logger.info(f"Results for {version}")
    pprint(results_dict, sort_dicts=False)
    return results_dict


def create_table(
    version_dict: Dict[str, Dict[str, Dict[str, int]]],
    data_subset: str = None,
    grpmax_counts: bool = False,
    non_syn_only: bool = False,
    eur_filter: bool = False,
):
    """
    Create tables of grpmax stats.

    Tables have a column for each version as well as the difference for each grpmax group.

    :param version_dict: Dictionary of grpmax stats
    :param data_subset: Data subset
    :param grpmax_counts: Results calculated using only grpmax counts for diverse genetic ancestry groups
    :param non_syn_only: Keep only non-synonymous variants
    :param eur_filter: Filter to variants where eur AF is < threshold while grpmax > threshold
    """
    headers = ["genetic ancestry group", "v2", "v4", "difference(v4-v2)"]
    v2_dict_index = (
        version_dict["v2"][data_subset]
        if grpmax_counts and eur_filter
        else version_dict["v2"]
    )
    v4_dict_index = (
        version_dict["v4"][data_subset]
        if grpmax_counts and eur_filter
        else version_dict["v4"]
    )
    grps = deepcopy(DIVERSE_GRPS)
    grps = grps + ["all"]

    if not eur_filter:
        grps.append("nfe")

    for threshold in AF_THRESHOLDS:
        table = []
        for grp in grps:
            v2_val = v2_dict_index[threshold].get(grp, 0)
            v4_val = v4_dict_index[threshold].get(grp, 0)
            diff = v4_val - v2_val
            if diff < 0:
                diff = f"\033[91m{diff}\033[0m"
            else:
                diff = f"\033[92m{diff}\033[0m"
            table.append([grp, v2_val, v4_val, diff])
        logger.info(
            f"\n{'non-synonymous' if non_syn_only else ''} variant count by genetic "
            f"ancestry group where the{' grpmax' if grpmax_counts else ''} AF is above"
            f" {threshold}{' and the ' + str(data_subset) + ' AF is below it' if grpmax_counts and eur_filter else ''}...\n"
            f"{tabulate(table, headers=headers, tablefmt='fancy_grid')}"
        )


def main(args):
    """Find grpmax stats for gnomAD v4 and v2."""
    version_dict = {}
    grpmax_counts = args.grpmax_counts
    non_syn_only = args.non_syn_only
    csq_terms = args.csq_terms
    eur_filter = args.eur_filter
    csq_terms = csq_terms if csq_terms else non_syn_only

    for version in ["v4", "v2"]:
        if version == "v4":
            ht = release_sites().ht()
        elif args.use_v2_liftover:
            ht = hl.read_table(get_gnomad_liftover_data_path("exomes"))
            vep_ht = vep_context.versions["105"].ht()
            ht = ht.annotate(vep=vep_ht[ht.key].vep)
        else:
            ht = get_gnomad_public_data("exomes")

        if args.test:
            # Filter to two partitions
            ht = ht._filter_partitions([0, 1])
        t_variants = ht.count()

        logger.info(
            "Total number of variants in %s before filtering: %s", version, t_variants
        )
        msg = ""
        # All MANE Select transcripts are canonical
        if args.canonical:
            logger.info("Filtering to only MANE Select and canonical transcripts...")
            ht = process_consequences(ht)
            ht = filter_vep_to_canonical_transcripts(ht, filter_empty_csq=True)
            msg = "canonical transcript filtering"

        if csq_terms:
            logger.info(
                "Filtering to keep only %s ...",
                (
                    str(csq_terms) + " variants"
                    if csq_terms and not non_syn_only
                    else "non-synonymous variants"
                ),
            )
            # NOTE: There is no guarantee the most severe consequence is from the
            # canonical transcript if table is filtered to canonical transcripts
            ht = ht.filter(
                hl.literal(csq_terms).contains(
                    ht.vep.worst_csq_for_variant_canonical
                )  # Need to determine field, process_consequences produces this one, should that function move up?
                # Should we just run on canonical alwways?
                # IF so move and keep this field, otherwise add more conditionals to the
                # vep field that we access
            )
            ns_variants = ht.count()

            if msg:
                msg += " and "
            msg += (
                f"{'non-synonymous consequence' if non_syn_only else str(csq_terms) + ' filtering'}"
            )

            logger.info(
                "Total number of variants in %s %s: %s (%.2f%% of total variants)",
                version,
                "after " + msg if msg else "",
                ns_variants,
                (ns_variants / t_variants) * 100,
            )

        ht = ht.checkpoint(
            f"gs://gnomad-tmp-4day/grp_comps_{version}_canonical_non_syn.ht",
            overwrite=True,
        )

        if args.no_variant_qc_filters:
            msg = "of variants, regardless of variant QC status"
        elif args.drop_hard_filtered:
            # Filter out AC0 and InbreedingCoeff variants
            ht = ht.filter(~ht.filters.any(lambda x: FILTER_VALUES_TO_DROP.contains(x)))
            msg = "of variants passing ONLY hard filters (RF/VQSR are retained)"
        else:
            # Filter to only PASS
            ht = ht.filter(hl.len(ht.filters) == 0)
            msg = "of variants passing all filters"

        p_ns_variants = ht.count()
        logger.info(
            "Total number %s: %s (%.2f%% of total variants)",
            msg,
            p_ns_variants,
            (p_ns_variants / t_variants) * 100,
        )

        version_dict[version] = version_stats(
            ht,
            t_variants=t_variants,
            version=version,
            grpmax_counts=grpmax_counts,
            eur_filter=eur_filter,
        )

    if grpmax_counts and eur_filter:
        for eur_grp in EUR_GRPS.keys():
            create_table(
                version_dict,
                data_subset=eur_grp,
                grpmax_counts=grpmax_counts,
                non_syn_only=non_syn_only,
                eur_filter=eur_filter,
            )
    else:
        create_table(
            version_dict, grpmax_counts=grpmax_counts, non_syn_only=non_syn_only
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--canonical",
        action="store_true",
        help="Only consider MANE Select and canonical transcripts",
    )
    parser.add_argument(
        "--drop-hard-filtered",
        action="store_true",
        help="Drop variants that have been hard filtered",
    )
    parser.add_argument(
        "--no-variant-qc-filters",
        action="store_true",
        help="Keep all variants regardless of variant QC status",
    )
    parser.add_argument(
        "--non-syn-only",
        action="store_true",
        help="Keep only non-synonymous variants",
    )
    parser.add_argument(
        "--grpmax-counts",
        action="store_true",
        help="Calculate using only grpmax counts for diverse genetic ancestry groups",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Filter to two partitions of each HT for testing purposes",
    )
    parser.add_argument(
        "--csq-terms",
        nargs="+",
        help="Consequence terms to consider",
        default=NS_CONSEQ_TERMS,
    )
    parser.add_argument(
        "--eur-filter",
        action="store_true",
        help="Filter to variants where eur AF is < threshold while grpmax > threshold",
    )
    parser.add_argument(
        "--use-v2-liftover",
        action="store_true",
        help="Use v2 liftover annotated with vep 105 for comparison",
    )

    args = parser.parse_args()
    if args.drop_hard_filtered and args.no_variant_qc_filters:
        raise ValueError(
            "Cannot use --drop-hard-filtered-variants and --all-variants together"
        )
    if args.non_syn_only and (args.csq_terms != NS_CONSEQ_TERMS):
        raise ValueError("Cannot use --non-syn-only and --csq-terms together")
    if args.eur_filter and not args.grpmax_counts:
        raise ValueError("Cannot use --eur-filter without --grpmax-counts")
    main(args)
