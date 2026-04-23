import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad.utils.filtering import filter_to_autosomes

from gnomad_qc.v5.resources.basics import get_logging_path

from typing import List
import argparse
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hom_depletion_check")
logger.setLevel(logging.INFO)


def get_field_name(group: str) -> str:
    """
    Get the field name for a group (for clarity, "adj" -> "overall").

    :param group: Group key (e.g., "adj", "afr_adj")
    :return: Field name (e.g., "overall", "afr_adj")
    """
    return "overall" if group == "adj" else group.replace("_adj", "")


def annotate_group_freqs_and_cutoffs(
    ht: hl.Table,
    groups: List[str],
    min_expected: int = 3,
):
    """
    Annotate a Hail Table with group-specific allele frequencies and AF cutoffs.

    For each group, this function:
      - Retrieves the corresponding frequency index from `ht.freq_index_dict`
      - Annotates `<field_name>_freq` with the group's frequency
      - Annotates `<field_name>_af_cutoff` using:
            sqrt(min_expected / (AN / 2))

    # If a group's index is not defined, the annotations are set to missing.

    :param ht: Hail Table containing `freq` (array of frequency structs) and`freq_index_dict` (dict mapping group names to indices in `freq`).
    :param groups: List of group keys to use (e.g., ["adj", "afr_adj"]).
    :param min_expected: Minimum expected allele count used to compute AF cutoff
    :return: Annotated Hail Table with added frequency and cutoff fields.
    """
    group_freqs = {}
    cutoffs = {}

    for group in groups:
        idx = ht.freq_index_dict.get(group)
        field_name = get_field_name(group)

        defined_idx = hl.is_defined(idx)

        freq_expr = hl.or_missing(
            defined_idx,
            ht.freq[idx],
        )

        group_freqs[f"{field_name}_freq"] = freq_expr

        cutoffs[f"{field_name}_af_cutoff"] = hl.or_missing(
            defined_idx,
            hl.sqrt(min_expected / (freq_expr.AN / 2)),
        )

    return ht.annotate(**group_freqs, **cutoffs)


def flag_hom_depletion(ht, groups):
    """
    Flag variants with depleted homozygous counts based on a specified metric and AF cutoff.

    A variant is flagged as hom-depleted if its frequency is greater than the AF cutoff and it has zero homozygous counts in the corresponding group.

    :param ht: Hail Table containing group-specific frequency and cutoff annotations.
    :param groups: List of group keys to use (e.g., ["adj", "afr_adj"]).
    :return: Hail Table with added boolean fields indicating hom-depletion for each group (e.g., "overall_flagged", "afr_adj_flagged").
    """
    # Dynamic threshold: AF must be greater than the square root of (desired_homs / N).
    flagged_groups = {}
    for group in groups:
        field_name = get_field_name(group)
        freq_field = f"{field_name}_freq"
        flagged_groups[f"{field_name}_hom_depletion"] = (
            hl.is_defined(ht[freq_field])
            & (ht[freq_field].AF > ht[f"{field_name}_af_cutoff"])
            & (ht[freq_field].homozygote_count == 0)
        )

    return ht.annotate(**flagged_groups)


def annotate_most_severe_gene(ht: hl.Table) -> hl.Table:
    most_severe_tx = hl.find(
        lambda tx: tx.consequence_terms.contains(ht.most_severe_csq),
        ht.vep.transcript_consequences,
    )

    return ht.annotate(most_severe_gene=most_severe_tx.gene_symbol)


def select_group_output_fields(ht, groups):
    """
    Select and flatten group-specific annotations.

    :param ht: Hail Table containing:
        - `<field_name>_freq` structs with AC, AN, AF, homozygote_count
        - `<field_name>_hom_depletion` boolean fields
        - `<field_name>_af_cutoff` fields
        - `most_severe_csq`, `most_severe_gene`, `lcr_or_segdup`
    :param groups: List of group keys (e.g., ["adj", "afr_adj"]).
    :return: Hail Table with selected and flattened fields.
    """
    output_fields = {}

    # per-group fields.
    for group in groups:
        field_name = get_field_name(group)
        freq = ht[f"{field_name}_freq"]

        output_fields.update(
            {
                f"{field_name}_AC": freq.AC,
                f"{field_name}_AN": freq.AN,
                f"{field_name}_AF": freq.AF,
                f"{field_name}_homozygote_count": freq.homozygote_count,
                f"{field_name}_hom_depletion": ht[f"{field_name}_hom_depletion"],
                f"{field_name}_af_cutoff": ht[f"{field_name}_af_cutoff"],
            }
        )

    # Per variant fields.
    output_fields.update(
        {
            "most_severe_csq": ht.most_severe_csq,
            "most_severe_gene": ht.most_severe_gene,
            "lcr_or_segdup": ht.lcr_or_segdup,
        }
    )

    return ht.select(**output_fields)


def main(args):
    """Perform checks for homozygous depletions."""
    hl.init(
        log="/hom_depletion_check.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")

    try:

        if args.adj_only:
            logger.info("Running hom depletion check with adj-only groups...")
            groups = ["adj"]
        else:
            logger.info("Running hom depletion check with all groups...")
            groups = [
                "adj",
                "afr_adj",
                "amr_adj",
                "asj_adj",
                "eas_adj",
                "fin_adj",
                "nfe_adj",
                "sas_adj",
            ]

        ht = public_release(args.data_type).ht()

        logger.info(
            "Filter to PASS variants in autosomes and annotate lcr/segdup/most_severe_consequence..."
        )
        ht = ht.filter(hl.len(ht.filters) == 0)
        ht = filter_to_autosomes(ht)
        ht = ht.annotate(
            lcr_or_segdup=(ht.region_flags.lcr | ht.region_flags.segdup),
            most_severe_csq=ht.vep.most_severe_consequence,
        )

        # Keep only fields needed for downstream depletion annotations.
        ht = ht.select(
            "freq",
            "most_severe_csq",
            "vep",
            "lcr_or_segdup",
        )

        logger.info("Annotate group-specific frequencies and AF cutoffs...")
        ht = annotate_group_freqs_and_cutoffs(ht, groups, args.min_expected)

        logger.info("Flag variants with depleted homozygous counts...")
        ht = flag_hom_depletion(ht, groups)

        logger.info("Pull out first gene with most severe consequence...")
        ht = annotate_most_severe_gene(ht)

        logger.info("Select and flatten output fields...")
        ht = select_group_output_fields(ht, groups)

        logger.info("Writing hom depletion check results...")
        ht = ht.naive_coalesce(50)
        ht.write("gs://gnomad-kristen/hom_depletion_af_check_all.ht", overwrite=True)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("hom_depletion_check"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--data_type",
        help="Dataset to run on, one of 'exomes', 'genomes', or 'joint'",
        default="exomes",
        choices=["exomes", "genomes", "joint"],
        type=str,
    )
    parser.add_argument(
        "--min_expected",
        help="Minimum expected allele count used to compute AF cutoff.",
        default=3,
        type=int,
    )
    parser.add_argument(
        "--adj_only",
        help="Whether to run hom depletion check using only the adj group (overall).",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
