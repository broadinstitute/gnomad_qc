"""Script to create release de novo HT for v4 exomes."""

import argparse
import logging
from typing import Tuple

import hail as hl
from gnomad.resources.grch38.gnomad import public_release
from gnomad.sample_qc.relatedness import default_get_de_novo_expr
from gnomad.utils.annotations import annotate_adj, get_copy_state_by_sex
from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import (
    CSQ_CODING_HIGH_IMPACT,
    CSQ_CODING_LOW_IMPACT,
    CSQ_CODING_MEDIUM_IMPACT,
    process_consequences,
)
from hail.utils import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_vep
from gnomad_qc.v4.resources.release import release_de_novo
from gnomad_qc.v4.resources.sample_qc import dense_trio_mt, pedigree, trio_denovo_ht
from gnomad_qc.v4.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_de_novo_release_ht")
logger.setLevel(logging.INFO)


def get_releasable_de_novo_calls_ht(
    mt: hl.MatrixTable,
    priors_ht: hl.Table,
    ped: hl.Pedigree,
    test: bool = False,
    n_partitions: int = 1000,
) -> hl.Table:
    """
    Get de novo calls Hail Table.

    :param mt: Dense MatrixTable of the releasable trios.
    :param priors_ht: Table with AFs used as population frequency priors in de novo calculations.
    :param ped: Pedigree.
    :param test: Whether to filter to chr20 for testing. Default is False.
    :param n_partitions: Number of partitions to use for the output Table. Default is 1000.
    :return: Hail Table with de novo calls.
    """
    if test:
        logger.info("Filtering to chr20 for testing...")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20")])

    # The split MT was made with these extra annotations.
    mt = mt.annotate_rows(
        n_unsplit_alleles=hl.len(mt.alleles),
        mixed_site=(
            (hl.len(mt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:])
        ),
    )
    mt = hl.experimental.sparse_split_multi(mt).checkpoint(
        dense_trio_mt(releasable=True, split=True, test=test).path, _read_if_exists=True
    )
    mt = annotate_adj(mt)

    # v3 and v4 VDS have the PL and AD fields for homref genotypes intentionally
    # removed to save storage space and costs. We need to approximate the
    # AD and PL fields when missing.
    mt = mt.annotate_entries(
        AD=hl.or_else(mt.AD, [mt.DP, 0]),
        PL=hl.or_else(mt.PL, [0, mt.GQ, 2 * mt.GQ]),
    )
    mt = mt.annotate_rows(prior=priors_ht[mt.row_key].freq[0].AF)
    mt = mt.select_entries("GT", "AD", "DP", "GQ", "PL", "adj")

    tm = hl.trio_matrix(mt, ped, complete_trios=True)
    tm = tm.transmute_cols(is_xx=tm.is_female)
    tm = tm.checkpoint(new_temp_file("trio_matrix", "mt"))

    ht = tm.entries()

    ht = ht.annotate(
        de_novo_call_info=default_get_de_novo_expr(
            locus_expr=ht.locus,
            alleles_expr=ht.alleles,
            proband_expr=ht.proband_entry,
            father_expr=ht.father_entry,
            mother_expr=ht.mother_entry,
            is_xx_expr=ht.is_xx,
            freq_prior_expr=ht.prior,
        )
    )
    final_partitions = n_partitions // 10 if test else n_partitions
    ht = ht.filter(ht.de_novo_call_info.is_de_novo).naive_coalesce(final_partitions)
    return ht


def aggregate_and_annotate_de_novos(ht: hl.Table) -> hl.Table:
    """
    Aggregate and annotate de novo calls.

    This step produces high-quality coding de novos for the release, including variants
    that are not in low-confidence regions, do not have a `*` alt allele, and are not
    excluded by variant QC A subset of filtered de novos, further restricted to HIGH
    confidence or MEDIUM confidence with a P-value >= 0.9, coding consequence,
    and gnomAD v4.1 exomes allele frequency (AF) and callset allele count (AC) filters.

    :param ht: De novo calls Table.
    :return: Aggregated and annotated Table.
    """
    # Get needed annotation Tables.
    filters_ht = final_filter(all_variants=True, only_filters=True).ht()
    freq_ht = public_release("exomes").ht()
    vep_ht = get_vep(data_type="exomes").ht()

    # Filter to variants with a confidence level
    ht = ht.filter(~hl.is_missing(ht.de_novo_call_info.confidence))
    logger.info(f"Getting {ht.count()} de novos with a confidence level...")

    ht = filter_low_conf_regions(
        ht, filter_lcr=True, filter_decoy=False, filter_segdup=True
    )
    logger.info(
        f"Getting {ht.count()} de novos after removing low confidence regions..."
    )

    diploid_expr, hemi_x_expr, hemi_y_expr = get_copy_state_by_sex(
        ht.locus,
        ht.is_xx,
    )

    adj_trio = (
        hl.case()
        .when(
            diploid_expr,
            ht.proband_entry.adj & ht.father_entry.adj & ht.mother_entry.adj,
        )
        .when(hemi_x_expr, ht.proband_entry.adj & ht.mother_entry.adj)
        .when(hemi_y_expr, ht.proband_entry.adj & ht.father_entry.adj)
        .or_missing()
    )

    ht = ht.annotate(adj_trio=adj_trio)

    ht = (
        ht.group_by("locus", "alleles")
        .aggregate(
            de_novo_AC=hl.struct(
                AC_all=hl.agg.count_where(ht.adj_trio),
                AC_all_raw=hl.agg.count(),
                AC_high_conf=hl.agg.count_where(
                    (ht.de_novo_call_info.confidence == "HIGH") & ht.adj_trio
                ),
                AC_medium_conf=hl.agg.count_where(
                    (ht.de_novo_call_info.confidence == "MEDIUM") & ht.adj_trio
                ),
                # We're adding some MEDIUM confidence calls to high-quality because
                # we found the rate of some coding consequences didn't meet the
                # expectations, based on the distribution of p_de_novo values at
                # MEDIUM level, we chose a reasonable threshold of 0.9.
                AC_medium_conf_P_0_9=hl.agg.count_where(
                    (ht.de_novo_call_info.confidence == "MEDIUM")
                    & (ht.de_novo_call_info.p_de_novo >= 0.9)
                    & ht.adj_trio
                ),
                AC_low_conf=hl.agg.count_where(
                    (ht.de_novo_call_info.confidence == "LOW") & ht.adj_trio
                ),
            ),
            p_de_novo_stats=hl.agg.stats(ht.de_novo_call_info.p_de_novo),
            # The mixed_site info should stay the same for each variant.
            mixed_site=hl.agg.take(ht.mixed_site, 1)[0],
        )
        .key_by("locus", "alleles")
    )
    ht = ht.annotate(p_de_novo_stats=ht.p_de_novo_stats.drop("n", "sum")).checkpoint(
        new_temp_file("denovo_ac", "ht")
    )
    logger.info(f"Getting {ht.count()} unique de novo variants after computing AC...")

    ht = ht.annotate(vep=vep_ht[ht.key].vep)
    ht = process_consequences(ht, has_polyphen=False)

    ht = ht.annotate(
        pass_filters=filters_ht[ht.key]
        # AC0 is not used as a criteria in de novos because AC0 was used for
        # release sites and some of the probands weren't included in the release.
        .filters.intersection(hl.set(["AS_VQSR", "InbreedingCoeff"])).length() == 0,
        worst_csq_for_variant=ht.vep.worst_csq_for_variant.most_severe_consequence,
        worst_csq_gene_id=ht.vep.worst_csq_for_variant.gene_id,
        worst_csq_transcript_id=ht.vep.worst_csq_for_variant.transcript_id,
        gnomAD_v4_exomes_AF=freq_ht[ht.key].freq[0].AF,
    ).drop("vep")

    ht = ht.filter(ht.alleles[1] != "*").checkpoint(
        new_temp_file("denovo_no_star", "ht")
    )
    logger.info(f"Getting {ht.count()} de novos after filtering out '*' alts...")

    ht = ht.filter(ht.pass_filters).checkpoint(
        new_temp_file("denovo_pass_filters", "ht")
    )
    logger.info(
        f"Getting {ht.count()} de novos after AS_VQSR and InbreedingCoeff "
        f"filters..."
    )

    ht = ht.filter(
        (ht.de_novo_AC.AC_high_conf > 0) | (ht.de_novo_AC.AC_medium_conf_P_0_9 > 0)
    ).checkpoint(new_temp_file("denovo_high_conf", "ht"))
    logger.info(f"Getting {ht.count()} high-quality de novos...")

    # Get the set of coding variant except for the splice region variants for
    # calibration.
    coding_csqs = hl.literal(
        set(CSQ_CODING_HIGH_IMPACT + CSQ_CODING_MEDIUM_IMPACT + CSQ_CODING_LOW_IMPACT)
        - {
            "splice_polypyrimidine_tract_variant",
            "splice_donor_5th_base_variant",
            "splice_region_variant",
            "splice_donor_region_variant",
            "incomplete_terminal_codon_variant",
            "coding_transcript_variant",
        }
    )
    ht = ht.filter(coding_csqs.contains(ht.worst_csq_for_variant)).checkpoint(
        new_temp_file("denovo_high_coding", "ht")
    )
    logger.info(f"Getting {ht.count()} high-quality coding de novos...")

    ht = ht.filter(
        (ht.gnomAD_v4_exomes_AF < 0.001) | hl.is_missing(ht.gnomAD_v4_exomes_AF)
    ).checkpoint(new_temp_file("denovo_coding_af0001", "ht"))
    logger.info(
        f"Getting {ht.count()} high-quality coding de novos with AF < 0.001 in "
        f"gnomAD..."
    )

    ht = ht.filter(
        ht.de_novo_AC.AC_high_conf + ht.de_novo_AC.AC_medium_conf_P_0_9 <= 15
    )
    logger.info(
        f"Getting {ht.count()} high-quality coding de novos with AC <= 15 in the "
        f"callset..."
    )

    ht = ht.drop(ht.de_novo_AC.AC_medium_conf, ht.de_novo_AC.AC_low_conf)
    return ht


def restructure_for_tsv(ht: hl.Table) -> hl.Table:
    """
    Restructure the de novo release HT for TSV export.

    .. note::
       This function will only be used for the TSV of high-quality coding de novos.

    :param ht: De novo release HT.
    :return: Restructured HT.
    """
    return ht.flatten().rename(
        {
            "de_novo_AC.AC_all": "de_novo_AC_all",
            "de_novo_AC.AC_all_raw": "de_novo_AC_all_raw",
            "de_novo_AC.AC_high_conf": "de_novo_AC_high_conf",
            "de_novo_AC.AC_medium_conf_P_0_9": "de_novo_AC_medium_conf_P_0_9",
            "p_de_novo_stats.mean": "p_de_novo_mean",
            "p_de_novo_stats.stdev": "p_de_novo_stdev",
            "p_de_novo_stats.min": "p_de_novo_min",
            "p_de_novo_stats.max": "p_de_novo_max",
        }
    )


def main(args):
    """Create de novo release ht."""
    hl.init(
        log=f"/create_de_novo_release_ht.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")

    test = args.test
    overwrite = args.overwrite

    if args.generate_de_novo_calls:
        mt = dense_trio_mt(releasable=True).mt()
        priors_ht = public_release("exomes").ht()
        ped = hl.Pedigree.read(pedigree().path, delimiter="\t")

        ht = get_releasable_de_novo_calls_ht(
            mt, priors_ht, ped, test, n_partitions=args.n_partitions
        )
        ht.write(
            trio_denovo_ht(
                test=test,
            ).path,
            overwrite=overwrite,
        )

    if args.generate_final_hts:
        ht = hl.read_table(trio_denovo_ht(test=test).path)
        ht = aggregate_and_annotate_de_novos(ht)
        ht.naive_coalesce(20).write(
            release_de_novo(test=test).path,
            overwrite=overwrite,
        )
        logger.info("Final de novo release HT schema:")
        ht.describe()

    if args.generate_final_tsv:
        ht = release_de_novo(test=test).ht()
        ht = restructure_for_tsv(ht)
        ht.export(
            release_de_novo(test=test).path.replace(".ht", ".tsv.bgz"),
            header=True,
            delimiter="\t",
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "-t",
        "--test",
        help="Runs a test on chr20",
        action="store_true",
    )
    parser.add_argument(
        "--n_partitions",
        help="Number of partitions to use for the output Table.",
        type=int,
        default=1000,
    )
    parser.add_argument(
        "--generate-de-novo-calls",
        help="Generate de novo calls HT.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-final-hts",
        help="Generate the final HT for de novo release.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-final-tsv",
        help="Generate the final TSV for de novo release.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
