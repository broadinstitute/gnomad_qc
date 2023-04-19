"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl
from gnomad.sample_qc.relatedness import generate_trio_stats_expr
from gnomad.utils.annotations import (
    add_variant_type,
    annotate_adj,
    get_adj_expr,
    get_lowqual_expr,
)
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import (
    INFO_INT32_SUM_AGG_FIELDS,
    INFO_SUM_AGG_FIELDS,
    default_compute_info,
    get_as_info_expr,
    get_site_info_expr,
    split_info_annotation,
    split_lowqual_annotation,
)
from gnomad.utils.vcf import adjust_vcf_incompatible_types
from gnomad.utils.vep import vep_or_lookup_vep

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import (
    allele_data,
    fam_stats,
    get_info,
    get_transmitted_singleton_vcf_path,
    info_vcf_path,
    qc_ac,
    vep,
)
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import trios

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def split_info() -> hl.Table:
    """
    Generate an info table that splits multi-allelic sites from the multi-allelic info table.

    :return: Info table with split multi-allelics
    :rtype: Table
    """
    info_ht = get_info(split=False).ht()

    # Create split version
    info_ht = hl.split_multi(info_ht)

    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation(info_ht.info, info_ht.a_index),
        ),
        AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )
    return info_ht


def generate_allele_data(ht: hl.Table) -> hl.Table:
    """
    Return bi-allelic sites HT with an 'allele_data' annotation.

    'allele_data' is a struct with the following information:
        - nonsplit_alleles
        - has_star
        - variant_type
        - n_alt_alleles

    :param Table ht: Full unsplit HT
    :return: Table with allele data annotations
    :rtype: Table
    """
    ht = ht.select()
    allele_data = hl.struct(
        nonsplit_alleles=ht.alleles, has_star=hl.any(lambda a: a == "*", ht.alleles)
    )
    ht = ht.annotate(allele_data=allele_data.annotate(**add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    ht = ht.filter(hl.len(ht.alleles) > 1)
    allele_type = (
        hl.case()
        .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv")
        .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), "ins")
        .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), "del")
        .default("complex")
    )
    ht = ht.annotate(
        allele_data=ht.allele_data.annotate(
            allele_type=allele_type, was_mixed=ht.allele_data.variant_type == "mixed"
        )
    )
    return ht


def generate_ac(mt: hl.MatrixTable) -> hl.Table:
    """
    Create Table containing allele counts per variant.

    Returns table containing the following annotations:
        - `ac_qc_samples_raw`: Allele count of high quality samples
        - `ac_qc_samples_unrelated_raw`: Allele count of high quality unrelated samples
        - `ac_release_samples_raw`: Allele count of release samples
        - `ac_qc_samples_adj`: Allele count of high quality samples after adj filtering
        - `ac_qc_samples_unrelated_adj`: Allele count of high quality unrelated samples after adj filtering
        - `ac_release_samples_adj`: Allele count of release samples after adj filtering

    :param mt: Input MatrixTable
    :return: Table containing allele counts
    """
    mt = mt.filter_cols(mt.meta.high_quality)
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = annotate_adj(mt)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        ac_qc_samples_unrelated_raw=hl.agg.filter(
            ~mt.meta.sample_filters.all_samples_related,
            hl.agg.sum(mt.GT.n_alt_alleles()),
        ),
        ac_release_samples_raw=hl.agg.filter(
            mt.meta.release, hl.agg.sum(mt.GT.n_alt_alleles())
        ),
        ac_qc_samples_adj=hl.agg.filter(mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_unrelated_adj=hl.agg.filter(
            ~mt.meta.sample_filters.all_samples_related & mt.adj,
            hl.agg.sum(mt.GT.n_alt_alleles()),
        ),
        ac_release_samples_adj=hl.agg.filter(
            mt.meta.release & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())
        ),
    )
    return mt.rows()


def generate_fam_stats(mt: hl.MatrixTable, fam_file: str) -> hl.Table:
    """
    Calculate transmission and de novo mutation statistics using trios in the dataset.

    :param mt: Input MatrixTable
    :param fam_file: path to text file containing trio pedigree
    :return: Table containing trio stats
    """
    # Load Pedigree data and filter MT to samples present in any of the trios
    ped = hl.Pedigree.read(fam_file, delimiter="\t")
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    fam_ht = fam_ht.annotate(fam_members=[fam_ht.id, fam_ht.pat_id, fam_ht.mat_id])
    fam_ht = fam_ht.explode("fam_members", name="s")
    fam_ht = fam_ht.key_by("s").select().distinct()

    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.col_key]))
    logger.info(
        f"Generating family stats using {mt.count_cols()} samples from"
        f" {len(ped.trios)} trios."
    )

    mt = filter_to_autosomes(mt)
    mt = annotate_adj(mt)
    mt = mt.select_entries("GT", "GQ", "AD", "END", "adj")
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mt = hl.trio_matrix(mt, pedigree=ped, complete_trios=True)
    trio_adj = mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj

    ht = mt.select_rows(
        **generate_trio_stats_expr(
            mt,
            transmitted_strata={"raw": True, "adj": trio_adj},
            de_novo_strata={
                "raw": True,
                "adj": trio_adj,
            },
            proband_is_female_expr=mt.is_female,
        )
    ).rows()

    return ht.filter(
        ht.n_de_novos_raw + ht.n_transmitted_raw + ht.n_untransmitted_raw > 0
    )


def export_transmitted_singletons_vcf():
    """
    Export the transmitted singleton Table to a VCF.

    :return: None
    """
    qc_ac_ht = qc_ac.ht()

    for transmission_confidence in ["raw", "adj"]:
        ts_ht = qc_ac_ht.filter(
            (
                fam_stats.ht()[qc_ac_ht.key][f"n_transmitted_{transmission_confidence}"]
                == 1
            )
            & (qc_ac_ht.ac_qc_samples_raw == 2)
        )

        ts_ht = ts_ht.annotate(s=hl.null(hl.tstr))

        ts_mt = ts_ht.to_matrix_table_row_major(columns=["s"], entry_field_name="s")
        ts_mt = ts_mt.filter_cols(False)
        hl.export_vcf(
            ts_mt,
            get_transmitted_singleton_vcf_path(transmission_confidence),
            tabix=True,
        )


def run_vep(vep_version: str = "101") -> hl.Table:
    """
    Return a table with a VEP annotation for each variant in the raw MatrixTable.

    :param vep_version: Version of VEPed context Table to use in `vep_or_lookup_vep`
    :return: VEPed Table
    """
    ht = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, remove_hard_filtered_samples=False
    ).rows()
    ht = ht.filter(hl.len(ht.alleles) > 1)
    ht = hl.split_multi_hts(ht)
    ht = vep_or_lookup_vep(ht, vep_version=vep_version)
    ht = ht.annotate_globals(version=f"v{vep_version}")

    return ht


def main(args):  # noqa: D103
    hl.init(default_reference="GRCh38", log="/qc_annotations.log")
    mt = get_gnomad_v4_vds(test=args.test, high_quality_only=True).variant_data
    if args.compute_info:
        # TODO: is there any reason to also compute info per platform?
        default_compute_info(mt, site_annotations=True).write(
            get_info(split=False).path, overwrite=args.overwrite
        )

    if args.split_info:
        split_info().write(get_info(split=True).path, overwrite=args.overwrite)

    if args.export_info_vcf:
        info_ht = get_info(split=False).ht()
        hl.export_vcf(adjust_vcf_incompatible_types(info_ht), info_vcf_path())

    if args.generate_allele_data:
        mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
        generate_allele_data(mt.rows()).write(
            allele_data.path, overwrite=args.overwrite
        )

    if args.generate_ac:  # TODO: compute AC and qc_AC as part of compute_info
        mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True, samples_meta=True)
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

        ht = generate_ac(mt).checkpoint(
            "gs://gnomad-tmp/ac_tmp.ht",
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        ht.repartition(10000, shuffle=False).write(qc_ac.path, overwrite=args.overwrite)

    if args.generate_fam_stats:
        mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True, samples_meta=True)
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
        fam_stats_ht = generate_fam_stats(mt, trios.path)
        fam_stats_ht = fam_stats_ht.checkpoint(
            "gs://gnomad-tmp/fam_stats_tmp.ht",
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        fam_stats_ht = fam_stats_ht.repartition(10000, shuffle=False)
        fam_stats_ht.write(fam_stats.path, overwrite=args.overwrite)

    if args.export_transmitted_singletons_vcf:
        export_transmitted_singletons_vcf()

    if args.vep:
        run_vep(vep_version=args.vep_version).write(vep.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")

    parser.add_argument("--compute-info", help="Computes info HT.", action="store_true")
    parser.add_argument("--split-info", help="Splits info HT.", action="store_true")
    parser.add_argument(
        "--export-info-vcf", help="Export info as VCF.", action="store_true"
    )
    parser.add_argument(
        "--generate-allele-data", help="Calculates allele data.", action="store_true"
    )
    parser.add_argument(
        "--generate-ac",
        help=(
            "Creates a table with ACs for QC, unrelated QC and release samples (raw and"
            " adj)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--generate-fam-stats",
        help="Creates a table with transmitted allele counts and de novo counts.",
        action="store_true",
    )
    parser.add_argument(
        "--export-transmitted-singletons-vcf",
        help="Exports transmitted singletons to VCF files.",
        action="store_true",
    )
    parser.add_argument("--vep", help="Generates vep annotations.", action="store_true")
    parser.add_argument(
        "--vep-version",
        help="Version of VEPed context Table to use in vep_or_lookup_vep.",
        action="store_true",
        default="105",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
