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
    get_as_info_expr,
    get_site_info_expr,
    INFO_INT32_SUM_AGG_FIELDS,
    INFO_SUM_AGG_FIELDS,
    split_info_annotation,
    split_lowqual_annotation,
)
from gnomad.utils.vcf import adjust_vcf_incompatible_types
from gnomad.utils.vep import vep_or_lookup_vep

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import (
    allele_data,
    fam_stats,
    get_info,
    info_vcf_path,
    qc_ac,
    vep,
)
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import trios
from gnomad_qc.v3.resources.annotations import get_transmitted_singleton_vcf_path

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def compute_info() -> hl.Table:
    """
    Computes a HT with the typical GATK AS and site-level info fields as well as ACs and lowqual fields.

    Note that this table doesn't split multi-allelic sites.

    :return: Table with info fields
    :rtype: Table
    """
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, remove_hard_filtered_samples=False
    )

    mt = mt.filter_rows((hl.len(mt.alleles) > 1))
    mt = mt.transmute_entries(**mt.gvcf_info)
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    # Compute AS and site level info expr
    # Note that production defaults have changed:
    # For new releases, the `RAWMQ_andDP` field replaces the `RAW_MQ` and `MQ_DP` fields
    info_expr = get_site_info_expr(
        mt,
        sum_agg_fields=INFO_SUM_AGG_FIELDS + ["RAW_MQ"],
        int32_sum_agg_fields=INFO_INT32_SUM_AGG_FIELDS + ["MQ_DP"],
        array_sum_agg_fields=["SB"],
    )
    info_expr = info_expr.annotate(
        **get_as_info_expr(
            mt,
            sum_agg_fields=INFO_SUM_AGG_FIELDS + ["RAW_MQ"],
            int32_sum_agg_fields=INFO_INT32_SUM_AGG_FIELDS + ["MQ_DP"],
            array_sum_agg_fields=["SB"],
        )
    )

    # Add AC and AC_raw:
    # First compute ACs for each non-ref allele, grouped by adj
    grp_ac_expr = hl.agg.array_agg(
        lambda ai: hl.agg.filter(
            mt.LA.contains(ai),
            hl.agg.group_by(
                get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
                hl.agg.sum(
                    mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[
                        mt.LA.index(ai)
                    ]
                ),
            ),
        ),
        mt.alt_alleles_range_array,
    )

    # Then, for each non-ref allele, compute
    # AC as the adj group
    # AC_raw as the sum of adj and non-adj groups
    info_expr = info_expr.annotate(
        AC_raw=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))),
        AC=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0))),
    )

    # Annotating raw MT with pab max
    info_expr = info_expr.annotate(
        AS_pab_max=hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai) & mt.LGT.is_het(),
                hl.agg.max(
                    hl.binom_test(
                        mt.LAD[mt.LA.index(ai)], hl.sum(mt.LAD), 0.5, "two-sided"
                    )
                ),
            ),
            mt.alt_alleles_range_array,
        )
    )

    info_ht = mt.select_rows(info=info_expr).rows()

    # Add lowqual flag
    info_ht = info_ht.annotate(
        lowqual=get_lowqual_expr(
            info_ht.alleles,
            info_ht.info.QUALapprox,
            # The indel het prior used for gnomad v3 was 1/10k bases (phred=40).
            # This value is usually 1/8k bases (phred=39).
            indel_phred_het_prior=40,
        ),
        AS_lowqual=get_lowqual_expr(
            info_ht.alleles, info_ht.info.AS_QUALapprox, indel_phred_het_prior=40
        ),
    )

    return info_ht.naive_coalesce(7500)


def split_info() -> hl.Table:
    """
    Generates an info table that splits multi-allelic sites from the multi-allelic info table.

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
    Returns bi-allelic sites HT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)

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
    Creates Table containing allele counts per variant.

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
        f"Generating family stats using {mt.count_cols()} samples from {len(ped.trios)} trios."
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
            de_novo_strata={"raw": True, "adj": trio_adj,},
            proband_is_female_expr=mt.is_female,
        )
    ).rows()

    return ht.filter(
        ht.n_de_novos_raw + ht.n_transmitted_raw + ht.n_untransmitted_raw > 0
    )


def export_transmitted_singletons_vcf():
    """
    Exports the transmitted singleton Table to a VCF.

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
    Returns a table with a VEP annotation for each variant in the raw MatrixTable.

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


def main(args):
    hl.init(default_reference="GRCh38", log="/qc_annotations.log")

    if args.compute_info:
        compute_info().write(get_info(split=False).path, overwrite=args.overwrite)

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
    parser.add_argument("--compute_info", help="Computes info HT", action="store_true")
    parser.add_argument("--split_info", help="Splits info HT", action="store_true")
    parser.add_argument(
        "--export_info_vcf", help="Export info as VCF", action="store_true"
    )
    parser.add_argument(
        "--generate_allele_data", help="Calculates allele data", action="store_true"
    )
    parser.add_argument(
        "--generate_ac",
        help="Creates a table with ACs for QC, unrelated QC and release samples (raw and adj)",
        action="store_true",
    )
    parser.add_argument(
        "--generate_fam_stats",
        help="Creates a table with transmitted allele counts and de novo counts.",
        action="store_true",
    )
    parser.add_argument("--vep", help="Generates vep annotations.", action="store_true")
    parser.add_argument(
        "--vep_version",
        help="Version of VEPed context Table to use in vep_or_lookup_vep",
        action="store_true",
        default="101",
    )
    parser.add_argument(
        "--export_transmitted_singletons_vcf",
        help="Exports transmitted singletons to VCF files.",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
