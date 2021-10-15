import argparse
import sys
from gnomad.utils.annotations import unphase_call_expr, add_variant_type
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.slack import slack_notifications
from gnomad.utils.file_utils import write_temp_gcs
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources import *


def generate_allele_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - allele_data (nonsplit_alleles, has_star, variant_type, and n_alt_alleles)

    :param MatrixTable mt: Full unsplit MT
    :return: Table with allele data annotations
    :rtype: Table
    """
    ht = mt.rows().select()
    allele_data = hl.struct(nonsplit_alleles=ht.alleles,
                            has_star=hl.any(lambda a: a == '*', ht.alleles))
    ht = ht.annotate(allele_data=allele_data.annotate(**add_variant_type(ht.alleles)))

    ht = hl.split_multi_hts(ht)
    allele_type = (hl.case()
                   .when(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv')
                   .when(hl.is_insertion(ht.alleles[0], ht.alleles[1]), 'ins')
                   .when(hl.is_deletion(ht.alleles[0], ht.alleles[1]), 'del')
                   .default('complex')
                   )
    ht = ht.annotate(allele_data=ht.allele_data.annotate(allele_type=allele_type,
                                                              was_mixed=ht.allele_data.variant_type == 'mixed'))
    return ht


def generate_call_stats(mt: hl.MatrixTable) -> hl.Table:
    """
    Add qc_callstats for 3 categories (high-quality samples, release, and all in VDS)
    """
    sample_group_filters = {
        "qc_samples_raw": mt.meta.high_quality,
        "release_samples_raw": mt.meta.release,
        "all_samples_raw": True
    }
    mt = mt.select_cols(**sample_group_filters)
    mt = mt.annotate_entries(GT=unphase_call_expr(mt.GT))
    mt = mt.select_rows()
    call_stats_expression = []
    for group in sample_group_filters.keys():
        callstats = hl.agg.filter(mt[group], hl.agg.call_stats(mt.GT, mt.alleles))
        call_stats_expression.append(callstats.annotate(meta={'group': group}))

    return mt.annotate_rows(qc_callstats=call_stats_expression).rows()


def generate_qc_annotations(mt: hl.MatrixTable, all_annotations: bool = True, medians: bool = True) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - qc_stats (allele-specific statistics and call_stats for sample groups)
    Calculates stats for: gq, dp, nrq, ab, p(ab)
    Also, calculates best ab, nrdp, qual, combined p(ab), and qd
    Optionally, adds medians for gq, dp, nrq, AB and p(ab)

    :param MatrixTable mt: Full MT
    :param bool all_annotations: Whether to compute all annotations possible
    :param bool medians: Whether to compute the allele-specific stats medians
    :return: Table with qc annotations
    :rtype: Table
    """
    sample_group_filters = {
        "qc_samples_raw": mt.meta.high_quality,
    }
    mt = mt.select_cols(**sample_group_filters)
    mt = mt.select_rows()

    lin = 10 ** (-mt.PL / 10.0)
    gp0 = lin[0] / hl.sum(lin)
    # gp0 = hl.bind(lambda lin: lin[0] / hl.sum(lin), 10 ** (-mt.PL / 10.0))
    pab = hl.binom_test(mt.AD[1], hl.sum(mt.AD), 0.5, "two.sided")

    qc_stats_expression = []
    for group in sample_group_filters.keys():
        criteria = mt[group]
        hets = criteria & mt.GT.is_het()
        non_refs = criteria & mt.GT.is_non_ref()
        nrq = -1*hl.log10(gp0)

        het_length = hl.agg.count_where(hets)

        qual_agg = -10 * hl.agg.filter(non_refs, hl.agg.sum(hl.cond(
            mt.PL[0] > 3000, -300, hl.log10(gp0)
        )))

        stats_expression = hl.struct(
            pab=hl.agg.filter(hets, hl.agg.stats(pab)),
            qd=hl.min(35, qual_agg / hl.agg.filter(non_refs, hl.agg.sum(mt.DP))),
            meta={'group': group})

        if all_annotations:
            stats_expression = stats_expression.annotate(
                gq=hl.agg.filter(non_refs, hl.agg.stats(mt.GQ)),
                dp=hl.agg.filter(non_refs, hl.agg.stats(mt.DP)),
                nrq=hl.agg.filter(non_refs, hl.agg.stats(nrq)),
                ab=hl.agg.filter(hets, hl.agg.stats(mt.AD[1] / hl.sum(mt.AD))),
                best_ab=hl.agg.filter(hets, hl.agg.min(hl.abs(mt.AD[1] / mt.DP - 0.5))),
                nrdp=hl.agg.filter(non_refs, hl.agg.sum(mt.DP)),
                qual=qual_agg,
                combined_pab=hl.or_missing(
                    het_length == 0,
                    -10 * hl.log10(hl.pchisqtail(-2 * hl.agg.filter(hets, hl.agg.sum(hl.log10(pab))), 2 * het_length))
                )
            )
        if medians:
            stats_expression = stats_expression.annotate(
                gq_median=hl.median(hl.agg.filter(non_refs, hl.agg.collect(mt.GQ))),
                dp_median=hl.median(hl.agg.filter(non_refs, hl.agg.collect(mt.DP))),
                nrq_median=hl.median(hl.agg.filter(non_refs, hl.agg.collect(nrq))),
                ab_median=hl.median(hl.agg.filter(hets, hl.agg.collect(mt.AD[1] / hl.sum(mt.AD))))
                # pab_median=hl.median(hl.agg.filter(hets, hl.agg.collect(pab)))
            )
        qc_stats_expression.append(stats_expression)
    return mt.annotate_rows(qc_stats=qc_stats_expression).rows()


def generate_qual_hists(mt: hl.MatrixTable) -> hl.Table:
    mt = hl.split_multi_hts(mt)
    return mt.annotate_rows(
        gq_hist_alt=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.GQ, 0, 100, 20)),
        gq_hist_all=hl.agg.hist(mt.GQ, 0, 100, 20),
        dp_hist_alt=hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.DP, 0, 100, 20)),
        dp_hist_all=hl.agg.hist(mt.DP, 0, 100, 20),
        ab_hist_alt=hl.agg.filter(mt.GT.is_het(), hl.agg.hist(mt.AD[1] / hl.sum(mt.AD), 0, 1, 20))
    ).rows()


def generate_family_stats(mt: hl.MatrixTable, fam_file: str, calculate_adj: bool = False) -> Tuple[hl.Table, hl.Table]:
    """
    Writes bi-allelic sites MT with the following annotations:
     - family_stats (TDT, Mendel Errors, AC_unrelated_qc)
     - truth_data (presence in Omni, HapMap, 1KG/TGP high conf SNVs, Mills)

    :param MatrixTable mt: Full MT
    :param str fam_file: Fam pedigree file location
    :param bool calculate_adj: Whether to also calculate family metrics for adj genotypes
    :return: Table with qc annotations
    :rtype: Table
    """
    mt = mt.select_cols(high_quality=mt.meta.high_quality)
    mt = mt.select_rows()
    mt = annotate_unrelated_sample(mt, fam_file)

    # Unphased for now, since mendel_errors does not support phased alleles
    mt = mt.annotate_entries(GT=unphase_call_expr(mt.GT))
    ped = hl.Pedigree.read(fam_file, delimiter='\\t')
    family_stats_struct, family_stats_sample_ht = family_stats(mt, ped, 'raw')
    mt = mt.annotate_rows(family_stats=[family_stats_struct])

    if calculate_adj:
        mt = filter_to_adj(mt)
        adj_family_stats_struct, adj_family_stats_sample_ht = family_stats(mt, ped, 'adj')

        family_stats_sample_ht = family_stats_sample_ht.annotate(
            adj=adj_family_stats_sample_ht[family_stats_sample_ht.s])

        mt = mt.annotate_rows(family_stats=mt.family_stats.append(adj_family_stats_struct))

    return mt.rows(), family_stats_sample_ht


def family_stats(mt: hl.MatrixTable, ped: hl.Pedigree, group_name: str) -> Tuple[hl.expr.StructExpression, hl.Table]:
    tdt_table = hl.transmission_disequilibrium_test(mt, ped)
    _, _, per_sample, per_variant = hl.mendel_errors(mt.GT, ped)
    family_stats_struct = hl.struct(mendel=per_variant[mt.row_key],
                                    tdt=tdt_table[mt.row_key],
                                    unrelated_qc_callstats=hl.agg.filter(mt.high_quality & mt.unrelated_sample,
                                                                         hl.agg.call_stats(mt.GT, mt.alleles)),
                                    meta={'group': group_name})
    return family_stats_struct, per_sample


def read_fam(fam_file: str) -> hl.Table:
    columns = ['fam_id', 's', 'pat_id', 'mat_id', 'is_female']
    return hl.import_table(fam_file, no_header=True).rename({f'f{i}': c for i, c in enumerate(columns)}).key_by('s')


def annotate_unrelated_sample(mt: hl.MatrixTable, fam_file: str) -> hl.MatrixTable:
    fam_ht = read_fam(fam_file)
    return mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))


def generate_de_novos(mt: hl.MatrixTable, fam_file: str, freq_data: hl.Table) -> hl.Table:
    mt = mt.select_cols()
    fam_ht = read_fam(fam_file).key_by()
    fam_ht = fam_ht.select(s=[fam_ht.s, fam_ht.pat_id, fam_ht.mat_id]).explode('s').key_by('s')
    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.s]))
    mt = mt.select_rows()
    mt = hl.split_multi_hts(mt)
    mt = mt.annotate_rows(family_stats=freq_data[mt.row_key].family_stats)
    ped = hl.Pedigree.read(fam_file, delimiter='\\t')

    de_novo_table = hl.de_novo(mt, ped, mt.family_stats[0].unrelated_qc_callstats.AF[1])
    de_novo_table = de_novo_table.key_by('locus', 'alleles').collect_by_key('de_novo_data')

    return de_novo_table


def annotate_truth_data(mt: hl.MatrixTable) -> hl.Table:
    """
    Writes bi-allelic sites MT with the following annotations:
     - truth_data (presence in Omni, HapMap, 1KG/TGP high conf SNVs, Mills)

    :param MatrixTable mt: Full MT
    :return: Table with qc annotations
    :rtype: Table
    """
    mt = mt.select_rows()

    truth_mtes = {
        'hapmap': hapmap_mt_path(),
        'omni': omni_mt_path(),
        'mills': mills_mt_path(),
        'kgp_high_conf_snvs': kgp_high_conf_snvs_mt_path()
    }
    truth_mtes = {key: hl.split_multi_hts(hl.read_matrix_table(path).repartition(1000), left_aligned=False)
                  for key, path in truth_mtes.items()}

    return mt.annotate_rows(truth_data=hl.struct(**{root: hl.is_defined(truth_mt[mt.row_key, :])
                                                    for root, truth_mt in truth_mtes.items()})).rows()


def main(args):
    hl.init(log='/qc_annotations.log')

    data_type = 'genomes' if args.genomes else 'exomes'

    if args.vep:  # CPU-hours: 250 (E), 600 (G)
        mt = get_gnomad_data(data_type).rows().select()
        hl.vep(mt, vep_config).write(annotations_ht_path(data_type, 'vep'), args.overwrite)

        mt = get_gnomad_data(data_type).rows().select()
        hl.vep(mt, vep_config, csq=True).write(annotations_ht_path(data_type, 'vep_csq'), args.overwrite)

    if args.generate_allele_data:  # CPU-hours: 100 (E), 200 (G)
        mt = get_gnomad_data(data_type, split=False)
        generate_allele_data(mt).write(annotations_ht_path(data_type, 'allele_data'), args.overwrite)

    if args.generate_qc_annotations:  # CPU-hours: 750 (E), 2400 (G) - ~1 individual task each runs long
        # Turn on spark speculation
        mt = get_gnomad_data(data_type, non_refs_only=True)
        mt = generate_qc_annotations(mt, all_annotations=args.calculate_all_annotations, medians=args.calculate_medians)
        write_temp_gcs(mt, annotations_ht_path(data_type, 'qc_stats'), args.overwrite)

    if args.generate_qual_hists:  # CPU-hours: 4000 (E), 8000 (G)
        mt = get_gnomad_data(data_type, raw=True, split=False, release_samples=True)
        ht = generate_qual_hists(mt)
        write_temp_gcs(ht, annotations_ht_path(data_type, 'qual_hists'), args.overwrite)

    if args.generate_call_stats:  # CPU-hours: 750 (E), 1500 (G)
        mt = get_gnomad_data(data_type)
        generate_call_stats(mt).write(annotations_ht_path(data_type, 'call_stats'), args.overwrite)

    if args.generate_family_stats:  # CPU-hours: 8K (E), 13K (G)
        mt = get_gnomad_data(data_type)
        mt, sample_table = generate_family_stats(mt, fam_path(data_type, true_trios=True), args.include_adj_family_stats)
        write_temp_gcs(mt, annotations_ht_path(data_type, 'family_stats'), args.overwrite)
        write_temp_gcs(sample_table, sample_annotations_table_path(data_type, 'family_stats'), args.overwrite)

    if args.generate_de_novos:  # (2.2 min/part @ 100K = 3K CPU-hours) + (7.4 m/p = 12K) + (34 m/p = ~44K) = 59K
        # Turn on spark speculation?
        mt = get_gnomad_data(data_type, raw=True, split=False)
        freq_data = hl.read_table(annotations_ht_path(data_type, 'family_stats'))
        mt = generate_de_novos(mt, fam_path(data_type), freq_data)
        mt.write(annotations_ht_path(data_type, 'de_novos'), args.overwrite)

    if args.annotate_truth_data:
        mt = get_gnomad_data(data_type, meta_root=None)
        annotate_truth_data(mt).write(annotations_ht_path(data_type, 'truth_data'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--vep', help='Runs VEP', action='store_true')
    parser.add_argument('--generate_allele_data', help='Calculates allele data', action='store_true')
    parser.add_argument('--generate_qc_annotations', help='Calculates QC annotations', action='store_true')
    parser.add_argument('--generate_qual_hists', help='Calculates GQ, DP, AB histograms per variant', action='store_true')
    parser.add_argument('--generate_call_stats', help='Calculates call stats', action='store_true')
    parser.add_argument('--generate_family_stats', help='Calculates family stats', action='store_true')
    parser.add_argument('--include_adj_family_stats', help='Also calculate family stats for adj genotypes', action='store_true')
    parser.add_argument('--generate_de_novos', help='Calculates de novo data', action='store_true')
    parser.add_argument('--annotate_truth_data', help='Annotates MT with truth data', action='store_true')
    parser.add_argument('--calculate_medians', help='Calculate metric medians (warning: slow)', action='store_true')
    parser.add_argument('--calculate_all_annotations', help='Calculation many more annotations (warning: slow)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
