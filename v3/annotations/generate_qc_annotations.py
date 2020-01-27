from gnomad_hail import *
from gnomad_hail.utils.sparse_mt import *
from gnomad_hail.utils.variant_qc import *
from v3.resources import *
import argparse


def compute_info() -> hl.Table:
    """
    Computes a HT with the typical GATK AS and site-level info fields
    as well as ACs and lowqual fields.
    Note that this table doesn't split multi-allelic sites.

    :return: Table with info fields
    :rtype: Table
    """
    mt = get_full_mt(split=False, key_by_locus_and_alleles=True, remove_hard_filtered_samples=False)
    mt = mt.filter_rows((hl.len(mt.alleles) > 1))
    mt = mt.transmute_entries(**mt.gvcf_info)

    # Compute AS and site level info expr
    # Note that production defaults have changed:
    # For new releases, the `RAWMQ_andDP` field replaces the `RAW_MQ` and `MQ_DP` fields
    info_expr = get_site_info_expr(
        mt,
        sum_agg_fields=INFO_SUM_AGG_FIELDS + ['RAW_MQ'],
        int32_sum_agg_fields=INFO_INT32_SUM_AGG_FIELDS  + ['MQ_DP'],
        array_sum_agg_fields=['SB']
    )
    info_expr = info_expr.annotate(
        **get_as_info_expr(
            mt,
            sum_agg_fields=INFO_SUM_AGG_FIELDS + ['RAW_MQ'],
            int32_sum_agg_fields=INFO_INT32_SUM_AGG_FIELDS + ['MQ_DP'],
            array_sum_agg_fields=['SB']
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
                    mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[mt.LA.index(ai)]
                )
            )
        ),
        hl.range(1, hl.len(mt.alleles))
    )

    # Then, for each non-ref allele, compute
    # AC as the adj group
    # AC_raw as the sum of adj and non-adj groups
    info_expr = info_expr.annotate(
        AC_raw=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0) + i.get(False, 0))),
        AC=grp_ac_expr.map(lambda i: hl.int32(i.get(True, 0)))
    )

    info_ht = mt.select_rows(
        info=info_expr
    ).rows()

    # Add lowqual flag
    info_ht = info_ht.annotate(
        lowqual=get_lowqual_expr(
            info_ht.alleles,
            info_ht.info.QUALapprox,
            # The indel het prior used for gnomad v3 was 1/10k bases (phred=40).
            # This value is usually 1/8k bases (phred=39).
            indel_phred_het_prior=40
        )
    )

    return info_ht.naive_coalesce(5000)


def split_info() -> hl.Table:
    """
    Generates an info table that splits multi-allelic sites from
    the multi-allelic info table.

    :return: Info table with split multi-allelics
    :rtype: Table
    """
    info_ht = hl.read_table(get_info_ht_path(split=False))

    # Create split version
    info_ht = hl.split_multi(info_ht)

    # Index AS annotations
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **{f: info_ht.info[f][info_ht.a_index - 1] for f in info_ht.info if f.startswith("AC") or (f.startswith("AS_") and not f == 'AS_SB_TABLE')},
            AS_SB_TABLE=info_ht.info.AS_SB_TABLE[0].extend(info_ht.info.AS_SB_TABLE[info_ht.a_index])
        ),
        lowqual=info_ht.lowqual[info_ht.a_index - 1]
    )
    return info_ht


def generate_ac(mt: hl.MatrixTable, fam_file: str) -> hl.Table:
    """
    Creates Table with QC samples, QC samples removing children and release samples raw and adj ACs.
    """
    mt = mt.filter_cols(mt.meta.high_quality)
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    mt = mt.annotate_cols(unrelated_sample=hl.is_missing(fam_ht[mt.s]))
    mt = mt.filter_rows(hl.len(mt.alleles)>1)
    mt = annotate_adj(mt)
    mt = mt.annotate_rows(
        ac_qc_samples_raw=hl.agg.sum(mt.GT.n_alt_alleles()),
        ac_qc_samples_unrelated_raw=hl.agg.filter(~mt.meta.all_samples_related, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_release_samples_raw=hl.agg.filter(mt.meta.release, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_adj=hl.agg.filter(mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_qc_samples_unrelated_adj=hl.agg.filter(~mt.meta.all_samples_related & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
        ac_release_samples_adj=hl.agg.filter(mt.meta.release & mt.adj, hl.agg.sum(mt.GT.n_alt_alleles())),
    )
    return mt.rows()


def generate_fam_stats(
        mt: hl.MatrixTable,
        fam_file: str
) -> hl.Table:
    # Load Pedigree data and filter MT to samples present in any of the trios
    ped = hl.Pedigree.read(fam_file, delimiter="\t")
    fam_ht = hl.import_fam(fam_file, delimiter="\t")
    fam_ht = fam_ht.annotate(
        fam_members=[fam_ht.id, fam_ht.pat_id, fam_ht.mat_id]
    )
    fam_ht = fam_ht.explode('fam_members', name='s')
    fam_ht = fam_ht.key_by('s').select().distinct()

    mt = mt.filter_cols(hl.is_defined(fam_ht[mt.col_key]))
    logger.info(f"Generating family stats using {mt.count_cols()} samples from {len(ped.trios)} trios.")

    mt = filter_to_autosomes(mt)
    mt = annotate_adj(mt)
    mt = mt.select_entries('GT', 'GQ', 'AD', 'END', 'adj')
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows(hl.len(mt.alleles)==2)
    mt = hl.trio_matrix(mt, pedigree=ped, complete_trios=True)
    trio_adj = (mt.proband_entry.adj & mt.father_entry.adj & mt.mother_entry.adj)
    parents_no_alt = (mt.mother_entry.AD[1] == 0) & (mt.father_entry.AD[1] == 0)
    parents_high_depth = (mt.mother_entry.AD[0] + mt.mother_entry.AD[1] > 20) & (mt.father_entry.AD[0] + mt.father_entry.AD[1] > 20)
    parents_high_gq = (mt.mother_entry.GQ >= 30) & (mt.father_entry.GQ >= 30)

    ht = mt.select_rows(
        **generate_fam_stats_expr(
            mt,
            transmitted_strata={
                'raw':  None,
                'adj': trio_adj
            },
            de_novo_strata={
                'raw': None,
                'adj': trio_adj,
                'hq':  trio_adj &  parents_high_gq & parents_high_depth & parents_no_alt
            },
            proband_is_female_expr=mt.is_female
        )
    ).rows()

    return ht.filter(
        ht.n_de_novos_raw + ht.n_transmitted_raw + ht.n_untransmitted_raw > 0
    )


def export_transmitted_singletons_vcf():
    qc_ac_ht = hl.read_table(ac_ht_path)
    fam_stats = hl.read_table(fam_stats_ht_path)

    for transmission_confidence in ['raw', 'adj']:
        ts_ht = qc_ac_ht.filter(
            (fam_stats[qc_ac_ht.key][f'n_transmitted_{transmission_confidence}'] == 1) &
            (qc_ac_ht.ac_qc_samples_raw == 2)
        )

        ts_ht = ts_ht.annotate(
            s=hl.null(hl.tstr)
        )

        ts_mt = ts_ht.to_matrix_table_row_major(columns=['s'], entry_field_name='s')
        ts_mt = ts_mt.filter_cols(False)
        hl.export_vcf(ts_mt, get_transmitted_singleton_vcf_path(transmission_confidence))


def run_vep() -> hl.Table:

    def get_mt_partitions(mt_path: str) -> List[hl.Interval]:
        """
        This function loads the partitioning from a given MT.
        Note that because it relies on hardcoded paths within the MT that are still in flux,
        it isn't guaranteed to work on future versions of the MT format.

        :param str mt_path: MT path
        :return: MT partitions
        :rtype: List of Interval
        """
        logger.info(f'Reading partitions for {mt_path}')
        import json
        from os import path
        mt = hl.read_matrix_table(mt_path)
        with hl.hadoop_open(path.join(mt_path, 'rows', 'rows', 'metadata.json.gz')) as f:
            intervals_json = json.load(f)['jRangeBounds']
            return hl.tarray(hl.tinterval(hl.tstruct(locus=mt.locus.dtype)))._convert_from_json(intervals_json)

    # Load VEP table with the same partitioning as the gnomAD v3 MT
    intervals = get_mt_partitions(get_full_mt_path())
    vep_ht = hl.read_table(synthetic_grch38_vep_ht_path, _intervals=intervals)

    ht = get_full_mt(key_by_locus_and_alleles=True).rows()
    ht = ht.filter(hl.len(ht.alleles)>1)

    return vep_or_lookup_vep(ht, reference_vep_ht=vep_ht, reference='GRCh38')


def main(args):
    hl.init(default_reference='GRCh38', log='/qc_annotations.log')

    if args.compute_info:
        compute_info().write(get_info_ht_path(split=False), overwrite=args.overwrite)

    if args.split_info:
        split_info().write(get_info_ht_path(split=True), overwrite=args.overwrite)

    if args.export_info_vcf:
        info_ht = hl.read_table(get_info_ht_path(split=False))
        hl.export_vcf(ht_to_vcf_mt(info_ht), info_vcf_path)

    if args.generate_ac:
        mt = get_full_mt(split=True, key_by_locus_and_alleles=True)
        meta_ht = hl.read_table(meta_ht_path)
        mt = mt.annotate_cols(meta=meta_ht[mt.col_key])
        ht = generate_ac(mt, ped_path).checkpoint('gs://gnomad-tmp/v3_ac_tmp.ht', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        ht.repartition(10000, shuffle=False).write(ac_ht_path, overwrite=args.overwrite)

    if args.generate_fam_stats:
        mt = get_full_mt(split=True, key_by_locus_and_alleles=True)
        meta_ht = hl.read_table(meta_ht_path)
        mt = mt.annotate_cols(meta=meta_ht[mt.col_key])
        fam_stats_ht = generate_fam_stats(mt, trios_path)
        fam_stats_ht = fam_stats_ht.checkpoint('gs://gnomad-tmp/v3_fam_stats_tmp.ht', overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        fam_stats_ht = fam_stats_ht.repartition(10000, shuffle=False)
        fam_stats_ht.write(fam_stats_ht_path, overwrite=args.overwrite)

    if args.export_transmitted_singletons_vcf:
        export_transmitted_singletons_vcf()

    if args.vep:
        run_vep().write(vep_ht_path, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--compute_info', help='Computes info HT', action='store_true')
    parser.add_argument('--split_info', help='Splits info HT', action='store_true')
    parser.add_argument('--export_info_vcf', help='Export info as VCF', action='store_true')
    parser.add_argument('--generate_ac', help='Creates a table with ACs for QC, unrelated QC and release samples (raw and adj)', action='store_true')
    parser.add_argument('--generate_fam_stats', help='Creates a table with transmitted allele counts and de novo counts.', action='store_true')
    parser.add_argument('--vep', help='Generates vep annotations.', action='store_true')
    parser.add_argument('--export_transmitted_singletons_vcf', help='Exports transmitted singletons to VCF files.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
