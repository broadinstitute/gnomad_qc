from gnomad_hail import *
from gnomad_hail.resources.sample_qc import *

def get_f_ht(mt: hl.MatrixTable, by_data_type: bool = False, others='oth', min_af=0.001):
    if by_data_type:
        mt = mt.annotate_rows(
            call_stats=hl.agg.group_by(
                mt.data_type,
                hl.agg.group_by(
                    mt.pop,
                    hl.agg.filter(~mt.related & (hl.len(mt.pop_platform_filters) == 0), hl.agg.call_stats(mt.GT, mt.alleles))
                )
            )
        )

        global_af = mt.call_stats.map_values(lambda x: hl.sum(x.values().map(lambda y: y.AC[1])) / hl.sum(x.values().map(lambda y: y.AN)))

        return mt.select_cols(
            pop_f=hl.agg.filter(
                (mt.pop == others) | (mt.call_stats[mt.data_type][mt.pop].AF[1] > min_af),
                hl.agg.inbreeding(mt.GT, hl.cond(mt.pop == others, global_af[mt.data_type], mt.call_stats[mt.data_type][mt.pop].AF[1]))
            ),
            global_f=hl.agg.inbreeding(mt.GT, global_af[mt.data_type])
        ).cols()

    else:
        mt = mt.annotate_rows(
            call_stats= hl.agg.group_by(
                    mt.pop,
                    hl.agg.filter(~mt.related & (hl.len(mt.pop_platform_filters) == 0), hl.agg.call_stats(mt.GT, mt.alleles))
                )
            )

        global_af = hl.sum(mt.call_stats.values().map(lambda y: y.AC[1])) / hl.sum(mt.call_stats.values().map(lambda y: y.AN))

        return mt.select_cols(
            pop_f = hl.agg.filter(
                (mt.pop == 'oth') | (mt.call_stats[mt.pop].AF[1] > min_af),
                hl.agg.inbreeding(mt.GT, hl.cond(mt.pop == 'oth', global_af, mt.call_stats[mt.pop].AF[1]))
            ),
            global_f = hl.agg.inbreeding(mt.GT, global_af)
        ).cols().repartition(10, shuffle=False)


def main(args):

    if args.adj:
        hts = []
        for data_type in ['exomes', 'genomes']:
            mt = hl.read_matrix_table(get_gnomad_data_path(data_type, hardcalls=True, split=True))
            if args.purcell_5k:
                p5k = hl.import_locus_intervals(purcell5k_intervals_path)
                p5k = p5k.key_by(interval=hl.interval(p5k.interval.start, p5k.interval.end, includes_end=True))
                mt = mt.filter_rows(hl.is_defined(p5k[mt.locus]))
                print(f"Found {mt.count_rows()} sites overlapping with purcell 5k")
            else:
                sites = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True)).rows().select()
                mt = mt.filter_rows(hl.is_defined(sites[mt.row_key]))

            hard_filters_ht = hl.read_table(qc_ht_path(data_type, 'hard_filters'))
            meta = hl.read_table(qc_ht_path(data_type, 'pop_platform')).select('pop', 'related', 'pop_platform_filters')
            mt = mt.key_cols_by(data_type=data_type, s=mt.s)
            mt = mt.filter_cols(hl.len(hard_filters_ht[mt.col_key].hard_filters) == 0)
            mt = mt.annotate_cols(**meta[mt.col_key])

            if args.filter_release:
                release = get_gnomad_public_data(data_type)
                release = release.filter(hl.len(release.filters) == 0)
                mt = mt.filter_rows(hl.is_defined(release[mt.row_key]))
            elif args.filter_inbreeding:
                mt = mt.filter_rows(mt.info.InbreedingCoeff > -0.3)

            if args.eur_subpops:
                subpop_ht = hl.read_table(subpop_ht_path('eur')).select('subpop')
                mt = mt.annotate_cols(pop=subpop_ht[mt.col_key].subpop)
                mt = mt.filter_cols(hl.is_defined(mt.pop))

            # mt = mt.filter_rows(hl.agg.fraction(mt.adj) > 0.8)
            mt = mt.filter_entries(mt.adj)

            mt = unphase_mt(mt)
            hts.append(get_f_ht(mt, False, others='oeu' if args.eur_subpops else 'oth', min_af=args.min_af).persist())

        ht = hts[0].union(hts[1])


    else:
        mt = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True))
        exomes_meta = hl.read_table(qc_ht_path('exomes', 'pop_platform'))
        exomes_meta = exomes_meta.key_by(data_type='exomes', s=exomes_meta.s).select('pop', 'related', 'pop_platform_filters')
        genomes_meta = hl.read_table(qc_ht_path('genomes', 'pop_platform'))
        genomes_meta = genomes_meta.key_by(data_type='genomes', s=genomes_meta.s).select('pop', 'related', 'pop_platform_filters')
        meta = exomes_meta.union(genomes_meta)
        mt = mt.annotate_cols(**meta[mt.col_key])

        if args.purcell_5k:
            p5k = hl.import_locus_intervals(purcell5k_intervals_path)
            mt = mt.filter_rows(hl.is_defined(p5k[mt.locus]))

        if args.filter_release:
            exomes_release = get_gnomad_public_data('exomes')
            exomes_release = exomes_release.filter(hl.len(exomes_release.filters)==0)
            genomes_release = get_gnomad_public_data('genomes')
            genomes_release = genomes_release.filter(hl.len(genomes_release.filters) == 0)
            release = exomes_release.select().union(genomes_release.select())
            mt = mt.filter_rows(hl.is_defined(release[mt.row_key]))
        elif args.filter_inbreeding:
            mt = mt.filter_rows(mt.info.InbreedingCoeff > -0.3)

        if args.eur_subpops:
            subpop_ht = hl.read_table(subpop_ht_path('eur')).select('subpop')
            mt = mt.annotate_cols(pop=subpop_ht[mt.col_key].subpop)
            mt = mt.filter_cols(hl.is_defined(mt.pop))

        ht = get_f_ht(mt, args.af_by_data_type, others='oeu' if args.eur_subpops else 'oth', min_af=args.min_af)

    ht.write(
        'gs://gnomad-tmp/test_f_{}{}{}{}{}{}.ht'.format(
            args.min_af,
            '_p5k' if args.purcell_5k else '',
            '_adj' if args.adj else '',
            '_by_data_type' if args.af_by_data_type else '',
            '_release' if args.filter_release else '_filtered' if args.filter_inbreeding else '',
            '_subpops' if args.eur_subpops else ''
        ),
        overwrite=args.overwrite
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--af_by_data_type', help='Compute AF on exomes and genomes separately (always the case if --adj).', action='store_true')
    parser.add_argument('--filter_inbreeding', help='Applies inbreeding coeff sites filtering.', action='store_true')
    parser.add_argument('--filter_release', help='Only compute on release sites.', action='store_true')
    parser.add_argument('--adj', help='Only compute on adj genotypes.', action='store_true')
    parser.add_argument('--min_af', help='Minimum AF.', default = 0.001, type=float)
    parser.add_argument('--purcell_5k', help='Uses purcell 5k SNPs only.', action='store_true')
    parser.add_argument('--eur_subpops', help='Generates F for eur subpops.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)

