from gnomad.utils.file_utils import write_temp_gcs
from gnomad.utils.slack import slack_notifications
from gnomad_qc.v2.resources import *
from gnomad_qc.slack_creds import slack_token
from collections import Counter
from gnomad.utils.annotations import pop_max_expr, project_max_expr
import argparse
import sys
import logging

logger = logging.getLogger("generate_frequency_data")

DOWNSAMPLINGS = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000, 50000,
                 55000, 60000, 65000, 70000, 75000, 80000, 85000, 90000, 95000, 100000, 110000, 120000]
POPS_TO_REMOVE_FOR_POPMAX = ['asj', 'fin', 'oth']
F_CUTOFF = 0.05


def generate_downsamplings_cumulative(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, List[int]]:
    pop_data = mt.meta.pop.collect()
    pops = Counter(pop_data)
    downsamplings = DOWNSAMPLINGS + list(pops.values())
    downsamplings = sorted([x for x in downsamplings if x <= sum(pops.values())])
    kt = mt.cols()
    kt = kt.annotate(r=hl.rand_unif(0, 1))
    kt = kt.order_by(kt.r).add_index('global_idx')

    for i, pop in enumerate(pops):
        pop_kt = kt.filter(kt.meta.pop == pop).add_index('pop_idx')
        if not i:
            global_kt = pop_kt
        else:
            global_kt = global_kt.union(pop_kt)
    global_kt = global_kt.key_by('s')
    return mt.annotate_cols(downsampling=global_kt[mt.s]), downsamplings


def add_faf_expr(freq: hl.expr.ArrayExpression, freq_meta: hl.expr.ArrayExpression, locus: hl.expr.LocusExpression, populations: Set[str]) -> hl.expr.ArrayExpression:
    """
    Calculates popmax (add an additional entry into freq with popmax: pop)

    :param ArrayExpression freq: ArrayExpression of Structs with ['ac', 'an', 'hom']
    :param ArrayExpression freq_meta: ArrayExpression of meta dictionaries corresponding to freq
    :param LocusExpression locus: LocusExpression
    :param set of str populations: Set of populations over which to calculate popmax
    :return: Frequency data with annotated popmax
    :rtype: ArrayExpression
    """
    pops_to_use = hl.literal(populations)
    freq = hl.map(lambda x: x[0].annotate(meta=x[1]), hl.zip(freq, freq_meta))
    freqs_to_use = hl.filter(lambda f:
                             ((f.meta.size() == 1) & (f.meta.get('group') == 'adj')) |
                             ((f.meta.size() == 2) & (f.meta.get('group') == 'adj') & pops_to_use.contains(f.meta.get('pop'))) |
                             (~locus.in_autosome_or_par() & (
                                     ((f.meta.size() == 2) & (f.meta.get('group') == 'adj') & f.meta.contains('sex')) |
                                     ((f.meta.size() == 3) & (f.meta.get('group') == 'adj') & pops_to_use.contains(f.meta.get('pop')) & f.meta.contains('sex')))),
                             freq)
    return freqs_to_use.map(lambda f: hl.struct(
        meta=f.meta,
        faf95=hl.experimental.filtering_allele_frequency(f.AC, f.AN, 0.95),
        faf99=hl.experimental.filtering_allele_frequency(f.AC, f.AN, 0.99)
    ))


def generate_frequency_data(mt: hl.MatrixTable, calculate_downsampling: bool = False,
                            calculate_by_platform: bool = False) -> Tuple[hl.Table, hl.Table]:
    """
    :param MatrixTable mt: Input MatrixTable
    :param bool calculate_downsampling: Calculate frequencies for downsampled data
    :param bool calculate_by_platform: Calculate frequencies for PCR-free data
    """
    if calculate_downsampling:
        mt, downsamplings = generate_downsamplings_cumulative(mt)
        print(f'Got {len(downsamplings)} downsamplings: {downsamplings}')
    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.meta.pop), hl.agg.counter(mt.meta.pop)),
                'sex': hl.agg.filter(hl.is_defined(mt.meta.sex), hl.agg.collect_as_set(mt.meta.sex)),
                'subpop': hl.agg.filter(hl.is_defined(mt.meta.subpop) & hl.is_defined(mt.meta.pop),
                                        hl.agg.collect_as_set(hl.struct(subpop=mt.meta.subpop, pop=mt.meta.pop)))
                }
    if calculate_by_platform:
        cut_dict['platform'] = hl.agg.filter(hl.is_defined(mt.meta.qc_platform),
                                             hl.agg.collect_as_set(mt.meta.qc_platform))
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))

    sample_group_filters = [({}, True)]
    sample_group_filters.extend([
        ({'pop': pop}, mt.meta.pop == pop) for pop in cut_data.pop
    ] + [
        ({'sex': sex}, mt.meta.sex == sex) for sex in cut_data.sex
    ] + [
        ({'pop': pop, 'sex': sex}, (mt.meta.sex == sex) & (mt.meta.pop == pop))
        for sex in cut_data.sex for pop in cut_data.pop
    ] + [
        ({'subpop': subpop.subpop, 'pop': subpop.pop},
         mt.meta.subpop == subpop.subpop)
        for subpop in cut_data.subpop
    ])

    if calculate_by_platform:
        sample_group_filters.extend([
            ({'platform': str(platform)}, mt.meta.qc_platform == platform)
            for platform in cut_data.platform
        ])

    if calculate_downsampling:
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': 'global'},
             mt.downsampling.global_idx < ds) for ds in downsamplings
        ])
        sample_group_filters.extend([
            ({'downsampling': str(ds), 'pop': pop},
             (mt.downsampling.pop_idx < ds) & (mt.meta.pop == pop))
            for ds in downsamplings for pop, pop_count in cut_data.pop.items() if ds <= pop_count
        ])
    mt = mt.select_cols(group_membership=tuple(x[1] for x in sample_group_filters), project_id=mt.meta.project_id, age=mt.meta.age)
    mt = mt.select_rows()

    frequency_expression = []
    meta_expressions = []
    for i in range(len(sample_group_filters)):
        subgroup_dict = sample_group_filters[i][0]
        subgroup_dict['group'] = 'adj'
        call_stats = hl.agg.filter(mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles))
        call_stats_bind = hl.bind(lambda cs: cs.annotate(
            AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
        ), call_stats)
        frequency_expression.append(call_stats_bind)
        meta_expressions.append(subgroup_dict)

    raw_stats = hl.agg.call_stats(mt.GT, mt.alleles)
    raw_stats_bind = hl.bind(lambda cs: cs.annotate(
        AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
    ), raw_stats)
    frequency_expression.insert(1, raw_stats_bind)
    meta_expressions.insert(1, {'group': 'raw'})

    print(f'Calculating {len(frequency_expression)} aggregators...')
    global_expression = {
        'freq_meta': meta_expressions
    }
    mt = mt.annotate_rows(freq=frequency_expression,
                          age_hist_het=hl.agg.filter(mt.adj & mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
                          age_hist_hom=hl.agg.filter(mt.adj & mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10)))
    if calculate_downsampling: global_expression['downsamplings'] = downsamplings
    mt = mt.annotate_globals(**global_expression)
    sample_data = mt.cols()

    pops = set(cut_data.pop.keys())
    [pops.discard(x) for x in POPS_TO_REMOVE_FOR_POPMAX]

    mt = mt.annotate_rows(
        popmax=pop_max_expr(mt.freq, mt.freq_meta, pops_to_exclude=set(POPS_TO_REMOVE_FOR_POPMAX)),
        faf=add_faf_expr(mt.freq, mt.freq_meta, mt.locus, populations=pops),
        project_max=project_max_expr(mt.project_id, mt.GT, mt.alleles)

    )

    return mt.rows(), sample_data


def generate_consanguineous_frequency_data(ht: hl.Table, f_ht: hl.Table, mt: hl.MatrixTable) -> hl.Table:
    mt = mt.annotate_cols(meta=mt.meta.annotate(f=f_ht[mt.col_key].pop_f))
    projects_to_use = hl.literal({'PROMIS', 'T2D_28K_PROMIS_exomes', 'T2D_28K_Singapore_S.Asian_exomes', 'T2DGENES'})
    consang_call_stats = hl.agg.filter((mt.meta.f.f_stat >= F_CUTOFF) & mt.adj &
                                       (mt.meta.pop == 'sas') & projects_to_use.contains(mt.meta.project_description),
                                       hl.agg.call_stats(mt.GT, mt.alleles))
    consang_call_stats_bind = hl.bind(lambda cs: cs.annotate(
        AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
    ), consang_call_stats)
    mt = mt.annotate_rows(consang=consang_call_stats_bind)
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.append({'group': 'adj', 'pop': 'sas', 'sample_set': 'consanguineous'})
    )
    return ht.annotate(
        freq=ht.freq.append(mt.rows()[ht.key].consang)
    )


def main(args):
    hl.init(log='/frequency_annotations.log')

    data_type = 'genomes' if args.genomes else 'exomes'

    mt = get_gnomad_data(data_type, release_samples=True)
    if args.subset_samples_by_field:
        mt = mt.filter_cols(mt.meta[args.subset_samples_by_field], keep=not args.invert)
    logger.info(f'Calculating frequencies for {mt.count_cols()} samples with {args.subset_samples_by_field} == {not args.invert}')

    location = f'frequencies_{args.subset_samples_by_field}' if args.subset_samples_by_field else 'frequencies'
    if args.calculate_frequencies:
        ht, sample_table = generate_frequency_data(mt, args.downsampling, args.by_platform)

        write_temp_gcs(ht, annotations_ht_path(data_type, location), args.overwrite)
        if args.downsampling:
            sample_table.write(sample_annotations_table_path(data_type, 'downsampling'), args.overwrite)

    if args.add_consanguineous_frequencies and not args.genomes:
        f_ht = hl.read_table('gs://gnomad/sample_qc/f_stat/test_f_0.05_adj_filtered.ht')
        f_ht = f_ht.filter(f_ht.data_type == 'exomes').key_by('s')
        ht = hl.read_table(annotations_ht_path(data_type, location))
        ht = generate_consanguineous_frequency_data(ht, f_ht, mt)
        write_temp_gcs(ht, annotations_ht_path(data_type, f'{location}_with_consanguineous'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--downsampling', help='Also calculate downsampling frequency data', action='store_true')
    parser.add_argument('--calculate_frequencies', help='Calcualte most frequency data', action='store_true')
    parser.add_argument('--add_consanguineous_frequencies', help='Add consanguineous frequencies to existing file', action='store_true')
    parser.add_argument('--by_platform', help='Also calculate frequencies by platform', action='store_true')
    parser.add_argument('--subset_samples_by_field', help='Annotation field containing boolean values describing samples to keep in subset')
    parser.add_argument('--invert', help='Remove samples in --subset_samples_by_field instead of keeping (e.g. neuro)', action='store_true')
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
