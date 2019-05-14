from gnomad_hail import *
from hail.utils import new_temp_file
from hail.utils.java import Env
from hail.linalg import BlockMatrix


def ld_matrix_path(data_type: str, pop: str, version: str = CURRENT_RELEASE):
    return f'gs://gnomad-resources/ld/gnomad.{data_type}.r{version}.{pop}.ld.bm'


def ld_index_path(data_type: str, pop: str, version: str = CURRENT_RELEASE):
    return f'gs://gnomad-resources/ld/gnomad.{data_type}.r{version}.{pop}.ld.variant_indices.ht'


def ld_pruned_path(data_type: str, pop: str, r2: str, version: str = CURRENT_RELEASE):
    return f'gs://gnomad-resources/ld/gnomad.{data_type}.r{version}.{pop}.ld.pruned_set.r2_{r2}.ht'


def ld_scores_path(data_type: str, pop: str, version: str = CURRENT_RELEASE):
    return f'gs://gnomad-resources/ld/gnomad.{data_type}.r{version}.{pop}.ld_scores.ht'


def get_pop_and_subpop_counters(mt):
    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.meta.pop) & (mt.meta.pop != 'oth'), hl.agg.counter(mt.meta.pop)),
                'subpop': hl.agg.filter(hl.is_defined(mt.meta.subpop) & (mt.meta.subpop != 'oea') &
                                        (mt.meta.subpopt != 'onf'), hl.agg.counter(mt.meta.subpop))
                }
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    return cut_data


def filter_mt_for_ld(mt, label, pop):
    pop_mt = mt.filter_cols(mt.meta[label] == pop)
    meta_index = hl.zip_with_index(pop_mt.freq_meta).find(lambda f: (f[1].get(label) == pop)).collect()[0][0]
    pop_mt = pop_mt.annotate_rows(pop_freq=pop_mt.freq[meta_index])
    pop_mt = pop_mt.filter_rows((hl.len(pop_mt.filters) == 0) &
                                (pop_mt.pop_freq.AC > 1) &
                                (pop_mt.pop_freq.AN - pop_mt.pop_freq.AC > 1) &
                                #  the following filter is needed for "het-only" monomorphic sites
                                #  as otherwise variance is 0, and so, row_correlation errors out
                                ~((pop_mt.pop_freq.AF == 0.5) &
                                  (pop_mt.pop_freq.homozygote_count == 0)))
    return pop_mt


def generate_ld_pruned_set(mt: hl.MatrixTable, pop_data: dict, data_type: str, r2: str = '0.2',
                           radius: int = 1000000, overwrite: bool = False):
    for label, pops in dict(pop_data).items():
        for pop in pops:
            pop_mt = filter_mt_for_ld(mt, label, pop)

            pruned_ht = hl.ld_prune(pop_mt.GT, r2=float(r2), bp_window_size=radius)

            ht = pop_mt.rows().select('pop_freq')
            ht = ht.filter(hl.is_defined(pruned_ht[ht.key]))
            ht.write(ld_pruned_path(data_type, pop, r2), overwrite)


def generate_ld_matrix(mt, pop_data, data_type, radius: int = 1000000, overwrite: bool = False):
    for label, pops in dict(pop_data).items():
        for pop in pops:
            pop_mt = filter_mt_for_ld(mt, label, pop)

            pop_mt.rows().select('pop_freq').add_index().write(ld_index_path(data_type, pop), overwrite)
            ld = hl.ld_matrix(pop_mt.GT.n_alt_alleles(), pop_mt.locus, radius).sparsify_triangle()
            ld.write(ld_matrix_path(data_type, pop), overwrite)


def generate_ld_scores_from_ld_matrix(pop_data, data_type, min_frequency=0.01, call_rate_cutoff=0.8,
                                      radius: int = 1000000,
                                      overwrite=False):
    for label, pops in dict(pop_data).items():
        for pop, n in pops.items():
            ht = hl.read_table(ld_index_path(data_type, pop))
            ht = ht.filter((ht.pop_freq.AF >= min_frequency) &
                           (ht.pop_freq.AF <= 1 - min_frequency) &
                           (ht.pop_freq.AN / n >= 2 * call_rate_cutoff)).add_index(name='new_idx')

            indices = ht.idx.collect()

            r2 = BlockMatrix.read(ld_matrix_path(data_type, pop))
            r2 = r2.filter(indices, indices) ** 2
            r2_adj = ((n - 1.0) / (n - 2.0)) * r2 - (1.0 / (n - 2.0))

            starts_and_stops = hl.linalg.utils.locus_windows(ht.locus, radius, _localize=False)

            # Lifted directly from https://github.com/hail-is/hail/blob/555e02d6c792263db2c3ed97db8002b489e2dacb/hail/python/hail/methods/statgen.py#L2595
            # for the time being, until efficient BlockMatrix filtering gets an easier interface
            r2_adj = BlockMatrix._from_java(r2_adj._jbm.filterRowIntervalsIR(
                Env.backend()._to_java_ir(starts_and_stops._ir),
                False))

            l2row = r2_adj.sum(axis=0).T
            l2col = r2_adj.sum(axis=1)
            l2 = l2row + l2col + 1

            l2_bm_tmp = new_temp_file()
            l2_tsv_tmp = new_temp_file()
            l2.write(l2_bm_tmp, force_row_major=True)
            BlockMatrix.export(l2_bm_tmp, l2_tsv_tmp)

            ht_scores = hl.import_table(l2_tsv_tmp, no_header=True, impute=True)
            ht_scores = ht_scores.add_index().rename({'f0': 'ld_score'})
            ht_scores = ht_scores.key_by('idx')

            ht = ht.annotate(**ht_scores[ht.new_idx])
            ht.filter(hl.is_defined(ht.ld_score)).write(ld_scores_path(data_type, pop), overwrite)


def main(args):
    hl.init(log='/ld_annotations.log')

    data_type = 'genomes' if args.genomes else 'exomes'

    mt = get_gnomad_data(data_type, release_samples=True, release_annotations=True)
    pop_data = get_pop_and_subpop_counters(mt)

    if args.generate_ld_pruned_set:
        generate_ld_pruned_set(mt, pop_data, data_type, args.r2, args.radius, args.overwrite)

    if args.generate_ld_matrix:
        generate_ld_matrix(mt, pop_data, data_type, args.radius, args.overwrite)

    if args.generate_ld_scores:
        generate_ld_scores_from_ld_matrix(pop_data, data_type, args.min_frequency, args.min_call_rate, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--generate_ld_pruned_set', help='Calculates LD pruned set of variants', action='store_true')
    parser.add_argument('--generate_ld_matrix', help='Calculates LD matrix', action='store_true')
    parser.add_argument('--generate_ld_scores', help='Calculates LD scores from LD matrix', action='store_true')
    parser.add_argument('--min_frequency', help='Minimum allele frequency to compute LD scores (default 0.01)', default=0.01, type=float)
    parser.add_argument('--min_call_rate', help='Minimum call rate to compute LD scores (default 0.8)', default=0.8, type=float)
    parser.add_argument('--r2', help='r-squared to which to prune LD (default 0.2)', default="0.2")
    parser.add_argument('--radius', help='Radius at which to calculate LD information (bp; default 1e6)', default=1000000, type=int)
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
