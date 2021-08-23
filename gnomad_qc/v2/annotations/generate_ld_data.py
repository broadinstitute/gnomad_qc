import gnomad.resources.grch37.gnomad_ld as ld_resources
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources import *
from hail.utils import new_temp_file
from hail.utils.java import Env
from hail.linalg import BlockMatrix
import sys
import argparse

COMMON_FREQ = 0.005
RARE_FREQ = 0.0005
SNV_SV_POPS = ('nfe', 'afr')


def ld_pruned_path(data_type: str, pop: str, r2: str, version: str = CURRENT_RELEASE):
    return f'gs://gnomad-resources/ld/gnomad.{data_type}.r{version}.{pop}.ld.pruned_set.r2_{r2}.ht'


def get_pop_and_subpop_counters(mt):
    cut_dict = {'pop': hl.agg.filter(hl.is_defined(mt.meta.pop) & (mt.meta.pop != 'oth'), hl.agg.counter(mt.meta.pop)),
                'subpop': hl.agg.filter(hl.is_defined(mt.meta.subpop) & (mt.meta.subpop != 'oea') &
                                        (mt.meta.subpop != 'onf'), hl.agg.counter(mt.meta.subpop))
                }
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    return cut_data


def filter_mt_for_ld(mt, label, pop, common_only: bool = True, sv_dataset: bool = False):
    frequency = COMMON_FREQ if common_only else RARE_FREQ
    # At the moment, ld_matrix OOMs above ~30M variants, so filtering all populations to 0.05% to cut down on variants
    # All work on standard machines except AFR
    # Also creating `common_only` ones for most practical uses (above 0.5%)
    pop_mt = mt.filter_cols(mt.meta[label] == pop)
    meta_index = hl.zip_with_index(pop_mt.freq_meta).find(lambda f: (f[1].get(label) == pop)).collect()[0][0]
    if sv_dataset:
        call_stats = hl.agg.call_stats(pop_mt.GT, pop_mt.alleles)
        call_stats_bind = hl.bind(lambda cs: cs.annotate(
            AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
        ), call_stats)
        pop_freq = call_stats_bind
    else:
        pop_mt = pop_mt.filter_rows((hl.len(pop_mt.filters) == 0))
        pop_freq = pop_mt.freq[meta_index]
    pop_mt = pop_mt.annotate_rows(pop_freq=pop_freq)
    pop_mt = pop_mt.filter_rows((pop_mt.pop_freq.AC > 1) &
                                (pop_mt.pop_freq.AN - pop_mt.pop_freq.AC > 1) &
                                (pop_mt.pop_freq.AF > frequency) &
                                #  the following filter is needed for "het-only" monomorphic sites
                                #  as otherwise variance is 0, and so, row_correlation errors out
                                ~((pop_mt.pop_freq.AF == 0.5) &
                                  (pop_mt.pop_freq.homozygote_count == 0)))
    if sv_dataset:
        pop_mt = pop_mt.filter_rows(pop_mt.pop_freq.AN >= 2000)
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


def generate_ld_matrix(mt, pop_data, data_type, radius: int = 1000000, common_only: bool = True,
                       adj: bool = False, overwrite: bool = False):
    # Takes about 4 hours on 20 n1-standard-8 nodes (with SSD - not sure if necessary) per population
    # Total of ~37 hours ($400)
    for label, pops in dict(pop_data).items():
        for pop in pops:
            if data_type == 'genomes_snv_sv' and pop not in SNV_SV_POPS: continue
            pop_mt = filter_mt_for_ld(mt, label, pop, common_only, sv_dataset=data_type == 'genomes_snv_sv')

            pop_mt.rows().select('pop_freq').add_index().write(ld_resources._ld_index_path(data_type, pop, common_only, adj), overwrite)
            ld = hl.ld_matrix(pop_mt.GT.n_alt_alleles(), pop_mt.locus, radius)
            if data_type != 'genomes_snv_sv':
                ld = ld.sparsify_triangle()
            ld.write(ld_resources._ld_matrix_path(data_type, pop, common_only, adj), overwrite)


def export_snv_sv_ld_matrix(pop_data, data_type, common_only: bool = True, adj: bool = False, overwrite: bool = False):
    for label, pops in dict(pop_data).items():
        for pop in pops:
            if pop not in SNV_SV_POPS: continue
            bm = BlockMatrix.read(ld_resources._ld_matrix_path(data_type, pop, common_only, adj))
            ld_index = hl.read_table(ld_resources._ld_index_path(data_type, pop, common_only, adj))
            snvs = ld_index.filter(ld_index.alleles[0] != "N")
            svs = ld_index.filter(ld_index.alleles[0] == "N")
            snv_indices = snvs.idx.collect()
            sv_indices = svs.idx.collect()
            ht = bm.filter(snv_indices, sv_indices).entries(keyed=False)
            ht.filter(ht.entry != 0).write(ld_resources._ld_snv_sv_path(pop), overwrite)

            hl.read_table(ld_resources._ld_snv_sv_path(pop)).export(ld_resources._ld_snv_sv_path(pop).replace('.ht', '.txt.bgz'))
            snvs = snvs.add_index().key_by()
            svs = svs.add_index().key_by()
            snvs.select(chrom=snvs.locus.contig, pos=snvs.locus.position, ref=snvs.alleles[0], alt=snvs.alleles[1],
                        i=snvs.idx).export(ld_resources._ld_snv_sv_index_path(pop, 'snv'))
            svs.select(chrom=svs.locus.contig, pos=svs.locus.position, ref=svs.alleles[0], alt=svs.alleles[1],
                       j=svs.idx).export(ld_resources._ld_snv_sv_index_path(pop, 'sv'))


def generate_ld_scores_from_ld_matrix(pop_data, data_type, min_frequency=0.01, call_rate_cutoff=0.8,
                                      adj: bool = False, radius: int = 1000000,
                                      overwrite=False):
    # This function required a decent number of high-mem machines (with an SSD for good measure) to complete the AFR
    # For the rest, on 20 n1-standard-8's, 1h15m to export block matrix, 15 mins to compute LD scores per population (~$150 total)
    for label, pops in dict(pop_data).items():
        for pop, n in pops.items():
            ht = hl.read_table(ld_resources._ld_index_path(data_type, pop, adj=adj))
            ht = ht.filter((ht.pop_freq.AF >= min_frequency) &
                           (ht.pop_freq.AF <= 1 - min_frequency) &
                           (ht.pop_freq.AN / n >= 2 * call_rate_cutoff)).add_index(name='new_idx')

            indices = ht.idx.collect()

            r2 = BlockMatrix.read(ld_resources._ld_matrix_path(data_type, pop, min_frequency >= COMMON_FREQ, adj=adj))
            r2 = r2.filter(indices, indices) ** 2
            r2_adj = ((n - 1.0) / (n - 2.0)) * r2 - (1.0 / (n - 2.0))

            out_name = ld_resources._ld_scores_path(data_type, pop, adj)
            compute_and_annotate_ld_score(ht, r2_adj, radius, out_name, overwrite)


def generate_all_cross_pop_ld_scores(pop_data, data_type, min_frequency=0.01, call_rate_cutoff=0.8,
                                     adj: bool = False, radius: int = 1000000,
                                     overwrite=False):
    for label, pops in dict(pop_data).items():
        for pop1, n in pops.items():
            for label, pops in dict(pop_data).items():
                for pop2, n in pops.items():
                    if (pop1, pop2) in (
                            ('eas', 'nfe'), ('afr', 'nfe'), ('eas', 'afr'),
                            ('fin', 'nfe'), ('fin', 'eas'),
                    ):
                        generate_cross_pop_ld_scores_from_ld_matrices(pop1, pop2, data_type, pop_data, min_frequency,
                                                                      call_rate_cutoff, adj, radius, overwrite)


def generate_cross_pop_ld_scores_from_ld_matrices(pop1, pop2, data_type, pop_data, min_frequency=0.01, call_rate_cutoff=0.8,
                                                  adj: bool = False, radius: int = 1000000, overwrite=False,
                                                  temp_bucket='gs://gnomad-tmp/ld'):
    n1 = pop_data.pop[pop1]
    n2 = pop_data.pop[pop2]
    ht1 = hl.read_table(ld_resources._ld_index_path(data_type, pop1, adj=adj))
    ht1 = ht1.filter((ht1.pop_freq.AF >= min_frequency) &
                     (ht1.pop_freq.AF <= 1 - min_frequency) &
                     (ht1.pop_freq.AN / n1 >= 2 * call_rate_cutoff))

    ht2 = hl.read_table(ld_resources._ld_index_path(data_type, pop2, adj=adj))
    ht2 = ht2.filter((ht2.pop_freq.AF >= min_frequency) &
                     (ht2.pop_freq.AF <= 1 - min_frequency) &
                     (ht2.pop_freq.AN / n2 >= 2 * call_rate_cutoff))

    ht1 = ht1.filter(hl.is_defined(ht2[ht1.key])).add_index(name='new_idx').checkpoint(f'{temp_bucket}/{pop1}_{pop2}.ht', overwrite=overwrite, _read_if_exists=not overwrite)
    ht2 = ht2.filter(hl.is_defined(ht1[ht2.key])).add_index(name='new_idx').checkpoint(f'{temp_bucket}/{pop2}_{pop1}.ht', overwrite=overwrite, _read_if_exists=not overwrite)
    indices1 = ht1.idx.collect()
    indices2 = ht2.idx.collect()
    assert len(indices1) == len(indices2)

    r1 = BlockMatrix.read(ld_resources._ld_matrix_path(data_type, pop1, min_frequency >= COMMON_FREQ, adj=adj)).filter(indices1, indices1)
    r2 = BlockMatrix.read(ld_resources._ld_matrix_path(data_type, pop2, min_frequency >= COMMON_FREQ, adj=adj)).filter(indices2, indices2)
    r_bm = r1 * r2

    # TODO: is a bias adjustment needed?
    # r2_adj = ((n - 1.0) / (n - 2.0)) * r2 - (1.0 / (n - 2.0))

    out_name = ld_resources._cross_pop_ld_scores_path(data_type, pop1, pop2, adj)
    compute_and_annotate_ld_score(ht1, r_bm, radius, out_name, overwrite)


def compute_and_annotate_ld_score(ht, r2_adj, radius, out_name, overwrite):
    starts_and_stops = hl.linalg.utils.locus_windows(ht.locus, radius, _localize=False)

    # Lifted directly from https://github.com/hail-is/hail/blob/555e02d6c792263db2c3ed97db8002b489e2dacb/hail/python/hail/methods/statgen.py#L2595
    # for the time being, until efficient BlockMatrix filtering gets an easier interface
    # This is required, as the squaring/multiplication densifies, so this re-sparsifies.
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
    ht = ht.annotate(**ht_scores[ht.new_idx]).select_globals()
    ht.filter(hl.is_defined(ht.ld_score)).write(out_name, overwrite)


def main(args):
    hl.init(log='/ld_annotations.log')

    data_type = 'genomes' if args.genomes else 'exomes'
    # Only run on genomes so far

    mt = get_gnomad_data(data_type, release_samples=True, release_annotations=True, adj=args.adj)
    pop_data = get_pop_and_subpop_counters(mt)

    if args.snv_sv_comparison:
        mt = hl.read_matrix_table('gs://gnomad/projects/genomes_sv_ld/genomes_union_sv.mt')
        mt = mt.annotate_cols(meta=get_gnomad_meta('genomes')[mt.genomes_s])
        if args.adj:
            mt = mt.filter_entries(hl.is_missing(mt.GQ) | (mt.GQ >= 20))
        mt = mt.annotate_rows(**get_gnomad_public_data('genomes')[mt.row_key]).annotate_globals(
            **get_gnomad_public_data('genomes').index_globals())
        pop_data = get_pop_and_subpop_counters(mt)
        generate_ld_matrix(mt, pop_data, 'genomes_snv_sv', args.radius, args.common_only, args.adj, args.overwrite)
        export_snv_sv_ld_matrix(pop_data, 'genomes_snv_sv', args.common_only, args.adj, args.overwrite)


    if args.generate_ld_pruned_set:
        generate_ld_pruned_set(mt, pop_data, data_type, args.r2, args.radius, args.overwrite)

    if args.generate_ld_matrix:
        generate_ld_matrix(mt, pop_data, data_type, args.radius, args.common_only, args.adj, args.overwrite)

    if args.generate_ld_scores:
        generate_ld_scores_from_ld_matrix(pop_data, data_type, args.min_frequency, args.min_call_rate, args.adj, overwrite=args.overwrite)

    if args.cross_pop_ld_scores:
        generate_all_cross_pop_ld_scores(pop_data, data_type, args.min_frequency, args.min_call_rate, args.adj, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--snv_sv_comparison', help='Compute LD between SNVs and SVs.', action='store_true')
    parser.add_argument('--generate_ld_pruned_set', help='Calculates LD pruned set of variants', action='store_true')
    parser.add_argument('--generate_ld_matrix', help='Calculates LD matrix', action='store_true')
    parser.add_argument('--common_only', help='Calculates LD matrix only on common variants (above 0.5%)', action='store_true')
    parser.add_argument('--adj', help='Calculates LD matrix only on using high-quality genotypes', action='store_true')
    parser.add_argument('--generate_ld_scores', help='Calculates LD scores from LD matrix', action='store_true')
    parser.add_argument('--cross_pop_ld_scores', help='Calculates cross-pop LD scores from LD matrix', action='store_true')
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
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
