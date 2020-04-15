import pickle
from gnomad.utils.annotations import get_adj_expr
from gnomad_qc.v3.resources import get_gnomad_v3_mt
import argparse
import hail as hl


def print_ref_block_stats(path: str):
    import numpy as np
    def _print_block_stats(stats: hl.Struct):
        def get_quantile(cum_prop, quantile):
            return [i for i, x in enumerate(cum_prop) if x > quantile][0]

        for strat, ref_block_stats in [('all', stats.ref_block_stats), ('adj', stats.adj_ref_block_stats)]:
            n_blocks = np.sum(ref_block_stats.hist.bin_freq) + ref_block_stats.hist.n_larger
            cum_blocks = np.cumsum(ref_block_stats.hist.bin_freq)
            cum_prop = [x / n_blocks for x in cum_blocks]

            print(f"Stats for {strat}")
            print(f"Number of blocks: {n_blocks}")
            print(f"Largest block size: {ref_block_stats.stats.max}")
            print(f"95% block size: {get_quantile(cum_prop, 0.95)}")
            print(f"99% block size: {get_quantile(cum_prop, 0.99)}")
            print(f"99.9% block size: {get_quantile(cum_prop, 0.999)}")
            print(f"99.95% block size: {get_quantile(cum_prop, 0.9995)}")
            print(f"Percentage blocks below 10k: {1-(ref_block_stats.hist.n_larger/n_blocks)}")

    if path.startswith('gs://'):
        with hl.hadoop_open(path, 'rb') as f:
            _print_block_stats(pickle.load(f))
    else:
        with open(path, 'rb') as f:
            _print_block_stats(pickle.load(f))

def compute_stats(stats_path: str):
    mt = get_gnomad_v3_mt()
    mt = mt.filter_entries(hl.is_defined(mt.END))
    ref_block_stats = mt.aggregate_entries(
        hl.struct(
            ref_block_stats=hl.struct(
                stats=hl.agg.stats(mt.END - mt.locus.position),
                hist=hl.agg.hist(mt.END - mt.locus.position, 0, 9999, 10000),
                hist_log=hl.agg.hist(hl.log10(1 + mt.END - mt.locus.position), 0, 5, 100)
            ),
            adj_ref_block_stats=hl.agg.filter(
                get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
                hl.struct(
                    stats=hl.agg.stats(mt.END - mt.locus.position),
                    hist=hl.agg.hist(mt.END - mt.locus.position, 0, 9999, 10000),
                    hist_log=hl.agg.hist(hl.log10(1 + mt.END - mt.locus.position), 0, 5, 100)
                )
            )
        )
    )

    with hl.hadoop_open(stats_path, 'wb') as f:
        pickle.dump(ref_block_stats, f)


def main(args):
    if args.compute_stats:
        compute_stats(args.stats_path)
        print_ref_block_stats('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3.repartitioned.ref_block_stats.pickle')

    elif args.print_stats:
        print_ref_block_stats(args.stats_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--compute_stats', help='Lift-over Purcell5k intervals', action='store_true')
    parser.add_argument('--stats_path', help='Path to stats file', default='gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3.repartitioned.ref_block_stats.pickle')
    parser.add_argument('--print_stats', help='Print existing stats', action='store_true')
    main(parser.parse_args())
