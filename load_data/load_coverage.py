
from gnomad_hail import *


def main(args):
    data_type = 'exomes' if args.exomes else 'genomes'
    num_partitions = 1000  # if args.exomes else 10000
    hl.init(min_block_size=0, log='/load_coverage.log')

    source_root = f'gs://gnomad/coverage/source/{data_type}'
    root = f'gs://gnomad/coverage/hail-0.2/{data_type}'

    all_file_data = []
    with hl.hadoop_open(f'{source_root}/coverage_files.txt', 'r') as f:
        for file_data in f:
            fname, anchored, sample_data = file_data.strip().split('\t', 2)
            base = fname.split('.txt')[0]
            anchored = anchored == 'anchored'
            sample_data = sample_data.split('\t')
            all_file_data.append([base, anchored, sample_data])

    sample_count = 9733 if args.exomes else 1279
    assert sum([len(x[2]) for x in all_file_data]) == sample_count

    meta_kt = get_gnomad_meta(data_type, full_meta=True)
    bam_dict = dict(get_sample_data(meta_kt, [meta_kt.bam, meta_kt.s]))

    assert all([all([y in bam_dict for y in x[2]]) for x in all_file_data])

    # Modify all_file_data in place
    for file_data in all_file_data:
        for i in range(len(file_data[2])):
            file_data[2][i] = bam_dict[file_data[2][i]]

    if args.read_coverage_files:
        # 20-30 seconds each for exomes, 2m45s for genomes
        for base, anchored, sample_data in all_file_data:
            fname = f'{source_root}/parts/full_{base}.gz'
            print('Loading:', fname)
            if anchored:
                mt = hl.import_matrix_table(fname, no_header=True, row_fields={'f0': hl.tstr, 'f1': hl.tint}, min_partitions=num_partitions, force_bgz=True)
                mt = mt.transmute_rows(locus=hl.locus(mt.f0, mt.f1))
            else:
                mt = hl.import_matrix_table(fname, no_header=True, min_partitions=num_partitions, force_bgz=True)
            mt = mt.key_cols_by(s=hl.literal(sample_data)[mt.col_id]).drop('col_id')
            mt.transmute_entries(coverage=mt.x).write(f'{root}/parts/part_{base}.mt', args.overwrite)

    if args.merge_coverage_mts:
        # Exomes: first merges are ~7.5 mins each, final one is ~1.5 hours (on 40 n1-standard-8s)
        chunks = int(len(all_file_data) ** 0.5) + 1
        for i in range(chunks):
            if i * chunks >= len(all_file_data): break
            base, anchored, sample_data = all_file_data[i * chunks]
            mt = hl.read_matrix_table(f'{root}/parts/part_{base}.mt')
            if i:
                mt = mt.select_rows()
            else:
                assert anchored
            for j in range(1, chunks):
                if i * chunks + j >= len(all_file_data): break
                base, anchored, sample_data = all_file_data[i * chunks + j]
                next_mt = hl.read_matrix_table(f'{root}/parts/part_{base}.mt').select_rows()
                mt = mt.union_cols(next_mt)
            mt.write(f'{root}/intermediates/intermediate_{i}.mt', args.overwrite)

        mt = hl.read_matrix_table(f'{root}/intermediates/intermediate_0.mt')
        for i in range(1, chunks):
            try:
                next_mt = hl.read_matrix_table(f'{root}/intermediates/intermediate_{i}.mt')
                mt = mt.union_cols(next_mt)
            except Exception:
                pass
        # This part has some trouble for genomes
        if args.exomes:
            mt.write(f'{root}/intermediates/final_intermediate.mt', args.overwrite)
            mt = hl.read_matrix_table(f'{root}/intermediates/final_intermediate.mt')
            mt = mt.key_rows_by('locus')
        mt.write(coverage_mt_path(data_type), args.overwrite)

    if args.aggregate_coverage:
        mt = hl.read_matrix_table(coverage_mt_path(data_type))
        meta_ht = get_gnomad_meta(data_type)
        mt = mt.filter_cols(meta_ht[mt.s].release)
        mt = mt.annotate_rows(mean=hl.agg.mean(mt.coverage),
                              median=hl.median(hl.agg.collect(mt.coverage)),
                              over_1=hl.agg.fraction(mt.coverage >= 1),
                              over_5=hl.agg.fraction(mt.coverage >= 5),
                              over_10=hl.agg.fraction(mt.coverage >= 10),
                              over_15=hl.agg.fraction(mt.coverage >= 15),
                              over_20=hl.agg.fraction(mt.coverage >= 20),
                              over_25=hl.agg.fraction(mt.coverage >= 25),
                              over_30=hl.agg.fraction(mt.coverage >= 30),
                              over_50=hl.agg.fraction(mt.coverage >= 50),
                              over_100=hl.agg.fraction(mt.coverage >= 100))
        ht = mt.rows()
        if args.exomes:
            ht.write(coverage_ht_path(data_type), args.overwrite)
        else:
            ht.write(f'{root}/intermediates/final_aggregated.ht', args.overwrite)
            ht = hl.read_table(f'{root}/intermediates/final_aggregated.ht')
            ht.key_by('locus').write(coverage_ht_path(data_type), args.overwrite)

    if args.aggregate_coverage_pops:
        mt = hl.read_matrix_table(coverage_mt_path(data_type))
        meta_ht = get_gnomad_meta(data_type)
        mt = mt.annotate_cols(meta=meta_ht[mt.s])
        mt = mt.filter_cols(mt.meta.release)
        agg_ds = (mt
                  .group_cols_by(mt.meta.pop, mt.meta.qc_platform)
                  .aggregate(mean=hl.agg.mean(mt.coverage),
                             median=hl.median(hl.agg.collect(mt.coverage)),
                             over_1=hl.agg.fraction(mt.coverage >= 1),
                             over_5=hl.agg.fraction(mt.coverage >= 5),
                             over_10=hl.agg.fraction(mt.coverage >= 10),
                             over_15=hl.agg.fraction(mt.coverage >= 15),
                             over_20=hl.agg.fraction(mt.coverage >= 20),
                             over_25=hl.agg.fraction(mt.coverage >= 25),
                             over_30=hl.agg.fraction(mt.coverage >= 30),
                             over_50=hl.agg.fraction(mt.coverage >= 50),
                             over_100=hl.agg.fraction(mt.coverage >= 100)))
        agg_ds.write(coverage_mt_path(data_type, by_population=True, by_platform=True), args.overwrite)

    if args.export_coverage:
        ht = hl.read_table(coverage_ht_path(data_type)).key_by()
        ht = ht.transmute(chrom=ht.locus.contig, pos=ht.locus.position).select(
            'chrom', 'pos', *list(ht.drop('locus', 'row_id').row))
        ht.export(coverage_ht_path(data_type).replace('.ht', '.tsv.bgz'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--read_coverage_files', help='Read raw coverage .gz files', action='store_true')
    parser.add_argument('--merge_coverage_mts', help='Merge individual coverage .mt files', action='store_true')
    parser.add_argument('--aggregate_coverage', help='Aggregate coverage data', action='store_true')
    parser.add_argument('--aggregate_coverage_pops', help='Aggregate coverage data', action='store_true')
    parser.add_argument('--export_coverage', help='Export coverage data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
