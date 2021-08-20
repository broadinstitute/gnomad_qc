from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources import *
import sys
import argparse


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

    meta_ht = get_gnomad_meta(data_type, full_meta=True)
    bam_dict = dict(hl.tuple([meta_ht.bam, meta_ht.s]).collect())
    # bam_dict = dict(get_sample_data(meta_ht, [meta_ht.bam, meta_ht.s]))

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
        ht = mt.annotate_rows(
            mean=hl.agg.mean(mt.coverage),
            median=hl.agg.approx_quantiles(mt.coverage, 0.5),
            count_array=hl.rbind(
                hl.agg.counter(hl.min(100, mt.coverage)),
                lambda c: hl.range(0, 100).map(lambda i: c.get(i, 0))
            )
        ).rows()

        ht = ht.transmute(
            **{f'over_{x}': hl.sum(mt.count_array[x:]) for x in [1, 5, 10, 15, 20, 25, 30, 50, 100]}
        )

        if args.exomes:
            ht.write(coverage_ht_path(data_type), args.overwrite)
        else:
            ht = ht.checkpoint('gs://gnomad-tmp/coverage/genomes_agg.ht', overwrite=args.overwrite)
            ht.key_by('locus').write(coverage_ht_path(data_type), args.overwrite)

    if args.aggregate_coverage_platforms:
        meta_ht = get_gnomad_meta(data_type)
        if data_type == 'genomes':
            meta_ht = meta_ht.annotate(qc_platform=hl.cond(meta_ht.pcr_free, 'pcr_free', 'pcr+'))

        mt = hl.read_matrix_table(coverage_mt_path(data_type))
        mt = mt.annotate_cols(meta=meta_ht.select('release', 'qc_platform', 'sex')[mt.s])
        mt = mt.filter_cols(mt.meta.release)

        n = mt.aggregate_cols(hl.agg.counter((mt.meta.qc_platform, mt.meta.sex)), _localize=False)

        grp_mt = (
            mt
                .group_cols_by(mt.meta.qc_platform, mt.meta.sex)
                .aggregate(
                mean=hl.agg.mean(mt.coverage),
                median=hl.agg.approx_quantiles(mt.coverage, 0.5),
                count_array=hl.rbind(
                    hl.agg.counter(hl.min(100, mt.coverage)),
                    lambda c: hl.range(0, 100).map(lambda i: c.get(i, 0))
                )
            )
        )

        grp_mt = grp_mt.transmute_entries(
            **{f'over_{x}': hl.sum(grp_mt.count_array[x:]) for x in [1, 5, 10, 15, 20, 25, 30, 50, 100]}
        )

        grp_mt = grp_mt.annotate_cols(
            n=n[(grp_mt.qc_platform, grp_mt.sex)]
        )
        if args.exomes:
            grp_mt.write(coverage_mt_path(data_type, grouped=True), args.overwrite)
        else:
            grp_mt.checkpoint('gs://gnomad-tmp/coverage/genomes_agg_grouped.mt', overwrite=args.overwrite)
            grp_mt.key_rows_by('locus').write(coverage_mt_path(data_type, grouped=True), args.overwrite)

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
    parser.add_argument('--aggregate_coverage_platforms', help='Aggregate coverage data by platform and sex', action='store_true')
    parser.add_argument('--export_coverage', help='Export coverage data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
