from gnomad.utils.annotations import annotate_adj
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import get_gnomad_data, get_gnomad_data_path
import sys
import argparse
import hail as hl


def main(args):
    hl.init()

    data_type = 'genomes' if args.genomes else 'exomes'

    if args.write_hardcalls:
        mt = get_gnomad_data(data_type, split=False, raw=True, meta_root=None)
        ht = hl.read_table(qc_ht_path(data_type, 'hard_filters'))
        mt = annotate_adj(mt.select_cols(sex=ht[hl.literal(data_type), mt.s].sex))
        mt = mt.select_entries(GT=hl.case(missing_false=True).when(hl.call(mt.PGT[0], mt.PGT[1]) == mt.GT, mt.PGT).default(mt.GT),
                               PID=mt.PID, adj=mt.adj)
        mt = adjust_sex_ploidy(mt, mt.sex)
        mt = mt.select_cols().naive_coalesce(10000)
        mt.write(get_gnomad_data_path(data_type, hardcalls=True, split=False), args.overwrite)

    if args.split_hardcalls:
        mt = get_gnomad_data(data_type, split=False, meta_root=None)
        mt = hl.split_multi_hts(mt)
        mt.write(get_gnomad_data_path(data_type, hardcalls=True, split=True), args.overwrite)

    if args.write_nonrefs:  # CPU-hours: 600 (E)
        mt = get_gnomad_data(data_type, split=False, raw=True, meta_root=None).select_cols()
        mt = mt.annotate_entries(is_missing=hl.is_missing(mt.GT))
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt = annotate_adj(mt)
        if args.exomes:
            mt = mt.naive_coalesce(10000)
        mt.write(get_gnomad_data_path(data_type, split=False, non_refs_only=True), args.overwrite)

    if args.split_nonrefs:  # CPU-hours: 300 (E)
        mt = get_gnomad_data(data_type, split=False, non_refs_only=True)
        mt = hl.split_multi_hts(mt)
        mt = mt.filter_entries(mt.is_missing | mt.GT.is_non_ref())
        mt.write(get_gnomad_data_path(data_type, split=True, non_refs_only=True), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--write_hardcalls', help='Creates a hardcalls mt', action='store_true')
    parser.add_argument('--split_hardcalls', help='Creates a split hardcalls mt from the hardcalls mt', action='store_true')
    parser.add_argument('--write_nonrefs', help='Creates a sparse mt with only non-ref genotypes', action='store_true')
    parser.add_argument('--split_nonrefs', help='Creates a split sparse mt with only non-ref genotypes', action='store_true')
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
