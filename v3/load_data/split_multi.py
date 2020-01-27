import hail as hl
import argparse
from v3.resources import get_full_mt, get_full_mt_path

def main(args):
    mt = get_full_mt(False)
    mt = hl.MatrixTable(
        hl.ir.MatrixKeyRowsBy(mt._mir, ['locus', 'alleles'], is_sorted=True))
    mt = mt.annotate_entries(
        gvcf_info=mt.gvcf_info.drop('ClippingRankSum', 'ReadPosRankSum'))
    mt = mt.annotate_rows(
        n_unsplit_alleles=hl.len(mt.alleles),
        mixed_site=(hl.len(mt.alleles) > 2) & hl.any(
            lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:]) & hl.any(
            lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]))
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt.write(get_full_mt_path(), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrites existing files', action='store_true')
    main(parser.parse_args())
