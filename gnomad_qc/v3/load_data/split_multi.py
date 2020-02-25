import hail as hl
import argparse
from gnomad_qc.v3.resources import get_gnomad_v3_mt

def main(args):
    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
    mt = mt.annotate_entries(
        gvcf_info=mt.gvcf_info.drop('ClippingRankSum', 'ReadPosRankSum')
    )
    mt = mt.annotate_rows(
        n_unsplit_alleles=hl.len(mt.alleles),
        mixed_site=(hl.len(mt.alleles) > 2) & hl.any(
            lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:]) & hl.any(
            lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]))
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt.write(args.split_mt_location, overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--create_split_for_real', help='Needed to create the split MT. IMPORTANT: The cost of storing the split MT is likely larger than recomputing it!', action='store_true')
    parser.add_argument('--split_mt_location', help='Specifies where to write the split MT. This should likely be a temporary location as storing the split MT is not a cost-effective long-term solution.', action='store_true')
    parser.add_argument('--overwrite', help='Overwrites existing files', action='store_true')
    args = parser.parse_args()
    if args.create_split_for_real:
        main(args)
    else:
        print("ERROR: If you really want to create the split MT, please add --create_split_for_real. IMPORTANT: The cost of storing the split MT is likely larger than recomputing it and this is therefore not recommended!")
