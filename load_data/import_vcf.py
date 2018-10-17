from gnomad_hail import *
import argparse


def main(args):
    data_type = 'exomes' if args.exomes else 'genomes'
    hl.init(min_block_size=args.min_block_size)

    mt = hl.import_vcf(args.vcf, force_bgz=args.force_bgz, call_fields=['GT','PGT'],
                        header_file=args.header if args.header else None)

    mt = hl.min_rep(mt, left_aligned=True)

    mt.write(get_gnomad_data_path(data_type, hail_version=0.2),
             overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input VCFs are exomes',
                        action='store_true')
    parser.add_argument('--genomes', help='Input VCFs are genomes',
                        action='store_true')
    parser.add_argument('--vcf', help='Location of VCF files.')
    parser.add_argument('--header', help='Header file if not using VCF header.')
    parser.add_argument('--min_block_size', help='Minimum file split size in MB.', type=int, default=1536)
    parser.add_argument('--force_bgz', help='Forces BGZ decoding for VCF files with .gz extension.',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite if MT exists already.',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.exomes:
        sys.exit("Exome VCFs aren't cloudable :(")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
