from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.variant_qc import *
from gnomad_qc.v2.resources.sample_qc import *
import argparse
import sys
import logging

logger = logging.getLogger("topmed_dups")


def create_shared_sites_table(data_type: str, overwrite: bool):
    freq_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    topmed_ht = hl.read_matrix_table('gs://gnomad-public/resources/hail-0.2/topmed.b37.mt').rows()

    methyation_ht = hl.read_table(methylation_sites_ht_path())

    gnomad = get_gnomad_data(data_type, non_refs_only=True, release_samples=True)
    gnomad = gnomad.select_cols(**gnomad.meta)
    gnomad = filter_to_autosomes(gnomad)

    rf_path = annotations_ht_path(data_type, 'rf')
    if hl.hadoop_exists(annotations_ht_path(data_type, 'rf')):
        logger.info(f"Filtering sites based on {rf_path}")
        filter_ht = hl.read_table(rf_path)
        gnomad = gnomad.filter_rows(hl.len(filter_ht[gnomad.row_key].filters) == 0)
    else:
        logger.warn(f"Could not find filtering table {rf_path}. Not filtering poor quality sites leads to lower performance.")

    gnomad = gnomad.filter_rows(hl.is_snp(gnomad.alleles[0], gnomad.alleles[1]) &
                                hl.or_else(methyation_ht[gnomad.locus].MEAN < 0.6, True)
                                )

    gnomad = gnomad.annotate_rows(topmed_ac=topmed_ht[gnomad.row_key].info.AC[0],
                                  gnomad_ac=freq_ht[gnomad.row_key].freq[0].AC[1])
    gnomad.annotate_cols(n_singletons=hl.agg.count_where((gnomad.gnomad_ac == 1) & gnomad.GT.is_het()),
                                n_doubletons=hl.agg.count_where((gnomad.gnomad_ac < 3) & gnomad.GT.is_het()),
                                n_tripletons=hl.agg.count_where((gnomad.gnomad_ac < 4) & gnomad.GT.is_het()),
                                n_topmed_singletons=hl.agg.count_where((gnomad.gnomad_ac == 1) & gnomad.GT.is_het() & hl.is_defined(gnomad.topmed_ac)),
                                n_topmed_doubletons=hl.agg.count_where((gnomad.gnomad_ac < 3) & gnomad.GT.is_het() & hl.is_defined(gnomad.topmed_ac)),
                                n_topmed_tripletons=hl.agg.count_where((gnomad.gnomad_ac < 4) & gnomad.GT.is_het() & hl.is_defined(gnomad.topmed_ac)),
                                n_both_singletons=hl.agg.count_where((gnomad.gnomad_ac == 1) & gnomad.GT.is_het() & (gnomad.topmed_ac == 1))
                                ).cols().write(get_topmed_shared_sites_ht_path(data_type), overwrite=overwrite)


def main(args):
    data_type = 'exomes' if args.exomes else 'genomes'

    create_shared_sites_table(data_type, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
