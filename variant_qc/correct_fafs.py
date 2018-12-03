from gnomad_qc import *
from gnomad_qc.variant_qc.prepare_data_release import *
from gnomad_hail import *


def make_faf_ac_struct(freq_ht, freq_index_dict, faf_index_dict):
    faf_ac_dict = {}
    for x in faf_index_dict.keys():
        faf_ac_dict.update({x: freq_ht.freq[freq_index_dict[x]].AC})
    return faf_ac_dict

def make_faf_struct(freq_ht, faf_index_dict):
    faf_dict = {}
    for x in faf_index_dict.keys():
        faf_dict.update({x: freq_ht.faf[faf_index_dict[x]].faf95})
    return faf_dict


def main(args):
    hl.init(log='/release.log')

    data_type = 'genomes' if args.genomes else 'exomes'
    freq_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies')).drop('project_max')

    freq_index_dict = make_index_dict(freq_ht)
    faf_index_dict = make_faf_index_dict(freq_ht) # {'adj_afr': 2, 'adj_amr': 5, 'adj_eas': 3, 'adj_nfe': 1, 'adj_sas': 4, 'adj': 0}

    # For any population where there is only a singleton, the faf for that population should be 0
    # Get pop names where ac_pop == 1; then match faf entry metas to these pop names and set fafs to 0
    freq_ht = freq_ht.annotate(faf_ac=hl.Struct(**make_faf_ac_struct(freq_ht, freq_index_dict, faf_index_dict)))
    freq_ht = freq_ht.annotate(ac1_pops=hl.zip_with_index(freq_ht.freq.map(lambda x: x.AC == 1))
                               .filter(lambda x: x[1]).map(lambda x: x[0])
                               .map(lambda x: hl.set(freq_ht.freq_meta[x].values())))
    freq_ht = freq_ht.annotate(orig_faf95=hl.Struct(**make_faf_struct(freq_ht, faf_index_dict)))
    freq_ht = freq_ht.annotate(faf=freq_ht.faf.map(lambda x: hl.cond(freq_ht.ac1_pops.contains(hl.set(x.meta.values())), x.annotate(faf95=0, faf99=0), x)))
    freq_ht = freq_ht.annotate(new_faf95=hl.Struct(**make_faf_struct(freq_ht, faf_index_dict)))

    # Check annotations
    freq_ht.select('faf_ac', 'orig_faf95', 'new_faf95').show()

    # Write and overwrite
    freq_ht.drop('faf_ac', 'ac1_pops', 'orig_faf95', 'new_faf95').write(annotations_ht_path(data_type, 'frequencies_faf_corrected'), args.overwrite)
    new_freq_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies_faf_corrected'))
    new_freq_ht.write(annotations_ht_path(data_type, 'frequencies'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
