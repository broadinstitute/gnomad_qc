from gnomad_hail import *
from gnomad_hail.resources import sample_qc
import pandas as pd
from datetime import datetime


def write_duplicates(version: str, overwrite: bool) -> None:

    kin_ht = hl.read_table(sample_qc.relatedness_ht_path)
    kin_ht = kin_ht.filter((kin_ht.i.data_type != kin_ht.j.data_type))
    kin_ht = kin_ht.key_by(
        i=kin_ht.i.data_type + "_" + kin_ht.i.s,
        j=kin_ht.j.data_type + "_" + kin_ht.j.s
    )

    dups = get_duplicated_samples(kin_ht)
    dup_combinations = []
    dup_set_id = 0
    for dup_set in dups:
        dup_set_id += 1
        for exome in [x for x in dup_set if x.startswith("exome")]:
            for genome in [x for x in dup_set if x.startswith("genome")]:
                dup_combinations.append((dup_set_id, exome[7:], genome[8:]))

    dup_pd = pd.DataFrame(dup_combinations, columns=['dup_set', 'exomes_s', 'genomes_s'])

    rank_ht = hl.import_table(sample_qc.rank_annotations_path('joint'), impute=True)
    rank_pd = rank_ht.select(rank_ht.data_type, rank_ht.s, rank_ht.rank).to_pandas()
    dup_pd = pd.merge(dup_pd, rank_pd[rank_pd.data_type == 'exomes'], left_on=['exomes_s'], right_on=['s']).rename(columns={'rank': 'exomes_rank'})
    dup_pd = pd.merge(dup_pd, rank_pd[rank_pd.data_type == 'genomes'], left_on=['genomes_s'], right_on=['s']).rename(columns={'rank': 'genomes_rank'})

    dup_pd = dup_pd.groupby('dup_set').apply(lambda df: df.sort_values(['exomes_rank', 'genomes_rank']).reset_index(drop=True))
    dup_pd = dup_pd.reset_index(level=1).rename(columns={'level_1': 'dup_pair_rank'})
    dup_pd = dup_pd.reset_index(level=0, drop=True)

    with hl.hadoop_open(genomes_exomes_duplicate_ids_tsv_path(version), 'w' if overwrite else 'x') as out:
        dup_pd.to_csv(out, sep="\t", index=False)


def get_qc_samples_filtered_gnomad_data(data_type: str, autosomes_only: bool = True) -> hl.MatrixTable:
    mt = get_gnomad_data(data_type)
    mt = mt.filter_cols(mt.meta.high_quality)
    mt = mt.select_cols(meta=mt.meta.select('qc_platform'))
    mt = mt.select_rows(
        a_index=mt.a_index,
        was_split=mt.was_split
    )
    if autosomes_only:
        mt = filter_to_autosomes(mt)

    return mt


def compute_concordance(mt: hl.MatrixTable, other_mt: hl.MatrixTable, name: str) -> Tuple[hl.Table, hl.Table]:
    # Filter to sites present in mt samples
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    other_mt = other_mt.filter_rows(hl.agg.any(other_mt.GT.is_non_ref()))

    summary, sample_concordance_ht, sites_concordance_ht = hl.concordance(mt, other_mt)
    logger.info(f'{name} concordance summary: %s', pformat(summary))

    return sample_concordance_ht, sites_concordance_ht


def write_omes_concordance(data_type: str, dup_version: str, by_platform: bool, overwrite: bool) -> None:
    other_data_type = 'exomes' if data_type == 'genomes' else 'genomes'

    mt = get_qc_samples_filtered_gnomad_data(data_type)
    other_mt = get_qc_samples_filtered_gnomad_data(other_data_type)

    dup_ht = hl.import_table(genomes_exomes_duplicate_ids_tsv_path(dup_version), impute=True)
    dup_ht = dup_ht.filter(dup_ht.dup_pair_rank == 0)

    # Unify sample names based on inferred duplicates (from pc_relate)
    mt = unphase_mt(mt.filter_cols(hl.is_defined(dup_ht.key_by(f'{data_type}_s')[mt.s])))
    other_mt = other_mt.key_cols_by(s=dup_ht.key_by(f'{other_data_type}_s')[other_mt.s][f'{data_type}_s'])
    other_mt = unphase_mt(other_mt.filter_cols(hl.is_defined(other_mt.s)))

    exome_calling_intervals = hl.import_locus_intervals(exome_calling_intervals_path, skip_invalid_intervals=True)
    mt = mt.filter_rows(hl.is_defined(exome_calling_intervals[mt.locus]))
    other_mt = other_mt.filter_rows(hl.is_defined(exome_calling_intervals[other_mt.locus]))

    if by_platform:
        omes_conc = hl.read_table(annotations_ht_path(data_type, "omes_concordance"))
        omes_conc = omes_conc.transmute(
            concordance=[hl.struct(concordance_matrix=omes_conc.concordance,
                                   n_discordant=omes_conc.n_discordant,
                                   meta={'data_type': 'omes', 'platform': 'all'}
                                   )]
        )
        platforms = mt.aggregate_cols(hl.agg.collect_as_set(mt.meta.qc_platform))
        logger.info("Computing concordance by platform for platforms: {}".format(",".join([str(x) for x in platforms])))

        for platform in platforms:
            plat_mt = mt.filter_cols(mt.meta.qc_platform == platform)
            _, sites_concordance_ht = compute_concordance(plat_mt, other_mt, name=f'omes (platform: {platform})')
            omes_conc = omes_conc.annotate(
                concordance=omes_conc.concordance.append(
                    hl.struct(concordance_matrix=sites_concordance_ht[omes_conc.key].concordance,
                              n_discordant=sites_concordance_ht[omes_conc.key].n_discordant,
                              meta={'data_type': 'omes', 'platform': hl.str(platform)})
                )
            )
        omes_conc.write(annotations_ht_path(data_type, "omes_by_platform_concordance"), overwrite=overwrite)

    else:
        sample_concordance_ht, sites_concordance_ht = compute_concordance(mt, other_mt, name='omes')
        sites_concordance_ht.write(annotations_ht_path(data_type, "omes_concordance"), overwrite=overwrite)
        sample_concordance_ht.write(sample_annotations_table_path(data_type, "omes_concordance"), overwrite=overwrite)


def write_truth_concordance(data_type: str, truth_sample: str, overwrite: bool) -> None:
    sample_mapping = {
        'NA12878': {
            'exomes': 'C1975::NA12878',
            'genomes': 'G94982_NA12878'
        },
        'syndip': {
            'exomes': 'CHMI_CHMI3_Nex1',
            'genomes': 'CHMI_CHMI3_WGS1'
        }
    }

    mt = get_qc_samples_filtered_gnomad_data(data_type, autosomes_only=False)
    mt = mt.filter_cols(mt.s == sample_mapping[truth_sample][data_type])
    mt = unphase_mt(mt)
    mt = mt.key_cols_by(s=hl.str(truth_sample))
    mt = mt.repartition(1000 if data_type == 'genomes' else 100, shuffle=False)

    truth_mt = hl.read_matrix_table(NA12878_mt_path() if truth_sample == 'NA12878' else syndip_mt_path())
    truth_mt = truth_mt.key_cols_by(s=hl.str(truth_sample))
    if data_type == 'exomes':
        exome_calling_intervals = hl.import_locus_intervals(exome_calling_intervals_path, skip_invalid_intervals=True)
        truth_mt = truth_mt.filter_rows(hl.is_defined(exome_calling_intervals[truth_mt.locus]))
    truth_mt = unphase_mt(hl.split_multi_hts(truth_mt, left_aligned=False))

    sample_concordance_ht, sites_concordance_ht = compute_concordance(mt, truth_mt, name=truth_sample)
    sites_concordance_ht.write(annotations_ht_path(data_type, f'{truth_sample}_concordance'), overwrite=overwrite)
    sample_concordance_ht.write(sample_annotations_table_path(data_type, f'{truth_sample}_concordance'), overwrite=overwrite)


def main(args):

    data_type = 'genomes' if args.genomes else 'exomes'

    if args.debug:
        logger.setLevel(logging.DEBUG)

    if args.write_duplicates:
        write_duplicates(datetime.today().strftime('%Y-%m-%d'), args.overwrite)

    if args.na12878:
        write_truth_concordance(data_type, 'NA12878', overwrite=args.overwrite)

    if args.syndip:
        write_truth_concordance(data_type, 'syndip', args.overwrite)

    if args.omes:
        write_omes_concordance(data_type,
                               dup_version=datetime.today().strftime('%Y-%m-%d') if args.write_duplicates else args.dup_version,
                               by_platform=False,
                               overwrite=args.overwrite)

    if args.omes_by_platform:
        write_omes_concordance(data_type,
                               dup_version=datetime.today().strftime('%Y-%m-%d') if args.write_duplicates else args.dup_version,
                               by_platform=True,
                               overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    dup_grp = parser.add_argument_group('Create exomes/genomes duplicate file')
    dup_grp.add_argument('--write_duplicates', help='Write a new exomes/genomes duplicate file.', action='store_true')

    conc_grp = parser.add_argument_group('Compute concordance')
    conc_grp.add_argument('--omes', help='Computes omes concordance.', action='store_true')
    conc_grp.add_argument('--omes_by_platform', help='Whether to also run by exome platform', action='store_true')
    conc_grp.add_argument('--na12878', help='Computes NA12878 concordance.', action='store_true')
    conc_grp.add_argument('--syndip', help='Computes Syndip concordance.', action='store_true')
    conc_grp.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    conc_grp.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    conc_grp.add_argument('--dup_version', help="Version of the exomes/genomes duplicates. (default is current version). Note that if --write_duplicates is used, then the new file generated is used.", default=CURRENT_DUPS)

    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if (args.omes or args.omes_by_platform or args.na12878 or args.syndip) and int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified when running concordance.')

    if (args.omes_by_platform and
            not args.omes and
            not hl.utils.hadoop_exists(annotations_ht_path('genomes' if args.genomes else 'exomes', 'omes_concordance'))):
        sys.exit("Error: Computing --omes_by_platform requires computing --omes first. Please run --omes first (you can run with both options in a single run)")

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
