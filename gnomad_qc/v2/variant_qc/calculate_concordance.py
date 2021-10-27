from gnomad.sample_qc.relatedness import get_duplicated_samples
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.annotations import unphase_call_expr
from gnomad.variant_qc.evaluation import add_rank
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources import sample_qc
from gnomad_qc.v2.resources.variant_qc import *
import pandas as pd
from datetime import datetime
import argparse
import sys
import hail as hl
import logging
from pprint import pformat

logger = logging.getLogger("calculate_concordance")


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
    mt = mt.filter_cols(hl.is_defined(dup_ht.key_by(f'{data_type}_s')[mt.s]))
    mt = mt.annotate_entries(GT=unphase_call_expr(mt.GT))
    other_mt = other_mt.key_cols_by(s=dup_ht.key_by(f'{other_data_type}_s')[other_mt.s][f'{data_type}_s'])
    other_mt = other_mt.filter_cols(hl.is_defined(other_mt.s))
    other_mt = other_mt.annotate_entries(GT=unphase_call_expr(other_mt.GT))

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
    mt = mt.annotate_entries(GT=unphase_call_expr(mt.GT))
    mt = mt.key_cols_by(s=hl.str(truth_sample))
    mt = mt.repartition(1000 if data_type == 'genomes' else 100, shuffle=False)

    truth_mt = hl.read_matrix_table(NA12878_mt_path() if truth_sample == 'NA12878' else syndip_mt_path())
    truth_mt = truth_mt.key_cols_by(s=hl.str(truth_sample))
    if data_type == 'exomes':
        exome_calling_intervals = hl.import_locus_intervals(exome_calling_intervals_path, skip_invalid_intervals=True)
        truth_mt = truth_mt.filter_rows(hl.is_defined(exome_calling_intervals[truth_mt.locus]))
    truth_mt = hl.split_multi_hts(truth_mt, left_aligned=False)
    truth_mt = truth_mt.annotate_entries(GT=unphase_call_expr(truth_mt.GT))

    sample_concordance_ht, sites_concordance_ht = compute_concordance(mt, truth_mt, name=truth_sample)
    sites_concordance_ht.write(annotations_ht_path(data_type, f'{truth_sample}_concordance'), overwrite=overwrite)
    sample_concordance_ht.write(sample_annotations_table_path(data_type, f'{truth_sample}_concordance'), overwrite=overwrite)


def create_binned_concordance(data_type: str, truth_sample: str, metric: str, nbins: int, overwrite: bool) -> None:
    """
    Creates and writes a concordance table binned by rank (both absolute and relative) for a given data type, truth sample and metric.

    :param str data_type: One 'exomes' or 'genomes'
    :param str truth_sample: Which truth sample concordance to load
    :param str metric: One of the evaluation metrics (or a RF hash)
    :param int nbins: Number of bins for the rank
    :param bool overwrite: Whether to overwrite existing table
    :return: Nothing -- just writes the table
    :rtype: None
    """

    if hl.hadoop_exists(binned_concordance_path(data_type, truth_sample, metric)+'/_SUCCESS') and not overwrite:
        logger.warn(f"Skipping binned concordance creation as {binned_concordance_path(data_type, truth_sample, metric)} exists and overwrite=False")
    else:
        ht = hl.read_table(annotations_ht_path(data_type, f'{truth_sample}_concordance'))
        # Remove 1bp indels for syndip as cannot be trusted
        if truth_sample == 'syndip':
            ht = ht.filter(hl.is_indel(ht.alleles[0], ht.alleles[1]) & (hl.abs(hl.len(ht.alleles[0]) - hl.len(ht.alleles[1])) == 1), keep=False)
            high_conf_intervals = hl.import_locus_intervals(syndip_high_conf_regions_bed_path)
        else:
            high_conf_intervals = hl.import_locus_intervals(NA12878_high_conf_regions_bed_path)

        lcr = hl.import_locus_intervals(lcr_intervals_path)
        segdup = hl.import_locus_intervals(segdup_intervals_path)
        ht = ht.filter(
            hl.is_defined(high_conf_intervals[ht.locus]) &
            hl.is_missing(lcr[ht.locus]) &
            hl.is_missing(segdup[ht.locus])
        )

        if metric in ['vqsr', 'rf_2.0.2', 'rf_2.0.2_beta', 'cnn']:
            metric_ht = hl.read_table(score_ranking_path(data_type, metric))
        else:
            metric_ht = hl.read_table(rf_path(data_type, 'rf_result', run_hash=metric))

        metric_snvs, metrics_indels = metric_ht.aggregate([hl.agg.count_where(hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1])),
                                                           hl.agg.count_where(~hl.is_snp(metric_ht.alleles[0], metric_ht.alleles[1]))])

        snvs, indels = ht.aggregate([hl.agg.count_where(hl.is_snp(ht.alleles[0], ht.alleles[1])),
                                     hl.agg.count_where(~hl.is_snp(ht.alleles[0], ht.alleles[1]))])

        ht = ht.annotate_globals(global_counts=hl.struct(snvs=metric_snvs, indels=metrics_indels),
                                 counts=hl.struct(snvs=snvs, indels=indels))

        ht = ht.annotate(snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
                         score=metric_ht[ht.key].score,
                         global_rank=metric_ht[ht.key].rank,
                         # TP => allele is found in both data sets
                         n_tp=ht.concordance[3][3] + ht.concordance[3][4] + ht.concordance[4][3] + ht.concordance[4][4],
                         # FP => allele is found only in test data set
                         n_fp=hl.sum(ht.concordance[3][:2]) + hl.sum(ht.concordance[4][:2]),
                         # FN => allele is found only in truth data set
                         n_fn=hl.sum(ht.concordance[:2].map(lambda x: x[3] + x[4]))
                         )

        ht = add_rank(ht, -1.0*ht.score)

        ht = ht.annotate(rank=[
            hl.tuple(['global_rank', (ht.global_rank + 1) / hl.cond(
                             ht.snv,
                             ht.globals.global_counts.snvs,
                             ht.globals.global_counts.indels
                         )]),
            hl.tuple(['truth_sample_rank', (ht.rank + 1) / hl.cond(
                             ht.snv,
                             ht.globals.counts.snvs,
                             ht.globals.counts.indels
                         )])
        ])

        ht = ht.explode(ht.rank)
        ht = ht.annotate(
            rank_name=ht.rank[0],
            bin=hl.int(ht.rank[1] * nbins)
        )

        ht = ht.group_by(
            'rank_name',
            'snv',
            'bin'
        ).aggregate(
            # Look at site-level metrics -> tp > fp > fn -- only important for multi-sample comparisons
            tp=hl.agg.count_where(ht.n_tp > 0),
            fp=hl.agg.count_where((ht.n_tp == 0) & (ht.n_fp > 0)),
            fn=hl.agg.count_where((ht.n_tp == 0) & (ht.n_fp == 0) & (ht.n_fn > 0)),
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n_alleles=hl.agg.count()
        ).repartition(5)

        ht.write(binned_concordance_path(data_type, truth_sample, metric), overwrite=overwrite)


def main(args):

    if args.debug:
        logger.setLevel(logging.DEBUG)

    if args.write_duplicates:
        logger.info(f"Writing duplicates")
        write_duplicates(datetime.today().strftime('%Y-%m-%d'), args.overwrite)

    data_types = []
    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:

        if args.compute_concordance:
            if args.na12878:
                logger.info(f"Computing {data_type} / NA12878 concordance")
                write_truth_concordance(data_type, 'NA12878', overwrite=args.overwrite)

            if args.syndip:
                logger.info(f"Computing {data_type} / syndip concordance")
                write_truth_concordance(data_type, 'syndip', args.overwrite)

            if args.omes:
                logger.info("Computing {} / {} concordance".format(data_type, 'exomes' if data_type == 'genomes' else 'exomes'))
                write_omes_concordance(data_type,
                                       dup_version=datetime.today().strftime('%Y-%m-%d') if args.write_duplicates else args.dup_version,
                                       by_platform=False,
                                       overwrite=args.overwrite)

            if args.omes_by_platform:
                logger.info("Computing {} / {} by platform concordance".format(data_type, 'exomes' if data_type == 'genomes' else 'exomes'))
                write_omes_concordance(data_type,
                                       dup_version=datetime.today().strftime('%Y-%m-%d') if args.write_duplicates else args.dup_version,
                                       by_platform=True,
                                       overwrite=args.overwrite)

        if args.bin_concordance:

            truth_samples = []
            if args.na12878:
                truth_samples.append('NA12878')
            if args.syndip:
                truth_samples.append('syndip')

            metrics = [] if not args.run_hash else [args.run_hash] if isinstance(args.run_hash, str) else args.run_hash
            if args.vqsr:
                metrics.append('vqsr')
            if args.rf_2_0_2:
                metrics.append('rf_2.0.2')
            if args.rf_2_0_2_beta:
                metrics.append('rf_2.0.2_beta')

            for truth_sample in truth_samples:
                for metric in metrics:
                    logger.info(f"Creating binned concordance table for {data_type} / {truth_sample} for metric {metric}")
                    create_binned_concordance(data_type, truth_sample, metric, args.n_bins, args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')

    action_grp = parser.add_argument_group('Actions')
    action_grp.add_argument('--write_duplicates', help='Write a new exomes/genomes duplicate file.', action='store_true')
    action_grp.add_argument('--compute_concordance', help='Computes concordance for each of the dataset specified separately.', action='store_true')
    action_grp.add_argument('--bin_concordance', help='Merges individual concordance results with specified metrics ranks.', action='store_true')

    data_grp = parser.add_argument_group('Datasets available')
    data_grp.add_argument('--omes', help='Exomes/genomes (--compute_concordance only)', action='store_true')
    data_grp.add_argument('--omes_by_platform', help='Exomes/genomes by platform (--compute_concordance only)', action='store_true')
    data_grp.add_argument('--na12878', help='NA12878 ', action='store_true')
    data_grp.add_argument('--syndip', help='Syndip', action='store_true')

    conc_grp = parser.add_argument_group('Compute concordance options')
    conc_grp.add_argument('--dup_version', help="Version of the exomes/genomes duplicates. (default is current version). Note that if --write_duplicates is used, then the new file generated is used.", default=CURRENT_DUPS)

    bin_grp = parser.add_argument_group('Bin concordance options')
    bin_grp.add_argument('--run_hash', help='RF hash(es) for annotation.', nargs='+')
    bin_grp.add_argument('--vqsr', help='When set, annotates with VQSR rank file.', action='store_true')
    bin_grp.add_argument('--rf_2_0_2', help='When set, annotates with 2.0.2 RF rank file.', action='store_true')
    bin_grp.add_argument('--rf_2_0_2_beta', help='When set, annotates with 2.0.2 RF beta (no medians, QD / max(pAB)) rank file.', action='store_true')
    bin_grp.add_argument('--n_bins', help='Number of bins for the binned file (default: 100)', default=100, type=int)

    args = parser.parse_args()

    if not (args.write_duplicates or args.compute_concordance or args.bin_concordance):
        sys.exit('Error: You must specify at least one action.')

    if (args.compute_concordance or args.bin_concordance) and not (args.exomes or args.genomes):
        sys.exit('Error: At least one of --exomes or --genomes must be specified when using --compute_concordance or --bin_concordance.')

    if args.compute_concordance and not(args.na12878 or args.syndip or args.omes or args.omes_by_platform):
        sys.exit('Error At least one dataset should be specified when running --compute_concordance.')

    if args.bin_concordance and not(args.na12878 or args.syndip):
        sys.exit('Error: At least one dataset should be specified when running --bin_concordance. (At the moment --omes is not implemented)')

    if (args.omes_by_platform and
            not args.omes and
            not hl.utils.hadoop_exists(annotations_ht_path('genomes' if args.genomes else 'exomes', 'omes_concordance'))):
        sys.exit("Error: Computing --omes_by_platform requires computing --omes first. Please run --omes first (you can run with both options in a single run)")

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
