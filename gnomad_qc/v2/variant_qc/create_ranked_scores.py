from gnomad.variant_qc.evaluation import add_rank
from gnomad.utils.annotations import add_variant_type
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.variant_qc import *
import argparse
import sys
import logging
from pprint import pformat

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("rank_rf")
logger.setLevel(logging.INFO)


def run_sanity_checks(ht: hl.Table) -> None:
    """
    Runs and prints sanity checks on rank table.

    :param Table ht: input ranks Table
    :return: Nothing
    :rtype: None
    """
    print(ht.aggregate(hl.struct(was_split=hl.agg.counter(ht.was_split),
                                 has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_rank)),
                                 was_singleton=hl.agg.counter(ht._singleton),
                                 has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_rank)),
                                 was_split_singleton=hl.agg.counter(ht._singleton & ~ht.was_split),
                                 has_biallelic_singleton_rank=hl.agg.counter(hl.is_defined(ht.biallelic_singleton_rank)))))


def create_rf_rank(data_type: str, run_hash: str) -> None:
    """
    Creates a ranked table for a RF run and writes it to its correct location in annotations.

    :param str data_type: One of 'exomes' or 'genomes'
    :param str run_hash: RF run hash
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Creating rank file for {data_type} RF run {run_hash}")

    if not hl.hadoop_exists(f'gs://gnomad-tmp/gnomad_{data_type}_rf_{run_hash}.ht/_SUCCESS'):
        gnomad_ht = get_gnomad_annotations(data_type)
        ht = hl.read_table(rf_path(data_type, 'rf_result', run_hash=run_hash))
        ht = ht.annotate(**gnomad_ht[ht.key],
                         score=ht.rf_probability['TP'])

        # Write to temp location as result will be overwritten
        ht.write(f'gs://gnomad-tmp/gnomad_{data_type}_rf_{run_hash}.ht', overwrite=True)
    ht = hl.read_table(f'gs://gnomad-tmp/gnomad_{data_type}_rf_{run_hash}.ht')

    ht = add_rank(ht,
                  score_expr=1-ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac > 0)
                  }
                  )
    ht.write(rf_path(data_type, 'rf_result', run_hash=run_hash), overwrite=True)


def create_vqsr_rank_ht(data_type: str):
    """
    Creates a rank table for VQSR and writes it to its correct location in annotations.

    :param str data_type: One of 'exomes' or 'genomes'
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Creating rank file for {data_type} VQSR")
    if not hl.utils.hadoop_exists(f'gs://gnomad-tmp/gnomad_{data_type}_vqsr_pre_ranking.ht/_SUCCESS'):
        ht = get_gnomad_data(data_type).rows()
        gnomad_ht = get_gnomad_annotations(data_type)
        logger.info('Filtering to high_quality samples and n_nonref==1...')
        ht = ht.select(**gnomad_ht[ht.key],
                            score=ht.info.VQSLOD,
                            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
                            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
                            culprit=ht.info.culprit,
                            )
        ht = ht.filter(ht.n_nonref > 0)
        ht = ht.repartition(1000, shuffle=False)
        ht.write(f'gs://gnomad-tmp/gnomad_{data_type}_vqsr_pre_ranking.ht', overwrite=True)
    ht = hl.read_table(f'gs://gnomad-tmp/gnomad_{data_type}_vqsr_pre_ranking.ht')

    ht = add_rank(ht,
                  score_expr=-1*ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac > 0)
                  }
                  )
    ht.write(score_ranking_path(data_type, 'vqsr'), overwrite=True)


def get_gnomad_annotations(data_type: str) -> hl.Table:
    """
    Gets a Table with gnomad annotations needed to create ranked tables.
    Because this is computationally expensive, the table is written the first time it is created and then loaded after that.

    :param str data_type: One of 'exomes' or 'genomes'
    :return: Table with annotations
    :rtype: Table
    """
    if not hl.utils.hadoop_exists(f'gs://gnomad/variant_qc/temp/gnomad_{data_type}_annotations.ht'):
        ht = get_gnomad_data(data_type).rows()
        call_stats_data = hl.read_table(annotations_ht_path(data_type, 'call_stats'))
        freq_data = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
        ht = ht.select(
            n_nonref=call_stats_data[ht.key].qc_callstats[0].AC[1],
            singleton=ht.info.AC[ht.a_index - 1] == 1,
            info_ac=ht.info.AC[ht.a_index - 1],
            was_split=ht.was_split,
            ac=freq_data[ht.key].freq[0].AC[1],
            ac_raw=freq_data[ht.key].freq[1].AC[1]
        )
        ht.repartition(1000, shuffle=False).write(f'gs://gnomad/variant_qc/temp/gnomad_{data_type}_annotations.ht')
    return hl.read_table(f'gs://gnomad/variant_qc/temp/gnomad_{data_type}_annotations.ht')


def create_cnn_rank_file() -> None:
    """
    Creates a rank file for the CNN data and writes it to its correct location.

    :return: Nothing
    :rtype: None
    """

    logger.info("Creating CNN rank file.")

    if not hl.utils.hadoop_exists('gs://gnomad/variant_qc/temp/friedman_cnn_scores.ht/_SUCCESS'):
        logger.info(f"Importing CNN scores")
        ht = hl.import_table('gs://gnomad/variant_qc/temp/friedman_cnn_scores.tsv.bgz', min_partitions=1000, impute=True)
        ht.write('gs://gnomad/variant_qc/temp/friedman_cnn_scores.ht', overwrite=True)

    logger.info('Formatting CNN scores...')
    ht = hl.read_table('gs://gnomad/variant_qc/temp/friedman_cnn_scores.ht')
    ht = ht.annotate(alt_alleles=ht.Alt.split(','), was_split=ht.Alt.contains(','))  # This transforms to a list
    ht = ht.explode('alt_alleles')
    ht = ht.annotate(locus=hl.locus(hl.str(ht.Contig), ht.Pos))

    # Minrep
    ht = ht.annotate(**hl.min_rep(ht.locus, [ht.Ref, ht.alt_alleles]))
    ht = ht.annotate(vartype=add_variant_type(ht.alleles))
    ht = ht.key_by('locus', 'alleles')

    # Select relevant annotations and add singleton / n_nonref, was_split
    gnomad_ht = get_gnomad_annotations('genomes')
    ht = ht.select(**gnomad_ht[ht.key],
                   variant_type=ht.vartype.variant_type,
                   n_alt_alleles=ht.vartype.n_alt_alleles,
                   score=ht.CNN_1D_Score,
                   )
    ht = ht.filter(ht.n_nonref > 0)
    ht = add_rank(ht,
                  score_expr=-1*ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac > 0)
                  }
                  )

    ht.write(score_ranking_path('genomes', 'cnn'), overwrite=True)


def create_rf_2_0_2_rank(data_type: str, beta: bool) -> None:
    """
    Creates a rank file for 2.0.2 RF and writes it to its correct location.

    :param str data_type: One of 'exomes' or 'genomes'
    :param bool beta: If set, then creates the table for the "beta" 2.0.2 RF with QD / max(p(AB))
    :return: Nothing
    :rtype: None
    """
    logger.info(f"Creating rank file for {data_type} RF 2.0.2{'beta' if beta else ''}")

    if not hl.hadoop_exists(f'gs://gnomad-tmp/gnomad_rf_2_0_2_{data_type}_{str(beta)}_tmp.ht'):
        ht = hl.import_table(get_2_0_2_rf_path(data_type, beta), types={'chrom': hl.tstr}, impute=True, min_partitions=1000)
        if 'chrom' in ht.row:
            ht = ht.transmute(
                locus=hl.locus(ht.chrom, ht.pos),
                alleles=[ht.ref, ht.alt]
            )
        else:
            ht = ht.transmute(
                v=hl.parse_variant(ht.v),
                rfprob=ht.rf_rpob_tp # Yes, this is awful
            )
            ht = ht.transmute(
                locus=ht.v.locus,
                alleles=ht.v.alleles
            )

        ht = ht.key_by('locus', 'alleles')

        gnomad_ht = get_gnomad_annotations(data_type)
        ht = ht.annotate(
            **gnomad_ht[ht.key],
            score=ht.rfprob
        )

        ht.write(f'gs://gnomad-tmp/gnomad_rf_2_0_2_{data_type}_{str(beta)}_tmp.ht')
    ht = hl.read_table(f'gs://gnomad-tmp/gnomad_rf_2_0_2_{data_type}_{str(beta)}_tmp.ht')
    ht = add_rank(ht,
                  score_expr=1-ht.score,
                  subrank_expr={
                      'singleton_rank': ht.singleton,
                      'biallelic_rank': ~ht.was_split,
                      'biallelic_singleton_rank': ~ht.was_split & ht.singleton,
                      'adj_rank': ht.ac > 0,
                      'adj_biallelic_rank': ~ht.was_split & (ht.ac > 0),
                      'adj_singleton_rank': ht.singleton & (ht.ac > 0),
                      'adj_biallelic_singleton_rank': ~ht.was_split & ht.singleton & (ht.ac > 0)
                  }
                  )

    ht.write(score_ranking_path(data_type, 'rf_2.0.2{}'.format('_beta' if beta else '')), overwrite=True)


def create_binned_data(ht: hl.Table, data: str, data_type: str, n_bins: int) -> hl.Table:
    """
    Creates binned data from a rank Table grouped by rank_id (rank, biallelic, etc.), contig, snv, bi_allelic and singleton
    containing the information needed for evaluation plots.

    :param Table ht: Input rank table
    :param str data: Which data/run hash is being created
    :param str data_type: one of 'exomes' or 'genomes'
    :param int n_bins: Number of bins.
    :return: Binned Table
    :rtype: Table
    """

    # Count variants for ranking
    count_expr = {x: hl.agg.filter(hl.is_defined(ht[x]), hl.agg.counter(hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel'))) for x in ht.row if x.endswith('rank')}
    rank_variant_counts = ht.aggregate(hl.Struct(**count_expr))
    logger.info(f"Found the following variant counts:\n {pformat(rank_variant_counts)}")
    ht = ht.annotate_globals(rank_variant_counts=rank_variant_counts)

    # Load external evaluation data
    clinvar_ht = hl.read_table(clinvar_ht_path)
    denovo_ht = get_validated_denovos_ht()
    if data_type == 'exomes':
        denovo_ht = denovo_ht.filter(denovo_ht.gnomad_exomes.high_quality)
    else:
        denovo_ht = denovo_ht.filter(denovo_ht.gnomad_genomes.high_quality)
    denovo_ht = denovo_ht.select(validated_denovo=denovo_ht.validated, high_confidence_denovo=denovo_ht.Confidence == 'HIGH')
    ht_truth_data = hl.read_table(annotations_ht_path(data_type, 'truth_data'))
    fam_ht = hl.read_table(annotations_ht_path(data_type, 'family_stats'))
    fam_ht = fam_ht.select(
        family_stats=fam_ht.family_stats[0]
    )
    gnomad_ht = get_gnomad_data(data_type).rows()
    gnomad_ht = gnomad_ht.select(
        vqsr_negative_train_site=gnomad_ht.info.NEGATIVE_TRAIN_SITE,
        vqsr_positive_train_site=gnomad_ht.info.POSITIVE_TRAIN_SITE,
        fail_hard_filters=(gnomad_ht.info.QD < 2) | (gnomad_ht.info.FS > 60) | (gnomad_ht.info.MQ < 30)
    )
    lcr_intervals = hl.import_locus_intervals(lcr_intervals_path)

    ht = ht.annotate(
        **ht_truth_data[ht.key],
        **fam_ht[ht.key],
        **gnomad_ht[ht.key],
        **denovo_ht[ht.key],
        clinvar=hl.is_defined(clinvar_ht[ht.key]),
        indel_length=hl.abs(ht.alleles[0].length()-ht.alleles[1].length()),
        rank_bins=hl.array(
            [hl.Struct(
                rank_id=rank_name,
                bin=hl.int(hl.ceil(hl.float(ht[rank_name] + 1) / hl.floor(ht.globals.rank_variant_counts[rank_name][hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), 'snv', 'indel')] / n_bins)))
            )
                for rank_name in rank_variant_counts]
        ),
        lcr=hl.is_defined(lcr_intervals[ht.locus])
    )

    ht = ht.explode(ht.rank_bins)
    ht = ht.transmute(
        rank_id=ht.rank_bins.rank_id,
        bin=ht.rank_bins.bin
    )
    ht = ht.filter(hl.is_defined(ht.bin))

    ht = ht.checkpoint(f'gs://gnomad-tmp/gnomad_score_binning_{data_type}_tmp_{data}.ht', overwrite=True)

    # Create binned data
    return (
        ht
        .group_by(
            rank_id=ht.rank_id,
            contig=ht.locus.contig,
            snv=hl.is_snp(ht.alleles[0], ht.alleles[1]),
            bi_allelic=hl.is_defined(ht.biallelic_rank),
            singleton=ht.singleton,
            release_adj=ht.ac>0,
            bin=ht.bin
        )._set_buffer_size(20000)
        .aggregate(
            min_score=hl.agg.min(ht.score),
            max_score=hl.agg.max(ht.score),
            n=hl.agg.count(),
            n_ins=hl.agg.count_where(hl.is_insertion(ht.alleles[0], ht.alleles[1])),
            n_del=hl.agg.count_where(hl.is_deletion(ht.alleles[0], ht.alleles[1])),
            n_ti=hl.agg.count_where(hl.is_transition(ht.alleles[0], ht.alleles[1])),
            n_tv=hl.agg.count_where(hl.is_transversion(ht.alleles[0], ht.alleles[1])),
            n_1bp_indel=hl.agg.count_where(ht.indel_length == 1),
            n_mod3bp_indel=hl.agg.count_where((ht.indel_length % 3) == 0),
            n_clinvar=hl.agg.count_where(ht.clinvar),
            n_singleton=hl.agg.count_where(ht.singleton),
            n_validated_de_novos=hl.agg.count_where(ht.validated_denovo),
            n_high_confidence_de_novos=hl.agg.count_where(ht.high_confidence_denovo),
            n_de_novo=hl.agg.filter(ht.family_stats.unrelated_qc_callstats.AC[1] == 0, hl.agg.sum(ht.family_stats.mendel.errors)),
            n_de_novo_no_lcr=hl.agg.filter(~ht.lcr & (ht.family_stats.unrelated_qc_callstats.AC[1] == 0), hl.agg.sum(ht.family_stats.mendel.errors)),
            n_de_novo_sites=hl.agg.filter(ht.family_stats.unrelated_qc_callstats.AC[1] == 0, hl.agg.count_where(ht.family_stats.mendel.errors > 0)),
            n_de_novo_sites_no_lcr=hl.agg.filter(~ht.lcr & (ht.family_stats.unrelated_qc_callstats.AC[1] == 0), hl.agg.count_where(ht.family_stats.mendel.errors > 0)),
            n_trans_singletons=hl.agg.filter((ht.info_ac < 3) & (ht.family_stats.unrelated_qc_callstats.AC[1] == 1), hl.agg.sum(ht.family_stats.tdt.t)),
            n_untrans_singletons=hl.agg.filter((ht.info_ac < 3) & (ht.family_stats.unrelated_qc_callstats.AC[1] == 1), hl.agg.sum(ht.family_stats.tdt.u)),
            n_train_trans_singletons=hl.agg.count_where((ht.family_stats.unrelated_qc_callstats.AC[1] == 1) & (ht.family_stats.tdt.t == 1)),
            n_omni=hl.agg.count_where(ht.truth_data.omni),
            n_mills=hl.agg.count_where(ht.truth_data.mills),
            n_hapmap=hl.agg.count_where(ht.truth_data.hapmap),
            n_kgp_high_conf_snvs=hl.agg.count_where(ht.truth_data.kgp_high_conf_snvs),
            fail_hard_filters=hl.agg.count_where(ht.fail_hard_filters),
            n_vqsr_pos_train=hl.agg.count_where(ht.vqsr_positive_train_site),
            n_vqsr_neg_train=hl.agg.count_where(ht.vqsr_negative_train_site)
        )
    )


def main(args):

    def run_data(rank_func: Callable[..., None], rank_func_args: List[Any], score_type: str, data_type: str) -> None:
        """
        Wrapper for running script actions on given data.

        :param callable rank_func: Function for creating ranking file
        :param list of Any rank_func_args: Arguments to pass to the ranking function
        :param str score_type: Score being processed
        :param str data_type: One of 'exomes' or 'genomes'
        :return: Nothing -- this runs the script actions
        :rtype:  None
        """
        if score_type in ['vqsr', 'cnn', 'rf_2.0.2', 'rf_2.0.2_beta']:
            rank_file_path = score_ranking_path(data_type, score_type)
        else:
            rank_file_path = rf_path(data_type, 'rf_result', run_hash=args.run_hash)

        if args.create_rank_file:
            rank_func(*rank_func_args)
        if args.run_sanity_checks:
            run_sanity_checks(hl.read_table(rank_file_path))
        if args.create_binned_file:
            ht = hl.read_table(rank_file_path)
            binned_ht = create_binned_data(ht, score_type, data_type, args.n_bins)
            binned_ht.write(score_ranking_path(data_type, score_type, binned=True), overwrite=True)

    hl.init(log='/create_rank.log')
    data_types = []
    if args.exomes:
        data_types.append('exomes')
    if args.genomes:
        data_types.append('genomes')

    for data_type in data_types:

        if args.run_hash:
            run_data(create_rf_rank, [data_type, args.run_hash], args.run_hash, data_type)

        if args.vqsr:
            run_data(create_vqsr_rank_ht, [data_type], 'vqsr', data_type)

        if args.cnn and data_type == 'genomes':
            run_data(create_cnn_rank_file, [], 'cnn', data_type)

        if args.rf_2_0_2:
            run_data(create_rf_2_0_2_rank, [data_type, False], 'rf_2.0.2', data_type)

        if args.rf_2_0_2_beta:
            run_data(create_rf_2_0_2_rank, [data_type, True], 'rf_2.0.2_beta', data_type)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--exomes', help='Run on exomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. At least one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--run_hash', help='Run hash for RF results to be ranked.')
    parser.add_argument('--vqsr', help='When set, creates the VQSR rank file.', action='store_true')
    parser.add_argument('--cnn', help='When set, creates the CNN rank file. Note that there is no CNN file for exomes.', action='store_true')
    parser.add_argument('--rf_2_0_2', help='When set, creates the 2.0.2 RF rank file.', action='store_true')
    parser.add_argument('--rf_2_0_2_beta', help='When set, creates the 2.0.2 RF beta (no medians, QD / max(pAB)) rank file.', action='store_true')
    parser.add_argument('--create_rank_file', help='When set, creates ranking file.', action='store_true')
    parser.add_argument('--run_sanity_checks', help='When set, runs ranking sanity checks.', action='store_true')
    parser.add_argument('--create_binned_file', help='When set, creates binned ranked file.', action='store_true')
    parser.add_argument('--n_bins', help='Number of bins for the binned file (default: 100)', default=100, type=int)
    args = parser.parse_args()

    if not args.exomes and not args.genomes:
        sys.exit('Error: At least one of --exomes or --genomes must be specified.')

    if not args.create_rank_file and not args.run_sanity_checks and not args.create_binned_file:
        sys.exit('Error: At least one of --create_rank_file, --run_sanity_checks or --create_binned_file must be specified.')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
