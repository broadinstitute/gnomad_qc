from gnomad_hail import try_slack
import hail as hl
import argparse
from gnomad_qc.v3.resources import qual_hists_json_path, internal_ht_path


def define_hist_ranges(ht):
    hist_dict = {
        'FS': hl.agg.hist(ht.info.FS, 0, 50, 50),
        'InbreedingCoeff': hl.agg.hist(ht.info.InbreedingCoeff, -0.25, 0.25, 100),
        'MQ': hl.agg.hist(ht.info.MQ, 20, 60, 40),
        'RAW_MQ': hl.agg.hist(hl.log10(ht.info.RAW_MQ), 2, 13, 33),
        'MQRankSum': hl.agg.hist(ht.info.MQRankSum, -15, 15, 60),
        'QD': hl.agg.hist(ht.info.QD, 0, 40, 40),
        'ReadPosRankSum': hl.agg.hist(ht.info.ReadPosRankSum, -15, 15, 60),
        'SOR': hl.agg.hist(ht.info.SOR, 0, 10, 50),
        'VarDP': hl.agg.hist(hl.log10(ht.info.VarDP), 1, 9, 32),
        'AS_VQSLOD': hl.agg.hist(ht.info.AS_VQSLOD, -30, 30, 60)
    }
    return hist_dict


def aggregate_qual_stats_by_bin(ht):
    ht = ht.annotate(metric=(
        hl.case()
            .when(ht.freq[1].AC == 1, "binned_singleton")
            .when(ht.freq[1].AC == 2, "binned_doubleton")
            .when(ht.freq[1].AF < 0.00005, "binned_0.00005")
            .when(ht.freq[1].AF < 0.0001, "binned_0.0001")
            .when(ht.freq[1].AF < 0.0002, "binned_0.0002")
            .when(ht.freq[1].AF < 0.0005, "binned_0.0005")
            .when(ht.freq[1].AF < 0.001, "binned_0.001")
            .when(ht.freq[1].AF < 0.002, "binned_0.002")
            .when(ht.freq[1].AF < 0.005, "binned_0.005")
            .when(ht.freq[1].AF < 0.01, "binned_0.01")
            .when(ht.freq[1].AF < 0.02, "binned_0.02")
            .when(ht.freq[1].AF < 0.05, "binned_0.05")
            .when(ht.freq[1].AF < 0.1, "binned_0.1")
            .when(ht.freq[1].AF < 0.2, "binned_0.2")
            .when(ht.freq[1].AF < 0.5, "binned_0.5")
            .when(ht.freq[1].AF <= 1, "binned_1")
            .default(hl.null(hl.tstr))
    ))
    return ht


def main(args):
    hl.init(default_reference='GRCh38', log='/variant_histograms.log')

    metrics = ['FS', 'InbreedingCoeff', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR', 'BaseQRankSum',
               'ClippingRankSum', 'DP', 'VQSLOD', 'rf_tp_probability', 'pab_max']

    ht = hl.read_table(internal_ht_path)
    # NOTE: histogram aggregations are done on the entire callset (not just PASS variants), on raw data

    # NOTE: run the following code in a first pass to determine bounds for metrics
    # Evaluate minimum and maximum values for each metric of interest
    if args.first_pass:
        minmax_dict = {}
        for metric in metrics:
            minmax_dict[metric] = hl.struct(min=hl.agg.min(ht[metric]), max=hl.cond(hl.agg.max(ht[metric])<1e10, hl.agg.max(ht[metric]), 1e10))
        minmax = ht.aggregate(hl.struct(**minmax_dict))
        print(minmax)
    else:
        # Aggregate hists over hand-tooled ranges
        hists = ht.aggregate(hl.struct(**define_hist_ranges(ht)))
        hist_out = hl.array([hists[f'{metric}'].annotate(metric=metric) for metric in hists])

        # Aggregate QUAL stats by bin:
        ht = aggregate_qual_stats_by_bin(ht)
        ht = ht.group_by('metric').aggregate(hist=hl.agg.hist(hl.log10(ht.info.QUALapprox), 1, 10, 36))
        hists = ht.collect()
        hist_out = hist_out.extend(hl.array([x.hist.annotate(metric=x.metric) for x in hists]))

        with hl.hadoop_open(qual_hists_json_path, 'w') as f:
            f.write(hl.eval(hl.json(hist_out)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_pass', help='Determine min/max values for each variant metric (to be used in hand-tooled histogram ranges', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
