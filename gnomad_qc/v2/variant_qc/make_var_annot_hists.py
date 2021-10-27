from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.variant_qc import *
import argparse
import sys
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_histograms")
logger.setLevel(logging.INFO)


def define_hist_ranges(ht):
    hist_dict = {
        'FS': hl.agg.hist(ht.FS, 0, 50, 50),  # NOTE: in 2.0.2 release this was on (0,20)
        'InbreedingCoeff': hl.agg.hist(ht.InbreedingCoeff, -0.25, 0.25, 50),
        'MQ': hl.agg.hist(ht.MQ, 0, 80, 40),
        'MQRankSum': hl.agg.hist(ht.MQRankSum, -15, 15, 60),
        'QD': hl.agg.hist(ht.QD, 0, 40, 40),
        'ReadPosRankSum': hl.agg.hist(ht.ReadPosRankSum, -15, 15, 60),
        'SOR': hl.agg.hist(ht.SOR, 0, 10, 50),
        'BaseQRankSum': hl.agg.hist(ht.BaseQRankSum, -15, 15, 60),
        'ClippingRankSum': hl.agg.hist(ht.ClippingRankSum, -5, 5, 40),
        'DP': hl.agg.hist(hl.log(ht.DP, base=10), 1, 9, 32),  # NOTE: in 2.0.2 release this was on (0,8)
        'VQSLOD': hl.agg.hist(ht.VQSLOD, -30, 30, 60),  # NOTE: in 2.0.2 release this was on (-20,20)
        'rf_tp_probability': hl.agg.hist(ht.rf_tp_probability, 0, 1, 50),
        'pab_max': hl.agg.hist(ht.pab_max, 0, 1, 50)
    }
    return hist_dict


def aggregate_qual_stats_by_bin(ht):
    ht = ht.annotate(metric=(hl.case()
                             .when(ht.gnomad_AC_raw == 1, "binned_singleton")
                             .when(ht.gnomad_AC_raw == 2, "binned_doubleton")
                             .when((ht.gnomad_AC_raw > 2) & (ht.gnomad_AF_raw < 0.00005), "binned_0.00005")
                             .when((ht.gnomad_AF_raw >= 0.00005) & (ht.gnomad_AF_raw < 0.0001), "binned_0.0001")
                             .when((ht.gnomad_AF_raw >= 0.0001) & (ht.gnomad_AF_raw < 0.0002), "binned_0.0002")
                             .when((ht.gnomad_AF_raw >= 0.0002) & (ht.gnomad_AF_raw < 0.0005), "binned_0.0005")
                             .when((ht.gnomad_AF_raw >= 0.0005) & (ht.gnomad_AF_raw < 0.001), "binned_0.001")
                             .when((ht.gnomad_AF_raw >= 0.001) & (ht.gnomad_AF_raw < 0.002), "binned_0.002")
                             .when((ht.gnomad_AF_raw >= 0.002) & (ht.gnomad_AF_raw < 0.005), "binned_0.005")
                             .when((ht.gnomad_AF_raw >= 0.005) & (ht.gnomad_AF_raw < 0.01), "binned_0.01")
                             .when((ht.gnomad_AF_raw >= 0.01) & (ht.gnomad_AF_raw < 0.02), "binned_0.02")
                             .when((ht.gnomad_AF_raw >= 0.02) & (ht.gnomad_AF_raw < 0.05), "binned_0.05")
                             .when((ht.gnomad_AF_raw >= 0.05) & (ht.gnomad_AF_raw < 0.1), "binned_0.1")
                             .when((ht.gnomad_AF_raw >= 0.1) & (ht.gnomad_AF_raw < 0.2), "binned_0.2")
                             .when((ht.gnomad_AF_raw >= 0.2) & (ht.gnomad_AF_raw < 0.5), "binned_0.5")
                             .when((ht.gnomad_AF_raw >= 0.5) & (ht.gnomad_AF_raw <= 1), "binned_1")
                             .default(hl.null(hl.tstr))))
    return ht

def main(args):
    hl.init(log='/variant_histograms.log')
    data_type = 'genomes' if args.genomes else 'exomes'

    metrics = ['FS', 'InbreedingCoeff', 'MQ', 'MQRankSum', 'QD', 'ReadPosRankSum', 'SOR', 'BaseQRankSum',
               'ClippingRankSum', 'DP', 'VQSLOD', 'rf_tp_probability', 'pab_max']

    ht = hl.read_table(release_ht_path(data_type, nested=False))
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
        hist_out = hl.array([hists[f'{metric}'].annotate(metric=metric) for metric in metrics])

        # Aggregate QUAL stats by bin:
        ht = aggregate_qual_stats_by_bin(ht)
        ht = ht.group_by('metric').aggregate(hist=hl.agg.hist(hl.log(ht.qual, base=10), 1, 10, 36))
        hists = ht.collect()
        hist_out = hist_out.extend(hl.array([x.hist.annotate(metric=x.metric) for x in hists]))

        with hl.hadoop_open(release_var_hist_path(data_type), 'w') as f:
            f.write(hl.eval(hl.json(hist_out)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Input MT is exomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--genomes', help='Input MT is genomes. One of --exomes or --genomes is required.', action='store_true')
    parser.add_argument('--first_pass', help='Determine min/max values for each variant metric (to be used in hand-tooled histogram ranges', action='store_true')
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
