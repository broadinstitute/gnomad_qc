from gnomad_hail import try_slack
from gnomad_hail.utils.annotations import get_annotations_hists, create_frequency_bins_expr
from gnomad_hail.utils.constants import ANNOTATIONS_HISTS
import hail as hl
import argparse
from gnomad_qc.v3.resources import qual_hists_json_path, release_ht_path


def main(args):
    hl.init(default_reference='GRCh38', log='/variant_histograms.log')

    ht = hl.read_table(release_ht_path())
    # NOTE: histogram aggregations are done on the entire callset (not just PASS variants), on raw data

    hist_dict = ANNOTATIONS_HISTS
    hist_dict['MQ'] = (20, 60, 40) # Boundaries changed for v3, but could be a good idea to settle on a standard
    hist_ranges_expr = get_annotations_hists(
        ht,
        ANNOTATIONS_HISTS
    )

    # NOTE: run the following code in a first pass to determine bounds for metrics
    # Evaluate minimum and maximum values for each metric of interest
    if args.first_pass:
        minmax_dict = {}
        for metric in hist_ranges_expr.keys():
            minmax_dict[metric] = hl.struct(min=hl.agg.min(ht[metric]), max=hl.if_else(hl.agg.max(ht[metric])<1e10, hl.agg.max(ht[metric]), 1e10))
        minmax = ht.aggregate(hl.struct(**minmax_dict))
        print(minmax)
    else:
        # Aggregate hists over hand-tooled ranges
        hists = ht.aggregate(
            hl.array(
                [hist_expr.annotate(metric=hist_metric) for hist_metric, hist_expr in hist_ranges_expr.items()]
            ).extend(
                hl.array(
                    hl.agg.group_by(
                        create_frequency_bins_expr(
                            AC=ht.freq[1].AC,
                            AF=ht.freq[1].AF
                        ),
                        hl.agg.hist(hl.log10(ht.info.QUALapprox), 1, 10, 36)
                    )
                ).map(
                    lambda x: x[1].annotate(metric=x[0])
                )
            ),
            _localize=False
        )

        with hl.hadoop_open(qual_hists_json_path, 'w') as f:
            f.write(hl.eval(hl.json(hists)))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_pass', help='Determine min/max values for each variant metric and prints them to stdout (to be used in hand-tooled histogram ranges)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
