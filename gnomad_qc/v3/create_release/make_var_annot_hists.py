import argparse
import json
import logging

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import (create_frequency_bins_expr,
                                      get_annotations_hists)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.constants import CURRENT_RELEASE
from gnomad_qc.v3.resources.release import (annotation_hists_path,
                                            qual_hists_json_path,
                                            release_ht_path)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


LOG10_ANNOTATIONS = ["AS_VarDP", "QUALapprox", "AS_QUALapprox"]
"""
List of annotations to log scale when creating histograms. 
"""


def create_frequency_bins_expr_inbreeding(
        AF: hl.expr.NumericExpression
) -> hl.expr.StringExpression:
    """
    Creates bins for frequencies in preparation for aggregating QUAL by frequency bin.

    Uses bins of < 0.0005 and >= 0.0005

    NOTE: Frequencies should be frequencies from raw data.
    Used when creating site quality distribution json files.

    :param AC: Field in input that contains the allele count information
    :param AF: Field in input that contains the allele frequency information
    :return: Expression containing bin name
    :rtype: hl.expr.StringExpression
    """
    bin_expr = (
        hl.case()
            .when(AF < 0.0005, "under_0.0005")
            .when((AF >= 0.0005) & (AF <= 1), "over_0.0005")
            .default(hl.null(hl.tstr))
    )
    return bin_expr


def main(args):
    hl.init(default_reference='GRCh38', log='/variant_histograms.log')

    logger.info("Loading ANNOTATIONS_HISTS dictionary...")
    if not file_exists(annotation_hists_path('3.1')):
        raise DataException(
            "Annotation hists JSON file not found. Need to create this JSON before running script!"
        )

    with hl.hadoop_open(annotation_hists_path('3.1')) as a:
        ANNOTATIONS_HISTS = json.loads(a.read())

    # NOTE: histogram aggregations on these metrics are done on the entire callset (not just PASS variants), on raw data
    metrics = list(ANNOTATIONS_HISTS.keys())


    ht = hl.read_table('gs://gnomad/release/v3.1/ht/genomes/gnomad.genomes.v3.1.sites.ht')  # TODO: update path to use resources

    # NOTE: histogram aggregations are done on the entire callset (not just PASS variants), on raw data
    ht = ht.select(freq=ht.freq, info=ht.info.select(*metrics))

    logger.info("Getting info annotation histograms...")
    hist_ranges_expr = get_annotations_hists(
        ht, ANNOTATIONS_HISTS, LOG10_ANNOTATIONS
    )

    # NOTE: run the following code in a first pass to determine bounds for metrics
    # Evaluate minimum and maximum values for each metric of interest
    # This doesn't need to be run unless the defaults do not result in nice-looking histograms.
    if args.first_pass:
        logger.info(
            "Evaluating minimum and maximum values for each metric of interest, and capping maximum value at 1e10"
        )
        minmax_dict = {}
        for metric in metrics:
            minmax_dict[metric] = hl.struct(
                min=hl.agg.min(ht.info[metric]),
                max=hl.if_else(
                    hl.agg.max(ht.info[metric]) < 1e10,
                    hl.agg.max(ht.info[metric]),
                    1e10,
                ),
            )
        minmax = ht.aggregate(hl.struct(**minmax_dict))
        logger.info(f"Metrics bounds: {minmax}")
    else:
        logger.info(
            "Aggregating hists over ranges determined using first pass run..."
        )
        hists = ht.aggregate(
            hl.array(
                [
                    hist_expr.annotate(metric=hist_metric)
                    for hist_metric, hist_expr in hist_ranges_expr.items()
                ]
            ).extend(
                hl.array(
                    hl.agg.group_by(
                        create_frequency_bins_expr(
                            AC=ht.freq[1].AC, AF=ht.freq[1].AF
                        ),
                        # Decided to use QUALapprox because its formula is easier to interpret than QUAL's
                        hl.agg.hist(
                            hl.log10(ht.info.QUALapprox),
                            *ANNOTATIONS_HISTS["QUALapprox"],
                        ),
                    )
                ).map(lambda x: x[1].annotate(metric="QUALapprox-"+x[0]))
            ).extend(
                hl.array(
                    hl.agg.group_by(
                        create_frequency_bins_expr(
                            AC=ht.freq[1].AC, AF=ht.freq[1].AF
                        ),
                        hl.agg.hist(
                            hl.log10(ht.info.AS_QUALapprox),
                            *ANNOTATIONS_HISTS["AS_QUALapprox"],
                        ),
                    )
                ).map(lambda x: x[1].annotate(metric="AS_QUALapprox-"+x[0]))
            ),
            _localize=False,
        )

        testing = {'under_0.0005': [-0.0005, 0.0005, 50], 'over_0.0005': [-0.25, 1, 50]}
        ht = ht.annotate(_af_bin=create_frequency_bins_expr_inbreeding(AF=ht.freq[1].AF))
        inbreeding_hists = [ht.aggregate(hl.agg.filter(ht._af_bin == x,
                                    hl.agg.hist(
                                        ht.info.InbreedingCoeff,
                                        *testing[x],
                                    ))
                      ).annotate(metric="InbreedingCoeff" + "-" + x) for x in testing
                            ]

        hists = hl.eval(hl.json(hists))
        inbreeding_hists = hl.eval(hl.json(inbreeding_hists))
        hists = hists[:-1] + "," + inbreeding_hists[1:]
        print(hists)
        logger.info("Writing output")
        with hl.hadoop_open(qual_hists_json_path(), "w") as f:
            f.write(hists)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--first_pass', help='Determine min/max values for each variant metric and prints them to stdout (to be used in hand-tooled histogram ranges). Note that this should only be run if the defaults do not result in well-behaved histograms.', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
