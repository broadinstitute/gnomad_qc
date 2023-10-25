# noqa: D100

import argparse
import json
import logging

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import create_frequency_bins_expr, get_annotations_hists
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v4.resources.release import (
    annotation_hists_params_path,
    qual_hists_json_path,
    release_ht_path,
)

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
    AF: hl.expr.NumericExpression,
) -> hl.expr.StringExpression:
    """
    Create bins for frequencies in preparation for aggregating Inbreedind Coefficient by frequency bin.

    Uses bins of < 0.0005 and >= 0.0005

    NOTE: Frequencies should be frequencies from raw data.
    Used when creating site quality distribution json files.

    :param AF: Field in input that contains the allele frequency information
    :return: Expression containing bin name
    :rtype: hl.expr.StringExpression
    """
    bin_expr = (
        hl.case()
        .when(AF < 0.0005, "under_0.0005")
        .when(AF >= 0.0005, "over_0.0005")
        .or_missing()
    )
    return bin_expr


def main(args):  # noqa: D103
    hl.init(
        default_reference="GRCh38",
        log="/variant_histograms.log",
        tmp_dir="gs://gnomad-tmp-4day/variant_histograms",
    )
    data_type = args.data_type
    test = args.test

    logger.info(
        "Loading ANNOTATIONS_HISTS dictionary with metric histogram parameters..."
    )
    if not file_exists(annotation_hists_params_path(data_type=data_type)):
        raise DataException(
            "Annotation hists JSON file not found. Need to create this JSON before"
            "running script. Use previous release's JSON as a starting point for "
            "determining metrics histogram bounds."
        )

    with hl.hadoop_open(annotation_hists_params_path(data_type=data_type)) as a:
        ANNOTATIONS_HISTS = json.loads(a.read())

    # NOTE: histogram aggregations on these metrics are done on the entire
    # callset (not just PASS variants), on raw data # noqa
    logger.info("Loading release Hail Table...")
    ht = hl.read_table(release_ht_path(data_type=data_type, public=False))

    if test:
        # Keep only PCSK9.
        logger.info("Filtering to PCSK9 region...")
        ht = hl.filter_intervals(
            ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )

    ht = ht.select(freq=ht.freq, info=ht.info.select(*ANNOTATIONS_HISTS))

    inbreeding_bin_ranges = ANNOTATIONS_HISTS["inbreeding_coeff"]

    # Remove InbreedingCoeff from ANNOTATIONS_HISTS. It requires different ranges
    # by allele frequency and needs to be handled differently. It is stored as a
    # dictionary in annotation_hists_path
    ANNOTATIONS_HISTS = {
        k: v for k, v in ANNOTATIONS_HISTS.items() if k != "inbreeding_coeff"
    }

    logger.info("Getting info annotation histograms...")
    hist_ranges_expr = get_annotations_hists(ht, ANNOTATIONS_HISTS, LOG10_ANNOTATIONS)

    # Note: This evaluates the minimum and maximum values for each metric of interest
    # to help determine the bounds of the hists. Run this first, then update values in
    # annotation_hists_path JSON if necessary after examining data in a notebook
    if args.determine_bounds:
        logger.info(
            "Evaluating minimum and maximum values for each metric of interest. Maximum"
            " values capped at 1e10."
        )
        minmax_dict = {}
        for metric in ANNOTATIONS_HISTS:
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
            "Aggregating hists over ranges defined in the annotation_hists_path JSON"
            " file. --determine_bounds can be used to help define these ranges if"
            " initial results seem off..."
        )
        # NOTE: QUALApprox and AS_QUALApprox are log10 scaled and also split into
        # frequency bins
        hists = ht.aggregate(
            hl.array(
                [
                    hist_expr.annotate(metric=hist_metric)
                    for hist_metric, hist_expr in hist_ranges_expr.items()
                ]
            )
            .extend(
                hl.array(
                    hl.agg.group_by(
                        create_frequency_bins_expr(AC=ht.freq[1].AC, AF=ht.freq[1].AF),
                        hl.agg.hist(
                            hl.log10(ht.info.QUALapprox),
                            *ANNOTATIONS_HISTS["QUALapprox"],
                        ),
                    )
                ).map(lambda x: x[1].annotate(metric="QUALapprox-" + x[0]))
            )
            .extend(
                hl.array(
                    hl.agg.group_by(
                        create_frequency_bins_expr(AC=ht.freq[1].AC, AF=ht.freq[1].AF),
                        hl.agg.hist(
                            hl.log10(ht.info.AS_QUALapprox),
                            *ANNOTATIONS_HISTS["AS_QUALapprox"],
                        ),
                    )
                ).map(lambda x: x[1].annotate(metric="AS_QUALapprox-" + x[0]))
            ),
            _localize=False,
        )

        # # Defining hist range for inbreeding coefficient and bins for allele frequency
        # groups because they needed different ranges
        ht = ht.annotate(af_bin=create_frequency_bins_expr_inbreeding(AF=ht.freq[1].AF))
        inbreeding_hists = [
            ht.aggregate(
                hl.agg.filter(
                    ht.af_bin == x,
                    hl.agg.hist(
                        ht.info.inbreeding_coeff,
                        *inbreeding_bin_ranges[x],
                    ),
                )
            ).annotate(metric="inbreeding_coeff" + "-" + x)
            for x in inbreeding_bin_ranges
        ]

        hists = hl.eval(hl.json(hists))
        inbreeding_hists = hl.eval(hl.json(inbreeding_hists))

        # Note: The following removes '}' from the JSON stored in hists and '{' from the JSON stored in
        # inbreeding_hists then joins them together to be written out as a single JSON
        hists = hists[:-1] + "," + inbreeding_hists[1:]

        output_path = qual_hists_json_path(data_type=data_type, test=test)
        logger.info("Writing output to %s ...", output_path)
        with hl.hadoop_open(output_path, "w") as f:
            f.write(hists)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--determine_bounds",
        help=(
            "Determine min/max values for each variant metric and prints them to stdout"
            " (to be used in hand-tooled histogram ranges). Note that this should only"
            " be run if the defaults do not result in well-behaved histograms."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--data-type", help="Data type", default="exomes", choices=["exomes", "genomes"]
    )
    parser.add_argument(
        "--test",
        help="Use test resources and test on PCSK9 region",
        action="store_true",
    )
    args = parser.parse_args()
    main(args)
