"""Run summary stats on gnomAD v4 data."""

import argparse
import logging

import hail as hl
from gnomad.assessment.summary_stats import get_summary_counts
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.assessment import release_summary_stats
from gnomad_qc.v4.resources.basics import get_logging_path
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("summary_stats")
logger.setLevel(logging.INFO)


def main(args):
    """Collect summary stats of gnomAD v4 data."""
    hl.init(
        log="/summary_stats.log",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    hl.default_reference("GRCh38")

    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    test = args.test
    overwrite = args.overwrite
    data_type = args.data_type
    public = args.get_public_sites

    try:
        if args.get_summary_counts:
            ht = release_sites(data_type=data_type, public=public).ht()
            filter_name = None
            if test:
                ht = ht._filter_partitions(range(2))

            if args.interval_qc_pass_only:
                logger.info("Removing regions that fail interval QC...")
                ht = ht.filter(~ht.region_flags.fail_interval_qc)
                filter_name = "interval_qc_pass_only"

            elif args.union_calling_interval_only:
                logger.info(
                    "Removing regions that are outside the union of UKB and Broad "
                    "intervals..."
                )
                ht = ht.filter(
                    ~ht.region_flags.outside_ukb_capture_region
                    | ~ht.region_flags.outside_broad_capture_region
                )
                filter_name = "union_calling_interval_only"

            elif args.intersection_calling_interval_only:
                logger.info(
                    "Removing regions that are outside the intersection of UKB and "
                    "Broad intervals..."
                )
                ht = ht.filter(
                    ~(
                        ht.region_flags.outside_ukb_capture_region
                        | ht.region_flags.outside_broad_capture_region
                    )
                )
                filter_name = "intersection_calling_interval_only"

            freq_index = 0
            if args.use_100k_downsampling:
                freq_index = hl.eval(
                    ht.freq_meta.index(
                        {"downsampling": "100000", "gen_anc": "global", "group": "adj"}
                    )
                )
                filter_name = filter_name + ".100k_downsampling"

            logger.info("Getting summary counts per variant category...")
            ht = get_summary_counts(
                ht.select_globals(),
                canonical_only=False,
                mane_select_only=True,
                index=freq_index,
            )
            meta_ht = meta(data_type=data_type).ht()
            meta_ht = meta_ht.filter(meta_ht.release)

            logger.info(f"Number of release samples: {meta_ht.count()}")
            ht = ht.annotate_globals(num_release_samples=meta_ht.count())
            ht.write(
                release_summary_stats(
                    test=test, data_type=data_type, filter_name=filter_name
                ).path,
                overwrite,
            )

        if args.generate_gene_lof_matrix:
            raise NotImplementedError(
                "Gene LoF matrix generation is not yet implemented."
            )

        if args.summarize_gene_lof_matrix:
            raise NotImplementedError("Gene LoF matrix summary is not yet implemented.")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("summary_stats"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite data (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--test",
        help="Test on a small number of variants",
        action="store_true",
    )
    parser.add_argument(
        "--data-type",
        default="exomes",
        choices=["exomes", "genomes"],
        help="Data type (exomes or genomes) to produce summary stats for.",
    )
    parser.add_argument(
        "--get-summary-counts",
        help="Get summary counts per variant category.",
        action="store_true",
    )
    parser.add_argument(
        "--interval-qc-pass-only",
        help="Generate summary stats on interval QC pass regions only.",
        action="store_true",
    )
    parser.add_argument(
        "--union-calling-interval-only",
        help=(
            "Generate summary stats on regions in the union between UKB and Broad "
            "intervals only."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--intersection-calling-interval-only",
        help=(
            "Generate summary stats on regions in the intersection between UKB and "
            "Broad intervals only."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--use-100k-downsampling",
        help="Use the 100K downsampling for summary stats.",
        action="store_true",
    )

    parser.add_argument(
        "--generate-gene-lof-matrix",
        help="Generate gene LoF matrix.",
        action="store_true",
    )
    parser.add_argument(
        "--summarize-gene-lof-matrix",
        help="Creates gene LoF matrix summary Table.",
        action="store_true",
    )
    parser.add_argument(
        "--get-public-sites",
        help="Get sites table from public bucket.",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
