"""Produce per-exome stats for v4 exomes."""
import argparse
import logging

import hail as hl
from gnomad.assessment.summary_stats import (
    default_generate_gene_lof_matrix,
    default_generate_gene_lof_summary,
    get_summary_counts,
)
from gnomad.utils.filtering import filter_to_adj
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds, get_logging_path
import gnomad_qc.v4.resources.meta as meta
from gnomad_qc.v4.resources.release import (
    release_lof,
    release_sites,
    release_summary_stats,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("summary_stats")
logger.setLevel(logging.INFO)


def main(args):
    """Collect summary stats of gnomAD v4 data."""
    hl.init(
        log="/per_exome.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    test = args.test
    overwrite = args.overwrite
    data_type = args.data_type

    ht_release_exomes = release_sites(data_type).ht()
    vds_exomes = get_gnomad_v4_vds(test=test,release_only=True)
    vds_exomes_vd = vds_exomes.variant_data
    vds_exomes_ref = vds_exomes.reference_data

    vds_exomes_vd_filtered = vds_exomes_vd.filter_rows(~hl.is_missing(ht_release_exomes[vds_exomes_vd.row_key])) # woo this method is right!

    meta_ht = meta.meta(data_type=data_type).ht()
 
    vds_exomes_vd_filtered_s = vds_exomes_vd_filtered.filter_cols(meta_ht[vds_exomes_vd_filtered.s].release)

    vds_exomes_vd_filtered_s = vds_exomes_vd_filtered_s.annotate_entries(GT = vds_exomes_vd_filtered_s.LGT)

    vds_exomes_vd_filtered_s_rowselect = vds_exomes_vd_filtered_s.select_rows()
    vds_exomes_vd_filtered_s_rowselect_colselect = vds_exomes_vd_filtered_s_rowselect.select_cols()

    vds_exomes_vd_filtered_s_rowselect_colselect = vds_exomes_vd_filtered_s_rowselect_colselect.select_entries('GT')

    sample_qc_mt = hl.sample_qc(vds_exomes_vd_filtered_s_rowselect_colselect)

    sample_qc_mt = sample_qc_mt.checkpoint('gs://gnomad-marten/sample_qc_mt.mt',overwrite=overwrite)

    mean_called_per_exome = sample_qc_mt.aggregate_entries(hl.agg.mean(sample_qc_mt.sample_qc.n_called))

    print('mean called per exome: ',mean_called_per_exome)

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
        "--data-type",
        default="exomes",
        choices=["exomes","genomes"],
        help="Data type (exomes or genomes) to produce summary stats for."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
