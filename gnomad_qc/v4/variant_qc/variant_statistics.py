"""Produce per-exome stats, per-genome stats, and variant type counts for v4."""
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
    # release_lof,
    release_sites,
    # release_summary_stats
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("summary_stats")
logger.setLevel(logging.INFO)


def per_sample_datatype(
    data_type: str = "exomes",
    test: bool = False,
    overwrite: bool = False,
) -> float:
    """
    Return a float of the mean number of variants called per sample for a chosen data type .

    :param data_type: String of either "exomes" or "genomes" for the sample type.
    :param test: Boolean for if you would like to use a small test set.
    :param overwrite: Boolean to overwrite checkpoint files if requested.
    """

    # Read in release data and metadata and VDS
    ht_release_data = release_sites(data_type).ht()
    if test:
        ht_release_data = ht_release_data.filter(
            ht_release_data.locus.contig == "chr22"
        )
    meta_ht = meta(data_type=data_type).ht()
    vds_data = get_gnomad_v4_vds(test=test, release_only=True)
    vds_data_mt = vds_data.variant_data

    # Filter to only variants and samples in release data
    vds_data_filtered = vds_data_mt.filter_rows(
        ~hl.is_missing(ht_release_data[vds_data_mt.row_key])
    )

    vds_data_filtered = vds_data_filtered.filter_cols(
        meta_ht[vds_data_filtered.s].release
    )

    # Select data down to: only variant locus & alleles and only if samples do or don't have it
    # Create then select down to GT
    vds_data_filtered = vds_data_filtered.annotate_entries(
        GT=hl.vds.lgt_to_gt(vds_data_filtered.LGT, vds_data_filtered.LA)
    )
    vds_data_filtered = vds_data_filtered.select_entries("GT")
    vds_data_filtered = vds_data_filtered.select_rows()
    vds_data_filtered = vds_data_filtered.select_cols()

    # Perform Hail's Sample QC module
    sample_qc_ht = hl.sample_qc(vds_data_filtered).cols().key_by("s")

    sample_qc_ht = sample_qc_ht.checkpoint(
        f"gs://gnomad-tmp-4day/sample_qc_mt_per_{data_type}_sample.ht",
        overwrite=overwrite,
    )

    # Column 'n_called' is a per-sample metric of the number of variants called
    # Is what we want to report and return
    mean_called_per = sample_qc_ht.aggregate(
        hl.agg.mean(sample_qc_ht.sample_qc.n_called)
    )

    logger.info(f"mean called per exome: {mean_called_per}")

    return mean_called_per


def variant_types():
    # what is needed: stats from your notebook for: 
    # Indels and Snps per data type
    # How many (number and proportion) pass VQSR and pass All Filters
    pass

def vep_consequences():
    # what is needed: stats from notebook for: 
    # VEP consequence breakdown , in those two different returns
    # let users specify: pass all filters, pass VQSR, or all variants
    pass


def main(args):
    """Collect summary stats of gnomAD v4 data."""
    hl.init(
        log="/per_exome.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    data_type = args.data_type
    test = args.test
    overwrite = args.overwrite

    if args.get_per_sample:
        per_sample_datatype(data_type, test, overwrite)


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
        "--get-per-sample",
        help="Get variants per sample data type.",
        action="store_true",
    )
    parser.add_argument(
        "--data-type",
        default="exomes",
        choices=["exomes", "genomes"],
        help="Data type (exomes or genomes) to produce summary stats for.",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
