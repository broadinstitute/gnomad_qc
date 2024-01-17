"""Produce per-exome stats, per-genome stats, and variant type counts for v4."""
import argparse
import logging

import hail as hl
from gnomad.assessment.summary_stats import (
    default_generate_gene_lof_matrix,
    default_generate_gene_lof_summary,
    get_summary_counts,
)
from gnomad.utils import vep
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
) -> None:
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
    vds_data_filtered = vds_data_filtered.annotate_entries(GT=hl.vds.lgt_to_gt(vds_data_filtered.LGT, vds_data_filtered.LA))
    vds_data_filtered = vds_data_filtered.select_entries("GT")
    vds_data_filtered = vds_data_filtered.select_rows()
    vds_data_filtered = (
        vds_data_filtered.select_cols()
    )

    # Perform Hail's Sample QC module
    sample_qc_ht = hl.sample_qc(vds_data_filtered).cols().key_by("s")

    sample_qc_ht = sample_qc_ht.checkpoint(
        f"gs://gnomad-tmp-4day/sample_qc_mt_per_{data_type}_sample.ht", overwrite=overwrite
    )

    mean_called_per = sample_qc_ht.aggregate(
        hl.agg.mean(sample_qc_ht.sample_qc.n_called)
    )

    print("mean called per exome: ", mean_called_per)


def count_variant_types(
    data_type: str = "exomes",
    exclusive: bool = False,
    do_snp: bool = True,
    check_snp_filters: bool = True,
    do_vep_lof: bool = False,
    synonymous: bool = False,
    filter_empty_csq: bool = True,
    canonical: bool = False,
    mane_select: bool = True,
    check_is_snp: bool = False,
    check_iqc: bool = False,
    check_calling: bool = False,
    check_filters: bool = False,
) -> dict:
    """
    Function to return variant count information, either by snp/indel or vep/lof consequence terms.

    :param data_type: String of either "exomes" or "genomes" for the sample type.
    :param exclusive: Boolean to anti-join chosen data type with opposite, to have a HT without variants found in the other
    :param do_snp: Boolean to calculate amount of snp and indel variants.
    :param check_snp_filters: Boolean to count snp and indel variants by all 'filters' fields, instead of the default of passing all filters or not
    :param do_vep_lof: Boolean to generate vep consequence and lof term counts 
    """
    # Read in release sites table to calculate numbers for
    ht_release = release_sites(data_type).ht()
    if exclusive:
        # Do check if what we're looking for already exists, should add
        opposite = list(set(["exomes", "genomes"]) - set([data_type]))[0]
        ht_opposite = release_sites(opposite).ht()
        ht_release = ht_opposite.anti_join(ht_opposite)
    # ht_release = ht_release.select() there's much to maybe select from, so I'll hold off

    return_tables = {}

    # Generate snp and indel counts
    if do_snp:
        snp_args = []
        if check_snp_filters:
            snp_args.append("filters")

        snp_aggregate_table = ht_release.group_by(
            is_snp=hl.is_snp(ht_release.alleles[0], ht_release.alleles[1]),
            all_pass=hl.len(ht_release.filters) == 0,
            *snp_args,
        ).aggregate(n=hl.agg.count())

        # Print and return snp results
        snp_aggregate_table.show()

        return_tables["snp"] = snp_aggregate_table

    # Generate vep and lof counts
    if do_vep_lof:
        # Generate and appropriately filter vep counts
        ht_vep = vep.filter_vep_transcript_csqs(
            ht_release,
            synonymous=synonymous,
            filter_empty_csq=filter_empty_csq,
            canonical=canonical,
            mane_select=mane_select,
        )

        ht_vep = vep.get_most_severe_consequence_for_summary(ht_vep)

        # For all vep terms
        vep_aggregate_table = ht_vep.group_by(csq=ht_vep.most_severe_csq).aggregate(
            n=hl.agg.count()
        )

        vep_aggregate_table.show(n=100)

        # Filter to only loss-of-function and appropriately return these 
        ht_lof = ht_vep.filter(~hl.is_missing(ht_vep.lof))

        lof_group_by_list = []
        for argument in [check_is_snp, check_iqc, check_calling, check_filters]:
            if argument is True:
                lof_group_by_list.append(argument)

        lof_aggregate_table = ht_lof.group_by(
            *lof_group_by_list,
            all_pass=hl.len(ht_release.filters) == 0,
            csq=ht_lof.most_severe_csq,
        ).aggregate(n=hl.agg.count())

        lof_aggregate_table.show(n=100)

        # Return all outputs
        return_tables["vep"] = vep_aggregate_table
        return_tables["lof"] = lof_aggregate_table

        return return_tables


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

    returns = []

    if args.get_per_sample:
        r1 = per_sample_datatype(data_type, test, overwrite)
        returns.append(r1)
    if args.get_snp_indel_counts or args.get_vep_counts:
        r2 = count_variant_types(
            data_type, do_snp=args.get_snp_indel_counts, do_vep_lof=args.get_vep_counts
        )
        returns.append(r2)

    
    


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
        "--get-snp-indel-counts",
        help="Get counts for snps and indels.",
        action="store_true",
    )
    parser.add_argument(
        "--get-vep-counts",
        help="Get counts for VEP consequences.",
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
