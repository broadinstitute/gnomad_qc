"""
Script to remove duplicate UKB samples from split VDS for frequency calculations.

Final sample count should be 730947, count with duplicate UKB samples is 730970.
"""
import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.annotations.generate_freq import get_downsampling_ht
from gnomad_qc.v4.resources.basics import gnomad_v4_genotypes, meta, ukb_dups_idx_path
from gnomad_qc.v4.resources.sample_qc import hard_filtered_samples

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("remove_ukb_dups_from_frequency_split_vds")
logger.setLevel(logging.INFO)


def main(args):
    """Remove duplicate UKB samples from split VDS for frequency calculations."""
    hl.init(
        log="/generate_frequency_data.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )

    test = args.test

    # Get count of release samples.
    meta_ht = meta.ht()
    release_meta_ht = meta_ht.filter(meta_ht.release)
    logger.info("Count of samples in release: %d", release_meta_ht.count())

    # Count of samples in split VDS.
    split_vds = hl.vds.read_vds(
        "gs://gnomad/v4.0/annotations/exomes/temp/gnomad.exomes.v4.0.split.vds"
    )
    logger.info(
        "Count of samples in frequency split VDS: %d",
        split_vds.variant_data.count_cols(),
    )

    # Make a HT indicating which samples to remove from the split VDS.
    # Read in full VDS
    vds = gnomad_v4_genotypes.vds()

    ############################
    # THIS IS THE ORIGINAL REMOVAL THAT SHOULD HAVE HAPPENED, WITH ANNOTATION INSTEAD
    # OF REMOVAL.
    ############################

    # Remove 27 fully duplicated IDs (same exact name for 's' in the VDS). See:
    # https://github.com/broadinstitute/ukbb_qc/blob/70c4268ab32e9efa948fe72f3887e1b81d8acb46/ukbb_qc/resources/basics.py#L286.
    # Confirmed that the col_idx of the 27 dup samples in the original UKB MT match
    # the col_idx of the dup UKB samples in the VDS.
    dup_ids = []
    with hl.hadoop_open(ukb_dups_idx_path, "r") as d:
        for line in d:
            line = line.strip().split("\t")
            dup_ids.append(f"{line[0]}_{line[1]}")

    def _annotate_ukb_dup_by_index(
        mt: hl.MatrixTable, dup_ids: hl.expr.ArrayExpression
    ) -> hl.MatrixTable:
        mt = mt.add_col_index()
        mt = mt.annotate_cols(new_s=hl.format("%s_%s", mt.s, mt.col_idx))
        mt = mt.annotate_cols(keep_sample=~dup_ids.contains(mt.new_s))
        return mt

    dup_ids = hl.literal(dup_ids)
    vd = _annotate_ukb_dup_by_index(vds.variant_data, dup_ids)
    rd = _annotate_ukb_dup_by_index(vds.reference_data, dup_ids)
    vds = hl.vds.VariantDataset(rd, vd)

    ############################
    # Remove hard filtered samples. Only adding because this was done prior to writing
    # out the split VDS, and I want to make sure we are using the correct index. I
    # don't think this impacts that, but it's better to add for consistency.
    ############################
    filter_ht = hard_filtered_samples.ht()
    vds = hl.vds.filter_samples(vds, filter_ht, keep=False)

    # Filter VDS to release samples.
    vds = hl.vds.filter_samples(vds, release_meta_ht)

    ############################
    # THIS IS NEW IN ORDER TO GET THE SPLIT VDS COL IDX.
    ############################
    def _annotate_ukb_dup_by_release_index(
        mt: hl.MatrixTable,
    ) -> hl.MatrixTable:
        mt = mt.add_col_index("col_idx2")
        mt = mt.annotate_cols(new_s2=hl.format("%s_%s", mt.s, mt.col_idx2))
        return mt

    vd = _annotate_ukb_dup_by_release_index(vds.variant_data)
    rd = _annotate_ukb_dup_by_release_index(vds.reference_data)
    vds = hl.vds.VariantDataset(rd, vd)

    dup_ht = vds.variant_data.cols()
    dup_ht = dup_ht.checkpoint(
        "gs://gnomad-tmp/ukb_duplicate_samples_to_remove_from_split_vds.ht",
        overwrite=True,
    )
    dup_remove_ht = dup_ht.filter(~dup_ht.keep_sample)

    ############################
    # Read in split VDS and remove UKB samples and write to new location.
    ############################
    split_vds = hl.vds.read_vds(
        "gs://gnomad/v4.0/annotations/exomes/temp/gnomad.exomes.v4.0.split.vds"
    )

    if test:
        split_vds = hl.vds.VariantDataset(
            split_vds.reference_data._filter_partitions(range(2)),
            split_vds.variant_data._filter_partitions(range(2)),
        )

    def _remove_ukb_dup_by_index(
        mt: hl.MatrixTable, dup_ids: hl.expr.ArrayExpression
    ) -> hl.MatrixTable:
        mt = mt.add_col_index()
        mt = mt.annotate_cols(new_s=hl.format("%s_%s", mt.s, mt.col_idx))
        mt = mt.filter_cols(~dup_ids.contains(mt.new_s)).drop("new_s", "col_idx")
        return mt

    dup_ids = hl.set(dup_remove_ht.new_s2.collect())
    vd = _remove_ukb_dup_by_index(split_vds.variant_data, dup_ids)
    rd = _remove_ukb_dup_by_index(split_vds.reference_data, dup_ids)
    split_vds = hl.vds.VariantDataset(rd, vd)

    if test:
        split_vds = split_vds.checkpoint(
            "gs://gnomad-tmp/gnomad.exomes.v4.0.split_multi.test.vds", overwrite=True
        )
    else:
        split_vds = split_vds.checkpoint(
            "gs://gnomad/v4.0/annotations/exomes/temp/gnomad.exomes.v4.0.split_multi.vds",
            overwrite=True,
        )
    logger.info(
        "Count of samples in variant data of new split VDS: %d",
        split_vds.variant_data.count_cols(),
    )
    logger.info(
        "Count of samples in reference data of new split VDS: %d",
        split_vds.reference_data.count_cols(),
    )

    # Fix downsamplings HT.
    ds_ht = get_downsampling_ht(split_vds.variant_data)

    if test:
        ds_ht = ds_ht.checkpoint(
            "gs://gnomad-tmp/gnomad.exomes.v4.0.downsampling.ht", overwrite=True
        )
    else:
        ds_ht = ds_ht.checkpoint(
            "gs://gnomad/v4.0/annotations/exomes/gnomad.exomes.v4.0.downsampling.ht",
            overwrite=True,
        )

    logger.info("Count of samples in updated downsampling Table: %d", ds_ht.count())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Use only 2 partitions of the VDS as input for testing purposes.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
