# NOTE:
# This script is needed for the v3.1.2 release to correct a problem discovered in the gnomAD v3.1 release.
# The fix we implemented to correct a technical artifact creating the depletion of homozygous alternate genotypes
# was not performed correctly in situations where samples that were individually heterozygous non-reference had both
# a heterozygous and homozygous alternate call at the same locus.

import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import densify_sites
from gnomad.utils.vcf import SPARSE_ENTRIES

from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.annotations import last_END_position
from gnomad_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("get_impacted_variants")
logger.setLevel(logging.INFO)


def main(args):
    """
    Script used to get variants impacted by homalt hotfix.
    """

    hl.init(log="/get_impacted_variants.log", default_reference="GRCh38")

    logger.info("Reading sparse MT and metadata table with release only samples...")
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, samples_meta=True, release_only=True,
    ).select_entries(*SPARSE_ENTRIES)

    if args.test:
        logger.info("Filtering to two partitions...")
        mt = mt._filter_partitions(range(2))

    logger.info(
        "Adding annotation for whether original genotype (pre-splitting multiallelics) is het nonref..."
    )
    # Adding a Boolean for whether a sample had a heterozygous non-reference genotype
    # Need to add this prior to splitting MT to make sure these genotypes
    # are not adjusted by the homalt hotfix downstream
    mt = mt.annotate_entries(_het_non_ref=mt.LGT.is_het_non_ref())

    logger.info("Splitting multiallelics...")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Reading in v3.0 release HT (to get frequency information)...")
    freq_ht = (
        release_sites(public=True).versions["3.0"].ht().select_globals().select("freq")
    )
    mt = mt.annotate_rows(AF=freq_ht[mt.row_key].freq[0].AF)

    logger.info(
        "Filtering to common (AF > 0.01) variants with at least one het nonref call..."
    )
    sites_ht = mt.filter_rows(
        (mt.AF > 0.01) & hl.agg.any(mt._het_non_ref & ((mt.AD[1] / mt.DP) > 0.9))
    ).rows()

    logger.info("Densifying MT to het nonref sites only...")
    # NOTE: densify_sites operates on locus only (doesn't check alleles)
    # https://github.com/broadinstitute/gnomad_methods/blob/master/gnomad/utils/sparse_mt.py#L102
    # Thus, it doesn't matter that `sites_ht` has already been split
    # NOTE: set semi_join_rows to False here because the sites_HT is small
    mt = densify_sites(
        mt, sites_ht, last_END_position.versions["3.1"].ht(), semi_join_rows=False,
    )
    mt = mt.filter_rows(
        (hl.len(mt.alleles) > 1) & (hl.is_defined(sites_ht[mt.row_key]))
    )

    logger.info(
        "Writing out dense MT with only het nonref calls that may have been incorrectly adjusted with the homalt hotfix..."
    )
    if args.test:
        mt.write(
            "gs://gnomad-tmp/release_3.1.2/het_nonref_fix_sites_3.1.2_test.mt",
            overwrite=True,
        )
    else:
        mt.write(
            "gs://gnomad-tmp/release_3.1.2/het_nonref_fix_sites.mt", overwrite=True
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--test", help="Runs a test on two partitions.", action="store_true"
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
