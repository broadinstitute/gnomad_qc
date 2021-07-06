import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import densify_sites
from gnomad.utils.vcf import SPARSE_ENTRIES


from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.annotations import get_info, last_END_position
from gnomad_qc.slack_creds import slack_token


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("get_impacted_variants")
logger.setLevel(logging.INFO)


def get_het_non_ref_impacted_var(
    mt: hl.MatrixTable, info_ht: hl.Table, freq_ht: hl.Table
) -> None:
    """
    Filter to variants where homalt hotfix incorrectly adjusts het nonref genotype calls.

    :param hl.MatrixTable mt: Raw, split MatrixTable annotated with original genotype het nonref status.
    :param hl.Table info_ht: Info Hail Table containing AS_lowqual information.
    :param hl.Table freq_ht: Hail Table containing 455K frequency information.
    :return: None
    """

    logger.info("Removing low qual variants and star alleles...")
    mt = mt.filter_rows(
        (~info_ht[mt.row_key].AS_lowqual)
        & ((hl.len(mt.alleles) > 1) & (mt.alleles[1] != "*"))
    )

    logger.info("Filtering to common (AF > 0.01) variants...")
    mt = mt.annotate_rows(AF=freq_ht[mt.row_key].freq[0].AF)
    mt = mt.filter_rows(mt.AF > 0.01)

    logger.info(
        "Filtering to variants with at least one het nonref call and checkpointing..."
    )
    mt = mt.filter_rows(hl.agg.sum(mt.het_non_ref) >= 1)

    logger.info(
        "Adding new genotype entries with original homalt hotfix and het nonref-aware hotfix..."
    )
    mt = mt.annotate_entries(
        GT_hotfix=hl.if_else(
            mt.GT.is_het() & (mt.AF > 0.01) & (mt.AD[1] / mt.DP > 0.9),
            hl.call(1, 1),
            mt.GT,
        ),
        GT_hetnonref_fix=hl.if_else(
            mt.GT.is_het()
            & ~mt.het_non_ref
            & (mt.AF > 0.01)
            & (mt.AD[1] / mt.DP > 0.9),
            hl.call(1, 1),
            mt.GT,
        ),
    )
    logger.info("Checking homozygote counts using the different genotypes...")
    # Taking index 1 of gt_stats because it contains an entry for the reference allele at index 0
    # homozygote_count (tarray of tint32) -
    # Homozygote genotype counts for each allele, including the reference. Only diploid genotype calls are counted.
    # (from hail docs, https://hail.is/docs/0.2/aggregators.html#hail.expr.aggregators.call_stats)

    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT_hotfix, mt.alleles))
    mt = mt.transmute_rows(hom_count_hotfix=mt.gt_stats.homozygote_count[1])

    mt = mt.annotate_rows(gt_stats=hl.agg.call_stats(mt.GT_hetnonref_fix, mt.alleles))
    mt = mt.transmute_rows(hom_count_hetnonref_fix=mt.gt_stats.homozygote_count[1])

    logger.info(
        "Filtering to rows where homalt hotfix erroneously adjusts het nonref calls..."
    )
    mt = mt.filter_rows(mt.hom_count_hotfix != mt.hom_count_hetnonref_fix)

    return mt.rows()


def main(args):
    """
    Script used to get variants impacted by homalt hotfix.
    Most of the work was done in notebooks, and any results created in notebooks are noted with comments.
    """

    hl.init(log="/get_impacted_variants.log", default_reference="GRCh38")

    logger.info("Reading sparse MT and metadata table with release only samples...")
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, samples_meta=True, release_only=True,
    ).select_entries(*SPARSE_ENTRIES)

    if args.test:
        logger.info("Filtering to two partitions on chr20")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20:1-1000000")])
        mt = mt._filter_partitions(range(2))

    logger.info(
        "Adding annotation for whether original genotype (pre-splitting multiallelics) is het nonref..."
    )
    # Adding a Boolean for whether a sample had a heterozygous non-reference genotype
    # Need to add this prior to splitting MT to make sure these genotypes
    # are not adjusted by the homalt hotfix downstream
    mt = mt.annotate_entries(het_non_ref=mt.LGT.is_het_non_ref())
    # TODO: Need to fix main code to add this annotation and use it also, where to add?

    logger.info("Splitting multiallelics...")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Reading in v3.0 release HT (to get frequency information)...")
    freq_ht = release_sites().versions["3"].ht().select_globals().select("freq")

    logger.info("Reading in info HT...")
    info_ht = get_info().ht()

    logger.info(
        "Checking for variants with het nonref calls that were incorrectly adjusted with homalt hotfix..."
    )
    sites_ht = get_het_non_ref_impacted_var(mt, info_ht, freq_ht)

    logger.info("Densifying MT to het non ref sites only...")
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
        mt = mt.checkpoint(
            "gs://gnomad-tmp/release_3.1.2/het_nonref_fix_sites_3.1.2_test.mt",
            overwrite=True,
        )
    else:
        mt = mt.checkpoint(
            "gs://gnomad-tmp/release_3.1.2/het_nonref_fix_sites.mt", overwrite=True
        )

    logger.info(
        "Writing out variants with het nonref calls that may have been incorrectly adjusted with the homalt hotfix..."
    )
    if args.test:
        mt.rows().write(
            "gs://gnomad-tmp/release_3.1.2/het_nonref_fix_sites_3.1.2_test.ht",
            overwrite=True,
        )
    else:
        mt.rows().write(
            "gs://gnomad/release/3.1.2/ht/het_nonref_fix_sites.ht", overwrite=True
        )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--test", help="Runs a test on two partitions of chr20.", action="store_true"
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
