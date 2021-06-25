import argparse
import logging

import hail as hl

from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SPARSE_ENTRIES

from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.annotations import get_info
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
    mt = mt.checkpoint("gs://gnomad-tmp/gnomad_temp_het_nonref_sites.mt")

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
    mt.write("gs://gnomad/release/3.1.2/mt/het_nonref_fix_sites.mt", overwrite=True)


def main(args):
    """
    Script used to get variants impacted by homalt hotfix.
    Most of the work was done in notebooks, and any results created in notebooks are noted with comments.
    """
    data_source = "broad"
    freeze = args.freeze
    tranche_data = (data_source, freeze)

    hl.init(log="/get_impacted_variants.log", default_reference="GRCh38")

    logger.info("Reading in raw MT...")
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True,
        release_only=not args.include_non_release,
        samples_meta=True,
    ).select_entries(
        *SPARSE_ENTRIES
    )

    logger.info(
        "Adding annotation for whether original genotype (pre-splitting multiallelics) is het nonref..."
    )
    mt = mt.annotate_entries(het_non_ref=mt.LGT.is_het_non_ref())

    logger.info("Splitting multiallelics...")
    mt = hl.experimental.sparse_split_multi(mt)

    logger.info("Reading in v3.0 release HT (to get frequency information)...")
    freq_ht = release_sites().versions['3'].ht().select_globals().select("freq")

    if args.get_het_non_ref_impacted_var:
        logger.info("Reading in info HT...")
        info_ht = get_info().ht()

        logger.info(
            "Checking for variants with het nonref calls that were incorrectly adjusted with homalt hotfix..."
        )
        get_het_non_ref_impacted_var(mt, info_ht, freq_ht)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--get-het-non-ref-impacted-var",
        help="Filter raw MT to variants with het nonref genotypes impacted by homalt hotfix",
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
