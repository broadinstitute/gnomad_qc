import argparse
import logging
from pprint import pformat

import hail as hl

from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.resources.grch38.reference_data import clinvar_pathogenic
from gnomad.utils.slack import slack_notifications
from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.variant_qc.evaluation import (
    compute_binned_truth_sample_concordance,
    compute_grouped_binned_ht,
    create_truth_sample_ht,
)
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg

from gnomad_qc.v3.resources import (
    fam_stats,
    get_binned_concordance,
    get_callset_truth_data,
    get_checkpoint_path,
    get_filters,
    get_gnomad_v3_mt,
    get_info,
    get_rf_result,
    get_rf_annotated,
    get_score_quantile_bins,
    TRUTH_SAMPLES,
)
from gnomad_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_evaluation")
logger.setLevel(logging.INFO)


def create_quantile_bin_ht(
    model_id: str, n_bins: int, vqsr: bool = False, overwrite: bool = False
) -> None:
    """
    Creates a table with quantile bin annotations added for a RF run and writes it to its correct location in
    annotations.

    :param model_id: Which variant QC model (rf or vqsr model ID) to quantile bin
    :param n_bins: Number of bins to bin the data into
    :param vqsr: Set True is `model_id` refers to a VQSR filtering model
    :param overwrite: Should output files be overwritten if present
    :return: Nothing
    """
    logger.info(f"Annotating {model_id} HT with quantile bins using {n_bins}")
    info_ht = get_info(split=True).ht()
    if vqsr:
        rf_ht = get_rf_annotated().ht()
        ht = get_filters(model_id, split=True).ht()

        ht = ht.filter(
            ~info_ht[ht.key].AS_lowqual
        )
        ht = ht.annotate(
            **rf_ht[ht.key],
            info=info_ht[ht.key].info,
            score=ht.info.AS_VQSLOD,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            AS_culprit=ht.info.AS_culprit,
        )

    else:
        ht = get_rf_result(model_id=model_id).ht()
        ht = ht.annotate(
            info=info_ht[ht.key].info,
            positive_train_site=ht.tp,
            negative_train_site=ht.fp,
            score=ht.rf_probability["TP"],
        )

    ht = ht.filter(ht.ac_raw > 0)
    ht = ht.filter(
        ~info_ht[ht.key].AS_lowqual & ~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus])
    )
    ht_lcr = filter_low_conf_regions(
        ht,
        filter_lcr=True,
        # TODO: Uncomment when we have decoy path
        filter_decoy=False,  # True,
        filter_segdup=True,
    )
    ht = ht.annotate(non_lcr=hl.is_defined(ht_lcr[ht.key]))
    bin_ht = create_binned_ht(ht, n_bins, add_substrat={"non_lcr": ht.non_lcr})
    bin_ht.write(
        get_score_quantile_bins(model_id, aggregated=False).path, overwrite=overwrite
    )


def create_grouped_bin_ht(model_id: str, overwrite: bool = False) -> None:
    """
    Creates binned data from a quantile bin annotated Table grouped by bin_id (rank, bi-allelic, etc.), contig, snv,
    bi_allelic and singleton containing the information needed for evaluation plots.

    :param str model_id: Which variant QC model (rf or vqsr model ID) to group
    :param bool overwrite: Should output files be overwritten if present
    :return: None
    :rtype: None
    """

    ht = get_score_quantile_bins(model_id, aggregated=False).ht()

    # Count variants for ranking
    count_expr = {
        x: hl.agg.filter(
            hl.is_defined(ht[x]),
            hl.agg.counter(
                hl.cond(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv", "indel")
            ),
        )
        for x in ht.row
        if x.endswith("bin")
    }
    bin_variant_counts = ht.aggregate(hl.struct(**count_expr))
    logger.info(f"Found the following variant counts:\n {pformat(bin_variant_counts)}")
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)

    # Load ClinVar pathogenic data
    ht = ht.annotate(clinvar_path=hl.is_defined(clinvar_pathogenic.ht()[ht.key]))
    trio_stats_ht = fam_stats.ht()

    logger.info(f"Creating grouped bin table...")
    grouped_binned_ht = compute_grouped_binned_ht(
        ht, checkpoint_path=get_checkpoint_path(f"grouped_bin_{model_id}"),
    )

    logger.info(f"Aggregating grouped bin table...")
    parent_ht = grouped_binned_ht._parent
    agg_ht = grouped_binned_ht.aggregate(
        n_clinvar_path=hl.agg.count_where(parent_ht.clinvar_path),
        **score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht)
    )

    agg_ht.write(
        get_score_quantile_bins(model_id, aggregated=True).path, overwrite=overwrite,
    )


def main(args):
    hl.init(log="/variant_qc_evaluation.log")

    vqsr = args.model_id.startswith("vqsr_")
    if args.create_quantile_bin_ht:
        create_quantile_bin_ht(args.model_id, args.n_bins, vqsr, args.overwrite)

    if args.run_sanity_checks:
        ht = get_score_quantile_bins(args.model_id, aggregated=False).ht()
        logger.info("Running sanity checks...")
        print(
            ht.aggregate(
                hl.struct(
                    was_split=hl.agg.counter(ht.was_split),
                    has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_bin)),
                    was_singleton=hl.agg.counter(ht.singleton),
                    has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_bin)),
                    was_split_singleton=hl.agg.counter(ht.singleton & ~ht.was_split),
                    has_biallelic_singleton_rank=hl.agg.counter(
                        hl.is_defined(ht.biallelic_singleton_bin)
                    ),
                )
            )
        )

    if args.create_aggregated_bin_ht:
        logger.warning("Use only workers, it typically crashes with preemptibles")
        create_grouped_bin_ht(args.model_id, args.overwrite)

    if args.extract_truth_samples:
        logger.info(f"Extracting truth samples from MT...")
        mt = get_gnomad_v3_mt(
            key_by_locus_and_alleles=True, remove_hard_filtered_samples=False
        )

        mt = mt.filter_cols(
            hl.literal([v["s"] for k, v in TRUTH_SAMPLES.items()]).contains(mt.s)
        )
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

        # Checkpoint to prevent needing to go through the large table a second time
        mt = mt.checkpoint(
            get_checkpoint_path("truth_samples", mt=True), overwrite=args.overwrite,
        )

        for truth_sample in TRUTH_SAMPLES:
            truth_sample_mt = mt.filter_cols(
                mt.s == TRUTH_SAMPLES[truth_sample]["s"]
            )
            # Filter to variants in truth data
            truth_sample_mt = truth_sample_mt.filter_rows(
                hl.agg.any(truth_sample_mt.GT.is_non_ref())
            )
            truth_sample_mt.naive_coalesce(args.n_partitions).write(
                get_callset_truth_data(truth_sample).path, overwrite=args.overwrite,
            )

    if args.merge_with_truth_data:
        for truth_sample in TRUTH_SAMPLES:
            logger.info(
                f"Creating a merged table with callset truth sample and truth data for {truth_sample}..."
            )

            # Load truth data
            mt = get_callset_truth_data(truth_sample).mt()
            truth_hc_intervals = TRUTH_SAMPLES[truth_sample]["hc_intervals"].ht()
            truth_mt = TRUTH_SAMPLES[truth_sample]["truth_mt"].mt()
            truth_mt = truth_mt.key_cols_by(
                s=hl.str(TRUTH_SAMPLES[truth_sample]["s"])
            )

            # remove low quality sites
            info_ht = get_info(split=True).ht()
            mt = mt.filter_rows(
                ~info_ht[mt.row_key].AS_lowqual
            )

            ht = create_truth_sample_ht(mt, truth_mt, truth_hc_intervals)
            ht.write(
                get_callset_truth_data(truth_sample, mt=False).path,
                overwrite=args.overwrite,
            )

    if args.bin_truth_sample_concordance:
        for truth_sample in TRUTH_SAMPLES:
            logger.info(
                f"Creating binned concordance table for {truth_sample} for metric {model_id}"
            )
            ht = get_callset_truth_data(truth_sample, mt=False).ht()

            info_ht = get_info(split=True).ht()
            ht = ht.filter(
                ~info_ht[ht.key].AS_lowqual & ~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus])
            )

            logger.info("Filtering out low confidence regions and segdups...")
            ht = filter_low_conf_regions(
                ht,
                filter_lcr=True,
                # TODO: Uncomment when we have decoy path
                filter_decoy=False,  # True,
                filter_segdup=True,
            )

            logger.info(
                "Loading HT containing RF or VQSR scores annotated with a bin based on metric quantiles..."
            )
            metric_ht = get_score_quantile_bins(args.model_id, aggregated=False).ht()
            ht = ht.filter(hl.is_defined(metric_ht[ht.key]))

            ht = ht.annotate(score=metric_ht[ht.key].score)

            ht = compute_binned_truth_sample_concordance(ht, metric_ht, args.n_bins)
            ht.write(
                get_binned_concordance(args.model_id, truth_sample).path,
                overwrite=args.overwrite,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--model_id",
        help="Model ID.",
        required=False,
    )

    parser.add_argument(
        "--create_quantile_bin_ht",
        help="When set, creates file annotated with quantile bin based on vqsr/RF score.",
        action="store_true",
    )
    parser.add_argument(
        "--run_sanity_checks",
        help="When set, runs ranking sanity checks.",
        action="store_true",
    )
    parser.add_argument(
        "--create_aggregated_bin_ht",
        help="When set, creates a file with aggregate counts of variants based on quantile bins.",
        action="store_true",
    )
    parser.add_argument(
        "--n_bins",
        help="Number of bins for the binned file (default: 100)",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--extract_truth_samples",
        help="Extract truth samples from matrix table",
        action="store_true",
    )
    parser.add_argument(
        "--n_partitions",
        help="Desired number of partitions for output Table/MatrixTable",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--merge_with_truth_data",
        help="Computes a table for each truth sample comparing the truth sample in the callset vs the truth.",
        action="store_true",
    )
    parser.add_argument(
        "--bin_truth_sample_concordance",
        help="Merges individual concordance results with specified metric binned files.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
