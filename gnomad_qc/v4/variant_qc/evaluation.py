# noqa: D100

import argparse
import logging
from pprint import pformat

import hail as hl
from gnomad.resources.grch38.reference_data import clinvar, telomeres_and_centromeres
from gnomad.utils.filtering import filter_low_conf_regions, filter_to_clinvar_pathogenic
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.evaluation import (
    compute_binned_truth_sample_concordance,
    compute_grouped_binned_ht,
    create_truth_sample_ht,
)
from gnomad.variant_qc.pipeline import create_binned_ht, score_bin_agg
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources import (
    allele_data,
    fam_stats,
    get_binned_concordance,
    get_checkpoint_path,
    get_info,
    get_rf_annotations,
    get_rf_result,
    get_score_bins,
    get_vqsr_filters,
)
from gnomad_qc.v3.resources.variant_qc import TRUTH_SAMPLES, get_callset_truth_data
from gnomad_qc.v4.resources.basics import TRUTH_SAMPLES_S, get_gnomad_v4_vds

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_evaluation")
logger.setLevel(logging.INFO)


def create_bin_ht(model_id: str, n_bins: int, hgdp_tgp_subset: bool) -> hl.Table:
    """
    Create a table with bin annotations added for a RF or VQSR run and writes it to its correct location in annotations.

    :param model_id: Which variant QC model (RF or VQSR model ID) to annotate with bin
    :param n_bins: Number of bins to bin the data into
    :param hgdp_tgp_subset: Whether this is for the HGDP + 1KG/TGP subset, if True this will filter the HT using the subset's raw allele count instead of gnomAD's high quality raw allele count
    :return: Table with bin annotations
    """
    logger.info(f"Annotating {model_id} HT with bins using {n_bins} bins")
    info_ht = get_info(split=True).ht()
    if model_id.startswith("vqsr"):
        rf_ht = get_rf_annotations().ht()
        ht = get_vqsr_filters(model_id, split=True).ht()

        ht = ht.annotate(
            **rf_ht[ht.key],
            info=info_ht[ht.key].info,
            score=ht.info.AS_VQSLOD,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            AS_culprit=ht.info.AS_culprit,
        )

        if hgdp_tgp_subset:
            allele_data_ht = allele_data.ht()
            ht = ht.annotate(**allele_data_ht[ht.key].allele_data)
        else:
            # Remove all variants with an undefined ac_raw because ac_raw was calculated on the high quality samples only
            # and VQSR was run before sample filtering
            ht = ht.filter(hl.is_defined(ht.ac_raw))

    else:
        if hgdp_tgp_subset:
            raise ValueError(
                "The hgdp_tgp_subset option is only relevant for VQSR models in all v3"
                " releases."
            )

        ht = get_rf_result(model_id=model_id).ht()
        ht = ht.annotate(
            info=info_ht[ht.key].info,
            positive_train_site=ht.tp,
            negative_train_site=ht.fp,
            score=ht.rf_probability["TP"],
        )

    ht = ht.filter(
        ~info_ht[ht.key].AS_lowqual
        & ~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus])
    )
    ht_non_lcr = filter_low_conf_regions(
        ht,
        filter_lcr=True,
        # TODO: Uncomment when we have decoy path
        filter_decoy=False,  # True,
        filter_segdup=True,
    )
    ht = ht.annotate(non_lcr=hl.is_defined(ht_non_lcr[ht.key]))
    bin_ht = create_binned_ht(ht, n_bins, add_substrat={"non_lcr": ht.non_lcr})

    return bin_ht


def create_aggregated_bin_ht(model_id: str) -> hl.Table:
    """
    Aggregate variants into bins using previously annotated bin information.

    Variants are grouped by:
        -'bin_id' (rank, bi-allelic, etc.)
        -'contig'
        -'snv'
        -'bi_allelic'
        -'singleton'

    For each bin, aggregates statistics needed for evaluation plots.

    :param str model_id: Which variant QC model (RF or VQSR model ID) to group
    :return: Table of aggregate statistics by bin
    """
    ht = get_score_bins(model_id, aggregated=False).ht()

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
    clinvar_pathogenic_ht = filter_to_clinvar_pathogenic(clinvar.ht())
    ht = ht.annotate(clinvar_path=hl.is_defined(clinvar_pathogenic_ht[ht.key]))
    trio_stats_ht = fam_stats.ht()

    logger.info(f"Creating grouped bin table...")
    grouped_binned_ht = compute_grouped_binned_ht(
        ht,
        checkpoint_path=get_checkpoint_path(f"grouped_bin_{model_id}"),
    )

    logger.info(f"Aggregating grouped bin table...")
    parent_ht = grouped_binned_ht._parent
    agg_ht = grouped_binned_ht.aggregate(
        n_clinvar_path=hl.agg.count_where(parent_ht.clinvar_path),
        **score_bin_agg(grouped_binned_ht, fam_stats_ht=trio_stats_ht),
    )

    return agg_ht


def get_evaluation_resources(
    test: bool,
    overwrite: bool,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC evaluation pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC evaluation pipeline.
    """
    evaluation_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc_evaluation",
        overwrite=overwrite,
        pipeline_resources={
            "annotations/generate_variant_qc_annotations.py --compute-info --split-info": {
                "info_ht": get_info(split=True)
            },
        },
    )

    extract_truth_samples = PipelineStepResourceCollection(
        "--extract-truth-samples",
        output_resources={
            f"{s}_mt": get_callset_truth_data(s) for s in TRUTH_SAMPLES_S
        },
    )
    merge_with_truth_data = PipelineStepResourceCollection(
        "--merge-with-truth-data",
        pipeline_input_steps=[extract_truth_samples],
        add_input_resources={
            "truth data high confidence intervals": {
                f"{s}_hc_intervals": TRUTH_SAMPLES[s]["hc_intervals"]
                for s in TRUTH_SAMPLES_S
            },
            "truth data matrix tables": {
                f"{s}_truth_mt": TRUTH_SAMPLES[s]["truth_mt"] for s in TRUTH_SAMPLES_S
            },
        },
        output_resources={
            f"{s}_ht": get_callset_truth_data(s, mt=False, test=test)
            for s in TRUTH_SAMPLES_S
        },
    )

    # Add all steps to the variant QC evaluation pipeline resource collection.
    evaluation_pipeline.add_steps(
        {
            "extract_truth_samples": extract_truth_samples,
            "merge_with_truth_data": merge_with_truth_data,
        }
    )

    return evaluation_pipeline


def main(args):  # noqa: D103
    hl.init(
        default_reference="GRCh38",
        log="/variant_qc_evaluation.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    overwrite = args.overwrite
    test = args.test

    evaluation_resources = get_evaluation_resources(
        test=test,
        overwrite=overwrite,
    )

    if args.create_bin_ht:
        create_bin_ht(args.model_id, args.n_bins, args.hgdp_tgp_subset).write(
            get_score_bins(
                args.model_id, aggregated=False, hgdp_tgp_subset=args.hgdp_tgp_subset
            ).path,
            overwrite=args.overwrite,
        )

    if args.run_sanity_checks:
        ht = get_score_bins(args.model_id, aggregated=False).ht()
        logger.info("Running sanity checks...")
        print(
            ht.aggregate(
                hl.struct(
                    was_biallelic=hl.agg.counter(~ht.was_split),
                    has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_bin)),
                    was_singleton=hl.agg.counter(ht.singleton),
                    has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_bin)),
                    was_biallelic_singleton=hl.agg.counter(
                        ht.singleton & ~ht.was_split
                    ),
                    has_biallelic_singleton_rank=hl.agg.counter(
                        hl.is_defined(ht.biallelic_singleton_bin)
                    ),
                )
            )
        )

    if args.create_aggregated_bin_ht:
        logger.warning("Use only workers, it typically crashes with preemptibles")
        create_aggregated_bin_ht(args.model_id).write(
            get_score_bins(args.model_id, aggregated=True).path,
            overwrite=args.overwrite,
        )

    if args.extract_truth_samples:
        logger.info(f"Extracting truth samples from VDS...")
        res = evaluation_resources.extract_truth_samples
        res.check_resource_existence()
        vds = get_gnomad_v4_vds(controls_only=True, split=True)
        # Checkpoint to prevent needing to go through the large VDS multiple times.
        vds = vds.checkpoint(new_temp_file("truth_samples", "vds"), overwrite=True)

        for s in TRUTH_SAMPLES_S:
            ts_mt = vds.variant_data
            ts_mt = ts_mt.filter_cols(ts_mt.s == s)
            # Filter to variants in truth data.
            ts_mt = ts_mt.filter_rows(hl.agg.any(ts_mt.GT.is_non_ref()))
            ts_mt.naive_coalesce(args.n_partitions)
            ts_mt.write(getattr(res, f"{s}_mt").path, overwrite=overwrite)

    if args.merge_with_truth_data:
        res = evaluation_resources.merge_with_truth_data
        res.check_resource_existence()
        for s in TRUTH_SAMPLES_S:
            logger.info(
                "Creating a merged table with callset truth sample and truth data for"
                f" {s}..."
            )
            # Remove low quality sites.
            mt = getattr(res, f"{s}_mt").mt()
            mt = mt.filter_rows(~info_ht[mt.row_key].AS_lowqual)

            # Load truth data.
            truth_mt = getattr(res, f"{s}_truth_mt").mt()
            truth_hc_intervals = getattr(res, f"{s}_hc_intervals").ht()
            truth_mt = truth_mt.key_cols_by(s=s)

            ht = create_truth_sample_ht(mt, truth_mt, truth_hc_intervals)
            ht.write(getattr(res, f"{s}_ht").path, overwrite=overwrite)

    if args.bin_truth_sample_concordance:
        for truth_sample in TRUTH_SAMPLES:
            logger.info(
                f"Creating binned concordance table for {truth_sample} for model"
                f" {args.model_id}"
            )
            ht = get_callset_truth_data(truth_sample, mt=False).ht()

            info_ht = get_info(split=True).ht()
            ht = ht.filter(
                ~info_ht[ht.key].AS_lowqual
                & ~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus])
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
                "Loading HT containing RF or VQSR scores annotated with a bin based on"
                " the rank of score..."
            )
            metric_ht = get_score_bins(args.model_id, aggregated=False).ht()
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
        help="Overwrite all data from this subset. (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--model_id",
        help="Model ID.",
        required=False,
    )
    parser.add_argument(
        "--hgdp_tgp_subset",
        help=(
            "Use HGDP + 1KG/TGP subset raw allele count as a filter in `create_bin_ht`"
            " instead of the full gnomAD raw allele count. Only relevant to VQSR"
            " because VQSR was used for filtering in all v3 releases. This is"
            " neededonly for the HGDP + 1KG/TGP subset release and is not used in"
            " evaluation."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--create_bin_ht",
        help=(
            "When set, creates file annotated with bin based on rank of VQSR/RF score."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--run_sanity_checks",
        help="When set, runs ranking sanity checks.",
        action="store_true",
    )
    parser.add_argument(
        "--create_aggregated_bin_ht",
        help=(
            "When set, creates a file with aggregate counts of variants based on bins."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--n_bins",
        help="Number of bins for the binned file (default: 100).",
        default=100,
        type=int,
    )
    parser.add_argument(
        "--extract-truth-samples",
        help="Extract truth samples from callset MatrixTable.",
        action="store_true",
    )
    parser.add_argument(
        "--n-partitions",
        help="Desired number of partitions for output Table. Default is 5000.",
        default=5000,
        type=int,
    )
    parser.add_argument(
        "--merge-with-truth-data",
        help=(
            "Computes a table for each truth sample comparing the truth sample in the"
            " callset vs the truth."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--bin_truth_sample_concordance",
        help=(
            "Merges concordance results (callset vs. truth) for a given truth sample"
            " with bins from specified model"
        ),
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
