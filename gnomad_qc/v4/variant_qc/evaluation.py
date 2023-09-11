# noqa: D100

import argparse
import logging
from pprint import pformat

import hail as hl
from gnomad.resources.grch38.reference_data import clinvar
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
from gnomad_qc.v4.resources.annotations import (
    get_info,
    get_trio_stats,
    get_variant_qc_annotations,
)
from gnomad_qc.v4.resources.basics import TRUTH_SAMPLES_S, get_gnomad_v4_vds
from gnomad_qc.v4.resources.variant_qc import (
    TRUTH_SAMPLES,
    get_binned_concordance,
    get_callset_truth_data,
    get_rf_result,
    get_score_bins,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_evaluation")
logger.setLevel(logging.INFO)


def create_bin_ht(
    ht: hl.Table,
    info_ht: hl.Table,
    rf_annotations_ht: hl.Table,
    n_bins: int,
    vqsr: bool = False,
) -> hl.Table:
    """
    Create a table with bin annotations added for a RF or VQSR run and writes it to its correct location in annotations.

    :param n_bins: Number of bins to bin the data into.
    :return: Table with bin annotations.
    """
    if vqsr:
        ht = ht.annotate(
            **rf_annotations_ht[ht.key],
            info=info_ht[ht.key].info,
            score=ht.info.AS_VQSLOD,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            AS_culprit=ht.info.AS_culprit,
        )
    else:
        ht = ht.annotate(
            score=ht.rf_probability["TP"],
            info=info_ht[ht.key].info,
            positive_train_site=ht.tp,
            negative_train_site=ht.fp,
        )

    ht = ht.filter(~info_ht[ht.key].AS_lowqual)
    ht_non_lcr = filter_low_conf_regions(
        ht,
        filter_decoy=False,
        filter_telomeres_and_centromeres=True,
    )
    ht = ht.annotate(non_lcr=hl.is_defined(ht_non_lcr[ht.key]))
    bin_ht = create_binned_ht(ht, n_bins, add_substrat={"non_lcr": ht.non_lcr})

    return bin_ht


def score_bin_validity_check(ht: hl.Table) -> None:
    """
    Check that the bin scoring looks correct by printing so aggregate numbers.

    :param ht: Table with bin annotations.
    :return: None
    """
    print(
        ht.aggregate(
            hl.struct(
                was_biallelic=hl.agg.counter(~ht.was_split),
                has_biallelic_rank=hl.agg.counter(hl.is_defined(ht.biallelic_bin)),
                was_singleton=hl.agg.counter(ht.singleton),
                has_singleton_rank=hl.agg.counter(hl.is_defined(ht.singleton_bin)),
                was_biallelic_singleton=hl.agg.counter(ht.singleton & ~ht.was_split),
                has_biallelic_singleton_rank=hl.agg.counter(
                    hl.is_defined(ht.biallelic_singleton_bin)
                ),
            )
        )
    )


def create_aggregated_bin_ht(ht: hl.Table, trio_stats_ht: hl.Table) -> hl.Table:
    """
    Aggregate variants into bins using previously annotated bin information.

    Variants are grouped by:
        -'bin_id' (rank, bi-allelic, etc.)
        -'contig'
        -'snv'
        -'bi_allelic'
        -'singleton'

    For each bin, aggregates statistics needed for evaluation plots.

    :param ht: Table with bin annotations.
    :return: Table of aggregate statistics by bin.
    """
    # Count variants for ranking.
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
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)
    logger.info(f"Found the following variant counts:\n {pformat(bin_variant_counts)}")

    # Load ClinVar pathogenic data.
    clinvar_pathogenic_ht = filter_to_clinvar_pathogenic(clinvar.ht())
    ht = ht.annotate(clinvar_path=hl.is_defined(clinvar_pathogenic_ht[ht.key]))

    logger.info(f"Creating grouped bin table...")
    grouped_binned_ht = compute_grouped_binned_ht(
        ht,
        checkpoint_path=new_temp_file(f"grouped_bin", "ht"),
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
    model_id: str = None,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC evaluation pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :param model_id: Model ID of RF model results to use.
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC evaluation pipeline.
    """
    # if vqsr:
    #    ht = get_vqsr_filters(model_id, split=True).ht()
    # else:
    #    ht = get_rf_result(model_id=model_id).ht()
    # Initialize variant QC evaluation pipeline resource collection.
    evaluation_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc_evaluation",
        overwrite=overwrite,
        pipeline_resources={
            "annotations/generate_variant_qc_annotations.py --compute-info --split-info": {
                "info_ht": get_info(split=True)
            },
        },
    )

    # Create resource collection for each step of the variant QC evaluation pipeline.
    create_bin_table = PipelineStepResourceCollection(
        "--create-bin-ht",
        input_resources={
            "annotations/generate_variant_qc_annotations.py --create-variant-qc-annotation-ht": {
                "rf_ht": get_variant_qc_annotations()
            },
            "random_forest.py --apply-rf": {
                "rf_result_ht": get_rf_result(model_id=model_id, test=test)
            },
        },
        output_resources={
            "bin_ht": get_score_bins(model_id, aggregated=False, test=test),
        },
    )
    validity_check = PipelineStepResourceCollection(
        "--score-bin-validity-check",
        pipeline_input_steps=[create_bin_table],
    )
    create_aggregated_bin_table = PipelineStepResourceCollection(
        "--create-aggregated-bin-ht",
        output_resources={
            "agg_bin_ht": get_score_bins(model_id, aggregated=True, test=test)
        },
        pipeline_input_steps=[create_bin_table],
        add_input_resources={
            "annotations/generate_variant_qc_annotations.py --generate-trio-stats": {
                "trio_stats_ht": get_trio_stats()
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
    bin_truth_sample_concordance = PipelineStepResourceCollection(
        "--bin-truth-sample-concordance",
        pipeline_input_steps=[create_bin_table, merge_with_truth_data],
        output_resources={
            f"{s}_bin_ht": get_binned_concordance(model_id, s, test=test)
            for s in TRUTH_SAMPLES_S
        },
    )

    # Add all steps to the variant QC evaluation pipeline resource collection.
    evaluation_pipeline.add_steps(
        {
            "create_bin_table": create_bin_table,
            "validity_check": validity_check,
            "create_aggregated_bin_table": create_aggregated_bin_table,
            "extract_truth_samples": extract_truth_samples,
            "merge_with_truth_data": merge_with_truth_data,
            "bin_truth_sample_concordance": bin_truth_sample_concordance,
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
    model_id = args.model_id
    n_bins = args.n_bins

    evaluation_resources = get_evaluation_resources(
        test=test,
        overwrite=overwrite,
        model_id=model_id,
    )
    info_ht = evaluation_resources.info_ht.ht()

    if args.create_bin_ht:
        logger.info(f"Annotating {model_id} HT with bins using {n_bins} bins...")
        res = evaluation_resources.create_bin_table
        res.check_resource_existence()
        ht = create_bin_ht(res.rf_result_ht.ht(), info_ht, res.rf_ht.ht(), n_bins)
        ht.write(res.bin_ht.path, overwrite=overwrite)

    if args.score_bin_validity_check:
        logger.info("Running validity checks on score bins Table...")
        score_bin_validity_check(res.bin_ht.ht())

    if args.create_aggregated_bin_ht:
        logger.warning("Use only workers, it typically crashes with preemptibles")
        res = evaluation_resources.create_aggregated_bin_table
        res.check_resource_existence()
        ht = create_aggregated_bin_ht(res.bin_ht.ht(), res.trio_stats_ht.ht())
        ht.write(res.agg_bin_ht.path, overwrite=overwrite)

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
        res = evaluation_resources.bin_truth_sample_concordance
        res.check_resource_existence()
        for s in TRUTH_SAMPLES_S:
            logger.info(
                f"Creating binned concordance table for {s} for model {model_id}..."
            )
            metric_ht = res.bin_ht.ht()
            ht = getattr(res, f"{s}_ht").ht()
            ht = ht.annotate(score=metric_ht[ht.key].score)

            logger.info("Filtering out low confidence regions and segdups...")
            ht = filter_low_conf_regions(
                ht,
                filter_decoy=False,
                filter_telomeres_and_centromeres=True,
            )

            ht = ht.filter(hl.is_defined(ht.score) & ~info_ht[ht.key].AS_lowqual)
            ht = compute_binned_truth_sample_concordance(ht, metric_ht, n_bins)
            ht.write(getattr(res, f"{s}_bin_ht").path, overwrite=overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset. (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--model-id",
        help="Model ID.",
        required=False,
    )
    parser.add_argument(
        "--create-bin-ht",
        help=(
            "When set, creates file annotated with bin based on rank of VQSR/RF score."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--score-bin-validity-check",
        help="When set, runs ranking validity checks.",
        action="store_true",
    )
    parser.add_argument(
        "--create-aggregated-bin-ht",
        help=(
            "When set, creates a file with aggregate counts of variants based on bins."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--n-bins",
        help="Number of bins for the binned file. Default is 100.",
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
        "--bin-truth-sample-concordance",
        help=(
            "Merges concordance results (callset vs. truth) for a given truth sample"
            " with bins from specified model."
        ),
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
