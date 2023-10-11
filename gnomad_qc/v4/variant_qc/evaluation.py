"""Script to create Tables with aggregate variant statistics by variant QC score bins needed for evaluation plots."""
import argparse
import logging
from pprint import pformat

import hail as hl
from gnomad.utils.filtering import filter_low_conf_regions
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
from gnomad_qc.v4.resources.basics import calling_intervals, get_gnomad_v4_vds
from gnomad_qc.v4.resources.sample_qc import interval_qc_pass
from gnomad_qc.v4.resources.variant_qc import (
    TRUTH_SAMPLES,
    get_binned_concordance,
    get_callset_truth_data,
    get_score_bins,
    get_variant_qc_result,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_evaluation")
logger.setLevel(logging.INFO)


def create_bin_ht(
    ht: hl.Table,
    info_ht: hl.Table,
    rf_annotations_ht: hl.Table,
    n_bins: int,
    model_type: bool = False,
) -> hl.Table:
    """
    Create a table with bin annotations added for a variant QC run.

    :param ht: Table with variant QC result annotations.
    :param info_ht: Table with info annotations.
    :param rf_annotations_ht: Table with RF annotations.
    :param n_bins: Number of bins to bin the data into.
    :param model_type: Type of variant QC model used to annotate the data. Must be one
        of 'vqsr', 'rf', or 'if'.
    :return: Table with bin annotations.
    """
    ht = ht.filter(~info_ht[ht.key].AS_lowqual)
    non_lcr_ht = filter_low_conf_regions(
        ht,
        filter_decoy=False,
        filter_telomeres_and_centromeres=True,
    )
    struct_expr = hl.struct(**rf_annotations_ht[ht.key])
    if model_type == "vqsr":
        struct_expr = struct_expr.annotate(
            score=ht.info.AS_VQSLOD,
            positive_train_site=ht.info.POSITIVE_TRAIN_SITE,
            negative_train_site=ht.info.NEGATIVE_TRAIN_SITE,
            AS_culprit=ht.info.AS_culprit,
        )
    elif model_type == "rf":
        struct_expr = struct_expr.annotate(
            score=ht.rf_probability["TP"],
            positive_train_site=ht.tp,
            negative_train_site=ht.fp,
        )
    elif model_type == "if":
        struct_expr = struct_expr.annotate(
            score=ht.info.SCORE,
            calibration_sensitivity=ht.info.CALIBRATION_SENSITIVITY,
            calibration=ht.info.calibration,
            extracted=ht.info.extracted,
            snp=ht.info.snp,
            positive_train_site=ht.info.training,
            # We only use positive training sites for IF.
            negative_train_site=hl.missing(hl.tbool),
        )
    else:
        raise ValueError(
            f"Model type {model_type} not recognized. Must be one of 'vqsr', 'rf', or"
            " 'if'."
        )

    ht = ht.annotate(**struct_expr, non_lcr=hl.is_defined(non_lcr_ht[ht.key]))
    ht = ht.checkpoint(new_temp_file("pre_bin", "ht"), overwrite=True)
    bin_ht = create_binned_ht(
        ht,
        n_bins,
        add_substrat={
            "non_lcr": ht.non_lcr,
            "in_calling_intervals": ht.in_calling_intervals,
            "pass_interval": ht.pass_interval_qc,
        },
    )

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
    :param trio_stats_ht: Table with trio statistics.
    :return: Table of aggregate statistics by bin.
    """
    # Count variants for ranking.
    count_expr = {
        x: hl.agg.filter(
            hl.is_defined(ht[x]),
            hl.agg.counter(
                hl.if_else(hl.is_snp(ht.alleles[0], ht.alleles[1]), "snv", "indel")
            ),
        )
        for x in ht.row
        if x.endswith("bin")
    }
    bin_variant_counts = ht.aggregate(hl.struct(**count_expr))
    ht = ht.annotate_globals(bin_variant_counts=bin_variant_counts)
    logger.info(f"Found the following variant counts:\n {pformat(bin_variant_counts)}")

    logger.info(f"Creating grouped bin table...")
    # score_bin_agg expects the unrelated AC to be named ac_qc_samples_unrelated_raw.
    ht = ht.transmute(ac_qc_samples_unrelated_raw=ht.ac_unrelated_raw)
    grouped_binned_ht = compute_grouped_binned_ht(
        ht, checkpoint_path=new_temp_file(f"grouped_bin", "ht")
    )

    logger.info(f"Aggregating grouped bin table...")
    agg_ht = grouped_binned_ht.aggregate(
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
    # Initialize variant QC evaluation pipeline resource collection.
    evaluation_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc_evaluation",
        overwrite=overwrite,
        pipeline_resources={
            "annotations/generate_variant_qc_annotations.py --compute-info --split-info": {
                "info_ht": get_info(split=True)
            },
            "Interval Tables": {
                "interval_ht": calling_intervals(
                    interval_name="intersection",
                    calling_interval_padding=50,
                ),
                "interval_qc_pass_ht": interval_qc_pass(all_platforms=True),
            },
        },
    )

    # Create resource collection for each step of the variant QC evaluation pipeline.
    model_type = model_id.split("_")[0]
    if model_type == "rf":
        model_res_key = "random_forest.py --apply-rf"
    elif model_type in ["vqsr", "if"]:
        model_res_key = "GATK batch output"
    else:
        raise ValueError(
            f"Model ID {model_id} not recognized. Must start with 'rf', 'vqsr', or "
            "'if'."
        )

    create_bin_table = PipelineStepResourceCollection(
        "--create-bin-ht",
        input_resources={
            "annotations/generate_variant_qc_annotations.py --create-variant-qc-annotation-ht": {
                "vqc_annotations_ht": get_variant_qc_annotations()
            },
            model_res_key: {
                "vqc_result_ht": get_variant_qc_result(model_id=model_id, test=test)
            },
        },
        output_resources={
            "bin_ht": get_score_bins(model_id, aggregated=False, test=test),
        },
    )
    validity_check = PipelineStepResourceCollection(
        "--score-bin-validity-check",
        pipeline_input_steps=[create_bin_table],
        output_resources={},  # No output resources.
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
            f"{s}_mt": get_callset_truth_data(s, test=test) for s in TRUTH_SAMPLES
        },
    )
    merge_with_truth_data = PipelineStepResourceCollection(
        "--merge-with-truth-data",
        pipeline_input_steps=[extract_truth_samples],
        add_input_resources={
            "truth data high confidence intervals": {
                f"{s}_hc_intervals": TRUTH_SAMPLES[s]["hc_intervals"]
                for s in TRUTH_SAMPLES
            },
            "truth data matrix tables": {
                f"{s}_truth_mt": TRUTH_SAMPLES[s]["truth_mt"] for s in TRUTH_SAMPLES
            },
        },
        output_resources={
            f"{s}_ht": get_callset_truth_data(s, mt=False, test=test)
            for s in TRUTH_SAMPLES
        },
    )
    bin_truth_sample_concordance = PipelineStepResourceCollection(
        "--bin-truth-sample-concordance",
        pipeline_input_steps=[create_bin_table, merge_with_truth_data],
        output_resources={
            f"{s}_bin_ht": get_binned_concordance(model_id, s, test=test)
            for s in TRUTH_SAMPLES
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


def main(args):
    """Script to create Tables with aggregate variant statistics by variant QC score bins needed for evaluation plots."""
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
    interval_ht = evaluation_resources.interval_ht.ht()
    interval_qc_pass_ht = evaluation_resources.interval_qc_pass_ht.ht()

    if args.create_bin_ht:
        logger.info(f"Annotating {model_id} HT with bins using {n_bins} bins...")
        res = evaluation_resources.create_bin_table
        res.check_resource_existence()
        ht = res.vqc_result_ht.ht()
        ht = ht.annotate(
            in_calling_intervals=hl.is_defined(interval_ht[ht.locus]),
            pass_interval_qc=interval_qc_pass_ht[ht.locus].pass_interval_qc,
        )
        ht = create_bin_ht(
            ht,
            info_ht,
            res.vqc_annotations_ht.ht(),
            n_bins,
            model_id.split("_")[0],
        )
        ht.write(res.bin_ht.path, overwrite=overwrite)

    if args.score_bin_validity_check:
        res = evaluation_resources.validity_check
        res.check_resource_existence()
        logger.info("Running validity checks on score bins Table...")
        score_bin_validity_check(res.bin_ht.ht())

    if args.create_aggregated_bin_ht:
        res = evaluation_resources.create_aggregated_bin_table
        res.check_resource_existence()
        ht = create_aggregated_bin_ht(res.bin_ht.ht(), res.trio_stats_ht.ht())
        ht.write(res.agg_bin_ht.path, overwrite=overwrite)

    if args.extract_truth_samples:
        logger.info(f"Extracting truth samples from VDS...")
        res = evaluation_resources.extract_truth_samples
        res.check_resource_existence()
        vds = get_gnomad_v4_vds(
            controls_only=True, split=True, filter_partitions=range(2) if test else None
        )
        # Checkpoint to prevent needing to go through the large VDS multiple times.
        vds = vds.checkpoint(new_temp_file("truth_samples", "vds"))

        for ts in TRUTH_SAMPLES:
            s = TRUTH_SAMPLES[ts]["s"]
            ts_mt = vds.variant_data
            ts_mt = ts_mt.filter_cols(ts_mt.s == s)
            # Filter to variants in truth data.
            ts_mt = ts_mt.filter_rows(hl.agg.any(ts_mt.GT.is_non_ref()))
            ts_mt.naive_coalesce(args.truth_sample_mt_n_partitions)
            ts_mt.write(getattr(res, f"{ts}_mt").path, overwrite=overwrite)

    if args.merge_with_truth_data:
        res = evaluation_resources.merge_with_truth_data
        res.check_resource_existence()
        for ts in TRUTH_SAMPLES:
            s = TRUTH_SAMPLES[ts]["s"]
            logger.info(
                "Creating a merged table with callset truth sample and truth data for"
                f" {ts}..."
            )
            # Remove low quality sites.
            mt = getattr(res, f"{ts}_mt").mt()
            mt = mt.filter_rows(~info_ht[mt.row_key].AS_lowqual)

            # Load truth data.
            truth_mt = getattr(res, f"{ts}_truth_mt").mt()
            truth_hc_intervals = getattr(res, f"{ts}_hc_intervals").ht()
            truth_mt = truth_mt.key_cols_by(s=s)

            ht = create_truth_sample_ht(mt, truth_mt, truth_hc_intervals)
            ht.write(getattr(res, f"{ts}_ht").path, overwrite=overwrite)

    if args.bin_truth_sample_concordance:
        res = evaluation_resources.bin_truth_sample_concordance
        res.check_resource_existence()
        for ts in TRUTH_SAMPLES:
            logger.info(
                f"Creating binned concordance table for {ts} for model {model_id}..."
            )
            metric_ht = res.bin_ht.ht()
            ht = getattr(res, f"{ts}_ht").ht()
            ht = ht.annotate(score=metric_ht[ht.key].score)

            logger.info("Filtering out low confidence regions and segdups...")
            ht = filter_low_conf_regions(
                ht,
                filter_decoy=False,
                filter_telomeres_and_centromeres=True,
            )

            ht = ht.filter(
                hl.is_defined(interval_ht[ht.locus])
                & hl.is_defined(ht.score)
                & ~info_ht[ht.key].AS_lowqual
            )
            ht = compute_binned_truth_sample_concordance(
                ht,
                metric_ht,
                n_bins,
                add_bins={
                    "pass_interval_bin": interval_qc_pass_ht[ht.locus].pass_interval_qc,
                    "calling_interval_bin": hl.is_defined(interval_ht[ht.locus]),
                },
            )
            ht.write(getattr(res, f"{ts}_bin_ht").path, overwrite=overwrite)


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
        "--test",
        help="If the model being evaluated is a test model.",
        action="store_true",
    )
    parser.add_argument(
        "--model-id",
        help="Model ID.",
        required=False,
    )
    parser.add_argument(
        "--n-bins",
        help="Number of bins for the binned file. Default is 100.",
        default=100,
        type=int,
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
        "--extract-truth-samples",
        help="Extract truth samples from callset MatrixTable.",
        action="store_true",
    )
    parser.add_argument(
        "--truth-sample-mt-n-partitions",
        help=(
            "Desired number of partitions for th truth sample MatrixTable. Default is "
            "5000."
        ),
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
