"""Script to create final filter Table for release."""
import argparse
import logging
from typing import Dict, Optional

import hail as hl
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_freq, get_info, get_vqsr_filters
from gnomad_qc.v4.resources.variant_qc import final_filter, get_score_bins

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("final_filter")
logger.setLevel(logging.INFO)

REGION_INFO_FIELDS = ["non_lcr", "in_calling_intervals", "pass_interval_qc"]
"""Annotations to keep in the 'region_info' field of the final filter Table."""
TRUTH_SET_FIELDS = [
    "hapmap",
    "omni",
    "mills",
    "kgp_phase1_hc",
    "transmitted_singleton_raw",
    "transmitted_singleton_adj",
    "sibling_singleton_adj",
    "sibling_singleton_raw",
]
"""Annotations to keep in the 'truth_sets' field of the final filter Table."""
TRAINING_INFO_FIELDS = {
    "RF": [
        "tp",
        "fp",
        "rf_train",
        "rf_label",
        "positive_train_site",
        "negative_train_site",
    ],
    "AS_VQSR": ["positive_train_site", "negative_train_site"],
    "IF": ["positive_train_site", "training", "calibration", "extracted"],
}
"""Annotations to keep in the 'training_info' field of the final filter Table."""
VARIANT_QC_RESULT_FIELDS = {
    "RF": ["rf_tp_probability", "rf_prediction"],
    "AS_VQSR": ["AS_VQSLOD", "AS_culprit"],
    "IF": ["SCORE", "CALIBRATION_SENSITIVITY"],
}
"""Annotations to keep in the 'results' field of the final filter Table."""
FINAL_FILTER_FIELDS = [
    "a_index",
    "was_split",
    "region_info",
    "features",
    "fail_hard_filters",
    "truth_sets",
    "training_info",
    "results",
    "evaluation_score_bins",
    "singleton",
    "transmitted_singleton",
    "monoallelic",
    "only_het",
    "filters",
]
"""Top level annotations to keep in the final filter Table."""
VARIANT_QC_GLOBAL_FIELDS = {
    "standard": [
        "transmitted_singletons",
        "sibling_singletons",
        "adj",
        "interval_qc_filter",
        "compute_info_method",
    ]
}
"""Variant QC global annotations to keep in the final filter Table."""
VARIANT_QC_GLOBAL_FIELDS["RF"] = VARIANT_QC_GLOBAL_FIELDS["standard"] + [
    "feature_medians",
    "features_importance",
    "test_results",
    "fp_to_tp",
    "num_trees",
    "max_depth",
    "test_intervals",
]
VARIANT_QC_GLOBAL_FIELDS["AS_VQSR"] = VARIANT_QC_GLOBAL_FIELDS["standard"]
VARIANT_QC_GLOBAL_FIELDS["IF"] = VARIANT_QC_GLOBAL_FIELDS["standard"]
VQSR_FEATURES = {
    "snv": [
        "AS_QD",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_FS",
        "AS_MQ",
    ],
    "indel": [
        "AS_QD",
        "AS_MQRankSum",
        "AS_ReadPosRankSum",
        "AS_FS",
    ],
}
"""List of features used in the VQSR model."""
IF_FEATURES = [
    "AS_MQRankSum",
    "AS_pab_max",
    "AS_MQ",
    "AS_QD",
    "AS_ReadPosRankSum",
    "AS_FS",
    "AS_SOR",
]
"""List of features used in the isolation forest model."""


def process_score_cutoffs(
    ht: hl.Table,
    snv_bin_cutoff: int = None,
    indel_bin_cutoff: int = None,
    snv_score_cutoff: float = None,
    indel_score_cutoff: float = None,
    aggregated_bin_ht: Optional[hl.Table] = None,
    snv_bin_id: Optional[str] = None,
    indel_bin_id: Optional[str] = None,
) -> Dict[str, hl.expr.StructExpression]:
    """
    Determine SNP and indel score cutoffs if given bin instead of score.

    .. note::

        - `snv_bin_cutoff` and `snv_score_cutoff` are mutually exclusive, and one must
          be supplied.
        - `indel_bin_cutoff` and `indel_score_cutoff` are mutually exclusive, and one
          must be supplied.
        - If a `snv_bin_cutoff` or `indel_bin_cutoff` cutoff is supplied then an
          `aggregated_bin_ht` and `bin_id` must also be supplied to determine the SNP
          and indel scores to use as cutoffs from an aggregated bin Table like one
          created by `compute_grouped_binned_ht` in combination with `score_bin_agg`.

    :param ht: Filtering Table to prepare as the final filter Table.
    :param snv_bin_cutoff: Bin cutoff to use for SNP variant QC filter. Can't be used
        with `snv_score_cutoff`.
    :param indel_bin_cutoff: Bin cutoff to use for indel variant QC filter. Can't be
        used with `indel_score_cutoff`.
    :param snv_score_cutoff: Score cutoff (e.g. RF probability or AS_VQSLOD) to use
        for SNP variant QC filter. Can't be used with `snv_bin_cutoff`.
    :param indel_score_cutoff: Score cutoff (e.g. RF probability or AS_VQSLOD) to use
        for indel variant QC filter. Can't be used with `indel_bin_cutoff`.
    :param aggregated_bin_ht: Table with aggregate counts of variants based on bins
    :param snv_bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to
        determine the SNP score cutoff.
    :param indel_bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to
        determine the indel score cutoff.
    :return: Finalized random forest Table annotated with variant filters.
    """
    if snv_bin_cutoff is not None and snv_score_cutoff is not None:
        raise ValueError(
            "snv_bin_cutoff and snv_score_cutoff are mutually exclusive, please only"
            " supply one SNP filtering cutoff."
        )

    if indel_bin_cutoff is not None and indel_score_cutoff is not None:
        raise ValueError(
            "indel_bin_cutoff and indel_score_cutoff are mutually exclusive, please"
            " only supply one indel filtering cutoff."
        )

    if snv_bin_cutoff is None and snv_score_cutoff is None:
        raise ValueError(
            "One (and only one) of the parameters snv_bin_cutoff and snv_score_cutoff"
            " must be supplied."
        )

    if indel_bin_cutoff is None and indel_score_cutoff is None:
        raise ValueError(
            "One (and only one) of the parameters indel_bin_cutoff and"
            " indel_score_cutoff must be supplied."
        )

    if (
        (snv_bin_cutoff is not None and snv_bin_id is None)
        or (indel_bin_cutoff is not None and indel_bin_id is None)
    ) and (aggregated_bin_ht is None):
        raise ValueError(
            "If using snv_bin_cutoff or indel_bin_cutoff, both aggregated_bin_ht and"
            " snv_bin_id/indel_bin_id must be supplied"
        )

    cutoffs = {
        "snv": {"score": snv_score_cutoff, "bin": snv_bin_cutoff, "bin_id": snv_bin_id},
        "indel": {
            "score": indel_score_cutoff,
            "bin": indel_bin_cutoff,
            "bin_id": indel_bin_id,
        },
    }

    min_score = ht.aggregate(hl.agg.min(ht.score))
    max_score = ht.aggregate(hl.agg.max(ht.score))
    aggregated_bin_ht = aggregated_bin_ht.annotate(indel=~aggregated_bin_ht.snv)

    cutoff_globals = {}
    for variant_type, cutoff in cutoffs.items():
        if cutoff["bin"] is not None:
            score_cutoff = aggregated_bin_ht.aggregate(
                hl.agg.filter(
                    aggregated_bin_ht[variant_type]
                    & (aggregated_bin_ht.bin_id == cutoff["bin_id"])
                    & (aggregated_bin_ht.bin == cutoff["bin"]),
                    hl.agg.min(aggregated_bin_ht.min_score),
                )
            )
            cutoff_globals[variant_type] = hl.struct(
                bin=cutoff["bin"], min_score=score_cutoff, bin_id=cutoff["bin_id"]
            )
        else:
            cutoff_globals[variant_type] = hl.struct(min_score=cutoff["score"])

        score_cutoff = hl.eval(cutoff_globals[variant_type].min_score)
        if score_cutoff < min_score or score_cutoff > max_score:
            raise ValueError(
                f"{variant_type}_score_cutoff is not within the range of score ("
                f"{min_score, max_score})."
            )

    logger.info(
        f"Using a SNP score cutoff of {hl.eval(cutoff_globals['snv'].min_score)} and an"
        f" indel score cutoff of {hl.eval(cutoff_globals['indel'].min_score)}."
    )

    return cutoff_globals


def generate_final_filter_ht(
    ht: hl.Table,
    filter_name: str,
    score_name: str,
    ac0_filter_expr: hl.expr.BooleanExpression,
    ts_ac_filter_expr: hl.expr.BooleanExpression,
    mono_allelic_flag_expr: hl.expr.BooleanExpression,
    only_het_flag_expr: hl.expr.BooleanExpression,
    score_cutoff_globals: Dict[str, hl.expr.StructExpression],
    inbreeding_coeff_cutoff: float = INBREEDING_COEFF_HARD_CUTOFF,
) -> hl.Table:
    """
    Prepare finalized filtering model given a filtering HT and cutoffs for filtering.

    :param ht: Filtering Table to prepare as the final filter Table.
    :param filter_name: Filtering model name to use in the 'filters' field (AS_VQSR or
        RF).
    :param score_name: Name to use for the filtering score annotation. This will be
        used in place of 'score' in the release HT info struct and the INFO field of
        the VCF (e.g. RF or AS_VQSLOD).
    :param ac0_filter_expr: Expression that indicates if a variant should be filtered
        as allele count 0 (AC0).
    :param ts_ac_filter_expr: Allele count expression in `ht` to use as a filter for
        determining a transmitted singleton.
    :param mono_allelic_flag_expr: Expression indicating if a variant is mono-allelic.
    :param only_het_flag_expr: Expression indicating if all carriers of a variant are
        het.
    :param score_cutoff_globals: Dictionary of score cutoffs to use for filtering.
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter to use for variants.
    :return: Finalized random forest Table annotated with variant filters.
    """
    if ht.any(hl.is_missing(ht.score)):
        logger.warning("Missing Score!")
        ht.filter(hl.is_missing(ht.score)).show()

    # Construct dictionary of filters for final 'filters' annotation.
    filters = {
        "InbreedingCoeff": hl.or_else(
            ht.inbreeding_coeff < inbreeding_coeff_cutoff, False
        ),
        "AC0": ac0_filter_expr,
        filter_name: hl.is_missing(ht.score),
    }
    snv_indel_expr = {"snv": hl.is_snp(ht.alleles[0], ht.alleles[1])}
    snv_indel_expr["indel"] = ~snv_indel_expr["snv"]
    for var_type, score_cut in score_cutoff_globals.items():
        filters[filter_name] = filters[filter_name] | (
            snv_indel_expr[var_type] & (ht.score < score_cut.min_score)
        )

    variant_qc_globals = ht.index_globals()
    compute_info_method = f"{hl.eval(ht.compute_info_method)}_info"
    bin_stats_expr = variant_qc_globals.bin_group_variant_counts
    if filter_name == "RF":
        # Fix RF annotations for release.
        vqc_expr = hl.struct(
            **ht[compute_info_method],
            positive_train_site=hl.or_else(ht.positive_train_site, False),
            rf_tp_probability=ht.rf_probability["TP"],
        )
        snv_training_variables = variant_qc_globals.features
        indel_training_variables = variant_qc_globals.features
        variant_qc_globals = variant_qc_globals.annotate(
            feature_medians=variant_qc_globals.feature_medians.map_values(
                lambda x: x[compute_info_method]
            )
        )
    elif filter_name == "AS_VQSR":
        vqc_expr = hl.struct(**ht[compute_info_method])
        snv_training_variables = VQSR_FEATURES["snv"]
        indel_training_variables = VQSR_FEATURES["indel"]
    else:
        vqc_expr = ht.info.annotate(**ht[compute_info_method])
        vqc_expr = vqc_expr.drop("SCORE")
        snv_training_variables = IF_FEATURES
        indel_training_variables = IF_FEATURES

    variant_qc_globals = variant_qc_globals.select(
        filtering_model_specific_info=variant_qc_globals.select(
            *VARIANT_QC_GLOBAL_FIELDS[filter_name]
        )
    )
    keep_features = hl.eval(
        hl.set(snv_training_variables).union(hl.set(indel_training_variables))
    )

    # Restructure bin annotations as a struct with a struct for 'adj' and 'raw' bins.
    bin_names = [x for x in ht.row if x.endswith("bin")]
    bin_names = {
        "adj": {x: x.replace("adj_", "") for x in bin_names if "adj_" in x},
        "raw": {x: x for x in bin_names if "adj_" not in x},
    }
    bin_expr = {g: hl.struct(**{n[k]: ht[k] for k in n}) for g, n in bin_names.items()}

    ht = ht.annotate(
        **vqc_expr,
        **{score_name: ht.score},
        evaluation_score_bins=hl.struct(**bin_expr),
        transmitted_singleton=hl.or_missing(
            ts_ac_filter_expr, ht.transmitted_singleton_raw
        ),
        monoallelic=mono_allelic_flag_expr,
        only_het=only_het_flag_expr,
        filters=add_filters_expr(filters=filters),
    ).select_globals()

    # Restructure annotations into groups or related annotations.
    ann_groups = {
        "region_info": REGION_INFO_FIELDS,
        "truth_sets": TRUTH_SET_FIELDS,
        "features": keep_features,
        "training_info": TRAINING_INFO_FIELDS[filter_name],
        "results": VARIANT_QC_RESULT_FIELDS[filter_name],
    }
    ann_groups = {k: hl.struct(**{x: ht[x] for x in v}) for k, v in ann_groups.items()}
    ht = ht.annotate(**ann_groups)

    # Select only the fields we want to keep in the final HT in the order we want them.
    ht = ht.select(*FINAL_FILTER_FIELDS, score_name)

    # Restructure final filter Table global annotations.
    ht = ht.select_globals(
        **variant_qc_globals,
        bin_stats=hl.struct(
            **{
                g: hl.struct(**{n[k]: bin_stats_expr[k] for k in n})
                for g, n in bin_names.items()
            }
        ),
        filtering_model=hl.struct(
            filter_name=filter_name,
            score_name=score_name,
            snv_cutoff=score_cutoff_globals["snv"],
            indel_cutoff=score_cutoff_globals["indel"],
            snv_training_variables=snv_training_variables,
            indel_training_variables=indel_training_variables,
        ),
        inbreeding_coeff_cutoff=inbreeding_coeff_cutoff,
    )

    return ht


def get_final_variant_qc_resources(
    test: bool,
    overwrite: bool,
    model_id: str,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the finalizing variant QC pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :param model_id: Model ID to use for final variant QC.
    :return: PipelineResourceCollection containing resources for all steps of finalizing
        variant QC pipeline.
    """
    # Initialize finalizing variant QC pipeline resource collection.
    variant_qc_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc",
        overwrite=overwrite,
    )
    input_resources = {
        "variant_qc/evaluation.py --create-bin-ht": {
            "bin_ht": get_score_bins(model_id, aggregated=False)
        },
        "variant_qc/evaluation.py --create-aggregated-bin-ht": {
            "agg_bin_ht": get_score_bins(model_id, aggregated=True)
        },
        "annotations/generate_variant_qc_annotations.py --split-info": {
            "info_ht": get_info(split=True)
        },
        "annotations/generate_freq.py --finalize-freq-ht": {
            "freq_ht": get_freq(finalized=True)
        },
    }
    if model_id.startswith("vqsr"):
        input_resources["variant_qc/load_vqsr.py"] = {
            "vqsr_ht": get_vqsr_filters(model_id, split=True)
        }

    # Create resource collection for the finalizing variant QC pipeline.
    finalize_variant_qc = PipelineStepResourceCollection(
        "final_filter.py",
        input_resources=input_resources,
        output_resources={"final_ht": final_filter(test=test)},
    )

    # Add all steps to the finalizing variant QC pipeline resource collection.
    variant_qc_pipeline.add_steps({"finalize_variant_qc": finalize_variant_qc})

    return variant_qc_pipeline


def main(args):
    """Create final filter Table for release."""
    hl.init(
        default_reference="GRCh38",
        log="/variant_qc_finalize.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test = args.test
    overwrite = args.overwrite
    final_vqc_resources = get_final_variant_qc_resources(
        test=test,
        overwrite=overwrite,
        model_id=args.model_id,
    )
    res = final_vqc_resources.finalize_variant_qc
    res.check_resource_existence()

    bin_ht = res.bin_ht.ht()
    freq_ht = res.freq_ht.ht()

    if test:
        bin_ht = bin_ht._filter_partitions(range(5))

    bin_ht = bin_ht.filter(~res.info_ht.ht()[bin_ht.key].AS_lowqual)

    if args.model_id.startswith("vqsr_"):
        bin_ht = bin_ht.drop("info")
        filter_name = "AS_VQSR"
        score_name = "AS_VQSLOD"
    elif args.model_id.startswith("if_"):
        filter_name = "IF"  # TODO: Ask group what they want this named.
        score_name = "SCORE"  # TODO: Ask group what they want this named.
    else:
        filter_name = score_name = "RF"

    snv_bin_id = "bin"
    indel_bin_id = "bin"
    if args.snv_use_interval_qc_bin:
        snv_bin_id = "pass_interval_bin"
    elif args.snv_use_calling_interval_bin:
        snv_bin_id = "in_calling_intervals_bin"
    if args.indel_use_interval_qc_bin:
        indel_bin_id = "pass_interval_bin"
    elif args.indel_use_calling_interval_bin:
        indel_bin_id = "in_calling_intervals_bin"

    score_cutoff_globals = process_score_cutoffs(
        bin_ht,
        snv_bin_cutoff=args.snv_bin_cutoff,
        indel_bin_cutoff=args.indel_bin_cutoff,
        snv_score_cutoff=args.snv_score_cutoff,
        indel_score_cutoff=args.indel_score_cutoff,
        aggregated_bin_ht=res.agg_bin_ht.ht(),
        snv_bin_id=snv_bin_id,
        indel_bin_id=indel_bin_id,
    )

    bin_ht = bin_ht.annotate(inbreeding_coeff=freq_ht[bin_ht.key].inbreeding_coeff)
    freq_idx = freq_ht[bin_ht.key]

    mono_allelic_flag_expr = (freq_idx.freq[1].AF == 1) | (freq_idx.freq[1].AF == 0)
    only_het_flag_expr = ((freq_idx.freq[0].AC * 2) == freq_idx.freq[0].AN) & (
        freq_idx.freq[0].homozygote_count == 0
    )

    ht = generate_final_filter_ht(
        bin_ht,
        filter_name,
        score_name,
        ac0_filter_expr=freq_idx.freq[0].AC == 0,
        ts_ac_filter_expr=freq_idx.freq[1].AC == 1,
        mono_allelic_flag_expr=mono_allelic_flag_expr,
        only_het_flag_expr=only_het_flag_expr,
        score_cutoff_globals=score_cutoff_globals,
        inbreeding_coeff_cutoff=args.inbreeding_coeff_threshold,
    )
    ht = ht.annotate_globals(
        filtering_model=ht.filtering_model.annotate(model_id=args.model_id)
    )

    ht.write(res.final_ht.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            "Whether to run a test using only the first 2 partitions of the variant QC "
            "Table."
        ),
        action="store_true",
    )
    parser.add_argument("--model-id", help="Filtering model ID to use.")
    parser.add_argument(
        "--inbreeding-coeff-threshold",
        help="InbreedingCoeff hard filter to use for variants.",
        type=float,
        default=INBREEDING_COEFF_HARD_CUTOFF,
    )
    snv_cutoff = parser.add_mutually_exclusive_group(required=True)
    snv_cutoff.add_argument(
        "--snv-bin-cutoff",
        help=(
            "RF or VQSR score bin to use as cutoff for SNVs. Value should be between 1"
            " and 100."
        ),
        type=int,
    )
    snv_cutoff.add_argument(
        "--snv-score-cutoff",
        help="RF or VQSR score to use as cutoff for SNVs.",
        type=float,
    )
    parser.add_argument(
        "--snv-use-interval-qc-bin",
        help=(
            "Whether the '--snv-bin-cutoff' should be applied to the interval QC bin "
            "instead of overall bin."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--snv-use-calling-interval-bin",
        help=(
            "Whether the '--snv-bin-cutoff' should be applied to the calling interval "
            "bin instead of overall bin."
        ),
        action="store_true",
    )
    indel_cutoff = parser.add_mutually_exclusive_group(required=True)
    indel_cutoff.add_argument(
        "--indel-bin-cutoff",
        help=(
            "RF or VQSR score bin to use as cutoff for indels. Value should be between"
            " 1 and 100."
        ),
        type=int,
    )
    indel_cutoff.add_argument(
        "--indel-score-cutoff",
        help="RF or VQSR score to use as cutoff for indels.",
        type=float,
    )
    parser.add_argument(
        "--indel-use-interval-qc-bin",
        help=(
            "Whether the '--indel-bin-cutoff' should be applied to the interval QC bin "
            "instead of overall bin."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--indel-use-calling-interval-bin",
        help=(
            "Whether the '--indel-bin-cutoff' should be applied to the calling interval"
            " bin instead of overall bin."
        ),
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
