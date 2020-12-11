import argparse
import logging
import sys
from typing import Optional, Union

import hail as hl

from gnomad.resources.grch38.reference_data import telomeres_and_centromeres, get_truth_ht
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import add_filters_expr
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources import (
    get_score_quantile_bins,
    freq,
    get_filters,
    get_info,
    final_filter,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_filtering")
logger.setLevel(logging.INFO)

INBREEDING_COEFF_HARD_CUTOFF = -0.3


def generate_final_rf_ht(
    ht: hl.Table,
    model_name: str,
    score_name: str,
    ac0_filter_expr: hl.expr.BooleanExpression,
    ts_ac_filter_expr: hl.expr.BooleanExpression,
    mono_allelic_flag_expr: hl.expr.BooleanExpression,
    snp_cutoff: Union[int, float],
    indel_cutoff: Union[int, float],
    inbreeding_coeff_cutoff: float = INBREEDING_COEFF_HARD_CUTOFF,
    determine_cutoff_from_bin: bool = False,
    aggregated_bin_ht: Optional[hl.Table] = None,
    bin_id: Optional[str] = None,
    vqsr_ht: hl.Table = None,
) -> hl.Table:
    """
    Prepares finalized RF model given an RF result table from `rf.apply_rf_model` and cutoffs for filtering.
    
    If `determine_cutoff_from_bin` is True, `aggregated_bin_ht` must be supplied to determine the SNP and indel RF
    probabilities to use as cutoffs from an aggregated quantile bin Table like one created by
    `compute_grouped_binned_ht` in combination with `score_bin_agg`.

    :param ht: RF result table from `rf.apply_rf_model` to prepare as the final RF Table
    :param ac0_filter_expr: Expression that indicates if a variant should be filtered as allele count 0 (AC0)
    :param ts_ac_filter_expr: Expression in `ht` that indicates if a variant is a transmitted singleton
    :param mono_allelic_flag_expr: Expression indicating if a variant is mono-allelic
    :param snp_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for SNP variant QC filter
    :param indel_cutoff: RF probability or bin (if `determine_cutoff_from_bin` True) to use for indel variant QC filter
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter to use for variants
    :param determine_cutoff_from_bin: If True RF probability will be determined using bin info in `aggregated_bin_ht`
    :param aggregated_bin_ht: File with aggregate counts of variants based on quantile bins
    :param bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to use to determine probability cutoff
    :return: Finalized random forest Table annotated with variant filters
    """
    # Determine SNP and indel RF cutoffs if given bin instead of RF probability
    if determine_cutoff_from_bin:
        snp_rf_cutoff, indel_rf_cutoff = aggregated_bin_ht.aggregate(
            [
                hl.agg.filter(
                    snv
                    & (aggregated_bin_ht.bin_id == bin_id)
                    & (aggregated_bin_ht.bin == cutoff),
                    hl.agg.min(aggregated_bin_ht.min_score),
                )
                for snv, cutoff in [
                    (aggregated_bin_ht.snv, snp_cutoff),
                    (~aggregated_bin_ht.snv, indel_cutoff),
                ]
            ]
        )
        snp_cutoff_global = hl.struct(bin=snp_cutoff, min_score=snp_rf_cutoff)
        indel_cutoff_global = hl.struct(bin=indel_cutoff, min_score=indel_rf_cutoff)

        logger.info(
            f"Using a SNP RF probability cutoff of {snp_rf_cutoff} and an indel RF probability cutoff of {indel_rf_cutoff}."
        )
    else:
        snp_cutoff_global = hl.struct(min_score=snp_cutoff)
        indel_cutoff_global = hl.struct(min_score=indel_cutoff)

    # Add filters to RF HT
    filters = dict()

    if "VQSR" in model_name:
        ht2 = ht.filter(hl.is_missing(ht.filters) | ~ht.filters.contains("LowQual"))
    if ht2.any(hl.is_missing(ht2.score)):
        ht2.filter(hl.is_missing(ht2.score)).show()
        raise ValueError("Missing Score!")

    filters[model_name] = hl.is_missing(ht.score) | (
        hl.is_snp(ht.alleles[0], ht.alleles[1])
        & (ht.score < snp_cutoff_global.min_score)
    ) | (
        ~hl.is_snp(ht.alleles[0], ht.alleles[1])
        & (ht.score < indel_cutoff_global.min_score)
    )

    filters["InbreedingCoeff"] = hl.or_else(
        ht.InbreedingCoeff < inbreeding_coeff_cutoff, False
    )
    filters["AC0"] = ac0_filter_expr

    annotations_expr = dict()
    if model_name == "RF":
        # Fix annotations for release
        annotations_expr = annotations_expr.update(
            {
                "positive_train_site": hl.or_else(ht.positive_train_site, False),
                "rf_tp_probability": ht.rf_probability["TP"],
            }
        )
    annotations_expr.update(
        {
            "transmitted_singleton": hl.or_missing(
                ts_ac_filter_expr, ht.transmitted_singleton
            )
        }
    )
    if "feature_imputed" in ht.row:
        annotations_expr.update(
            {
                x: hl.or_missing(~ht.feature_imputed[x], ht[x])
                for x in [f for f in ht.row.feature_imputed]
            }
        )

    ht = ht.transmute(
        filters=add_filters_expr(filters=filters),
        monoallelic=mono_allelic_flag_expr,
        **{score_name: ht.score},
        **annotations_expr
    )

    bin_names = [x for x in ht.row if x.endswith("bin")]
    bin_names = [(x, x.split("adj_")[0] + x.split("adj_")[1] if len(x.split("adj_")) == 2 else "raw_" + x) for x in bin_names]
    ht = ht.transmute(
        **{j: ht[i] for i, j in bin_names}
    )

    ht = ht.annotate_globals(
        bin_stats=hl.struct(**{j: ht.bin_stats[i] for i, j in bin_names}),
        filtering_model=hl.struct(
            model_name=model_name,
            score_name=score_name,
            snv_cutoff=snp_cutoff_global,
            indel_cutoff=indel_cutoff_global,
        )
    )
    vqsr = vqsr_ht[ht.key]
    ht = ht.annotate(
        vqsr=hl.struct(
            AS_VQSLOD=vqsr.info.AS_VQSLOD,
            AS_culprit=vqsr.info.AS_culprit,
            NEGATIVE_TRAIN_SITE=vqsr.info.NEGATIVE_TRAIN_SITE,
            POSITIVE_TRAIN_SITE=vqsr.info.POSITIVE_TRAIN_SITE,
        ),
        SOR=vqsr.info.SOR,
    )

    ht = ht.drop('AS_culprit')

    return ht


def main(args):
    hl.init(log="/variant_qc_finalize.log")

    ht = get_score_quantile_bins(args.model_id, aggregated=False).ht()
    if args.filter_centromere_telomere:
        ht = ht.filter(~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus]))

    info_ht = get_info(split=True).ht()
    ht = ht.filter(~info_ht[ht.key].AS_lowqual)

    if args.model_id.startswith('vqsr_'):
        ht = ht.drop("info")

    freq_ht = freq.ht()
    ht = ht.annotate(InbreedingCoeff=freq_ht[ht.key].InbreedingCoeff)
    freq_idx = freq_ht[ht.key]
    aggregated_bin_path = get_score_quantile_bins(args.model_id, aggregated=True).path
    if not file_exists(aggregated_bin_path):
        sys.exit(
            f"Could not find binned HT for RF  run {args.model_id} ({aggregated_bin_path}). Please run create_ranked_scores.py for that hash."
        )
    aggregated_bin_ht = get_score_quantile_bins(args.model_id, aggregated=True).ht()

    ht = generate_final_rf_ht(
        ht,
        args.model_name,
        args.score_name,
        ac0_filter_expr=freq_idx.freq[0].AC == 0,
        ts_ac_filter_expr=freq_idx.freq[1].AC == 1,
        mono_allelic_flag_expr=(freq_idx.freq[1].AF == 1) | (freq_idx.freq[1].AF == 0),
        snp_cutoff=args.snp_cutoff,
        indel_cutoff=args.indel_cutoff,
        determine_cutoff_from_bin=not args.treat_cutoff_as_prob,
        aggregated_bin_ht=aggregated_bin_ht,
        bin_id='bin',
        inbreeding_coeff_cutoff=args.inbreeding_coeff_threshold,
        vqsr_ht=get_filters('vqsr_alleleSpecificTrans', split=True).ht(),
    )
    ht = ht.annotate_globals(
        filtering_model=ht.filtering_model.annotate(
            model_id=args.model_id,
        )
    )
    if args.model_id.startswith('vqsr_'):
        ht = ht.annotate_globals(
            filtering_model=ht.filtering_model.annotate(
                snv_training_variables=["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_SOR", "AS_MQ"],
                indel_training_variables=["AS_QD", "AS_MQRankSum", "AS_ReadPosRankSum", "AS_FS", "AS_SOR"],
        )
        )
    else:
        ht = ht.annotate_globals(
            filtering_model=ht.filtering_model.annotate(
                snv_training_variables=ht.features,
                indel_training_variables=ht.features,
        )
        )

    ht.write(final_filter.path, args.overwrite)

    final_filter_ht = final_filter.ht()
    final_filter_ht.summarize()


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
    parser.add_argument("--model_id", help="Filtering model ID to use")
    parser.add_argument(
        "--model_name",
        help="Filtering model name to use in the filters. (e.g. VQSR or RF)",
    )
    parser.add_argument(
        "--score_name",
        help="Name for the score in the info struct in HT / INFO field inVCF (e.g. RF or VQSLOD",
    )
    parser.add_argument(
        "--inbreeding_coeff_threshold",
        help="Name for the score in the info struct in HT / INFO field inVCF (e.g. RF or VQSLOD",
        type=float,
        default=INBREEDING_COEFF_HARD_CUTOFF,
    )
    parser.add_argument(
        "--snp_cutoff", help="Percentile to set RF cutoff", type=float, default=90.0
    )
    parser.add_argument(
        "--indel_cutoff", help="Percentile to set RF cutoff", type=float, default=80.0
    )
    parser.add_argument(
        "--treat_cutoff_as_prob",
        help="If set snp_cutoff and indel_cutoff will be probability rather than percentile ",
        action="store_true",
    )
    parser.add_argument(
        "--filter_centromere_telomere",
        help="Train RF without centromeres and telomeres.",
        action="store_true",
    )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)