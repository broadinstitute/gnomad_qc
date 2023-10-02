"""Script to create final filter Table for release."""
import argparse
import logging
from typing import Optional

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
logger = logging.getLogger("variant_qc_filtering")
logger.setLevel(logging.INFO)


def generate_final_filter_ht(
    ht: hl.Table,
    model_name: str,
    score_name: str,
    ac0_filter_expr: hl.expr.BooleanExpression,
    ts_ac_filter_expr: hl.expr.BooleanExpression,
    mono_allelic_flag_expr: hl.expr.BooleanExpression,
    only_het_flag_expr: hl.expr.BooleanExpression,
    inbreeding_coeff_cutoff: float = INBREEDING_COEFF_HARD_CUTOFF,
    snp_bin_cutoff: int = None,
    indel_bin_cutoff: int = None,
    snp_score_cutoff: float = None,
    indel_score_cutoff: float = None,
    aggregated_bin_ht: Optional[hl.Table] = None,
    bin_id: Optional[str] = None,
    vqsr_ht: hl.Table = None,
) -> hl.Table:
    """
    Prepare finalized filtering model given a filtering HT and cutoffs for filtering.

    .. note::

        - `snp_bin_cutoff` and `snp_score_cutoff` are mutually exclusive, and one must
          be supplied.
        - `indel_bin_cutoff` and `indel_score_cutoff` are mutually exclusive, and one
          must be supplied.
        - If a `snp_bin_cutoff` or `indel_bin_cutoff` cutoff is supplied then an
          `aggregated_bin_ht` and `bin_id` must also be supplied to determine the SNP
          and indel scores to use as cutoffs from an aggregated bin Table like one
          created by `compute_grouped_binned_ht` in combination with `score_bin_agg`.

    :param ht: Filtering Table to prepare as the final filter Table.
    :param model_name: Filtering model name to use in the 'filters' field (AS_VQSR or
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
    :param inbreeding_coeff_cutoff: InbreedingCoeff hard filter to use for variants.
    :param snp_bin_cutoff: Bin cutoff to use for SNP variant QC filter. Can't be used
        with `snp_score_cutoff`.
    :param indel_bin_cutoff: Bin cutoff to use for indel variant QC filter. Can't be
        used with `indel_score_cutoff`.
    :param snp_score_cutoff: Score cutoff (e.g. RF probability or AS_VQSLOD) to use
        for SNP variant QC filter. Can't be used with `snp_bin_cutoff`.
    :param indel_score_cutoff: Score cutoff (e.g. RF probability or AS_VQSLOD) to use
        for indel variant QC filter. Can't be used with `indel_bin_cutoff`.
    :param aggregated_bin_ht: Table with aggregate counts of variants based on bins
    :param bin_id: Name of bin to use in 'bin_id' column of `aggregated_bin_ht` to use
        to determine probability cutoff.
    :param vqsr_ht: If a VQSR HT is supplied a 'vqsr' annotation containing AS_VQSLOD,
        AS_culprit, NEGATIVE_TRAIN_SITE, and POSITIVE_TRAIN_SITE will be included in the
        returned Table.
    :return: Finalized random forest Table annotated with variant filters.
    """
    if snp_bin_cutoff is not None and snp_score_cutoff is not None:
        raise ValueError(
            "snp_bin_cutoff and snp_score_cutoff are mutually exclusive, please only"
            " supply one SNP filtering cutoff."
        )

    if indel_bin_cutoff is not None and indel_score_cutoff is not None:
        raise ValueError(
            "indel_bin_cutoff and indel_score_cutoff are mutually exclusive, please"
            " only supply one indel filtering cutoff."
        )

    if snp_bin_cutoff is None and snp_score_cutoff is None:
        raise ValueError(
            "One (and only one) of the parameters snp_bin_cutoff and snp_score_cutoff"
            " must be supplied."
        )

    if indel_bin_cutoff is None and indel_score_cutoff is None:
        raise ValueError(
            "One (and only one) of the parameters indel_bin_cutoff and"
            " indel_score_cutoff must be supplied."
        )

    if (snp_bin_cutoff is not None or indel_bin_cutoff is not None) and (
        aggregated_bin_ht is None or bin_id is None
    ):
        raise ValueError(
            "If using snp_bin_cutoff or indel_bin_cutoff, both aggregated_bin_ht and"
            " bin_id must be supplied"
        )

    # Determine SNP and indel score cutoffs if given bin instead of score
    if snp_bin_cutoff:
        snp_score_cutoff = aggregated_bin_ht.aggregate(
            hl.agg.filter(
                aggregated_bin_ht.snv
                & (aggregated_bin_ht.bin_id == bin_id)
                & (aggregated_bin_ht.bin == snp_bin_cutoff),
                hl.agg.min(aggregated_bin_ht.min_score),
            )
        )
        snp_cutoff_global = hl.struct(bin=snp_bin_cutoff, min_score=snp_score_cutoff)

    if indel_bin_cutoff:
        indel_score_cutoff = aggregated_bin_ht.aggregate(
            hl.agg.filter(
                ~aggregated_bin_ht.snv
                & (aggregated_bin_ht.bin_id == bin_id)
                & (aggregated_bin_ht.bin == indel_bin_cutoff),
                hl.agg.min(aggregated_bin_ht.min_score),
            )
        )
        indel_cutoff_global = hl.struct(
            bin=indel_bin_cutoff, min_score=indel_score_cutoff
        )

    min_score = ht.aggregate(hl.agg.min(ht.score))
    max_score = ht.aggregate(hl.agg.max(ht.score))

    if snp_score_cutoff:
        if snp_score_cutoff < min_score or snp_score_cutoff > max_score:
            raise ValueError("snp_score_cutoff is not within the range of score.")
        snp_cutoff_global = hl.struct(min_score=snp_score_cutoff)

    if indel_score_cutoff:
        if indel_score_cutoff < min_score or indel_score_cutoff > max_score:
            raise ValueError("indel_score_cutoff is not within the range of score.")
        indel_cutoff_global = hl.struct(min_score=indel_score_cutoff)

    logger.info(
        f"Using a SNP score cutoff of {snp_score_cutoff} and an indel score cutoff of"
        f" {indel_score_cutoff}."
    )

    # Add filters to HT
    filters = dict()

    if ht.any(hl.is_missing(ht.score)):
        ht.filter(hl.is_missing(ht.score)).show()
        raise ValueError("Missing Score!")

    filters[model_name] = (
        hl.is_missing(ht.score)
        | (
            hl.is_snp(ht.alleles[0], ht.alleles[1])
            & (ht.score < snp_cutoff_global.min_score)
        )
        | (
            ~hl.is_snp(ht.alleles[0], ht.alleles[1])
            & (ht.score < indel_cutoff_global.min_score)
        )
    )

    filters["InbreedingCoeff"] = hl.or_else(
        ht.inbreeding_coeff < inbreeding_coeff_cutoff, False
    )
    filters["AC0"] = ac0_filter_expr

    annotations_expr = dict()
    if model_name == "RF":
        # Fix annotations for release.
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
        only_het=only_het_flag_expr,
        **{score_name: ht.score},
        **annotations_expr,
    )

    bin_names = [x for x in ht.row if x.endswith("bin")]
    bin_names = [
        (
            x,
            (
                x.split("adj_")[0] + x.split("adj_")[1]
                if len(x.split("adj_")) == 2
                else "raw_" + x
            ),
        )
        for x in bin_names
    ]
    ht = ht.transmute(**{j: ht[i] for i, j in bin_names})

    ht = ht.annotate_globals(
        bin_stats=hl.struct(**{j: ht.bin_stats[i] for i, j in bin_names}),
        filtering_model=hl.struct(
            model_name=model_name,
            score_name=score_name,
            snv_cutoff=snp_cutoff_global,
            indel_cutoff=indel_cutoff_global,
        ),
        inbreeding_coeff_cutoff=inbreeding_coeff_cutoff,
    )
    if vqsr_ht:
        vqsr = vqsr_ht[ht.key]
        ht = ht.annotate(
            vqsr=hl.struct(
                AS_VQSLOD=vqsr.info.AS_VQSLOD,
                AS_culprit=vqsr.info.AS_culprit,
                NEGATIVE_TRAIN_SITE=vqsr.info.NEGATIVE_TRAIN_SITE,
                POSITIVE_TRAIN_SITE=vqsr.info.POSITIVE_TRAIN_SITE,
            ),
        )

    ht = ht.drop("AS_culprit")

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

    bin_ht = bin_ht.annotate(inbreeding_coeff=freq_ht[bin_ht.key].inbreeding_coeff)
    freq_idx = freq_ht[bin_ht.key]

    ac0_filter_expr = freq_idx.freq[0].AC == 0
    mono_allelic_flag_expr = (freq_idx.freq[1].AF == 1) | (freq_idx.freq[1].AF == 0)
    only_het_flag_expr = ((freq_idx.freq[0].AC * 2) == freq_idx.freq[0].AN) & (
        freq_idx.freq[0].homozygote_count == 0
    )

    ht = generate_final_filter_ht(
        bin_ht,
        args.model_name,
        args.score_name,
        ac0_filter_expr=ac0_filter_expr,
        ts_ac_filter_expr=freq_idx.freq[1].AC == 1,
        mono_allelic_flag_expr=mono_allelic_flag_expr,
        only_het_flag_expr=only_het_flag_expr,
        snp_bin_cutoff=args.snp_bin_cutoff,
        indel_bin_cutoff=args.indel_bin_cutoff,
        snp_score_cutoff=args.snp_score_cutoff,
        indel_score_cutoff=args.indel_score_cutoff,
        inbreeding_coeff_cutoff=args.inbreeding_coeff_threshold,
        aggregated_bin_ht=res.agg_bin_ht.ht(),
        bin_id="bin",
        vqsr_ht=res.vqsr_ht.ht() if args.model_id.startswith("vqsr_") else None,
    )
    ht = ht.annotate_globals(
        filtering_model=ht.filtering_model.annotate(model_id=args.model_id)
    )
    if args.model_id.startswith("vqsr_"):
        ht = ht.annotate_globals(
            filtering_model=ht.filtering_model.annotate(
                snv_training_variables=[
                    "AS_QD",
                    "AS_MQRankSum",
                    "AS_ReadPosRankSum",
                    "AS_FS",
                    "AS_SOR",
                    "AS_MQ",
                ],
                indel_training_variables=[
                    "AS_QD",
                    "AS_MQRankSum",
                    "AS_ReadPosRankSum",
                    "AS_FS",
                    "AS_SOR",
                ],
            )
        )
    else:
        ht = ht.annotate_globals(
            filtering_model=ht.filtering_model.annotate(
                snv_training_variables=ht.features,
                indel_training_variables=ht.features,
            )
        )

    final_filter_ht = ht.checkpoint(res.final_ht.path, overwrite=args.overwrite)
    final_filter_ht.summarize()


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
        "--model-name",
        help="Filtering model name to use in the filters field ('AS_VQSR', or 'RF').",
        choices=["AS_VQSR", "RF"],
    )
    parser.add_argument(
        "--score-name",
        help=(
            "What to rename the filtering score annotation. This will be used in place"
            " of 'score' in the release HT info struct and the INFO field of the VCF"
            " (e.g. 'RF', 'AS_VQSLOD')."
        ),
    )
    parser.add_argument(
        "--inbreeding-coeff-threshold",
        help="InbreedingCoeff hard filter to use for variants.",
        type=float,
        default=INBREEDING_COEFF_HARD_CUTOFF,  # TODO: Is this still a good cutoff?
    )
    snp_cutoff = parser.add_mutually_exclusive_group(required=True)
    snp_cutoff.add_argument(
        "--snp-bin-cutoff",
        help=(
            "RF or VQSR score bin to use as cutoff for SNPs. Value should be between 1"
            " and 100."
        ),
        type=float,
    )
    snp_cutoff.add_argument(
        "--snp-score-cutoff",
        help="RF or VQSR score to use as cutoff for SNPs.",
        type=float,
    )
    indel_cutoff = parser.add_mutually_exclusive_group(required=True)
    indel_cutoff.add_argument(
        "--indel-bin-cutoff",
        help=(
            "RF or VQSR score bin to use as cutoff for indels. Value should be between"
            " 1 and 100."
        ),
        type=float,
    )
    indel_cutoff.add_argument(
        "--indel-score-cutoff",
        help="RF or VQSR score to use as cutoff for indels.",
        type=float,
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
