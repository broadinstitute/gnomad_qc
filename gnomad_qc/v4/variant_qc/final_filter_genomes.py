"""Script to create final filter Table for v4 genomes release."""
import argparse
import logging
from copy import deepcopy

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.utils.filtering import add_filters_expr, remove_fields_from_constant
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import ALLELE_TYPE_FIELDS
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_info, get_vqsr_filters
from gnomad_qc.v3.resources.variant_qc import get_score_bins
from gnomad_qc.v4.resources.annotations import get_freq
from gnomad_qc.v4.resources.variant_qc import VQSR_FEATURES, final_filter
from gnomad_qc.v4.variant_qc.final_filter import (
    FINAL_FILTER_FIELDS,
    TRAINING_INFO_FIELDS,
    VARIANT_QC_RESULT_FIELDS,
    process_score_cutoffs,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("final_filter_genomes")
logger.setLevel(logging.INFO)

ALLELE_TYPE_FIELDS = deepcopy(ALLELE_TYPE_FIELDS)
# Remove original_alleles for containing non-releasable alleles.
ALLELE_TYPE_FIELDS = remove_fields_from_constant(
    ALLELE_TYPE_FIELDS, ["original_alleles", "has_star"]
)
"""Annotations for allele info to keep in the final filter Table."""

REGION_INFO_FIELDS = ["non_lcr"]
"""Annotations to keep in the 'region_info' field of the final filter Table."""

TRUTH_SET_FIELDS = [
    "hapmap",
    "omni",
    "mills",
    "kgp_phase1_hc",
    "transmitted_singleton",
]
"""Annotations to keep in the 'truth_sets' field of the final filter Table."""

FINAL_FILTER_FIELDS = deepcopy(FINAL_FILTER_FIELDS)
FINAL_FILTER_FIELDS = remove_fields_from_constant(
    FINAL_FILTER_FIELDS, ["sibling_singleton"]
)
FINAL_FILTER_FIELDS = ["allele_info"] + FINAL_FILTER_FIELDS + ["AS_VQSLOD"]
"""Top level annotations to keep in the final filter Table."""

VARIANT_QC_GLOBAL_FIELDS = {
    "transmitted_singletons",
    "adj",
}
"""Variant QC global annotations to keep in the final filter Table."""


def get_pipeline_resources(
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
            "info_ht": get_info()
        },
        "annotations/generate_freq_genomes.py --finalize-freq-ht": {
            "freq_ht": get_freq(data_type="genomes", finalized=True, test=test)
        },
        "VQSR result HT": {"vqsr_ht": get_vqsr_filters(model_id, split=True)},
    }

    # Create resource collection for the finalizing variant QC pipeline.
    finalize_variant_qc = PipelineStepResourceCollection(
        "final_filter.py",
        input_resources=input_resources,
        output_resources={"final_ht": final_filter(data_type="genomes", test=test)},
    )

    # Add all steps to the finalizing variant QC pipeline resource collection.
    variant_qc_pipeline.add_steps({"finalize_variant_qc": finalize_variant_qc})

    return variant_qc_pipeline


def main(args):
    """Create final filter Table for v4 genome release."""
    hl.init(
        default_reference="GRCh38",
        log="/genome_variant_qc_finalize.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test = args.test
    overwrite = args.overwrite
    inbreeding_coeff_cutoff = args.inbreeding_coeff_threshold

    # Call method to return final VQC resources from evaluation runs.
    final_vqc_resources = get_pipeline_resources(
        test=test,
        overwrite=overwrite,
        model_id=args.model_id,
    )
    res = final_vqc_resources.finalize_variant_qc
    res.check_resource_existence()

    # Get Bin and Freq HTs, where Bin contains info about ordered bins arranged by
    # quality and freq contains frequency info for each variant.
    ht = res.bin_ht.ht()
    freq_ht = res.freq_ht.ht()

    if test:
        logger.info("Filtering to a region in PCSK9 for testing purposes...")
        ht = hl.filter_intervals(
            ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )

    # Filter out bins which fail hard cutoffs.
    ht = ht.filter(~res.info_ht.ht()[ht.key].AS_lowqual)

    if args.filter_centromere_telomere:
        ht = ht.filter(~hl.is_defined(telomeres_and_centromeres.ht()[ht.locus]))

    if ht.any(hl.is_missing(ht.score)):
        logger.warning("Missing Score!")
        ht.filter(hl.is_missing(ht.score)).show()

    # Return a dictionary containing cutoffs by variant type.
    score_cutoff_globals = process_score_cutoffs(
        ht,
        snv_bin_cutoff=args.snv_bin_cutoff,
        indel_bin_cutoff=args.indel_bin_cutoff,
        snv_score_cutoff=args.snv_score_cutoff,
        indel_score_cutoff=args.indel_score_cutoff,
        aggregated_bin_ht=res.agg_bin_ht.ht(),
        snv_bin_id="bin",
        indel_bin_id="bin",
    )

    # Append with frequency information.
    ht = ht.annotate(inbreeding_coeff=freq_ht[ht.key].InbreedingCoeff)
    freq_idx = freq_ht[ht.key]

    # Expression that indicates if a variant should be filtered as allele count 0 (AC0)
    ac0_filter_expr = freq_idx.freq[0].AC == 0

    # Allele count expression in `ht` to use as a filter for determining a transmitted
    # singleton.
    ts_ac_filter_expr = freq_idx.freq[1].AC == 1

    # Generate expressions for mono-allelic and only-het status.
    mono_allelic_flag_expr = (freq_idx.freq[0].AC > 0) & (freq_idx.freq[1].AF == 1)
    only_het_flag_expr = (
        (freq_idx.freq[0].AC > 0)
        & ((freq_idx.freq[0].AC * 2) == freq_idx.freq[0].AN)
        & (freq_idx.freq[0].homozygote_count == 0)
    )

    # Construct dictionary of filters for final 'filters' annotation.
    vqsr_filter = hl.is_missing(ht.score)
    snv_indel_expr = {"snv": hl.is_snp(ht.alleles[0], ht.alleles[1])}
    snv_indel_expr["indel"] = ~snv_indel_expr["snv"]
    for var_type, score_cut in score_cutoff_globals.items():
        vqsr_filter = vqsr_filter | (
            snv_indel_expr[var_type] & (ht.score < score_cut.min_score)
        )

    filters = {
        "InbreedingCoeff": hl.or_else(
            ht.inbreeding_coeff < inbreeding_coeff_cutoff, False
        ),
        "AC0": ac0_filter_expr,
        "AS_VQSR": vqsr_filter,
    }

    # Restructure bin annotations as a struct with a struct for 'adj' and 'raw' bins.
    bin_names = [x for x in ht.row if x.endswith("bin")]
    bin_names = {
        "adj": {x: x.replace("adj_", "") for x in bin_names if "adj_" in x},
        "raw": {x: x for x in bin_names if "adj_" not in x},
    }

    ht = ht.annotate(
        **ht.info,
        AS_VQSLOD=ht.score,
        evaluation_score_bins=hl.struct(
            **{g: hl.struct(**{n[k]: ht[k] for k in n}) for g, n in bin_names.items()}
        ),
        transmitted_singleton=hl.or_missing(
            ts_ac_filter_expr, ht.transmitted_singleton
        ),
        monoallelic=mono_allelic_flag_expr,
        only_het=only_het_flag_expr,
        filters=add_filters_expr(filters=filters),
    )

    # Restructure annotations into groups of related annotations.
    snv_features = VQSR_FEATURES["genomes"]["snv"]
    indel_features = VQSR_FEATURES["genomes"]["indel"]
    keep_features = set(snv_features) | set(indel_features)
    ann_groups = {
        "allele_info": ALLELE_TYPE_FIELDS,
        "region_info": REGION_INFO_FIELDS,
        "truth_sets": TRUTH_SET_FIELDS,
        "features": keep_features,
        "training_info": TRAINING_INFO_FIELDS["AS_VQSR"],
        "results": VARIANT_QC_RESULT_FIELDS["AS_VQSR"],
    }
    ht = ht.annotate(
        **{k: hl.struct(**{x: ht[x] for x in v}) for k, v in ann_groups.items()}
    )

    # Select only the fields we want to keep in the final HT in the order we want them.
    ht = ht.select(*FINAL_FILTER_FIELDS)

    # NOTE: This was required for v4.0 genomes to grab the SOR annotation, we now
    # compute this in `get_site_info_expr`.
    ht = ht.annotate(SOR=res.vqsr_ht.ht()[ht.key].info.SOR)

    # Restructure final filter Table global annotations.
    ht = ht.select_globals(
        filtering_model_specific_info=hl.struct(
            transmitted_singletons=True,
            adj=False,
        ),
        bin_stats=hl.struct(
            **{
                g: hl.struct(
                    **{n[k]: ht.index_globals().bin_group_variant_counts[k] for k in n}
                )
                for g, n in bin_names.items()
            }
        ),
        filtering_model=hl.struct(
            model_id=args.model_id,
            filter_name="AS_VQSR",
            score_name="AS_VQSLOD",
            snv_cutoff=score_cutoff_globals["snv"],
            indel_cutoff=score_cutoff_globals["indel"],
            snv_training_variables=snv_features,
            indel_training_variables=indel_features,
        ),
        inbreeding_coeff_cutoff=inbreeding_coeff_cutoff,
    )

    # Write out final filtered table to path defined above in resources.
    ht.write(res.final_ht.path, overwrite=args.overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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
        help="Whether to run a test using only a region in PCSK9.",
        action="store_true",
    )
    parser.add_argument(
        "--model-id",
        help="Filtering model ID to use.",
        default="vqsr_alleleSpecificTrans",
    )
    parser.add_argument(
        "--inbreeding-coeff-threshold",
        help="InbreedingCoeff hard filter to use for variants.",
        type=float,
        default=INBREEDING_COEFF_HARD_CUTOFF,
    )
    parser.add_argument(
        "--filter-centromere-telomere",
        help="Filter centromeres and telomeres from final filter Table.",
        action="store_true",
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
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
