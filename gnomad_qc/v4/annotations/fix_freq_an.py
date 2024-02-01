"""Updates the freq HT with the the correct AN."""

import argparse
import logging
from copy import deepcopy
from typing import Optional

import hail as hl
from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    POPS_TO_REMOVE_FOR_POPMAX,
    release_all_sites_an,
)
from gnomad.utils.annotations import merge_freq_arrays, merge_histograms, qual_hist_expr
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SORT_ORDER

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.annotations.generate_freq import (
    compute_inbreeding_coeff,
    correct_for_high_ab_hets,
    create_final_freq_ht,
    generate_faf_grpmax,
)
from gnomad_qc.v4.resources.annotations import get_freq
from gnomad_qc.v4.resources.basics import get_logging_path

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_exomes_AN_fix")
logger.setLevel(logging.INFO)


AGE_HISTS = [
    "age_hist_het",
    "age_hist_hom",
]
"""Age histograms to compute and keep on the frequency Table."""
QUAL_HISTS = [
    "gq_hist_all",
    "dp_hist_all",
    "gq_hist_alt",
    "dp_hist_alt",
    "ab_hist_alt",
]
"""Quality histograms to compute and keep on the frequency Table."""
FREQ_HIGH_AB_HET_ROW_FIELDS = [
    "high_ab_hets_by_group",
    "high_ab_het_adjusted_ab_hists",
    "high_ab_het_adjusted_age_hists",
]
"""
List of top level row and global annotations relating to the high allele balance
heterozygote correction that we want on the frequency HT before deciding on the AF
cutoff.
"""
FREQ_ROW_FIELDS = [
    "freq",
    "qual_hists",
    "raw_qual_hists",
    "age_hists",
]
"""
List of top level row and global annotations with no high allele balance heterozygote
correction that we want on the frequency HT.
"""
ALL_FREQ_ROW_FIELDS = FREQ_ROW_FIELDS + FREQ_HIGH_AB_HET_ROW_FIELDS
"""
List of final top level row and global annotations created from dense data that we
want on the frequency HT before deciding on the AF cutoff.
"""
FREQ_GLOBAL_FIELDS = [
    "downsamplings",
    "freq_meta",
    "age_distribution",
    "freq_index_dict",
    "freq_meta_sample_count",
]
"""
List of final global annotations created from dense data that we want on the frequency
HT before deciding on the AF cutoff.
"""
SUBSET_DICT = {"gnomad": 0, "non_ukb": 1}
"""
Dictionary for accessing the annotations with subset specific annotations such as
age_hists, popmax, and faf.
"""


def get_freq_an_resources(
    overwrite: bool = False,
    test: Optional[bool] = False,
    chrom: Optional[str] = None,
) -> PipelineResourceCollection:
    """
    Get update frequency HT AN resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :param chrom: Chromosome used in freq calculations.
    :return: Update frequency HT AN resources.
    """
    update_freq_an_pipeline = PipelineResourceCollection(
        pipeline_name="frequency_an_fix",
        overwrite=overwrite,
    )

    update_freq_an = PipelineStepResourceCollection(
        "--update-freq-an",
        pipeline_input_steps={
            "annotations/generate_freq.py --combine-freq-hts": {
                "freq_ht": get_freq(
                    version="4.0",
                    test=test,
                    hom_alt_adjusted=False,
                    chrom=chrom,
                    finalized=False,
                )
            },
            "annotations/compute_coverage.py --compute-allele-number-ht": {
                "an_ht": release_all_sites_an(public=False, test=test)
            },
        },
        output_resources={
            "updated_an_freq_ht": get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=False,
                chrom=chrom,
                finalized=False,
            )
        },
    )
    update_af_dependent_annotations = PipelineStepResourceCollection(
        "--update-af-dependent-annotations",
        pipeline_input_steps=[update_freq_an],
        output_resources={
            "updated_freq_ht": get_freq(
                version="4.1",
                test=test,
                hom_alt_adjusted=True,
                chrom=chrom,
                finalized=False,
            )
        },
    )
    finalize_freq_ht = PipelineStepResourceCollection(
        "--finalize-freq-ht",
        pipeline_input_steps=[update_af_dependent_annotations],
        output_resources={
            "final_freq_ht": get_freq(version="4.1", test=test, finalized=True)
        },
    )

    update_freq_an_pipeline.add_steps(
        {
            "update_freq_an": update_freq_an,
            "update_af_dependent_annotations": update_af_dependent_annotations,
            "finalize_freq_ht": finalize_freq_ht,
        }
    )
    return update_freq_an_pipeline


def update_freq_an(freq_ht: hl.Table, an_ht: hl.Table) -> hl.Table:
    """
    Update frequency HT AN field with the correct AN from the AN HT.

    This module updates the freq_ht AN field which is an array of ints with the correct
    AN from the AN HT integer array annotation. It uses the an_ht dictionary and the
    freq_ht dictionary to match the indexing of array annotations so each group is
    correctly updated. It then recomputes AF, AC/AN.

    :param freq_ht: Frequency HT to update.
    :param an_ht: AN HT to update from.
    :return: Updated frequency HT.
    """
    freq_meta = freq_ht.index_globals().freq_meta
    an_meta = an_ht.index_globals().strata_meta

    # Set all freq_ht's ANs to 0 so can use the merge_freq function to update
    # the ANs (expects int64s)
    freq_ht = freq_ht.annotate(
        freq=freq_ht.freq.map(
            lambda x: x.select(
                AC=x.AC, AN=hl.int64(0), homozygote_count=x.homozygote_count, AF=x.AF
            )
        )
    )

    # Create freq array annotation of structs with AC and nhomozygotes set to 0 to
    # use merge_freq function (expects int64s)
    an_ht = an_ht.transmute(
        freq=an_ht.AN.map(
            lambda x: hl.struct(
                AC=0, AN=x, homozygote_count=0, AF=hl.missing(hl.tfloat64)
            )
        )
    )

    # Add an_ht to freq_ht to call merge_freq_arrays
    freq_ht = freq_ht.annotate(an_freq=an_ht[freq_ht.locus].freq)
    freq_ht = freq_ht.annotate_globals(an_strata_meta=an_meta)
    freq, freq_meta = merge_freq_arrays(
        [freq_ht.freq, freq_ht.an_freq], [freq_ht.freq_meta, freq_ht.an_strata_meta]
    )
    freq_ht = freq_ht.annotate(freq=freq)
    freq_ht = freq_ht.annotate_globals(freq_meta=freq_meta)
    return freq_ht


def drop_gatk_groupings(ht: hl.Table) -> hl.Table:
    """
    Drop GATK groupings from the frequency HT.

    :param ht: Frequency HT to drop GATK groupings from.
    :return: Frequency HT with GATK groupings dropped.
    """
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            "freq": ht.freq,
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=["gatk_version"],
        keep=False,
    )
    ht = ht.annotate(freq=array_exprs["freq"])
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )
    return ht


def main(args):
    """Script to generate frequency and dense dependent annotations on v4 exomes."""
    overwrite = args.overwrite
    use_test_dataset = args.use_test_dataset
    test_n_partitions = args.test_n_partitions
    test_gene = args.test_gene
    test = use_test_dataset or test_n_partitions or test_gene
    chrom = args.chrom
    af_threshold = args.af_threshold

    hl.init(
        log="/frequency_an_update.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")
    resources = get_freq_an_resources(overwrite, test, chrom)

    try:
        if args.update_freq_an:
            logger.info("Updating AN in freq HT...")
            res = resources.update_freq_an
            res.check_resource_existence()
            freq_ht = res.freq_ht.ht()
            freq_ht = drop_gatk_groupings(freq_ht)
            an_ht = res.an_ht.ht()
            ht = update_freq_an(freq_ht, an_ht)
            ht.write(res.updated_an_freq_ht.path, overwrite=overwrite)

        if args.update_af_dependent_annotations:
            logger.info("Updatings AF dependent annotations...")
            res = resources.update_af_dependent_annotations
            res.check_resource_existence()
            ht = res.freq_ht.ht()

            logger.info("Correcting call stats, qual AB hists, and age hists...")
            ht = correct_for_high_ab_hets(ht, af_threshold=af_threshold)

            logger.info("Computing FAF & grpmax...")
            ht = generate_faf_grpmax(ht)

            logger.info("Calculating InbreedingCoeff...")
            ht = compute_inbreeding_coeff(ht)

            logger.info("High AB het corrected frequency HT schema...")
            ht.describe()

            logger.info("Writing corrected frequency Table...")
            ht.write(res.corrected_freq_ht.path, overwrite=overwrite)

        # Fill this in, can we do this on VDS, run hists over MT then merge hists?
        if args.regenerate_gq_dp_hists:
            pass

        if args.finalize_freq_ht:
            logger.info("Writing final frequency Table...")
            res = resources.finalize_freq_ht
            res.check_resource_existence()
            ht = create_final_freq_ht(res.corrected_freq_ht.ht())

            logger.info("Final frequency HT schema...")
            ht.describe()
            ht.write(res.final_freq_ht.path, overwrite=overwrite)
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("frequency_data"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--use-test-dataset",
        help="Runs a test on the gnomad test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--test-gene",
        help="Runs a test on the DRD2 gene in the gnomad test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-partitions",
        help=(
            "Use only N partitions of the VDS as input for testing purposes. Defaults"
            "to 2 if passed without a value."
        ),
        nargs="?",
        const=2,
        type=int,
    )
    parser.add_argument(
        "--chrom",
        help="If passed, script will only run on passed chromosome.",
        type=str,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--run-freq-and-dense-annotations",
        help=(
            "Calculate frequencies, histograms, and high AB sites per sample grouping."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--combine-freq-hts",
        help=(
            "Combine frequency and histogram Tables for UKB and non-UKB samples into a"
            " single Table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--correct-for-high-ab-hets",
        help=(
            "Correct each frequency entry to account for homozygous alternate depletion"
            " present in GATK versions released prior to 4.1.4.1 and run chosen"
            " downstream annotations."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--af-threshold",
        help=(
            "Threshold at which to adjust site group frequencies at sites for"
            " homozygous alternate depletion present in GATK versions released prior to"
            " 4.1.4.1."
        ),
        type=float,
        default=0.01,
    )
    parser.add_argument(
        "--finalize-freq-ht",
        help=(
            "Finalize frequency Table by dropping unnecessary fields and renaming"
            " remaining fields."
        ),
        action="store_true",
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
