# noqa: D100

import argparse
import logging
from copy import deepcopy
from typing import Dict, List, Optional, Set, Union

import hail as hl
from gnomad.assessment.validity_checks import validate_release_t
from gnomad.resources.grch38.gnomad import POPS, SEXES, SUBSETS
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.file_utils import file_exists
from gnomad.utils.filtering import remove_fields_from_constant
from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    AS_VQSR_FIELDS,
    ENTRIES,
    FAF_POPS,
    FORMAT_DICT,
    HISTS,
    IN_SILICO_ANNOTATIONS_INFO_DICT,
    INFO_DICT,
    REGION_FLAG_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
    VRS_FIELDS_DICT,
    add_as_info_dict,
    adjust_vcf_incompatible_types,
    build_vcf_export_reference,
    create_label_groups,
    make_hist_bin_edges_expr,
    make_hist_dict,
    make_info_dict,
    make_vcf_filter_dict,
    rekey_new_reference,
)
from gnomad.utils.vep import VEP_CSQ_HEADER, vep_struct_to_csq
from gnomad.variant_qc.pipeline import INBREEDING_COEFF_HARD_CUTOFF

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v4.resources.release import (
    append_to_vcf_header_path,
    release_header_path,
    release_sites,
    release_vcf_path,
    validated_release_ht,
)

POPS = deepcopy(POPS["v4"])
SUBSETS = deepcopy(SUBSETS["v4"])

# Remove unnecessary pop names from POP_NAMES dict
POPS = {pop: POP_NAMES[pop] for pop in POPS}

# Remove unnecessary pop names from FAF_POPS dict
FAF_POPS = {pop: POP_NAMES[pop] for pop in FAF_POPS}

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("vcf_release")
logger.setLevel(logging.INFO)


def get_export_resources(
    overwrite: bool = False,
    test: Optional[bool] = False,
) -> PipelineResourceCollection:
    """
    Get export resources.

    :param overwrite: Whether to overwrite existing files.
    :param test: Whether to use test resources.
    :return: Export resources.
    """
    export_pipeline = PipelineResourceCollection(
        pipeline_name="validate_and_export_vcf",
        overwrite=overwrite,
    )
    validate_release = PipelineStepResourceCollection(
        "--validate-release-ht",
        output_resources={"validated_ht": validated_release_ht(test=test)},
    )
    export_pipeline.add_steps(
        {
            "validate_ht": validate_release,
        }
    )
    return export_pipeline


def filter_to_test(ht: hl.Table, num_partitions: int = 2) -> hl.Table:
    """
    Filter Table to `num_partitions` partitions on chr20, chrX, and chrY for testing.

    :param t: Input Table to filter.
    :param num_partitions: Number of partitions to grab from each chromosome.
    :return: Input Table filtered to `num_partitions` on chr20, chrX, and chrY.
    """
    logger.info(
        "Filtering to %d partitions on chr20, chrX, and chrY (for tests only)...",
        num_partitions,
    )
    ht_chr20 = hl.filter_intervals(ht, [hl.parse_locus_interval("chr20")])
    ht_chr20 = ht_chr20._filter_partitions(range(num_partitions))

    ht_chrx = hl.filter_intervals(ht, [hl.parse_locus_interval("chrX")])
    ht_chrx = ht_chrx._filter_partitions(range(num_partitions))

    ht_chry = hl.filter_intervals(ht, [hl.parse_locus_interval("chrY")])
    ht_chry = ht_chry._filter_partitions(range(num_partitions))
    ht = ht_chr20.union(ht_chrx, ht_chry)

    return ht


def prepare_ht_for_validation(
    ht: hl.Table, freq_entries_to_remove: List[str], vcf_info_reorder: List[str]
) -> hl.Table:
    """
    Prepare HT for validity checks and export.

    :param ht: Release Hail Table
    :param freq_entries_to_remove: List of entries to remove from freq
    :param vcf_info_reorder: Order of VCF INFO fields
    :return: Hail Table prepared for validity checks and export
    """


def main(args):  # noqa: D103
    hl.init(
        log="/validate_and_export_vcf.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    overwrite = args.overwrite
    test = args.test
    resources = get_export_resources(overwrite, test)

    try:
        if args.validate_release_ht:
            logger.info("Validating release HT...")
            res = resources.validate_release
            res.check_resource_existence()

            ht = release_sites().ht()
            if test:
                logger.info("Filtering to test partitions...")
                ht = filter_to_test(ht)

            ht = prepare_ht_for_validation(ht)

            validated_ht = validate_release_t(ht)

            # Note: Checkpoint saves time for the final export by not needing to run the
            # VCF HT prep on each chromosome
            logger.info(
                "Checkpointing prepared VCF HT for validity checks and export..."
            )
            validated_ht = validated_ht.checkpoint(
                res.validated_ht.path, overwrite=True
            )
            validated_ht.describe()

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("validity_checks_and_export"))
