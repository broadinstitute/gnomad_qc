# noqa: D100

import argparse
import logging
import pickle
from typing import Dict, List, Optional, Set, Union

import hail as hl
from gnomad.assessment.validity_checks import validate_release_t, vcf_field_check
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
from gnomad_qc.v4.resources.basics import get_checkpoint_path, qc_temp_prefix
from gnomad_qc.v4.resources.release import (
    append_to_vcf_header_path,
    release_header_path,
    release_sites,
    release_vcf_path,
)

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
        output_resources={
            "validated_release_ht": get_validated_release_ht(
                test=test, overwrite=overwrite
            )
        },
    )
    export_pipeline.add_steps(
        {
            "validate_release_ht": validate_release,
        }
    )
    return export_pipeline
