"""Script to run VQSR on an AS-annotated Sites VCF."""

import argparse
import json
import logging
from typing import Dict, List, Optional

import hail as hl
import hailtop.batch as hb
from hailtop.batch.job import Job

from gnomad_qc.v5.resources.annotations import get_true_positive_vcf_path, info_vcf_path
from gnomad_qc.v5.resources.basics import calling_intervals
from gnomad_qc.v5.resources.sample_qc import interval_qc_pass
from gnomad_qc.v5.resources.variant_qc import VQSR_FEATURES, get_variant_qc_result
from gnomad_qc.v5.variant_qc.import_variant_qc_vcf import (
    import_variant_qc_vcf as import_vqsr,
)

logger = logging.getLogger(__file__)
logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger.setLevel(logging.INFO)


SNP_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.8,
    99.6,
    99.5,
    99.4,
    99.3,
    99.0,
    98.0,
    97.0,
    90.0,
]
INDEL_RECALIBRATION_TRANCHE_VALUES = [
    100.0,
    99.95,
    99.9,
    99.5,
    99.0,
    97.0,
    96.0,
    95.0,
    94.0,
    93.5,
    93.0,
    92.0,
    91.0,
    90.0,
]
