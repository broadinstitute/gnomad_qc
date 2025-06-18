"""Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels or spike ins."""

import argparse
import json
import logging
import pickle
from typing import Any, Dict, List, Optional, Tuple

import hail as hl
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.sample_qc import hgdp_tgp_pop_outliers
from gnomad_qc.v4.sample_qc.assign_ancestry import V3_SPIKE_PROJECTS, V4_POP_SPIKE_DICT
from gnomad_qc.v5.resources.basics import get_checkpoint_path
from gnomad_qc.v5.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    get_joint_qc,
    get_pop_ht,
    get_pop_pr_ht,
    # joint_qc_meta,
    per_pop_min_rf_probs_json_path,
    pop_rf_path,
    # related_samples_to_drop,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("ancestry_assignment")
logger.setLevel(logging.INFO)
