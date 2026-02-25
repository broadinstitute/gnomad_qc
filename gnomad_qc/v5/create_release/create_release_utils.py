"""Common utils for creating gnomAD v5 genome release."""

import hail as hl

from gnomad_qc.v5.resources.constants import DEFAULT_VEP_VERSION

# Tool versions used in the release that were not on the source tables.
DBSNP_VERSION = "b156"
SIFT_VERSION = "5.2.2"
POLYPHEN_VERSION = "2.2.2"
VRS_SCHEMA_VERSION = "2.0.1"
VRS_PYTHON_VERSION = "2.2.0"
SEQREPO_VERSION = "2024-12-20"

# GENCODE versions by VEP version.
GENCODE_VERSIONS = {
    "105": "Release 39",
    "115": "Release 49",
}
# MANE Select versions by VEP version.
MANE_SELECT_VERSIONS = {
    "105": "v0.95",
    "115": "v1.4",
}
# Backward compatibility: default to DEFAULT_VEP_VERSION.
GENCODE_VERSION = GENCODE_VERSIONS[DEFAULT_VEP_VERSION]
MANE_SELECT_VERSION = MANE_SELECT_VERSIONS[DEFAULT_VEP_VERSION]

# NOTE: VEP 115 annotations were added in v4.1.1.
VEP_VERSIONS_TO_ADD = {"5.0": ["115"]}
