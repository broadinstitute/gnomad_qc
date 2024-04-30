"""Script containing version and release constants."""
DATA_TYPES = ["exomes", "genomes"]
RELEASE_DATA_TYPES = DATA_TYPES + ["joint"]

VERSIONS = ["4.0", "4.1"]
CURRENT_VERSION = "4.1"

# NOTE: There were no changes to the genomes frequency data for v4.1.
FREQ_VERSIONS = {"exomes": ["4.0", "4.1"], "genomes": ["4.0"], "joint": ["4.0", "4.1"]}
CURRENT_FREQ_VERSION = {"exomes": "4.1", "genomes": "4.0", "joint": "4.1"}

# NOTE: There were no changes to raw data for v4.1.
RAW_VERSIONS = ["4.0"]
CURRENT_RAW_VERSION = "4.0"

# NOTE: There were no changes to project metadata for v4.1.
PROJECT_META_VERSIONS = ["4.0"]
CURRENT_PROJECT_META_VERSION = "4.0"

# NOTE: There were no changes to sample QC for v4.1.
SAMPLE_QC_VERSIONS = ["4.0"]
CURRENT_SAMPLE_QC_VERSION = "4.0"

# NOTE: There were no changes to most variant annotation files for v4.1.
ANNOTATION_VERSIONS = ["4.0"]
CURRENT_ANNOTATION_VERSION = "4.0"

# NOTE: There were no changes to most variant QC files for v4.1.
VARIANT_QC_VERSIONS = ["4.0"]
CURRENT_VARIANT_QC_VERSION = "4.0"

# NOTE: There were fixes put in to the exomes VQSR loading for v4.1.
VARIANT_QC_RESULT_VERSIONS = {"exomes": ["4.0", "4.1"], "genomes": ["4.0"]}
CURRENT_VARIANT_QC_RESULT_VERSION = {"exomes": "4.1", "genomes": "4.0"}

RELEASES = ["4.0", "4.1"]
CURRENT_RELEASE = "4.1"

COVERAGE_RELEASES = {"exomes": ["4.0"], "genomes": ["3.0"]}
CURRENT_COVERAGE_RELEASE = {"exomes": "4.0", "genomes": "3.0"}

ALL_SITES_AN_RELEASES = {"exomes": ["4.1"], "genomes": ["4.1"]}
CURRENT_ALL_SITES_AN_RELEASE = {"exomes": "4.1", "genomes": "4.1"}

HGDP_TGP_RELEASES = ["4.0"]
CURRENT_HGDP_TGP_RELEASE = "4.0"
