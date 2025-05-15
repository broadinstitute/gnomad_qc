"""Script containing version and release constants."""

VERSIONS = ["5.0"]
CURRENT_VERSION = VERSIONS[-1]
AOU_VERSIONS = ["8"]
CURRENT_AOU_VERSION = AOU_VERSIONS[-1]
CURRENT_PROJECT_META_VERSION = "5.0"

SAMPLE_QC_VERSIONS = ["5.0"]
CURRENT_SAMPLE_QC_VERSION = "5.0"

# Storage constants.
WORKSPACE_BUCKET = "fc-secure-b25d1307-7763-48b8-8045-fcae9caadfa1"
GNOMAD_BUCKET = "gnomad"
GNOMAD_TMP_BUCKET = "gnomad-tmp"
AOU_BUCKET = "gs://fc-aou-datasets-controlled/v8"
AOU_WGS_BUCKET = f"{AOU_BUCKET}/wgs/short_read/snpindel"
AOU_WGS_AUX_BUCKET = f"{AOU_WGS_BUCKET}/aux"
AOU_LOW_QUALITY_PATH = f"{AOU_BUCKET}/known_issues/wgs_v8_known_issue_1.txt"
AOU_GENOMIC_METRICS_PATH = f"{AOU_WGS_AUX_BUCKET}/qc/genomic_metrics.tsv"
