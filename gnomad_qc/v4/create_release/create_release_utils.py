"""Common utils for creating gnomAD v4 exomes and genomes releases."""

import hail as hl

from gnomad_qc.v4.resources.constants import DEFAULT_VEP_VERSION

# Tool versions used in the release that were not on the source tables
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
    "115": "v1.2",
}
# Backward compatibility: default to DEFAULT_VEP_VERSION.
GENCODE_VERSION = GENCODE_VERSIONS[DEFAULT_VEP_VERSION]
MANE_SELECT_VERSION = MANE_SELECT_VERSIONS[DEFAULT_VEP_VERSION]

# NOTE: VEP 115 annotations were added in v4.1.1.
VEP_VERSIONS_TO_ADD = {"4.1.1": ["115"]}


def remove_missing_vep_fields(vep_expr: hl.StructExpression) -> hl.StructExpression:
    """
    Remove fields from VEP 105 annotations that have been excluded in past releases or are missing in all rows.

    See https://github.com/broadinstitute/gnomad_production/issues/936 for more
    information.

    :param vep_expr: StructExpression containing VEP 105 annotations.
    :return: StructExpression containing VEP 105 annotations with missing fields
        removed.
    """
    vep_expr = vep_expr.drop("colocated_variants", "context")

    vep_expr = vep_expr.annotate(
        transcript_consequences=vep_expr.transcript_consequences.map(
            lambda x: x.drop("minimised", "swissprot", "trembl", "uniparc")
        )
    )

    for consequence in [
        "intergenic_consequences",
        "motif_feature_consequences",
        "regulatory_feature_consequences",
    ]:
        vep_expr = vep_expr.annotate(
            **{consequence: vep_expr[consequence].map(lambda x: x.drop("minimised"))}
        )
    return vep_expr
