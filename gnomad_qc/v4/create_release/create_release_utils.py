"""Common utils for creating gnomAD v4 exomes and genomes releases."""
import hail as hl

# Tool versions used in the release that were not on the source tables
GENCODE_VERSION = "Release 39"
MANE_SELECT_VERSION = "v0.95"
DBSNP_VERSION = "b156"
SIFT_VERSION = "5.2.2"
POLYPHEN_VERSION = "2.2.2"
VRS_SCHEMA_VERSION = "1.3.0"
VRS_PYTHON_VERSION = "0.8.4"
SEQREPO_VERSION = "2018-11-26"


def remove_missing_vep_fields(vep_expr: hl.StructExpression) -> hl.StructExpression:
    """
    Remove fields from VEP 105 annotations that have been excluded in past releases or are missing in all rows.

    See https://github.com/broadinstitute/gnomad_production/issues/936 for more information.

    :param vep_expr: StructExpression containing VEP 105 annotations.
    :return: StructExpression containing VEP 105 annotations with missing fields removed.
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
