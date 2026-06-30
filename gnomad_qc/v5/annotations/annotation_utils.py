"""AoU-specific annotation utilities."""

import logging
from typing import Union

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("annotation_utils")
logger.setLevel(logging.INFO)


# AoU-specific adj annotation utilities.
# Adapted from https://github.com/broadinstitute/gatk/pull/8772/files.
# After discussion, we decided to use GQ 30 threshold for both haploid and diploid genotypes.
# See thread: https://atgu.slack.com/archives/CRA2TKTV0/p1762443074200349
def annotate_adj_no_dp(
    mt: hl.MatrixTable,
    adj_gq: int = 30,
    adj_ab: float = 0.2,
    gq_only: bool = False,
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria.

    Defaults are similar to gnomAD values, but GQ >= 20 changed to GQ >= 30 to make up for lack of DP filter
    (raised GQ threshold to 30 to match DP 10 threshold in gnomAD defaults).

    :param mt: Input MatrixTable.
    :param adj_gq: Minimum GQ. Default is 30.
    :param adj_ab: Minimum allele balance. Default is 0.2.
    :param gq_only: If True, annotate adj from GQ alone (``adj = GQ >= adj_gq``),
        bypassing the allele-balance test (which needs AD/LAD). Use for reference
        data: hom-ref reference blocks carry no AD/LAD, and the AB test is vacuously
        true for non-het calls, so GQ-only is exactly ``get_adj_expr`` restricted to
        reference blocks. Default is False.
    :return: MatrixTable with adj annotation.
    """
    if gq_only:
        return mt.annotate_entries(adj=mt.GQ >= adj_gq)
    if "LGT" in mt.entry and "LAD" in mt.entry:
        gt_expr = mt.LGT
        ad_expr = mt.LAD
    else:
        assert "GT" in mt.entry and "AD" in mt.entry
        gt_expr = mt.GT
        ad_expr = mt.AD
    return mt.annotate_entries(
        adj=get_adj_expr(gt_expr, mt.GQ, ad_expr, adj_gq, adj_ab)
    )


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 30,
    adj_ab: float = 0.2,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.

    Defaults are similar to gnomAD values, but GQ >= 20 changed to GQ >= 30 to make up for lack of DP filter.

    .. note::

        Assumes that the genotype expression is already adjusted for sex ploidy.

    :param gt_expr: Genotype expression.
    :param gq_expr: GQ expression.
    :param ad_expr: Allele depth expression.
    :param adj_gq: Minimum GQ. Default is 30.
    :param adj_ab: Minimum allele balance. Default is 0.2.
    :return: Expression for adj genotype annotation.
    """
    total_ad = hl.sum(ad_expr)
    return (gq_expr >= adj_gq) & (
        hl.case()
        .when(~gt_expr.is_het(), True)
        .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / total_ad >= adj_ab)
        .default(
            (ad_expr[gt_expr[0]] / total_ad >= adj_ab)
            & (ad_expr[gt_expr[1]] / total_ad >= adj_ab)
        )
    )
