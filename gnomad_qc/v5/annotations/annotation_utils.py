"""AoU-specific annotation utilities."""

import logging
from typing import Union

import hail as hl
from gnomad.utils.annotations import get_adj_expr as get_gnomad_adj_expr

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("annotation_utils")
logger.setLevel(logging.INFO)


# AoU-specific adj annotation utilities.
# Adapted from https://github.com/broadinstitute/gatk/pull/8772/files.
#
# The AoU VDS has no DP entry field, so DP is approximated as sum(LAD). Hom-ref calls
# (reference blocks, which carry GQ but no LAD) are filtered on GQ only. All other
# calls use the usual gnomAD cutoffs: GQ, DP (as sum(LAD)), and AB for het calls.
# After discussion, we decided to use a GQ 20 threshold for both haploid and diploid
# genotypes. See thread: https://atgu.slack.com/archives/CRA2TKTV0/p1787333884174549
def annotate_adj_no_dp(
    mt: hl.MatrixTable,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.MatrixTable:
    """
    Annotate genotypes with adj criteria.

    Defaults correspond to gnomAD values. DP is approximated as the sum of the allele
    depths and is only applied to non-hom-ref calls; hom-ref calls are filtered on GQ
    only. See the module comment for details.

    :param mt: Input MatrixTable.
    :param adj_gq: Minimum GQ. Default is 20.
    :param adj_dp: Minimum DP (sum of allele depths) for non-hom-ref calls. Default is
        10.
    :param adj_ab: Minimum allele balance for het calls. Default is 0.2.
    :param haploid_adj_dp: Minimum DP (sum of allele depths) for haploid non-ref
        calls. Default is 5.
    :return: MatrixTable with adj annotation.
    """
    if "LGT" in mt.entry and "LAD" in mt.entry:
        gt_expr = mt.LGT
        ad_expr = mt.LAD
    else:
        assert "GT" in mt.entry and "AD" in mt.entry
        gt_expr = mt.GT
        ad_expr = mt.AD
    return mt.annotate_entries(
        adj=get_adj_expr(
            gt_expr, mt.GQ, ad_expr, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    )


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: Union[hl.expr.Int32Expression, hl.expr.Int64Expression],
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.

    Defaults correspond to gnomAD values. Hom-ref calls are filtered on GQ only. All
    other calls use the standard gnomAD adj criteria (GQ, DP, and AB for het calls)
    with DP approximated as the sum of `ad_expr`.

    .. note::

        Assumes that the genotype expression is already adjusted for sex ploidy.

    :param gt_expr: Genotype expression.
    :param gq_expr: GQ expression.
    :param ad_expr: Allele depth expression.
    :param adj_gq: Minimum GQ. Default is 20.
    :param adj_dp: Minimum DP (sum of allele depths) for non-hom-ref calls. Default is
        10.
    :param adj_ab: Minimum allele balance for het calls. Default is 0.2.
    :param haploid_adj_dp: Minimum DP (sum of allele depths) for haploid non-ref
        calls. Default is 5.
    :return: Expression for adj genotype annotation.
    """
    return hl.if_else(
        gt_expr.is_hom_ref(),
        gq_expr >= adj_gq,
        get_gnomad_adj_expr(
            gt_expr,
            gq_expr,
            hl.sum(ad_expr),
            ad_expr,
            adj_gq=adj_gq,
            adj_dp=adj_dp,
            adj_ab=adj_ab,
            haploid_adj_dp=haploid_adj_dp,
        ),
    )
