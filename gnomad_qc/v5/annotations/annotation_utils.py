"""AoU-specific annotation utilities."""

import logging

import hail as hl
from gnomad.utils.annotations import get_adj_expr as get_gnomad_adj_expr

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("annotation_utils")
logger.setLevel(logging.INFO)


# AoU-specific adj annotation utilities.
# Adapted from https://github.com/broadinstitute/gatk/pull/8772/files.
#
# The AoU VDS has no DP entry field, so DP is approximated as sum(LAD). Hom-ref calls
# with no allele depths (reference blocks, which carry GQ but no LAD) are filtered on
# GQ only. All other calls, including hom-ref calls that do carry allele depths (e.g.
# calls downcoded to hom-ref by split_multi), use the usual gnomAD cutoffs: GQ, DP (as
# sum(LAD)), and AB for het calls.
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
    depths. Hom-ref calls with no allele depths (reference blocks) are filtered on GQ
    only; every other call uses the gnomAD cutoffs. See the module comment for details.

    Accepts three entry layouts:

        - ``LGT`` and ``LAD`` (VDS variant data): full adj criteria.
        - ``GT`` and ``AD`` (split or dense data): full adj criteria.
        - ``GT`` without allele depths (VDS reference data, e.g. the AoU reference
          blocks with only ``GT``, ``GQ``, ``END`` and ``LEN``): every call must be
          hom-ref and adj is ``GQ >= adj_gq``. A non-hom-ref call raises at runtime
          rather than being silently marked non-adj.

    ``GQ`` is required in all three cases.

    :param mt: Input MatrixTable.
    :param adj_gq: Minimum GQ. Default is 20.
    :param adj_dp: Minimum DP (sum of allele depths) for calls with allele depths.
        Default is 10.
    :param adj_ab: Minimum allele balance for het calls. Default is 0.2.
    :param haploid_adj_dp: Minimum DP (sum of allele depths) for haploid calls with
        allele depths. Default is 5.
    :return: MatrixTable with adj annotation.
    """
    entry_fields = set(mt.entry)
    if "GQ" not in entry_fields:
        raise ValueError("annotate_adj_no_dp requires a 'GQ' entry field.")

    if "LGT" in entry_fields and "LAD" in entry_fields:
        adj_expr = get_adj_expr(
            mt.LGT, mt.GQ, mt.LAD, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    elif "GT" in entry_fields and "AD" in entry_fields:
        adj_expr = get_adj_expr(
            mt.GT, mt.GQ, mt.AD, adj_gq, adj_dp, adj_ab, haploid_adj_dp
        )
    elif "GT" in entry_fields:
        # Reference data: hom-ref blocks that carry GQ but no allele depths, so
        # adj is GQ only. A missing genotype (e.g. after a sex-ploidy adjustment)
        # gets a missing adj, matching get_adj_expr. Error on a non-hom-ref call
        # rather than silently marking it non-adj.
        adj_expr = (
            hl.case()
            .when(hl.is_missing(mt.GT), hl.missing(hl.tbool))
            .when(mt.GT.is_hom_ref(), mt.GQ >= adj_gq)
            .or_error(
                "Found a non-hom-ref genotype in a MatrixTable with no AD or LAD "
                "entry field. Allele depths are required to compute adj for "
                "non-hom-ref calls."
            )
        )
    else:
        raise ValueError(
            "annotate_adj_no_dp requires 'LGT' and 'LAD', 'GT' and 'AD', or 'GT' "
            "alone (hom-ref reference data) in the entry fields."
        )

    return mt.annotate_entries(adj=adj_expr)


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: hl.expr.Int32Expression | hl.expr.Int64Expression,
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.expr.BooleanExpression:
    """
    Get adj genotype annotation.

    Defaults correspond to gnomAD values. Hom-ref calls with a missing `ad_expr`
    (reference blocks) are filtered on GQ only. All other calls, including hom-ref
    calls with defined allele depths, use the standard gnomAD adj criteria (GQ, DP,
    and AB for het calls) with DP approximated as the sum of `ad_expr`.

    .. note::

        Assumes that the genotype expression is already adjusted for sex ploidy.

    :param gt_expr: Genotype expression.
    :param gq_expr: GQ expression.
    :param ad_expr: Allele depth expression.
    :param adj_gq: Minimum GQ. Default is 20.
    :param adj_dp: Minimum DP (sum of allele depths) for calls with allele depths.
        Default is 10.
    :param adj_ab: Minimum allele balance for het calls. Default is 0.2.
    :param haploid_adj_dp: Minimum DP (sum of allele depths) for haploid calls with
        allele depths. Default is 5.
    :return: Expression for adj genotype annotation.
    """
    return hl.if_else(
        gt_expr.is_hom_ref() & hl.is_missing(ad_expr),
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
