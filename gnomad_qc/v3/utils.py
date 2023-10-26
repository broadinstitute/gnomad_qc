# noqa: D100

import hail as hl


def hom_alt_depletion_fix(
    gt_expr: hl.expr.CallExpression,
    het_non_ref_expr: hl.expr.BooleanExpression,
    af_expr: hl.expr.Float64Expression,
    ab_expr: hl.expr.Float64Expression,
    af_cutoff: float = 0.01,
    ab_cutoff: float = 0.9,
    use_v3_1_correction: bool = False,
) -> hl.expr.CallExpression:
    """
    Get expression for genotypes with temporary fix for the depletion of homozygous alternate genotypes.

    More details about the problem can be found on the gnomAD blog:
    https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#tweaks-and-updates

    :param gt_expr: Genotype expression to be adjusted.
    :param het_non_ref_expr: Expression indicating whether the original genotype (pre-
        split multi) is het non ref.
    :param af_expr: Allele frequency expression to determine which variants need the
        hom alt fix.
    :param ab_expr: Allele balance expression to determine which genotypes need the
        hom alt fix.
    :param af_cutoff: Allele frequency cutoff for variants that need the hom alt fix.
        Default is 0.01.
    :param ab_cutoff: Allele balance cutoff to determine which genotypes need the hom
        alt fix. Default is 0.9.
    :param use_v3_1_correction: Use the version of the correction that was used for
        the v3.1.2 release. This version was missing the `missing_false=True` in the
        if_else statement. Default is False.
    :return: Expression for genotypes adjusted for the hom alt depletion fix.
    """
    return hl.if_else(
        (af_expr > af_cutoff) & gt_expr.is_het()
        # Skip adjusting genotypes if sample originally had a het nonref genotype.
        & ~het_non_ref_expr & (ab_expr > ab_cutoff),
        hl.call(1, 1),
        gt_expr.GT,
        missing_false=False if use_v3_1_correction else True,
    )
