import hail as hl


def hom_alt_depletion_fix(
    mt: hl.MatrixTable,
    het_non_ref_expr: hl.expr.BooleanExpression,
    af_expr: hl.expr.Float64Expression,
    af_cutoff: float = 0.01,
    ab_cutoff: float = 0.9,
) -> hl.MatrixTable:
    """
    Adjust MT genotypes with temporary fix for the depletion of homozygous alternate genotypes.
    
    More details about the problem can be found on the gnomAD blog:
    https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#tweaks-and-updates
    
    :param mt: Input MT that needs hom alt genotype fix
    :param het_non_ref_expr: Expression indicating whether the original genotype (pre split multi) is het non ref
    :param af_expr: Allele frequency expression to determine which variants need the hom alt fix
    :param af_cutoff: Allele frequency cutoff for variants that need the hom alt fix. Default is 0.01
    :param ab_cutoff: Allele balance cutoff to determine which genotypes need the hom alt fix. Default is 0.9
    :return: MatrixTable with genotypes adjusted for the hom alt depletion fix
    """
    return mt.annotate_entries(
        GT=hl.if_else(
            mt.GT.is_het()
            # Skip adjusting genotypes if sample originally had a het nonref genotype
            & ~het_non_ref_expr
            & (af_expr > af_cutoff)
            & (mt.AD[1] / mt.DP > ab_cutoff),
            hl.call(1, 1),
            mt.GT,
        )
    )
