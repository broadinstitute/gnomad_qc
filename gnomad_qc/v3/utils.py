import hail as hl


def hom_alt_depletion_fix(
        mt: hl.MatrixTable,
        af_expr: hl.expr.Float32Expression,
        af_cutoff: float = 0.01,
        ab_cutoff: float = 0.9,
) -> hl.MatrixTable:
    """
    Adjust MT genotypes with temporary fix for the depletion of homozygous alternate genotypes.

    :param mt: Input MT to add hom alt depletion genotype fix to
    :param af_expr: Allele frequency expression to determine variants that need the hom alt fix
    :param af_cutoff: Allele frequency cutoff for variants that need the hom alt fix
    :param ab_cutoff: Allele balance cutoff to determine which genotypes need the hom alt fix
    :return: MatrixTable with genotypes adjusted for the hom alt depletion fix
    """
    return mt.annotate_entries(
        GT=hl.cond(
            (af_expr > af_cutoff) & mt.GT.is_het() & ((mt.AD[1] / mt.DP) > ab_cutoff),
            hl.call(1, 1),
            mt.GT,
        )
    )
