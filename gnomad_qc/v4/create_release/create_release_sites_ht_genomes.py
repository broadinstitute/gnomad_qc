"""Script to create release sites HT for v4 genomes."""
import hail as hl


def remove_missing_vep_fields(vep_expr: hl.StructExpression) -> hl.StructExpression:
    """
    Remove fields from VEP 105 annotations that have been excluded in past releases or are missing in all rows.

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


def replace_oth_with_remaining(ht: hl.Table) -> hl.Table:
    """
    Replace 'oth' with 'remaining' in Global fields of a Table.

    .. note::
        This function is used to rename ancestry group fields to match v4 exomes.
        Which means that all "oth" will be changed to "remaining".
        2 Global fields are renamed: freq_meta (in the value of the key "pop") and freq_index_dict (in the keys).

    :param ht: release sites Table to be modified
    :return: release sites Table with 'oth' replaced with 'remaining'
    """
    # in freq_meta, type: array<dict<str, str>>, the value of key "pop" is
    # changed from 'oth' to 'remaining'.
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: x.replace("oth", "remaining"))
        )
    )

    # in freq_index_dict, type: dict<str, int>, the keys containing 'oth'
    # inside are changed to 'remaining'.
    oldkeys = hl.eval(ht.freq_index_dict.keys())
    newkeys = [s.replace("oth", "remaining") for s in oldkeys]
    vals = hl.eval(ht.freq_index_dict.values())
    freq_index_dict_new = {k: v for k, v in zip(newkeys, vals)}
    ht = ht.annotate_globals(freq_index_dict=freq_index_dict_new)
    return ht
