"""Temporary script to set XX samples chrY faf stats to missing in final freq HT."""

import hail as hl
from Typing import Union

from gnomad_qc.v4.resources.annotations import get_freq


def missing_faf_expr() -> hl.expr.StructExpression:
    """
    Create a missing faf stats struct for insertion into faf annotation arrays when data is missing.

    :return: Hail Struct with missing values for each callstats element
    """
    return hl.struct(
        faf95=hl.missing(hl.tfloat64),
        faf99=hl.missing(hl.tfloat64),
    )


def set_female_y_metrics_to_na_expr(
    t: Union[hl.Table, hl.MatrixTable],
    faf_expr: Union[hl.expr.ArrayExpression, str] = "faf",
    faf_meta_expr: Union[hl.expr.ArrayExpression, str] = "faf_meta",
    faf_index_dict_expr: Union[hl.expr.DictExpression, str] = "faf_index_dict",
) -> hl.expr.ArrayExpression:
    """
    Set Y-variant faf stats for female-specific metrics to missing structs.

    :param t: Table or MatrixTable for which to adjust female metrics.
    :param faf_expr: Array expression or string annotation name for the faf
        array. Default is "faf".
    :param faf_meta_expr: Array expression or string annotation name for the faf
        metadata. Default is "faf_meta".
    :param freq_index_dict_expr: Dict expression or string annotation name for the
        faf metadata index dictionary. Default is "faf_index_dict".
    :return: Hail array expression to set female Y-variant metrics to missing values.
    """
    if isinstance(faf_expr, str):
        faf_expr = t[faf_expr]
    if isinstance(faf_meta_expr, str):
        faf_meta_expr = t[faf_meta_expr]
    if isinstance(faf_index_dict_expr, str):
        faf_index_dict_expr = t[faf_index_dict_expr]

    female_idx = hl.map(
        lambda x: faf_index_dict_expr[x],
        hl.filter(lambda x: x.contains("XX"), faf_index_dict_expr.keys()),
    )
    faf_idx_range = hl.range(hl.len(faf_meta_expr))

    new_faf_expr = hl.if_else(
        (t.locus.in_y_nonpar() | t.locus.in_y_par()),
        hl.map(
            lambda x: hl.if_else(
                female_idx.contains(x), missing_faf_expr(), faf_expr[x]
            ),
            faf_idx_range,
        ),
        faf_expr,
    )

    return new_faf_expr


hl.init(
    log="/fix_chrY_missing_in_freq.log",
    default_reference="GRCh38",
    tmp_dir="gs://gnomad-tmp-30day",
)

ht = get_freq().ht()
ht = ht.checkpoint(hl.utils.new_temp_file("exomes_freq_4_1", "ht"))

ht = ht.annotate(faf=set_female_y_metrics_to_na_expr(ht))
ht.write(get_freq().path, overwrite=True)
