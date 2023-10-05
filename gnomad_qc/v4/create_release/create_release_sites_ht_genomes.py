"""Script to create release sites HT for v4 genomes."""
from typing import Tuple

import hail as hl
from gnomad.utils.vep import filter_vep_transcript_csqs


def get_sift_polyphen_from_vep(
    ht: hl.Table,
) -> Tuple[hl.ArrayExpression, hl.ArrayExpression]:
    """
    Get the max SIFT and PolyPhen scores from VEP 105 annotations.

     This retrieves the max of SIFT and PolyPhen scores for a variant's MANE Select
     transcript or, if MANE Select does not exist, canonical transcript. This also
     drops SIFT and PolyPhen scores and predictions from VEP's transcript_consequences
     struct.

    :param ht: VEP 105 annotated Hail Table.
    :return: Tuple of max SIFT and PolyPhen scores.
    """
    mane = filter_vep_transcript_csqs(
        ht, synonymous=False, canonical=False, mane_select=True
    )
    canonical = filter_vep_transcript_csqs(ht, synonymous=False, canonical=True)
    ht = ht.annotate(
        sift_mane=mane[ht.key].vep.transcript_consequences.sift_score,
        polyphen_mane=mane[ht.key].vep.transcript_consequences.polyphen_score,
        sift_canonical=canonical[ht.key].vep.transcript_consequences.sift_score,
        polyphen_canonical=canonical[ht.key].vep.transcript_consequences.polyphen_score,
    )

    sift_max = hl.or_else(hl.max(ht.sift_mane), hl.max(ht.sift_canonical))
    polyphen_max = hl.or_else(hl.max(ht.polyphen_mane), hl.max(ht.polyphen_canonical))

    return sift_max, polyphen_max


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
        This function renames v3 ancestry groups to match v4 exomes' ancestry groups.
        The value of the key "pop" within `freq_meta` and the keys of `freq_index_dict`
        that contain "oth" are replaced with "remaining".

    :param ht: release sites Table to be modified.
    :return: release sites Table with 'oth' replaced with 'remaining'.
    """
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: x.replace("oth", "remaining"))
        ),
        freq_index_dict=hl.dict(
            hl.zip(
                ht.freq_index_dict.keys().map(lambda k: k.replace("oth", "remaining")),
                ht.freq_index_dict.values(),
            )
        ),
    )
    return ht
