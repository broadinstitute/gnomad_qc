"""Script to create release sites HT for v4 genomes."""
import hail as hl

from gnomad.utils.vep import filter_vep_transcript_csqs


# TODO: steps to create release sites HT
# input: VRS-annotated sites HT v3.1.4
# move SIFT and PolyPhen to insilico_predictors in vep_ht
# remove missing fields in vep_ht
# drop vep from site HT and merge with vep_ht
# remove old insilico predictors and replace with new ones (the left predictors)
# update dbSNP to dbSNP 156
# update freq, freq_meta, freq_index_dict by integrating updates in HGDP/TGP
# rerun inbreeding_coefficient with the callstats from the new freq
# replace oth with remaining in global fields
# update global fields


def get_sift_polyphen_from_vep(ht: hl.Table) -> hl.Table:
    """
    Get SIFT and PolyPhen scores from VEP 105 annotations.

    .. note::
     This function is to get a max of SIFT and PolyPhen scores from mane_select
     transcripts, otherwise get a max of SIFT and PolyPhen scores from canonical
     transcripts. It will also drop the SIFT and PolyPhen scores and predictions
     from the transcript_consequences struct.

    :param ht: VEP105 annotated Hail Table.
    :return: Hail Table with VEP105 annotations and  SIFT and PolyPhen scores
    extracted and stored in insilico_predictors struct.
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

    ht = ht.annotate(
        vep=ht.vep.annotate(
            transcript_consequences=ht.vep.transcript_consequences.map(
                lambda x: x.drop(
                    "sift_prediction",
                    "sift_score",
                    "polyphen_prediction",
                    "polyphen_score",
                )
            )
        )
    )
    ht = ht.annotate(
        insilico_predictors=hl.struct(
            sift_max=hl.or_else(hl.max(ht.sift_mane), hl.max(ht.sift_canonical)),
            polyphen_max=hl.or_else(
                hl.max(ht.polyphen_mane), hl.max(ht.polyphen_canonical)
            ),
        )
    )
    ht = ht.drop("sift_mane", "polyphen_mane", "sift_canonical", "polyphen_canonical")

    ht = ht.annotate_globals(
        tool_versions=hl.struct(sift_version="5.2.2", polyphen_version="2.2.2")
    )
    return ht


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
