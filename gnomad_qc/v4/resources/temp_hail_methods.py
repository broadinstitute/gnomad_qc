"""THIS IS A TEMPORARY FILE TO BE USED UNTIL THESE HAIL METHODS ARE RELEASED. 
CODE WAS COPIED DIRECTLY FROM A PR IN GITHUB BROWSER ON 02/20 (iirc) AND SHOULD NOT BE RAN ON RELEASE."""
import hail as hl
from typing import Optional, Sequence
from hail.expr.expressions import Expression
from hail.expr.functions import _allele_types, _num_allele_type
from hail.utils.misc import divide_null, new_temp_file, wrap_to_list


def vmt_sample_qc(
    *,
    global_gt: "Expression",
    gq: "Expression",
    variant_ac: "Expression",
    variant_atypes: "Expression",
    dp: Optional["Expression"] = None,
    gq_bins: "Sequence[int]" = (0, 20, 60),
    dp_bins: "Sequence[int]" = (0, 1, 10, 20, 30),
    dp_field: Optional[str] = None,
) -> "Expression":
    """"""
    allele_types = _allele_types[:]
    allele_types.extend(["Transition", "Transversion"])
    allele_enum = {i: v for i, v in enumerate(allele_types)}
    allele_ints = {v: k for k, v in allele_enum.items()}

    bound_exprs = {}

    bound_exprs["n_het"] = hl.agg.count_where(global_gt.is_het())
    bound_exprs["n_hom_var"] = hl.agg.count_where(global_gt.is_hom_var())
    bound_exprs["n_singleton"] = hl.agg.sum(  # returns the sum per sample
        hl.rbind(
            global_gt,  # for each GT entry in GT
            lambda global_gt: hl.sum(  # taking the sum of:
                # either 0 or 1
                hl.range(
                    0, global_gt.ploidy
                ).map(  # [0,1].map() # -> [False, True] or [False, False]
                    # for every allele (ref or alt, 0,1):
                    lambda i: hl.rbind(
                        global_gt[i], lambda gti: (gti != 0) & (variant_ac[gti] == 1)
                    )
                )
            ),
        )
    )
    bound_exprs["n_singleton_ti"] = hl.agg.sum(
        hl.rbind(
            global_gt,
            lambda global_gt: hl.sum(
                hl.range(0, global_gt.ploidy).map(
                    lambda i: hl.rbind(
                        global_gt[i],
                        lambda gti: (gti != 0)
                        & (variant_ac[gti] == 1)
                        & (variant_atypes[gti - 1] == allele_ints["Transition"]),
                    )
                )
            ),
        )
    )
    bound_exprs["n_singleton_tv"] = hl.agg.sum(
        hl.rbind(
            global_gt,
            lambda global_gt: hl.sum(
                hl.range(0, global_gt.ploidy).map(
                    lambda i: hl.rbind(
                        global_gt[i],
                        lambda gti: (gti != 0)
                        & (variant_ac[gti] == 1)
                        & (variant_atypes[gti - 1] == allele_ints["Transversion"]),
                    )
                )
            ),
        )
    )

    bound_exprs["allele_type_counts"] = hl.agg.explode(
        lambda allele_type: hl.tuple(
            hl.agg.count_where(allele_type == i) for i in range(len(allele_ints))
        ),
        (
            hl.range(0, global_gt.ploidy)
            .map(lambda i: global_gt[i])
            .filter(lambda allele_idx: allele_idx > 0)
            .map(lambda allele_idx: variant_atypes[allele_idx - 1])
        ),
    )

    dp_exprs = {}
    if dp is not None:
        dp_exprs["dp"] = hl.tuple(hl.agg.count_where(dp >= x) for x in dp_bins)

    gq_dp_exprs = hl.struct(
        **{"gq": hl.tuple(hl.agg.count_where(gq >= x) for x in gq_bins)}, **dp_exprs
    )

    return hl.rbind(
        hl.struct(**bound_exprs),
        lambda x: hl.rbind(
            hl.struct(
                **{
                    "gq_dp_exprs": gq_dp_exprs,
                    "n_het": x.n_het,
                    "n_hom_var": x.n_hom_var,
                    "n_non_ref": x.n_het + x.n_hom_var,
                    "n_singleton": x.n_singleton,
                    "n_singleton_ti": x.n_singleton_ti,
                    "n_singleton_tv": x.n_singleton_tv,
                    "n_snp": (
                        x.allele_type_counts[allele_ints["Transition"]]
                        + x.allele_type_counts[allele_ints["Transversion"]]
                    ),
                    "n_insertion": x.allele_type_counts[allele_ints["Insertion"]],
                    "n_deletion": x.allele_type_counts[allele_ints["Deletion"]],
                    "n_transition": x.allele_type_counts[allele_ints["Transition"]],
                    "n_transversion": x.allele_type_counts[allele_ints["Transversion"]],
                    "n_star": x.allele_type_counts[allele_ints["Star"]],
                }
            ),
            lambda s: s.annotate(
                r_ti_tv=divide_null(hl.float64(s.n_transition), s.n_transversion),
                r_ti_tv_singleton=divide_null(
                    hl.float64(s.n_singleton_ti), s.n_singleton_tv
                ),
                r_het_hom_var=divide_null(hl.float64(s.n_het), s.n_hom_var),
                r_insertion_deletion=divide_null(
                    hl.float64(s.n_insertion), s.n_deletion
                ),
            ),
        ),
    )


def vmt_sample_qc_variant_annotations(
    *,
    global_gt: "Expression",
    alleles: "Expression",
) -> "Expression":
    """"""

    allele_types = _allele_types[:]
    allele_types.extend(["Transition", "Transversion"])
    allele_enum = {i: v for i, v in enumerate(allele_types)}
    allele_ints = {v: k for k, v in allele_enum.items()}

    def allele_type(ref, alt):
        return hl.bind(
            lambda at: hl.if_else(
                at == allele_ints["SNP"],
                hl.if_else(
                    hl.is_transition(ref, alt),
                    allele_ints["Transition"],
                    allele_ints["Transversion"],
                ),
                at,
            ),
            _num_allele_type(ref, alt),
        )

    return hl.struct(
        variant_ac=hl.agg.call_stats(global_gt, alleles).AC,
        variant_atypes=alleles[1:].map(lambda alt: allele_type(alleles[0], alt)),
    )
