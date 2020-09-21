from gnomad_qc.v3.resources import score_ranking_path, get_filtering_model, freq
from gnomad.utils.annotations import bi_allelic_site_inbreeding_expr
from gnomad.utils.filtering import add_filters_expr
import argparse
import hail as hl
import logging

logger = logging.getLogger("create_filters")


def main(args):

    # Get snv and indel cutoffs based on bin cutoffs
    binned_ht = hl.read_table(score_ranking_path(args.model_id, binned=True))
    snp_rf_cutoff, indel_rf_cutoff = binned_ht.aggregate(
        [
            hl.agg.filter(
                binned_ht.snv & (binned_ht.bin == args.snv_bin_cutoff),
                hl.agg.min(binned_ht.min_score),
            ),
            hl.agg.filter(
                ~binned_ht.snv & (binned_ht.bin == args.indel_bin_cutoff),
                hl.agg.min(binned_ht.min_score),
            ),
        ]
    )

    # Add filters to the HT
    ht = hl.read_table(score_ranking_path(args.model_id, binned=False))
    model_ht = get_filtering_model(args.model_id, split=True, finalized=False).ht()
    ht = ht.annotate(info=model_ht[ht.key].info)

    ht = ht.annotate_globals(
        filtering_model=hl.struct(
            model_id=args.model_id,
            model_name=args.model_name,
            score_name=args.score_name,
            snv_cutoff=hl.struct(bin=args.snv_bin_cutoff, min_score=snp_rf_cutoff),
            indel_cutoff=hl.struct(
                bin=args.indel_bin_cutoff, min_score=indel_rf_cutoff
            ),
        )
    )

    n_missing = ht.aggregate(hl.agg.count_where(hl.is_missing((ht.score))))
    if n_missing:
        logger.warn(
            f"Found {n_missing} sites where {args.model_id} score is missing. These sites will be filtered."
        )

    indexed_freq_ht = freq.ht()[ht.key]
    ht = ht.annotate(
        # InbreedingCoeff=bi_allelic_site_inbreeding_from_call_stats_expr(indexed_freq_ht.freq[1]),
        AC0=indexed_freq_ht.freq[0].AC
        < 1
    )

    filters = {
        args.model_name: hl.or_else(
            hl.cond(
                hl.is_snp(ht.alleles[0], ht.alleles[1]),
                ht.score < ht.filtering_model.snv_cutoff.min_score,
                ht.score < ht.filtering_model.indel_cutoff.min_score,
            ),
            True,
        ),
        "InbreedingCoeff": ht.InbreedingCoeff < args.inbreeding_coeff_threshold,
        "AC0": ht.AC0,
    }

    ht = ht.annotate(filters=add_filters_expr(filters))

    ht = ht.select(
        "filters",
        info=hl.struct(
            InbreedingCoeff=ht.InbreedingCoeff,
            **{args.score_name: ht.score},
            **{field: ht.info[field] for field in args.info_fields},
        ),
    )

    ht.describe()

    # Fix annotations for release
    # annotations_expr = {
    #     'tp': hl.or_else(ht.tp, False),
    #     'transmitted_singleton': hl.or_missing(freq_ht[ht.key].freq[1].AC[1] == 1, ht.transmitted_singleton)
    #     'rf_probability': ht.rf_probability["TP"]
    # }
    # if 'feature_imputed' in ht.row:
    #     annotations_expr.update(
    #         {x: hl.or_missing(~ht.feature_imputed[x], ht[x]) for x in [f for f in ht.row.feature_imputed]}
    #     )
    #
    # ht = ht.transmute(
    #     **annotations_expr
    # )

    # This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
    # ht = ht.drop(ht.rf_prediction)
    ht.write(
        get_filtering_model(args.model_id, split=True, finalized=True).path,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument('--write_filters_ht', help='Creates a table with VQSR rank', action='store_true')
    parser.add_argument("--model_id", help="Filtering model ID to use")
    parser.add_argument(
        "--model_name",
        help="Filtering model name to use in the filters. (e.g. VQSR or RF)",
    )
    parser.add_argument(
        "--score_name",
        help="Name for the score in the info struct in HT / INFO field inVCF (e.g. RF or VQSLOD",
    )
    parser.add_argument(
        "--info_fields",
        help="Fields other than the score to include in the finalized table info / INFO",
        nargs="+",
        type=str,
    )
    parser.add_argument(
        "--inbreeding_coeff_threshold",
        help="Name for the score in the info struct in HT / INFO field inVCF (e.g. RF or VQSLOD",
        type=float,
        default=-0.3,
    )
    parser.add_argument(
        "--snv_bin_cutoff", help="Cutoff for SNVs", default=90, type=int
    )
    parser.add_argument(
        "--indel_bin_cutoff", help="Cutoff for Indels", default=80, type=int
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    main(args)
