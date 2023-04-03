# noqa: D100
import argparse
import logging

from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import sv_age_and_gq_hists
from gnomad_qc.v3.resources.basics import (
    gnomad_sv_release_samples_list_path,
    gnomad_sv_vcf_path,
    temp_gnomad_sv_mt_path,
)
from gnomad_qc.v3.resources.meta import get_gnomad_meta

logger = logging.getLogger("generate_gnomad_sv_histograms")


def generate_hists(mt: hl.MatrixTable) -> hl.Table:
    """
    Generate histograms for age and GQ for gnomAD SVs.

    :param mt: gnomAD SV MatrixTable.
    :return: Table with histograms.
    """
    meta = get_gnomad_meta("genomes")
    meta = meta.key_by(s=(meta.project_id + "_" + meta.s).replace("\W+", "_"))
    mt = mt.annotate_cols(
        age=meta[mt.col_key].age, in_v2=hl.is_defined(meta[mt.col_key])
    )

    logger.info(
        "Found %i samples with age data.",
        mt.aggregate_cols(hl.agg.count_where(hl.is_defined(mt.age))),
    )

    hists = mt.select_rows(
        age_hist_het=hl.or_missing(
            ~mt.filters.contains("MULTIALLELIC"),
            hl.agg.filter(mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
        ),
        age_hist_hom=hl.or_missing(
            ~mt.filters.contains("MULTIALLELIC"),
            hl.agg.filter(mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10)),
        ),
        gq_hist_alt=hl.or_missing(
            ~mt.filters.contains("MULTIALLELIC"),
            hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.GQ, 0, 100, 20)),
        ),
        gq_hist_all=hl.or_missing(
            ~mt.filters.contains("MULTIALLELIC"), hl.agg.hist(mt.GQ, 0, 100, 20)
        ),
    ).rows()

    hist_expr = {}
    for hist in ["age_hist_het", "age_hist_hom", "gq_hist_alt", "gq_hist_all"]:
        hist_expr.update(
            {
                f"{hist}_bin_freq": hl.delimit(hists[hist].bin_freq, delimiter="|"),
                f"{hist}_bin_edges": hl.delimit(hists[hist].bin_edges, delimiter="|"),
                f"{hist}_n_smaller": hists[hist].n_smaller,
                f"{hist}_n_larger": hists[hist].n_larger,
            }
        )

    hists = hists.select(**hist_expr)

    return hists


def main(args):
    """Generate age and GQ histograms for gnomAD SVs."""
    hl.init(
        log="/generate_gnomad_sv_histograms.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day/",
    )
    mt = hl.import_vcf(gnomad_sv_vcf_path, force_bgz=True, min_partitions=300)
    meta = hl.import_table(gnomad_sv_release_samples_list_path, force=True, key="s")
    mt = mt.annotate_cols(release=hl.is_defined(meta[mt.col_key]))
    mt = mt.checkpoint(temp_gnomad_sv_mt_path, overwrite=args.overwrite)

    hists_ht = generate_hists(mt)
    hists_ht.write(
        sv_age_and_gq_hists,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--import_vcf",
        help="Imports gnomAD SV VCF and writes it as MT.",
        action="store_true",
    )
    parser.add_argument(
        "--generate_hists", help="Generate age and GQ histograms", action="store_true"
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
