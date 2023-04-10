# noqa: D100
import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import sv_age_and_gq_hists
from gnomad_qc.v3.resources.basics import (
    gnomad_sv_release_samples_list_path,
    gnomad_sv_vcf_path,
    temp_gnomad_sv_mt_path,
)
from gnomad_qc.v3.resources.meta import meta

logger = logging.getLogger("generate_gnomad_sv_histograms")


def generate_hists(mt: hl.MatrixTable) -> hl.Table:
    """
    Generate histograms for age and GQ for gnomAD SVs.

    :param mt: gnomAD SV MatrixTable.
    :return: Table with histograms.
    """
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


def get_sample_age(sv_list: hl.Table) -> hl.Table:
    """
    Annotate sample list Table with age for generating histograms.

    :return: Prepared metadata Table.
    """
    sample_meta = meta.ht().key_by()

    # NOTE: some 1KG samples were already in v3.0 (SV data) and were given a new prefix in v3.1(meta)
    # To access the correct metadata using the SV sample list, we need to update the v3.1 meta IDs
    # to match the SV list Table.
    s_updates = hl.dict(
        {
            "v3.1::HG00512": "HG00512",
            "v3.1::HG00513": "HG00513",
            "v3.1::HG00731": "HG00731",
            "v3.1::HG00732": "HG00732",
        }
    )
    sample_meta = sample_meta.transmute(
        s=hl.if_else(
            s_updates.contains(sample_meta.s), s_updates[sample_meta.s], sample_meta.s
        )
    ).key_by("s")

    # NOTE: Add age to sample list. Most age data is stored as integers in 'age' annotation, # noqa
    #  but for a select number of samples, age is stored as a bin range and 'age_alt' # noqa
    #  corresponds to an integer in the middle of the bin # noqa
    sv_list = sv_list.annotate(
        age=hl.if_else(
            hl.is_defined(sample_meta[sv_list.s].project_meta.age),
            sample_meta[sv_list.s].project_meta.age,
            sample_meta[sv_list.s].project_meta.age_alt,
        ),
        release=sample_meta[sv_list.s].release,
    )
    logger.info(
        "%i out of %i samples in the SV sample list have a defined age.",
        sv_list.filter(hl.is_defined(sv_list.age)).count(),
        sv_list.count(),
    )
    logger.info(
        "%i out of %i samples in the SV sample list are in the core release",
        sv_list.filter(hl.is_defined(sv_list.release)).count(),
        sv_list.count(),
    )
    return sv_list


def main(args):
    """Generate age and GQ histograms for gnomAD SVs."""
    hl.init(
        log="/generate_gnomad_sv_histograms.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day/",
    )
    mt = hl.import_vcf(gnomad_sv_vcf_path, force_bgz=True, min_partitions=300)

    if args.test:
        logger.info("Running on chr22 only for testing purposes.")
        test_interval = [hl.parse_locus_interval("chr22")]
        mt = hl.filter_intervals(mt, test_interval)

    sv_list = hl.import_table(
        gnomad_sv_release_samples_list_path, force=True, no_header=True
    )
    sv_list = sv_list.transmute(s=sv_list.f0).key_by("s")
    sv_list = get_sample_age(sv_list)
    mt = mt.annotate_cols(**sv_list[mt.col_key])
    mt = mt.checkpoint(temp_gnomad_sv_mt_path, overwrite=args.overwrite)

    hists_ht = generate_hists(mt)
    hists_ht.write(
        temp_gnomad_sv_mt_path.replace(".mt", ".hists.ht")
        if args.test
        else sv_age_and_gq_hists.path,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--test", help="Test run", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
