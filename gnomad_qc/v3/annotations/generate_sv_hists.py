# noqa: D100
import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import sv_age_and_gq_hists
from gnomad_qc.v3.resources.basics import (
    gnomad_sv_bucket_path,
    gnomad_sv_release_samples_list_path,
    temp_gnomad_sv_mt_path,
)
from gnomad_qc.v3.resources.meta import meta

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("generate_gnomad_sv_histograms")
logger.setLevel(logging.INFO)


def generate_hists(mt: hl.MatrixTable) -> hl.Table:
    """
    Generate histograms for age and GQ for gnomAD SVs.

    :param mt: gnomAD SV MatrixTable.
    :return: Table with histograms.
    """
    logger.info("Generating histograms for age and GQ...")
    hists = mt.select_rows(
        age_hist_het=hl.or_missing(
            ~mt.info.MULTIALLELIC,
            hl.agg.filter(mt.GT.is_het(), hl.agg.hist(mt.age, 30, 80, 10)),
        ),
        age_hist_hom=hl.or_missing(
            ~mt.info.MULTIALLELIC,
            hl.agg.filter(mt.GT.is_hom_var(), hl.agg.hist(mt.age, 30, 80, 10)),
        ),
        gq_hist_alt=hl.or_missing(
            ~mt.info.MULTIALLELIC,
            hl.agg.filter(mt.GT.is_non_ref(), hl.agg.hist(mt.GQ, 0, 100, 20)),
        ),
        gq_hist_all=hl.or_missing(
            ~mt.info.MULTIALLELIC, hl.agg.hist(mt.GQ, 0, 100, 20)
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
    # TODO: Check with browser team if they want to keep the format of the
    # histograms as is or if we can store bin edges as globals
    return hists.select(**hist_expr)


def get_sample_age(sv_list: hl.Table) -> hl.Table:
    """
    Annotate sample list Table with age for generating histograms.

    :return: Prepared metadata Table.
    """
    sample_meta = meta.ht().key_by()

    # NOTE: some 1KG samples were already in v3.0 (SV data) and were given a new prefix in v3.1 (meta)
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
    #  but for a select number of samples, age was recorded as a bin range and store in  # noqa
    #  the 'age_bin' annotation.'age_alt' corresponds to the average of the 'age_bin' edges # noqa
    sv_list = sv_list.annotate(
        age=hl.if_else(
            hl.is_defined(sample_meta[sv_list.s].project_meta.age),
            sample_meta[sv_list.s].project_meta.age,
            sample_meta[sv_list.s].project_meta.age_alt,
        ),
        release=sample_meta[sv_list.s].release,
    )
    return sv_list


def get_sex_and_autosome_mt() -> hl.MatrixTable:
    """
    Read in and union autosome and sex chromosome VCFs.

    :return: MatrixTable containing calls from autosome and sex chromosome VCFs.
    """
    logger.info("Importing VCFs...")
    # NOTE: The sex chromosome VCFs have an extra "PAR" field that needs to
    # be dropped in order to union with the autosome VCFs
    sex_chr_vcf_paths = [
        f"{gnomad_sv_bucket_path}/gnomAD.v3.SV.chr{s}.vcf.gz" for s in ("X", "Y")
    ]
    s_mt = hl.import_vcf(
        sex_chr_vcf_paths,
        reference_genome="GRCh38",
        min_partitions=300,
        force_bgz=True,
    )
    s_mt = s_mt.annotate_rows(info=s_mt.info.drop("PAR"))

    autosomes_vcf_paths = [
        f"{gnomad_sv_bucket_path}/gnomAD.v3.SV.chr{i}.vcf.gz" for i in range(1, 23)
    ]
    a_mt = hl.import_vcf(
        autosomes_vcf_paths,
        reference_genome="GRCh38",
        min_partitions=300,
        force_bgz=True,
    )
    mt = s_mt.union_rows(a_mt)
    return mt


def main(args):
    """Generate age and GQ histograms for gnomAD SVs."""
    hl.init(
        log="/generate_gnomad_sv_histograms.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day/",
    )
    logger.info("Running generate_gnomad_sv_histograms.py...")
    mt = get_sex_and_autosome_mt()
    if args.test:
        logger.info("Running on chr22 only for testing purposes.")
        mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr22")])
        logger.info(
            "Imported VCF with %i variants and %i samples",
            mt.count_rows(),
            mt.count_cols(),
        )

    logger.info("Importing sample list...")
    sv_list = hl.import_table(
        gnomad_sv_release_samples_list_path, force=True, no_header=True
    )
    sv_list = sv_list.transmute(s=sv_list.f0).key_by("s")

    logger.info("Annotating sample list with age...")
    sv_list = get_sample_age(sv_list)

    logger.info("Checkpoint SV MT...")
    mt = mt.annotate_cols(**sv_list[mt.col_key])
    mt = mt.checkpoint(temp_gnomad_sv_mt_path, overwrite=args.overwrite)
    logger.info(
        "Generating age and GQ histograms for %s variants and %s samples...",
        mt.count_rows(),
        mt.count_cols(),
    )

    hists_ht = generate_hists(mt)
    hists_ht.write(
        (
            temp_gnomad_sv_mt_path.replace(".mt", ".hists.ht")
            if args.test
            else sv_age_and_gq_hists.path
        ),
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
