from gnomad.utils.slack import slack_notifications
from gnomad.resources.grch38.gnomad import (
    SUBSETS,
    DOWNSAMPLINGS,
    POPS,
    POPS_TO_REMOVE_FOR_POPMAX,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
    age_hists_expr,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.vcf import make_label_combos, index_globals

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta

from typing import Tuple
import argparse
import logging
import hail as hl

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency_data")
logger.setLevel(logging.INFO)


def make_faf_index_dict(ht):
    """
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array
    :param Table ht: Table containing filter allele frequency annotations to be indexed
    :return: Dictionary of faf annotation population groupings, where values are the corresponding 0-based indices for the
        groupings in the faf_meta array
    :rtype: Dict of str: int
    """
    faf_meta = hl.eval(ht.faf_meta)

    index_dict = index_globals(faf_meta, dict(group=["adj"]))
    index_dict.update(index_globals(faf_meta, dict(group=["adj"], pop=POPS)))
    index_dict.update(index_globals(faf_meta, dict(group=["adj"], sex=SEXES)))
    index_dict.update(index_globals(faf_meta, dict(group=["adj"], pop=POPS, sex=SEXES)))

    return index_dict


def main(args):
    subset = args.subset
    log_tag = "." + subset if subset else ""
    hl.init(log=f"/generate_frequency_data{log_tag}.log", default_reference="GRCh38")

    try:
        logger.info("Reading full sparse MT and metadata table...")
        mt = get_gnomad_v3_mt(
            key_by_locus_and_alleles=True, release_only=True, samples_meta=True
        )

        if args.test:
            logger.info("Filtering to two partitions on chr20")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20:1-1000000")])
            mt = mt._filter_partitions(range(2))

        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

        if subset:
            mt = mt.filter_cols(mt.meta.subsets[subset], keep=True)
            logger.info(
                f"Running frequency generation pipeline on {mt.count_cols()} samples in {subset} subset..."
            )
        else:
            logger.info(
                f"Running frequency generation pipeline on {mt.count_cols()} samples..."
            )

        logger.info("Computing adj and sex adjusted genotypes...")
        mt = mt.annotate_entries(
            GT=adjusted_sex_ploidy_expr(
                mt.locus, mt.GT, mt.meta.sex_imputation.sex_karyotype
            ),
            adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
        )

        logger.info("Densify-ing...")
        mt = hl.experimental.densify(mt)
        mt = mt.filter_rows(hl.len(mt.alleles) > 1)

        # Temporary hotfix for depletion of homozygous alternate genotypes
        logger.info(
            "Setting het genotypes at sites with >1% AF (using v3.0 frequencies) and > 0.9 AB to homalt..."
        )
        # Load v3.0 allele frequencies to avoid an extra frequency calculation
        # NOTE: Using previous callset AF works for small incremental changes to a callset, but we will need to revisit for large increments
        freq_ht = get_freq(version="3").ht()
        freq_ht = freq_ht.select(AF=freq_ht.freq[0].AF)

        mt = mt.annotate_entries(
            GT=hl.cond(
                (freq_ht[mt.row_key].AF > 0.01)
                & mt.GT.is_het()
                & (mt.AD[1] / mt.DP > 0.9),
                hl.call(1, 1),
                mt.GT,
            )
        )

        logger.info("Generating frequency data...")
        if subset:
            mt = annotate_freq(
                mt,
                sex_expr=mt.meta.sex_imputation.sex_karyotype,
                pop_expr=mt.meta.population_inference.pop
                if subset not in ["tgp", "hgdp"]
                else mt.meta.project_meta.project_subpop,
            )

            # NOTE: no FAFs or popmax needed for subsets
            mt = mt.select_rows("freq")

            logger.info(f"Writing out frequency data for {subset} subset...")
            if args.test:
                mt.rows().write(
                    f"gs://gnomad-tmp/gnomad_freq/chr20_test_freq.{subset}.ht",
                    overwrite=True,
                )
            else:
                mt.rows().write(get_freq(subset=subset).path, overwrite=args.overwrite)

        else:
            logger.info("Computing age histograms for each variant...")
            mt = mt.annotate_rows(
                **age_hists_expr(
                    mt.adj,
                    mt.GT,
                    hl.if_else(
                        hl.is_defined(mt.meta.project_meta.age),
                        mt.meta.project_meta.age,
                        mt.meta.project_meta.age_alt,
                    ),
                )
            )

            # Compute callset-wide age histogram global
            meta_ht = meta.ht()
            meta_ht = meta_ht.filter(meta_ht.release)
            mt = mt.annotate_globals(
                age_distribution=meta_ht.aggregate(
                    hl.agg.hist(
                        hl.if_else(
                            hl.is_defined(meta_ht.project_meta.age),
                            meta_ht.project_meta.age,
                            meta_ht.project_meta.age_alt,
                        ),
                        30,
                        80,
                        10,
                    )
                )
            )

            mt = annotate_freq(
                mt,
                sex_expr=mt.meta.sex_imputation.sex_karyotype,
                pop_expr=mt.meta.population_inference.pop,
                downsamplings=DOWNSAMPLINGS,
            )
            # Remove all loci with raw AC=0
            mt = mt.filter_rows(mt.freq[1].AC > 0, keep=True)

            logger.info("Calculating InbreedingCoefficient...")
            # NOTE: This is not the ideal location to calculate this, but added here to avoid another densify
            mt = mt.annotate_rows(
                InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT)
            )

            # Select freq, FAF and popmax
            faf, faf_meta = faf_expr(
                mt.freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX
            )
            mt = mt.select_rows(
                "InbreedingCoeff",
                "freq",
                faf=faf,
                popmax=pop_max_expr(mt.freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX),
            )
            mt = mt.annotate_globals(faf_meta=faf_meta)
            mt = mt.annotate_globals(faf_index_dict=make_faf_index_dict(mt))
            mt = mt.annotate_rows(
                popmax=mt.popmax.annotate(
                    faf95=mt.faf[
                        mt.faf_meta.index(
                            lambda x: x.values() == ["adj", mt.popmax.pop]
                        )
                    ].faf95
                )
            )

            # Annotate quality metrics histograms, as these also require densifying
            mt = mt.annotate_rows(**qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD, mt.adj))
            ht = mt.rows()

            qual_hists = [field for field in ht.row_value if "_hist_" in field]
            adj_hists = [field for field in qual_hists if "_adj" in field]
            raw_hists = [field for field in qual_hists if "_adj" not in field]
            ht = ht.annotate(
                qual_hists=hl.Struct(
                    **{i.replace("_adj", ""): ht[i] for i in adj_hists}
                ),
                raw_qual_hists=hl.Struct(**{i: ht[i] for i in raw_hists}),
            )
            ht = ht.drop(*qual_hists)

            logger.info("Writing out frequency data...")
            if args.test:
                ht.write(
                    "gs://gnomad-tmp/gnomad_freq/chr20_test_freq.ht", overwrite=True
                )
            else:
                ht.write(get_freq().path, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log("gs://gnomad-tmp/logs/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of chr20", action="store_true"
    )
    parser.add_argument(
        "--subset",
        help="Name of subset for which to generate frequency data",
        choices=SUBSETS,
    )

    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
