import argparse
import logging
from typing import Dict, List

import hail as hl

from gnomad.resources.grch38.gnomad import (
    COHORTS_WITH_POP_STORED_AS_SUBPOP,
    DOWNSAMPLINGS,
    POPS,
    POPS_TO_REMOVE_FOR_POPMAX,
    SEXES,
    SUBSETS,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
)
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import index_globals

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq
from gnomad_qc.v3.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v3_mt,
    qc_temp_prefix,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency_data")
logger.setLevel(logging.INFO)


def make_faf_index_dict(faf_meta: List[Dict[str, str]]) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array

    :param List of Dict faf_meta: Global annotation containing the set of groupings for each element of the faf array (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :return: Dictionary of faf annotation population groupings, where values are the corresponding 0-based indices for the
        groupings in the faf_meta array
    :rtype: Dict of str: int
    """

    index_dict = index_globals(faf_meta, dict(group=["adj"]))
    index_dict.update(index_globals(faf_meta, dict(group=["adj"], pop=POPS)))
    index_dict.update(index_globals(faf_meta, dict(group=["adj"], sex=SEXES)))
    index_dict.update(index_globals(faf_meta, dict(group=["adj"], pop=POPS, sex=SEXES)))

    return index_dict


def main(args):
    subsets = args.subsets
    hl.init(
        log=f"/generate_frequency_data{'.' + '_'.join(subsets) if subsets else ''}.log",
        default_reference="GRCh38",
    )

    invalid_subsets = []
    n_subsets_in_use_subpops = 0
    for s in subsets:
        if s not in SUBSETS:
            invalid_subsets.append(s)
        if s in COHORTS_WITH_POP_STORED_AS_SUBPOP:
            n_subsets_in_use_subpops += 1

    if len(invalid_subsets) > 0:
        raise ValueError(
            f"{', '.join(invalid_subsets)} subset(s) are not one of the following official subsets: {SUBSETS}"
        )
    if (n_subsets_in_use_subpops > 0) & (n_subsets_in_use_subpops != len(subsets)):
        raise ValueError(
            f"All or none of the supplied subset(s) should be in the list of cohorts that need to use subpops instead "
            f"of pops in frequency calculations: {COHORTS_WITH_POP_STORED_AS_SUBPOP}"
        )

    try:
        logger.info("Reading full sparse MT and metadata table...")
        mt = get_gnomad_v3_mt(
            key_by_locus_and_alleles=True,
            release_only=not args.include_non_release,
            samples_meta=True,
        )

        logger.info("Filtering MT columns to high quality and HGDP + TGP samples")
        mt = mt.filter_cols(
            mt.meta.high_quality & (mt.meta.subsets.hgdp | mt.meta.subsets.tgp)
        )

        print("Numcols: ", mt.count_cols())
        mt.describe()

        if args.test:
            logger.info("Filtering to two partitions on chr20")
            mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20:1-1000000")])
            mt = mt._filter_partitions(range(2))

        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

        if subsets:
            mt = mt.filter_cols(hl.any([mt.meta.subsets[s] for s in subsets]))
            logger.info(
                f"Running frequency generation pipeline on {mt.count_cols()} samples in {', '.join(subsets)} subset(s)..."
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
        if subsets:
            mt = annotate_freq(
                mt,
                sex_expr=mt.meta.sex_imputation.sex_karyotype,
                pop_expr=mt.meta.population_inference.pop
                if n_subsets_in_use_subpops == 0
                else mt.meta.project_meta.project_subpop,
                # NOTE: TGP and HGDP labeled populations are highly specific and are stored in the project_subpop meta field
            )

            # NOTE: no FAFs or popmax needed for subsets
            mt = mt.select_rows("freq")
            mt = mt.filter_rows(
                mt.freq[1].AC > 0, keep=True
            )  # TODO: Not in master do we need?

            logger.info(
                f"Writing out frequency data for {', '.join(subsets)} subset(s)..."
            )
            if args.test:
                mt.rows().write(
                    get_checkpoint_path(f"chr20_test_freq.{'_'.join(subsets)}.ht"),
                    overwrite=True,
                )
            else:
                mt.rows().write(
                    get_freq(subset="_".join(subsets)).path, overwrite=args.overwrite
                )

        else:
            logger.info("Computing age histograms for each variant...")
            mt = mt.annotate_cols(
                age=hl.if_else(
                    hl.is_defined(mt.meta.project_meta.age),
                    mt.meta.project_meta.age,
                    mt.meta.project_meta.age_alt,
                    # NOTE: most age data is stored as integers in 'age' annotation, but for a select number of samples, age is stored as a bin range and 'age_alt' corresponds to an integer in the middle of the bin
                )
            )
            mt = mt.annotate_rows(**age_hists_expr(mt.adj, mt.GT, mt.age,))

            # Compute callset-wide age histogram global
            mt = mt.annotate_globals(
                age_distribution=mt.aggregate_cols(hl.agg.hist(mt.age, 30, 80, 10,))
            )

            mt = annotate_freq(
                mt,
                sex_expr=mt.meta.sex_imputation.sex_karyotype,
                pop_expr=mt.meta.population_inference.pop,
                downsamplings=DOWNSAMPLINGS,
            )
            # Remove all loci with raw AC=0
            mt = mt.filter_rows(mt.freq[1].AC > 0)

            logger.info("Calculating InbreedingCoeff...")
            # NOTE: This is not the ideal location to calculate this, but added here to avoid another densify
            mt = mt.annotate_rows(
                InbreedingCoeff=bi_allelic_site_inbreeding_expr(mt.GT)
            )

            logger.info("Computing filtering allele frequencies and popmax...")
            faf, faf_meta = faf_expr(
                mt.freq, mt.freq_meta, mt.locus, POPS_TO_REMOVE_FOR_POPMAX
            )
            mt = mt.select_rows(
                "InbreedingCoeff",
                "freq",
                faf=faf,
                popmax=pop_max_expr(mt.freq, mt.freq_meta, POPS_TO_REMOVE_FOR_POPMAX),
            )
            mt = mt.annotate_globals(
                faf_meta=faf_meta, faf_index_dict=make_faf_index_dict(faf_meta)
            )
            mt = mt.annotate_rows(
                popmax=mt.popmax.annotate(
                    faf95=mt.faf[
                        mt.faf_meta.index(
                            lambda x: x.values() == ["adj", mt.popmax.pop]
                        )
                    ].faf95
                )
            )

            logger.info("Annotating quality metrics histograms...")
            # NOTE: these are performed here as the quality metrics histograms also require densifying
            mt = mt.annotate_rows(
                qual_hists=qual_hist_expr(mt.GT, mt.GQ, mt.DP, mt.AD, mt.adj)
            )
            ht = mt.rows()
            ht = ht.annotate(
                qual_hists=hl.Struct(
                    **{
                        i.replace("_adj", ""): ht.qual_hists[i]
                        for i in ht.qual_hists
                        if "_adj" in i
                    }
                ),
                raw_qual_hists=hl.Struct(
                    **{i: ht.qual_hists[i] for i in ht.qual_hists if "_adj" not in i}
                ),
            )

            logger.info("Writing out frequency data...")
            if args.test:
                ht.write(get_checkpoint_path("chr20_test_freq"), overwrite=True)
            else:
                ht.write(get_freq().path, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(f"gs://{qc_temp_prefix()}/logs/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of chr20.", action="store_true"
    )
    parser.add_argument(
        "--include_non_release",
        help="Includes un-releasable samples in the frequency calculations.",
        action="store_true",
    )
    parser.add_argument(
        "--subsets",
        help="Comma separated list of subsets for which to generate combined frequency data.",
        nargs="+",
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
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
