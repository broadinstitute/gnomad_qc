import argparse
import logging

import hail as hl

from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.resources.grch37.gnomad import liftover
from gnomad.resources.grch38.gnomad import public_release as v3_public_release
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("compare_genome_exome_af")
logger.setLevel(logging.INFO)

gnomad_public_resource_configuration.source = (
    GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
)

POPS = ["afr", "sas", "amr", "eas", "nfe"]
POP_FORMAT = {"v2": "gnomad_{}", "v3": "{}-adj"}


def main(args):
    hl.init(log="/compare_genome_exome_af.log")

    v3_ht = (
        v3_public_release("genomes")
        .ht()
        .select("freq", "filters")
        .select_globals("freq_index_dict")
    )
    v2_ht = (
        liftover("exomes")
        .ht()
        .select("freq", "filters")
        .select_globals("freq_index_dict")
    )

    if args.test:
        v3_ht = v3_ht._filter_partitions(range(2))
        v2_ht = v2_ht._filter_partitions(range(2))

    # Keep only variants present in both the genomes and exomes
    ht = v3_ht.join(v2_ht, how="inner",)
    ht = ht.select(
        freq_v3=ht.freq,
        freq_v2=ht.freq_1,
        filters_v3=ht.filters,
        filters_v2=ht.filters_1,
    )
    ht = ht.select_globals(
        freq_index_dict_v3=ht.freq_index_dict, freq_index_dict_v2=ht.freq_index_dict_1,
    )

    # Only keep sites that are PASS in both datasets
    ht = ht.filter((hl.len(ht.filters_v3) == 0) & (hl.len(ht.filters_v2) == 0))

    ht = ht.annotate(
        **{
            f"n_ref_{version}_{pop}": ht[f"freq_{version}"][
                ht[f"freq_index_dict_{version}"][POP_FORMAT[version].format(pop)]
            ].AN
            - ht[f"freq_{version}"][
                ht[f"freq_index_dict_{version}"][POP_FORMAT[version].format(pop)]
            ].AC
            for pop in POPS
            for version in ["v2", "v3"]
        },
        **{
            f"n_alt_{version}_{pop}": ht[f"freq_{version}"][
                ht[f"freq_index_dict_{version}"][POP_FORMAT[version].format(pop)]
            ].AC
            for pop in POPS
            for version in ["v2", "v3"]
        },
    )

    # Filtered to minimum AC in each population
    # min_ac = 5

    ht = ht.annotate(
        **{
            f"chi_squared_{pop}": hl.chi_squared_test(
                ht.n_alt_v3, ht.n_ref_v3, ht.n_alt_v2, ht.n_ref_v2
            )
            for pop in POPS
        }
    )

    # Biallelic sites


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset. (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Create release files using only 2 partitions on chr20, chrX, and chrY for testing purposes",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
