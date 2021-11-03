"""
Compare frequencies for two gnomAD versions.

Generate a Hail Table containing frequencies for two gnomAD versions (specified using `args.versions_to_compare`),
and the following tests comparing the two frequencies:
    - Chi squared test
    - Fisher's exact test

The Table is written out to the location of the more recent of the two gnomAD versions being compared.
"""
import argparse
import logging

import hail as hl

from gnomad.resources.config import (
    gnomad_public_resource_configuration,
    GnomadPublicResourceSource,
)
from gnomad.resources.resource_utils import ResourceNotAvailable
from gnomad.resources.grch37.gnomad import (
    CURRENT_EXOME_RELEASE,
    CURRENT_GENOME_RELEASE,
    EXOME_POPS,
    GENOME_POPS,
    liftover,
)
from gnomad.resources.grch37.gnomad import public_release as v2_public_release
from gnomad.resources.grch38.gnomad import POPS
from gnomad.resources.grch38.gnomad import public_release as v3_public_release
from gnomad.resources.grch38.gnomad import CURRENT_GENOME_RELEASE as V3_CURRENT_RELEASE
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq_comparison

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("compare_freq")
logger.setLevel(logging.INFO)

POPS_MAP = {
    "v2_exomes": {pop.lower() for pop in EXOME_POPS},
    "v2_genomes": {pop.lower() for pop in GENOME_POPS},
    "v3_genomes": {pop.lower() for pop in POPS},
}
POP_FORMAT = {
    "v2_exomes": "gnomad_{}",
    "v2_genomes": "gnomad_{}",
    "v3_genomes": "{}-adj",
}
CURRENT_RELEASE_MAP = {
    "v2_exomes": CURRENT_EXOME_RELEASE,
    "v2_genomes": CURRENT_GENOME_RELEASE,
    "v3_genomes": V3_CURRENT_RELEASE,
}


def main(args):
    hl.init(log="/compare_freq.log")

    # Reorder so that v3_genomes is first in version. Forces the output location to be in more recent release location
    versions = [v for v in ["v3_genomes", "v2_genomes", "v2_exomes"] if v in args.versions_to_compare]
    v1, v2 = versions
    pops = list(POPS_MAP[v1] & POPS_MAP[v2])
    release_resource_map = {"v2": v2_public_release, "v3": v3_public_release}

    if "v3_genomes" in versions:
        release_resource_map["v2"] = liftover

    logger.info("Reading release HTs and extracting frequency info...")
    hts = []
    for version in versions:
        v, d = version.split("_")
        try:
            gnomad_public_resource_configuration.source = (
                GnomadPublicResourceSource.GOOGLE_CLOUD_PUBLIC_DATASETS
            )
            ht = release_resource_map[v](d).ht()
        except ResourceNotAvailable:
            logger.warning(
                "Using data in the gnomAD requester-pays bucket because the current datasets are not available from "
                "the Google Cloud Public Datasets. If the current machine is outside the US, egress charges will apply!"
            )
            gnomad_public_resource_configuration.source = (
                GnomadPublicResourceSource.GNOMAD
            )
            ht = release_resource_map[v](d).ht()

        logger.info("Filtering %s to only variants that are PASS...", version)
        ht = ht.filter(hl.len(ht.filters) == 0)
        freq_index_dict = hl.eval(ht.freq_index_dict)
        # Keep full gnomAD callset adj [0], raw [1], and ancestry-specific adj frequencies
        freq_idx = [0, 1] + [
            freq_index_dict[POP_FORMAT[version].format(pop)] for pop in pops
        ]

        logger.info(
            "Keeping only frequencies that are needed (adj, raw, adj by pop)..."
        )
        ht = ht.select(**{f"freq_{version}": [ht.freq[i] for i in freq_idx]})
        ht = ht.select_globals(
            **{f"freq_meta_{version}": [ht.freq_meta[i] for i in freq_idx]}
        )

        if args.test:
            ht = ht._filter_partitions(range(2))

        hts.append(ht)

    logger.info("Performing an inner join on frequency HTs...")
    ht = hts[0].join(hts[1], how="inner")

    ht = ht.annotate(
        n_ref=hl.struct(
            **{
                version: ht[f"freq_{version}"][2:].AN - ht[f"freq_{version}"][2:].AC
                for version in versions
            }
        ),
        n_alt=hl.struct(
            **{version: ht[f"freq_{version}"][2:].AC for version in versions}
        ),
    )

    logger.info(
        "Computing chi squared and fisher exact tests on per population frequencies..."
    )
    ht = ht.transmute(
        **{
            test: hl.struct(
                **{
                    pop: f(
                        ht.n_alt[v1][i],
                        ht.n_ref[v1][i],
                        ht.n_alt[v2][i],
                        ht.n_ref[v2][i],
                    )
                    for i, pop in enumerate(pops)
                }
            )
            for test, f in [
                ("chi_squared", hl.chi_squared_test),
                ("fisher_exact", hl.fisher_exact_test),
            ]
        }
    )

    ht.write(
        f"gs://gnomad-tmp/{v1}_{v2}.compare_freq.test.ht"
        if args.test
        else get_freq_comparison(
            CURRENT_RELEASE_MAP[v1],
            v1.split("_")[1],
            CURRENT_RELEASE_MAP[v2],
            v2.split("_")[1],
        ).path,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset. (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Filter Tables to only the first 2 partitions for testing.",
        action="store_true",
    )
    parser.add_argument(
        "--versions-to-compare",
        help="Two gnomAD versions to compare allele frequencies.",
        nargs=2,
        choices=["v3_genomes", "v2_exomes", "v2_genomes"],
        default=["v3_genomes", "v2_exomes"],
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
