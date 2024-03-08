"""Create v4.1 genomes release HT by dropping joint annotations from v4.0 genomes release HT."""

import argparse
import logging
from datetime import datetime

import hail as hl

from gnomad_qc.v4.resources.basics import qc_temp_prefix
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_genomes_v4_1_release_ht")
logger.setLevel(logging.INFO)

# Joint annotations to drop
joint_globals = [
    "joint_freq_meta",
    "joint_freq_index_dict",
    "joint_freq_meta_sample_count",
    "joint_faf_meta",
    "joint_faf_index_dict",
]
joint_rows = [
    "joint_freq",
    "joint_grpmax",
    "joint_faf",
    "joint_fafmax",
]


def main(args):
    """Create v4.1 genome release ht."""
    test = args.test

    logger.info("Loading v4.0 genomes release HT...")
    ht = release_sites(data_type="genomes").versions["4.0"].ht()

    if test:
        # Keep only PCSK9.
        logger.info("Filtering to PCSK9 region for test...")
        ht = hl.filter_intervals(
            ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )

    logger.info("Dropping joint annotations and writing v4.1 genomes release HT...")
    ht = ht.drop(*joint_rows + joint_globals)
    ht = ht.annotate_globals(
        date=datetime.now().isoformat(),
        version="4.1",
    )
    ht.describe()
    output_path = (
        f"{qc_temp_prefix(data_type='genomes')}release/gnomad.genomes.sites.test.{datetime.today().strftime('%Y-%m-%d')}.ht"
        if test
        else release_sites(data_type="genomes").versions["4.1"].path
    )
    ht.write(
        output_path,
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite script output.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Filter to PCSK9 region for testing purposes.",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
