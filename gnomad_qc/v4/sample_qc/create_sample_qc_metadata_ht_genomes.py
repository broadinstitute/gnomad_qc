"""
Script to create sample QC metadata HT for genomes.

The main updates from v3.1 to v4.0 are due to the following updates in the HGDP + TGP subset meta:
- 169 samples added
- 110 samples removed

The 169 samples added were diverse HGDP/TGP samples in v3.1 that were unintentionally removed during sample outlier detection.

The 110 samples were removed due to the following reasons:
- Contamination flagged by CHARR (CHARR score stored in `freemix` column).
- Subcontinental PCA outliers within the subset.
- Updated relatedness results within the subset, re-evaluation of relatedness
  inference of the subset to the rest of gnomAD release samples, and removal of any
  samples that are duplicates of v4.0 exomes release samples.

All samples impacted (169 added, 110 released) are releasable; the updates we made here were to their sample QC status (high quality vs filtered).
"""


import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import TableResource
from gnomad.utils.slack import slack_notifications
from hail.utils import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.basics import meta as v3_meta
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v4.resources.meta import meta_tsv_path as v4_meta_tsv_path
from gnomad_qc.v4.resources.sample_qc import hgdp_tgp_meta_updated

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_v4_genomes_sample_meta")
logger.setLevel(logging.INFO)

N_DIFF_SAMPLES = 169 - 110
"""
The number of samples that have been updated between v3.1 and v4.0.
"""


def import_updated_annotations(ht: hl.Table, subset_ht: hl.Table) -> hl.Table:
    """
    Import updated annotations from HGDP/TGP subset meta HT and annotate full meta HT.

    :param ht: v3.1 genomes meta HT
    :param subset_ht: Updated HGDP + TGP subset meta HT
    :return: Updated v4.0 genomes meta HT
    """
    # Some HGDP/TGP samples in v3.1 meta HT have a prefix of "v3.1::" in the key 's',
    # but the project_meta.sample_id field does not have this prefix. This is used as
    # a temp key to join the subset meta HT.
    release_v3 = ht.filter(ht.release).count()
    ht1 = ht.filter(ht.subsets.hgdp | ht.subsets.tgp)
    ht1 = ht1.checkpoint(new_temp_file("genomes_meta", "ht"))
    ht = ht.filter(~(ht.subsets.hgdp | ht.subsets.tgp))
    ht = ht.checkpoint(new_temp_file("genomes_meta", "ht"))

    ht1 = ht1.key_by(ht1.project_meta.sample_id)
    annot = subset_ht[ht1.key]
    subpop = annot.hgdp_tgp_meta.population
    freemix = annot.bam_metrics.freemix
    hard_filters = annot.gnomad_sample_filters.hard_filters
    hard_filtered = annot.gnomad_sample_filters.hard_filtered
    release_related = (
        annot.relatedness_inference.related
        | annot.gnomad_sample_filters.related_to_nonsubset
        | annot.gnomad_sample_filters.v4_exome_duplicate
    )
    high_quality = annot.gnomad_high_quality
    release = annot.gnomad_release

    ht1 = ht1.annotate(
        project_meta=ht1.project_meta.annotate(project_subpop=subpop),
        bam_metrics=ht1.bam_metrics.annotate(freemix=freemix),
        sample_filters=ht1.sample_filters.annotate(
            hard_filters=hard_filters,
            hard_filtered=hard_filtered,
            release_related=release_related,
        ),
        high_quality=high_quality,
        release=release,
    )

    ht1 = ht1.key_by("s").drop("sample_id")
    ht = ht1.union(ht)

    # Update gnomAD inferred population 'oth' to be 'remaining'.
    ht = ht.annotate(
        population_inference=ht.population_inference.annotate(
            pop=hl.if_else(
                ht.population_inference.pop == "oth",
                "remaining",
                ht.population_inference.pop,
            )
        )
    )
    release_v4 = ht.filter(ht.release).count()
    logger.info(
        "Number of genome release samples will change from %s in v3.1 to %s in v4.0",
        release_v3,
        release_v4,
    )
    if release_v4 <= release_v3 or release_v4 != release_v3 + N_DIFF_SAMPLES:
        raise ValueError("Number of release samples in v4.0 is not as expected.")

    return ht


def main(args):
    """Create updated v4 genomes sample QC metadata HT."""
    hl.init(
        log="/create_v4_genomes_meta.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    logger.info("Importing v3.1 genomes meta HT...")
    ht = v3_meta.ht()
    logger.info("Importing updated HGDP + TGP meta HT...")
    subset_ht = hgdp_tgp_meta_updated.ht()
    logger.info("Updating annotations...")
    ht = import_updated_annotations(ht, subset_ht)
    logger.info("Writing HT and exporting TSV...")
    ht.write(
        v4_meta(data_type="genomes").path,
        overwrite=args.overwrite,
    )
    ht.export(v4_meta_tsv_path(data_type="genomes"), delimiter="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        """This script is used to create the sample QC metadata HT for genomes."""
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite the existing meta HT.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
