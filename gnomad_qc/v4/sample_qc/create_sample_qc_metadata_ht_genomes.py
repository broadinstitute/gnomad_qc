"""
Script to create sample QC metadata HT for genomes.

This script is used to create the sample QC metadata HT for genomes. The main updates
from v3.1 to v4.0 are:
Sample QC annotations for HGDP + TGP subset, this script directly merges the updated
annotations from the subset meta HT to the full meta HT.
"""


import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import TableResource
from hail.utils import new_temp_file

from gnomad_qc.v3.resources.basics import meta as v3_meta
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v4.resources.meta import meta_tsv_path as v4_meta_tsv_path
from gnomad_qc.v4.resources.sample_qc import hgdp_tgp_meta_updated

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("Create v4.0 genomes sample QC metadata HT.")
logger.setLevel(logging.INFO)


_meta_versions = {
    "4.0": TableResource(
        path="gs://gnomad/v4.0/metadata/exomes/gnomad.exomes.v4.0.sample_qc_metadata.ht"
    ),
}


def import_updated_annotations(ht: hl.Table, subset_ht: hl.Table) -> hl.Table:
    """
    Import updated annotations from subset meta HT to full meta HT.

    :param ht: v3.1 genomes meta HT of genomes
    :param subset_ht: Updated meta HT of HGDP + TGP subset
    :return: Updated v4 genomes meta HT
    """
    # Some HGDP/TGP samples in v3.1 meta HT have a prefix of "v3.1::" in the key 's',
    # but the project_meta.sample_id field does not have this prefix. This is used as
    # a temp key to join the subset meta HT.
    release_v3 = ht.filter(ht.release).count()
    ht1 = ht.filter(ht.subsets.hgdp | ht.subsets.tgp)
    ht1 = ht1.naive_coalesce(5).checkpoint(new_temp_file("v4.0_genomes_meta", "ht"))
    ht2 = ht.filter(~(ht.subsets.hgdp | ht.subsets.tgp))
    ht2 = ht2.checkpoint(new_temp_file("v4.0_genomes_meta", "ht"))

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
    ht = ht1.union(ht2)

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
        "Number of genome release samples will change from %s in v3.1 to %s in v4.0.",
        release_v3,
        release_v4,
    )

    return ht


def main(args):
    """Create v4.0 genomes sample QC metadata HT."""
    hl.init(
        log="/create_v4.0_genomes_meta.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    logger.info("Importing v3.1 genomes meta HT...")
    ht = v3_meta.ht()
    logger.info("Importing updated HGDP + TGP meta HT...")
    subset_ht = hgdp_tgp_meta_updated.ht()
    logger.info("Updating annotations...")
    ht = import_updated_annotations(ht, subset_ht)
    ht = ht.naive_coalesce(160)
    ht.write(v4_meta(data_type="genomes").path, overwrite=args.overwrite)
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
    args = parser.parse_args()
    main(args)
