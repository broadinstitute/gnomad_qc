"""script to create release sites ht for v4 genomes."""
import argparse
import logging

import hail as hl

from gnomad_qc.v4.resources.annotations import get_vep


def add_release_annotations(freq_ht: hl.Table) -> hl.Table:
    """
    Load and join all Tables with variant annotations.

    :param freq_ht: Table with frequency annotations
    :return: Table containing joined annotations
    """
    logger.info("Loading annotation tables...")

    vep_ht = get_vep(data_type="genomes").ht()

    vep_ht = vep_ht.annotate(vep=vep_ht.vep.drop("colocated_variants", "context"))

    vep_ht = vep_ht.annotate(
        vep=vep_ht.vep.annotate(
            transcript_consequences=vep_ht.vep.transcript_consequences.map(
                lambda x: x.drop("minimised", "swissprot", "trembl", "uniparc")
            )
        )
    )

    for consequence in [
        "intergenic_consequences",
        "motif_feature_consequences",
        "regulatory_feature_consequences",
    ]:
        vep_ht = vep_ht.annotate(
            vep=vep_ht.vep.annotate(
                **{
                    consequence: vep_ht.vep[consequence].map(
                        lambda x: x.drop("minimised")
                    )
                }
            )
        )
