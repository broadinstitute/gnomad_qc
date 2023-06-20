"""script to create release sites ht for v4 genomes."""
import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import dbsnp

from gnomad_qc.v4.resources.annotations import get_vep
from gnomad_qc.v4.resources.constants import CURRENT_VERSION

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)


def remove_missing_vep_fields(vep_ht: hl.Table) -> hl.Table:
    """
    Remove fields from VEP 105 annotations that are missing in all rows.

    :param vep_ht: Table containing VEP 105 annotations
    :return: Table containing VEP 105 annotations with missing fields removed
    """
    logger.info("Loading annotation tables...")

    vep_ht = get_vep(version=CURRENT_VERSION, data_type="genomes").ht()

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
    return vep_ht
