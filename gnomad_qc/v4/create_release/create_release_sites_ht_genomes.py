"""script to create release sites ht for v4 genomes."""
import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import dbsnp

from gnomad_qc.v3.resources.release import release_sites
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


def replace_oth_with_remaining(ht: hl.Table) -> hl.Table:
    """
    Replace 'oth' with 'remaining' in Global fields of a Table.

    .. note::
        This function is used to rename ancestry group fields to match v4 exomes.
        Which means that all "oth" will be changed to "remaining".
        2 Global fields are renamed: freq_meta (in the value of the key "pop") and freq_index_dict (in the keys).

    :param ht: release sites Table to be modified
    :return: release sites Table with 'oth' replaced with 'remaining'
    """
    # in freq_meta, type: array<dict<str, str>>, the value of key "pop" is
    # changed from 'oth' to 'remaining'.
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: x.replace("oth", "remaining"))
        )
    )

    # in freq_index_dict, type: dict<str, int>, the keys containing 'oth'
    # inside are changed to 'remaining'.
    oldkeys = hl.eval(ht.freq_index_dict.keys())
    newkeys = [s.replace("oth", "remaining") for s in oldkeys]
    vals = hl.eval(ht.freq_index_dict.values())
    freq_index_dict_new = {k: v for k, v in zip(newkeys, vals)}
    ht = ht.annotate_globals(freq_index_dict=freq_index_dict_new)
    return ht
