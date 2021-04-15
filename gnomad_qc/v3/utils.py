import logging
from typing import List, Union

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def remove_fields_from_globals(global_field: List[str], fields_to_remove: List[str]):
    """
    Remove fields from the pre-defined global field variables.

    :param global_field: Global list of fields
    :param fields_to_remove: List of fields to remove from global (they must be in the global list)
    """
    for field in fields_to_remove:
        if field in global_field:
            global_field.remove(field)
        else:
            logger.info(f"'{field}'' missing from {global_field}")


def build_export_reference() -> hl.ReferenceGenome:
    """
    Create export reference based on GRCh38. Eliminates all non-standard contigs

    :return: Reference for VCF export containing chr1-22,X,Y, and M
    """
    ref = hl.get_reference("GRCh38")
    my_contigs = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

    export_reference = hl.ReferenceGenome(
        name="gnomAD_GRCh38",
        contigs=my_contigs,
        lengths={my_contig: ref.lengths[my_contig] for my_contig in my_contigs},
        x_contigs=ref.x_contigs,
        y_contigs=ref.y_contigs,
        par=[
            (interval.start.contig, interval.start.position, interval.end.position)
            for interval in ref.par
        ],
        mt_contigs=ref.mt_contigs,
    )

    return export_reference


def rekey_new_reference(
    t: Union[hl.Table, hl.MatrixTable], reference: hl.ReferenceGenome
):
    """
    Re-key Table or MatrixTable with a new reference genome.

    :param t: Input Table/MatrixTable.
    :param reference: Reference genome to re-key with.
    :return: Re-keyed Table/MatrixTable
    """
    t = t.rename({"locus": "locus_original"})
    locus_expr = hl.locus(
        t.locus_original.contig, t.locus_original.position, reference_genome=reference,
    )

    if isinstance(t, hl.MatrixTable):
        t = t.annotate_rows(locus=locus_expr)
        t = t.key_rows_by("locus", "alleles").drop("locus_original")
    else:
        t = t.annotate(locus=locus_expr)
        t = t.key_by("locus", "alleles").drop("locus_original")

    return t
