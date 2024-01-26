"Script to create custom liftover file for genes CBS, KCNE1, and CRYAA."
import argparse
import logging

import hail as hl

from gnomad.utils.annotations import merge_freq_arrays
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad_qc.v2.resources.basics import get_gnomad_liftover_data_path
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.release import _release_root

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("false_dup_genes")
logger.setLevel(logging.INFO)

FALSE_DUP_GENES = ["KCNE1", "CBS", "CRYAA"]


def get_false_dup_genes_path(
    release_version: str = CURRENT_RELEASE,
    test: bool = False,
    data_type: str = "genomes",
) -> str:
    """
    Retrieve path for the liftover table containing three genes of interest within a false duplication in GRCh38.

    :return: Combined custom liftover table path for the three genes in false duplication.
    """
    return (
        f"{_release_root(version=release_version, test=test, data_type=data_type)}/gnomad.{data_type}.v{release_version}_three_false_dup_genes_liftover.ht"
    )


def read_and_filter(
    data_type: str,
) -> hl.Table:
    """
    Read in gnomAD v2 liftover table, filter to genes of interest, and return the Table.

    :param data_type: String of either exome of genome
    :return: filtered Hail Table
    """
    ht = hl.read_table(
        get_gnomad_liftover_data_path(data_type=data_type, version="2.1.1")
    )

    # Filter to chr21, since it is known that is where all 3 of the genes are.
    ht = ht.filter(ht.locus.contig == "chr21")

    # Filter using lambda statement for any of the 3 genes of interest in the gene_symbol array.
    ht = ht.filter(
        hl.set(FALSE_DUP_GENES).any(
            lambda gene_symbol: ht.vep.transcript_consequences.gene_symbol.contains(
                gene_symbol
            )
        )
    )

    return ht


def main(args):
    """
    These are the three clinically relevant genes impacted by false duplications in the GRCh38 reference.
    To make it easier for our users to switch to gnomAD v4 and the new reference build,
    we should create a new release file that combines information from the v2 exomes and genomes only for these three genes.
    To create this file, filter the v2 liftover files (https://gnomad.broadinstitute.org/downloads#v2-liftover) to these genes,
    and then merge frequency information across the exomes and genomes
    """
    # Read in liftover files.
    exome_ht = read_and_filter("exomes")
    genome_ht = read_and_filter("genomes")

    # Union and merge to contain any variant present in either.
    # Note: only about ~550 overlap, with ~187k exome variants and ~12k genome variants.
    # This makes sense in exonic regions.
    ht = exome_ht.select().join(genome_ht.select(), how="outer").distinct()

    # Annotate with information from both exomes and genomes - many will be NaN.
    ht = ht.annotate(
        v2_exomes=exome_ht[ht.locus, ht.alleles],
        v2_genomes=genome_ht[ht.locus, ht.alleles],
    )

    # Annotate with globals from v2 exomes and v2 genomes.
    ht = ht.annotate_globals(
        v2_exome_globals=exome_ht.globals.collect(),
        v2_genome_globals=genome_ht.globals.collect(),
    )

    # Checkpoint output to created resource.
    ht = ht.checkpoint(hl.utils.new_temp_file("three_dup_genes_liftover", "ht"))

    # Merge exomes and genomes frequencies.
    # Code adapted from
    # gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L155
    # as of 01-26-2025
    freq, freq_meta, count_arrays_dict = merge_freq_arrays(
        farrays=[ht.v2_exomes.freq, ht.v2_genomes.freq],
        fmeta=[ht.index_globals().freq_meta, ht.index_globals().freq_meta_1],
        count_arrays={
            "counts": [
                ht.index_globals().freq_index_dict.values(),
                ht.index_globals().freq_index_dict_1.values(),
            ],
        },
    )

    ht = ht.annotate(v2_joint_freq=freq)
    ht = ht.annotate_globals(
        joint_freq_meta=freq_meta,
        joint_freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta)),
    )

    # Checkpoint output to created resource.
    ht = ht.checkpoint(
        get_false_dup_genes_path(test=args.test), overwrite=args.overwrite
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Option to overwrite existing custom liftover table",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Option to store output at test paths",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
