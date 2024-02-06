"Script to create custom liftover file for genes CBS, KCNE1, and CRYAA."
import argparse
import logging

import hail as hl
from gnomad.utils.annotations import (
    faf_expr,
    gen_anc_faf_max_expr,
    merge_freq_arrays,
    pop_max_expr,
)
from gnomad.utils.release import make_freq_index_dict_from_meta

from gnomad_qc.v2.annotations.generate_frequency_data import POPS_TO_REMOVE_FOR_POPMAX
from gnomad_qc.v2.resources.basics import get_gnomad_liftover_data_path
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.release import get_false_dup_genes_path

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("false_dup_genes")
logger.setLevel(logging.INFO)

FALSE_DUP_GENES = ["KCNE1", "CBS", "CRYAA"]


def filter_liftover_to_false_dups(
    data_type: str,
) -> hl.Table:
    """
    Read in gnomAD v2 liftover table, filter to genes of interest, and return the Table.

    :param data_type: Data type 'exomes' or 'genomes'. Default is 'exomes'.
    :return: Filtered Hail Table.
    """
    ht = hl.read_table(
        get_gnomad_liftover_data_path(data_type=data_type, version="2.1.1")
    )

    # Filter to chr21, since it is known that is where all 3 of the genes are.
    ht = ht.filter(ht.locus.contig == "chr21")

    # Filter to any variant located in one of the 3 genes of interest.
    ht = ht.filter(
        ht.vep.transcript_consequences.gene_symbol.any(
            lambda x: hl.literal(FALSE_DUP_GENES).contains(x)
        )
    )

    return ht


def main(args):
    """
    Writes merged v2 exomes and genomes liftover table with joint information for three clinically relevant genes.

    There are three clinically relevant genes impacted by false duplications in the GRCh38 reference (KCNE1, CBS, and CRYAA).
    To make it easier for our users to switch to gnomAD v4 and the new reference build,
    this script creates a new release file that combines information from the v2 exomes and genomes only for these three genes.
    To create this file, filter the v2 liftover files (https://gnomad.broadinstitute.org/downloads#v2-liftover) to these genes,
    and then merge frequency information across the exomes and genomes
    """
    # Read in liftover files.
    exome_ht = filter_liftover_to_false_dups("exomes")
    genome_ht = filter_liftover_to_false_dups("genomes")

    # Union and merge to contain any variant present in either.
    # Note: only about ~550 overlap, with ~187k exome variants and ~12k genome variants.
    # Genome global annotations stored with _1.
    ht = exome_ht.select().join(genome_ht.select(), how="outer")

    # Clear naming for global annotations.
    # Globals from exomes and genomes end with _exomes or _genomes.
    renaming_dictionary = {}
    for global_name in sorted(list(ht.globals)):
        if "_1" in global_name:
            renaming_dictionary[global_name] = global_name.replace("_1", "_genomes")
        else:
            renaming_dictionary[global_name] = f"{global_name}_exomes"

    ht = ht.rename(renaming_dictionary)

    # Annotate with information from both exomes and genomes.
    # Many will be NaN from the lack of overlap.
    ht = ht.annotate(
        v2_exomes=exome_ht[ht.locus, ht.alleles],
        v2_genomes=genome_ht[ht.locus, ht.alleles],
    )

    # Checkpoint output to new temp file.
    ht = ht.checkpoint(hl.utils.new_temp_file("three_dup_genes_liftover", "ht"))

    # Merge exomes and genomes frequencies.
    # Code adapted from
    # gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L155
    # as of 01-26-2025 with changing 'gen_anc' -> pop.
    freq, freq_meta = merge_freq_arrays(
        farrays=[ht.v2_exomes.freq, ht.v2_genomes.freq],
        fmeta=[
            ht.index_globals().freq_meta_exomes,
            ht.index_globals().freq_meta_genomes,
        ],
    )

    # Select to first 28 groups, containing:
    # Raw, adj, adj male & female, and adj sex-split pops for 8 pops.
    # Hail sets this to first 28 of each entry in the array, not first 28 of array.
    freq = freq[:28]
    freq_meta = freq_meta[:28]

    # Create a new struct to hold all joint information.
    ht = ht.annotate(v2_joint=hl.struct(joint_freq=freq))

    ht = ht.annotate_globals(
        joint_freq_meta=freq_meta,
        joint_freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta)),
    )

    # Checkpoint output to new temp file.
    ht = ht.checkpoint(hl.utils.new_temp_file("three_dup_genes_liftover", "ht"))

    # Compute FAF on the merged exomes + genomes frequencies.
    logger.info("Generating faf metrics...")
    faf, faf_meta = faf_expr(
        ht.v2_joint.joint_freq,
        ht.joint_freq_meta,
        ht.locus,
        pops_to_exclude=POPS_TO_REMOVE_FOR_POPMAX,
        pop_label="pop",
    )
    faf_meta_by_pop = {
        m.get("pop"): i for i, m in enumerate(faf_meta) if m.get("pop") and len(m) == 2
    }
    faf_meta_by_pop = hl.literal(faf_meta_by_pop)

    # Compute group max (popmax) on the merged exomes + genomes frequencies.
    logger.info("Computing grpmax...")
    grpmax = pop_max_expr(
        ht.v2_joint.joint_freq,
        ht.joint_freq_meta,
        pops_to_exclude=POPS_TO_REMOVE_FOR_POPMAX,
        pop_label="pop",
    )

    # Annotate Table with all joint exomes + genomes computations.
    ht = ht.annotate(
        v2_joint=ht.v2_joint.annotate(
            joint_faf=faf,
            joint_fafmax=gen_anc_faf_max_expr(
                faf, hl.literal(faf_meta), pop_label="pop"
            ),
            joint_grpmax=grpmax,
        )
    )

    ht = ht.annotate_globals(
        joint_faf_meta=faf_meta,
        joint_faf_index_dict=make_freq_index_dict_from_meta(hl.literal(faf_meta)),
    )

    # Checkpoint output to created resource.
    ht = ht.write(get_false_dup_genes_path(test=args.test), overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Option to overwrite existing custom liftover table.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Option to store output at test paths.",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
