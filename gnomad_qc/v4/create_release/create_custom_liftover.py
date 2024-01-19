" Script to create custom liftover file for genes CBS, KCNE1, and CRYAA"
import argparse
import logging

import hail as hl

from gnomad.resources.resource_utils import (
    MatrixTableResource,
    PedigreeResource,
    TableResource,
    VersionedMatrixTableResource,
    VersionedPedigreeResource,
    VersionedTableResource,
)
from gnomad_qc.v2.resources.basics import get_gnomad_liftover_data_path
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE, RELEASES
from gnomad_qc.v4.resources.release import (
    # release_lof,
    _release_root,
    release_sites,
    # release_summary_stats
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)


def get_custom_liftover_three_genes_path(
    release_version: str = CURRENT_RELEASE,
    test: bool = False,
    data_type: str = "genomes",
) -> str:
    """
    Retrieve path for the liftover table containing three genes of interest.

    :return: Combined custom liftover table path.
    """
    return (
        f"{_release_root(version=release_version, test=test, data_type=data_type)}/gnomad.{data_type}.v{release_version}_custom_liftover_three_genes.ht"
    )


def main(args):
    # Read in liftover file
    liftover_ht = hl.read_table(
        get_gnomad_liftover_data_path(data_type="exomes", version="2.1.1").path
    )

    # Create and filter to genes of interest
    # Current coordinates grabbed from gnomAD v4 browser, have reached out to browser team about something more professional.
    # CBS: 21:43053191-43076943
    cbs_interval = hl.locus_interval(
        contig="chr21", start=43053191, end=43076943, reference_genome="GRCh38"
    )
    # KCNE1: 21:34446688-34512214
    kcne1_interval = hl.locus_interval(
        contig="chr21", start=34446688, end=34512214, reference_genome="GRCh38"
    )
    # CRYAA interval: 21:43169008-43172805
    cryaa_interval = hl.locus_interval(
        contig="chr21", start=43169008, end=43172805, reference_genome="GRCh38"
    )
    interval_list = [kcne1_interval, cbs_interval, cryaa_interval]

    # Create liftover ht at only three relevant genes
    lifted_three_genes = hl.filter_intervals(liftover_ht, interval_list)

    # Read in v4 information
    genomes_frequency_info = hl.filter_intervals(
        hl.read_table(release_sites("genomes").path),
        interval_list,
    )

    exomes_frequency_info = hl.filter_intervals(
        hl.read_table(release_sites("genomes").path),
        interval_list,
    )

    # Create a new struct with v2, v4 exomes, and v4 genomes frequency annotation
    lifted_three_genes_combined_freq = lifted_three_genes.annotate(
        combined_frequency=hl.struct(
            v2_exomes_freq=lifted_three_genes.freq,
            v4_exomes_freq=exomes_frequency_info[
                lifted_three_genes.locus, lifted_three_genes.alleles
            ].freq,
            v4_genomes_freq=genomes_frequency_info[
                lifted_three_genes.locus, lifted_three_genes.alleles
            ].freq,
        )
    )

    # Annotate table with the names of the relevant genes - useful for clinicians possibly
    lifted_three_genes_combined_freq = lifted_three_genes_combined_freq.annotate(
        gene_symbol="No Gene Name"
    )
    for pair in [
        ["KCNE1", kcne1_interval],
        ["CBS", cbs_interval],
        ["CRYAA", cryaa_interval],
    ]:
        gene_interval = pair[1]
        new_gene_symbol = pair[0]
        lifted_three_genes_combined_freq = lifted_three_genes_combined_freq.annotate(
            gene_symbol=hl.if_else(
                gene_interval.contains(lifted_three_genes_combined_freq.locus),
                new_gene_symbol,
                lifted_three_genes_combined_freq.gene_symbol,
            )
        )

    # Write output to created resource with release root
    lifted_three_genes_combined_freq = lifted_three_genes_combined_freq.checkpoint(
        get_custom_liftover_three_genes_path(), overwrite=args.overwrite
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Option to overwrite existing custom liftover table",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
