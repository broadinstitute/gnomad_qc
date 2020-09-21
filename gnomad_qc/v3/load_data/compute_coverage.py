import argparse
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import coverage, coverage_tsv_path, CURRENT_GENOME_COVERAGE_RELEASE
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres

from gnomad.utils.reference_genome import get_reference_ht
from gnomad.utils.sparse_mt import compute_coverage_stats

from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta


def main(args):
    hl.init(default_reference='GRCh38')
    coverage_version = args.coverage_version if args.coverage_version else CURRENT_GENOME_COVERAGE_RELEASE
    logger = logging.getLogger("gnomad_qc.v3.load_data.compute_coverage")

    logger.warning("Last time this was run (July 2020), this script required high-mem machines.")

    if args.compute_coverage_ht:

        print("Building reference context HT")
        ref_ht = get_reference_ht(
            hl.get_reference('GRCh38'),
            contigs=[f'chr{x}' for x in range(1,23)] + ['chrX', 'chrY'],
            excluded_intervals=telomeres_and_centromeres.ht().interval.collect()
        )
        ref_ht = ref_ht.checkpoint("gs://gnomad-tmp/ref_context.ht", overwrite=True)
        logger.info("Done building reference context HT")

        mt = get_gnomad_v3_mt()
        mt = mt.filter_cols(meta.ht()[mt.col_key].release)

        coverage_ht = compute_coverage_stats(
            mt,
            ref_ht
        )

        coverage_ht = coverage_ht.checkpoint('gs://gnomad-tmp/gnomad.genomes_v3.coverage.summary.ht', overwrite=True)
        coverage_ht = coverage_ht.naive_coalesce(5000)

        coverage_ht.write(coverage('genomes').versions[coverage_version].path, overwrite=args.overwrite)

    if args.export_coverage:
        ht = coverage('genomes').versions[coverage_version].ht()
        if 'count_array' in ht.row_value: # Note that count_array isn't computed any more, so this is v3.0-specific
            ht = ht.drop('count_array')
        ht.export(coverage_tsv_path('genomes', coverage_version))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrites existing hail tables', action='store_true')
    parser.add_argument('--compute_coverage_ht', help='Computes the coverage HT', action='store_true')
    parser.add_argument('--export_coverage', help='Exports coverage TSV file', action='store_true')
    parser.add_argument('--coverage_version', help="Specifies coverage version to read/write. If not set, gnomad.resources.grch38.gnomad.CURRENT_GENOME_COVERAGE_RELEASE is used.")

    main(parser.parse_args())
