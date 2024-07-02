"""Script to combine gVCFs into a single VDS."""

import argparse
import logging
from typing import List, Optional

import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_vds")
logger.setLevel(logging.INFO)


def create_vds(
    gvcfs: List[str],
    output_path: str,
    temp_path: str = "gs://gnomad-tmp-4day/combiner/",
    save_path: str = "gs://gnomad-tmp-4day/combiner/",
    use_genome_default_intervals: bool = False,
    use_exome_default_intervals: bool = False,
    intervals: Optional[List[str]] = None,
) -> hl.VariantDataset:
    """
    Combine gVCFs into a single VDS.

    :param List[str] gvcfs: List of gVCF paths.
    :param str output_path: Path to write output VDS.
    :param str temp_path: Path to write temporary files.
    :param str save_path: Path to write save path on failure.
    :param bool use_genome_default_intervals: Use the default genome intervals.
    :param bool use_exome_default_intervals: Use the default exome intervals.
    :param List[str] intervals: List of intervals to use.
    :return: Combined VDS.
    """
    logger.info("Combining %s gVCFs into a single VDS", len(gvcfs))
    combiner = hl.vds.new_combiner(
        temp_path=temp_path,
        output_path=output_path,
        save_path=save_path,
        gvcf_paths=gvcfs,
        use_genome_default_intervals=use_genome_default_intervals,
        use_exome_default_intervals=use_exome_default_intervals,
        intervals=intervals,
    )
    combiner.run()
    vds = hl.vds.read_vds(output_path)
    return vds


def main(args):
    """Create VDS from table of GVCFs paths and write to output path."""
    hl.init(
        log="/tmp/gvcf_combiner.log", default_reference="GRCh38", tmp_dir=args.temp_path
    )

    if args.use_genome_default_intervals and args.use_exome_default_intervals:
        raise ValueError(
            "Cannot use both genome and exome default intervals. Please choose one."
        )

    gvcfs = args.gvcfs
    output_path = args.output_path
    temp_path = args.temp_path
    save_path = args.temp_path

    logger.info("Reading gVCFs from %s", gvcfs)
    gvcfs = hl.import_table(gvcfs).collect()
    vds = create_vds(
        gvcfs=gvcfs,
        output_path=output_path,
        temp_path=temp_path,
        save_path=save_path,
        use_genome_default_intervals=args.use_genome_default_intervals,
        use_exome_default_intervals=args.use_exome_default_intervals,
        intervals=args.intervals,
    )
    logger.info("VDS variant data schema: %s", vds.variant_data.describe())
    logger.info("VDS sample data schema: %s", vds.sample_data.describe())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gvcfs", help="Path to file with only gVCF paths", required=True
    )
    parser.add_argument(
        "--output-path", help="Path to write output VDS to", required=True
    )
    parser.add_argument(
        "--temp-path",
        help="Path to write temporary files, including save path on failure",
        default="gs://gnomad-tmp-4day/combiner/",
    )
    parser.add_argument(
        "--use-genome-default-intervals",
        help="Use the default genome intervals",
        action="store_true",
    )
    parser.add_argument(
        "--use-exome-default-intervals",
        help="Use the default exome intervals",
        action="store_true",
    )
    parser.add_argument(
        "--intervals", help="Optional list of intervals to use", nargs="+"
    )
    args = parser.parse_args()
    main(args)
