"""Script to combine GVCFs into a single VDS."""

import argparse
import logging

import hail as hl
from gnomad.utils.file_utils import create_vds

from gnomad_qc.v5.resources.basics import dragen_tgp_vds

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_vds")
logger.setLevel(logging.INFO)


def main(args):
    """Create VDS from text file of GVCFs paths and write to output path."""
    hl.init(
        log="/tmp/gvcf_combiner.log", default_reference="GRCh38", tmp_dir=args.temp_path
    )
    gvcfs = args.gvcfs
    temp_path = args.temp_path
    save_path = temp_path + "combiner_plan.json"

    if args.output_path:
        output_path = args.output_path
    else:
        logger.info(
            "No output path provided, using DRAGEN TGP VDS path as output path."
        )
        output_path = dragen_tgp_vds.path

    vds = create_vds(
        gvcfs=gvcfs,
        output_path=output_path,
        temp_path=temp_path,
        save_path=save_path,
        use_genome_default_intervals=args.use_genome_default_intervals,
        intervals=args.intervals,
        gvcf_batch_size=args.gvcf_batch_size,
    )
    logger.info("VDS variant data schema: %s", vds.variant_data.describe())
    logger.info("VDS reference data schema: %s", vds.reference_data.describe())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--gvcfs",
        help="Path to text file with only gVCF paths.",
        required=True,
        type=str,
    )
    parser.add_argument("--output-path", help="Path to write output VDS to.", type=str)
    parser.add_argument(
        "--temp-path",
        help="Directory path to write temporary files, including save path on failure.",
        default="gs://gnomad-tmp-4day/combiner/",
        type=str,
    )
    parser.add_argument(
        "--use-genome-default-intervals",
        help="Use hail's default genome intervals.",
        action="store_true",
    )
    parser.add_argument(
        "--intervals",
        help="Text file of intervals to use for VDS creation.",
        type=str,
    )
    parser.add_argument(
        "--gvcf-batch-size",
        help="Number of GVCFs to combine into a Variant Dataset at once.",
        type=int,
    )
    args = parser.parse_args()
    main(args)
