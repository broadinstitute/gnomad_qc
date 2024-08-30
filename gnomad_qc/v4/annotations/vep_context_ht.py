"""Script to add VEP annotations to the GRCh38 context Table."""

import argparse
import logging
import os

import hail as hl
from gnomad.resources.config import (
    GnomadPublicResourceSource,
    gnomad_public_resource_configuration,
)
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import VEP_CONFIG_PATH, get_vep_help

from gnomad_qc.slack_creds import slack_token

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def main(args):
    """Add VEP annotations to the GRCh38 context Table."""
    hl.init(default_reference="GRCh38", log="/vep.log", tmp_dir="gs://gnomad-tmp-4day")

    vep_version = args.vep_version
    ht = vep_context.versions[args.input_context_ht_vep_version].ht()
    if args.test:
        ht = ht._filter_partitions(range(2))

    # Drop the 'vep' annotation which will be replaced by the new vep version.
    # In the vep 101 context HT we stored the 'vep_proc_id' annotation, so we need to
    # drop it if it is present.
    if "vep_proc_id" in list(ht.row):
        ht = ht.drop("vep", "vep_proc_id")
    else:
        ht = ht.drop("vep")

    ht = hl.vep(ht, VEP_CONFIG_PATH)

    vep_help = get_vep_help(VEP_CONFIG_PATH)
    with hl.hadoop_open(VEP_CONFIG_PATH) as vep_config_file:
        vep_config = vep_config_file.read()

    ht = ht.annotate_globals(
        version=f"v{vep_version}", vep_help=vep_help, vep_config=vep_config
    )

    # Switch gnomAD public resource source to get path in
    # gs://gnomad-public-requester-pays to write to.
    gnomad_public_resource_configuration.source = GnomadPublicResourceSource.GNOMAD
    out_vep_path = vep_context.versions[vep_version].path
    if args.test:
        ht.write(
            "gs://gnomad-tmp/resources/context/" + os.path.split(out_vep_path)[-1],
            overwrite=args.overwrite,
        )
    else:
        ht.write(out_vep_path, overwrite=args.overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-context-ht-vep-version",
        help=(
            "Version of VEPed context Table to annotate onto. Drops 'vep' and "
            "'vep_proc_id' annotations before re-annotating."
        ),
        action="store_true",
        default="101",
    )
    parser.add_argument(
        "--vep-version",
        help="Version of VEP the context Table is being annotated with.",
        action="store_true",
        default="105",
    )
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the MT.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--overwrite", help="Overwrite data if it exists.", action="store_true"
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
