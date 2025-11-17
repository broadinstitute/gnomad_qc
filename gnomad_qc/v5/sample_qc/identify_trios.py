"""Script to identify trios from relatedness data and filter based on Mendel errors and de novos."""

import argparse
import logging

import hail as hl
from gnomad.sample_qc.relatedness import (
    get_duplicated_samples,
    get_duplicated_samples_ht,
)

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.resources.basics import get_logging_path
from gnomad_qc.v5.resources.sample_qc import (
    duplicates,
    finalized_outlier_filtering,
    relatedness,
    sample_rankings,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("identify_trios")
logger.setLevel(logging.INFO)


def filter_relatedness_ht(ht: hl.Table, filter_ht: hl.Table) -> hl.Table:
    """
    Filter relatedness Table to only include pairs of samples that are both AoU genomes and not QC-filtered.

    :param ht: Relatedness Table.
    :param filter_ht: Outlier filtering Table.
    :return: Filtered relatedness Table.
    """
    ht = ht.filter(
        ((ht.i.data_type == "genomes") & (ht.i.project == "aou"))
        & ((ht.j.data_type == "genomes") & (ht.j.project == "aou"))
    )
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)

    # Remove all pairs with a QC-filtered sample
    ht = ht.filter(
        filter_ht[ht.i].outlier_filtered | filter_ht[ht.j].outlier_filtered,
        keep=False,
    )

    return ht


def main(args):
    """Identify trios and filter based on Mendel errors and de novos."""
    hl.init(
        log="/identify_trios.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test = args.test

    try:

        rel_ht = relatedness().ht()
        rel_ht = filter_relatedness_ht(rel_ht, finalized_outlier_filtering().ht())

        if args.identify_duplicates:
            logger.info("Selecting best duplicate per duplicated sample set...")
            dup_ht_path = duplicates().path
            check_resource_existence(
                output_step_resources={"duplicates_ht": [dup_ht_path]},
                overwrite=overwrite,
            )
            ht = get_duplicated_samples_ht(
                get_duplicated_samples(rel_ht),
                sample_rankings(release=True).ht(),
            )
            ht.write(dup_ht_path, overwrite=overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("identify_trios", environment=args.environment))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing data.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Run mendel errors on only five partitions of the MT.",
        action="store_true",
    )

    identify_dup_args = parser.add_argument_group("Duplicate identification")
    identify_dup_args.add_argument(
        "--identify-duplicates",
        help=(
            "Create a table with duplicate samples indicating which one is the best to"
            " use based on the ranking of all samples after sample QC metric outlier "
            "filtering (ranking used to determine related samples to drop for the "
            "release)."
        ),
        action="store_true",
    )

    inter_fam_args = parser.add_argument_group("Pedigree inference")
    inter_fam_args.add_argument(
        "--infer-families",
        help=(
            "Infer families and trios using the relatedness Table and duplicate Table."
        ),
        action="store_true",
    )

    fake_ped_args = parser.add_argument_group("Fake Pedigree creation")
    fake_ped_args.add_argument(
        "--create-fake-pedigree",
        help=(
            "Create a fake Pedigree from unrelated samples in the data for comparison "
            "to the inferred Pedigree."
        ),
        action="store_true",
    )
    fake_ped_args.add_argument(
        "--fake-fam-prop",
        help=(
            "Number of fake trios to generate as a proportion of the total number of"
            " trios found in the data. Default is 0.1."
        ),
        default=0.1,
        type=float,
    )

    mendel_err_args = parser.add_argument_group("Mendel error calculation")
    mendel_err_args.add_argument(
        "--run-mendel-errors",
        help="Calculate mendel errors for the inferred and fake Pedigrees on chr20.",
        action="store_true",
    )
    finalize_ped_args = parser.add_argument_group(
        "Pedigree filtering for final Pedigree generation"
    )
    finalize_ped_args.add_argument(
        "--finalize-ped",
        help=(
            "Create final families/trios ped files by excluding trios where the number "
            "of Mendel errors or de novos are outliers. Outliers can be defined as "
            "trios that have Mendel errors or de novos higher than those specified in "
            "--max-mendel and --max-de-novo respectively. They can also be defined as "
            "trios that have Mendel errors or de novos higher than --max-mendel-z or "
            "--max-de-novo-z standard deviations above the mean across inferred trios."
        ),
        action="store_true",
    )
    finalize_ped_args.add_argument(
        "--max-mendel-z",
        help=(
            "Max number of standard deviations above the mean Mendel errors across "
            "inferred trios to keep a trio. If flag is set, default is 3."
        ),
        nargs="?",
        const=3,
        type=int,
    )
    finalize_ped_args.add_argument(
        "--max-de-novo-z",
        help=(
            "Max number of standard deviations above the mean de novos across inferred "
            "trios to keep a trio. If flag is set, default is 3."
        ),
        nargs="?",
        const=3,
        type=int,
    )
    finalize_ped_args.add_argument(
        "--max-mendel",
        help=(
            "Maximum number of raw Mendel errors for real trios. If specified and "
            "--ukb-max-mendel is not, --max-mendel will be used for all samples. If "
            "both --max-mendel and --ukb-max-mendel are specified, --max-mendel will "
            "be used for non-UKB samples and --ukb-max-mendel will be used for UKB "
            "samples."
        ),
        type=int,
    )
    finalize_ped_args.add_argument(
        "--max-de-novo",
        help=(
            "Maximum number of raw de novo mutations for real trios. If specified and"
            " --ukb-max-de-novo is not, --max-de-novo will be used for all samples. If"
            " both --max-de-novo and --ukb-max-de-novo are specified, --max-de-novo"
            " will be used for non-UKB samples and --ukb-max-de-novo will be used for"
            " UKB samples."
        ),
        type=int,
    )
    finalize_ped_args.add_argument(
        "--seed",
        help=(
            "Random seed for choosing one random trio per family to keep after "
            "filtering."
        ),
        type=int,
        default=24,
    )

    dense_trio_mt_args = parser.add_argument_group("Create dense trio MT.")
    dense_trio_mt_args.add_argument(
        "--create-dense-trio-mt",
        help=("Create a dense MT for high quality trios."),
        action="store_true",
    )
    dense_trio_mt_args.add_argument(
        "--naive-coalesce-partitions",
        help=("Number of partitions to coalesce the VDS to."),
        type=int,
        default=5000,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
