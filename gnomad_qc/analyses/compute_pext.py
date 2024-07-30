"""Script to get release pext annotation files."""

import argparse
import logging
from typing import Tuple

import hail as hl
from gnomad.resources import grch37 as gnomad_grch37
from gnomad.resources import grch38 as gnomad_grch38
from gnomad.utils.slack import slack_notifications
from gnomad.utils.transcript_annotation import (
    TISSUES_TO_EXCLUDE,
    clean_tissue_name_for_browser,
    create_tx_annotation_by_region,
    perform_tx_annotation_pipeline,
    summarize_transcript_expression,
)
from hail.utils import new_temp_file

from gnomad_qc.analyses.resources import get_pext
from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("pext")
logger.setLevel(logging.INFO)


def filter_to_test(
    ht: hl.Table,
    tx_mt: hl.MatrixTable,
) -> Tuple[hl.Table, hl.MatrixTable]:
    """
    Filter `ht` and `tx_mt` to a test set of genes.

    Genes to filter to:

        - APOE
        - HTR2A
        - C1orf170/PERM1
        - TNFRSF25
        - TNFRSF18

    :param ht: VEP context HT to filter.
    :param tx_mt: Transcript MT to filter.
    :return: Filtered HT and MT.
    """
    test_intervals = [
        hl.parse_locus_interval("19:45409011-45412650"),  # ENSG00000130203, APOE
        hl.parse_locus_interval("13:47405685-47471169"),  # ENSG00000102468, HTR2A
        hl.parse_locus_interval("1:910579-917497"),  # ENSG00000187642, C1orf170/PERM1
        hl.parse_locus_interval("1:6521211-6526255"),  # ENSG00000215788, TNFRSF25
        hl.parse_locus_interval("1:1138888-1142071"),  # ENSG00000186891, TNFRSF18
    ]
    test_genes = [
        "ENSG00000130203",
        "ENSG00000102468",
        "ENSG00000187642",
        "ENSG00000215788",
        "ENSG00000186891",
    ]

    logger.info(
        "Filtering VEP context Table to the following genes: %s...",
        ", ".join(test_genes),
    )
    ht = hl.filter_intervals(ht, test_intervals)

    tx_mt = tx_mt.filter_rows(hl.literal(test_genes).contains(tx_mt.gene_id))
    logger.info(
        "Filter GTEx RSEM MatrixTable to %d transcripts in %d genes...",
        tx_mt.count()[0],
        len(test_genes),
    )

    return ht, tx_mt


def get_pipeline_resources(
    test: bool,
    gtex_version: str,
    vep_version: str,
    overwrite: bool,
) -> PipelineResourceCollection:
    """
    Get pipeline resources for pext annotation pipeline.

    :param test: Whether to get resources for a test run.
    :param gtex_version: GTEx version to get resources for.
    :param vep_version: VEP version to get resources for.
    :param overwrite: Whether to overwrite existing data.
    :return: Pipeline resource collection for pext annotation pipeline.
    """
    if gtex_version == "v7":
        gnomad_res = gnomad_grch37
    elif gtex_version == "v10":
        gnomad_res = gnomad_grch38
    else:
        raise ValueError("Invalid GTEx version")

    # Initialize pipeline resource collection.
    per_sample_pipeline = PipelineResourceCollection(
        pipeline_name="compute_pext",
        pipeline_resources={
            "context HT": {"context_ht": gnomad_res.vep_context.versions[vep_version]},
            "GTEx RSEM MT": {"tx_mt": gnomad_res.gtex_rsem.versions[gtex_version]},
        },
        overwrite=overwrite,
    )

    # Create resource collection for each step of the pipeline.
    annotation_level = PipelineStepResourceCollection(
        "--annotation-level",
        output_resources={
            "annotation_pext_ht": get_pext(
                "annotation_level",
                gtex_version,
                vep_version,
                test=test,
            ),
        },
    )
    base_level = PipelineStepResourceCollection(
        "--base-level",
        output_resources={
            "base_pext_ht": get_pext(
                "base_level",
                gtex_version,
                vep_version,
                test=test,
            ),
        },
    )
    browser_ht = PipelineStepResourceCollection(
        "--browser-ht",
        pipeline_input_steps=[base_level],
        output_resources={
            "browser_pext_ht": get_pext(
                "browser",
                gtex_version,
                vep_version,
                test=test,
            ),
        },
    )

    # Add all steps to the pipeline resource collection.
    per_sample_pipeline.add_steps(
        {
            "annotation_level": annotation_level,
            "base_level": base_level,
            "browser_ht": browser_ht,
        }
    )

    return per_sample_pipeline


def main(args):
    """Script to get release pext annotation files."""
    test = args.test
    gtex_version = args.gtex_version
    vep_version = args.vep_version
    overwrite = args.overwrite

    hl.init(
        log="/pext.log",
        default_reference="GRCh37" if gtex_version == "v7" else "GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    pext_res = get_pipeline_resources(test, gtex_version, vep_version, overwrite)
    ht = pext_res.context_ht.ht().naive_coalesce(5000)
    tx_mt = pext_res.tx_mt.mt()

    if test:
        ht, tx_mt = filter_to_test(ht, tx_mt)

    logger.info(
        "Summarizing transcript expression for GTEx version %s...",
        gtex_version,
    )
    tx_ht = summarize_transcript_expression(tx_mt).checkpoint(
        new_temp_file(f"{f'test_' if test else ''}tx_ht", "ht"), overwrite=overwrite
    )

    if args.annotation_level:
        logger.info("Creating annotation-level pext scores...")
        res = pext_res.annotation_level
        res.check_resource_existence()

        perform_tx_annotation_pipeline(
            ht, tx_ht, tissues_to_exclude_from_mean=TISSUES_TO_EXCLUDE[gtex_version]
        ).write(res.annotation_pext_ht.path, overwrite=overwrite)

    if args.base_level:
        logger.info("Creating base-level pext scores...")
        res = pext_res.base_level
        res.check_resource_existence()

        perform_tx_annotation_pipeline(
            ht,
            tx_ht,
            additional_group_by=["gene_symbol"],
            tissues_to_exclude_from_mean=TISSUES_TO_EXCLUDE[gtex_version],
        ).write(res.base_pext_ht.path, overwrite=overwrite)

    if args.browser_ht:
        logger.info("Creating pext Table for browser loading...")
        res = pext_res.browser_ht
        res.check_resource_existence()
        base_pext_ht = res.base_pext_ht.ht()

        tissue_map = {
            t: clean_tissue_name_for_browser(t) for t in hl.eval(base_pext_ht.tissues)
        }
        base_pext_ht = base_pext_ht.annotate(
            **{k: base_pext_ht[k].expression_proportion for k in tissue_map}
        )
        base_pext_ht = base_pext_ht.rename(tissue_map)
        base_pext_ht = base_pext_ht.annotate_globals(
            tissues=base_pext_ht.tissues.map(lambda x: hl.literal(tissue_map).get(x))
        )
        create_tx_annotation_by_region(base_pext_ht).write(
            res.browser_pext_ht.path, overwrite=overwrite
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script creates the pext annotation files for release."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--slack-token", help="Slack API token.")
    parser.add_argument(
        "--test", help="Test on a small number of genes.", action="store_true"
    )
    parser.add_argument(
        "--gtex-version",
        help="GTEx version to use for pext scores.",
        default="v10",
        choices=["v7", "v10"],
    )
    parser.add_argument(
        "--vep-version",
        help="VEP version to use for pext scores.",
        default="105",
        choices=["85", "95", "101", "105"],
    )
    parser.add_argument(
        "--annotation-level",
        help="Get annotation-level pext scores.",
        action="store_true",
    )
    parser.add_argument(
        "--base-level", help="Get base-level pext scores.", action="store_true"
    )
    parser.add_argument(
        "--browser-ht",
        help="Get pext Table for browser loading.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
