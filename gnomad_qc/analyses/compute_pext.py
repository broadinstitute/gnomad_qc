"""Script to get release pext annotation files."""

import argparse
import logging
from typing import Tuple

import hail as hl
from gnomad.resources import grch37 as gnomad_grch37
from gnomad.resources import grch38 as gnomad_grch38
from gnomad.resources.grch37.gnomad import public_release as gnomad_release_grch37
from gnomad.resources.grch38.gnomad import public_release as gnomad_release_grch38
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.utils.slack import slack_notifications
from gnomad.utils.transcript_annotation import (
    clean_tissue_name_for_browser,
    create_tx_annotation_by_region,
    get_tissues_to_exclude,
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
    test_genes = {
        "ENSG00000130203": {  # APOE
            "GRCh37": "19:45409011-45412650",
            "GRCh38": "chr19:44905791-44909393",
        },
        "ENSG00000102468": {  # HTR2A
            "GRCh37": "13:47405685-47471169",
            "GRCh38": "chr13:46831546-46897076",
        },
        "ENSG00000187642": {  # C1orf170/PERM1
            "GRCh37": "1:910579-917497",
            "GRCh38": "chr1:975204-982093",
        },
        "ENSG00000215788": {  # TNFRSF25
            "GRCh37": "1:6521211-6526255",
            "GRCh38": "chr1:6460786-6466175",
        },
        "ENSG00000186891": {  # TNFRSF18
            "GRCh37": "1:1138888-1142071",
            "GRCh38": "chr1:1203508-1206592",
        },
    }

    logger.info(
        "Filtering VEP context Table to the following genes: %s...",
        ", ".join(list(test_genes.keys())),
    )

    # Get the reference genome build from the HT.
    build = get_reference_genome(ht.locus).name

    test_intervals = [hl.parse_locus_interval(g[build]) for g in test_genes.values()]
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
    overwrite: bool,
) -> PipelineResourceCollection:
    """
    Get pipeline resources for pext annotation pipeline.

    :param test: Whether to get resources for a test run.
    :param gtex_version: GTEx version to get resources for. Must be one of 'v7' or 'v10'.
    :param overwrite: Whether to overwrite existing data.
    :return: Pipeline resource collection for pext annotation pipeline.
    """
    if gtex_version == "v7":
        gnomad_res = gnomad_grch37
        vep_version = "85"
        gnomad_exomes = gnomad_release_grch37("exomes").versions["2.1.1"]
        gnomad_genomes = gnomad_release_grch37("genomes").versions["2.1.1"]

        # TODO: Remove this once the GTEx v10 MT public resource is available.
        tx_mt = gnomad_res.gtex_rsem.versions[gtex_version]
    elif gtex_version == "v10":
        gnomad_res = gnomad_grch38
        vep_version = "105"
        gnomad_exomes = gnomad_release_grch38("exomes").versions["4.1"]
        gnomad_genomes = gnomad_release_grch38("genomes").versions["4.1"]

        # TODO: Update to use the GTEx v10 MT public resource once it is available.
        from gnomad.resources.resource_utils import MatrixTableResource

        tx_mt = MatrixTableResource(
            "gs://gnomad/resources/grch38/gtex/gtex_rsem_v10.mt"
        )
    else:
        raise ValueError("Invalid GTEx version")

    # Initialize pipeline resource collection.
    per_sample_pipeline = PipelineResourceCollection(
        pipeline_name="compute_pext",
        pipeline_resources={
            "context HT": {"context_ht": gnomad_res.vep_context.versions[vep_version]},
            # TODO: Change to use gnomad_res.gtex_rsem.versions[gtex_version], once
            #  the public version is available.
            "GTEx RSEM MT": {"tx_mt": tx_mt},
        },
        overwrite=overwrite,
    )

    # Create resource collection for each step of the pipeline.
    annotation_level = PipelineStepResourceCollection(
        "--annotation-level",
        output_resources={
            "annotation_level_pext_ht": get_pext(
                "annotation_level", gtex_version, test=test
            ),
        },
    )
    base_level = PipelineStepResourceCollection(
        "--base-level",
        output_resources={
            "base_pext_ht": get_pext("base_level", gtex_version, test=test),
        },
    )
    gnomad_exomes = PipelineStepResourceCollection(
        "--gnomad-exomes-ht",
        input_resources={
            "gnomad exomes HT": {"gnomad_exomes_ht": gnomad_exomes},
        },
        output_resources={
            "gnomad_exomes_pext_ht": get_pext("exomes", gtex_version, test=test),
        },
    )
    gnomad_genomes = PipelineStepResourceCollection(
        "--gnomad-genomes-ht",
        input_resources={
            "gnomad genomes HT": {"gnomad_genomes_ht": gnomad_genomes},
        },
        output_resources={
            "gnomad_genomes_pext_ht": get_pext("genomes", gtex_version, test=test),
        },
    )
    browser_ht = PipelineStepResourceCollection(
        "--browser-ht",
        pipeline_input_steps=[base_level],
        output_resources={
            "browser_pext_ht": get_pext("browser", gtex_version, test=test),
        },
    )
    export_tsv = PipelineStepResourceCollection(
        "--export-tsv",
        pipeline_input_steps=[annotation_level, base_level],
        output_resources={
            "pext_base_tsv": get_pext(
                "base_level", gtex_version, suffix="tsv.gz", test=test
            ),
            "pext_annotation_tsv": get_pext(
                "annotation_level", gtex_version, suffix="tsv.gz", test=test
            ),
        },
    )

    # Add all steps to the pipeline resource collection.
    per_sample_pipeline.add_steps(
        {
            "annotation_level": annotation_level,
            "base_level": base_level,
            "gnomad_exomes": gnomad_exomes,
            "gnomad_genomes": gnomad_genomes,
            "browser_ht": browser_ht,
            "export_tsv": export_tsv,
        }
    )

    return per_sample_pipeline


def main(args):
    """Script to get release pext annotation files."""
    test = args.test
    gtex_version = args.gtex_version
    overwrite = args.overwrite
    min_samples = args.min_samples
    browser_min_samples = args.browser_min_samples
    if browser_min_samples < min_samples:
        raise ValueError(
            "The --browser-min-samples threshold must be greater than or equal to the "
            "--min-samples threshold."
        )

    annotation_level = ["annotation_level"] if args.annotation_level else []
    annotation_level += ["gnomad_exomes"] if args.gnomad_exomes_ht else []
    annotation_level += ["gnomad_genomes"] if args.gnomad_genomes_ht else []

    hl.init(
        log="/pext.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh37" if gtex_version == "v7" else "GRCh38")
    hl._set_flags(use_new_shuffle="1")

    pext_res = get_pipeline_resources(test, gtex_version, overwrite)

    # Reduce the number of partitions in the context HT to avoid long run times.
    # The context HTs are >30K partitions, which led to long run times for the pipeline,
    # reducing the number of partitions to 5K helped reduce the run time.
    ht = pext_res.context_ht.ht().naive_coalesce(5000)
    tx_mt = pext_res.tx_mt.mt()

    # Get tissues to completely exclude from the pext calculations. These are tissues
    # that have less than `min_samples` samples and should not be included in the
    # individual tissue pext values or in the overall pext.
    tissues_to_exclude = get_tissues_to_exclude(
        tx_mt, reproductive=False, cell_lines=False, min_samples=min_samples
    )

    # Get tissues to exclude from the overall mean pext calculation. This includes
    # reproductive tissues and cell lines, as well as any tissues that have less than
    # `min_samples` samples (`tissues_to_exclude`).
    tissues_to_exclude_from_mean = (
        get_tissues_to_exclude(tx_mt, min_samples=None) + tissues_to_exclude
    )

    # If `keep_tissues_under_min_samples` is set, tissues that have less than
    # `min_samples` samples will be excluded from the overall mean pext calculation,
    # but will still be included in the individual tissue pext values.
    if args.keep_tissues_under_min_samples:
        tissues_to_exclude = []

    if test:
        ht, tx_mt = filter_to_test(ht, tx_mt)

    if annotation_level or args.base_level:
        logger.info(
            "Summarizing transcript expression for GTEx version %s...",
            gtex_version,
        )
        tx_ht = summarize_transcript_expression(tx_mt).checkpoint(
            new_temp_file(f"{f'test_' if test else ''}tx_ht", "ht"), overwrite=overwrite
        )

    if annotation_level:
        logger.info("Creating annotation-level pext scores...")

        for dataset in annotation_level:
            dataset_name = "context" if dataset == "annotation_level" else dataset
            logger.info("Annotating the %s HT...", dataset_name.replace("_", " "))
            res = getattr(pext_res, dataset)
            res.check_resource_existence()

            perform_tx_annotation_pipeline(
                (
                    ht
                    if dataset == "annotation_level"
                    else getattr(res, f"{dataset}_ht").ht()
                ),
                tx_ht,
                tissues_to_exclude=tissues_to_exclude,
                tissues_to_exclude_from_mean=tissues_to_exclude_from_mean,
            ).write(getattr(res, f"{dataset}_pext_ht").path, overwrite=overwrite)

    if args.base_level:
        logger.info("Creating base-level pext scores...")
        res = pext_res.base_level
        res.check_resource_existence()

        perform_tx_annotation_pipeline(
            ht,
            tx_ht,
            additional_group_by=["gene_symbol"],
            tissues_to_exclude=tissues_to_exclude,
            tissues_to_exclude_from_mean=tissues_to_exclude_from_mean,
        ).write(res.base_pext_ht.path, overwrite=overwrite)

    if args.browser_ht:
        logger.info("Creating pext Table for browser loading...")
        res = pext_res.browser_ht
        res.check_resource_existence()
        base_pext_ht = res.base_pext_ht.ht()

        tissue_map = {
            t: clean_tissue_name_for_browser(t) for t in hl.eval(base_pext_ht.tissues)
        }

        browser_tissues_to_exclude = tissues_to_exclude
        if browser_min_samples > min_samples:
            # Get tissues that have less than `browser_min_samples` samples and should
            # not be displayed as a tissue pext track on the browser.
            browser_tissues_to_exclude = get_tissues_to_exclude(
                tx_mt,
                reproductive=False,
                cell_lines=False,
                min_samples=browser_min_samples,
            )
            logger.info(
                "Excluding tissues with less than %d samples from the browser: %s",
                browser_min_samples,
                ", ".join(
                    list(set(browser_tissues_to_exclude) - set(tissues_to_exclude))
                ),
            )

        base_pext_ht = base_pext_ht.annotate(
            **{k: base_pext_ht[k].expression_proportion for k in tissue_map}
        )
        base_pext_ht = base_pext_ht.rename(tissue_map)
        base_pext_ht = base_pext_ht.annotate_globals(
            tissues=base_pext_ht.tissues.filter(
                lambda x: ~hl.literal(browser_tissues_to_exclude).contains(x)
            ).map(lambda x: hl.literal(tissue_map).get(x))
        )
        create_tx_annotation_by_region(base_pext_ht).write(
            res.browser_pext_ht.path, overwrite=overwrite
        )

    if args.export_tsv:
        logger.info("Exporting pext Hail Tables to TSVs...")
        res = pext_res.export_tsv
        res.check_resource_existence()

        hts = [
            (res.base_pext_ht.ht(), res.pext_base_tsv),
            (res.annotation_level_pext_ht.ht(), res.pext_annotation_tsv),
        ]

        for ht, tsv_path in hts:
            ht = ht.annotate(
                **{k: ht[k].expression_proportion for k in hl.eval(ht.tissues)}
            )
            ht.export(tsv_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script creates the pext annotation files for release."
    )
    parser.add_argument("--overwrite", help="Overwrite data.", action="store_true")
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
        "--min-samples",
        help=(
            "Minimum number of samples to include a tissue in the overall mean pext. "
            "If keep-tissues-under-min-samples is not set, these tissues will also be "
            "excluded from the individual tissue pext values."
        ),
        type=int,
        default=50,
    )
    parser.add_argument(
        "--keep-tissues-under-min-samples",
        help=(
            "Keep tissues under the min-samples threshold in the individual tissue "
            "pext values."
        ),
        action="store_true",
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
        "--gnomad-exomes-ht",
        help="Get pext Table for gnomAD exomes sites.",
        action="store_true",
    )
    parser.add_argument(
        "--gnomad-genomes-ht",
        help="Get pext Table for gnomAD genomes sites.",
        action="store_true",
    )
    parser.add_argument(
        "--browser-ht",
        help="Get pext Table for browser loading.",
        action="store_true",
    )
    parser.add_argument(
        "--browser-min-samples",
        help=(
            "Minimum number of samples to include a tissue in the individual tissue "
            "pext values displayed as a track on the browser."
        ),
        type=int,
        default=100,
    )
    parser.add_argument(
        "--export-tsv",
        help="Export pext Hail Tables to TSVs.",
        action="store_true",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
