"""Script to get base-level pext annotations for gnomAD v2.1.1."""

import argparse
import logging

import hail as hl
from gnomad.resources.grch37.reference_data import gtex_rsem, vep_context
from gnomad.utils.slack import slack_notifications
from gnomad.utils.transcript_annotation import (
    get_max_pext_per_gene,
    get_worst_csq_per_variant,
    perform_tx_annotation_pipeline,
    summarize_transcript_expression,
)
from hail.utils import new_temp_file

from gnomad_qc.slack_creds import slack_token

# List of tissues to exclude from the analysis, including reproductive tissues,
# cell lines and any tissue with less than 100 samples in GTEx v7.
TISSUES_TO_EXCLUDE = [
    "Bladder",
    "Brain_Spinalcord_cervicalc_1",
    "Brain_Substantianigra",
    "Cells_EBV_transformedlymphocytes",
    "Cells_Transformedfibroblasts",
    "Cervix_Ectocervix",
    "Cervix_Endocervix",
    "FallopianTube",
    "Kidney_Cortex",
    "MinorSalivaryGland",
    "Ovary",
    "Prostate",
    "Testis",
    "Uterus",
    "Vagina",
]

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("pext_v2_pipeline")
logger.setLevel(logging.INFO)


def main(args):
    """
    Script to get pext annotations for gnomAD v2.1.1.

    :param str args: Command-line arguments passed to the script.
    """
    test = args.test
    max_pext = args.max_pext
    overwrite = args.overwrite
    annotation_level = args.annotation_level
    base_level = args.base_level
    worst_csq_per_variant = args.worst_csq_per_variant

    hl.init(
        log="/v2_pext.log",
        default_reference="GRCh37",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    ht = vep_context.ht()
    tx = gtex_rsem.mt()

    if test:
        test_intervals = [
            hl.parse_locus_interval("19:45409011-45412650"),  # ENSG00000130203, APOE
            hl.parse_locus_interval("13:47405685-47471169"),  # ENSG00000102468, HTR2A
            hl.parse_locus_interval(
                "1:910579-917497"
            ),  # ENSG00000187642, C1orf170/PERM1
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
            "Filtering vep context_ht to the following genes: %s...",
            ", ".join(test_genes),
        )
        ht = hl.filter_intervals(ht, test_intervals)
        logger.info("Filteringg vep context_ht to %d variants...", ht.count())
        tx = tx.filter_rows(hl.literal(test_genes).contains(tx.gene_id))
        logger.info(
            "Filter tx_ht to %d transcripts in %d genes...",
            tx.count()[0],
            len(test_genes),
        )

    tx = summarize_transcript_expression(tx)

    # TODO: Whether to annotate the result Table back to original variant HT.
    if annotation_level:
        ht = perform_tx_annotation_pipeline(ht, tx)
        ht = ht.checkpoint(
            new_temp_file(f"{f'test_' if test else ''}pext_annotation-level_v2", "ht")
        )
        if worst_csq_per_variant:
            ht = ht.collect_by_key("tx_annotation")
            worst_ht = get_worst_csq_per_variant(ht)
            worst_ht.checkpoint(
                new_temp_file(f"{f'test_' if test else ''}worst_pext_per_variant", "ht")
            )
        if max_pext:
            ht = ht.collect_by_key("tx_annotation")
            max_ht = get_max_pext_per_gene(ht, tissues_to_filter=TISSUES_TO_EXCLUDE)
            max_ht.checkpoint(
                new_temp_file(
                    f"{f'test_' if test else ''}max_pext_per_gene if test", "ht"
                ),
            )

    if base_level:
        ht = perform_tx_annotation_pipeline(
            ht,
            tx,
            additional_group_by=["gene_symbol"],
        )
        ht = ht.checkpoint(
            new_temp_file(f"{f'test_' if test else ''}pext_base-level_v2", "ht")
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script creates the v2 pext scores on GTEx v7 RSEM."
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
        "--annotation-level",
        help="Get annotation-level pext scores.",
        action="store_true",
    )
    parser.add_argument(
        "--base-level", help="Get base-level pext scores.", action="store_true"
    )
    parser.add_argument(
        "--worst-csq-per-variant", help="Pull worst pext per gene.", action="store_true"
    )
    parser.add_argument(
        "--max-pext", help="Get maximum pext per gene.", action="store_true"
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
