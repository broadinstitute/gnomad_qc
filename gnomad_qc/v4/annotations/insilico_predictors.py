"""Script to generate Hail Tables with in silico predictors."""

import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def create_cadd_grch38_ht() -> hl.Table:
    """
    Create a Hail Table with CADD scores for GRCh38.

    .. note:: The CADD scores are from the following sources:
        - all SNVs: `cadd.v1.6.whole_genome_SNVs.tsv.bgz` (81G) downloaded from `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz`. It contains 8,812,917,339 SNVs.
        - gnomad 3.0 indels: `cadd.v1.6.indels.tsv.bgz` (1.1G) downloaded from `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz`. It contains 100,546,109 indels in gnomaD v3.0.
        - gnomad 3.1 indels: `CADD-indels-gnomad.3.1.ht` was run by William Phu with CADD v1.6 in 2020. It contains 166,122,720 indels in gnomAD v3.1.
        - gnomad 3.1 complex indels: `CADD-1.6-gnomad-complex-variants.ht` was run by William Phu with CADD v1.6 in 2020. It contains 2,307 complex variants that were not filtered out by Hail is.indels().
        - gnomAD v4 indels: `cadd.v1.6.gnomAD_v4_new_indels.tsv.bgz` (368M) was run by Qin He with CADD v1.6 in 2023. It contains 32,561,253 indel that are new in gnomAD v4.
        - 1,972,208 indels were duplicated in gnomAD v3.0 and v4.0, in gnomAD v3.1 and v4.0. However, CADD only generate a score per loci. We use collect_by_key() to keep only one of the duplicated indels.
    :return: Hail Table with CADD scores for GRCh38.
    """

    def load_cadd_raw(cadd_tsv) -> hl.Table:
        column_names = {
            "f0": "chr",
            "f1": "pos",
            "f2": "ref",
            "f3": "alt",
            "f4": "RawScore",
            "f5": "PHRED",
        }
        types = {"f0": hl.tstr, "f1": hl.tint32, "f4": hl.tfloat32, "f5": hl.tfloat32}
        ht = hl.import_table(
            cadd_tsv,
            delimiter="\t",
            types=types,
            no_header=True,
            force_bgz=True,
            comment="#",
            min_partitions=1000,
            missing="NA",
        )
        ht = ht.rename(column_names)
        chr = hl.if_else(ht.chr.startswith("chr"), ht.chr, hl.format("chr%s", ht.chr))
        ht = ht.annotate(
            locus=hl.locus(chr, ht.pos, reference_genome="GRCh38"),
            alleles=hl.array([ht.ref, ht.alt]),
        )
        ht = ht.select("locus", "alleles", "RawScore", "PHRED")
        ht = ht.key_by("locus", "alleles")

        return ht

    snvs = load_cadd_raw(
        "gs://gnomad-qin/v4_annotations/cadd.v1.6.whole_genome_SNVs.tsv.bgz"
    )
    indel3_0 = load_cadd_raw(
        "gs://gnomad-qin/v4_annotations/cadd.v1.6.gnomad.genomes.r3.0.indel.tsv.bgz"
    )
    indel3_1 = hl.read_table("gs://gnomad-wphu/CADD-indels-gnomad.3.1.ht")
    indel3_1_complex = hl.read_table(
        "gs://gnomad-wphu/CADD-1.6-gnomad-complex-variants.ht"
    )
    indel4 = load_cadd_raw(
        "gs://gnomad-qin/v4_annotations/gnomAD_v4_new_indels.tsv.bgz"
    )

    # Merge the CADD predictions run for v3 versions.
    indel3 = indel3_0.union(indel3_1, indel3_1_complex)

    # This will avoid duplicated indels in gnomAD v3 and v4.
    indel3 = indel3.anti_join(indel4)

    ht = snvs.union(indel3, indel4)
    return ht


def main(args):
    """Generate Hail Tables with in silico predictors."""
    hl.init(
        default_reference="GRCh38",
        log="/insilico_predictors.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    if args.cadd:
        logger.info("Creating CADD Hail Table for GRCh38...")
        ht = create_cadd_grch38_ht()
        ht.write(
            "gs://gnomad/v4.0/annotations/in_silico_predictors/gnomad.v4.0.cadd.grch38.ht",
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
