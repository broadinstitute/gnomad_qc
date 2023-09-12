"""Script to generate Hail Tables with in silico predictors."""

import argparse
import logging

import hail as hl
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_insilico_predictors, get_insilico_raw

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def create_cadd_grch38_ht() -> hl.Table:
    """
    Create a Hail Table with CADD scores for GRCh38.

    The combined CADD scores in the returned table are from the following sources:
        - all SNVs: `cadd.v1.6.whole_genome_SNVs.tsv.bgz` (81G) downloaded from
        `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz`.
        It contains 8,812,917,339 SNVs.
        - gnomad 3.0 indels: `cadd.v1.6.indels.tsv.bgz` (1.1G) downloaded from `
        https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz`.
        It contains 100,546,109 indels from gnomaD v3.0.
        - gnomad 3.1 indels: `CADD-indels-gnomad.3.1.ht` was run on gnomAD v3.1
        with CADD v1.6 in 2020. It contains 166,122,720 indels from gnomAD v3.1.
        - gnomad 3.1 complex indels: `CADD-1.6-gnomad-complex-variants.ht` was
        run on gnomAD v3.1 with CADD v1.6 in 2020. It contains 2,307 complex
        variants that do not fit Hail's criteria for an indel and thus exist in
        a separate table than the gnomad 3.1 indels.
        - gnomAD v4 indels: `cadd.v1.6.gnomAD_v4_new_indels.tsv.bgz` (368M) was
        run on gnomAD v4 with CADD v1.6 in 2023. It contains 32,561,253 indels
        that are new in gnomAD v4.

         .. note::
         1,972,208 indels were duplicated in gnomAD v3.0 and v4.0 or in gnomAD
         v3.1 and v4.0. However, CADD only generates a score per loci.
         We keep only the latest prediction, v4.0, for these loci.
    :return: Hail Table with CADD scores for GRCh38.
    """

    def _load_cadd_raw(cadd_tsv) -> hl.Table:
        """Functions to load CADD raw data in TSV format to Hail Table."""
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
            types=types,
            no_header=True,
            force_bgz=True,
            comment="#",
            min_partitions=1000,
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

    snvs = _load_cadd_raw(
        "gs://gnomad-qin/v4_annotations/cadd.v1.6.whole_genome_SNVs.tsv.bgz"
    )
    indel3_0 = _load_cadd_raw(
        "gs://gnomad-qin/v4_annotations/cadd.v1.6.gnomad.genomes.r3.0.indel.tsv.bgz"
    )
    indel3_1 = hl.read_table("gs://gnomad-wphu/CADD-indels-gnomad.3.1.ht")
    indel3_1_complex = hl.read_table(
        "gs://gnomad-wphu/CADD-1.6-gnomad-complex-variants.ht"
    )
    indel4 = _load_cadd_raw(
        "gs://gnomad-qin/v4_annotations/cadd.v1.6.gnomAD_v4_new_indels.tsv.bgz"
    )

    # Merge the CADD predictions run for v3 versions.
    indel3 = indel3_0.union(indel3_1, indel3_1_complex)

    # This will avoid duplicated indels in gnomAD v3 and v4.
    indel3 = indel3.anti_join(indel4)

    ht = snvs.union(indel3, indel4)
    ht = ht.annotate_globals(cadd_version="v1.6")
    return ht


def create_revel_grch38_ht(revel_csv: str, ensembl_ids: str) -> hl.Table:
    """
    Create a Hail Table with REVEL scores for GRCh38.

    .. note::
    Starting with gnomAD v4, we use REVEL scores for only MANE Select and
    canonical transcripts. Even when a variant falls on multiple MANE/canonical
    transcripts of different genes, the scores are equal.
    REVEL's Ensembl ID is not from Ensembl 105, so we filter to transcripts
    that are in Ensembl 105.
    This deprecates the `has_duplicate` field present in gnomAD v3.

    :param revel_csv: Path to REVEL CSV file, downloaded from
       https://www.google.com/url?q=https%3A%2F%2Frothsj06.dmz.hpc.mssm.edu%2Frevel-v1.3_all_chromosomes.zip&sa=D&sntz=1&usg=AOvVaw2DS2TWUYl__0vqijzzxp5M
       size ~648M, ~82,100,677 variants
    :param ensembl_ids: Path to Ensembl 105 ID file,
       downloaded from Ensembl 105 archive. It contains the following columns:
       Transcript stable ID, Ensembl Canonical, MANE Select
    :return: Hail Table with REVEL scores for GRCh38.
    """
    ht = hl.import_table(
        revel_csv,
        delimiter=",",
        min_partitions=1000,
        types={"grch38_pos": hl.tstr, "REVEL": hl.tfloat64},
    )

    logger.info("Annotating REVEL table...")
    ht = ht.drop("hg19_pos", "aaref", "aaalt")
    # drop variants that have no position in GRCh38 when lifted over from GRCh37
    ht = ht.filter(ht.grch38_pos.contains("."), keep=False)
    ht = ht.transmute(chr="chr" + ht.chr)
    ht = ht.annotate(
        locus=hl.locus(ht.chr, hl.int(ht.grch38_pos), reference_genome="GRCh38"),
        alleles=hl.array([ht.ref, ht.alt]),
        Transcript_stable_ID=ht.Ensembl_transcriptid.strip().split(";"),
    )
    ht = ht.select("locus", "alleles", "REVEL", "Transcript_stable_ID")
    ht = ht.explode("Transcript_stable_ID")
    ht = ht.key_by("Transcript_stable_ID")
    print(ht.count())

    logger.info("Import and process the ensembl ID file...")
    ensembl_ids = hl.import_table(
        ensembl_ids, min_partitions=200, impute=True, key="Transcript_stable_ID"
    )
    ensembl_ids = ensembl_ids.select("Ensembl_Canonical", "MANE_Select")
    logger.info("Annotating revel HT with canonical and MANE Select transcripts...")
    ht = ht.annotate(**ensembl_ids[ht.Transcript_stable_ID])
    print(ht.count())

    logger.info(
        "Annotating REVEL scores for MANE Select transcripts and canonical"
        " transcripts..."
    )
    ht = ht.key_by("locus", "alleles")
    ht = ht.annotate(
        revel_mane=hl.or_missing(ht.MANE_Select != "", ht.REVEL),
        revel_canonical=hl.or_missing(ht.Ensembl_Canonical == "1", ht.REVEL),
        revel_noncanonical=hl.or_missing(ht.Ensembl_Canonical != "1", ht.REVEL),
    )

    logger.info("Taking max REVEL scores for MANE Select transcripts...")
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/revel-v1.3_in_Ensembl105.ht", _read_if_exists=True
    )
    mane_ht = ht.filter(hl.is_defined(ht.revel_mane), keep=True)
    max_revel_mane = mane_ht.group_by(*mane_ht.key).aggregate(
        max_revel=hl.agg.max(mane_ht.revel_mane),
    )

    logger.info("Taking max REVEL scores for canonical transcripts...")
    canonical_ht = ht.filter(
        (~hl.is_defined(ht.revel_mane)) & (hl.is_defined(ht.revel_canonical)), keep=True
    )
    max_revel_canonical = canonical_ht.group_by(*canonical_ht.key).aggregate(
        max_revel=hl.agg.max(canonical_ht.revel_canonical),
    )
    logger.info(
        "Merge max REVEL scores for MANE Select transcripts and canonical transcripts"
        " to one HT..."
    )
    final_ht = max_revel_mane.union(max_revel_canonical)
    return final_ht


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
            get_insilico_predictors(predictor="cadd").path,
            overwrite=args.overwrite,
        )
        logger.info("CADD Hail Table for GRCh38 created.")

    if args.revel:
        logger.info("Creating REVEL Hail Table for GRCh38...")
        revel_csv = get_insilico_raw(predictor="revel-v1.3", postfix="csv.bgz")
        ensembl_ids = get_insilico_raw(predictor="ensembl105id", postfix="tsv.bgz")

        ht = create_revel_grch38_ht(revel_csv, ensembl_ids)
        ht.write(
            get_insilico_predictors(predictor="revel").path,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT", action="store_true")
    parser.add_argument("--revel", help="Create REVEL HT", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)