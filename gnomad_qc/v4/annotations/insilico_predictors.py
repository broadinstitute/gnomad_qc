"""Script to generate Hail Tables with in silico predictors."""

import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import NO_CHR_TO_CHR_CONTIG_RECODING
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_insilico_predictors

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
        - gnomad 3.0 indels: `cadd.v1.6.indels.tsv.bgz` (1.1G) downloaded from
        `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz`.
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
        "gs://gnomad-insilico/cadd/cadd.v1.6.whole_genome_SNVs.tsv.bgz"
    )
    indel3_0 = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.r3.0.indel.tsv.bgz"
    )
    indel3_1 = hl.read_table("gs://gnomad-insilico/cadd/CADD-indels-gnomad.3.1.ht")
    indel3_1_complex = hl.read_table(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.v3.1.indels.complex.ht"
    )
    indel4 = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.v4.0.indels.new.tsv.bgz"
    )

    # Merge the CADD predictions run for v3 versions.
    indel3 = indel3_0.union(indel3_1, indel3_1_complex)

    # This will avoid duplicated indels in gnomAD v3 and v4.
    indel3 = indel3.anti_join(indel4)

    ht = snvs.union(indel3, indel4)
    ht = ht.annotate_globals(cadd_version="v1.6")
    return ht


def create_spliceai_grch38_ht() -> hl.Table:
    """
    Create a Hail Table with SpliceAI scores for GRCh38.

    SpliceAI scores are from the following resources:
        - Precomputed SNVs: spliceai_scores.masked.snv.hg38.vcf.bgz,
          downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
        - Precomputed indels: spliceai_scores.masked.indel.hg38.vcf.bgz,
          downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
        - gnomAD v4 indels: gnomad_v4_new_indels.spliceai_masked.vcf.bgz,
          computed on v4 by Illumina.

    :return: Hail Table with SpliceAI scores for GRCh38.
    """
    snvs_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.snv.hg38.vcf.bgz"
    indels_path = (
        "gs://gnomad-insilico/spliceai/spliceai_scores.masked.indel.hg38.vcf.bgz"
    )
    new_indels_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v4_new_indels.hg38.vcf.bgz"
    header_file_path = "gs://gnomad-insilico/spliceai/spliceai.vcf.header"

    def import_spliceai_vcf(path: str) -> hl.Table:
        """
        Import SpliceAI VCF into Hail Table.

        :param str path: Path to SpliceAI VCF.
        :return: Hail MatrixTable with SpliceAI scores.
        """
        ht = hl.import_vcf(
            path,
            force_bgz=True,
            header_file=header_file_path,
            reference_genome="GRCh38",
            contig_recoding=NO_CHR_TO_CHR_CONTIG_RECODING,
            skip_invalid_loci=True,
            min_partitions=1000,
        ).rows()
        return ht

    logger.info("Importing vcf of SpliceAI scores into HT...")
    spliceai_snvs = import_spliceai_vcf(snvs_path)
    spliceai_indels = import_spliceai_vcf(indels_path)
    spliceai_new_indels = import_spliceai_vcf(new_indels_path)

    ht = spliceai_snvs.union(spliceai_indels, spliceai_new_indels)

    # `explode` because some variants fall on multiple genes and have a score per gene.
    # All rows without a SpliceAI score, an empty array, are removed through explode.
    logger.info("Exploding SpliceAI scores...")
    ht = ht.explode(ht.info.SpliceAI)

    # delta_score array for 4 splicing consequences: DS_AG|DS_AL|DS_DG|DS_DL
    logger.info("Annotating SpliceAI scores...")
    delta_scores = ht.info.SpliceAI.split(delim="\\|")[2:6]
    ht = ht.annotate(delta_scores=hl.map(lambda x: hl.float32(x), delta_scores))

    logger.info(
        "Getting the max SpliceAI score across consequences for each variant per"
        " gene..."
    )
    ht = ht.select(ds_max=hl.max(ht.delta_scores))

    logger.info("Getting the max SpliceAI score for each variant across genes...")
    ht = ht.collect_by_key()
    ht = ht.select(splice_ai=hl.struct(ds_max=hl.max(ht.values.ds_max)))
    ht = ht.annotate_globals(spliceai_version="v1.3")
    return ht


def create_revel_grch38_ht() -> hl.Table:
    """
    Create a Hail Table with REVEL scores for GRCh38.

    .. note::
    Starting with gnomAD v4, we use REVEL scores for only MANE Select and
    canonical transcripts. Even when a variant falls on multiple MANE/canonical
    transcripts of different genes, the scores are equal.
    REVEL scores were downloaded from:
       https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
       size ~648M, ~82,100,677 variants
    REVEL's Ensembl ID is not from Ensembl 105, so we filter to transcripts
    that are in Ensembl 105. The Ensembl 105 ID file was downloaded from Ensembl
    105 archive. It contains the following columns:
       Transcript stable ID, Ensembl Canonical, MANE Select
    This deprecates the `has_duplicate` field present in gnomAD v3.

    :return: Hail Table with REVEL scores for GRCh38.
    """
    revel_csv = "gs://gnomad-insilico/revel/revel-v1.3_all_chromosomes_with_transcript_ids.csv.bgz"
    ensembl_path = "gs://gnomad-insilico/ensembl/ensembl105id.grch38.tsv.bgz"

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
    ht = ht.select(
        locus=hl.locus(ht.chr, hl.int(ht.grch38_pos), reference_genome="GRCh38"),
        alleles=hl.array([ht.ref, ht.alt]),
        REVEL=ht.REVEL,
        Transcript_stable_ID=ht.Ensembl_transcriptid.strip().split(";"),
    )
    ht = ht.explode("Transcript_stable_ID")
    ht = ht.key_by("Transcript_stable_ID")
    logger.info("Number of rows in REVEL table: %s", ht.count())

    logger.info("Import and process the ensembl ID file...")
    ensembl_ids_ht = hl.import_table(
        ensembl_path, min_partitions=200, impute=True, key="Transcript_stable_ID"
    )
    ensembl_ids_ht = ensembl_ids_ht.select("Ensembl_Canonical", "MANE_Select")

    logger.info("Annotating REVEL HT with canonical and MANE Select transcripts...")
    ht = ht.annotate(**ensembl_ids_ht[ht.key])

    logger.info(
        "Annotating REVEL scores for MANE Select transcripts and canonical"
        " transcripts..."
    )
    ht = ht.key_by("locus", "alleles")
    ht = ht.annotate(
        revel_mane=hl.or_missing(ht.MANE_Select != "", ht.REVEL),
        revel_canonical=hl.or_missing(ht.Ensembl_Canonical == "1", ht.REVEL),
    )
    ht = ht.checkpoint(
        "gs://gnomad-tmp-4day/revel-v1.3_in_Ensembl105.ht", overwrite=not args.overwrite
    )

    # Since the REVEL score for each variant is transcript-specific, we
    # prioritize the scores predicted on MANE Select and canonical transcripts,
    # and take the max if a variants falls on multiple MANE Select or canonical
    # transcripts. Normally, the score should be equal across MANE Select
    # and canonical.
    logger.info("Taking max REVEL scores for MANE Select transcripts...")
    mane_ht = ht.filter(hl.is_defined(ht.revel_mane), keep=True)
    max_revel_mane = mane_ht.group_by(*mane_ht.key).aggregate(
        revel_max=hl.agg.max(mane_ht.revel_mane),
    )

    logger.info("Taking max REVEL scores for canonical transcripts...")
    canonical_ht = ht.filter(
        (~hl.is_defined(ht.revel_mane)) & (hl.is_defined(ht.revel_canonical)), keep=True
    )
    max_revel_canonical = canonical_ht.group_by(*canonical_ht.key).aggregate(
        revel_max=hl.agg.max(canonical_ht.revel_canonical),
    )
    logger.info(
        "Merge max REVEL scores for MANE Select transcripts and canonical transcripts"
        " to one HT..."
    )
    final_ht = max_revel_mane.union(max_revel_canonical)
    logger.info("Number of rows in final REVEL HT: %s", final_ht.count())
    final_ht = final_ht.annotate_globals(revel_version="v1.3")
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

    if args.spliceai:
        logger.info("Creating SpliceAI Hail Table for GRCh38...")
        ht = create_spliceai_grch38_ht()
        ht.write(
            get_insilico_predictors(predictor="spliceai").path,
            overwrite=args.overwrite,
        )
        logger.info("SpliceAI Hail Table for GRCh38 created.")

    if args.revel:
        logger.info("Creating REVEL Hail Table for GRCh38...")

        ht = create_revel_grch38_ht()
        ht.write(
            get_insilico_predictors(predictor="revel").path,
            overwrite=args.overwrite,
        )
        logger.info("REVEL Hail Table for GRCh38 created.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT", action="store_true")
    parser.add_argument("--spliceai", help="Create SpliceAI HT", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT.", action="store_true")
    parser.add_argument("--revel", help="Create REVEL HT.", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
