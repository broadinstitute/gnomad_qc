"""Script to generate Hail Tables with in silico predictors."""

import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import NO_CHR_TO_CHR_CONTIG_RECODING
from gnomad.utils.slack import slack_notifications

from gnomad_qc.resource_utils import check_resource_existence
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
        - all SNVs: `cadd.v1.6.whole_genome_SNVs.tsv.bgz` (81G) downloaded from `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz`. It contains 8,812,917,339 SNVs.
        - gnomad 3.0 indels: `cadd.v1.6.indels.tsv.bgz` (1.1G) downloaded from `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz`. It contains 100,546,109 indels from gnomaD v3.0.
        - gnomad 3.1 indels: `CADD-indels-gnomad.3.1.ht` was run on gnomAD v3.1 with CADD v1.6 in 2020. It contains 166,122,720 indels from gnomAD v3.1.
        - gnomad 3.1 complex indels: `CADD-1.6-gnomad-complex-variants.ht` was run on gnomAD v3.1 with CADD v1.6 in 2020. It contains 2,307 complex variants that do not fit Hail's criteria for an indel and thus exist in a separate table than the gnomad 3.1 indels.
        - gnomAD v4 indels: `cadd.v1.6.gnomAD_v4_new_indels.tsv.bgz` (368M) was run on gnomAD v4 with CADD v1.6 in 2023. It contains 32,561,253 indels that are new in gnomAD v4.

         .. note::
         1,972,208 indels were duplicated in gnomAD v3.0 and v4.0 or in gnomAD v3.1 and v4.0. However, CADD only generates a score per loci. We keep only the latest prediction, v4.0, for these loci.
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

    return ht


def create_pangolin_grch38_ht() -> hl.Table:
    """
    Create a Hail Table with Pangolin score for splicing for GRCh38.

    .. note::
    The score was based on the splicing prediction tool Pangolin:
    Zeng, T., Li, Y.I. Predicting RNA splicing from DNA sequence using Pangolin.
     Genome Biol 23, 103 (2022). https://doi.org/10.1186/s13059-022-02664-4

    There's no precomputed score for all possible variants, the scores were
    generated for gnomAD v4 genomes (=v3 genomes) and v4 exomes variants in
    gene body only.

    :return: Hail Table with Pangolin score for splicing for GRCh38.
    """
    v4_genomes = (
        "gs://gnomad-insilico/pangolin/gnomad.v4.0.genomes.pangolin.vcf.bgz/*.bgz"
    )
    v4_exomes = "gs://gnomad-insilico/pangolin/gnomad.v4.0.exomes.pangolin.vcf.bgz/*.gz"

    def import_pangolin_vcf(vcf_path: str) -> hl.Table:
        """
        Import Pangolin VCF to Hail Table.

        :param vcf_path: Path to Pangolin VCF.
        :return: Hail Table with Pangolin scores.
        """
        ht = hl.import_vcf(
            vcf_path,
            min_partitions=1000,
            force_bgz=True,
            reference_genome="GRCh38",
            skip_invalid_loci=True,
        ).rows()
        return ht

    ht_g = import_pangolin_vcf(v4_genomes)
    ht_e = import_pangolin_vcf(v4_exomes)

    logger.info(
        "Number of rows in genomes Pangolin HT: %s; "
        "Number of rows in raw exomes Pangolin HT: %s. ",
        ht_g.count(),
        ht_e.count(),
    )

    ht = ht_g.union(ht_e)

    logger.info("Exploding Pangolin scores...")
    # `explode` will eliminate rows with empty array
    # The Pangolin annotation is imported as an array of strings containing
    # one element with the following format:
    # gene1|pos_splice_gain:largest_increase|pos_splice_loss:largest_decrease|Warnings:||gene2...
    # for example:
    # Pangolin=ENSG00000121005.9|-86:0.25|38:-0.49|Warnings:||ENSG00000254238.1|-40:0.01|30:-0.17|Warnings:
    ht = ht.annotate(pangolin=ht.info.Pangolin[0].split(delim="\|\|"))
    ht = ht.explode(ht.pangolin)
    logger.info("Number of rows in exploded Pangolin Hail Table: %s", ht.count())

    logger.info("Annotating Pangolin scores...")
    # The Pangolin score is the delta score of splice gain and splice loss,
    # which is the second and fourth element in the array after splitting.
    ht = ht.transmute(
        pangolin=hl.struct(
            delta_scores=(
                hl.empty_array(hl.tfloat64).append(
                    hl.float(ht.pangolin.split(delim=":|\\|")[2])
                )
            ).append(hl.float(ht.pangolin.split(delim=":|\\|")[4])),
        )
    )

    # Using > instead of >= in case where delta score of splice gain and splice
    # loss are the same but different sign, we will keep the positive score of
    # splice gain
    ht = ht.annotate(
        largest_ds_gene=hl.if_else(
            hl.abs(hl.min(ht.pangolin.delta_scores))
            > hl.abs(hl.max(ht.pangolin.delta_scores)),
            hl.min(ht.pangolin.delta_scores),
            hl.max(ht.pangolin.delta_scores),
        )
    )
    ht = ht.select(ht.largest_ds_gene)
    ht = ht.collect_by_key()
    logger.info(
        "Number of rows in Pangolin Hail Table after collect_by_key: %s", ht.count()
    )
    ht = ht.select(
        pangolin=hl.struct(
            largest_ds=hl.if_else(
                hl.abs(hl.min(ht.values.largest_ds_gene))
                > hl.abs(hl.max(ht.values.largest_ds_gene)),
                hl.min(ht.values.largest_ds_gene),
                hl.max(ht.values.largest_ds_gene),
            )
        )
    )
    logger.info(
        "\nNumber of variants indicating splice gain: %s;\n"
        "Number of variants indicating splice loss: %s; \n"
        "Number of variants with no splicing consequence: %s \n",
        hl.agg.count_where(ht.pangolin.largest_ds > 0),
        hl.agg.count_where(ht.pangolin.largest_ds < 0),
        hl.agg.count_where(ht.pangolin.largest_ds == 0),
    )
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

    if args.pangolin:
        logger.info("Creating Pangolin Hail Table for GRCh38...")

        ht = create_pangolin_grch38_ht()
        ht.write(
            get_insilico_predictors(predictor="pangolin").path,
            overwrite=args.overwrite,
        )
        logger.info("Pangolin Hail Table for GRCh38 created.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT", action="store_true")
    parser.add_argument("--spliceai", help="Create SpliceAI HT", action="store_true")
    parser.add_argument("--pangolin", help="Create Pangolin HT", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
