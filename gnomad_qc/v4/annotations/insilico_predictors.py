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

    snvs = hl.read_table(
        "gs://gcp-public-data--gnomad/resources/grch38/in_silico_predictors/CADD_v1.6_SNVs.ht"
    )
    indel3_0 = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.r3.0.indel.tsv.bgz"
    )
    indel3_1 = hl.read_table("gs://gnomad-insilico/cadd/CADD-indels-gnomad.3.1.ht")
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


def create_spliceai_grch38_ht(
    snvs_path: str, indels_path: str, new_indels_path: str, header_file_path: str
) -> hl.Table:
    """
    Create a Hail Table with SpliceAI scores for GRCh38.

    SpliceAI scores are from the following resources:
        - Precomputed SNVs: spliceai_scores.masked.snv.hg38.vcf.bgz,
          downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
        - Precomputed indels: spliceai_scores.masked.indel.hg38.vcf.bgz,
          downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
        - gnomAD v4 indels: gnomad_v4_new_indels.spliceai_masked.vcf.bgz,
          computed on v4 by Illumina.

    :param str snvs_path: Path to the precomputed SNVs.
    :param str indels_path: Path to the precomputed indels.
    :param str new_indels_path: Path to the gnomAD v4 indels.
    :return: Hail Table with SpliceAI scores for GRCh38.
    """

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

    logger.info("Importing vcf of SpliceAI scores into MT...")

    spliceai_snvs = import_spliceai_vcf(snvs_path)
    spliceai_indels = import_spliceai_vcf(indels_path)
    spliceai_new_indels = import_spliceai_vcf(new_indels_path)

    ht = spliceai_snvs.union(spliceai_indels).union(spliceai_new_indels)

    logger.info("Exploding SpliceAI scores...")
    # `explode` will eliminate rows with empty array, but the varaints with
    # multiple genes will extend the number of rows.
    ht = ht.explode(ht.info.SpliceAI)

    logger.info("Annotating SpliceAI scores...")
    # there will only be one gene in the array after exploding
    gene_symbol = ht.info.SpliceAI.split(delim="\\|")[1:2][0]
    delta_scores = ht.info.SpliceAI.split(delim="\\|")[2:6]
    positions = ht.info.SpliceAI.split(delim="\\|")[6:10]
    ht = ht.annotate(
        splice_ai=hl.struct(
            gene_symbol=gene_symbol,
            delta_scores=hl.map(lambda x: hl.float32(x), delta_scores),
            positions=hl.map(lambda x: hl.int32(x), positions),
        )
    )

    # Annotate info.max_DS with the max of DS_AG, DS_AL, DS_DG, DS_DL in info.
    # delta_score array is |DS_AG|DS_AL|DS_DG|DS_DL
    logger.info(
        "Getting the max SpliceAI score for each variant across consequences per"
        " gene..."
    )
    consequences = hl.literal(
        ["acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss"]
    )
    ht = ht.annotate(
        splice_ai=ht.splice_ai.annotate(ds_max=hl.max(ht.splice_ai.delta_scores))
    )
    ht = ht.annotate(
        splice_ai=ht.splice_ai.annotate(
            consequence_max=hl.if_else(
                ht.splice_ai.ds_max > 0,
                consequences[ht.splice_ai.delta_scores.index(ht.splice_ai.ds_max)],
                "no_consequence",
            ),
            position_max=hl.if_else(
                ht.splice_ai.ds_max > 0,
                ht.splice_ai.positions[
                    ht.splice_ai.delta_scores.index(ht.splice_ai.ds_max)
                ],
                hl.missing(hl.tint32),
            ),
        )
    )

    logger.info("Getting the max SpliceAI score for each variant across genes...")
    ht2 = ht.group_by(*ht.key).aggregate(
        staging=hl.agg.take(
            hl.struct(
                ds_max=ht.splice_ai.ds_max,
                position_max=ht.splice_ai.position_max,
                consequence_max=ht.splice_ai.consequence_max,
                gene=ht.splice_ai.gene_symbol,
            ),
            1,
            ordering=-ht.splice_ai.ds_max,
        )
    )

    logger.info("Annotating SpliceAI scores in right format...")
    # `aggregate` put everything in an array, we need to extract the struct from the array.
    ht2 = ht2.annotate(
        splice_ai=hl.struct(
            ds_max=hl.float32(ht2.staging.ds_max[0]),
            position_max=hl.int32(ht2.staging.position_max[0]),
            consequence_max=hl.str(ht2.staging.consequence_max[0]),
            gene=hl.str(ht2.staging.gene[0]),
        )
    )

    ht2 = ht2.select("splice_ai")
    return ht2


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
        snvs_path = (
            "gs://gnomad-insilico/spliceai/spliceai_scores.masked.snv.hg38.vcf.bgz"
        )
        indels_path = (
            "gs://gnomad-insilico/spliceai/spliceai_scores.masked.indel.hg38.vcf.bgz"
        )
        new_indels_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v4_new_indels.hg38.vcf.bgz"
        header_file_path = "gs://gnomad-insilico/spliceai/spliceai.vcf.header"
        ht = create_spliceai_grch38_ht(
            snvs_path, indels_path, new_indels_path, header_file_path
        )
        ht.write(
            get_insilico_predictors(predictor="spliceai").path,
            overwrite=args.overwrite,
        )
        logger.info("SpliceAI Hail Table for GRCh38 created.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT", action="store_true")
    parser.add_argument("--spliceai", help="Create SpliceAI HT", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
