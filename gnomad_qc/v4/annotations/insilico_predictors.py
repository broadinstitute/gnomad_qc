"""Script to generate Hail Tables with in silico predictors."""

import argparse
import logging

import hail as hl
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


def create_spliceai_grch38_ht(
    snvs_path: str, indels_path: str, new_indels_path: str, header_file_path: str
) -> hl.Table:
    """
    Create a Hail Table with SpliceAI scores for GRCh38.

    ..Note::
    The SpliceAI scores are from the following resources:
    all SNVs precomputed: spliceai_scores.masked.snv.hg38.vcf.bgz, downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
    all indels precomputed: spliceai_scores.masked.indel.hg38.vcf.bgz, downloaded from https://basespace.illumina.com/s/5u6ThOblecrh
    new indels from gnomAD v4: gnomad_v4_new_indels.spliceai_masked.vcf.bgz, requested Illumina to compute for us.
    We may have lost the SpliceAI score for gnomAD v3 indels, which makes it hard to process the raw scores to get a maximum and retain the position information.
    Before we didn't `explode` the SpliceAI scores, so we only had scores for the first gene for each variant.

    :return: Hail Table with SpliceAI scores for GRCh38.
    """
    logger.info("Importing vcf of SpliceAI scores into MT...")
    recode = {f"{i}": f"chr{i}" for i in (list(range(1, 23)) + ["X", "Y"])}
    spliceai_snvs = hl.import_vcf(
        snvs_path,
        force_bgz=True,
        min_partitions=1000,
        reference_genome="GRCh38",
        contig_recoding=recode,
        skip_invalid_loci=True,
    )

    spliceai_indels = hl.import_vcf(
        indels_path,
        force_bgz=True,
        reference_genome="GRCh38",
        contig_recoding=recode,
        skip_invalid_loci=True,
        min_partitions=1000,
    )

    spliceai_new_indels = hl.import_vcf(
        new_indels_path,
        header_file=header_file_path,
        force_bgz=True,
        reference_genome="GRCh38",
        contig_recoding=recode,
        skip_invalid_loci=True,
        min_partitions=1000,
    )

    mt = spliceai_snvs.union_rows(spliceai_indels).union_rows(spliceai_new_indels)
    ht = mt.rows().repartition(3000)
    logger.info("Number of rows in original SpliceAI Hail Table: %s", ht.count())

    logger.info("Exploding SpliceAI scores...")
    # `explode` will eliminate rows with empty array, so we expect to have fewer rows after exploding.
    ht = ht.explode(ht.info.SpliceAI)
    logger.info("Number of rows in exploded SpliceAI Hail Table: %s", ht.count())

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
    ht = ht.checkpoint("gs://gnomad-tmp-4day/spliceai.ht", overwrite=True)

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
    logger.info(
        "Number of rows in SpliceAI Hail Table after aggregation: %s", ht2.count()
    )

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
            "gs://gnomad-qin/v4_annotations/spliceai_scores.masked.snv.hg38.vcf.bgz"
        )
        indels_path = (
            "gs://gnomad-qin/v4_annotations/spliceai_scores.masked.indel.hg38.vcf.bgz"
        )
        new_indels_path = "gs://gnomad-qin/v4_annotations/gnomad_v4_new_indels.spliceai_masked.vcf.bgz"
        header_file_path = "gs://gnomad-qin/v4_annotations/spliceai.vcf.header"
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
