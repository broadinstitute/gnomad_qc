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


def create_pangolin_grch38_ht(vcf_path: str) -> hl.Table:
    """
    Create a Hail Table with Pangolin score for splicing for GRCh38.

    .. note::
    The score was based on the splicing prediction tool Pangolin:
    Zeng, T., Li, Y.I. Predicting RNA splicing from DNA sequence using Pangolin.
     Genome Biol 23, 103 (2022). https://doi.org/10.1186/s13059-022-02664-4

    There's no precomputed for all variants, the scores were generated for
    gnomAD v3 variants in gene body only.

    :param vcf_path: path to the VCF files with Pangolin scores for splicing.
    :return: Hail Table with Pangolin score for splicing for GRCh38.
    """
    vcf_path = (
        "gs://gnomad-insilico/pangolin/gnomad.v4.0.genomes.pangolin.vcf.bgz/*.bgz"
    )
    ht = hl.import_vcf(vcf_path, min_partitions=1000, reference_genome="GRCh38").rows()
    logger.info("Number of rows in original Pangolin Hail Table: %s", ht.count())

    logger.info("Exploding Pangolin scores...")
    # `explode` will eliminate rows with empty array
    ht = ht.annotate(pango=ht.info.Pangolin[0].split(delim="\|\|"))
    ht = ht.explode(ht.pango)
    logger.info("Number of rows in exploded Pangolin Hail Table: %s", ht.count())

    logger.info("Annotating Pangolin scores...")
    ht = ht.annotate(
        pangolin=hl.struct(
            gene=ht.pango.split(delim=":|\\|")[0],
            positions=(
                hl.empty_array(hl.tint32).append(
                    hl.int(ht.pango.split(delim=":|\\|")[1])
                )
            ).append(hl.int(ht.pango.split(delim=":|\\|")[3])),
            delta_scores=(
                hl.empty_array(hl.tfloat64).append(
                    hl.float(ht.pango.split(delim=":|\\|")[2])
                )
            ).append(hl.abs(hl.float(ht.pango.split(delim=":|\\|")[4]))),
        )
    )

    logger.info(
        "Getting the max Pangolin score across consequences for each variant  per"
        " gene..."
    )
    consequences = hl.literal(["splice_gain", "splice_loss"])
    ht = ht.annotate(
        pangolin=ht.pangolin.annotate(ds_max=hl.max(ht.pangolin.delta_scores))
    )
    ht = ht.annotate(
        pangolin=ht.pangolin.annotate(
            consequence_max=hl.if_else(
                ht.pangolin.ds_max > 0,
                consequences[ht.pangolin.delta_scores.index(ht.pangolin.ds_max)],
                "no_consequence",
            ),
            position_max=hl.if_else(
                ht.pangolin.ds_max > 0,
                ht.pangolin.positions[
                    ht.pangolin.delta_scores.index(ht.pangolin.ds_max)
                ],
                hl.missing(hl.tint32),
            ),
        )
    )
    ht = ht.checkpoint("gs://gnomad-tmp-4day/pangolin.ht", overwrite=True)

    logger.info("Getting the max Pangolin score for each variant across genes...")
    ht2 = ht.group_by(*ht.key).aggregate(
        staging=hl.agg.take(
            hl.struct(
                ds_max=ht.pangolin.ds_max,
                position_max=ht.pangolin.position_max,
                consequence_max=ht.pangolin.consequence_max,
                gene=ht.pangolin.gene,
            ),
            1,
            ordering=-ht.pangolin.ds_max,
        )
    )
    logger.info(
        "Number of rows in Pangolin Hail Table after aggregation: %s", ht2.count()
    )

    # `aggregate` put everything in an array, we need to extract the struct from the array.
    ht2 = ht2.annotate(
        pangolin=hl.struct(
            ds_max=hl.float32(ht2.staging.ds_max[0]),
            position_max=hl.int32(ht2.staging.position_max[0]),
            consequence_max=hl.str(ht2.staging.consequence_max[0]),
            gene=hl.str(ht2.staging.gene[0]),
        )
    )

    ht2 = ht2.select("pangolin")
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
    parser.add_argument("--pangolin", help="Create Pangolin HT", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
