"""Script to generate Hail Tables with in silico predictors."""
import argparse
import logging

import hail as hl
from gnomad.resources.resource_utils import NO_CHR_TO_CHR_CONTIG_RECODING
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import filter_vep_transcript_csqs
from hail.utils import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_insilico_predictors

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_sift_polyphen_from_vep(ht: hl.Table) -> hl.Table:
    """
    Get the max SIFT and PolyPhen scores from VEP 105 annotations.

    This retrieves the max of SIFT and PolyPhen scores for a variant's MANE Select
    transcript or, if MANE Select does not exist, canonical transcript.

    :param ht: VEP 105 annotated Hail Table.
    :return: Table annotated with max SIFT and PolyPhen scores.
    """
    mane = filter_vep_transcript_csqs(
        ht, synonymous=False, canonical=False, mane_select=True
    )
    canonical = filter_vep_transcript_csqs(ht, synonymous=False, canonical=True)
    ht = ht.annotate(
        sift_mane=mane[ht.key].vep.transcript_consequences.sift_score,
        polyphen_mane=mane[ht.key].vep.transcript_consequences.polyphen_score,
        sift_canonical=canonical[ht.key].vep.transcript_consequences.sift_score,
        polyphen_canonical=canonical[ht.key].vep.transcript_consequences.polyphen_score,
    )

    ht = ht.select(
        sift_max=hl.or_else(hl.max(ht.sift_mane), hl.max(ht.sift_canonical)),
        polyphen_max=hl.or_else(
            hl.max(ht.polyphen_mane), hl.max(ht.polyphen_canonical)
        ),
    )

    return ht


def create_cadd_grch38_ht() -> hl.Table:
    """
    Create a Hail Table with CADD scores for GRCh38.

    The combined CADD scores in the returned table are from the following sources:
        - all SNVs: `cadd.v1.6.whole_genome_SNVs.tsv.bgz` (81G) downloaded from
          `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/whole_genome_SNVs.tsv.gz`.
          It contains 8,812,917,339 SNVs.
        - gnomad 3.0 indels: `cadd.v1.6.gnomad.genomes.v3.0.indel.tsv.bgz` (1.1G)
          downloaded from `https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz`.
          It contains 100,546,109 indels from gnomaD v3.0.
        - gnomad 3.1 indels: `cadd.v1.6.gnomad.genomes.v3.1.indels.new.ht` was run
          on gnomAD v3.1 with CADD v1.6 in 2020. It contains 166,122,720 new indels from
          gnomAD v3.1 compared to v3.0.
        - gnomad 3.1 complex indels: `cadd.v1.6.gnomad.genomes.v3.1.indels.complex.ht`
          was run on gnomAD v3.1 with CADD v1.6 in 2020. It contains 2,307 complex
          variants that do not fit Hail's criteria for an indel and thus exist in
          a separate table than the gnomad 3.1 indels.
        - gnomAD v4 exomes indels: `cadd.v1.6.gnomad.exomes.v4.0.indels.new.tsv.bgz`
          (368M) was run on gnomAD v4 with CADD v1.6 in 2023. It contains 32,561,
          253 indels that are new in gnomAD v4.
        - gnomAD v4 genomes indels: `cadd.v1.6.gnomad.genomes.v4.0.indels.new.tsv.bgz`
          (13M) was run on gnomAD v4 with CADD v1.6 in 2023. It contains 904,906 indels
          that are new in gnomAD v4 genomes because of the addition of HGDP/TGP samples.

         .. note::
         ~1,9M indels were duplicated in gnomAD v3.0 and v4.0 or in gnomAD v3.1 and
         v4.0. However, CADD only generates a score per loci. We keep only the latest
         prediction, v4.0, for these loci.
         The output generated a CADD HT with 9,110,177,520 rows.
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
        ht = ht.checkpoint(new_temp_file("cadd", "ht"))
        return ht

    snvs = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.whole_genome_SNVs.tsv.bgz"
    )
    indel3_0 = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.v3.0.indel.tsv.bgz"
    )
    indel3_1 = hl.read_table(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.v3.1.indels.new.ht"
    )
    indel3_1_complex = hl.read_table(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.v3.1.indels.complex.ht"
    )
    indel4_e = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.exomes.v4.0.indels.new.tsv.bgz"
    )
    indel4_g = _load_cadd_raw(
        "gs://gnomad-insilico/cadd/cadd.v1.6.gnomad.genomes.v4.0.indels.new.tsv.bgz"
    )

    # Merge the CADD predictions run for v3 versions.
    indel3 = indel3_0.union(indel3_1, indel3_1_complex)

    # Merge the CADD predictions run for v4 versions.
    indel4 = indel4_e.union(indel4_g).distinct()

    # This will avoid duplicated indels in gnomAD v3 and v4.
    indel3 = indel3.anti_join(indel4)

    ht = snvs.union(indel3, indel4)

    ht = ht.select(cadd=hl.struct(phred=ht.PHRED, raw_score=ht.RawScore))
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
        - gnomAD v3 indels: gnomad_v3_indel.spliceai_masked.vcf.bgz,
          computed on v3.1 indels by Illumina in 2020.
        - gnomAD v4 indels: gnomad_v4_new_indels.spliceai_masked.vcf.bgz,
          computed on v4 indels that are new compared to v3 indels by Illumina in
          February 2023.
        - gnomAD v3 and v4 unscored indels:
          spliceai_scores.masked.gnomad_v3_v4_unscored_indels.hg38.vcf.bgz,
          another set of indels were not scored in v3 or v4 but computed by Illumina in
          September 2023.

    :return: Hail Table with SpliceAI scores for GRCh38.
    """
    snvs_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.snv.hg38.vcf.bgz"
    indels_path = (
        "gs://gnomad-insilico/spliceai/spliceai_scores.masked.indel.hg38.vcf.bgz"
    )
    v3_indels_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v3_indel.hg38.vcf.bgz"
    v4_indels_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v4_new_indels.hg38.vcf.bgz"
    v3_v4_unscored_path = "gs://gnomad-insilico/spliceai/spliceai_scores.masked.gnomad_v3_v4_unscored_indels.hg38.vcf.bgz"
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
    v3_indels = import_spliceai_vcf(v3_indels_path)
    v4_indels = import_spliceai_vcf(v4_indels_path)
    v3_v4_unscored_indels = import_spliceai_vcf(v3_v4_unscored_path)

    ht = spliceai_snvs.union(
        spliceai_indels, v3_indels, v4_indels, v3_v4_unscored_indels
    )

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
    ht = ht.select(spliceai_ds_max=hl.max(ht.values.ds_max))
    ht = ht.annotate_globals(spliceai_version="v1.3")
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
    gene body only with code from developers at Invitae:
    https://github.com/invitae/pangolin. All v4 genomes variants (except ~20M
    bug-affected and ~3M new variants from HGDP/TGP samples, noted below) were run on
    Pangolin v1.3.12, the others were run on Pangolin v1.4.4.

    :return: Hail Table with Pangolin score for splicing for GRCh38.
    """
    v4_genomes = (
        "gs://gnomad-insilico/pangolin/gnomad.v4.0.genomes.pangolin.vcf.bgz/*.gz"
    )
    # ~3 million new variants in gnomAD v4 genomes from the addition of HGDP/TGP
    # samples, they were run with the fixed code.
    v4_genomes_new = "gs://gnomad-insilico/pangolin/gnomad.v4.0.genomes.new_variants.pangolin.vcf.bgz/*.gz"
    v4_exomes = "gs://gnomad-insilico/pangolin/gnomad.v4.0.exomes.pangolin.vcf.bgz/*.gz"

    # There was a bug in original Pangolin code, that would generate incorrect
    # scores when a variant falls on multiple genes on the same strand. This bug
    # impacted v4 genomes but was fixed before running v4 exomes, the affected
    #  ~20M v4 genome variants were rerun.
    v4_bugfix = (
        "gs://gnomad-insilico/pangolin/gnomad.v4.0.genomes.bugfix.pangolin.vcf.bgz/*.gz"
    )

    header_file_path = "gs://gnomad-insilico/pangolin/pangolin.vcf.header"

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
            header_file=header_file_path,
        ).rows()
        ht = ht.checkpoint(hl.utils.new_temp_file("pangolin", "ht"))
        return ht

    ht_g = import_pangolin_vcf(v4_genomes)
    ht_g_new = import_pangolin_vcf(v4_genomes_new)
    ht_e = import_pangolin_vcf(v4_exomes)
    ht_bugfix = import_pangolin_vcf(v4_bugfix)
    ht_g_correct = ht_g.anti_join(ht_bugfix)
    ht_g = ht_g_correct.union(ht_bugfix)

    logger.info(
        f"Number of rows in genomes Pangolin HT: {ht_g.count()};\n Number of rows in"
        f" raw exomes Pangolin HT: {ht_e.count()};\n Number of rows in Pangolin HT for"
        f" new variants of genomes: {ht_g_new.count()}.\n"
    )

    ht = ht_g.union(ht_e, ht_g_new)

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
        pangolin_largest_ds=hl.if_else(
            hl.abs(hl.min(ht.values.largest_ds_gene))
            > hl.abs(hl.max(ht.values.largest_ds_gene)),
            hl.min(ht.values.largest_ds_gene),
            hl.max(ht.values.largest_ds_gene),
        )
    )
    # TODO: get the version for genomes run in May 2023.
    ht = ht.annotate_globals(pangolin_version=["v1.3.12", "v1.4.4"])
    logger.info(
        "\nNumber of variants indicating splice gain: %s;\n"
        "Number of variants indicating splice loss: %s; \n"
        "Number of variants with no splicing consequence: %s \n",
        hl.agg.count_where(ht.pangolin_largest_ds > 0),
        hl.agg.count_where(ht.pangolin_largest_ds < 0),
        hl.agg.count_where(ht.pangolin_largest_ds == 0),
    )
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
        "gs://gnomad-tmp-4day/revel-v1.3_in_Ensembl105.ht", overwrite=True
    )

    # Since the REVEL score for each variant is transcript-specific, we
    # prioritize the scores predicted on MANE Select and canonical transcripts,
    # and take the max if a variants falls on multiple MANE Select or canonical
    # transcripts. Normally, the score should be equal across MANE Select
    # and canonical.
    logger.info("Taking max REVEL scores for MANE Select transcripts...")
    mane_ht = ht.filter(hl.is_defined(ht.revel_mane), keep=True)
    max_revel_mane = mane_ht.group_by(*mane_ht.key).aggregate(
        revel_mane_max=hl.agg.max(mane_ht.revel_mane),
    )

    logger.info("Taking max REVEL scores for canonical transcripts...")
    canonical_ht = ht.filter(
        (~hl.is_defined(ht.revel_mane)) & (hl.is_defined(ht.revel_canonical)), keep=True
    )
    max_revel_canonical = canonical_ht.group_by(*canonical_ht.key).aggregate(
        revel_canonical_max=hl.agg.max(canonical_ht.revel_canonical),
    )
    logger.info(
        "Merge max REVEL scores for MANE Select transcripts and canonical transcripts"
        " to one HT..."
    )
    final_ht = max_revel_mane.union(max_revel_canonical)
    final_ht = final_ht.select(
        revel_max=hl.or_else(final_ht.revel_mane_max, final_ht.revel_canonical_max)
    )
    logger.info("Number of rows in final REVEL HT: %s", final_ht.count())
    final_ht = final_ht.annotate_globals(revel_version="v1.3")
    return final_ht


def create_phylop_grch38_ht() -> hl.Table:
    """
    Convert PhyloP scores to Hail Table.

    .. note::
    BigWig format of Phylop was download from here:
    https://cgl.gi.ucsc.edu/data/cactus/241-mammalian-2020v2-hub/Homo_sapiens/241-mammalian-2020v2.bigWig
    and converted it to bedGraph format with bigWigToBedGraph from the kent packages
    of UCSC (https://hgdownload.cse.ucsc.edu/admin/exe/) with the following command:
    `./bigWigToBedGraph ~/Downloads/241-mammalian-2020v2.bigWig ~/Downloads/241-mammalian-2020v2.bedGraph`
    The bedGraph file is bigzipped before importing to Hail.
    Different to other in silico predictors, the Phylop HT is keyed by locus only. Since
     the PhyloP scores have one value per position, we exploded the interval to store
     the HT by locus. In result, we have Phylop scores for 2,852,623,265 locus from
     2,648,607,958 intervals.

    :return: Hail Table with Phylop Scores for GRCh38
    """
    bg_path = "gs://gnomad-insilico/phylop/Human-GRCh38-Phylop-241-mammalian-2020v2.bedGraph.bgz"
    columns = ["chr", "start", "end", "phylop"]
    ht = hl.import_table(
        bg_path,
        min_partitions=1000,
        impute=True,
        no_header=True,
    ).rename({f"f{i}": c for i, c in enumerate(columns)})

    # We add 1 to both start and end because input bedGraph is 0-indexed interval and
    # the interval is end exclusive
    ht = ht.annotate(pos=hl.range(ht.start + 1, ht.end + 1))
    ht = ht.explode("pos")
    ht = ht.annotate(locus=hl.locus(ht.chr, ht.pos, reference_genome="GRCh38"))
    ht = ht.select("locus", "phylop")
    ht = ht.key_by("locus")
    ht = ht.annotate_globals(phylop_version="v2")

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

    if args.revel:
        logger.info("Creating REVEL Hail Table for GRCh38...")

        ht = create_revel_grch38_ht()
        ht.write(
            get_insilico_predictors(predictor="revel").path,
            overwrite=args.overwrite,
        )
        logger.info("REVEL Hail Table for GRCh38 created.")
    if args.phylop:
        logger.info("Creating PhyloP Hail Table for GRCh38...")
        ht = create_phylop_grch38_ht()
        ht.write(
            get_insilico_predictors(predictor="phylop").path,
            overwrite=args.overwrite,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument("--cadd", help="Create CADD HT", action="store_true")
    parser.add_argument("--spliceai", help="Create SpliceAI HT", action="store_true")
    parser.add_argument("--pangolin", help="Create Pangolin HT", action="store_true")
    parser.add_argument("--revel", help="Create REVEL HT.", action="store_true")
    parser.add_argument("--phylop", help="Create PhyloP HT.", action="store_true")
    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
