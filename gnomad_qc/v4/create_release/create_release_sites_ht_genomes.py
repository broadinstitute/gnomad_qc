"""Script to create release sites HT for v4 genomes."""
import hail as hl


def create_revel_grch38_ht(revel_csv, ensembl_ids):
    """
    Create a Hail Table with REVEL scores for GRCh38.

    :param revel_csv: Path to REVEL CSV file.
    :param ensembl_ids: Path to Ensembl 105 ID file, downloaded from Ensembl 105 archive.
    """
    # import and process the revel csv
    # TODO: update path
    ht = hl.import_table(
        revel_csv,
        delimiter=",",
        min_partitions=200,
        types={"grch38_pos": hl.tstr, "REVEL": hl.tfloat64},
    )
    ht = ht.drop("hg19_pos", "aaref", "aaalt")
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

    # import the ensembl 105 id file
    # TODO: update path
    ensembl_ids = hl.import_table(ensembl_ids, min_partitions=200)
    ensembl_ids = ensembl_ids.select("Transcript_stable_ID", "Ensembl_Canonical")
    ensembl_ids = ensembl_ids.key_by("Transcript_stable_ID")

    # join the two tables
    ht = ht.join(ensembl_ids, how="inner")

    # get the canonical and non-canonical revel scores
    ht = ht.key_by("locus", "alleles")
    ht = ht.annotate(
        revel_canonical=hl.if_else(
            ht.Ensembl_Canonical == "1", ht.REVEL, hl.missing(hl.tfloat)
        ),
        revel_noncanonical=hl.if_else(
            ht.Ensembl_Canonical != "1", ht.REVEL, hl.missing(hl.tfloat)
        ),
    )

    # get the max revel score for each variant
    max_revel_canonical = ht.group_by(*ht.key).aggregate(
        max_revel_canonical=hl.agg.max(ht.revel_canonical)
    )
    max_revel_noncanonical = ht.group_by(*ht.key).aggregate(
        max_revel_noncanonical=hl.agg.max(ht.revel_noncanonical)
    )
    revel = max_revel_canonical.join(max_revel_noncanonical, how="inner")
    revel.checkpoint(
        "gs://gnomad-qin/v4_annotations/revel-v1.3_processed_max_canonical.ht",
        overwrite=True,
    )
    # TODO: update path


revel_csv = "gs://gnomad-qin/v4_annotations/revel-v1.3_all_chromosomes_with_transcript_ids.csv.bgz"
ensembl_ids = (
    "gs://gnomad-qin/ensembl/ensembl105_chr1-22XY_MANE_canonical_ucscid.tsv.bgz"
)
create_revel_grch38_ht(revel_csv, ensembl_ids)


def create_primateai_grch38_ht(primateai_tsv, ucsc_ids, ensembl_ids):
    """
    Create a Hail Table with PrimateAI scores for GRCh38.

    :param primateai_tsv: Path to PrimateAI TSV file. The raw had the locus lifted over from GRCh37 to GRCh38, but the uscs ids are still in GRCh37.
    :param ucsc_ids: Path to UCSC ID file, downloaded from UCSC Table Browser, in version GRCh37, known2Ensembl.
    :param ensembl_ids: Path to Ensembl 105 ID file, downloaded from Ensembl 105 archive.
    """
    # import and process the primateAI tsv
    ht = hl.import_table(
        primateai_tsv,
        comment="#",
        skip_blank_lines=True,
        types={
            "pos": hl.tint32,
            "primateDL_score": hl.tfloat32,
            "ExAC_coverage": hl.tfloat32,
        },
        min_partitions=200,
    )
    ht = ht.annotate(locus=hl.locus(ht.chr, ht.pos, reference_genome="GRCh38"))
    ht = ht.annotate(alleles=hl.array([ht.ref, ht.alt]))
    ht = ht.select("locus", "alleles", "primateDL_score", "UCSC_gene")
    ht = ht.key_by("UCSC_gene")

    # import the ucsc id file
    ucsc = hl.import_table(ucsc_ids, min_partitions=200)
    ucsc = ucsc.annotate(
        UCSC_gene=ucsc.hg19_knownToEnsembl_name,
        Transcript_stable_ID=ucsc.hg19_ensGene_transcript,
    )
    ucsc = ucsc.select("UCSC_gene", "Transcript_stable_ID")
    ucsc = ucsc.key_by("UCSC_gene")

    # join the two tables
    ht = ht.join(ucsc, how="left")
    ht = ht.filter(hl.is_defined(ht.Transcript_stable_ID))
    ht = ht.key_by("Transcript_stable_ID")

    # import the Ensembl 105 id file
    canon = hl.import_table(ensembl_ids, min_partitions=200)
    canon = canon.select("Transcript_stable_ID", "Ensembl_Canonical")
    canon = canon.key_by("Transcript_stable_ID")

    # join the two tables
    ht = ht.join(canon, how="inner")

    # get the canonical and non-canonical primateAI scores
    ht = ht.key_by("locus", "alleles")
    ht = ht.annotate(
        primateai_canonical=hl.if_else(
            ht.Ensembl_Canonical == "1", ht.primateDL_score, hl.missing(hl.tfloat)
        ),
        primateai_noncanonical=hl.if_else(
            ht.Ensembl_Canonical != "1", ht.primateDL_score, hl.missing(hl.tfloat)
        ),
    )

    # get the max primateAI score for each variant
    max_primateai_canonical = ht.group_by(*ht.key).aggregate(
        max_primateai_canonical=hl.agg.max(ht.primateai_canonical)
    )
    max_primateai_noncanonical = ht.group_by(*ht.key).aggregate(
        max_primateai_noncanonical=hl.agg.max(ht.primateai_noncanonical)
    )
    primateai = max_primateai_canonical.join(max_primateai_noncanonical, how="inner")
    primateai.checkpoint(
        "gs://gnomad-qin/v4_annotations/primateai-v0.2_processed_max_canonical.ht",
        overwrite=True,
    )


primateai_tsv = "gs://gnomad-qin/v4_annotations/PrimateAI_scores_v0.2_hg38.tsv.bgz"
ucsc_ids = "gs://gnomad-qin/ensembl/ucsc_hg19_knownToEnsembl.tsv.bgz"
ensembl_ids = (
    "gs://gnomad-qin/ensembl/ensembl105_chr1-22XY_MANE_canonical_ucscid.tsv.bgz"
)
create_primateai_grch38_ht(primateai_tsv, ucsc_ids, ensembl_ids)


def remove_missing_vep_fields(vep_expr: hl.StructExpression) -> hl.StructExpression:
    """
    Remove fields from VEP 105 annotations that have been excluded in past releases or are missing in all rows.

    :param vep_expr: StructExpression containing VEP 105 annotations.
    :return: StructExpression containing VEP 105 annotations with missing fields removed.
    """
    vep_expr = vep_expr.drop("colocated_variants", "context")

    vep_expr = vep_expr.annotate(
        transcript_consequences=vep_expr.transcript_consequences.map(
            lambda x: x.drop("minimised", "swissprot", "trembl", "uniparc")
        )
    )

    for consequence in [
        "intergenic_consequences",
        "motif_feature_consequences",
        "regulatory_feature_consequences",
    ]:
        vep_expr = vep_expr.annotate(
            **{consequence: vep_expr[consequence].map(lambda x: x.drop("minimised"))}
        )
    return vep_expr


def replace_oth_with_remaining(ht: hl.Table) -> hl.Table:
    """
    Replace 'oth' with 'remaining' in Global fields of a Table.

    .. note::
        This function renames v3 ancestry groups to match v4 exomes' ancestry groups.
        The value of the key "pop" within `freq_meta` and the keys of `freq_index_dict`
        that contain "oth" are replaced with "remaining".

    :param ht: release sites Table to be modified.
    :return: release sites Table with 'oth' replaced with 'remaining'.
    """
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: x.replace("oth", "remaining"))
        ),
        freq_index_dict=hl.zip(
            ht.freq_index_dict.keys().map(lambda k: k.replace("oth", "remaining")),
            ht.freq_index_dict.values(),
        ),
    )
    return ht


# TODO: drop old_locus, old_alleles, seems to be already dropped in v3.1.4 release table
# TODO: rename ancestry group fields to match v4 exomes
# TODO: drop missing fields from VEP 105 annotations
# TODO: merge vep_ht with v3.1.4 release table
# TODO: rerun inbreeding coefficient with callstats instead of GT calls
# TODO: annotate with newest dbsnp release
# TODO: add Pangolin prediction
# TODO: add zoonomia constraint scores


def main(args):
    """Script to generate release sites ht for v4 genomes."""
    logger.info("Loading release table of v3.1.4...")
    # ht = release_sites(public=True).versions["3.1.4"].ht()
    ht = hl.read_table(
        "gs://gnomad/release/3.1.4/ht/genomes/gnomad.genomes.v3.1.4.sites.ht"
    )

    logger.info("Replacing ancestry group 'oth' with 'remaining'...")
    ht = replace_oth_with_remaining(ht)

    logger.info("Loading VEP-annotated table...")
    vep_ht = get_vep(version=CURRENT_VERSION, data_type="genomes").ht()

    logger.info("Removing missing VEP fields...")
    vep_ht = ht.annotate(vep=remove_missing_vep_fields(vep_ht.vep))

    logger.info("Loading dbsnp table...")
    dbsnp_ht = dbsnp.ht().select("rsid")

    logger.info("Loading pangolin table...")
    pangolin_ht = hl.read_table(
        "gs://gnomad/v4.0/annotations/genomes/gnomad.genomes.v4.0.pangolin.ht"
    )  # TODO: make Pangolin a versioned resource
