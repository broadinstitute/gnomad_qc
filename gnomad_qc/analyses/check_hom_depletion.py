"""
Script to flag and characterize variants depleted of homozygous individuals.

A variant is flagged as hom-depleted when its allele frequency exceeds a dynamic
cutoff -- sqrt(`min_expected` / N) -- yet no homozygotes were observed, overall or
within a genetic ancestry group. Depletion of that kind can reflect true biological
intolerance of biallelic genotypes, but it can equally reflect a technical artifact,
so variants are annotated with external evidence that helps tell the two apart:
constraint, ClinVar, OMIM recessive labels, mouse viability (IMPC/MGI),
structural-variant overlap, and published human knockout and recessive-association
results (Lassen et al. 2026, Heyne et al. 2023, and the BRaVa cohort-level counts).

The output keeps every PASS autosomal variant, flagged or not; filter on the
`*_hom_depletion` fields to get the candidates.
"""

import argparse
import logging
from typing import List

import hail as hl
from gnomad.resources.grch38.gnomad import (
    GEN_ANC_GROUPS,
    constraint,
    public_release,
    sv_sites_vcf_path,
)
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.vep import filter_vep_transcript_csqs_expr
from hail.utils import new_temp_file

from gnomad_qc.analyses.resources import get_hom_depletion
from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.basics import get_logging_path

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hom_depletion_check")
logger.setLevel(logging.INFO)


# TODO: Move resources to more stable path.
OMIM_PATH = "gs://gnomad-kristen/hom_depletion/omim_2023.tsv"
IMPC_PATH = "gs://gnomad-kristen/hom_depletion/viability.csv.gz"  # 2026-03-16 13:49
# accessed 07/07/2026 https://www.informatics.jax.org/homology.shtml
MGI_PATH = "gs://gnomad-kristen/hom_depletion/HOM_AllOrganism.txt"
# Table S4 of Lassen et al. 2026 (AJHG 113:1330-1346, doi:10.1016/j.ajhg.2026.04.005):
# genes with >=1 observed biallelic pLoF genotype ("knockout") in their six-biobank
# meta-analysis (948,690 individuals), compared against previously published
# knockout surveys.
LASSEN_S4_PATH = (
    "gs://gnomad-kristen/hom_depletion/lassen2026_meta_analysis_supp_s4.txt"
)
# Per-survey knockout counts in LASSEN_S4_PATH, read by column name in
# `annotate_lassen_biallelic_ko`. The gnomAD, Sulem, Narasimhan, and Saleheen counts
# are secondhand, taken from the Oddsson et al. 2023 (doi:10.1038/s41467-023-38951-2)
# compilation.
LASSEN_S4_DATASETS = [
    "lassen2026",
    "sulem",
    "narasimhan",
    "saleheen",
    "gnomad_via_oddsson",
    "oddsson",
    "sun",
]
# Supplementary Table 1 of Heyne et al. 2023 (Nature 613:519-525,
# doi:10.1038/s41586-022-05420-7): every coding variant association with P < 1e-7
# under the additive *or* the recessive model, across all FinnGen (release 4)
# disease endpoints. Read with its own header.
HEYNE_S1_PATH = "gs://gnomad-kristen/hom_depletion/heyne_supp1.txt"
# S1 symbols that GENCODE v39 (VEP 105) calls by a newer name, found as MGI human
# synonyms and checked against v39. Not mapped through MGI synonyms at run time:
# some S1 names are also an old synonym of a different gene (e.g. FLG of FGFR1).
HEYNE_S1_RENAMED_GENES = {
    "ALS2CR12": "FLACC1",
    "BHMG1": "MEIOSIN",
    "FAM103A1": "RAMAC",
    "FAM109A": "PHETA1",
    "KDELC2": "POGLUT3",
    "MARCH8": "MARCHF8",
    "OBFC1": "STN1",
    "SEPN1": "SELENON",
}
# Table S5 of Lassen et al. 2026: the 58 significant gene-trait recessive
# associations (39 genes), each classified as Recessive / Additive / Ambiguous.
# Export the S5 sheet to this path as tab-delimited text, as for S4.
LASSEN_S5_PATH = (
    "gs://gnomad-kristen/hom_depletion/lassen2026_meta_analysis_supp_s5.txt"
)
# Cohort-level (pre-meta-analysis) recessive burden results behind Lassen et al.
# 2026, from the BRaVa recessive sumstats browser repo
# (https://github.com/giorkala/recessive-sumstats-browser, data/chunks/pre_meta/):
# one row per (biobank, ancestry, sex, trait, encoding, annotation mask, gene).
# Mirror the 128 gene-bucket .tsv.gz chunks to this prefix.
BRAVA_PRE_META_PATH = "gs://gnomad-kristen/hom_depletion/brava_pre_meta/*.tsv.gz"
# Genetic ancestry labels used in the BRaVa cohort-level results.
BRAVA_ANCESTRIES = ["afr", "amr", "eas", "eur", "sas"]


def get_readable_field_name(group: str) -> str:
    """
    Get the field name for a group (for clarity, "adj" -> "overall").

    :param group: Group key (e.g., "adj", "afr_adj")
    :return: Field name (e.g., "overall", "afr")
    """
    return "overall" if group == "adj" else group.replace("_adj", "")


def table_to_dict_expr(ht: hl.Table) -> hl.expr.DictExpression:
    """
    Collect a small single-key Table into a dict expression of key to row value.

    Looking up a gene table this way avoids the re-key of the variant table by gene
    (and back) that a join on a non-key field requires.

    :param ht: Table keyed by a single field, small enough to broadcast.
    :return: Dict expression mapping each key to its row value (non-key fields).
    """
    return ht.aggregate(
        hl.dict(hl.agg.collect((ht.key[0], ht.row_value))), _localize=False
    )


def annotate_group_freqs_and_cutoffs(
    ht: hl.Table,
    groups: List[str],
    min_expected: int = 3,
) -> hl.Table:
    """
    Annotate a Hail Table with group-specific allele frequencies and AF cutoffs.

    For each group, this function:
      - Retrieves the corresponding frequency index from `ht.freq_index_dict`
      - Annotates `<field_name>_freq` with the group's frequency
      - Annotates `<field_name>_af_cutoff` using:
            sqrt(min_expected / (AN / 2))

    Every group in `groups`, and "raw", must be a key in `freq_index_dict`;
    a missing group raises a KeyError.

    Also annotates `raw_freq` from the "raw" index. The group frequencies above are
    adj-filtered, so their homozygote counts exclude genotypes failing genotype QC
    (GQ/DP/allele balance). Carrying raw alongside them lets a zero adj homozygote
    count be separated into "no homozygous genotypes were called" versus "homozygous
    genotypes were called and then filtered", which are different findings.

    Note that `freq` carries a single all-samples raw group, not one per genetic
    ancestry group, so `raw_freq` is only comparable to the overall (adj) group.

    :param ht: Hail Table containing `freq` (array of frequency structs) and
        `freq_index_dict` (dict mapping group names to indices in `freq`).
    :param groups: List of group keys to use (e.g., ["adj", "afr_adj"]).
    :param min_expected: Minimum expected number of homozygous individuals used to
        compute the AF cutoff. Default is 3, the threshold Arnadottir et al. 2022
        (https://doi.org/10.1038/s41467-022-28330-8) used to define a complete
        homozygote deficit.
    :return: Annotated Hail Table with added frequency and cutoff fields.
    """
    freq_index_dict = hl.eval(ht.freq_index_dict)
    group_freqs = {}
    cutoffs = {}
    metrics = {}

    for group in groups:
        freq_expr = ht.freq[freq_index_dict[group]]
        field_name = get_readable_field_name(group)

        group_freqs[f"{field_name}_freq"] = freq_expr
        cutoffs[f"{field_name}_af_cutoff"] = hl.sqrt(min_expected / (freq_expr.AN / 2))
        metrics[f"{field_name}_expected_hom"] = freq_expr.AF**2 * (freq_expr.AN / 2)

    group_freqs["raw_freq"] = ht.freq[freq_index_dict["raw"]]

    return ht.annotate(**group_freqs, **cutoffs, **metrics)


def flag_hom_depletion(ht: hl.Table, groups: List[str]) -> hl.Table:
    """
    Flag variants with depleted homozygous counts based on a metric and AF cutoff.

    A variant is flagged as hom-depleted if its frequency is greater than the AF
    cutoff and it has zero homozygous counts in the corresponding group.

    :param ht: Hail Table containing group-specific frequency and cutoff annotations.
    :param groups: List of group keys to use (e.g., ["adj", "afr_adj"]).
    :return: Hail Table with added boolean fields indicating hom-depletion for each
        group (e.g., "overall_hom_depletion", "afr_hom_depletion").
    """
    # Dynamic threshold: AF must be greater than the square root of (desired_homs / N).
    flagged_groups = {}
    for group in groups:
        field_name = get_readable_field_name(group)
        freq_field = f"{field_name}_freq"
        flagged_groups[f"{field_name}_hom_depletion"] = (
            hl.is_defined(ht[freq_field])
            & (ht[freq_field].AF > ht[f"{field_name}_af_cutoff"])
            & (ht[freq_field].homozygote_count == 0)
        )

    return ht.annotate(**flagged_groups)


def annotate_most_severe_gene(ht: hl.Table) -> hl.Table:
    """
    Annotate the gene carrying `most_severe_csq`, preferring the best transcript.

    The most severe consequence can fall on transcripts of more than one gene
    (overlapping genes, readthroughs, antisense overlaps), so the transcript is chosen
    by quality rather than by Ensembl VEP's array order: MANE Select, then Ensembl
    canonical, then any Ensembl transcript. This is the same preference
    `annotate_constraint` applies, so the two annotations describe the same transcript.

    Ensembl transcripts are preferred over RefSeq ones because `most_severe_gene_id`
    is joined against BRaVa's Ensembl gene IDs. The last fallback drops that
    restriction, so no variant that resolved to a gene before becomes missing.

    All three fields are missing when no transcript carries the most severe
    consequence.

    :param ht: Hail Table containing `vep.transcript_consequences` and
        `most_severe_csq`.
    :return: Hail Table with `most_severe_gene`, `most_severe_gene_id`, and
        `most_severe_hgnc_id` annotations.
    """

    def _find(**kwargs: bool) -> hl.expr.StructExpression:
        """Find the first Ensembl transcript with the most severe consequence."""
        return filter_vep_transcript_csqs_expr(
            ht.vep.transcript_consequences,
            ensembl_only=True,
            additional_filtering_criteria=[
                lambda tx: tx.consequence_terms.contains(ht.most_severe_csq)
            ],
            **kwargs,
        ).first()

    most_severe_tx = hl.coalesce(
        _find(mane_select=True),
        _find(canonical=True),
        _find(),
        hl.find(
            lambda tx: tx.consequence_terms.contains(ht.most_severe_csq),
            ht.vep.transcript_consequences,
        ),
    )

    return ht.annotate(
        most_severe_gene=most_severe_tx.gene_symbol,
        most_severe_gene_id=most_severe_tx.gene_id,
        most_severe_hgnc_id=most_severe_tx.hgnc_id,
    )


def annotate_omim_recessive_gene_label(
    ht: hl.Table,
    omim_path: str = OMIM_PATH,
) -> hl.Table:
    """
    Annotate whether `most_severe_gene_id` is OMIM-labeled as recessive.

    :param ht: Hail Table containing `most_severe_gene_id`.
    :param omim_path: Path to a tab-delimited OMIM table with `gene_id` and
        `phenotype_inheritance` columns. Default is OMIM_PATH.
    :return: Hail Table with `is_omim_recessive_gene` boolean and `omim` array
        annotations.
    """
    omim_ht = hl.import_table(omim_path, impute=True)

    omim_ht = omim_ht.filter(
        omim_ht.phenotype_inheritance.lower().contains("recessive")
    )

    omim_ht = omim_ht.group_by(gene_id=omim_ht.gene_id).aggregate(
        omim=hl.agg.collect(
            hl.struct(
                phenotype_inheritance=omim_ht.phenotype_inheritance,
                phenotype_description=omim_ht.phenotype_description,
                gene_symbols=omim_ht.gene_symbols,
                gene_description=omim_ht.gene_description,
                comments=omim_ht.comments,
            )
        )
    )

    o = table_to_dict_expr(omim_ht).get(ht.most_severe_gene_id)
    return ht.annotate(
        omim=o.omim,
        is_omim_recessive_gene=hl.is_defined(o),
    )


def annotate_gerp(ht: hl.Table) -> hl.Table:
    """Annotate GERP scores by locus."""
    gerp_ht = hl.experimental.load_dataset(
        name="gerp_scores", version="hg19", reference_genome="GRCh38"
    )

    # -12.3 to 6.17 is the range of GERP values where 6.17 is the most conserved.
    return ht.annotate(gerp=gerp_ht[ht.locus].S)


def annotate_constraint(ht: hl.Table) -> hl.Table:
    """
    Annotate gene-level constraint metrics for `most_severe_gene_id`.

    Pulls LoF, missense, and synonymous z-scores, observed/expected (oe) ratios,
    and oe confidence intervals, plus constraint and gene flags. One transcript is
    selected per gene, preferring MANE Select and falling back to the canonical
    transcript. The join is keyed on Ensembl gene ID.

    :param ht: Hail Table containing `most_severe_gene_id` (Ensembl gene ID).
    :return: Hail Table with added `constraint` struct annotation.
    """
    constraint_ht = constraint().ht()

    # Restrict to Ensembl genes.
    constraint_ht = constraint_ht.filter(constraint_ht.gene_id.startswith("ENSG"))

    # One row per gene: prefer MANE Select, fall back to the canonical transcript.
    mane_ht = constraint_ht.filter(constraint_ht.mane_select).key_by("gene_id")
    canonical_ht = constraint_ht.filter(constraint_ht.canonical).key_by("gene_id")
    canonical_ht = canonical_ht.anti_join(mane_ht)
    constraint_ht = mane_ht.union(canonical_ht)

    # Keep only the fields read below, so the broadcast lookup dict stays small.
    constraint_ht = constraint_ht.select(
        "lof",
        "mis",
        "syn",
        "constraint_flags",
        "gene_flags",
        "canonical",
        "mane_select",
    )
    c = table_to_dict_expr(constraint_ht).get(ht.most_severe_gene_id)

    def _metrics(x: hl.expr.StructExpression) -> hl.expr.StructExpression:
        """Select z-score, o/e, and o/e CI bounds from a constraint metric struct."""
        return hl.struct(
            z_score=x.z_score,
            oe=x.oe,
            oe_ci_lower=x.oe_ci.lower,
            oe_ci_upper=x.oe_ci.upper,
        )

    return ht.annotate(
        constraint=hl.struct(
            lof=_metrics(c.lof),
            mis=_metrics(c.mis),
            syn=_metrics(c.syn),
            constraint_flags=c.constraint_flags,
            gene_flags=c.gene_flags,
            canonical=c.canonical,
            mane_select=c.mane_select,
        )
    )


def annotate_clinvar(ht: hl.Table) -> hl.Table:
    """
    Annotate per-variant ClinVar clinical significance and phenotype (disease) info.

    Joins on locus/alleles and pulls disease names (CLNDN), clinical significance
    (CLNSIG), and review status (CLNREVSTAT) for any ClinVar variant. Variants not
    present in ClinVar receive a missing `clinvar` struct.

    :param ht: Hail Table keyed by locus/alleles.
    :return: Hail Table with a `clinvar` struct annotation.
    """
    clinvar_ht = clinvar.ht()

    clinvar_row = clinvar_ht[ht.key]
    info = clinvar_row.info
    return ht.annotate(
        clinvar=hl.or_missing(
            hl.is_defined(clinvar_row),
            hl.struct(
                disease=info.CLNDN,  # phenotype / disease name(s)
                clinical_significance=info.CLNSIG,  # benign, pathogenic, VUS, etc.
                review_status=info.CLNREVSTAT,
            ),
        )
    )


def get_mouse_to_human_ortholog_ht(mgi_path: str = MGI_PATH) -> hl.Table:
    """
    Build a mouse->human ortholog map from MGI's HomoloGene HOM_AllOrganism report.

    Genes sharing a `DB Class Key` are homologs across organisms. Grouping by that key
    and separating human (taxon 9606) and mouse (taxon 10090) members yields, per class,
    the human HGNC IDs and mouse MGI IDs. Result is keyed by mouse MGI ID (`mouse_id`),
    matching IMPC's `Gene Accession Id`, and records the human ortholog(s) plus an
    `ortholog_type` flag ("one2one" vs "one2many") so non-1:1 mappings can be handled
    downstream.

    :param mgi_path: Path to MGI HOM_AllOrganism.txt (tab-delimited).
        Default is MGI_PATH.
    :return: Table keyed by mouse MGI ID (`mouse_id`) with `human_hgnc_ids` (set) and
        `ortholog_type`.
    """
    homer = hl.import_table(mgi_path, impute=True)

    by_class = homer.group_by(class_key=homer["DB Class Key"]).aggregate(
        human_hgnc_ids=hl.agg.filter(
            homer["NCBI Taxon ID"] == 9606, hl.agg.collect_as_set(homer["HGNC ID"])
        ),
        mouse_ids=hl.agg.filter(
            homer["NCBI Taxon ID"] == 10090,
            hl.agg.collect_as_set(homer["Mouse MGI ID"]),
        ),
    )

    # Keep classes with both a human and a mouse member.
    by_class = by_class.filter(
        (hl.len(by_class.human_hgnc_ids) > 0) & (hl.len(by_class.mouse_ids) > 0)
    )
    by_class = by_class.annotate(
        ortholog_type=hl.if_else(
            (hl.len(by_class.human_hgnc_ids) == 1) & (hl.len(by_class.mouse_ids) == 1),
            "one2one",
            "one2many",
        )
    )

    # One row per mouse MGI ID, keyed for joining to the IMPC report.
    ortholog_ht = by_class.annotate(mouse_id=hl.array(by_class.mouse_ids)).explode(
        "mouse_id"
    )

    return ortholog_ht.key_by("mouse_id").select("human_hgnc_ids", "ortholog_type")


def annotate_impc_viability(
    ht: hl.Table, impc_path: str = IMPC_PATH, mgi_path: str = MGI_PATH
) -> hl.Table:
    """
    Annotate mouse-knockout homozygous viability (IMPC) on `most_severe_hgnc_id`.

    Viability calls (viable / subviable / lethal) are mapped from mouse to human genes
    via MGI HomoloGene orthologs. `ortholog_type` records whether the mapping was 1:1;
    non-one2one genes (duplicated families) should be interpreted with care.

    :param ht: Hail Table containing `most_severe_hgnc_id` (HGNC ID).
    :param impc_path: Path to the IMPC viability.csv.gz report. Default is IMPC_PATH.
    :param mgi_path: Path to MGI HOM_AllOrganism.txt ortholog report.
        Default is MGI_PATH.
    :return: Hail Table with an `impc` struct annotation.
    """
    impc_ht = hl.import_table(
        impc_path, impute=True, force=True, delimiter=",", quote='"'
    )
    ortholog_ht = get_mouse_to_human_ortholog_ht(mgi_path)

    # Match on MGI ID: IMPC still uses some old symbols (e.g. Car12, now Ca12).
    ortholog = ortholog_ht[impc_ht["Gene Accession Id"]]
    impc_ht = impc_ht.annotate(
        human_hgnc_ids=ortholog.human_hgnc_ids,
        ortholog_type=ortholog.ortholog_type,
        viability_call=impc_ht["Viability Phenotype HOMs/HEMIs"],
    )

    # Keep genes with a human ortholog, then explode to one row per human HGNC ID.
    impc_ht = impc_ht.filter(hl.is_defined(impc_ht.human_hgnc_ids))
    impc_ht = impc_ht.annotate(human_hgnc_id=hl.array(impc_ht.human_hgnc_ids)).explode(
        "human_hgnc_id"
    )

    impc_ht = impc_ht.group_by(hgnc_id=impc_ht.human_hgnc_id).aggregate(
        viability=hl.agg.collect_as_set(impc_ht.viability_call),
        ortholog_type=hl.agg.collect_as_set(impc_ht.ortholog_type),
    )

    # Match on HGNC ID: MGI has current human symbols, some renamed since VEP 105
    # (e.g. MAPDA, which VEP 105 calls ADAL).
    impc = table_to_dict_expr(impc_ht).get(ht.most_severe_hgnc_id)
    return ht.annotate(
        impc=hl.struct(
            viability=impc.viability,
            ortholog_type=impc.ortholog_type,
            is_mouse_ko_lethal=hl.any(
                lambda x: x.lower().contains("lethal"), hl.array(impc.viability)
            ),
        )
    )


def annotate_lassen_biallelic_ko(
    ht: hl.Table,
    s4_path: str = LASSEN_S4_PATH,
) -> hl.Table:
    """
    Annotate whether a biallelic pLoF knockout of `most_severe_gene` is known.

    Uses Lassen et al. 2026 table S4 (see LASSEN_S4_PATH). An observed knockout
    suggests the gene tolerates biallelic LoF, so depletion is more likely
    allele-specific or technical. `ko_novel_to_lassen2026` marks knockouts first
    seen in that study.

    Absence from S4 conflates intolerance with low power; pair it with constraint
    and mouse viability. S4 is pLoF-only (says nothing about missense alleles).
    Compound hets are phased (AOU: homozygotes only). Genes with a Lassen 2026
    knockout are matched on Ensembl ID and the rest on gene symbol, so older-survey
    aliases appear absent.

    :param ht: Hail Table containing `most_severe_gene` and `most_severe_gene_id`.
    :param s4_path: Path to the tab-delimited export of Lassen et al. 2026 table S4.
        Default is LASSEN_S4_PATH.
    :return: Hail Table with a `lassen2026_ko` struct annotation.
    """
    # Read by column name, so a reordered export can't shift counts between surveys
    # and a renamed column fails here.
    s4 = hl.import_table(s4_path, delimiter="\t")
    # Drop the TOTAL row and the export's blank trailing rows. Rows with Gene "NA"
    # (read as missing) are known only by Ensembl ID and are kept.
    s4 = s4.filter(~hl.literal({"TOTAL", ""}).contains(s4.Gene))

    def _count(field: hl.expr.StringExpression) -> hl.expr.Int32Expression:
        """Parse an S4 count, treating blank or non-numeric cells as 0."""
        return hl.or_else(hl.parse_int32(field), 0)

    s4 = s4.select(
        gene_symbol=s4.Gene,
        # Set only for genes with a Lassen 2026 knockout.
        gene_id=s4["Ensembl ID"],
        lassen2026=_count(s4["This study"]),
        sulem=_count(s4["Sulem et al"]),
        narasimhan=_count(s4["Narasimhan et al"]),
        saleheen=_count(s4["Saleheen et al"]),
        # The header has a trailing space here.
        gnomad_via_oddsson=_count(s4["gnomAD "]),
        oddsson=_count(s4["Oddsson et al"]),
        sun=_count(s4["Sun et al"]),
    )

    # Look up both keys (import_table reads S4's "NA" as missing). S4 lists 69
    # genes twice: older surveys under the old symbol, without an ID, and
    # Lassen 2026 under a newer symbol (or none), with one. Their counts are merged.
    # A symbol row with its own ID is used only when no row matches by ID, since it
    # may belong to another Ensembl gene with the same symbol (e.g. two TBCE genes).
    by_id = table_to_dict_expr(
        s4.filter(hl.is_defined(s4.gene_id))
        .key_by("gene_id")
        .select(*LASSEN_S4_DATASETS)
    )
    by_symbol = table_to_dict_expr(
        s4.filter(hl.is_defined(s4.gene_symbol))
        .key_by("gene_symbol")
        .select("gene_id", *LASSEN_S4_DATASETS)
    )
    id_rec = by_id.get(ht.most_severe_gene_id)
    symbol_rec = by_symbol.get(ht.most_severe_gene)
    symbol_rec = hl.or_missing(
        hl.is_missing(id_rec) | hl.is_missing(symbol_rec.gene_id), symbol_rec
    )
    # A split gene's two rows count different surveys, so the max is their union.
    counts = hl.or_missing(
        hl.is_defined(id_rec) | hl.is_defined(symbol_rec),
        hl.struct(
            **{
                d: hl.max(hl.or_else(id_rec[d], 0), hl.or_else(symbol_rec[d], 0))
                for d in LASSEN_S4_DATASETS
            }
        ),
    )

    ko_datasets = hl.array(
        [hl.or_missing(counts[dataset] > 0, dataset) for dataset in LASSEN_S4_DATASETS]
    ).filter(hl.is_defined)

    return ht.annotate(
        lassen2026_ko=hl.struct(
            gene_in_lassen_s4=hl.is_defined(counts),
            biallelic_plof_ko_observed=hl.or_else(hl.len(ko_datasets) > 0, False),
            n_datasets_observing_ko=hl.or_else(hl.len(ko_datasets), 0),
            datasets_observing_ko=ko_datasets,
            # Per-survey observed-knockout counts, missing when the gene is absent
            # from S4. A single observed knockout is weaker evidence of tolerated
            # biallelic loss than a recurrent one.
            ko_counts=counts,
            ko_observed_in_gnomad_per_oddsson=hl.or_else(
                counts.gnomad_via_oddsson > 0, False
            ),
            # Same as S4's "Novel" column on every S4 row, but derived from the merged
            # counts so a split gene isn't novel when older surveys saw its knockout.
            ko_novel_to_lassen2026=hl.or_else(
                (counts.lassen2026 > 0) & (hl.len(ko_datasets) == 1), False
            ),
        )
    )


def annotate_lassen_recessive_assoc(
    ht: hl.Table,
    s5_path: str = LASSEN_S5_PATH,
) -> hl.Table:
    """
    Annotate published recessive gene-trait associations for `most_severe_gene_id`.

    Uses table S5 of Lassen et al. 2026: the significant associations from their
    six-biobank recessive analysis, each labeled Recessive, Additive, or Ambiguous
    according to whether the recessive model fit better than the additive one.

    An association is evidence that biallelic carriers of the gene both exist and
    reach adulthood with a measurable phenotype, so -- like the S4 knockout tally --
    it argues against absolute intolerance of biallelic loss while marking the gene
    as disease relevant. The tables are small (58 associations over 39 genes), so
    most genes get `has_recessive_assoc=False`; this is a lookup for the occasional
    hit, not a broadly informative annotation. False there means only that S5 holds
    no significant *recessive* association for the gene -- S5 never contained
    additive-only hits, so it is not evidence the gene is unassociated.

    :param ht: Hail Table containing `most_severe_gene_id`.
    :param s5_path: Path to the tab-delimited export of Lassen et al. 2026 table S5.
        Default is LASSEN_S5_PATH.
    :return: Hail Table with a `lassen2026_recessive_assoc` struct annotation.
    """
    # Read by column name, as for S4.
    s5 = hl.import_table(s5_path, delimiter="\t")
    # Drop the export's blank trailing rows.
    s5 = s5.filter(s5["ENSEMBL ID"] != "")
    s5 = s5.select(
        # Match on Ensembl ID, which S5 gives for every association.
        gene_id=s5["ENSEMBL ID"],
        phenotype=s5.Phenotype,
        mode=s5["Mode (Sup Note 7)"],
        p_rec=hl.parse_float64(s5["P rec"]),
        p_add=hl.parse_float64(s5["P add"]),
    )
    s5 = s5.group_by(gene_id=s5.gene_id).aggregate(
        associations=hl.agg.collect(
            hl.struct(
                phenotype=s5.phenotype,
                mode=s5.mode,
                p_rec=s5.p_rec,
                p_add=s5.p_add,
            )
        )
    )

    rec = table_to_dict_expr(s5).get(ht.most_severe_gene_id)
    associations = hl.or_else(
        rec.associations, hl.empty_array(rec.associations.dtype.element_type)
    )

    return ht.annotate(
        lassen2026_recessive_assoc=hl.struct(
            # S5 holds only recessive-model hits, so False means "no significant
            # recessive association", not "no association": a gene with an
            # additive-only signal is absent from S5 entirely.
            has_recessive_assoc=hl.len(associations) > 0,
            n_associations=hl.len(associations),
            # True if any S5 trait is labeled "Recessive" (Lassen Note S7:
            # p_rec < p_add/100 and either beta_rec/beta_add > 1.5 or p_add > 0.05).
            # Stricter than Heyne's `has_recessive_specific`, which has no beta check
            # and passes when p_add is missing.
            has_recessive_specific=associations.any(lambda a: a.mode == "Recessive"),
            modes=hl.set(associations.map(lambda a: a.mode)),
            phenotypes=associations.map(lambda a: a.phenotype),
            associations=associations,
        )
    )


def annotate_brava_cohort_carriers(
    ht: hl.Table,
    pre_meta_path: str = BRAVA_PRE_META_PATH,
) -> hl.Table:
    """
    Annotate per-cohort biallelic-carrier counts from the BRaVa recessive analysis.

    Uses the cohort-level (pre-meta) results behind Lassen et al. 2026 (UKB, AOU, GEL,
    GNH, BioMe, BBJ) to show which ancestries have biallelic carriers of the gene. GNH
    has high parental relatedness, so no carriers there is stronger evidence than
    gnomAD alone.

    Recessive encoding is 0/0/1, so AC = number of biallelic individuals (max over
    traits). Counts are summed across biobanks within an ancestry, never across
    masks, which are nested (pLoF is a subset of nonsynonymous).

    A missing row means fewer than 5 carriers, not none. GEL reports no AC, and AOU
    reports "<40" (5-39); neither is summed. `cohorts_missing_counts` lists which of
    `cohorts_tested` lack usable counts, and `ac_suppressed` in `cohorts` separates
    AOU suppression from GEL blanks. Check `ancestries_counted` before reading
    `<anc>_n_biallelic`. Cohort p-values are deliberately not surfaced.

    :param ht: Hail Table containing `most_severe_gene_id` (Ensembl gene ID).
    :param pre_meta_path: Glob of the BRaVa pre-meta gene-bucket .tsv.gz chunks.
        Default is BRAVA_PRE_META_PATH.
    :return: Hail Table with a `brava_carriers` struct annotation.
    """
    pm = hl.import_table(pre_meta_path, delimiter="\t", force=True)

    # Recessive-encoded, all-sex rows for damaging masks. Synonymous is excluded
    # here (it is a calibration mask, not carrier evidence).
    pm = pm.filter(
        (pm.ENCODING == "R") & (pm.SEX == "ALL") & (pm.ANNOTATION != "synonymous")
    )
    # Two ways a count can be unusable: AOU suppresses small counts as "<40" (5-39
    # carriers), and GEL leaves AC blank. Both add 0 to the sums below, and
    # `counts_missing` flags them so a gene tested there isn't read as carrier-free.
    pm = pm.select(
        gene_id=pm.ID,
        biobank=pm.BIOBANK,
        ancestry=pm.ANCESTRY,
        annotation=pm.ANNOTATION,
        ac=hl.parse_float64(pm.AC),
        ac_suppressed=pm.AC.startswith("<"),
        ac_reported=(pm.AC != "") & (pm.AC != "NA"),
        n=hl.parse_int32(pm.N),
    )

    # One row per (gene, biobank, ancestry, mask): max over traits dedupes the
    # per-trait N wobble.
    pm = pm.group_by(pm.gene_id, pm.biobank, pm.ancestry, pm.annotation).aggregate(
        ac=hl.agg.max(pm.ac),
        n=hl.agg.max(pm.n),
        ac_suppressed=hl.agg.any(pm.ac_suppressed),
        ac_reported=hl.agg.any(pm.ac_reported),
    )
    pm = pm.annotate(
        n_biallelic=hl.int32(hl.or_else(pm.ac, 0.0)),
        counts_missing=~pm.ac_reported | pm.ac_suppressed,
    )

    def _ancestry_sums(mask_filter) -> dict:
        """Per-ancestry biallelic-carrier sums across biobanks for one mask."""
        return {
            f"{anc}_n_biallelic": hl.agg.filter(
                mask_filter & (pm.ancestry == anc.upper()),
                hl.agg.sum(pm.n_biallelic),
            )
            for anc in BRAVA_ANCESTRIES
        }

    pm = pm.group_by(pm.gene_id).aggregate(
        cohorts=hl.agg.collect(
            hl.struct(
                biobank=pm.biobank,
                ancestry=pm.ancestry,
                annotation=pm.annotation,
                n_biallelic=pm.n_biallelic,
                n=pm.n,
                counts_missing=pm.counts_missing,
                ac_suppressed=pm.ac_suppressed,
            )
        ),
        # (biobank, ancestry) pairs that tested this gene, and the subset of those
        # whose carrier counts are unavailable. If a candidate's driving ancestry
        # appears only in the latter, the carrier sums are uninformative for it.
        cohorts_tested=hl.agg.collect_as_set(hl.tuple([pm.biobank, pm.ancestry])),
        cohorts_missing_counts=hl.agg.filter(
            pm.counts_missing,
            hl.agg.collect_as_set(hl.tuple([pm.biobank, pm.ancestry])),
        ),
        # Broadest damaging mask: any two damaging alleles.
        **_ancestry_sums(pm.annotation == "nonsynonymous"),
        # pLoF-only, comparable to the S4 knockout tally but ancestry-resolved.
        **{f"plof_{k}": v for k, v in _ancestry_sums(pm.annotation == "pLoF").items()},
    )

    rec = table_to_dict_expr(pm).get(ht.most_severe_gene_id)

    return ht.annotate(
        brava_carriers=hl.struct(
            gene_in_brava=hl.is_defined(rec),
            cohorts=rec.cohorts,
            cohorts_tested=rec.cohorts_tested,
            cohorts_missing_counts=rec.cohorts_missing_counts,
            # Ancestries tested in at least one cohort that reports carrier counts,
            # i.e. those for which `<anc>_n_biallelic` is interpretable.
            ancestries_counted=hl.or_missing(
                hl.is_defined(rec),
                hl.set(
                    rec.cohorts.filter(lambda c: ~c.counts_missing).map(
                        lambda c: c.ancestry
                    )
                ),
            ),
            **{
                f"{anc}_n_biallelic": rec[f"{anc}_n_biallelic"]
                for anc in BRAVA_ANCESTRIES
            },
            **{
                f"plof_{anc}_n_biallelic": rec[f"plof_{anc}_n_biallelic"]
                for anc in BRAVA_ANCESTRIES
            },
        )
    )


def annotate_heyne_assoc(
    ht: hl.Table,
    s1_path: str = HEYNE_S1_PATH,
    max_p: float = 1e-7,
) -> hl.Table:
    """
    Annotate FinnGen coding-variant associations for `most_severe_gene`.

    Uses Heyne et al. 2023 table S1 (see HEYNE_S1_PATH): every coding variant
    reaching P < 1e-7 under the additive or recessive model across all FinnGen
    disease endpoints. The join is on gene symbol, so a hit usually shows the gene is
    disease relevant, not this exact variant (see `this_variant_associated`). A large
    `max_ac_hom` shows the gene tolerates homozygosity for some allele, though not
    necessarily this one.

    S1 uses some gene symbols older than VEP 105's: the 8 in
    `HEYNE_S1_RENAMED_GENES` are switched to their newer names before the join. 5
    more (RP11-/RP5- clone names and FAM206BP) have no known mapping and never match.

    `has_recessive_specific` is the strict tier: the gene has an association passing
    Heyne's 100x rule (recessive p at least two orders of magnitude below additive),
    the criterion behind the paper's 124 recessive-specific associations.

    :param ht: Hail Table containing `most_severe_gene`.
    :param s1_path: Path to the tab-delimited export of Heyne et al. 2023 table S1.
        Default is HEYNE_S1_PATH.
    :param max_p: Significance threshold applied to both models; rows passing
        neither are dropped. Default is 1e-7, S1's stated pre-filter, which 41 of
        its 2,167 rows still fail. The paper's genome-wide threshold is 5e-8.
    :return: Hail Table with a `heyne_results` struct annotation.
    """
    # Both p-value columns can be empty (the additive-or-recessive filter), and
    # `parse_float64` maps those to missing, which is what the aggregations expect.
    s1 = hl.import_table(s1_path, delimiter="\t")
    s1 = s1.select(
        gene_symbol=hl.literal(HEYNE_S1_RENAMED_GENES).get(s1.gene, s1.gene),
        phenotype=s1.phenotype,
        variant=s1.variant,
        category=s1["variant category"],
        p_add=hl.parse_float64(s1["additive pvalue"]),
        p_rec=hl.parse_float64(s1["recessive pvalue"]),
        ac_hom=hl.parse_int32(s1.AC_Hom),
    )
    # Keep only rows passing `max_p`; every association downstream passed it in at
    # least one model.
    s1 = s1.annotate(
        add_significant=hl.is_defined(s1.p_add) & (s1.p_add < max_p),
        rec_significant=hl.is_defined(s1.p_rec) & (s1.p_rec < max_p),
    )
    # The paper's recessive-specific ("more-likely recessive") criterion: recessive
    # significance plus a recessive p at least 100-fold below the additive p (the
    # Heyne 100x rule). S1 blanks the additive p when it was not extracted (its
    # additive results only cover P < 1e-4), so a missing additive p is treated as
    # passing -- the permissive reading. Checked against table S3, which records a
    # 1e-4 ceiling instead of blanks: same qualifying variant set (39), a handful
    # of borderline extra associations genome-wide.
    s1 = s1.annotate(
        rec_specific=s1.rec_significant
        & (hl.is_missing(s1.p_add) | (s1.p_rec < s1.p_add / 100)),
    )
    s1 = s1.filter(s1.add_significant | s1.rec_significant)
    s1 = s1.group_by(gene_symbol=s1.gene_symbol).aggregate(
        associations=hl.agg.collect(
            hl.struct(
                phenotype=s1.phenotype,
                variant=s1.variant,
                category=s1.category,
                p_add=s1.p_add,
                p_rec=s1.p_rec,
                add_significant=s1.add_significant,
                rec_significant=s1.rec_significant,
                rec_specific=s1.rec_specific,
                ac_hom=s1.ac_hom,
            )
        )
    )

    hey = table_to_dict_expr(s1).get(ht.most_severe_gene)
    associations = hl.or_else(
        hey.associations, hl.empty_array(hey.associations.dtype.element_type)
    )

    # This row's variant in S1's "chrom:pos:ref:alt" form (GRCh38, no "chr").
    variant_str = hl.delimit(
        [
            ht.locus.contig.replace("chr", ""),
            hl.str(ht.locus.position),
            ht.alleles[0],
            ht.alleles[1],
        ],
        ":",
    )

    has_additive = associations.any(lambda a: a.add_significant)
    has_recessive = associations.any(lambda a: a.rec_significant)

    return ht.annotate(
        heyne_results=hl.struct(
            has_association=hl.len(associations) > 0,
            n_associations=hl.len(associations),
            has_additive_assoc=has_additive,
            has_recessive_assoc=has_recessive,
            # The strict tier: at least one association passing the 100x rule --
            # Heyne's "recessive-specific" set (the paper's 124 associations over
            # 39 variants). Looser than Lassen's `has_recessive_specific`; see
            # `annotate_lassen_recessive_assoc`.
            has_recessive_specific=associations.any(lambda a: a.rec_specific),
            # Every kept row passed `max_p` in at least one model, so a gene with
            # any association is always one of the three labels; genes with none
            # are missing here (and False above).
            assoc_models=hl.or_missing(
                hl.len(associations) > 0,
                hl.case()
                .when(has_additive & has_recessive, "both")
                .when(has_additive, "additive")
                .default("recessive"),
            ),
            this_variant_associated=associations.any(
                lambda a: a.variant == variant_str
            ),
            min_p_rec=hl.min(associations.map(lambda a: a.p_rec)),
            min_p_add=hl.min(associations.map(lambda a: a.p_add)),
            # Homozygotes seen for any associated variant in the gene.
            max_ac_hom=hl.max(associations.map(lambda a: a.ac_hom)),
            phenotypes=hl.set(associations.map(lambda a: a.phenotype)),
            associations=associations,
        )
    )


def import_sv_sites_ht(sv_vcf_path: str = sv_sites_vcf_path()) -> hl.Table:
    """
    Import the gnomAD SV sites VCF as a rows-only Hail Table.

    :param sv_vcf_path: Path to the gnomAD SV sites VCF (bgzipped).
        Default is the current release from `sv_sites_vcf_path`.
    :return: Hail Table of SV sites (one row per SV), for use with
        `annotate_sv_cnv_overlap`.
    """
    mt = hl.import_vcf(
        sv_vcf_path,
        reference_genome="GRCh38",
        force_bgz=True,
    )
    return mt.rows()


def annotate_sv_cnv_overlap(
    ht: hl.Table,
    sv_ht: hl.Table,
    use_grpmax: bool = True,  # DUP frequency: prefer max-over-ancestries AF
    require_pass_sv: bool = False,  # if True, only consider PASS SVs
) -> hl.Table:
    """
    Annotate `ht.sv_cnv` with raw per-variant SV/CNV overlap evidence (no cutoffs).

    For each SNV, aggregates over every overlapping copy-number-gain-capable SV
    (biallelic DUP or multiallelic CNV) and records the raw signals needed to judge
    whether an apparent homozygote deficit is a copy-number artifact rather than
    biology. A DUP/CNV inflates local copy number, so true homozygotes are rarely
    observed -- a fake hom deficit. DEL is excluded on purpose (hemizygous loss is a
    different mechanism).

    No thresholds are applied and no rows are removed: this is annotation only. The
    fields are maxima/minima over the overlapping SVs, so a downstream caller can
    reproduce any per-SV `any(...)` test exactly -- `max_x >= t` equals
    `any(x >= t)` and `min_x < t` equals `any(x < t)`. Fallbacks for no-overlap
    variants are chosen to fail every such test (max -> 0, min -> 1).

    `sv_cnv` struct fields:
      - n_overlaps            number of overlapping DUP/CNV SVs
      - sv_ids, svtypes       their ids and SVTYPEs
      - max_dup_af            biallelic DUP: max GRPMAX_AF (fallback AF[0]) over hits
      - max_mcnv_gain_freq    mCNV: max frequency of gain states (CN>2)
      - max_mcnv_modal_cn     mCNV: max modal copy number
      - min_mcnv_cn2_freq     mCNV: min diploid (CN=2) frequency
      - any_lowconf_dup       any hit flagged LOW_CONFIDENCE_REPETITIVE_LARGE_DUP
      - any_gt_overdispersed  any hit flagged PESR_GT_OVERDISPERSION (informational)

    :param ht: Hail Table keyed by locus (the SNVs to annotate).
    :param sv_ht: SV sites Table, e.g. from `import_sv_sites_ht`.
    :param use_grpmax: DUP frequency uses GRPMAX_AF (fallback AF[0]) if True, else
        AF[0].
    :param require_pass_sv: if True, restrict to PASS SVs before computing overlaps.
    :return: `ht` with an added `sv_cnv` struct of raw overlap evidence.
    """
    rg = ht.locus.dtype.reference_genome

    # gain-carrying + multiallelic classes; drop cross-chrom / undefined spans (BND/CTX)
    i0 = sv_ht.info
    sv = sv_ht.filter((i0.SVTYPE == "DUP") | (i0.SVTYPE == "CNV") | i0.MULTIALLELIC)
    sv = sv.filter(hl.is_defined(sv.info.END) & (sv.info.END > sv.locus.position))
    if require_pass_sv:
        sv = sv.filter(hl.len(sv.filters) == 0)
    info = sv.info
    is_mcnv = info.MULTIALLELIC

    # --- biallelic DUP frequency (grpmax, falling back to AF[0]); safe on empty AF ---
    af0 = hl.if_else(hl.len(info.AF) > 0, info.AF[0], hl.missing(hl.tfloat64))
    dup_af_raw = hl.or_else(info.GRPMAX_AF, af0) if use_grpmax else af0
    dup_af = hl.if_else(~is_mcnv, dup_af_raw, hl.missing(hl.tfloat64))

    # --- mCNV: CN_STATUS holds the CN value at each index, CN_FREQ its frequency ---
    cn_pairs = hl.zip(info.CN_STATUS, info.CN_FREQ).filter(
        lambda p: hl.is_defined(p[0]) & hl.is_defined(p[1])
    )
    gain_freq = hl.if_else(
        is_mcnv,
        hl.sum(cn_pairs.filter(lambda p: p[0] > 2).map(lambda p: p[1])),
        hl.missing(hl.tfloat64),
    )
    cn2_freq = hl.if_else(
        is_mcnv,
        hl.sum(cn_pairs.filter(lambda p: p[0] == 2).map(lambda p: p[1])),
        hl.missing(hl.tfloat64),
    )
    modal_cn = hl.if_else(
        is_mcnv & (hl.len(cn_pairs) > 0),
        hl.sorted(cn_pairs, key=lambda p: -p[1])[0][0],
        hl.missing(hl.tint32),
    )

    sv = sv.annotate(
        _is_mcnv=is_mcnv,
        _dup_af=dup_af,
        _gain_freq=gain_freq,
        _cn2_freq=cn2_freq,
        _modal_cn=modal_cn,
        _lowconf_dup=hl.or_else(info.LOW_CONFIDENCE_REPETITIVE_LARGE_DUP, False),
        _gt_overdisp=hl.or_else(info.PESR_GT_OVERDISPERSION, False),
        _sv_id=sv.rsid,
        _svtype=info.SVTYPE,
    )

    # --- interval-key for one-to-many overlap, then annotate every SNV ---
    sv = sv.annotate(
        _interval=hl.locus_interval(
            sv.locus.contig,
            sv.locus.position,
            sv.info.END,
            includes_start=True,
            includes_end=True,
            reference_genome=rg,
        )
    )
    sv = sv.key_by("_interval").select(
        "_is_mcnv",
        "_dup_af",
        "_gain_freq",
        "_cn2_freq",
        "_modal_cn",
        "_lowconf_dup",
        "_gt_overdisp",
        "_sv_id",
        "_svtype",
    )

    # Checkpoint to avoid a decode NullPointerException in the interval join.
    sv = sv.checkpoint(new_temp_file("sv_cnv_intervals", "ht"))

    ht = ht.annotate(_hits=sv.index(ht.locus, all_matches=True))
    hits = ht._hits
    is_mcnv_s = lambda s: hl.or_else(s._is_mcnv, False)
    dup_hits = hits.filter(lambda s: ~is_mcnv_s(s))  # biallelic DUP hits
    mcnv_hits = hits.filter(is_mcnv_s)  # multiallelic CNV hits

    # Raw aggregates only -- all thresholding is left to the caller. Fallbacks make
    # a no-overlap variant fail every downstream cutoff (max -> 0, min -> 1); e.g.
    # `max_dup_af >= t` reproduces the old `any(dup_af >= t)` exactly.
    ht = ht.annotate(
        sv_cnv=hl.struct(
            n_overlaps=hl.len(hits),
            sv_ids=hits.map(lambda s: s._sv_id),
            svtypes=hits.map(lambda s: s._svtype),
            max_dup_af=hl.or_else(hl.max(dup_hits.map(lambda s: s._dup_af)), 0.0),
            max_mcnv_gain_freq=hl.or_else(
                hl.max(mcnv_hits.map(lambda s: s._gain_freq)), 0.0
            ),
            max_mcnv_modal_cn=hl.or_else(
                hl.max(mcnv_hits.map(lambda s: s._modal_cn)), 0
            ),
            min_mcnv_cn2_freq=hl.or_else(
                hl.min(mcnv_hits.map(lambda s: s._cn2_freq)), 1.0
            ),
            any_lowconf_dup=hl.or_else(hits.any(lambda s: s._lowconf_dup), False),
            any_gt_overdispersed=hl.or_else(hits.any(lambda s: s._gt_overdisp), False),
        )
    ).drop("_hits")

    return ht


def select_group_output_fields(ht: hl.Table, groups: List[str]) -> hl.Table:
    """
    Select and flatten group-specific annotations.

    :param ht: Hail Table containing:

        - `<field_name>_freq` structs with AC, AN, AF, homozygote_count
        - `<field_name>_hom_depletion` boolean fields
        - `<field_name>_af_cutoff` fields
        - `raw_freq` struct with AC, AN, AF, homozygote_count
        - `most_severe_csq`, `most_severe_gene`, `lcr_or_segdup`
        - `lassen2026_ko` struct of observed-human-knockout evidence
        - `brava_carriers` struct of per-cohort biallelic-carrier counts
        - `heyne_results` and `lassen2026_recessive_assoc` structs of gene-level
          associations

        All four are written to the output nested under `association_studies`:
        `lassen2026_ko` and `brava_carriers` under `biallelic_carriers`,
        `heyne_results` and `lassen_results` (from `lassen2026_recessive_assoc`)
        under `phenotype_assoc`.
    :param groups: List of group keys (e.g., ["adj", "afr_adj"]).
    :return: Hail Table with selected and flattened fields.
    """
    output_fields = {}

    # per-group fields.
    for group in groups:
        field_name = get_readable_field_name(group)
        freq = ht[f"{field_name}_freq"]

        output_fields.update(
            {
                f"{field_name}_AC": freq.AC,
                f"{field_name}_AN": freq.AN,
                f"{field_name}_AF": freq.AF,
                f"{field_name}_homozygote_count": freq.homozygote_count,
                f"{field_name}_hom_depletion": ht[f"{field_name}_hom_depletion"],
                f"{field_name}_af_cutoff": ht[f"{field_name}_af_cutoff"],
                f"{field_name}_expected_hom": ht[f"{field_name}_expected_hom"],
            }
        )

    # Sum of expected homozygotes across all genetic ancestry groups
    # (excluding adj/overall).
    ancestry_expected_hom_fields = [
        ht[f"{get_readable_field_name(group)}_expected_hom"]
        for group in groups
        if group != "adj"
    ]
    output_fields["group_total_expected_hom"] = (
        hl.sum(hl.array(ancestry_expected_hom_fields))
        if ancestry_expected_hom_fields
        else hl.missing(hl.tfloat64)
    )

    # Raw is all-samples only, so these pair with the overall flag.
    # `hom_depleted_adj_only` marks overall flags whose raw homozygotes all failed
    # adj filtering: worth a closer look, not necessarily an artifact.
    raw_freq = ht.raw_freq
    zero_hom_raw = raw_freq.homozygote_count == 0
    output_fields.update(
        {
            "raw_AC": raw_freq.AC,
            "raw_AN": raw_freq.AN,
            "raw_AF": raw_freq.AF,
            "raw_homozygote_count": raw_freq.homozygote_count,
            "zero_hom_raw": zero_hom_raw,
        }
    )
    if "adj" in groups:
        output_fields["hom_depleted_adj_only"] = (
            ht.overall_hom_depletion & ~zero_hom_raw
        )

    # Per variant fields.
    output_fields.update(
        {
            "most_severe_csq": ht.most_severe_csq,
            "most_severe_gene": ht.most_severe_gene,
            "most_severe_gene_id": ht.most_severe_gene_id,
            "is_omim_recessive_gene": ht.is_omim_recessive_gene,
            "omim": ht.omim,
            "gerp": ht.gerp,
            "lcr_or_segdup": ht.lcr_or_segdup,
            "phylop_score": ht.phylop_score,
            "constraint": ht.constraint,
            "clinvar": ht.clinvar,
            "impc": ht.impc,
            # Gene-level lookups against published association studies, grouped by
            # what the evidence speaks to: whether biallelic humans have been
            # observed, and whether the gene shows a phenotype association. One
            # child struct per study within each group.
            "association_studies": hl.struct(
                biallelic_carriers=hl.struct(
                    lassen2026_ko=ht.lassen2026_ko,
                    brava_carriers=ht.brava_carriers,
                ),
                phenotype_assoc=hl.struct(
                    heyne_results=ht.heyne_results,
                    lassen_results=ht.lassen2026_recessive_assoc,
                ),
            ),
            "sv_cnv": ht.sv_cnv,
        }
    )

    return ht.select(**output_fields)


def main(args: argparse.Namespace) -> None:
    """Perform checks for homozygous depletions."""
    if args.test_n_partitions is not None and not args.test:
        raise ValueError("--test-n-partitions requires --test.")

    hl.init(
        log="/hom_depletion_check.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")

    try:
        out = get_hom_depletion(data_type=args.data_type, test=args.test)
        check_resource_existence(
            output_step_resources={"check_hom_depletion.py": [out]},
            overwrite=args.overwrite,
        )

        logger.info("Running hom depletion check with all groups...")
        # Obtain list of genetic ancestry groups, `remaining` is excluded: it's a
        # grab-bag of samples that didn't assign to any group, so an HWE expectation
        # computed within it isn't meaningful.
        gen_anc_groups = [
            g for g in GEN_ANC_GROUPS["v4"][args.data_type] if g != "remaining"
        ]
        groups = ["adj"] + [f"{g}_adj" for g in gen_anc_groups]

        ht = public_release(args.data_type).ht()
        if args.test_n_partitions is not None:
            logger.info(
                "Filtering to the first %d partitions for testing...",
                args.test_n_partitions,
            )
            ht = ht._filter_partitions(range(args.test_n_partitions))

        logger.info(
            "Filter to PASS variants in autosomes and annotate "
            "lcr/segdup/most_severe_consequence..."
        )
        ht = ht.filter(hl.len(ht.filters) == 0)
        ht = filter_to_autosomes(ht)
        ht = ht.annotate(
            lcr_or_segdup=(ht.region_flags.lcr | ht.region_flags.segdup),
            most_severe_csq=ht.vep.most_severe_consequence,
            phylop_score=ht.in_silico_predictors.phylop,
        )

        # Keep only fields needed for downstream depletion annotations.
        ht = ht.select(
            "freq",
            "most_severe_csq",
            "vep",
            "lcr_or_segdup",
            "phylop_score",
        )

        logger.info("Annotate group-specific frequencies and AF cutoffs...")
        ht = annotate_group_freqs_and_cutoffs(ht, groups, args.min_expected)

        logger.info("Flag variants with depleted homozygous counts...")
        ht = flag_hom_depletion(ht, groups)

        logger.info("Pull out first gene with most severe consequence...")
        ht = annotate_most_severe_gene(ht)

        # vep and freq are not used downstream and are large; drop them before the
        # gene annotations to cut serialization volume.
        ht = ht.drop("vep", "freq")

        logger.info("Annotating OMIM recessive-gene label...")
        ht = annotate_omim_recessive_gene_label(ht, OMIM_PATH)

        logger.info("Annotating GERP scores...")
        ht = annotate_gerp(ht)

        logger.info("Annotating gene constraint metrics...")
        ht = annotate_constraint(ht)

        logger.info("Annotating ClinVar clinical significance and phenotype...")
        ht = annotate_clinvar(ht)

        logger.info("Annotating IMPC mouse-knockout viability...")
        ht = annotate_impc_viability(ht, IMPC_PATH, MGI_PATH)

        logger.info("Annotating Lassen et al. 2026 biallelic pLoF knockout evidence...")
        ht = annotate_lassen_biallelic_ko(ht, LASSEN_S4_PATH)

        logger.info("Annotating Lassen et al. 2026 recessive associations...")
        ht = annotate_lassen_recessive_assoc(ht, LASSEN_S5_PATH)

        logger.info("Annotating BRaVa per-cohort biallelic carriers...")
        ht = annotate_brava_cohort_carriers(ht, BRAVA_PRE_META_PATH)

        logger.info("Annotating FinnGen coding associations (Heyne et al. 2023)...")
        ht = annotate_heyne_assoc(ht, HEYNE_S1_PATH)

        # Break lineage before the SV interval join. Only a select sits between SV
        # and the final write, so checkpointing here isolates the long upstream
        # annotation DAG (freq + OMIM/GERP/constraint/ClinVar/IMPC joins) from the
        # SV shuffle that was failing under shuffle-fetch loss.
        logger.info("Checkpointing annotations before SV/CNV overlap join...")
        ht = ht.checkpoint(new_temp_file("hom_depletion_pre_sv", "ht"))

        logger.info("Annotating raw SV/CNV overlap evidence...")
        sv_ht = import_sv_sites_ht()
        ht = annotate_sv_cnv_overlap(ht, sv_ht)

        logger.info("Select and flatten output fields...")
        ht = select_group_output_fields(ht, groups)

        logger.info("Checkpointing annotated table before final write...")
        ht = ht.checkpoint(new_temp_file("hom_depletion_pre_write", "ht"))

        logger.info("Writing hom depletion check results to %s...", out.path)
        ht = ht.naive_coalesce(50)
        ht.write(out.path, overwrite=args.overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("hom_depletion_check"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--data-type",
        help="Dataset to run on, one of 'exomes' or 'genomes'",
        default="exomes",
        choices=["exomes", "genomes"],
        type=str,
    )
    parser.add_argument(
        "--min-expected",
        help=(
            "Minimum expected number of homozygous individuals used to compute the"
            " AF cutoff. Default is 3, following Arnadottir et al. 2022"
            " (https://doi.org/10.1038/s41467-022-28330-8)."
        ),
        default=3,
        type=int,
    )
    parser.add_argument(
        "--test",
        help="Whether to write results to a tmp path instead of the analyses bucket.",
        action="store_true",
    )
    parser.add_argument(
        "--test-n-partitions",
        help=(
            "Number of partitions of the release HT to run on. Requires --test. "
            "Default is None, which runs on the full dataset."
        ),
        type=int,
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite an existing hom depletion output table.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
