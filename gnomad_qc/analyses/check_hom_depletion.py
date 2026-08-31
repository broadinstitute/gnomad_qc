import hail as hl
from hail.utils import new_temp_file
from gnomad.resources.grch38.gnomad import constraint, public_release
from gnomad.resources.grch38.reference_data import clinvar
from gnomad.utils.filtering import filter_to_autosomes

from gnomad_qc.v5.resources.basics import get_logging_path

from typing import List
import argparse
import logging


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hom_depletion_check")
logger.setLevel(logging.INFO)


SV_SITES_VCF_PATH = (
    "gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz"
)
OMIM_PATH = "gs://gnomad-kristen/hom_depletion/omim_2023.tsv"
IMPC_PATH = "gs://gnomad-kristen/hom_depletion/viability.csv.gz"  # 2026-03-16 13:49
MGI_PATH = "gs://gnomad-kristen/hom_depletion/HOM_AllOrganism.txt"  # accessed 07/07/2026 https://www.informatics.jax.org/homology.shtml
# Table S4 of Lassen et al. 2026 (AJHG 113:1330-1346, doi:10.1016/j.ajhg.2026.04.005):
# genes with >=1 observed biallelic pLoF genotype ("knockout") in their six-biobank
# meta-analysis (948,690 individuals), compared against previously published
# knockout surveys.
LASSEN_S4_PATH = (
    "gs://gnomad-kristen/hom_depletion/lassen2026_meta_analysis_supp_s4.txt"
)
# Column order of LASSEN_S4_PATH, which is read positionally (see
# `annotate_lassen_biallelic_ko`). "Novel" is their flag for genes whose knockout
# was first observed in the 2026 meta-analysis. gnomAD/Sulem/Narasimhan/Saleheen
# counts are secondhand, taken from the Oddsson et al. 2023 compilation.
# Table S5 of the same paper: the 58 significant gene-trait recessive associations
# (39 genes), each classified as Recessive / Additive / Ambiguous. Export the S5
# sheet to this path as tab-delimited text, as for S4.
# Supplementary Table 10 of Heyne et al. 2023 (Nature 613:519-525,
# doi:10.1038/s41586-022-05420-7): per-variant FinnGen release 4 summary
# statistics, including homozygote counts, for coding variants with any additive
# association P < 1e-4. Upload the unzipped r4_coding_variants_web_table_v0.txt
# to this path.
HEYNE_S10_PATH = (
    "gs://gnomad-kristen/hom_depletion/heyne2023_finngen_r4_coding_variants.txt"
)
# Number of FinnGen individuals analysed in Heyne et al. (release 4), used as the
# denominator when computing how many homozygotes FinnGen was expected to see.
FINNGEN_R4_N = 176899
LASSEN_S5_PATH = (
    "gs://gnomad-kristen/hom_depletion/lassen2026_meta_analysis_supp_s5.txt"
)
LASSEN_S4_DATASETS = [
    "lassen2026",
    "sulem",
    "narasimhan",
    "saleheen",
    "gnomad_via_oddsson",
    "oddsson",
    "sun",
]


def get_field_name(group: str) -> str:
    """
    Get the field name for a group (for clarity, "adj" -> "overall").

    :param group: Group key (e.g., "adj", "afr_adj")
    :return: Field name (e.g., "overall", "afr_adj")
    """
    return "overall" if group == "adj" else group.replace("_adj", "")


def annotate_group_freqs_and_cutoffs(
    ht: hl.Table,
    groups: List[str],
    min_expected: int = 3,
):
    """
    Annotate a Hail Table with group-specific allele frequencies and AF cutoffs.

    For each group, this function:
      - Retrieves the corresponding frequency index from `ht.freq_index_dict`
      - Annotates `<field_name>_freq` with the group's frequency
      - Annotates `<field_name>_af_cutoff` using:
            sqrt(min_expected / (AN / 2))

    # If a group's index is not defined, the annotations are set to missing.

    Also annotates `raw_freq` from the "raw" index. The group frequencies above are
    adj-filtered, so their homozygote counts exclude genotypes failing genotype QC
    (GQ/DP/allele balance). Carrying raw alongside them lets a zero adj homozygote
    count be separated into "no homozygous genotypes were called" versus "homozygous
    genotypes were called and then filtered", which are different findings.

    Note that `freq` carries a single all-samples raw group, not one per genetic
    ancestry group, so `raw_freq` is only comparable to the overall (adj) group.

    :param ht: Hail Table containing `freq` (array of frequency structs) and`freq_index_dict` (dict mapping group names to indices in `freq`).
    :param groups: List of group keys to use (e.g., ["adj", "afr_adj"]).
    :param min_expected: Minimum expected allele count used to compute AF cutoff
    :return: Annotated Hail Table with added frequency and cutoff fields.
    """
    group_freqs = {}
    cutoffs = {}
    metrics = {}

    for group in groups:
        idx = ht.freq_index_dict.get(group)
        field_name = get_field_name(group)

        defined_idx = hl.is_defined(idx)

        freq_expr = hl.or_missing(
            defined_idx,
            ht.freq[idx],
        )

        group_freqs[f"{field_name}_freq"] = freq_expr

        cutoffs[f"{field_name}_af_cutoff"] = hl.or_missing(
            defined_idx,
            hl.sqrt(min_expected / (freq_expr.AN / 2)),
        )

        metrics[f"{field_name}_expected_hom"] = hl.or_missing(
            defined_idx,
            freq_expr.AF**2 * (freq_expr.AN / 2),
        )

    raw_idx = ht.freq_index_dict.get("raw")
    group_freqs["raw_freq"] = hl.or_missing(hl.is_defined(raw_idx), ht.freq[raw_idx])

    return ht.annotate(**group_freqs, **cutoffs, **metrics)


def flag_hom_depletion(ht, groups):
    """
    Flag variants with depleted homozygous counts based on a specified metric and AF cutoff.

    A variant is flagged as hom-depleted if its frequency is greater than the AF cutoff and it has zero homozygous counts in the corresponding group.

    :param ht: Hail Table containing group-specific frequency and cutoff annotations.
    :param groups: List of group keys to use (e.g., ["adj", "afr_adj"]).
    :return: Hail Table with added boolean fields indicating hom-depletion for each group (e.g., "overall_flagged", "afr_adj_flagged").
    """
    # Dynamic threshold: AF must be greater than the square root of (desired_homs / N).
    flagged_groups = {}
    for group in groups:
        field_name = get_field_name(group)
        freq_field = f"{field_name}_freq"
        flagged_groups[f"{field_name}_hom_depletion"] = (
            hl.is_defined(ht[freq_field])
            & (ht[freq_field].AF > ht[f"{field_name}_af_cutoff"])
            & (ht[freq_field].homozygote_count == 0)
        )

    return ht.annotate(**flagged_groups)


def annotate_most_severe_gene(ht: hl.Table) -> hl.Table:
    most_severe_tx = hl.find(
        lambda tx: tx.consequence_terms.contains(ht.most_severe_csq),
        ht.vep.transcript_consequences,
    )

    return ht.annotate(
        most_severe_gene=most_severe_tx.gene_symbol,
        most_severe_gene_id=most_severe_tx.gene_id,
    )


def annotate_omim_recessive_gene_label(
    ht: hl.Table,
    omim_path: str,
) -> hl.Table:
    """
    Annotate whether `most_severe_gene` is OMIM-labeled as recessive.

    :param ht: Hail Table containing `most_severe_gene`.
    :param omim_path: Path to OMIM table (TSV/CSV) with gene + inheritance columns.
    :return: Hail Table with `is_omim_recessive_gene` boolean annotation.
    """
    omim_ht = hl.import_table(omim_path, impute=True)

    omim_ht = omim_ht.filter(
        hl.or_else(hl.str(omim_ht["phenotype_inheritance"]), "")
        .lower()
        .contains("recessive")
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

    o = omim_ht[ht.most_severe_gene_id]
    return ht.annotate(
        omim=o.omim,
        is_omim_recessive_gene=hl.is_defined(o),
    )


def annotate_gerp(ht: hl.Table) -> hl.Table:
    """Annotate GERP scores by locus."""
    gerp_ht = hl.experimental.load_dataset(
        name="gerp_scores", version="hg19", reference_genome="GRCh38"
    )

    #### -12.3 to 6.17 is the range of GERP values where 6.17 is the most conserved.
    return ht.annotate(gerp=hl.or_else(gerp_ht[ht.locus].S, 0))


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

    c = constraint_ht[ht.most_severe_gene_id]

    def _metrics(x: hl.expr.StructExpression) -> hl.expr.StructExpression:
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

    info = clinvar_ht[ht.key].info
    return ht.annotate(
        clinvar=hl.struct(
            disease=info.CLNDN,  # phenotype / disease name(s)
            clinical_significance=info.CLNSIG,  # benign, pathogenic, VUS, etc.
            review_status=info.CLNREVSTAT,
        )
    )


def get_mouse_to_human_ortholog_ht(mgi_path: str) -> hl.Table:
    """
    Build a mouse->human ortholog map from MGI's HomoloGene HOM_AllOrganism report.

    Genes sharing a `DB Class Key` are homologs across organisms. Grouping by that key
    and separating human (taxon 9606) and mouse (taxon 10090) members yields, per class,
    the human and mouse symbol sets. Result is keyed by mouse `Symbol` and records the
    human ortholog(s) plus an `ortholog_type` flag ("one2one" vs "one2many") so non-1:1
    mappings can be handled downstream.

    :param mgi_path: Path to MGI HOM_AllOrganism.txt (tab-delimited).
    :return: Table keyed by mouse `Symbol` with `human_symbols` (set) and `ortholog_type`.
    """
    homer = hl.import_table(mgi_path, impute=True)

    by_class = homer.group_by(class_key=homer["DB Class Key"]).aggregate(
        human_symbols=hl.agg.filter(
            homer["NCBI Taxon ID"] == 9606, hl.agg.collect_as_set(homer["Symbol"])
        ),
        mouse_symbols=hl.agg.filter(
            homer["NCBI Taxon ID"] == 10090, hl.agg.collect_as_set(homer["Symbol"])
        ),
    )

    # Keep classes with both a human and a mouse member.
    by_class = by_class.filter(
        (hl.len(by_class.human_symbols) > 0) & (hl.len(by_class.mouse_symbols) > 0)
    )
    by_class = by_class.annotate(
        ortholog_type=hl.if_else(
            (hl.len(by_class.human_symbols) == 1)
            & (hl.len(by_class.mouse_symbols) == 1),
            "one2one",
            "one2many",
        )
    )

    # One row per mouse symbol, keyed for joining to the IMPC report.
    ortholog_ht = by_class.annotate(
        mouse_symbol=hl.array(by_class.mouse_symbols)
    ).explode("mouse_symbol")

    return ortholog_ht.key_by("mouse_symbol").select("human_symbols", "ortholog_type")


def annotate_impc_viability(ht: hl.Table, impc_path: str, mgi_path: str) -> hl.Table:
    """
    Annotate mouse-knockout homozygous viability (IMPC) on `most_severe_gene`.

    Viability calls (viable / subviable / lethal) are mapped from mouse to human genes
    via MGI HomoloGene orthologs. `ortholog_type` records whether the mapping was 1:1;
    non-one2one genes (duplicated families) should be interpreted with care.

    :param ht: Hail Table containing `most_severe_gene` (human gene symbol).
    :param impc_path: Path to the IMPC viability.csv.gz report.
    :param mgi_path: Path to MGI HOM_AllOrganism.txt ortholog report.
    :return: Hail Table with an `impc` struct annotation.
    """
    impc_ht = hl.import_table(
        impc_path, impute=True, force=True, delimiter=",", quote='"'
    )
    ortholog_ht = get_mouse_to_human_ortholog_ht(mgi_path)

    ortholog = ortholog_ht[impc_ht["Gene Symbol"]]
    impc_ht = impc_ht.annotate(
        human_symbols=ortholog.human_symbols,
        ortholog_type=ortholog.ortholog_type,
        viability_call=impc_ht["Viability Phenotype HOMs/HEMIs"],
    )

    # Keep genes with a human ortholog, then explode to one row per human symbol.
    impc_ht = impc_ht.filter(hl.is_defined(impc_ht.human_symbols))
    impc_ht = impc_ht.annotate(human_symbol=hl.array(impc_ht.human_symbols)).explode(
        "human_symbol"
    )

    impc_ht = impc_ht.group_by(human_gene=impc_ht.human_symbol).aggregate(
        viability=hl.agg.collect_as_set(impc_ht.viability_call),
        ortholog_type=hl.agg.collect_as_set(impc_ht.ortholog_type),
    )

    impc = impc_ht[ht.most_severe_gene]
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
    Annotate whether a human biallelic pLoF knockout of `most_severe_gene` has been observed.

    Uses table S4 of Lassen et al. 2026, which reports genes carrying at least one
    biallelic pLoF genotype in their six-biobank meta-analysis (948,690 individuals)
    alongside the same tally from earlier knockout surveys (Sun et al. 2024 directly;
    Sulem, Narasimhan, Saleheen, gnomAD, and Oddsson via the Oddsson et al. 2023
    compilation). A gene that is depleted of homozygotes here but has an observed
    knockout there is unlikely to be intolerant of biallelic loss, so the depletion is
    more plausibly allele-specific or technical. Conversely, `ko_novel_to_lassen2026`
    marks genes whose first observed knockout is this 2026 meta-analysis -- evidence
    that post-dates gnomAD's own knockout tallies.

    Genes absent from S4 get `biallelic_plof_ko_observed=False` with
    `gene_in_lassen_s4=False`; absence conflates true intolerance with lack of power
    to observe a knockout, so pair it with constraint and mouse-viability annotations.

    S4 is restricted to pLoF genotypes, so it speaks to whether the gene tolerates
    biallelic loss-of-function, not whether a specific missense or in-frame allele is
    tolerated.

    Matching is by gene symbol. Symbols that VEP and S4 spell differently (aliases)
    will not match and are reported as absent, which biases toward "no knockout
    observed"; treat a lone absence as weaker evidence than a positive hit.

    :param ht: Hail Table containing `most_severe_gene`.
    :param s4_path: Path to the tab-delimited export of Lassen et al. 2026 table S4.
    :return: Hail Table with a `lassen2026_ko` struct annotation.
    """
    # Read positionally: the header carries empty trailing column names and a
    # trailing space in "gnomAD ", both of which break named import.
    s4 = hl.import_table(s4_path, no_header=True, delimiter="\t")
    s4 = s4.filter(~hl.literal({"Gene", "TOTAL"}).contains(s4.f0))

    def _count(field: hl.expr.StringExpression) -> hl.expr.Int32Expression:
        return hl.or_else(hl.parse_int32(field), 0)

    s4 = s4.select(
        gene_symbol=s4.f0,
        lassen2026=_count(s4.f1),
        novel=_count(s4.f2),
        sulem=_count(s4.f3),
        narasimhan=_count(s4.f4),
        saleheen=_count(s4.f5),
        gnomad_via_oddsson=_count(s4.f6),
        oddsson=_count(s4.f7),
        sun=_count(s4.f8),
    )

    # Join on gene symbol: it reaches 9,092 of the 9,176 S4 rows, whereas the
    # Ensembl ID column is "NA" for ~36% of them. Drop rows whose symbol is unset
    # ("NA"), the sole repeated symbol key, which would otherwise match any variant
    # lacking a gene symbol.
    value_fields = LASSEN_S4_DATASETS + ["novel"]
    s4 = (
        s4.filter(hl.is_defined(s4.gene_symbol) & (s4.gene_symbol != "NA"))
        .key_by("gene_symbol")
        .select(*value_fields)
    )

    rec = s4[ht.most_severe_gene]

    observing = hl.array(
        [
            hl.or_missing(rec[dataset] > 0, dataset)
            for dataset in LASSEN_S4_DATASETS
        ]
    ).filter(hl.is_defined)

    return ht.annotate(
        lassen2026_ko=hl.struct(
            gene_in_lassen_s4=hl.is_defined(rec),
            biallelic_plof_ko_observed=hl.or_else(hl.len(observing) > 0, False),
            n_datasets_observing_ko=hl.or_else(hl.len(observing), 0),
            datasets_observing_ko=observing,
            ko_observed_in_gnomad_per_oddsson=hl.or_else(
                rec.gnomad_via_oddsson > 0, False
            ),
            ko_novel_to_lassen2026=hl.or_else(rec.novel > 0, False),
        )
    )


def annotate_lassen_recessive_assoc(
    ht: hl.Table,
    s5_path: str = LASSEN_S5_PATH,
) -> hl.Table:
    """
    Annotate published recessive gene-trait associations for `most_severe_gene`.

    Uses table S5 of Lassen et al. 2026: the significant associations from their
    six-biobank recessive analysis, each labeled Recessive, Additive, or Ambiguous
    according to whether the recessive model fit better than the additive one.

    An association is evidence that biallelic carriers of the gene both exist and
    reach adulthood with a measurable phenotype, so -- like the S4 knockout tally --
    it argues against absolute intolerance of biallelic loss while marking the gene
    as disease relevant. The tables are small (58 associations over 39 genes), so
    most genes get `has_association=False`; this is a lookup for the occasional hit,
    not a broadly informative annotation.

    :param ht: Hail Table containing `most_severe_gene`.
    :param s5_path: Path to the tab-delimited export of Lassen et al. 2026 table S5.
    :return: Hail Table with a `lassen2026_recessive_assoc` struct annotation.
    """
    # Positional read, as for S4: the exported header carries empty trailing column
    # names and a trailing space in "P rec pLoF ".
    s5 = hl.import_table(s5_path, no_header=True, delimiter="\t")
    s5 = s5.filter(s5.f0 != "Phenotype")
    s5 = s5.select(
        gene_symbol=s5.f1,
        phenotype=s5.f0,
        mode=s5.f2,
        p_rec=hl.parse_float64(s5.f3),
        p_add=hl.parse_float64(s5.f4),
    )
    s5 = s5.group_by(gene_symbol=s5.gene_symbol).aggregate(
        associations=hl.agg.collect(
            hl.struct(
                phenotype=s5.phenotype,
                mode=s5.mode,
                p_rec=s5.p_rec,
                p_add=s5.p_add,
            )
        )
    )

    rec = s5[ht.most_severe_gene]
    associations = hl.or_else(
        rec.associations, hl.empty_array(rec.associations.dtype.element_type)
    )

    return ht.annotate(
        lassen2026_recessive_assoc=hl.struct(
            has_association=hl.len(associations) > 0,
            n_associations=hl.len(associations),
            # True only when the recessive model fit better than the additive one
            # for at least one trait, i.e. the paper's "more likely recessive" set.
            has_recessive_mode=associations.any(lambda a: a.mode == "Recessive"),
            modes=hl.set(associations.map(lambda a: a.mode)),
            phenotypes=associations.map(lambda a: a.phenotype),
            associations=associations,
        )
    )


def annotate_finngen_homozygotes(
    ht: hl.Table,
    s10_path: str = HEYNE_S10_PATH,
    n_finngen: int = FINNGEN_R4_N,
    min_info: float = 0.9,
    min_expected_hom: float = 5.0,
) -> hl.Table:
    """
    Annotate FinnGen homozygote evidence for each variant from Heyne et al. 2023.

    FinnGen is an independent cohort on an independent technology (imputed array
    data rather than sequencing), so a homozygote deficit seen in both datasets is
    unlikely to share a mapping or paralog artifact. Founder-population enrichment
    also means FinnGen can be better powered than its sample size suggests for
    Finnish-enriched alleles, and worse for alleles that are rarer in Finns.

    A zero homozygote count is only meaningful relative to what FinnGen was powered
    to observe, so this computes the expected count under Hardy-Weinberg,
    `expected_hom = AF^2 * n_finngen`, and the Poisson probability of seeing none,
    `p_zero_hom = exp(-expected_hom)`. `verdict` collapses this into:

      - `homozygotes_observed`   FinnGen has homozygotes; the depletion is refuted.
      - `depletion_replicated`   well imputed, `expected_hom >= min_expected_hom`,
                                 and none observed.
      - `underpowered`           none observed but too few expected to conclude
                                 anything (the common case, easily misread as
                                 corroboration).
      - `low_imputation_quality` INFO below `min_info`; counts are unreliable.
      - `not_in_finngen`         variant absent from the table, which reflects
                                 FinnGen imputation coverage (and the table's own
                                 additive P < 1e-4 filter), not biology.

    :param ht: Hail Table keyed by locus and alleles (GRCh38).
    :param s10_path: Path to the Heyne et al. table S10 text file.
    :param n_finngen: Number of FinnGen individuals behind the homozygote counts.
    :param min_info: Minimum imputation INFO score to trust a count.
    :param min_expected_hom: Expected homozygotes needed to call a zero count a
        replication rather than a lack of power.
    :return: Hail Table with a `finngen` struct annotation.
    """
    fg = hl.import_table(s10_path, delimiter="\t")

    # `variant` is GRCh38 "chrom:pos:ref:alt" without the "chr" prefix.
    fg = fg.annotate(_v=fg.variant.split(":"))
    fg = fg.annotate(
        locus=hl.locus(
            "chr" + fg._v[0], hl.parse_int32(fg._v[1]), reference_genome="GRCh38"
        ),
        alleles=[fg._v[2], fg._v[3]],
        af=hl.parse_float64(fg.AF),
        ac_hom=hl.parse_int32(fg.AC_Hom),
        info_score=hl.parse_float64(fg.INFO),
        enrichment_nfsee=hl.parse_float64(fg.enrichment_nfsee),
        pval_recessive=hl.parse_float64(fg.pval_recessive),
    )

    # The table carries one row per variant-phenotype pair; collapse to one row
    # per variant, keeping the per-variant fields and the best recessive p-value.
    fg = fg.group_by("locus", "alleles").aggregate(
        af=hl.agg.take(fg.af, 1)[0],
        ac_hom=hl.agg.take(fg.ac_hom, 1)[0],
        info_score=hl.agg.take(fg.info_score, 1)[0],
        enrichment_nfsee=hl.agg.take(fg.enrichment_nfsee, 1)[0],
        min_pval_recessive=hl.agg.min(fg.pval_recessive),
        n_phenotypes=hl.agg.count(),
    )

    f = fg[ht.key]
    expected_hom = f.af**2 * n_finngen

    return ht.annotate(
        finngen=hl.struct(
            in_finngen=hl.is_defined(f.af),
            af=f.af,
            ac_hom=f.ac_hom,
            info_score=f.info_score,
            enrichment_nfsee=f.enrichment_nfsee,
            expected_hom=expected_hom,
            # Poisson P(0 | expected); only meaningful when no homozygote was seen.
            p_zero_hom=hl.or_missing(f.ac_hom == 0, hl.exp(-expected_hom)),
            min_pval_recessive=f.min_pval_recessive,
            verdict=(
                hl.case()
                .when(hl.is_missing(f.af), "not_in_finngen")
                .when(f.ac_hom > 0, "homozygotes_observed")
                .when(f.info_score < min_info, "low_imputation_quality")
                .when(expected_hom >= min_expected_hom, "depletion_replicated")
                .default("underpowered")
            ),
        )
    )


def import_sv_sites_ht(sv_vcf_path: str = SV_SITES_VCF_PATH) -> hl.Table:
    """
    Import the gnomAD SV sites VCF as a rows-only Hail Table.

    :param sv_vcf_path: Path to the gnomAD SV sites VCF (bgzipped).
    :return: Hail Table of SV sites (one row per SV), for use with `annotate_sv_cnv_overlap`.
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
    :param use_grpmax: DUP frequency uses GRPMAX_AF (fallback AF[0]) if True, else AF[0].
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

    # Materialize the slim interval table to native Hail format before the join:
    # the VCF-derived rows are decoded once here rather than re-decoded inside the
    # interval join (which was throwing a decode NullPointerException).
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


def select_group_output_fields(ht, groups):
    """
    Select and flatten group-specific annotations.

    :param ht: Hail Table containing:
        - `<field_name>_freq` structs with AC, AN, AF, homozygote_count
        - `<field_name>_hom_depletion` boolean fields
        - `<field_name>_af_cutoff` fields
        - `raw_freq` struct with AC, AN, AF, homozygote_count
        - `most_severe_csq`, `most_severe_gene`, `lcr_or_segdup`
        - `lassen2026_ko` struct of observed-human-knockout evidence
    :param groups: List of group keys (e.g., ["adj", "afr_adj"]).
    :return: Hail Table with selected and flattened fields.
    """
    output_fields = {}

    # per-group fields.
    for group in groups:
        field_name = get_field_name(group)
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

    # Sum of expected homozygotes across all genetic ancestry groups (excluding adj/overall).
    ancestry_expected_hom_fields = [
        ht[f"{get_field_name(group)}_expected_hom"]
        for group in groups
        if group != "adj"
    ]
    output_fields["group_total_expected_hom"] = (
        hl.sum(hl.array(ancestry_expected_hom_fields))
        if ancestry_expected_hom_fields
        else hl.missing(hl.tfloat64)
    )

    # Raw (pre-genotype-QC) homozygote evidence. `hom_depleted_adj_only` marks
    # variants whose homozygotes were called but removed by adj filtering, so the
    # apparent depletion is a genotype-quality artifact rather than a real absence.
    #
    # Raw exists only as a single all-samples group, so these columns pair with the
    # overall flag only. They cannot adjudicate a per-ancestry flag: raw homozygotes
    # in one ancestry group would make `zero_hom_raw` False for every group.
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
            "lassen2026_ko": ht.lassen2026_ko,
            "lassen2026_recessive_assoc": ht.lassen2026_recessive_assoc,
            "finngen": ht.finngen,
            "sv_cnv": ht.sv_cnv,
        }
    )

    return ht.select(**output_fields)


def main(args):
    """Perform checks for homozygous depletions."""
    hl.init(
        log="/hom_depletion_check.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    hl.default_reference("GRCh38")

    try:

        if args.adj_only:
            logger.info("Running hom depletion check with adj-only groups...")
            groups = ["adj"]
        else:
            logger.info("Running hom depletion check with all groups...")
            groups = [
                "adj",
                "afr_adj",
                "amr_adj",
                "asj_adj",
                "eas_adj",
                "fin_adj",
                "nfe_adj",
                "sas_adj",
            ]

        ht = public_release(args.data_type).ht()

        logger.info(
            "Filter to PASS variants in autosomes and annotate lcr/segdup/most_severe_consequence..."
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

        # vep and freq are not used downstream and are large; drop before the
        # gene-keyed joins and checkpoint to cut shuffle/serialization volume.
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

        logger.info("Annotating FinnGen homozygote evidence (Heyne et al. 2023)...")
        ht = annotate_finngen_homozygotes(ht, HEYNE_S10_PATH)

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

        logger.info("Writing hom depletion check results...")
        ht = ht.naive_coalesce(50)
        ht.write("gs://gnomad-kristen/hom_depletion_af_check_all_2.ht", overwrite=True)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("hom_depletion_check"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument(
        "--data-type",
        help="Dataset to run on, one of 'exomes', 'genomes', or 'joint'",
        default="exomes",
        choices=["exomes", "genomes", "joint"],
        type=str,
    )
    parser.add_argument(
        "--min-expected",
        help="Minimum expected allele count used to compute AF cutoff.",
        default=3,
        type=int,
    )
    parser.add_argument(
        "--adj-only",
        help="Whether to run hom depletion check using only the adj group (overall).",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)
