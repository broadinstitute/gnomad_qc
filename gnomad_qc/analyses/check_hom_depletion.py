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


OMIM_PATH = "gs://gnomad-kristen/hom_depletion/omim_2023.tsv"
IMPC_PATH="gs://gnomad-kristen/hom_depletion/viability.csv.gz" 	#2026-03-16 13:49
MGI_PATH= "gs://gnomad-kristen/hom_depletion/HOM_AllOrganism.txt"  # accessed 07/07/2026 https://www.informatics.jax.org/homology.shtml	


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


def annotate_impc_viability(
    ht: hl.Table, impc_path: str, mgi_path: str
) -> hl.Table:
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
    impc_ht = impc_ht.annotate(
        human_symbol=hl.array(impc_ht.human_symbols)
    ).explode("human_symbol")

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


def select_group_output_fields(ht, groups):
    """
    Select and flatten group-specific annotations.

    :param ht: Hail Table containing:
        - `<field_name>_freq` structs with AC, AN, AF, homozygote_count
        - `<field_name>_hom_depletion` boolean fields
        - `<field_name>_af_cutoff` fields
        - `most_severe_csq`, `most_severe_gene`, `lcr_or_segdup`
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

        logger.info("Select and flatten output fields...")
        ht = select_group_output_fields(ht, groups)

        logger.info("Checkpointing annotated table before final write...")
        ht = ht.checkpoint(new_temp_file("hom_depletion_pre_write", "ht"))

        logger.info("Writing hom depletion check results...")
        ht = ht.naive_coalesce(50)
        ht.write("gs://gnomad-kristen/hom_depletion_af_check_all_1.ht", overwrite=True)

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
