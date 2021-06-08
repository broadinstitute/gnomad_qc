import argparse
import logging

import hail as hl

from gnomad.resources.grch38.gnomad import POPS_STORED_AS_SUBPOPS
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
    telomeres_and_centromeres,
)
from gnomad.sample_qc.relatedness import UNRELATED
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import get_adj_expr, region_flag_expr
from gnomad.utils.release import make_freq_index_dict
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    AS_FIELDS,
    SITE_FIELDS,
    SPARSE_ENTRIES,
)

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import (
    analyst_annotations,
    get_freq,
    get_info,
    vep,
)
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.release import (
    release_sites,
    hgdp_1kg_subset,
    hgdp_1kg_subset_annotations,
    hgdp_1kg_subset_sample_tsv,
)
from gnomad_qc.v3.resources.sample_qc import relatedness
from gnomad_qc.v3.resources.variant_qc import final_filter, SYNDIP
from gnomad_qc.v3.utils import hom_alt_depletion_fix

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_subset")
logger.setLevel(logging.INFO)

AS_FIELDS.remove("InbreedingCoeff")
SITE_FIELDS.remove("BaseQRankSum")

GLOBAL_SAMPLE_ANNOTATION_DICT = {
    "sex_imputation_ploidy_cutoffs": {
        "Description": (
            "Contains sex chromosome ploidy cutoffs used when determining sex chromosome karyotypes. Format: "
            "(upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for "
            "triple X) and (lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)."
        )
    },
    "population_inference_pca_metrics": {
        "Description": (
            "Contains the number of principal components (PCs) used when running PC-project and the minimum cutoff "
            "probability of belonging to a given population."
        )
    },
    "hard_filter_cutoffs": {
        "Description": (
            "Contains the cutoffs used for hard-filtering samples prior to sample QC. Sample QC metrics are "
            "computed using the Hail sample_qc module on all autosomal bi-allelic SNVs. Samples are removed if "
            "they are clear outliers for any of the following metrics: number of snps (n_snp), ratio of heterozygous "
            "variants to homozygous variants (r_het_hom_var), number of singletons (n_singleton), and mean coverage on "
            "chromosome 20 (cov). Additionally, we filter based on outliers of the following Picard metrics: % "
            "contamination (freemix), % chimera, and median insert size."
        )
    },
}
GLOBAL_VARIANT_ANNOTATION_DICT = {
    "cohort_freq_meta": {
        "Description": (
            "HGDP and 1KG frequency metadata. An ordered list containing the frequency aggregation group"
            "for each element of the cohort_freq array row annotation."
        )
    },
    "gnomad_freq_meta": {
        "Description": (
            "gnomAD frequency metadata. An ordered list containing the frequency aggregation group"
            "for each element of the gnomad_freq array row annotation."
        )
    },
    "cohort_freq_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: HGDP or 1KG "
            "subpopulation, sex: sex karyotype), with values describing the corresponding index of each grouping "
            "entry in the HGDP + 1KG frequency array annotation."
        )
    },
    "gnomad_freq_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD "
            "inferred global population sex: sex karyotype), with values describing the corresponding index "
            "of each grouping entry in the gnomAD frequency array annotation."
        )
    },
    "gnomad_faf_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD "
            "inferred global population sex: sex karyotype), with values describing the corresponding index "
            "of each grouping entry in the filtering allele frequency (using Poisson 99% CI) annotation."
        )
    },
    "gnomad_faf_meta": {
        "Description": (
            "gnomAD filtering allele frequency (using Poisson 99% CI) metadata. An ordered list "
            "containing the frequency aggregation group for each element of the gnomad_faf array row annotation."
        )
    },
    "vep_version": {"Description": "VEP version."},
    "vep_csq_header": {"Description": "VEP header for VCF export."},
    "dbsnp_version": {"Description": "dbSNP version."},
    "filtering_model": {
        "Description": "The variant filtering model used and its specific cutoffs.",
        "sub_globals": {
            "model_name": {
                "Description": (
                    "Variant filtering model name used in the 'filters'  row annotation to indicate"
                    "the variant was filtered by the model during variant QC."
                )
            },
            "score_name": {"Description": "Name of score used for variant filtering."},
            "snv_cutoff": {
                "Description": "SNV filtering cutoff information.",
                "sub_globals": {
                    "bin": {"Description": "Filtering percentile cutoff for SNVs."},
                    "min_score": {
                        "Description": "Minimum score (score_name) at SNV filtering percentile cutoff."
                    },
                },
            },
            "indel_cutoff": {
                "Description": "Information about cutoff used for indel filtering.",
                "sub_globals": {
                    "bin": {"Description": "Filtering percentile cutoff for indels."},
                    "min_score": {
                        "Description": "Minimum score (score_name) at indel filtering percentile cutoff."
                    },
                },
            },
            "snv_training_variables": {
                "Description": "Variant annotations used as features in SNV filtering model."
            },
            "indel_training_variables": {
                "Description": "Variant annotations used as features in indel filtering model."
            },
        },
    },
    "inbreeding_coeff_cutoff": {
        "Description": "Hard-filter cutoff for InbreedingCoeff on variants."
    },
}
GLOBAL_ANNOTATION_DICT = hl.struct(
    **GLOBAL_SAMPLE_ANNOTATION_DICT, **GLOBAL_VARIANT_ANNOTATION_DICT
)

SAMPLE_ANNOTATION_DICT = {
    "s": {"Description": "Sample ID."},
    "bam_metrics": {
        "Description": "Sample level metrics obtained from BAMs/CRAMs.",
        "sub_annotations": {
            "pct_bases_20x": {
                "Description": "The fraction of bases that attained at least 20X sequence coverage in post-filtering bases."
            },
            "pct_chimeras": {
                "Description": (
                    "The fraction of reads that map outside of a maximum insert size (usually 100kb) or that have the "
                    "two ends mapping to different chromosomes."
                )
            },
            "freemix": {"Description": "Estimate of contamination (0-100 scale)."},
            "mean_coverage": {
                "Description": "The mean coverage in bases of the genome territory, after all filters are applied, https://broadinstitute.github.io/picard/picard-metric-definitions.html."
            },
            "median_coverage": {
                "Description": "The median coverage in bases of the genome territory, after all filters are applied, https://broadinstitute.github.io/picard/picard-metric-definitions.html."
            },
            "mean_insert_size": {
                "Description": (
                    "The mean insert size of the 'core' of the distribution. Artefactual outliers in the distribution "
                    "often cause calculation of nonsensical mean and stdev values. To avoid this the distribution is "
                    "first trimmed to a 'core' distribution of +/- N median absolute deviations around the median "
                    "insert size."
                )
            },
            "median_insert_size": {
                "Description": "The median insert size of all paired end reads where both ends mapped to the same chromosome."
            },
            "pct_bases_10x": {
                "Description": "The fraction of bases that attained at least 10X sequence coverage in post-filtering bases."
            },
        },
    },
    "subsets": {
        "Description": "Struct containing information from the sample's cohort: HGDP or 1KG (tgp).",
        "sub_annotations": {
            "tgp": {"Description": "True if sample is from the 1KG dataset."},
            "hgdp": {"Description": "True if the sample is from the HGDP dataset."},
        },
    },
    "sex_imputation": {
        "Description": "Struct containing sex imputation information.",
        "sub_annotations": {
            "f_stat": {
                "Description": "Inbreeding coefficient (excess heterozygosity) on chromosome X."
            },
            "n_called": {"Description": "Number of variants with a genotype call."},
            "expected_homs": {"Description": "Expected number of homozygotes."},
            "observed_homs": {"Description": "Observed number of homozygotes."},
            "chr20_mean_dp": {
                "Description": "Sample's mean depth across chromosome 20."
            },
            "chrX_mean_dp": {"Description": "Sample's mean depth across chromosome X."},
            "chrY_mean_dp": {"Description": "Sample's mean depth across chromosome Y."},
            "chrX_ploidy": {
                "Description": "Sample's chromosome X ploidy (chrX_mean_dp normalized using chr20_mean_dp)."
            },
            "chrY_ploidy": {
                "Description": "Sample's chromosome Y ploidy (chrY_mean_dp normalized using chr20_mean_dp)."
            },
            "X_karyotype": {"Description": "Sample's chromosome X karyotype."},
            "Y_karyotype": {"Description": "Sample's chromosome Y karyotype."},
            "sex_karyotype": {
                "Description": "Sample's sex karyotype (combined X and Y karyotype)."
            },
        },
    },
    "sample_qc": {
        "Description": "Struct containing sample QC metrics calculated using hl.sample_qc().",
        "sub_annotations": {
            "n_hom_ref": {"Description": "Number of homozygous reference calls."},
            "n_het": {"Description": "Number of heterozygous calls."},
            "n_hom_var": {"Description": "Number of homozygous alternate calls."},
            "n_non_ref": {"Description": "Sum of n_het and n_hom_var."},
            "n_snp": {"Description": "Number of SNP alternate alleles."},
            "n_insertion": {"Description": "Number of insertion alternate alleles."},
            "n_deletion": {"Description": "Number of deletion alternate alleles."},
            "n_transition": {
                "Description": "Number of transition (A-G, C-T) alternate alleles."
            },
            "n_transversion": {
                "Description": "Number of transversion alternate alleles."
            },
            "r_ti_tv": {"Description": "Transition/Transversion ratio."},
            "r_het_hom_var": {"Description": "Het/HomVar call ratio."},
            "r_insertion_deletion": {"Description": "Insertion/Deletion allele ratio."},
        },
    },
    "population_inference": {
        "Description": (
            "Struct containing ancestry information assigned by applying a principal component analysis (PCA) on "
            "gnomAD samples and using those PCs in a random forest classifier trained on known gnomAD ancestry labels."
        ),
        "sub_annotations": {
            "pca_scores": {
                "Description": "Sample's scores for each gnomAD population PC."
            },
            "pop": {"Description": "Sample's inferred gnomAD population label."},
            "prob_afr": {
                "Description": "Random forest probability that the sample is of African/African-American ancestry."
            },
            "prob_ami": {
                "Description": "Random forest probability that the sample is of Amish ancestry."
            },
            "prob_amr": {
                "Description": "Random forest probability that the sample is of Latino ancestry."
            },
            "prob_asj": {
                "Description": "Random forest probability that the sample is of Ashkenazi Jewish ancestry."
            },
            "prob_eas": {
                "Description": "Random forest probability that the sample is of East Asian ancestry."
            },
            "prob_fin": {
                "Description": "Random forest probability that the sample is of Finnish ancestry."
            },
            "prob_mid": {
                "Description": "Random forest probability that the sample is of Middle Eastern ancestry."
            },
            "prob_nfe": {
                "Description": "Random forest probability that the sample is of Non-Finnish European ancestry."
            },
            "prob_oth": {
                "Description": "Random forest probability that the sample is of Other ancestry."
            },
            "prob_sas": {
                "Description": "Random forest probability that the sample is of South Asian ancestry."
            },
        },
    },
    "labeled_subpop": {
        "Description": "The sample's population label supplied by HGDP or 1KG."
    },
    "gnomad_release": {
        "Description": (
            "Indicates whether the sample was included in the gnomAD release dataset. For the full gnomAD release, "
            "relatedness inference is performed on the full dataset, and release samples are chosen in a way that "
            "maximizes the number of samples retained while filtering the dataset to include only samples with less "
            "than second-degree relatedness. For the HGDP + 1KG subset, samples passing all other sample QC "
            "metrics are retained."
        )
    },
    "high_quality": {
        "Description": (
            "Indicates whether a sample has passed all sample QC metrics except for relatedness."
        )
    },
}

SAMPLE_QC_METRICS = [
    "n_deletion",
    "n_het",
    "n_hom_ref",
    "n_hom_var",
    "n_insertion",
    "n_non_ref",
    "n_snp",
    "n_transition",
    "n_transversion",
    "r_het_hom_var",
    "r_insertion_deletion",
    "r_ti_tv",
]


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Parse relatedness Table to get every relationship (except UNRELATED) per sample.

    Return Table keyed by sample with all sample relationships, kin, ibd0, ibd1, and ibd2 in a struct.

    :param relatedness_ht: Table with inferred relationship information output by pc_relate.
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationship information annotated as a struct.
    """
    relatedness_ht = relatedness_ht.filter(relatedness_ht.relationship != UNRELATED)
    relationship_struct = hl.struct(
        kin=relatedness_ht.kin,
        ibd0=relatedness_ht.ibd0,
        ibd1=relatedness_ht.ibd1,
        ibd2=relatedness_ht.ibd2,
        relationship=relatedness_ht.relationship,
    )

    relatedness_ht_i = relatedness_ht.group_by(s=relatedness_ht.i.s).aggregate(
        relationships=hl.agg.collect_as_set(
            hl.struct(s=relatedness_ht.j.s, **relationship_struct)
        )
    )

    relatedness_ht_j = relatedness_ht.group_by(s=relatedness_ht.j.s).aggregate(
        relationships=hl.agg.collect_as_set(
            hl.struct(s=relatedness_ht.i.s, **relationship_struct)
        )
    )

    relatedness_ht = relatedness_ht_i.union(relatedness_ht_j)

    return relatedness_ht


def prepare_sample_annotations() -> hl.Table:
    """
    Load meta HT and select row and global annotations for HGDP + TGP subset.

    .. note::

        Expects that `meta.ht()` and `relatedness.ht()` exist. Relatedness pair information will be subset to only
        samples within HGDP + TGP and stored as the `relatedness_inference` annotation of the returned HT.

    :return: Table containing sample metadata for the subset
    """

    logger.info(
        "Subsetting and modifying sample QC metadata to desired globals and annotations"
    )
    meta_ht = meta.ht()
    meta_ht = meta_ht.filter(
        meta_ht.subsets.hgdp | meta_ht.subsets.tgp | (meta_ht.s == SYNDIP)
    )
    meta_ht = meta_ht.select_globals(
        global_annotation_descriptions=hl.literal(GLOBAL_SAMPLE_ANNOTATION_DICT),
        sample_annotation_descriptions=hl.literal(SAMPLE_ANNOTATION_DICT),
        sex_imputation_ploidy_cutoffs=meta_ht.sex_imputation_ploidy_cutoffs,
        population_inference_pca_metrics=hl.struct(
            n_pcs=meta_ht.population_inference_pca_metrics.n_pcs,
            min_prob=meta_ht.population_inference_pca_metrics.min_prob,
        ),
        hard_filter_cutoffs=meta_ht.hard_filter_cutoffs,
    )

    relatedness_ht = relatedness.ht()
    subset_samples = meta_ht.s.collect()
    relatedness_ht = relatedness_ht.filter(
        subset_samples.contains(relatedness_ht.i.s)
        & subset_samples.contains(relatedness_ht.j.s)
    )

    relatedness_ht = get_relatedness_set_ht(relatedness_ht)
    meta_ht = meta_ht.select(
        bam_metrics=meta_ht.bam_metrics,
        subsets=meta_ht.subsets.select("tgp", "hgdp"),
        sex_imputation=meta_ht.sex_imputation.drop("is_female"),
        sample_qc=meta_ht.sample_qc.select(*SAMPLE_QC_METRICS),
        relatedness_inference=relatedness_ht[meta_ht.key],
        population_inference=meta_ht.population_inference.drop(
            "training_pop", "training_pop_all"
        ),
        labeled_subpop=meta_ht.project_meta.project_subpop,
        gnomad_release=meta_ht.release,
        high_quality=meta_ht.high_quality,
    )

    return meta_ht


# TODO: Might be good to generalize this because a similar function is used in creating the release sites HT.
def prepare_variant_annotations(ht: hl.Table, filter_lowqual: bool = True) -> hl.Table:
    """
    Load and join all Tables with variant annotations.

    :param ht: Input HT to add variant annotations to.
    :param filter_lowqual: If True, filter out lowqual variants using the info HT's AS_lowqual.
    :return: Table containing joined annotations.
    """
    logger.info("Loading annotation tables...")
    filters_ht = final_filter.ht()
    vep_ht = vep.ht()
    dbsnp_ht = dbsnp.ht().select("rsid")
    info_ht = get_info().ht()
    analyst_ht = analyst_annotations.ht()
    freq_ht = get_freq().ht()
    score_name = hl.eval(filters_ht.filtering_model.score_name)

    if filter_lowqual:
        logger.info("Filtering lowqual variants...")
        ht = ht.filter(info_ht[ht.key].AS_lowqual, keep=False)

    logger.info("Assembling 'info' field...")
    info_fields = SITE_FIELDS + AS_FIELDS
    missing_info_fields = set(info_fields).difference(info_ht.info.keys())
    select_info_fields = set(info_fields).intersection(info_ht.info.keys())
    logger.info(
        "The following fields are not found in the info HT: %s", missing_info_fields,
    )

    # NOTE: SOR and AS_SOR annotations are now added to the info HT by default with get_as_info_expr and
    # get_site_info_expr in gnomad_methods, but they were not for v3 or v3.1
    filters = filters_ht[info_ht.key]
    info_ht = info_ht.transmute(
        info=info_ht.info.select(
            *select_info_fields,
            AS_SOR=filters.AS_SOR,
            SOR=filters.SOR,
            singleton=filters.singleton,
            transmitted_singleton=filters.transmitted_singleton,
            omni=filters.omni,
            mills=filters.mills,
            monoallelic=filters.monoallelic,
            InbreedingCoeff=freq_ht[info_ht.key].InbreedingCoeff,
            **{f"{score_name}": filters[f"{score_name}"]},
        )
    )

    logger.info("Adding annotations...")
    filters_ht = filters_ht.annotate(
        allele_info=hl.struct(
            variant_type=filters_ht.variant_type,
            allele_type=filters_ht.allele_type,
            n_alt_alleles=filters_ht.n_alt_alleles,
            was_mixed=filters_ht.was_mixed,
        ),
    )

    keyed_filters = filters_ht[ht.key]
    info = info_ht[ht.key]
    ht = ht.annotate(
        a_index=info.a_index,
        was_split=info.was_split,
        rsid=dbsnp_ht[ht.key].rsid,
        filters=keyed_filters.filters,
        info=info.info,
        vep=vep_ht[ht.key].vep.drop("colocated_variants"),
        vqsr=keyed_filters.vqsr,
        region_flag=region_flag_expr(
            ht,
            non_par=False,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
        ),
        allele_info=keyed_filters.allele_data,
        **analyst_ht[ht.key],
    )

    return ht


def adjust_subset_alleles(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Modeled after Hail's `filter_alleles` module to adjust the allele annotation to include only alleles present in the MT.

    .. note::

        Should be used only on sparse Matrix Tables

    Uses `hl.agg.any` to determine if an allele if found in the MT. The alleles annotations will only include reference
    alleles and alternate alleles that are in MT. `mt.LA` will be adjusted to the new alleles annotation.

    :param mt: MatrixTable to subset locus alleles
    :return: MatrixTable with alleles adjusted to only those with a sample containing a non reference allele
    """
    mt = mt.annotate_rows(
        _keep_allele=hl.agg.array_agg(
            lambda i: hl.agg.any(mt.LA.contains(i[0])), hl.enumerate(mt.alleles)
        )
    )
    new_to_old = (
        hl.enumerate(mt._keep_allele).filter(lambda elt: elt[1]).map(lambda elt: elt[0])
    )
    old_to_new_dict = hl.dict(
        hl.enumerate(
            hl.enumerate(mt.alleles).filter(lambda elt: mt._keep_allele[elt[0]])
        ).map(lambda elt: (elt[1][1], elt[0]))
    )
    mt = mt.annotate_rows(
        _old_to_new=hl.bind(
            lambda d: mt.alleles.map(lambda a: d.get(a)), old_to_new_dict
        ),
        _new_to_old=new_to_old,
    )
    new_locus_alleles = hl.min_rep(
        mt.locus, mt._new_to_old.map(lambda i: mt.alleles[i])
    )
    mt = mt.key_rows_by(
        locus=new_locus_alleles.locus, alleles=new_locus_alleles.alleles
    )
    mt = mt.annotate_entries(LA=mt.LA.map(lambda x: mt._old_to_new[x]))

    return mt.drop("_keep_allele", "_new_to_old", "_old_to_new")


def create_full_subset_dense_mt(mt: hl.MatrixTable, meta_ht: hl.Table):
    """
    Create the subset dense release MatrixTable with multi-allelic variants and all sample and variant annotations.

    .. note::

        This function uses the sparse subset MT and assumes that centromeres and telomeres have already been filtered
        if necessary.

    :param mt: Sparse subset release MatrixTable
    :param meta_ht: Metadata HT to use for sample (column) annotations
    :return: Dense release MatrixTable with all row, column, and global annotations
    """
    release_ht = release_sites().ht()
    subset_freq = get_freq(subset="hgdp_tgp").ht()
    info_ht = get_info(split=True).ht()
    filters_ht = final_filter.ht()

    logger.info(
        "Adding subset's sample QC metadata to MT columns and global annotations to MT globals"
    )
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.annotate_globals(
        global_annotation_descriptions=hl.literal(GLOBAL_ANNOTATION_DICT),
        **meta_ht.drop("global_annotation_descriptions").index_globals(),
    )

    logger.info(
        "Annotate entries with het non ref status for use in the homozygous alternate depletion fix..."
    )
    mt = mt.annotate_entries(_het_non_ref=mt.LGT.is_het_non_ref())

    logger.info("Splitting multi-allelics")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    variant_annotation_ht = prepare_variant_annotations(
        mt.rows().select().select_globals()
    )

    logger.info("Computing adj and sex adjusted genotypes...")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus, mt.GT, mt.meta.sex_imputation.sex_karyotype
        ),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    logger.info(
        "Setting het genotypes at sites with > 1% AF (using precomputed v3.0 frequencies) and > 0.9 AB to homalt..."
    )
    # NOTE: Using v3.0 frequencies here and not v3.1 frequencies because 
    # the frequency code adjusted genotypes (homalt depletion fix) using v3.0 frequencies
    # https://github.com/broadinstitute/gnomad_qc/blob/efea6851a421f4bc66b73db588c0eeeb7cd27539/gnomad_qc/v3/annotations/generate_freq_data_hgdp_tgp.py#L129
    freq_ht = get_freq(version="3").ht()
    freq_ht = freq_ht.select(AF=freq_ht.freq[0].AF)

    mt = hom_alt_depletion_fix(
        mt, het_non_ref_expr=mt._het_non_ref, af_expr=freq_ht[mt.row_key].freq[0].AF
    )
    mt = mt.drop("_het_non_ref")

    logger.info(
        "Add gnomad freq from release HT, remove downsampling and subset info from freq, freq_meta, and freq_index_dict"
    )
    full_release_freq_meta = release_ht.freq_meta.collect()[0]
    freq_meta = [
        x
        for x in full_release_freq_meta
        if "downsampling" not in x and "subset" not in x
    ]
    index_keep = [
        i
        for i, x in enumerate(full_release_freq_meta)
        if "downsampling" not in x and "subset" not in x
    ]
    freq_index_dict = release_ht.freq_index_dict.collect()[0]
    freq_index_dict = {k: v for k, v in freq_index_dict.items() if v in index_keep}

    release_struct = release_ht[mt.row_key]
    mt = mt.annotate_rows(
        cohort_freq=subset_freq[mt.row_key].freq,
        gnomad_freq=release_struct.freq[: len(freq_meta)],
        gnomad_raw_qual_hists=release_struct.raw_qual_hists,
        gnomad_popmax=release_struct.popmax,
        gnomad_qual_hists=release_struct.qual_hists,
        gnomad_faf=release_struct.faf,
    )
    mt = mt.annotate_globals(
        cohort_freq_meta=subset_freq.index_globals().freq_meta,
        cohort_freq_index_dict=make_freq_index_dict(
            subset_freq.index_globals().freq_meta,
            pops=POPS_STORED_AS_SUBPOPS,
            subsets=["hgdp|tgp"],
            label_delimiter="-",
        ),
        gnomad_freq_meta=freq_meta,
        gnomad_freq_index_dict=freq_index_dict,
        gnomad_faf_index_dict=release_ht.index_globals().faf_index_dict,
        gnomad_faf_meta=release_ht.index_globals().faf_meta,
        vep_version=release_ht.index_globals().vep_version,
        vep_csq_header=release_ht.index_globals().vep_csq_header,
        dbsnp_version=release_ht.index_globals().dbsnp_version,
        filtering_model=release_ht.index_globals().filtering_model.drop("model_id"),
        inbreeding_coeff_cutoff=filters_ht.index_globals().inbreeding_coeff_cutoff,
    )

    logger.info("Add all other variant annotations from release HT...")
    mt = mt.annotate_rows(**variant_annotation_ht[mt.row_key])
    mt = mt.drop("was_split", "a_index")

    # TODO: lowqual variants how to handle them on both sparse and dense, dense is easy because I can
    #  remove the variants easily, sparse is more difficult, maybe just make sure annotation is on the variant HT,
    #  which at the moment is not the case because we just remove them.
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows((~info_ht[mt.row_key].AS_lowqual) & (hl.len(mt.alleles) > 1))

    return mt


def main(args):
    hl.init(log="/hgdp_1kg_subset.log", default_reference="GRCh38")

    if args.create_sample_meta:
        meta_ht = prepare_sample_annotations()
        meta_ht.write(hgdp_1kg_subset_annotations().path)

    if args.export_meta_txt:
        hgdp_1kg_subset_annotations().ht().export(hgdp_1kg_subset_sample_tsv())

    if args.create_subset_sparse_mt:
        # NOTE: no longer filtering to high_quality by request from Alicia Martin, but we do filter to variants in
        # high_quality samples, so how to handle that in the future?
        meta_ht = hgdp_1kg_subset_annotations().ht()
        mt = get_gnomad_v3_mt(
            key_by_locus_and_alleles=True, remove_hard_filtered_samples=False
        )

        logger.info(
            "Filtering MT columns to HGDP + TGP samples and the CHMI haploid sample (syndip)"
        )
        mt = mt.filter_cols(
            (
                meta_ht[mt.col_key].subsets.hgdp
                | meta_ht[mt.col_key].subsets.tgp
                | (mt.s == SYNDIP)
            )
        )
        mt = mt.annotate_rows(
            _telomere_or_centromere=hl.is_defined(
                telomeres_and_centromeres.ht()[mt.locus]
            ),
        )
        # Filter to rows with a non ref GT outside of telomeres and centomeres, or if there is END info
        mt = mt.filter_rows(
            (~mt._telomere_or_centromere & hl.agg.any(mt.LGT.is_non_ref()))
            | hl.agg.any(hl.is_defined(mt.END))
        )

        # On rows with END info that fall in telomeres and centomeres, keep only END info, set other entries to missing
        mt = mt.filter_entries(
            mt._telomere_or_centromere & hl.is_missing(mt.END), keep=False
        )

        # Adjust alleles and LA to include only alleles present in the subset
        mt = adjust_subset_alleles(mt)

        logger.info(
            "Note: for the finalized HGDP + TGP subset frequency, dense MT, VCFs we adjust the sex genotypes and add "
            "a fix for older GATK gVCFs with a known depletion of homozygous alternate alleles, and remove standard "
            "GATK LowQual variants"
        )

        mt = mt.drop("_telomere_or_centromere")
        mt.write(
            hgdp_1kg_subset(dense=False).path, overwrite=args.overwrite,
        )

    if args.create_subset_dense_mt:
        mt = hgdp_1kg_subset(dense=False).mt()
        mt = mt.select_entries(*SPARSE_ENTRIES)
        meta_ht = hgdp_1kg_subset_annotations().ht()
        mt = create_full_subset_dense_mt(mt, meta_ht)

        logger.info("Writing HGDP + TGP MT")
        mt.write(hgdp_1kg_subset(dense=True).path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script subsets the gnomAD v3.1 release to only HGDP and 1KG samples"
    )
    parser.add_argument(
        "--export-meta",
        help="Pull sample subset metadata and export to a .tsv",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
