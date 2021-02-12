import argparse
import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch38.gnomad import (
    GROUPS,
    HGDP_POPS,
    KG_POPS,
    SEXES,
)
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    telomeres_and_centromeres,
    seg_dup_intervals,
)
from gnomad.utils.filtering import filter_low_conf_regions
from gnomad.utils.vcf import (
    AS_FIELDS,
    index_globals,
    SITE_FIELDS,
    SPARSE_ENTRIES,
)

from gnomad_qc.v3.resources.annotations import (
    analyst_annotations,
    get_freq,
    get_info,
    vep,
)
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta as metadata
from gnomad_qc.v3.resources.release import release
from gnomad_qc.v3.resources.variant_qc import final_filter

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_subset")
logger.setLevel(logging.INFO)

AS_FIELDS.remove("InbreedingCoeff")
AS_FIELDS.extend(["AS_QUALapprox", "AS_SB_TABLE"])
SITE_FIELDS.remove("BaseQRankSum")
SITE_FIELDS.extend(["SB", "QUALapprox"])

GLOBAL_ANNOTATION_DICT = hl.struct(
    sex_imputation_ploidy_cutoffs=hl.struct(
        Description=(
            "Contains sex chromosome ploidy cutoffs used when determining sex chromosome karyotypes. Format: "
            "(upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for "
            "triple X) and (lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)."
        )
    ),
    population_inference_pca_metrics=hl.struct(
        Description=(
            "Contains the number of principal components (PCs) used when running PC-project and the minimum cutoff "
            "probability of belonging to a given population for assignment."
        )
    ),
    hard_filter_cutoffs=hl.struct(
        Description=(
            "Contains the cutoffs used for hard-filtering of samples prior to sample QC. Sample QC metrics are "
            "computed using the Hail sample_qc module on all autosomal bi-allelic SNVs and samples are removed if "
            "they are clear outliers for number of snps (n_snp), ratio of heterozygous variants to homozygous "
            "variants (r_het_hom_var), number of singletons (n_singleton), and mean coverage on chromosome 20 (cov). "
            "Additionally, we filter based on outliers of some Picard metrics: % contamination (freemix), % chimera, "
            "and median insert size."
        )
    ),
    cohort_freq_meta=hl.struct(
        Description=(
            "HGDP and 1KG frequency metadata. Contains the information about what frequency aggregation is in each "
            "element of the cohort_freq array row annotation."
        )
    ),
    gnomad_freq_meta=hl.struct(
        Description=(
            "gnomAD release frequency metadata. Contains the information about what frequency aggregation is in each "
            "element of the gnomad_freq array row annotation."
        )
    ),
    cohort_freq_index_dict=hl.struct(
        Description=(
            "Dictionary where keys are combinations of all possible groupings used in the HGDP + 1KG frequency "
            "aggregations (group: adj/raw, pop: HGDP or 1KG subpopulation, sex: sex karyotype) and values are the "
            "index into the cohort_freq array row annotation."
        )
    ),
    gnomad_freq_index_dict=hl.struct(
        Description=(
            "Dictionary where keys are combinations of all possible groupings used in gnomAD frequency aggregations "
            "(group: adj/raw, pop: gnomAD inferred global population, sex: sex karyotype) and values are the index "
            "into the gnomad_freq array row annotation."
        )
    ),
    gnomad_faf_index_dict=hl.struct(
        Description=(
            "Dictionary where keys are combinations of all possible groupings used in gnomAD filtering allele "
            "frequency (using Poisson 99% CI) aggregations (group: adj/raw, pop: gnomAD inferred global population, "
            "sex: sex karyotype) and values are the index into the gnomad_faf array row annotation."
        )
    ),
    gnomad_faf_meta=hl.struct(
        Description=(
            "gnomAD release filtering allele frequency (using Poisson 99% CI) metadata. Contains the information "
            "about what faf aggregation is in each element of the gnomad_faf array row annotation."
        )
    ),
    vep_version=hl.struct(Description="VEP version."),
    vep_csq_header=hl.struct(Description="VEP header for VCF export."),
    dbsnp_version=hl.struct(Description="dbSNP version."),
    filtering_model=hl.struct(
        Description="Contains information about the filtering model used and specific cutoffs.",
        sub_globals=hl.struct(
            model_name=hl.struct(
                Description=(
                    "Variant filtering model name. This is used as the filter in the filters row annotation to "
                    "indicate the variant was filtered during variant QC."
                )
            ),
            score_name=hl.struct(
                Description="Name of score used for variant filtering."
            ),
            snv_cutoff=hl.struct(
                Description="Information about cutoff used for SNV filtering.",
                sub_globals=hl.struct(
                    bin=hl.struct(
                        Description="Percentile cutoff for SNVs for filtering."
                    ),
                    min_score=hl.struct(
                        Description="Min score (score_name) at percentile cutoff for SNV filtering."
                    ),
                ),
            ),
            indel_cutoff=hl.struct(
                Description="Information about cutoff used for indel filtering.",
                sub_globals=hl.struct(
                    bin=hl.struct(
                        Description="Percentile cutoff for indels for filtering."
                    ),
                    min_score=hl.struct(
                        Description="Min score (score_name) at percentile cutoff for indel filtering."
                    ),
                ),
            ),
            snv_training_variables=hl.struct(
                Description="Variant annotations used as features in SNV filtering model."
            ),
            indel_training_variables=hl.struct(
                Description="Variant annotations used as features in indel filtering model."
            ),
        ),
    ),
    inbreeding_coeff_cutoff=hl.struct(
        Description="Hard-filter cutoff for InbreedingCoeff on variants (rows)."
    ),
)

SAMPLE_ANNOTATION_DICT = hl.struct(
    s=hl.struct(Description="Sample ID."),
    bam_metrics=hl.struct(
        Description="Sample level metrics obtained from BAMs/CRAMs.",
        sub_annotations=hl.struct(
            pct_bases_20x=hl.struct(
                Description="The fraction of bases that attained at least 20X sequence coverage in post-filtering bases."
            ),
            pct_chimeras=hl.struct(
                Description=(
                    "The fraction of reads that map outside of a maximum insert size (usually 100kb) or that have the "
                    "two ends mapping to different chromosomes."
                )
            ),
            freemix=hl.struct(Description="Estimate of contamination (0-100 scale)."),
            mean_coverage=hl.struct(
                Description="The mean coverage in bases of the genome territory, after all filters are applied."
            ),
            median_coverage=hl.struct(
                Description="The median coverage in bases of the genome territory, after all filters are applied."
            ),
            mean_insert_size=hl.struct(
                Description=(
                    "The mean insert size of the 'core' of the distribution. Artefactual outliers in the distribution "
                    "often cause calculation of nonsensical mean and stdev values. To avoid this the distribution is "
                    "first trimmed to a 'core' distribution of +/- N median absolute deviations around the median "
                    "insert size. By default N=10, but this is configurable."
                )
            ),
            median_insert_size=hl.struct(
                Description="The MEDIAN insert size of all paired end reads where both ends mapped to the same chromosome."
            ),
            pct_bases_10x=hl.struct(
                Description="The fraction of bases that attained at least 10X sequence coverage in post-filtering bases."
            ),
        ),
    ),
    subsets=hl.struct(
        Description="Contains information about the specific cohort the sample is from: HGDP or 1KG (tgp).",
        sub_annotations=hl.struct(
            tgp=hl.struct(Description="True if sample is from the 1KG dataset."),
            hgdp=hl.struct(Description="True if the sample is from the HGDP dataset."),
        ),
    ),
    sex_imputation=hl.struct(
        Description="Struct containing sex imputation information.",
        sub_annotations=hl.struct(
            f_stat=hl.struct(
                Description="Inbreeding coefficient (excess heterozygosity) on chromosome X."
            ),
            n_called=hl.struct(Description="Number of variants with a genotype call."),
            expected_homs=hl.struct(Description="Expected number of homozygotes."),
            observed_homs=hl.struct(Description="Observed number of homozygotes."),
            chr20_mean_dp=hl.struct(
                Description="Sample's mean depth across chromosome 20."
            ),
            chrX_mean_dp=hl.struct(
                Description="Sample's mean depth across chromosome X."
            ),
            chrY_mean_dp=hl.struct(
                Description="Sample's mean depth across chromosome Y."
            ),
            chrX_ploidy=hl.struct(
                Description="Sample's chromosome X ploidy (chrX_mean_dp normalized using chr20_mean_dp)."
            ),
            chrY_ploidy=hl.struct(
                Description="Sample's chromosome Y ploidy (chrY_mean_dp normalized using chr20_mean_dp)."
            ),
            X_karyotype=hl.struct(Description="Sample's chromosome X karyotype."),
            Y_karyotype=hl.struct(Description="Sample's chromosome Y karyotype."),
            sex_karyotype=hl.struct(
                Description="Sample's sex karyotype (combined X and Y karyotype)."
            ),
        ),
    ),
    sample_qc=hl.struct(
        Description="Struct containing sample QC metrics calculated using hl.sample_qc().",
        sub_annotations=hl.struct(
            n_hom_ref=hl.struct(Description="Number of homozygous reference calls."),
            n_het=hl.struct(Description="Number of heterozygous calls."),
            n_hom_var=hl.struct(Description="Number of homozygous alternate calls."),
            n_non_ref=hl.struct(Description="Sum of n_het and n_hom_var."),
            n_snp=hl.struct(Description="Number of SNP alternate alleles."),
            n_insertion=hl.struct(Description="Number of insertion alternate alleles."),
            n_deletion=hl.struct(Description="Number of deletion alternate alleles."),
            n_transition=hl.struct(
                Description="Number of transition (A-G, C-T) alternate alleles."
            ),
            n_transversion=hl.struct(
                Description="Number of transversion alternate alleles."
            ),
            r_ti_tv=hl.struct(Description="Transition/Transversion ratio."),
            r_het_hom_var=hl.struct(Description="Het/HomVar call ratio."),
            r_insertion_deletion=hl.struct(
                Description="Insertion/Deletion allele ratio."
            ),
        ),
    ),
    population_inference=hl.struct(
        Description=(
            "Struct containing ancestry information assigned by applying a Principal component analysis (PCA) on "
            "gnomAD samples and using those PCs in a random forest classifier trained on known gnomAD ancestry labels."
        ),
        sub_annotations=hl.struct(
            pca_scores=hl.struct(
                Description="Sample's scores for each gnomAD population PC."
            ),
            pop=hl.struct(Description="Sample's gnomAD PC project population label."),
            prob_afr=hl.struct(
                Description="Random forest probability that the sample is of African/African-American ancestry."
            ),
            prob_ami=hl.struct(
                Description="Random forest probability that the sample is of Amish ancestry."
            ),
            prob_amr=hl.struct(
                Description="Random forest probability that the sample is of Latino ancestry."
            ),
            prob_asj=hl.struct(
                Description="Random forest probability that the sample is of Ashkenazi Jewish ancestry."
            ),
            prob_eas=hl.struct(
                Description="Random forest probability that the sample is of East Asian ancestry."
            ),
            prob_fin=hl.struct(
                Description="Random forest probability that the sample is of Finnish ancestry."
            ),
            prob_mid=hl.struct(
                Description="Random forest probability that the sample is of Middle Eastern ancestry."
            ),
            prob_nfe=hl.struct(
                Description="Random forest probability that the sample is of Non-Finnish European ancestry."
            ),
            prob_oth=hl.struct(
                Description="Random forest probability that the sample is of Other ancestry."
            ),
            prob_sas=hl.struct(
                Description="Random forest probability that the sample is of South Asian ancestry."
            ),
        ),
    ),
    labeled_subpop=hl.struct(
        Description="The sample's population label supplied by HGDP or 1KG."
    ),
    gnomad_release=hl.struct(
        Description=(
            "Indicates whether the sample was included in the gnomAD release dataset. For the full gnomAD release, "
            "relatedness inference is performed on the full dataset, and release samples are chosen in a way that "
            "maximizes the number of samples retained while filtering the dataset to include only samples with less "
            "than second-degree relatedness. For the HGDP + 1KG subset, all samples passing all other sample QC "
            "metrics are retained."
        )
    ),
)

POPS = KG_POPS + HGDP_POPS

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


def make_freq_index_dict(freq_meta: List[Dict[str, str]]) -> Dict[str, int]:
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array

    :param List of Dict freq_meta: Global annotation containing the set of groupings for each element of the freq array
        (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    """
    return {
        **index_globals(freq_meta, dict(group=GROUPS), label_delimiter="-"),
        **index_globals(freq_meta, dict(group=GROUPS, pop=POPS), label_delimiter="-"),
        **index_globals(freq_meta, dict(group=GROUPS, sex=SEXES), label_delimiter="-"),
        **index_globals(
            freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES), label_delimiter="-"
        ),
    }


def region_flag_expr(
    t: Union[hl.Table, hl.MatrixTable],
    non_par: bool = True,
    prob_regions: Dict[str, hl.Table] = None,
) -> hl.expr.StructExpression:
    """
    Creates `region_flag` struct.
    Struct contains flags for problematic regions (i.e., LCR, decoy, segdup, and nonpar regions).
    .. note::
        No hg38 resources for decoy or self chain available yet.
    :param Table/MatrixTable t: Input Table/MatrixTable.
    :return: `region_flag` struct row annotation.
    :rtype: hl.expr.StructExpression
    """

    prob_flags_expr = (
        {"non_par": (t.locus.in_x_nonpar() | t.locus.in_y_nonpar())} if non_par else {}
    )

    if prob_regions is not None:
        prob_flags_expr.update(
            {
                region_name: hl.is_defined(region_table[t.locus])
                for region_name, region_table in prob_regions.items()
            }
        )

    return hl.struct(**prob_flags_expr)


def prepare_sample_annotations() -> hl.Table:
    """


    :return:
    """
    logger.info(
        "Sub-setting and modifying sample QC metadata to desired globals and annotations"
    )
    meta_ht = metadata.ht()
    meta_ht = meta_ht.select_globals(
        sex_imputation_ploidy_cutoffs=meta_ht.sex_imputation_ploidy_cutoffs,
        population_inference_pca_metrics=hl.struct(
            n_pcs=meta_ht.population_inference_pca_metrics.n_pcs,
            min_prob=meta_ht.population_inference_pca_metrics.min_prob,
        ),
        hard_filter_cutoffs=hl.struct(
            min_cov=15,
            max_n_snp=3.75e6,
            min_n_snp=2.4e6,
            max_n_singleton=1e5,
            max_r_het_hom_var=3.3,
            max_pct_contamination=5.00,
            max_pct_chimera=5.00,
            min_median_insert_size=250,
        ),
    )
    meta_ht = meta_ht.select(
        bam_metrics=meta_ht.bam_metrics,
        subsets=meta_ht.subsets.select("tgp", "hgdp"),
        sex_imputation=meta_ht.sex_imputation.drop("is_female"),
        sample_qc=meta_ht.sample_qc.select(*SAMPLE_QC_METRICS),
        population_inference=meta_ht.population_inference.drop(
            "training_pop", "training_pop_all"
        ),
        labeled_subpop=meta_ht.project_meta.project_subpop,
        gnomad_release=meta_ht.release,
        high_quality=meta_ht.high_quality,
    )

    return meta_ht


def prepare_variant_annotations(ht: hl.Table, filter_lowqual: bool = True) -> hl.Table:
    """
    Load and join all Tables with variant annotations.

    :param Table ht:
    :param bool filter_lowqual: If true the returned Table will filter out lowqual variants.
    :return: Table containing joined annotations.
    :rtype: hl.Table
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
        f"The following fields are not found in the info HT: {missing_info_fields}"
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

    filters = filters_ht[ht.key]
    info = info_ht[ht.key]
    ht = ht.annotate(
        a_index=info.a_index,
        was_split=info.was_split,
        rsid=dbsnp_ht[ht.key].rsid,
        filters=filters.filters,
        info=info.info,
        vep=vep_ht[ht.key].vep.drop("colocated_variants"),
        vqsr=filters.vqsr,
        region_flag=region_flag_expr(
            ht,
            non_par=False,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht(),},
        ),
        allele_info=filters.allele_data,
        **analyst_ht[ht.key],
    )

    return ht


# TODO: Should we put this hotfix into a common function somewhere? I'm not sure the best place to put it,
#  but we use it a couple times
def hom_alt_depletion_fix(mt, af_expr, af_cutoff=0.01, ab_cutoff=0.9):
    """
    Hotfix for the depletion of homozygous alternate genotypes

    :param mt:
    :param af_expr:
    :param af_cutoff:
    :param ab_cutoff:
    :return:
    """
    return mt.annotate_entries(
        GT=hl.cond(
            (af_expr > af_cutoff) & mt.GT.is_het() & ((mt.AD[1] / mt.DP) > ab_cutoff),
            hl.call(1, 1),
            mt.GT,
        )
    )


def main(args):
    hl.init(log="/subset.log", default_reference="GRCh38")

    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
    mt = mt.select_entries(*SPARSE_ENTRIES)
    meta_ht = prepare_sample_annotations()
    release_ht = release(public=False).ht()
    subset_freq = get_freq(subset="hgdp_tgp")
    info_ht = get_info(split=True).ht()
    filters_ht = final_filter.ht()

    logger.info("Filtering MT columns to high quality and HGDP + TGP samples")
    mt = mt.filter_cols(
        meta_ht[mt.col_key].high_quality
        & (meta_ht[mt.col_key].subsets.hgdp | meta_ht[mt.col_key].subsets.tgp)
    )

    logger.info(
        "Adding desired sample QC metadata to MT columns and global annotations to MT globals"
    )
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.annotate_globals(
        global_annotation_descriptions=hl.literal(GLOBAL_ANNOTATION_DICT),
        sample_annotation_descriptions=hl.literal(SAMPLE_ANNOTATION_DICT),
        **meta_ht.index_globals(),
    )

    if args.export_meta:
        meta = meta_ht.semi_join(mt.cols())
        meta.export(
            f"{args.output_path}/metadata.tsv.bgz"
        )  # TODO: add path to resources

    # TODO: probably checkpoint filtered mt before changing the alleles

    mt = mt.annotate_rows(
        _site_has_non_ref=hl.agg.any(mt.LGT.is_non_ref()),
        _site_has_end=hl.agg.any(hl.is_defined(mt.END)),
    )

    mt = mt.filter_rows(mt._site_has_non_ref | mt._site_has_end)

    # For those sites with no non ref variants only END info, remove the alt allele,
    # because it isn't actually present in this subset
    mt = mt.key_rows_by(
        "locus",
        alleles=hl.if_else(
            ~mt._site_has_non_ref & mt._site_has_end, [mt.alleles[0]], mt.alleles
        ),
    )

    logger.info("Splitting multi-allelics")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Add HGDP + TGP subset frequency info")
    mt = mt.annotate_rows(
        AS_lowqual=info_ht[mt.row_key].AS_lowqual,
        telomere_or_centromere=hl.is_defined(telomeres_and_centromeres.ht()[mt.locus]),
        cohort_freq=subset_freq[mt.row_key].freq,
    )

    # If the site has some non-ref variants, but this specific allele has none, remove this allele.
    # This removes alleles not in subset, but will keep and END info at the site
    mt = mt.filter_rows(
        mt._site_has_non_ref & hl.is_missing(mt.cohort_freq), keep=False
    )
    variant_annotation_ht = prepare_variant_annotations(
        mt.rows().select().select_globals(), release_ht.select("freq")
    )

    logger.info("Removing centromeres and telomeres")
    variant_annotation_ht = filter_low_conf_regions(
        variant_annotation_ht,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False,
        filter_telomeres_and_centromeres=True,
    )

    logger.info(
        "Setting het genotypes at sites with >1% AF (using precomputed v3.1 frequencies) and > 0.9 AB to homalt..."
    )
    hom_alt_depletion_fix(mt, af_expr=release_ht[mt.row_key].freq[0].AF)

    logger.info(
        "Add gnomad freq from release HT, remove downsampling and subset info from freq, freq_meta, and freq_index_dict"
    )
    freq_meta = [
        x
        for x in release_ht.freq_meta.collect()[0]
        if "downsampling" not in x and "subset" not in x
    ]
    index_keep = [
        i
        for i, x in enumerate(release_ht.freq_meta.collect()[0])
        if "downsampling" not in x and "subset" not in x
    ]
    freq_index_dict = release_ht.freq_index_dict.collect()[0]
    freq_index_dict = {k: v for k, v in freq_index_dict.items() if v in index_keep}

    logger.info("Filtering freq HT to only variants found in subset...")
    mt = mt.annotate_rows(_is_subset_variant=hl.agg.any(mt.GT.is_non_ref()))
    ht_variants = mt.rows()
    release_ht = release_ht.filter(ht_variants[release_ht.key]._is_subset_variant)
    release = release_ht[mt.row_key]
    variant_annotation_ht = variant_annotation_ht.filter(
        ht_variants[variant_annotation_ht.key]._is_subset_variant
    )
    mt = mt.annotate_rows(
        cohort_freq=hl.or_missing(
            ~mt.AS_lowqual & ~mt.telomere_or_centromere, mt.cohort_freq
        ),
        gnomad_freq=release.freq[:len(freq_meta)],
        gnomad_raw_qual_hists=release.raw_qual_hists,
        gnomad_popmax=release.popmax,
        gnomad_qual_hists=release.qual_hists,
        gnomad_faf=release.faf,
    )
    mt = mt.annotate_globals(cohort_freq_meta=subset_freq.index_globals().freq_meta)
    mt = mt.annotate_globals(
        cohort_freq_index_dict=make_freq_index_dict(mt),
        gnomad_freq_meta=freq_meta,
        gnomad_freq_index_dict=freq_index_dict,
        gnomad_faf_index_dict=release_ht.index_globals().faf_index_dict,
        gnomad_faf_meta=release_ht.index_globals().faf_meta,
        vep_version=release_ht.index_globals().vep_version,
        vep_csq_header=release_ht.index_globals().vep_csq_header,
        dbsnp_version=release_ht.index_globals().dbsnp_version,
        filtering_model=release_ht.index_globals().filtering_model.drop("model_id"),
        inbreeding_coeff_cutoff=filters_ht.index_globals().inbreeding_coeff_cutoff,  # TODO: Should have this added onto the final_filter Table
    )

    logger.info("Add all other variant annotations from release HT (dropping freq)")
    mt = mt.annotate_rows(**variant_annotation_ht[mt.row_key])

    mt = mt.drop(
        "_is_subset_variant",
        "_site_has_non_ref",
        "_site_has_end",
        "was_split",
        "a_index",
    )

    logger.info("Writing HGDP + TGP MT")
    mt.write(
        "gs://gnomad/release/3.1/mt/genomes/hgdp_1kg.gnomad.v3.1.mt", overwrite=True
    )  # TODO replace path


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script subsets gnomAD using a list of samples or population"
    )
    parser.add_argument(
        "--export-meta",
        help="Pull sample subset metadata and export to a .tsv",
        action="store_true",
    )

    parser.add_argument(
        "--num-vcf-shards", help="Number of shards in output VCF", type=int
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
