import argparse
import logging
from typing import Dict

import hail as hl

from gnomad.resources.grch38.gnomad import POPS_STORED_AS_SUBPOPS
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
    telomeres_and_centromeres,
)
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import get_adj_expr, region_flag_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.release import make_freq_index_dict
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vep import vep_or_lookup_vep
from gnomad.utils.vcf import (
    AS_FIELDS,
    SITE_FIELDS,
    SPARSE_ENTRIES,
)

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.create_release.hgdp_tgp_constants import (
    GLOBAL_ANNOTATIONS,
    GLOBAL_SAMPLE_ANNOTATIONS,
    GLOBAL_VARIANT_ANNOTATIONS,
    SAMPLE_ANNOTATIONS,
    VARIANT_ANNOTATIONS,
)
from gnomad_qc.v3.resources.annotations import (
    get_freq,
    get_info,
)
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.release import (
    release_sites,
    hgdp_tgp_subset,
    hgdp_tgp_subset_annotations,
    hgdp_tgp_subset_sample_tsv,
)
from gnomad_qc.v3.resources.sample_qc import (
    hgdp_tgp_meta,
    hgdp_tgp_pcs,
    hgdp_tgp_pop_outliers,
    hgdp_tgp_relatedness,
    hgdp_tgp_related_samples_to_drop,
)
from gnomad_qc.v3.resources.variant_qc import final_filter, SYNDIP
from gnomad_qc.v3.utils import hom_alt_depletion_fix

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("create_subset")
logger.setLevel(logging.INFO)

AS_FIELDS.remove("InbreedingCoeff")
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


def convert_heterogeneous_dict_to_struct(global_dict: Dict) -> hl.struct:
    """
    Convert heterogeneous dictionary (one with unspecified levels of nested dicts) into a multi-level Hail struct.

    :param global_dict: Heterogeneous dictionary to convert into a Hail struct.
    :return: Global dictionary converted to struct.
    """
    if isinstance(global_dict, dict):
        return hl.struct(
            **{
                k: convert_heterogeneous_dict_to_struct(global_dict[k])
                for k in global_dict
            }
        )
    else:
        return global_dict


def get_sample_qc_filter_struct_expr(ht: hl.Table) -> hl.struct:
    """
    Create expression for the final filter struct indicating gnomAD hard filtered samples and HGDP + 1KG/TGP subcontinental PCA outliers.

    .. note::

        This assumes that the input Hail Table contains a `gnomad_sample_filters` annotation with
        `hard_filters` and `hard_filtered` as sub-annotations.

    :param ht: Input Table containing hard filter information.
    :return: Struct expression for sample QC filters.
    """
    logger.info("Read in population-specific PCA outliers...")
    hgdp_tgp_pop_outliers_ht = hgdp_tgp_pop_outliers.ht()
    set_to_remove = hgdp_tgp_pop_outliers_ht.s.collect(_localize=False)

    num_outliers = hl.eval(hl.len(set_to_remove))
    num_outliers_found = ht.filter(set_to_remove.contains(ht["s"])).count()
    if hl.eval(hl.len(set_to_remove)) != num_outliers:
        raise ValueError(
            f"Expected {num_outliers} samples to be labeled as population PCA outliers, but found {num_outliers_found}"
        )

    return hl.struct(
        hard_filters=ht.gnomad_sample_filters.hard_filters,
        hard_filtered=ht.gnomad_sample_filters.hard_filtered,
        pop_outlier=hl.if_else(
            ht.s == SYNDIP, hl.missing(hl.tbool), set_to_remove.contains(ht["s"])
        ),
    )


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Create Table of all related samples and the relatedness information for all samples they are related to.

    Return Table keyed by sample with a `related_samples` annotation that is a set containing a struct of relatedness
    information for each sample it is related to. Each struct has the following: kin, ibd0, ibd1, and ibd2.

    :param relatedness_ht: Table with inferred relationship information output by pc_relate.
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationship information annotated as a struct.
    """
    relationship_struct = hl.struct(
        kin=relatedness_ht.kin,
        ibd0=relatedness_ht.ibd0,
        ibd1=relatedness_ht.ibd1,
        ibd2=relatedness_ht.ibd2,
    )

    relatedness_ht_i = relatedness_ht.group_by(s=relatedness_ht.i.s).aggregate(
        related_samples=hl.agg.collect_as_set(
            hl.struct(s=relatedness_ht.j.s, **relationship_struct)
        )
    )

    relatedness_ht_j = relatedness_ht.group_by(s=relatedness_ht.j.s).aggregate(
        related_samples=hl.agg.collect_as_set(
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
        (meta_ht.subsets.hgdp | meta_ht.subsets.tgp | (meta_ht.s == SYNDIP))
    )

    meta_ht = meta_ht.select_globals(
        global_annotation_descriptions=convert_heterogeneous_dict_to_struct(
            GLOBAL_SAMPLE_ANNOTATIONS
        ),
        sample_annotation_descriptions=convert_heterogeneous_dict_to_struct(
            SAMPLE_ANNOTATIONS
        ),
        gnomad_sex_imputation_ploidy_cutoffs=meta_ht.sex_imputation_ploidy_cutoffs,
        gnomad_population_inference_pca_metrics=hl.struct(
            n_pcs=meta_ht.population_inference_pca_metrics.n_pcs,
            min_prob=meta_ht.population_inference_pca_metrics.min_prob,
        ),
        sample_hard_filter_cutoffs=meta_ht.hard_filter_cutoffs,
        gnomad_sample_qc_metric_outlier_cutoffs=meta_ht.outlier_detection_metrics,
        gnomad_age_distribution=release_sites(public=True)
        .versions["3.1.1"]
        .ht()
        .index_globals()
        .age_distribution,
    )

    # Use a pre-computed relatedness HT from the Martin group - details of its creation are
    # here: https://github.com/atgu/hgdp_tgp/blob/master/pca_subcont.ipynb
    relatedness_ht = hgdp_tgp_relatedness.ht()
    subset_samples = meta_ht.s.collect(_localize=False)
    relatedness_ht = relatedness_ht.filter(
        subset_samples.contains(relatedness_ht.i.s)
        & subset_samples.contains(relatedness_ht.j.s)
    )

    relatedness_ht = get_relatedness_set_ht(relatedness_ht)

    # Note: Needs to be done before adding the relatedness info because the relatedness HT doesn't have the prefix
    logger.info(
        "Removing 'v3.1::' from the sample names, these were added because there are duplicates of some 1KG/TGP samples"
        " in the full gnomAD dataset..."
    )
    meta_ht = meta_ht.key_by(s=meta_ht.s.replace("v3.1::", ""))
    meta_ht = meta_ht.select(
        bam_metrics=meta_ht.bam_metrics,
        sample_qc=meta_ht.sample_qc.select(*SAMPLE_QC_METRICS),
        gnomad_sex_imputation=meta_ht.sex_imputation.annotate(
            **meta_ht.sex_imputation.impute_sex_stats
        ).drop("is_female", "impute_sex_stats"),
        gnomad_population_inference=meta_ht.population_inference.drop(
            "training_pop", "training_pop_all"
        ),
        gnomad_sample_qc_residuals=meta_ht.sample_qc.select(
            *[
                k
                for k in meta_ht.sample_qc.keys()
                if "_residual" in k and k.replace("_residual", "") in SAMPLE_QC_METRICS
            ]
        ),
        gnomad_sample_filters=meta_ht.sample_filters.select(
            "hard_filters", "hard_filtered", "release_related", "qc_metrics_filters"
        ),
        gnomad_high_quality=meta_ht.high_quality,
        gnomad_release=meta_ht.release,
        relatedness_inference=hl.struct(
            related_samples=hl.or_else(
                relatedness_ht[meta_ht.key].related_samples,
                hl.empty_set(
                    hl.dtype(
                        "struct{s: str, kin: float64, ibd0: float64, ibd1: float64, ibd2: float64}"
                    )
                ),
            ),
            related=hl.is_defined(hgdp_tgp_related_samples_to_drop.ht()[meta_ht.key]),
        ),
        gnomad_labeled_subpop=meta_ht.project_meta.project_subpop,
    )

    logger.info("Loading additional sample metadata from Martin group...")
    hgdp_tgp_meta_ht = hgdp_tgp_meta.ht()
    hgdp_tgp_meta_ht = hgdp_tgp_meta_ht.select(
        project=hgdp_tgp_meta_ht.hgdp_tgp_meta.Project,
        study_region=hgdp_tgp_meta_ht.hgdp_tgp_meta.Study.region,
        population=hgdp_tgp_meta_ht.hgdp_tgp_meta.Population,
        genetic_region=hgdp_tgp_meta_ht.hgdp_tgp_meta.Genetic.region,
        latitude=hgdp_tgp_meta_ht.hgdp_tgp_meta.Latitude,
        longitude=hgdp_tgp_meta_ht.hgdp_tgp_meta.Longitude,
        hgdp_technical_meta=hgdp_tgp_meta_ht.bergstrom.select("source", "library_type"),
    )
    hgdp_tgp_meta_ht = hgdp_tgp_meta_ht.union(
        hl.Table.parallelize(
            [hl.struct(s=SYNDIP, project="synthetic_diploid_truth_sample")]
        ).key_by("s"),
        unify=True,
    )

    logger.info("Adding sample QC struct and sample metadata from Martin group...")
    meta_ht = meta_ht.annotate(sample_filters=get_sample_qc_filter_struct_expr(meta_ht))
    hgdp_tgp_pcs_indexed = hgdp_tgp_pcs.ht()[meta_ht.key]

    meta_ht = meta_ht.transmute(
        hgdp_tgp_meta=hl.struct(
            **hgdp_tgp_meta_ht[meta_ht.key],
            global_pca_scores=hgdp_tgp_pcs_indexed.pca_preoutlier_global_scores,
            subcontinental_pca=hl.struct(
                pca_scores=hgdp_tgp_pcs_indexed.pca_scores,
                pca_scores_outliers_removed=hgdp_tgp_pcs_indexed.pca_scores_outliers_removed,
                outlier=meta_ht.sample_filters.pop_outlier,
            ),
            gnomad_labeled_subpop=meta_ht.gnomad_labeled_subpop,
        ),
        high_quality=~meta_ht.sample_filters.hard_filtered
        & ~meta_ht.sample_filters.pop_outlier,
    )

    return meta_ht


# TODO: Might be good to generalize this because a similar function is used in creating the release sites HT.
def prepare_variant_annotations(
    ht: hl.Table, filter_lowqual: bool = True, vep_version: str = "101"
) -> hl.Table:
    """
    Load and join all Tables with variant annotations.

    :param ht: Input HT to add variant annotations to.
    :param filter_lowqual: If True, filter out lowqual variants using the info HT's AS_lowqual.
    :param vep_version: Version of VEPed context Table to use (if None, the default vep_context resource will be used).
    :return: Table containing joined annotations.
    """
    logger.info("Loading annotation tables...")
    filters_ht = final_filter(hgdp_tgp_subset=True).ht()
    # vep_ht = vep.ht()  # Commented out for v3.1.2 release because annotation file has been removed
    dbsnp_ht = dbsnp.ht().select("rsid")
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    subset_freq = get_freq(subset="hgdp-tgp").ht()
    release_ht = release_sites(public=True).versions["3.1.1"].ht()

    # NOTE: Added for v3.1.2 release because this annotation was removed and not a full duplicate of variants in the release HT
    vep_ht = vep_or_lookup_vep(ht, vep_version=vep_version)
    vep_ht = vep_ht.annotate_globals(version=f"v{vep_version}")

    if file_exists(get_info().path):
        info_ht = get_info().ht()
    elif not file_exists(get_info().path) and file_exists(get_info(split=False).path):
        from gnomad_qc.v3.annotations.generate_qc_annotations import split_info

        info_ht = split_info()
    else:
        raise DataException("There is no available split or unsplit info HT for use!")
    info_ht = info_ht.drop("old_locus", "old_alleles")

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
    keyed_filters = filters_ht[info_ht.key]
    info_ht = info_ht.transmute(
        info=info_ht.info.select(
            *select_info_fields,
            AS_SOR=keyed_filters.AS_SOR,
            SOR=keyed_filters.SOR,
            transmitted_singleton=keyed_filters.transmitted_singleton,
            omni=keyed_filters.omni,
            mills=keyed_filters.mills,
            monoallelic=keyed_filters.monoallelic,
            InbreedingCoeff=release_ht[
                info_ht.key
            ].info.InbreedingCoeff,  # NOTE: Changed to use release HT instead of freq
            **{f"{score_name}": keyed_filters[f"{score_name}"]},
        )
    )

    logger.info(
        "Preparing gnomad freq information from the release HT: removing downsampling and subset info from freq, "
        "freq_meta, and freq_index_dict"
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

    logger.info("Assembling all variant annotations...")
    filters_ht = filters_ht.annotate(
        allele_info=hl.struct(
            variant_type=filters_ht.variant_type,
            allele_type=filters_ht.allele_type,
            n_alt_alleles=filters_ht.n_alt_alleles,
            was_mixed=filters_ht.was_mixed,
        ),
    )

    keyed_filters = filters_ht[ht.key]
    keyed_release = release_ht[ht.key]
    keyed_info = info_ht[ht.key]
    ht = ht.annotate(
        a_index=keyed_info.a_index,
        was_split=keyed_info.was_split,
        rsid=dbsnp_ht[ht.key].rsid,
        filters=keyed_filters.filters,
        info=keyed_info.info,
        vep=vep_ht[ht.key].vep.drop("colocated_variants"),
        vqsr=keyed_filters.vqsr,
        region_flag=region_flag_expr(
            ht,
            non_par=False,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
        ),
        allele_info=keyed_filters.allele_info,
        hgdp_tgp_freq=subset_freq[ht.key].freq,
        gnomad_freq=[keyed_release.freq[i] for i in index_keep],
        gnomad_popmax=keyed_release.popmax,
        gnomad_faf=keyed_release.faf,
        gnomad_raw_qual_hists=keyed_release.raw_qual_hists,
        gnomad_qual_hists=keyed_release.qual_hists,
        gnomad_age_hist_het=keyed_release.age_hist_het,
        gnomad_age_hist_hom=keyed_release.age_hist_hom,
        cadd=keyed_release.cadd,
        revel=keyed_release.revel,
        splice_ai=keyed_release.splice_ai,
        primate_ai=keyed_release.primate_ai,
        AS_lowqual=keyed_info.AS_lowqual,
        telomere_or_centromere=hl.is_defined(telomeres_and_centromeres.ht()[ht.locus]),
    )

    logger.info("Adding global variant annotations...")
    # The `freq_meta` global annotation on subset frequency Tables include a "subset" key in each element. 
    # E.g., `{'group': 'adj', 'pop': 'cdx', 'subset': 'hgdp|tgp'}`
    # This needs to be removed for the HGDP + 1KG variant annotations
    hgdp_tgp_freq_meta = [
        {k: v for k, v in x.items() if k != "subset"}
        for x in hl.eval(subset_freq.freq_meta)
    ]

    ht = ht.annotate_globals(
        global_annotation_descriptions=convert_heterogeneous_dict_to_struct(
            GLOBAL_VARIANT_ANNOTATIONS
        ),
        variant_annotation_descriptions=convert_heterogeneous_dict_to_struct(
            VARIANT_ANNOTATIONS
        ),
        hgdp_tgp_freq_meta=hgdp_tgp_freq_meta,
        hgdp_tgp_freq_index_dict=make_freq_index_dict(
            hgdp_tgp_freq_meta, pops=POPS_STORED_AS_SUBPOPS, label_delimiter="-",
        ),
        gnomad_freq_meta=freq_meta,
        gnomad_freq_index_dict=freq_index_dict,
        gnomad_faf_index_dict=release_ht.index_globals().faf_index_dict,
        gnomad_faf_meta=release_ht.index_globals().faf_meta,
        vep_version=release_ht.index_globals().vep_version,
        vep_csq_header=release_ht.index_globals().vep_csq_header,
        dbsnp_version=release_ht.index_globals().dbsnp_version,
        variant_filtering_model=release_ht.index_globals().filtering_model.drop(
            "model_id"
        ),
        variant_inbreeding_coeff_cutoff=filters_ht.index_globals().inbreeding_coeff_cutoff,
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

    def split_shuffle(mt: hl.MatrixTable) -> hl.MatrixTable:
        """
        Split rows that do and do not have the same new and old locus prior to the row rekey.
        
        Most rows will have the same new and old locus annotations.

        Re-keying rows triggers a shuffle in hail.
        This function is needed to force the shuffle to only happen on the small number of rows that have a different new and
        old locus. Fixing a problem discussed in more detail here: https://discuss.hail.is/t/preventing-a-shuffle-error/2261/6

        :param mt: MatrixTable to rekey
        :return: Rekeyed MatrixTable
        """
        # Filter to rows that have the same old and new locus annotations
        # These rows do not trigger a shuffle when re-keyed
        mt1 = mt.filter_rows(mt.locus == mt.new_locus)
        mt1 = mt1.key_rows_by(locus=mt1.new_locus, alleles=mt1.new_alleles)
        # Filter to row with different old and new locus annotations
        # These rows will trigger a shuffle when re-keyed
        mt2 = mt.filter_rows(mt.locus != mt.new_locus)
        mt2 = mt2.key_rows_by(locus=mt2.new_locus, alleles=mt2.new_alleles)

        return mt1.union_rows(mt2)

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
    mt = mt.annotate_rows(
        new_locus=new_locus_alleles.locus, new_alleles=new_locus_alleles.alleles
    )

    mt = split_shuffle(mt)
    mt = mt.annotate_entries(LA=mt.LA.map(lambda x: mt._old_to_new[x]))

    return mt.drop(
        "_keep_allele", "_new_to_old", "_old_to_new", "new_locus", "new_alleles"
    )


def create_full_subset_dense_mt(
    mt: hl.MatrixTable, meta_ht: hl.Table, variant_annotation_ht: hl.Table
):
    """
    Create the subset dense release MatrixTable with multi-allelic variants and all sample and variant annotations.

    .. note::

        This function uses the sparse subset MT and filters out LowQual variants and centromeres and telomeres.

    :param mt: Sparse subset release MatrixTable
    :param meta_ht: Metadata HT to use for sample (column) annotations
    :param variant_annotation_ht: Metadata HT to use for variant (row) annotations
    :return: Dense release MatrixTable with all row, column, and global annotations
    """
    logger.info(
        "Adding subset's sample QC metadata to MT columns and global annotations to MT globals..."
    )
    mt = mt.annotate_cols(**meta_ht[mt.col_key])
    mt = mt.annotate_globals(
        global_annotation_descriptions=convert_heterogeneous_dict_to_struct(
            GLOBAL_ANNOTATIONS
        ),
        **meta_ht.drop("global_annotation_descriptions").index_globals(),
    )

    logger.info(
        "Annotate entries with het non ref status for use in the homozygous alternate depletion fix..."
    )
    mt = mt.annotate_entries(_het_non_ref=mt.LGT.is_het_non_ref())

    logger.info("Splitting multi-allelics...")
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Computing adj and sex adjusted genotypes...")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus, mt.GT, mt.gnomad_sex_imputation.sex_karyotype
        ),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    logger.info(
        "Setting het genotypes at sites with > 1% AF (using precomputed v3.0 frequencies) and > 0.9 AB to homalt..."
    )
    # NOTE: Using v3.0 frequencies here and not v3.1 frequencies because
    # the frequency code adjusted genotypes (homalt depletion fix) using v3.0 frequencies
    # https://github.com/broadinstitute/gnomad_qc/blob/efea6851a421f4bc66b73db588c0eeeb7cd27539/gnomad_qc/v3/annotations/generate_freq_data_hgdp_tgp.py#L129
    freq_ht = release_sites(public=True).versions["3.0"].ht().select("freq")
    mt = hom_alt_depletion_fix(
        mt, het_non_ref_expr=mt._het_non_ref, af_expr=freq_ht[mt.row_key].freq[0].AF
    )
    mt = mt.drop("_het_non_ref")

    logger.info("Add all variant annotations and variant global annotations...")
    mt = mt.annotate_rows(**variant_annotation_ht[mt.row_key])
    mt = mt.annotate_globals(
        **variant_annotation_ht.drop("global_annotation_descriptions").index_globals()
    )

    logger.info("Removing chrM...")
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chrM")], keep=False)

    logger.info("Densify MT...")
    mt = hl.experimental.densify(mt)

    logger.info(
        "Filter out LowQual variants (using allele-specific annotation) and variants within centromere and telomere "
        "regions..."
    )
    mt = mt.filter_rows(
        ~mt.AS_lowqual & ~mt.telomere_or_centromere & (hl.len(mt.alleles) > 1)
    )
    mt = mt.drop("AS_lowqual", "telomere_or_centromere")

    return mt


def main(args):
    hl.init(log="/hgdp_tgp_subset.log", default_reference="GRCh38")

    test = args.test
    sample_annotation_resource = hgdp_tgp_subset_annotations(test=test)
    variant_annotation_resource = hgdp_tgp_subset_annotations(sample=False, test=test)
    sparse_mt_resource = hgdp_tgp_subset(test=test)
    dense_mt_resource = hgdp_tgp_subset(dense=True, test=test)

    if args.create_sample_annotation_ht:
        meta_ht = prepare_sample_annotations()
        meta_ht = meta_ht.checkpoint(
            sample_annotation_resource.path, overwrite=args.overwrite
        )

        meta_ht.relatedness_inference.summarize()
        logger.info("Number of samples in meta HT: %d", meta_ht.count())
        meta_ht = meta_ht.filter(
            meta_ht.high_quality & ~meta_ht.relatedness_inference.related
        )
        logger.info(
            "Number of high quality unrelated samples in meta HT: %d", meta_ht.count()
        )

    if (
        test
        and (
            args.export_sample_annotation_tsv
            or args.create_subset_sparse_mt
            or args.create_subset_dense_mt
        )
        and not file_exists(sample_annotation_resource.path)
    ):
        raise DataException(
            "There is currently no sample meta HT for the HGDP + TGP subset written to temp for testing. "
            "Run '--create_sample_annotation_ht' with '--test' to create one."
        )

    if args.export_sample_annotation_tsv:
        meta_ht = sample_annotation_resource.ht()
        meta_ht.export(hgdp_tgp_subset_sample_tsv(test=test))

    if args.create_subset_sparse_mt:
        # NOTE: We no longer remove samples that fail the gnomAD-wide sample QC high_quality filter. However, for frequency calculations we still remove samples failing hard filters and subcontinental PCA outliers.
        mt = get_gnomad_v3_mt(
            key_by_locus_and_alleles=True, remove_hard_filtered_samples=False
        )
        if test:
            logger.info(
                "Filtering MT to first %d partitions for testing",
                args.test_n_partitions,
            )
            mt = mt._filter_partitions(range(args.test_n_partitions))

        logger.info(
            "Filtering MT columns to HGDP + TGP samples and the CHMI haploid sample (syndip)"
        )
        # Note: Need to use sample names with the v3.1:: prefix
        meta_ht = meta.ht()
        meta_ht = meta_ht.filter(
            (meta_ht.subsets.hgdp | meta_ht.subsets.tgp | (meta_ht.s == SYNDIP))
        )
        mt = mt.filter_cols(hl.is_defined(meta_ht[mt.col_key]))
        logger.info("Number of samples in sparse MT: %d", mt.count_cols())

        logger.info(
            "Removing 'v3.1::' from the column names, these were added because there are duplicates of some 1KG/TGP samples"
            " in the full gnomAD dataset..."
        )
        mt = mt.key_cols_by(s=mt.s.replace("v3.1::", ""))

        # Adjust alleles and LA to include only alleles present in the subset
        mt = adjust_subset_alleles(mt)

        logger.info(
            "NOTE: This sparse MT is the raw subset MT and therefore it does not have adjusted sex genotypes and "
            "the fix for older GATK gVCFs with a known depletion of homozygous alternate alleles. This also does not "
            "remove standard GATK LowQual variants and variants in centromeres and telomeres (to preserve the ref "
            "block END annotation), which we recommend ultimately removing (and is removed for the dense MT release)."
        )
        mt.write(sparse_mt_resource.path, overwrite=args.overwrite)

    if (
        test
        and (args.create_variant_annotation_ht or args.create_subset_dense_mt)
        and not file_exists(sparse_mt_resource.path)
    ):
        raise DataException(
            "There is currently no sparse test MT for the HGDP + TGP subset. Run '--create_subset_sparse_mt' "
            "with '--test' to create one."
        )

    if args.create_variant_annotation_ht:
        logger.info("Creating variant annotation Hail Table")
        ht = sparse_mt_resource.mt().rows().select().select_globals()

        logger.info("Splitting multi-allelics and filtering out ref block variants...")
        ht = hl.split_multi(ht)
        ht = ht.filter(hl.len(ht.alleles) > 1)

        ht = prepare_variant_annotations(
            ht, filter_lowqual=False, vep_version=args.vep_version
        )
        ht.write(variant_annotation_resource.path, overwrite=args.overwrite)

    if args.create_subset_dense_mt:
        meta_ht = sample_annotation_resource.ht()
        if test and not file_exists(variant_annotation_resource.path):
            raise DataException(
                "There is currently no variant annotation HT for the HGDP + TGP subset written to temp for testing. "
                "Run '--create_variant_annotation_ht' with '--test' to create one."
            )
        variant_annotation_ht = variant_annotation_resource.ht()

        mt = sparse_mt_resource.mt()
        mt = mt.select_entries(*SPARSE_ENTRIES)
        mt = create_full_subset_dense_mt(mt, meta_ht, variant_annotation_ht)

        logger.info(
            "Writing dense HGDP + TGP MT with all sample and variant annotations"
        )
        mt.write(dense_mt_resource.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script subsets the gnomAD v3.1 release to only HGDP and 1KG/TGP samples."
    )
    parser.add_argument(
        "--create_sample_annotation_ht",
        help="Create the HGDP + 1KG/TGP subset sample metadata Hail Table.",
        action="store_true",
    )
    parser.add_argument(
        "--export_sample_annotation_tsv",
        help="Pull sample subset metadata and export to a .tsv.",
        action="store_true",
    )
    parser.add_argument(
        "--create_subset_sparse_mt",
        help="Create the HGDP + 1KG/TGP subset sparse MT. NOTE: This needs to be run without preemptibles because the allele adjustment requires a shuffle!",
        action="store_true",
    )
    parser.add_argument(
        "--create_variant_annotation_ht",
        help="Create the HGDP + 1KG/TGP subset variant annotation Hail Table.",
        action="store_true",
    )
    parser.add_argument(
        "--vep_version",
        help="Version of VEPed context Table to use in vep_or_lookup_vep",
        action="store_true",
        default="101",
    )
    parser.add_argument(
        "--create_subset_dense_mt",
        help="Create the HGDP + 1KG/TGP subset dense MT.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            "Run small test export on a subset of partitions of the MT. Writes to temp rather than writing to the "
            "main bucket."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test_n_partitions",
        default=5,
        type=int,
        help="Number of partitions to use for testing.",
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
