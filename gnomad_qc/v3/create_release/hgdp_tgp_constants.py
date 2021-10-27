GLOBAL_SAMPLE_ANNOTATIONS = {
    "gnomad_sex_imputation_ploidy_cutoffs": {
        "Description": (
            "Contains sex chromosome ploidy cutoffs used when determining sex chromosome karyotypes for the gnomAD sex "
            "imputation. Format: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), "
            "lower cutoff for triple X) and (lower cutoff for single Y, upper cutoff for single Y), lower cutoff for "
            "double Y)."
        )
    },
    "gnomad_population_inference_pca_metrics": {
        "Description": (
            "Contains the number of principal components (PCs) used when running PC-project and the minimum cutoff "
            "probability of belonging to a given population for the gnomAD population inference."
        )
    },
    "sample_hard_filter_cutoffs": {
        "Description": (
            "Contains the cutoffs used for hard-filtering samples prior to sample QC. Sample QC metrics are computed "
            "using the Hail sample_qc module on all autosomal bi-allelic SNVs. Samples are removed if they are clear "
            "outliers for any of the following metrics: number of snps (n_snp), ratio of heterozygous variants to "
            "homozygous variants (r_het_hom_var), number of singletons (n_singleton), and mean coverage on chromosome "
            "20 (cov). Additionally, we filter based on outliers of the following BAM/CRAM-derived metrics: "
            "% contamination (freemix), % chimera, and median insert size."
        )
    },
    "gnomad_sample_qc_metric_outlier_cutoffs": {
        "Description": (
            "Contains the cutoffs used for filtering outlier samples based on QC metrics (reported in the sample_qc "
            "and gnomad_sample_qc_residuals annotations). The first eight PCs computed during the gnomAD ancestry "
            "assignment were regressed out and the sample filter cutoffs were determined based on the residuals for "
            "each of the sample QC metrics. Samples were filtered if they fell outside four median absolute deviations "
            "(MADs) from the median for the following sample QC metrics: n_snp, r_ti_tv, r_insertion_deletion, "
            "n_insertion, n_deletion, n_het, n_hom_var, n_transition, and n_transversion. Samples over 8 MADs above "
            "the median n_singleton metric and over 4 MADs above the median r_het_hom_var metric were also filtered."
        )
    },
    "gnomad_age_distribution": {
        "Description": "gnomAD callset-wide age histogram calculated on release samples.",
        "sub_globals": {
            "bin_edges": {"Description": "Bin edges for the age histogram."},
            "bin_freq": {
                "Description": "Bin frequencies for the age histogram. This is the number of records found in each bin."
            },
            "n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge."
            },
            "n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge."
            },
        },
    },
}
GLOBAL_VARIANT_ANNOTATIONS = {
    "hgdp_tgp_freq_meta": {
        "Description": (
            "HGDP and 1KG frequency metadata. An ordered list containing the frequency aggregation group for each "
            "element of the hgdp_tgp_freq array row annotation."
        )
    },
    "gnomad_freq_meta": {
        "Description": (
            "gnomAD frequency metadata. An ordered list containing the frequency aggregation group for each element of "
            "the gnomad_freq array row annotation."
        )
    },
    "hgdp_tgp_freq_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: HGDP or 1KG subpopulation, "
            "sex: sex karyotype), with values describing the corresponding index of each grouping entry in the "
            "HGDP + 1KG frequency array annotation."
        )
    },
    "gnomad_freq_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD inferred global "
            "population sex: sex karyotype), with values describing the corresponding index of each grouping entry in "
            "the gnomAD frequency array annotation."
        )
    },
    "gnomad_faf_meta": {
        "Description": (
            "gnomAD filtering allele frequency metadata. An ordered list containing the frequency aggregation group "
            "for each element of the gnomad_faf array row annotation."
        )
    },
    "gnomad_faf_index_dict": {
        "Description": (
            "Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD inferred global "
            "population sex: sex karyotype), with values describing the corresponding index of each grouping entry in "
            "the filtering allele frequency (using Poisson 99% CI) annotation."
        )
    },
    "variant_filtering_model": {
        "Description": {"The variant filtering model used and its specific cutoffs."},
        "sub_globals": {
            "model_name": {
                "Description": (
                    "Variant filtering model name used in the 'filters' row annotation to indicate the variant was "
                    "filtered by the model during variant QC."
                )
            },
            "score_name": {"Description": "Name of score used for variant filtering."},
            "snv_cutoff": {
                "Description": "SNV filtering cutoff information.",
                "sub_globals": {
                    "bin": {"Description": "Filtering percentile cutoff for SNVs."},
                    "min_score": {
                        "Description": "Minimum score at SNV filtering percentile cutoff."
                    },
                },
            },
            "indel_cutoff": {
                "Description": "Information about cutoff used for indel filtering.",
                "sub_globals": {
                    "bin": {"Description": "Filtering percentile cutoff for indels."},
                    "min_score": {
                        "Description": "Minimum score at indel filtering percentile cutoff."
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
    "variant_inbreeding_coeff_cutoff": {
        "Description": "Hard-filter cutoff for InbreedingCoeff on variants."
    },
    "vep_version": {"Description": "VEP version."},
    "vep_csq_header": {"Description": "VEP header for VCF export."},
    "dbsnp_version": {"Description": "dbSNP version."},
}
GLOBAL_ANNOTATIONS = {**GLOBAL_SAMPLE_ANNOTATIONS, **GLOBAL_VARIANT_ANNOTATIONS}

SAMPLE_ANNOTATIONS = {
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
                "Description": (
                    "The mean coverage in bases of the genome territory after all filters are applied; see: "
                    "https://broadinstitute.github.io/picard/picard-metric-definitions.html."
                )
            },
            "median_coverage": {
                "Description": (
                    "The median coverage in bases of the genome territory after all filters are applied; see: "
                    "https://broadinstitute.github.io/picard/picard-metric-definitions.html."
                )
            },
            "mean_insert_size": {
                "Description": (
                    "The mean insert size of the 'core' of the distribution. Artefactual outliers in the distribution "
                    "often cause calculation of nonsensical mean and stdev values. To avoid this, the distribution is "
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
    "sample_qc": {
        "Description": "Struct containing sample QC metrics calculated using hl.sample_qc().",
        "sub_annotations": {
            "n_deletion": {"Description": "Number of deletion alternate alleles."},
            "n_het": {"Description": "Number of heterozygous calls."},
            "n_hom_ref": {"Description": "Number of homozygous reference calls."},
            "n_hom_var": {"Description": "Number of homozygous alternate calls."},
            "n_insertion": {"Description": "Number of insertion alternate alleles."},
            "n_non_ref": {"Description": "Sum of n_het and n_hom_var."},
            "n_snp": {"Description": "Number of SNP alternate alleles."},
            "n_transition": {
                "Description": "Number of transition (A-G, C-T) alternate alleles."
            },
            "n_transversion": {
                "Description": "Number of transversion alternate alleles."
            },
            "r_het_hom_var": {"Description": "Het/HomVar call ratio."},
            "r_insertion_deletion": {"Description": "Insertion/Deletion allele ratio."},
            "r_ti_tv": {"Description": "Transition/Transversion ratio."},
        },
    },
    "gnomad_sex_imputation": {
        "Description": "Struct containing sex imputation information.",
        "sub_annotations": {
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
            "f_stat": {
                "Description": "Inbreeding coefficient (excess heterozygosity) on chromosome X."
            },
            "n_called": {"Description": "Number of variants with a genotype call."},
            "expected_homs": {"Description": "Expected number of homozygotes."},
            "observed_homs": {"Description": "Observed number of homozygotes."},
        },
    },
    "gnomad_population_inference": {
        "Description": (
            "Struct containing ancestry information assigned by applying a principal components analysis (PCA) on "
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
    "gnomad_sample_qc_residuals": {
        "Description": (
            "Struct containing the residuals after regressing out the first eight PCs computed during the gnomAD "
            "ancestry assignment from each sample QC metric calculated using hl.sample_qc().",
        ),
        "sub_annotations": {
            "n_snp_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the number of SNP alternate "
                    "alleles."
                ),
            },
            "r_ti_tv_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the Transition/Transversion ratio."
                )
            },
            "r_insertion_deletion_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the Insertion/Deletion allele ratio."
                )
            },
            "n_insertion_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the number of insertion alternate "
                    "alleles."
                ),
            },
            "n_deletion_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the number of deletion alternate "
                    "alleles."
                )
            },
            "r_het_hom_var_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the Het/HomVar call ratio."
                ),
            },
            "n_transition_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the number of transition "
                    "(A-G, C-T) alternate alleles."
                )
            },
            "n_transversion_residual": {
                "Description": (
                    "Residuals after regressing out the first eight ancestry PCs from the number of transversion "
                    "alternate alleles."
                )
            },
        },
    },
    "gnomad_sample_filters": {
        "Description": "Sample QC filter annotations used for the gnomAD release.",
        "sub_annotations": {
            "hard_filters": {
                "Description": (
                    "Set of hard filters applied to each sample prior to additional sample QC. Samples are hard "
                    "filtered if they are extreme outliers for any of the following metrics: number of snps (n_snp), "
                    "ratio of heterozygous variants to homozygous variants (r_het_hom_var), number of singletons "
                    "(n_singleton), and mean coverage on chromosome 20 (cov). Additionally, we filter based on outliers "
                    "of the following Picard metrics: % contamination (freemix), % chimera, and median insert size."
                )
            },
            "hard_filtered": {
                "Description": (
                    "Whether a sample was hard filtered. The gnomad_sample_filters.hard_filters set is empty if this "
                    "annotation is True."
                )
            },
            "release_related": {
                "Description": (
                    "Whether a sample had a second-degree or greater relatedness to another sample in the gnomAD release."
                )
            },
            "qc_metrics_filters": {
                "Description": (
                    "Set of all sample QC metrics for which each sample was found to be an outlier after computing "
                    "sample QC metrics using the Hail sample_qc() module and regressing out the first 8 gnomAD ancestry "
                    "assignment PCs."
                )
            },
        },
    },
    "gnomad_high_quality": {
        "Description": (
            "Whether a sample has passed gnomAD sample QC metrics except for relatedness "
            "(i.e., gnomad_sample_filters.hard_filters and gnomad_sample_filters.qc_metrics_filters are empty sets)."
        )
    },
    "gnomad_release": {
        "Description": (
            "Whether the sample was included in the gnomAD release dataset. For the full gnomAD release, relatedness "
            "inference is performed on the full dataset, and release samples are chosen in a way that maximizes the "
            "number of samples retained while filtering the dataset to include only samples with less than "
            "second-degree relatedness. For the HGDP + 1KG subset, samples passing all other sample QC metrics "
            "are retained."
        )
    },
    "relatedness_inference": {
        "Description": "Information about the sample’s relatedness to other samples within the callset.",
        "sub_annotations": {
            "related_samples": {
                "Description": (
                    "Set of all HGDP or 1KG samples that have a kinship estimate (kin) > 0.05 determined using Hail’s "
                    "pc_relate module (https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) "
                    "with this sample. More details on the relatedness inference can be found here: "
                    "https://github.com/atgu/hgdp_tgp/blob/master/pca_subcont.ipynb. This set is empty if the sample "
                    "has no such relationships within this HGDP + 1KG subset. Each entry of the set consists of a "
                    "struct containing the following information about this relationship:"
                ),
                "sub_annotations": {
                    "s": {"Description": "Sample ID."},
                    "kin": {"Description": "Kinship estimate."},
                    "ibd0": {"Description": "IBD0 estimate."},
                    "ibd1": {"Description": "IBD1 estimate."},
                    "ibd2": {"Description": "IBD2 estimate."},
                },
            },
            "related": {
                "Description": (
                    "Indicates whether this sample is excluded from variant frequency calculations (hgdp_tgp_freq) "
                    "because of relatedness to other samples. Closely related individuals are pruned from the dataset "
                    "to enhance the accuracy of population frequency estimates. Pruning using Hail’s "
                    "maximal_independent_set module (https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set) "
                    "maintains the maximal number of individuals in the dataset."
                )
            },
        },
    },
    "hgdp_tgp_meta": {
        "Description": "Sample metadata specific to the HGDP + 1KG subset.",
        "sub_annotations": {
            "project": {
                "Description": (
                    "Indicates if the sample is part of the Human Genome Diversity Project (‘HGDP’), the ‘1000 Genomes’ "
                    "project, or is the synthetic-diploid sample (a mixture of DNA from two haploid CHM cell lines; "
                    "https://www.nature.com/articles/s41592-018-0054-7?WT.feed_name=subjects_standards; "
                    "https://github.com/lh3/CHM-eval)."
                )
            },
            "study_region": {"Description": "Study-specific global population labels."},
            "population": {
                "Description": (
                    "Population label for HGDP or 1KG. The HGDP populations are detailed in "
                    "https://science.sciencemag.org/content/367/6484/eaay5012. 1KG populations are described here: "
                    "https://www.internationalgenome.org/category/population."
                )
            },
            "genetic_region": {
                "Description": "Global population labels harmonized across both studies."
            },
            "latitude": {
                "Description": "Approximate latitude of the geographical place of origin of the population."
            },
            "longitude": {
                "Description": "Approximate longitude of the geographical place of origin of the population."
            },
            "hgdp_technical_meta": {
                "Description": (
                    "Technical considerations for HGDP detailed in https://science.sciencemag.org/content/367/6484/eaay5012. "
                    "This struct will be missing for 1KG samples and syndip."
                ),
                "sub_annotations": {
                    "source": {
                        "Description": (
                            "Which batch/project these HGDP samples were sequenced in (Sanger vs Simons Genome "
                            "Diversity Project)."
                        )
                    },
                    "library_type": {
                        "Description": "Whether library prep was PCR-free or used PCR."
                    },
                },
            },
            "global_pca_scores": {
                "Description": (
                    "Array of the first 20 principal components analysis (PCA) scores on the full HGDP + 1KG subset. "
                    "Obtained by first dividing the sample set into relateds and unrelateds (using "
                    "relatedness_inference.related). PCA was run on the unrelated samples using Hail’s "
                    "hwe_normalized_pca module (https://hail.is/docs/0.2/methods/genetics.html#hail.methods.hwe_normalized_pca), "
                    "producing a global PCA score for each of the unrelated samples and loadings for each variant. "
                    "The related samples were then projected (https://hail.is/docs/0.2/experimental/index.html#hail.experimental.pc_project) "
                    "onto the predefined PCA space using the variant loadings from the unrelated sample PCA to produce "
                    "the PC scores for the related samples. Code used to obtain these scores can be found under "
                    "'global pca' here: https://github.com/atgu/hgdp_tgp/blob/master/pca_subcont.ipynb"
                )
            },
            "subcontinental_pca": {
                "Description": (
                    "The subcontinental PCAs were obtained in a similar manner as the global PCA scores "
                    "(hgdp_tgp_meta.global_pca_scores). The full HGDP + 1KG subset was split by genetic region "
                    "(hgdp_tgp_meta.genetic_region) prior to performing the same PCA and PC projection steps "
                    "described for the global PCA scores. The code used to obtain these scores can be found under "
                    "subcontinental pca here: https://github.com/atgu/hgdp_tgp/blob/master/pca_subcont.ipynb "
                ),
                "sub_annotations": {
                    "pca_scores": {
                        "Description": (
                            "Array of the first 20 subcontinental PCA scores for the sample based on its value in the "
                            "hgdp_tgp_meta.genetic_region annotation (one of: AFR, AMR, CSA, EAS, EUR, MID, or OCE)."
                        )
                    },
                    "pca_scores_outliers_removed": {
                        "Description": (
                            "Array of the first 20 subcontinental PCA scores following the removal of samples labeled "
                            "as subcontinental PCA outliers (hgdp_tgp_meta.subcontinental_pca.outlier)."
                        )
                    },
                    "outlier": {
                        "Description": "Whether the sample was an outlier in the subcontinental PCAs."
                    },
                },
            },
            "gnomad_labeled_subpop": {
                "Description": (
                    "Similar to the 'hgdp_tgp_meta.population' annotation, this is the sample's population label "
                    "supplied by HGDP or 1KG with slight modifications that were used to harmonize labels with other "
                    "gnomAD samples."
                )
            },
        },
    },
    "high_quality": {
        "Description": (
            "Samples that pass all ‘gnomad_sample_filters.hard_filters’ and were not found to be outliers in global "
            "population-specific principal component analysis hgdp_tgp_meta.subcontinental_pca.outlier"
        )
    },
}

VARIANT_ANNOTATIONS = {
    "locus": {
        "Description": "Variant locus. Contains contig and position information.",
    },
    "alleles": {"Description": "Variant alleles.",},
    "rsid": {"Description": "dbSNP reference SNP identification (rsID) numbers.",},
    "a_index": {
        "Description": (
            "The original index of this alternate allele in the multiallelic representation (1 is the first alternate "
            "allele or the only alternate allele in a biallelic variant)."
        )
    },
    "was_split": {
        "Description": "True if this variant was originally multiallelic, otherwise False."
    },
    "hgdp_tgp_freq": {
        "Description": (
            "Allele frequency information (AC, AN, AF, homozygote count) in HGDP + 1KG samples that pass the "
            "high_quality sample annotation and are inferred as unrelated (False in relatedness_inference.related "
            "annotation)."
        ),
        "sub_annotations": {
            "AC": {
                "Description": "Alternate allele count  in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
            "AF": {
                "Description": "Alternate allele frequency  in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
            "AN": {
                "Description": "Total number of alleles in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
            "homozygote_count": {
                "Description": "Count of homozygous individuals in HGDP + 1KG samples that pass the high_quality sample annotation."
            },
        },
    },
    "gnomad_freq": {
        "Description": "Allele frequency information (AC, AN, AF, homozygote count) in gnomAD release.",
        "sub_annotations": {
            "AC": {"Description": "Alternate allele count in gnomAD release."},
            "AF": {"Description": "Alternate allele frequency in gnomAD release."},
            "AN": {"Description": "Total number of alleles in gnomAD release."},
            "homozygote_count": {
                "Description": "Count of homozygous individuals in gnomAD release."
            },
        },
    },
    "gnomad_popmax": {
        "Description": "Allele frequency information (AC, AN, AF, homozygote count) for the population with maximum AF in gnomAD.",
        "sub_annotations": {
            "AC": {
                "Description": "Allele count in the population with the maximum AF in gnomAD."
            },
            "AF": {
                "Description": "Maximum allele frequency across populations in gnomAD."
            },
            "AN": {
                "Description": "Total number of alleles in the population with the maximum AF in gnomAD."
            },
            "homozygote_count": {
                "Description": "Count of homozygous individuals in the population with the maximum allele frequency in gnomAD."
            },
            "pop": {"Description": "Population with maximum AF in gnomAD."},
            "faf95": {
                "Description": (
                    "Filtering allele frequency (using Poisson 95% CI) for the population with the maximum allele "
                    "frequency in gnomAD."
                )
            },
        },
    },
    "gnomad_faf": {
        "Description": "Filtering allele frequency in gnomAD release.",
        "sub_annotations": {
            "faf95": {
                "Description": "Filtering allele frequency in gnomAD release (using Poisson 95% CI)."
            },
            "faf99": {
                "Description": "Filtering allele frequency in gnomAD release (using Poisson 99% CI)."
            },
        },
    },
    "gnomad_qual_hists": {
        "Description": "gnomAD genotype quality metric histograms for high quality genotypes.",
        "sub_annotations": {
            "gq_hist_all": {
                "Description": "Histogram for GQ calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the GQ histogram calculated on high quality genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                        ),
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the GQ histogram calculated on high quality genotypes. The number of "
                            "records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ calculated on high "
                            "quality genotypes"
                        ),
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ calculated on high "
                            "quality genotypes"
                        ),
                    },
                },
            },
            "dp_hist_all": {
                "Description": "Histogram for DP calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the DP histogram calculated on high quality genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the DP histogram calculated on high quality genotypes. The number of "
                            "records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of DP values falling below lowest histogram bin edge, for DP calculated on high "
                            "quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of DP values falling above highest histogram bin edge, for DP calculated on high "
                            "quality genotypes."
                        )
                    },
                },
            },
            "gq_hist_alt": {
                "Description": "Histogram for GQ in heterozygous individuals calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of GQ in heterozygous individuals calculated on high quality "
                            "genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        ),
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of GQ in heterozygous individuals calculated on high "
                            "quality genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                },
            },
            "dp_hist_alt": {
                "Description": "Histogram for DP in heterozygous individuals calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of DP in heterozygous individuals calculated on high quality "
                            "genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of DP in heterozygous individuals calculated on high "
                            "quality genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of DP values falling below lowest histogram bin edge, for DP in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of DP values falling above highest histogram bin edge, for DP in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                },
            },
            "ab_hist_alt": {
                "Description": "Histogram for AB in heterozygous individuals calculated on high quality genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of AB in heterozygous individuals calculated on high quality "
                            "genotypes are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of AB in heterozygous individuals calculated on high "
                            "quality genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of AB values falling below lowest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of AB values falling above highest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on high quality genotypes."
                        )
                    },
                },
            },
        },
    },
    "gnomad_raw_qual_hists": {
        "Description": "gnomAD genotype quality metric histograms.",
        "sub_annotations": {
            "gq_hist_all": {
                "Description": "Histogram for GQ calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the GQ histogram calculated on all genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the GQ histogram calculated on all genotypes. The number of records "
                            "found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ calculated on all "
                            "genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ calculated on all "
                            "genotypes."
                        )
                    },
                },
            },
            "dp_hist_all": {
                "Description": "Histogram for DP calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the DP histogram calculated on all genotypes are: "
                            "0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100"
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the DP histogram calculated on all genotypes. The number of records "
                            "found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of DP values falling below lowest histogram bin edge, for DP calculated on all "
                            "genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of DP values falling above highest histogram bin edge, for DP calculated on all "
                            "genotypes."
                        )
                    },
                },
            },
            "gq_hist_alt": {
                "Description": "Histogram for GQ in heterozygous individuals calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of GQ in heterozygous individuals calculated on all genotypes "
                            "are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of GQ in heterozygous individuals calculated on all "
                            "genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of GQ values falling below lowest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of GQ values falling above highest histogram bin edge, for GQ in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                },
            },
            "dp_hist_alt": {
                "Description": "Histogram for DP in heterozygous individuals calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of DP in heterozygous individuals calculated on all genotypes "
                            "are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of DP in heterozygous individuals calculated on all "
                            "genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": "Count of DP values falling below lowest histogram bin edge, for DP in heterozygous individuals calculated on all genotypes."
                    },
                    "n_larger": {
                        "Description": "Count of DP values falling above highest histogram bin edge, for DP in heterozygous individuals calculated on all genotypes."
                    },
                },
            },
            "ab_hist_alt": {
                "Description": "Histogram for AB in heterozygous individuals calculated on all genotypes.",
                "sub_annotations": {
                    "bin_edges": {
                        "Description": (
                            "Bin edges for the histogram of AB in heterozygous individuals calculated on all genotypes "
                            "are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00."
                        )
                    },
                    "bin_freq": {
                        "Description": (
                            "Bin frequencies for the histogram of AB in heterozygous individuals calculated on all "
                            "genotypes. The number of records found in each bin."
                        )
                    },
                    "n_smaller": {
                        "Description": (
                            "Count of AB values falling below lowest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                    "n_larger": {
                        "Description": (
                            "Count of AB values falling above highest histogram bin edge, for AB in heterozygous "
                            "individuals calculated on all genotypes."
                        )
                    },
                },
            },
        },
    },
    "gnomad_age_hist_het": {
        "Description": "Histogram for age in all heterozygous gnomAD release samples calculated on high quality genotypes.",
        "sub_annotations": {
            "bin_edges": {"Description": "Bin edges for the age histogram."},
            "bin_freq": {
                "Description": "Bin frequencies for the age histogram. This is the number of records found in each bin."
            },
            "n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge."
            },
            "n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge."
            },
        },
    },
    "gnomad_age_hist_hom": {
        "Description": "Histogram for age in all homozygous gnomAD release samples calculated on high quality genotypes.",
        "sub_annotations": {
            "bin_edges": {"Description": "Bin edges for the age histogram."},
            "bin_freq": {
                "Description": "Bin frequencies for the age histogram. This is the number of records found in each bin."
            },
            "n_smaller": {
                "Description": "Count of age values falling below lowest histogram bin edge."
            },
            "n_larger": {
                "Description": "Count of age values falling above highest histogram bin edge."
            },
        },
    },
    "filters": {
        "Description": (
            "Variant filters; AC0: Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10;"
            " and AB < 0.2 for het calls), AS_VQSR: Failed VQSR filtering thresholds of -2.7739 for SNPs and -1.0606 "
            "for indels, InbreedingCoeff: GATK InbreedingCoeff < -0.3, PASS: Passed all variant filters."
        )
    },
    "info": {
        "Description": "Struct containing typical GATK allele-specific (AS) info fields and additional variant QC fields.",
        "sub_annotations": {
            "QUALapprox": {
                "Description": "Sum of PL[0] values; used to approximate the QUAL score."
            },
            "SB": {
                "Description": (
                    "Per-sample component statistics which comprise the Fisher's exact test to detect strand bias. "
                    "Values are: depth of reference allele on forward strand, depth of reference allele on reverse "
                    "strand, depth of alternate allele on forward strand, depth of alternate allele on reverse strand."
                )
            },
            "MQ": {
                "Description": "Root mean square of the mapping quality of reads across all samples."
            },
            "MQRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities."
            },
            "VarDP": {
                "Description": "Depth over variant genotypes (does not include depth of reference samples)."
            },
            "AS_ReadPosRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read position bias."
            },
            "AS_pab_max": {
                "Description": (
                    "Maximum p-value over callset for binomial test of observed allele balance for a heterozygous "
                    "genotype, given expectation of 0.5."
                )
            },
            "AS_QD": {
                "Description": "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant."
            },
            "AS_MQ": {
                "Description": "Allele-specific root mean square of the mapping quality of reads across all samples."
            },
            "QD": {
                "Description": "Variant call confidence normalized by depth of sample reads supporting a variant."
            },
            "AS_MQRankSum": {
                "Description": "Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities."
            },
            "FS": {
                "Description": "Phred-scaled p-value of Fisher's exact test for strand bias."
            },
            "AS_FS": {
                "Description": "Allele-specific phred-scaled p-value of Fisher's exact test for strand bias."
            },
            "ReadPosRankSum": {
                "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias."
            },
            "AS_QUALapprox": {
                "Description": "Allele-specific sum of PL[0] values; used to approximate the QUAL score."
            },
            "AS_SB_TABLE": {
                "Description": "Allele-specific forward/reverse read counts for strand bias tests."
            },
            "AS_VarDP": {
                "Description": "Allele-specific depth over variant genotypes (does not include depth of reference samples)."
            },
            "AS_SOR": {
                "Description": "Allele-specific strand bias estimated by the symmetric odds ratio test."
            },
            "SOR": {
                "Description": "Strand bias estimated by the symmetric odds ratio test."
            },
            "transmitted_singleton": {
                "Description": (
                    "Variant was a callset-wide doubleton that was transmitted within a family from a parent to a "
                    "child (i.e., a singleton amongst unrelated samples in cohort)."
                )
            },
            "omni": {
                "Description": "Variant is present on the Omni 2.5 genotyping array and found in 1000 Genomes data."
            },
            "mills": {"Description": "Indel is present in the Mills and Devine data."},
            "monoallelic": {
                "Description": "All samples are all homozygous alternate for the variant."
            },
            "InbreedingCoeff": {
                "Description": (
                    "Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number "
                    "of heterozygous genotypes)/(the number of heterozygous genotypes expected under Hardy-Weinberg "
                    "equilibrium)."
                )
            },
        },
    },
    "vep": {
        "Description": (
            "Consequence annotations from Ensembl VEP. More details about VEP output is described here: "
            "https://uswest.ensembl.org/info/docs/tools/vep/vep_formats.html#output. VEP was run using the LOFTEE "
            "plugin and information about the additional LOFTEE annotations can be found here: "
            "https://github.com/konradjk/loftee."
        )
    },
    "vqsr": {
        "Description": "VQSR related variant annotations.",
        "sub_annotations": {
            "AS_VQSLOD": {
                "Description": (
                    "Allele-specific log-odds ratio of being a true variant versus being a false positive under the "
                    "trained VQSR Gaussian mixture model."
                )
            },
            "AS_culprit": {
                "Description": "Allele-specific worst-performing annotation in the VQSR Gaussian mixture model."
            },
            "NEGATIVE_TRAIN_SITE": {
                "Description": "Variant was used to build the negative training set of low-quality variants for VQSR."
            },
            "POSITIVE_TRAIN_SITE": {
                "Description": "Variant was used to build the positive training set of high-quality variants for VQSR."
            },
        },
    },
    "region_flag": {
        "Description": "Struct containing flags for problematic regions.",
        "sub_annotations": {
            "lcr": {"Description": "Variant falls within a low complexity region."},
            "segdup": {
                "Description": "Variant falls within a segmental duplication region."
            },
        },
    },
    "allele_info": {
        "Description": "Allele information.",
        "sub_annotations": {
            "variant_type": {
                "Description": "Variant type (snv, indel, multi-snv, multi-indel, or mixed).",
            },
            "allele_type": {
                "Description": "Allele type (snv, insertion, deletion, or mixed).",
            },
            "n_alt_alleles": {
                "Description": "Total number of alternate alleles observed at variant locus.",
            },
        },
    },
    "was_mixed": {"Description": "Variant type was mixed."},
    "cadd": {
        "sub_annotations": {
            "raw_score": {
                "Description": (
                    "Raw CADD scores are interpretable as the extent to which the annotation profile for a given "
                    "variant suggests that the variant is likely to be 'observed' (negative values) vs 'simulated' "
                    "(positive values); higher values indicate that a variant is more likely to be simulated (or 'not "
                    "observed') and therefore more likely to have deleterious effects. More information can be found "
                    "on the CADD website: https://cadd.gs.washington.edu/info."
                )
            },
            "phred": {
                "Description": (
                    "CADD Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on the rank of each "
                    "variant relative to all possible 8.6 billion substitutions in the human reference genome. Larger "
                    "values are more deleterious. More information can be found on the CADD website: "
                    "https://cadd.gs.washington.edu/info."
                )
            },
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one CADD score associated with "
                    "it. For a small set of variants, the in silico predictors calculated multiple scores per variant "
                    "based on additional information. For example, if a variant is found in multiple transcripts or if "
                    "it has multiple trinucleotide contexts, an in silico predictor may report scores for multiple "
                    "scenarios. The highest score was taken for each variant in the cases where the in silico predictor "
                    "calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "revel": {
        "Description": (
            "dbNSFP's Revel score, ranging from 0 to 1. Variants with higher scores are predicted to be more likely to "
            "be deleterious."
        ),
        "sub_annotations": {
            "revel_score": {"Description": "Revel’s numerical score from 0 to 1."},
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one revel_score associated "
                    "with it. For a small set of variants, the in silico predictors calculated multiple scores per "
                    "variant based on additional information. For example, if a variant is found in multiple "
                    "transcripts or if it has multiple trinucleotide contexts, an in silico predictor may report "
                    "scores for multiple scenarios. The highest score was taken for each variant in the cases where "
                    "the in silico predictor calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "splice_ai": {
        "sub_annotations": {
            "splice_ai": {
                "Description": "The maximum delta score, interpreted as the probability of the variant being splice-altering."
            },
            "splice_consequence": {
                "Description": "The consequence term associated with the max delta score in 'splice_ai’."
            },
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one splice_ai score associated "
                    "with it. For a small set of variants, the in silico predictors calculated multiple scores per "
                    "variant based on additional information. For example, if a variant is found in multiple "
                    "transcripts or if it has multiple trinucleotide contexts, an in silico predictor may report "
                    "scores for multiple scenarios. The highest score was taken for each variant in the cases where "
                    "the in silico predictor calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "primate_ai": {
        "sub_annotations": {
            "primate_ai_score": {
                "Description": "PrimateAI's deleteriousness score from 0 (less deleterious) to 1 (more deleterious)."
            },
            "has_duplicate": {
                "Description": (
                    "A True/False flag that indicates whether the variant has more than one primate_ai_score associated "
                    "with it. For a small set of variants, the in silico predictors calculated multiple scores per "
                    "variant based on additional information. For example, if a variant is found in multiple "
                    "transcripts or if it has multiple trinucleotide contexts, an in silico predictor may report "
                    "scores for multiple scenarios. The highest score was taken for each variant in the cases where "
                    "the in silico predictor calculates multiple scores, and we flag variants with multiple scores."
                )
            },
        },
    },
    "AS_lowqual": {
        "Description": (
            "Whether the variant falls below a low quality threshold and was excluded from the gnomAD dataset. We "
            "recommend filtering all such variants. This is similar to the GATK LowQual filter, but is allele-specific. "
            "GATK computes this annotation at the site level, which uses the least stringent prior for mixed sites."
        )
    },
    "telomere_or_centromere": {
        "Description": (
            "Whether the variant falls within a telomere or centromere region. These variants were excluded from the "
            "gnomAD dataset. We recommend filtering all such variants."
        )
    },
}
