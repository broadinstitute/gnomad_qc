## Sample QC annotations

Below are gnomAD sample QC metadata annotation definitions. These annotations are aggregated into a single metadata Hail Table using the `--generate_metadata` flag in the `sample_qc.py` script.

## Global annotations

* **sex_imputation_ploidy_cutoffs**: Contains sex chromosome ploidy cutoffs used when determining sex chromosome karyotypes for sex imputation.
    * **x_ploidy_cutoffs**: Cutoffs for X ploidy.
        * upper_cutoff_X: Upper cutoff for single X.
        * lower_cutoff_XX: Lower cutoff for double X.
        * upper_cutoff_XX: Upper cutoff for double X.
        * lower_cutoff_XXX: Lower cutoff for triple X. 
    * **y_ploidy_cutoffs**: Cutoffs for Y ploidy.
        * lower_cutoff_Y: Lower cutoff for single Y.
        * upper_cutoff_Y: Upper cutoff for single Y.
        * lower_cutoff_YY: Lower cutoff for double Y.
    * **f_stat_cutoff**: F-statistic cutoff to roughly divide 'XX' from 'XY' samples. XX samples are below cutoff and XY are above cutoff. Not used in final ploidy annotation.

 * **population_inference_pca_metrics**: Contains parameters from population assignment.
    * **min_prob**: The minimum cutoff probability of belonging to a given population at which the population was set.	
    * **include_unreleasable_samples**: Whether unreleasable samples were included in the ancestry PCA.
    * **max_mislabeled_training_samples**: The max number of mislabeled training samples that were allowed in population assignment.
    * **n_pcs**: The number of principal components (PCs) used in [PC-project]([https://github.com/broadinstitute/gnomad_methods/blob/0b660daf23854efa5f6baa233ecbedce40d06fd3/gnomad/sample_qc/ancestry.py#L80](https://github.com/broadinstitute/gnomad_methods/blob/0b660daf23854efa5f6baa233ecbedce40d06fd3/gnomad/sample_qc/ancestry.py#L80)).
    * **pop_assignment_iterations**: How many population assignment iterations ran.

* **relatedness_inference_cutoffs**: Contains relatedness inference cutoffs used when running Hail's [pc-relate]([https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate)).
    * **min_individual_maf**: The minimum individual-specific minor allele frequency used to estimate relatedness for a pair of samples.
    * **min_emission_kinship**: The minimum kinship coefficient cutoff in the results.
    * **ibd0_0_max**: The IBD0 cutoff used to determine parent-offspring vs full sibling relationships.
    * **second_degree_kin_cutoff**: The minimum kinship threshold used for filtering a pair of samples with a second degree relationship.
    * **first_degree_kin_thresholds**: The first degree kinship threshold used for filtering a pair of samples with a first degree relationship. 

* **outlier_detection_metrics**: Contains the linear regression statistics and cutoffs used for filtering outlier samples based on QC metrics. 
    * **lms**: Linear regression statistics for QC metrics.
        * **n_snp**: SNP regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **n_singleton**: Singleton regression statistics.
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **r_ti_tv**: Transition/Transversion ratio regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **r_insertion_deletion**: Insertion/Deletion ratio regression statistics.
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **n_insertion**: Insertion regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **n_deletion**: Deletion regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **r_het_hom_var**: Het/HomVar regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **n_transition**: Transition regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
        * **n_tranversion**: Transversion regression statistics
            * **beta**: Estimated regression coefficient for each PC.
            * **standard_error**: Estimated standard error for each PC.
            * **t_stat**: t-statistic for each PC.
            * **p_value**: p-value for each PC.
            * **multiple_standard_error**: Estimated standard deviation of the random error.
            * **multiple_r_squared**: Coefficient of determination for nested models.
            * **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
            * **f_stat**: F-statistic for nested models.
            * **multiple_p_value**: (p-value for the F-test of nested models.
            * **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
    * **qc_metric_stats**: The statistics generated for outlier filter cutoffs based on the residuals for each of the sample QC metrics.
        * **n_snp_residual**:
            * **median**: The n_snp median.
            * **mad**: The n_snp median absolute deviation (MAD).
            * **lower**: The n_snp lower median absolute deviation threshold.
            * **upper**: The n_snp upper median absolute deviation threshold.
        * **n_singleton_residual**:
            * **median**: The n_singleton median.
            * **mad**: The n_singleton median absolute deviation (MAD).
            * **lower**: The n_singleton lower median absolute deviation threshold.
            * **upper**: The n_singleton upper median absolute deviation threshold.
        * **r_ti_tv_residual**:
            * **median**: The r_ti_tv median.
            * **mad**: The r_ti_tv median absolute deviation (MAD).
            * **lower**: The r_ti_tv lower median absolute deviation threshold.
            * **upper**: The r_ti_tv upper median absolute deviation threshold.
        * **r_insertition_deletion_residual**:
            * **median**: The r_insertition_deletion median.
            * **mad**: The r_insertition_deletion median absolute deviation (MAD).
            * **lower**: The r_insertition_deletion lower median absolute deviation threshold.
            * **upper**: The r_insertition_deletion upper median absolute deviation threshold.
        * **n_insertion_residual**:
            * **median**: The n_insertion median.
            * **mad**: The n_insertion median absolute deviation (MAD).
            * **lower**: The n_insertion lower median absolute deviation threshold.
            * **upper**: The n_insertion upper median absolute deviation threshold.
        * **n_deletion_residual**:
            * **median**: The n_deletion median.
            * **mad**: The n_deletion median absolute deviation (MAD).
            * **lower**: The n_deletion lower median absolute deviation threshold.
            * **upper**: The n_deletion upper median absolute deviation threshold.
        * **r_het_hom_residual**:
            * **median**: The r_het_hom median.
            * **mad**: The r_het_hom median absolute deviation (MAD).
            * **lower**: The r_het_hom lower median absolute deviation threshold.
            * **upper**: The r_het_hom upper median absolute deviation threshold.
        * **n_transition_residual**:
            * **median**: The n_transition median.
            * **mad**: The n_transition median absolute deviation (MAD).
            * **lower**: The n_transition lower median absolute deviation threshold.
            * **upper**: The n_transition upper median absolute deviation threshold.
        * **n_transversion_residual**:
            * **median**: The n_transversion median.
            * **mad**: The n_transversion median absolute deviation (MAD).
            * **lower**: The n_transversion lower median absolute deviation threshold.
            * **upper**: The n_transversion upper median absolute deviation threshold.
    * **n_pcs**: The number of PCs computed during the ancestry assignment that were regressed out.
    * **used_regressed_metrics**: Whether the regression outlier Hail Table was used for outlier detection rather than the population stratified Hail Table.

* **hard_filter_cutoffs**: Contains the cutoffs used for hard-filtering samples prior to sample QC. 
    * **min_cov**: Filtering threshold to use for chr20 coverage.
    * **max_n_snp**: Filtering threshold to use for the max number of SNPs.
    * **min_n_snp**: Filtering threshold to use for the min number of SNPs.
    * **max_n_singleton**: Filtering threshold to use for the max number of singletons.
    * **max_r_het_hom_var**: Filtering threshold to use for the max ratio of heterozygotes to alternate homozygotes.
    * **max_pct_contamination**: Filtering threshold to use for max percent contamination (this is a percent not a proportion, e.g. 5% == 5.00, %5 != 0.05).
    * **max_pct_chimera**: Filtering threshold to use for max percent chimera (this is a percent not a proportion, e.g. 5% == 5.00, %5 != 0.05).
    * **min_median_insert_size**: Filtering threshold to use for min median insert size.

**Row annotations**:

* **s**: Unique sample ID
* **project_meta**: Metadata collected from collaborators or previous versions of gnomAD. The description of annotations in this struct can be found in the [gnomad_meta repo](https://github.com/broadinstitute/gnomad_meta/v3.1/README.md) (for internal use only).

* **subsets**: Subsets of the gnomAD release HT that the sample belongs to
    * **non_topmed**: Whether the sample was included in the ‘non_topmed’ subset. Subset contains variants from samples that are not the TOPMED cohort. Generated using the negation of the ‘project_meta.topmed’ annotation
    * **controls_and_biobanks**: Whether the sample was included in the ‘controls_and_biobanks’ subset. Subset contains variants from samples that are either a control or from a biobank project. Generated matching on ‘case_control’ annotation values of ‘control' or 'biobank'
    * **non_neuro**: Whether the sample was included in the ‘non_neuro’ subset. Subset contains variants from samples that are either 1. not cases in neuro cohorts or 2. neuro cohort samples missing ‘case_control’ information. Generated using the ‘neuro_case’ annotation which is calculated using’ neuro_cohort’ and ‘case_control’.
    * **non_v2**: Whether the sample was included in the ‘non_v2’ subset. Subset contains variants from samples that are not in gnomAD v2, exomes or genomes. Generated using ‘v2_release’ annotation. 
    * **non_cancer**: Whether the sample was included in the ‘non_cancer’ subset. Subset contains variants from samples that are not from TCGA, regardless of germline or somatic status. 
    * **tgp**: Whether the sample was included in the 1000 Genomes Project subset.
    * **hgdp**: Whether the sample was included in the Human Genome Diversity Project subset.

* **bam_metrics**: Sample level metrics obtained from BAMs/CRAMs.
    * **pct_bases_20x**: The fraction of bases that attained at least 20X sequence coverage in post-filtering bases.
    * **pct_chimeras**: The fraction of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes.
    * **freemix**: Estimate of contamination (0-100 scale).
    * **mean_coverage**: The mean coverage in bases of the genome territory after all filters are applied; see [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html).
    * **median_coverage**: The median coverage in bases of the genome territory after all filters are applied; see [here](https://broadinstitute.github.io/picard/picard-metric-definitions.html).
    * **mean_insert_size**: The mean insert size of the 'core' of the distribution. Artefactual outliers in the distribution often cause calculation of nonsensical mean and stdev values. To avoid this, the distribution is first trimmed to a 'core' distribution of +/- N median absolute deviations around the median insert size.
    * **median_insert_size**: The median insert size of all paired end reads where both ends mapped to the same chromosome.
    * **pct_bases_10x**: The fraction of bases that attained at least 10X sequence coverage in post-filtering bases.

* **sex_imputation**: Struct containing sex imputation information.
    * **is_female**: Whether the sample was imputed female.
    * **chr20_mean_dp**: Sample's mean depth across chromosome 20.
    * **chrX_mean_dp**: Sample's mean depth across chromosome X.
    * **chrY_mean_dp**: Sample's mean depth across chromosome Y.
    * **chrX_ploidy**: Sample's chromosome X ploidy (chrX_mean_dp normalized using chr20_mean_dp).
    * **chrY_ploidy**: Sample's chromosome Y ploidy (chrY_mean_dp normalized using chr20_mean_dp).
    * **X_karyotype**: Sample's chromosome X karyotype.
    * **Y_karyotype**: Sample's chromosome Y karyotype.
    * **sex_karyotype**: Sample's sex karyotype (combined X and Y karyotype).
    * **impute_sex_stats**:
        * **f_stat**: Inbreeding coefficient (excess heterozygosity) on chromosome X.
        * **n_called**: Number of variants with a genotype call.
        * **expected_homs**: Expected number of homozygotes.
        * **observed_homs**: Observed number of homozygotes.

* **sample_qc**: Struct containing sample QC metrics calculated using Hail’s  [sample_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc).
    * **n_hom_ref**: Number of homozygous reference calls.
    * **n_het**: Number of heterozygous calls.
    * **n_hom_var**: Number of homozygous alternate calls.
    * **n_non_ref**: Sum of n_het and n_hom_var.
    * **n_singleton**: Number of singletons.
    * **n_snp**: Number of SNP alternate alleles.
    * **n_insertion**: Number of insertion alternate alleles.
    * **n_deletion**: Number of deletion alternate alleles.
    * **n_transition**: Number of transition (A-G, C-T) alternate alleles.
    * **n_transversion**: Number of transversion alternate alleles.
    * **n_star**: Number of star alleles.
    * **r_ti_tv**: Transition/Transversion ratio.
    * **r_het_hom_var**: Het/HomVar call ratio.
    * **r_insertion_deletion**: Insertion/Deletion allele ratio.
    * **n_snp_residual**: Residuals after regressing out the first eight ancestry PCs from the number of SNP alternate alleles.
    * **n_singleton_residual**: Residuals after regressing out the first eight ancestry PCs from the number of singletons.
    * **r_ti_tv_residual**: Residuals after regressing out the first eight ancestry PCs from the Transition/Transversion ratio.
    * **r_insertion_deletion_residual**: Residuals after regressing out the first eight ancestry PCs from the Insertion/Deletion allele ratio.
    * **n_insertion_residual**: Residuals after regressing out the first eight ancestry PCs from the number of insertion alternate alleles.
    * **n_deletion_residual**: Residuals after regressing out the first eight ancestry PCs from the number of deletion alternate alleles.
    * **r_het_hom_var_residual**: Residuals after regressing out the first eight ancestry PCs from the Het/HomVar call ratio.
    * **n_transition_residual**: Residuals after regressing out the first eight ancestry PCs from the number of transition (A-G, C-T) alternate alleles.
    * **n_transversion_residual**: Residuals after regressing out the first eight ancestry PCs from the number of transversion alternate alleles.

* **population_inference**: Struct containing ancestry information assigned by applying a principal components analysis (PCA) on samples and using those PCs in a random forest classifier trained on known ancestry labels.
    * **training_pop**: Sample's project_pop or v2_pop if it is not 'oth'.
    * **pca_scores**: Sample's scores for each population PC.
    * **pop**: Sample's inferred population label.
    * **prob_afr**: Random forest probability that the sample is of African/African American ancestry.
    * **prob_ami**: Random forest probability that the sample is of Amish ancestry.
    * **prob_amr**: Random forest probability that the sample is of Latino ancestry.
    * **prob_asj**: Random forest probability that the sample is of Ashkenazi Jewish ancestry.
    * **prob_eas**: Random forest probability that the sample is of East Asian ancestry.
    * **prob_fin**: Random forest probability that the sample is of Finnish ancestry.
    * **prob_mid**: Random forest probability that the sample is of Middle Eastern ancestry.
    * **prob_nfe**: Random forest probability that the sample is of Non-Finnish European ancestry.
    * **prob_oth**: Random forest probability that the sample is of Other ancestry.
    * **prob_sas**: Random forest probability that the sample is of South Asian ancestry.
    * **training_pop_all**: Sample's project_pop or v2_pop if it is not ‘oth'.

* **sample_filters**: Sample QC filter annotations used for the release.
    * **sex_aneuploidy**: Whether the sample was flagged for sex aneuploidy.
    * **insert_size**: Whether the sample was flagged for insert size.
    * **chimera**: Whether the sample was flagged for chimera.
    * contamination: Whether the sample was flagged for contamination.
    * **bad_qc_metrics**: Whether the sample was flagged for hard filters from hard_filter_cutoffs.
    * **low_coverage**: Whether the sample was flagged for low coverage.
    * **ambiguous_sex**: Whether the sample was flagged for ambiguous sex.
    * **failed_fingerprinting**: Whether the sample was flagged for failing fingerprinting.
    * **TCGA_tumor_sample**: Whether the sample was flagged for being a TCGA tumor sample.
    * **hard_filters**: The sample's set of failed hard filters.
    * **hard_filtered**: Whether the sample was hard filtered.
    * **release_related**: Whether the sample had a second-degree or greater relatedness to a sample in the release.
    * **release_duplicate**: Whether the sample was a duplicate of a sample in the release.
    * **release_parent_child**: Whether the sample had a parent child relationship to a sample in the release.
    * **release_sibling**: Whether the sample had a sibling relationship to another sample in the release.
    * **all_samples_related**: Whether the sample had second-degree or greater relatedness to another sample in the raw callset.
    * **all_samples_duplicate**: Whether the sample was a duplicate of another sample in the raw callset.
    * **all_samples_parent_child**: Whether the sample had a parent child relationship to another sample in the raw callset.
    * **all_samples_sibling**: Whether the sample had a sibling relationship to another sample in the raw callset.
    * **fail_n_snp_residual**: Whether the sample failed the number of SNP alternate alleles cutoff based on the metric residuals.
    * **fail_n_singleton_residual**: Whether the sample failed the number of singletons cutoff based on the metric residuals.
    * **fail_r_ti_tv_residual**: Whether the sample failed the Transition/Transversion ratio cutoff based on the metric residuals.
    * **fail_r_insertion_deletion_residual**: Whether the sample failed the Insertion/Deletion allele ratio cutoff based on the metric residuals.
    * **fail_n_insertion_residual**: Whether the sample failed the number of insertion alternate alleles cutoff based on the metric residuals.
    * **fail_n_deletion_residual**: Whether the sample failed the number of deletion alternate alleles cutoff based on the metric residuals.
    * **fail_r_het_hom_var_residual**: Whether the sample failed the Het/HomVar call ratio cutoff based on the metric residuals.
    * **fail_n_transition_residual**: Whether the sample failed the number of transition alternate alleles.
    * **fail_n_transversion_residual**: Whether the sample failed the number of transversion alternate alleles.

* **qc_metrics_filters**: Set of all sample QC metrics for which each sample was found to be an outlier after computing sample QC metrics using Hail’s [sample_qc]([https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.sample_qc)) and regressing out the first 8 ancestry assignment PCs.
* **relatedness_inference**: Information about the sample’s relatedness to other samples within the callset.
    * relationships: Samples that have a kinship estimate (kin) > 0.05 determined using Hail’s [pc_relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) with this sample.

* **high_quality**: Whether the sample is considered high quality (no hard filter or qc_metric filter flags)
* **release**: Whether the sample meets criteria for release (Permission to release, high quality, not to be excluded, and not related to any other sample in the release).
