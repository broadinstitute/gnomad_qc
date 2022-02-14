## Contents and ordering of scripts
The scripts below are run approximately in the order they are listed. We begin quality control with a sparse MT created using Hail's [gVCF combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)

### Quality Control
**load_data/**
* `create_last_END_positions.py` — Compute the genomic position of the most upstream reference block overlapping each row on the raw sparse MatrixTable. This computation is used in downstream steps to filter to only relevant rows before a densification to only a subset of sites using (densify_sites)[https://github.com/broadinstitute/gnomad_methods/blob/32b3c1d50d3faf300597528272f8efd6973b37e4/gnomad/utils/sparse_mt.py#L84].

**sample_qc/**
* `sample_qc.py`
  * `--sample_qc` - Compute Hail's sample QC metrics on the raw split MatrixTable stratified by bi-allelic and multi-allelic variants.
  * `--impute_sex` - Impute chromosomal sex karyotype annotation.
  * `--compute_hard_filters` - Determine the samples that fail hard filtering thresholds.
  * `--compute_qc_mt` - Filter the full sparse MatrixTable to a smaller dense QC MatrixTable based on liftover of gnomAD v2 QC sites and Purcell 5k sites.
  * `--run_pc_relate` - Run Hail's implementation of [PC-relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) to compute relatedness estimates among pairs of samples in the callset.
  * `--run_pca` - Perform genotype PCA on unrelated samples and project related samples onto PCs.
  * `--assign_pops` - Fit random forest (RF) model on PCs for samples with known ancestry labels and apply that RF model to assign ancestry to remaining samples.
  * `--calculate_inbreeding` Calculate sample level inbreeding coefficient. This is not currently recommended for use because of ancestral differences that will impact this calculation.
  * `--calculate_clinvar` Calculate counts of ClinVar and ClinVar Pathogenic/Likely Pathogenic variants per sample. Used to investigate potential project specific differences.
  * `--apply_stratified_filters` - Evaluate sample quality metrics by population and flag outlier samples.
  * `--apply_regressed_filters` - Regress out the PCs used for the ancestry assignment (`--run_pca`) and flag samples that are outliers based on the residuals for each of the QC metrics.
  * `--compute_related_samples_to_drop` - Determine related samples to drop by computing global sample ranking (based on hard-filters, releasable and coverage) and then using Hail's [maximal_independent_set](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set).
  * `--generate_metadata` - Combine project specific metadata, sample_qc (`--sample_qc`), sex imputation (`--impute_sex`), hard filter (`--compute_hard_filters`), relatedness (`--run_pc_relate` and `--compute_related_samples_to_drop`) and ancestry inference (`--run_pca` and `--assign_pops`) into a unified metadata Table. Define the release sample set.

* `create_fam.py`
  * `--find_dups` - Create a table with duplicate samples indicating which one is the best to use.
  * `--infer_families` - Infer all complete trios from kinship coefficients (`sample_qc.py --run_pc_relate`) and sex imputation annotations (`sample_qc.py --impute_sex`), including duplicates.
  * `--run_mendel_errors` - Calculate Mendelian violations on the inferred complete trios (`--infer_families`) and a random set of trios.
  * `--finalize_ped` - Create a final ped file by excluding families where the number of Mendel errors or de novos are higher than those specified in `--max_dnm` and `--max_mendel`.

**annotations/**
* `generate_qc_annotations.py`
  * `--compute_info` - Compute a Table with the typical GATK allele-specific (AS) and site-level info fields as well as ACs and lowqual fields.
  * `--split_info` - Split the alleles of the info Table (`--compute_info`).
  * `--export_info_vcf` - Export the split info Table (`--split_info`) as a VCF.

**load_data/**
* `load_vqsr.py` - Import VQSR VCF output.

**variant_qc/**
* `random_forest.py`
  * `--annotate_for_rf` - Gather variant- and allele-level annotations used as features for training the variant QC random forests model, impute any missing entries.
  * `--train_rf` - Select variants for training examples and train random forests model using specified parameters for depth and number of trees. If specified, test resulting model on a pre-selected region. Save the training data with metadata describing random forest parameters used.
  * `--apply_rf` - Apply random forest model to full variant set.

* `evaluation.py`
  * `--create_bin_ht` - Create Table with bin annotations based on the ranking of random forest and/or VQSR variant quality scores, with SNPs and indels handled separately. Additional bin annotations are added for the following stratifications: bi-allelic variants, singletons, and bi-allelic singletons.
  * `--create_aggregated_bin_ht` - Compute aggregated metrics for each bin that are useful for comparison of variant filtering performance across multiple random forest models and VQSR.
  * `--extract_truth_samples` - Extract truth samples (NA12878 and synthetic diploid sample) from the full callset MatrixTable for comparison to their truth data.
  * `--merge_with_truth_data` - Compute a table for each truth sample (NA12878 and synthetic diploid sample) comparing the truth sample in the callset to the truth data.
  * `--bin_truth_sample_concordance` - Create a concordance Table of the filtering model (e.g., VQSR, random forest) against truth data (NA12878 and synthetic diploid sample) binned by rank (both absolute and relative. Used for evaluating the variant filtering models.

* `final_filter.py` - Create Table containing the variant filter status based on SNP and indel quality cutoffs and an inbreeding coefficient threshold.

### Making the Release
**annotations/**
* `generate_freq_data.py` - Compute frequencies of variants in gnomAD for various sample groupings (e.g., ancestry, sex, subsets) and for downsampled sets and add filtering allele frequency annotations for ancestry sample groupings.
* `generate_qc_annotations.py`
  * `--generate_allele_data` - Determine the following annotations for each variant in the split info Table: variant type (SNV, indel, multi-SNV, multi-indel, mixed), allele type (SNV, insertion, deletion, complex), and the total number of alleles present at the site.
  * `--generate_ac` - Allele count per variant (raw and adj filtered genotypes) of high quality samples, high quality unrelated samples, and release samples.
  * `--generate_fam_stats` - Calculate transmission and de novo statistics using trios (`sample_qc.py --finalize_ped`).
  * `--vep` - Get the Ensembl Variant Effect Predictor (VEP) annotation for each variant.
  * `--export_transmitted_singletons_vcf` - Export transmitted singletons (`--generate_fam_stats`) to VCF files, used for running VQSR with transmitted singletons in the training.

**create_release/**
* `create_release_sites_ht.py` - Combine frequency, filtering allele frequency, variant QC, VEP, dbSNP, in silico scores, and variant QC metric histogram annotations into unified table. Add frequency annotations for each sample subset (e.g., controls/biobanks, non-neuro, non-cancer, non-TOPMed, non-v2, 1KG/TGP, HGDP).
* `create_hgdp_tgp_subset.py` - Subset raw sparse MatrixTable to only samples in the Human Genome Diversity Project (HGDP) and 1000 Genomes Project (1KG/TGP) to create an unsplit sparse MatrixTable and split dense MatrixTable (with full sample and variant annotations) that can be made publicly available.
* `make_var_annot_hists.py` - Aggregate background distributions of specified variant QC metrics into histograms. The first-pass run determines minimum and maximum values per metric and the second-pass run aggregates QUAL metrics binned by allele frequency and aggregates other metrics over manually set bins/ranges. Writes results out to json files (primarily for browser loading purposes).
* `prepare_vcf_data_release.py` - Select and reformat release annotations for VCF export, create VCF header text, and perform validity checks on release annotations.

### Resources
**load_data/resources/**
* `dbsnp.header.fixed` - Header that fixes the dbSNP version 151 VCF header for import into Hail.
* `vqsr.header.fixed` - Header supplied to `load_vqsr.py --header_path` to fix incorrect Number on AS_SOR description for the classic VQSR method (using site annotations rather than allele-specific annotations).
* `vqsr.alleleSpecific.header.fixed` - Header supplied to `load_vqsr.py --header_path` to fix incorrect Number on AS_SOR description for the allele-specific VQSR method.

**resources/** - Resources for all gnomAD input and output files at each step of the QC and release pipeline.


## Others
**load_data/**
* `compute_coverage.py` — Compute coverage statistics for every base of the GRCh38 reference excluding centromeres and telomeres. The following coverage stats are calculated: mean, median, total DP, and fraction of samples with coverage above x, for each x in [1, 5, 10, 15, 20, 25, 30, 50, 100].
* `compute_ref_block_stats.py` — Compute stats on the length of reference blocks in the raw MatrixTable.
* `split_multi.py` - Script used in v3 to split the raw MatrixTable. v3.1 switched to split the raw MatrixTable on the fly where needed after determining that storage costs were more than rerunning the split when needed.

**sample_qc/**
* `v2_pc_relate.py` - Join the liftover v2 exomes QC MatrixTable (release samples) and the v3 QC MatrixTable then run Hail's [PC-relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) on the v2 exomes/v3 combined QC MatrixTable to obtain relatedness estimates for all pairs of samples in the combined MatrixTable.

**utils.py** - Contains a function that is used in a few scripts within the v3 code. `hom_alt_depletion_fix` adjusts MatrixTable genotypes with a temporary fix for the depletion of homozygous alternate genotypes described in more detail [here](https://gnomad.broadinstitute.org/blog/2020-10-gnomad-v3-1-new-content-methods-annotations-and-data-availability/#tweaks-and-updates).

## gnomAD v3 release HT

The gnomAD release HT annotations are defined below:

**Global fields**:

* **freq_meta**: Allele frequency metadata. An ordered list containing the frequency aggregation group for each element of the ‘freq’ array row annotation.
* **freq_index_dict**: Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD inferred global population, sex: sex karyotype), with values describing the corresponding index of each grouping entry in the ‘freq’ array row annotation.
* **faf_index_dict**: Dictionary keyed by specified label grouping combinations (group: adj/raw, pop: gnomAD inferred global population, sex: sex karyotype), with values describing the corresponding index of each grouping entry in the filtering allele frequency (‘faf’) row annotation.
* **faf_meta**: Filtering allele frequency metadata. An ordered list containing the frequency aggregation group for each element of the ‘faf’ array row annotation.
* **VEP version**: VEP version that was run on the callset. 
* **vep_csq_header**: VEP header for VCF export.
* **dbsnp_version**: dbSNP version used in the callset.
* **filtering_model**: The variant filtering model used and its specific cutoffs.
    * **model_name**: Variant filtering model name used in the 'filters' row annotation, indicating the variant was filtered by this model during variant QC.
    * **score_name**: Annotation name of the score used for variant filtering.
    * **snv_cutoff**: SNV filtering cutoff information.
        * **bin**: Filtering percentile cutoff for SNVs.
        * **min_score**: Minimum score at SNV filtering percentile cutoff.
    * **indel_cutoff**: Indel filtering cutoff information.
        * **bin**: Filtering percentile cutoff for indels.
        * **min_score**: Minimum score at indel filtering percentile cutoff.
    * **model_id**: Variant filtering model ID for score data (used for internal specification of the model).
    * **snv_training_variables**: Variant annotations used as features in the SNV filtering model.
    * **indel_training_variables**: Variant annotations used as features in the indel filtering model.
* **age_distribution**: Callset-wide age histogram calculated on release samples.
    * **bin_edges**: Bin edges for the age histogram.
    * **bin_freq**: Bin frequencies for the age histogram. This is the number of records found in each bin.
    * **n_smaller**: Count of age values falling below lowest histogram bin edge.
    * **n_larger**: Count of age values falling above highest histogram bin edge.
* **freq_sample_count**: A sample count per sample grouping defined in the** ‘**freq_meta’ global annotation.

**Row fields**

* **locus**: Variant locus. Contains contig and position information.
* **alleles**: Variant alleles.
* **freq**: Array of allele frequency information  (AC, AN, AF, homozygote count) for each frequency aggregation group in the gnomAD release. 
    * **AC**: Alternate allele count in release.
    * **AF**: Alternate allele frequency, (AC/AN), in release.
    * **AN**: Total number of alleles in release.
    * **homozygote_count**: Count of homozygous alternate individuals in release.
* **raw_qual_hists**: Genotype quality metric histograms for all genotypes as opposed to high quality genotypes.
    * **gq_hist_all**: Histogram for GQ calculated on all genotypes.
        * **bin_edges**: Bin edges for the GQ histogram calculated on all genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the GQ histogram calculated on all genotypes. The number of records found in each bin.
        * **n_smaller**: Count of GQ values falling below lowest histogram bin edge, for GQ calculated on all genotypes.
        * **n_larger**: Count of GQ values falling above highest histogram bin edge, for GQ calculated on all genotypes.
    * **dp_hist_all**: Histogram for DP calculated on all genotypes.
        * **bin_edges**: Bin edges for the DP histogram calculated on all genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100
        * **bin_freq**: Bin frequencies for the DP histogram calculated on all genotypes. The number of records found in each bin.
        * **n_smaller**: Count of DP values falling below lowest histogram bin edge, for DP calculated on all genotypes.
        * **n_larger**: Count of DP values falling above highest histogram bin edge, for DP calculated on all genotypes.
    * **gq_hist_alt**: Histogram for GQ in heterozygous individuals calculated on all genotypes.
        * **bin_edges**: Bin edges for the histogram of GQ in heterozygous individuals calculated on all genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the histogram of GQ in heterozygous individuals calculated on all genotypes. The number of records found in each bin.
        * **n_smaller**: Count of GQ values in heterozygous individuals falling below lowest histogram bin edge, calculated on all genotypes.
        * **n_larger**: Count of GQ values in heterozygous individuals falling above highest histogram bin edge, calculated on all genotypes.
    * **dp_hist_alt**: Histogram for DP in heterozygous individuals calculated on all genotypes.
        *  **bin_edges**: Bin edges for the histogram of DP in heterozygous individuals calculated on all genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the histogram of DP in heterozygous individuals calculated on all genotypes. The number of records found in each bin.
        * **n_smaller**: Count of DP values in heterozygous individuals falling below lowest histogram bin edge, calculated on all genotypes.
        * **n_larger**: Count of DP values  in heterozygous individuals  falling above highest histogram bin edge, calculated on all genotypes.
    * **ab_hist_alt**: Histogram for AB in heterozygous individuals calculated on all genotypes.
        * **bin_edges**: Bin edges for the histogram of AB in heterozygous individuals calculated on all genotypes are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00.
        * **bin_freq**: Bin frequencies for the histogram of AB in heterozygous individuals calculated on all genotypes. The number of records found in each bin.
        * **n_smaller**: Count of AB values in heterozygous individuals falling below lowest histogram bin edge, calculated on all genotypes.
        * **n_larger**: Count of AB values in heterozygous individuals falling above highest histogram bin edge, calculated on all genotypes.
* **popmax**: Allele frequency information (AC, AN, AF, homozygote count) for the population with maximum allele frequency.
    * **AC**: Alternate allele count in the population with the maximum allele frequency.
    * **AF**: Maximum alternate allele frequency, (AC/AN), across populations in gnomAD.
    * **AN**: Total number of alleles in the population with the maximum allele frequency.
    * **homozygote_count**: Count of homozygous individuals in the population with the maximum allele frequency.
    * **pop**: Population with maximum allele frequency
    * **faf95**: Filtering allele frequency (using Poisson 95% CI) for the population with the maximum allele frequency. 
* **qual_hists**: Genotype quality metric histograms for high quality genotypes.
    * **gq_hist_all**: Histogram for GQ calculated on high quality genotypes.
        * **bin_edges**: Bin edges for the GQ histogram calculated on high quality genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the GQ histogram calculated on high quality genotypes. The number of records found in each bin.
        * **n_smaller**: Count of GQ values falling below the lowest histogram bin edge, calculated on high quality genotypes.
        * **n_larger**: Count of GQ values falling above the highest histogram bin edge, calculated on high quality genotypes.
    * **dp_hist_all**: Histogram for DP calculated on high quality genotypes.
        * **bin_edges**: Bin edges for the DP histogram calculated on high quality genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the DP histogram calculated on high quality genotypes. The number of records found in each bin.
        * **n_smaller**: Count of DP values falling below the lowest histogram bin edge, calculated on high quality genotypes.
        * **n_larger**: Count of DP values falling above the highest histogram bin edge, calculated on high quality genotypes.
    * **gq_hist_alt**: Histogram for GQ in heterozygous individuals calculated on high quality genotypes.
        * **bin_edges**: Bin edges for the histogram of GQ in heterozygous individuals calculated on high quality genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the histogram of GQ in heterozygous individuals calculated on high quality genotypes. The number of records found in each bin.
        * **n_smaller**: Count of GQ values in heterozygous individuals falling below the lowest histogram bin edge, calculated on high quality genotypes.
        * **n_larger**: Count of GQ values in heterozygous individuals falling above the highest histogram bin edge, calculated on high quality genotypes.
    * **dp_hist_alt**: Histogram for DP in heterozygous individuals calculated on high quality genotypes.
        * **bin_edges**: Bin edges for the histogram of DP in heterozygous individuals calculated on high quality genotypes are: 0|5|10|15|20|25|30|35|40|45|50|55|60|65|70|75|80|85|90|95|100.
        * **bin_freq**: Bin frequencies for the histogram of DP in heterozygous individuals calculated on high quality genotypes. The number of records found in each bin.
        * **n_smaller**: Count of DP values in heterozygous individuals falling below the lowest histogram bin edge, calculated on high quality genotypes.
        * **n_larger**: Count of DP values in heterozygous individuals falling above highest histogram bin edge, calculated on high quality genotypes.
    * **ab_hist_alt**: Histogram for AB in heterozygous individuals calculated on high quality genotypes.
        * **bin_edges**: Bin edges for the histogram of AB in heterozygous individuals calculated on high quality genotypes are: 0.00|0.05|0.10|0.15|0.20|0.25|0.30|0.35|0.40|0.45|0.50|0.55|0.60|0.65|0.70|0.75|0.80|0.85|0.90|0.95|1.00.
        * **bin_freq**: Bin frequencies for the histogram of AB in heterozygous individuals calculated on high quality genotypes. The number of records found in each bin.
        * **n_smaller**: Count of AB values in heterozygous individuals falling below the lowest histogram bin edge, calculated on high quality genotypes.
        * **n_larger**: Count of AB values in heterozygous individuals falling above the highest histogram bin edge, calculated on high quality genotypes.
* **faf**: Filtering allele frequency.
    * **faf95**: Filtering allele frequency (using Poisson 95% CI).
    * **faf99**: Filtering allele frequency (using Poisson 99% CI).
* **a_index**: The original index of this alternate allele in the multiallelic representation (1 is the first alternate allele or the only alternate allele in a biallelic variant).
* **was_split**: True if this variant was originally multiallelic, otherwise False.
* **rsid**: dbSNP reference SNP identification (rsID) numbers.
* **filters**: Variant filters; AC0: Allele count is zero after filtering out low-confidence genotypes (GQ &lt; 20; DP &lt; 10; and AB &lt; 0.2 for het calls), AS_VQSR: Failed allele-specific VQSR filtering thresholds of -2.7739 for SNPs and -1.0606 for indels, InbreedingCoeff: GATK InbreedingCoeff &lt; -0.3, PASS: Passed all variant filters.
* **info**: Struct containing typical GATK allele-specific (AS) info fields and additional variant QC fields.
    * **QUALapprox**: Sum of PL[0] values; used to approximate the QUAL score.
    * **SB**: Per-sample component statistics which comprise the Fisher's exact test to detect strand bias. Values are: depth of reference allele on forward strand, depth of reference allele on reverse strand, depth of alternate allele on forward strand, depth of alternate allele on reverse strand.
    * **MQ**: Root mean square of the mapping quality of reads across all samples.
    * **MQRankSum**: Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities.
    * **VarDP**: Depth over variant genotypes (does not include depth of reference samples).
    * **AS_ReadPosRankSum**: Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read position bias.
    * **AS_pab_max**: Maximum p-value over callset for binomial test of observed allele balance for a heterozygous genotype, given expectation of 0.5.
    * **AS_QD**: Allele-specific variant call confidence normalized by depth of sample reads supporting a variant.
    * **AS_MQ**: Allele-specific root mean square of the mapping quality of reads across all samples.
    * **QD**: Variant call confidence normalized by depth of sample reads supporting a variant.
    * **AS_MQRankSum**: Allele-specific z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities.
    * **FS**: Phred-scaled p-value of Fisher's exact test for strand bias.
    * **AS_FS**: Allele-specific phred-scaled p-value of Fisher's exact test for strand bias.
    * **ReadPosRankSum**: Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias.
    * **AS_QUALapprox**: Allele-specific sum of PL[0] values; used to approximate the QUAL score.
    * **AS_SB_TABLE**: Allele-specific forward/reverse read counts for strand bias tests.
    * **AS_VarDP**: Allele-specific depth over variant genotypes (does not include depth of reference samples).
    * **AS_SOR**: Allele-specific strand bias estimated by the symmetric odds ratio test.
    * **SOR**: Strand bias estimated by the symmetric odds ratio test.
    * **singleton**: Variant is seen once in the callset.
    * **transmitted_singleton**: Variant was a callset-wide doubleton that was transmitted within a family from a parent to a child (i.e., a singleton amongst unrelated samples in cohort).
    * **omni**: Variant is present on the Omni 2.5 genotyping array and found in 1000 Genomes data.
    * **mills**: Indel is present in the Mills and Devine data.
    * **monoallelic**: All samples are homozygous alternate for the variant.
    * **AS_VQSLOD**: Allele-specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model.
    * **InbreedingCoeff**: Inbreeding coefficient, the excess heterozygosity at a variant site, computed as 1 - (the number of heterozygous genotypes) / (the number of heterozygous genotypes expected under Hardy-Weinberg equilibrium).
* **vep**: Consequence annotations from Ensembl VEP. More details about VEP output is described [here](https://ensembl.org/info/docs/tools/vep/vep_formats.html#output). VEP was run using the LOFTEE plugin and information about the additional LOFTEE annotations can be found [here](https://github.com/konradjk/loftee).
* **vqsr**: VQSR related variant annotations.
    * **AS_VQSLOD**: Allele-specific log-odds ratio of being a true variant versus being a false positive under the trained VQSR Gaussian mixture model.
    * **AS_culprit**: Allele-specific worst-performing annotation in the VQSR Gaussian mixture model.
    * **NEGATIVE_TRAIN_SITE**: Variant was used to build the negative training set of low-quality variants for VQSR.
    * **POSITIVE_TRAIN_SITE**: Variant was used to build the positive training set of high-quality variants for VQSR.
* **region_flag**: Struct containing flags for problematic regions.
    * **lcr**: Variant falls within a low complexity region.
    * **segdup**: Variant falls within a segmental duplication region.
* **allele_info**: Allele information.
    * **variant_type**: Variant type (snv, indel, multi-snv, multi-indel, or mixed).
    * **allele_type**: Allele type (snv, insertion, deletion, or mixed).
    * **n_alt_alleles**: Total number of alternate alleles observed at variant locus.
    * **was_mixed**: Variant type was mixed.
* **age_hist_het**: Histogram for age in all heterozygous release samples calculated on high quality genotypes.
    * **bin_edges**: Bin edges for the age histogram.
    * **bin_freq**: Bin frequencies for the age histogram. This is the number of records found in each bin.
    * **n_smaller**: Count of age values falling below lowest histogram bin edge.
    * **n_larger**: Count of age values falling above highest histogram bin edge.
* **age_hist_hom**: Histogram for age in all homozygous release samples calculated on high quality genotypes.
    * **bin_edges**: Bin edges for the age histogram.
    * **bin_freq**: Bin frequencies for the age histogram. This is the number of records found in each bin.
    * **n_smaller**: Count of age values falling below lowest histogram bin edge.
    * **n_larger**: Count of age values falling above highest histogram bin edge.
* **cadd**:
    * **phred**: Cadd Phred-like scores ('scaled C-scores') ranging from 1 to 99, based on the rank of each variant relative to all possible 8.6 billion substitutions in the human reference genome. Larger values are more deleterious.
    * **raw_score**: Raw CADD scores are interpretable as the extent to which the annotation profile for a given variant suggests that the variant is likely to be 'observed' (negative values) vs 'simulated' (positive values). Larger values are more deleterious.
    * **has_duplicate**:Whether the variant has more than one CADD score associated with it.
* **revel**:
    * **revel_score**: dbNSFP's Revel score from 0 to 1. Variants with higher scores are predicted to be more likely to be deleterious.
    * **has_duplicate**: Whether the variant has more than one revel_score associated with it.
* **splice_ai**:
    * **splice_ai**: The maximum delta score, interpreted as the probability of the variant being splice-altering.
    * **splice_consequence**: The consequence term associated with the max delta score in 'splice_ai_max_ds'.
    * **has_duplicate**: Whether the variant has more than one splice_ai score associated with it.
* **primate_ai**:
    * **primate_ai_score**: PrimateAI's deleteriousness score from 0 (less deleterious) to 1 (more deleterious).
    * **has_duplicate**: Whether the variant has more than one primate_ai_score associated with it.