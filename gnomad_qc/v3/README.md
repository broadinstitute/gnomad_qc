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
