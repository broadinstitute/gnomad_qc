## Contents and ordering of scripts
The scripts below are run approximately in the order they are listed. We begin quality control with a sparse MT created using Hail's [gVCF combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)

### Quality Control
**load_data/**
* `create_last_END_positions.py` — Compute the genomic position of the most upstream reference block overlapping each row on the raw sparse MatrixTable.

**sample_qc/**
* `sample_qc.py`
  * `--sample_qc` - Compute Hail's sample QC metrics on the raw split MatrixTable stratified by bi-allelic and multi-allelic variants.
  * `--impute_sex` - Impute chromosomal sex karyotype annotation.
  * `--compute_hard_filters` - Determine the samples that fail hard filtering thresholds.
  * `--compute_samples_ranking` - Compute global sample ranking based on hard-filters, releasable and coverage. This is used by `--compute_related_samples_to_drop` to determine sample removal due to relatedness.
  * `--compute_qc_mt` - Filter the full sparse MatrixTable to a smaller dense QC MatrixTable based on liftover of gnomAD v2 QC sites and Purcell 5k sites.
  * `--run_pc_relate` - Run Hail's implementation of [PC-relate](https://hail.is/docs/0.2/methods/relatedness.html#hail.methods.pc_relate) to compute relatedness estimates among pairs of samples in the callset.
  * `--compute_related_samples_to_drop` - Determines related samples to drop using the ranking computed with `--compute_samples_ranking` and Hail's [maximal_independent_set](https://hail.is/docs/0.2/methods/misc.html#hail.methods.maximal_independent_set).
  * `--run_pca` - Perform genotype PCA on unrelated samples; project related samples onto PCs.
  * `--assign_pops` - Fit random forest (RF) model on PCs for samples with known ancestry labels and apply that RF model to assign ancestry to remaining samples.
  * `--calculate_inbreeding` Calculate sample level inbreeding coefficient. This is not currently recommended for use because of ancestral differences that will impact this calculation.
  * `--calculate_clinvar` Calculate counts of ClinVar and ClinVar Pathogenic/Likely Pathogenic variants per sample. Used to investigate potential project specific differences.
  * `--apply_stratified_filters` - Evaluate sample quality metrics by population and flag outlier samples.
  * `--apply_regressed_filters` - Regress out the PCs used for the ancestry assignment (`--run_pca`) and flag samples that are outliers based on the residuals for each of the QC metrics.
  * `--generate_metadata` - Combine project specific metadata, sample_qc (`--sample_qc`), sex imputation (`--impute_sex`), hard filter (`--compute_hard_filters`), relatedness (`--run_pc_relate` and `--compute_related_samples_to_drop`), ancestry inference (`--run_pca` and `--assign_pops`), into a unified metadata Table. Define the release sample set.

* `create_fam.py`
  * `--infer_families` - Infer all complete trios from kinship coefficients (`sample_qc.py --run_pc_relate`) and sex imputation annotations (`--impute_sex`), including duplicates.
  * `-run_mendel_errors` - Calculate Mendelian violations on the inferred complete trios (`--infer_families`) and a random set of trios.
  * `--finalize_ped` - Create a final ped file by excluding families where the number of Mendel errors or de novos are higher than those specified in `--max_dnm` and `--max_mendel`.

**annotations/**
* `generate_qc_annotations.py`
  * `--compute_info` - Compute a Table with the typical GATK AS and site-level info fields as well as ACs and lowqual fields.
  * `--split_info` - Split the alleles of the info Table (`--compute_info`).
  * `--export_info_vcf` - Export the split info Table (`--split_info`) as a VCF.

**load_data/**
* `load_vqsr.py` - Import the VQSR VCF output.

**variant_qc/**
* `random_forest.py`
  * `--annotate_for_rf` - Gather variant- and allele-level annotations used as features for training the variant QC random forests model, impute any missing entries.
  * `--train_rf` - Select variants for training examples and train random forests model using specified parameters for depth and number of trees. If specified, test resulting model on a pre-selected region. Save the training data with metadata describing random forest parameters used.
  * `--apply_rf` - Apply random forest model to full variant set.

* `evaluation.py` -annotate results with binned scores (dependent on `create_ranked_scores.py`) and 
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
* `create_release_sites_ht.py` - combine frequency, filtering alelle frequency, variant QC, VEP, dbSNP, and variant QC metric histogram annotations into unified table; add frequency annotations for each sample subset (e.g., controls, non-neuro, non-cancer, non-TOPMED); select and reformat release annotations for VCF export; create VCF header text; perform sanity check on release annotations.
* `create_hgdp_tgp_subset.py` -
* `prepare_vcf_data_release.py` -
* `make_var_annot_hists.py` - aggregate background distributions of specified variant QC metrics into histograms; first-pass run determines minimum and maximum values per metric; second-pass run aggregates QUAL metrics binned by allele frequency and aggregates other metrics over manually set bins/ranges, writes results out to json files (primarily for browser loading purposes).

### Resources -- 
**load_data/resources/**
* `dbsnp.header.fixed` -
* `vqsr.alleleSpecific.header.fixed` -
* `vqsr.header.fixed` -

**resources/**
* `annotations.py` -
* `basics.py` -
* `constants.py` -
* `meta.py` -
* `release.py` -
* `sample_qc.py` -
* `variant_qc.py` -

## Others
**load_data/**
* `compute_coverage.py` — 
* `compute_ref_block_stats.py` — 
* `split_multi.py` - 

**sample_qc/**
* `v2_pc_relate.py` - 
* `utils.py` - 