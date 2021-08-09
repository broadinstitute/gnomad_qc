## Contents and ordering of scripts
The scripts below are run approximately in the order they are listed.

### Quality Control
**load_data/**
* `import_vcf.py` — import raw VCF files into Hail, left-align indel representation, and write to MatrixTable format 
* `import_resources.py` — import reference datasets (e.g., Clinvar annotations, methylation and CpG annotations, ExAC site annotations, truth sets) into Hail and write to MatrixTable or Table format
* `load_coverage.py` — import individual-level coverage files, write to MatrixTable and Table formats, compute summary metrics, and export summary data for release
* `import_gnomad_sv.py` - import the gnomAD SV VCF and write it as a MT; merge gnomAD SVs and short variants into a single MT using an allele frequency cutoff for short variants; generate age and GQ histograms

**sample_qc/**
* `apply_hard_filters.py` — filter dataset to bi-allelic, high-call rate, common SNPs; import and annotate relevant metadata; infer chromosomal sex; annotate samples failing hard filtering thresholds; export relevant annotations for ranking samples for removal due to relatedness (see `joint_sample_qc.py`)
* `generate_hardcalls.py` — generate “hardcalls” version of dataset (dropping extraneous GT fields); adjust for sex ploidy; split multi-allelic sites; generate dataset version with non-reference genotypes only
* `exomes_platform_pca.py` — transform exome genotype data into missingness data over exome intervals; run PCA on missingness data; cluster and assign generic platform labels
* `joint_sample_qc.py` — join exome and genome data; filter joint callset to high-call rate variants and remove rare variants; LD-prune joint callset; evaluate relatedness and drop related samples; perform genotype PCA on unrelated samples; project related samples onto PCs; fit random forest model on PCs for samples with known ancestry labels; apply RF model to assign ancestry to remaining samples; evaluate sample quality metrics on a joint platform-and-population basis and flag outlier samples
* `assign_subpops.py` — annotate LD-pruned joint exome and genome data with hard filters and previously inferred platform and population labels; split joint data by relatedness status; filter to samples in continental ancestry group of interest; filter relateds/unrelateds to high-call rate variants; run genotype PCA on unrelateds; project relateds onto subpop PCs; fit random forest model on subpop PCs for samples with known subpopulation labels; apply RF model to assign subpop ancestry to remaining samples
* `finalize_sample_qc.py` — combine hard filter, platform inference, relatedness, ancestry inference, subpopulation inference, and TOPMED annotations into unified metadata tables for both exomes and genomes; define release sample set; collapse small population sizes into “other” category
* `create_fam.py` — depends on relatedness output (`joint_sample_qc.py`) and hard filter flags (`apply_hard_filters.py`); infer all complete trios from kinship coefficients and sex imputation annotations, including duplicates; select representative set of unique trios in call set; evaluate inferred trios by generating random trios and comparing Mendelian violations for random trios against inferred trios; generate final list of inferred trios by filtering out trios with Mendelian violation counts exceeding standard deviation cutoff
* `get_topmed_dups.py` — depends on frequency table annotations (`annotations/generate_frequency_data.py`); create table of high-quality autosomal SNPs shared across gnomAD and TOPMED along with allele count and singleton/doubleton/tripleton annotations (on variants and samples, respectively) for identifying samples present in both gnomAD and TOPMED

**variant_qc/**
* `variantqc.py` — gather variant- and allele-level annotations used as features for training random forests model, impute missing entries; select variants for training examples; train random forests model using specified parameters for depth and number of trees; if specified, test resulting model on pre-selected region; save training data with metadata describing random forest parameters used; apply random forest model to full variant set; annotate results with binned scores (dependent on `create_ranked_scores.py`) and create variant filter status based on SNP and indel quality cutoffs, inbreeding coefficient thresholds
* `create_ranked_scores.py` — create table with rank annotations based on random forest and/or VQSR variant quality scores, with SNPs and indels handled separately, along with optional additional rankings for more specific variant types (e.g., bi-allelic variants, singletons, bi-allelic singletons); bin the variants in the rank tables by all possible rank groupings and add variant annotations from external validation datasets; compute evaluation metrics for each bin for eventual comparison of variant filtering performance across multiple random forest models and VQSR
* `calculate_concordance.py` — compute concordance against truth data (NA12878 and synthetic diploid sample) and/or for duplicate samples across exomes and genomes; create a concordance table binned by rank (both absolute and relative) for a given data type (exomes/genomes), truth sample, and filtering model (e.g., VQSR, random forest) for the purposes of evaluating variant filtering models

### Making the Release
**annotations/**
* `generate_qc_annotations.py` — compute and collate all variant- and allele-level annotations used for variant QC, including VEP; allele type; variant type; call stats (AC, AN, AF, homozygote count) for high-quality, release, and all samples; histograms of GQ, DP, and AB  metric distributions; and trio statistics (Mendelian errors and transmission disequilibrium test results)
* `generate_frequency_data.py` — compute frequencies of variants in gnomAD for various sample groupings (e.g., ancestry, sex, and subpopulation) and for downsampled sets; add filtering allele frequency annotations for ancestry sample groupings
* `generate_ld_data.py` - generate a linkage disequilibrium (LD) pruned set of variants using Hail's [ld_prune](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_prune) method; calculate LD matrices by population using Hail's [ld_matrix](https://hail.is/docs/0.2/methods/genetics.html#hail.methods.ld_matrix) method; calculate LD scores for each population from the LD matrices; calculate cross-pop LD scores from the LD matricies; compute LD between SNVs and SVs

**variant_qc/**
* `correct_fafs.py` — temporary fix to address 2.1 bug; filter allele frequency annotations so that populations in which AC=1 have filtering allele frequencies set to 0
* `make_var_annot_hists.py` — aggregate background distributions of specified variant QC metrics into histograms; first-pass run determines minimum and maximum values per metric; second-pass run aggregates QUAL metrics binned by allele frequency and aggregates other metrics over manually set bins/ranges, writes results out to json files (primarily for browser loading purposes)
* `prepare_data_release.py` — combine frequency, filtering alelle frequency, variant QC, VEP, dbSNP, and variant QC metric histogram annotations into unified table; add frequency annotations for each sample subset (e.g., controls, non-neuro, non-cancer, non-TOPMED); select and reformat release annotations for VCF export; create VCF header text; perform sanity check on release annotations
* `exomes_genomes_coverage.py` - kept only for archival purposes. This script was used to creates coverage files for the [gnomAD v2 flagship LOF paper](https://www.nature.com/articles/s41586-020-2308-7), but is not used as a regular part of gnomAD production


### Others
**load_data/**
* `import_exomes_vcf_on_prem.py` - kept only for archival purposes. It has no practical utility and is outdated

**variant_qc/**
* `select_qc_set.py` - kept only for archival purposes. This script was used for a one-time analysis to assess variant QC, but is not used as a regular part of gnomAD production
* `variant_qc_plots.py` - helper functions to create Bokeh plots for assessment of variant QC results in a Jupyter notebook
