## Contents and ordering of scripts
The scripts below are run approximately in the order they are listed. We begin quality control with a sparse MT created using Hail's [gVCF combiner](https://hail.is/docs/0.2/experimental/vcf_combiner.html)

### Quality Control
**load_data/**
* `create_last_END_positions.py` — 

**sample_qc/**
* `sample_qc.py`
  * Main arguments: 
    * `--impute_sex` Runs sex imputation. Also runs sex karyotyping annotation
    * `--compute_hard_filters` Computes samples to be hard-filtered
    * `--compute_samples_ranking` Computes global samples ranking based on hard-filters, releasable and coverage.
    * `--compute_qc_mt` Creates the QC MT based on liftover of v2 QC and Purcell 5k sites
    * `--run_pc_relate` Run PC-relate
    * `--compute_related_samples_to_drop` Flags related samples to drop
    * `--run_pca` Compute PCA
    * `--assign_pops` Assigns pops from PCA
    * `--calculate_inbreeding` Calculate sample level inbreeding
    * `--calculate_clinvar` Calculate counts of ClinVar and ClinVar P/LP variants per sample
    * `--apply_stratified_filters` Compute per pop filtering
    * `--apply_regressed_filters` Computes qc_metrics adjusted for pop
    * `--generate_metadata` Generates the metadata HT

* `create_fam.py` -

**annotations/**
* `generate_qc_annotations.py` -

**load_data/**
* `load_vqsr.py` -

**variant_qc/**
* `random_forest.py` -
* `evaluation.py` -
* `final_filter.py` -

### Making the Release
**annotations/**
* `generate_freq_data.py` -
* `generate_qc_annotations.py` -

**create_release/**
* `create_release_sites_ht.py` -
* `create_hgdp_tgp_subset.py` -
* `prepare_vcf_data_release.py` -
* `make_var_annot_hists.py` -

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