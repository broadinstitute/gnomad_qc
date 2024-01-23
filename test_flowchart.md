## annotations:
### compute_coverage: Script to compute coverage statistics on gnomAD v4 exomes.
```mermaid
flowchart TD;
  --compute-coverage-ht --> --export-release-files
```

### generate_freq: Script to generate the frequency data annotations across v4 exomes.
```mermaid
flowchart TD;
  --write-split-vds-and-downsampling-ht --> --run-freq-and-dense-annotations
  --run-freq-and-dense-annotations --> --combine-freq-hts
  --combine-freq-hts --> --correct-for-high-ab-hets
  --correct-for-high-ab-hets --> --finalize-freq-ht
```

### generate_freq_genomes: Script to create frequencies HT for v4.0 genomes.
```mermaid
flowchart TD;
  --update-annotations --> --get-callstats-for-updated-samples
  --get-callstats-for-updated-samples --> --join-callstats-for-update
  --join-callstats-for-update --> --compute-allele-number-for-new-variants
  --compute-allele-number-for-new-variants --> --update-release-callstats
  --update-release-callstats --> --apply-patch-to-freq-ht
  --join-callstats-for-update --> --compute-allele-number-for-pop-diff
  --compute-allele-number-for-pop-diff --> --update-release-callstats
  --update-release-callstats --> --apply-patch-to-freq-ht
  --join-callstats-for-update --> --update-release-callstats
  --update-release-callstats --> --apply-patch-to-freq-ht
  --update-annotations --> --compute-allele-number-for-new-variants
  --compute-allele-number-for-new-variants --> --update-release-callstats
  --update-release-callstats --> --apply-patch-to-freq-ht
  --update-annotations --> --compute-allele-number-for-pop-diff
  --compute-allele-number-for-pop-diff --> --update-release-callstats
  --update-release-callstats --> --apply-patch-to-freq-ht
  --update-annotations --> --update-release-callstats
  --update-release-callstats --> --apply-patch-to-freq-ht
  --update-annotations --> --apply-patch-to-freq-ht
```

### generate_variant_qc_annotations: Script to generate annotations for variant QC on gnomAD v4.
```mermaid
flowchart TD;
  --compute-info --> --split-info
  --split-info --> --create-variant-qc-annotation-ht
  --compute-info --> --export-info-vcf
  --compute-info --> --export-split-info-vcf
  --run-vep --> --validate-vep
  --generate-trio-stats --> --create-variant-qc-annotation-ht
  --generate-sib-stats --> --create-variant-qc-annotation-ht
```

### insilico_predictors: Script to generate Hail Tables with in silico predictors.
### vep_context_ht: Script to add VEP annotations to the GRCh38 context Table.
### vrs_annotation_batch: This is a batch script which adds VRS IDs to a Hail Table by creating sharded VCFs, running a vrs-annotation script on each shard, and merging the results into the original Hail Table.
## create_release:
### create_combined_faf_release_ht: Create a joint gnomAD v4 exome and genome frequency and FAF.
```mermaid
flowchart TD;
  --create-combined-frequency-table --> --perform-contingency-table-test
  --create-combined-frequency-table --> --perform-cochran-mantel-haenszel-test
  --create-combined-frequency-table --> --finalize-combined-faf-release
```

### create_release_sites_ht: Script to create release sites HT for v4.0 exomes and genomes.
### create_release_utils: Common utils for creating gnomAD v4 exomes and genomes releases.
### make_var_annot_hists:
### validate_and_export_vcf:
```mermaid
flowchart TD;
  --validate-release-ht --> --prepare-vcf-header
  --prepare-vcf-header --> --export-vcf
  --validate-release-ht --> --export-vcf
```

## resources:
### annotations: Script containing annotation related resources.
### basics: Script containing generic resources.
### constants: Script containing version and release constants.
### meta: Script containing metadata related resources.
### release: Script containing release related resources.
### sample_qc: Script containing sample QC related resources.
### variant_qc: Script containing variant QC related resources.
## sample_qc:
### assign_ancestry: Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels.
### create_sample_qc_metadata_ht: Script to merge the output of all sample QC modules into a single Table.
### create_sample_qc_metadata_ht_genomes: Script to create sample QC metadata HT for genomes.
### cuKING:
#### cloud_batch_submit:
#### cuking_outputs_to_ht: Converts cuKING outputs in Parquet format to a Hail table.
#### mt_to_cuking_inputs: Converts a Hail MatrixTable to the Parquet format suitable for cuKING.
### generate_qc_mt: Script to create a dense MatrixTable filtered to a diverse set of variants for relatedness/ancestry PCA using CCDG, gnomAD v3, and UK Biobank.
### hard_filters: Script to determine samples that fail hard filtering thresholds.
### identify_trios: Script to identify trios from relatedness data and filter based on Mendel errors and de novos.
```mermaid
flowchart TD;
  --identify-duplicates --> --infer-families
  --infer-families --> --create-fake-pedigree
  --create-fake-pedigree --> --run-mendel-errors
  --run-mendel-errors --> --finalize-ped
  --infer-families --> --run-mendel-errors
  --run-mendel-errors --> --finalize-ped
  --infer-families --> --finalize-ped
```

### interval_qc: Script to define high quality intervals based on per interval aggregate statistics over samples.
### outlier_filtering: Script to determine sample QC metric outliers that should be filtered.
```mermaid
flowchart TD;
  --apply-regressed-filters --> --create-finalized-outlier-filter
  --determine-nearest-neighbors --> --apply-nearest-neighbor-filters
  --apply-nearest-neighbor-filters --> --create-finalized-outlier-filter
```

### platform_inference: Script to assign platforms based on per interval fraction of bases over DP 0 PCA results using HDBSCAN.
### relatedness: Script to compute relatedness estimates among pairs of samples in the callset.
```mermaid
flowchart TD;
  --prepare-cuking-inputs --> --print-cuking-command
  --print-cuking-command --> --create-cuking-relatedness-table
  --create-cuking-relatedness-table --> --run-ibd-on-cuking-pairs
  --create-cuking-relatedness-table --> --finalize-relatedness-ht
  --finalize-relatedness-ht --> --compute-related-samples-to-drop
  --run-pc-relate-pca --> --create-pc-relate-relatedness-table
```

### sex_inference: Script to impute chromosomal sex karyotype annotation.
## cuKING:
### cloud_batch_submit:
### cuking_outputs_to_ht: Converts cuKING outputs in Parquet format to a Hail table.
### mt_to_cuking_inputs: Converts a Hail MatrixTable to the Parquet format suitable for cuKING.
## subset: Script to filter the gnomAD v4 VariantDataset to a subset of specified samples.
## variant_qc:
### evaluation: Script to create Tables with aggregate variant statistics by variant QC score bins needed for evaluation plots.
```mermaid
flowchart TD;
  --create-bin-ht --> --score-bin-validity-check
  --create-bin-ht --> --create-aggregated-bin-ht
  --create-bin-ht --> --bin-truth-sample-concordance
  --extract-truth-samples --> --merge-with-truth-data
  --merge-with-truth-data --> --bin-truth-sample-concordance
```

### final_filter: Script to create final filter Table for release.
```mermaid
flowchart TD;
```

### final_filter_genomes: Script to create final filter Table for v4 genomes release.
```mermaid
flowchart TD;
```

### import_variant_qc_vcf: Script to load variant QC result VCF into a Hail Table.
### random_forest: Script for running random forest model on gnomAD v4 variant QC data.
```mermaid
flowchart TD;
  --train-rf --> --apply-rf
```

### vqsr: Script to run VQSR on an AS-annotated Sites VCF.
