# gnomAD v4 production overview:
```mermaid
flowchart TD;
  sample_qc --> annotations;
  click sample_qc "#sampleqc-"
  annotations --> variant_qc;
  click annotations "#annotations-"
  variant_qc --> create_release;
  click variant_qc "#variantqc-"
  click create_release "#createrelease-"
```
## sample_qc:
```mermaid
flowchart TD;
  prepare-cuking-inputs["relatedness.py --prepare-cuking-inputs"] --> print-cuking-command["relatedness.py --print-cuking-command"];
  print-cuking-command["relatedness.py --print-cuking-command"] --> create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"];
  create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"] --> run-ibd-on-cuking-pairs["relatedness.py --run-ibd-on-cuking-pairs"];
  create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"] --> finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"];
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"] --> compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"];
  create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"] --> compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"];
  run-pc-relate-pca["relatedness.py --run-pc-relate-pca"] --> create-pc-relate-relatedness-table["relatedness.py --create-pc-relate-relatedness-table"];
  apply-regressed-filters["outlier_filtering.py --apply-regressed-filters"] --> create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"];
  apply-nearest-neighbor-filters["outlier_filtering.py --apply-nearest-neighbor-filters"] --> create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"];
  determine-nearest-neighbors["outlier_filtering.py --determine-nearest-neighbors"] --> apply-nearest-neighbor-filters["outlier_filtering.py --apply-nearest-neighbor-filters"];
  compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"] --> identify-duplicates["identify_trios.py --identify-duplicates"];
  identify-duplicates["identify_trios.py --identify-duplicates"] --> infer-families["identify_trios.py --infer-families"];
  infer-families["identify_trios.py --infer-families"] --> create-fake-pedigree["identify_trios.py --create-fake-pedigree"];
  infer-families["identify_trios.py --infer-families"] --> run-mendel-errors["identify_trios.py --run-mendel-errors"];
  create-fake-pedigree["identify_trios.py --create-fake-pedigree"] --> run-mendel-errors["identify_trios.py --run-mendel-errors"];
  infer-families["identify_trios.py --infer-families"] --> finalize-ped["identify_trios.py --finalize-ped"];
  run-mendel-errors["identify_trios.py --run-mendel-errors"] --> finalize-ped["identify_trios.py --finalize-ped"];
```
## annotations:
```mermaid
flowchart TD;
  compute-info["generate_variant_qc_annotations.py --compute-info"] --> split-info["generate_variant_qc_annotations.py --split-info"];
  split-info["generate_variant_qc_annotations.py --split-info"] --> create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"];
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"] --> create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"];
  generate-sib-stats["generate_variant_qc_annotations.py --generate-sib-stats"] --> create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"];
  compute-info["generate_variant_qc_annotations.py --compute-info"] --> export-info-vcf["generate_variant_qc_annotations.py --export-info-vcf"];
  compute-info["generate_variant_qc_annotations.py --compute-info"] --> export-split-info-vcf["generate_variant_qc_annotations.py --export-split-info-vcf"];
  run-vep["generate_variant_qc_annotations.py --run-vep"] --> validate-vep["generate_variant_qc_annotations.py --validate-vep"];
  finalize-ped["identify_trios.py --finalize-ped"] --> generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"];
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"] --> generate-sib-stats["generate_variant_qc_annotations.py --generate-sib-stats"];
  write-split-vds-and-downsampling-ht["generate_freq.py --write-split-vds-and-downsampling-ht"] --> run-freq-and-dense-annotations["generate_freq.py --run-freq-and-dense-annotations"];
  combine-freq-hts["generate_freq.py --combine-freq-hts"] --> correct-for-high-ab-hets["generate_freq.py --correct-for-high-ab-hets"];
  correct-for-high-ab-hets["generate_freq.py --correct-for-high-ab-hets"] --> finalize-freq-ht["generate_freq.py --finalize-freq-ht"];
  run-freq-and-dense-annotations["generate_freq.py --run-freq-and-dense-annotations"] --> combine-freq-hts["generate_freq.py --combine-freq-hts"];
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"] --> get-duplicated-to-exomes["generate_freq_genomes.py --get-duplicated-to-exomes"];
  update-annotations["generate_freq_genomes.py --update-annotations"] --> get-callstats-for-updated-samples["generate_freq_genomes.py --get-callstats-for-updated-samples"];
  get-callstats-for-updated-samples["generate_freq_genomes.py --get-callstats-for-updated-samples"] --> join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"];
  update-annotations["generate_freq_genomes.py --update-annotations"] --> compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"];
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"] --> compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"];
  update-annotations["generate_freq_genomes.py --update-annotations"] --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"];
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"] --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"];
  compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"] --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"];
  compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"] --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"];
  update-annotations["generate_freq_genomes.py --update-annotations"] --> apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"];
  update-release-callstats["generate_freq_genomes.py --update-release-callstats"] --> apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"];
  update-annotations["generate_freq_genomes.py --update-annotations"] --> compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"];
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"] --> compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"];
  compute-coverage-ht["compute_coverage.py --compute-coverage-ht"] --> export-release-files["compute_coverage.py --export-release-files"];
```
## variant_qc:
```mermaid
flowchart TD;
  train-rf["random_forest.py --train-rf"] --> apply-rf["random_forest.py --apply-rf"];
  create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"] --> create-bin-ht["evaluation.py --create-bin-ht"];
  create-bin-ht["evaluation.py --create-bin-ht"] --> score-bin-validity-check["evaluation.py --score-bin-validity-check"];
  create-bin-ht["evaluation.py --create-bin-ht"] --> create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"];
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"] --> create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"];
  create-bin-ht["evaluation.py --create-bin-ht"] --> bin-truth-sample-concordance["evaluation.py --bin-truth-sample-concordance"];
  merge-with-truth-data["evaluation.py --merge-with-truth-data"] --> bin-truth-sample-concordance["evaluation.py --bin-truth-sample-concordance"];
  extract-truth-samples["evaluation.py --extract-truth-samples"] --> merge-with-truth-data["evaluation.py --merge-with-truth-data"];
  create-bin-ht["evaluation.py --create-bin-ht"] --> final_filter.py["final_filter.py"];
  create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"] --> final_filter.py["final_filter.py"];
  split-info["generate_variant_qc_annotations.py --split-info"] --> final_filter.py["final_filter.py"];
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"] --> final_filter.py["final_filter.py"];
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"] --> final_filter.py["final_filter.py"];
```
## create_release:
```mermaid
flowchart TD;
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"] --> create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"];
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"] --> create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"];
  final_filter.py["final_filter.py"] --> create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"];
  final_filter.py["final_filter.py"] --> create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"];
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"] --> perform-contingency-table-test["create_combined_faf_release_ht.py --perform-contingency-table-test"];
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"] --> perform-cochran-mantel-haenszel-test["create_combined_faf_release_ht.py --perform-cochran-mantel-haenszel-test"];
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"] --> finalize-combined-faf-release["create_combined_faf_release_ht.py --finalize-combined-faf-release"];
  validate-release-ht["validate_and_export_vcf.py --validate-release-ht"] --> prepare-vcf-header["validate_and_export_vcf.py --prepare-vcf-header"];
  validate-release-ht["validate_and_export_vcf.py --validate-release-ht"] --> export-vcf["validate_and_export_vcf.py --export-vcf"];
  prepare-vcf-header["validate_and_export_vcf.py --prepare-vcf-header"] --> export-vcf["validate_and_export_vcf.py --export-vcf"];
```
## sample_qc:
### [generate_qc_mt.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/generate_qc_mt.py):     Generate combined gnomAD v3 and v4 QC MatrixTable for use in relatedness and ancestry inference.
```mermaid
flowchart TD;

```
### [hard_filters.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/hard_filters.py): Script to determine samples that fail hard filtering thresholds.
```mermaid
flowchart TD;

```
### [platform_inference.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/platform_inference.py): Script to assign platforms based on per interval fraction of bases over DP 0 PCA results using HDBSCAN.
```mermaid
flowchart TD;

```
### [interval_qc.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/interval_qc.py): Script to define high quality intervals based on per interval aggregate statistics over samples.
```mermaid
flowchart TD;

```
### [sex_inference.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/sex_inference.py): Script to impute chromosomal sex karyotype annotation.
```mermaid
flowchart TD;

```
### [hard_filters.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/hard_filters.py): Script to determine samples that fail hard filtering thresholds.
```mermaid
flowchart TD;

```
### [relatedness.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py): Script to compute relatedness estimates among pairs of samples in the callset.
```mermaid
flowchart TD;
  prepare-cuking-inputs["--prepare-cuking-inputs"] --> print-cuking-command["--print-cuking-command"];
  print-cuking-command["--print-cuking-command"] --> create-cuking-relatedness-table["--create-cuking-relatedness-table"];
  create-cuking-relatedness-table["--create-cuking-relatedness-table"] --> run-ibd-on-cuking-pairs["--run-ibd-on-cuking-pairs"];
  create-cuking-relatedness-table["--create-cuking-relatedness-table"] --> finalize-relatedness-ht["--finalize-relatedness-ht"];
  finalize-relatedness-ht["--finalize-relatedness-ht"] --> compute-related-samples-to-drop["--compute-related-samples-to-drop"];
  create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"] --> compute-related-samples-to-drop["--compute-related-samples-to-drop"];
  run-pc-relate-pca["--run-pc-relate-pca"] --> create-pc-relate-relatedness-table["--create-pc-relate-relatedness-table"];
```
### [assign_ancestry.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/assign_ancestry.py): Script to assign global ancestry labels to samples using known v3 population labels or TGP and HGDP labels.
```mermaid
flowchart TD;

```
### [outlier_filtering.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/outlier_filtering.py): Script to determine sample QC metric outliers that should be filtered.
```mermaid
flowchart TD;
  apply-regressed-filters["--apply-regressed-filters"] --> create-finalized-outlier-filter["--create-finalized-outlier-filter"];
  apply-nearest-neighbor-filters["--apply-nearest-neighbor-filters"] --> create-finalized-outlier-filter["--create-finalized-outlier-filter"];
  determine-nearest-neighbors["--determine-nearest-neighbors"] --> apply-nearest-neighbor-filters["--apply-nearest-neighbor-filters"];
```
### [identify_trios.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py): Script to identify trios from relatedness data and filter based on Mendel errors and de novos.
```mermaid
flowchart TD;
  compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"] --> identify-duplicates["--identify-duplicates"];
  identify-duplicates["--identify-duplicates"] --> infer-families["--infer-families"];
  infer-families["--infer-families"] --> create-fake-pedigree["--create-fake-pedigree"];
  infer-families["--infer-families"] --> run-mendel-errors["--run-mendel-errors"];
  create-fake-pedigree["--create-fake-pedigree"] --> run-mendel-errors["--run-mendel-errors"];
  infer-families["--infer-families"] --> finalize-ped["--finalize-ped"];
  run-mendel-errors["--run-mendel-errors"] --> finalize-ped["--finalize-ped"];
```
### [create_sample_qc_metadata_ht.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/create_sample_qc_metadata_ht.py): Script to merge the output of all sample QC modules into a single Table.
```mermaid
flowchart TD;

```
### [create_sample_qc_metadata_ht_genomes.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/create_sample_qc_metadata_ht_genomes.py): Script to create sample QC metadata HT for genomes.
```mermaid
flowchart TD;

```
## annotations:
### [generate_variant_qc_annotations.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py): Script to generate annotations for variant QC on gnomAD v4.
```mermaid
flowchart TD;
  compute-info["--compute-info"] --> split-info["--split-info"];
  split-info["--split-info"] --> create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"];
  generate-trio-stats["--generate-trio-stats"] --> create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"];
  generate-sib-stats["--generate-sib-stats"] --> create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"];
  compute-info["--compute-info"] --> export-info-vcf["--export-info-vcf"];
  compute-info["--compute-info"] --> export-split-info-vcf["--export-split-info-vcf"];
  run-vep["--run-vep"] --> validate-vep["--validate-vep"];
  finalize-ped["identify_trios.py --finalize-ped"] --> generate-trio-stats["--generate-trio-stats"];
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"] --> generate-sib-stats["--generate-sib-stats"];
```
### [generate_freq.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py): Script to generate the frequency data annotations across v4 exomes.
```mermaid
flowchart TD;
  write-split-vds-and-downsampling-ht["--write-split-vds-and-downsampling-ht"] --> run-freq-and-dense-annotations["--run-freq-and-dense-annotations"];
  combine-freq-hts["--combine-freq-hts"] --> correct-for-high-ab-hets["--correct-for-high-ab-hets"];
  correct-for-high-ab-hets["--correct-for-high-ab-hets"] --> finalize-freq-ht["--finalize-freq-ht"];
  run-freq-and-dense-annotations["--run-freq-and-dense-annotations"] --> combine-freq-hts["--combine-freq-hts"];
```
### [generate_freq_genomes.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py): Script to create frequencies HT for v4.0 genomes.
```mermaid
flowchart TD;
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"] --> get-duplicated-to-exomes["--get-duplicated-to-exomes"];
  update-annotations["--update-annotations"] --> get-callstats-for-updated-samples["--get-callstats-for-updated-samples"];
  get-callstats-for-updated-samples["--get-callstats-for-updated-samples"] --> join-callstats-for-update["--join-callstats-for-update"];
  update-annotations["--update-annotations"] --> compute-allele-number-for-new-variants["--compute-allele-number-for-new-variants"];
  join-callstats-for-update["--join-callstats-for-update"] --> compute-allele-number-for-new-variants["--compute-allele-number-for-new-variants"];
  update-annotations["--update-annotations"] --> update-release-callstats["--update-release-callstats"];
  join-callstats-for-update["--join-callstats-for-update"] --> update-release-callstats["--update-release-callstats"];
  compute-allele-number-for-new-variants["--compute-allele-number-for-new-variants"] --> update-release-callstats["--update-release-callstats"];
  compute-allele-number-for-pop-diff["--compute-allele-number-for-pop-diff"] --> update-release-callstats["--update-release-callstats"];
  update-annotations["--update-annotations"] --> apply-patch-to-freq-ht["--apply-patch-to-freq-ht"];
  update-release-callstats["--update-release-callstats"] --> apply-patch-to-freq-ht["--apply-patch-to-freq-ht"];
  update-annotations["--update-annotations"] --> compute-allele-number-for-pop-diff["--compute-allele-number-for-pop-diff"];
  join-callstats-for-update["--join-callstats-for-update"] --> compute-allele-number-for-pop-diff["--compute-allele-number-for-pop-diff"];
```
### [vep_context_ht.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/vep_context_ht.py): Script to add VEP annotations to the GRCh38 context Table.
```mermaid
flowchart TD;

```
### [vrs_annotation_batch.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/vrs_annotation_batch.py): This is a batch script which adds VRS IDs to a Hail Table by creating sharded VCFs, running a vrs-annotation script on each shard, and merging the results into the original Hail Table.
```mermaid
flowchart TD;

```
### [insilico_predictors.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/insilico_predictors.py): Script to generate Hail Tables with in silico predictors.
```mermaid
flowchart TD;

```
### [compute_coverage.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/compute_coverage.py): Script to compute coverage statistics on gnomAD v4 exomes.
```mermaid
flowchart TD;
  compute-coverage-ht["--compute-coverage-ht"] --> export-release-files["--export-release-files"];
```
## variant_qc:
### [random_forest.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py): Script for running random forest model on gnomAD v4 variant QC data.
```mermaid
flowchart TD;
  train-rf["--train-rf"] --> apply-rf["--apply-rf"];
```
### [vqsr.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/vqsr.py): Script to run VQSR on an AS-annotated Sites VCF.
```mermaid
flowchart TD;

```
### [import_variant_qc_vcf.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/import_variant_qc_vcf.py):     Import variant QC result site VCF into a HT.
```mermaid
flowchart TD;

```
### [evaluation.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py): Script to create Tables with aggregate variant statistics by variant QC score bins needed for evaluation plots.
```mermaid
flowchart TD;
  create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"] --> create-bin-ht["--create-bin-ht"];
  create-bin-ht["--create-bin-ht"] --> score-bin-validity-check["--score-bin-validity-check"];
  create-bin-ht["--create-bin-ht"] --> create-aggregated-bin-ht["--create-aggregated-bin-ht"];
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"] --> create-aggregated-bin-ht["--create-aggregated-bin-ht"];
  create-bin-ht["--create-bin-ht"] --> bin-truth-sample-concordance["--bin-truth-sample-concordance"];
  merge-with-truth-data["--merge-with-truth-data"] --> bin-truth-sample-concordance["--bin-truth-sample-concordance"];
  extract-truth-samples["--extract-truth-samples"] --> merge-with-truth-data["--merge-with-truth-data"];
```
### [final_filter.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/final_filter.py): Script to create final filter Table for release.
```mermaid
flowchart TD;
  create-bin-ht["evaluation.py --create-bin-ht"] --> final_filter.py["final_filter.py"];
  create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"] --> final_filter.py["final_filter.py"];
  split-info["generate_variant_qc_annotations.py --split-info"] --> final_filter.py["final_filter.py"];
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"] --> final_filter.py["final_filter.py"];
```
### [final_filter_genomes.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/final_filter_genomes.py): Script to create final filter Table for v4 genomes release.
```mermaid
flowchart TD;
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"] --> final_filter.py["final_filter.py"];
```
## create_release:
### [create_combined_faf_release_ht.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py): Create a joint gnomAD v4 exome and genome frequency and FAF.
```mermaid
flowchart TD;
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"] --> create-combined-frequency-table["--create-combined-frequency-table"];
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"] --> create-combined-frequency-table["--create-combined-frequency-table"];
  final_filter.py["final_filter.py"] --> create-combined-frequency-table["--create-combined-frequency-table"];
  final_filter.py["final_filter.py"] --> create-combined-frequency-table["--create-combined-frequency-table"];
  create-combined-frequency-table["--create-combined-frequency-table"] --> perform-contingency-table-test["--perform-contingency-table-test"];
  create-combined-frequency-table["--create-combined-frequency-table"] --> perform-cochran-mantel-haenszel-test["--perform-cochran-mantel-haenszel-test"];
  create-combined-frequency-table["--create-combined-frequency-table"] --> finalize-combined-faf-release["--finalize-combined-faf-release"];
```
### [create_release_sites_ht.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_release_sites_ht.py): Script to create release sites HT for v4.0 exomes and genomes.
```mermaid
flowchart TD;

```
### [validate_and_export_vcf.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py):
```mermaid
flowchart TD;
  validate-release-ht["--validate-release-ht"] --> prepare-vcf-header["--prepare-vcf-header"];
  validate-release-ht["--validate-release-ht"] --> export-vcf["--export-vcf"];
  prepare-vcf-header["--prepare-vcf-header"] --> export-vcf["--export-vcf"];
```
### [make_var_annot_hists.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/make_var_annot_hists.py):
```mermaid
flowchart TD;

```
