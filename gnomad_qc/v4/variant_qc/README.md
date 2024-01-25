# gnomAD v4 variantQC overview:
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  train-rf["random_forest.py --train-rf"]:::step_color --> apply_rf["apply_rf"]:::step_color;
  train-rf["random_forest.py --train-rf"]:::step_color --> apply-rf["random_forest.py --apply-rf"]:::step_color;
  create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"]:::step_color --> create_bin_ht["create_bin_ht"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> score_bin_validity_check["score_bin_validity_check"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> create_aggregated_bin_ht["create_aggregated_bin_ht"]:::step_color;
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"]:::step_color --> create_aggregated_bin_ht["create_aggregated_bin_ht"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> bin_truth_sample_concordance["bin_truth_sample_concordance"]:::step_color;
  merge-with-truth-data["evaluation.py --merge-with-truth-data"]:::step_color --> bin_truth_sample_concordance["bin_truth_sample_concordance"]:::step_color;
  extract-truth-samples["evaluation.py --extract-truth-samples"]:::step_color --> merge_with_truth_data["merge_with_truth_data"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> score-bin-validity-check["evaluation.py --score-bin-validity-check"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> bin-truth-sample-concordance["evaluation.py --bin-truth-sample-concordance"]:::step_color;
  extract-truth-samples["evaluation.py --extract-truth-samples"]:::step_color --> merge-with-truth-data["evaluation.py --merge-with-truth-data"]:::step_color;
  merge-with-truth-data["evaluation.py --merge-with-truth-data"]:::step_color --> bin-truth-sample-concordance["evaluation.py --bin-truth-sample-concordance"]:::step_color;
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  split-info["generate_variant_qc_annotations.py --split-info"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
```
### [random_forest.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py): Script for running random forest model on gnomAD v4 variant QC data.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  train-rf["--train-rf"]:::step_color --> train_rf[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py#L57'>train_rf</a>"]]:::func_color;
  train_rf[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py#L57'>train_rf</a>"]]:::func_color --> add_model_to_run_list[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py#L172'>add_model_to_run_list</a>"]]:::func_color;
  add_model_to_run_list[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py#L172'>add_model_to_run_list</a>"]]:::func_color --> apply_rf["apply_rf"]:::step_color;
  add_model_to_run_list[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/random_forest.py#L172'>add_model_to_run_list</a>"]]:::func_color --> apply-rf["--apply-rf"]:::step_color;
```
### [evaluation.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py): Script to create Tables with aggregate variant statistics by variant QC score bins needed for evaluation plots.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"]:::step_color --> create_bin_ht["create_bin_ht"]:::step_color;
  create_bin_ht["create_bin_ht"]:::step_color --> create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color;
  create-bin-ht["--create-bin-ht"]:::step_color --> create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color;
  create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color --> score_bin_validity_check["score_bin_validity_check"]:::step_color;
  score_bin_validity_check["score_bin_validity_check"]:::step_color --> score_bin_validity_check[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L112'>score_bin_validity_check</a>"]]:::func_color;
  create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color --> create_aggregated_bin_ht["create_aggregated_bin_ht"]:::step_color;
  create_aggregated_bin_ht["create_aggregated_bin_ht"]:::step_color --> create_aggregated_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L135'>create_aggregated_bin_ht</a>"]]:::func_color;
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"]:::step_color --> create_aggregated_bin_ht["create_aggregated_bin_ht"]:::step_color;
  create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color --> bin_truth_sample_concordance["bin_truth_sample_concordance"]:::step_color;
  merge-with-truth-data["--merge-with-truth-data"]:::step_color --> bin_truth_sample_concordance["bin_truth_sample_concordance"]:::step_color;
  extract-truth-samples["--extract-truth-samples"]:::step_color --> merge_with_truth_data["merge_with_truth_data"]:::step_color;
  create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color --> score-bin-validity-check["--score-bin-validity-check"]:::step_color;
  score-bin-validity-check["--score-bin-validity-check"]:::step_color --> score_bin_validity_check[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L112'>score_bin_validity_check</a>"]]:::func_color;
  create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color --> create-aggregated-bin-ht["--create-aggregated-bin-ht"]:::step_color;
  create-aggregated-bin-ht["--create-aggregated-bin-ht"]:::step_color --> create_aggregated_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L135'>create_aggregated_bin_ht</a>"]]:::func_color;
  create_bin_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/evaluation.py#L42'>create_bin_ht</a>"]]:::func_color --> bin-truth-sample-concordance["--bin-truth-sample-concordance"]:::step_color;
  extract-truth-samples["--extract-truth-samples"]:::step_color --> merge-with-truth-data["--merge-with-truth-data"]:::step_color;
  merge-with-truth-data["--merge-with-truth-data"]:::step_color --> bin-truth-sample-concordance["--bin-truth-sample-concordance"]:::step_color;
```
### [final_filter.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/final_filter.py): Script to create final filter Table for release.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  create-bin-ht["evaluation.py --create-bin-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  create-aggregated-bin-ht["evaluation.py --create-aggregated-bin-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  split-info["generate_variant_qc_annotations.py --split-info"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
```
### [final_filter_genomes.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/variant_qc/final_filter_genomes.py): Script to create final filter Table for v4 genomes release.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"]:::step_color --> final_filter.py["final_filter.py"]:::step_color;
```
