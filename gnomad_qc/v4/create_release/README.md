# gnomAD v4 createrelease overview:
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  final_filter.py["final_filter.py"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  final_filter.py["final_filter.py"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"]:::step_color --> perform_contingency_table_test["perform_contingency_table_test"]:::step_color;
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"]:::step_color --> perform_cochran_mantel_haenszel_test["perform_cochran_mantel_haenszel_test"]:::step_color;
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"]:::step_color --> finalize_combined_faf_release["finalize_combined_faf_release"]:::step_color;
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"]:::step_color --> perform-contingency-table-test["create_combined_faf_release_ht.py --perform-contingency-table-test"]:::step_color;
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"]:::step_color --> perform-cochran-mantel-haenszel-test["create_combined_faf_release_ht.py --perform-cochran-mantel-haenszel-test"]:::step_color;
  create-combined-frequency-table["create_combined_faf_release_ht.py --create-combined-frequency-table"]:::step_color --> finalize-combined-faf-release["create_combined_faf_release_ht.py --finalize-combined-faf-release"]:::step_color;
  validate-release-ht["validate_and_export_vcf.py --validate-release-ht"]:::step_color --> prepare_vcf_header["prepare_vcf_header"]:::step_color;
  validate-release-ht["validate_and_export_vcf.py --validate-release-ht"]:::step_color --> export_vcf["export_vcf"]:::step_color;
  prepare-vcf-header["validate_and_export_vcf.py --prepare-vcf-header"]:::step_color --> export_vcf["export_vcf"]:::step_color;
  validate-release-ht["validate_and_export_vcf.py --validate-release-ht"]:::step_color --> prepare-vcf-header["validate_and_export_vcf.py --prepare-vcf-header"]:::step_color;
  prepare-vcf-header["validate_and_export_vcf.py --prepare-vcf-header"]:::step_color --> export-vcf["validate_and_export_vcf.py --export-vcf"]:::step_color;
  validate-release-ht["validate_and_export_vcf.py --validate-release-ht"]:::step_color --> export-vcf["validate_and_export_vcf.py --export-vcf"]:::step_color;
```
### [create_combined_faf_release_ht.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py): Create a joint gnomAD v4 exome and genome frequency and FAF.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  finalize-freq-ht["generate_freq.py --finalize-freq-ht"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  create_combined_frequency_table["create_combined_frequency_table"]:::step_color --> extract_freq_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L63'>extract_freq_info</a>"]]:::func_color;
  extract_freq_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L63'>extract_freq_info</a>"]]:::func_color --> extract_freq_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L63'>extract_freq_info</a>"]]:::func_color;
  extract_freq_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L63'>extract_freq_info</a>"]]:::func_color --> get_joint_freq_and_faf[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L154'>get_joint_freq_and_faf</a>"]]:::func_color;
  get_joint_freq_and_faf[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L154'>get_joint_freq_and_faf</a>"]]:::func_color --> filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color;
  apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  final_filter.py["final_filter.py"]:::step_color --> create_combined_frequency_table["create_combined_frequency_table"]:::step_color;
  create-combined-frequency-table["--create-combined-frequency-table"]:::step_color --> extract_freq_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L63'>extract_freq_info</a>"]]:::func_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> perform_contingency_table_test["perform_contingency_table_test"]:::step_color;
  perform_contingency_table_test["perform_contingency_table_test"]:::step_color --> perform_contingency_table_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L244'>perform_contingency_table_test</a>"]]:::func_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> perform_cochran_mantel_haenszel_test["perform_cochran_mantel_haenszel_test"]:::step_color;
  perform_cochran_mantel_haenszel_test["perform_cochran_mantel_haenszel_test"]:::step_color --> perform_cmh_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L305'>perform_cmh_test</a>"]]:::func_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> finalize_combined_faf_release["finalize_combined_faf_release"]:::step_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> perform-contingency-table-test["--perform-contingency-table-test"]:::step_color;
  perform-contingency-table-test["--perform-contingency-table-test"]:::step_color --> perform_contingency_table_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L244'>perform_contingency_table_test</a>"]]:::func_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> perform-cochran-mantel-haenszel-test["--perform-cochran-mantel-haenszel-test"]:::step_color;
  perform-cochran-mantel-haenszel-test["--perform-cochran-mantel-haenszel-test"]:::step_color --> perform_cmh_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L305'>perform_cmh_test</a>"]]:::func_color;
  filter_gene_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/create_combined_faf_release_ht.py#L50'>filter_gene_to_test</a>"]]:::func_color --> finalize-combined-faf-release["--finalize-combined-faf-release"]:::step_color;
```
### [validate_and_export_vcf.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py):
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  validate-release-ht["--validate-release-ht"]:::step_color --> check_globals_for_retired_terms[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L892'>check_globals_for_retired_terms</a>"]]:::func_color;
  check_globals_for_retired_terms[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L892'>check_globals_for_retired_terms</a>"]]:::func_color --> prepare_ht_for_validation[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L481'>prepare_ht_for_validation</a>"]]:::func_color;
  prepare_ht_for_validation[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L481'>prepare_ht_for_validation</a>"]]:::func_color --> filter_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/interval_qc.py#L92'>filter_to_test</a>"]]:::func_color;
  filter_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/interval_qc.py#L92'>filter_to_test</a>"]]:::func_color --> prepare_vcf_header["prepare_vcf_header"]:::step_color;
  prepare_vcf_header["prepare_vcf_header"]:::step_color --> prepare_vcf_header_dict[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L730'>prepare_vcf_header_dict</a>"]]:::func_color;
  filter_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/interval_qc.py#L92'>filter_to_test</a>"]]:::func_color --> export_vcf["export_vcf"]:::step_color;
  export_vcf["export_vcf"]:::step_color --> format_validated_ht_for_export[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L808'>format_validated_ht_for_export</a>"]]:::func_color;
  prepare-vcf-header["--prepare-vcf-header"]:::step_color --> prepare_vcf_header_dict[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L730'>prepare_vcf_header_dict</a>"]]:::func_color;
  prepare_vcf_header_dict[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L730'>prepare_vcf_header_dict</a>"]]:::func_color --> export_vcf["export_vcf"]:::step_color;
  filter_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/interval_qc.py#L92'>filter_to_test</a>"]]:::func_color --> prepare-vcf-header["--prepare-vcf-header"]:::step_color;
  prepare_vcf_header_dict[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L730'>prepare_vcf_header_dict</a>"]]:::func_color --> export-vcf["--export-vcf"]:::step_color;
  export-vcf["--export-vcf"]:::step_color --> format_validated_ht_for_export[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/create_release/validate_and_export_vcf.py#L808'>format_validated_ht_for_export</a>"]]:::func_color;
  filter_to_test[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/interval_qc.py#L92'>filter_to_test</a>"]]:::func_color --> export-vcf["--export-vcf"]:::step_color;
```
