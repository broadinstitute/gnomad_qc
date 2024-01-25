# gnomAD v4 annotations overview:
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  compute-info["generate_variant_qc_annotations.py --compute-info"]:::step_color --> split_info["split_info"]:::step_color;
  split-info["generate_variant_qc_annotations.py --split-info"]:::step_color --> create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color;
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"]:::step_color --> create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color;
  generate-sib-stats["generate_variant_qc_annotations.py --generate-sib-stats"]:::step_color --> create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color;
  compute-info["generate_variant_qc_annotations.py --compute-info"]:::step_color --> export_info_vcf["export_info_vcf"]:::step_color;
  compute-info["generate_variant_qc_annotations.py --compute-info"]:::step_color --> export_split_info_vcf["export_split_info_vcf"]:::step_color;
  run-vep["generate_variant_qc_annotations.py --run-vep"]:::step_color --> validate_vep["validate_vep"]:::step_color;
  finalize-ped["identify_trios.py --finalize-ped"]:::step_color --> generate_trio_stats["generate_trio_stats"]:::step_color;
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color --> generate_sib_stats["generate_sib_stats"]:::step_color;
  compute-info["generate_variant_qc_annotations.py --compute-info"]:::step_color --> split-info["generate_variant_qc_annotations.py --split-info"]:::step_color;
  split-info["generate_variant_qc_annotations.py --split-info"]:::step_color --> create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"]:::step_color;
  compute-info["generate_variant_qc_annotations.py --compute-info"]:::step_color --> export-info-vcf["generate_variant_qc_annotations.py --export-info-vcf"]:::step_color;
  compute-info["generate_variant_qc_annotations.py --compute-info"]:::step_color --> export-split-info-vcf["generate_variant_qc_annotations.py --export-split-info-vcf"]:::step_color;
  run-vep["generate_variant_qc_annotations.py --run-vep"]:::step_color --> validate-vep["generate_variant_qc_annotations.py --validate-vep"]:::step_color;
  generate-trio-stats["generate_variant_qc_annotations.py --generate-trio-stats"]:::step_color --> create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"]:::step_color;
  generate-sib-stats["generate_variant_qc_annotations.py --generate-sib-stats"]:::step_color --> create-variant-qc-annotation-ht["generate_variant_qc_annotations.py --create-variant-qc-annotation-ht"]:::step_color;
  write-split-vds-and-downsampling-ht["generate_freq.py --write-split-vds-and-downsampling-ht"]:::step_color --> run_freq_and_dense_annotations["run_freq_and_dense_annotations"]:::step_color;
  combine-freq-hts["generate_freq.py --combine-freq-hts"]:::step_color --> correct_for_high_ab_hets["correct_for_high_ab_hets"]:::step_color;
  correct-for-high-ab-hets["generate_freq.py --correct-for-high-ab-hets"]:::step_color --> finalize_freq_ht["finalize_freq_ht"]:::step_color;
  write-split-vds-and-downsampling-ht["generate_freq.py --write-split-vds-and-downsampling-ht"]:::step_color --> run-freq-and-dense-annotations["generate_freq.py --run-freq-and-dense-annotations"]:::step_color;
  run-freq-and-dense-annotations["generate_freq.py --run-freq-and-dense-annotations"]:::step_color --> combine-freq-hts["generate_freq.py --combine-freq-hts"]:::step_color;
  combine-freq-hts["generate_freq.py --combine-freq-hts"]:::step_color --> correct-for-high-ab-hets["generate_freq.py --correct-for-high-ab-hets"]:::step_color;
  correct-for-high-ab-hets["generate_freq.py --correct-for-high-ab-hets"]:::step_color --> finalize-freq-ht["generate_freq.py --finalize-freq-ht"]:::step_color;
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color --> get_duplicated_to_exomes["get_duplicated_to_exomes"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> get_callstats_for_updated_samples["get_callstats_for_updated_samples"]:::step_color;
  get-callstats-for-updated-samples["generate_freq_genomes.py --get-callstats-for-updated-samples"]:::step_color --> join_callstats_for_update["join_callstats_for_update"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> compute_allele_number_for_new_variants["compute_allele_number_for_new_variants"]:::step_color;
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color --> compute_allele_number_for_new_variants["compute_allele_number_for_new_variants"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> update_release_callstats["update_release_callstats"]:::step_color;
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color --> update_release_callstats["update_release_callstats"]:::step_color;
  compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"]:::step_color --> update_release_callstats["update_release_callstats"]:::step_color;
  compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"]:::step_color --> update_release_callstats["update_release_callstats"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> apply_patch_to_freq_ht["apply_patch_to_freq_ht"]:::step_color;
  update-release-callstats["generate_freq_genomes.py --update-release-callstats"]:::step_color --> apply_patch_to_freq_ht["apply_patch_to_freq_ht"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> compute_allele_number_for_pop_diff["compute_allele_number_for_pop_diff"]:::step_color;
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color --> compute_allele_number_for_pop_diff["compute_allele_number_for_pop_diff"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> get-callstats-for-updated-samples["generate_freq_genomes.py --get-callstats-for-updated-samples"]:::step_color;
  get-callstats-for-updated-samples["generate_freq_genomes.py --get-callstats-for-updated-samples"]:::step_color --> join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color;
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color --> compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"]:::step_color;
  compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"]:::step_color --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"]:::step_color;
  update-release-callstats["generate_freq_genomes.py --update-release-callstats"]:::step_color --> apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"]:::step_color;
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color --> compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"]:::step_color;
  compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"]:::step_color --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"]:::step_color;
  join-callstats-for-update["generate_freq_genomes.py --join-callstats-for-update"]:::step_color --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> compute-allele-number-for-new-variants["generate_freq_genomes.py --compute-allele-number-for-new-variants"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> compute-allele-number-for-pop-diff["generate_freq_genomes.py --compute-allele-number-for-pop-diff"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> update-release-callstats["generate_freq_genomes.py --update-release-callstats"]:::step_color;
  update-annotations["generate_freq_genomes.py --update-annotations"]:::step_color --> apply-patch-to-freq-ht["generate_freq_genomes.py --apply-patch-to-freq-ht"]:::step_color;
  compute-coverage-ht["compute_coverage.py --compute-coverage-ht"]:::step_color --> export_release_files["export_release_files"]:::step_color;
  compute-coverage-ht["compute_coverage.py --compute-coverage-ht"]:::step_color --> export-release-files["compute_coverage.py --export-release-files"]:::step_color;
```
### [generate_variant_qc_annotations.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py): Script to generate annotations for variant QC on gnomAD v4.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  compute-info["--compute-info"]:::step_color --> run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color;
  run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color --> split_info["split_info"]:::step_color;
  split_info["split_info"]:::step_color --> split_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L399'>split_info</a>"]]:::func_color;
  split-info["--split-info"]:::step_color --> split_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L399'>split_info</a>"]]:::func_color;
  split_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L399'>split_info</a>"]]:::func_color --> create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color;
  create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color --> create_variant_qc_annotation_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L478'>create_variant_qc_annotation_ht</a>"]]:::func_color;
  generate-trio-stats["--generate-trio-stats"]:::step_color --> run_generate_trio_stats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L427'>run_generate_trio_stats</a>"]]:::func_color;
  run_generate_trio_stats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L427'>run_generate_trio_stats</a>"]]:::func_color --> create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color;
  generate-sib-stats["--generate-sib-stats"]:::step_color --> create_variant_qc_annotation_ht["create_variant_qc_annotation_ht"]:::step_color;
  run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color --> export_info_vcf["export_info_vcf"]:::step_color;
  run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color --> export_split_info_vcf["export_split_info_vcf"]:::step_color;
  run-vep["--run-vep"]:::step_color --> validate_vep["validate_vep"]:::step_color;
  finalize-ped["identify_trios.py --finalize-ped"]:::step_color --> generate_trio_stats["generate_trio_stats"]:::step_color;
  generate_trio_stats["generate_trio_stats"]:::step_color --> run_generate_trio_stats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L427'>run_generate_trio_stats</a>"]]:::func_color;
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color --> generate_sib_stats["generate_sib_stats"]:::step_color;
  run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color --> split-info["--split-info"]:::step_color;
  split_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L399'>split_info</a>"]]:::func_color --> create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"]:::step_color;
  create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"]:::step_color --> create_variant_qc_annotation_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L478'>create_variant_qc_annotation_ht</a>"]]:::func_color;
  run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color --> export-info-vcf["--export-info-vcf"]:::step_color;
  run_compute_info[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L221'>run_compute_info</a>"]]:::func_color --> export-split-info-vcf["--export-split-info-vcf"]:::step_color;
  run-vep["--run-vep"]:::step_color --> validate-vep["--validate-vep"]:::step_color;
  run_generate_trio_stats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_variant_qc_annotations.py#L427'>run_generate_trio_stats</a>"]]:::func_color --> create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"]:::step_color;
  generate-sib-stats["--generate-sib-stats"]:::step_color --> create-variant-qc-annotation-ht["--create-variant-qc-annotation-ht"]:::step_color;
```
### [generate_freq.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py): Script to generate the frequency data annotations across v4 exomes.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  write-split-vds-and-downsampling-ht["--write-split-vds-and-downsampling-ht"]:::step_color --> get_vds_for_freq[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L202'>get_vds_for_freq</a>"]]:::func_color;
  get_vds_for_freq[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L202'>get_vds_for_freq</a>"]]:::func_color --> get_downsampling_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L548'>get_downsampling_ht</a>"]]:::func_color;
  get_downsampling_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L548'>get_downsampling_ht</a>"]]:::func_color --> run_freq_and_dense_annotations["run_freq_and_dense_annotations"]:::step_color;
  run_freq_and_dense_annotations["run_freq_and_dense_annotations"]:::step_color --> densify_and_prep_vds_for_freq[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L505'>densify_and_prep_vds_for_freq</a>"]]:::func_color;
  densify_and_prep_vds_for_freq[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L505'>densify_and_prep_vds_for_freq</a>"]]:::func_color --> generate_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L395'>generate_freq_ht</a>"]]:::func_color;
  generate_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L395'>generate_freq_ht</a>"]]:::func_color --> get_downsampling_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L548'>get_downsampling_ht</a>"]]:::func_color;
  get_downsampling_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L548'>get_downsampling_ht</a>"]]:::func_color --> mt_hists_fields[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L636'>mt_hists_fields</a>"]]:::func_color;
  combine-freq-hts["--combine-freq-hts"]:::step_color --> combine_freq_hts[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L667'>combine_freq_hts</a>"]]:::func_color;
  combine_freq_hts[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L667'>combine_freq_hts</a>"]]:::func_color --> update_non_ukb_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L570'>update_non_ukb_freq_ht</a>"]]:::func_color;
  update_non_ukb_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L570'>update_non_ukb_freq_ht</a>"]]:::func_color --> correct_for_high_ab_hets["correct_for_high_ab_hets"]:::step_color;
  correct_for_high_ab_hets["correct_for_high_ab_hets"]:::step_color --> correct_for_high_ab_hets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L754'>correct_for_high_ab_hets</a>"]]:::func_color;
  correct_for_high_ab_hets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L754'>correct_for_high_ab_hets</a>"]]:::func_color --> generate_faf_grpmax[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L815'>generate_faf_grpmax</a>"]]:::func_color;
  generate_faf_grpmax[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L815'>generate_faf_grpmax</a>"]]:::func_color --> compute_inbreeding_coeff[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L882'>compute_inbreeding_coeff</a>"]]:::func_color;
  correct-for-high-ab-hets["--correct-for-high-ab-hets"]:::step_color --> correct_for_high_ab_hets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L754'>correct_for_high_ab_hets</a>"]]:::func_color;
  compute_inbreeding_coeff[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L882'>compute_inbreeding_coeff</a>"]]:::func_color --> finalize_freq_ht["finalize_freq_ht"]:::step_color;
  finalize_freq_ht["finalize_freq_ht"]:::step_color --> create_final_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L896'>create_final_freq_ht</a>"]]:::func_color;
  get_downsampling_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L548'>get_downsampling_ht</a>"]]:::func_color --> run-freq-and-dense-annotations["--run-freq-and-dense-annotations"]:::step_color;
  run-freq-and-dense-annotations["--run-freq-and-dense-annotations"]:::step_color --> densify_and_prep_vds_for_freq[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L505'>densify_and_prep_vds_for_freq</a>"]]:::func_color;
  mt_hists_fields[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L636'>mt_hists_fields</a>"]]:::func_color --> combine-freq-hts["--combine-freq-hts"]:::step_color;
  update_non_ukb_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L570'>update_non_ukb_freq_ht</a>"]]:::func_color --> correct-for-high-ab-hets["--correct-for-high-ab-hets"]:::step_color;
  compute_inbreeding_coeff[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L882'>compute_inbreeding_coeff</a>"]]:::func_color --> finalize-freq-ht["--finalize-freq-ht"]:::step_color;
  finalize-freq-ht["--finalize-freq-ht"]:::step_color --> create_final_freq_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq.py#L896'>create_final_freq_ht</a>"]]:::func_color;
```
### [generate_freq_genomes.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py): Script to create frequencies HT for v4.0 genomes.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color --> get_duplicated_to_exomes["get_duplicated_to_exomes"]:::step_color;
  update-annotations["--update-annotations"]:::step_color --> add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> get_callstats_for_updated_samples["get_callstats_for_updated_samples"]:::step_color;
  get_callstats_for_updated_samples["get_callstats_for_updated_samples"]:::step_color --> get_updated_release_samples[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L398'>get_updated_release_samples</a>"]]:::func_color;
  get_updated_release_samples[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L398'>get_updated_release_samples</a>"]]:::func_color --> get_hgdp_tgp_callstats_for_selected_samples[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L518'>get_hgdp_tgp_callstats_for_selected_samples</a>"]]:::func_color;
  get_hgdp_tgp_callstats_for_selected_samples[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L518'>get_hgdp_tgp_callstats_for_selected_samples</a>"]]:::func_color --> filter_freq_arrays[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L566'>filter_freq_arrays</a>"]]:::func_color;
  get-callstats-for-updated-samples["--get-callstats-for-updated-samples"]:::step_color --> get_updated_release_samples[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L398'>get_updated_release_samples</a>"]]:::func_color;
  filter_freq_arrays[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L566'>filter_freq_arrays</a>"]]:::func_color --> join_callstats_for_update["join_callstats_for_update"]:::step_color;
  join_callstats_for_update["join_callstats_for_update"]:::step_color --> join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> compute_allele_number_for_new_variants["compute_allele_number_for_new_variants"]:::step_color;
  compute_allele_number_for_new_variants["compute_allele_number_for_new_variants"]:::step_color --> get_group_membership_ht_for_an[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L688'>get_group_membership_ht_for_an</a>"]]:::func_color;
  get_group_membership_ht_for_an[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L688'>get_group_membership_ht_for_an</a>"]]:::func_color --> compute_an_by_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L814'>compute_an_by_group_membership</a>"]]:::func_color;
  join-callstats-for-update["--join-callstats-for-update"]:::step_color --> join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color;
  join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color --> compute_allele_number_for_new_variants["compute_allele_number_for_new_variants"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> update_release_callstats["update_release_callstats"]:::step_color;
  update_release_callstats["update_release_callstats"]:::step_color --> generate_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L948'>generate_v4_genomes_callstats</a>"]]:::func_color;
  generate_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L948'>generate_v4_genomes_callstats</a>"]]:::func_color --> finalize_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1046'>finalize_v4_genomes_callstats</a>"]]:::func_color;
  join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color --> update_release_callstats["update_release_callstats"]:::step_color;
  compute-allele-number-for-new-variants["--compute-allele-number-for-new-variants"]:::step_color --> get_group_membership_ht_for_an[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L688'>get_group_membership_ht_for_an</a>"]]:::func_color;
  compute_an_by_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L814'>compute_an_by_group_membership</a>"]]:::func_color --> update_release_callstats["update_release_callstats"]:::step_color;
  compute-allele-number-for-pop-diff["--compute-allele-number-for-pop-diff"]:::step_color --> compute_an_by_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L814'>compute_an_by_group_membership</a>"]]:::func_color;
  compute_an_by_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L814'>compute_an_by_group_membership</a>"]]:::func_color --> get_pop_diff_v3_vds_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1221'>get_pop_diff_v3_vds_group_membership</a>"]]:::func_color;
  get_pop_diff_v3_vds_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1221'>get_pop_diff_v3_vds_group_membership</a>"]]:::func_color --> update_release_callstats["update_release_callstats"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> apply_patch_to_freq_ht["apply_patch_to_freq_ht"]:::step_color;
  apply_patch_to_freq_ht["apply_patch_to_freq_ht"]:::step_color --> patch_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1259'>patch_v4_genomes_callstats</a>"]]:::func_color;
  patch_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1259'>patch_v4_genomes_callstats</a>"]]:::func_color --> finalize_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1046'>finalize_v4_genomes_callstats</a>"]]:::func_color;
  update-release-callstats["--update-release-callstats"]:::step_color --> generate_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L948'>generate_v4_genomes_callstats</a>"]]:::func_color;
  finalize_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1046'>finalize_v4_genomes_callstats</a>"]]:::func_color --> apply_patch_to_freq_ht["apply_patch_to_freq_ht"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> compute_allele_number_for_pop_diff["compute_allele_number_for_pop_diff"]:::step_color;
  compute_allele_number_for_pop_diff["compute_allele_number_for_pop_diff"]:::step_color --> compute_an_by_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L814'>compute_an_by_group_membership</a>"]]:::func_color;
  join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color --> compute_allele_number_for_pop_diff["compute_allele_number_for_pop_diff"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> get-callstats-for-updated-samples["--get-callstats-for-updated-samples"]:::step_color;
  filter_freq_arrays[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L566'>filter_freq_arrays</a>"]]:::func_color --> join-callstats-for-update["--join-callstats-for-update"]:::step_color;
  join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color --> compute-allele-number-for-new-variants["--compute-allele-number-for-new-variants"]:::step_color;
  compute_an_by_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L814'>compute_an_by_group_membership</a>"]]:::func_color --> update-release-callstats["--update-release-callstats"]:::step_color;
  finalize_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1046'>finalize_v4_genomes_callstats</a>"]]:::func_color --> apply-patch-to-freq-ht["--apply-patch-to-freq-ht"]:::step_color;
  apply-patch-to-freq-ht["--apply-patch-to-freq-ht"]:::step_color --> patch_v4_genomes_callstats[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1259'>patch_v4_genomes_callstats</a>"]]:::func_color;
  join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color --> compute-allele-number-for-pop-diff["--compute-allele-number-for-pop-diff"]:::step_color;
  get_pop_diff_v3_vds_group_membership[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1221'>get_pop_diff_v3_vds_group_membership</a>"]]:::func_color --> update-release-callstats["--update-release-callstats"]:::step_color;
  join_release_ht_with_subsets[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L652'>join_release_ht_with_subsets</a>"]]:::func_color --> update-release-callstats["--update-release-callstats"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> compute-allele-number-for-new-variants["--compute-allele-number-for-new-variants"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> compute-allele-number-for-pop-diff["--compute-allele-number-for-pop-diff"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> update-release-callstats["--update-release-callstats"]:::step_color;
  add_updated_sample_qc_annotations[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/generate_freq_genomes.py#L253'>add_updated_sample_qc_annotations</a>"]]:::func_color --> apply-patch-to-freq-ht["--apply-patch-to-freq-ht"]:::step_color;
```
### [compute_coverage.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/annotations/compute_coverage.py): Script to compute coverage statistics on gnomAD v4 exomes.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  compute-coverage-ht["--compute-coverage-ht"]:::step_color --> export_release_files["export_release_files"]:::step_color;
  compute-coverage-ht["--compute-coverage-ht"]:::step_color --> export-release-files["--export-release-files"]:::step_color;
```
