# gnomAD v4 sampleQC overview:
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  prepare-cuking-inputs["relatedness.py --prepare-cuking-inputs"]:::step_color --> print_cuking_command["print_cuking_command"]:::step_color;
  print-cuking-command["relatedness.py --print-cuking-command"]:::step_color --> create_cuking_relatedness_table["create_cuking_relatedness_table"]:::step_color;
  create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"]:::step_color --> run_ibd_on_cuking_pairs["run_ibd_on_cuking_pairs"]:::step_color;
  create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"]:::step_color --> finalize_relatedness_ht["finalize_relatedness_ht"]:::step_color;
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color --> compute_related_samples_to_drop["compute_related_samples_to_drop"]:::step_color;
  create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"]:::step_color --> compute_related_samples_to_drop["compute_related_samples_to_drop"]:::step_color;
  run-pc-relate-pca["relatedness.py --run-pc-relate-pca"]:::step_color --> create_pc_relate_relatedness_table["create_pc_relate_relatedness_table"]:::step_color;
  prepare-cuking-inputs["relatedness.py --prepare-cuking-inputs"]:::step_color --> print-cuking-command["relatedness.py --print-cuking-command"]:::step_color;
  print-cuking-command["relatedness.py --print-cuking-command"]:::step_color --> create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"]:::step_color;
  create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"]:::step_color --> run-ibd-on-cuking-pairs["relatedness.py --run-ibd-on-cuking-pairs"]:::step_color;
  create-cuking-relatedness-table["relatedness.py --create-cuking-relatedness-table"]:::step_color --> finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color;
  finalize-relatedness-ht["relatedness.py --finalize-relatedness-ht"]:::step_color --> compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"]:::step_color;
  run-pc-relate-pca["relatedness.py --run-pc-relate-pca"]:::step_color --> create-pc-relate-relatedness-table["relatedness.py --create-pc-relate-relatedness-table"]:::step_color;
  apply-regressed-filters["outlier_filtering.py --apply-regressed-filters"]:::step_color --> create_finalized_outlier_filter["create_finalized_outlier_filter"]:::step_color;
  apply-nearest-neighbor-filters["outlier_filtering.py --apply-nearest-neighbor-filters"]:::step_color --> create_finalized_outlier_filter["create_finalized_outlier_filter"]:::step_color;
  determine-nearest-neighbors["outlier_filtering.py --determine-nearest-neighbors"]:::step_color --> apply_nearest_neighbor_filters["apply_nearest_neighbor_filters"]:::step_color;
  apply-regressed-filters["outlier_filtering.py --apply-regressed-filters"]:::step_color --> create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"]:::step_color;
  determine-nearest-neighbors["outlier_filtering.py --determine-nearest-neighbors"]:::step_color --> apply-nearest-neighbor-filters["outlier_filtering.py --apply-nearest-neighbor-filters"]:::step_color;
  apply-nearest-neighbor-filters["outlier_filtering.py --apply-nearest-neighbor-filters"]:::step_color --> create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"]:::step_color;
  compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"]:::step_color --> identify_duplicates["identify_duplicates"]:::step_color;
  identify-duplicates["identify_trios.py --identify-duplicates"]:::step_color --> infer_families["infer_families"]:::step_color;
  infer-families["identify_trios.py --infer-families"]:::step_color --> create_fake_pedigree["create_fake_pedigree"]:::step_color;
  infer-families["identify_trios.py --infer-families"]:::step_color --> run_mendel_errors["run_mendel_errors"]:::step_color;
  create-fake-pedigree["identify_trios.py --create-fake-pedigree"]:::step_color --> run_mendel_errors["run_mendel_errors"]:::step_color;
  infer-families["identify_trios.py --infer-families"]:::step_color --> finalize_ped["finalize_ped"]:::step_color;
  run-mendel-errors["identify_trios.py --run-mendel-errors"]:::step_color --> finalize_ped["finalize_ped"]:::step_color;
  identify-duplicates["identify_trios.py --identify-duplicates"]:::step_color --> infer-families["identify_trios.py --infer-families"]:::step_color;
  infer-families["identify_trios.py --infer-families"]:::step_color --> create-fake-pedigree["identify_trios.py --create-fake-pedigree"]:::step_color;
  create-fake-pedigree["identify_trios.py --create-fake-pedigree"]:::step_color --> run-mendel-errors["identify_trios.py --run-mendel-errors"]:::step_color;
  run-mendel-errors["identify_trios.py --run-mendel-errors"]:::step_color --> finalize-ped["identify_trios.py --finalize-ped"]:::step_color;
  infer-families["identify_trios.py --infer-families"]:::step_color --> run-mendel-errors["identify_trios.py --run-mendel-errors"]:::step_color;
  infer-families["identify_trios.py --infer-families"]:::step_color --> finalize-ped["identify_trios.py --finalize-ped"]:::step_color;
```
### [relatedness.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py): Script to compute relatedness estimates among pairs of samples in the callset.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  prepare-cuking-inputs["--prepare-cuking-inputs"]:::step_color --> print_cuking_command["print_cuking_command"]:::step_color;
  print_cuking_command["print_cuking_command"]:::step_color --> print_cuking_command[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L46'>print_cuking_command</a>"]]:::func_color;
  print-cuking-command["--print-cuking-command"]:::step_color --> print_cuking_command[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L46'>print_cuking_command</a>"]]:::func_color;
  print_cuking_command[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L46'>print_cuking_command</a>"]]:::func_color --> create_cuking_relatedness_table["create_cuking_relatedness_table"]:::step_color;
  create-cuking-relatedness-table["--create-cuking-relatedness-table"]:::step_color --> run_ibd_on_cuking_pairs["run_ibd_on_cuking_pairs"]:::step_color;
  run_ibd_on_cuking_pairs["run_ibd_on_cuking_pairs"]:::step_color --> compute_ibd_on_cuking_pair_subset[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L93'>compute_ibd_on_cuking_pair_subset</a>"]]:::func_color;
  create-cuking-relatedness-table["--create-cuking-relatedness-table"]:::step_color --> finalize_relatedness_ht["finalize_relatedness_ht"]:::step_color;
  finalize_relatedness_ht["finalize_relatedness_ht"]:::step_color --> finalize_relatedness_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L175'>finalize_relatedness_ht</a>"]]:::func_color;
  finalize-relatedness-ht["--finalize-relatedness-ht"]:::step_color --> finalize_relatedness_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L175'>finalize_relatedness_ht</a>"]]:::func_color;
  finalize_relatedness_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L175'>finalize_relatedness_ht</a>"]]:::func_color --> compute_related_samples_to_drop["compute_related_samples_to_drop"]:::step_color;
  compute_related_samples_to_drop["compute_related_samples_to_drop"]:::step_color --> run_compute_related_samples_to_drop[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L315'>run_compute_related_samples_to_drop</a>"]]:::func_color;
  create-finalized-outlier-filter["outlier_filtering.py --create-finalized-outlier-filter"]:::step_color --> compute_related_samples_to_drop["compute_related_samples_to_drop"]:::step_color;
  run-pc-relate-pca["--run-pc-relate-pca"]:::step_color --> create_pc_relate_relatedness_table["create_pc_relate_relatedness_table"]:::step_color;
  prepare-cuking-inputs["--prepare-cuking-inputs"]:::step_color --> print-cuking-command["--print-cuking-command"]:::step_color;
  print_cuking_command[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L46'>print_cuking_command</a>"]]:::func_color --> create-cuking-relatedness-table["--create-cuking-relatedness-table"]:::step_color;
  create-cuking-relatedness-table["--create-cuking-relatedness-table"]:::step_color --> run-ibd-on-cuking-pairs["--run-ibd-on-cuking-pairs"]:::step_color;
  run-ibd-on-cuking-pairs["--run-ibd-on-cuking-pairs"]:::step_color --> compute_ibd_on_cuking_pair_subset[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L93'>compute_ibd_on_cuking_pair_subset</a>"]]:::func_color;
  create-cuking-relatedness-table["--create-cuking-relatedness-table"]:::step_color --> finalize-relatedness-ht["--finalize-relatedness-ht"]:::step_color;
  finalize_relatedness_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L175'>finalize_relatedness_ht</a>"]]:::func_color --> compute-related-samples-to-drop["--compute-related-samples-to-drop"]:::step_color;
  compute-related-samples-to-drop["--compute-related-samples-to-drop"]:::step_color --> run_compute_related_samples_to_drop[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/relatedness.py#L315'>run_compute_related_samples_to_drop</a>"]]:::func_color;
  run-pc-relate-pca["--run-pc-relate-pca"]:::step_color --> create-pc-relate-relatedness-table["--create-pc-relate-relatedness-table"]:::step_color;
```

## Relatedness

We use [`cuKING`](https://github.com/populationgenomics/cuKING) for pruning related samples for v4. Make sure that you've cloned the `cuKING` git submodule in this directory (`git submodule update --init`). Then follow these steps:

### Preparing the production environment

Below, `$PROJECT_ID` refers to the GCP project ID for gnomAD production.

```sh
gcloud config set project $PROJECT_ID
```

Enable the necessary services:

```sh
gcloud services enable \
    artifactregistry.googleapis.com \
    batch.googleapis.com \
    cloudbuild.googleapis.com
```

Create an Artifact Registry repository for Docker images:

```sh
gcloud artifacts repositories create images --location=us-central1 --repository-format=docker
```

Build the Docker image:

```sh
cd cuKING
gcloud builds submit --region=us-central1 --config cloudbuild.yaml --substitutions=TAG_NAME=$(git describe --tags) .
```

This should take about half an hour.

Create a service account for running `cuKING`:

```sh
gcloud iam service-accounts create cuking
```

Grant the service account access permissions for:

* the Artifact Registry repository,

  ```sh
  gcloud artifacts repositories add-iam-policy-binding images \
      --location=us-central1 \
      --member=serviceAccount:cuking@$PROJECT_ID.iam.gserviceaccount.com \
      --role=roles/artifactregistry.reader
  ```

* Cloud Batch related roles on the project level, and

  ```sh
  for role in role=roles/batch.agentReporter roles/serviceusage.serviceUsageConsumer roles/logging.logWriter roles/monitoring.metricWriter
  do
      gcloud projects add-iam-policy-binding $PROJECT_ID --member=serviceAccount:cuking@$PROJECT_ID.iam.gserviceaccount.com --role=$role
  done
  ```

* the `gnomad` and `gnomad-tmp` buckets.

  ```sh
  for bucket in gnomad gnomad-tmp
  do
      gsutil iam ch serviceAccount:cuking@$PROJECT_ID.iam.gserviceaccount.com:objectAdmin gs://$bucket
  done
  ```

Create a VM instance template:

```sh
gcloud compute instance-templates create cuking-instance-template \
    --project=$PROJECT_ID \
    --machine-type=a2-highgpu-1g \
    --network-interface=network=default,network-tier=PREMIUM,address="" \
    --metadata-from-file=startup-script=instance_startup_script.sh \
    --maintenance-policy=TERMINATE \
    --provisioning-model=STANDARD \
    --service-account=cuking@$PROJECT_ID.iam.gserviceaccount.com \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --accelerator=count=1,type=nvidia-tesla-a100 \
    --create-disk=auto-delete=yes,boot=yes,device-name=cuking-instance-template,image=projects/ubuntu-os-cloud/global/images/ubuntu-1804-bionic-v20220712,mode=rw,size=10,type=pd-balanced \
    --no-shielded-secure-boot \
    --shielded-vtpm \
    --shielded-integrity-monitoring \
    --reservation-affinity=any
```

The above instance template uses A2 machines with NVIDIA A100 GPUs. Make sure that your [quotas](https://console.cloud.google.com/iam-admin/quotas) are set sufficiently high for `NVIDIA_A100_GPUS` and `A2_CPUS` (e.g. 16 and 192, respectively).

Make sure that the [firewall configuration](https://console.cloud.google.com/networking/firewalls/list) contains a rule to allow internal traffic (like [`default-allow-internal`](https://cloud.google.com/vpc/docs/firewalls#default_firewall_rules) for the default network).

### Prepare inputs

Start an autoscaling Hail Dataproc cluster with 20 secondary workers and submit `relatedness.py --prepare-inputs`. This should take about an hour.

### Run `cuKING`

This step doesn't need a Dataproc cluster. Run the command printed by `relatedness.py --print-cuking-command` to submit a Cloud Batch job. This should take about two hours.

### Convert outputs

`cuKING` writes Parquet files. To convert these to a standard Hail Table, run `relatedness.py --create-relatedness-table` on a small Dataproc cluster. This should only take a minute.

### Determine samples to drop

To prune samples up to second degree relatedness, run `relatedness.py --compute-related-samples-to-drop` on a small Dataproc cluster. This should only take a few minutes.

### Inferring relationships for closely related samples

`get_relatedness_annotated_ht()` annotates the relatedness Table for closely related samples based on the kinship and IBS values computed by cuKING.

### [outlier_filtering.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/outlier_filtering.py): Script to determine sample QC metric outliers that should be filtered.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  apply-regressed-filters["--apply-regressed-filters"]:::step_color --> create_finalized_outlier_filter["create_finalized_outlier_filter"]:::step_color;
  create_finalized_outlier_filter["create_finalized_outlier_filter"]:::step_color --> create_finalized_outlier_filter_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/outlier_filtering.py#L587'>create_finalized_outlier_filter_ht</a>"]]:::func_color;
  apply-nearest-neighbor-filters["--apply-nearest-neighbor-filters"]:::step_color --> create_finalized_outlier_filter["create_finalized_outlier_filter"]:::step_color;
  determine-nearest-neighbors["--determine-nearest-neighbors"]:::step_color --> apply_nearest_neighbor_filters["apply_nearest_neighbor_filters"]:::step_color;
  apply-regressed-filters["--apply-regressed-filters"]:::step_color --> create-finalized-outlier-filter["--create-finalized-outlier-filter"]:::step_color;
  create-finalized-outlier-filter["--create-finalized-outlier-filter"]:::step_color --> create_finalized_outlier_filter_ht[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/outlier_filtering.py#L587'>create_finalized_outlier_filter_ht</a>"]]:::func_color;
  determine-nearest-neighbors["--determine-nearest-neighbors"]:::step_color --> apply-nearest-neighbor-filters["--apply-nearest-neighbor-filters"]:::step_color;
  apply-nearest-neighbor-filters["--apply-nearest-neighbor-filters"]:::step_color --> create-finalized-outlier-filter["--create-finalized-outlier-filter"]:::step_color;
```
### [identify_trios.py](https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py): Script to identify trios from relatedness data and filter based on Mendel errors and de novos.
```mermaid
flowchart TB;
  classDef script_color fill:#2C5D4A
  classDef step_color fill:#2C5D4A
  classDef func_color fill:#2C5D4A
  classDef gnomad_methods_color fill:#2C5D4A
  classDef hail_color fill:#2C5D4A
  classDef resource_color fill:#2C5D4A
  classDef validity_check_color fill:#2C5D4A
  compute-related-samples-to-drop["relatedness.py --compute-related-samples-to-drop"]:::step_color --> identify_duplicates["identify_duplicates"]:::step_color;
  identify-duplicates["--identify-duplicates"]:::step_color --> infer_families["infer_families"]:::step_color;
  infer_families["infer_families"]:::step_color --> filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color;
  infer-families["--infer-families"]:::step_color --> filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color;
  filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color --> create_fake_pedigree["create_fake_pedigree"]:::step_color;
  create_fake_pedigree["create_fake_pedigree"]:::step_color --> run_create_fake_pedigree[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L105'>run_create_fake_pedigree</a>"]]:::func_color;
  filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color --> run_mendel_errors["run_mendel_errors"]:::step_color;
  run_mendel_errors["run_mendel_errors"]:::step_color --> run_mendel_errors[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L140'>run_mendel_errors</a>"]]:::func_color;
  create-fake-pedigree["--create-fake-pedigree"]:::step_color --> run_create_fake_pedigree[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L105'>run_create_fake_pedigree</a>"]]:::func_color;
  run_create_fake_pedigree[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L105'>run_create_fake_pedigree</a>"]]:::func_color --> run_mendel_errors["run_mendel_errors"]:::step_color;
  filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color --> finalize_ped["finalize_ped"]:::step_color;
  finalize_ped["finalize_ped"]:::step_color --> filter_ped[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L218'>filter_ped</a>"]]:::func_color;
  filter_ped[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L218'>filter_ped</a>"]]:::func_color --> families_to_trios[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L46'>families_to_trios</a>"]]:::func_color;
  run-mendel-errors["--run-mendel-errors"]:::step_color --> run_mendel_errors[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L140'>run_mendel_errors</a>"]]:::func_color;
  run_mendel_errors[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L140'>run_mendel_errors</a>"]]:::func_color --> finalize_ped["finalize_ped"]:::step_color;
  identify-duplicates["--identify-duplicates"]:::step_color --> infer-families["--infer-families"]:::step_color;
  filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color --> create-fake-pedigree["--create-fake-pedigree"]:::step_color;
  run_create_fake_pedigree[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L105'>run_create_fake_pedigree</a>"]]:::func_color --> run-mendel-errors["--run-mendel-errors"]:::step_color;
  run_mendel_errors[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L140'>run_mendel_errors</a>"]]:::func_color --> finalize-ped["--finalize-ped"]:::step_color;
  finalize-ped["--finalize-ped"]:::step_color --> filter_ped[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L218'>filter_ped</a>"]]:::func_color;
  filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color --> run-mendel-errors["--run-mendel-errors"]:::step_color;
  filter_ped_to_same_platform[["<a href='https://github.com/broadinstitute/gnomad_qc/tree/main/gnomad_qc/v4/sample_qc/identify_trios.py#L188'>filter_ped_to_same_platform</a>"]]:::func_color --> finalize-ped["--finalize-ped"]:::step_color;
```
