# Sample QC

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


## Sample QC metadata Table annotation definitions

### Global annotations:

* **sample_qc_parameters**: Parameters passed to Hail's `sample_qc` module.
	* **gq_bins**: Tuple containing cutoffs for genotype quality (GQ) scores.
	* **dp_bins**: Tuple containing cutoffs for depth (DP) scores.

* **platform_inference_parameters**: Parameters used during platform assignment.
	* **hdbscan_min_cluster_size**: HDBSCAN `min_cluster_size` parameter.
	* **hdbscan_min_samples**: HDBSCAN `min_samples` parameter.
	* **n_pcs**: Number of platform PCs used for platform assignment.
	* **calling_interval_name**: Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'.
	* **calling_interval_padding**: Number of base pair padding to use on the calling intervals. One of 0 or 50 bp.

* **sex_imputation_parameters**: Parameters used during chromosomal sex imputation.
	* **normalization_contig**: Which autosomal chromosome was used for normalizing the coverage of chromosomes X and Y.
	* **variants_only_x_ploidy**: Whether depth of variant data within calling intervals was used instead of reference data for chrX ploidy estimation.
	* **variants_only_y_ploidy**: Whether depth of variant data within calling intervals was used instead of reference data for chrY ploidy estimation.
	* **variants_filter_lcr**: Whether variants in LCR regions were filtered out of the variants only ploidy estimation and fraction of homozygous alternate variants on chromosome X.
	* **variants_segdup**: Whether variants in segdup regions were filtered out of the variants only ploidy estimation and fraction of homozygous alternate variants on chromosome X.
	* **variants_filter_decoy**: Whether variants in decoy regions were filtered out for variants only ploidy estimation and fraction of homozygous alternate variants on chromosome X.
	* **variants_snv_only**: Whether variants were filtered to only single nucleotide variants for variants only ploidy estimation and fraction of homozygous alternate variants on chromosome X.
	* **f_stat_cutoff**: F-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
	* **f_stat_min_af**: Minimum alternate allele frequency cutoff used to filter sites.
	* **f_stat_ukb_var**: Whether UK Biobank high callrate (0.99) and common variants (UKB allele frequency > value specified by '--min-af') were used for f-stat computation instead of the sites determined by '--determine-fstat-sites'.
	* **x_ploidy_cutoffs**: X ploidy cutoffs for karyotype assignment per platform.
		* **platform_-1**: Platform -1 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_13**: Platform 13 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_1**: Platform 1 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_12**: Platform 12 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_4**: Platform 4 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_19**: Platform 19 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_18**: Platform 18 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_8**: Platform 8 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_10**: Platform 10 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_17**: Platform 17 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_0**: Platform 0 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_7**: Platform 7 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_2**: Platform 2 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_3**: Platform 3 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_9**: Platform 9 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_15**: Platform 15 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_14**: Platform 14 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_6**: Platform 6 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_5**: Platform 5 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_16**: Platform 16 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
		* **platform_11**: Platform 11 X ploidies cutoffs.
			* **upper_cutoff_X**: Upper cutoff for X ploidy.
			* **lower_cutoff_XX**: Lower cutoff for XX ploidy.
			* **upper_cutoff_XX**: Upper cutoff for XX ploidy.
			* **lower_cutoff_XXX**: Lower cutoff for XXX ploidy.
	* **y_ploidy_cutoffs**: Y ploidy cutoffs for karyotype assignment per platform.
		* **platform_-1**: Platform -1 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_13**: Platform 13 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_1**: Platform 1 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_12**: Platform 12 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_4**: Platform 4 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_19**: Platform 19 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_18**: Platform 18 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_8**: Platform 8 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_10**: Platform 10 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_17**: Platform 17 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_0**: Platform 0 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_7**: Platform 7 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_2**: Platform 2 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_3**: Platform 3 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_9**: Platform 9 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_15**: Platform 15 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_14**: Platform 14 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_6**: Platform 6 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_5**: Platform 5 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_16**: Platform 16 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
		* **platform_11**: Platform 11 Y ploidies cutoffs.
			* **lower_cutoff_Y**: Lower cutoff for Y ploidy.
			* **upper_cutoff_Y**: Upper cutoff for Y ploidy.
			* **lower_cutoff_YY**: Lower cutoff for YY ploidy.
	* **x_frac_hom_alt_cutoffs**: Cutoffs for the fraction homozygous alternate genotypes (hom-alt/(hom-alt + het)) on chromosome X per platform.
		* **platform_-1**: Platform -1 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_13**: Platform 13 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_1**: Platform 1 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_12**: Platform 12 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_4**: Platform 4 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_19**: Platform 19 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_18**: Platform 18 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_8**: Platform 8 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_10**: Platform 10 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_17**: Platform 17 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_0**: Platform 0 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_7**: Platform 7 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_2**: Platform 2 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_3**: Platform 3 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_9**: Platform 9 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_15**: Platform 15 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_14**: Platform 14 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_6**: Platform 6 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_5**: Platform 5 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_16**: Platform 16 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.
		* **platform_11**: Platform 11 fraction chromosome X homozygous alternate genotypes (hom-alt/(hom-alt + het)) cutoffs.
			* **lower_cutoff_more_than_one_X**: Lower cutoff for more than one X.
			* **upper_cutoff_more_than_one_X**: Upper cutoff for more than one X.
			* **lower_cutoff_single_X**: Lower cutoff for a single X.

* **hard_filter_metrics_parameters**: Parameters used during hard filtering.
	* **contamination_approximation_parameters**: Parameters used to compute contamination estimate as the mean of reference allele balances of high-quality (DP >= contam-dp-cutoff), autosomal, bi-allelic homozygous SNVs per sample.
		* **dp_cutoff**: Minimum genotype depth to be included in contamination estimate calculation.
	* **chr20_sample_mean_dp_parameters**: Parameters used to compute per sample mean DP on chromosome 20 using interval coverage results.
		* **calling_interval_name**: Name of calling intervals used for interval coverage. One of: 'ukb', 'broad', or 'intersection'.
		* **calling_interval_padding**: Number of base pair padding used on the calling intervals. One of 0 or 50 bp.
	* **sample_qc_mt_callrate_parameters**: Parameters used to compute sample QC MatrixTable callrate.
		* **min_af**: Minimum AF for inclusion in sample callrate computation.
		* **min_site_callrate**: Filtering threshold used for sample callrate computed on only predetermined QC variants (predetermined using CCDG genomes/exomes, gnomAD v3.1 genomes, and UKB exomes) after ADJ filtering.
	* **bi_allelic_sample_qc_parameters**: Parameters used to compute sample QC metrics on bi-allelic SNVs.
		* **gq_bins**: Tuple containing cutoffs for genotype quality (GQ) scores.
		* **dp_bins**: Tuple containing cutoffs for depth (DP) scores.

* **population_inference_parameters**: Parameters used to assign ancestry groups.
	* **assign_pops_from_pc_params**: Parameters used in a random forest model to assign ancestry group labels based on the results of PCA.
		* **min_assignment_prob**: Minimum probability of belonging to a given ancestry group for the ancestry group to be set (otherwise set to 'remaining').
		* **error_rate**: Random forests error rate on the training data for acestry group assignment.
	* **min_prob**: Minimum probability of belonging to a given ancestry group for the ancestry group to be set (otherwise set to 'remaining').
	* **include_unreleasable_samples**: Whether unreleasable samples were included in PCA.
	* **pcs**: List of PCs used in the RF.
	* **include_v2_known_in_training**: Whether v2 known populations were included in the RF training data.
	* **v3_population_spike**: List of v3 populations spiked into the RF training data.
	* **v4_population_spike**: List of v4 populations spiked into the RF training data.
	* **min_prob_cutoffs**: The minimum probability of belonging to a given ancestry group for the ancestry group to be set (otherwise set to 'remaining').

* **relatedness_inference_parameters**: Parameters used for relatedness inference.
	* **relationship_inference_method**: Which relatedness method was used for finalized relatedness Table. Options are 'cuking' and 'pc_relate'.
	* **relationship_cutoffs**: Cutoffs for relatedness inference.
		* **second_degree_sibling_lower_cutoff_slope**: Slope of the line used as a lower cutoff for second degree relatives and siblings from parent-child pairs.
		* **second_degree_sibling_lower_cutoff_intercept**: Intercept of the line used as a lower cutoff for second degree relatives and siblings from parent-child pairs.
		* **second_degree_upper_sibling_lower_cutoff_slope**: Slope of the line used as an upper cutoff for second degree relatives and a lower cutoff for siblings.
		* **second_degree_upper_sibling_lower_cutoff_intercept**: Intercept of the line used as an upper cutoff for second degree relatives and a lower cutoff for siblings.
		* **duplicate_twin_min_kin**: Minimum kinship for duplicate or twin pairs.
		* **second_degree_min_kin**: Minimum kinship threshold used to filter a pair of samples with a second degree relationship when filtering related individuals. This cutoff was determined by evaluation of the kinship distribution.
		* **duplicate_twin_ibd1_min**: Minimum IBD1 cutoff for duplicate or twin pairs if `relatedness_method` is 'pc_relate'. Note: the min is used because pc_relate can output large negative values in some corner cases.
		* **duplicate_twin_ibd1_max**: Maximum IBD1 cutoff for duplicate or twin pairs if `relatedness_method` is 'pc_relate'.
		* **parent_child_max_ibs0_over_ibs2**: Maximum value of IBD0 (if `relatedness_method` is 'pc_relate') or IBS0/IBS2 (if `relatedness_method` is 'cuking') for a parent-child pair.

* **sample_filters_parameters**: Parameters used for sample filtering.
	* **hard_filter_cutoffs**: Cutoffs used for hard-filtering samples prior to sample QC.
		* **max_n_singleton**: Maximum number of singletons.
		* **max_r_het_hom_var**: Maximum ratio of heterozygous variants to homozygous variants.
		* **min_bases_dp_over_1**: Minimum number of bases with a depth over 1.
		* **min_bases_dp_over_20**: Minimum number of bases with a depth over 20.
		* **max_chimera**: Maximum chimera (this is a proportion not a percent, e.g. 5% == 0.05, %5 != 5).
		* **max_contamination_estimate**: Maximum contamination estimate.
		* **chr_20_dp_threshold**: Minimum mean depth over chromosome 20.
	* **outlier_detection_parameters**: The parameters used for filtering outlier samples by QC metrics.
		* **regressed_filters_globals**: The sample filter cutoffs determined by the residuals for each of the sample QC metrics after PCs computed during ancestry assignment were regressed out.
			* **lms**: Linear regression statistics for QC metrics.
				* **r_ti_tv**: Transition : transversion (TiTv) ratio regression statistics.
					* **beta**: Estimated regression coefficient for each PC.
					* **standard_error**: Estimated standard error for each PC.
					* **t_stat**: The t-statistic for each PC.
					* **p_value**: The p-value for each PC.
					* **multiple_standard_error**: Estimated standard deviation of the random error.
					* **multiple_r_squared**: Coefficient of determination for nested models.
					* **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
					* **f_stat**: F-statistic for nested models.
					* **multiple_p_value**: The p-value for the F-test of nested models.
					* **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
				* **r_ti_tv_singleton**: Singleton transition : transversion (TiTv) ratio regression statistics.
					* **beta**: Estimated regression coefficient for each PC.
					* **standard_error**: Estimated standard error for each PC.
					* **t_stat**: The t-statistic for each PC.
					* **p_value**: The p-value for each PC.
					* **multiple_standard_error**: Estimated standard deviation of the random error.
					* **multiple_r_squared**: Coefficient of determination for nested models.
					* **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
					* **f_stat**: F-statistic for nested models.
					* **multiple_p_value**: The p-value for the F-test of nested models.
					* **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
				* **r_het_hom_var**: Heterozygous : homozygous ratio regression statistics.
					* **beta**: Estimated regression coefficient for each PC.
					* **standard_error**: Estimated standard error for each PC.
					* **t_stat**: The t-statistic for each PC.
					* **p_value**: The p-value for each PC.
					* **multiple_standard_error**: Estimated standard deviation of the random error.
					* **multiple_r_squared**: Coefficient of determination for nested models.
					* **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
					* **f_stat**: F-statistic for nested models.
					* **multiple_p_value**: The p-value for the F-test of nested models.
					* **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
				* **r_insertion_deletion**: Insertion : deletion ratio regression statistics.
					* **beta**: Estimated regression coefficient for each PC.
					* **standard_error**: Estimated standard error for each PC.
					* **t_stat**: The t-statistic for each PC.
					* **p_value**: The p-value for each PC.
					* **multiple_standard_error**: Estimated standard deviation of the random error.
					* **multiple_r_squared**: Coefficient of determination for nested models.
					* **adjusted_r_squared**: Adjusted multiple_r_squared taking into account degrees of freedom.
					* **f_stat**: F-statistic for nested models.
					* **multiple_p_value**: The p-value for the F-test of nested models.
					* **n**: Number of samples included in the regression. A sample is included if and only if y, all elements of x, and weight (if set) are non-missing.
			* **qc_metrics_stats**: The statistics generated for outlier filter cutoffs based on the residuals for each of the Hail sample QC metrics.
				* **r_ti_tv_residual**: Transition : transversion (TiTv) ratio residuals statistics.
					* **median**: The transition : transversion (TiTv) ratio median.
					* **mad**: The transition : transversion (TiTv) ratio median aboslute deviation(MAD).
					* **lower**: The transition : transversion (TiTv) ratio lower median absolute deviation threshold.
					* **upper**: The transition : transversion (TiTv) ratio upper median absolute deviation threshold.
				* **r_insertion_deletion_residual**: Insertion : deletion ratio residuals statistics.
					* **median**: The insertion : deletion ratio median.
					* **mad**: The insertion : deletion ratio median aboslute deviation(MAD).
					* **lower**: The insertion : deletion ratio lower median absolute deviation threshold.
					* **upper**: The insertion : deletion ratio upper median absolute deviation threshold.
				* **r_het_hom_var_residual**: Heterozygous : homozygous ratio residuals statistics.
					* **median**: The heterozygous : homozygous ratio median.
					* **mad**: The heterozygous : homozygous ratio median aboslute deviation(MAD).
					* **lower**: The heterozygous : homozygous ratio lower median absolute deviation threshold.
					* **upper**: The heterozygous : homozygous ratio upper median absolute deviation threshold.
				* **r_ti_tv_singleton_residual**: Singleton transition : transversion (TiTv) ratio residuals statistics.
					* **median**: The singleton transition : transversion (TiTv) ratio median.
					* **mad**: The singleton transition : transversion (TiTv) ratio median aboslute deviation(MAD).
					* **lower**: The singleton transition : transversion (TiTv) ratio lower median absolute deviation threshold.
					* **upper**: The singleton transition : transversion (TiTv) ratio upper median absolute deviation threshold.
			* **strata**: Strata used in outlier filtering.
			* **qc_metrics**: QC Metrics used in outlier filtering.
			* **regress_pop_n_pcs**: Number of population PCA scores used in regression.
			* **regress_per_platform**: Whether the regression per platform ran.
			* **apply_r_ti_tv_singleton_filter**: Whether the TiTv singleton filter was applied.
			* **include_unreleasable_samples**: Whether unreleasable samples were included in the regression and filtering.
		* **nearest_neighbor_filters_globals**: Filters used in nearest neighbor filterings.
			* **include_unreleasable_samples**: Whether unreleasable samples were included in the regression and filtering.
			* **apply_r_ti_tv_singleton_filter**: Whether the TiTv singleton filter was applied.
			* **nearest_neighbors_platform_stratified**: Whether the nearest neighbor filtering was stratified by platform.
			* **nearest_neighbors_approximation**: Whether the nearest neighbors approximation used the Annoy nearest neighbors algorithm.
		* **filter_qc_metrics**: List of sample QC metrics used for the finalized outlier filtering.
		* **ensemble_combination_operator**: Logical operator used for combining filtering methods for the ensemble method.
		* **include_unreleasable_samples**: Whether unreleasable samples were included in filtering.

* **date**: Date sample QC metadata HT was generated.

### Row annotations:

* **s**: Sample ID.

* **bam_metrics**: Sample metrics generated from the BAM.
	* **contam_rate**: Freemix value (verifyBamID defintion: Sequence-only estimate of contamination (0-1 scale)). This value can be from either version 1 or version 2 of verifyBamID (version used has not been recorded).
	* **target_bases_20x_rate**: Picard definition: The fraction of all target bases achieving 20X or greater coverage.
	* **median_target_coverage**: Picard definition: The median coverage of a target region.
	* **mean_target_coverage**: Picard definition: The mean coverage of a target region.
	* **target_bases_10x_rate**: Picard definition: The fraction of all target bases achieving 10X or greater coverage.
	* **mean_insert_size**: Picard definition: The mean insert size of the core of the distribution.
	* **chimeras_rate**: Picard definition: The fraction of reads that map outside of a maximum insert size (usually 100kb) or that have the two ends mapping to different chromosomes.
	* **median_insert_size**: Picard definition: The median insert size of all paired end reads where both ends mapped to the same chromosome. Note: the median insert size is not reliably reported in internal metadata - it is sometimes obtained from a coordinator-supplied file, from the Broad database via API request, or from a terra workspace.

* **project_meta**: Sample project metadata.
	* **terra_workspace**: Name of workspace to which sample belongs on Terra.
	* **investigator**: Principal investigator who provided a sample to the Broad.
	* **cohort**: Abbreviation/description of the cohort to which a sample belongs.
	* **orsp_cg**: Consent group code (as assigned by Office of Sponsored Research).
	* **consent**: Data use category (description of abbreviations can be found at https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4721915/).
	* **seq_id**: Combination of project and subject_id connected by '_'.
	* **project**: Research project ID. For internal samples this will be the RP project ID, C project ID (antiquated), or G project ID (antiquated). C and G project IDs exist for older projects that began on-prem. For external samples, 'project' will typically be a short description of the cohort.
	* **subject_id**: Collaborator-assigned sample ID.
	* **sm_id**: Sample ID (assigned by Broad Genomics Platform).
	* **pdo**: Genomics platform's product order in JIRA. Used to track a sample through the lab (e.g. PDO-13180).
	* **capture**: Technology used for sequencing a sample (e.g. Nextera).
	* **sex**: Project provided chromosomal sex.
	* **cram**: Path to cram at the time of metadata collection, if available.
	* **gvcf**: Path to the original gVCF.
	* **disease_description**: Abbreviated description of the phenotype associated with the cohort to which a sample belongs (eg 'IBD').
	* **primary_disease**: Short description of the phenotype associated with the cohort to which a sample belongs, typically more descriptive than 'disease_description'.
	* **race_ethnicity**: Collaborator provided ethnicity of a sample.
	* **delivered_gvcf**: Path to which the reblocked gVCF was delivered. Note: this file may no longer exist.
	* **controls**: Collaborator provided case or control status. Sometimes also determined by presence of 'case' or 'control' in the terra workspace name.
	* **age**: Collaborator provided age in years (e.g. 25).
	* **v2_meta**: Metadata from gnomAD v2.
		* **v2_age**: Last known sample age in gnomAD v2 metadata.
		* **v2_sex**: Sample's chromosomal sex in gnomAD v2 metadata.
		* **v2_hard_filters**: Sample's hard filters in gnomAD v2 ('contamination', 'callrate', 'chimera', 'insert_size', 'coverage', 'ambiguous_sex', and 'sex_aneuploidy').
		* **v2_perm_filters**: Permission filters in gnomad v2 ('tcga_tumor', 'tcga_barcode', 'tcga_below_30', 'specific_exclusion', 'esp', 'not_releasable', and 'syndip').
		* **v2_pop_platform_filters**: Population- and platform-specific outlier filters from gnomAD v2.
		* **v2_related**: Whether a sample was found to be related to another sample within gnomAD v2.
		* **v2_data_type**: Sample datatype, exomes or genomes. Metadata from only gnomAD v2 exomes was merged into this v4 exome metadata Table.
		* **v2_product**: Technology used within the Genomics Platform (e.g. PCR-Free Human WGS - 20x v1.1) from gnomAD v2 metadata.
		* **v2_product_simplified**: Reduced list of sample products ('v2_product'), where similar products are combined.
		* **v2_qc_platform**: GnomAD v2 PCA imputed platform based on sample call rate.
		* **v2_project_id**: GnomAD v2 project ID: RP project, C project ID (antiquated), or G project ID (antiquated). Will match 'seq_project' and 'research_project_key' if the project is in gnomAD v3 and began after cloud processing. C and G project IDs exist for older projects that began on-prem.
		* **v2_project_description**: Text description of 'v2_project_id' (e.g. PROMIS, GTEx).
		* **v2_internal**: Whether or not a sample was processed at the Broad from gnomAD v2 metadata. If False, a sample was processed externally.
		* **v2_investigator**: GnomAD v2 Collaborator, similar to 'contact_pi' which was used in gnomAD v3.
		* **v2_known_pop**: Collaborator provided global population from gnomAD v2.
		* **v2_known_subpop**: Collaborator provided subcontinental population from gnomAD v2.
		* **v2_pop**: Sample's inferred global population label in gnomAD v2.
		* **v2_subpop**: Sample's inferred subcontinental population label in gnomAD v2.
		* **v2_neuro**: Whether a sample is in a neurology cohort and thus excluded from the v2 non-neuro subset. See gnomAD's dataset selection for more information.
		* **v2_control**: Whether a sample is a control in gnomAD v2. See gnomAD's dataset selection for more information.
		* **v2_topmed**: Whether a sample is a TOPMED sample in gnomAD v2. See gnomAD's dataset selection for more information.
		* **v2_high_quality**: Whether a sample is considered high-quality in gnomAD v2, i.e. sample does not have any hard filters or population and platform outlier filters.
		* **v2_release**: Whether a sample is in the gnomAD v2 release (sample was high-quality and not related to any other sample in the v2 release dataset). The field equates to True if a sample is found in the v2/v3 relatedness HT with kinship coefficient > 0.4 (considered a duplicate), the v3 sample ID ('s' in this v3 dataset)) matches the previous v2 genome sample ID ('s' in the in the v2 dataset), or the newly added 1000 genome sample ID ('s' in the v3 dataset) matches a v2 exome sample ID ('s' in the v2 dataset) and the matching v2 exome was released.
		* **v2_release_2_0_2**: Sample was included in the gnomAD v2.0.2 release.
		* **v2_f_stat**: Inbreeding coefficient (excess heterozygosity) on chromosome X from the gnomAD v2 callset.
		* **v2_ambiguous_sex**: Whether a sample was flagged for ambiguous sex in the gnomAD v2 callset.
		* **v2_sex_aneuploidy**: Whether a sample was flagged for sex aneuploidy in the gnomAD v2 callset.
		* **v2_normalized_y_coverage**: Normalized coverage on chromosome Y (normalized using chromosome 20) from the gnomAD v2 metadata.
		* **v2_freemix**: Freemix value (verifyBamID defintion: Sequence-only estimate of contamination (0-1 scale)) from the gnomAD v2 callset.
	* **remapped_primary_disease**: Remapping of 'primary_disease' into more general categories.
	* **ukb_meta**: UKB metadata.
		* **ukb_projected_pop**: Inferred population after projecting onto gnomAD population PCs for UKB samples.
		* **ukb_array_sex**: Sample sex as inferred from genotype data for UKB samples.
		* **ukb_ambiguous_sex**: Whether a sample was flagged for ambiguous sex in the UKB callset.
		* **ukb_sex_aneuploidy**: Whether a sample was flagged for sex aneuploidy in the UKB callset.
		* **ukb_sex_imputation**: UKB sex imputation results.
			* **is_female**: Whether a sample was imputed female in the UKB callset (based on default values of Hail's impute_sex module).
			* **f_stat**: Inbreeding coefficient (excess heterozygosity) on chromosome X from the UKB callset.
			* **n_called**: Number of variants with a genotype call from the UKB callset.
			* **expected_homs**: Expected number of homozygotes from the UKB callset.
			* **observed_homs**: Observed number of homozygotes from the UKB callset.
			* **chr20_mean_dp**: Sample's mean depth across chromosome 20 from the UKB callset.
			* **chrX_mean_dp**: Sample's mean depth across chromosome X from the UKB callset.
			* **chrY_mean_dp**: Sample's mean depth across chromosome Y from the UKB callset.
			* **chrX_ploidy**: Sample's chromosome X ploidy (chrX_mean_dp normalized using chr20_mean_dp) from the UKB callset.
			* **chrY_ploidy**: Sample's chromosome Y ploidy (chrY_mean_dp normalized using chr20_mean_dp) from the UKB callset.
			* **X_karyotype**: Sample's chromosome X karyotype from the UKB callset.
			* **Y_karyotype**: Sample's chromosome Y karyotype from the UKB callset.
			* **sex_karyotype**: Sample's sex karyotype (combined X and Y karyotype) from the UKB callset.
		* **eid_26041**: UK Biobank (UKB) sample ID from application 26041.
		* **eid_31063**: UK Biobank (UKB) sample ID from application 31063.
		* **eid_48511**: UK Biobank (UKB) sample ID from application 48511.
		* **ukb_withdraw**: True if the UKB sample should be withdrawn from the VDS and any analyses based on withdraw list for application 31063 on 02/22/2022.
		* **ukb_batch**: Sample's data tranche (150K, 100K, 200K, 300K).
		* **ukb_batch_num**: Numeric representation of batch. 0: 150K, 1: 100K, 2: 200K, 3: 300K.
	* **releasable**: Whether a sample has permissions for release. Releasable is not the same as 'release' ('release' means a sample was actually included in an official gnomAD release, the 'release' annotation is not included in this Table - it is added after sample QC).
	* **broad_external**: Whether a sample was processed at the Broad or externally.
	* **ukb_sample**: Whether a sample is from the UK Biobank.
	* **gatk_version**: GATK version used to call the sample.
	* **fixed_homalt_model**: Whether a sample is depleted for homozygous alternate calls.

* **sample_qc**: Sample QC metrics generated by Hail's `sample_qc` module on variants with fewer than 3 alternate alleles.
	* **n_het**: Number of heterozygous calls.
	* **n_hom_var**: Number of homozygous alternate calls.
	* **n_non_ref**: Number of homozygous reference calls.
	* **n_singleton**: Number of private alleles. Reference alleles are never counted as singletons, even if every other allele at a site is non-reference.
	* **n_singleton_ti**: Number of singletons that are transitions.
	* **n_singleton_tv**: Number of singletons that are transversions.
	* **n_snp**: Number of SNP alternate alleles.
	* **n_insertion**: Number of insertion alternate alleles.
	* **n_deletion**: Number of deletion alternate alleles.
	* **n_transition**: Number of transition (A-G, C-T) alternate alleles.
	* **n_transversion**: Number of transversion alternate alleles.
	* **n_star**: Number of star (upstream deletion) alleles.
	* **r_ti_tv**: Transition/Transversion ratio.
	* **r_ti_tv_singleton**: Singleton Transition/Transversion ratio.
	* **r_het_hom_var**: Het/HomVar call ratio.
	* **r_insertion_deletion**: Insertion/Deletion allele ratio.
	* **bases_over_gq_threshold**: Number of bases in the interval over passed GQ thresholds.
	* **bases_over_dp_threshold**: Number of bases in the interval over passed DP thresholds.

* **platform_inference**: Struct containing platform inference information. Platforms were assigned by applying a principal components analysis (PCA) on the per-sample per-interval fraction of bases over DP 0 and using HBDSCAN on those PCs to assign a platform.
	* **scores**: Array of principal components analysis scores used when assigning platforms with HBDSCAN.
	* **qc_platform**: Inferred sequencing platform label.

* **sex_imputation**: Annotations generated or used in chromosomal sex inference.
	* **var_data_chr20_mean_dp**: Sample's mean depth of non-reference genotypes in calling intervals on chromosome 20.
	* **chrX_mean_dp**: Sample's mean depth across calling intervals on chromosome X. If 'sex_imputation_parameters.variants_only_x_ploidy' is True, only the depth of non-reference genotypes is used. Otherwise depth of reference blocks is also included in the mean.
	* **chrX_ploidy**: Sample's chromosome X ploidy ('chrX_mean_dp' normalized using 'autosomal_mean_dp' or 'var_data_chr20_mean_dp' if 'sex_imputation_parameters.variants_only_x_ploidy').
	* **chrY_mean_dp**: Sample's mean depth across calling intervals on chromosome Y. If 'sex_imputation_parameters.variants_only_y_ploidy' is True, only the depth of non-reference genotypes is used. Otherwise depth of reference blocks is also included in the mean.
	* **chrY_ploidy**: Sample's chromosome Y ploidy ('chrY_mean_dp' normalized using 'autosomal_mean_dp' or 'var_data_chr20_mean_dp' if 'sex_imputation_parameters.variants_only_Y_ploidy').
	* **autosomal_mean_dp**: Sample's mean depth across calling intervals on the normalization contig (defined by 'sex_imputation_parameters.normalization_contig'). If 'sex_imputation_parameters.variants_only_y_ploidy' is True, only the depth of non-reference genotypes is used. Otherwise depth of reference blocks is also included in the mean.
	* **chrx_frac_hom_alt**: The fraction of homozygous alternate variants on chromosome X.
	* **chrx_frac_hom_alt_adj**: Sample's fraction of high-quality (GQ >= 20; DP >= 10; and AB >= 0.2 for het calls) homozygous alternate genotypes on chromosome X.
	* **platform**: Platform assigned during platform inference.
	* **X_karyotype**: Sample's inferred chromosome X karyotype.
	* **Y_karyotype**: Sample's inferred chromosome Y karyotype.
	* **sex_karyotype**: Sample's inferred sex karyotype (combined X and Y karyotype).
	* **impute_sex_stats**: Statistics calculated by Hail's `impute_sex` module.
		* **f_stat**: Inbreeding coefficient (excess heterozygosity) on chromosome X.
		* **n_called**: Number of variants with a genotype call.
		* **expected_homs**: Expected number of homozygotes.
		* **observed_homs**: Observed number of homozygotes.

* **hard_filter_metrics**: Annotations generated or used during hard-filtering.
	* **mean_AB_snp_biallelic**: Mean bi-allelic SNP allele balance used as a contamination estimate.
	* **chr20_mean_dp**: Chromosome 20 mean DP.
	* **sample_qc_mt_callrate**: Sample callrate in the v4 precomputed QC MatrixTable.
	* **sample_qc_mt_callrate_adj**: Sample callrate for high-quality (GQ >= 20; DP >= 10; and AB >= 0.2 for het calls) genotypes in the v4 precomputed QC MatrixTable.
	* **bi_allelic_sample_qc**: Sample QC metrics calculated from bi-allelic sites only.
		* **n_het**: Number of heterozygous calls.
		* **n_hom_var**: Number of homozygous alternate calls.
		* **n_non_ref**: Number of homozygous reference calls.
		* **n_singleton**: Number of private alleles. Reference alleles are never counted as singletons, even if every other allele at a site is non-reference.
		* **n_singleton_ti**: Number of singletons that are transitions.
		* **n_singleton_tv**: Number of singletons that are transversions.
		* **n_snp**: Number of SNP alternate alleles.
		* **n_insertion**: Number of insertion alternate alleles.
		* **n_deletion**: Number of deletion alternate alleles.
		* **n_transition**: Number of transition (A-G, C-T) alternate alleles.
		* **n_transversion**: Number of transversion alternate alleles.
		* **n_star**: Number of star (upstream deletion) alleles.
		* **r_ti_tv**: Transition/Transversion ratio.
		* **r_ti_tv_singleton**: Singleton Transition/Transversion ratio.
		* **r_het_hom_var**: Het/HomVar call ratio.
		* **r_insertion_deletion**: Insertion/Deletion allele ratio.
		* **bases_over_gq_threshold**: Number of bases in the interval over passed GQ thresholds.
		* **bases_over_dp_threshold**: Number of bases in the interval over passed DP thresholds.

* **population_inference**: Struct containing genetic ancestry information assigned by applying a principal components analysis (PCA) on samples and using those PCs in a random forest classifier trained on known ancestry labels.
	* **training_pop**: Sample's genetic ancestry group if used in training the Random Forest.
	* **pca_scores**: Sample's scores for each genetic ancestry PC.
	* **pop**: Sample's inferred genetic ancestry group label.
	* **prob_afr**: Random forest probability that the sample clusters with samples in the African/African-American genetic ancestry group.
	* **prob_ami**: Random forest probability that the sample clusters with samples in the Amish genetic ancestry group.
	* **prob_amr**: Random forest probability that the sample clusters with samples in the Admixed American genetic ancestry group.
	* **prob_asj**: Random forest probability that the sample clusters with samples in the Ashkenazi Jewish genetic ancestry group.
	* **prob_eas**: Random forest probability that the sample clusters with samples in the East Asian genetic ancestry group.
	* **prob_fin**: Random forest probability that the sample clusters with samples in the Finnish genetic ancestry group.
	* **prob_mid**: Random forest probability that the sample clusters with samples in the Middle Eastern genetic ancestry group.
	* **prob_nfe**: Random forest probability that the sample clusters with samples in the Non-Finnish European genetic ancestry group.
	* **prob_sas**: Random forest probability that the sample clusters with samples in the South Asian genetic ancestry group.
	* **evaluation_sample**: Whether a sample was used in the Random Forest's evaluation.
	* **training_sample**: Whether a sample was used in the Random Forest's training.

* **relatedness_inference**: Struct containing information about a sampleâs relatedness to other samples within the callset.
	* **relationships**: Dictionary of all relationships (second-degree or closer) a sample has with other samples in the dataset. The key is the relationship and the value is a set of all samples with that key relationship to the given sample.
	* **relationships_high_quality**: Dictionary of all relationships (second-degree or closer) a high-quality sample has with other high-quality samples in the dataset. The key is the relationship and the value is a set of all high-quality samples with that key relationship to the given sample.
	* **gnomad_v3_duplicate**: Whether sample is a duplicate of one of all v3.1 samples passing hard-filtering.
	* **gnomad_v3_release_duplicate**: Whether sample is a duplicate of one of the v3.1 release samples.

* **sample_filters**: Sample QC filter annotations used for the release.
	* **chimera**: Whether a sample was flagged for a high percent chimera: 'pct_chimeras' > 5%.
	* **contamination**: Whether a sample was flagged for high contamination: 'mean_AB_snp_biallelic' > 0.015.
	* **low_coverage**: Whether a sample was flagged for low Chromosome 20 coverage: 'chr20_mean_dp' < 10.
	* **sample_qc_metrics**: Whether a sample was flagged for failing a Hail's `vds.sample_qc` module metric on bi-allelic variants: 'n_singleton' > 5000, 'r_het_hom_var' > 10, 'bases_dp_over_1' < 5e7, or 'bases_dp_over_20' < 4e7.
	* **low_adj_callrate**: Whether a sample was flagged for low callrate using high-quality genotypes (GQ >= 20; DP >= 10; and AB >= 0.2 for het calls): 'sample_qc_mt_callrate_adj' < 0.8.
	* **failed_fingerprinting**: Whether a sample was flagged for failing fingerprinting.
	* **sex_aneuploidy**: Whether a sample was flagged for sex aneuploidy.
	* **ambiguous_sex**: Whether a sample was flagged for ambiguous sex.
	* **hard_filters**: The sample's set of hard filters it failed.
	* **hard_filtered**: Whether a sample failed hard filtering.
	* **regressed_filters**: Sample QC metric outlier filter information for the per platform ancestry PC regression/residual filtering approach.
		* **qc_metrics_residuals**: Sample QC metric residuals after regressing out ancestry PCs per platform.
			* **r_ti_tv_residual**: Transition : transversion (TiTv) ratio residual.
			* **r_het_hom_var_residual**: Heterozygous : homozygous ratio residual.
			* **r_insertion_deletion_residual**: Insertion : deletion ratio residual.
			* **r_ti_tv_singleton_residual**: Singleton transition : transversion (TiTv) ratio residual.
		* **qc_metrics_fail**: Struct of booleans for each sample QC metric indicating whether a sample was found to be an outlier after regressing out ancestry PCs per platform.
			* **fail_r_ti_tv_residual**: Whether a sample failed the transition : transversion (TiTv) ratio residual filter.
			* **fail_r_insertion_deletion_residual**: Whether a sample failed the insertion : deletion ratio residual filter.
			* **fail_r_het_hom_var_residual**: Whether a sample failed the heterozygous : homozygous ratio residual filter.
			* **fail_r_ti_tv_singleton_residual**: Whether a sample failed the singleton transition : transversion (TiTv) ratio residual filter.
		* **qc_metrics_filters**: Set of all sample QC metrics for which a sample was found to be an outlier after regressing out the ancestry assignment PCs per platform.
	* **nearest_neighbor_filters**: Sample QC metric outlier filter information for the nearest neighbors filtering approach.
		* **qc_metrics_stats**: Sample QC metric statistics for a sample's nearest neighbors in the same platform.
			* **r_ti_tv**: Transition : transversion (TiTv) ratio.
				* **median**: Median transition : transversion (TiTv) ratio of nearest neighbors.
				* **mad**: MAD of transition : transversion (TiTv) ratio of nearest neighbors.
				* **lower**: Lower cutoff threshold of transition : transversion (TiTv) ratio of nearest neighbors.
				* **upper**: Upper cutoff threshold of transition : transversion (TiTv) ratio of nearest neighbors.
			* **r_ti_tv_singleton**: Singleton transition : transversion (TiTv) ratio.
				* **median**: Median singleton transition : transversion (TiTv) ratio of nearest neighbors.
				* **mad**: MAD of singleton transition : transversion (TiTv) ratio of nearest neighbors.
				* **lower**: Lower cutoff threshold of singleton transition : transversion (TiTv) ratio of nearest neighbors.
				* **upper**: Upper cutoff threshold of singleton transition : transversion (TiTv) ratio of nearest neighbors.
			* **r_het_hom_var**: Heterozygous : homozygous ratio.
				* **median**: Median heterozygous : homozygous ratio of nearest neighbors.
				* **mad**: MAD of heterozygous : homozygous ratio of nearest neighbors.
				* **lower**: Lower cutoff threshold of heterozygous : homozygous ratio of nearest neighbors.
				* **upper**: Upper cutoff threshold of heterozygous : homozygous ratio of nearest neighbors.
			* **r_insertion_deletion**: Insertion : deletion ratio.
				* **median**: Median insertion : deletion ratio of nearest neighbors.
				* **mad**: MAD of insertion : deletion ratio of nearest neighbors.
				* **lower**: Lower cutoff threshold of insertion : deletion ratio of nearest neighbors.
				* **upper**: Upper cutoff threshold of insertion : deletion ratio of nearest neighbors.
		* **qc_metrics_fail**: Struct of booleans for each sample QC metric indicating whether a sample was found to be an outlier in the per platform nearest neighbors approach.
			* **fail_r_het_hom_var**: Whether a sample failed the heterozygous : homozygous ratio filter in the per platform nearest neighbors approach.
			* **fail_r_ti_tv_singleton**: Whether a sample failed the singleton transition : transversion (TiTv) ratio filter in the per platform nearest neighbors approach.
			* **fail_r_ti_tv**: Whether a sample failed the transition : transversion (TiTv) ratio filter in the per platform nearest neighbors approach.
			* **fail_r_insertion_deletion**: Whether a sample failed the insertion : deletion ratio filter in the per platform nearest neighbors approach.
		* **qc_metrics_filters**: Set of all sample QC metrics for which a sample was found to be an outlier in the nearest neighbors filtering approach.
	* **qc_metrics_fail**: Struct of booleans for each sample QC metric indicating whether a sample was found to be an outlier in 'nearest_neighbor_filters' and/or (depends on value of 'ensemble_combination_operator' global) 'regressed_filters'.
		* **r_ti_tv**: Whether a sample failed the transition : transversion (TiTv) ratio ensemble filter.
		* **r_insertion_deletion**: Whether a sample failed the insertion : deletion ratio ensemble filter.
		* **r_het_hom_var**: Whether a sample failed the heterozygous : homozygous ratio ensemble filter.
		* **r_ti_tv_singleton**: Whether a sample failed the singleton transition : transversion (TiTv) ratio ensemble filter.
	* **qc_metrics_filters**: Combined set of all sample QC metrics for which a sample was found to be an outlier in 'nearest_neighbor_filters' and/or (depends on value of 'ensemble_combination_operator' global) 'regressed_filters'.
	* **outlier_filtered**: Whether a sample was filtered as an outlier for QC metrics ('qc_metrics_filters' is not empty).
	* **relatedness_filters**: Struct of relatedness filters. Any sample that is hard-filtered will have a missing value for these annotations.
		* **related**: Whether a sample was filtered for second-degree (or closer) relatedness in the ancestry inference PCA.
		* **duplicate_or_twin**: Whether a sample has a duplicate or twin among all samples that are not hard-filtered.
		* **parent_child**: Whether a sample has a parent or child among all samples that are not hard-filtered.
		* **sibling**: Whether a sample has a sibling among all samples that are not hard-filtered.
	* **release_relatedness_filters**: Struct of release sample relatedness filters. Any sample that is hard-filtered or outlier-filtered will have a missing value for these annotations.
		* **related**: Whether a sample was filtered for second-degree (or closer) relatedness in the final release.
		* **duplicate_or_twin**: Whether a release filtered sample has a duplicate or twin among all samples that are not hard-filtered or outlier-filtered.
		* **parent_child**: Whether a release filtered sample has a parent or child among all samples that are not hard-filtered or outlier-filtered.
		* **sibling**: Whether a release filtered sample has a sibling among all samples that are not hard-filtered or outlier-filtered.
	* **control**: Whether a sample is used as a control sample in variant QC assesment. This includes a synthetic-diploid (syndip) sample (a mixture of DNA from two haploid CHM cell lines) and two NA12878 samples.

* **high_quality**: Whether a sample is considered high-quality (no hard filter or sample QC metric outlier filter flags).

* **release**: Whether a sample meets criteria for release (Permission to release, high-quality, not to be excluded, and not related to any other sample in the release).
