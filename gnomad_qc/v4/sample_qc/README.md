# Sample QC

## Relatedness

We use [`cuKING`](https://github.com/populationgenomics/cuKING) for pruning related samples for v4. It involves the following steps:

### Preparing the production environment

Below, `$GCP_PROJECT` refers to the GCP project name for gnomAD production.

```sh
gcloud config set project $GCP_PROJECT
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
gcloud artifacts repositories create images --location=us-central1
```

Fetch a `cuKING` release and replace the GCP project ID:

```sh
curl -sSL https://github.com/populationgenomics/cuKING/archive/refs/tags/v1.0.0.tar.gz | tar xz -

cd cuKING-1.0.0

perl -pi -e "s/cpg-gnomad-production-27bb/$GCP_PROJECT/g" *
```

Build the Docker image:

```sh
gcloud builds submit --config cloudbuild.yaml .
```

This should take about half an hour.

Create a service account for running `cuKING`:

```sh
gcloud iam service-accounts create cuking
```

Grant the service account write access to the `gnomad` bucket:

```sh
gsutil iam ch serviceAccount:cuking@$GCP_PROJECT.iam.gserviceaccount.com:objectAdmin gs://gnomad
```

Create a VM instance template:

```sh
gcloud compute instance-templates create cuking-instance-template \
    --project=$GCP_PROJECT \
    --machine-type=a2-highgpu-1g \
    --network-interface=network=default,network-tier=PREMIUM,address="" \
    --metadata-from-file=startup-script=instance_startup_script.sh \
    --maintenance-policy=TERMINATE \
    --provisioning-model=STANDARD \
    --service-account=cuking@$GCP_PROJECT.iam.gserviceaccount.com \
    --scopes=https://www.googleapis.com/auth/cloud-platform \
    --accelerator=count=1,type=nvidia-tesla-a100 \
    --create-disk=auto-delete=yes,boot=yes,device-name=cuking-instance-template,image=projects/ubuntu-os-cloud/global/images/ubuntu-1804-bionic-v20220712,mode=rw,size=10,type=pd-balanced \
    --no-shielded-secure-boot \
    --shielded-vtpm \
    --shielded-integrity-monitoring \
    --reservation-affinity=any
```

### Prepare inputs

Start an autoscaling Hail Dataproc cluster with 20 secondary workers and submit `relatedness.py --prepare-inputs`. This should take about an hour.

### Run `cuKING`

This step doesn't need a Dataproc cluster. Run `relatedness.py --run-cuking` to submit a Cloud Batch job. The script will print the command that can be used to check progress. This should take about two hours.

### Convert outputs

`cuKING` writes Parquet files. To convert these to a standard Hail table, run `relatedness.py --create-relatedness-table` on a small Dataproc cluster. This should only take a minute.

### Determine samples to drop

To prune samples up to second degree relatedness, run `relatedness.py --compute-related-samples-to-drop` on a small Dataproc cluster. This should only take a few minutes.

### Inferring relationships for closely related samples

To create the annotations for distinguishing between siblings, parent-child pairs, etc., run `relatedness.py --create-relatedness-annotated-table` on a small Dataproc cluster. This should only take a minute.
