"""Script containing generic resources."""

import logging
import os
from typing import List, Optional, Set, Union

import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VariantDatasetResource,
    VersionedVariantDatasetResource,
)
from gnomad.utils.file_utils import file_exists
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import _split_and_filter_variant_data_for_loading
from gnomad_qc.v5.resources.constants import (
    AOU_GENOMIC_METRICS_PATH,
    AOU_LOW_QUALITY_PATH,
    AOU_WGS_BUCKET,
    BATCH_BUCKET,
    BATCH_READ_ONLY_BUCKET,
    BATCH_TMP_BUCKET,
    CURRENT_AOU_VERSION,
    CURRENT_VERSION,
    GNOMAD_BUCKET,
    GNOMAD_TMP_BUCKET,
    WORKSPACE_BUCKET,
)

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)

# Constants for environment validation.
_SAMPLE_DATA_ENVIRONMENTS = frozenset({"rwb", "batch"})
_ALL_ENVIRONMENTS = frozenset({"rwb", "batch", "dataproc"})


def _file_exists_for_env(path: str, environment: str) -> bool:
    """Check if a file exists, handling permission errors in batch mode.

    When running with the batch ServiceBackend, ``file_exists`` executes
    locally on the driver via hailtop's GCS client, which may lack access to
    Terra workspace buckets that batch workers can reach.  In that case, assume
    the resource already exists so that downstream ``read_expression`` calls
    (which run on batch) succeed.

    :param path: GCS path to check.
    :param environment: Compute environment (``"rwb"``, ``"batch"``, or
        ``"dataproc"``).
    :return: True if the file exists (or existence cannot be verified in batch
        mode), False otherwise.
    """
    try:
        return file_exists(path)
    except Exception:
        if environment != "batch":
            raise
        logger.warning(
            "Unable to verify existence of '%s' from the local driver in "
            "batch mode. Assuming the resource exists; use --overwrite to "
            "force recomputation.",
            path,
        )
        return True


def _validate_environment(
    environment: str,
    allowed: frozenset = _ALL_ENVIRONMENTS,
) -> None:
    """
    Raise ValueError if `environment` is not in `allowed`.

    :param environment: Environment string to validate.
    :param allowed: Set of permitted environment strings. Default is
        `_ALL_ENVIRONMENTS`.
    """
    if environment not in allowed:
        raise ValueError(
            f"Environment '{environment}' is not allowed for this resource. "
            f"Must be one of: {sorted(allowed)}"
        )


def _get_base_bucket(environment: str = "batch", read_only: bool = False) -> str:
    """
    Return the top-level GCS bucket for the given environment.

    :param environment: Environment to use. Default is "batch". Must be one of "rwb", "batch", or
        "dataproc".
    :param read_only: If True and environment is "batch", return the read-only
        bucket within the AoU authorization Domain instead of the primary batch bucket.
        Default is False.
    :return: Bucket name string (without gs:// prefix).
    """
    _validate_environment(environment)
    if read_only and environment != "batch":
        raise ValueError(
            f"read_only=True is only supported when environment is 'batch', "
            f"got '{environment}'."
        )
    if environment == "rwb":
        return WORKSPACE_BUCKET
    elif environment == "batch":
        return BATCH_READ_ONLY_BUCKET if read_only else BATCH_BUCKET
    else:
        return GNOMAD_BUCKET


# v5 DRAGEN TGP test VDS.
dragen_tgp_vds = VariantDatasetResource(
    "gs://gnomad/v5.0/testing/genomes/dragen_tgp_v5.0_test.vds"
)

_aou_genotypes = {
    "8": VariantDatasetResource(
        "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
    )
}

aou_genotypes = VersionedVariantDatasetResource(
    CURRENT_AOU_VERSION,
    _aou_genotypes,
)

# This test VDS was initially created to count the number of singletons in AoU and compare the results to those
# provided in the sample QC metrics. It contains 10 randomly chosen samples based on if they
# `fail_singleton_residual` (5 passed, 5 failed) from the AoU sample QC data and is also used for testing our code.
# 'hl.vds.filter_samples()' is used to filter the VDS to these samples, hence it includes all the variants present in
# these samples on all the contigs as the original VDS.
aou_test_dataset = VariantDatasetResource(
    f"gs://{WORKSPACE_BUCKET}/v5.0/hard_filtering/10sample_for_singleton_test.vds"
)


def _init_hail(
    log_name: str,
    environment: str = "batch",
    billing_project: Optional[str] = None,
    tmp_dir_days: Optional[int] = 4,
    **kwargs,
) -> None:
    """
    Initialize Hail with environment-appropriate settings and set GRCh38 as default reference.

    :param log_name: Base name for the log file (without path or extension).
    :param environment: Compute environment. One of "rwb", "batch", or "dataproc". Default is "batch".
    :param billing_project: GCP billing project for requester-pays buckets (batch only).
        Default is None. When None, uses "broad-mpg-gnomad".
    :param tmp_dir_days: Retention days for the tmp directory passed to qc_temp_prefix.
        Must be None, 4, or 30. Default is 4.
    :param kwargs: Additional keyword arguments forwarded to hl.init() in all
        environments. None values are silently dropped, so optional params (e.g.
        batch resource params from :func:`_get_batch_resource_kwargs`, or
        ``spark_conf`` for dataproc) can be passed unconditionally.
    """
    if environment == "rwb":
        log = f"/home/jupyter/workspaces/gnomadproduction/{log_name}.log"
    elif environment == "batch":
        log = f"/tmp/{log_name}.log"
    else:
        log = get_logging_path(
            log_name, environment=environment, tmp_dir_days=tmp_dir_days
        )
    init_kwargs = {
        "log": log,
        "tmp_dir": qc_temp_prefix(environment=environment, days=tmp_dir_days),
        **{k: v for k, v in kwargs.items() if v is not None},
    }
    if environment == "batch":
        os.environ["HAIL_BATCH_BILLING_PROJECT"] = (
            billing_project or "gnomad-production"
        )
        os.environ["HAIL_BATCH_REMOTE_TMPDIR"] = f"gs://{BATCH_TMP_BUCKET}"
        init_kwargs.update(
            {
                "backend": "batch",
                "gcs_requester_pays_configuration": "broad-mpg-gnomad",
                "regions": ["us-central1"],
            }
        )
    hl.init(**init_kwargs, default_reference="GRCh38")


def _init_hail_local_spark(
    log_name: str,
    tmp_dir: Optional[str] = None,
    driver_memory: str = "10g",
    local_cores: int = 2,
) -> None:
    """
    Initialize Hail in single-node local-Spark mode.

    Intended for use inside a Hail Batch BashJob container so each job runs
    its Hail Query work locally on its own cores (no recursive QoB).

    :param log_name: Base name for the log file (without path or extension).
    :param tmp_dir: Optional temporary directory prefix for Hail.
    :param driver_memory: JVM heap size for the local Spark driver.
        Default ``"10g"``.
    :param local_cores: Number of Spark threads. Fewer threads means less
        concurrent memory pressure. Default ``2``.
    """
    init_kwargs = {
        "log": f"/tmp/{log_name}.log",
        "default_reference": "GRCh38",
        "gcs_requester_pays_configuration": "broad-mpg-gnomad",
        "master": f"local[{local_cores}]",
        "spark_conf": {
            "spark.driver.memory": driver_memory,
        },
    }
    if tmp_dir is not None:
        init_kwargs["tmp_dir"] = tmp_dir
    hl.init(**init_kwargs)


# hl.default_reference("GRCh38")


def qc_temp_prefix(
    version: str = CURRENT_VERSION,
    environment: str = "dataproc",
    days: Optional[int] = None,
) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return.
    :param environment: Compute environment, either 'dataproc','rwb', or 'batch'. Default is 'dataproc'.
    :param days: Number of days to keep temporary data. Default is None.
    :return: Path to bucket with temporary QC data.
    """
    if days not in [None, 4, 30]:
        raise ValueError("Days must be either None, 4, or 30.")

    if environment == "rwb":
        env_bucket = (
            f"{WORKSPACE_BUCKET}/tmp{f'/{days}_day' if days is not None else ''}"
        )
    elif environment in ("dataproc", "batch"):
        tmp_bucket = (
            GNOMAD_TMP_BUCKET if environment == "dataproc" else BATCH_TMP_BUCKET
        )
        env_bucket = f"{tmp_bucket}{f'-{days}day' if days is not None else ''}"
    else:
        raise ValueError(
            f"Environment {environment} not recognized. Choose 'rwb', 'dataproc', or 'batch'."
        )

    return f"gs://{env_bucket}/gnomad.genomes.v{version}.qc_data/"


def get_checkpoint_path(
    name: str,
    version: str = CURRENT_VERSION,
    mt: bool = False,
    environment: str = "dataproc",
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable.
    :param version: Version of annotation path to return.
    :param bool mt: Whether path is for a MatrixTable, default is False.
    :param environment: Compute environment, either 'dataproc','rwb', or 'batch'. Default is 'dataproc'.
    :return: Output checkpoint path.
    """
    return f'{qc_temp_prefix(version, environment)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(
    name: str,
    version: str = CURRENT_VERSION,
    environment: str = "dataproc",
    tmp_dir_days: Optional[int] = None,
) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :param version: Version of annotation path to return.
    :param environment: Compute environment, 'dataproc', 'rwb', or 'batch'. Default is 'dataproc'.
    :param tmp_dir_days: Number of days to keep temporary data. Default is None.
    :return: Output log path.
    """
    return f"{qc_temp_prefix(version, environment, tmp_dir_days)}{name}.log"


_BATCH_RESOURCE_PARAMS = [
    "app_name",
    "driver_cores",
    "driver_memory",
    "worker_cores",
    "worker_memory",
]


def _get_batch_resource_kwargs(args) -> dict:
    """
    Extract optional Hail Batch resource parameters from parsed args, omitting None values.

    Intended for use with scripts that expose ``--app-name``, ``--driver-cores``,
    ``--driver-memory``, ``--worker-cores``, and ``--worker-memory`` arguments. The
    result can be unpacked directly into :func:`_init_hail`.

    :param args: Parsed command-line arguments.
    :return: Dict of non-None batch resource kwargs.
    """
    return {
        p: getattr(args, p)
        for p in _BATCH_RESOURCE_PARAMS
        if getattr(args, p, None) is not None
    }


def get_aou_vds(
    split: bool = False,
    remove_hard_filtered_samples: bool = True,
    high_quality_only: bool = False,
    release_only: bool = False,
    filter_samples: Optional[Union[List[str], hl.Table]] = None,
    test: bool = False,
    n_partitions_on_read: Optional[int] = None,
    filter_partitions: Optional[List[int]] = None,
    repartition_after_filter: Optional[int] = None,
    chrom: Optional[Union[str, List[str], Set[str]]] = None,
    autosomes_only: bool = False,
    sex_chr_only: bool = False,
    filter_variant_ht: Optional[hl.Table] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    remove_dead_alleles: bool = True,
    annotate_meta: bool = False,
    entries_to_keep: Optional[List[str]] = None,
    checkpoint_variant_data: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
    add_project_prefix: bool = False,
    environment: str = "batch",
) -> hl.vds.VariantDataset:
    """
    Load the AOU VDS.

    :param split: Whether to split multi-allelic variants in the VDS. Note: this will perform a split on the VDS
        rather than grab an already split VDS. Default is False.
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard
        filters (only relevant after hard filtering is complete). Default is True.
    :param high_quality_only: Whether to filter the VDS to only high quality samples. Default is False.
    :param release_only: Whether to filter the VDS to only samples available for
        release (can only be used if metadata is present).
    :param filter_samples: Optional samples to filter the VDS to. Can be a list of sample IDs or a Table with sample IDs.
        If `filter_samples` contains sample IDs with collisions with gnomAD samples,
        `add_project_prefix` must be set to True to filter properly. Default is None.
    :param test: Whether to load the test VDS instead of the full VDS. The test VDS includes 10 samples selected from the full dataset for testing purposes. Default is False.
    :param n_partitions_on_read: Optional number of partitions to repartition the VDS
        into on read. This increases parallelism so individual tasks are small enough
        to complete before preemption on GCP preemptible VMs. Default is None.
    :param filter_partitions: Optional argument to filter the VDS to a list of specific partitions.
    :param repartition_after_filter: Optional number of partitions to repartition the
        filtered VDS into after applying ``filter_partitions``. Useful for testing
        when original partitions are too large for parallelism. Requires
        ``filter_partitions`` to be set. Default is None.
    :param chrom: Optional argument to filter the VDS to a specific chromosome(s).
    :param autosomes_only: Whether to include only autosomes. Default is False.
    :param sex_chr_only: Whether to include only sex chromosomes. Default is False.
    :param filter_variant_ht: Optional argument to filter the VDS to a specific set of variants. Only supported when splitting the VDS.
    :param filter_intervals: Optional argument to filter the VDS to specific intervals.
    :param split_reference_blocks: Whether to split the reference data at the edges of the intervals defined by `filter_intervals`. Default is True.
    :param remove_dead_alleles: Whether to remove dead alleles when removing samples. Default is True.
    :param annotate_meta: Whether to annotate the VDS with the sample QC metadata. Default is False.
    :param entries_to_keep: Optional list of entries to keep in the variant data. If splitting the VDS, use the global entries (e.g. 'GT') instead of the local entries (e.g. 'LGT') to keep.
    :param checkpoint_variant_data: Whether to checkpoint the variant data MT after splitting and filtering. Default is False.
    :param naive_coalesce_partitions: Optional number of partitions to coalesce the VDS to. Default is None.
    :param add_project_prefix: Whether to prefix sample IDs (e.g., ``'aou_'``) for samples that exist in multiple projects to avoid ID collisions. Default is False.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb" or "batch".
    :return: AoU v8 VDS.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    aou_v8_resource = aou_test_dataset if test else aou_genotypes

    if isinstance(chrom, str):
        chrom = [chrom]

    if n_partitions_on_read and (
        chrom or autosomes_only or sex_chr_only or filter_intervals
    ):
        logger.warning(
            "n_partitions_on_read is set with chromosome or interval filtering. "
            "Partition boundaries will be computed over the entire genome before "
            "filtering, which may result in many empty partitions."
        )

    vds = aou_v8_resource.vds(
        read_args=(
            {"n_partitions": n_partitions_on_read} if n_partitions_on_read else None
        )
    )

    if autosomes_only and sex_chr_only:
        raise ValueError(
            "Only one of 'autosomes_only' or 'sex_chr_only' can be set to True."
        )

    # Apply chromosome filtering.
    if sex_chr_only:
        vds = hl.vds.filter_chromosomes(vds, keep=["chrX", "chrY"])
    elif autosomes_only:
        vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)
    elif chrom and len(chrom) > 0:
        logger.info("Filtering to chromosome(s) %s...", chrom)
        vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    # --- Row filtering first, before any sample operations. ---
    # Narrowing rows early means downstream sample filtering (especially
    # filter_samples with remove_dead_alleles) operates on less data.

    # Apply partition filtering.
    if filter_partitions and len(filter_partitions) > 0:
        logger.info("Filtering to %s partitions...", len(filter_partitions))
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(filter_partitions),
            vds.variant_data._filter_partitions(filter_partitions),
        )

    # Repartition filtered VDS into smaller partitions for parallelism.
    if repartition_after_filter:
        if not filter_partitions:
            raise ValueError(
                "'repartition_after_filter' requires 'filter_partitions' to be set."
            )
        logger.info(
            "Repartitioning filtered VDS into %s partitions...",
            repartition_after_filter,
        )
        vds = hl.vds.VariantDataset(
            vds.reference_data.repartition(repartition_after_filter),
            vds.variant_data.repartition(repartition_after_filter),
        )

    # Remove the chr4 site with excessive numbers of alleles (n=22233) to
    # avoid memory issues with `split_multi`.
    logger.info("Dropping excessively multi-allelic site at chr4:12237652...")
    vds = hl.vds.filter_intervals(
        vds,
        [hl.parse_locus_interval("chr4:12237652-12237653", reference_genome="GRCh38")],
        keep=False,
    )

    # Apply interval filtering.
    if filter_intervals and len(filter_intervals) > 0:
        logger.info("Filtering to %s intervals...", len(filter_intervals))
        if isinstance(filter_intervals[0], str):
            filter_intervals = [
                hl.parse_locus_interval(x, reference_genome="GRCh38")
                for x in filter_intervals
            ]
        vds = hl.vds.filter_intervals(
            vds, filter_intervals, split_reference_blocks=split_reference_blocks
        )

    if naive_coalesce_partitions:
        vds = hl.vds.VariantDataset(
            vds.reference_data.naive_coalesce(naive_coalesce_partitions),
            vds.variant_data.naive_coalesce(naive_coalesce_partitions),
        )

    # --- Sample filtering, now on the row-narrowed VDS. ---

    # Count initial number of samples.
    n_samples_before = vds.variant_data.count_cols()

    # Remove samples that should have been excluded from the AoU v8 release
    # and samples with non-XX/XY ploidies.
    hard_filtered_samples_ht = None
    if remove_hard_filtered_samples:
        from gnomad_qc.v5.resources.sample_qc import get_hard_filtered_samples

        logger.info("Removing hard filtered samples from AoU VDS...")
        hard_filtered_samples_ht = get_hard_filtered_samples(
            environment=environment
        ).ht()
    s_to_exclude = list(
        hl.eval(
            get_samples_to_exclude(hard_filtered_samples_ht, environment=environment)
        )
    )
    vds = hl.vds.filter_samples(
        vds, s_to_exclude, keep=False, remove_dead_alleles=remove_dead_alleles
    )

    # Report final sample exclusion count.
    n_samples_after = vds.variant_data.count_cols()
    logger.info("Removed %d samples from VDS.", n_samples_before - n_samples_after)

    vmt = vds.variant_data
    rmt = vds.reference_data

    if release_only or high_quality_only or annotate_meta or add_project_prefix:
        # Import here to avoid circular imports.
        from gnomad_qc.v5.resources.meta import get_sample_id_collisions, meta

        meta_ht = meta(data_type="genomes", environment=environment).ht()

        logger.warning(
            "Adding 'aou_' prefix to samples that had ID collisions with gnomAD samples..."
        )
        sample_collisions = get_sample_id_collisions(environment=environment).ht()
        vmt = add_project_prefix_to_sample_collisions(
            t=vmt, sample_collisions=sample_collisions, project="aou"
        )
        rmt = add_project_prefix_to_sample_collisions(
            t=rmt, sample_collisions=sample_collisions, project="aou"
        )

        if annotate_meta:
            logger.info("Annotating VDS variant_data with metadata...")
            vmt = vmt.annotate_cols(meta=meta_ht[vmt.col_key])

    if filter_variant_ht is not None and split is False:
        raise ValueError(
            "Filtering to a specific set of variants is only supported when splitting"
            " the VDS."
        )
    if split:
        vmt = _split_and_filter_variant_data_for_loading(
            vmt, filter_variant_ht, entries_to_keep, checkpoint_variant_data
        )
        # Don't need to filter entries again if already done in split function
        entries_filtered_during_split = entries_to_keep is not None
    else:
        entries_filtered_during_split = False

    if entries_to_keep is not None and not entries_filtered_during_split:
        vmt = vmt.select_entries(*entries_to_keep)

    if checkpoint_variant_data:
        vmt = vmt.checkpoint(new_temp_file("vds_loading.variant_data", "mt"))

    vds = hl.vds.VariantDataset(rmt, vmt)

    if release_only:
        logger.info("Filtering VDS to release samples only...")
        filter_expr = meta_ht.release
        vds = hl.vds.filter_samples(vds, meta_ht.filter(filter_expr))

    if high_quality_only:
        logger.info("Filtering VDS to high quality samples only...")
        filter_expr = meta_ht.high_quality
        vds = hl.vds.filter_samples(vds, meta_ht.filter(filter_expr))

    if filter_samples:
        logger.info(
            "Filtering to %s samples...",
            (
                len(filter_samples)
                if isinstance(filter_samples, list)
                else filter_samples.count()
            ),
        )
        logger.warning(
            "This filter will not work correctly if `filter_samples` contains sample"
            " IDs with collisions with gnomAD samples and `add_project_prefix` is not"
            " set."
        )
        vds = hl.vds.filter_samples(vds, filter_samples)

    return vds


def get_gnomad_v5_genomes_vds(
    split: bool = False,
    release_only: bool = False,
    consent_drop_only: bool = False,
    annotate_meta: bool = False,
    test: bool = False,
    filter_partitions: Optional[List[int]] = None,
    chrom: Optional[Union[str, List[str], Set[str]]] = None,
    autosomes_only: bool = False,
    sex_chr_only: bool = False,
    filter_variant_ht: Optional[hl.Table] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    entries_to_keep: Optional[List[str]] = None,
    annotate_het_non_ref: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
    filter_samples_ht: Optional[hl.Table] = None,
) -> hl.vds.VariantDataset:
    """
    Get gnomAD v5 genomes VariantDataset with desired filtering and metadata annotations.

    :param split: Perform split on VDS - Note: this will perform a split on the VDS
        rather than grab an already split VDS.
    :param release_only: Whether to filter the VDS to only samples available for
        v5 release (distinct from v4 release due to samples to drop for consent reasons).
        Requires that v5 sample metadata has been computed.
    :param consent_drop_only: Whether to filter the VDS to only consent drop samples.
    :param annotate_meta: Whether to add v4 genomes metadata to VDS variant_data in
        'meta' column.
    :param test: Whether to use the test VDS instead of the full v4 genomes VDS.
    :param filter_partitions: Optional argument to filter the VDS to specific partitions
        in the provided list.
    :param chrom: Optional argument to filter the VDS to a specific chromosome(s).
    :param autosomes_only: Whether to filter the VDS to autosomes only. Default is
        False.
    :param sex_chr_only: Whether to filter the VDS to sex chromosomes only. Default is
        False.
    :param filter_variant_ht: Optional argument to filter the VDS to a specific set of
        variants. Only supported when splitting the VDS.
    :param filter_intervals: Optional argument to filter the VDS to specific intervals.
    :param split_reference_blocks: Whether to split the reference data at the edges of
        the intervals defined by `filter_intervals`. Default is True.
    :param entries_to_keep: Optional argument to keep only specific entries in the
        returned VDS. If splitting the VDS, use the global entries (e.g. 'GT') instead
        of the local entries (e.g. 'LGT') to keep.
    :param annotate_het_non_ref: Whether to annotate non reference heterozygotes (as
        '_het_non_ref') to the variant data. Default is False.
    :param naive_coalesce_partitions: Optional argument to coalesce the VDS to a
        specific number of partitions using naive coalesce.
    :param filter_samples_ht: Optional Table of samples to filter the VDS to.
    :return: gnomAD v4 genomes VariantDataset with chosen annotations and filters.
    """
    # Import v3 basics and v4 meta here to avoid circular imports.
    from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds
    from gnomad_qc.v4.resources.meta import meta as v4_meta

    vds = get_gnomad_v3_vds(
        split=split,
        # False because v3 hard filtered samples HT no longer exists.
        remove_hard_filtered_samples=False,
        release_only=False,
        samples_meta=False,
        test=test,
        filter_partitions=filter_partitions,
        chrom=chrom,
        autosomes_only=autosomes_only,
        sex_chr_only=sex_chr_only,
        filter_variant_ht=filter_variant_ht,
        filter_intervals=filter_intervals,
        split_reference_blocks=split_reference_blocks,
        entries_to_keep=entries_to_keep,
        annotate_het_non_ref=annotate_het_non_ref,
        naive_coalesce_partitions=naive_coalesce_partitions,
        filter_samples_ht=filter_samples_ht,
    )

    # NOTE: Not using `project_meta` here to allow all gnomAD steps to be run
    # in Dataproc.
    meta_ht = v4_meta(data_type="genomes").ht()
    meta_ht = meta_ht.annotate(
        consent_drop=hl.is_defined(meta_ht.project_meta.research_project_key)
        & (
            (meta_ht.project_meta.research_project_key == "RP-1061")
            | (meta_ht.project_meta.research_project_key == "RP-1411")
        )
    )

    if release_only or consent_drop_only or annotate_meta:
        filter_expr = True
        if release_only:
            if not consent_drop_only:
                meta_ht = meta_ht.annotate(
                    release=hl.if_else(
                        meta_ht.consent_drop,
                        False,
                        meta_ht.release,
                    )
                )
            filter_expr &= meta_ht.release

        if consent_drop_only:
            filter_expr &= meta_ht.consent_drop

        if annotate_meta:
            vd = vds.variant_data
            vds = hl.vds.VariantDataset(
                vds.reference_data, vd.annotate_cols(meta=meta_ht[vd.col_key])
            )
        vds = hl.vds.filter_samples(vds, meta_ht.filter(filter_expr))

    return vds


aou_acaf_mt = MatrixTableResource(
    path=f"gs://{AOU_WGS_BUCKET}/acaf_threshold/splitMT/hail.mt"
)
"""
AoU v8 ACAF (Allele Count/Allele Frequency threshold) MatrixTable.

MatrixTable contains only variants with AF > 1% or AC > 100 in any genetic ancestry group.

See https://support.researchallofus.org/hc/en-us/articles/29475228181908-How-the-All-of-Us-Genomic-data-are-organized#01JJK0HH53FX9XQRDQ5HQFZW9B
and https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK
for more information.
"""


aou_exome_mt = MatrixTableResource(path=f"gs://{AOU_WGS_BUCKET}/exome/splitMT/hail.mt")
"""
AoU v8 Exome MatrixTable.

MatrixTable contains only variants in exons (with 15 bp padding on either side) as defined by GENCODE v42 basic.

See same links as above (in `acaf_mt`) for more information.
"""


def get_aou_failing_genomic_metrics_samples() -> Set[str]:
    """
    Import AoU genomic metrics and filter to samples that fail specific quality criteria, including low coverage and ambiguous sex ploidy.

    .. note::

        Samples with low mean coverage (<30x), genome coverage (<90% at 20x), All of Us Hereditary Disease Risk gene
        (AoUHDR) coverage (<95% at 20x), or aligned_q30_bases (<8e10) were expected to be excluded from the AoU
        callset. However, some such samples are still present. AoU is preparing to publish a known issue in their
        quality report related to this. This note will be updated with a link once the issue is published.

        In addition, we exclude samples with ambiguous sex ploidy (i.e., not "XX" or "XY") from the callset.

    :return: Set of sample IDs failing coverage filters or with non-XX-XY sex ploidies.
    """
    types = {
        "research_id": hl.tstr,
        "sample_source": hl.tstr,
        "site_id": hl.tstr,
        "sex_at_birth": hl.tstr,
        "dragen_sex_ploidy": hl.tstr,
        "mean_coverage": hl.tfloat,
        "genome_coverage": hl.tfloat,
        "aou_hdr_coverage": hl.tfloat,
        "dragen_contamination": hl.tfloat,
        "aligned_q30_bases": hl.tint64,
        "verify_bam_id2_contamination": hl.tfloat,
        "biosample_collection_date": hl.tstr,
    }
    ht = hl.import_table(AOU_GENOMIC_METRICS_PATH, types=types).key_by("research_id")

    low_cov_samples = ht.filter(
        (ht.mean_coverage < 30)
        | (ht.genome_coverage < 0.9)
        | (ht.aou_hdr_coverage < 0.95)
        | (ht.aligned_q30_bases < 8e10)  # AoU threshold; no v8 samples failed this.
    ).select()
    logger.info("%s samples with low coverage...", low_cov_samples.count())

    ambiguous_sex_samples = ht.filter(
        (ht.dragen_sex_ploidy != "XX") & (ht.dragen_sex_ploidy != "XY")
        | hl.is_missing(ht.dragen_sex_ploidy)
    ).select()
    logger.info(
        "%s samples with ambiguous sex ploidy...", ambiguous_sex_samples.count()
    )

    ht = low_cov_samples.union(ambiguous_sex_samples).distinct()
    return ht.aggregate(hl.agg.collect_as_set(ht.research_id))


def get_samples_to_exclude(
    filter_samples: Optional[Union[List[str], hl.Table]] = None,
    overwrite: bool = False,
    environment: str = "batch",
) -> hl.expr.SetExpression:
    """
    Get set of AoU sample IDs to exclude.

    .. note::

        If `filter_samples` is a Hail Table, it must contain a field named 's' with sample IDs.

    :param filter_samples: Optional additional samples to remove. Can be a list of sample IDs or a Table with sample IDs.
    :param overwrite: Whether to overwrite the existing `samples_to_exclude` resource. Default is False.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: SetExpression containing IDs of samples to exclude from v5 analysis.
    """
    _validate_environment(environment, _SAMPLE_DATA_ENVIRONMENTS)
    # Import here to avoid circular imports.
    from gnomad_qc.v5.resources.meta import (
        get_failing_metrics_samples,
        get_low_quality_samples,
        get_samples_to_exclude_resource,
    )

    lq_resource = get_low_quality_samples(environment=environment)
    fm_resource = get_failing_metrics_samples(environment=environment)
    ste_resource = get_samples_to_exclude_resource(environment=environment)

    if overwrite or not _file_exists_for_env(ste_resource.path, environment):

        if not _file_exists_for_env(lq_resource.path, environment):
            # Load samples flagged in AoU Known Issues #1.
            logger.info("Removing 3 known low-quality samples (Known Issues #1)...")
            low_quality_ht = hl.import_table(AOU_LOW_QUALITY_PATH).key_by("research_id")
            low_quality_sample_ids = low_quality_ht.aggregate(
                hl.agg.collect_as_set(low_quality_ht.research_id)
            )
            hl.experimental.write_expression(
                hl.set(low_quality_sample_ids), lq_resource.path
            )
        if not _file_exists_for_env(fm_resource.path, environment):
            # Load and count samples failing genomic metrics filters.
            failing_genomic_metrics_samples = get_aou_failing_genomic_metrics_samples()
            logger.info(
                "Removing %d samples failing genomic QC (low coverage or ambiguous sex)...",
                len(failing_genomic_metrics_samples),
            )
            hl.experimental.write_expression(
                hl.set(failing_genomic_metrics_samples), fm_resource.path
            )

        # Union all samples to exclude and write out.
        low_quality_sample_ids = hl.experimental.read_expression(lq_resource.path)
        failing_genomic_metrics_samples = hl.experimental.read_expression(
            fm_resource.path
        )
        s_to_exclude = low_quality_sample_ids.union(failing_genomic_metrics_samples)
        hl.experimental.write_expression(
            s_to_exclude, ste_resource.path, overwrite=True
        )

    s_to_exclude = hl.experimental.read_expression(ste_resource.path)

    if filter_samples is None:
        additional_samples = hl.empty_set(hl.tstr)
    elif isinstance(filter_samples, hl.Table):
        if "s" not in filter_samples.row:
            raise ValueError(
                "Hail Table must contain a field named 's' with sample IDs."
            )
        additional_samples = filter_samples.aggregate(
            hl.agg.collect_as_set(filter_samples.s)
        )
    elif isinstance(filter_samples, list):
        additional_samples = hl.literal(set(filter_samples))
    else:
        raise ValueError(
            "`filter_samples` must be a list of sample IDs or a Hail Table with sample IDs."
        )

    return s_to_exclude.union(additional_samples)


def add_project_prefix_to_sample_collisions(
    t: Union[hl.Table, hl.MatrixTable],
    sample_collisions: hl.Table,
    project: Optional[str] = None,
    sample_id_field: str = "s",
) -> hl.Table:
    """
    Add project prefix to sample IDs that exist in multiple projects.

    :param t: Table/MatrixTable to add project prefix to sample IDs.
    :param sample_collisions: Table of sample IDs that exist in multiple projects.
    :param project: Optional project name to prepend to sample collisions. If not set, will use 'ht.project' annotation. Default is None.
    :param sample_id_field: Field name for sample IDs in the table.
    :return: Table with project prefix added to sample IDs.
    """
    logger.info(
        "Adding project prefix to sample IDs that exists in multiple projects..."
    )
    collisions = sample_collisions.aggregate(hl.agg.collect_as_set(sample_collisions.s))

    if project:
        prefix_expr = hl.literal(project)
    else:
        try:
            prefix_expr = t["project"]
        except KeyError:
            raise ValueError(
                "No project name provided and 't' does not contain a 'project' field."
            )

    renaming_expr = {
        f"{sample_id_field}": hl.if_else(
            hl.literal(collisions).contains(t[sample_id_field]),
            hl.delimit([prefix_expr, t[sample_id_field]], "_"),
            t[sample_id_field],
        )
    }
    if isinstance(t, hl.Table):
        return t.key_by(**renaming_expr)

    logger.warning("Input is a MatrixTable, rekeying columns...")
    return t.key_cols_by(**renaming_expr)
