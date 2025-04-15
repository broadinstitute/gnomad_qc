"""Script containing generic resources."""

import logging
from typing import List, Optional, Set, Union

import hail as hl
from gnomad.resources.resource_utils import (
    VariantDatasetResource,
    VersionedVariantDatasetResource,
)
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import _split_and_filter_variant_data_for_loading
from gnomad_qc.v5.resources.constants import (
    CURRENT_AOU_VERSION,
    CURRENT_VERSION,
    WORKSPACE_BUCKET,
)

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)

# v5 DRAGEN TGP test VDS.
dragen_tgp_vds = VariantDatasetResource(
    "gs://gnomad/v5.0/testing/genomes/dragen_tgp_v5.0_test.vds"
)


def qc_temp_prefix(
    version: str = CURRENT_VERSION, environment: str = "dataproc"
) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return.
    :param environment: Compute environment, either 'dataproc' or 'rwb'. Defaults to 'dataproc'.
    :return: Path to bucket with temporary QC data.
    """
    if environment == "rwb":
        env_bucket = f"{WORKSPACE_BUCKET}/tmp"
    elif environment == "dataproc":
        env_bucket = "gnomad-tmp"
    else:
        raise ValueError(
            f"Environment {environment} not recognized. Choose 'rwb' or 'dataproc'."
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
    :param environment: Compute environment, either 'dataproc' or 'rwb'. Defaults to 'dataproc'.
    :return: Output checkpoint path.
    """
    return f'{qc_temp_prefix(version, environment)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(
    name: str, version: str = CURRENT_VERSION, environment: str = "dataproc"
) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :param version: Version of annotation path to return.
    :param environment: Compute environment, either 'dataproc' or 'rwb'. Defaults to 'dataproc'.
    :return: Output log path.
    """
    return f"{qc_temp_prefix(version, environment)}{name}.log"


_aou_genotypes = {
    "v8": VariantDatasetResource(
        "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
    )
}

aou_genotypes = VersionedVariantDatasetResource(
    CURRENT_AOU_VERSION,
    _aou_genotypes,
)


def get_aou_vds(
    split: bool = False,
    chrom: Optional[Union[str, List[str], Set[str]]] = None,
    autosomes_only: bool = False,
    sex_chr_only: bool = False,
    filter_samples: Optional[Union[List[str], hl.Table]] = None,
    filter_partitions: Optional[List[int]] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    filter_variant_ht: Optional[hl.Table] = None,
    entries_to_keep: Optional[List[str]] = None,
    checkpoint_variant_data: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
) -> hl.vds.VariantDataset:
    """
    Load the AOU VDS.

    :param split: Whether to split the VDS into separate datasets for each chromosome. Default is False.
    :param chrom: Chromosome(s) to filter the VDS to. Can be a single chromosome or a list of chromosomes.
    :param autosomes_only: Whether to include only autosomes.
    :param sex_chr_only: whether to include only sex chromosomes.
    :param filter_samples: List of samples to filter the VDS to. Can be a list of sample IDs or a Table with sample IDs.
    :param filter_intervals: List of intervals to filter the VDS.
    :param split_reference_blocks: Whether to split the reference blocks. Default is True.
    :param filter_variant_ht: Table to filter the variant data.
    :param filter_partitions: List of partitions to filter the VDS.
    :param entries_to_keep: List of entries to keep in the variant data.
    :param checkpoint_variant_data: Whether to checkpoint the variant data. Default is False.
    :param naive_coalesce_partitions: Number of partitions to coalesce the VDS to. Default is None.
    :return: The AOU VDS.
    """
    aou_v8_resource = aou_genotypes
    vds = aou_v8_resource.vds()

    if isinstance(chrom, str):
        chrom = [chrom]

    if autosomes_only and sex_chr_only:
        raise ValueError(
            "Only one of 'autosomes_only' or 'sex_chr_only' can be set to True."
        )

        # Apply chromosome filtering
    if sex_chr_only:
        vds = hl.vds.filter_chromosomes(vds, keep=["chrX", "chrY"])
    elif autosomes_only:
        vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)
    elif chrom and len(chrom) > 0:
        logger.info("Filtering to chromosome(s) %s...", chrom)
        vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    # Apply sample filtering
    if filter_samples:
        logger.info(
            "Filtering to %s samples...",
            (
                len(filter_samples)
                if isinstance(filter_samples, list)
                else filter_samples.count()
            ),
        )
        vds = hl.vds.filter_samples(vds, filter_samples)

    if naive_coalesce_partitions:
        vds = hl.vds.VariantDataset(
            vds.reference_data.naive_coalesce(naive_coalesce_partitions),
            vds.variant_data.naive_coalesce(naive_coalesce_partitions),
        )

    # Apply interval filtering
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

        # Apply partition filtering
    if filter_partitions and len(filter_partitions) > 0:
        logger.info("Filtering to %s partitions...", len(filter_partitions))
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(filter_partitions),
            vds.variant_data._filter_partitions(filter_partitions),
        )

    vmt = vds.variant_data

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

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    return vds
