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
    n_partitions: Optional[int] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    filter_variant_ht: Optional[hl.Table] = None,
    entries_to_keep: Optional[List[str]] = None,
    checkpoint_variant_data: bool = False,
) -> hl.vds.VariantDataset:
    """
    Load the AOU VDS.

    :param split: Whether to split the VDS into separate datasets for each chromosome. Default is False.
    :param chrom: Chromosome(s) to filter the VDS to. Can be a single chromosome or a list of chromosomes.
    :param autosomes_only: Whether to include only autosomes.
    :param sex_chr_only: whether to include only sex chromosomes.
    :param n_partitions: Number of partitions to use for the VDS.
    :param filter_intervals: List of intervals to filter the VDS.
    :param split_reference_blocks: Whether to split the reference blocks. Default is True.
    :param filter_variant_ht: Table to filter the variant data.
    :param entries_to_keep: List of entries to keep in the variant data.
    :param checkpoint_variant_data: Whether to checkpoint the variant data. Default is False.
    :return: The AOU VDS.
    """
    aou_v8_resource = aou_genotypes

    if isinstance(chrom, str):
        chrom = [chrom]

    if autosomes_only or sex_chr_only:
        if sex_chr_only:
            vds = hl.vds.filter_chromosomes(aou_v8_resource.vds, keep=["chrX", "chrY"])
        else:
            vds = hl.vds.filter_chromosomes(aou_v8_resource.vds, keep_autosomes=True)
    elif autosomes_only and sex_chr_only:
        raise ValueError(
            "Only one of 'autosomes_only' or 'sex_chr_only' can be set to True."
        )

    if n_partitions and chrom:
        logger.info(
            "Filtering to chromosome(s) %s with %s partitions...", chrom, n_partitions
        )
        reference_data = hl.read_matrix_table(
            hl.vds.VariantDataset._reference_path(aou_v8_resource.path)
        )
        reference_data = hl.filter_intervals(
            reference_data,
            [hl.parse_locus_interval(x, reference_genome="GRCh38") for x in chrom],
        )
        intervals = reference_data._calculate_new_partitions(n_partitions)
        reference_data = hl.read_matrix_table(
            hl.vds.VariantDataset._reference_path(aou_v8_resource.path),
            _intervals=intervals,
        )
        variant_data = hl.read_matrix_table(
            hl.vds.VariantDataset._variants_path(aou_v8_resource.path),
            _intervals=intervals,
        )
        vds = hl.vds.VariantDataset(reference_data, variant_data)
    elif n_partitions:
        vds = hl.vds.read_vds(aou_v8_resource.path, n_partitions=n_partitions)
    else:
        vds = aou_v8_resource.vds()
        if chrom:
            logger.info("Filtering to chromosome %s...", chrom)
            vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    if filter_intervals:
        logger.info("Filtering to %s intervals...", len(filter_intervals))
        if isinstance(filter_intervals[0], str):
            filter_intervals = [
                hl.parse_locus_interval(x, reference_genome="GRCh38")
                for x in filter_intervals
            ]
        vds = hl.vds.filter_intervals(
            vds, filter_intervals, split_reference_blocks=split_reference_blocks
        )

    vmt = vds.variant_data

    if split:
        vmt = _split_and_filter_variant_data_for_loading(
            vmt, filter_variant_ht, entries_to_keep, checkpoint_variant_data
        )

    if entries_to_keep is not None:
        vmt = vmt.select_entries(*entries_to_keep)

    if checkpoint_variant_data:
        vmt = vmt.checkpoint(new_temp_file("vds_loading.variant_data", "mt"))

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    return vds
