# noqa: D100
import logging
from typing import List, Optional

import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VariantDatasetResource,
    VersionedMatrixTableResource,
    VersionedVariantDatasetResource,
)

from gnomad_qc.v3.resources.constants import CURRENT_RELEASE, CURRENT_VERSION
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)


def get_gnomad_v3_vds(
    split: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
    test: bool = False,
    filter_partitions: Optional[List[int]] = None,
) -> hl.vds.VariantDataset:
    """
    Get gnomAD VariantDataset with desired filtering and metadata annotations.

    :param split: Perform split on VDS - Note: this will perform a split on the VDS
        rather than grab an already split VDS.
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard
        filters (only relevant after sample QC).
    :param release_only: Whether to filter the VDS to only samples available for
        release (can only be used if metadata is present).
    :param samples_meta: Whether to add metadata to VDS variant_data in 'meta' column.
    :param test: Whether to use the test VDS instead of the full v3 VDS.
    :param filter_partitions: Optional argument to filter the VDS to specific partitions in the provided list.
    :return: gnomAD v3 dataset with chosen annotations and filters.
    """
    if test:
        vds = gnomad_v3_testset_vds.vds()
    else:
        vds = gnomad_v3_genotypes_vds.vds()

    if filter_partitions:
        logger.info("Filtering to %s partitions...", len(filter_partitions))
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(filter_partitions),
            vds.variant_data._filter_partitions(filter_partitions),
        )

    if remove_hard_filtered_samples:
        vds = hl.vds.filter_samples(
            vds,
            hard_filtered_samples.ht(),
            keep=False,
        )

    if samples_meta or release_only:
        meta_ht = meta.ht()
        if release_only:
            vds = hl.vds.filter_samples(
                vds,
                meta_ht.filter(meta_ht.release),
            )

        if samples_meta:
            vd = vds.variant_data
            vds = hl.vds.VariantDataset(
                vds.reference_data, vd.annotate_cols(meta=meta_ht[vd.col_key])
            )

    if split:
        vd = vds.variant_data
        vd = vd.annotate_rows(
            n_unsplit_alleles=hl.len(vd.alleles),
            mixed_site=(hl.len(vd.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(vd.alleles[0], a), vd.alleles[1:])
            & hl.any(lambda a: hl.is_snp(vd.alleles[0], a), vd.alleles[1:]),
        )
        vd = hl.experimental.sparse_split_multi(vd, filter_changed_loci=True)
        vds = hl.vds.VariantDataset(
            vds.reference_data,
            vd,
        )

    return vds


def get_gnomad_v3_mt(
    split: bool = False,
    key_by_locus_and_alleles: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
    test: bool = False,
) -> hl.MatrixTable:
    """
    Get gnomAD data with desired filtering and metadata annotations.

    :param split: Perform split on MT - Note: this will perform a split on the MT rather than grab an already split MT
    :param key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles (only needed for v3)
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :param samples_meta: Whether to add metadata to MT in 'meta' column
    :param test: Whether to use the test MT instead of the full v3 MT
    :return: gnomAD v3 dataset with chosen annotations and filters
    """
    logger.warning(
        "The raw v3.1 sparse MatrixTable has been deleted. This function now calls "
        "'get_gnomad_v3_vds' and turns the VariantDataset into a sparse MatrixTable."
        " We recommend direct use of the VariantDataset by calling 'get_gnomad_v3_vds' "
        "instead."
    )
    vds = get_gnomad_v3_vds(
        split=split,
        remove_hard_filtered_samples=remove_hard_filtered_samples,
        release_only=release_only,
        samples_meta=samples_meta,
        test=test,
    )
    mt = hl.vds.to_merged_sparse_mt(vds)
    if not key_by_locus_and_alleles:
        mt = mt.key_rows_by(mt.locus)

    return mt


# NOTE: These have been removed, use the VariantDataset!
_gnomad_v3_genotypes = {
    "3": MatrixTableResource(
        "gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned.mt"
    ),
    "3.1": MatrixTableResource(
        "gs://gnomad/raw/genomes/3.1/gnomad_v3.1_sparse_unsplit.repartitioned.mt"
    ),
}


# The same raw VDS is used for v3.1.1, v3.1.2, and v3.1
gnomad_v3_genotypes_vds = VersionedVariantDatasetResource(
    CURRENT_VERSION,
    {"3.1": VariantDatasetResource("gs://gnomad/v3.1/raw/gnomad_v3.1.vds")},
)
# v3 test dataset VDS
gnomad_v3_testset_vds = VariantDatasetResource(
    "gs://gnomad/v3.1/raw/gnomad_v3.1.test.vds"
)


# NOTE: These have been removed, use the VariantDataset!
# The same raw MT is used for v3.1.1, v3.1.2, and v3.1
gnomad_v3_genotypes = VersionedMatrixTableResource(
    CURRENT_VERSION,
    _gnomad_v3_genotypes,
)

# v3 test dataset sparse MT
gnomad_v3_testset = MatrixTableResource(
    "gs://gnomad/raw/genomes/3.1/gnomad_v3.1_sparse.test.mt"
)


def qc_temp_prefix(version: str = CURRENT_RELEASE) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return
    :return: Path to bucket with temporary QC data
    """
    return f"gs://gnomad-tmp-4day/gnomad_v{version}_qc_data/"


def get_checkpoint_path(
    name: str, version: str = CURRENT_RELEASE, mt: bool = False
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable
    :param version: Version of annotation path to return
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    """
    return f'{qc_temp_prefix(version)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(name: str, version: str = CURRENT_RELEASE) -> str:
    """
    Create a path for Hail log files.

    :param str name: Name of log file
    :param version: Version of annotation path to return
    :return: Output log path
    """
    return f"{qc_temp_prefix(version)}{name}.log"


# Resources need for SV histogram generation
gnomad_sv_bucket_path = "gs://talkowski-sv-gnomad-v3-release/Releasable_freeze1_202303"
"""
Path to bucket with gnomAD SV results (VCFs per chromosome).

Bucket also contains list of releasable samples (see below).
"""

gnomad_sv_release_samples_list_path = (
    f"{gnomad_sv_bucket_path}/gnomAD.v3.SV.releasable.samples"
)
