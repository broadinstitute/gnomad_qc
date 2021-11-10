import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VersionedMatrixTableResource,
)

from gnomad_qc.v4.resources.constants import (
    CURRENT_RELEASE,
    CURRENT_VERSION,
)
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.sample_qc import hard_filtered_samples


def get_gnomad_v4_mt(
    split=False,
    key_by_locus_and_alleles: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
) -> hl.MatrixTable:
    """
    Wrapper function to get gnomAD data with desired filtering and metadata annotations.

    :param split: Perform split on MT - Note: this will perform a split on the MT rather than grab an already split MT
    :param key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :param samples_meta: Whether to add metadata to MT in 'meta' column
    :return: gnomAD v4 dataset with chosen annotations and filters
    """
    mt = gnomad_v4_genotypes.mt()
    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(
            hl.ir.MatrixKeyRowsBy(
                mt._mir, ["locus", "alleles"], is_sorted=True
            )  # Prevents hail from running sort on genotype MT which is already sorted by a unique locus
        )

    if remove_hard_filtered_samples:
        mt = mt.filter_cols(
            hl.is_missing(
                hard_filtered_samples.versions[CURRENT_VERSION].ht()[mt.col_key]
            )
        )

    if samples_meta:
        mt = mt.annotate_cols(meta=meta.versions[CURRENT_VERSION].ht()[mt.col_key])

        if release_only:
            mt = mt.filter_cols(mt.meta.release)

    elif release_only:
        mt = mt.filter_cols(meta.versions[CURRENT_VERSION].ht()[mt.col_key].release)

    if split:
        mt = mt.annotate_rows(
            n_unsplit_alleles=hl.len(mt.alleles),
            mixed_site=(hl.len(mt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(mt.alleles[0], a), mt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(mt.alleles[0], a), mt.alleles[1:]),
        )
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    return mt


# TODO: Will sparse_split_multi still work in v4?

_gnomad_v4_genotypes = {
    "4": MatrixTableResource(
        "gs://gnomad/v4/raw/exomes/gnomad_v4_sparse_unsplit.repartitioned.mt"
    ),
}


gnomad_v4_genotypes = VersionedMatrixTableResource(
    CURRENT_VERSION, _gnomad_v4_genotypes,
)


def qc_temp_prefix(version: str = CURRENT_RELEASE) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return
    :return: Path to bucket with temporary QC data
    """
    return f"gs://gnomad-tmp/gnomad_v{version}_qc_data/"


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
