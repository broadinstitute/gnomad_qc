import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VersionedMatrixTableResource,
)

from gnomad_qc.v3.resources.constants import (
    CURRENT_RELEASE,
    CURRENT_VERSION,
)
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples


def get_gnomad_v3_mt(
    split=False,
    key_by_locus_and_alleles: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
) -> hl.MatrixTable:
    """
    Wrapper function to get gnomAD data with desired filtering and metadata annotations

    :param split: Perform split on MT - Note: this will perform a split on the MT rather than grab an already split MT
    :param key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles (only needed for v3)
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :param samples_meta: Whether to add metadata to MT in 'meta' column
    :return: gnomAD v3 dataset with chosen annotations and filters
    """
    mt = gnomad_v3_genotypes.mt()
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


_gnomad_v3_genotypes = {
    "3": MatrixTableResource(
        "gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned.mt"
    ),
    "3.1": MatrixTableResource(
        "gs://gnomad/raw/genomes/3.1/gnomad_v3.1_sparse_unsplit.repartitioned.mt"
    ),
}

# The same raw MT is used for v3.1.1, v3.1.2, and v3.1
gnomad_v3_genotypes = VersionedMatrixTableResource(
    CURRENT_VERSION, _gnomad_v3_genotypes,
)


def qc_temp_prefix(version: str = CURRENT_RELEASE) -> str:
    """
    Returns path to temporary QC bucket.

    :param version: Version of annotation path to return
    :return: Path to bucket with temporary QC data
    """
    return f"gs://gnomad-tmp/gnomad_v{version}_qc_data/"


def get_checkpoint_path(
    name: str, version: str = CURRENT_RELEASE, mt: bool = False
) -> str:
    """
    Creates a checkpoint path for Table or MatrixTable

    :param str name: Name of intermediate Table/MatrixTable
    :param version: Version of annotation path to return
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    :rtype: str
    """
    return f'{qc_temp_prefix(version)}{name}.{"mt" if mt else "ht"}'
