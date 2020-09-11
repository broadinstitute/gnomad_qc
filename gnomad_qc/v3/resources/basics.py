import hail as hl
from gnomad.resources import DataException, MatrixTableResource
from gnomad_qc.v3.resources import (
    CURRENT_META_VERSION,
    CURRENT_RELEASE,
    hard_filtered_samples,
    meta,
    RELEASES,
)


def get_gnomad_v3_mt(
    version: str = CURRENT_RELEASE,
    key_by_locus_and_alleles: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
    meta_version: str = CURRENT_META_VERSION,
) -> hl.MatrixTable:
    """
    Wrapper function to get gnomAD data with desired filtering and metadata annotations

    :param version: Version of MT to return
    :param key_by_locus_and_alleles: Whether to key the MatrixTable by locus and alleles (only needed for v3)
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :param samples_meta: Whether to add metadata to MT in 'meta' column
    :param meta_version: Version of metadata (Default is CURRENT_META_VERSION)
    :return: gnomAD v3 dataset with chosen annotations and filters
    """
    mt = gnomad_v3_genotypes(version).mt()
    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(
            hl.ir.MatrixKeyRowsBy(mt._mir, ["locus", "alleles"], is_sorted=True)
        )

    if remove_hard_filtered_samples:
        mt = mt.filter_cols(hl.is_missing(hard_filtered_samples().ht()[mt.col_key]))

    if samples_meta:
        mt = mt.annotate_cols(meta=meta(version, meta_version).ht()[mt.col_key])

        if release_only:
            mt = mt.filter_cols(mt.meta.release)

    elif release_only:
        mt = mt.filter_cols(meta(version, meta_version).ht()[mt.col_key].release)

    return mt


# V3 genotype data
def gnomad_v3_genotypes(version: str = CURRENT_RELEASE) -> MatrixTableResource:
    """
    Get gnomAD v3 raw MatrixTable

    :param version: Version of raw MT to return
    :return: gnomAD raw v3 MatrixTable
    """
    if version not in RELEASES:
        return DataException("Select version as one of: {}".format(",".join(RELEASES)))

    if version == "3":
        return MatrixTableResource(
            "gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned.mt"
        )
    if version == "3.1":
        return MatrixTableResource(
            "gs://gnomad/raw/genomes/3.1/gnomad_v3.1_sparse_unsplit.repartitioned.mt"
        )
