import hail as hl
from gnomad.resources.resource_utils import (
    VariantDatasetResource,
    VersionedVariantDatasetResource,
)

from gnomad_qc.v4.resources.constants import CURRENT_VERSION
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.sample_qc import hard_filtered_samples

from ukbb_qc.resources.basics import excluded_samples_path
from ukbb_qc.resources.resource_utils import CURRENT_FREEZE as CURRENT_UKBB_FREEZE


def get_gnomad_v4_vds(
    split=False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
) -> hl.vds.VariantDataset:
    """
    Wrapper function to get gnomAD v4 data with desired filtering and metadata annotations.

    :param split: Perform split on MT - Note: this will perform a split on the MT rather than grab an already split MT
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :return: gnomAD v4 dataset with chosen annotations and filters
    """
    vds = gnomad_v4_genotypes.vds()
    if remove_hard_filtered_samples:
        vds = hl.vds.filter_samples(
            vds, hard_filtered_samples.versions[CURRENT_VERSION].ht(), keep=False
        )

    if release_only:
        meta_ht = meta.versions[CURRENT_VERSION].ht()
        meta_ht = meta_ht.filter(meta_ht.release)
        vds = hl.vds.filter_samples(vds, meta_ht)

    if split:
        vmt = vds.variant_data
        vmt = vmt.annotate_rows(
            n_unsplit_alleles=hl.len(vmt.alleles),
            mixed_site=(hl.len(vmt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(vmt.alleles[0], a), vmt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(vmt.alleles[0], a), vmt.alleles[1:]),
        )
        vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
        vds = VariantDataset(vds.reference_data, vmt)

    # Remove withdrawn UKBB samples
    excluded_ukbb_samples_ht = hl.import_table(
        excluded_samples_path(CURRENT_UKBB_FREEZE)
    )
    vds = hl.vds.filter_samples(vds, excluded_ukbb_samples_ht)

    return vds


_gnomad_v4_genotypes = {
    "4": VariantDatasetResource("gs://gnomad/raw/exomes/4.0/gnomad_v4.0.vds"),
}


gnomad_v4_genotypes = VersionedVariantDatasetResource(
    CURRENT_VERSION, _gnomad_v4_genotypes,
)


def qc_temp_prefix(version: str = CURRENT_VERSION) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return
    :return: Path to bucket with temporary QC data
    """
    return f"gs://gnomad-tmp/gnomad_v{version}_qc_data/"


def get_checkpoint_path(
    name: str, version: str = CURRENT_VERSION, mt: bool = False
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable
    :param version: Version of annotation path to return
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    """
    return f'{qc_temp_prefix(version)}{name}.{"mt" if mt else "ht"}'


def testset_vds(version: str = CURRENT_VERSION) -> str:
    """
    Return path to the testset VDS.

    :param version: Version of annotation path to return
    :return: Path to VDS with testset
    """
    vds = hl.vds.read_vds(
        "gs://gnomad/raw/exomes/{version}/testing/gnomad_v{version}_test.vds"
    )
    # Will there only be one version of this?

    return vds
