import hail as hl
from gnomad.resources import MatrixTableResource
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples
from gnomad_qc.v3.resources.meta import meta

RELEASES = ["3", "3.1"]
CURRENT_RELEASE = "3.1"
CURRENT_META_VERSION = "2019-09-27"
META_VERSIONS = ["2019-09-27"]
CURRENT_PROJECT_META_VERSION = "2020-09-11"
PROJECT_META_VERSIONS = ["09-09-2019"]


def get_gnomad_v3_mt(
    version: str = CURRENT_RELEASE,
    key_by_locus_and_alleles: bool = False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
    samples_meta: bool = False,
    meta_version: str = CURRENT_METADATA_VERSION,
) -> hl.MatrixTable:
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
        mt = mt.filter_cols(meta.ht()[mt.col_key].release)

    return mt


# V3 genotype data
def gnomad_v3_genotypes(version) -> MatrixTableResource:
    if version not in RELEASES:
        return DataException("Select version as one of: {}".format(",".join(RELEASES)))

    if version == 3.0:
        return MatrixTableResource(
            "gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned.mt"
        )
    if version == 3.1:
        return MatrixTableResource(
            "gs://gnomad/raw/genomes/3.1/gnomad_v3.1_sparse_unsplit.repartitioned.mt"
        )
