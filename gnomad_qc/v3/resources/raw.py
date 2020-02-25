import hail as hl
from gnomad_hail.resources import MatrixTableResource
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples
from gnomad_qc.v3.resources.meta import meta


def get_gnomad_v3_mt(
        key_by_locus_and_alleles: bool = False,
        remove_hard_filtered_samples: bool = True,
        release_only: bool = False,
        samples_meta: bool = False
) -> hl.MatrixTable:
    mt = gnomad_v3_genotypes.mt()
    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(hl.ir.MatrixKeyRowsBy(mt._mir, ['locus', 'alleles'], is_sorted=True))

    if remove_hard_filtered_samples:
        mt = mt.filter_cols(hl.is_missing(hard_filtered_samples.ht()[mt.col_key]))

    if samples_meta:
        mt = mt.annotate_cols(meta=meta.ht()[mt.col_key])

        if release_only:
            mt = mt.filter_cols(mt.meta.release)

    elif release_only:
        mt = mt.filter_cols(meta.ht()[mt.col_key].release)

    return mt


# V3 genotype data
gnomad_v3_genotypes = MatrixTableResource("gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned.mt")
