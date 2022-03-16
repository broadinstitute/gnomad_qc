import hail as hl
​
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
​
hl.init(log='hail.log', default_reference='GRCh38')
vds = get_gnomad_v4_vds(remove_hard_filtered_samples=False, split=False)
vmt = hl.read_matrix_table('gs://gnomad-tmp/cdv/v4-fingerprint-checkpoint-filter-loci.mt')
# from get_gnomad_v4_vds(split=True)
vmt = vmt.annotate_rows(
    n_unsplit_alleles=hl.len(vmt.alleles),
    mixed_site=(hl.len(vmt.alleles) > 2)
        & hl.any(lambda a: hl.is_indel(vmt.alleles[0], a), vmt.alleles[1:])
        & hl.any(lambda a: hl.is_snp(vmt.alleles[0], a), vmt.alleles[1:]),
)
vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
vds = hl.vds.VariantDataset(vds.reference_data, vmt)
ht = hl.import_table(
    "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.haplotype_database.txt",
     comment=("@"),
)
variants = ht.aggregate(hl.agg.collect(hl.struct(locus=hl.locus(ht['#CHROMOSOME'], hl.int(ht['POSITION'])), alleles=hl.array([ht.MAJOR_ALLELE, ht.MINOR_ALLELE]))))
vmt = vds.variant_data.filter_rows(hl.set(variants).contains(vds.variant_data.row_key))
vds = hl.vds.VariantDataset(reference_data=vds.reference_data, variant_data=vmt)
mt = hl.vds.to_dense_mt(vds)
mt = mt.checkpoint("gs://gnomad/v4.0/sample_qc/exomes/gnomad.exomes.v4.0.fingerprinting_variants.mt", overwrite=True)
mt = mt.naive_coalesce(1)


mt = hl.read_matrix_table('gs://gnomad-tmp/cdv/gnomad.exomes.v4.0.fingerprinting_variants.mt')
mt = mt.select_entries("GT", "PL")
hl.export_vcf(mt, 'gs://gnomad/v4.0/sample_qc/exomes/gnomad.exomes.v4.0.fingerprinting_variants.vcf.bgz')