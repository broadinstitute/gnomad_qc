
# NOTE
# This script is kept here only for archiving purpose.
# It has no practical utility and is outdated

# import hail as hl
#
# # hl.init(sc=sc)
#
# # Note: all files were moved to HDFS on spark-c002 to speed along imports (now at 18 hours on 100 cores)
#
# meta_path = 'gnomad_meta/metadata_import_table_intermediate_oct_09_2017.txt'
# vcfs_path = 'gnomad_exomes/scattered/*/genotypes.unfiltered.vcf.gz'
# out_vds_path = 'gnomad.vds'
# vqsr_file = 'gnomad_meta/ExAC.merged.sites_only.vcf.ICfiltered.recalibrated.vcf.bgz'
#
# meta_kt = hl.import_table(meta_path, impute=True).key_by('sample')
# vqsr_vds = hl.import_vcf(vqsr_file)
#
# vds = hl.import_vcf(vcfs_path, call_fields=['GT', 'PGT'], force_bgz=True, header_file='gnomad_exomes/gnomad_exomes_header_fixed.vcf')
#
# vds = vds.filter_cols(meta_kt[vds.s].cloudable)
# vds = vds.annotate_rows(info=vqsr_vds[(vds.locus, vds.alleles), :].info)
# vds = hl.min_rep(vds, left_aligned=True)
# vds.write(out_vds_path)

# To copy to cloud:
# cd /humgen/atgu1/fs03/konradk/gnomad/2.1/hail
# hadoop distcp -libjars gcs-connector-latest-hadoop2.jar -conf core-site.xml -m 50 hdfs:///user/konradk/gnomad.vds gs://gnomad/raw/hail-0.2/vds/exomes/gnomad.exomes.vds