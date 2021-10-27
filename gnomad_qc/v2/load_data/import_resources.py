from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources import *
import argparse


def import_clinvar(overwrite: bool = False):
    from datetime import datetime
    clinvar_ht = hl.import_vcf(clinvar_vcf_path, min_partitions=500, skip_invalid_loci=True).rows()
    clinvar_ht = clinvar_ht.annotate_globals(imported_on=datetime.now().strftime('%Y-%m-%d'))
    clinvar_ht = hl.vep(clinvar_ht, vep_config)
    clinvar_ht.write(clinvar_ht_path, overwrite)


def import_de_novos(overwrite: bool = False):
    denovo_ht = hl.import_table('gs://gnomad/resources/validated_de_novos.txt.bgz', types={'pos': hl.tint32})
    denovo_ht = denovo_ht.transmute(locus=hl.locus(hl.str(denovo_ht.chrom), denovo_ht.pos))
    denovo_ht = denovo_ht.annotate(alleles=hl.min_rep(denovo_ht.locus, [denovo_ht.ref, denovo_ht.alt])[1])
    denovo_ht = denovo_ht.drop('ref', 'alt').key_by('locus', 'alleles')
    denovo_ht.write('', overwrite)


def import_methylation(overwrite: bool = False):
    data_path = 'gs://gnomad-resources/methylation/source/all_methylation.txt.bgz'
    kt = hl.import_table(data_path, impute=True, min_partitions=100)
    kt = kt.transmute(CHROM=hl.cond(kt['#CHROM'] == 'M', 'MT', kt['#CHROM']))
    kt = kt.transmute(locus=hl.locus(kt.CHROM, kt.POS))
    kt.key_by('locus').write(methylation_sites_ht_path(), overwrite)

    ht = hl.read_table(methylation_sites_ht_path())
    ref_37 = hl.get_reference('GRCh37')
    ref_38 = hl.get_reference('GRCh38')
    ref_37.add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', ref_38)
    ht = ht.annotate(new_locus=hl.liftover(ht.locus, 'GRCh38', include_strand=True),
                     old_locus=ht.locus)
    ht = ht.key_by(locus=ht.new_locus.result)
    ht.write(methylation_sites_ht_path(ref='GRCh38'), overwrite=overwrite)


def import_exac_data(overwrite: bool = False):
    vcf_path = "gs://gnomad/raw/source/ExAC.r1.sites.vep.vcf.gz"
    vds = hl.import_vcf(vcf_path, force_bgz=True, min_partitions=5000).rows()
    vds = hl.split_multi_hts(vds)
    vds = hl.vep(vds, vep_config)
    vds.write(exac_release_sites_ht_path(), overwrite)


def import_cpgs(overwrite: bool = False):
    hl.import_vcf('gs://gnomad-public/resources/cpg.vcf.bgz', min_partitions=20).rows().write(cpg_sites_ht_path(), overwrite)


def import_truth_sets(overwrite: bool = False):
    root = 'gs://gnomad-public/truth-sets'
    truth_sets = [
        '1000G_omni2.5.b37.vcf.bgz',
        'hapmap_3.3.b37.vcf.bgz',
        'Mills_and_1000G_gold_standard.indels.b37.vcf.bgz',
        'NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.bgz',
        'hybrid.m37m.vcf.bgz',
        '1000G_phase1.snps.high_confidence.b37.vcf.bgz'
    ]
    for truth_vcf in truth_sets:
        vds_path = truth_vcf.replace('.vcf.bgz', '.mt')
        vds = hl.import_vcf('{}/source/{}'.format(root, truth_vcf), min_partitions=10)
        hl.split_multi_hts(vds).write('{}/hail-{}/{}'.format(root, CURRENT_HAIL_VERSION, vds_path), overwrite)


def main(args):
    hl.init(min_block_size=0)

    if args.import_truth_sets:
        import_truth_sets(args.overwrite)

    if args.import_cpgs:
        import_cpgs(args.overwrite)

    if args.import_methylation:
        import_methylation(args.overwrite)

    if args.import_exac_data:
        import_exac_data(args.overwrite)

    if args.import_clinvar:
        import_clinvar(args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--import_truth_sets', help='Import truth data', action='store_true')
    parser.add_argument('--import_cpgs', help='Import CpG data', action='store_true')
    parser.add_argument('--import_methylation', help='Import methylation data', action='store_true')
    parser.add_argument('--import_exac_data', help='Import ExAC data', action='store_true')
    parser.add_argument('--import_clinvar', help='Import clinvar data', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)

# # v1
# hc = HailContext(min_block_size=32)
# hc.import_vcf('gs://exac2/ExAC.r1.sites.vep.vcf.gz', sites_only=True, force_bgz=True).vep(config=vep_config).write(final_exac_sites_vds)

# exac = hc.import_vcf("gs://gnomad/raw/hail-0.1/vds/exac/exac_all.vcf.gz", header_file="gs://gnomad/raw/hail-0.1/vds/exac/header.vcf", force_bgz=True).min_rep().write("gs://gnomad/raw/hail-0.1/vds/exac/exac.vds")

# Raw counts:
# print hc.import_vcf('gs://exac2/variantqc/ExAC.merged.sites_only.vcf.ICfiltered.recalibrated.vcf.bgz').query_variants('variants.count()')[0]  # 18444471 sites
# print hc.read('gs://exac2/variantqc/exac2_vqsr.vds').query_variants('variants.count()')[0]  # 21352671 variants
