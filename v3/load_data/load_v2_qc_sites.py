from gnomad_hail import *
from v2.resources.sample_qc import qc_mt_path
import argparse


def gnomad_v2_qc_sites(overwrite: bool):
    v2_qc_sites = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True)).rows().select()
    v2_qc_sites = v2_qc_sites.key_by()
    v2_qc_sites = v2_qc_sites.transmute(
        locus=hl.liftover(v2_qc_sites.locus, 'GRCh38'),
        locus37=v2_qc_sites.locus
    )
    v2_qc_sites = v2_qc_sites.filter(v2_qc_sites.locus.contig == 'chr' + v2_qc_sites.locus37.contig).key_by('locus', 'alleles')
    v2_qc_sites.repartition(100, shuffle=False).write(gnomad_v2_qc_sites_path, overwrite=overwrite)


def main(args):

    hl.init(log='/load_resources.hail.log', default_reference='GRCh38')
    hl.get_reference('GRCh37').add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', hl.get_reference('GRCh38'))

    if args.gnomad_v2_qc_sites:
        gnomad_v2_qc_sites(args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--gnomad_v2_qc_sites', help='Lift-over gnomAD v2 High quality QC sites', action='store_true')
    parser.add_argument('--overwrite', help='Overwrites existing files', action='store_true')
    main(parser.parse_args())