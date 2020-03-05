from gnomad_qc.v2.resources.sample_qc import qc_mt_path
from gnomad_qc.v3.resources import gnomad_v2_qc_sites
import argparse
import hail as hl


def main(args):

    hl.init(log='/load_resources.hail.log', default_reference='GRCh38')
    hl.get_reference('GRCh37').add_liftover('gs://hail-common/references/grch37_to_grch38.over.chain.gz', hl.get_reference('GRCh38'))

    v2_qc_sites = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True)).rows().select()
    v2_qc_sites = v2_qc_sites.key_by()
    v2_qc_sites = v2_qc_sites.transmute(
        locus=hl.liftover(v2_qc_sites.locus, 'GRCh38'),
        locus37=v2_qc_sites.locus
    )
    v2_qc_sites = v2_qc_sites.filter(v2_qc_sites.locus.contig == 'chr' + v2_qc_sites.locus37.contig).key_by('locus', 'alleles')
    v2_qc_sites.repartition(100, shuffle=False).write(gnomad_v2_qc_sites.path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrites existing files', action='store_true')
    main(parser.parse_args())