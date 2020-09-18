import argparse
import logging

import hail as hl

from gnomad_qc.v2.resources.sample_qc import get_liftover_v2_qc_mt
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import (qc, v2_v3_pc_relate_pca_scores,
                                              v2_v3_relatedness)

logger = logging.getLogger("v2_pc_relate")


def main(args):
    if args.join_qc_mt:
        v2_qc_mt_liftover = get_liftover_v2_qc_mt('exomes', ld_pruned=True, release_only=True)
        v2_qc_mt_liftover = v2_qc_mt_liftover.key_cols_by(s=v2_qc_mt_liftover.s, data_type="v2_exomes")
        v3_qc_mt = qc.mt()
        v3_qc_mt = v3_qc_mt.filter_cols(meta.ht()[v3_qc_mt.col_key].release)
        v3_qc_mt = v3_qc_mt.select_rows().select_cols()
        v3_qc_mt = v3_qc_mt.key_cols_by(s=v3_qc_mt.s, data_type="v3_genomes")
        joint_qc_mt = v2_qc_mt_liftover.union_cols(v3_qc_mt)
        joint_qc_mt.write("gs://gnomad-tmp/v2_exomes_v3_joint_qc.mt", overwrite=args.overwrite)

    if args.run_pc_relate:
        logger.info('Running PC-Relate')
        logger.warning("PC-relate requires SSDs and doesn't work with preemptible workers!")
        joint_qc_mt = hl.read_matrix_table("gs://gnomad-tmp/v2_exomes_v3_joint_qc.mt")
        joint_qc_mt = joint_qc_mt.sample_rows(0.1)
        eig, scores, _ = hl.hwe_normalized_pca(joint_qc_mt.GT, k=10, compute_loadings=False)
        scores = scores.checkpoint(v2_v3_pc_relate_pca_scores.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        relatedness_ht = hl.pc_relate(joint_qc_mt.GT, min_individual_maf=0.01, scores_expr=scores[joint_qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=0.1, statistics='all')
        relatedness_ht.write(v2_v3_relatedness.path, args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--join_qc_mt', help="Join liftover v2 exomes release and v3 release QC MT", action='store_true')
    parser.add_argument('--run_pc_relate', help='Run PC-relate on v2 exomes/v3 combined release QC MT', action='store_true')

    main(parser.parse_args())
