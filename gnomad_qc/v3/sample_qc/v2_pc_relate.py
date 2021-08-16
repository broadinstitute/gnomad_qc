import argparse
import logging

import hail as hl

from gnomad_qc.v2.resources.basics import get_gnomad_meta
from gnomad_qc.v2.resources.sample_qc import get_liftover_v2_qc_mt
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import (
    qc,
    v2_v3_pc_relate_pca_scores,
    v2_v3_1_new_pc_relate_pca_scores,
    v2_v3_relatedness,
    v2_v3_1_new_relatedness,
)

logger = logging.getLogger("v2_pc_relate")


def main(args):
    if args.join_qc_mt:
        v2_qc_mt_liftover = get_liftover_v2_qc_mt('joint', ld_pruned=True)
        meta_ht = get_gnomad_meta('joint')
        v2_qc_mt_liftover = v2_qc_mt_liftover.key_cols_by("s", "data_type")
        v2_qc_mt_liftover = v2_qc_mt_liftover.filter_cols(meta_ht[v2_qc_mt_liftover.col_key].release)

        # Filter to only exomes
        v2_qc_mt_liftover = v2_qc_mt_liftover.filter_cols(v2_qc_mt_liftover.data_type == "exomes")
        v2_qc_mt_liftover = v2_qc_mt_liftover.key_cols_by(s=v2_qc_mt_liftover.s, data_type="v2_exomes").select_cols()
        v3_qc_mt = qc.mt()
        if args.filter_to_new_3_1_samples:
            v3_qc_mt = v3_qc_mt.filter_cols(
                meta.ht()[v3_qc_mt.col_key].release
                & hl.is_missing(meta.versions["3"].ht()[v3_qc_mt.col_key])
            )
        else:
            v3_qc_mt = v3_qc_mt.filter_cols(meta.ht()[v3_qc_mt.col_key].release)
        v3_qc_mt = v3_qc_mt.select_rows().select_cols()
        v3_qc_mt = v3_qc_mt.key_cols_by(s=v3_qc_mt.s, data_type="v3_genomes").select_cols()
        joint_qc_mt = v2_qc_mt_liftover.union_cols(v3_qc_mt)
        if args.filter_to_new_3_1_samples:
            joint_qc_mt.write("gs://gnomad-tmp/v2_exomes_v3.1_new_samples_joint_qc.mt", overwrite=args.overwrite)
        else:
            joint_qc_mt.write("gs://gnomad-tmp/v2_exomes_v3_joint_qc.mt", overwrite=args.overwrite)

    if args.run_pc_relate:
        scores_path = v2_v3_1_new_pc_relate_pca_scores.path if args.filter_to_new_3_1_samples else v2_v3_pc_relate_pca_scores.path
        relatedness_path = v2_v3_1_new_relatedness.path if args.filter_to_new_3_1_samples else v2_v3_relatedness.path

        logger.info('Running PC-Relate')
        logger.warning("PC-relate requires SSDs and doesn't work with preemptible workers!")
        if args.filter_to_new_3_1_samples:
            joint_qc_mt = hl.read_matrix_table("gs://gnomad-tmp/v2_exomes_v3.1_new_samples_joint_qc.mt")
        else:
            joint_qc_mt = hl.read_matrix_table("gs://gnomad-tmp/v2_exomes_v3_joint_qc.mt")
        #joint_qc_mt = joint_qc_mt.sample_rows(0.1)
        eig, scores, _ = hl.hwe_normalized_pca(joint_qc_mt.GT, k=10, compute_loadings=False)
        scores = scores.checkpoint(scores_path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        relatedness_ht = hl.pc_relate(joint_qc_mt.GT, min_individual_maf=0.01, scores_expr=scores[joint_qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=0.1, statistics='all')
        relatedness_ht.write(relatedness_path, args.overwrite)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--join_qc_mt', help="Join liftover v2 exomes release and v3 release QC MT", action='store_true')
    parser.add_argument('--run_pc_relate', help='Run PC-relate on v2 exomes/v3 combined release QC MT', action='store_true')
    parser.add_argument('--filter_to_new_3_1_samples', help='Only run the join/PC_relate on v2 exomes combined with v3.1 new samples.', action='store_true')


    main(parser.parse_args())
