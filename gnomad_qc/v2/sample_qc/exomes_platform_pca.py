from gnomad.utils.filtering import filter_to_autosomes
from gnomad.utils.slack import slack_notifications
import numpy as np
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import get_gnomad_data, evaluation_intervals_path
import argparse
import hdbscan
import logging
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_b")
logger.setLevel(logging.INFO)


def assign_platform_pcs(platform_pc_table: hl.Table) -> hl.Table:
    """
    Function assumes that platform_pc_table contains columns named 'combined_sample', 'gross_platform' (known labels), 'callratePC<n>'

    :param Table platform_pc_table: Table containing samples and callrate PCs
    :return: Table containing samples, callrate PCs, and imputed platform labels
    :rtype: Table
    """
    # Read and format data for clustering
    data = platform_pc_table.to_pandas()
    callrate_data = np.matrix(data['scores'].tolist())
    logger.info('Assigning platforms to {} exome samples in MT...'.format(len(callrate_data)))

    # Cluster data
    clusterer = hdbscan.HDBSCAN(min_cluster_size=100)
    cluster_labels = clusterer.fit_predict(callrate_data)
    n_clusters = len(set(cluster_labels)) - (-1 in cluster_labels)  # NOTE: -1 is the label for noisy (un-classifiable) data points
    logger.info('Found {} unique platforms during platform imputation...'.format(n_clusters))

    data['qc_platform'] = cluster_labels
    ht = hl.Table.from_pandas(data, key=['data_type', 's'])
    ht = ht.annotate(qc_platform=hl.str(ht.qc_platform))
    return ht


def main(args):
    hl.init(log='/platform_pca.log')

    if not args.skip_prepare_data_for_platform_pca:
        # ~1 hour on 800 cores (3/8/18)
        logger.info('Preparing data for platform PCA...')
        mt = get_gnomad_data('exomes', adj=True, raw=False, meta_root=None, fam_root=None, split=False)
        mt = filter_to_autosomes(mt)
        intervals = hl.import_locus_intervals(evaluation_intervals_path)
        mt = mt.annotate_rows(interval=intervals[mt.locus].target)
        mt = mt.filter_rows(hl.is_defined(mt.interval) & (hl.len(mt.alleles) == 2))
        mt = mt.select_entries(GT=hl.or_missing(hl.is_defined(mt.GT), hl.struct()))
        callrate_mt = mt.group_rows_by(mt.interval).aggregate(callrate=hl.agg.fraction(hl.is_defined(mt.GT)))
        callrate_mt.write(exome_callrate_mt_path, args.overwrite)

    if not args.skip_run_platform_pca:
        logger.info('Running platform PCA...')
        qc_ht = hl.read_table(qc_ht_path('exomes', 'hard_filters')).key_by('s')
        callrate_mt = hl.read_matrix_table(exome_callrate_mt_path)
        callrate_mt = callrate_mt.filter_cols(hl.len(qc_ht[callrate_mt.col_key].hard_filters) == 0)
        callrate_mt = callrate_mt.annotate_entries(callrate=hl.int(callrate_mt.callrate > 0.25))
        # Center until Hail's PCA does it for you
        callrate_mt = callrate_mt.annotate_rows(mean_callrate=hl.agg.mean(callrate_mt.callrate))
        callrate_mt = callrate_mt.annotate_entries(callrate=callrate_mt.callrate - callrate_mt.mean_callrate)
        eigenvalues, scores, _ = hl.pca(callrate_mt.callrate, compute_loadings=False)
        logger.info('Eigenvalues: {}'.format(eigenvalues))
        # [731282566.2824697, 78687228.90071851, 43837650.51729764, 33969298.61827205, 26308703.539534636, 21102437.512725923, 16949828.555817757, 12994894.187041137, 8372332.274295175, 8128326.814388647]
        scores.write(exome_callrate_scores_ht_path)

    logger.info('Annotating with platform PCs and known platform annotations...')
    scores = hl.read_table(exome_callrate_scores_ht_path).annotate(data_type='exomes')
    if args.pc_scores_in_separate_fields:
        scores = scores.transmute(scores=[
            scores[ann] for ann in sorted(
                [ann for ann in scores.row if ann.startswith("PC")],
                key=lambda x: int(x[2:])
            )
        ])
    platform_pcs = assign_platform_pcs(scores)
    platform_pcs.write(qc_ht_path('exomes', 'platforms'), overwrite=args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--skip_prepare_data_for_platform_pca', help='Skip prepping data for platform imputation (assuming already done)', action='store_true')
    parser.add_argument('--skip_run_platform_pca', help='Skip platform PCA (assuming already done)', action='store_true')
    parser.add_argument('--pc_scores_in_separate_fields', help='This option was added to deal with legacy scores HT, where the PC scores where stored in multiple annotations (PC1, ... PCn)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
