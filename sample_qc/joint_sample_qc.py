from gnomad_hail import *
from gnomad_qc.resources.sample_qc import *
from create_fam import *
from sklearn.ensemble import RandomForestClassifier
import pandas as pd
import random
import pickle

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_c")
logger.setLevel(logging.INFO)


def read_and_pre_process_data(mt_path: str, ht_path: str) -> hl.MatrixTable:
    """
    :param str mt_path: Path to MT to be formatted for joining
    :param str ht_path: Path to HT used to annotate MT
    :return: MatrixTable with uniquified (prefixed) sample IDs that retains previously annotated permissions and hard filters
    :rtype: MatrixTable
    """
    ht = hl.read_table(ht_path).select('data_type', 's', 'hard_filters', 'perm_filters').key_by('s')
    mt = hl.read_matrix_table(mt_path)
    mt = mt.annotate_cols(**ht[mt.s]).key_cols_by('data_type', 's')
    mt = filter_to_autosomes(mt.filter_cols(hl.len(mt.hard_filters) == 0))
    return mt.select_entries('GT')


def make_rank_file(outfile: str) -> hl.Table:
    """
    Assign a rank describing retention preference in de-duplication and familial pruning to each exome and genome sample
    # NOTE: use order_by() method on Tables only for big datasets -- order_by() is distributed across nodes; whereas pandas sort is all local
    # NOTE: missing annotations (e.g., in `pct_bases_20x`) are prioritized last

    :param str outfile: filepath to tsv containing ranks assigned across genome and exome samples jointly
    :return: Table organized by uid (e.g., 'exome_ABC1234' or 'genome_DEF5678') and global rank
    :rtype: Table
    """

    # Load data
    exome = hl.import_table(rank_annotations_path('exomes'), impute=True).to_pandas()
    genome = hl.import_table(rank_annotations_path('genomes'), impute=True).to_pandas()
    exome_trio = hl.import_table(dup_pedigree_tsv_path('exomes'), impute=True, no_header=True).to_pandas()
    genome_trio = hl.import_table(dup_pedigree_tsv_path('genomes'), impute=True, no_header=True).to_pandas()

    # Select complete trios
    exome_trio = exome_trio[(exome_trio.pat_id != '0') & (exome_trio.mat_id != '0')]
    genome_trio = genome_trio[(genome_trio.pat_id != '0') & (genome_trio.mat_id != '0')]

    # Sort genomes
    genome['data_type'] = 'genomes'
    # 1 is "parent", 2 is "NA", 3 is "child"
    genome['parent_child'] = [1 if x else 2 for x in pd.Series(list(genome['s'])).isin(list(genome_trio['pat_id']) + list(genome_trio['mat_id']))]
    genome.loc[pd.Series(list(genome['s'])).isin(list(genome_trio['s'])), 'parent_child'] = 3
    sorted_genome = genome.sort_values(['pcr_free', 'parent_child', 'mean_dp'], ascending=[False, True, False])

    # Harmonize column names with exomes
    sorted_genome['internal'] = 'NaN'
    sorted_genome['pct_bases_20x'] = 'NaN'
    sorted_genome['project_id'] = 'NaN'

    # Sort exomes by internal vs external status
    exome['data_type'] = 'exomes'
    exome['parent_child'] = [1 if x else 2 for x in pd.Series(list(exome['s'])).isin(list(exome_trio['pat_id']) + list(exome_trio['mat_id']))]
    exome.loc[pd.Series(list(exome['s'])).isin(list(exome_trio['s'])), 'parent_child'] = 3

    exome_internal = exome.loc[exome.internal, ].copy()
    exome_internal['project_id'].replace(to_replace="^RP-", value="", regex=True, inplace=True)
    exome_internal['project_id'].replace(to_replace="^C", value="1000", regex=True, inplace=True)  # NOTE: C-projects in exomes are more recent/desirable than RP-projects, so boosting these
    exome_internal['project_id'] = pd.to_numeric(exome_internal['project_id'])
    sorted_exome_internal = exome_internal.sort_values(['project_id', 'parent_child', 'pct_bases_20x'], ascending=[False, True, False])

    exome_external = exome.loc[~exome.internal, ].copy()
    sorted_exome_external = exome_external.sort_values(['parent_child', 'pct_bases_20x'], ascending=[True, False])
    sorted_exome = pd.concat([sorted_exome_internal, sorted_exome_external])

    # Harmonize column names with genomes
    sorted_exome['pcr_free'] = 'NaN'
    sorted_exome['mean_dp'] = 'NaN'
    sorted_exome = sorted_exome[sorted_genome.columns]

    # Combine and rearrange by permissions
    sorted_data = pd.concat([sorted_genome, sorted_exome])

    releasable = sorted_data.loc[sorted_data.releasable, ].copy()
    nonreleasable = sorted_data.loc[~sorted_data.releasable, ].copy()
    sorted_data = pd.concat([releasable, nonreleasable])
    sorted_data.rename(index=str, columns={'project_id': 'rank_project_id'})
    sorted_data['rank'] = range(1, len(sorted_data) + 1)

    with hl.hadoop_open(outfile, 'w') as out:
        sorted_data.to_csv(out, sep="\t", index=False)
    new_data = hl.import_table(outfile, impute=True).key_by('data_type', 's')
    return new_data


def run_assign_population_pcs(pop_pc_table: hl.Table, outfile: str, picklefile: str, pcs: List[int],
                              fit: RandomForestClassifier = None, seed: int = 42) -> Tuple[hl.Table, RandomForestClassifier]:
    """
    :param Table pop_pc_table: Table containing population PCs ('PC<n>') as well as a column 'known_pop' with population labels
    :param str outfile: filepath to tsv with input samples and imputed population labels
    :param str picklefile: filepath to which the pickled random forest model is written
    :param list of int pcs: 1-based list of PCs to train the model on
    :param RandomForestClassifier fit: fit from a previously trained random forest model (i.e., the output from a previous RandomForestClassifier() call)
    :param int seed: Random seed
    :return: Table containing sample IDs and imputed population labels, trained random forest model
    :rtype: Table, RandomForestClassifier
    """
    data = pop_pc_table.to_pandas()
    data = expand_pd_array_col(data, 'scores', max(pcs), 'PC')
    new_data, pop_clf = assign_population_pcs(
        data,
        pc_cols=['PC{}'.format(pc) for pc in pcs],
        fit=fit,
        seed=seed
    )

    if not fit:
        # Pickle RF
        with hl.hadoop_open(picklefile, 'wb') as out:
            pickle.dump(pop_clf, out)

    with hl.hadoop_open(outfile, 'w') as out:
        new_data.to_csv(out, sep="\t", na_rep="NA", index=False)
    return hl.import_table(outfile, impute=True).key_by('data_type', 's'), pop_clf


def gnomad_sample_qc(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Filter MTs to bi-allelic sites and remove problematic intervals, and performs sample QC
    # TODO: consider reinstating inbreeding coefficient filter

    :param MatrixTable mt: MT on which sample QC metrics need to be computed
    :return: MT filtered to autosomes and high-confidence regions, with computed sample QC column annotations
    :rtype: MatrixTable
    """
    mt = filter_to_autosomes(filter_low_conf_regions(mt))
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)  # NOTE: this does not work on a split VDS!
    mt = hl.sample_qc(mt)
    mt = mt.annotate_rows(variant_qc=hl.struct(af=hl.agg.mean(mt.GT.n_alt_alleles()) / 2))
    mt = mt.annotate_cols(sample_qc=mt.sample_qc.annotate(f_inbreeding=hl.agg.inbreeding(mt.GT, mt.variant_qc.af)))
    return mt


def make_pop_filters_expr(mt: hl.MatrixTable, qc_metrics: List[str]) -> hl.expr.SetExpression:
    """
    :param MatrixTable mt: MT to which the pop-/platform-specific filtering should apply
    :param list qc_metrics: list of sample QC metrics for which to collect pop-/platform-specific filtering status
    :return: SetExpression value used to create the pop_platform_filters sample annotation
    :rtype: SetExpression
    """
    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(mt['fail_{}'.format(metric)], '{}'.format(metric)) for metric in qc_metrics]))


def get_related_samples_to_drop(rank_table: hl.Table, relatedness_ht: hl.Table) -> hl.Table:
    """
    Use the maximal independence function in Hail to intelligently prune clusters of related individuals, removing
    less desirable samples while maximizing the number of unrelated individuals kept in the sample set

    :param Table rank_table: Table with ranking annotations across exomes and genomes, computed via make_rank_file()
    :param Table relatedness_ht: Table with kinship coefficient annotations computed via pc_relate()
    :return: Table containing sample IDs ('s') to be pruned from the combined exome and genome sample set
    :rtype: Table
    """
    # Define maximal independent set, using rank list
    related_pairs = relatedness_ht.filter(relatedness_ht.kin > 0.08838835).select('i', 'j')
    related_samples = related_pairs.aggregate(hl.agg.collect_as_set(hl.agg.explode([related_pairs.i, related_pairs.j])))
    logger.info('{} samples with at least 2nd-degree relatedness found in callset'.format(len(related_samples)))
    max_rank = rank_table.count()
    related_pairs = related_pairs.annotate(id1_rank=hl.struct(id=related_pairs.i, rank=rank_table[related_pairs.i].rank),
                                           id2_rank=hl.struct(id=related_pairs.j, rank=rank_table[related_pairs.j].rank)
                                           ).select('id1_rank', 'id2_rank')

    def tie_breaker(l, r):
        return hl.or_else(l.rank, max_rank + 1) - hl.or_else(r.rank, max_rank + 1)

    related_samples_to_drop_ranked = hl.maximal_independent_set(related_pairs.id1_rank, related_pairs.id2_rank,
                                                                keep=False, tie_breaker=tie_breaker)
    return related_samples_to_drop_ranked.select(**related_samples_to_drop_ranked.node.id).key_by('data_type', 's')


def calculate_metrics(mt: hl.MatrixTable, qc_metrics: List[str], data_type: str) -> hl.Table:
    """
    Compute median, MAD, and upper and lower thresholds for each metric used in pop- and platform-specific outlier filtering

    :param MatrixTable mt: MT containing relevant sample QC metric annotations
    :param list qc_metrics: list of metrics for which to compute the critical values for filtering outliers
    :param str data_type: 'exomes' or 'genomes'
    :return: Table grouped by pop and platform, with upper and lower threshold values computed for each sample QC metric
    :rtype: Table
    """
    cols = ['qc_pop', 'qc_platform']
    key = 's'
    colnames = cols + ['sample_qc.{}'.format(metric) for metric in qc_metrics]
    new_colnames = cols + ['{}'.format(metric) for metric in qc_metrics]
    ht = mt.cols().flatten().select(*colnames).rename(dict(zip(colnames, new_colnames))).key_by(key)

    key_expr = ['{0}_median'.format(metric) for metric in qc_metrics] + ['{0}_mad'.format(metric) for metric in qc_metrics]
    value_expr = [hl.median(hl.agg.collect(ht['{}'.format(metric)])) for metric in qc_metrics] + \
                 [1.4826 * hl.median(hl.abs(hl.agg.collect(ht['{0}'.format(metric)]) -
                                            hl.median(hl.agg.collect(ht['{}'.format(metric)])))) for metric in qc_metrics]
    ht = ht.group_by(ht.qc_pop, ht.qc_platform).aggregate(**dict(zip(key_expr, value_expr)))

    upper_key_expr = ['{0}_upper'.format(metric) for metric in qc_metrics]
    upper_value_expr = [ht['{0}_median'.format(metric)] + 4 * ht['{0}_mad'.format(metric)] if metric != 'callrate' else 1 for
                        metric in qc_metrics]

    lower_key_expr = ['{0}_lower'.format(metric) for metric in qc_metrics]
    lower_value_expr = [ht['{0}_median'.format(metric)] - 4 * ht['{0}_mad'.format(metric)] if metric != 'callrate' else 0.99 for
                        metric in qc_metrics]

    annotation = dict(zip(upper_key_expr, upper_value_expr))
    annotation.update(dict(zip(lower_key_expr, lower_value_expr)))
    ht = ht.annotate(**annotation)
    ht = ht.annotate(idx=ht.qc_pop + "_" + hl.str(ht.qc_platform))
    ht = ht.filter(hl.is_defined(ht.idx))
    ht.export(qc_temp_data_prefix(data_type) + '.platform_pop_MAD.txt.bgz')
    return ht
# TODO: investigate why some idx categories are NA (i.e., why some pop/platform information is missing)


def apply_pop_platform_filters(ht: hl.Table, mt: hl.MatrixTable, qc_metrics: List[str]) -> hl.MatrixTable:
    """
    Flag samples that fall outside pop- and platform-specific QC metric thresholds

    :param Table ht: Table containing upper- and lower-threshold values for each pop- and platform-specific cohort of samples
    :param MatrixTable mt: MT containing samples to which pop- and platform-specific outlier filtering is applied
    :param list qc_metrics: list of metrics on which to flag outlier samples
    :return: MatrixTable containing the pop_platform_filter annotation summarizing metrics on which a given sample is an outlier
    :rtype: MatrixTable
    """
    thresh_colnames = ['idx'] + ['{}_upper'.format(metric) for metric in qc_metrics] + [
        '{}_lower'.format(metric) for metric in qc_metrics]
    ht = ht.select(*thresh_colnames).key_by('idx')
    ht.write('/tmp.h', True)  # temporary write/read for hail assert bug
    ht = hl.read_table('/tmp.h')
    mt = mt.annotate_cols(_thresholds=ht[mt.idx])
    fail_exprs = {
        'fail_{}'.format(metric):
            (mt.sample_qc['{}'.format(metric)] >= mt._thresholds['{}_upper'.format(metric)]) |
            (mt.sample_qc['{}'.format(metric)] <= mt._thresholds['{}_lower'.format(metric)]) for metric in qc_metrics}
    mt = mt.annotate_cols(**fail_exprs).drop('_thresholds')
    pop_platform_filters = make_pop_filters_expr(mt, qc_metrics)
    return mt.annotate_cols(pop_platform_filters=pop_platform_filters).drop('idx')


def split_mt_by_relatedness(pruned_mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    """
    Split gnomAD Mt into unrelated and related MTs (based on hard-coded data paths)
    """
    related_samples_to_drop_ranked = hl.read_table(qc_temp_data_prefix('joint') + '.related_samples_to_drop.ht')
    rank_table = hl.import_table(rank_annotations_path('joint'), impute=True).key_by('data_type', 's')
    pruned_mt = pruned_mt.annotate_cols(
        drop_unrelated_ranked=hl.is_defined(related_samples_to_drop_ranked[pruned_mt.col_key]),
        **rank_table[pruned_mt.col_key])
    pca_mt = pruned_mt.filter_cols(pruned_mt.drop_unrelated_ranked, keep=False).annotate_cols(related=False)
    related_mt = pruned_mt.filter_cols(pruned_mt.drop_unrelated_ranked, keep=True).annotate_cols(related=True)
    return pca_mt, related_mt


def main(args):
    hl.init(log='/sample_qc.log', tmp_dir='hdfs:///pc_relate.tmp/')

    if not args.load_joint_pruned_qc_mt:
        logger.info('Joining exomes and genomes...')
        exome_qc_mt = read_and_pre_process_data(qc_mt_path('exomes'), qc_ht_path('exomes', 'hard_filters'))
        genome_qc_mt = read_and_pre_process_data(qc_mt_path('genomes'), qc_ht_path('genomes', 'hard_filters'))

        joint_qc_mt = exome_qc_mt.union_cols(genome_qc_mt)  # NOTE: this is an inner join on rows
        joint_qc_mt = joint_qc_mt.filter_rows((hl.agg.mean(joint_qc_mt.GT.n_alt_alleles()) / 2 > 0.001) &
                                              (hl.agg.fraction(hl.is_defined(joint_qc_mt.GT)) > 0.99))
        joint_qc_mt.write(qc_mt_path('joint'), args.overwrite)

        logger.info('LD-pruning joint mt of exomes and genomes...')
        joint_qc_mt = hl.read_matrix_table(qc_mt_path('joint'))
        variants, samples = joint_qc_mt.count()
        logger.info('Pruning {0} variants in {1} samples'.format(variants, samples))
        joint_qc_pruned_mt = hl.ld_prune(joint_qc_mt, r2=0.1, n_cores=args.num_cores)
        joint_qc_pruned_mt.write(qc_mt_path('joint', ld_pruned=True), args.overwrite)

    pruned_mt = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True))
    variants, samples = pruned_mt.count()
    logger.info('{0} samples, {1} variants found in LD-pruned joint MT'.format(samples, variants))

    if not args.skip_pc_relate:
        logger.info('Running PCA for PC-Relate...')
        eig, scores, _ = hl.hwe_normalized_pca(pruned_mt.GT, k=10, compute_loadings=False)
        scores.write(qc_temp_data_prefix('joint') + '.pruned.pca_scores.ht', args.overwrite)

        logger.info('Running PC-Relate...')
        scores = hl.read_table(qc_temp_data_prefix('joint') + '.pruned.pca_scores.ht')
        # NOTE: This needs SSDs on your workers (for the temp files) and no pre-emptibles while the BlockMatrix writes
        relatedness_ht = hl.pc_relate(pruned_mt.GT, min_individual_maf=0.05, scores_expr=scores[pruned_mt.col_key].scores,
                                      block_size=4096, min_kinship=0.05, statistics='kin2')
        relatedness_ht.write(relatedness_ht_path, args.overwrite)

    relatedness_ht = hl.read_table(relatedness_ht_path)

    if not args.skip_relatedness:
        infer_ped(GnomADRelatedData('exomes'))
        infer_ped(GnomADRelatedData('genomes'))

        logger.info('Making rank file...')
        rank_table = make_rank_file(rank_annotations_path('joint'))
        logger.info('Finished making rank file...')

        related_samples_to_drop_ranked = get_related_samples_to_drop(rank_table, relatedness_ht)
        related_samples_to_drop_ranked.write(qc_temp_data_prefix('joint') + '.related_samples_to_drop.ht', args.overwrite)

    pca_mt, related_mt = split_mt_by_relatedness(pruned_mt)

    if not args.skip_pop_pca:
        variants, samples = pca_mt.count()
        logger.info('{} samples after removing relateds'.format(samples))
        # TODO: Check that there are no longer any 2nd-degree relateds in the callset by running KING on the output file below
        plink_mt = pca_mt.annotate_cols(uid=pca_mt.data_type + '_' + pca_mt.s.replace(" ", "_")).replace("/", "_").key_cols_by('uid')
        hl.export_plink(plink_mt, qc_temp_data_prefix('joint') + '.unrelated.plink', fam_id=plink_mt.uid, ind_id=plink_mt.uid)

        logger.info('Computing population PCs and annotating with known population labels...')
        pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(pca_mt.GT, k=20, compute_loadings=True)
        pca_mt = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2)
        pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)
        pca_scores.write(ancestry_pca_scores_ht_path(), args.overwrite)
        pca_loadings.write(ancestry_pca_loadings_ht_path(), args.overwrite)

    pca_scores = hl.read_table(ancestry_pca_scores_ht_path())
    pca_loadings = hl.read_table(ancestry_pca_loadings_ht_path())
    pca_mt = pca_mt.annotate_cols(scores=pca_scores[pca_mt.col_key].scores)

    variants, samples = related_mt.count()
    logger.info('Projecting population PCs for {} related samples...'.format(samples))
    related_scores = pc_project(related_mt, pca_loadings)
    relateds = related_mt.cols()
    relateds = relateds.annotate(scores=related_scores[relateds.key].scores)

    logger.info('Assigning population annotations...')
    pop_colnames = ['related', 'known_pop', 'scores']
    pop_annots_ht = hl.import_table(known_population_annotations, impute=True).key_by('combined_sample')

    joint_ht = pca_mt.cols().union(relateds)
    joint_ht = joint_ht.annotate(known_pop=pop_annots_ht[joint_ht.data_type.replace('s', '') + '_' + joint_ht.s.replace(' ', '_')].known_pop) # FIXME: temporarily doing the underscore thing until known_population_annotations is fixed
    joint_pca_ht = joint_ht.select(*pop_colnames)
    joint_pca_ht, joint_pca_fit = run_assign_population_pcs(joint_pca_ht, qc_temp_data_prefix('joint') + '.RF_pop_assignments.txt.bgz', qc_temp_data_prefix('joint') + '.RF_fit.pkl', pcs=list(range(1,7)))
    joint_ht = joint_ht.annotate(pop=joint_pca_ht[joint_ht.key].pop).select('pop', *pop_colnames)

    # Add special Estonian pop category for genomes
    estonian_ht = (hl.import_table(estonian_batches, impute=True)
                   .annotate(data_type='genomes').key_by('data_type', 'sample'))
    joint_ht = joint_ht.annotate(batch=estonian_ht[joint_ht.key].batch)
    joint_ht = joint_ht.annotate(qc_pop=hl.case(missing_false=True)
                                 .when(hl.is_defined(joint_ht.pop) & (joint_ht.batch == 1), 'est_b1')
                                 .when(hl.is_defined(joint_ht.pop) & (joint_ht.batch == 2), 'est_b2')
                                 .default(joint_ht.pop))

    # These are keyed by only `s`
    genome_mt = get_gnomad_data('genomes', adj=False, split=False, meta_root=None).select_cols()
    exome_mt = get_gnomad_data('exomes', adj=False, split=False, meta_root=None).select_cols()

    # Population-specific filtering
    if not args.skip_calculate_sample_metrics:
        logger.info('Running mini sample QC for platform- and population-specific filtering...')
        gnomad_sample_qc(exome_mt).cols().select('sample_qc').write(qc_temp_data_prefix('exomes') + '.sample_qc.ht', args.overwrite)
        gnomad_sample_qc(genome_mt).cols().select('sample_qc').write(qc_temp_data_prefix('genomes') + '.sample_qc.ht', args.overwrite)
        # TODO: check that the pcr_free annotations are complete once samples are updated from Jessica's spreadsheet

    logger.info('Annotating population and platform assignments...')
    platform_ht = hl.read_table(qc_ht_path('exomes', 'platforms'))
    exome_mt = exome_mt.annotate_cols(qc_platform=platform_ht.key_by('s')[exome_mt.s].qc_platform,
                                      **joint_ht.filter(joint_ht.data_type == 'exomes').key_by('s')[exome_mt.s])

    genome_meta_ht = hl.read_table(qc_ht_path('genomes', 'hard_filters'))
    genome_mt = genome_mt.annotate_cols(qc_platform=genome_meta_ht.key_by('s')[genome_mt.s].qc_platform,
                                        **joint_ht.filter(joint_ht.data_type == 'genomes').key_by('s')[genome_mt.s])

    exome_sample_qc_ht = hl.read_table(qc_temp_data_prefix('exomes') + '.sample_qc.ht')
    genome_sample_qc_ht = hl.read_table(qc_temp_data_prefix('genomes') + '.sample_qc.ht')

    exome_mt = exome_mt.annotate_cols(**exome_sample_qc_ht[exome_mt.s])
    genome_mt = genome_mt.annotate_cols(**genome_sample_qc_ht[genome_mt.s])
    exome_mt = exome_mt.annotate_cols(idx=exome_mt.qc_pop + "_" + hl.str(exome_mt.qc_platform))
    genome_mt = genome_mt.annotate_cols(idx=genome_mt.qc_pop + "_" + hl.str(genome_mt.qc_platform))

    # For each population, aggregate sample QC metrics and calculate the MAD/mean/stdev
    logger.info('Calculating platform- and population-specific sample QC thresholds...')
    exome_qc_metrics = ['n_snp', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var']
    exome_ht = calculate_metrics(exome_mt, exome_qc_metrics, 'exomes')

    genome_qc_metrics = exome_qc_metrics + ['call_rate']
    genome_ht = calculate_metrics(genome_mt, genome_qc_metrics, 'genomes')

    # Annotate samples that fail their respective filters
    logger.info('Flagging samples failing pop/platform-specific sample qc thresholds...')
    exome_mt = apply_pop_platform_filters(exome_ht, exome_mt, exome_qc_metrics)
    checkpoint = exome_mt.aggregate_cols(hl.agg.count_where(hl.len(exome_mt.pop_platform_filters) == 0))
    logger.info('{0} exome samples found passing pop/platform-specific filtering'.format(checkpoint))
    exome_mt.cols().annotate(data_type='exomes').key_by('data_type', 's').write(qc_ht_path('exomes', 'pop_platform'), args.overwrite)

    genome_mt = apply_pop_platform_filters(genome_ht, genome_mt, genome_qc_metrics)
    checkpoint = genome_mt.aggregate_cols(hl.agg.count_where(hl.len(genome_mt.pop_platform_filters) == 0))
    logger.info('{0} genome samples found passing pop/platform-specific filtering'.format(checkpoint))
    genome_mt.cols().annotate(data_type='genomes').key_by('data_type', 's').write(qc_ht_path('genomes', 'pop_platform'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--load_joint_pruned_qc_mt', help='Load a pre-computed LD-pruned joint MT instead of rewriting', action='store_true')
    parser.add_argument('--num_cores', help='Required for LD Prune', default=32)
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')

    parser.add_argument('--skip_platform_pca', help='Skip annotating and assigning platform PCs', action='store_true')
    parser.add_argument('--skip_pc_relate', help='Skip PC relate calculations', action='store_true')
    parser.add_argument('--skip_relatedness', help='Skip calculating relatedness', action='store_true')
    parser.add_argument('--skip_pop_pca', help='Skip calculating population PCs on unrelated samples', action='store_true')
    parser.add_argument('--skip_calculate_sample_metrics', help='Skip calculating sample metrics (sample_qc)', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        try_slack(args.slack_channel, main, args)
    else:
        main(args)
