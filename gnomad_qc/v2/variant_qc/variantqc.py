from gnomad.variant_qc.evaluation import add_rank
from gnomad.utils.slack import slack_notifications
from gnomad.variant_qc import random_forest
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.variant_qc import *
from pprint import pformat
import json
import uuid
import sys
import argparse
import logging

logger = logging.getLogger("variant_qc")

LABEL_COL = 'rf_label'
TRAIN_COL = 'rf_train'
prediction_col_name = 'rf_prediction'

SITES_FEATURES = [
    'info_MQRankSum',
    'info_SOR',
    'info_InbreedingCoeff',
    'info_ReadPosRankSum'
]

VQSR_FEATURES = [
    'info_FS',
    'info_QD',
    'info_MQ',
    'info_DP'
]

ALLELE_FEATURES = [
    'variant_type',
    'allele_type',
    'n_alt_alleles',
    'was_mixed',
    'has_star',
    'qd',
    'pab_max'
]

MEDIAN_FEATURES = [
    'gq_median',
    'dp_median',
    'nrq_median',
    'ab_median'
]

INBREEDING_COEFF_HARD_CUTOFF = -0.3


def get_features_list(sites_features: bool, allele_features: bool, vqsr_features: bool, median_features: bool = False) -> List[str]:
    """
    Returns the list of features to use based on desired arguments (currently only VQSR / alleles)

    :param bool sites_features: Whether to use site-level features
    :param bool vqsr_features: Whether to use VQSR features
    :param bool allele_features: Whether to use Allele-specific features
    :param bool median_features: Whether to use 2.0.2 features
    :return: List of features to use
    """

    features = []
    if sites_features:
        features.extend(SITES_FEATURES)
    if allele_features:
        features.extend(ALLELE_FEATURES)
    if vqsr_features:
        features.extend(VQSR_FEATURES)
    if median_features:
        features.extend(MEDIAN_FEATURES)

    return features


def sample_rf_training_examples(
        ht: hl.Table,
        tp_col: str,
        fp_col: str,
        fp_to_tp: float = 1.0,
        label_col: str = LABEL_COL,
        train_col: str = TRAIN_COL
) -> hl.Table:
    """

    Adds an annotation to a Hail Table that defines which rows should be used for RF training, given
    a set of positive and negative training examples and a TP to FP ratio.
    Examples that are both TP and FP are never kept.

    Note that it expects this Table to be from a split MT.

    :param Table ht: Input Table
    :param str tp_col: TP examples column
    :param str fp_col: FP examples column
    :param int fp_to_tp: FP to TP ratio. If set to <= 0, all training examples are used.
    :param str label_col: Name of column to store the training label
    :param str train_col: Name of column to store whether the site is used for training or not
    :return: Table subset with corresponding TP and FP examples with desired FP to TP ratio.
    :rtype: Table
    """

    def get_train_counts(ht: hl.Table) -> Tuple[int, int]:

        if 'test_intervals' in ht.globals:
            interval_list = hl.eval(ht.globals.test_intervals)
            ht = hl.filter_intervals(ht, interval_list, keep=False)

        # Get stats about TP / FP sets
        train_stats = hl.struct(
            n=hl.agg.count()
        )

        if 'alleles' in ht.row and ht.row.alleles.dtype == hl.tarray(hl.tstr):
            train_stats = train_stats.annotate(
                ti=hl.agg.count_where(hl.expr.is_transition(ht.alleles[0], ht.alleles[1])),
                tv=hl.agg.count_where(hl.expr.is_transversion(ht.alleles[0], ht.alleles[1])),
                indel=hl.agg.count_where(hl.expr.is_indel(ht.alleles[0], ht.alleles[1]))
            )

        # Sample training examples
        pd_stats = ht.group_by(**{'contig': ht.locus.contig, tp_col: ht[tp_col], fp_col: ht[fp_col]}).aggregate(**train_stats).to_pandas()
        logger.info(pformat(pd_stats))
        pd_stats = pd_stats.fillna(False)
        # logger.info(pformat(pd_stats))

        n_tp = pd_stats[pd_stats[tp_col] & ~pd_stats[fp_col]]['n'].sum()
        n_fp = pd_stats[~pd_stats[tp_col] & pd_stats[fp_col]]['n'].sum()

        return n_tp, n_fp

    n_tp, n_fp = get_train_counts(ht)

    if fp_to_tp > 0:

        desired_fp = fp_to_tp * n_tp
        if desired_fp < n_fp:
            prob_tp = 1.0
            prob_fp = desired_fp / n_fp
        else:
            prob_tp = n_fp / desired_fp
            prob_fp = 1.0

        logger.info(f"Training examples sampling: tp={prob_tp}*{n_tp}, fp={prob_fp}*{n_fp}")

        # Due to hail bug, run in two passes!

        if prob_fp < 1.0:
            ht = ht.annotate(**{train_col: hl.cond(hl.or_else(ht[tp_col], False), hl.or_else(~ht[fp_col], True), ht[fp_col])})
            train_expr = hl.cond(ht[fp_col] & hl.or_else(~ht[tp_col], True), hl.rand_bool(prob_fp),  ht[train_col])
            #train_expr = hl.cond(hl.or_else(ht[tp_col], False), hl.or_else(~ht[fp_col], True), ht[fp_col] & hl.rand_bool(prob_fp))  # Note: Hail propagates missing values with invert operator
        elif prob_tp < 1.0:
            ht = ht.annotate(**{train_col: hl.cond(hl.or_else(ht[fp_col], False), hl.or_else(~ht[tp_col], True), ht[tp_col])})
            train_expr = hl.cond(ht[tp_col] & hl.or_else(~ht[fp_col], True), hl.rand_bool(prob_tp), ht[train_col])
            #train_expr = hl.cond(hl.or_else(ht[fp_col], False), hl.or_else(~ht[tp_col], True), ht[tp_col] & hl.rand_bool(prob_tp))
        else:
            train_expr = (ht[tp_col] ^ ht[fp_col])
    else:
        logger.info(f"Using all {n_tp} TP and {n_fp} FP training examples.")
        train_expr = (ht[tp_col] ^ ht[fp_col])

    label_expr = (hl.case(missing_false=True)
                  .when(ht[tp_col] & hl.or_else(~ht[fp_col], True), "TP")
                  .when(ht[fp_col] & hl.or_else(~ht[tp_col], True), "FP")
                  .default(hl.null(hl.tstr)))

    if 'test_intervals' in ht.globals:
        train_expr = ~ht.globals.test_intervals.any(lambda interval: interval.contains(ht.locus)) & train_expr

    return ht.annotate(**{label_col: label_expr,
                          train_col: train_expr})


def create_rf_ht(
        data_type: str,
        n_variants_median: int = 50000,
        impute_features_by_variant_type: bool = True,
        group: str = 'qc_samples_raw'
) -> hl.Table:
    """
    Creates a Table with all necessary columns for RF:
    - Features columns (both VQSR and RF)
    - Training criteria columns

    Numerical features are median-imputed. If impute_features_by_variant_type is set, imputation is done based on the median of the variant type.

    :param str data_type: One of 'exomes' or 'genomes'
    :param int n_variants_median: Number of variants to use for median computation
    :param bool impute_features_by_variant_type: Whether to impute features median by variant type
    :param str group: Whether to use 'raw' or 'adj' genotypes
    :return: Hail Table ready for RF
    :rtype: Table
    """

    def get_site_features_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to annotate site-level features used by both VQSR and RF.

        :param Table ht: Table to create annotation expression for.
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return dict(zip(SITES_FEATURES, [
            ht.info.MQRankSum,
            ht.info.SOR,
            ht.info.InbreedingCoeff,
            ht.info.ReadPosRankSum
        ]))

    def get_vqsr_features_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to annotate site-level features only used by VQSRs

        :param Table ht: Table to create annotation expression for.
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return dict(zip(VQSR_FEATURES, [
            ht.info.FS,
            ht.info.QD,
            ht.info.MQ,
            ht.info.DP
        ]))

    def get_allele_features_expr(ht: hl.Table, qc_stats_group_index: int) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to annotate allele-level features (RF only)

        :param Table ht: Table to create annotation expression for.
        :param int qc_stats_group_index: Index of group to get stats for
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return dict(zip(ALLELE_FEATURES, [
            ht.allele_data.variant_type,
            ht.allele_data.allele_type,
            ht.allele_data.n_alt_alleles,
            ht.allele_data.was_mixed,
            ht.allele_data.has_star,
            ht.qc_stats[qc_stats_group_index].qd,
            ht.qc_stats[qc_stats_group_index].pab.max
        ]))

    def get_median_features_expr(ht: hl.Table, qc_stats_group_index: int) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to annotate 2.0.2 allele-specific features (RF only)

        :param Table ht: Table to create annotation expression for.
        :param int qc_stats_group_index: Index of group to get stats for
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return dict(zip(MEDIAN_FEATURES, [
            ht.qc_stats[qc_stats_group_index].gq_median,
            ht.qc_stats[qc_stats_group_index].dp_median,
            ht.qc_stats[qc_stats_group_index].nrq_median,
            ht.qc_stats[qc_stats_group_index].ab_median
        ]))

    def get_training_sites_expr(ht: hl.Table, family_stats_group_index: int) -> Dict[str, hl.expr.Expression]:
        """
        Returns expressions to columns to select training examples

        :param Table ht: Table to create annotation expression for.
        :param int family_stats_group_index: Index of group to get stats for
        :return: Dict with keys containing column names and values expressions
        :rtype: Dict of str: Expression
        """
        return {
            'transmitted_singleton': (ht.family_stats[family_stats_group_index].tdt.t == 1) &
                                     (ht.family_stats[family_stats_group_index].unrelated_qc_callstats.AC[1] == 1),

            'fail_hard_filters': (ht.info.QD < 2) | (ht.info.FS > 60) | (ht.info.MQ < 30),
            'info_POSITIVE_TRAIN_SITE': ht.info.POSITIVE_TRAIN_SITE,
            'info_NEGATIVE_TRAIN_SITE': ht.info.NEGATIVE_TRAIN_SITE
        }

    mt = get_gnomad_data(data_type)
    ht_family_stats = hl.read_table(annotations_ht_path(data_type, 'family_stats'))
    ht_qc_stats = hl.read_table(annotations_ht_path(data_type, 'qc_stats'))
    ht_call_stats = hl.read_table(annotations_ht_path(data_type, 'call_stats'))
    ht_truth_data = hl.read_table(annotations_ht_path(data_type, 'truth_data'))
    ht_allele_data = hl.read_table(annotations_ht_path(data_type, 'allele_data'))
    call_stats_data = hl.read_table(annotations_ht_path(data_type, 'call_stats'))

    ht = mt.annotate_rows(
        n_nonref=call_stats_data[mt.row_key].qc_callstats[0].AC[1],
        singleton=mt.info.AC[mt.a_index - 1] == 1,
        info_ac=mt.info.AC[mt.a_index - 1],
        **ht_call_stats[mt.row_key],
        **ht_family_stats[mt.row_key],
        **ht_qc_stats[mt.row_key],
        **ht_truth_data[mt.row_key],
        **ht_allele_data[mt.row_key]
    ).rows()

    ht = ht.filter(ht.qc_callstats.find(lambda x: x.meta.get('group') == group).AC[1] > 0)

    qc_stats_groups = ht_qc_stats.aggregate(hl.agg.take(ht_qc_stats.qc_stats.map(lambda x: x.meta['group']), 1))[0]
    qc_stats_group_index = next(x[0] for x in enumerate(qc_stats_groups) if x[1] == group)
    family_stats_groups = ht_family_stats.aggregate(hl.agg.take(ht_family_stats.family_stats.map(lambda x: x.meta['group']), 1))[0]
    family_stats_group_index = next(x[0] for x in enumerate(family_stats_groups) if x[1] == 'raw')
    ht = ht.select(
        **get_allele_features_expr(ht, qc_stats_group_index),
        **get_site_features_expr(ht),
        **get_vqsr_features_expr(ht),
        **get_training_sites_expr(ht, family_stats_group_index),
        **get_median_features_expr(ht, qc_stats_group_index),
        omni=ht.truth_data.omni,
        mills=ht.truth_data.mills,
        n_nonref=ht.n_nonref,
        singleton=ht.singleton,
        was_split=ht.was_split,
        info_ac=ht.info_ac
    )

    # Annotate variants by type (e.g. for downsampling purpose)
    variants_by_type = ht.aggregate(hl.agg.counter(ht.variant_type))
    logger.info('Variants by type:\n{}'.format('\n'.join(['{}: {}'.format(k, v) for k, v in variants_by_type.items()])))

    ht = ht.annotate_globals(variants_by_type=variants_by_type)
    ht = ht.repartition(1000, shuffle=False)
    ht = ht.persist()

    # Compute medians
    # numerical_features = [k for k, v in ht.row.dtype.items() if annotation_type_is_numeric(v)]
    numerical_features = [k for k, v in ht.row.dtype.items() if v == hl.tint or v == hl.tfloat]
    # print(numerical_features)

    if impute_features_by_variant_type:
        prob_sample = hl.literal({v: min([n, n_variants_median])/n for v, n in variants_by_type.items()})
        ht_by_variant = ht.group_by(ht.variant_type).partition_hint(1)
        medians = ht_by_variant.aggregate(
            **{feature: hl.median(
                hl.agg.filter(
                    hl.rand_bool(prob_sample[ht.variant_type]),
                    hl.agg.collect(ht[feature])
                )
            ) for feature in numerical_features}
        ).collect()

        ht = ht.annotate_globals(feature_medians={x.variant_type: x for x in medians})
        ht = ht.annotate(
            **{f: hl.or_else(ht[f], ht.globals.feature_medians[ht.variant_type][f]) for f in numerical_features},
            feature_imputed=hl.struct(
                **{f: hl.is_missing(ht[f]) for f in numerical_features}
            )
        )
    else:
        prob_sample = min([sum(variants_by_type.values()), n_variants_median]) / sum(variants_by_type.values())
        medians = ht.aggregate(
            hl.struct(
                **{feature: hl.median(
                    hl.agg.filter(
                        hl.rand_bool(prob_sample),
                        hl.agg.collect(ht[feature])
                    )
                ) for feature in numerical_features}
            )
        )
        ht = ht.annotate_globals(features_median=medians)
        ht = ht.annotate(
            **{f: hl.or_else(ht[f], ht.globals.features_median[f]) for f in numerical_features},
            feature_imputed=hl.struct(
                **{f: hl.is_missing(ht[f]) for f in numerical_features}
            )
        )

    return ht


def get_run_data(
        data_type: str,
        vqsr_training: bool,
        no_transmitted_singletons: bool,
        adj: bool,
        test_intervals: List[str],
        features_importance: Dict[str, float],
        test_results: List[hl.tstruct] = None
) -> Dict:
    run_data = {
        'data': data_type,
        'input_args': {
            'vqsr_training': vqsr_training,
            'no_transmitted_singletons': no_transmitted_singletons,
            'adj': adj
        },
        'features_importance': features_importance,
        'test_intervals': test_intervals
    }

    if test_results:
        tps = 0
        total = 0
        for row in test_results:
            values = list(row.values())
            if values[0] == values[1]:
                tps += values[2]
            total += values[2]
        run_data['test_results'] = [dict(x) for x in test_results]
        run_data['test_accuracy'] = tps / total

    return run_data


def train_rf(data_type, args):

    # Get unique hash for run and load previous runs
    run_hash = str(uuid.uuid4())[:8]
    rf_runs = get_rf_runs(data_type)
    while run_hash in rf_runs:
        run_hash = str(uuid.uuid4())[:8]

    ht = hl.read_table(rf_annotated_path(data_type, args.adj))

    ht = ht.repartition(500, shuffle=False)

    if not args.vqsr_training:
        if args.no_transmitted_singletons:
            tp_expr = ht.omni | ht.mills | ht.info_POSITIVE_TRAIN_SITE
        else:
            tp_expr = ht.omni | ht.mills | ht.info_POSITIVE_TRAIN_SITE | ht.transmitted_singleton

        ht = ht.annotate(
            tp=tp_expr
        )

    test_intervals_str = [] if not args.test_intervals else [args.test_intervals] if isinstance(args.test_intervals, str) else args.test_intervals
    test_intervals_locus = [hl.parse_locus_interval(x) for x in test_intervals_str]

    if test_intervals_locus:
        ht = ht.annotate_globals(
            test_intervals=test_intervals_locus
        )

    ht = sample_rf_training_examples(ht,
                                     tp_col='info_POSITIVE_TRAIN_SITE' if args.vqsr_training else 'tp',
                                     fp_col='info_NEGATIVE_TRAIN_SITE' if args.vqsr_training else 'fail_hard_filters',
                                     fp_to_tp=args.fp_to_tp)

    ht = ht.persist()

    rf_ht = ht.filter(ht[TRAIN_COL])
    rf_features = get_features_list(True,
                                    not (args.vqsr_features or args.median_features),
                                    args.vqsr_features,
                                    args.median_features)

    logger.info("Training RF model:\nfeatures: {}\nnum_tree: {}\nmax_depth:{}\nTest intervals: {}".format(
        ",".join(rf_features),
        args.num_trees,
        args.max_depth,
        ",".join(test_intervals_str)))

    rf_model = random_forest.train_rf(rf_ht,
                           features=rf_features,
                           label=LABEL_COL,
                           num_trees=args.num_trees,
                           max_depth=args.max_depth)

    logger.info("Saving RF model")
    random_forest.save_model(rf_model,
                  rf_path(data_type, data='model', run_hash=run_hash),
                  overwrite=args.overwrite
                  )

    test_results = None
    if args.test_intervals:
        logger.info("Testing model {} on intervals {}".format(run_hash, ",".join(test_intervals_str)))
        test_ht = hl.filter_intervals(ht, test_intervals_locus, keep=True)
        test_ht = test_ht.checkpoint('gs://gnomad-tmp/test_random_forest.ht', overwrite=True)
        test_ht = test_ht.filter(hl.is_defined(test_ht[LABEL_COL]))
        test_results = random_forest.test_model(test_ht,
                                     rf_model,
                                     features=get_features_list(True, not args.vqsr_features, args.vqsr_features),
                                     label=LABEL_COL)
        ht = ht.annotate_globals(
            test_results=test_results
        )

    logger.info("Writing RF training HT")
    features_importance = random_forest.get_features_importance(rf_model)
    ht = ht.annotate_globals(
        features_importance=features_importance,
        features=get_features_list(True, not args.vqsr_features, args.vqsr_features),
        vqsr_training=args.vqsr_training,
        no_transmitted_singletons=args.no_transmitted_singletons,
        adj=args.adj
    )
    ht.write(rf_path(data_type, data='training', run_hash=run_hash), overwrite=args.overwrite)

    logger.info("Adding run to RF run list")
    rf_runs[run_hash] = get_run_data(
        data_type,
        args.vqsr_training,
        args.no_transmitted_singletons,
        args.adj,
        test_intervals_str,
        features_importance,
        test_results
    )
    with hl.hadoop_open(rf_run_hash_path(data_type), 'w') as f:
        json.dump(rf_runs, f)

    return run_hash


def prepare_final_ht(data_type: str, run_hash: str, snp_bin_cutoff: int, indel_bin_cutoff: int) -> hl.Table:

    # Get snv and indel RF cutoffs based on bin
    binned_ht_path = score_ranking_path(data_type, run_hash, binned=True)
    if not hl.hadoop_exists(score_ranking_path(data_type, run_hash, binned=True)):
        sys.exit(f"Could not find binned HT for RF  run {run_hash} ({binned_ht_path}). Please run create_ranked_scores.py for that hash.")
    binned_ht = hl.read_table(binned_ht_path)
    snp_rf_cutoff, indel_rf_cutoff = binned_ht.aggregate([hl.agg.filter(binned_ht.snv & (binned_ht.bin == snp_bin_cutoff), hl.agg.min(binned_ht.min_score)),
                                                          hl.agg.filter(~binned_ht.snv & (binned_ht.bin == indel_bin_cutoff), hl.agg.min(binned_ht.min_score))])

    # Add filters to RF HT
    ht = hl.read_table(rf_path(data_type, 'rf_result', run_hash=run_hash))
    ht = ht.annotate_globals(rf_hash=run_hash,
                             rf_snv_cutoff=hl.struct(bin=snp_bin_cutoff, min_score=snp_rf_cutoff),
                             rf_indel_cutoff=hl.struct(bin=indel_bin_cutoff, min_score=indel_rf_cutoff))
    rf_filter_criteria = (hl.is_snp(ht.alleles[0], ht.alleles[1]) & (ht.rf_probability['TP'] < ht.rf_snv_cutoff.min_score)) | (
            ~hl.is_snp(ht.alleles[0], ht.alleles[1]) & (ht.rf_probability['TP'] < ht.rf_indel_cutoff.min_score))
    ht = ht.annotate(filters=hl.case()
                     .when(rf_filter_criteria, {'RF'})
                     .when(~rf_filter_criteria, hl.empty_set(hl.tstr))
                     .or_error('Missing RF probability!'))

    inbreeding_coeff_filter_criteria = hl.is_defined(ht.info_InbreedingCoeff) & (
            ht.info_InbreedingCoeff < INBREEDING_COEFF_HARD_CUTOFF)
    ht = ht.annotate(filters=hl.cond(inbreeding_coeff_filter_criteria,
                                     ht.filters.add('InbreedingCoeff'), ht.filters))

    freq_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies'))
    ac0_filter_criteria = freq_ht[ht.key].freq[0].AC[1] == 0

    ht = ht.annotate(
        filters=hl.cond(ac0_filter_criteria, ht.filters.add('AC0'), ht.filters)
    )

    # Fix annotations for release
    annotations_expr = {
        'tp': hl.or_else(ht.tp, False),
        'transmitted_singleton': hl.or_missing(freq_ht[ht.key].freq[1].AC[1] == 1, ht.transmitted_singleton),
        'rf_probability': ht.rf_probability["TP"]
    }
    if 'feature_imputed' in ht.row:
        annotations_expr.update(
            {x: hl.or_missing(~ht.feature_imputed[x], ht[x]) for x in [f for f in ht.row.feature_imputed]}
        )

    ht = ht.transmute(
        **annotations_expr
    )

    #This column is added by the RF module based on a 0.5 threshold which doesn't correspond to what we use
    ht = ht.drop(ht.rf_prediction)

    return ht


def main(args):
    hl.init(log='/variantqc.log')

    data_type = 'exomes' if args.exomes else 'genomes'

    if args.debug:
        logger.setLevel(logging.DEBUG)

    if args.list_rf_runs:
        logger.info(f"RF runs for {data_type}:")
        random_forest.pretty_print_runs(get_rf_runs(data_type))

    if args.annotate_for_rf:
        ht = create_rf_ht(data_type,
                          n_variants_median=args.n_variants_median,
                          impute_features_by_variant_type=not args.impute_features_no_variant_type,
                          group='adj' if args.adj else 'raw')
        ht.write(rf_annotated_path(data_type, args.adj), overwrite=args.overwrite)

    run_hash = train_rf(data_type, args) if args.train_rf else args.run_hash

    if args.apply_rf:
        logger.info(f"Applying RF model {run_hash} to {data_type}.")

        rf_model = random_forest.load_model(rf_path(data_type, data='model', run_hash=run_hash))
        ht = hl.read_table(rf_path(data_type, data='training', run_hash=run_hash))

        ht = random_forest.apply_rf_model(ht, rf_model, get_features_list(True, not args.vqsr_features, args.vqsr_features), label=LABEL_COL)

        if 'singleton' in ht.row and 'was_split' in ht.row:  # Needed for backwards compatibility for RF runs that happened prior to updating annotations
            ht = add_rank(ht,
                          score_expr=ht.rf_probability['FP'],
                          subrank_expr={
                              'singleton_rank': ht.singleton,
                              'biallelic_rank': ~ht.was_split,
                              'biallelic_singleton_rank': ~ht.was_split & ht.singleton
                          }
                      )
        else:
            logger.warn("Ranking was not added  because of missing annotations -- please run 'create_ranked_scores.py' to add rank.")

        ht.write(rf_path(data_type, 'rf_result', run_hash=run_hash), overwrite=args.overwrite)

    if args.finalize:
        ht = prepare_final_ht(data_type, args.run_hash, args.snp_bin_cutoff, args.indel_bin_cutoff)
        ht.write(annotations_ht_path(data_type, 'rf'), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--exomes', help='Run on exomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Run on genomes. One and only one of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--run_hash', help='Run hash. Created by --train_rf and only needed for --apply_rf without running --train_rf', required=False)
    parser.add_argument('--snp_bin_cutoff', help='Percentile to set RF cutoff', type=int, default=90)
    parser.add_argument('--indel_bin_cutoff', help='Percentile to set RF cutoff', type=int, default=80)
    parser.add_argument('--debug', help='Prints debug statements', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')

    actions = parser.add_argument_group('Actions')
    actions.add_argument('--list_rf_runs', help='Lists all previous RF runs, along with their hash,  parameters and testing results.', action='store_true')
    actions.add_argument('--annotate_for_rf', help='Creates an annotated ht with features for RF', action='store_true')
    actions.add_argument('--train_rf', help='Trains RF model', action='store_true')
    actions.add_argument('--apply_rf', help='Applies RF model to the data', action='store_true')
    actions.add_argument('--finalize', help='Write final RF model', action='store_true')

    annotate_params = parser.add_argument_group('Annotate features params')
    annotate_params.add_argument('--n_variants_median', help='Number of variants to use for median computation.', type=int, default=20000)
    annotate_params.add_argument('--impute_features_no_variant_type', help='If set, feature imputation is NOT split by variant type.', action='store_true')

    rf_params = parser.add_argument_group('Random Forest parameters')
    rf_params.add_argument('--fp_to_tp', help='Ratio of FPs to TPs for creating the RF model. If set to 0, all training examples are used. (default=1.0)', default=1.0,
                           type=float)
    rf_params.add_argument('--test_intervals', help='The specified interval(s) will be held out for testing and used for evaluation only. (default to "20")', nargs='+', type=str, default='20')
    rf_params.add_argument('--num_trees', help='Number of trees in the RF model. (default=500)', default=500, type=int)
    rf_params.add_argument('--max_depth', help='Maxmimum tree depth in the RF model. (default=5)', default=5, type=int)

    training_params = parser.add_argument_group('Training data parameters')
    training_params.add_argument('--adj', help='Use adj genotypes.', action='store_true')
    training_params.add_argument('--vqsr_features', help='Use VQSR features only (+ snv/indel variant type)',
                                 action='store_true')
    training_params.add_argument('--median_features', help='Use gnomAD 2.0.2 features',
                                 action='store_true')
    training_params.add_argument('--vqsr_training', help='Use VQSR training examples', action='store_true')
    training_params.add_argument('--no_transmitted_singletons', help='Do not use transmitted singletons for training.',
                                 action='store_true')

    args = parser.parse_args()
    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified.')

    if not args.run_hash and not args.train_rf and args.apply_rf:
        sys.exit('Error: --run_hash is required when running --apply_rf without running --train_rf too.')

    if args.run_hash and args.train_rf:
        sys.exit('Error: --run_hash and --train_rf are mutually exclusive. --train_rf will generate a run hash.')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
