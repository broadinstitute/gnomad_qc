"""Script for running random forest model on gnomAD v5 variant QC data."""

import argparse
import json
import logging
import math
import pickle
import pprint
import uuid
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
import numpy as np
import pandas as pd
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.variant_qc.pipeline import INFO_FEATURES, NON_INFO_FEATURES, train_rf_model
from gnomad.variant_qc.random_forest import (
    apply_rf_model,
    check_ht_fields_for_spark,
    get_rf_runs,
    get_run_data,
    load_model,
    pretty_print_runs,
    save_model,
)
from gnomad.variant_qc.training import sample_training_examples
from sklearn.ensemble import RandomForestClassifier

from gnomad_qc.v5.resources.annotations import get_variant_qc_annotations
from gnomad_qc.v5.resources.basics import (
    _check_resource_existence,
    _get_batch_resource_kwargs,
    _init_hail,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.sample_qc import DENSE_TRIO_TEST_CHROMS
from gnomad_qc.v5.resources.variant_qc import (
    _validate_model_id,
    get_rf_model_path,
    get_rf_run_path,
    get_rf_training,
    get_variant_qc_result,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("variant_qc_random_forest")
logger.setLevel(logging.INFO)

FEATURES = INFO_FEATURES + NON_INFO_FEATURES
"""Combined list of all features used in the RF model."""

RF_COLS = {
    "label": "rf_label",
    "prediction": "rf_prediction",
    "probability": "rf_probability",
    "train": "rf_train",
}
"""RF column names for training labels, predictions, probabilities, and train/test split."""


def load_variant_qc_annotations_ht(
    test: bool = False,
    environment: str = "rwb",
) -> hl.Table:
    """
    Load the pre-computed variant QC annotations Table for RF training.

    This Table is created by the generate_variant_qc_annotations.py script and contains
    all annotations needed for RF training:
    - Info features (in `info` struct): AS_MQRankSum, AS_pab_max, AS_MQ, AS_QD,
      AS_ReadPosRankSum, AS_SOR, AS_FS
    - Non-info features: variant_type, allele_type, n_alt_alleles, was_mixed, has_star
    - Truth data: hapmap, omni, mills, kgp_phase1_hc
    - Training labels: transmitted_singleton_raw, transmitted_singleton_adj,
      sibling_singleton_raw, sibling_singleton_adj, fail_hard_filters

    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "rwb".
    :return: Table with variant QC annotations for RF training.
    """
    logger.info("Loading variant QC annotations HT...")
    ht = get_variant_qc_annotations(test=test, environment=environment).ht()

    # Extract info features from the info struct for RF training.
    ht = ht.transmute(**{f: ht.info[f] for f in INFO_FEATURES})

    return ht


def _get_rf_features_defined_expr(
    ht: hl.Table, features: List[str]
) -> hl.expr.BooleanExpression:
    """
    Get whether all RF features of a row are defined and not NaN.

    Rows with a missing or NaN feature are dropped before RF training and scoring, as
    in ``ht_to_rf_df``.

    :param ht: Input Table.
    :param features: List of feature columns.
    :return: Expression that is True if no feature is missing or NaN.
    """
    return hl.all(
        [
            (
                hl.is_defined(ht[f]) & ~hl.is_nan(ht[f])
                if ht[f].dtype in (hl.tfloat32, hl.tfloat64)
                else hl.is_defined(ht[f])
            )
            for f in features
        ]
    )


def test_model_sklearn(
    ht: hl.Table,
    rf_model: Dict[str, Any],
    features: List[str],
    label: str,
    prediction_col_name: str = "rf_prediction",
) -> List[hl.Struct]:
    """
    Test a scikit-learn RF model on a set of examples with known labels.

    Mirrors ``test_model`` from gnomad_methods: runs the model on the data, logs the
    confusion matrix and accuracy, and returns the confusion matrix.

    :param ht: Input Table.
    :param rf_model: scikit-learn RF model.
    :param features: Columns containing features that were used in the model.
    :param label: Column containing label to be predicted.
    :param prediction_col_name: Where to store the prediction. Default is
        "rf_prediction".
    :return: A list containing structs with {label, prediction, n}.
    """
    ht = apply_rf_model_sklearn(
        ht.filter(hl.is_defined(ht[label])),
        rf_model,
        features,
        label,
        prediction_col_name=prediction_col_name,
    )

    test_results = (
        ht.group_by(ht[prediction_col_name], ht[label])
        .aggregate(n=hl.agg.count())
        .collect()
    )

    # Print results.
    df = pd.DataFrame(test_results)
    df = df.pivot(index=label, columns=prediction_col_name, values="n")
    logger.info("Testing results:\n%s", pprint.pformat(df))
    logger.info(
        "Accuracy: %f",
        sum([x.n for x in test_results if x[label] == x[prediction_col_name]])
        / sum([x.n for x in test_results]),
    )

    return test_results


def apply_rf_model_sklearn(
    ht: hl.Table,
    rf_model: Dict[str, Any],
    features: List[str],
    label: Optional[str] = None,
    probability_col_name: str = "rf_probability",
    prediction_col_name: str = "rf_prediction",
) -> hl.Table:
    """
    Apply a scikit-learn RF model to a Table and annotate the RF probabilities and predictions.

    Mirrors ``apply_rf_model`` from gnomad_methods without Spark, so it runs on the
    Hail Batch backend. The trees are passed to Hail as a literal and evaluated in a
    Hail expression, so the Table is never collected.

    As in ``apply_rf_model``, rows with a missing or NaN feature get missing
    probabilities and predictions, and string feature values not seen in training
    get the next index, like StringIndexer's "keep" option.

    :param ht: Input Table.
    :param rf_model: scikit-learn RF model: a dict with the fitted forest (``rf``),
        ``features``, ``labels``, and the StringIndexer labels of each string feature
        (``string_feature_labels``).
    :param features: List of feature columns in the model. Should match the model
        list of features.
    :param label: Optional column containing labels. Only used to check its type, as
        in ``apply_rf_model``. Default is None.
    :param probability_col_name: Name of the column that will store the RF
        probabilities. Default is "rf_probability".
    :param prediction_col_name: Name of the column that will store the RF predictions.
        Default is "rf_prediction".
    :return: Table with RF columns.
    """
    logger.info("Applying RF model.")

    check_fields = features[:]
    if label:
        check_fields.append(label)
    check_ht_fields_for_spark(ht, check_fields)

    # Flatten all trees into one set of node arrays. Leaves keep -1 as their
    # children; other children are shifted to forest-wide node indexes.
    roots, feature, threshold, left, right, value = [], [], [], [], [], []
    for estimator in rf_model["rf"].estimators_:
        tree = estimator.tree_
        offset = len(feature)
        roots.append(offset)
        feature.extend(tree.feature.tolist())
        threshold.extend(tree.threshold.tolist())
        left.extend([c + offset if c >= 0 else -1 for c in tree.children_left.tolist()])
        right.extend(
            [c + offset if c >= 0 else -1 for c in tree.children_right.tolist()]
        )
        # Each tree votes with the class proportions of its leaf.
        node_value = tree.value[:, 0, :]
        value.extend((node_value / node_value.sum(axis=1, keepdims=True)).tolist())

    forest = hl.literal(
        hl.Struct(
            roots=roots,
            feature=feature,
            threshold=threshold,
            left=left,
            right=right,
            value=value,
        ),
        dtype=hl.tstruct(
            roots=hl.tarray(hl.tint32),
            feature=hl.tarray(hl.tint32),
            threshold=hl.tarray(hl.tfloat64),
            left=hl.tarray(hl.tint32),
            right=hl.tarray(hl.tint32),
            value=hl.tarray(hl.tarray(hl.tfloat64)),
        ),
    )
    depth = max(estimator.tree_.max_depth for estimator in rf_model["rf"].estimators_)

    # Build the feature vector in model order, indexing string features like
    # StringIndexer. scikit-learn trains on float32 values, so features are cast to
    # float32 for the split thresholds to match.
    string_feature_labels = rf_model["string_feature_labels"]
    x = []
    for f in rf_model["features"]:
        expr = ht[f]
        if f in string_feature_labels:
            levels = string_feature_labels[f]
            expr = hl.literal({v: i for i, v in enumerate(levels)}).get(
                expr, len(levels)
            )
        x.append(hl.float32(expr))
    x = hl.array(x)

    # Rows with a missing or NaN feature are not scored, like in ``ht_to_rf_df``.
    defined = _get_rf_features_defined_expr(ht, features)

    # Walk each tree from its root to a leaf. A leaf has no children and stays put.
    leaves = forest.roots.map(
        lambda root: hl.fold(
            lambda node, _: hl.if_else(
                forest.left[node] < 0,
                node,
                hl.if_else(
                    x[forest.feature[node]] <= forest.threshold[node],
                    forest.left[node],
                    forest.right[node],
                ),
            ),
            root,
            hl.range(depth),
        )
    )
    # Average the leaf class proportions over trees, like ``predict_proba``.
    labels = rf_model["labels"]
    probability = hl.bind(
        lambda tree_leaves: hl.range(len(labels)).map(
            lambda i: hl.sum(tree_leaves.map(lambda leaf: forest.value[leaf][i]))
            / len(roots)
        ),
        leaves,
    )

    ht = ht.annotate(**{probability_col_name: hl.or_missing(defined, probability)})
    ht = ht.annotate(
        **{
            probability_col_name: {
                c: ht[probability_col_name][i] for i, c in enumerate(labels)
            },
            prediction_col_name: hl.literal(labels)[
                hl.argmax(ht[probability_col_name])
            ],
        }
    )

    return ht


def save_model_sklearn(
    rf_model: Dict[str, Any], out_path: str, overwrite: bool = False
) -> None:
    """
    Save a scikit-learn RF model as a pickle.

    Load it with the same scikit-learn version it was trained with.

    :param rf_model: scikit-learn RF model to save.
    :param out_path: Output path.
    :param overwrite: If set, will overwrite an existing file at the output location.
        Default is False.
    :return: None.
    """
    logger.info("Saving model to %s", out_path)
    if not overwrite and hl.hadoop_exists(out_path):
        raise FileExistsError(f"{out_path} already exists and overwrite is not set.")
    with hl.hadoop_open(out_path, "wb") as f:
        pickle.dump(rf_model, f)


def load_model_sklearn(input_path: str) -> Dict[str, Any]:
    """
    Load a scikit-learn RF model saved by ``save_model_sklearn``.

    :param input_path: Location of model to load.
    :return: scikit-learn RF model.
    """
    logger.info("Loading model from %s", input_path)
    with hl.hadoop_open(input_path, "rb") as f:
        return pickle.load(f)


def _get_string_indexer_labels(values: pd.Series) -> List[str]:
    """
    Get the distinct values of a column in StringIndexer order.

    StringIndexer orders values by frequency, most frequent first, and breaks ties
    alphabetically.

    :param values: Column to index.
    :return: Distinct values in index order.
    """
    counts = values.value_counts()
    return sorted(counts.index, key=lambda v: (-counts[v], v))


def train_rf_sklearn(
    ht: hl.Table,
    features: List[str],
    label: str,
    num_trees: int = 500,
    max_depth: int = 5,
) -> Dict[str, Any]:
    """
    Train a scikit-learn RF model.

    Mirrors ``train_rf`` from gnomad_methods without Spark, so it runs on the Hail
    Batch backend. As in ``train_rf``, rows with a missing or NaN feature are dropped,
    missing labels are replaced with "NA", labels and string features are indexed
    like StringIndexer, and each split considers the square root of the number of
    features, rounded up.

    .. note::

        The training rows are collected into memory with ``to_pandas``. Unlike Spark,
        scikit-learn trains on float32 values, finds exact splits instead of binning
        each feature into 32 bins, and treats indexed string features as ordered
        rather than categorical.

    :param ht: Input Table.
    :param features: List of columns to be used as features.
    :param label: Column containing the label to predict.
    :param num_trees: Number of trees to use. Default is 500.
    :param max_depth: Maximum tree depth. Default is 5.
    :return: scikit-learn RF model: a dict with the fitted forest (``rf``),
        ``features``, ``labels``, and the StringIndexer labels of each string feature
        (``string_feature_labels``).
    """
    logger.info(
        "Training RF model using:\nfeatures: %s\nlabels: %s\nnum_trees:"
        " %d\nmax_depth: %d",
        ",".join(features),
        label,
        num_trees,
        max_depth,
    )

    check_ht_fields_for_spark(ht, features + [label])

    # Like ``ht_to_rf_df``, drop rows with a missing or NaN feature (in Hail, so they
    # are not collected) and replace missing labels with "NA".
    ht = ht.key_by().select(*features, label)
    df = ht.filter(_get_rf_features_defined_expr(ht, features)).to_pandas()
    df[label] = df[label].fillna("NA")

    labels = _get_string_indexer_labels(df[label])
    logger.info("Found labels: %s", labels)

    string_features = [f for f in features if ht[f].dtype == hl.tstr]
    if string_features:
        logger.info("Indexing string features: %s", ",".join(string_features))
    string_feature_labels = {
        f: _get_string_indexer_labels(df[f]) for f in string_features
    }

    # Assemble the features in order, with string features replaced by their index.
    x = np.column_stack(
        [
            (
                df[f].map({v: i for i, v in enumerate(string_feature_labels[f])})
                if f in string_feature_labels
                else df[f]
            ).to_numpy(dtype=np.float32)
            for f in features
        ]
    )
    y = df[label].map({v: i for i, v in enumerate(labels)}).to_numpy(dtype=np.int64)

    rf = RandomForestClassifier(
        n_estimators=num_trees,
        max_depth=max_depth,
        max_features=math.ceil(math.sqrt(len(features))),
        n_jobs=-1,
    )

    # Train model.
    logger.info("Training RF model")
    rf.fit(x, y)

    feature_importance = dict(zip(features, rf.feature_importances_.tolist()))
    logger.info(
        "RF features importance:\n%s",
        "\n".join([f"{f}: {i}" for f, i in feature_importance.items()]),
    )

    return {
        "rf": rf,
        "features": features,
        "labels": labels,
        "string_feature_labels": string_feature_labels,
    }


def train_rf_model_sklearn(
    ht: hl.Table,
    rf_features: List[str],
    tp_expr: hl.expr.BooleanExpression,
    fp_expr: hl.expr.BooleanExpression,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    test_expr: hl.expr.BooleanExpression = False,
) -> Tuple[hl.Table, Dict[str, Any]]:
    """
    Perform scikit-learn random forest (RF) training using a Table annotated with features and training data.

    Mirrors ``train_rf_model`` from gnomad_methods, but trains and tests the model
    with ``train_rf_sklearn`` and ``test_model_sklearn``, so it runs on the Hail Batch
    backend. The returned Table has the same annotations as in ``train_rf_model``.

    :param ht: Table annotated with features for the RF model and the positive and
        negative training data.
    :param rf_features: List of column names to use as features in the RF training.
    :param tp_expr: TP training expression.
    :param fp_expr: FP training expression.
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all
        training examples are used. Default is 1.0.
    :param num_trees: Number of trees in the RF model. Default is 500.
    :param max_depth: Maximum tree depth in the RF model. Default is 5.
    :param test_expr: An expression specifying variants to hold out for testing and
        use for evaluation only. Default is False.
    :return: Table with TP and FP training sets used in the RF training and the
        resulting scikit-learn RF model.
    """
    ht = ht.annotate(_tp=tp_expr, _fp=fp_expr, rf_test=test_expr)

    rf_ht = sample_training_examples(
        ht, tp_expr=ht._tp, fp_expr=ht._fp, fp_to_tp=fp_to_tp, test_expr=ht.rf_test
    )
    ht = ht.annotate(rf_train=rf_ht[ht.key].train, rf_label=rf_ht[ht.key].label)

    summary = ht.group_by("_tp", "_fp", "rf_train", "rf_label", "rf_test").aggregate(
        n=hl.agg.count()
    )
    logger.info("Summary of TP/FP and RF training data:")
    summary.show(n=20)

    logger.info(
        "Training RF model:\nfeatures: %s\nnum_tree: %d\nmax_depth:%d",
        ",".join(rf_features),
        num_trees,
        max_depth,
    )

    rf_model = train_rf_sklearn(
        ht.filter(ht.rf_train),
        features=rf_features,
        label="rf_label",
        num_trees=num_trees,
        max_depth=max_depth,
    )

    test_results = None
    if test_expr is not None:
        logger.info("Testing model on specified variants or intervals...")
        test_ht = ht.filter(hl.is_defined(ht.rf_label) & ht.rf_test)
        test_results = test_model_sklearn(
            test_ht, rf_model, features=rf_features, label="rf_label"
        )

    features_importance = dict(
        zip(rf_model["features"], rf_model["rf"].feature_importances_.tolist())
    )
    ht = ht.select_globals(
        features_importance=features_importance,
        features=rf_features,
        test_results=test_results,
    )

    return ht.select("rf_train", "rf_label", "rf_test"), rf_model


def train_rf(
    ht: hl.Table,
    test: bool = False,
    features: List[str] = None,
    fp_to_tp: float = 1.0,
    num_trees: int = 500,
    max_depth: int = 5,
    transmitted_singletons: bool = False,
    sibling_singletons: bool = False,
    adj: bool = False,
    filter_centromere_telomere: bool = False,
    test_intervals: Union[str, List[str]] = "chr20",
    environment: str = "rwb",
) -> Tuple[hl.Table, Any]:
    """
    Train random forest model using `train_rf_model`, or `train_rf_model_sklearn` on Batch.

    :param ht: Table containing annotations needed for RF training.
    :param test: Whether to filter the input Table to `DENSE_TRIO_TEST_CHROMS` and
        `test_intervals` for test purposes. Default is False.
    :param features: List of features to use in the random forests model. When no
        `features` list is provided, the global FEATURES is used.
    :param fp_to_tp: Ratio of FPs to TPs for creating the RF model. If set to 0, all
        training examples are used. Default is 1.0.
    :param num_trees: Number of trees in the RF model. Default is 500.
    :param max_depth: Maximum tree depth in the RF model. Default is 5.
    :param transmitted_singletons: Whether to use transmitted singletons for training.
        Default is False.
    :param sibling_singletons: Whether to use sibling singletons for training. Default
        is False.
    :param adj: Whether to use adj genotypes for transmitted/sibling singletons instead
        of raw. Default is False and raw is used.
    :param filter_centromere_telomere: Filter centromeres and telomeres before training.
        Default is False.
    :param test_intervals: Specified interval(s) will be held out for testing and
        evaluation only. Default is "chr20".
    :param environment: Environment to run in. On "batch", which has no Spark, the RF
        is trained with scikit-learn using `train_rf_model_sklearn`. Default is "rwb".
    :return: Input `ht` annotated with training information and the RF model.
    """
    if features is None:
        features = FEATURES

    if test_intervals:
        if isinstance(test_intervals, str):
            test_intervals = [test_intervals]
        test_intervals = [
            hl.parse_locus_interval(x, reference_genome="GRCh38")
            for x in test_intervals
        ]

    if test:
        # Test trio stats (singletons) exist only for DENSE_TRIO_TEST_CHROMS.
        logger.info(
            "Filtering to %s and evaluation intervals for testing...",
            ", ".join(DENSE_TRIO_TEST_CHROMS),
        )
        test_chrom_intervals = [
            hl.parse_locus_interval(c, reference_genome="GRCh38")
            for c in DENSE_TRIO_TEST_CHROMS
        ]
        ht = hl.filter_intervals(ht, test_chrom_intervals + test_intervals)

    logger.info("Annotating true positives and false positives in HT...")
    fp_expr = ht.fail_hard_filters
    tp_expr = ht.omni | ht.mills
    transmit_expr = ht.transmitted_singleton_raw
    sibling_expr = ht.sibling_singleton_raw

    if adj:
        transmit_expr = ht.transmitted_singleton_adj
        sibling_expr = ht.sibling_singleton_adj

    if transmitted_singletons:
        tp_expr |= transmit_expr
    if sibling_singletons:
        tp_expr |= sibling_expr

    ht = ht.annotate(tp=tp_expr, fp=fp_expr)

    rf_ht = ht
    if filter_centromere_telomere:
        logger.info("Filtering centromeres and telomeres from HT...")
        rf_ht = rf_ht.filter(
            ~hl.is_defined(telomeres_and_centromeres.ht()[rf_ht.locus])
        )

    # The Hail Batch backend has no Spark, so the RF is trained with scikit-learn.
    train_fn = train_rf_model_sklearn if environment == "batch" else train_rf_model
    rf_ht, rf_model = train_fn(
        rf_ht,
        rf_features=features,
        tp_expr=rf_ht.tp,
        fp_expr=rf_ht.fp,
        fp_to_tp=fp_to_tp,
        num_trees=num_trees,
        max_depth=max_depth,
        test_expr=hl.literal(test_intervals).any(
            lambda interval: interval.contains(rf_ht.locus)
        ),
    )

    rf_ht = rf_ht.annotate_globals(
        fp_to_tp=fp_to_tp,
        num_trees=num_trees,
        max_depth=max_depth,
        transmitted_singletons=transmitted_singletons,
        sibling_singletons=sibling_singletons,
        adj=adj,
        filter_centromere_telomere=filter_centromere_telomere,
        test_intervals=test_intervals,
    )

    ht = ht.select("tp", "fp").join(rf_ht, how="left")

    return ht, rf_model


def add_model_to_run_list(
    ht: hl.Table,
    model_id: str,
    rf_runs: Dict[str, Any],
    rf_run_path: str,
) -> None:
    """
    Add RF model run to RF run list.

    :param ht: Table containing RF model run information as globals.
    :param model_id: ID of RF model run.
    :param rf_runs: Dictionary containing current RF run information.
    :param rf_run_path: Path to RF run list.
    :return: None
    """
    ht = ht.annotate_globals(test_intervals=hl.str(ht.test_intervals))
    ht_globals = hl.eval(ht.globals)
    input_args = [
        "transmitted_singletons",
        "sibling_singletons",
        "adj",
        "filter_centromere_telomere",
    ]
    rf_output = ["test_intervals", "features_importance", "test_results"]
    rf_runs[model_id] = get_run_data(
        input_args={k: ht_globals[k] for k in input_args},
        **{k: ht_globals[k] for k in rf_output},
    )

    with hl.hadoop_open(rf_run_path, "w") as f:
        json.dump(rf_runs, f)


def get_or_create_model_id(
    model_id: str = None,
    test: bool = False,
    environment: str = "rwb",
) -> Tuple[str, str]:
    """
    Get or create a unique model ID for the RF run.

    :param model_id: Existing model ID to use. If None, generates a new unique ID.
    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "rwb".
    :return: Tuple of (model_id, rf_run_path).
    """
    rf_run_path = get_rf_run_path(test=test, environment=environment)
    rf_runs = get_rf_runs(rf_run_path)

    if model_id is None:
        model_id = f"rf_{str(uuid.uuid4())[:8]}"
        while model_id in rf_runs:
            model_id = f"rf_{str(uuid.uuid4())[:8]}"

    return model_id, rf_run_path


def main(args):
    """Run random forest variant QC pipeline."""
    overwrite = args.overwrite
    test = args.test
    environment = args.environment
    if args.model_id is not None and _validate_model_id(args.model_id) != "rf":
        raise ValueError(f"model_id must start with 'rf_', but got {args.model_id}")

    _init_hail(
        "random_forest",
        environment=environment,
        billing_project=args.gcp_billing_project,
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{qc_temp_prefix(environment=environment, days=args.tmp_dir_days)}random_forest",
        **_get_batch_resource_kwargs(args),
    )

    model_id, rf_run_path = get_or_create_model_id(
        model_id=args.model_id,
        test=test,
        environment=environment,
    )

    # Get resource paths.
    vqc_annotations_ht_resource = get_variant_qc_annotations(
        test=test, environment=environment
    )
    rf_training_ht_resource = get_rf_training(
        model_id=model_id, test=test, environment=environment
    )
    rf_model_path = get_rf_model_path(
        model_id=model_id, test=test, environment=environment
    )
    rf_result_ht_resource = get_variant_qc_result(
        model_id=model_id, test=test, environment=environment
    )

    try:
        if args.list_rf_runs:
            logger.info("RF runs:")
            pretty_print_runs(get_rf_runs(rf_run_path))

        if args.train_rf:
            # Check input exists.
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "variant_qc_annotations_ht": [vqc_annotations_ht_resource]
                },
                output_step_resources={
                    "rf_training_ht": [rf_training_ht_resource],
                    "rf_model": [rf_model_path],
                },
                overwrite=overwrite,
            )

            logger.info("Loading variant QC annotations...")
            vqc_annotation_ht = load_variant_qc_annotations_ht(
                test=test, environment=environment
            )

            logger.info("Training RF model...")
            ht, rf_model = train_rf(
                vqc_annotation_ht,
                test=test,
                features=args.features,
                fp_to_tp=args.fp_to_tp,
                num_trees=args.num_trees,
                max_depth=args.max_depth,
                transmitted_singletons=args.transmitted_singletons,
                sibling_singletons=args.sibling_singletons,
                adj=args.adj,
                filter_centromere_telomere=args.filter_centromere_telomere,
                test_intervals=args.test_intervals,
                environment=environment,
            )
            ht = ht.checkpoint(rf_training_ht_resource.path, overwrite=overwrite)

            logger.info("Adding run to RF run list")
            add_model_to_run_list(ht, model_id, get_rf_runs(rf_run_path), rf_run_path)

            logger.info("Saving RF model")
            if environment == "batch":
                save_model_sklearn(rf_model, rf_model_path, overwrite=overwrite)
            else:
                save_model(rf_model, rf_model_path, overwrite=overwrite)

        if args.apply_rf:
            # Check inputs exist.
            _check_resource_existence(
                environment=environment,
                input_step_resources={
                    "rf_training_ht": [rf_training_ht_resource],
                    "rf_model": [rf_model_path],
                },
                output_step_resources={"rf_result_ht": [rf_result_ht_resource]},
                overwrite=overwrite,
            )

            logger.info("Applying RF model %s...", model_id)
            rf_ht = rf_training_ht_resource.ht()
            rf_features = hl.eval(rf_ht.features)

            logger.info("Loading variant QC annotations...")
            vqc_annotation_ht = load_variant_qc_annotations_ht(
                test=test, environment=environment
            )
            ht = rf_ht.annotate(**vqc_annotation_ht[rf_ht.key].select(*rf_features))
            if environment == "batch":
                ht = apply_rf_model_sklearn(
                    ht,
                    rf_model=load_model_sklearn(rf_model_path),
                    features=rf_features,
                )
            else:
                ht = apply_rf_model(
                    ht,
                    rf_model=load_model(rf_model_path),
                    features=rf_features,
                )

            logger.info("Finished applying RF model...")
            summary_cols = [
                "tp",
                "fp",
                RF_COLS["train"],
                RF_COLS["label"],
                RF_COLS["prediction"],
            ]
            ht = ht.select(*summary_cols, RF_COLS["probability"])
            ht = ht.annotate_globals(model_id=model_id)
            ht = ht.checkpoint(rf_result_ht_resource.path, overwrite=overwrite)
            ht.group_by(*summary_cols).aggregate(n=hl.agg.count()).show(-1)

    finally:
        if environment == "rwb":
            logger.info("Copying log to logging bucket...")
            hl.copy_log(get_logging_path("random_forest", environment=environment))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            f"If the dataset should be filtered to {', '.join(DENSE_TRIO_TEST_CHROMS)}"
            " for testing (also filtered to evaluation interval specified by"
            " --test-intervals)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--model-id",
        help=(
            "Model ID. Created by --train-rf and only needed for --apply-rf without"
            " running --train-rf."
        ),
        required=False,
    )
    parser.add_argument(
        "--environment",
        help="Environment to run in.",
        choices=["rwb", "batch", "dataproc"],
        default="rwb",
    )
    parser.add_argument(
        "--tmp-dir-days",
        help="Number of days to keep temporary files.",
        default=4,
        type=int,
    )
    batch_args = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
    )

    batch_args.add_argument(
        "--gcp-billing-project",
        help="GCP billing project to use for requester pays buckets (required for --environment=batch).",
        default=None,
    )
    batch_args.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_args.add_argument(
        "--driver-cores",
        type=int,
        default=None,
        help="Number of cores for driver node.",
    )
    batch_args.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_args.add_argument(
        "--worker-cores",
        type=int,
        default=None,
        help="Number of cores for worker nodes.",
    )
    batch_args.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
    )

    actions = parser.add_argument_group("Actions")
    actions.add_argument(
        "--list-rf-runs",
        help=(
            "Lists all previous RF runs, along with their model ID, parameters and"
            " testing results."
        ),
        action="store_true",
    )
    actions.add_argument("--train-rf", help="Trains RF model.", action="store_true")
    actions.add_argument(
        "--apply-rf", help="Applies RF model to the data.", action="store_true"
    )

    rf_params = parser.add_argument_group("Random Forest parameters")
    rf_params.add_argument(
        "--features",
        help="Features to use in the random forests model.",
        default=FEATURES,
        type=str,
        nargs="+",
    )
    rf_params.add_argument(
        "--fp-to-tp",
        help=(
            "Ratio of FPs to TPs for training the RF model. If 0, all training examples"
            " are used. Default is 1.0."
        ),
        default=1.0,
        type=float,
    )
    rf_params.add_argument(
        "--test-intervals",
        help=(
            "The specified interval(s) will be held out for testing and evaluation"
            ' only. Default is "chr20".'
        ),
        nargs="+",
        type=str,
        default="chr20",
    )
    rf_params.add_argument(
        "--num-trees",
        help="Number of trees in the RF model. Default is 500.",
        default=500,
        type=int,
    )
    rf_params.add_argument(
        "--max-depth",
        help="Maximum tree depth in the RF model. Default is 5.",
        default=5,
        type=int,
    )

    training_params = parser.add_argument_group("Training data parameters")
    training_params.add_argument(
        "--adj",
        help="Use adj genotypes for transmitted/sibling singletons.",
        action="store_true",
    )
    training_params.add_argument(
        "--transmitted-singletons",
        help="Include transmitted singletons in training.",
        action="store_true",
    )
    training_params.add_argument(
        "--sibling-singletons",
        help="Include sibling singletons in training.",
        action="store_true",
    )
    training_params.add_argument(
        "--filter-centromere-telomere",
        help="Train RF without centromeres and telomeres.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    batch_args = [
        "app_name",
        "driver_cores",
        "driver_memory",
        "worker_cores",
        "worker_memory",
    ]
    provided_batch_args = [arg for arg in batch_args if getattr(args, arg) is not None]
    if provided_batch_args and args.environment != "batch":
        parser.error(
            f"Batch configuration arguments ({', '.join('--' + a.replace('_', '-') for a in provided_batch_args)}) "
            f"require --environment=batch"
        )

    if args.environment == "batch" and args.gcp_billing_project is None:
        parser.error("--gcp-billing-project is required when --environment=batch")

    if not args.model_id and not args.train_rf and args.apply_rf:
        parser.error(
            "Error: --model-id is required when running --apply-rf without running"
            " --train-rf too."
        )

    if args.model_id and args.train_rf:
        parser.error(
            "Error: --model-id and --train-rf are mutually exclusive. --train-rf will"
            " generate a run model ID."
        )

    main(args)
