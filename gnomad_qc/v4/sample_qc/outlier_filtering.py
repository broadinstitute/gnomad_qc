"""Add script docstring."""
import argparse
import logging
import math
from typing import Dict, List, Optional, Tuple

import hail as hl
import pandas as pd
from annoy import AnnoyIndex
from gnomad.sample_qc.filtering import (
    compute_qc_metrics_residuals,
    compute_stratified_metrics_filter,
)
from gnomad.utils.gen_stats import get_median_and_mad_expr
from gnomad.utils.slack import slack_notifications
from sklearn.neighbors import NearestNeighbors

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import check_resource_existence
from gnomad_qc.v4.resources.sample_qc import (
    ancestry_pca_scores,
    get_pop_ht,
    get_sample_qc,
    hard_filtered_samples,
    joint_qc_meta,
    nearest_neighbors,
    nearest_neighbors_filtering,
    platform,
    regressed_filtering,
    stratified_filtering,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("outlier_filtering")
logger.setLevel(logging.INFO)


def apply_stratified_filters(
    sample_qc_ht: hl.Table,
    filtering_qc_metrics: List[str],
    pop_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Use population stratified QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to run hl.sample_qc and return the metrics stratified by assigned
    population with a `qc_metrics_filters` annotation indicating if the sample falls a certain number of MAD
    outside the distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4, we modify this for n_singleton to be
        (4.0, 8.0)

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    """
    logger.info(
        "Computing stratified QC metrics filters using metrics: "
        + ", ".join(filtering_qc_metrics)
    )

    strata = {}
    if pop_ht is not None:
        strata["pop"] = pop_ht[sample_qc_ht.key].qc_pop
    if platform_ht is not None:
        strata["platform"] = platform_ht[sample_qc_ht.key].qc_platform

    stratified_metrics_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        strata=strata,
        metric_threshold={"n_singleton": (4.0, 8.0)},
    )
    return stratified_metrics_ht


def apply_regressed_filters(
    ht: hl.Table,
    filtering_qc_metrics: List[str],
    pop_scores_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
    project_meta_ht: Optional[hl.Table] = None,
    regress_pop_n_pcs: int = 16,
    regress_platform_n_pcs: int = 9,
    include_unreleasable_samples: bool = False,
) -> hl.Table:
    """
    Compute sample QC metrics residuals after regressing out population PCs and determine what samples are outliers that should be filtered.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4. We
        modify this for n_singleton_residual
        to be (math.inf, 8.0) and for r_het_hom_var_residual to be (math.inf,
        4.0). The math.inf is used to prevent a
        lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :param include_unreleasable_samples: Should unreleasable samples be included in
    the regression
    :param n_pcs: Number of population PCs to use in regression
    :return: Table with stratified metrics and filters
    """
    # TODO: add option to only regress out population, and do this regression
    # per platform
    if pop_scores_ht is None and platform_ht is None:
        ValueError(
            "At least one of 'pop_scores_ht' or 'platform_ht' must be specified!"
        )

    pc_scores = None
    ann_args = {}
    if pop_scores_ht is not None:
        ann_args["pop_scores"] = pc_scores = pop_scores_ht[ht.key].scores[
            :regress_pop_n_pcs
        ]
    if platform_ht is not None:
        ann_args["platform_scores"] = platform_ht[ht.key].scores[
            :regress_platform_n_pcs
        ]
        if pc_scores is None:
            pc_scores = ann_args["platform_scores"]
        else:
            pc_scores = pc_scores.extend(ann_args["platform_scores"])
    if not include_unreleasable_samples:
        ann_args["releasable"] = project_meta_ht[ht.key].releasable

    ann_args["pc_scores"] = pc_scores
    ht = ht.annotate(**ann_args)
    ht = compute_qc_metrics_residuals(
        ht,
        pc_scores=ht.pc_scores,
        qc_metrics={metric: ht[metric] for metric in filtering_qc_metrics},
        regression_sample_inclusion_expr=None
        if include_unreleasable_samples
        else ht.releasable,
    )
    stratified_metrics_ht = compute_stratified_metrics_filter(
        ht,
        qc_metrics=dict(ht.row_value),
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
        },
    )

    ht = ht.annotate(**stratified_metrics_ht[ht.key])
    ht = ht.annotate_globals(
        **stratified_metrics_ht.index_globals(),
        regress_pop_n_pcs=regress_pop_n_pcs,
        regress_platform_n_pcs=regress_platform_n_pcs,
    )

    return ht


def determine_nearest_neighbors(
    ht,
    pop_scores_ht: Optional[hl.Table] = None,
    # platform_ht: Optional[hl.Table] = None,
    n_pcs: int = 16,
    n_neighbors: int = 50,
    n_jobs: int = -2,
    add_neighbor_distances: bool = False,
    distance_metric: str = "euclidean",
    use_approximation: bool = False,
    n_trees: int = 50,
):
    """
    Add module docstring.

    :param ht:
    :param pop_scores_ht:
    :param n_pcs:
    :param n_neighbors:
    :param n_jobs:
    :param add_neighbor_distances:
    :param distance_metric:
    :param use_approximation:
    :param n_trees:
    :return:
    """
    # TODO: add option to stratify nearest neighbors by platform

    ht = ht.annotate(
        pop_scores=pop_scores_ht[ht.key].scores,
        # platform=platform_ht[ht.key].platform,
    )
    ht = ht.filter(hl.is_defined(ht.pop_scores))
    ht = ht.select(**{f"PC{i}": ht.pop_scores[i - 1] for i in range(1, n_pcs + 1)})
    pop_scores_pd = ht.to_pandas()
    pop_scores_pd_s = pop_scores_pd.s
    pop_scores_pd = pop_scores_pd[[f"PC{i}" for i in range(1, n_pcs + 1)]]

    if use_approximation:
        nbrs = AnnoyIndex(n_pcs, distance_metric)
        for i, row in pop_scores_pd.iterrows():
            nbrs.add_item(i, row)
        nbrs.build(n_trees, n_jobs=n_jobs)

        indexes = []
        for i in range(pop_scores_pd.shape[0]):
            indexes.append(
                nbrs.get_nns_by_item(
                    i, n_neighbors, include_distances=add_neighbor_distances
                )
            )

        if add_neighbor_distances:
            distances = [d for i, d in indexes]
            indexes = [i for i, d in indexes]
    else:
        vals = pop_scores_pd.values
        nbrs = NearestNeighbors(
            n_neighbors=n_neighbors, n_jobs=n_jobs, metric=distance_metric
        )
        nbrs.fit(vals)
        indexes = nbrs.kneighbors(vals, return_distance=add_neighbor_distances)
        if add_neighbor_distances:
            distances, indexes = indexes

    # Format neighbor indexes as a Hail Table
    indexes_pd = pd.DataFrame(indexes)
    indexes_pd = pd.concat([pop_scores_pd_s, indexes_pd], axis=1)
    indexes_pd = indexes_pd.rename(
        columns={i: f"nbrs_index_{i}" for i in range(n_neighbors)}
    )
    indexes_ht = hl.Table.from_pandas(indexes_pd, key=["s"])
    indexes_ht = indexes_ht.transmute(
        nearest_neighbor_idxs=hl.array(
            [indexes_ht[f"nbrs_index_{i}"] for i in range(n_neighbors)]
        )
    )

    if add_neighbor_distances:
        # Format neighbor distances as a Hail Table
        distances_pd = pd.DataFrame(distances)
        print(distances_pd)
        distances_pd = distances_pd.rename(
            columns={i: f"nbrs_{i}" for i in range(n_neighbors)}
        )
        distances_pd = pd.concat([pop_scores_pd_s, distances_pd], axis=1)
        distances_ht = hl.Table.from_pandas(distances_pd, key=["s"])
        distances_ht = distances_ht.transmute(
            nearest_neighbor_dists=hl.array(
                [distances_ht[f"nbrs_{str(i)}"] for i in range(n_neighbors)]
            )
        )

        # Join neighbor distances Table and indexes Table
        nbrs_ht = indexes_ht.join(distances_ht)
    else:
        nbrs_ht = indexes_ht

    # Add nearest_neighbors annotation
    nbrs_ht = nbrs_ht.add_index()
    explode_nbrs_ht = nbrs_ht.key_by("idx").explode("nearest_neighbor_idxs")
    explode_nbrs_ht = explode_nbrs_ht.annotate(
        nbr=explode_nbrs_ht[hl.int64(explode_nbrs_ht.nearest_neighbor_idxs)].s
    )
    explode_nbrs_ht = explode_nbrs_ht.group_by("s").aggregate(
        nearest_neighbors=hl.agg.collect_as_set(explode_nbrs_ht.nbr)
    )
    nbrs_ht = nbrs_ht.annotate(
        nearest_neighbors=explode_nbrs_ht[nbrs_ht.key].nearest_neighbors
    )
    nbrs_ht = nbrs_ht.annotate_globals(n_pcs=n_pcs, n_neighbors=n_neighbors)

    return nbrs_ht


def get_nearest_neighbors_median_and_mad(
    ht: hl.Table,
    nn_ht: hl.Table,
    filtering_qc_metrics: List[str],
):
    """
    Add module docstring.

    :param ht:
    :param nn_ht:
    :param filtering_qc_metrics:
    :return:
    """
    nn_ht = nn_ht.explode("nearest_neighbors", name="nbr")
    nn_ht = nn_ht.annotate(**ht[nn_ht.nbr])

    agg_expr = hl.struct(
        **{
            metric: get_median_and_mad_expr(nn_ht[metric])
            for metric in filtering_qc_metrics
        }
    )

    return nn_ht.group_by("s").aggregate(**agg_expr)


def compute_nearest_neighbor_metrics_filter(
    ht: hl.Table,
    nn_ht: hl.Table,
    qc_metrics: Dict[str, hl.expr.NumericExpression],
    lower_threshold: float = 4.0,
    upper_threshold: float = 4.0,
    metric_threshold: Optional[Dict[str, Tuple[float, float]]] = None,
    filter_name: str = "qc_metrics_filters",
) -> hl.Table:
    """
    Compute median, MAD, and upper and lower thresholds for each metric used in outlier filtering.

    :param ht: HT containing relevant sample QC metric annotations
    :param qc_metrics: list of metrics (name and expr) for which to compute the critical values for filtering outliers
    :param lower_threshold: Lower MAD threshold
    :param upper_threshold: Upper MAD threshold
    :param metric_threshold: Can be used to specify different (lower, upper) thresholds for one or more metrics
    :param filter_name: Name of resulting filters annotation
    :return: Table grouped by strata, with upper and lower threshold values computed for each sample QC metric
    """
    _metric_threshold = {
        metric: (lower_threshold, upper_threshold) for metric in qc_metrics
    }
    if metric_threshold is not None:
        _metric_threshold.update(metric_threshold)

    ht = ht.select(**qc_metrics)
    median_mad_ht = get_nearest_neighbors_median_and_mad(ht, nn_ht, qc_metrics)

    ht = ht.annotate(
        qc_metrics_stats=hl.struct(
            **{
                metric: hl.bind(
                    lambda x: x.annotate(
                        lower=x.median - _metric_threshold[metric][0] * x.mad,
                        upper=x.median + _metric_threshold[metric][1] * x.mad,
                    ),
                    median_mad_ht[ht.key][metric],
                )
                for metric in qc_metrics
            }
        )
    )

    ht = ht.transmute(
        **{
            f"fail_{metric}": (ht[metric] <= ht.qc_metrics_stats[metric].lower)
            | (ht[metric] >= ht.qc_metrics_stats[metric].upper)
            for metric in qc_metrics
        }
    )
    ht = ht.annotate(
        **{
            filter_name: hl.set(
                hl.filter(
                    lambda x: hl.is_defined(x),
                    [
                        hl.or_missing(ht[f"fail_{metric}"], metric)
                        for metric in qc_metrics
                    ],
                )
            )
        }
    )
    return ht


def apply_nearest_neighbors_filters(
    sample_qc_ht: hl.Table,
    nn_ht: hl.Table,
    filtering_qc_metrics: List[str],
) -> hl.Table:
    """
    Compute sample QC metrics residuals after regressing out population PCs and determine what samples are outliers that should be filtered.

    :param sample_qc_ht:
    :param nn_ht:
    :param filtering_qc_metrics:
    :return:
    """
    logger.info(
        "Computing stratified QC metrics filters using metrics: "
        + ", ".join(filtering_qc_metrics)
    )

    return compute_nearest_neighbor_metrics_filter(
        sample_qc_ht,
        nn_ht,
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        metric_threshold={"n_singleton": (4.0, 8.0)},
    )


def main(args):
    """Add main docstring."""
    hl.init(
        log="/outlier_filtering.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test = args.test
    overwrite = args.overwrite
    regress_pop_n_pcs = args.regress_pop_n_pcs
    regress_platform_n_pcs = args.regress_platform_n_pcs
    filtering_qc_metrics = args.filtering_qc_metrics.split(",")
    population_correction = args.population_correction
    platform_correction = args.platform_correction
    # pop_ht = get_pop_ht(test=test).ht()
    pop_ht = hl.read_table(
        "gs://gnomad/v4.0/sample_qc/joint/ancestry_inference/gnomad.joint.v4.0.hgdp_tgp_training.pop.rf_w_16pcs_0.75_pop_prob.ht"
    )
    pop_scores_ht = ancestry_pca_scores(
        include_unreleasable_samples=args.include_unreleasable_samples,
    ).ht()
    platform_ht = platform.ht()  # get_platform_ht(test=test).ht()

    sample_qc_ht = get_sample_qc("bi_allelic").ht()

    # Convert tuples to lists so we can find the index of the passed threshold.
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            f"bases_dp_over_{hl.eval(sample_qc_ht.dp_bins[i])}": sample_qc_ht.bases_over_dp_threshold[
                i
            ]
            for i in range(len(sample_qc_ht.dp_bins))
        },
    )
    sample_qc_ht = sample_qc_ht.filter(
        hl.is_missing(hard_filtered_samples.ht()[sample_qc_ht.key])
    )
    sample_qc_ht = sample_qc_ht.annotate(
        r_snp_indel=sample_qc_ht.n_snp
        / (sample_qc_ht.n_insertion + sample_qc_ht.n_deletion)
    )

    if test:
        sample_qc_ht = sample_qc_ht.sample(0.01)

    if args.apply_regressed_filters:
        if not population_correction:
            regress_pop_n_pcs = None
        if not platform_correction:
            regress_platform_n_pcs = None
        apply_regressed_filters(
            sample_qc_ht,
            filtering_qc_metrics,
            pop_scores_ht,
            platform_ht,
            joint_qc_meta.ht(),
            regress_pop_n_pcs=regress_pop_n_pcs,
            regress_platform_n_pcs=regress_platform_n_pcs,
        ).write(regressed_filtering(test=test).path, overwrite=overwrite)

    if args.apply_stratified_filters:
        if not population_correction:
            pop_ht = None
        if not platform_correction:
            platform_ht = None
        apply_stratified_filters(
            sample_qc_ht,
            filtering_qc_metrics,
            pop_ht,
            platform_ht,
        ).write(stratified_filtering(test=test).path, overwrite=overwrite)

    if args.determine_nearest_neighbors:
        determine_nearest_neighbors(
            sample_qc_ht,
            pop_scores_ht,
            # platform_ht,
            n_pcs=args.nearest_neighbors_pop_n_pcs,
            n_neighbors=args.n_nearest_neighbors,
            n_jobs=args.n_jobs,
            add_neighbor_distances=args.get_neighbor_distances,
            distance_metric=args.distance_metric,
            use_approximation=args.use_nearest_neighbors_approximation,
            n_trees=args.n_trees,
        ).write(nearest_neighbors(test=test).path, overwrite=overwrite)

    if args.apply_nearest_neighbor_filters:
        ht = apply_nearest_neighbors_filters(
            sample_qc_ht, nearest_neighbors(test=test).ht(), filtering_qc_metrics
        )
        ht.show()
        ht.describe()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test", help="Use a test MatrixTableResource as input.", action="store_true"
    )
    parser.add_argument(
        "--filtering-qc-metrics",
        help="List of QC metrics for filtering.",
        default=",".join(
            [
                "n_snp",
                "n_singleton",
                "r_ti_tv",
                "r_insertion_deletion",
                "n_insertion",
                "n_deletion",
                "r_het_hom_var",
                "n_transition",
                "n_transversion",
            ]
        ),
    )
    parser.add_argument(
        "--population-correction",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--platform-correction",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--include-unreleasable-samples",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--apply-stratified-filters",
        help="Compute per pop filtering.",
        action="store_true",
    )
    parser.add_argument(
        "--apply-regressed-filters",
        help="Computes qc_metrics adjusted for pop.",
        action="store_true",
    )
    parser.add_argument(
        "--regress-pop-n-pcs",
        help="Number of population PCs to use for qc metric regressions",
        default=16,
        type=int,
    )
    parser.add_argument(
        "--regress-platform-n-pcs",
        help="Number of platform PCs to use for qc metric regressions",
        default=9,
        type=int,
    )
    parser.add_argument(
        "--determine-nearest-neighbors",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--nearest-neighbors-pop-n-pcs",
        help="",
        default=16,
        type=int,
    )
    parser.add_argument(
        "--n-nearest-neighbors",
        help="",
        default=50,
        type=int,
    )
    parser.add_argument(
        "--n-jobs",
        help="",
        default=-1,
        type=int,
    )
    parser.add_argument(
        "--get-neighbor-distances",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--distance-metric",
        help="",
        default="euclidean",
        type=str,
        # choices=[],
    )
    parser.add_argument(
        "--use-nearest-neighbors-approximation",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--n-trees",
        help="",
        default=10,
        type=int,
    )

    parser.add_argument(
        "--apply-nearest-neighbor-filters",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
# Add platforms or platform PCs somehow
# Add interval QC options? Redo sample QC after interval filtering?
