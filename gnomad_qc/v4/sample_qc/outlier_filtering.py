"""Add script docstring."""
import argparse
import logging
import math
from typing import List, Optional

import hail as hl
import pandas as pd
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
    platform,
    regressed_metrics,
    stratified_metrics,
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
    :rtype: hl.Table
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
        qc_metrics={
            metric: sample_qc_ht.sample_qc[metric] for metric in filtering_qc_metrics
        },
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
    pop_n_pcs: int = 16,
    platform_n_pcs: int = 9,
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
    :rtype: hl.Table
    """
    if pop_scores_ht is None and platform_ht is None:
        ValueError(
            "At least one of 'pop_scores_ht' or 'platform_ht' must be specified!"
        )

    ann_args = {}
    if pop_scores_ht is not None:
        ann_args["pop_scores"] = pop_scores_ht[ht.key].scores
    if platform_ht is not None:
        ann_args["platform_scores"] = platform_ht[ht.key].scores
    if not include_unreleasable_samples:
        ann_args["releasable"] = project_meta_ht[ht.key].project_meta.releasable

    ht = ht.annotate(**ann_args)

    ht = compute_qc_metrics_residuals(
        ht,
        pc_scores=ht.pop_scores[:pop_n_pcs] + ht.platform_scores[:platform_n_pcs],
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
        pop_n_pcs=pop_n_pcs,
        platform_n_pcs=platform_n_pcs,
    )

    return ht


def determine_nearest_neighbors(
    ht,
    platform_ht,
    project_meta_ht,
    n_pcs: int = 16,
    n_neighbors: int = 50,
    n_jobs: int = -2,
    include_unreleasable_samples: bool = False,
):
    """
    Add module docstring.

    :param ht:
    :param platform_ht:
    :param project_meta_ht:
    :param n_pcs:
    :param n_neighbors:
    :param n_jobs:
    :param include_unreleasable_samples:
    :return:
    """
    # TODO: Filter HT to what we need before pd export
    ht = ht.annotate(
        pop_scores=ancestry_pca_scores(include_unreleasable_samples)
        .ht()[ht.key]
        .scores,
        platform=platform_ht[ht.key].platform,
        releasable=project_meta_ht[ht.key].project_meta.releasable,
    )
    ht = ht.filter(hl.is_defined(ht.pop_scores))
    pop_scores_pd = ht.to_pandas()
    pop_scores_pd_s = pop_scores_pd.s
    pop_scores_pd = pop_scores_pd[[f"PC{i}" for i in range(1, n_pcs + 1)]]

    vals = pop_scores_pd.values
    nbrs = NearestNeighbors(n_neighbors=n_neighbors, n_jobs=n_jobs)
    nbrs.fit(vals)
    distances, indexes = nbrs.kneighbors(vals)

    # Format neighbor distances as a Hail Table
    distances_pd = pd.DataFrame(distances)
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

    # Join neighbor distances Table and indexes Table
    nbrs_ht = distances_ht.join(indexes_ht)
    nbrs_ht = nbrs_ht.add_index()

    # Add nearest_neighbors annotation
    explode_nbrs_ht = nbrs_ht.key_by("idx").explode("nearest_neighbor_idxs")
    explode_nbrs_ht = explode_nbrs_ht.annotate(
        nbr=explode_nbrs_ht[hl.int64(explode_nbrs_ht.nbrs_indexes)].s
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
    nn_ht: hl.Table,
    sample_qc_ht: hl.Table,
    filtering_qc_metrics: List[str],
):
    """
    Add module docstring.

    :param nn_ht:
    :param sample_qc_ht:
    :param filtering_qc_metrics:
    :return:
    """
    explode_nbrs_ht = nn_ht.explode("nearest_neighbors", name="nbr")
    explode_nbrs_ht = explode_nbrs_ht.annotate(**sample_qc_ht[explode_nbrs_ht.nbr])

    agg_expr = hl.struct(
        **{
            metric: get_median_and_mad_expr(explode_nbrs_ht[metric])
            for metric in filtering_qc_metrics
        }
    )

    return explode_nbrs_ht.group_by("s").aggregate(**agg_expr)


def main(args):
    """Add main docstring."""
    hl.init(
        log="/outlier_filtering.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test = args.test
    overwrite = args.overwrite
    filtering_qc_metrics = args.filtering_qc_metrics.split(",")
    population_correction = args.population_correction
    platform_correction = args.platform_correction
    pop_ht = get_pop_ht(test=test).ht()
    pop_scores_ht = ancestry_pca_scores(
        include_unreleasable_samples=args.include_unreleasable_samples,
        test=test,
    ).ht()
    platform_ht = platform.ht()  # get_platform_ht(test=test).ht()

    sample_qc_ht = get_sample_qc("bi_allelic", test=test).ht()

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

    if args.apply_regressed_filters:
        if not population_correction:
            regress_pop_n_pcs = None
        if not platform_correction:
            regress_platform_n_pcs = None
        apply_regressed_filters(
            sample_qc_ht,
            pop_scores_ht,
            filtering_qc_metrics,
            regress_pop_n_pcs=regress_pop_n_pcs,
            regress_platform_n_pcs=regress_platform_n_pcs,
        ).write(regressed_metrics.path, overwrite=overwrite)

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
        ).write(stratified_metrics.path, overwrite=overwrite)

    # if args.apply_nearest_neighbor_filters:


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
        "--filtering_qc_metrics",
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
        "--apply_stratified_filters",
        help="Compute per pop filtering.",
        action="store_true",
    )
    parser.add_argument(
        "--apply_regressed_filters",
        help="Computes qc_metrics adjusted for pop.",
        action="store_true",
    )
    parser.add_argument(
        "--regress_n_pcs",
        help="Number of PCs to use for qc metric regressions",
        default=10,
        type=int,
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
