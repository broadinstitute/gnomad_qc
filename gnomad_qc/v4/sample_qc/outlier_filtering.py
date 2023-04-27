"""Script to determine sample QC metric outliers that should be filtered."""
import argparse
import logging
import math
from collections import defaultdict
from typing import Any, Dict, List, Optional, Union

import hail as hl
from gnomad.sample_qc.filtering import (
    compute_qc_metrics_residuals,
    compute_stratified_metrics_filter,
    determine_nearest_neighbors,
)
from gnomad.utils.slack import slack_notifications
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.sample_qc import (
    ancestry_pca_scores,
    finalized_outlier_filtering,
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


def get_sample_qc_ht(
    sample_qc_ht: hl.Table, test: bool = False, seed: int = 24
) -> hl.Table:
    """
    Get sample QC Table with modifications needed for outlier filtering.

    The following modifications are made:
        - Add 'r_snp_indel' metric
        - Convert each element of `bases_over_dp_threshold` to a top level annotation
          with 'bases_dp_over_{dp_threshold}'
        - Exclude hard filtered samples
        - Sample 1% of the dataset if `test` is True

    :param sample_qc_ht: Sample QC Table.
    :param test: Whether to filter the input Table to a random sample of 1% of the
        dataset. Default is False.
    :param seed: Random seed for making test dataset. Default is 24.
    :return: Sample QC Table for outlier filtering.
    """
    if test:
        sample_qc_ht = sample_qc_ht.sample(0.01, seed=seed)

    # Convert each element of `bases_over_dp_threshold` to an annotation with a name
    # that contains the DP threshold.
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            f"bases_dp_over_{hl.eval(sample_qc_ht.dp_bins[i])}": sample_qc_ht.bases_over_dp_threshold[
                i
            ]
            for i in range(len(sample_qc_ht.dp_bins))
        },
    ).drop("bases_over_gq_threshold")
    # Exclude hard filtered samples from the sample QC Table.
    # They should not be included in the metric distribution stats.
    sample_qc_ht = sample_qc_ht.filter(
        hl.is_missing(hard_filtered_samples.ht()[sample_qc_ht.key])
    )
    # Add 'r_snp_indel' annotation the sample QC HT.
    sample_qc_ht = sample_qc_ht.annotate(
        r_snp_indel=sample_qc_ht.n_snp
        / (sample_qc_ht.n_insertion + sample_qc_ht.n_deletion)
    )

    return sample_qc_ht.select_globals()


def apply_filter(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    filtering_method: str,
    apply_r_ti_tv_singleton_filter: bool = False,
    pop_scores_ht: Optional[hl.Table] = None,
    pop_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
    **kwargs: Any,
) -> hl.Table:
    """
    Apply stratified, regressed, or nearest neighbors outlier filtering method.

    :param sample_qc_ht: Sample QC HT.
    :param qc_metrics: Specific metrics to use for outlier detection.
    :param filtering_method: The filtering method to apply to `sample_qc_ht`. One of
        'stratified', 'regressed', or 'nearest_neighbors'.
    :param apply_r_ti_tv_singleton_filter: Whether to apply the filtering method to
        only samples with: Number of singletons (or n_singleton residuals) >
        median(comparison group number of singletons/n_singleton residuals).
    :param pop_scores_ht: Optional Table with population PCA scores.
    :param pop_ht: Optional Table with population assignment.
    :param platform_ht: Optional Table with platform assignment.
    :param kwargs: Additional parameters to pass to the requested filtering method
        function.
    :return: Table with outlier filter annotations.
    """
    # Create dict of expression parameters needed in filtering method.
    ann_exprs = {}
    if pop_ht is not None:
        ann_exprs["pop_expr"] = pop_ht[sample_qc_ht.key].pop
    if platform_ht is not None:
        ann_exprs["platform_expr"] = platform_ht[sample_qc_ht.key].qc_platform
        if filtering_method == "regressed":
            ann_exprs["platform_scores_expr"] = platform_ht[sample_qc_ht.key].scores
    if pop_scores_ht is not None:
        ann_exprs["pop_scores_expr"] = pop_scores_ht[sample_qc_ht.key].scores

    # Run filtering method using defined expressions, and any other passed parameters.
    if filtering_method == "stratified":
        ht = apply_stratified_filtering_method(
            sample_qc_ht, qc_metrics, **ann_exprs, **kwargs
        )
    elif filtering_method == "regressed":
        ht = apply_regressed_filtering_method(
            sample_qc_ht, qc_metrics, **ann_exprs, **kwargs
        )
    elif filtering_method == "nearest_neighbors":
        ht = apply_nearest_neighbor_filtering_method(sample_qc_ht, qc_metrics, **kwargs)
    else:
        raise ValueError(
            "'filtering_method' must be one of: 'stratified', 'regressed', "
            "'nearest_neighbors'!"
        )

    # Apply the n_singleton median filter for the r_ti_tv_singleton filter.
    if apply_r_ti_tv_singleton_filter:
        ht = ht.checkpoint(new_temp_file("outlier_filtering", extension="ht"))
        ann_exprs = {ann.split("_expr")[0]: expr for ann, expr in ann_exprs.items()}
        ht = apply_n_singleton_filter_to_r_ti_tv_singleton(
            ht, sample_qc_ht, filtering_method, ann_exprs, **kwargs
        )
    ht = ht.annotate_globals(
        apply_r_ti_tv_singleton_filter=apply_r_ti_tv_singleton_filter
    )

    return ht


def apply_stratified_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    pop_expr: Optional[hl.expr.StringExpression] = None,
    platform_expr: Optional[hl.expr.StringExpression] = None,
) -> hl.Table:
    """
    Use population stratified QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to compute the median, MAD, and upper and
    lower thresholds for each `qc_metrics` stratified by assigned population and/or
    platform and return a `qc_metrics_filters` annotation indicating if the sample
    falls a certain number of MAD outside the distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4. We
        modify this for 'n_singleton' to be (math.inf, 8.0) and for
        'r_het_hom_var' to be (math.inf, 4.0). The math.inf (infinity) is used to
        prevent a lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT.
    :param qc_metrics: Specific metrics to use for outlier detection.
    :param pop_expr: Optional expression with population assignment.
    :param platform_expr: Optional expression with platform assignment.
    :return: Table with stratified metrics and filters.
    """
    if pop_expr is None and platform_expr is None:
        raise ValueError(
            "At least one of 'pop_expr' or 'platform_expr' must be specified!"
        )
    logger.info(
        "Computing QC metrics outlier filters with stratification, using metrics: %s",
        ", ".join(qc_metrics),
    )
    strata = {}
    if pop_expr is not None:
        strata["pop"] = pop_expr
    if platform_expr is not None:
        strata["platform"] = platform_expr

    logger.info("Outlier filtering is stratified by: %s", ", ".join(strata.keys()))
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        strata=strata,
        metric_threshold={
            "n_singleton": (math.inf, 8.0),
            "r_het_hom_var": (math.inf, 4.0),
        },
    )

    return filter_ht


def apply_regressed_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    pop_scores_expr: Optional[hl.expr.ArrayExpression] = None,
    platform_scores_expr: Optional[hl.expr.ArrayExpression] = None,
    platform_expr: Optional[hl.expr.StringExpression] = None,
    regress_pop_n_pcs: Optional[int] = 30,
    regress_platform_n_pcs: Optional[int] = 9,
    regress_per_platform: bool = False,
) -> hl.Table:
    """
    Compute sample QC metrics residuals after regressing out specified PCs and determine what samples are outliers that should be filtered.

    The following are all possible filtering options:
        - Include `pop_scores_expr`: Regress population PCs only and determine outliers
          for each metric in `qc_metrics`.
        - Include `platform_scores_expr`: Regress platform PCs only and determine
          outliers for each metric in `qc_metrics`.
        - Include `pop_scores_expr` and `platform_scores_expr`: Regress both population
          PCs and platform PCs and determine outliers for each metric in `qc_metrics`.
        - For all of the above filtering options, there is also a `regress_per_platform`
          option, which will stratify the dataset by platform before performing the
          regression and outlier filtering.

    After regression, `compute_stratified_metrics_filter` is used to compute the
    median, MAD, and upper and lower thresholds for each of the `qc_metrics` residuals,
    optionally stratified by assigned platform, and return a `qc_metrics_filters`
    annotation indicating if the sample falls a certain number of MAD outside the
    distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4. We
        modify this for 'n_singleton_residual' to be (math.inf, 8.0) and for
        'r_het_hom_var_residual' to be (math.inf, 4.0). The math.inf (infinity) is used
        to prevent a lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT.
    :param qc_metrics: Specific metrics to use for outlier detection.
    :param pop_scores_expr: Optional expression with population PCA scores.
    :param platform_scores_expr: Optional expression with platform PCA scores.
    :param platform_expr: Optional expression with platform assignment.
    :param regress_pop_n_pcs: Number of population PCA scores to use in regression.
        Default is 30.
    :param regress_platform_n_pcs: Number of platform PCA scores to use in regression.
        Default is 9.
    :param regress_per_platform: Whether to perform the regression per platform.
        Default is False.
    :return: Table with regression residuals and outlier filters.
    """
    if pop_scores_expr is None and platform_expr is None:
        raise ValueError(
            "At least one of 'pop_scores_expr' or 'platform_scores_expr' must be "
            "specified!"
        )
    if regress_per_platform and platform_expr is None:
        raise ValueError(
            "When using 'regress_per_platform', 'platform_expr' must be specified!"
        )
    logger.info(
        "Computing QC metrics outlier filters with PC regression, using metrics: %s",
        ", ".join(qc_metrics),
    )
    ann_expr = {"scores": hl.empty_array(hl.tfloat64)}
    global_expr = {}
    log_str = []
    if pop_scores_expr is not None:
        ann_expr["scores"] = ann_expr["scores"].extend(
            pop_scores_expr[:regress_pop_n_pcs]
        )
        log_str.append("pop PCs")
        global_expr["regress_pop_n_pcs"] = regress_pop_n_pcs
    if regress_per_platform:
        ann_expr["platform"] = platform_expr
        log_str.append("is stratified by platform")
    elif platform_scores_expr is not None:
        ann_expr["scores"] = ann_expr["scores"].extend(
            platform_scores_expr[:regress_platform_n_pcs]
        )
        log_str.append("platform PCs")
        global_expr["regress_platform_n_pcs"] = regress_platform_n_pcs

    sample_qc_ht = sample_qc_ht.annotate(**ann_expr)

    logger.info(
        "QC metric residuals are being computed using a regression with %s.",
        " and ".join(log_str),
    )
    sample_qc_res_ht = compute_qc_metrics_residuals(
        sample_qc_ht,
        pc_scores=sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        strata={"platform": sample_qc_ht.platform} if regress_per_platform else None,
    )
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_res_ht,
        qc_metrics=dict(sample_qc_res_ht.row_value),
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
        },
        strata={"platform": sample_qc_ht[sample_qc_res_ht.key].platform}
        if regress_per_platform
        else None,
    )
    sample_qc_res_ht = sample_qc_res_ht.annotate(**filter_ht[sample_qc_res_ht.key])
    filter_ht = sample_qc_res_ht.select_globals(
        **filter_ht.index_globals(),
        **global_expr,
        regress_per_platform=regress_per_platform,
    )

    return filter_ht


def apply_nearest_neighbor_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    nn_ht: hl.Table,
) -> hl.Table:
    """
    Use nearest neighbor QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to compute the median, MAD, and upper and
    lower thresholds for each `qc_metrics` of the samples nearest neighbors and return
    a `qc_metrics_filters` annotation indicating if the sample falls a certain number
    of MAD outside the nearest neighbors metric distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4. We
        modify this for 'n_singleton' to be (math.inf, 8.0) and for
        'r_het_hom_var' and 'r_ti_tv_singleton' to be (math.inf, 4.0). The
        math.inf is used to prevent a lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT.
    :param qc_metrics: Specific metrics to use for outlier detection.
    :param nn_ht: Table with nearest neighbors annotation 'nearest_neighbors'. This
        annotation should be a List of samples that are the nearest neighbors of each
        sample.
    :return: Table with nearest neighbor outlier filters.
    """
    logger.info(
        "Computing QC metrics outlier filters by nearest neighbor comparison, using "
        "metrics: %s",
        ", ".join(qc_metrics),
    )
    sample_qc_ht = sample_qc_ht.annotate(
        nearest_neighbors=nn_ht[sample_qc_ht.key].nearest_neighbors
    )
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
        },
        comparison_sample_expr=sample_qc_ht.nearest_neighbors,
    )
    filter_ht = filter_ht.select_globals(**nn_ht.index_globals())

    return filter_ht


def apply_n_singleton_filter_to_r_ti_tv_singleton(
    ht: hl.Table,
    sample_qc_ht: hl.Table,
    filtering_method: str,
    ann_exprs: Dict[str, hl.expr.Expression],
    **kwargs: Any,
) -> hl.Table:
    """
    Apply n_singleton median sample filter and update the r_ti_tv_singleton outlier filtering.

    This function takes a Table (`ht`) returned from one of the filtering functions:
    `apply_stratified_filtering_method`, `apply_regressed_filtering_method`, or
    `apply_nearest_neighbor_filtering_method`, which must include filtering information
    for 'n_singleton' or 'n_singleton_residuals'. Median values for 'n_singleton' or
    'n_singleton_residuals' are extracted from 'qc_metrics_stats' in `ht`.
    'qc_metrics_stats' takes different forms depending on `filtering_method` and the
    use of 'strata' (defined as a global on `ht` if used) in filtering (described
    below). Then the `sample_qc_ht` is filtered to only samples with 'n_singleton'
    (or 'n_singleton_residuals') > the median in 'qc_metrics_stats' that is relevant to
    the sample. This filtered `sample_qc_ht` is rerun through the filtering function
    for `filtering_method` using only 'r_ti_tv_singleton', and `ht` is updated with new
    global and row values relevant to 'r_ti_tv_singleton' filtering.

    If `filtering_method` is 'stratified':
        - 'strata' is defined as a global annotation on `ht` in the form:
          tuple<str, ...> where the number of string elements in the tuple depends on
          the number of strata.
        - 'qc_metrics_stats' is defined as a global annotation on `ht` in the form:
          dict<strata, struct {
            metric: struct {
              'median': float64, 'mad': float64, 'lower': float64, 'upper': float64
            }
          }>
    If `filtering_method` is 'regressed':
        - 'qc_metrics_stats' is defined as a global annotation on `ht.
        - If platform stratification was performed, the 'strata' global annotation is
          tuple<'platform'>, and 'qc_metrics_stats' is in the same form as described
          for 'stratified', except metric -> metric_residual.
        - If there was no platform stratification, 'qc_metrics_stats' is in the form:
          struct {
            metric_residual: struct {
              'median': float64, 'mad': float64, 'lower': float64, 'upper': float64
            }
          }
    If `filtering_method` is 'nearest_neighbor':
        - 'qc_metrics_stats' is defined as a row annotation on `ht` in the form:
          struct {
            metric: struct {
              'median': float64, 'mad': float64, 'lower': float64, 'upper': float64
            }
          }

    :param ht: Input filtering Table with n_singleton or n_singleton_residuals medians
        computed.
    :param sample_qc_ht: Sample QC HT.
    :param filtering_method: The filtering method to apply to `sample_qc_ht`. One of
        'stratified', 'regressed', or 'nearest_neighbors'.
    :param ann_exprs: Dictionary of the expressions for population assignment and/or
        platform assignment on `sample_qc_ht` to be used if filtering is stratified.
    :return: Table with outlier filter annotations updated for r_ti_tv_singleton or
        r_ti_tv_singleton_residuals.
    """
    filtering_methods = {"stratified", "regressed", "nearest_neighbors"}
    if filtering_method not in filtering_methods:
        raise ValueError(
            f"Filtering method must be one of: {', '.join(filtering_methods)}"
        )

    # Setup repeatedly used variables.
    median_filter_metric = "n_singleton"
    update_metric = "r_ti_tv_singleton"
    filtering_qc_metrics = [update_metric]
    empty_stats_struct = hl.struct(
        median=hl.missing(hl.tfloat64),
        mad=hl.missing(hl.tfloat64),
        lower=hl.missing(hl.tfloat64),
        upper=hl.missing(hl.tfloat64),
    )

    # Annotate the sample QC Table with annotations provided in 'ann_expr'.
    sample_qc_ht = sample_qc_ht.annotate(**ann_exprs)

    # Extract the strata annotations from input Table globals (pop and/or platform
    # if present).
    ht_globals = hl.eval(ht.index_globals())
    if "strata" in ht_globals:
        strata = ht_globals["strata"]
    else:
        strata = None

    # If filtering method is 'nearest_neighbors', the 'qc_metrics_stats' is per sample,
    # so it is stored in the rows of 'ht' instead of the globals.
    if filtering_method == "nearest_neighbors":
        qc_metrics_stats = ht[sample_qc_ht.key].qc_metrics_stats
    else:
        qc_metrics_stats = ht_globals["qc_metrics_stats"]

    # If the filtering method is 'regressed', the metrics will have '_residual' on the
    # end of the annotation name and be stored in 'ht' instead of 'sample_qc_ht'.
    if filtering_method == "regressed":
        median_filter_metric += "_residual"
        median_filter_metric_expr = ht[sample_qc_ht.key][median_filter_metric]
        update_metric += "_residual"
    else:
        median_filter_metric_expr = sample_qc_ht[median_filter_metric]

    # Determine the expression for median cutoffs of 'median_filter_metric'.
    # When stratification is applied, 'qc_metrics_stats' has a median value per strata.
    if median_filter_metric in qc_metrics_stats:
        median_expr = qc_metrics_stats[median_filter_metric].median
    else:
        medians = hl.literal(
            {
                strat: cutoff_dict[median_filter_metric].median
                for strat, cutoff_dict in qc_metrics_stats.items()
                if cutoff_dict[median_filter_metric].median is not None
            }
        )
        median_expr = medians.get(hl.tuple([sample_qc_ht[x] for x in strata]))

    # Filter 'sample_qc_ht' to only samples that have 'median_filter_metric' greater
    # than the 'median_expr'. These are the only samples that should be included in the
    # updated outlier filtering.
    sample_qc_ht = sample_qc_ht.filter(
        median_filter_metric_expr > median_expr
    ).select_globals()

    filter_ht = apply_filter(
        filtering_method=filtering_method,
        sample_qc_ht=sample_qc_ht,
        qc_metrics=filtering_qc_metrics,
        apply_r_ti_tv_singleton_filter=False,
        **{f"{ann}_expr": sample_qc_ht[ann] for ann in ann_exprs},
        **kwargs,
    )

    # For all filtering methods other than nearest_neighbors, the global
    # 'qc_metrics_stats' needs to be updated for the 'update_metric'.
    if filtering_method != "nearest_neighbors":
        updated_stats = hl.eval(filter_ht.qc_metrics_stats)
        if strata is None:
            qc_metrics_stats.annotate(**{update_metric: updated_stats[update_metric]})
        else:
            for strat in qc_metrics_stats:
                update_stats_struct = empty_stats_struct
                if strat in updated_stats:
                    update_stats_struct = updated_stats[strat].get(
                        update_metric, empty_stats_struct
                    )
                qc_metrics_stats[strat].annotate(**{update_metric: update_stats_struct})

        ht = ht.annotate_globals(qc_metrics_stats=qc_metrics_stats)

    # For all filtering methods: add a sample annotation indicating whether the sample
    # was included in the updated outlier filtering, update the 'fail_{update_metric}'
    # annotation, and remove 'update_metric' from 'qc_metrics_filters' if no longer
    # filtered or add it if the sample would now be filtered.
    ht_idx = filter_ht[ht.key]
    ann_expr = {
        f"{median_filter_metric}_median_filtered": hl.is_defined(ht_idx),
        f"fail_{update_metric}": hl.coalesce(ht_idx[f"fail_{update_metric}"], False),
        "qc_metrics_filters": (
            (ht.qc_metrics_filters - hl.set({update_metric}))
            | hl.coalesce(ht_idx.qc_metrics_filters, hl.empty_set(hl.tstr))
        ),
    }

    # For the 'regressed' filtering method, residuals were recomputed, so update them.
    if filtering_method == "regressed":
        ann_expr[update_metric] = ht_idx[update_metric]

    # For the 'nearest_neighbors' filtering method, the per sample 'qc_metrics_stats'
    # annotation needs to be updated with the new medians and cutoffs for
    # 'update_metric'.
    if filtering_method == "nearest_neighbors":
        ann_expr["qc_metrics_stats"] = ht.qc_metrics_stats.annotate(
            **{
                update_metric: hl.coalesce(
                    ht_idx.qc_metrics_stats[update_metric],
                    ht.qc_metrics_stats[update_metric],
                )
            }
        )

    ht = ht.annotate(**ann_expr)

    return ht


def create_finalized_outlier_filter_ht(
    finalized_outlier_hts: Dict[str, hl.Table],
    qc_metrics: List[str],
    ensemble_operator: str = "and",
) -> hl.Table:
    """
    Create the finalized outlier filtering Table.

    .. note::

        `qc_metrics` must be the same as, or a subset of the original QC metrics used
        to create each Table in `finalized_outlier_hts`.

    Tables included in `finalized_outlier_hts` are reformatted to include only
    information from metrics in`qc_metrics`, and the 'qc_metrics_filters' annotation
    is modified to include only metrics in `qc_metrics`.

    If `finalized_outlier_hts` includes multiple Tables, the reformatted Table
    annotations are grouped under top level annotations specified by the keys of
    `finalized_outlier_hts`. The following annotations are also added:
        - `qc_metrics_fail` - includes the combined (using `ensemble_operator`) outlier
          fail boolean for each metric.
        - `qc_metrics_filters` - The combined set of QC metrics that a sample is an
          outlier for using `ensemble_operator` for the combination ('and' uses set
          intersection, 'or' uses set union).

    Finally, an `outlier_filtered` annotation is added indicating if a sample has any
    metrics in `qc_metrics_filters`.

    :param finalized_outlier_hts: Dictionary containing Tables to use for the finalized
        outlier filtering, keyed by a string to use as a top level annotation to group
        the annotations and global annotations under.
    :param qc_metrics: List of sample QC metrics to use for the finalized outlier
        filtering.
    :param ensemble_operator: Logical operator to use for combining filtering methods
        when multiple Tables are in `finalized_outlier_hts`. Options are ['or', 'and'].
         Default is 'and'.
    :return: Finalized outlier filtering Table.
    """
    if ensemble_operator == "and":
        fail_combine_func = hl.all
    elif ensemble_operator == "or":
        fail_combine_func = hl.any
    else:
        raise ValueError("'ensemble_operator' must be one of: 'and' or 'or'!")

    def _ensemble_func(x: hl.SetExpression, y: hl.SetExpression) -> hl.SetExpression:
        """
        Combine sets `x` and `y`.

        If x is missing return `y`.

        Logical operator for combination is determined by external `ensemble_operator`
        variable. Uses set intersection when `ensemble_operator` is 'and', and set
        union when `ensemble_operator` is 'or'.

        :param x: Input set to combine with `y`.
        :param y: Input set to combine with `x`.
        :return: Combined SetExpression.
        """
        return (
            hl.case()
            .when(hl.is_missing(x), y)
            .when(ensemble_operator == "and", x.intersection(y))
            .when(ensemble_operator == "or", x.union(y))
            .or_missing()
        )

    num_methods = len(finalized_outlier_hts)
    hts = []
    fail_map = defaultdict(dict)
    for filter_method, ht in finalized_outlier_hts.items():
        logger.info(
            "Reformatting Hail Table for %s outlier filtering method...", filter_method
        )

        # Determine the annotations that should be kept in the final filter Table
        annotations = list(ht.row.keys())
        fail_annotations_keep = set([])
        residual_annotations_keep = set([])
        qc_metrics_filters_keep = set([])
        for m in qc_metrics:
            if f"fail_{m}" in annotations:
                fail_map[filter_method][m] = f"fail_{m}"
                fail_annotations_keep.add(f"fail_{m}")
                qc_metrics_filters_keep.add(m)
            elif f"fail_{m}_residual" in annotations:
                fail_map[filter_method][m] = f"fail_{m}_residual"
                residual_annotations_keep.add(f"{m}_residual")
                fail_annotations_keep.add(f"fail_{m}_residual")
                qc_metrics_filters_keep.add(m)
            else:
                raise ValueError(
                    "The following metric is not found in the Table(s) provided for "
                    f"outlier filtering finalization: {m}."
                )

        # If residuals are annotated on the filtering Table, re-annotate under
        # 'qc_metrics_residuals' rather than keeping them top level
        select_expr = {}
        if residual_annotations_keep:
            select_expr["qc_metrics_residuals"] = ht[ht.key].select(
                *residual_annotations_keep
            )

        # Nearest neighbors filtering has 'qc_metrics_stats' as a row annotation
        # instead of a global annotation. This adds it to the final Table.
        if "qc_metrics_stats" in ht.row:
            select_expr["qc_metrics_stats"] = ht.qc_metrics_stats.select(
                *qc_metrics_filters_keep
            )

        # Group all filter fail annotations under 'qc_metrics_fail' rather than
        # keeping them top level
        select_expr.update(
            {
                "qc_metrics_fail": ht[ht.key].select(*fail_annotations_keep),
                "qc_metrics_filters": (
                    hl.literal(qc_metrics_filters_keep)
                    & ht.qc_metrics_filters.map(lambda x: x.replace("_residual", ""))
                ),
            }
        )
        ht = ht.select(**select_expr)

        def _update_globals(
            global_ann: Union[hl.struct, hl.dict],
            metrics_keep: List[str],
        ) -> Union[hl.struct, hl.dict]:
            """
            Update a global annotation to have values for only the requested metrics.

            :param global_ann: Dict or Struct for global annotation to modify.
            :param metrics_keep: Metrics to keep in the global annotation.
            :return: Updated global annotation Dict or Struct.
            """
            if isinstance(global_ann, hl.tstruct):
                global_ann = global_ann.select(*metrics_keep)
            else:
                global_ann = {
                    s: hl.struct(**{x: global_ann[s][x] for x in metrics_keep})
                    for s in hl.eval(global_ann.keys())
                }
            return global_ann

        # If multiple filtering Tables are provided, group all Table annotations under a
        # filtering method annotation rather than keeping it top level. This is needed
        # for joining requested filtering method Tables.
        if num_methods > 1:
            ht = ht.select(**{filter_method: ht[ht.key]})
            filter_globals = ht.index_globals()
            updated_globals = {}
            for g in filter_globals:
                if g == "qc_metrics_stats":
                    update_g = _update_globals(
                        filter_globals[g],
                        [x.split("fail_")[1] for x in fail_annotations_keep],
                    )
                elif g == "lms":
                    update_g = _update_globals(
                        filter_globals[g],
                        qc_metrics_filters_keep,
                    )
                else:
                    update_g = filter_globals[g]
                updated_globals[g] = update_g

            ht = ht.select_globals(
                **{f"{filter_method}_globals": hl.struct(**updated_globals)}
            )
            hts.append(ht)

    # If multiple filtering Tables are provided, join all filtering Tables, and make
    # top level annotations for 'qc_metrics_fail' and 'qc_metrics_filters' based on the
    # combination of all methods using 'ensemble_operator' as the combination method.
    if num_methods > 1:
        ht = hts[0]
        for ht2 in hts[1:]:
            ht = ht.join(ht2, how="outer")

        ht = ht.annotate(
            qc_metrics_fail=hl.struct(
                **{
                    m: fail_combine_func(
                        [
                            ht[f]["qc_metrics_fail"][fail_map[f][m]]
                            for f in finalized_outlier_hts
                        ]
                    )
                    for m in qc_metrics
                }
            ),
            qc_metrics_filters=hl.fold(
                _ensemble_func,
                hl.missing(hl.tset(hl.tstr)),
                [ht[f]["qc_metrics_filters"] for f in finalized_outlier_hts],
            ),
        )

    ht = ht.annotate(outlier_filtered=hl.len(ht.qc_metrics_filters) > 0)
    ht = ht.annotate_globals(
        filter_qc_metrics=qc_metrics, ensemble_combination_operator=ensemble_operator
    )

    return ht


def get_outlier_filtering_resources(
    args: argparse.Namespace,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the outlier filtering pipeline.

    :param args: argparse Namespace for arguments passed to the outlier_filtering.py
        script.
    :return: PipelineResourceCollection containing resources for all steps of the
        outlier filtering pipeline.
    """
    test = args.test
    overwrite = args.overwrite
    include_unreleasable_samples = args.include_unreleasable_samples

    # Adding resources from previous scripts that are used by multiple steps in the
    # outlier filtering pipeline.
    sample_qc_ht = get_sample_qc("under_three_alt_alleles")

    pop_ht = get_pop_ht()
    pop_scores_ht = ancestry_pca_scores(
        include_unreleasable_samples=include_unreleasable_samples
    )

    sample_qc_input = {
        "hard_filters.py --sample-qc": {"sample_qc_ht": sample_qc_ht},
    }
    pop_scores_input = {
        "assign_ancestry.py --run-pca": {"pop_scores_ht": pop_scores_ht}
    }
    pop_assign_input = {"assign_ancestry.py --assign-pops": {"pop_ht": pop_ht}}
    platform_input = {
        "platform_inference.py --assign-platforms": {"platform_ht": platform}
    }
    joint_qc_meta_input = {
        "generate_qc_mt.py --generate-qc-meta": {"joint_qc_meta": joint_qc_meta}
    }

    # Initialize outlier filtering pipeline resource collection.
    outlier_filtering_pipeline = PipelineResourceCollection(
        pipeline_name="outlier_filtering",
        pipeline_resources={
            **sample_qc_input,
            **pop_scores_input,
            **pop_assign_input,
            **platform_input,
            **joint_qc_meta_input,
        },
        overwrite=overwrite,
    )

    # Create resource collection for each step of the outlier filtering pipeline.
    apply_regressed_filters = PipelineStepResourceCollection(
        "--apply-regressed-filters",
        output_resources={
            "regressed_filter_ht": regressed_filtering(
                test=test,
                pop_pc_regressed=args.regress_population,
                platform_pc_regressed=args.regress_platform,
                platform_stratified=args.regress_per_platform,
            )
        },
        input_resources={
            **sample_qc_input,
            **pop_scores_input,
            **pop_assign_input,
            **platform_input,
            **joint_qc_meta_input,
        },
    )
    apply_stratified_filters = PipelineStepResourceCollection(
        "--apply-stratified-filters",
        output_resources={
            "stratified_filter_ht": stratified_filtering(
                test=test,
                pop_stratified=args.stratify_population,
                platform_stratified=args.stratify_platform,
            )
        },
        input_resources={**sample_qc_input, **pop_assign_input, **platform_input},
    )
    determine_nearest_neighbors = PipelineStepResourceCollection(
        "--determine-nearest-neighbors",
        output_resources={
            "nn_ht": nearest_neighbors(
                test=test,
                platform_stratified=args.nearest_neighbors_per_platform,
                approximation=args.use_nearest_neighbors_approximation,
                include_unreleasable_samples=include_unreleasable_samples,
            )
        },
        input_resources={**sample_qc_input, **pop_assign_input, **platform_input},
    )
    apply_nearest_neighbor_filters = PipelineStepResourceCollection(
        "--apply-nearest-neighbor-filters",
        output_resources={"nn_filter_ht": nearest_neighbors_filtering(test=test)},
        pipeline_input_steps=[determine_nearest_neighbors],
        add_input_resources=sample_qc_input,
    )

    finalized_input_steps = []
    if args.apply_regressed_filters:
        finalized_input_steps.append(apply_regressed_filters)
    if args.apply_stratified_filters:
        finalized_input_steps.append(apply_stratified_filters)
    if args.apply_nearest_neighbor_filters:
        finalized_input_steps.append(apply_nearest_neighbor_filters)

    if args.create_finalized_outlier_filter and len(finalized_input_steps) == 0:
        raise ValueError(
            "At least one filtering method and relevant options must be supplied "
            "when using '--create-finalized-outlier-filter'"
        )

    create_finalized_outlier_filter = PipelineStepResourceCollection(
        "--create-finalized-outlier-filter",
        output_resources={"finalized_ht": finalized_outlier_filtering(test=test)},
        pipeline_input_steps=finalized_input_steps,
    )

    # Add all steps to the outlier filtering pipeline resource collection.
    outlier_filtering_pipeline.add_steps(
        {
            "apply_regressed_filters": apply_regressed_filters,
            "apply_stratified_filters": apply_stratified_filters,
            "determine_nearest_neighbors": determine_nearest_neighbors,
            "apply_nearest_neighbor_filters": apply_nearest_neighbor_filters,
            "create_finalized_outlier_filter": create_finalized_outlier_filter,
        }
    )

    return outlier_filtering_pipeline


def main(args):
    """Determine sample QC metric outliers that should be filtered."""
    hl.init(
        log="/outlier_filtering.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    overwrite = args.overwrite
    filtering_qc_metrics = args.filtering_qc_metrics
    apply_r_ti_tv_singleton_filter = args.apply_n_singleton_filter_to_r_ti_tv_singleton
    include_unreleasable_samples = args.include_unreleasable_samples
    nn_platform_stratified = args.nearest_neighbors_per_platform
    nn_approximation = args.use_nearest_neighbors_approximation

    if args.apply_n_singleton_filter_to_r_ti_tv_singleton:
        err_msg = ""
        for metric in {"n_singleton", "r_ti_tv_singleton"}:
            if metric not in filtering_qc_metrics:
                err_msg += (
                    "'--apply-n-singleton-filter-to-r-ti-tv-singleton' flag is set, but"
                    f" {metric} is not in requested 'filtering_qc_metrics'!\n"
                )
        if err_msg:
            raise ValueError(err_msg)

    outlier_resources = get_outlier_filtering_resources(args)
    outlier_resources.check_resource_existence()
    pop_ht = outlier_resources.pop_ht.ht()
    platform_ht = outlier_resources.platform_ht.ht()
    pop_scores_ht = outlier_resources.pop_scores_ht.ht()
    joint_qc_meta_ht = outlier_resources.joint_qc_meta.ht()
    sample_qc_ht = get_sample_qc_ht(
        outlier_resources.sample_qc_ht.ht(), test=args.test, seed=args.seed
    )

    if not include_unreleasable_samples:
        sample_qc_ht = sample_qc_ht.filter(
            joint_qc_meta_ht[sample_qc_ht.key].releasable
        )

    if args.create_finalized_outlier_filter and args.use_existing_filter_tables:
        rerun_filtering = False
    else:
        rerun_filtering = True

    if args.apply_regressed_filters and rerun_filtering:
        res = outlier_resources.apply_regressed_filters
        res.check_resource_existence()

        ht = apply_filter(
            filtering_method="regressed",
            apply_r_ti_tv_singleton_filter=apply_r_ti_tv_singleton_filter,
            sample_qc_ht=sample_qc_ht,
            qc_metrics=filtering_qc_metrics,
            pop_scores_ht=pop_scores_ht,
            platform_ht=platform_ht,
            regress_pop_n_pcs=(
                args.regress_pop_n_pcs if args.regress_population else None
            ),
            regress_platform_n_pcs=(
                args.regress_platform_n_pcs if args.regress_platform else None
            ),
            regress_per_platform=args.regress_per_platform,
        )
        ht.annotate_globals(
            include_unreleasable_samples=include_unreleasable_samples
        ).write(res.regressed_filter_ht.path, overwrite=overwrite)

    if args.apply_stratified_filters and rerun_filtering:
        res = outlier_resources.apply_stratified_filters
        res.check_resource_existence()

        ht = apply_filter(
            filtering_method="stratified",
            apply_r_ti_tv_singleton_filter=apply_r_ti_tv_singleton_filter,
            sample_qc_ht=sample_qc_ht,
            qc_metrics=filtering_qc_metrics,
            pop_ht=pop_ht if args.stratify_population else None,
            platform_ht=platform_ht if args.stratify_platform else None,
        )
        ht.annotate_globals(
            include_unreleasable_samples=include_unreleasable_samples
        ).write(res.stratified_filter_ht.path, overwrite=overwrite)

    if args.determine_nearest_neighbors:
        res = outlier_resources.determine_nearest_neighbors
        res.check_resource_existence()

        if nn_platform_stratified:
            strata = {"platform": res.platform_ht.ht()[sample_qc_ht.key].qc_platform}
        else:
            strata = None
        ht = determine_nearest_neighbors(
            sample_qc_ht,
            pop_scores_ht[sample_qc_ht.key].scores,
            strata=strata,
            n_pcs=args.nearest_neighbors_pop_n_pcs,
            n_neighbors=args.n_nearest_neighbors,
            n_jobs=args.n_jobs,
            add_neighbor_distances=args.get_neighbor_distances,
            distance_metric=args.distance_metric,
            use_approximation=nn_approximation,
            n_trees=args.n_trees,
        )
        ht.annotate_globals(
            include_unreleasable_samples=include_unreleasable_samples
        ).write(res.nn_ht.path, overwrite=overwrite)

    if args.apply_nearest_neighbor_filters and rerun_filtering:
        res = outlier_resources.apply_nearest_neighbor_filters
        res.check_resource_existence()

        ht = apply_filter(
            filtering_method="nearest_neighbors",
            apply_r_ti_tv_singleton_filter=apply_r_ti_tv_singleton_filter,
            sample_qc_ht=sample_qc_ht,
            qc_metrics=filtering_qc_metrics,
            nn_ht=res.nn_ht.ht(),
        )
        ht.annotate_globals(
            include_unreleasable_samples=include_unreleasable_samples,
            nearest_neighbors_platform_stratified=nn_platform_stratified,
            nearest_neighbors_approximation=nn_approximation,
        ).write(res.nn_filter_ht.path, overwrite=overwrite)

    if args.create_finalized_outlier_filter:
        res = outlier_resources.create_finalized_outlier_filter
        res.check_resource_existence()
        # Reformat input step names for use as annotation labels.
        ht = create_finalized_outlier_filter_ht(
            {
                k.split("--apply-")[1].replace("-", "_"): v[0].ht()
                for k, v in res.input_resources.items()
            },
            qc_metrics=args.final_filtering_qc_metrics,
            ensemble_operator=args.ensemble_method_logical_operator,
        )
        ht.annotate_globals(
            include_unreleasable_samples=include_unreleasable_samples
        ).write(res.finalized_ht.path, overwrite=overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Test filtering using a random sample of 1%% of the dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--seed",
        help="Random seed for making random test dataset.",
        type=int,
        default=24,
    )
    parser.add_argument(
        "--include-unreleasable-samples",
        help=(
            "Include unreleasable samples in the nearest_neighbors determination, "
            "sample QC metric regressions, and/or sample QC filtering."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--filtering-qc-metrics",
        help="List of sample QC metrics to use for filtering.",
        default=[
            "n_snp",
            "n_singleton",
            "r_ti_tv",
            "r_insertion_deletion",
            "n_insertion",
            "n_deletion",
            "r_het_hom_var",
            "n_transition",
            "n_transversion",
            "r_ti_tv_singleton",
            "r_snp_indel",
        ],
        type=list,
        nargs="+",
    )
    parser.add_argument(
        "--apply-n-singleton-filter-to-r-ti-tv-singleton",
        help=(
            "Whether to apply the filtering method to only samples with: Number of "
            "singletons (or n_singleton residuals) > median(number of "
            "singletons/n_singleton residuals) where the median is computed on only "
            "samples within the same strata or on only the 50 nearest neighbors."
        ),
        action="store_true",
    )

    stratified_args = parser.add_argument_group(
        "Apply stratified filtering method.",
        "Arguments specific to stratified outlier filtering.",
    )
    stratified_args.add_argument(
        "--apply-stratified-filters",
        help="Apply stratified outlier filtering method.",
        action="store_true",
    )
    stratified_args.add_argument(
        "--stratify-population",
        help="Stratify by population for filtering.",
        action="store_true",
    )
    stratified_args.add_argument(
        "--stratify-platform",
        help="Stratify by platform for filtering.",
        action="store_true",
    )

    regressed_args = parser.add_argument_group(
        "Apply regression filtering method.",
        "Arguments specific to regression outlier filtering.",
    )
    regressed_args.add_argument(
        "--apply-regressed-filters",
        help="Apply regression outlier filtering method.",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-population",
        help="Use population PCs in sample QC metric regressions.",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-platform",
        help="Use platform PCs in sample QC metric regressions.",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-pop-n-pcs",
        help="Number of population PCs to use for sample QC metric regressions.",
        default=20,
        type=int,
    )
    regressed_args.add_argument(
        "--regress-platform-n-pcs",
        help="Number of platform PCs to use for sample QC metric regressions.",
        default=9,
        type=int,
    )
    regressed_args.add_argument(
        "--regress-per-platform",
        help="Perform sample QC metric regressions and outlier filtering per platform.",
        action="store_true",
    )

    nn_args = parser.add_argument_group(
        "Determine nearest neighbors for each sample.",
        "Arguments specific to determining nearest neighbors.",
    )
    nn_args.add_argument(
        "--determine-nearest-neighbors",
        help="Determine nearest neighbors for each sample.",
        action="store_true",
    )
    nn_args.add_argument(
        "--nearest-neighbors-pop-n-pcs",
        help="Number of population PCs to use for nearest neighbor determination.",
        default=20,
        type=int,
    )
    nn_per_platform = nn_args.add_argument(
        "--nearest-neighbors-per-platform",
        help=(
            "Stratify samples by platform assignment when determining the population "
            "PC nearest neighbors."
        ),
        action="store_true",
    )
    nn_args.add_argument(
        "--n-nearest-neighbors",
        help="Number of nearest neighbors to report for each sample.",
        default=50,
        type=int,
    )
    nn_args.add_argument(
        "--n-jobs",
        help=(
            "Number of threads to use when finding the nearest neighbors. Default "
            "is -1 which uses the number of CPUs on the head node -1."
        ),
        default=-1,
        type=int,
    )
    nn_args.add_argument(
        "--get-neighbor-distances",
        help="Return distances of the nearest neighbors too.",
        action="store_true",
    )
    nn_args.add_argument(
        "--distance-metric",
        help="Distance metric to use for nearest neighbor determination.",
        default="euclidean",
        type=str,
        choices=[
            "euclidean",
            "cityblock",
            "cosine",
            "haversine",
            "l1",
            "l2",
            "manhattan",
            "angular",
            "hamming",
            "dot",
        ],
    )
    nn_approx = nn_args.add_argument(
        "--use-nearest-neighbors-approximation",
        help=(
            "Whether to use the package Annoy to determine approximate nearest "
            "neighbors instead of using scikit-learn's `NearestNeighbors`. This method "
            "is faster, but only needed for very large datasets, for instance "
            "> 500,000 samples."
        ),
        action="store_true",
    )
    nn_args.add_argument(
        "--n-trees",
        help=(
            "Number of trees to build in the Annoy approximate nearest neighbors "
            "method. Only used if 'use-nearest-neighbors-approximation' is set. "
            "'n_trees' is provided during build time and affects the build time and "
            "the index size. A larger value will give more accurate results, but "
            "larger indexes."
        ),
        default=10,
        type=int,
    )

    nn_filter_args = parser.add_argument_group(
        "Apply nearest neighbors filtering method.",
        "Arguments specific to nearest neighbors outlier filtering.",
    )
    nn_filter_args.add_argument(
        "--apply-nearest-neighbor-filters",
        help="Apply nearest neighbors outlier filtering method.",
        action="store_true",
    )
    # Indicate that the --nearest-neighbors-per-platform and
    # --use-nearest-neighbors-approximation options apply to the
    # "nn_filter_args" argument group as well as the "nn_args" group.
    nn_filter_args._group_actions.append(nn_per_platform)
    nn_filter_args._group_actions.append(nn_approx)

    final_filter_args = parser.add_argument_group(
        "Create finalized outlier filtering Table.",
        "Arguments for finalizing outlier filtering.",
    )
    final_filter_args.add_argument(
        "--create-finalized-outlier-filter",
        help=(
            "Create the finalized outlier filtering Table. To specify filtering method "
            "to use for the finalized outlier filtering Table, use the parameter "
            "options from that method. e.g. To use the stratified outlier filtering "
            "with population stratification include the following arguments: "
            "--apply-stratified-filters --stratify-population. If an ensemble between "
            "two methods is desired, include arguments for both. The filtering Table "
            "for each method requested must already exist, the outlier filtering will "
            "not be rerun."
        ),
        action="store_true",
    )
    final_filter_args.add_argument(
        "--use-existing-filter-tables",
        help=(
            "Whether to use existing filter Tables if they exist instead of re-running "
            "any specified filtering methods before creating the finalized outlier "
            "filtering Table."
        ),
        action="store_true",
    )
    final_filter_args.add_argument(
        "--final-filtering-qc-metrics",
        help=(
            "List of sample QC metrics to use for the finalized outlier filtering. "
            "Note: This must be the same as, or a subset of `--filtering-qc-metrics` "
            "used when generating the filtering Table(s) desired for final filtering."
        ),
        default=[
            "r_ti_tv",
            "r_insertion_deletion",
            "r_het_hom_var",
            "r_ti_tv_singleton",
        ],
        type=list,
        nargs="+",
    )
    final_filter_args.add_argument(
        "--ensemble-method-logical-operator",
        help=(
            "Logical operator to use for combining filtering methods for the ensemble "
            "method."
        ),
        type=str,
        default="and",
        choices=["and", "or"],
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
