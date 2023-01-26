"""Script to determine sample QC metric outliers that should be filtered."""
import argparse
import logging
import math
from collections import defaultdict
from typing import Dict, List, Optional

import hail as hl
from gnomad.sample_qc.filtering import (
    compute_qc_metrics_residuals,
    compute_stratified_metrics_filter,
    determine_nearest_neighbors,
)
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import check_resource_existence
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
        - Sample 1% of the dataset id `test` is True

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


def apply_stratified_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    pop_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
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
        'r_het_hom_var' and 'r_ti_tv_singleton' to be (math.inf, 4.0). The
        math.inf (infinity) is used to prevent a lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT.
    :param qc_metrics: Specific metrics to use for outlier detection.
    :param pop_ht: Optional Table with population assignment annotation 'pop'.
    :param platform_ht: Optional Table with platform annotation 'qc_platform'.
    :return: Table with stratified metrics and filters.
    """
    if pop_ht is None and platform_ht is None:
        raise ValueError(
            "At least one of 'pop_scores_ht' or 'platform_ht' must be specified!"
        )
    logger.info(
        "Computing QC metrics outlier filters with stratification, using metrics: %s",
        ", ".join(qc_metrics),
    )
    strata = {}
    if pop_ht is not None:
        strata["pop"] = pop_ht[sample_qc_ht.key].pop
    if platform_ht is not None:
        strata["platform"] = platform_ht[sample_qc_ht.key].qc_platform
    logger.info("Outlier filtering is stratified by: %s", ", ".join(strata.keys()))
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        strata=strata,
        metric_threshold={
            "n_singleton": (math.inf, 8.0),
            "r_het_hom_var": (math.inf, 4.0),
            "r_ti_tv_singleton": (math.inf, 4.0),
        },
    )

    return filter_ht


def apply_regressed_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    pop_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
    regress_pop_n_pcs: int = 30,
    regress_platform_n_pcs: int = 9,
    regress_per_platform: bool = False,
    regression_include_unreleasable: bool = False,
    project_meta_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Compute sample QC metrics residuals after regressing out specified PCs and determine what samples are outliers that should be filtered.

    The following are all possible filtering options:
        - Include `pop_ht`: Regress population PCs only and determine outliers for each
          metric in `qc_metrics`.
        - Include `platform_ht`: Regress platform PCs only and determine outliers for
          each metric in `qc_metrics`.
        - Include `pop_ht` and `platform_ht`: Regress both population PCs and platform
          PCs and determine outliers for each metric in `qc_metrics`.
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
        'r_het_hom_var_residual' and 'r_ti_tv_singleton_residual' to be (math.inf, 4.0).
        The math.inf (infinity) is used to prevent a lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT.
    :param qc_metrics: Specific metrics to use for outlier detection.
    :param pop_ht: Optional Table with population PCA scores annotation 'scores'.
    :param platform_ht: Optional Table with platform annotation 'qc_platform' and
        platform PCA scores annotation 'scores'.
    :param regress_pop_n_pcs: Number of population PCA scores to use in regression.
        Default is 30.
    :param regress_platform_n_pcs: Number of platform PCA scores to use in regression.
        Default is 9.
    :param regress_per_platform: Whether to perform the regression per platform.
        Default is False.
    :param regression_include_unreleasable: Whether to include unreleasable samples in
        the regression. Default is False.
    :param project_meta_ht: Optional project meta Table that must be supplied if
        'regression_include_unreleasable' is False. By default, this Table should be
        included (because the default for `regression_include_unreleasable` is False).
    :return: Table with regression residuals and outlier filters.
    """
    if pop_ht is None and platform_ht is None:
        raise ValueError("At least one of 'pop_ht' or 'platform_ht' must be specified!")
    if regress_per_platform and (pop_ht is None or platform_ht is None):
        raise ValueError(
            "When using 'per_platform', 'pop_ht' and 'platform_ht' must "
            "both be specified!"
        )
    if not regression_include_unreleasable and project_meta_ht is None:
        raise ValueError(
            "When 'regression_include_unreleasable' is False, 'project_meta_ht' must "
            "be specified!"
        )
    logger.info(
        "Computing QC metrics outlier filters with PC regression, using metrics: %s",
        ", ".join(qc_metrics),
    )
    ann_expr = {"scores": hl.empty_array(hl.tfloat64)}
    global_expr = {}
    log_str = []
    if pop_ht is not None:
        ann_expr["scores"] = ann_expr["scores"].extend(
            pop_ht[sample_qc_ht.key].scores[:regress_pop_n_pcs]
        )
        log_str.append("pop PCs")
        global_expr["regress_pop_n_pcs"] = regress_pop_n_pcs
    if platform_ht is not None:
        if regress_per_platform:
            ann_expr["platform"] = platform_ht[sample_qc_ht.key].qc_platform
            log_str.append("is stratified by platform")
        else:
            ann_expr["scores"] = ann_expr["scores"].extend(
                platform_ht[sample_qc_ht.key].scores[:regress_platform_n_pcs]
            )
            log_str.append("platform PCs")
            global_expr["regress_platform_n_pcs"] = regress_platform_n_pcs

    if not regression_include_unreleasable:
        ann_expr["releasable"] = project_meta_ht[sample_qc_ht.key].releasable

    sample_qc_ht = sample_qc_ht.annotate(**ann_expr)

    logger.info(
        "QC metric residuals are being computed using a regression with %s.",
        " and ".join(log_str),
    )
    sample_qc_res_ht = compute_qc_metrics_residuals(
        sample_qc_ht,
        pc_scores=sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        regression_sample_inclusion_expr=None
        if regression_include_unreleasable
        else sample_qc_ht.releasable,
        strata={"platform": sample_qc_ht.platform} if regress_per_platform else None,
    )
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_res_ht,
        qc_metrics=dict(sample_qc_res_ht.row_value),
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
            "r_ti_tv_singleton_residual": (math.inf, 4.0),
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
        regression_include_unreleasable=regression_include_unreleasable,
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
            "r_ti_tv_singleton_residual": (math.inf, 4.0),
        },
        comparison_sample_expr=sample_qc_ht.nearest_neighbors,
    )
    filter_ht = filter_ht.select_globals(**nn_ht.index_globals())

    return filter_ht


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
        - `qc_metrics_filters` - The combined set of qc metrics that a sample is an
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
            select_expr["qc_metrics_stats"] = ht.qc_metrics_stats

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

        # If multiple filtering Tables are provided, group all Table annotations under a
        # filtering method annotation rather than keeping it top level. This is needed
        # for joining requested filtering method Tables.
        if num_methods > 1:
            ht = ht.select(**{filter_method: ht[ht.key]})
            ht = ht.select_globals(**{f"{filter_method}_globals": ht.index_globals()})
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


def main(args):
    """Determine sample QC metric outliers that should be filtered."""
    hl.init(
        log="/outlier_filtering.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test = args.test
    overwrite = args.overwrite
    filtering_qc_metrics = args.filtering_qc_metrics
    create_finalized_outlier_filter = args.create_finalized_outlier_filter
    sample_qc_ht = get_sample_qc("bi_allelic")
    pop_ht = get_pop_ht()
    include_unreleasable_samples = args.regression_include_unreleasable
    pop_scores_ht = ancestry_pca_scores(
        include_unreleasable_samples=include_unreleasable_samples
    )
    platform_ht = platform
    nn_per_platform = args.nearest_neighbors_per_platform
    nn_approximation = args.use_nearest_neighbors_approximation
    nn_ht = nearest_neighbors(
        test=test,
        platform_stratified=nn_per_platform,
        approximation=nn_approximation,
    )

    check_resource_existence(
        input_step_resources={
            "hard_filters.py --sample-qc": [sample_qc_ht],
        }
    )
    sample_qc_ht = get_sample_qc_ht(sample_qc_ht.ht(), test=test, seed=args.seed)
    finalized_outlier_input_resources = {}

    if args.apply_regressed_filters:
        regress_pop_n_pcs = args.regress_pop_n_pcs
        regress_platform_n_pcs = args.regress_platform_n_pcs
        regress_per_platform = args.regress_per_platform
        regress_population = args.regress_population
        regress_platform = args.regress_platform
        output = regressed_filtering(
            test=test,
            pop_pc_regressed=regress_population,
            platform_pc_regressed=regress_platform,
            platform_stratified=regress_per_platform,
            include_unreleasable_samples=include_unreleasable_samples,
        )
        if create_finalized_outlier_filter:
            finalized_outlier_input_resources["--apply-regressed-filters"] = [output]
        else:
            check_resource_existence(
                input_step_resources={
                    "assign_ancestry.py --run-pca": [pop_scores_ht],
                    "assign_ancestry.py --assign-pops": [platform_ht],
                    "platform_inference.py --assign-platforms": [joint_qc_meta],
                },
                output_step_resources={"--apply-regressed-filters": [output]},
                overwrite=overwrite,
            )
            if not regress_population:
                regress_pop_n_pcs = None
            if not regress_platform:
                regress_platform_n_pcs = None

            apply_regressed_filtering_method(
                sample_qc_ht,
                filtering_qc_metrics,
                pop_ht=pop_scores_ht.ht(),
                platform_ht=platform_ht.ht(),
                regress_pop_n_pcs=regress_pop_n_pcs,
                regress_platform_n_pcs=regress_platform_n_pcs,
                regress_per_platform=regress_per_platform,
                regression_include_unreleasable=include_unreleasable_samples,
                project_meta_ht=joint_qc_meta.ht(),
            ).write(output.path, overwrite=overwrite)

    if args.apply_stratified_filters:
        stratify_population = args.stratify_population
        stratify_platform = args.stratify_platform
        output = stratified_filtering(
            test=test,
            pop_stratified=stratify_population,
            platform_stratified=stratify_platform,
        )
        if create_finalized_outlier_filter:
            finalized_outlier_input_resources["--apply-stratified-filters"] = [output]
        else:
            check_resource_existence(
                input_step_resources={
                    "assign_ancestry.py --assign-pops": [pop_ht],
                    "platform_inference.py --assign-platforms": [platform_ht],
                },
                output_step_resources={"--apply-stratified-filters": [output]},
                overwrite=overwrite,
            )
            if not stratify_population:
                pop_ht = None
            if not stratify_platform:
                platform_ht = None

            apply_stratified_filtering_method(
                sample_qc_ht,
                filtering_qc_metrics,
                pop_ht=pop_ht.ht(),
                platform_ht=platform_ht.ht(),
            ).write(output.path, overwrite=overwrite)

    if args.determine_nearest_neighbors:
        check_resource_existence(
            input_step_resources={
                "assign_ancestry.py --run-pca": [pop_scores_ht],
                "platform_inference.py --assign-platforms": [platform_ht],
            },
            output_step_resources={"--determine-nearest-neighbors": [nn_ht]},
            overwrite=overwrite,
        )
        if nn_per_platform:
            strata = {"platform": platform_ht.ht()[sample_qc_ht.key].qc_platform}
        else:
            strata = None
        determine_nearest_neighbors(
            sample_qc_ht,
            pop_scores_ht.ht()[sample_qc_ht.key].scores,
            strata=strata,
            n_pcs=args.nearest_neighbors_pop_n_pcs,
            n_neighbors=args.n_nearest_neighbors,
            n_jobs=args.n_jobs,
            add_neighbor_distances=args.get_neighbor_distances,
            distance_metric=args.distance_metric,
            use_approximation=nn_approximation,
            n_trees=args.n_trees,
        ).write(nn_ht.path, overwrite=overwrite)

    if args.apply_nearest_neighbor_filters:
        output = nearest_neighbors_filtering(test=test)
        if create_finalized_outlier_filter:
            finalized_outlier_input_resources["--apply-nearest-neighbor-filters"] = [
                output
            ]
        else:
            check_resource_existence(
                input_step_resources={"--determine-nearest-neighbors": [nn_ht]},
                output_step_resources={"--apply-nearest-neighbor-filters": [output]},
                overwrite=overwrite,
            )
            apply_nearest_neighbor_filtering_method(
                sample_qc_ht,
                filtering_qc_metrics,
                nn_ht=nn_ht.ht(),
            ).write(output.path, overwrite=overwrite)

    if create_finalized_outlier_filter:
        if len(finalized_outlier_input_resources) == 0:
            raise ValueError(
                "At least one filtering method and relevant options must be supplied "
                "when using '--create-finalized-outlier-filter'"
            )
        output = finalized_outlier_filtering(test=test)
        check_resource_existence(
            input_step_resources=finalized_outlier_input_resources,
            output_step_resources={"--create-finalized-outlier-filter": [output]},
            overwrite=overwrite,
        )
        create_finalized_outlier_filter_ht(
            {
                k.replace("--apply-", "").replace("-", "_"): v[0].ht()
                for k, v in finalized_outlier_input_resources.items()
            },
            qc_metrics=args.final_filtering_qc_metrics,
            ensemble_operator=args.ensemble_method_logical_operator,
        ).write(output.path, overwrite=overwrite)


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
        help="Test filtering using a random sample of 1% of the dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--seed",
        help="Random seed for making random test dataset.",
        type=int,
        default=24,
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
        "--regression-include-unreleasable",
        help="Include unreleasable samples in the sample QC metric regressions.",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-pop-n-pcs",
        help="Number of population PCs to use for sample QC metric regressions.",
        default=30,
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
        default=30,
        type=int,
    )
    nn_args.add_argument(
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
            "Number of threads to use when finding the nearest neighbors. Default"
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
    nn_args.add_argument(
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
