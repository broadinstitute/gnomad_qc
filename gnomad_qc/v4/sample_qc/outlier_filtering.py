"""Add script docstring."""
import argparse
import logging
import math
from typing import List, Optional

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


def apply_stratified_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
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
        the regression.
    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    """
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
        metric_threshold={"n_singleton": (4.0, 8.0)},
    )

    return filter_ht


def apply_regressed_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    pop_scores_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
    regress_pop_n_pcs: int = 16,
    regress_platform_n_pcs: int = 9,
    regression_include_unreleasable: bool = False,
    regress_per_platform: bool = False,
    project_meta_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Use population stratified QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to run hl.sample_qc and return the metrics stratified by assigned
    population with a `qc_metrics_filters` annotation indicating if the sample falls a certain number of MAD
    outside the distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4, we modify this for n_singleton to be
        (4.0, 8.0)

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
        the regression.
    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    """
    if pop_scores_ht is None and platform_ht is None:
        ValueError(
            "At least one of 'pop_scores_ht' or 'platform_ht' must be specified!"
        )
    if regress_per_platform and (pop_scores_ht is None or platform_ht is None):
        ValueError(
            "When using 'per_platform', 'pop_scores_ht' and 'platform_ht' must "
            "both be specified!"
        )
    logger.info(
        "Computing QC metrics outlier filters with PC regression, using metrics: %s",
        ", ".join(qc_metrics),
    )
    ann_args = {"scores": hl.empty_array(hl.tfloat64)}
    log_str = []
    strata = None
    if pop_scores_ht is not None:
        ann_args["scores"] = ann_args["scores"].extend(
            pop_scores_ht[sample_qc_ht.key].scores[:regress_pop_n_pcs]
        )
        log_str.append("pop PCs")
    if platform_ht is not None:
        if regress_per_platform:
            strata["platform"] = platform_ht[sample_qc_ht.key].qc_platform
            log_str.append("is stratified by platform")
        else:
            ann_args["scores"] = ann_args["scores"].extend(
                platform_ht[sample_qc_ht.key].scores[:regress_platform_n_pcs]
            )
            log_str.append("platform PCs")
    logger.info(
        "QC metric residuals are being computed using a regression with %s.",
        "and ".join(log_str),
    )

    if not regression_include_unreleasable:
        ann_args["releasable"] = project_meta_ht[sample_qc_ht.key].releasable

    sample_qc_ht = sample_qc_ht.annotate(**ann_args)
    sample_qc_ht = compute_qc_metrics_residuals(
        sample_qc_ht,
        pc_scores=sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        regression_sample_inclusion_expr=None
        if regression_include_unreleasable
        else sample_qc_ht.releasable,
        strata=strata,
    )
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics=dict(sample_qc_ht.row_value),
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
        },
    )

    sample_qc_ht = sample_qc_ht.annotate(**filter_ht[sample_qc_ht.key])
    filter_ht = sample_qc_ht.annotate_globals(
        **filter_ht.index_globals(),
        regress_pop_n_pcs=regress_pop_n_pcs,
        regress_platform_n_pcs=regress_platform_n_pcs,
        regress_per_platform=regress_per_platform,
    )

    return filter_ht


def apply_nearest_neighbor_filtering_method(
    sample_qc_ht: hl.Table,
    qc_metrics: List[str],
    nn_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Use population stratified QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to run hl.sample_qc and return the metrics stratified by assigned
    population with a `qc_metrics_filters` annotation indicating if the sample falls a certain number of MAD
    outside the distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4, we modify this for n_singleton to be
        (4.0, 8.0)

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
        the regression.
    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    """
    logger.info(
        "Computing QC metrics outlier filters by nearest neighbor comparison, using "
        "metrics: %s",
        ", ".join(qc_metrics),
    )
    if nn_ht is None:
        ValueError(
            "When filtering 'method' is 'nearest neighbors', 'nn_ht' must be supplied!"
        )
    sample_qc_ht = sample_qc_ht.annotate(
        nearest_neighbors=nn_ht[sample_qc_ht.key].nearest_neighbors
    )
    filter_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht[metric] for metric in qc_metrics},
        metric_threshold={"n_singleton": (4.0, 8.0)},
        comparison_sample_expr=sample_qc_ht.nearest_neighbors,
    )

    return filter_ht


def main(args):
    """Add main docstring."""
    hl.init(
        log="/outlier_filtering.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # TODO: Add interval QC options? Redo sample QC after interval filtering?
    test = args.test
    overwrite = args.overwrite
    filtering_qc_metrics = args.filtering_qc_metrics.split(",")
    sample_qc_ht = get_sample_qc("bi_allelic")
    # pop_ht = get_pop_ht(test=test)
    from gnomad.resources.resource_utils import TableResource

    pop_ht = TableResource(
        "gs://gnomad/v4.0/sample_qc/joint/ancestry_inference/gnomad.joint.v4.0.hgdp_tgp_training.pop.rf_w_16pcs_0.75_pop_prob.ht"
    )
    pop_scores_ht = ancestry_pca_scores(
        include_unreleasable_samples=args.regression_include_unreleasable
    )
    platform_ht = platform
    nn_per_platform = args.nn_per_platform
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

    # Convert tuples to lists, so we can find the index of the passed threshold.
    sample_qc_ht = sample_qc_ht.ht()
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            f"bases_dp_over_{hl.eval(sample_qc_ht.dp_bins[i])}": sample_qc_ht.bases_over_dp_threshold[
                i
            ]
            for i in range(len(sample_qc_ht.dp_bins))
        },
    )
    # Exclude hard filtered samples from the sample QC MT.
    # They should not be included in the metric distribution stats.
    sample_qc_ht = sample_qc_ht.filter(
        hl.is_missing(hard_filtered_samples.ht()[sample_qc_ht.key])
    )
    # Add 'r_snp_indel' annotation the the sample QC HT.
    sample_qc_ht = sample_qc_ht.annotate(
        r_snp_indel=sample_qc_ht.n_snp
        / (sample_qc_ht.n_insertion + sample_qc_ht.n_deletion)
    )

    if test:
        sample_qc_ht = sample_qc_ht.sample(0.01)

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
        )
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
            pop_scores_ht=pop_scores_ht.ht(),
            platform_ht=platform_ht.ht(),
            regress_pop_n_pcs=regress_pop_n_pcs,
            regress_platform_n_pcs=regress_platform_n_pcs,
            regress_per_platform=regress_per_platform,
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
        check_resource_existence(
            input_step_resources={"--determine-nearest-neighbors": [nn_ht]},
            output_step_resources={"--apply-nearest-neighbor": [output]},
            overwrite=overwrite,
        )
        apply_nearest_neighbor_filtering_method(
            sample_qc_ht,
            filtering_qc_metrics,
            nn_ht=nn_ht.ht(),
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
    stratified_args = parser.add_argument_group(
        "",
        "",
    )
    stratified_args.add_argument(
        "--apply-stratified-filters",
        help="Compute per pop filtering.",
        action="store_true",
    )
    stratified_args.add_argument(
        "--stratify-population",
        help="",
        action="store_true",
    )
    stratified_args.add_argument(
        "--stratify-platform",
        help="",
        action="store_true",
    )
    regressed_args = parser.add_argument_group(
        "",
        "",
    )
    regressed_args.add_argument(
        "--apply-regressed-filters",
        help="Computes qc_metrics adjusted for pop.",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-population",
        help="",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-platform",
        help="",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regression-include-unreleasable",
        help="",
        action="store_true",
    )
    regressed_args.add_argument(
        "--regress-pop-n-pcs",
        help="Number of population PCs to use for qc metric regressions",
        default=30,
        type=int,
    )
    regressed_args.add_argument(
        "--regress-platform-n-pcs",
        help="Number of platform PCs to use for qc metric regressions",
        default=9,
        type=int,
    )
    regressed_args.add_argument(
        "--regress-per-platform",
        help="Number of platform PCs to use for qc metric regressions",
        default=9,
        type=int,
    )
    nn_args = parser.add_argument_group(
        "",
        "",
    )
    nn_args.add_argument(
        "--determine-nearest-neighbors",
        help="",
        action="store_true",
    )
    nn_args.add_argument(
        "--nearest-neighbors-pop-n-pcs",
        help="",
        default=16,
        type=int,
    )
    regressed_args.add_argument(
        "--nn-per-platform",
        help="Number of platform PCs to use for qc metric regressions",
        default=9,
        type=int,
    )
    nn_args.add_argument(
        "--n-nearest-neighbors",
        help="",
        default=50,
        type=int,
    )
    nn_args.add_argument(
        "--n-jobs",
        help="",
        default=-1,
        type=int,
    )
    nn_args.add_argument(
        "--get-neighbor-distances",
        help="",
        action="store_true",
    )
    nn_args.add_argument(
        "--distance-metric",
        help="",
        default="euclidean",
        type=str,
        # choices=[],
    )
    nn_args.add_argument(
        "--use-nearest-neighbors-approximation",
        help="",
        action="store_true",
    )
    nn_args.add_argument(
        "--n-trees",
        help="",
        default=10,
        type=int,
    )
    nn_filter_args = parser.add_argument_group(
        "",
        "",
    )
    nn_filter_args.add_argument(
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
