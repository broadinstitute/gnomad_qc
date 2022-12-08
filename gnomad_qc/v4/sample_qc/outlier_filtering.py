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


def apply_filtering_method(
    sample_qc_ht: hl.Table,
    filtering_qc_metrics: List[str],
    method: str = "stratified",
    pop_ht: Optional[hl.Table] = None,
    pop_scores_ht: Optional[hl.Table] = None,
    platform_ht: Optional[hl.Table] = None,
    regress_pop_n_pcs: int = 16,
    regress_platform_n_pcs: int = 9,
    include_unreleasable_samples: bool = False,
    project_meta_ht: Optional[hl.Table] = None,
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
    the regression

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    """
    # Convert tuples to lists, so we can find the index of the passed threshold.
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

    if method == "stratified":
        logger.info(
            "Computing stratified QC metrics filters using metrics: "
            + ", ".join(filtering_qc_metrics)
        )

        strata = {}
        if pop_ht is not None:
            strata["pop"] = pop_ht[sample_qc_ht.key].pop
        if platform_ht is not None:
            strata["platform"] = platform_ht[sample_qc_ht.key].qc_platform

        filter_ht = compute_stratified_metrics_filter(
            sample_qc_ht,
            qc_metrics={
                metric: sample_qc_ht[metric] for metric in filtering_qc_metrics
            },
            strata=strata,
            metric_threshold={"n_singleton": (4.0, 8.0)},
        )
    elif method == "regressed":
        # TODO: add option to only regress out population, and do this regression
        # per platform
        if pop_scores_ht is None and platform_ht is None:
            ValueError(
                "At least one of 'pop_scores_ht' or 'platform_ht' must be specified!"
            )

        ann_args = {"pc_scores": hl.empty_array(hl.tfloat64)}
        if pop_scores_ht is not None:
            ann_args["pc_scores"] = ann_args["pc_scores"].extend(
                pop_scores_ht[sample_qc_ht.key].scores[:regress_pop_n_pcs]
            )
        if platform_ht is not None:
            ann_args["pc_scores"] = ann_args["pc_scores"].extend(
                platform_ht[sample_qc_ht.key].scores[:regress_platform_n_pcs]
            )

        if not include_unreleasable_samples:
            ann_args["releasable"] = project_meta_ht[sample_qc_ht.key].releasable

        sample_qc_ht = sample_qc_ht.annotate(**ann_args)
        sample_qc_ht = compute_qc_metrics_residuals(
            sample_qc_ht,
            pc_scores=sample_qc_ht.pc_scores,
            qc_metrics={
                metric: sample_qc_ht[metric] for metric in filtering_qc_metrics
            },
            regression_sample_inclusion_expr=None
            if include_unreleasable_samples
            else sample_qc_ht.releasable,
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
        )
    elif method == "nearest neighbors":
        if nn_ht is None:
            ValueError(
                "When filtering 'method' is 'nearest neighbors', 'nn_ht' must be "
                "supplied!"
            )
        logger.info(
            "Computing stratified QC metrics filters using metrics: "
            + ", ".join(filtering_qc_metrics)
        )
        sample_qc_ht = sample_qc_ht.annotate(
            nearest_neighbors=nn_ht[sample_qc_ht.key].nearest_neighbors
        )
        sex_ht = hl.read_table(
            "gs://gnomad/v4.0/sample_qc/exomes/sex_inference/gnomad.exomes.v4.0.sex.ht"
        )
        filter_ht = compute_stratified_metrics_filter(
            sample_qc_ht,
            qc_metrics={
                metric: sample_qc_ht[metric] for metric in filtering_qc_metrics
            },
            metric_threshold={"n_singleton": (4.0, 8.0)},
            comparison_sample_expr=sample_qc_ht.nearest_neighbors,
        )
    else:
        raise ValueError(
            f"Filtering method {method} is not valid, method must be one of "
            "'stratified', 'regressed', or 'nearest neighbors'!"
        )

    return filter_ht


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
    sample_qc_ht = get_sample_qc("bi_allelic")
    # pop_ht = get_pop_ht(test=test)
    from gnomad.resources.resource_utils import TableResource

    pop_ht = TableResource(
        "gs://gnomad/v4.0/sample_qc/joint/ancestry_inference/gnomad.joint.v4.0.hgdp_tgp_training.pop.rf_w_16pcs_0.75_pop_prob.ht"
    )
    pop_scores_ht = ancestry_pca_scores(
        include_unreleasable_samples=args.include_unreleasable_samples
    )
    platform_ht = platform
    nn_ht = nearest_neighbors(test=test)

    check_resource_existence(
        input_step_resources={
            "hard_filters.py --sample-qc": [sample_qc_ht],
        }
    )

    sample_qc_ht = sample_qc_ht.ht()
    if test:
        sample_qc_ht = sample_qc_ht.sample(0.01)

    if args.apply_regressed_filters:
        output = regressed_filtering(test=test)
        check_resource_existence(
            input_step_resources={
                "assign_ancestry.py --run-pca": [pop_scores_ht],
                "assign_ancestry.py --assign-pops": [platform_ht],
                "platform_inference.py --assign-platforms": [joint_qc_meta],
            },
            output_step_resources={"--apply-regressed-filters": [output]},
            overwrite=overwrite,
        )
        if not population_correction:
            regress_pop_n_pcs = None
        if not platform_correction:
            regress_platform_n_pcs = None

        apply_filtering_method(
            sample_qc_ht,
            filtering_qc_metrics,
            method="regressed",
            pop_scores_ht=pop_scores_ht.ht(),
            platform_ht=platform_ht.ht(),
            regress_pop_n_pcs=regress_pop_n_pcs,
            regress_platform_n_pcs=regress_platform_n_pcs,
            project_meta_ht=joint_qc_meta.ht(),
        ).write(output.path, overwrite=overwrite)

    if args.apply_stratified_filters:
        output = stratified_filtering(test=test)
        check_resource_existence(
            input_step_resources={
                "assign_ancestry.py --assign-pops": [pop_ht],
                "platform_inference.py --assign-platforms": [platform_ht],
            },
            output_step_resources={"--apply-stratified-filters": [output]},
            overwrite=overwrite,
        )
        if not population_correction:
            pop_ht = None
        if not platform_correction:
            platform_ht = None

        apply_filtering_method(
            sample_qc_ht,
            filtering_qc_metrics,
            method="stratified",
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
        if platform_correction:
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
            use_approximation=args.use_nearest_neighbors_approximation,
            n_trees=args.n_trees,
        ).write(nn_ht.path, overwrite=overwrite)

    if args.apply_nearest_neighbor_filters:
        output = nearest_neighbors_filtering(test=test)
        check_resource_existence(
            input_step_resources={"--determine-nearest-neighbors": [nn_ht]},
            output_step_resources={"--apply-nearest-neighbor": [output]},
            overwrite=overwrite,
        )
        apply_filtering_method(
            sample_qc_ht,
            filtering_qc_metrics,
            method="nearest neighbors",
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

# Add interval QC options? Redo sample QC after interval filtering?
