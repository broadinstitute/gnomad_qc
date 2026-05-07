"""Script to compute coverage, allele number, and quality histograms on all gnomAD v5 genomes (AoU v8 + updated gnomAD v4)."""

import argparse
import logging
import subprocess
from functools import reduce
from os import getenv
from typing import List, Optional, Tuple

import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.gnomad import CURRENT_GENOME_AN_RELEASE as v4_AN_RELEASE
from gnomad.resources.grch38.gnomad import (
    CURRENT_GENOME_COVERAGE_RELEASE as v4_COVERAGE_RELEASE,
)
from gnomad.resources.grch38.gnomad import DOWNSAMPLINGS, public_release
from gnomad.resources.grch38.reference_data import (
    telomeres_and_centromeres,
    vep_context,
)
from gnomad.utils.annotations import (
    annotate_downsamplings,
    build_freq_stratification_list,
    generate_freq_group_membership_array,
    merge_array_expressions,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.sparse_mt import (
    compute_stats_per_ref_site,
    get_allele_number_agg_func,
    get_coverage_agg_func,
)
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.meta import meta as v4_meta
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_aou_downsampling,
    group_membership,
    qual_hists,
)
from gnomad_qc.v5.resources.basics import (
    _get_batch_resource_kwargs,
    _init_hail,
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.meta import meta
from gnomad_qc.v5.resources.release import (
    release_all_sites_an_tsv_path,
    release_coverage_path,
    release_coverage_tsv_path,
)

logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("v5_coverage_and_an")
logger.setLevel(logging.INFO)


def get_downsampling_ht(ht: hl.Table) -> hl.Table:
    """
    Get Table with downsampling groups for all samples.

    v5 downsampling is only applied to the AoU dataset.
    Desired groups:
    - 10,000
    - 100,000
    - Genetic ancestry group sizes for AFR, AMR, NFE
    Note that the only desired genetic ancestry group sizes are AFR, AMR, and NFE,
    but code will also generate downsamplings for all other groups.

    :param ht: Input Table.
    :return: Table with downsampling groups.
    """
    logger.info(
        "Determining downsampling groups for AoU...",
    )
    downsamplings = DOWNSAMPLINGS["v5"]
    ht = annotate_downsamplings(
        ht, downsamplings, ht.genetic_ancestry_inference.gen_anc
    )
    return ht


def get_group_membership_ht(
    meta_ht: hl.Table,
    project: str,
    ds_ht: Optional[hl.Table] = None,
    reduce_min_aggs: bool = False,
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    :param meta_ht: Meta HT.
    :param project: Project name. Must be "aou" or "gnomad". If "gnomad", function will filter meta HT to only consent drop samples.
    :param ds_ht: Optional downsampling HT. Only used for AoU.
    :param reduce_min_aggs: If True, build the HT with `reduce_to_minimal_groups=True`. Downstream `compute_stats_per_ref_site` calls must then list every non-reducible annotation under `entry_agg_group_membership` and pass `reducible_aggs={"AN"}` (or similar) so AN gets expanded back to full strata.
    :return: Group membership HT.
    """
    if project == "aou":
        ht = generate_freq_group_membership_array(
            meta_ht,
            build_freq_stratification_list(
                sex_expr=meta_ht.sex_karyotype,
                gen_anc_expr=meta_ht.genetic_ancestry_inference.gen_anc,
                downsampling_expr=ds_ht[meta_ht.key].downsampling,
            ),
            reduce_to_minimal_groups=reduce_min_aggs,
            downsamplings=hl.eval(ds_ht.downsamplings),
            ds_gen_anc_counts=hl.eval(ds_ht.ds_gen_anc_counts),
        )

        ht = ht.annotate_globals(
            freq_meta=ht.freq_meta.map(
                lambda d: hl.dict(
                    d.items().map(
                        lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x)
                    )
                )
            ),
        )

    elif project == "gnomad":
        # Filter to v4 consent drop samples.
        # NOTE: Not using v5 project meta here because this part will be run in
        # Dataproc.
        ht = meta_ht.filter(
            meta_ht.release
            & (
                (meta_ht.project_meta.research_project_key == "RP-1061")
                | (meta_ht.project_meta.research_project_key == "RP-1411")
            )
        )
        ht = generate_freq_group_membership_array(
            ht,
            build_freq_stratification_list(
                sex_expr=ht.sex_imputation.sex_karyotype,
                gen_anc_expr=ht.population_inference.pop,
            ),
            reduce_to_minimal_groups=reduce_min_aggs,
        )

    return ht


def validate_vds(vds: hl.vds.VariantDataset) -> None:
    """
    Validate VDS before densify.

    Code is taken from https://github.com/hail-is/hail/blob/858f3ab30c2bcc46d6e57fdbfe408284b4b3de53/hail/python/hail/vds/variant_dataset.py#L271
    at suggestion from Chris Vittal.

    :param vds: Input VDS.
    :return: None; raises ValueError if VDS is not valid.
    """
    rd = vds.reference_data
    vd = vds.variant_data

    ref_cols = rd.col_key.collect()
    var_cols = vd.col_key.collect()

    if len(ref_cols) != len(var_cols):
        raise ValueError(
            f"mismatch in number of columns: reference data has {ref_cols} columns, variant data has {var_cols} columns"
        )

    if ref_cols != var_cols:
        first_mismatch = 0
        while ref_cols[first_mismatch] == var_cols[first_mismatch]:
            first_mismatch += 1
        raise ValueError(
            f"mismatch in columns keys: ref={ref_cols[first_mismatch]}, var={var_cols[first_mismatch]} at position {first_mismatch}"
        )


def _file_exists_for_env(path: str, environment: str) -> bool:
    """
    Check if a path exists, tolerant of permission errors in batch mode.

    On the batch backend, anonymous file probes against requester-pays buckets
    can raise permission errors before getting to "exists / does not exist."
    Treat those as "exists" so the chunk is skipped rather than re-run; the
    next stage will surface a real error if the file is actually broken.
    """
    from gnomad.utils.file_utils import file_exists

    try:
        return file_exists(path)
    except Exception as e:
        if environment == "batch":
            logger.warning(
                "file_exists check on %s raised %s; assuming exists.", path, e
            )
            return True
        raise


def _chunk_path(cov_and_an_ht_path: str, kind: str, idx: int) -> str:
    """
    Return a sibling per-chunk HT path under ``<cov_and_an_path>_chunks/``.

    ``kind`` is ``"chunk"`` for a per-chunk worker output or ``"group"``
    for an intermediate group-merge output.
    """
    base = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
    return f"{base}_chunks/{idx:08d}.{kind}.ht"


def _derive_chunk_locus_intervals(
    vds_filtered: hl.vds.VariantDataset,
    reference_genome: str = "GRCh38",
) -> List[hl.expr.IntervalExpression]:
    """
    Derive per-contig locus intervals covering the filtered VDS reference_data.

    The vep_context HT (used as ``ref_ht``) and the AoU VDS have independent
    partition layouts, so chunking each by its own partition index would not
    yield the same locus range. Instead, we derive the locus extent of the
    VDS chunk and filter ``ref_ht`` to that range so the locus-keyed join in
    ``compute_stats_per_ref_site`` produces correct output.
    """
    rd = vds_filtered.reference_data
    bounds = rd.aggregate_rows(
        hl.agg.group_by(
            rd.locus.contig,
            hl.struct(
                lo=hl.agg.min(rd.locus.position),
                hi=hl.agg.max(rd.locus.position),
            ),
        )
    )
    return [
        hl.parse_locus_interval(
            f"{contig}:{b.lo}-{b.hi + 1}",
            reference_genome=reference_genome,
        )
        for contig, b in bounds.items()
    ]


def compute_all_release_stats_per_ref_site(
    vds: hl.vds.VariantDataset,
    ref_ht: hl.Table,
    sex_karyotype_field: str,
    project: str,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    interval_ht: Optional[hl.Table] = None,
    group_membership_ht: Optional[hl.Table] = None,
    reduce_min_aggs: bool = False,
) -> hl.Table:
    """
    Compute coverage, allele number, and quality histograms per reference site.

    .. note::

        Running this function prior to calculating frequencies removes the need for an additional
        densify for frequency calculations.

    :param vds: Input VDS.
    :param ref_ht: Reference HT.
    :param sex_karyotype_field: Field name for sex karyotype.
    :param project: Project name.
    :param coverage_over_x_bins: List of boundaries for computing samples over X depth.
    :param interval_ht: Interval HT.
    :param group_membership_ht: Group membership HT.
    :param reduce_min_aggs: If True, pass `reducible_aggs={"AN"}` to `compute_stats_per_ref_site` so AN goes through the leaf-reduction path. Requires `group_membership_ht` to have been built with `reduce_to_minimal_groups=True`.
    :return: HT with allele number and quality histograms per reference site.
    """

    def _get_hists(qual_expr) -> hl.expr.Expression:
        return qual_hist_expr(
            gq_expr=qual_expr[0],
            dp_expr=qual_expr[1],
            adj_expr=qual_expr[2] == 1,
            split_adj_and_raw=True,
        )

    # Set up coverage bins.
    cov_bins = sorted(coverage_over_x_bins)
    rev_cov_bins = list(reversed(cov_bins))
    max_cov_bin = cov_bins[-1]
    cov_bins = hl.array(cov_bins)

    entry_agg_funcs = {
        "AN": get_allele_number_agg_func("LGT"),
        "coverage_stats": get_coverage_agg_func(dp_field="DP", max_cov_bin=max_cov_bin),
    }
    # Restrict `coverage_stats` to the global adj-filtered group only —
    # downstream code uses `coverage_stats[0]` exclusively (the global adj
    # group; transmuted into flat mean/sum/over_X fields), so per-strata
    # coverage aggregation is wasted work. We use `{"group": "adj"}` (and
    # not `"raw"`) because `compute_stats_per_ref_site` does NOT rewrite
    # group labels to "raw" when a pre-built `group_membership_ht` is
    # passed (the rewrite only happens on the `strata_expr` code path), so
    # `freq_meta[0]` is `{"group": "adj"}` and we need to match it.
    # AN is intentionally omitted so it still fans out across all strata,
    # which is what downstream consumers need.
    entry_agg_group_membership = {
        "coverage_stats": [{"group": "adj"}],
    }
    # Only compute qual hists for AoU.
    if project == "aou":
        entry_agg_funcs["qual_hists"] = (lambda t: [t.GQ, t.DP, t.adj], _get_hists)

        # Below we use just the raw group for qual hist computations because qual hists
        # has its own built-in adj filtering when adj is passed as an argument and will
        # produce both adj and raw histograms.
        entry_agg_group_membership["qual_hists"] = [{"group": "raw"}]

    logger.info(
        "Computing coverage, allele number, and optionally qual hists per reference site..."
    )

    vmt = vds.variant_data
    sex_expr = reduce(lambda x, field: x[field], sex_karyotype_field.split("."), vmt)
    vmt = vmt.annotate_cols(sex_karyotype=sex_expr)
    rmt = vds.reference_data
    if "LEN" not in rmt.entry:
        rmt = rmt.annotate_entries(LEN=rmt.END - rmt.locus.position + 1)
    vds = hl.vds.VariantDataset(rmt, vmt)

    # TODO: Save dense MT for gnomAD consent drop samples?
    ht = compute_stats_per_ref_site(
        vds,
        ref_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["GQ", "DP"],
        reducible_aggs={"AN"},
        entry_agg_group_membership=entry_agg_group_membership,
        sex_karyotype_field="sex_karyotype",
    )

    # This expression aggregates the DP counter in reverse order of the cov_bins and
    # computes the cumulative sum over them. It needs to be in reverse order because we
    # want the sum over samples covered by > X.
    def _cov_stats(
        cov_stat: hl.expr.StructExpression,
    ) -> hl.expr.StructExpression:
        # The coverage was already floored to the max_coverage_bin, so no more
        # aggregation is needed for the max bin.
        count_expr = cov_stat.coverage_counter
        max_bin_expr = hl.int32(count_expr.get(max_cov_bin, 0))

        # For each of the other bins, coverage is summed between the boundaries.
        bin_expr = hl.range(hl.len(cov_bins) - 1, 0, step=-1)
        bin_expr = bin_expr.map(
            lambda i: hl.sum(
                hl.range(cov_bins[i - 1], cov_bins[i]).map(
                    lambda j: hl.int32(count_expr.get(j, 0))
                )
            )
        )
        bin_expr = hl.cumulative_sum(hl.array([max_bin_expr]).extend(bin_expr))
        # NOTE: Keeping these as sample counts rather than fractions to join with
        # gnomAD v4 genomes.
        bin_expr = {f"over_{x}": bin_expr[i] for i, x in enumerate(rev_cov_bins)}
        return cov_stat.annotate(**bin_expr).drop("coverage_counter")

    # Keep coverage stats from global adj grouping (index 0) only.
    ht = ht.annotate_globals(
        coverage_stats_meta_sample_count=ht.strata_sample_count[0],
    )
    cov_stats_expr = _cov_stats(ht.coverage_stats[0])

    ht = ht.transmute(**cov_stats_expr)

    if project == "aou":
        # `qual_hists` as returned by `compute_stats_per_ref_site` is an array of length 1 so we drop the array here.
        ht = ht.annotate(qual_hists=ht.qual_hists[0])
    return ht


def _rename_cov_annotations(
    ht: hl.Table,
    project: str,
    sample_count: int,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
) -> hl.Table:
    """
    Rename coverage annotations prior to merging Tables.

    Function transforms mean back into sum using `sample_count` argument.

    :param ht: Input HT.
    :param project: Project name.
    :param sample_count: Number of samples in HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :return: Renamed HT.
    """
    # Transform mean back into sum.
    ht = ht.transmute(sum=ht.mean * sample_count)

    # Rename annotations to include project.
    row_fields = list(ht.row_value)
    rename_dict = {f: f"{f}_{project}" for f in row_fields}
    ht = ht.rename(rename_dict)
    if project == "gnomad" or project == "gnomad_release":
        # Revert v4 genomes fraction over X bins to sample count over X bins.
        ht = ht.transmute(
            **{
                f"over_{x}_{project}": (
                    hl.float64(ht[f"over_{x}_{project}"]) * sample_count
                )
                for x in coverage_over_x_bins
            }
        )
    return ht


def _merge_coverage_fields(
    ht: hl.Table,
    project_1: str,
    project_2: str,
    sample_count: int,
    operation: str,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
) -> hl.expr.DictExpression:
    """
    Merge coverage fields from two Tables.

    .. note::
        - Function does not merge `median_approx` fields.

    :param ht: Input HT. Must have annotations from both projects.
    :param project_1: First project name.
    :param project_2: Second project name.
    :param sample_count: Total sample count.
    :param operation: Operation to perform on the coverage fields. Must be "sum" or "diff".
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :return: Merged fields.
    """
    if operation == "diff":
        # Make sure over_x fields are float64s.
        merged_fields = {
            "sum_gnomad": ht[f"sum_{project_1}"] - ht[f"sum_{project_2}"],
            "total_DP_gnomad": (
                ht[f"total_DP_{project_1}"] - ht[f"total_DP_{project_2}"]
            ),
        }
        merged_fields.update(
            {
                f"over_{x}_gnomad": (
                    ht[f"over_{x}_{project_1}"] - ht[f"over_{x}_{project_2}"]
                )
                for x in coverage_over_x_bins
            }
        )
    else:
        merged_fields = {
            "mean": (ht[f"sum_{project_1}"] + ht[f"sum_{project_2}"]) / sample_count,
        }
        merged_fields.update(
            {
                f"over_{x}": (
                    (ht[f"over_{x}_{project_1}"] + ht[f"over_{x}_{project_2}"])
                    / sample_count
                )
                for x in coverage_over_x_bins
            }
        )
    return merged_fields


def merge_gnomad_coverage_hts(
    gnomad_ht: hl.Table,
    gnomad_release_ht: hl.Table,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    v4_count: int = 76215,
    consent_drop_count: int = 866,
) -> hl.Table:
    """
    Subtract consent drop samples from gnomAD v4 genomes release HT to create gnomAD v5 genomes coverage HT.

    :param gnomad_ht: gnomAD coverage HT (contains coverage for consent drop samples only).
    :param gnomad_release_ht: gnomAD v4 genomes coverage release HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :param v4_count: Number of release gnomAD v4 genome samples. Default is 76215.
    :param consent_drop_count: Number of consent drop gnomAD v4 genome samples. Default is 866.
    :return: gnomAD v5 genomes coverage HT.
    """
    logger.info(
        "Subtracting gnomAD v4 consent drop samples from gnomAD v4 genomes release HT..."
    )
    gnomad_ht = _rename_cov_annotations(
        gnomad_ht, "gnomad", consent_drop_count, coverage_over_x_bins
    )
    gnomad_release_ht = _rename_cov_annotations(
        gnomad_release_ht, "gnomad_release", v4_count, coverage_over_x_bins
    )
    gnomad_v5_count = v4_count - consent_drop_count
    logger.info("Total number of gnomAD v5 release genomes: %s", gnomad_v5_count)

    gnomad_ht = gnomad_ht.join(gnomad_release_ht, "right")
    merged_fields = _merge_coverage_fields(
        ht=gnomad_ht,
        project_1="gnomad_release",
        project_2="gnomad",
        sample_count=gnomad_v5_count,
        operation="diff",
    )
    gnomad_ht = gnomad_ht.transmute(**merged_fields)

    # Keep median_approx from v4 release.
    gnomad_ht = gnomad_ht.transmute(
        median_approx_gnomad=gnomad_ht.median_approx_gnomad_release
    )

    # Drop unnecessary globals and add back v5 gnomAD genomes count.
    gnomad_ht = gnomad_ht.select_globals()
    gnomad_ht = gnomad_ht.annotate_globals(
        coverage_stats_meta_sample_count=gnomad_v5_count,
    )
    return gnomad_ht


def join_aou_and_gnomad_coverage_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    gnomad_v5_count: int = 76215 - 866,
) -> hl.Table:
    """
    Join AoU and gnomAD coverage HTs for release.

    :param aou_ht: AoU coverage HT.
    :param gnomad_ht: gnomAD v5 genomes coverage HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :param gnomad_v5_count: Number of release gnomAD v5 genome samples. Default is 76215 - 866.
    :return: Joined HT.
    """
    aou_count = hl.eval(aou_ht.coverage_stats_meta_sample_count)
    logger.info("Total number of AoU v8 release samples: %s", aou_count)

    logger.info("Merging AoU and gnomAD v5 coverage HTs...")
    aou_ht = _rename_cov_annotations(aou_ht, "aou", aou_count, coverage_over_x_bins)
    v5_count = aou_count + gnomad_v5_count
    logger.info("Total number of AoU + gnomAD v5 release genomes: %s", v5_count)
    ht = aou_ht.join(gnomad_ht, "left")
    merged_fields = _merge_coverage_fields(
        ht=ht,
        project_1="aou",
        project_2="gnomad",
        sample_count=v5_count,
        operation="sum",
    )
    ht = ht.transmute(**merged_fields)
    return ht.select_globals()


def _rename_fields(
    ht: hl.Table, field_name: str, project: str, rename_globals: bool
) -> hl.Table:
    """
    Rename fields by adding project name prior to merging Tables.

    Used for AN and qual hists merging but not coverage because
    coverage does not have globals and also requires extra transformations
    to transform v4 mean back into sum.

    :param ht: Input HT.
    :param project: Project name.
    :param rename_globals: Whether to rename globals.
    :return: Renamed HT.
    """
    if rename_globals:
        rename_globals = {
            f"strata_meta_{project}": ht.strata_meta,
            f"strata_sample_count_{project}": ht.strata_sample_count,
        }
        ht = ht.transmute_globals(**rename_globals)
    rename_dict = {f"{field_name}": f"{field_name}_{project}"}
    return ht.rename(rename_dict)


def _merge_an_fields(
    ht: hl.Table, project_1: str, project_2: str, operation: str
) -> Tuple[
    hl.expr.ArrayExpression,
    hl.expr.ArrayExpression,
    hl.expr.DictExpression,
]:
    """
    Merge AN fields from two projects.

    :param ht: Input HT. Must have annotations from both projects.
    :param project_1: First project name.
    :param project_2: Second project name.
    :param operation: Operation to perform on the AN fields. Must be "sum" or "diff".
    :return: Joined AN, strata meta, and count arrays.
    """
    joint_an, joint_strata_meta, count_arrays_dict = merge_array_expressions(
        arrays=[ht[f"AN_{project_1}"], ht[f"AN_{project_2}"]],
        meta=[
            ht.index_globals()[f"strata_meta_{project_1}"],
            ht.index_globals()[f"strata_meta_{project_2}"],
        ],
        count_arrays={
            "counts": [
                ht.index_globals()[f"strata_sample_count_{project_1}"],
                ht.index_globals()[f"strata_sample_count_{project_2}"],
            ],
        },
        operation=operation,
    )
    return joint_an, joint_strata_meta, count_arrays_dict


def merge_gnomad_an_hts(
    gnomad_ht: hl.Table,
    gnomad_release_ht: hl.Table,
) -> hl.Table:
    """
    Subtract consent drop samples from gnomAD v4 genomes release HT to create gnomAD v5 genomes AN HT.

    :param gnomad_ht: gnomAD AN HT (contains AN for consent drop samples only).
    :param gnomad_release_ht: gnomAD v4 genomes release AN HT.
    :return: gnomAD v5 genomes AN HT.
    """
    logger.info(
        "Subtracting gnomAD v4 consent drop samples from gnomAD v4 genomes release HT..."
    )
    gnomad_ht = _rename_fields(gnomad_ht, "AN", "gnomad", rename_globals=True)
    gnomad_release_ht = _rename_fields(
        gnomad_release_ht, "AN", "gnomad_release", rename_globals=True
    )

    gnomad_ht = gnomad_ht.join(gnomad_release_ht, "right")
    joint_an, joint_strata_meta, count_arrays_dict = _merge_an_fields(
        ht=gnomad_ht,
        project_1="gnomad_release",
        project_2="gnomad",
        operation="diff",
    )
    gnomad_ht = gnomad_ht.annotate(AN_gnomad=joint_an)
    gnomad_ht = gnomad_ht.annotate_globals(
        strata_meta_gnomad=joint_strata_meta,
        strata_sample_count_gnomad=count_arrays_dict["counts"],
    )
    gnomad_ht = gnomad_ht.drop(
        "strata_meta_gnomad_release",
        "strata_sample_count_gnomad_release",
        "coverage_stats_meta_sample_count",
    )
    gnomad_ht = gnomad_ht.select("AN_gnomad")
    return gnomad_ht


def join_aou_and_gnomad_an_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
) -> hl.Table:
    """
    Join AoU and gnomAD AN HTs for release.

    :param aou_ht: AoU AN HT.
    :param gnomad_ht: gnomAD v5 genomes AN HT.
    :return: Joined HT.
    """
    aou_ht = _rename_fields(aou_ht, "AN", "aou", rename_globals=True)

    # TODO: should we only merge the overall adj AN (like coverage)?
    logger.info("Merging AoU and gnomAD v5 AN HTs...")
    ht = aou_ht.join(gnomad_ht, "left")
    ht = ht.checkpoint(new_temp_file("aou_and_gnomad_join", "ht"))

    joint_an, joint_strata_meta, count_arrays_dict = _merge_an_fields(
        ht=ht,
        project_1="aou",
        project_2="gnomad",
        operation="sum",
    )
    ht = ht.annotate(AN=joint_an)
    ht = ht.annotate_globals(
        strata_meta=joint_strata_meta,
        strata_sample_count=count_arrays_dict["counts"],
    )
    return ht


def join_aou_and_gnomad_qual_hists_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
) -> hl.Table:
    """
    Join AoU and gnomAD qual hists HTs for release.

    .. note::
        We did not compute qual hists for the gnomAD v4 genomes release
        (https://github.com/broadinstitute/gnomad_qc/blob/e65bdbb5768113c0129199a875d845da245690e2/gnomad_qc/v4/annotations/generate_freq_genomes.py#L1139).
        This means we will not also not recompute hists on the gnomAD v4 genomes for v5,
        which also means we will not subtract values from the samples to drop for consent reasons.

    :param aou_ht: AoU qual hists HT.
    :param gnomad_ht: gnomAD qual hists HT.
    :return: Joined HT.
    """
    aou_ht = _rename_fields(aou_ht, "qual_hists", "aou", rename_globals=False)
    gnomad_ht = _rename_fields(gnomad_ht, "qual_hists", "gnomad", rename_globals=False)
    ht = aou_ht.join(gnomad_ht, "left")
    qual_hists = [
        "gq_hist_all",
        "dp_hist_all",
    ]
    hist_structs = {
        "qual_hists": qual_hists,
        "raw_qual_hists": qual_hists,
    }
    hists_expr = {
        hist_struct: hl.struct(
            **{
                h: merge_histograms(
                    [
                        ht.qual_hists_aou[hist_struct][h],
                        ht.qual_hists_gnomad[hist_struct][h],
                    ],
                    operation="sum",
                )
                for h in hists
            }
        )
        for hist_struct, hists in hist_structs.items()
    }
    ht = ht.annotate(**hists_expr)
    return ht


def _filter_to_locus_bounds(target_ht: hl.Table, source_ht: hl.Table) -> hl.Table:
    """
    Filter `target_ht` to per-contig locus bounds derived from `source_ht`.

    Used in test runs where `target_ht` (e.g., a public release HT) and
    `source_ht` (e.g., a freshly-computed test HT) do not share partition
    layouts, so their first-N-partition slices do not overlap. Filtering to
    the per-contig min..max locus range of `source_ht` guarantees test
    join/merge overlap at low cost (one small aggregation + partition-pruned
    interval filter).

    :param target_ht: Table to filter.
    :param source_ht: Table whose locus ranges define the filter intervals.
    :return: `target_ht` filtered to the per-contig locus ranges of
        `source_ht`.
    """
    bounds = source_ht.aggregate(
        hl.agg.group_by(
            source_ht.locus.contig,
            hl.tuple(
                [
                    hl.agg.min(source_ht.locus.position),
                    hl.agg.max(source_ht.locus.position),
                ]
            ),
        )
    )
    intervals = [
        hl.parse_locus_interval(f"{contig}:{lo}-{hi + 1}", reference_genome="GRCh38")
        for contig, (lo, hi) in bounds.items()
    ]
    return hl.filter_intervals(target_ht, intervals)


def _resolve_partition_range(
    args: argparse.Namespace,
) -> Tuple[Optional[List[int]], int, int]:
    """
    Resolve the VDS partition range for a chunk-style invocation.

    Returns ``(partition_range, start, stop)`` where ``partition_range`` is a
    list usable as ``filter_partitions=`` (or ``None`` for the full VDS).
    Explicit ``--chunk-start`` / ``--chunk-stop`` bounds win when set;
    otherwise ``--test-2-partitions`` is treated as an alias for
    ``--chunk-start 0 --chunk-stop 2``.

    :param args: Parsed CLI args.
    :return: Tuple of (partition list or None, start, stop).
    """
    start = args.chunk_start
    stop = args.chunk_stop
    if stop > start:
        return list(range(start, stop)), start, stop
    if args.test_2_partitions:
        return list(range(2)), 0, 2
    return None, start, stop


def _build_chunk_ref_ht(
    vds_filtered: hl.vds.VariantDataset,
    project: str,
    repartition_factor: int,
    partition_count: int,
    test_chr22_chrx_chry: bool,
    chrom: Optional[List[str]],
) -> hl.Table:
    """
    Build the per-chunk ``ref_ht`` aligned to the VDS chunk's locus extent.

    Filters ``vep_context`` to the chunk's contigs+locus range (derived from
    the filtered VDS), strips to the locus key, removes telomeres/centromeres,
    and optionally explodes into ``partition_count * repartition_factor``
    sub-partitions to shrink per-task densify size.
    """
    ref_ht = vep_context.versions["105"].ht()
    if test_chr22_chrx_chry:
        ref_ht = hl.filter_intervals(
            ref_ht, [hl.parse_locus_interval(c) for c in chrom]
        )
    elif partition_count > 0:
        chunk_intervals = _derive_chunk_locus_intervals(vds_filtered)
        logger.info("Chunk locus intervals: %s", chunk_intervals)
        ref_ht = hl.filter_intervals(ref_ht, chunk_intervals)

    ref_ht = ref_ht.key_by("locus").select().distinct()
    ref_ht = hl.filter_intervals(
        ref_ht,
        telomeres_and_centromeres.ht().interval.collect(),
        keep=False,
    )
    if partition_count > 0 and repartition_factor > 1:
        target_n = partition_count * repartition_factor
        logger.info(
            "Repartitioning ref_ht into %d sub-partitions (%d * %d).",
            target_n,
            partition_count,
            repartition_factor,
        )
        ref_ht = ref_ht.repartition(target_n, shuffle=True)
    return ref_ht.checkpoint(hl.utils.new_temp_file("ref", "ht"))


def _run_coverage_chunk(args: argparse.Namespace) -> None:
    """
    Compute one chunk of the coverage/AN HT and write it to ``args.chunk_output``.

    Invoked by the per-chunk Hail Batch relay container (see
    ``_orchestrate_coverage_batch``). Hail must already be initialized
    (``main`` does this via ``_init_hail`` before dispatching here).

    Steps:
      1. Filter VDS via ``filter_partitions=range(chunk_start, chunk_stop)``.
      2. Derive locus intervals from the filtered VDS reference_data and
         ``filter_intervals`` ``vep_context`` to that range.
      3. Repartition ``ref_ht`` into ``partitions_per_chunk * repartition_factor``
         sub-partitions to shrink per-task densify size.
      4. Run ``compute_all_release_stats_per_ref_site``.
      5. Intermediate-checkpoint, ``naive_coalesce`` to ``partitions_per_chunk``,
         and write to ``args.chunk_output``.
    """
    project = args.project_name
    environment = args.environment
    partition_range, start, stop = _resolve_partition_range(args)
    n = stop - start

    test = args.test_2_partitions or args.test_chr22_chrx_chry or args.test
    chrom = ["chr22", "chrX", "chrY"] if args.test_chr22_chrx_chry else None

    group_membership_ht_path = group_membership(
        test=test, data_set=project, environment=environment
    ).path
    if args.reduce_min_aggs:
        group_membership_ht_path = (
            group_membership_ht_path.rstrip("/").removesuffix(".ht") + "_reduce.ht"
        )

    if project == "aou":
        sex_karyotype_field = "meta.sex_karyotype"
        meta_ht = hl.read_table(meta(data_type="genomes", environment=environment).path)
        vds = get_aou_vds(
            release_only=True,
            filter_partitions=partition_range,
            annotate_meta=True,
            chrom=chrom,
            environment=environment,
        )
        # AoU v8 VDS lacks DP; annotate adj without DP and synthesize DP from
        # LAD so `compute_stats_per_ref_site` doesn't error out.
        vmt = vds.variant_data
        vmt = annotate_adj_no_dp(vmt)
        vmt = vmt.annotate_entries(DP=hl.sum(vmt.LAD))
        vds = hl.vds.VariantDataset(vds.reference_data, vmt)
        if test and args.reduce_samples_for_test:
            meta_ht = meta_ht.filter(
                (meta_ht.project_meta.project == project) & (meta_ht.release)
            ).select()
            meta_ht = meta_ht.sample(0.001)
            vds = hl.vds.filter_samples(vds, meta_ht)
    else:
        sex_karyotype_field = "meta.sex_imputation.sex_karyotype"
        vds = get_gnomad_v5_genomes_vds(
            release_only=True,
            consent_drop_only=True,
            filter_partitions=partition_range,
            annotate_meta=True,
            chrom=chrom,
        )

    ref_ht = _build_chunk_ref_ht(
        vds_filtered=vds,
        project=project,
        repartition_factor=args.repartition_factor,
        partition_count=n,
        test_chr22_chrx_chry=args.test_chr22_chrx_chry,
        chrom=chrom,
    )

    validate_vds(vds)

    cov_and_an_ht = compute_all_release_stats_per_ref_site(
        vds,
        ref_ht,
        sex_karyotype_field=sex_karyotype_field,
        project=project,
        group_membership_ht=hl.read_table(group_membership_ht_path),
        reduce_min_aggs=args.reduce_min_aggs,
    )
    cov_and_an_ht = cov_and_an_ht.checkpoint(
        new_temp_file(f"{project}_cov_and_an_chunk", "ht")
    )
    target_n = max(n, 1)
    cov_and_an_ht = cov_and_an_ht.naive_coalesce(target_n)
    cov_and_an_ht.write(args.chunk_output, overwrite=True)
    logger.info("Wrote chunk [%d, %d) to %s", start, stop, args.chunk_output)


def _run_coverage_merge(input_paths: List[str], output_path: str) -> None:
    """
    Union per-chunk coverage HTs and write the merged HT to ``output_path``.

    Globals (``strata_meta``, etc.) are identical across all chunks (computed
    from the same group_membership HT), so the union inherits them from the
    first input. Used for both group-level and final merges in the
    ``--use-batch-fanout`` pipeline.
    """
    logger.info("Merging %d HTs -> %s", len(input_paths), output_path)
    hts = [hl.read_table(p) for p in input_paths]
    merged = hl.Table.union(*hts) if len(hts) > 1 else hts[0]
    merged.write(output_path, overwrite=True)
    logger.info("Wrote merged HT to %s", output_path)


def _init_hail_local_spark(
    log_name: str,
    jvm_heap: str,
    local_cores: int,
    gcs_requester_pays_project: str = "broad-mpg-gnomad",
) -> None:
    """
    Init Hail with a local-Spark backend (no separate QoB driver pod).

    Used when ``--chunk-backend local``: the chunk's container runs Hail's
    JVM in-process using `local[N]` Spark threads instead of spawning a
    separate QoB driver+worker pair on Hail Batch. Cheaper at scale (one
    spot container per chunk vs. nonpreempt driver + spot workers) but the
    container itself must be sized to hold the densify peak.

    :param log_name: Per-invocation log name (per-chunk so concurrent chunks
        don't clobber a shared GCS log).
    :param jvm_heap: Spark driver heap, e.g. ``"22g"``. Should be roughly
        ``container_memory_gb - 4`` to leave headroom for Python + OS.
    :param local_cores: Number of local Spark threads. Typically half the
        container CPU count to balance executor parallelism vs. driver work.
    :param gcs_requester_pays_project: Project to bill for requester-pays
        GCS reads (the AoU VDS bucket).
    """
    import pyspark

    conf = pyspark.SparkConf()
    conf.set("spark.driver.memory", jvm_heap)
    # Container-local log file; Hail's `copy_log` at script exit (or the
    # surrounding Hail Batch container's log capture) is the canonical
    # output for failure debugging.
    log_path = f"/tmp/hail-{log_name}.log"
    hl.init(
        master=f"local[{local_cores}]",
        log=log_path,
        tmp_dir=f"gs://fc-11093c2b-590e-424a-91ac-0cc040d562fc/batch-tmp-4day/hail/tmp",
        gcs_requester_pays_configuration=gcs_requester_pays_project,
        spark_conf={
            "spark.driver.memory": jvm_heap,
        },
        skip_logging_configuration=True,
        default_reference="GRCh38",
    )


def _build_setup_command(commit: str, methods_branch: str = "main") -> str:
    """
    Build shell commands to download gnomad_qc and gnomad_methods.

    Both repos are actively developed, so chunk containers pull them at
    runtime rather than relying on what's baked into the Docker image. The
    image provides ``hail`` and system dependencies. Mirrors the equivalent
    helper in ``generate_frequency.py`` from earlier this branch.

    :param commit: Git commit hash to pin gnomad_qc to.
    :param methods_branch: Branch/commit of gnomad_methods to pull.
    :return: Shell command string (terminated with newline).
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = (
        "https://github.com/broadinstitute/gnomad_methods/archive/"
        f"{methods_branch}.tar.gz"
    )
    methods_dir_suffix = methods_branch.replace("/", "-")
    # Hail config so `hl.init(backend="batch")` finds the Batch billing
    # project and remote_tmpdir, plus the GCS requester-pays project.
    # The canonical path is the XDG-style `~/.config/hail/config.ini`
    # (what `hailtop.config.get_user_config_path()` returns); also write
    # `~/.hail/config.ini` as a legacy fallback for older Hail versions.
    config_body = (
        "[batch]\n"
        "billing_project = gnomad-production\n"
        "remote_tmpdir = gs://fc-11093c2b-590e-424a-91ac-0cc040d562fc/batch-tmp\n"
        "[gcs_requester_pays]\n"
        "project = broad-mpg-gnomad\n"
    )
    return (
        "set -euxo pipefail\n"
        "mkdir -p ~/.config/hail ~/.hail\n"
        "cat > ~/.config/hail/config.ini <<'HAILCFG'\n"
        f"{config_body}"
        "HAILCFG\n"
        "cp ~/.config/hail/config.ini ~/.hail/config.ini\n"
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_dir_suffix} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def _build_chunk_common_flags(args: argparse.Namespace) -> str:
    """Build the CLI flag string shared by every per-chunk relay invocation."""
    flags = [
        f"--project-name {args.project_name}",
        "--environment batch",
        f"--gcp-billing-project {args.gcp_billing_project}",
        f"--tmp-dir-days {args.tmp_dir_days}",
        f"--n-partitions {args.n_partitions}",
        f"--repartition-factor {args.repartition_factor}",
        f"--chunk-backend {args.chunk_backend}",
    ]
    if args.chunk_backend == "local":
        # Pass through container sizing so the chunk worker can configure
        # its local Spark heap / threadcount to match what the orchestrator
        # actually allocated. These default to None and are resolved
        # lazily in main() if not explicitly set.
        if args.jvm_heap is not None:
            flags.append(f"--jvm-heap {args.jvm_heap}")
        if args.local_cores is not None:
            flags.append(f"--local-cores {args.local_cores}")
    if args.experimental:
        # Pass through only if user opted in. With --experimental the inner
        # QoB driver attaches to this outer batch (HAIL_BATCH_ID); without
        # it, each chunk's QoB creates its own Hail Batch.
        flags.append("--experimental")
    if args.reduce_min_aggs:
        flags.append("--reduce-min-aggs")
    if args.reduce_samples_for_test:
        flags.append("--reduce-samples-for-test")
    if args.cov_and_an_output_suffix:
        flags.append(f"--cov-and-an-output-suffix {args.cov_and_an_output_suffix}")
    if args.test_chr22_chrx_chry:
        flags.append("--test-chr22-chrx-chry")
    # Promote --test-2-partitions to --test for the chunk: --test-2-partitions
    # would override --chunk-start/--chunk-stop in `_resolve_partition_range`,
    # but we need test-mode resource paths + (optional) sample subsampling to
    # behave the same as a direct `--test` run.
    if args.test or args.test_2_partitions:
        flags.append("--test")
    if args.app_name:
        flags.append(f"--app-name {args.app_name}")
    if args.qob_driver_cores is not None:
        flags.append(f"--driver-cores {args.qob_driver_cores}")
    if args.qob_driver_memory:
        flags.append(f"--driver-memory {args.qob_driver_memory}")
    if args.qob_worker_cores is not None:
        flags.append(f"--worker-cores {args.qob_worker_cores}")
    if args.qob_worker_memory:
        flags.append(f"--worker-memory {args.qob_worker_memory}")
    return " ".join(flags)


def _submit_chunk_wave(
    args: argparse.Namespace,
    backend: hb.ServiceBackend,
    wave_idx: Optional[int],
    chunk_indices: List[int],
    cov_and_an_ht_path: str,
    setup_cmd: str,
    common_flags_str: str,
    script: str,
) -> None:
    """
    Build and submit a single Hail Batch for one wave of chunk jobs.

    Chunks-only: each chunk job is a relay container that runs
    ``--run-chunk`` (which spawns its own QoB driver). No merge jobs are
    added here; merging is a separate, user-triggered step
    (see ``--merge-cov-chunks`` once wired). ``batch.run()`` blocks until
    the wave completes when ``wait=True`` (the default for ServiceBackend).

    Chunks whose ``_SUCCESS`` already exists are not added to the pipeline.
    """
    project = args.project_name
    total = args.total_partitions
    pp = args.partitions_per_chunk
    regions = args.regions or ["us-central1"]
    wave_label = "single" if wave_idx is None else f"wave{wave_idx:03d}"
    batch_name = (
        f"v5_cov_{project}_{total}p_{pp}ppc_repart{args.repartition_factor}"
        f"_{wave_label}"
    )
    if args.cov_and_an_output_suffix:
        batch_name += f"_{args.cov_and_an_output_suffix}"

    batch = hb.Batch(name=batch_name, backend=backend)

    chunk_paths = [
        _chunk_path(cov_and_an_ht_path, "chunk", idx) for idx in chunk_indices
    ]

    n_added = 0
    for local_pos, idx in enumerate(chunk_indices):
        path = chunk_paths[local_pos]
        if not args.overwrite and _file_exists_for_env(path, "batch"):
            logger.info("  skip chunk %d (exists at %s)", idx, path)
            continue
        start = idx * pp
        stop = min(start + pp, total)
        j = batch.new_job(name=f"cov_chunk_{idx:06d}_{start}_{stop}")
        j.image(args.batch_image)
        j.cpu(args.chunk_cpu)
        j.memory(args.chunk_memory)
        j.storage(args.chunk_storage)
        j.regions(regions)
        j.spot(True)
        j.n_max_attempts(args.chunk_attempts)
        j.command(
            f"{setup_cmd}"
            f"{script} --run-chunk"
            f" --chunk-start {start} --chunk-stop {stop}"
            f" --chunk-output {path}"
            f" {common_flags_str}"
        )
        n_added += 1

    logger.info(
        "Submitting Hail Batch '%s': %d chunk jobs (dry_run=%s)",
        batch_name,
        n_added,
        args.batch_dry_run,
    )
    if n_added == 0 and not args.batch_dry_run:
        logger.info("  no pending chunks; skipping batch.run() for %s", batch_name)
        return
    batch.run(dry_run=args.batch_dry_run)


def _orchestrate_coverage_batch(
    args: argparse.Namespace, cov_and_an_ht_path: str
) -> None:
    """
    Fan per-partition coverage/AN compute out as relay chunk jobs in Hail Batch.

    Chunks-only: this orchestrator submits one relay job per partition chunk
    (each spawning its own QoB driver via ``hl.init(backend="batch")``) and
    nothing else. The merge step is intentionally separate and runs only
    when triggered by its own CLI flag, so we can validate per-chunk writes
    end-to-end before paying for any merge tree.

    By default each chunk's QoB creates a separate Hail Batch; pass
    ``--experimental`` to instead attach the inner QoB to this outer batch
    (one batch graph for relay + driver + workers).

    Two dispatch modes:
      * ``--chunks-per-wave`` unset: one Hail Batch contains every chunk job.
      * ``--chunks-per-wave N`` set: chunks are split into waves of N. Each
        wave is its own Hail Batch; the orchestrator blocks on each wave
        (``batch.run`` is sync when ``wait=True``, the default) before
        submitting the next.

    Idempotency: chunks whose ``_SUCCESS`` already exists are skipped.
    """
    project = args.project_name
    total = args.total_partitions
    pp = args.partitions_per_chunk
    if total <= 0 or pp <= 0:
        raise ValueError(
            "--use-batch-fanout requires --total-partitions and"
            " --partitions-per-chunk to be > 0."
        )

    n_chunks = (total + pp - 1) // pp
    all_chunk_indices = list(range(n_chunks))

    n_pending = 0
    for idx in all_chunk_indices:
        path = _chunk_path(cov_and_an_ht_path, "chunk", idx)
        if not args.overwrite and _file_exists_for_env(path, "batch"):
            logger.info("Skipping already-complete chunk %d at %s", idx, path)
        else:
            n_pending += 1
    logger.info(
        "Coverage fan-out: %d chunks total, %d pending, %d skipped"
        " (overwrite=%s, chunks_per_wave=%s, project=%s)",
        n_chunks,
        n_pending,
        n_chunks - n_pending,
        args.overwrite,
        args.chunks_per_wave,
        project,
    )

    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    setup_cmd = _build_setup_command(commit, methods_branch=args.methods_branch)

    backend_kwargs = {"billing_project": args.batch_billing_project}
    if args.batch_remote_tmpdir:
        backend_kwargs["remote_tmpdir"] = args.batch_remote_tmpdir
    backend = hb.ServiceBackend(**backend_kwargs)

    common_flags_str = _build_chunk_common_flags(args)
    script = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/compute_coverage.py"

    cpw = args.chunks_per_wave
    if cpw is None or cpw <= 0 or cpw >= n_chunks:
        # Single-batch path: one Hail Batch covers all chunks.
        _submit_chunk_wave(
            args=args,
            backend=backend,
            wave_idx=None,
            chunk_indices=all_chunk_indices,
            cov_and_an_ht_path=cov_and_an_ht_path,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
        )
        return

    # Wave path: split chunks into windows of cpw and submit one Hail Batch
    # per wave, blocking on each before submitting the next.
    waves = [all_chunk_indices[i : i + cpw] for i in range(0, n_chunks, cpw)]
    for wave_idx, wave_chunks in enumerate(waves):
        logger.info(
            "Submitting wave %d/%d (%d chunks: %d..%d)",
            wave_idx + 1,
            len(waves),
            len(wave_chunks),
            wave_chunks[0],
            wave_chunks[-1],
        )
        _submit_chunk_wave(
            args=args,
            backend=backend,
            wave_idx=wave_idx,
            chunk_indices=wave_chunks,
            cov_and_an_ht_path=cov_and_an_ht_path,
            setup_cmd=setup_cmd,
            common_flags_str=common_flags_str,
            script=script,
        )


def main(args):
    """Compute all sites coverage, allele number, and quality histograms for v5 genomes (AoU v8 + gnomAD v4)."""
    project = args.project_name
    environment = args.environment

    # Resolve backend-aware container defaults (ignore explicit overrides
    # from the user). For QoB the chunk container is just a tiny relay; for
    # local Spark the container has to fit the densify peak.
    if args.chunk_backend == "local":
        if args.chunk_cpu is None:
            args.chunk_cpu = 4
        if args.chunk_memory is None:
            args.chunk_memory = "highmem"
        if args.chunk_storage is None:
            args.chunk_storage = "25Gi"
    else:  # qob
        if args.chunk_cpu is None:
            args.chunk_cpu = 0.5
        if args.chunk_memory is None:
            args.chunk_memory = "standard"
        if args.chunk_storage is None:
            args.chunk_storage = "5Gi"

    # Auto-derive jvm_heap and local_cores from container shape if running
    # the local-Spark backend and the user didn't set them explicitly.
    if args.chunk_backend == "local":
        mem_per_core_gb = {"highmem": 6.5, "standard": 3.75, "lowmem": 0.9}.get(
            args.chunk_memory, 3.75
        )
        total_mem_gb = int(args.chunk_cpu * mem_per_core_gb)
        if args.jvm_heap is None:
            args.jvm_heap = f"{max(total_mem_gb - 4, 4)}g"
        if args.local_cores is None:
            args.local_cores = max(int(args.chunk_cpu) // 2, 1)

    # Orchestrator-only mode: build and submit a Hail Batch pipeline, then
    # exit. Skip _init_hail entirely so we don't pay for an unused QoB
    # driver here (the per-chunk relays each spawn their own).
    if args.use_batch_fanout:
        # `--test-2-partitions` in fanout mode is an alias for
        # `--total-partitions 2 --partitions-per-chunk 2` (one chunk covering
        # the first 2 VDS partitions). Override here so the orchestrator
        # dispatches one chunk instead of the default 145k-partition fan-out.
        if args.test_2_partitions:
            if args.total_partitions != 145000:
                logger.warning(
                    "--test-2-partitions overrides --total-partitions=%s -> 2",
                    args.total_partitions,
                )
            if args.partitions_per_chunk != 2:
                logger.warning(
                    "--test-2-partitions overrides --partitions-per-chunk=%s -> 2",
                    args.partitions_per_chunk,
                )
            args.total_partitions = 2
            args.partitions_per_chunk = 2

        test = args.test_2_partitions or args.test_chr22_chrx_chry or args.test
        cov_and_an_ht_path = coverage_and_an_path(
            test=test,
            data_set=project,
            environment=environment,
        ).path
        if args.cov_and_an_output_suffix:
            cov_and_an_ht_path = (
                cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
                + f"_{args.cov_and_an_output_suffix}.ht"
            )
        _orchestrate_coverage_batch(args, cov_and_an_ht_path)
        return

    # Per-worker log names so concurrent chunks don't clobber a single
    # `v5_coverage_and_an_generation.log` in the logging bucket. The
    # monolithic flow keeps the original name; chunk/merge workers get
    # invocation-specific names derived from the partition slice or output
    # path so each worker writes to its own log file.
    log_name = "v5_coverage_and_an_generation"
    if args.run_chunk:
        log_name = f"v5_cov_chunk_{args.chunk_start:08d}_{args.chunk_stop:08d}"
    elif args.run_merge:
        merge_slug = "merge"
        if args.merge_output:
            merge_slug = (
                args.merge_output.rstrip("/").split("/")[-1].removesuffix(".ht")
            )
        log_name = f"v5_cov_{merge_slug}"

    # Local-Spark chunk worker takes a different init path: no QoB driver
    # pod, just an in-container Hail JVM. Skip _init_hail entirely.
    if args.run_chunk and args.chunk_backend == "local":
        _init_hail_local_spark(
            log_name=log_name,
            jvm_heap=args.jvm_heap,
            local_cores=args.local_cores,
            gcs_requester_pays_project=args.gcp_billing_project,
        )
        _run_coverage_chunk(args)
        return

    # When --experimental is set, attach the QoB driver to an existing
    # Hail Batch via the HAIL_BATCH_ID env var (via the experimental init
    # path). Useful when this script is invoked from inside a
    # hailtop.batch pipeline that wants all QoB jobs added to its own
    # batch graph. Without --experimental, the env var is ignored.
    batch_id = None
    if args.experimental:
        batch_id_env = getenv("HAIL_BATCH_ID")
        batch_id = int(batch_id_env) if batch_id_env else None

    _init_hail(
        log_name,
        environment,
        billing_project=getattr(args, "gcp_billing_project", None),
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{qc_temp_prefix(environment=environment, days=args.tmp_dir_days)}coverage_and_an_generation",
        experimental=args.experimental,
        batch_id=batch_id,
        **_get_batch_resource_kwargs(args),
    )

    # Per-chunk worker entry-points invoked by the Hail Batch fan-out
    # pipeline. Both run after _init_hail so they share the same QoB init
    # and tmp_dir conventions; they exit before touching the existing
    # monolithic flow below.
    if args.run_chunk:
        _run_coverage_chunk(args)
        return
    if args.run_merge:
        _run_coverage_merge(args.merge_inputs, args.merge_output)
        return

    test_2_partitions = args.test_2_partitions
    test_chr22_chrx_chry = args.test_chr22_chrx_chry
    test = test_2_partitions or test_chr22_chrx_chry or args.test
    overwrite = args.overwrite
    n_partitions = args.n_partitions

    chrom = None
    if test_chr22_chrx_chry:
        chrom = ["chr22", "chrX", "chrY"]

    try:
        cov_and_an_ht_path = coverage_and_an_path(
            test=test,
            data_set=project,
            environment=environment,
        ).path
        if args.cov_and_an_output_suffix:
            cov_and_an_ht_path = cov_and_an_ht_path.rstrip("/").removesuffix(".ht")
            cov_and_an_ht_path = (
                f"{cov_and_an_ht_path}_{args.cov_and_an_output_suffix}.ht"
            )
        downsampling_ht_path = get_aou_downsampling(
            test=test, environment=environment
        ).path
        meta_ht_path = meta(data_type="genomes", environment=environment).path
        group_membership_ht_path = group_membership(
            test=test, data_set=project, environment=environment
        ).path
        if args.reduce_min_aggs:
            group_membership_ht_path = (
                group_membership_ht_path.rstrip("/").removesuffix(".ht") + "_reduce.ht"
            )
        meta_ht = None
        if project == "aou":
            meta_ht = hl.read_table(meta_ht_path)

        if args.write_aou_downsampling_ht:
            check_resource_existence(
                output_step_resources={"downsampling_ht": [downsampling_ht_path]},
                overwrite=overwrite,
            )

            ht = meta_ht.filter(
                (meta_ht.project_meta.project == project) & (meta_ht.release)
            )
            ds_ht = get_downsampling_ht(ht)
            ds_ht.write(downsampling_ht_path, overwrite=overwrite)

        if args.write_group_membership_ht:
            check_resource_existence(
                output_step_resources={
                    "group_membership_ht": [group_membership_ht_path]
                },
                overwrite=overwrite,
            )

            if project == "gnomad":
                logger.info("Writing gnomAD group membership HT...")
                v4_meta_ht_path = v4_meta(data_type="genomes").path

                group_membership_ht = get_group_membership_ht(
                    hl.read_table(v4_meta_ht_path),
                    project=project,
                    reduce_min_aggs=args.reduce_min_aggs,
                )
                group_membership_ht.write(group_membership_ht_path, overwrite=overwrite)
            else:
                logger.info("Writing AoU group membership HT...")
                meta_ht = meta_ht.filter(
                    (meta_ht.project_meta.project == project) & (meta_ht.release)
                )
                group_membership_ht = get_group_membership_ht(
                    meta_ht=meta_ht,
                    project=project,
                    ds_ht=hl.read_table(downsampling_ht_path),
                    reduce_min_aggs=args.reduce_min_aggs,
                )
                group_membership_ht.write(
                    group_membership_ht_path,
                    overwrite=overwrite,
                )

        if args.compute_all_cov_release_stats_ht:
            logger.info(
                "Computing coverage, all sites allele number, and optionally quality histograms HT for %s...",
                project,
            )
            check_resource_existence(
                output_step_resources={
                    "coverage_and_an_ht": [cov_and_an_ht_path],
                },
                overwrite=overwrite,
            )

            # Context Table is used because it contains every locus in the GRCh38
            # reference as opposed to a ref-blocked VDS reference dataset.
            ref_ht = vep_context.versions["105"].ht()
            if test_chr22_chrx_chry:
                ref_ht = hl.filter_intervals(
                    ref_ht, [hl.parse_locus_interval(c) for c in chrom]
                )
            elif test_2_partitions:
                ref_ht = ref_ht._filter_partitions(range(2))

            # Retain only 'locus' annotation from context Table.
            ref_ht = ref_ht.key_by("locus").select().distinct()

            # The AoU dataset should have removed centromeres and telomeres, but
            # this serves as an additional safeguard.
            ref_ht = hl.filter_intervals(
                ref_ht,
                telomeres_and_centromeres.ht().interval.collect(),
                keep=False,
            )
            ref_ht = ref_ht.checkpoint(hl.utils.new_temp_file("ref", "ht"))

            if project == "aou":
                sex_karyotype_field = "meta.sex_karyotype"
                vds = get_aou_vds(
                    release_only=True,
                    filter_partitions=range(2) if test_2_partitions else None,
                    annotate_meta=True,
                    chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
                    environment=environment,
                )
                # NOTE: AoU v8 VDS does not have DP annotation, so using custom function
                # to annotate adj.
                # Also adding DP here to make sure `compute_stats_per_ref_site`
                # doesn't throw an error.
                vmt = vds.variant_data
                vmt = annotate_adj_no_dp(vmt)
                vmt = vmt.annotate_entries(DP=hl.sum(vmt.LAD))
                vds = hl.vds.VariantDataset(vds.reference_data, vmt)

                if test and args.reduce_samples_for_test:
                    meta_ht = meta_ht.filter(
                        (meta_ht.project_meta.project == project) & (meta_ht.release)
                    ).select()
                    meta_ht = meta_ht.sample(0.001)
                    vds = hl.vds.filter_samples(vds, meta_ht)
            else:
                sex_karyotype_field = "meta.sex_imputation.sex_karyotype"
                vds = get_gnomad_v5_genomes_vds(
                    release_only=True,
                    consent_drop_only=True,
                    filter_partitions=range(2) if test_2_partitions else None,
                    annotate_meta=True,
                    chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
                )

            validate_vds(vds)

            cov_and_an_ht = compute_all_release_stats_per_ref_site(
                vds,
                ref_ht,
                sex_karyotype_field=sex_karyotype_field,
                project=project,
                group_membership_ht=hl.read_table(group_membership_ht_path),
                reduce_min_aggs=args.reduce_min_aggs,
            )
            cov_and_an_ht = cov_and_an_ht.checkpoint(
                new_temp_file(f"{project}_cov_and_an", "ht")
            )
            # Write out the intermediate HT.
            cov_and_an_ht = cov_and_an_ht.naive_coalesce(n_partitions)
            cov_and_an_ht.write(cov_and_an_ht_path, overwrite=overwrite)

        if args.merge_gnomad_coverage:
            merged_gnomad_coverage_ht_path = f"{qc_temp_prefix(environment=environment)}gnomad_v5_genomes_coverage.ht"
            check_resource_existence(
                output_step_resources={
                    "merged_gnomad_coverage_ht": [merged_gnomad_coverage_ht_path],
                },
                overwrite=overwrite,
            )

            gnomad_ht = hl.read_table(cov_and_an_ht_path).drop("AN")
            gnomad_release_ht = hl.read_table(
                release_coverage_path(
                    release_version=v4_COVERAGE_RELEASE,
                    public=True,
                )
            )
            if test:
                gnomad_release_ht = _filter_to_locus_bounds(
                    gnomad_release_ht, gnomad_ht
                )

            ht = merge_gnomad_coverage_hts(gnomad_ht, gnomad_release_ht)
            ht.write(merged_gnomad_coverage_ht_path, overwrite=overwrite)

        if args.merge_gnomad_an:
            merged_gnomad_an_ht_path = (
                f"{qc_temp_prefix(environment=environment)}gnomad_v5_genomes_an.ht"
            )
            check_resource_existence(
                output_step_resources={
                    "merged_gnomad_an_ht": [merged_gnomad_an_ht_path],
                },
                overwrite=overwrite,
            )

            gnomad_ht = hl.read_table(cov_and_an_ht_path).select("AN")
            gnomad_release_ht = hl.read_table(
                release_coverage_path(
                    release_version=v4_AN_RELEASE,
                    public=True,
                    coverage_type="allele_number",
                )
            )

            if test:
                gnomad_release_ht = _filter_to_locus_bounds(
                    gnomad_release_ht, gnomad_ht
                )

            ht = merge_gnomad_an_hts(gnomad_ht, gnomad_release_ht)
            ht.write(merged_gnomad_an_ht_path, overwrite=overwrite)

        if args.export_coverage_release_files:
            cov_ht_path = release_coverage_path(
                public=False,
                test=test,
                coverage_type="coverage",
                environment=environment,
            )
            cov_tsv_path = release_coverage_tsv_path(test=test, environment=environment)
            gnomad_coverage_ht_path = f"{qc_temp_prefix(environment=environment)}gnomad_v5_genomes_coverage.ht"
            check_resource_existence(
                input_step_resources={
                    "gnomad_coverage_ht": [gnomad_coverage_ht_path],
                },
                output_step_resources={
                    "cov_release_ht": [cov_ht_path],
                    "cov_tsv": [cov_tsv_path],
                },
                overwrite=overwrite,
            )

            logger.info("Exporting coverage HT and TSV...")
            aou_ht = hl.read_table(cov_and_an_ht_path)
            aou_ht = aou_ht.drop("AN", "qual_hists")
            gnomad_ht = hl.read_table(gnomad_coverage_ht_path)

            ht = join_aou_and_gnomad_coverage_ht(aou_ht, gnomad_ht)
            ht = ht.checkpoint(new_temp_file("aou_and_gnomad_cov_join", "ht"))
            ht = ht.naive_coalesce(n_partitions)
            ht = ht.checkpoint(cov_ht_path, overwrite=overwrite)
            ht.export(cov_tsv_path)

        if args.export_an_release_files:
            an_ht_path = release_coverage_path(
                public=False,
                test=test,
                coverage_type="allele_number",
                environment=environment,
            )
            an_tsv_path = release_all_sites_an_tsv_path(
                test=test, environment=environment
            )
            gnomad_an_ht_path = (
                f"{qc_temp_prefix(environment=environment)}gnomad_v5_genomes_an.ht"
            )
            check_resource_existence(
                input_step_resources={
                    "gnomad_an_ht": [gnomad_an_ht_path],
                },
                output_step_resources={
                    "an_release_ht": [an_ht_path],
                    "an_release_tsv": [an_tsv_path],
                },
                overwrite=overwrite,
            )

            logger.info("Exporting AN HT and TSV...")
            aou_ht = hl.read_table(cov_and_an_ht_path)
            aou_ht = aou_ht.select("AN")
            gnomad_ht = hl.read_table(gnomad_an_ht_path)

            ht = join_aou_and_gnomad_an_ht(aou_ht, gnomad_ht)
            ht = ht.checkpoint(new_temp_file("aou_and_gnomad_an_join", "ht"))
            ht = ht.select("AN")
            ht = ht.select_globals(
                strata_meta=ht.strata_meta,
                strata_sample_count=ht.strata_sample_count,
            )
            ht = ht.naive_coalesce(n_partitions)
            ht = ht.checkpoint(an_ht_path, overwrite=overwrite)

            ht = ht.transmute(AN=ht.AN[0])
            ht.export(an_tsv_path)

        if args.merge_qual_hists:
            qual_hists_path = qual_hists(test=test, environment=environment).path
            if args.qual_hists_output_suffix:
                qual_hists_path = qual_hists_path.rstrip("/").removesuffix(".ht")
                qual_hists_path = (
                    f"{qual_hists_path}_{args.qual_hists_output_suffix}.ht"
                )
            check_resource_existence(
                output_step_resources={"qual_hists_ht": [qual_hists_path]},
                overwrite=overwrite,
            )

            logger.info("Merging qual hists HTs...")
            aou_ht = (
                hl.read_table(cov_and_an_ht_path).select("qual_hists").select_globals()
            )
            gnomad_ht = (
                public_release(data_type="genomes")
                .ht()
                .select("histograms")
                .select_globals()
            )
            if test:
                gnomad_ht = _filter_to_locus_bounds(gnomad_ht, aou_ht)

            # Drop age hists because they are handled in the frequency script
            # and re-key by locus. `distinct()` deduplicates the multi-row-per-
            # locus rows that result from rekeying a (locus, alleles)-keyed HT
            # to locus-only; the `_all` hists are locus-level and identical
            # across split rows, so distinct keeps a representative row.
            gnomad_ht = gnomad_ht.select(
                qual_hists=gnomad_ht.histograms.drop("age_hists")
            )
            gnomad_ht = gnomad_ht.key_by("locus").distinct()

            ht = join_aou_and_gnomad_qual_hists_ht(aou_ht, gnomad_ht)
            ht.write(qual_hists_path, overwrite=overwrite)

    finally:
        if environment != "batch":
            logger.info("Copying hail log to logging bucket...")
            hl.copy_log(
                get_logging_path(
                    "v5_coverage_and_an_generation",
                    environment=environment,
                    tmp_dir_days=args.tmp_dir_days,
                )
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-name",
        help="Project name.",
        default="aou",
        type=str,
        choices=["aou", "gnomad"],
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )

    # Environment configuration.
    env_group = parser.add_argument_group("environment configuration")
    env_group.add_argument(
        "--environment",
        help="Compute environment.",
        choices=["rwb", "batch", "dataproc"],
        default="rwb",
    )
    env_group.add_argument(
        "--tmp-dir-days",
        type=int,
        default=4,
        help="Number of days for temp directory retention. Default is 4.",
    )
    env_group.add_argument(
        "--gcp-billing-project",
        type=str,
        default="broad-mpg-gnomad",
        help="Google Cloud billing project for reading requester pays buckets.",
    )

    # Batch-specific configuration.
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
    )
    batch_group.add_argument(
        "--experimental",
        action="store_true",
        help=(
            "Route the QoB init through `hl.experimental.init` instead of"
            " `hl.init`. When set, also reads HAIL_BATCH_ID from the env"
            " and (if present) attaches the QoB driver to that existing"
            " Hail Batch instead of creating a new one. Without this flag,"
            " HAIL_BATCH_ID is ignored."
        ),
    )
    batch_group.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_group.add_argument(
        "--driver-cores",
        type=int,
        default=None,
        help="Number of cores for driver node.",
    )
    batch_group.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_group.add_argument(
        "--worker-cores",
        type=int,
        default=None,
        help="Number of cores for worker nodes.",
    )
    batch_group.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of partitions to use for the output Table.",
        type=int,
        default=5000,
    )
    parser.add_argument(
        "--cov-and-an-output-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix appended to the cov_and_an HT path (before the"
            " .ht extension). Use to write the output to a sibling location"
            " for A/B comparison. Applies to both writes (step 3) and reads"
            " (downstream steps), so pass the same suffix consistently."
        ),
    )
    parser.add_argument(
        "--qual-hists-output-suffix",
        type=str,
        default=None,
        help=(
            "Optional suffix appended to the qual_hists HT path (before the"
            " .ht extension) — analogous to --cov-and-an-output-suffix but"
            " for the --merge-qual-hists output. Use for A/B comparison."
        ),
    )
    parser.add_argument(
        "--reduce-min-aggs",
        action="store_true",
        help=(
            "Use the leaf-only `reduce_to_minimal_groups` optimization: build"
            " the group_membership_ht with `reduce_to_minimal_groups=True`"
            " and pass `reducible_aggs={'AN'}` to `compute_stats_per_ref_site`."
            " Must be set consistently for both --write-group-membership-ht"
            " and --compute-all-cov-release-stats-ht runs (the gmh's"
            " `freq_reduced` global determines whether reduction takes effect"
            " at compute time)."
        ),
    )
    parser.add_argument(
        "--reduce-samples-for-test",
        action="store_true",
        help=(
            "When `--test` (or `--test-2-partitions` / `--test-chr22-chrx-chry`)"
            " is set on AoU, additionally subsample the AoU sample set to"
            " ~0.1%% via `meta_ht.sample(0.001)`. Default off: a `--test` run"
            " uses all samples but the partition / chrom subset, so AN values"
            " are stable and comparable across runs. Useful when you want a"
            " cheap tiny-cohort sanity check rather than a real-scale slice."
        ),
    )

    test_group = parser.add_mutually_exclusive_group()
    test_group.add_argument(
        "--test-2-partitions",
        help=(
            "Whether to run a test using only the first 2 partitions of the VDS test"
            " dataset."
        ),
        action="store_true",
    )
    test_group.add_argument(
        "--test-chr22-chrx-chry",
        help=(
            "Whether to run a test using only the chr22, chrX, and chrY chromosomes of"
            " the VDS test dataset."
        ),
        action="store_true",
    )

    group_membership_args = parser.add_argument_group(
        "Get gnomAD genomes group membership HT.",
    )
    group_membership_args.add_argument(
        "--write-group-membership-ht",
        help="Write group membership HT.",
        action="store_true",
    )
    group_membership_args.add_argument(
        "--test",
        help="Write test group membership HT to test path.",
        action="store_true",
    )

    parser.add_argument(
        "--write-aou-downsampling-ht",
        help="Write v5 downsampling HT.",
        action="store_true",
    )
    parser.add_argument(
        "--compute-all-cov-release-stats-ht",
        help="Compute the all sites coverage, allele number, and quality histogram HT.",
        action="store_true",
    )

    coverage_args = parser.add_argument_group(
        "Compute coverage release stats HT.",
    )
    coverage_args.add_argument(
        "--merge-gnomad-coverage",
        help="Subtract consent drop samples from v4 release HT to create gnomAD v5 genomes coverage HT.",
        action="store_true",
    )
    coverage_args.add_argument(
        "--export-coverage-release-files",
        help="Join and export AoU + gnomAD v4 coverage release HT and TSV file.",
        action="store_true",
    )

    an_args = parser.add_argument_group(
        "Compute AN release stats HT.",
    )
    an_args.add_argument(
        "--merge-gnomad-an",
        help="Subtract consent drop samples from v4 release HT to create gnomAD v5 genomes AN HT.",
        action="store_true",
    )
    an_args.add_argument(
        "--export-an-release-files",
        help="Exports joint AoU + gnomAD v4 AN release HT and TSV file.",
        action="store_true",
    )

    parser.add_argument(
        "--merge-qual-hists",
        help="Merge variant quality histograms from AoU v8 and gnomAD v4 genomes.",
        action="store_true",
    )

    fanout_group = parser.add_argument_group(
        "Hail Batch fan-out for `--compute-all-cov-release-stats-ht`.",
        "Submit one Hail Batch job per partition chunk; each chunk is a tiny"
        " relay container that spawns its own QoB driver to do the densify."
        " Idempotent: chunks whose `_SUCCESS` exists are skipped on rerun.",
    )
    fanout_group.add_argument(
        "--use-batch-fanout",
        action="store_true",
        help=(
            "Build and submit a Hail Batch pipeline that fans out the"
            " coverage/AN compute by partition chunks (QoB-from-container"
            " per chunk + tree-reduce merge). Mutually exclusive with"
            " --compute-all-cov-release-stats-ht and the per-chunk worker"
            " flags."
        ),
    )
    fanout_group.add_argument(
        "--partitions-per-chunk",
        type=int,
        default=2,
        help=(
            "Number of VDS partitions per chunk job. Default 2 (matches the"
            " driver-comfort load validated in the variance-experiment"
            " runs)."
        ),
    )
    fanout_group.add_argument(
        "--repartition-factor",
        type=int,
        default=50,
        help=(
            "Per-chunk: explode the filtered ref_ht into"
            " `partitions_per_chunk * repartition_factor` sub-partitions"
            " before densify, to shrink per-task densify size. Default 50."
        ),
    )
    fanout_group.add_argument(
        "--total-partitions",
        type=int,
        default=145000,
        help=(
            "Total VDS partition count to scatter across. Default 145000"
            " (prod AoU VDS partition count). Used by --use-batch-fanout."
        ),
    )
    fanout_group.add_argument(
        "--merge-group-size",
        type=int,
        default=500,
        help="Chunk HTs per group-merge job. Default 500.",
    )
    fanout_group.add_argument(
        "--chunks-per-wave",
        type=int,
        default=None,
        help=(
            "If set, split chunk dispatch into sequential Hail Batch waves"
            " of this size (each its own batch + intra-wave merge). After"
            " all waves finish, a final cross-wave merge unions the"
            " per-wave HTs into the final output. Default unset: one Hail"
            " Batch contains every chunk + merge job (fewer moving parts;"
            " untested at 72k+ scale)."
        ),
    )
    fanout_group.add_argument(
        "--methods-branch",
        type=str,
        default="main",
        help="Branch/commit of gnomad_methods to clone in chunk containers.",
    )
    fanout_group.add_argument(
        "--batch-image",
        type=str,
        default=(
            "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:latest"
        ),
        help="Docker image for the chunk + merge BashJob containers.",
    )
    fanout_group.add_argument(
        "--batch-billing-project",
        type=str,
        default="gnomad-production",
        help="Hail Batch billing project for the orchestrator pipeline.",
    )
    fanout_group.add_argument(
        "--batch-remote-tmpdir",
        type=str,
        default=None,
        help=(
            "GCS scratch path for hb.ServiceBackend. Defaults to the value"
            " encoded in `_build_setup_command` config."
        ),
    )
    fanout_group.add_argument(
        "--regions",
        type=str,
        nargs="+",
        default=None,
        help="GCP regions for chunk + merge jobs. Default ['us-central1'].",
    )
    fanout_group.add_argument(
        "--chunk-backend",
        type=str,
        default="qob",
        choices=["qob", "local"],
        help=(
            "Chunk-level compute backend.\n"
            "  qob   (default): the chunk container is a tiny relay that"
            " spawns a separate QoB driver+workers pair via"
            " `hl.init(backend='batch')`. Per-task spot retries; pays a"
            " nonpreempt n1-highmem-8 driver tax per chunk.\n"
            "  local: the chunk container runs Hail's JVM in-process with"
            " `local[N]` Spark threads. No separate driver pod; the"
            " container is everything. Cheaper at prod scale (no driver"
            " tax) but per-chunk preemption loses the whole chunk's work."
        ),
    )
    fanout_group.add_argument(
        "--chunk-cpu",
        type=float,
        default=None,
        help=(
            "CPU per chunk container. Default auto-resolves by backend:"
            " 0.5 for qob (relay only), 4 for local (must fit densify"
            " peak)."
        ),
    )
    fanout_group.add_argument(
        "--chunk-memory",
        type=str,
        default=None,
        choices=["lowmem", "standard", "highmem"],
        help=(
            "Hail Batch memory preset per chunk container. Default"
            " auto-resolves by backend: `standard` for qob, `highmem` for"
            " local."
        ),
    )
    fanout_group.add_argument(
        "--chunk-storage",
        type=str,
        default=None,
        help=(
            "/io storage per chunk container. Default auto-resolves by"
            " backend: 5Gi for qob (relay writes nothing), 25Gi for local"
            " (Spark shuffle spills + checkpoints)."
        ),
    )
    fanout_group.add_argument(
        "--qob-driver-cores",
        type=int,
        default=None,
        help=(
            "Cores for the QoB driver pod each chunk relay spawns"
            " (--chunk-backend=qob). Forwarded to the relay's"
            " --driver-cores; decoupled from the orchestrator/merge-side"
            " --driver-cores. Ignored when --chunk-backend=local."
        ),
    )
    fanout_group.add_argument(
        "--qob-driver-memory",
        type=str,
        default=None,
        help=(
            "Memory profile for the QoB driver pod each chunk relay"
            " spawns (--chunk-backend=qob), e.g. 'highmem'. Forwarded to"
            " the relay's --driver-memory; decoupled from the"
            " orchestrator/merge-side --driver-memory. Ignored when"
            " --chunk-backend=local."
        ),
    )
    fanout_group.add_argument(
        "--qob-worker-cores",
        type=int,
        default=None,
        help=(
            "Cores per QoB worker pod the chunk's QoB driver dispatches"
            " (--chunk-backend=qob). Forwarded to the relay's"
            " --worker-cores; decoupled from the orchestrator/merge-side"
            " --worker-cores. Ignored when --chunk-backend=local."
        ),
    )
    fanout_group.add_argument(
        "--qob-worker-memory",
        type=str,
        default=None,
        help=(
            "Memory profile per QoB worker pod the chunk's QoB driver"
            " dispatches (--chunk-backend=qob), e.g. 'lowmem'. Forwarded"
            " to the relay's --worker-memory; decoupled from the"
            " orchestrator/merge-side --worker-memory. Ignored when"
            " --chunk-backend=local."
        ),
    )
    fanout_group.add_argument(
        "--jvm-heap",
        type=str,
        default=None,
        help=(
            "Spark driver heap for `--chunk-backend local` (e.g. `22g`)."
            " Default auto-resolves to `container_memory_gb - 4` to leave"
            " headroom for Python + OS. Ignored when --chunk-backend=qob."
        ),
    )
    fanout_group.add_argument(
        "--local-cores",
        type=int,
        default=None,
        help=(
            "Number of local Spark threads inside `--chunk-backend local`"
            " containers. Default auto-resolves to `chunk_cpu // 2` to"
            " balance executor parallelism vs. driver work. Ignored when"
            " --chunk-backend=qob."
        ),
    )
    fanout_group.add_argument(
        "--chunk-attempts",
        type=int,
        default=5,
        help="Max retry attempts (n_max_attempts) per chunk job.",
    )
    fanout_group.add_argument(
        "--merge-cpu",
        type=int,
        default=4,
        help="CPU per merge container.",
    )
    fanout_group.add_argument(
        "--merge-memory",
        type=str,
        default="standard",
        help="Hail Batch memory preset per merge container.",
    )
    fanout_group.add_argument(
        "--merge-storage",
        type=str,
        default="50Gi",
        help="Extra /io storage per group-merge container.",
    )
    fanout_group.add_argument(
        "--final-merge-storage",
        type=str,
        default="100Gi",
        help="Extra /io storage for the final merge container.",
    )
    fanout_group.add_argument(
        "--batch-dry-run",
        action="store_true",
        help="If set, validate the Hail Batch DAG without submitting.",
    )

    worker_group = parser.add_argument_group(
        "Per-chunk worker entry-points for --use-batch-fanout (set by the"
        " orchestrator; not for direct use unless debugging a single chunk).",
    )
    worker_group.add_argument(
        "--run-chunk",
        action="store_true",
        help=(
            "Worker entry-point: compute one chunk and write to"
            " --chunk-output. Requires --chunk-start, --chunk-stop,"
            " --chunk-output."
        ),
    )
    worker_group.add_argument(
        "--chunk-start",
        type=int,
        default=0,
        help="First VDS partition index for --run-chunk (inclusive).",
    )
    worker_group.add_argument(
        "--chunk-stop",
        type=int,
        default=0,
        help="Last VDS partition index for --run-chunk (exclusive).",
    )
    worker_group.add_argument(
        "--chunk-output",
        type=str,
        default=None,
        help="GCS output path for --run-chunk's per-chunk HT.",
    )
    worker_group.add_argument(
        "--run-merge",
        action="store_true",
        help=(
            "Worker entry-point: union the HTs at --merge-inputs and write"
            " to --merge-output."
        ),
    )
    worker_group.add_argument(
        "--merge-inputs",
        type=str,
        nargs="+",
        default=None,
        help="GCS paths of per-chunk (or per-group) HTs to union.",
    )
    worker_group.add_argument(
        "--merge-output",
        type=str,
        default=None,
        help="GCS output path for --run-merge's merged HT.",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    batch_args = [
        "app_name",
        "driver_cores",
        "driver_memory",
        "qob_driver_cores",
        "qob_driver_memory",
        "worker_cores",
        "worker_memory",
        "qob_worker_cores",
        "qob_worker_memory",
    ]
    provided_batch_args = [arg for arg in batch_args if getattr(args, arg) is not None]
    if provided_batch_args and args.environment != "batch":
        parser.error(
            f"Batch configuration arguments ({', '.join('--' + a.replace('_', '-') for a in provided_batch_args)}) "
            f"require --environment=batch"
        )

    if args.project_name == "aou" and args.environment not in ("rwb", "batch"):
        parser.error(
            "--project-name=aou requires --environment to be 'rwb' or 'batch'."
        )
    if args.project_name == "gnomad" and args.environment != "dataproc":
        parser.error("--project-name=gnomad requires --environment=dataproc.")

    if args.write_aou_downsampling_ht and args.project_name != "aou":
        raise ValueError(
            "--write-aou-downsampling-ht requires --project-name to be 'aou'."
        )

    if args.merge_gnomad_coverage and args.project_name != "gnomad":
        raise ValueError(
            "--merge-gnomad-coverage requires --project-name to be 'gnomad'."
        )

    if args.merge_gnomad_an and args.project_name != "gnomad":
        raise ValueError("--merge-gnomad-an requires --project-name to be 'gnomad'.")

    if args.export_coverage_release_files and args.project_name != "aou":
        raise ValueError(
            "--export-coverage-release-files requires --project-name to be 'aou'."
        )

    if args.export_an_release_files and args.project_name != "aou":
        raise ValueError(
            "--export-an-release-files requires --project-name to be 'aou'."
        )

    if args.merge_qual_hists and args.project_name != "aou":
        raise ValueError("--merge-qual-hists requires --project-name to be 'aou'.")

    # --use-batch-fanout / --run-chunk / --run-merge are mutually exclusive
    # entry-points. The orchestrator sets --run-chunk / --run-merge inside
    # the per-chunk relay containers; users should not pass them directly
    # alongside --use-batch-fanout from the command line.
    fanout_modes = sum(
        bool(x) for x in (args.use_batch_fanout, args.run_chunk, args.run_merge)
    )
    if fanout_modes > 1:
        parser.error(
            "--use-batch-fanout, --run-chunk, and --run-merge are mutually"
            " exclusive."
        )
    if args.use_batch_fanout and args.compute_all_cov_release_stats_ht:
        parser.error(
            "--use-batch-fanout is an alternative to"
            " --compute-all-cov-release-stats-ht; do not pass both."
        )
    if args.run_chunk:
        if args.chunk_output is None:
            parser.error("--run-chunk requires --chunk-output.")
        if not args.test_2_partitions and args.chunk_stop <= args.chunk_start:
            parser.error(
                "--run-chunk requires --chunk-stop > --chunk-start (or"
                " --test-2-partitions as an alias for [0, 2))."
            )
    if args.run_merge:
        if args.merge_output is None or not args.merge_inputs:
            parser.error("--run-merge requires --merge-output and --merge-inputs.")
    if args.use_batch_fanout and args.environment != "batch":
        parser.error("--use-batch-fanout requires --environment=batch.")

    main(args)
