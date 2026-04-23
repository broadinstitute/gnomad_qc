"""
Script to generate frequency data for gnomAD v5.

This script calculates variant frequencies and histograms for:
1. gnomAD dataset - updating v4 frequencies by subtracting consent withdrawal samples
2. AoU dataset - using either pre-computed allele numbers or a densify approach

Processing Workflow:
--------------------
gnomAD (--process-gnomad):
1. Load v4 frequency table (contains frequencies and age histograms)
2. Prepare consent withdrawal VDS (split multiallelics, annotate metadata)
3. Calculate frequencies and age histograms for consent samples
4. Subtract from v4 frequencies to get updated gnomAD v5 frequencies

AoU (--process-aou):
  Densify path (default):
  1. Ensure split/prepared AoU VDS checkpoint exists (one-off prepare step)
  2. Fan out frequency computation via Hail Batch PythonJobs
     (one per partition slice, hierarchical merge)
  All-sites-AN path (--use-all-sites-ans):
  1. Load AoU VDS with metadata, single-job frequency computation

Merged dataset (--merge-datasets):
1. Merge frequency data and histograms from both gnomAD and AoU datasets.
2. Calculate FAF, grpmax, and other post-processing annotations on merged dataset.

Usage Examples:
---------------
# Process AoU (densify path, batch fan-out).
python generate_frequency.py --process-aou --environment batch

# Process AoU using all-sites ANs (single-job, no fan-out).
python generate_frequency.py --process-aou --use-all-sites-ans --environment rwb

# Process gnomAD consent withdrawals
python generate_frequency.py --process-gnomad --environment batch

# Run gnomAD in test mode
python generate_frequency.py --process-gnomad --test --test-partitions 2

# Merge both datasets
python generate_frequency.py --merge-datasets --environment batch --app-name "merged_freq" --driver-cores 8 --worker-memory highmem
"""

import argparse
import copy
import logging
from typing import List, Optional

import hail as hl
import hailtop.batch as hb
from gnomad.resources.grch38.gnomad import GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    agg_by_strata,
    bi_allelic_site_inbreeding_expr,
    compute_freq_by_strata,
    expand_strata_array_from_leaves,
    faf_expr,
    find_minimal_strata_groups,
    gen_anc_faf_max_expr,
    get_adj_expr,
    grpmax_expr,
    merge_freq_arrays,
    merge_histograms,
    qual_hist_expr,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.vcf import SORT_ORDER
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v5.annotations.annotation_utils import annotate_adj_no_dp
from gnomad_qc.v5.resources.annotations import (
    coverage_and_an_path,
    get_aou_freq_chunk_path,
    get_freq,
    get_split_aou_vds,
    group_membership,
)
from gnomad_qc.v5.resources.basics import (
    _file_exists_for_env,
    _get_batch_resource_kwargs,
    _init_hail,
    _init_hail_local_spark,
    get_aou_vds,
    get_gnomad_v5_genomes_vds,
    get_logging_path,
    qc_temp_prefix,
)
from gnomad_qc.v5.resources.meta import meta

# Use force=True so that our root handler wins over any handler that
# Hail / hailtop / absl / other deps installed during their imports above.
# Without force=True, basicConfig is a no-op when the root logger already has
# a handler, and all `logger.info` calls get silently dropped when running
# locally against the batch backend (Mode 1 QoB).
logging.basicConfig(
    format="%(levelname)s (%(name)s %(lineno)s): %(message)s",
    level=logging.INFO,
    force=True,
)
logger = logging.getLogger("v5_frequency")
logger.setLevel(logging.INFO)


def mt_hist_fields(mt: hl.MatrixTable) -> hl.StructExpression:
    """
    Annotate allele balance quality metrics histograms and age histograms onto MatrixTable.

    :param mt: Input MatrixTable.
    :return: Struct with allele balance, quality metrics histograms, and age histograms.
    """
    logger.info(
        "Computing quality metrics histograms and age histograms for each variant..."
    )
    qual_hists = qual_hist_expr(
        gt_expr=mt.GT,
        gq_expr=mt.GQ,
        dp_expr=hl.sum(mt.AD),
        adj_expr=mt.adj,
        ab_expr=(mt.AD[1] / hl.sum(mt.AD)),
        split_adj_and_raw=True,
    )
    return hl.struct(
        qual_hists=qual_hists,
        age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
    )


def _prepare_aou_vds(
    aou_vds: hl.vds.VariantDataset,
    use_all_sites_ans: bool = False,
    test: bool = False,
    environment: str = "batch",
) -> hl.vds.VariantDataset:
    """
    Prepare AoU VDS for frequency calculations.

    :param aou_vds: AoU VariantDataset.
    :param use_all_sites_ans: Whether to use all sites ANs for frequency calculations.
    :param test: Whether running in test mode.
    :param environment: Environment being used. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Prepared AoU VariantDataset.
    """
    aou_vmt = aou_vds.variant_data
    # Use existing AoU group membership table and filter to variant samples.
    logger.info(
        "Loading AoU group membership table for variant frequency stratification..."
    )
    group_membership_ht = group_membership(
        test=test, data_set="aou", environment=environment
    ).ht()

    logger.info("Selecting cols for frequency stratification...")
    aou_vmt = aou_vmt.select_cols(
        sex_karyotype=aou_vmt.meta.sex_karyotype,
        gen_anc=aou_vmt.meta.genetic_ancestry_inference.gen_anc,
        age=aou_vmt.meta.project_meta.age,
        group_membership=group_membership_ht[aou_vmt.col_key].group_membership,
    )
    # NOTE: At the time of writing this code, it was not yet decided if we would use
    # the all sites AN for the frequency calculcation, a cost-savings approach to
    # avoid a densify, or if we would densify within the frequency script to get
    # call stats. The split timing depends on whether we use all-sites ANs:
    # - With all-sites ANs: adjust ploidy -> adj -> split (preserves LGT for adj)
    # - Without: split -> adjust ploidy -> adj (consistent with v4 workflow)
    # Please see project notes for more details:
    # https://docs.google.com/document/d/109NuIQlZ3a8uw__iBMK7FpoO2bFFGECzEDnBpB1xW8E/edit?tab=t.0#heading=h.neezqe6t09ul
    if use_all_sites_ans:
        aou_vmt = aou_vmt.select_entries(
            LGT=adjusted_sex_ploidy_expr(
                aou_vmt.locus, aou_vmt.LGT, aou_vmt.sex_karyotype
            ),
            GQ=aou_vmt.GQ,
            LAD=aou_vmt.LAD,
            LA=aou_vmt.LA,
        )
        aou_vmt = annotate_adj_no_dp(aou_vmt)
        aou_vds = hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt)
        aou_vds = hl.vds.split_multi(aou_vds, filter_changed_loci=True)
        aou_vmt = aou_vds.variant_data
    else:
        aou_vds = hl.vds.VariantDataset(aou_vds.reference_data, aou_vmt)
        aou_vds = hl.vds.split_multi(aou_vds, filter_changed_loci=True)
        aou_vmt = aou_vds.variant_data
        # Defer sex-ploidy adjustment and adj computation to
        # `_calculate_aou_frequencies_and_hists_using_densify` so they run on
        # the dense MT and apply to ref-block fill-ins as well as variant
        # entries (matches the v4 convention in
        # gnomad_qc/v4/annotations/generate_freq.py:534-548). Doing it here
        # would be wasted work — the var-only `adj` field gets filled with
        # missing on ref-block samples by `to_dense_mt` and would have to be
        # recomputed anyway.
        aou_vmt = aou_vmt.select_entries(
            GT=aou_vmt.GT,
            GQ=aou_vmt.GQ,
            AD=aou_vmt.AD,
        )

    logger.info("Annotating globals...")
    group_membership_globals = group_membership_ht.index_globals()
    aou_vmt = aou_vmt.select_globals(
        freq_meta=group_membership_globals.freq_meta,
        freq_meta_sample_count=group_membership_globals.freq_meta_sample_count,
        age_distribution=aou_vmt.aggregate_cols(hl.agg.hist(aou_vmt.age, 30, 80, 10)),
        downsamplings=group_membership_globals.downsamplings,
    )

    # Drop GT from the reference data to work around a hail 0.2.130 bug in
    # `vds.to_dense_mt`: its `coalesce_join` builds shared_fields by prepending
    # the call_field and then re-adding any field present in both ref and var,
    # which crashes with "duplicate identifier 'GT'" when (as with the AoU VDS)
    # the reference data already carries GT. Fixed in 0.2.130.post1 (hail-is/hail
    # commit 9b8e322379-area; see PR #14695). The drop is also safe on fixed
    # versions: `_calculate_aou_frequencies_and_hists_using_densify` re-applies
    # `adjusted_sex_ploidy_expr` to the dense MT, which restores the same
    # sex-aware ploidy on ref-block fill-ins that the fixed `coalesce_join`
    # would have preserved from ref's GT.
    ref_data = aou_vds.reference_data
    if "GT" in ref_data.entry:
        ref_data = ref_data.drop("GT")
    aou_vds = hl.vds.VariantDataset(ref_data, aou_vmt)

    return aou_vds


def _calculate_aou_frequencies_and_hists_using_all_sites_ans(
    aou_variant_mt: hl.MatrixTable, test: bool = False, environment: str = "batch"
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using all sites ANs.

    :param aou_variant_mt: Prepared variant MatrixTable.
    :param test: Whether to use test resources.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    all_sites_an_ht = coverage_and_an_path(test=test, environment=environment).ht()
    aou_variant_mt = aou_variant_mt.annotate_rows(
        hist_fields=mt_hist_fields(aou_variant_mt)
    )

    logger.info("Annotating frequencies with all sites ANs...")
    group_membership_ht = group_membership(
        test=test, data_set="aou", environment=environment
    ).ht()
    aou_variant_freq_ht = agg_by_strata(
        aou_variant_mt.select_entries(
            "GT",
            "adj",
            n_alt_alleles=aou_variant_mt.GT.n_alt_alleles(),
            is_hom_var=aou_variant_mt.GT.is_hom_var(),
        ),
        {
            "AC": (lambda t: t.n_alt_alleles, hl.agg.sum),
            "homozygote_count": (lambda t: t.is_hom_var, hl.agg.count_where),
        },
        group_membership_ht=group_membership_ht,
        select_fields=["hist_fields"],
    )

    # Load AN values from all sites ANs table (calculated by another script but used
    # same group membership HT so same strata order).
    logger.info("Annotating AN values from all sites ANs...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        all_sites_an=all_sites_an_ht[aou_variant_freq_ht.locus].AN
    )

    logger.info("Building complete frequency struct with imported AN values...")
    aou_variant_freq_ht = aou_variant_freq_ht.annotate(
        freq=hl.map(
            lambda AC, hom_alt, AN: hl.struct(
                AC=hl.int32(AC),
                AF=hl.if_else(AN > 0, AC / AN, hl.missing(hl.tfloat64)),
                AN=hl.int32(AN),
                homozygote_count=hl.int32(hom_alt),
            ),
            aou_variant_freq_ht.AC,
            aou_variant_freq_ht.homozygote_count,
            aou_variant_freq_ht.all_sites_an,
        ),
    ).drop("all_sites_an")

    # Nest histograms to match gnomAD structure.
    # Note: hists_fields.qual_hists already contains raw_qual_hists and qual_hists
    # as nested fields due to split_adj_and_raw=True in qual_hist_expr.
    aou_variant_freq_ht = aou_variant_freq_ht.select(
        freq=aou_variant_freq_ht.freq,
        histograms=hl.struct(
            qual_hists=aou_variant_freq_ht.hist_fields.qual_hists.qual_hists,
            raw_qual_hists=aou_variant_freq_ht.hist_fields.qual_hists.raw_qual_hists,
            age_hists=aou_variant_freq_ht.hist_fields.age_hists,
        ),
    )

    return aou_variant_freq_ht


def _calculate_aou_frequencies_and_hists_using_densify(
    aou_vds: hl.vds.VariantDataset,
    reduce_to_minimal_groups: bool = False,
) -> hl.Table:
    """
    Calculate frequencies and age histograms for AoU variant data using densify.

    :param aou_vds: Prepared AoU VariantDataset.
    :param reduce_to_minimal_groups: When True, only compute call stats for the
        minimal "leaf" stratification groups and reconstruct the rest by
        element-wise summation. Reduces the cost of the per-variant
        aggregation proportionally to the reduction in group count.
        Default is False.
    :return: Table with freq and age_hists annotations.
    """
    logger.info("Annotating quality metrics histograms and age histograms...")
    aou_mt = hl.vds.to_dense_mt(aou_vds)
    # Apply sex-ploidy adjustment and (re)compute adj on the dense MT so they
    # cover both variant entries and ref-block fill-ins. This matches the
    # convention in gnomad_qc/v4/annotations/generate_freq.py:534-548.
    #
    # Why post-densify rather than on the variant data:
    #   - On hail 0.2.130, `to_dense_mt`'s buggy `coalesce_join` overwrites
    #     every ref-block GT with diploid `hl.call(0, 0)` regardless of
    #     karyotype; on 0.2.130.post1+ the ref's own GT is preserved.
    #     Adjusting here normalizes the two versions: chrX/chrY non-PAR
    #     hom-ref calls become haploid for XY, chrY calls become missing for
    #     XX, autosomes are unchanged.
    #   - `adj` is a variant-only entry field, so `to_dense_mt` fills it with
    #     missing on ref-block samples on every hail version. Without
    #     recomputing here, hom-ref calls inside ref blocks would be silently
    #     excluded from adj-filtered aggregations in `compute_freq_by_strata`
    #     (which uses `hl.agg.filter(adj[i], ...)`, and missing → False).
    #     `get_adj_expr` short-circuits via `~is_het()` for hom-ref calls and
    #     only gates on `gq_expr >= adj_gq`, so the missing AD on ref-block
    #     fill-ins is not a problem.
    aou_mt = aou_mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(aou_mt.locus, aou_mt.GT, aou_mt.sex_karyotype)
    )
    aou_mt = annotate_adj_no_dp(aou_mt)
    aou_mt = aou_mt.annotate_rows(hist_fields=mt_hist_fields(aou_mt))

    # Optionally reduce group_membership to only the leaf groups before the
    # expensive per-variant aggregation.
    if reduce_to_minimal_groups:
        freq_meta_full = hl.eval(aou_mt.freq_meta)
        freq_meta_sample_count_full = hl.eval(aou_mt.freq_meta_sample_count)
        leaf_indices, decomposition = find_minimal_strata_groups(freq_meta_full)
        n_full = len(freq_meta_full)

        logger.info(
            "Reducing freq_meta from %d to %d leaf groups for aggregation.",
            n_full,
            len(leaf_indices),
        )

        aou_mt = aou_mt.annotate_cols(
            group_membership=hl.array(
                [aou_mt.group_membership[i] for i in leaf_indices]
            )
        )
        aou_mt = aou_mt.annotate_globals(
            freq_meta=[freq_meta_full[i] for i in leaf_indices],
            freq_meta_sample_count=[
                freq_meta_sample_count_full[i] for i in leaf_indices
            ],
        )

    # hl.agg.call_stats is used within compute_freq_by_strata which returns int32s for
    # AC, AN, and homozygote_count so do not need to convert to int32 here.
    aou_freq_ht = compute_freq_by_strata(aou_mt, select_fields=["hist_fields"])

    # Expand the leaf-only freq array back to the full set of groups and
    # restore the original freq_meta globals.
    if reduce_to_minimal_groups:
        aou_freq_ht = aou_freq_ht.annotate(
            freq=expand_strata_array_from_leaves(
                aou_freq_ht.freq,
                leaf_indices,
                decomposition,
                n_full,
                is_freq_struct=True,
            )
        )
        aou_freq_ht = aou_freq_ht.annotate_globals(
            freq_meta=freq_meta_full,
            freq_meta_sample_count=freq_meta_sample_count_full,
        )

    # Nest histograms to match gnomAD structure.
    # Note: hist_fields.qual_hists already contains raw_qual_hists and qual_hists
    # as nested fields due to split_adj_and_raw=True in qual_hist_expr.
    aou_freq_ht = aou_freq_ht.transmute(
        histograms=hl.struct(
            qual_hists=aou_freq_ht.hist_fields.qual_hists.qual_hists,
            raw_qual_hists=aou_freq_ht.hist_fields.qual_hists.raw_qual_hists,
            age_hists=aou_freq_ht.hist_fields.age_hists,
        )
    )

    return aou_freq_ht


def process_aou_dataset(
    test_vds: bool = False,
    test_partitions: int = None,
    repartition_after_filter: Optional[int] = None,
    use_all_sites_ans: bool = False,
    environment: str = "batch",
    reduce_to_minimal_groups: bool = False,
    overwrite_split_aou_vds: bool = False,
) -> hl.Table:
    """
    Process All of Us dataset for frequency calculations and age histograms.

    This function efficiently processes the AoU VDS by:
    1. Computing complete frequency struct (uses imported AN from AoU all site ANs if requested)
    2. Generating age histograms within the frequency calculation

    :param test_vds: Whether to run in test mode on test VDS.
    :param test_partitions: Number of partitions to use in test mode. Default is None.
    :param repartition_after_filter: Optional number of partitions to repartition
        the filtered VDS into after applying ``filter_partitions``. Default is None.
    :param use_all_sites_ans: Whether to use all sites ANs for frequency calculations.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :param reduce_to_minimal_groups: When True, only compute call stats for
        the minimal "leaf" stratification groups and reconstruct the rest by
        summation. Reduces per-variant aggregation cost. Default is False.
    :param overwrite_split_aou_vds: Whether to overwrite the persistent split
        AoU VDS checkpoint (only used on the densify path). When False and
        the checkpoint already exists, it is read back instead of recomputed.
        This is independent of the general ``--overwrite`` flag that controls
        the final freq HT so you can rebuild the freq HT from an existing
        split VDS without paying to recompute it. Default is False.
    :return: Table with freq and age_hists annotations for AoU dataset.
    """
    test = test_vds or test_partitions is not None

    # On the densify path, check whether the persistent prepared/split AoU VDS
    # checkpoint already exists. If so, skip the raw load + _prepare_aou_vds
    # work entirely and read from the checkpoint. This survives multi-day runs
    # and script reruns (the tmp bucket used by `new_temp_file` expires every
    # few days, which is too short for a full AoU freq run).
    split_vds_resource = None
    use_cached_split_vds = False
    if not use_all_sites_ans:
        split_vds_resource = get_split_aou_vds(test=test, environment=environment)
        if (
            _file_exists_for_env(split_vds_resource.path, environment)
            and not overwrite_split_aou_vds
        ):
            use_cached_split_vds = True

    if use_cached_split_vds:
        logger.info(
            "Prepared split AoU VDS already exists at %s; reading from "
            "persistent checkpoint. Pass --overwrite-split-aou-vds to force "
            "recomputation.",
            split_vds_resource.path,
        )
        aou_vds = split_vds_resource.vds()
    else:
        aou_vds = get_aou_vds(
            annotate_meta=True,
            release_only=True,
            test=test_vds,
            filter_partitions=(
                list(range(test_partitions)) if test_partitions else None
            ),
            repartition_after_filter=repartition_after_filter,
            environment=environment,
        )
        aou_vds = _prepare_aou_vds(
            aou_vds,
            use_all_sites_ans=use_all_sites_ans,
            test=test,
            environment=environment,
        )

    # Calculate frequencies and age histograms together
    logger.info("Calculating AoU frequencies and age histograms...")
    if use_all_sites_ans:
        logger.info("Using all sites ANs for frequency calculations...")
        aou_freq_ht = _calculate_aou_frequencies_and_hists_using_all_sites_ans(
            aou_vds.variant_data, test=test, environment=environment
        )
    else:
        # Persist the prepared split VDS to a durable (non-tmp) resource so
        # `hl.vds.split_multi` and the pre-densify col/entry/ref transforms in
        # `_prepare_aou_vds` run exactly once. Without this checkpoint, every
        # downstream scan of the dense MT (hists + freq) would re-run
        # split_multi from the raw VDS. This is a sparse-VDS write (not a
        # dense MT), so the cost is bounded by the raw VDS size and is
        # one-time. Skipped when we already loaded the cached VDS above.
        if not use_cached_split_vds:
            logger.info(
                "Checkpointing prepared split AoU VDS to persistent resource "
                "at %s (one-time write; re-used on rerun unless "
                "--overwrite-split-aou-vds)...",
                split_vds_resource.path,
            )
            aou_vds = aou_vds.checkpoint(
                split_vds_resource.path, overwrite=overwrite_split_aou_vds
            )
        logger.info("Using densify for frequency calculations...")
        aou_freq_ht = _calculate_aou_frequencies_and_hists_using_densify(
            aou_vds, reduce_to_minimal_groups=reduce_to_minimal_groups
        )
    aou_freq_ht = select_final_dataset_fields(aou_freq_ht, dataset="aou")

    return aou_freq_ht


# ---------------------------------------------------------------------------
# Hail Batch fan-out for --process-aou
# ---------------------------------------------------------------------------


def _run_aou_freq_chunk(
    start: int,
    stop: int,
    output_path: str,
    test_vds: bool,
    test: bool,
    environment: str,
    reduce_to_minimal_groups: bool,
    repartition_after_filter: Optional[int],
) -> None:
    """
    Compute the AoU frequency HT for a single partition slice.

    Runs inside a Hail Batch PythonJob container. Loads the AoU VDS via
    ``get_aou_vds`` (which handles requester-pays and environment routing),
    filters to partitions ``[start, stop)``, prepares, densifies, and writes
    a partial frequency HT.

    :param start: First partition index (inclusive).
    :param stop: Last partition index (exclusive).
    :param output_path: GCS path for this chunk's partial freq HT.
    :param test_vds: Whether to use the test VDS (10-sample subset).
    :param test: Whether running in test mode (for group_membership lookup).
    :param environment: Resource environment string.
    :param reduce_to_minimal_groups: Forwarded to
        ``_calculate_aou_frequencies_and_hists_using_densify``.
    :param repartition_after_filter: Optional number of partitions to
        repartition this chunk's slice into for intra-container parallelism.
    """
    _init_hail_local_spark(f"v5_freq_aou_chunk_{start:06d}_{stop:06d}")

    vds = get_aou_vds(
        test=test_vds,
        filter_partitions=list(range(start, stop)),
        repartition_after_filter=repartition_after_filter,
        annotate_meta=True,
        release_only=True,
        environment=environment,
    )

    vds = _prepare_aou_vds(vds, test=test, environment=environment)

    freq_ht = _calculate_aou_frequencies_and_hists_using_densify(
        vds, reduce_to_minimal_groups=reduce_to_minimal_groups
    )
    freq_ht = select_final_dataset_fields(freq_ht, dataset="aou")
    freq_ht.write(output_path, overwrite=True)


def _run_aou_freq_merge(input_paths: List[str], output_path: str) -> None:
    """
    Union a list of partial AoU frequency HTs and write the result.

    Runs inside a Hail Batch PythonJob. Globals are identical across
    all inputs (produced from the same split VDS), so the union inherits
    them from the first HT.

    :param input_paths: GCS paths of the partial HTs to union.
    :param output_path: GCS path to write the merged HT.
    """
    _init_hail_local_spark("v5_freq_aou_merge")
    hts = [hl.read_table(p) for p in input_paths]
    merged = hl.Table.union(*hts) if len(hts) > 1 else hts[0]
    merged.write(output_path, overwrite=True)


def _build_setup_command(commit: str, methods_branch: str = "main") -> str:
    """
    Build shell commands to download gnomad_qc and gnomad_methods.

    Both repos are actively developed, so we pull them at runtime rather
    than relying on what's baked into the Docker image. The image provides
    ``hail`` and system dependencies (g++, curl).

    :param commit: Git commit hash to pin gnomad_qc to.
    :param methods_branch: Branch/commit of gnomad_methods to pull.
        Default is ``"main"``.
    :return: Shell command string.
    """
    qc_tarball = f"https://github.com/broadinstitute/gnomad_qc/archive/{commit}.tar.gz"
    methods_tarball = f"https://github.com/broadinstitute/gnomad_methods/archive/{methods_branch}.tar.gz"
    return (
        "set -euxo pipefail\n"
        f"curl -sSL {methods_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_methods-{methods_branch.replace('/', '-')} /tmp/gnomad_methods\n"
        f"curl -sSL {qc_tarball} | tar xz -C /tmp\n"
        f"mv /tmp/gnomad_qc-{commit} /tmp/gnomad_qc\n"
        "export PYTHONPATH=/tmp/gnomad_qc:/tmp/gnomad_methods:${PYTHONPATH:-}\n"
    )


def run_aou_freq_as_batch(
    final_output_path: str,
    n_partitions: int,
    test_vds: bool = False,
    test: bool = False,
    environment: str = "batch",
    partitions_per_job: int = 2,
    repartition_after_filter: Optional[int] = 100,
    merge_group_size: int = 500,
    reduce_to_minimal_groups: bool = False,
    chunk_cpu: int = 2,
    chunk_memory: str = "highmem",
    chunk_storage: str = "25Gi",
    merge_cpu: int = 4,
    merge_memory: str = "standard",
    merge_storage: str = "50Gi",
    final_merge_storage: str = "100Gi",
    regions: Optional[List[str]] = None,
    billing_project: str = "gnomad-production",
    remote_tmpdir: Optional[str] = None,
    image: str = "us-central1-docker.pkg.dev/broad-mpg-gnomad/images/v5_freq_batch:latest",
    methods_branch: str = "main",
    dry_run: bool = False,
) -> None:
    """
    Submit a Hail Batch that computes the AoU frequency HT via fan-out.

    Each chunk job clones ``gnomad_qc`` at the current commit, installs
    dependencies, then runs ``generate_frequency.py --run-chunk``. Chunks
    are merged hierarchically via ``--run-merge`` jobs.

    :param final_output_path: GCS path for the final merged AoU freq HT.
    :param n_partitions: Total number of VDS partitions to process.
    :param test_vds: Whether to use the test VDS (10-sample subset).
    :param test: Whether running in test mode.
    :param environment: Resource environment string.
    :param partitions_per_job: Partitions per chunk job. Default 2.
    :param repartition_after_filter: Partitions each chunk explodes into
        before densifying. Default 100.
    :param merge_group_size: Chunk HTs per group-merge job. Default 500.
    :param reduce_to_minimal_groups: Forwarded to chunk workers.
    :param chunk_cpu: CPU request per chunk job.
    :param chunk_memory: Hail Batch memory preset per chunk job.
    :param chunk_storage: Extra ``/io`` storage per chunk job.
    :param merge_cpu: CPU request per merge job.
    :param merge_memory: Memory preset for merge jobs.
    :param merge_storage: Extra storage per group-merge job.
    :param final_merge_storage: Extra storage for the final merge job.
    :param regions: GCP regions for jobs. Default ``["us-central1"]``.
    :param billing_project: Hail Batch billing project.
    :param remote_tmpdir: ``gs://`` scratch path for the ServiceBackend.
    :param image: Docker image for BashJobs.
    :param dry_run: If True, validate the DAG without submitting.
    """
    import subprocess

    if regions is None:
        regions = ["us-central1"]
    if remote_tmpdir is None:
        from gnomad_qc.v5.resources.constants import BATCH_TMP_BUCKET

        remote_tmpdir = f"gs://{BATCH_TMP_BUCKET}"
        logger.info("Using default remote_tmpdir: %s", remote_tmpdir)

    # Pin to the current commit so every job runs the same code.
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    setup_cmd = _build_setup_command(commit, methods_branch=methods_branch)

    # Base args shared by all chunk jobs.
    script = "python3 /tmp/gnomad_qc/gnomad_qc/v5/annotations/generate_frequency.py"
    chunk_base_flags = f"--environment {environment}"
    if test_vds:
        chunk_base_flags += " --test-vds"
    if reduce_to_minimal_groups:
        chunk_base_flags += " --reduce-to-minimal-groups"
    if repartition_after_filter:
        chunk_base_flags += f" --repartition-after-filter {repartition_after_filter}"

    logger.info(
        "Processing %s partitions; submitting %s-partition chunks (~%s jobs), "
        "pinned to commit %s...",
        n_partitions,
        partitions_per_job,
        (n_partitions + partitions_per_job - 1) // partitions_per_job,
        commit[:12],
    )

    backend = hb.ServiceBackend(
        billing_project=billing_project,
        remote_tmpdir=remote_tmpdir,
    )
    batch = hb.Batch(name="v5_freq_aou_fanout", backend=backend)

    def _configure(j, cpu, memory, storage):
        j.image(image)
        j.cpu(cpu)
        j.memory(memory)
        j.storage(storage)
        j.regions(regions)
        j.spot(True)
        j.n_max_attempts(5)
        return j

    # --- Chunk jobs ---
    chunk_jobs = []
    chunk_paths = []
    for chunk_idx, start in enumerate(range(0, n_partitions, partitions_per_job)):
        stop = min(start + partitions_per_job, n_partitions)
        chunk_path = get_aou_freq_chunk_path(
            chunk_idx, kind="chunk", test=test, environment=environment
        )
        j = batch.new_job(name=f"aou_freq_chunk_{chunk_idx:06d}")
        _configure(j, chunk_cpu, chunk_memory, chunk_storage)
        j.command(
            f"{setup_cmd}"
            f"{script} --run-chunk"
            f" --chunk-start {start} --chunk-stop {stop}"
            f" --chunk-output {chunk_path}"
            f" {chunk_base_flags}"
        )
        chunk_jobs.append(j)
        chunk_paths.append(chunk_path)

    # --- Group-merge jobs (hierarchical to avoid blowing up Spark plans) ---
    group_jobs = []
    group_paths = []
    for group_idx, group_start in enumerate(
        range(0, len(chunk_paths), merge_group_size)
    ):
        group_end = min(group_start + merge_group_size, len(chunk_paths))
        group_inputs = chunk_paths[group_start:group_end]
        group_path = get_aou_freq_chunk_path(
            group_idx, kind="group", test=test, environment=environment
        )
        j = batch.new_job(name=f"aou_freq_group_{group_idx:04d}")
        _configure(j, merge_cpu, merge_memory, merge_storage)
        j.depends_on(*chunk_jobs[group_start:group_end])
        paths_arg = " ".join(group_inputs)
        j.command(
            f"{setup_cmd}"
            f"{script} --run-merge"
            f" --merge-inputs {paths_arg}"
            f" --merge-output {group_path}"
        )
        group_jobs.append(j)
        group_paths.append(group_path)

    # --- Final merge ---
    final = batch.new_job(name="aou_freq_final_merge")
    _configure(final, merge_cpu, merge_memory, final_merge_storage)
    final.depends_on(*group_jobs)
    paths_arg = " ".join(group_paths)
    final.command(
        f"{setup_cmd}"
        f"{script} --run-merge"
        f" --merge-inputs {paths_arg}"
        f" --merge-output {final_output_path}"
    )

    logger.info(
        "Submitting Hail Batch: %s chunk jobs, %s group-merge jobs, "
        "1 final merge (dry_run=%s)...",
        len(chunk_jobs),
        len(group_jobs),
        dry_run,
    )
    batch.run(dry_run=dry_run)


def _prepare_consent_vds(
    v4_ht: hl.Table,
    test_vds: bool = False,
    test_partitions: int = 2,
) -> hl.vds.VariantDataset:
    """
    Load and prepare VDS for consent withdrawal sample processing.

    :param v4_ht: v4 release table for AF annotation.
    :param test: Whether running in test mode.
    :param test_partitions: Number of partitions to use in test mode. Default is 2.
    :return: Prepared VDS with consent samples, split multiallelics, and annotations.
    """
    logger.info("Loading and preparing VDS for consent withdrawal samples...")

    vds = get_gnomad_v5_genomes_vds(
        release_only=True,
        test=test_vds,
        consent_drop_only=True,
        annotate_meta=True,
        filter_partitions=list(range(test_partitions)),
    )

    logger.info(
        "VDS has been filtered to %s consent withdrawal samples...",
        vds.variant_data.count_cols(),
    )

    vmt = vds.variant_data
    vmt = vmt.select_cols(
        gen_anc=vmt.meta.population_inference.pop,
        sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
        age=vmt.meta.project_meta.age,
    )

    logger.info("Selecting entries and annotating non_ref hets pre-split...")
    vmt = vmt.select_entries(
        "LA", "LAD", "DP", "GQ", "LGT", _het_non_ref=vmt.LGT.is_het_non_ref()
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    vds = vds.checkpoint(new_temp_file("consent_samples_vds", "vds"))

    logger.info("Splitting multiallelics in gnomAD sample withdrawal VDS...")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    # Annotate with v4 frequencies for hom alt depletion fix
    vmt = vds.variant_data
    vmt = vmt.annotate_rows(v4_af=v4_ht[vmt.row_key].freq[0].AF)

    # This follows the v3/v4 genomes workflow for adj and sex adjusted genotypes which
    # were added before the hom alt depletion fix.
    # The correct order is to do the hom alt fix before adjusting sex ploidy and before
    # determining the adj annotation because haploid GTs have different adj filtering
    # criteria, but the option to adjust ploidy after adj is included for consistency
    # with v3.1, where we added the adj annotation before adjusting for sex ploidy.
    logger.info("Computing sex adjusted genotypes and quality annotations...")
    vmt = vmt.annotate_entries(
        adj=get_adj_expr(vmt.GT, vmt.GQ, vmt.DP, vmt.AD),
    )
    vmt = vmt.select_entries(
        "AD",
        "DP",
        "GQ",
        "_het_non_ref",
        "adj",
        GT=adjusted_sex_ploidy_expr(vmt.locus, vmt.GT, vmt.sex_karyotype),
    )

    # We set use_v3_1_correction to True to mimic the v4 genomes approach.
    logger.info("Applying v4 genomes hom alt depletion fix...")
    vmt = vmt.annotate_entries(
        GT=hom_alt_depletion_fix(
            vmt.GT,
            het_non_ref_expr=vmt._het_non_ref,
            af_expr=vmt.v4_af,
            ab_expr=vmt.AD[1] / vmt.DP,
            use_v3_1_correction=True,
        )
    )
    logger.info("Annotating age distribution...")
    vmt = vmt.annotate_globals(
        age_distribution=vmt.aggregate_cols(hl.agg.hist(vmt.age, 30, 80, 10))
    )

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)
    return vds.checkpoint(new_temp_file("consent_samples_vds_prepared", "vds"))


def _calculate_consent_frequencies_and_age_histograms(
    vds: hl.vds.VariantDataset,
    test_run: bool = False,
    environment: str = "batch",
) -> hl.Table:
    """
    Calculate frequencies and age histograms for consent withdrawal samples.

    :param vds: Prepared VDS with consent samples.
    :param test_run: Whether running in test mode.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Table with freq and age_hists annotations for consent samples.
    """
    logger.info("Densifying VDS for frequency calculations...")
    mt = hl.vds.to_dense_mt(vds)
    # Group membership table is already filtered to consent drop samples and is in GCS.
    group_membership_ht = group_membership(
        test=test_run, data_set="gnomad", environment=environment
    ).ht()

    mt = mt.annotate_cols(
        group_membership=group_membership_ht[mt.col_key].group_membership,
    )
    mt = mt.annotate_globals(
        freq_meta=group_membership_ht.index_globals().freq_meta,
        freq_meta_sample_count=group_membership_ht.index_globals().freq_meta_sample_count,
    )

    logger.info(
        "Calculating frequencies and age histograms using compute_freq_by_strata..."
    )

    # Annotate hists_fields on MatrixTable rows before calling compute_freq_by_strata so
    # to keep the age_hists annotation on the frequency table.
    logger.info("Annotating hists_fields on MatrixTable rows...")
    mt = mt.annotate_rows(
        hists_fields=hl.struct(
            age_hists=age_hists_expr(mt.adj, mt.GT, mt.age),
        )
    )

    logger.info(
        "Computing frequencies for consent samples using compute_freq_by_strata..."
    )
    consent_freq_ht = compute_freq_by_strata(
        mt,
        select_fields=["hists_fields"],
    )

    consent_freq_ht = consent_freq_ht.transmute(
        age_hists=consent_freq_ht.hists_fields.age_hists,
    )

    return consent_freq_ht.checkpoint(new_temp_file("consent_freq_and_hists", "ht"))


def _subtract_consent_frequencies_and_age_histograms(
    v4_ht: hl.Table,
    consent_freq_ht: hl.Table,
) -> hl.Table:
    """
    Subtract consent withdrawal frequencies and age histograms from v4 frequency table.

    :param v4_ht: v4 release table (contains both freq and histograms.age_hists).
    :param consent_freq_ht: Consent withdrawal table with freq and age_hists annotations.
    :return: Updated frequency table with consent frequencies and age histograms subtracted.
    """
    logger.info(
        "Subtracting consent withdrawal frequencies and age histograms from v4 release table..."
    )

    joined_freq_ht = v4_ht.annotate(
        consent_freq=consent_freq_ht[v4_ht.key].freq,
        consent_age_hists=consent_freq_ht[v4_ht.key].age_hists,
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        consent_freq_meta=consent_freq_ht.index_globals().freq_meta,
        consent_freq_meta_sample_count=consent_freq_ht.index_globals().freq_meta_sample_count,
    )

    logger.info("Subtracting consent frequencies...")
    updated_freq_expr, updated_freq_meta, updated_sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.consent_freq],
        [
            joined_freq_ht.index_globals().freq_meta,
            joined_freq_ht.index_globals().consent_freq_meta,
        ],
        operation="diff",
        count_arrays={
            "freq_meta_sample_count": [
                joined_freq_ht.index_globals().freq_meta_sample_count,
                joined_freq_ht.index_globals().consent_freq_meta_sample_count,
            ],
        },
    )
    # Update the frequency table with freq changes.
    joined_freq_ht = joined_freq_ht.annotate(freq=updated_freq_expr)
    joined_freq_ht = joined_freq_ht.annotate_globals(
        freq_meta=updated_freq_meta,
        freq_meta_sample_count=updated_sample_counts["freq_meta_sample_count"],
    )

    logger.info("Subtracting consent age histograms...")
    updated_age_hist_het = merge_histograms(
        [
            joined_freq_ht.histograms.age_hists.age_hist_het,
            joined_freq_ht.consent_age_hists.age_hist_het,
        ],
        operation="diff",
    )
    updated_age_hist_hom = merge_histograms(
        [
            joined_freq_ht.histograms.age_hists.age_hist_hom,
            joined_freq_ht.consent_age_hists.age_hist_hom,
        ],
        operation="diff",
    )

    # Update the frequency table with age hist changes.
    joined_freq_ht = joined_freq_ht.annotate(
        histograms=joined_freq_ht.histograms.annotate(
            age_hists=joined_freq_ht.histograms.age_hists.annotate(
                age_hist_het=updated_age_hist_het,
                age_hist_hom=updated_age_hist_hom,
            )
        ),
    )

    return joined_freq_ht.checkpoint(new_temp_file("merged_freq_and_hists", "ht"))


def select_final_dataset_fields(ht: hl.Table, dataset: str = "gnomad") -> hl.Table:
    """
    Create final dataset freq Table with only desired annotations.

    :param ht: Hail Table containing all annotations.
    :param dataset: Dataset to select final fields, either "gnomad", "aou" or "merged".
    :return: Hail Table with final annotations.
    """
    if dataset not in ["gnomad", "aou", "merged"]:
        raise ValueError(f"Invalid dataset: {dataset}")

    if dataset in ["gnomad", "aou"]:
        final_globals = ["freq_meta", "freq_meta_sample_count", "age_distribution"]
        final_fields = ["freq", "histograms"]

        if dataset == "aou":
            # AoU has one extra 'downsamplings' global field that is not present in
            # gnomAD.
            final_globals.append("downsamplings")

        # Convert all int64 annotations in the freq struct to int32s for merging type
        # compatibility.
        ht = ht.annotate(
            freq=ht.freq.map(
                lambda x: x.annotate(
                    **{k: hl.int32(v) for k, v in x.items() if v.dtype == hl.tint64}
                )
            )
        )
    else:
        sort_order = copy.deepcopy(SORT_ORDER) + ["aou-downsampling"]

        ht = ht.annotate_globals(
            freq_index_dict=make_freq_index_dict_from_meta(
                ht.freq_meta, sort_order=sort_order
            ),
            faf_index_dict=make_freq_index_dict_from_meta(
                ht.faf_meta, sort_order=sort_order
            ),
        )
        final_globals = [
            "freq_meta",
            "freq_index_dict",
            "freq_meta_sample_count",
            "faf_meta",
            "faf_index_dict",
            "age_distribution",
            "aou_downsamplings",
        ]
        final_fields = [
            "freq",
            "faf",
            "grpmax",
            "fafmax",
            "inbreeding_coeff",
            "histograms",
        ]

    return ht.select(*final_fields).select_globals(*final_globals)


def _fix_v4_global_age_distribution(
    freq_ht: hl.Table, environment: str = "batch"
) -> hl.Table:
    """
    Fix the age distribution global annotation in the frequency table.

    :param freq_ht: Frequency table to annotate with the age distribution.
    :param environment: Environment to use. Default is "batch".
    :return: Frequency table with the age distribution global annotation fixed.
    """
    # Use v5 meta as age is already fixed in the v5 project metadata as are the consent
    # withdrawal samples' releasable field.
    meta_ht = meta(environment=environment).ht()
    meta_ht = meta_ht.filter(
        (meta_ht.release) & (meta_ht.project_meta.project == "gnomad")
    )
    meta_ht = meta_ht.annotate_globals(
        age_distribution=meta_ht.aggregate(
            hl.agg.hist(meta_ht.project_meta.age, 30, 80, 10)
        )
    )
    freq_ht = freq_ht.annotate_globals(
        age_distribution=meta_ht.index_globals().age_distribution
    )

    return freq_ht


def _add_non_aou_subset_entries(freq_ht: hl.Table) -> hl.Table:
    """
    Add non-AoU subset fields to the frequency table.

    Duplicates gnomAD freq array fields (adj, raw, gen_anc, sex, gen_anc-sex) and adds
    "subset": "non_aou" to freq_meta for downstream non-AoU subset frequency reporting.

    :param freq_ht: Frequency table to add non-AoU subset fields to.
    :return: Frequency table with non-AoU subset fields added.
    """
    logger.info(
        "Filtering to non-AoU subset strata (adj, raw, gen_anc, sex, gen_anc-sex)..."
    )
    # Filter to only adj, raw, gen_anc, sex, and gen_anc-sex strata by excluding
    # entries with downsampling or subset keys.
    non_aou_freq_meta, non_aou_array_exprs = filter_arrays_by_meta(
        freq_ht.freq_meta,
        {
            "freq": freq_ht.freq,
            "freq_meta_sample_count": freq_ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=["downsampling", "subset"],
        keep=False,
        combine_operator="or",
    )

    logger.info("Adding non-aou subset data to freq and freq_meta...")
    non_aou_freq_meta = non_aou_freq_meta.map(
        lambda d: hl.dict(d.items().append(("subset", "non_aou")))
    )

    # Can extend freq and freq_meta_sample_count because there are no overlap of strata.
    freq_ht = freq_ht.annotate(
        freq=freq_ht.freq.extend(non_aou_array_exprs["freq"]),
    )
    freq_ht = freq_ht.annotate_globals(
        freq_meta=freq_ht.freq_meta.extend(non_aou_freq_meta),
        freq_meta_sample_count=freq_ht.freq_meta_sample_count.extend(
            non_aou_array_exprs["freq_meta_sample_count"]
        ),
    )

    return freq_ht


def process_gnomad_dataset(
    test_vds: bool = False,
    test_partitions: int = 2,
    environment: str = "batch",
) -> hl.Table:
    """
    Process gnomAD dataset to update v4 frequency HT by removing consent withdrawal samples.

    This function performs frequency adjustment by:
    1. Loading v4 frequency HT (contains both frequencies and age histograms)
    2. Loading consent withdrawal VDS
    3. Filtering to sites present in BOTH consent VDS AND v4 frequency table
    4. Calculating frequencies and age histograms for consent withdrawal samples
    5. Subtracting both frequencies and age histograms from v4 frequency HT
    6. Only overwriting fields that were actually updated in the final output

    :param test_vds: Whether to run on test vds. Default is False.
    :param test_partitions: Number of partitions to filter to in test mode. Default is 2.
    :param environment: Environment to use. Default is "batch". Must be one of "rwb"
        or "batch".
    :return: Updated frequency HT with updated frequencies and age histograms for gnomAD dataset.
    """
    test_run = test_vds or test_partitions is not None
    v4_ht = release_sites(data_type="genomes").ht()

    vds = _prepare_consent_vds(
        v4_ht,
        test=test_vds,
        test_partitions=test_partitions,
    )

    logger.info("Calculating frequencies and age histograms for consent samples...")
    consent_freq_ht = _calculate_consent_frequencies_and_age_histograms(
        vds, test_run=test_run, environment=environment
    )

    if test_run:
        v4_ht = v4_ht.filter(hl.is_defined(consent_freq_ht[v4_ht.key]))
        v4_ht = v4_ht.naive_coalesce(100).checkpoint(
            new_temp_file("v4_ht_filtered_test", "ht")
        )

    logger.info("Subtracting consent frequencies and age histograms from v4...")
    updated_freq_ht = _subtract_consent_frequencies_and_age_histograms(
        v4_ht, consent_freq_ht
    )

    logger.info("Merging updated frequency fields...")
    freq_ht = _merge_updated_frequency_fields(v4_ht, updated_freq_ht)

    logger.info("Reannotating gnomAD's age distribution global annotation...")
    freq_ht = _fix_v4_global_age_distribution(freq_ht, environment=environment)

    logger.info(
        "Duplicating gnomAD adj, raw, gen_anc, sex, and gen_anc-sex array entries with"
        " 'subset': 'non_aou' in freq_meta..."
    )
    freq_ht = _add_non_aou_subset_entries(freq_ht)

    # Select only the fields that were updated as FAF/grpmax/inbreeding_coeff annotations
    # will be calculated on the final merged dataset.
    logger.info("Selecting gnomAD freq HT final fields...")
    final_freq_ht = select_final_dataset_fields(freq_ht, dataset="gnomad")

    return final_freq_ht


def _merge_updated_frequency_fields(
    v4_release_ht: hl.Table, updated_freq_ht: hl.Table
) -> hl.Table:
    """
    Merge frequency tables, only overwriting fields that were actually updated.

    For sites that exist in updated_freq_ht, use the updated values.
    For sites that don't exist in updated_freq_ht, keep original values.

    Note: FAF/grpmax/inbreeding_coeff annotations are not calculated during consent
    withdrawal processing and will be calculated later on the final merged dataset.

    :param v4_release_ht: Original v4 release table.
    :param updated_freq_ht: Updated frequency table with consent withdrawals subtracted.
    :return: Final frequency table with selective field updates.
    """
    logger.info("Merging frequency tables with selective field updates...")

    # Bring in updated values with a single lookup.
    updated_row = updated_freq_ht[v4_release_ht.key]

    # Update freq and age_hists in a single annotate to avoid source mismatch:
    # - freq: use updated if present, otherwise keep original
    # - histograms.age_hists: update only age_hists, preserving qual_hists and raw_qual_hists
    final_freq_ht = v4_release_ht.annotate(
        freq=hl.coalesce(updated_row.freq, v4_release_ht.freq),
        histograms=v4_release_ht.histograms.annotate(
            age_hists=hl.coalesce(
                updated_row.histograms.age_hists,
                v4_release_ht.histograms.age_hists,
            )
        ),
    )

    # Update globals from updated table.
    updated_globals = {}
    for global_field in ["freq_meta", "freq_meta_sample_count"]:
        if global_field in updated_freq_ht.globals:
            updated_globals[global_field] = updated_freq_ht.index_globals()[
                global_field
            ]

    if updated_globals:
        final_freq_ht = final_freq_ht.annotate_globals(**updated_globals)

    return final_freq_ht


def merge_gnomad_and_aou_frequencies(
    gnomad_freq_ht: hl.Table,
    aou_freq_ht: hl.Table,
) -> hl.Table:
    """
    Merge frequency data and histograms from gnomAD and All of Us datasets.

    :param gnomad_freq_ht: Frequency Table for gnomAD.
    :param aou_freq_ht: Frequency Table for AoU.
    :return: Merged frequency Table with combined frequencies and histograms.
    """
    joined_freq_ht = gnomad_freq_ht.annotate(
        aou_freq=aou_freq_ht[gnomad_freq_ht.key].freq
    )

    joined_freq_ht = joined_freq_ht.annotate_globals(
        aou_freq_meta=aou_freq_ht.index_globals().freq_meta,
        aou_freq_meta_sample_count=aou_freq_ht.index_globals().freq_meta_sample_count,
        aou_age_distribution=aou_freq_ht.index_globals().age_distribution,
        aou_downsamplings=aou_freq_ht.index_globals().downsamplings,
    )

    merged_freq, merged_meta, sample_counts = merge_freq_arrays(
        [joined_freq_ht.freq, joined_freq_ht.aou_freq],
        [
            joined_freq_ht.index_globals().freq_meta,
            joined_freq_ht.index_globals().aou_freq_meta,
        ],
        operation="sum",
        count_arrays={
            "counts": [
                joined_freq_ht.index_globals().freq_meta_sample_count,
                joined_freq_ht.index_globals().aou_freq_meta_sample_count,
            ],
        },
    )

    joined_freq_ht = joined_freq_ht.annotate(freq=merged_freq).annotate_globals(
        freq_meta=merged_meta,
        freq_meta_sample_count=sample_counts["counts"],
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(merged_meta)),
    )

    # Rename the 'downsampling' group in freq meta list to 'aou_downsampling' as aou
    # is the source dataset for downsampling group. gnomAD downsamplings can
    # be retrieved from v3.
    logger.info(
        "Renaming 'downsampling' group in freq meta list to 'aou_downsampling'..."
    )
    renamed_freq_meta = hl.literal(
        [
            {("aou-downsampling" if k == "downsampling" else k): m[k] for k in m}
            for m in hl.eval(merged_meta)
        ]
    )
    joined_freq_ht = joined_freq_ht.annotate_globals(
        freq_meta=renamed_freq_meta,
    )

    logger.info("Merging quality histograms and age histograms from both datasets...")
    joined_freq_ht = joined_freq_ht.annotate(
        aou_histograms=aou_freq_ht[joined_freq_ht.key].histograms,
    )

    def _merge_hist_struct(hist1, hist2, operation="sum"):
        """Merge all fields of two histogram structs."""
        return hl.struct(
            **{
                field: merge_histograms(
                    [hist1[field], hist2[field]], operation=operation
                )
                for field in hist1.dtype.fields
            }
        )

    merged_histograms = hl.struct(
        qual_hists=_merge_hist_struct(
            joined_freq_ht.histograms.qual_hists,
            joined_freq_ht.aou_histograms.qual_hists,
        ),
        raw_qual_hists=_merge_hist_struct(
            joined_freq_ht.histograms.raw_qual_hists,
            joined_freq_ht.aou_histograms.raw_qual_hists,
        ),
        age_hists=_merge_hist_struct(
            joined_freq_ht.histograms.age_hists,
            joined_freq_ht.aou_histograms.age_hists,
        ),
    )

    joined_freq_ht = joined_freq_ht.annotate(histograms=merged_histograms)

    # Merge age_distribution global histograms (single histogram per dataset, not
    # per-strata like freq_meta_sample_count).
    merged_age_distribution = merge_histograms(
        [
            joined_freq_ht.index_globals().age_distribution,
            joined_freq_ht.index_globals().aou_age_distribution,
        ],
        operation="sum",
    )
    joined_freq_ht = joined_freq_ht.annotate_globals(
        age_distribution=merged_age_distribution
    )

    return joined_freq_ht


def calculate_faf_and_grpmax_annotations(
    ht: hl.Table,
) -> hl.Table:
    """
    Calculate FAF, grpmax, gen_anc_faf_max, and inbreeding coefficient annotations.

    Computes filtering allele frequencies and grpmax for both the full dataset
    (gnomad) and the non-AoU subset, similar to v4's non-UKB subset handling.

    :param ht: Merged frequency table for AoU and gnomAD data containing 'freq' and
        'freq_meta' annotations.
    :return: Table with 'faf', 'grpmax', 'gen_anc_faf_max', and 'inbreeding_coeff'.
    """
    logger.info(
        "Filtering frequencies to just 'non_aou' subset entries for 'faf' "
        "calculations..."
    )
    # Filter to non_aou subset and remove the "subset" key from freq_meta so faf_expr
    # pulls the correct indices.
    non_aou_freq_meta, non_aou_array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {"freq": ht.freq},
        items_to_filter={"subset": ["non_aou"]},
        keep=True,
        combine_operator="or",
    )
    non_aou_freq_meta = non_aou_freq_meta.map(
        lambda d: hl.dict(d.items().filter(lambda x: x[0] != "subset"))
    )
    non_aou_ht = ht.annotate(freq=non_aou_array_exprs["freq"])
    non_aou_ht = non_aou_ht.annotate_globals(freq_meta=non_aou_freq_meta)

    freq_metas = {
        "gnomad": (ht.freq, ht.index_globals().freq_meta),
        "non_aou": (
            non_aou_ht[ht.key].freq,
            non_aou_ht.index_globals().freq_meta,
        ),
    }

    faf_exprs = []
    faf_meta_exprs = []
    grpmax_exprs = {}
    gen_anc_faf_max_exprs = {}

    for dataset, (freq, meta) in freq_metas.items():
        faf, faf_meta = faf_expr(
            freq, meta, ht.locus, GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v5"]
        )
        grpmax = grpmax_expr(freq, meta, GEN_ANC_GROUPS_TO_REMOVE_FOR_GRPMAX["v5"])
        gen_anc_faf_max = gen_anc_faf_max_expr(faf, faf_meta)

        # Add subset back to non_aou faf meta.
        if dataset == "non_aou":
            faf_meta = [{**x, **{"subset": "non_aou"}} for x in faf_meta]

        faf_exprs.append(faf)
        faf_meta_exprs.append(faf_meta)
        grpmax_exprs[dataset] = grpmax
        gen_anc_faf_max_exprs[dataset] = gen_anc_faf_max

    logger.info(
        "Annotating 'faf', 'grpmax', 'gen_anc_faf_max', and 'inbreeding_coeff'..."
    )
    ht = ht.annotate(
        faf=hl.flatten(faf_exprs),
        grpmax=hl.struct(**grpmax_exprs),
        fafmax=hl.struct(**gen_anc_faf_max_exprs),
        inbreeding_coeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1]),
    )
    faf_meta_exprs = hl.flatten(faf_meta_exprs)
    ht = ht.annotate_globals(
        faf_meta=faf_meta_exprs,
        faf_index_dict=make_freq_index_dict_from_meta(faf_meta_exprs),
    )

    ht = ht.checkpoint(new_temp_file("freq_with_faf", "ht"))

    return ht


def _get_default_app_name(args) -> Optional[str]:
    """
    Derive a Hail Batch app name from the processing-step args.

    Combines the requested processing steps (``--process-gnomad``,
    ``--process-aou``, ``--merge-datasets``) into a single name so the run is
    easy to identify in the Hail Batch UI. Returns ``None`` when no processing
    step is selected.

    :param args: Parsed command-line arguments.
    :return: Default app name, or ``None`` if no steps were requested.
    """
    steps = []
    if args.process_gnomad:
        steps.append("gnomad")
    if args.process_aou:
        steps.append("aou")
    if args.merge_datasets:
        steps.append("merge")

    if not steps:
        return None

    return f"v5_freq_{'_'.join(steps)}"


def _initialize_hail(args) -> None:
    """
    Initialize Hail with appropriate configuration for the environment.

    :param args: Parsed command-line arguments.
    """
    # Auto-derive an app name from the requested processing steps if the user
    # did not pass one explicitly. Only relevant for the batch backend.
    if args.environment == "batch" and args.app_name is None:
        args.app_name = _get_default_app_name(args)

    _init_hail(
        "v5_frequency_generation",
        args.environment,
        billing_project=getattr(args, "gcp_billing_project", None),
        tmp_dir_days=args.tmp_dir_days,
        tmp_dir=f"{qc_temp_prefix(environment=args.environment, days=args.tmp_dir_days)}frequency_generation",
        **_get_batch_resource_kwargs(args),
    )


def main(args):
    """Generate v5 frequency data."""
    # --- Batch worker subcommands (run inside Hail Batch containers) ---
    if args.run_chunk:
        _run_aou_freq_chunk(
            start=args.chunk_start,
            stop=args.chunk_stop,
            output_path=args.chunk_output,
            test_vds=args.test_vds,
            test=args.test_vds or args.test_partitions is not None,
            environment=args.environment,
            reduce_to_minimal_groups=args.reduce_to_minimal_groups,
            repartition_after_filter=args.repartition_after_filter,
        )
        return

    if args.run_merge:
        _run_aou_freq_merge(
            input_paths=args.merge_inputs,
            output_path=args.merge_output,
        )
        return

    # --- Normal orchestrator flow ---
    environment = args.environment
    use_all_sites_ans = args.use_all_sites_ans
    test_vds = args.test_vds
    test_partitions = args.test_partitions
    overwrite = args.overwrite
    overwrite_split_aou_vds = args.overwrite_split_aou_vds
    aou_freq_suffix = args.aou_freq_ht_suffix
    tmp_dir_days = args.tmp_dir_days
    test_run = test_vds or test_partitions is not None

    _initialize_hail(args)

    try:
        logger.info("Running generate_frequency.py...")

        if args.process_gnomad:
            logger.info("Processing gnomAD dataset...")

            gnomad_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set="gnomad",
                environment=environment,
            )

            check_resource_existence(
                output_step_resources={"process-gnomad": [gnomad_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = process_gnomad_dataset(
                test=test_vds,
                test_partitions=test_partitions,
                environment=environment,
            )

            logger.info(
                "Writing gnomAD frequency HT (with embedded age histograms) to %s...",
                gnomad_freq.path,
            )
            gnomad_freq_ht.write(gnomad_freq.path, overwrite=overwrite)

        if args.process_aou:
            logger.info("Processing All of Us dataset...")
            aou_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set="aou",
                environment=environment,
                suffix=aou_freq_suffix,
            )

            check_resource_existence(
                output_step_resources={"process-aou": [aou_freq]},
                overwrite=overwrite,
            )

            if use_all_sites_ans:
                # All-sites-AN path: cheap single-job, no fan-out needed.
                aou_freq_ht = process_aou_dataset(
                    test_vds=test_vds,
                    test_partitions=test_partitions,
                    repartition_after_filter=args.repartition_after_filter,
                    use_all_sites_ans=True,
                    environment=environment,
                    reduce_to_minimal_groups=args.reduce_to_minimal_groups,
                    overwrite_split_aou_vds=overwrite_split_aou_vds,
                )
                logger.info("Writing AoU frequency HT to %s...", aou_freq.path)
                aou_freq_ht.write(aou_freq.path, overwrite=overwrite)
            else:
                # Densify path: fan out via Hail Batch PythonJobs. Each
                # chunk reads the raw VDS, filters to its partition slice,
                # prepares (split_multi etc.), densifies, and aggregates.
                if test_partitions:
                    n_partitions = test_partitions
                elif args.n_partitions:
                    n_partitions = args.n_partitions
                else:
                    raise ValueError(
                        "The densify fan-out path requires either "
                        "--test-partitions or --n-partitions to know how "
                        "many partitions to process."
                    )
                logger.info(
                    "Fanning out %s partitions via Hail Batch...",
                    n_partitions,
                )

                run_aou_freq_as_batch(
                    final_output_path=aou_freq.path,
                    n_partitions=n_partitions,
                    test_vds=test_vds,
                    test=test_run,
                    environment=environment,
                    partitions_per_job=args.partitions_per_job,
                    repartition_after_filter=args.repartition_after_filter,
                    merge_group_size=args.merge_group_size,
                    reduce_to_minimal_groups=args.reduce_to_minimal_groups,
                    chunk_cpu=args.chunk_cpu,
                    chunk_memory=args.chunk_memory,
                    chunk_storage=args.chunk_storage,
                    merge_cpu=args.merge_cpu,
                    merge_memory=args.merge_memory,
                    merge_storage=args.merge_storage,
                    final_merge_storage=args.final_merge_storage,
                    billing_project=args.billing_project,
                    remote_tmpdir=args.batch_remote_tmpdir,
                    **({"image": args.batch_image} if args.batch_image else {}),
                    methods_branch=args.methods_branch,
                    dry_run=args.batch_dry_run,
                )

        if args.merge_datasets:
            logger.info(
                "Merging frequency data and age histograms from both datasets..."
            )

            merged_freq = get_freq(
                test=test_run,
                data_type="genomes",
                data_set="merged",
                environment=environment,
            )

            check_resource_existence(
                output_step_resources={"merge-datasets": [merged_freq]},
                overwrite=overwrite,
            )

            gnomad_freq_ht = get_freq(
                data_type="genomes",
                test=test_run,
                data_set="gnomad",
                environment=environment,
            )
            aou_freq_ht = get_freq(
                data_type="genomes",
                test=test_run,
                data_set="aou",
                environment=environment,
                suffix=aou_freq_suffix,
            )

            check_resource_existence(
                input_step_resources={
                    "process-gnomad": [gnomad_freq_ht],
                    "process-aou": [aou_freq_ht],
                }
            )
            merged_freq_ht = merge_gnomad_and_aou_frequencies(
                gnomad_freq_ht.ht(),
                aou_freq_ht.ht(),
            )
            merged_freq_ht = merged_freq_ht.checkpoint(
                new_temp_file("merged_freq", "ht")
            )

            logger.info(
                "Calculating FAF, grpmax, and other annotations on merged dataset..."
            )
            merged_freq_ht = calculate_faf_and_grpmax_annotations(merged_freq_ht)

            merged_freq_ht = select_final_dataset_fields(
                merged_freq_ht, dataset="merged"
            )

            logger.info("Writing merged frequency HT to %s...", merged_freq.path)
            merged_freq_ht.write(merged_freq.path, overwrite=overwrite)

    finally:
        # Skip in batch mode: the JVM runs remotely so there is no local log
        # file to copy, and Hail Batch already retains the full driver/worker
        # logs (accessible via the Batch UI or `hailctl batch log`).
        if environment != "batch":
            hl.copy_log(
                get_logging_path(
                    "v5_frequency_run",
                    environment=environment,
                    tmp_dir_days=tmp_dir_days,
                )
            )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(
        description="Generate frequency data for gnomAD v5."
    )

    # General arguments.
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--overwrite-split-aou-vds",
        help=(
            "Overwrite the persistent split AoU VDS checkpoint (the prepared/"
            "split_multi VDS written by --process-aou on the densify path). "
            "Independent of --overwrite so you can rebuild the freq HT from "
            "an existing split VDS without paying to recompute it. Only "
            "relevant with --process-aou (without --use-all-sites-ans)."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--aou-freq-ht-suffix",
        type=str,
        default=None,
        help=(
            "Optional filename suffix inserted before the '.ht' extension on "
            "the AoU frequency HT, e.g. 'split_vds' -> "
            "'...frequencies.split_vds.ht'. Useful for tagging test / "
            "experimental runs so they don't collide with the default output "
            "path. Applied to both --process-aou and --merge-datasets reads "
            "of the AoU freq HT so the merge picks up the same file. Default "
            "is None (no suffix)."
        ),
    )

    # Test/debug arguments.
    test_group = parser.add_argument_group("testing options")
    test_group.add_argument(
        "--test-vds",
        help="Use the test VDS for given project.",
        action="store_true",
    )
    test_group.add_argument(
        "--test-partitions",
        type=int,
        default=2,
        help="If supplied, filter full or test VDS to N partitions in test mode. Default is 2.",
    )
    test_group.add_argument(
        "--repartition-after-filter",
        type=int,
        default=None,
        help=(
            "Number of partitions to repartition the VDS into after applying"
            " partition filtering. Useful when filtered partitions are too large"
            " for parallelism on preemptible VMs."
        ),
    )

    # Processing step arguments.
    processing_group = parser.add_argument_group("processing steps")
    processing_group.add_argument(
        "--process-gnomad",
        help="Process gnomAD dataset for frequency calculations.",
        action="store_true",
    )
    processing_group.add_argument(
        "--process-aou",
        help=(
            "Process AoU dataset for frequency calculations. Ensures the "
            "split/prepared VDS checkpoint exists, then fans out frequency "
            "computation via Hail Batch PythonJobs."
        ),
        action="store_true",
    )
    processing_group.add_argument(
        "--use-all-sites-ans",
        help="Use all sites ANs in frequency calculations to avoid a densify.",
        action="store_true",
    )
    processing_group.add_argument(
        "--reduce-to-minimal-groups",
        help=(
            "Only compute call stats for the minimal 'leaf' stratification "
            "groups and reconstruct the rest by summation. Reduces the cost "
            "of per-variant aggregation. Only applies to the densify path "
            "(--process-aou without --use-all-sites-ans)."
        ),
        action="store_true",
    )
    processing_group.add_argument(
        "--merge-datasets",
        help="Merge frequency data from both gnomAD and AoU datasets.",
        action="store_true",
    )

    # Batch worker subcommands (invoked by Hail Batch jobs, not by users).
    worker_group = parser.add_argument_group(
        "batch worker subcommands",
        "Internal subcommands run by Hail Batch chunk/merge jobs.",
    )
    worker_group.add_argument(
        "--run-chunk",
        help=argparse.SUPPRESS,
        action="store_true",
    )
    worker_group.add_argument(
        "--chunk-start",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--chunk-stop",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--chunk-output",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--run-merge",
        help=argparse.SUPPRESS,
        action="store_true",
    )
    worker_group.add_argument(
        "--merge-inputs",
        nargs="+",
        default=None,
        help=argparse.SUPPRESS,
    )
    worker_group.add_argument(
        "--merge-output",
        type=str,
        default=None,
        help=argparse.SUPPRESS,
    )

    # Environment configuration.
    env_group = parser.add_argument_group("environment configuration")
    env_group.add_argument(
        "--environment",
        help="Environment to run in.",
        choices=["rwb", "batch"],
        default="batch",
    )
    env_group.add_argument(
        "--tmp-dir-days",
        type=int,
        default=4,
        help="Number of days for temp directory retention. Default is 4.",
    )
    env_group.add_argument(
        "--billing-project",
        type=str,
        default="gnomad-production",
        help="Google Cloud billing project for reading requester pays buckets.",
    )

    # Batch-specific configuration.
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for batch/QoB backend (only used when --environment=batch).",
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

    # Batch fan-out configuration (used by --process-aou densify path).
    fanout_group = parser.add_argument_group(
        "batch fan-out configuration",
        "Parameters for --process-aou Hail Batch PythonJob fan-out (densify path).",
    )
    fanout_group.add_argument(
        "--n-partitions",
        type=int,
        default=None,
        help=(
            "Total number of VDS partitions to process. Required for "
            "production runs (the AoU VDS has ~145000). For test runs, "
            "--test-partitions is used instead."
        ),
    )
    fanout_group.add_argument(
        "--partitions-per-job",
        type=int,
        default=2,
        help=(
            "Partitions of the split VDS to process per chunk job. Default 2"
            " keeps each job well under the ~30min GCP preemption window."
        ),
    )
    fanout_group.add_argument(
        "--merge-group-size",
        type=int,
        default=500,
        help="Number of chunk HTs unioned per group-merge job. Default 500.",
    )
    fanout_group.add_argument(
        "--chunk-cpu",
        type=int,
        default=2,
        help="CPU request per chunk job. Default 2.",
    )
    fanout_group.add_argument(
        "--chunk-memory",
        type=str,
        default="highmem",
        help="Memory preset per chunk job. Default 'highmem'.",
    )
    fanout_group.add_argument(
        "--chunk-storage",
        type=str,
        default="25Gi",
        help="Extra /io storage per chunk job. Default '25Gi'.",
    )
    fanout_group.add_argument(
        "--merge-cpu",
        type=int,
        default=4,
        help="CPU request per merge job. Default 4.",
    )
    fanout_group.add_argument(
        "--merge-memory",
        type=str,
        default="standard",
        help="Memory preset for merge jobs. Default 'standard'.",
    )
    fanout_group.add_argument(
        "--merge-storage",
        type=str,
        default="50Gi",
        help="Extra storage per group-merge job. Default '50Gi'.",
    )
    fanout_group.add_argument(
        "--final-merge-storage",
        type=str,
        default="100Gi",
        help="Extra storage for the final merge job. Default '100Gi'.",
    )
    fanout_group.add_argument(
        "--batch-image",
        type=str,
        default=None,
        help=(
            "Docker image for chunk and merge jobs. Defaults to"
            " hailgenetics/hail matching the current Hail version."
        ),
    )
    fanout_group.add_argument(
        "--batch-remote-tmpdir",
        type=str,
        default=None,
        help=(
            "gs:// path used by the Hail Batch ServiceBackend for its"
            " scratch. Defaults to the standard batch tmp bucket."
        ),
    )
    fanout_group.add_argument(
        "--methods-branch",
        type=str,
        default="main",
        help=(
            "Branch or commit of gnomad_methods to pull at runtime." " Default 'main'."
        ),
    )
    fanout_group.add_argument(
        "--batch-dry-run",
        help="Validate the fan-out DAG without actually running jobs.",
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

    main(args)
