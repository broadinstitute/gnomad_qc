"""Script to compute coverage, allele number, and quality histograms on all gnomAD v5 genomes (AoU v8 + updated gnomAD v4)."""

import argparse
import logging
from typing import List, Optional

import hail as hl
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
from gnomad_qc.v4.resources.annotations import get_downsampling as get_v4_downsampling
from gnomad_qc.v4.resources.meta import meta as v4_meta

# TODO: Switch from v4>v5 once v5 sample QC is complete
from gnomad_qc.v5.resources.annotations import (  # get_aou_downsampling
    group_membership,
    qual_hists,
)
from gnomad_qc.v5.resources.basics import (  # get_aou_vds
    get_gnomad_v5_genomes_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.meta import project_meta  # meta
from gnomad_qc.v5.resources.release import (
    release_all_sites_an_tsv_path,
    release_coverage_path,
    release_coverage_tsv_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("aou_coverage_and_an")
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
    # TODO: Update to v5 downsamplings when that exists.
    downsamplings = DOWNSAMPLINGS["v3"]
    ht = annotate_downsamplings(ht, downsamplings, ht.gen_anc)
    return ht


def get_group_membership_ht(
    meta_ht: hl.Table, project: str, ds_ht: Optional[hl.Table] = None
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    :param meta_ht: Meta HT.
    :param project: Project name.
    :param ds_ht: Optional downsampling HT. Only used for AoU.
    :return: Group membership HT.
    """
    if project == "aou":
        ht = generate_freq_group_membership_array(
            meta_ht,
            build_freq_stratification_list(
                sex_expr=meta_ht.sex_karyotype,
                gen_anc_expr=meta_ht.gen_anc,
                downsampling_expr=ds_ht[meta_ht.key].downsampling,
            ),
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
        # Filter to v4 release samples and drop consent drop samples.
        # NOTE: Not using v5 project meta here because this part will be run in
        # Dataproc.
        ht = meta_ht.filter(
            meta_ht.release
            & (
                hl.is_missing(meta_ht.project_meta.research_project_key)
                | (
                    (meta_ht.project_meta.research_project_key != "RP-1061")
                    & (meta_ht.project_meta.research_project_key != "RP-1411")
                )
            )
        )
        ht = generate_freq_group_membership_array(
            ht,
            build_freq_stratification_list(
                sex_expr=ht.sex_imputation.sex_karyotype,
                gen_anc_expr=ht.population_inference.pop,
            ),
        )

    return ht


def compute_all_release_stats_per_ref_site(
    vds: hl.vds.VariantDataset,
    ref_ht: hl.Table,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    interval_ht: Optional[hl.Table] = None,
    group_membership_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Compute coverage, allele number, and quality histograms per reference site.

    .. note::

        Running this function prior to calculating frequencies removes the need for an additional
        densify for frequency calculations.

    :param vds: Input VDS.
    :param ref_ht: Reference HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X depth.
    :param interval_ht: Interval HT.
    :param group_membership_ht: Group membership HT.
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
        "qual_hists": (lambda t: [t.GQ, t.DP, t.adj], _get_hists),
        "coverage_stats": get_coverage_agg_func(dp_field="DP", max_cov_bin=max_cov_bin),
    }

    logger.info(
        "Computing coverage, allele number, and qual hists per reference site..."
    )
    # Below we use just the raw group for qual hist computations because qual hists
    # has its own built-in adj filtering when adj is passed as an argument and will
    # produce both adj and raw histograms.

    vmt = vds.variant_data
    vmt = vmt.annotate_cols(sex_karyotype=vmt.meta.sex_karyotype)
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    # TODO: Add adj annotation and sex ploidy adjustment here for AoU.
    # This is to avoid the issue we saw in v4 where multiallelic variants have
    # different adj ANs because adj and sex ploidy adjustments were not included in AN calculation.
    # For more details, see:
    # https://github.com/broadinstitute/gnomad_qc/blob/23b962a08a054df48aa4a29fc8bdca298cb26b7a/gnomad_qc/v4/annotations/fix_freq_an.py#L583
    # We should calculate AN correctly for AoU, but we should be consistent with v4
    # for gnomAD samples that we need to remove from v5 (i.e., samples to drop
    # for consent reasons).
    ht = compute_stats_per_ref_site(
        vds,
        ref_ht,
        entry_agg_funcs,
        interval_ht=interval_ht,
        group_membership_ht=group_membership_ht,
        entry_keep_fields=["GQ", "DP"],
        entry_agg_group_membership={"qual_hists": [{"group": "raw"}]},
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
    # `qual_hists` as returned by `compute_stats_per_ref_site` is an array of length 1 so we drop the array here.
    return ht.annotate(qual_hists=ht.qual_hists[0])


def join_aou_and_gnomad_coverage_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
    v4_count: int = 76215,
) -> hl.Table:
    """
    Join AoU and gnomAD coverage HTs for release.

    :param aou_ht: AoU coverage HT.
    :param gnomad_ht: gnomAD coverage HT.
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
        Default is [1, 5, 10, 15, 20, 25, 30, 50, 100].
    :param v4_count: Number of release gnomAD v4 genome samples. Default is 76215.
    :return: Joined HT.
    """
    # TODO: Add support for subtracting gnomAD v4 samples.
    # Get total number of AoU v8 release samples.
    aou_count = hl.eval(aou_ht.coverage_stats_meta_sample_count)

    def _rename_cov_annotations(ht: hl.Table, project: str) -> hl.Table:
        """
        Rename coverage annotations prior to merging Tables.

        :param ht: Input HT.
        :param project: Project name.
        :return: Renamed HT.
        """
        # Transform mean back into sum.
        sample_count = v4_count if project == "gnomad" else aou_count
        ht = ht.transmute(sum=ht.mean * sample_count)

        # Rename annotations to include project.
        row_fields = list(ht.row_value)
        rename_dict = {f: f"{f}_{project}" for f in row_fields}
        ht = ht.rename(rename_dict)
        if project == "gnomad":
            # Revert v4 genomes fraction over X bins to sample count over X bins.
            ht = ht.transmute(
                **{
                    f"over_{x}_{project}": ht[f"over_{x}_{project}"] * v4_count
                    for x in coverage_over_x_bins
                }
            )
        return ht

    aou_ht = _rename_cov_annotations(aou_ht, "aou")
    gnomad_ht = _rename_cov_annotations(gnomad_ht, "gnomad")

    total_count = aou_count + v4_count
    ht = aou_ht.join(gnomad_ht, "left")
    merged_fields = {
        "mean": (ht.sum_aou + ht.sum_gnomad) / total_count,
    }
    merged_fields.update(
        {
            f"over_{x}": (ht[f"over_{x}_aou"] + ht[f"over_{x}_gnomad"]) / total_count
            for x in coverage_over_x_bins
        }
    )
    ht = ht.transmute(**merged_fields)
    return ht.select_globals()


def _rename_fields(ht: hl.Table, field_name: str, project: str) -> hl.Table:
    """
    Rename fields by adding project name prior to merging Tables.

    Used for AN and qual hists merging but not coverage because
    coverage does not have globals and also requires extra transformations
    to transform v4 mean back into sum.

    :param ht: Input HT.
    :param project: Project name.
    :return: Renamed HT.
    """
    rename_globals = {
        f"strata_meta_{project}": ht.strata_meta,
        f"strata_sample_count_{project}": ht.strata_sample_count,
    }
    ht = ht.rename_globals(rename_globals)
    rename_dict = {f"{field_name}": f"{field_name}_{project}"}
    return ht.rename(rename_dict)


def join_aou_and_gnomad_an_ht(
    aou_ht: hl.Table,
    gnomad_ht: hl.Table,
) -> hl.Table:
    """
    Join AoU and gnomAD AN HTs for release.

    :param aou_ht: AoU AN HT.
    :param gnomad_ht: gnomAD AN HT.
    :return: Joined HT.
    """
    # TODO: Add support for subtracting gnomAD v4 samples.
    aou_ht = _rename_fields(aou_ht, "AN", "aou")
    gnomad_ht = _rename_fields(gnomad_ht, "AN", "gnomad")
    ht = aou_ht.join(gnomad_ht, "left")
    ht = ht.checkpoint(new_temp_file("aou_and_gnomad_join", "ht"))

    joint_an, joint_strata_meta, count_arrays_dict = merge_array_expressions(
        arrays=[ht.AN_aou, ht.AN_gnomad],
        meta=[
            ht.index_globals().strata_meta_aou,
            ht.index_globals().strata_meta_gnomad,
        ],
        count_arrays={
            "counts": [
                ht.index_globals().strata_sample_count_aou,
                ht.index_globals().strata_sample_count_gnomad,
            ],
        },
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
    aou_ht = _rename_fields(aou_ht, "qual_hists", "aou")
    gnomad_ht = _rename_fields(gnomad_ht, "qual_hists", "gnomad")
    ht = aou_ht.join(gnomad_ht, "left")
    ht = ht.annotate(
        qual_hists=merge_histograms(
            [ht.qual_hists_aou, ht.qual_hists_gnomad],
            operation="sum",
        ),
    )
    return ht


def main(args):
    """Compute all sites coverage, allele number, and quality histograms for v5 genomes (AoU v8 + gnomAD v4)."""
    if args.environment == "rwb":
        hl.init(
            log="/home/jupyter/workspaces/gnomadproduction/compute_coverage.log",
            tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
        )
    else:
        hl.init(
            log="compute_coverage.log",
            tmp_dir="gs://gnomad-tmp-4day",
        )
    hl.default_reference("GRCh38")

    test_2_partitions = args.test_2_partitions
    test_chr22_chrx_chry = args.test_chr22_chrx_chry
    test = test_2_partitions or test_chr22_chrx_chry
    overwrite = args.overwrite
    n_partitions = args.n_partitions

    try:
        # Retrieve raw coverage table path here because it is used in most of the
        # script's run options.
        cov_and_an_ht_path = release_coverage_path(
            public=False,
            test=test,
            raw=True,
            coverage_type="coverage",
            data_set="aou",
            environment="rwb",
        )
        # TODO: Update this to use get_aou_downsampling once sample QC is complete.
        downsampling_ht_path = get_v4_downsampling(test=test).path
        meta_ht_path = project_meta.path

        if args.write_gnomad_group_membership_ht:
            logger.info("Writing group membership HT...")
            v4_meta_ht_path = v4_meta(data_type="genomes").path
            group_membership_ht_path = group_membership(
                test=args.test, data_set="gnomad"
            ).path
            check_resource_existence(
                input_step_resources={"v4_meta_ht": [v4_meta_ht_path]},
                output_step_resources={
                    "gnomad_group_membership_ht": [group_membership_ht_path]
                },
                overwrite=overwrite,
            )

            ht = get_group_membership_ht(
                hl.read_table(v4_meta_ht_path), project="gnomad"
            )
            ht.write(group_membership_ht_path, overwrite=overwrite)

        if args.write_aou_downsampling_ht:
            check_resource_existence(
                output_step_resources={"downsampling_ht": downsampling_ht_path},
                overwrite=overwrite,
            )
            # TODO: Update this for v5 once sample QC is complete.
            ht = hl.read_table(meta_ht_path)
            ht = ht.filter(
                (ht.project == "gnomad") & (ht.data_type == "genomes") & (ht.release)
            )
            # ds_ht = get_downsampling_ht(ht)
            ds_ht = get_downsampling_ht(ht)
            ds_ht.write(downsampling_ht_path, overwrite=overwrite)

        if args.compute_all_cov_release_stats_ht:
            # Context Table is used because it contains every locus in the GRCh38
            # reference as opposed to a ref-blocked VDS reference dataset.
            ref_ht = vep_context.versions["105"].ht()
            if test_chr22_chrx_chry:
                chrom = ["chr22", "chrX", "chrY"]
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

            vds = get_gnomad_v5_genomes_vds(
                release_only=True,
                remove_hard_filtered_samples=True,
                test=test,
                filter_partitions=range(2) if test_2_partitions else None,
                annotate_meta=True,
                chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
            )
            # vds = get_aou_vds(
            #    release_only=True,
            #    test=test,
            #    filter_partitions=range(2) if test_2_partitions else None,
            #    annotate_meta=True,
            #    chrom=["chr22", "chrX", "chrY"] if test_chr22_chrx_chry else None,
            # )

            logger.info(
                "Computing coverage, all sites allele number, and quality histograms HT..."
            )
            aou_group_membership_ht_path = group_membership(test=test).path
            check_resource_existence(
                output_step_resources={
                    "aou_group_membership_ht": aou_group_membership_ht_path,
                    "coverage_and_an_ht": cov_and_an_ht_path,
                },
                overwrite=overwrite,
            )
            aou_group_membership_ht = get_group_membership_ht(
                meta_ht=hl.read_table(meta_ht_path),
                project="aou",
                ds_ht=hl.read_table(downsampling_ht_path),
            )
            aou_group_membership_ht = aou_group_membership_ht.checkpoint(
                aou_group_membership_ht_path, overwrite=overwrite
            )

            cov_and_an_ht = compute_all_release_stats_per_ref_site(
                vds,
                ref_ht,
                group_membership_ht=aou_group_membership_ht,
            )
            cov_and_an_ht = cov_and_an_ht.checkpoint(
                new_temp_file("aou_cov_and_an", "ht")
            )
            # Write out the intermediate HT.
            cov_and_an_ht = cov_and_an_ht.naive_coalesce(n_partitions)
            cov_and_an_ht.write(cov_and_an_ht_path, overwrite=overwrite)

        if args.export_coverage_release_files:
            cov_ht_path = release_coverage_path(
                public=False,
                test=test,
                raw=False,
                coverage_type="coverage",
                data_set="joint",
            )
            cov_tsv_path = release_coverage_tsv_path(test=test)
            check_resource_existence(
                output_step_resources={
                    "cov_release_ht": cov_ht_path,
                    "cov_tsv": cov_tsv_path,
                },
                overwrite=overwrite,
            )

            logger.info("Exporting coverage HT and TSV...")
            aou_ht = hl.read_table(cov_and_an_ht_path)
            aou_ht = aou_ht.drop("AN", "qual_hists")
            gnomad_ht = hl.read_table(
                release_coverage_path(
                    data_type="genomes",
                    release_version=v4_COVERAGE_RELEASE,
                    public=True,
                    raw=False,
                )
            )
            if test_chr22_chrx_chry:
                chrom = ["chr22", "chrX", "chrY"]
                gnomad_ht = hl.filter_intervals(
                    gnomad_ht, [hl.parse_locus_interval(c) for c in chrom]
                )
            elif test_2_partitions:
                gnomad_ht = gnomad_ht._filter_partitions(range(2))

            ht = join_aou_and_gnomad_coverage_ht(aou_ht, gnomad_ht)
            ht = ht.checkpoint(new_temp_file("aou_and_gnomad_cov_join", "ht"))
            ht = ht.naive_coalesce(n_partitions)
            ht = ht.checkpoint(cov_ht_path, overwrite=overwrite)
            ht.export(cov_tsv_path)

        if args.export_an_release_files:
            an_raw_ht_path = release_coverage_path(
                public=False,
                test=test,
                raw=True,
                coverage_type="allele_number",
                data_set="joint",
            )
            an_ht_path = release_coverage_path(
                public=False,
                test=test,
                raw=False,
                coverage_type="allele_number",
                data_set="joint",
            )
            an_tsv_path = release_all_sites_an_tsv_path(test=test)
            check_resource_existence(
                output_step_resources={
                    "an_raw_ht": an_raw_ht_path,
                    "an_release_ht": an_ht_path,
                    "an_release_tsv": an_tsv_path,
                },
                overwrite=overwrite,
            )

            logger.info("Exporting AN HT and TSV...")
            aou_ht = hl.read_table(cov_and_an_ht_path)
            aou_ht = aou_ht.select("AN")
            gnomad_ht = hl.read_table(
                release_coverage_path(
                    data_type="genomes",
                    release_version=v4_AN_RELEASE,
                    public=True,
                    raw=False,
                    coverage_type="allele_number",
                )
            )
            if test_chr22_chrx_chry:
                chrom = ["chr22", "chrX", "chrY"]
                gnomad_ht = hl.filter_intervals(
                    gnomad_ht, [hl.parse_locus_interval(c) for c in chrom]
                )
            elif test_2_partitions:
                gnomad_ht = gnomad_ht._filter_partitions(range(2))

            ht = join_aou_and_gnomad_an_ht(aou_ht, gnomad_ht)
            ht = ht.checkpoint(an_raw_ht_path, overwrite=overwrite)
            ht = ht.select("AN")
            ht = ht.select_globals(
                strata_meta=ht.strata_meta,
                strata_sample_count=ht.strata_sample_count,
            )
            ht = ht.naive_coalesce(n_partitions)
            ht = ht.checkpoint(an_ht_path, overwrite=overwrite)

            # Only export the adj AN for all release samples.
            ht = ht.transmute(AN=ht.AN[0])
            ht.export(an_tsv_path)

        if args.merge_qual_hists:
            qual_hists_path = qual_hists(test=test).path
            check_resource_existence(
                output_step_resources={"qual_hists_ht": qual_hists_path},
                overwrite=overwrite,
            )

            logger.info("Merging qual hists HTs...")
            aou_ht = hl.read_table(cov_and_an_ht_path).select("qual_hists")
            gnomad_ht = public_release(data_type="genomes").select("histograms")

            # Drop age hists because they are handled in the frequency script.
            gnomad_ht = gnomad_ht.annotate(
                qual_hists=gnomad_ht.histograms.drop("age_hists")
            )

            ht = join_aou_and_gnomad_qual_hists_ht(aou_ht, gnomad_ht)
            ht.write(qual_hists_path, overwrite=overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compute_coverage", environment=args.environment))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--environment",
        help="Environment where script will run.",
        default="rwb",
        type=str,
        choices=["rwb", "dataproc"],
    )
    parser.add_argument(
        "--project",
        help="Project name.",
        default="aou",
        type=str,
        choices=["aou", "gnomad"],
    )
    parser.add_argument(
        "--overwrite", help="Overwrite existing hail Tables.", action="store_true"
    )
    parser.add_argument(
        "--test-2-partitions",
        help=(
            "Whether to run a test using only the first 2 partitions of the VDS test"
            " dataset."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--test-chr22-chrx-chry",
        help=(
            "Whether to run a test using only the chr22, chrX, and chrY chromosomes of"
            " the VDS test dataset."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of partitions to use for the output Table.",
        type=int,
        default=5000,
    )

    gnomad_group_membership_args = parser.add_argument_group(
        "Get gnomAD genomes group membership HT.",
    )
    gnomad_group_membership_args.add_argument(
        "--write-gnomad-group-membership-ht",
        help="Write gnomAD genomes group membership HT.",
        action="store_true",
    )
    gnomad_group_membership_args.add_argument(
        "--test",
        help="Write test gnomAD genomes group membership HT to test path.",
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
    parser.add_argument(
        "--export-coverage-release-files",
        help="Exports joint AoU + gnomAD v4 coverage release HT and TSV file.",
        action="store_true",
    )
    parser.add_argument(
        "--export-an-release-files",
        help="Exports joint AoU + gnomAD v4 AN release HT and TSV file.",
        action="store_true",
    )
    parser.add_argument(
        "--merge-qual-hists",
        help="Merge variant quality histograms from AoU v8 and gnomAD v4 genomes.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.project == "aou" and args.environment != "rwb":
        parser.error("--project aou requires --environment rwb")
    if args.project == "gnomad" and args.environment != "dataproc":
        parser.error("--project gnomad requires --environment dataproc")

    gnomad_args = [arg for arg in vars(args) if "gnomad" in arg.lower()]
    if any(
        getattr(args, arg)
        for arg in gnomad_args
        if isinstance(getattr(args, arg), bool)
    ):
        if args.environment != "dataproc":
            parser.error(
                "Steps run on only gnomAD must be run in Dataproc (`--environment dataproc`)!"
            )

    main(args)
