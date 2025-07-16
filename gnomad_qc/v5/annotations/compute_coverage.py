"""Script to compute coverage statistics on gnomAD v5 genomes."""

import argparse
import logging
from typing import List, Optional

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.annotations import (
    build_freq_stratification_list,
    generate_freq_group_membership_array,
    qual_hist_expr,
)
from gnomad.utils.sparse_mt import (
    compute_stats_per_ref_site,
    get_allele_number_agg_func,
    get_coverage_agg_func,
)

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds
from gnomad_qc.v4.resources.meta import meta

# TODO: Switch from v4>v5 once v5 sample QC is complete
from gnomad_qc.v5.resources.basics import get_logging_path  # get_aou_vds
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET

# from gnomad_qc.resources.meta import meta
from gnomad_qc.v5.resources.release import (
    release_coverage_path,
    release_coverage_tsv_path,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("coverage")
logger.setLevel(logging.INFO)


def get_genomes_group_membership_ht(
    meta_ht: hl.Table,
) -> hl.Table:
    """
    Get genomes group membership HT for all sites allele number stratification.

    :param meta_ht: Metadata HT.
    :return: Group membership HT.
    """
    # Filter to release samples.
    meta_ht = meta_ht.filter(meta_ht.release)

    # Create group membership HT.
    ht = generate_freq_group_membership_array(
        meta_ht,
        build_freq_stratification_list(
            sex_expr=meta_ht.sex_imputation.sex_karyotype,
            pop_expr=meta_ht.gen_anc_inference.gen_anc,
        ),
    )

    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: hl.dict(
                d.items().map(lambda x: hl.if_else(x[0] == "pop", ("gen_anc", x[1]), x))
            )
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
    :param coverage_over_x_bins: List of boundaries for computing samples over X.
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
    vmt = vmt.annotate_cols(sex_karyotype=vmt.meta.sex_imputation.sex_karyotype)
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

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
        cov_stat: hl.expr.StructExpression, n: hl.expr.Int32Expression
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
        bin_expr = {f"over_{x}": bin_expr[i] / n for i, x in enumerate(rev_cov_bins)}
        return cov_stat.annotate(**bin_expr).drop("coverage_counter")

    # Keep coverage stats from global adj grouping (index 0) only.
    ht = ht.annotate_globals(
        coverage_stats_meta_sample_count=ht.strata_sample_count[0],
    )
    cov_stats_expr = _cov_stats(ht.coverage_stats[0], ht.sample_count)

    ht = ht.transmute(**cov_stats_expr)
    return ht.annotate(qual_hists=ht.qual_hists[0])


def main(args):
    """Compute coverage statistics, including mean, median_approx, and coverage over certain DPs."""
    hl.init(
        log="/home/jupyter/workspaces/gnomadproduction/compute_coverage.log",
        tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
    )
    hl.default_reference("GRCh38")

    # TODO: Remove this?
    # SSA Logs are easier to troubleshoot with.
    hl._set_flags(use_ssa_logs="1")

    test_2_partitions = args.test_2_partitions
    test_chr22_chrx_chry = args.test_chr22_chrx_chry
    test = test_2_partitions or test_chr22_chrx_chry
    overwrite = args.overwrite
    n_partitions = args.n_partitions

    try:
        # NOTE: v5 genomes coverage returns Table with coverage, all sites AN, and
        # qual hists.
        cov_and_an_ht_path = release_coverage_path(
            public=False,
            test=test,
            include_meta=True,
            coverage_type="coverage",
        )

        if args.compute_all_cov_release_stats_ht:
            # Read in context Table.
            ref_ht = vep_context.versions["105"].ht()
            if test_chr22_chrx_chry:
                chrom = ["chr22", "chrX", "chrY"]
                ref_ht = hl.filter_intervals(
                    ref_ht, [hl.parse_locus_interval(c) for c in chrom]
                )
            elif test_2_partitions:
                ref_ht = ref_ht._filter_partitions(range(5))

            # Retain only 'locus' annotation from context Table.
            ref_ht = ref_ht.key_by("locus").select().distinct()

            # Read in VDS.
            vds = get_gnomad_v4_genomes_vds(
                release_only=True,
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
                "Running compute coverage, all sites allele number, and quality histograms HT..."
            )
            check_resource_existence(
                output_step_resources={"coverage_and_an_ht": cov_and_an_ht_path}
            )
            group_membership_ht = get_genomes_group_membership_ht(
                meta().ht(),
            )
            group_membership_ht = group_membership_ht.checkpoint(
                hl.utils.new_temp_file("group_membership", "ht")
            )

            cov_and_an_ht = compute_all_release_stats_per_ref_site(
                vds,
                ref_ht,
                group_membership_ht=group_membership_ht,
            )
            cov_and_an_ht = cov_and_an_ht.checkpoint(
                hl.utils.new_temp_file("cov_and_an", "ht")
            )
            # Naive coalesce and write out the intermediate HT.
            cov_and_an_ht = cov_and_an_ht.naive_coalesce(n_partitions)
            cov_and_an_ht.write(cov_and_an_ht_path, overwrite=overwrite)

        if args.export_coverage_release_files:
            release_ht_path = release_coverage_path(
                public=False, test=test, include_meta=False, coverage_type="coverage"
            )
            release_tsv_path = release_coverage_tsv_path(test=test)
            check_resource_existence(
                input_step_resources={"cov_and_an_ht": cov_and_an_ht_path},
                output_step_resources={
                    "cov_release_ht": release_ht_path,
                    "cov_tsv": release_tsv_path,
                },
            )

            logger.info("Exporting coverage and AN HT...")
            # Select coverage and AN fields for release.
            ht = hl.read_table(cov_and_an_ht_path)
            ht = ht.select_globals(
                an_strata_meta=ht.strata_meta,
                an_strata_sample_count=ht.strata_sample_count,
            )
            ht = ht.select(
                "AN",
                **{k: ht.coverage_stats[0][k] for k in ht.coverage_stats[0]},
            )
            ht = ht.checkpoint(release_ht_path, overwrite=overwrite)

            logger.info("Exporting coverage and AN tsv...")
            # Only export the adj AN for all release samples.
            ht = ht.transmute(AN=ht.AN[0])
            ht.export(release_tsv_path)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("compute_coverage", environment="rwb"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
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
    parser.add_argument(
        "--compute-coverage-ht", help="Compute the coverage HT.", action="store_true"
    )
    parser.add_argument(
        "--compute-all-cov-release-stats-ht",
        help="Compute the all sites coverage, allele number, and quality histogram HT.",
        action="store_true",
    )
    parser.add_argument(
        "--export-coverage-release-files",
        help="Exports coverage release HT and TSV file.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
