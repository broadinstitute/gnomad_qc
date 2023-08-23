# noqa: D100

import argparse
import logging
from typing import List, Union

import hail as hl
from gnomad.resources.grch38.gnomad import (
    CURRENT_GENOME_COVERAGE_RELEASE,
    coverage,
    coverage_tsv_path,
)
from gnomad.resources.grch38.reference_data import (
    telomeres_and_centromeres,
    vep_context,
)
from gnomad.utils.sparse_mt import compute_coverage_stats

from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("coverage")
logger.setLevel(logging.INFO)


def compute_coverage_stats(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    reference_ht: hl.Table,
    coverage_over_x_bins: List[int] = [1, 5, 10, 15, 20, 25, 30, 50, 100],
) -> hl.Table:
    """
    Compute coverage statistics for every base of the `reference_ht` provided.

    The following coverage stats are calculated:
        - mean
        - median
        - total DP
        - fraction of samples with coverage above X, for each x in `coverage_over_x_bins`

    The `reference_ht` is a table that contains row for each locus coverage should be computed on.
    It needs to be keyed with the same keys as `mt`, typically either `locus` or `locus, alleles`.
    The `reference_ht` can e.g. be created using `get_reference_ht`

    :param mtds: Input sparse MT
    :param reference_ht: Input reference HT
    :param coverage_over_x_bins: List of boundaries for computing samples over X
    :return: Table with per-base coverage stats
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)
    if is_vds:
        n_samples = mtds.variant_data.count_cols()
    else:
        n_samples = mtds.count_cols()

    print(f"Computing coverage stats on {n_samples} samples.")

    # Create an outer join with the reference Table
    def join_with_ref(mt: hl.MatrixTable) -> hl.MatrixTable:
        keep_entries = ["DP"]
        if "END" in mt.entry:
            keep_entries.append("END")
        mt.select_entries(*keep_entries).select_cols().select_rows()
        col_key_fields = list(mt.col_key)
        t = mt._localize_entries("__entries", "__cols")
        t = t.join(reference_ht.key_by(*mt.row_key).select(_in_ref=True), how="outer")
        t = t.annotate(
            __entries=hl.or_else(
                t.__entries,
                hl.range(n_samples).map(
                    lambda x: hl.missing(t.__entries.dtype.element_type)
                ),
            )
        )
        return t._unlocalize_entries("__entries", "__cols", col_key_fields)

    if is_vds:
        rmt = mtds.reference_data
        vmt = mtds.variant_data
        mtds = hl.vds.VariantDataset(join_with_ref(rmt), join_with_ref(vmt))
        mt = hl.vds.to_dense_mt(mtds)
    else:
        mtds = join_with_ref(mtds)
        # Densify
        mt = hl.experimental.densify(mtds)

    # Filter rows where the reference is missing
    mt = mt.filter_rows(mt._in_ref)

    # Unfilter entries so that entries with no ref block overlap aren't null
    mt = mt.unfilter_entries()

    # Compute coverage stats
    coverage_over_x_bins = sorted(coverage_over_x_bins)
    max_coverage_bin = coverage_over_x_bins[-1]
    hl_coverage_over_x_bins = hl.array(coverage_over_x_bins)

    # This expression creates a counter DP -> number of samples for DP between
    # 0 and max_coverage_bin
    coverage_counter_expr = hl.agg.counter(
        hl.min(max_coverage_bin, hl.or_else(mt.DP, 0))
    )

    # This expression aggregates the DP counter in reverse order of the coverage_over_x_bins
    # and computes the cumulative sum over them.
    # It needs to be in reverse order because we want the sum over samples
    # covered by > X.
    count_array_expr = hl.cumulative_sum(
        hl.array(
            [
                hl.int32(coverage_counter_expr.get(max_coverage_bin, 0))
            ]  # The coverage was already floored to the max_coverage_bin, so no more aggregation is needed for the max bin
            # For each of the other bins, coverage needs to be summed between the
            # boundaries
        ).extend(
            hl.range(hl.len(hl_coverage_over_x_bins) - 1, 0, step=-1).map(
                lambda i: hl.sum(
                    hl.range(
                        hl_coverage_over_x_bins[i - 1], hl_coverage_over_x_bins[i]
                    ).map(lambda j: hl.int32(coverage_counter_expr.get(j, 0)))
                )
            )
        )
    )
    mean_expr = hl.agg.mean(hl.or_else(mt.DP, 0))

    # Annotate rows now
    return mt.select_rows(
        mean=hl.if_else(hl.is_nan(mean_expr), 0, mean_expr),
        median_approx=hl.or_else(hl.agg.approx_median(hl.or_else(mt.DP, 0)), 0),
        total_DP=hl.agg.sum(mt.DP),
        **{
            f"over_{x}": count_array_expr[i] / n_samples
            for i, x in zip(
                range(
                    len(coverage_over_x_bins) - 1, -1, -1
                ),  # Reverse the bin index as count_array_expr has the reverse order
                coverage_over_x_bins,
            )
        },
    ).rows()


def main(args):  # noqa: D103
    hl.init(
        log="/coverage.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # TODO: Determine coverage version from args
    coverage_version = (
        args.coverage_version
        if args.coverage_version
        else CURRENT_GENOME_COVERAGE_RELEASE
    )
    test = args.test

    if args.compute_coverage_ht:
        # TODO: See if the context HT looks OK and use that one. Might need to remove
        #  telomeres and centromeres.
        ref_ht = vep_context.versions["105"].ht()
        if test:
            ref_ht = ref_ht._filter_partitions(range(50))

        vds = get_gnomad_v4_vds(
            release_only=True, filter_partitions=range(2) if test else None
        )

        coverage_ht = compute_coverage_stats(vds, ref_ht)

        coverage_ht = coverage_ht.checkpoint(
            "gs://gnomad-tmp/gnomad.genomes_v3.coverage.summary.ht", overwrite=True
        )
        coverage_ht.describe()
        coverage_ht.show()
        # coverage_ht = coverage_ht.naive_coalesce(5000)

        # coverage_ht.write(
        #    coverage("genomes").versions[coverage_version].path,
        #    overwrite=args.overwrite,
        # )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite", help="Overwrites existing hail tables", action="store_true"
    )
    parser.add_argument(
        "--test",
        help="",
        action="store_true",
    )
    parser.add_argument(
        "--compute_coverage_ht", help="Computes the coverage HT", action="store_true"
    )
    parser.add_argument(
        "--coverage_version",
        help=(
            "Specifies coverage version to read/write. If not set,"
            " gnomad.resources.grch38.gnomad.CURRENT_GENOME_COVERAGE_RELEASE is used."
        ),
    )

    main(parser.parse_args())
