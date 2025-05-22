"""Script to determine samples that fail hard filtering thresholds."""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.utils.annotations import bi_allelic_expr

from gnomad_qc.v5.resources.basics import get_aou_vds, get_checkpoint_path
from gnomad_qc.v5.resources.constants import WORKSPACE_BUCKET
from gnomad_qc.v5.resources.sample_qc import get_sample_qc, hard_filtered_samples

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hard_filters")
logger.setLevel(logging.INFO)


def compute_aou_sample_qc(
    n_partitions: int = 500,
    test: bool = False,
) -> hl.Table:
    """
    Perform sample QC on AoU VDS.

    .. note::

        We are not including the `n_alt_alleles_strata` parameter in this function—as we did for v4 exomes—
        because the distribution of alternate alleles in whole genome sequencing data is not as skewed as in exomes.
        For example, in AoU v8 genomes, 77.06% of variants are bi-allelic, compared to 76.65% in v4 genomes and
        only 35.77% in v4 exomes.

    :param n_partitions: Number of partitions to use when writing the sample QC table.
    :param test: If true, test the function on a smaller subset of the data.
    :return: Table containing sample QC metrics
    """
    logger.info("Computing sample QC")
    logger.info("Loading test VDS..." if test else "Loading VDS...")

    if test:
        n_partitions = n_partitions // 100

    vds = get_aou_vds(
        test=test,
        autosomes_only=True,
        split=False,
    )

    logger.info(
        "Excluding telomeres and centromeres from VDS (redundant but acts as a safety check)..."
    )
    # AoU pipeline has already filtered out the centromeres and telomeres, but
    # this serves as an additional safeguard.
    vds = hl.vds.filter_intervals(
        vds, intervals=telomeres_and_centromeres.ht(), keep=False
    )

    logger.info("Excluding loci with more than 100 alternative alleles...")
    # NOTE: Filtering loci with >100 alleles has no effect on the test dataset,
    # since it contains no such loci. However, this filtering is still necessary
    # on the full dataset because the number of alleles can change after removing
    # samples from the raw VDS. As a result, it may not always match the ``EXCESS_ALLELES`` flag
    # in the ``filters`` field of the raw VDS. (AoU v8 defined this as a site with >100 alternate alleles
    # in their quality report: `All of Us Genomic Quality Report
    # <https://support.researchallofus.org/hc/en-us/articles/29390274413716-All-of-Us-Genomic-Quality-Report>`_)
    vmt = vds.variant_data
    vmt = vmt.annotate_rows(n_unsplit_alleles=hl.len(vmt.alleles))
    vmt = vmt.filter_rows(vmt.n_unsplit_alleles < 101)

    logger.info("Splitting multi-allelic variants...")
    vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    logger.info("Computing sample QC metrics...")
    sample_qc_ht = compute_stratified_sample_qc(
        vds,
        strata={
            "bi_allelic": bi_allelic_expr(vds.variant_data),
            "multi_allelic": ~bi_allelic_expr(vds.variant_data),
        },
        tmp_ht_prefix=get_sample_qc(test=test).path[:-3],
        gt_col="GT",
    )

    return sample_qc_ht.naive_coalesce(n_partitions)


def compute_hard_filters(
    max_n_singleton: float = 100000,
    max_r_het_hom_var: float = 10,
    max_r_insertion_deletion: float = 0.42,
    min_r_ti_tv: float = 2.9,
    max_r_ti_tv_singleton: float = 5.2,
    test: bool = False,
) -> hl.Table:
    """
    Apply hard filters to samples and return a Table with the filtered samples and the reason for filtering.

    Function only filters outliers based on sample QC metrics; samples were filtered upstream based on
    contamination, sex karyotype, and coverage.

    .. warning::
        The defaults used in this function are callset specific; these hardfilter
        cutoffs will need to be re-examined for each callset.

    :param max_n_singleton: Filtering threshold to use for the maximum number of
        singletons.
    :param max_r_het_hom_var: Filtering threshold to use for the maximum ratio of
        heterozygotes to alternate homozygotes.
    :param max_r_insertion_deletion: Filtering threshold to use for the maximum ratio of
        insertions to deletions.
    :param min_r_ti_tv: Filtering threshold to use for the minimum ratio of
        transitions to transverions.
    :param max_r_ti_tv_singleton: Filtering threshold to use for maximum ratio of
        transitions to tranversions in singletons.
    :param test: Whether to use the gnomAD v4 test dataset. Default is to use the full
        dataset.
    :return: Table of hard filtered samples.
    """
    ht = get_aou_vds(
        test=test,
    ).variant_data.cols()
    ht = ht.annotate_globals(
        hard_filter_cutoffs=hl.struct(
            max_n_singleton=max_n_singleton,
            max_r_het_hom_var=max_r_het_hom_var,
            max_r_insertion_deletion=max_r_insertion_deletion,
            min_r_ti_tv=min_r_ti_tv,
            max_r_ti_tv_singleton=max_r_ti_tv_singleton,
        ),
    )

    # Flag extreme raw bi-allelic sample QC outliers.
    sample_qc_metric_hard_filters = dict()
    bi_allelic_qc_ht = get_sample_qc("bi_allelic", test=test).ht()

    bi_allelic_qc_struct = bi_allelic_qc_ht[ht.key]
    sample_qc_metric_hard_filters["high_n_singleton"] = (
        bi_allelic_qc_struct.n_singleton > max_n_singleton
    )
    sample_qc_metric_hard_filters["high_r_het_hom_var"] = (
        bi_allelic_qc_struct.r_het_hom_var > max_r_het_hom_var
    )
    sample_qc_metric_hard_filters["high_r_ins_del"] = (
        bi_allelic_qc_struct.r_insertion_deletion < max_r_insertion_deletion
    )
    sample_qc_metric_hard_filters["low_r_ti_tv"] = (
        bi_allelic_qc_struct.r_ti_tv < min_r_ti_tv
    )
    sample_qc_metric_hard_filters["high_r_ti_tv_singleton"] = (
        bi_allelic_qc_struct.r_ti_tv_singleton < max_r_ti_tv_singleton
    )
    ht = ht.annotate(
        sample_qc_metric_hard_filters=add_filters_expr(
            filters=sample_qc_metric_hard_filters
        ),
    )

    # Keep samples failing hard filters.
    ht = ht.filter(hl.len(ht.sample_qc_metric_hard_filters) > 0)
    return ht


def main(args):
    """Determine samples that fail hard filtering thresholds."""
    hl.init(tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day")
    hl.default_reference("GRCh38")

    test = args.test
    overwrite = args.overwrite

    if args.aou_sample_qc:
        ht = compute_aou_sample_qc(
            n_partitions=args.sample_qc_n_partitions,
            test=test,
        )
        ht.write(get_sample_qc(test=test).path, overwrite=overwrite)

    if args.compute_hard_filters:
        hard_filter_path = (
            get_checkpoint_path("test_aou_hard_filters")
            if test
            else hard_filtered_samples.path
        )
        ht = compute_hard_filters(
            args.max_n_singleton,
            args.max_r_het_hom_var,
            args.max_r_insertion_deletion,
            args.min_r_ti_tv,
            args.max_r_ti_tv_singleton,
            test,
        )
        ht = ht.checkpoint(hard_filter_path, overwrite=overwrite)
        ht.group_by("sample_qc_metric_hard_filters").aggregate(n=hl.agg.count()).show(
            20
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all MatrixTables/Tables. (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Use the AoU test dataset instead of the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--aou-sample-qc",
        help="Compute Hail's VDS sample QC metrics on AoU.",
        action="store_true",
    )
    parser.add_argument(
        "--sample-qc-n-partitions",
        help="Number of partitions to use when writing the sample QC table.",
        type=int,
        default=500,
    )
    parser.add_argument(
        "--compute-hard-filters",
        help=(
            "Computes samples to be hard-filtered. NOTE: Cutoffs should be determined"
            " by visual inspection of the metrics."
        ),
        action="store_true",
    )
    parser.add_argument_group()
    hard_filter_args = parser.add_argument_group(
        "Hard-filter cutoffs", "Arguments used for hard-filter cutoffs."
    )
    hard_filter_args.add_argument(
        "--max-n-singleton",
        type=float,
        default=100000,
        help=(
            "Filtering threshold to use for the maximum number of singletons. Default"
            " is 100000."
        ),
    )
    hard_filter_args.add_argument(
        "--max-r-het-hom-var",
        type=float,
        default=2.5,
        help=(
            "Filtering threshold to use for the maximum ratio of heterozygotes to"
            " alternate homozygotes. Default is 2.5."
        ),
    )
    hard_filter_args.add_argument(
        "--max-r-insertion-deletion",
        type=float,
        default=0.42,
        help=(
            "Filtering threshold to use for the maximum ratio of insertions to"
            " deletions. Default is 0.42."
        ),
    )
    hard_filter_args.add_argument(
        "--min-r-ti-tv",
        type=float,
        default=2.9,
        help=(
            "Filtering threshold to use for the minimum ratio of transitions to"
            " transverions. Default is 2.9."
        ),
    )
    hard_filter_args.add_argument(
        "--max-r-ti-tv-singleton",
        type=float,
        default=5.2,
        help=(
            "Filtering threshold to use for maximum ratio of transitions to"
            " tranversions in singletons. Default is 5.2."
        ),
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
