"""Script to determine samples that fail hard filtering thresholds."""
import argparse
import logging
from typing import Optional

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.utils.annotations import annotate_adj, bi_allelic_expr
from gnomad.utils.filtering import add_filters_expr, filter_to_autosomes
from gnomad.utils.slack import slack_notifications

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    calling_intervals,
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
    gnomad_v4_testset_meta,
)
from gnomad_qc.v4.resources.meta import project_meta
from gnomad_qc.v4.resources.sample_qc import (
    contamination,
    fingerprinting_failed,
    get_predetermined_qc,
    get_sample_qc,
    hard_filtered_samples,
    interval_coverage,
    sample_chr20_mean_dp,
    sample_qc_mt_callrate,
    sex,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hard_filters")
logger.setLevel(logging.INFO)


def compute_sample_qc(
    vds: hl.vds.VariantDataset,
    n_partitions: int = 500,
    test: bool = False,
    n_alt_alleles_strata: int = 3,
    n_alt_alleles_strata_name: str = "three",
) -> hl.Table:
    """
    Perform sample QC on the raw split matrix table using `compute_stratified_sample_qc`.

    :param vds: Input VDS.
    :param n_partitions: Number of partitions to write the output sample QC HT to.
    :param test: Whether to use the gnomAD v4 test dataset. Default is to use the full
        dataset.
    :param n_alt_alleles_strata: Number of alternative alleles to stratify sample QC on.
        For example `n_alt_alleles_strata` = 2 will stratify 'bi-allelic' and
        'multi-allelic'. `n_alt_alleles_strata` = 3 will stratify under 3 alt alleles
        and 3 or more alt alleles. Default is 3.
    :param n_alt_alleles_strata_name: Name to use for number of alternative allele
        stratification, where strata names are
        'under_{n_alt_alleles_strata_name}_alt_alleles' and
        '{n_alt_alleles_strata_name}_or_more_alt_alleles'. This is not used when
        `n_alt_alleles_strata` = 2, and instead uses 'bi-allelic' and 'multi-allelic'.
        Default is 'three'.
    :return: Table containing sample QC metrics.
    """
    logger.info("Computing sample QC")
    vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)

    # Remove centromeres and telomeres in case they were included.
    vds = hl.vds.filter_intervals(
        vds, intervals=telomeres_and_centromeres.ht(), keep=False
    )
    if n_alt_alleles_strata < 2:
        raise ValueError("'n_alt_alleles_strata' must be greater than or equal to 2!")
    elif n_alt_alleles_strata == 2:
        strata = {
            "bi_allelic": bi_allelic_expr(vds.variant_data),
            "multi_allelic": ~bi_allelic_expr(vds.variant_data),
        }
    else:
        # n_unsplit_alleles includes the reference allele. So, for example, if
        # spliting at 3 alt alleles (n_alt_alleles_strata = 3) the stratification will
        # be:
        # - Under 3 alt alleles (under_three_alt_alleles), which is the same as under 4
        # total alleles (n_unsplit_alleles < 4), and the same as less than or equal to
        # 3 total alleles (n_unsplit_alleles <= 3).
        # - 3 or more alt alleles (three_or_more_alt_alleles), which is the same as 4
        # or more total alleles (n_unsplit_alleles >= 4), and the same as greater than
        # 3 total alleles (n_unsplit_alleles > 3).
        strata = {
            f"under_{n_alt_alleles_strata_name}_alt_alleles": (
                vds.variant_data.n_unsplit_alleles <= n_alt_alleles_strata
            ),
            f"{n_alt_alleles_strata_name}_or_more_alt_alleles": (
                vds.variant_data.n_unsplit_alleles > n_alt_alleles_strata
            ),
        }
    sample_qc_ht = compute_stratified_sample_qc(
        vds,
        strata=strata,
        tmp_ht_prefix=get_sample_qc(test=test).path[:-3],
        gt_col="GT",
    )

    return sample_qc_ht.repartition(n_partitions)


def compute_hard_filters(
    ht: hl.Table,
    fingerprinting_failed_ht: hl.Table,
    project_meta_ht: hl.Table,
    sample_qc_ht: hl.Table,
    contamination_ht: hl.Table,
    max_n_singleton: float = 5000,
    max_r_het_hom_var: float = 10,
    min_bases_dp_over_1: float = 5e7,
    min_bases_dp_over_20: float = 4e7,
    max_chimera: float = 0.05,
    max_contamination_estimate: float = 0.015,
    test: bool = False,
    chr20_mean_dp_ht: Optional[hl.Table] = None,
    min_cov: Optional[int] = None,
    qc_mt_callrate_ht: Optional[hl.Table] = None,
    min_qc_mt_adj_callrate: Optional[float] = None,
    sex_ht: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Apply hard filters to samples and return a Table with the filtered samples and the reason for filtering.

    If `include_sex_filter` is True, this function expects a sex inference Table
    generated by `sex_inference.py --impute-sex`.

    .. warning::
        The defaults used in this function are callset specific, these hardfilter
        cutoffs will need to be re-examined for each callset

    :param ht: Table of samples to add hard filter information to.
    :param fingerprinting_failed_ht: Table containing samples that failed fingerprinting.
    :param project_meta_ht: Table containing project metadata with 'bam_metrics'
        annotation.
    :param sample_qc_ht: Table containing sample QC metrics.
    :param contamination_ht: Table containing contamination estimates.
    :param max_n_singleton: Filtering threshold to use for the maximum number of
        singletons.
    :param max_r_het_hom_var: Filtering threshold to use for the maximum ratio of
        heterozygotes to alternate homozygotes.
    :param min_bases_dp_over_1: Filtering threshold to use for the minimum number of
        bases with a DP over one.
    :param min_bases_dp_over_20: Filtering threshold to use for the minimum number of
        bases with a DP over 20.
    :param max_chimera: Filtering threshold to use for maximum chimera (this is a
        proportion not a percent, e.g. 5% == 0.05, %5 != 5).
    :param max_contamination_estimate: Filtering threshold to use for maximum
        contamination estimate.
    :param test: Whether to use the gnomAD v4 test dataset. Default is to use the full
        dataset.
    :param chr20_mean_dp_ht: Table containing the per sample chromosome 20 mean DP.
    :param min_cov: Filtering threshold to use for chr20 coverage.
    :param qc_mt_callrate_ht: Table containing the per sample QC MatrixTable callrate.
    :param min_qc_mt_adj_callrate: Filtering threshold to use for sample callrate
        computed on only predetermined QC variants (predetermined using CCDG
        genomes/exomes, gnomAD v3.1 genomes, and UKB exomes) after ADJ filtering.
    :param sex_ht: Optional Table with sex karyotype annotations. If no Table is
        supplied, no sex filters will be applied.
    :return: Table of hard filtered samples.
    """
    ht = ht.annotate_globals(
        hard_filter_cutoffs=hl.struct(
            max_n_singleton=max_n_singleton,
            max_r_het_hom_var=max_r_het_hom_var,
            min_bases_dp_over_1=min_bases_dp_over_1,
            min_bases_dp_over_20=min_bases_dp_over_20,
            max_chimera=max_chimera,
            max_contamination_estimate=max_contamination_estimate,
        ),
    )
    if min_cov is not None:
        ht = ht.annotate_globals(
            hard_filter_cutoffs=ht.hard_filter_cutoffs.annotate(
                chr_20_dp_threshold=min_cov
            )
        )

    hard_filters = dict()
    sample_qc_metric_hard_filters = dict()

    # Flag samples failing fingerprinting.
    hard_filters["failed_fingerprinting"] = hl.is_defined(
        fingerprinting_failed_ht[ht.key]
    )

    # Flag extreme raw sample QC outliers.
    # Convert tuples to lists so we can find the index of the passed threshold.
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            f"bases_dp_over_{hl.eval(sample_qc_ht.dp_bins[i])}": (
                sample_qc_ht.bases_over_dp_threshold[i]
            )
            for i in range(len(sample_qc_ht.dp_bins))
        },
    )
    sample_qc_struct = sample_qc_ht[ht.key]
    sample_qc_metric_hard_filters["high_n_singleton"] = (
        sample_qc_struct.n_singleton > max_n_singleton
    )
    sample_qc_metric_hard_filters["high_r_het_hom_var"] = (
        sample_qc_struct.r_het_hom_var > max_r_het_hom_var
    )
    sample_qc_metric_hard_filters["low_bases_dp_over_1"] = (
        sample_qc_struct.bases_dp_over_1 < min_bases_dp_over_1
    )
    sample_qc_metric_hard_filters["low_bases_dp_over_20"] = (
        sample_qc_struct.bases_dp_over_20 < min_bases_dp_over_20
    )
    hard_filters["sample_qc_metrics"] = (
        sample_qc_metric_hard_filters["high_n_singleton"]
        | sample_qc_metric_hard_filters["high_r_het_hom_var"]
        | sample_qc_metric_hard_filters["low_bases_dp_over_1"]
        | sample_qc_metric_hard_filters["low_bases_dp_over_20"]
    )

    # Flag samples that fail bam metric thresholds.
    if test:
        project_meta_ht = gnomad_v4_testset_meta.ht()
        # Use the gnomAD v4 test dataset's `rand_sampling_meta` annotation to get the
        # bam metrics needed for hard filtering.
        # This annotation includes all of the metadata for the random samples
        # chosen for the test dataset.
        bam_metrics_struct = project_meta_ht[ht.key].rand_sampling_meta
        bam_metrics_struct = bam_metrics_struct.annotate(
            chimeras_rate=bam_metrics_struct.pct_chimeras
        )
        contamination_struct = hl.read_table(
            get_checkpoint_path("test_gnomad.exomes.contamination")
        )[ht.key]
    else:
        bam_metrics_struct = project_meta_ht[ht.key].bam_metrics
        contamination_struct = contamination_ht[ht.key]

    hard_filters["chimera"] = bam_metrics_struct.chimeras_rate > max_chimera
    hard_filters["contamination"] = (
        contamination_struct.mean_AB_snp_biallelic > max_contamination_estimate
    )

    # Flag low-coverage samples using mean coverage on chromosome 20.
    if min_cov is not None:
        if chr20_mean_dp_ht is None:
            raise ValueError(
                "If a chromosome 20 coverage threshold is supplied, a chr20 mean DP"
                " Table must be supplied too."
            )
        hard_filters["low_coverage"] = chr20_mean_dp_ht[ht.key].chr20_mean_dp < min_cov

    if min_qc_mt_adj_callrate is not None:
        if qc_mt_callrate_ht is None:
            raise ValueError(
                "If a QC MatrixTable adj sample callrate threshold is supplied, a QC "
                "MatrixTable sample callrate Table must be supplied too."
            )
        hard_filters["low_adj_callrate"] = (
            qc_mt_callrate_ht[ht.key].callrate_adj < min_qc_mt_adj_callrate
        )
        ht = ht.annotate_globals(
            hard_filter_cutoffs=ht.hard_filter_cutoffs.annotate(
                min_qc_mt_adj_callrate=min_qc_mt_adj_callrate
            )
        )

    if sex_ht is not None:
        sex_struct = sex_ht[ht.key]
        # Remove samples with ambiguous sex assignments.
        hard_filters["ambiguous_sex"] = sex_struct.sex_karyotype == "ambiguous"
        hard_filters["sex_aneuploidy"] = hl.is_defined(
            sex_struct.sex_karyotype
        ) & ~hl.set(  # pylint: disable=invalid-unary-operand-type
            {"ambiguous", "XX", "XY"}
        ).contains(
            sex_struct.sex_karyotype
        )

    ht = ht.annotate(
        hard_filters=add_filters_expr(filters=hard_filters),
        sample_qc_metric_hard_filters=add_filters_expr(
            filters=sample_qc_metric_hard_filters
        ),
    )

    # Keep samples failing hard filters.
    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    return ht


def get_pipeline_resources(
    test: bool,
    overwrite: bool,
    calling_interval_name: str,
    calling_interval_padding: int,
    include_sex_filter: bool,
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the hard filters pipeline.

    :param test: Whether to gather all resources for the test dataset.
    :param overwrite: Whether to overwrite resources if they exist.
    :param calling_interval_name: Name of calling intervals to use.
    :param calling_interval_padding: Padding to use for calling intervals.
    :param include_sex_filter: Whether sex filters should be included in hard filtering.
    :return: PipelineResourceCollection containing resources for all steps of the
        hard filters pipeline.
    """
    # Initialize hard filter pipeline resource collection.
    hard_filter_pipeline = PipelineResourceCollection(
        pipeline_name="hard_filters",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the relatedness pipeline.
    run_sample_qc = PipelineStepResourceCollection(
        "--sample-qc",
        output_resources={"sample_qc_ht": get_sample_qc(test=test)},
    )
    compute_coverage = PipelineStepResourceCollection(
        "--compute-coverage",
        input_resources={
            "Calling intervals": {
                "interval_ht": calling_intervals(
                    calling_interval_name, calling_interval_padding
                )
            }
        },
        output_resources={"interval_cov_mt": interval_coverage(test=test)},
    )
    compute_contamination_estimate = PipelineStepResourceCollection(
        "--compute-contamination-estimate",
        output_resources={"contamination_ht": contamination(test=test)},
    )
    compute_chr20_mean_dp = PipelineStepResourceCollection(
        "--compute-chr20-mean-dp",
        pipeline_input_steps=[compute_coverage],
        output_resources={"chr20_mean_dp_ht": sample_chr20_mean_dp(test=test)},
    )
    compute_qc_mt_callrate = PipelineStepResourceCollection(
        "--compute-qc-mt-callrate",
        input_resources={
            "generate_qc_mt.py --create-v4-filtered-dense-mt": {
                "v4_predetermined_qc_ht": get_predetermined_qc(test=test)
            }
        },
        output_resources={"sample_qc_mt_callrate_ht": sample_qc_mt_callrate(test=test)},
    )

    compute_hard_filters_ht = PipelineStepResourceCollection(
        "--compute-hard-filters",
        output_resources={
            "hard_filter_ht": hard_filtered_samples(
                include_sex_filter=include_sex_filter, test=test
            )
        },
        pipeline_input_steps=[
            run_sample_qc,
            compute_contamination_estimate,
            compute_chr20_mean_dp,
            compute_qc_mt_callrate,
        ],
        add_input_resources=(
            {"sex_inference.py --annotate-sex-karyotype": {"sex_ht": sex(test=test)}}
            if include_sex_filter
            else None
        ),
    )

    # Add all steps to the relatedness pipeline resource collection.
    hard_filter_pipeline.add_steps(
        {
            "run_sample_qc": run_sample_qc,
            "compute_coverage": compute_coverage,
            "compute_contamination_estimate": compute_contamination_estimate,
            "compute_chr20_mean_dp": compute_chr20_mean_dp,
            "compute_qc_mt_callrate": compute_qc_mt_callrate,
            "compute_hard_filters_ht": compute_hard_filters_ht,
        }
    )

    return hard_filter_pipeline


def main(args):
    """Determine samples that fail hard filtering thresholds."""
    hl.init(
        log="/gnomad_hard_filters.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default.
    hl._set_flags(use_new_shuffle="1")

    calling_interval_name = args.calling_interval_name
    calling_interval_padding = args.calling_interval_padding
    test = args.test
    overwrite = args.overwrite

    hard_filtering_resources = get_pipeline_resources(
        test=test,
        overwrite=overwrite,
        calling_interval_name=calling_interval_name,
        calling_interval_padding=calling_interval_padding,
        include_sex_filter=args.include_sex_filter,
    )

    try:
        if args.sample_qc:
            res = hard_filtering_resources.sample_qc
            res.check_resource_existence()
            vds = get_gnomad_v4_vds(
                split=True, remove_hard_filtered_samples=False, test=test
            )
            compute_sample_qc(
                vds,
                n_partitions=args.sample_qc_n_partitions,
                test=test,
                n_alt_alleles_strata=args.n_alt_alleles_strata,
                n_alt_alleles_strata_name=args.n_alt_alleles_strata_name,
            ).write(res.sample_qc_ht.path, overwrite=overwrite)

        if args.compute_coverage:
            logger.info("Loading v4 VDS...")
            res = hard_filtering_resources.compute_coverage
            res.check_resource_existence()
            vds = get_gnomad_v4_vds(remove_hard_filtered_samples=False, test=test)

            logger.info(
                "Loading calling intervals: %s with padding of %d...",
                calling_interval_name,
                calling_interval_padding,
            )
            mt = hl.vds.interval_coverage(vds, intervals=res.interval_ht.ht())
            mt = mt.annotate_globals(
                calling_interval_name=calling_interval_name,
                calling_interval_padding=calling_interval_padding,
            )
            mt.write(res.interval_cov_mt.path, overwrite=overwrite)

        if args.compute_contamination_estimate:
            logger.info(
                "Loading v4 VDS, filtering to high-quality (DP >= %d), autosomal,"
                " bi-allelic homozygous SNVs and computing the mean of reference allele"
                " balances per sample...",
                args.contam_dp_cutoff,
            )
            res = hard_filtering_resources.compute_contamination_estimate
            res.check_resource_existence()
            mt = filter_to_autosomes(
                get_gnomad_v4_vds(
                    remove_hard_filtered_samples=False, test=test
                ).variant_data
            )
            mt = mt.filter_rows(
                (hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1])
            )
            mt = mt.filter_entries(
                mt.LGT.is_hom_var() & (mt.DP >= args.contam_dp_cutoff)
            )
            mt = mt.annotate_cols(
                mean_AB_snp_biallelic=hl.agg.mean(mt.LAD[0] / (mt.LAD[0] + mt.LAD[1]))
            )
            mt = mt.cols().annotate_globals(dp_cutoff=args.contam_dp_cutoff)
            mt.write(res.contamination_ht.path, overwrite=overwrite)

        if args.compute_chr20_mean_dp:
            # TODO: Add drop of `gq_thresholds` for future versions.
            res = hard_filtering_resources.compute_chr20_mean_dp
            res.check_resource_existence()
            coverage_mt = res.interval_cov_mt.mt()
            coverage_mt = coverage_mt.filter_rows(
                coverage_mt.interval.start.contig == "chr20"
            )
            coverage_mt.select_cols(
                chr20_mean_dp=hl.agg.sum(coverage_mt.sum_dp)
                / hl.agg.sum(coverage_mt.interval_size)
            ).cols().write(res.chr20_mean_dp_ht.path, overwrite=overwrite)

        if args.compute_qc_mt_callrate:
            res = hard_filtering_resources.compute_qc_mt_callrate
            res.check_resource_existence()
            mt = res.v4_predetermined_qc_ht.ht()
            num_samples = mt.count_cols()

            # Filter predetermined QC variants to AF > qc_mt_callrate_min_af
            # with site callrate > qc_mt_callrate_min_site_callrate for ADJ genotypes.
            mt = annotate_adj(mt)
            mt = mt.filter_rows(
                hl.agg.filter(
                    hl.is_defined(mt.GT) & mt.adj,
                    (
                        (hl.agg.count() / num_samples)
                        > args.qc_mt_callrate_min_site_callrate
                    )
                    & (
                        (hl.agg.sum(mt.GT.n_alt_alleles()) / (hl.agg.count() * 2))
                        > args.qc_mt_callrate_min_af
                    ),
                )
            )
            num_variants = mt.count_rows()
            ht = mt.annotate_cols(
                sample_qc_mt_callrate=hl.agg.count_where(hl.is_defined(mt.GT))
                / num_variants,
                sample_qc_mt_callrate_adj=hl.agg.count_where(
                    hl.is_defined(mt.GT) & mt.adj
                )
                / num_variants,
            ).cols()
            ht = ht.annotate_globals(
                min_af=args.qc_mt_callrate_min_af,
                min_site_callrate=args.qc_mt_callrate_min_site_callrate,
            )
            ht.write(res.sample_qc_mt_callrate_ht.path, overwrite=overwrite)

        if args.compute_hard_filters:
            res = hard_filtering_resources.compute_hard_filters_ht
            res.check_resource_existence()
            ht = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False, test=test
            ).variant_data.cols()
            ht = compute_hard_filters(
                ht,
                fingerprinting_failed.ht(),
                get_sample_qc("bi_allelic", test=test).ht(),
                args.include_sex_filter,
                args.max_n_singleton,
                args.max_r_het_hom_var,
                args.min_bases_dp_over_1,
                args.min_bases_dp_over_20,
                args.max_chimera,
                args.max_contamination_estimate,
                test,
                res.chr20_mean_dp_ht.ht(),
                args.min_cov,
                res.sample_qc_mt_callrate_ht.ht(),
                args.min_qc_mt_adj_callrate,
            )
            ht = ht.checkpoint(res.hard_filter_ht.path, overwrite=overwrite)
            ht.group_by("hard_filters").aggregate(n=hl.agg.count()).show(20)
            ht.group_by("sample_qc_metric_hard_filters").aggregate(
                n=hl.agg.count()
            ).show(20)
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("hard_filters"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all Matrixtables/Tables. (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Use the v4 test dataset instead of the full dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--sample-qc", help="Compute Hail's VDS sample QC metrics.", action="store_true"
    )
    parser.add_argument(
        "--n-alt-alleles-strata",
        help=(
            "Number of alternative alleles to stratify sample QC on. For example "
            "'n_alt_alleles_strata' = 2 will stratify 'bi-allelic' and 'multi-allelic'."
            " 'n_alt_alleles_strata' = 3 will stratify under 3 alt alleles and 3 or "
            "more alt alleles. Default is 3."
        ),
        type=int,
        default=3,
    )
    parser.add_argument(
        "--n-alt-alleles-strata-name",
        help=(
            "Name to use for number of alternative allele stratification, where strata "
            "names are 'under_{n_alt_alleles_strata_name}_alt_alleles' and "
            "'{n_alt_alleles_strata_name}_or_more_alt_alleles'. This is not used when "
            "`n_alt_alleles_strata` = 2, and instead uses 'bi-allelic' and "
            "'multi-allelic'. Default is 'three'."
        ),
        type=str,
        default="three",
    )
    parser.add_argument(
        "--sample-qc-n-partitions",
        help="Number of desired partitions for the sample QC output Table.",
        default=500,
        type=int,
    )
    parser.add_argument(
        "--compute-coverage",
        help=(
            "Compute per interval coverage metrics using Hail's vds.interval_coverage"
            " method."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-contamination-estimate",
        help=(
            "Compute contamination estimate as the mean of reference allele balances of"
            " high-quality (DP >= contam-dp-cutoff), autosomal, bi-allelic homozygous"
            " SNVs per sample."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--contam-dp-cutoff",
        help=(
            "Minimum genotype depth to be included in contamination estimate"
            " calculation."
        ),
        type=int,
        default=10,
    )
    parser.add_argument(
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb',"
            " 'broad', or 'intersection'."
        ),
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    parser.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or"
            " 50 bp."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )
    parser.add_argument(
        "--compute-chr20-mean-dp",
        help=(
            "Compute per sample mean DP on chromosome 20 using interval coverage"
            " results."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-qc-mt-callrate",
        help=(
            "Compute per sample callrate on the predetermined QC variants after "
            "filtering to common variants (AF > --qc-mt-callrate-min-af) with high "
            "site callrate ( > --qc-mt-callrate-min-site-callrate) for ADJ genotypes."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--qc-mt-callrate-min-af",
        default=0.0001,
        type=float,
        help=(
            "Filtering threshold to use for minimum variant allele frequency when "
            "computing sample callrate on the predetermined QC variants "
            "(--compute-qc-mt-callrate). Default is 0.0001."
        ),
    )
    parser.add_argument(
        "--qc-mt-callrate-min-site-callrate",
        default=0.99,
        type=float,
        help=(
            "Filtering threshold to use for minimum site ADJ callrate when computing "
            "sample callrate on predetermined QC variants (--compute-qc-mt-callrate). "
            "Default is 0.99."
        ),
    )
    parser.add_argument(
        "--compute-hard-filters",
        help=(
            "Computes samples to be hard-filtered. NOTE: Cutoffs should be determined"
            " by visual inspection of the metrics."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--include-sex-filter",
        help="If sex filters should be included in hard filtering.",
        action="store_true",
    )
    parser.add_argument_group()
    hard_filter_args = parser.add_argument_group(
        "Hard-filter cutoffs", "Arguments used for hard-filter cutoffs."
    )
    hard_filter_args.add_argument(
        "--max-n-singleton",
        type=float,
        default=5000,
        help=(
            "Filtering threshold to use for the maximum number of singletons. Default"
            " is 5000."
        ),
    )
    hard_filter_args.add_argument(
        "--max-r-het-hom-var",
        type=float,
        default=10,
        help=(
            "Filtering threshold to use for the maximum ratio of heterozygotes to"
            " alternate homozygotes. Default is 10."
        ),
    )
    hard_filter_args.add_argument(
        "--min-bases-dp-over-1",
        type=float,
        help=(
            "Filtering threshold to use for the minimum number of bases with a DP over"
            " one. Default is 5e7."
        ),
        default=5e7,
    )
    hard_filter_args.add_argument(
        "--min-bases-dp-over-20",
        type=float,
        help=(
            "Filtering threshold to use for the minimum number of bases with a DP over"
            " 20. Default is 4e7."
        ),
        default=4e7,
    )
    hard_filter_args.add_argument(
        "--max-chimera",
        type=float,
        default=0.05,
        help=(
            "Filtering threshold to use for maximum chimera (this is a proportion not a"
            " percent, e.g. 5% == 0.05, %5 != 5). Default is 0.05."
        ),
    )
    hard_filter_args.add_argument(
        "--max-contamination-estimate",
        default=0.015,
        type=float,
        help=(
            "Filtering threshold to use for maximum contamination estimate (from"
            " --compute-contamination-estimate). Default is 0.015."
        ),
    )
    hard_filter_args.add_argument(
        "--min-cov",
        help=(
            "Minimum chromosome 20 coverage for inclusion when computing hard-filters."
        ),
        default=None,
        type=int,
    )
    hard_filter_args.add_argument(
        "--min-qc-mt-adj-callrate",
        help=(
            "Minimum sample callrate computed on only predetermined QC variants"
            " (predetermined using CCDG genomes/exomes, gnomAD v3.1 genomes, and UKB"
            " exomes) after ADJ filtering."
        ),
        default=None,
        type=float,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
