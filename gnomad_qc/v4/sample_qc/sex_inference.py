import argparse
import logging
from typing import Callable, List, Optional

import hail as hl

from gnomad.sample_qc.pipeline import annotate_sex, infer_sex_karyotype
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    calling_intervals,
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
    ukb_f_stat,
)
from gnomad_qc.v4.resources.sample_qc import (
    f_stat_sites,
    platform,
    ploidy,
    sex,
    sex_imputation_coverage,
    sex_imputation_platform_coverage,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sex_inference")
logger.setLevel(logging.INFO)


def determine_fstat_sites(
    vds: hl.vds.VariantDataset,
    approx_af_and_no_callrate: bool = False,
    min_af: float = 0.001,
    min_callrate: float = 0.99,
) -> hl.Table:
    """
    Write a Table with chromosome X SNPs that are bi-allelic, common, and high callrate by default.

    This Table is designed to be used as a variant filter in sex imputation for f-stat computation.

    .. warning::

        By default `approx_af_and_no_callrate` is False and the final Table will be filtered to high callrate (> value
        specified by `min_callrate`) variants. This requires a densify of chrX!"

    .. note::

        If `approx_af_and_no_callrate` is True, allele frequency is approximated with AC/(n_samples * 2) and no callrate
        filter is used.

    :param vds: Input VariantDataset.
    :param approx_af_and_no_callrate: Whether to approximate allele frequency with AC/(n_samples * 2) and use no
        callrate cutoff to filter sites.
    :param min_af: Minimum alternate allele frequency cutoff used to filter sites.
    :param min_callrate: Minimum callrate cutoff used to filter sites.
    :return: Table of chromosome X sites to be used for f-stat computation.
    """
    vds = hl.vds.filter_chromosomes(vds, keep=["chrX"])
    vd = vds.variant_data
    vd = vd.filter_rows(
        (hl.len(vd.alleles) == 2) & hl.is_snp(vd.alleles[0], vd.alleles[1])
    )
    vd = vd.transmute_entries(GT=hl.experimental.lgt_to_gt(vd.LGT, vd.LA))

    if approx_af_and_no_callrate:
        n_samples = vd.count_cols()
        logger.info("Number of samples: %d", n_samples)
        ht = vd.select_rows(
            AF=hl.agg.sum(vd.GT.n_alt_alleles()) / (n_samples * 2),
        ).rows()
        ht = ht.filter(ht.AF > min_af)
        ht = ht.annotate_globals(
            min_approx_af=min_af,
        )
    else:
        mt = hl.vds.to_dense_mt(hl.vds.VariantDataset(vds.reference_data, vd))
        ht = hl.variant_qc(mt).rows()
        ht = ht.filter(
            (ht.variant_qc.call_rate > min_callrate) & (ht.variant_qc.AF[1] > min_af)
        )
        ht = ht.annotate(AF=ht.variant_qc.AF[1])
        ht = ht.annotate_globals(
            min_af=min_af,
            min_callrate=min_callrate,
        )

    return ht


def load_platform_ht(
    test: bool = False,
    calling_interval_name: str = "intersection",
    calling_interval_padding: int = 50,
) -> hl.Table:
    """
    Load platform assignment Table or test Table and return an error if requested Table does not exist.

    .. note::

        If `test` is True and the test platform assignment Table does not exist, the function will load the final
        platform assignment Table instead if it already exists.

    :param test: Whether a test platform assignment Table should be loaded.
    :param calling_interval_name: Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or
        'intersection'. Only used if `test` is True.
    :param calling_interval_padding: Number of base pair padding to use on the calling intervals. One of 0 or 50 bp.
        Only used if `test` is True.
    :return: Platform assignment Table.
    """
    logger.info("Loading platform information...")
    test_platform_path = get_checkpoint_path(
        f"test_platform_assignment.{calling_interval_name}.pad{calling_interval_padding}"
    )
    if test and file_exists(test_platform_path):
        ht = hl.read_table(test_platform_path)
    elif file_exists(platform.path):
        ht = platform.ht()
        if test:
            logger.warning(
                "Test platform file does not exist for calling interval %s and interval padding %s, using final "
                "platform assignment Table instead. To use a test platform assignment please run "
                "platform_inference.py --assign-platforms with the --test argument and needed "
                "--calling-interval-name/--calling-interval-padding arguments.",
                calling_interval_name,
                calling_interval_padding,
            )
    elif test:
        raise FileNotFoundError(
            f"There is no test platform assignment Table written for calling interval {calling_interval_name} and "
            f"interval padding {calling_interval_padding} and a final platform assignment Table does not exist. "
            f"Please run platform_inference.py --assign-platforms with the --test argument and needed "
            f"--calling-interval-name/--calling-interval-padding arguments."
        )
    else:
        raise FileNotFoundError(
            f"There is no final platform assignment Table written. Please run: "
            f"platform_inference.py --assign-platforms to compute the platform assignment Table."
        )

    return ht


def generate_sex_imputation_interval_coverage_mt(
    vds: hl.vds.VariantDataset,
    calling_intervals_ht: hl.Table,
    contigs: List[str] = ["chrX", "chrY", "chr20"],
) -> hl.MatrixTable:
    """
    Create a MatrixTable of interval-by-sample coverage on a specified list of contigs with PAR regions excluded.

    :param vds: Input VariantDataset.
    :param calling_intervals_ht: Calling interval Table.
    :param contigs: Which contigs to compute interval coverage on for sex imputation. Default: 'chrX', 'chrY', 'chr20'.
    :return: MatrixTable with interval coverage per sample on specified contigs and PAR regions excluded.
    """
    calling_intervals_ht = calling_intervals_ht.filter(
        hl.literal(contigs).contains(calling_intervals_ht.interval.start.contig)
    )
    logger.info(
        "Filtering VariantDataset to the following contigs: %s...",
        ", ".join(contigs),
    )
    vds = hl.vds.filter_chromosomes(vds, keep=contigs)
    rg = vds.reference_data.locus.dtype.reference_genome

    par_boundaries = []
    for par_interval in rg.par:
        par_boundaries.append(par_interval.start)
        par_boundaries.append(par_interval.end)

    # Segment on PAR interval boundaries
    calling_intervals = hl.segment_intervals(calling_intervals_ht, par_boundaries)

    # Remove intervals overlapping PAR
    calling_intervals = calling_intervals.filter(
        hl.all(lambda x: ~x.overlaps(calling_intervals.interval), hl.literal(rg.par))
    )

    kept_contig_filter = hl.array(contigs).map(
        lambda x: hl.parse_locus_interval(x, reference_genome=rg)
    )
    vds = hl.vds.VariantDataset(
        hl.filter_intervals(vds.reference_data, kept_contig_filter),
        hl.filter_intervals(vds.variant_data, kept_contig_filter),
    )
    mt = hl.vds.interval_coverage(vds, calling_intervals, gq_thresholds=()).drop(
        "gq_thresholds"
    )

    return mt


def generate_sex_imputation_interval_qc_mt(
    mt: hl.MatrixTable,
    platform_ht: hl.Table,
    mean_dp_thresholds: List[int] = [5, 10, 15, 20, 25],
) -> hl.MatrixTable:
    """
    Create a Table of fraction of samples per interval and per platform with mean DP over specified thresholds.

    :param mt: Input sex interval coverage MatrixTable.
    :param platform_ht: Input platform assignment Table.
    :param mean_dp_thresholds: List of mean DP thresholds to use for computing the fraction of samples with mean
        interval DP >= the threshold.
    :return: MatrixTable with annotations for the fraction of samples per interval and per platform over DP thresholds.
    """
    # Note: Default `hl.vds.interval_coverage` will return a list for `fraction_over_dp_threshold` where the
    # second element is dp >= 1 (dp > 0)
    mt = mt.annotate_cols(platform=platform_ht[mt.col_key].qc_platform)
    mt = mt.annotate_rows(
        **{
            f"fraction_over_{dp}x": hl.agg.fraction(mt.mean_dp >= dp)
            for dp in mean_dp_thresholds
        },
        mean_fraction_over_dp_0=hl.agg.mean(mt.fraction_over_dp_threshold[1]),
    )
    mt = mt.select_globals(mean_dp_thresholds=mean_dp_thresholds)

    logger.info("Adding per platform aggregation...")
    platform_mt = mt.group_cols_by(mt.platform).aggregate(
        **{
            f"platform_fraction_over_{dp}x": hl.agg.fraction(mt.mean_dp >= dp)
            for dp in mean_dp_thresholds
        },
        platform_mean_fraction_over_dp_0=hl.agg.mean(mt.fraction_over_dp_threshold[1]),
    )

    platform_ht = platform_ht.group_by(platform_ht.qc_platform).aggregate(
        n_samples=hl.agg.count()
    )
    platform_mt = platform_mt.annotate_cols(
        n_samples=platform_ht[platform_mt.col_key].n_samples
    )

    return platform_mt


def compute_sex_ploidy(
    vds: hl.vds.VariantDataset,
    interval_qc_mt: Optional[hl.MatrixTable] = None,
    high_cov_intervals: bool = False,
    high_cov_per_platform: bool = False,
    high_cov_all_platforms: bool = False,
    platform_ht: Optional[hl.Table] = None,
    min_platform_size: bool = 100,
    normalization_contig: str = "chr20",
    variant_depth_only_x_ploidy: bool = False,
    variant_depth_only_y_ploidy: bool = False,
    x_cov: int = None,
    y_cov: int = None,
    norm_cov: int = None,
    prop_samples_x: float = None,
    prop_samples_y: float = None,
    prop_samples_norm: float = None,
    freq_ht: Optional[hl.Table] = None,
    min_af: float = 0.001,
    f_stat_cutoff: float = -1.0,
) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Within this function there are some critical parameters to consider (more details on each are described below):
        - Filtering to intervals with high coverage for sex chromosome ploidy estimation (`high_cov_intervals`, `x_cov`,
         `y_cov`, `norm_cov`, `prop_samples_x`, `prop_samples_y`, and `prop_samples_norm`).
        - Per platform computations: if `high_cov_per_platform` is used, per platform sex chromosome ploidy estimation
         is performed using per platform high coverage intervals (`x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
         `prop_samples_y`, `prop_samples_norm`, and `min_platform_size`).
        - Use per platform stats to determine which are high coverage across all platforms (depending on platform sample
         size). Uses parameters `high_cov_all_platforms`, `x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
         `prop_samples_y`, `prop_samples_norm`.
        - We use `annotate_sex`, which uses Hail's VDS imputation of sex ploidy `hail.vds.impute_sex_chromosome_ploidy`,
          and uses `variant_depth_only_x_ploidy` and `variant_depth_only_y_ploidy`.
        - Use of a `freq_ht` to filter variants for the X-chromosome heterozygosity computation. This is the f_stat
         annotation applied by Hail's `impute_sex` module. (`freq_ht`, `min_af`).

    The following are the available options for chrX and chrY relative sex ploidy estimation:
        - Ploidy estimation using all intervals found in `interval_qc_mt`.
        - Per platform ploidy estimation using all intervals found in `interval_qc_mt`.
        - Ploidy estimation using only high coverage intervals across the entire dataset, defined as: chrX intervals
         having a proportion of all samples (`prop_samples_x`) with over a specified mean DP (`x_cov`), chrY intervals
         having a proportion of all samples (`prop_samples_y`) with over a specified mean DP (`y_cov`), and
         `normalization_contig` intervals having a proportion of all samples (`prop_samples_norm`) with over a
         specified mean DP (`norm_cov`).
        - Per platform ploidy estimation using only per platform high coverage intervals defined by: chrX intervals with
         the specific platform having a proportion of samples (`prop_samples_x`) with over a specified mean DP (`x_cov`),
         chrY intervals with the specific platform having a proportion of samples (`prop_samples_y`) with over a
         specified mean DP (`y_cov`), and `normalization_contig` intervals with the specific platform having a
         proportion of samples (`prop_samples_norm`) with over a specified mean DP (`norm_cov`).
        - A combined ploidy estimation using only per platform high coverage intervals. Each interval is assessed as high
         coverage per platform on chrX, chrY, and `normalization_contig` as described above for "per platform ploidy
         estimation". Then this method uses only platforms that have # of samples > `min_platform_size` to determine
         intervals that have a high coverage across all platforms.

    For each of the options described above, there is also the possibility to use a ploidy estimation that uses only
    variants within the specified calling intervals:
        - This can be defined differently for chrX and chrY using `variant_depth_only_x_ploidy` and
        `variant_depth_only_y_ploidy`.
        - `variant_depth_only_x_ploidy` is the preferred method for chrX ploidy computation.
        - If not using only variants Hail's `hail.vds.impute_sex_chromosome_ploidy` method will only compute chromosome
         ploidy using reference block DP per calling interval. This method breaks up the reference blocks at the
         calling interval boundries, maintaining all reference block END information for the mean DP per interval
         computation. This is different from the method used on sparse MatrixTables in gnomad_methods
         `gnomad.utils.sparse_mt.impute_sex_ploidy`. In this the chromosome coverage will be computed using all
         reference blocks and variants found in the sparse MatrixTable, but it only uses specified calling intervals to
         determine the contig size and doesn't break up the reference blocks in the same way the Hail method does.

    :param vds: Input VDS for use in sex inference.
    :param interval_qc_mt: Optional interval QC MatrixTable. This is only needed if `high_cov_intervals`,
        `high_cov_per_platform` or `high_cov_all_platforms` are True.
    :param high_cov_intervals: Whether to filter to high coverage intervals for the sex ploidy imputation. Default
        is False.
    :param high_cov_per_platform: Whether to filter to per platform high coverage intervals for the sex ploidy imputation.
        Default is False.
    :param high_cov_all_platforms: Whether to filter to high coverage intervals for the sex ploidy imputation. Using
        only intervals that are considered high coverage across all platforms. Default is False.
    :param platform_ht: Input platform assignment Table. This is only needed if `high_cov_per_platform` or
        `high_cov_all_platforms` are True.
    :param min_platform_size: Required size of a platform to be considered when using `high_cov_all_platforms`. Only
        platforms that have # of samples > 'min_platform_size' are used to determine intervals that have a high
        coverage across all platforms. Default is 100.
    :param normalization_contig: Which autosomal chromosome to use for normalizing the coverage of chromosomes X and Y.
        Default is 'chr20'.
    :param variant_depth_only_x_ploidy: Whether to use depth of variant data within calling intervals instead of
        reference data for chrX ploidy estimation. Default will only use reference data.
    :param variant_depth_only_y_ploidy: Whether to use depth of variant data within calling intervals instead of
        reference data for chrY ploidy estimation. Default will only use reference data.
    :param x_cov: Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the
        interval_coverage MatrixTable.
    :param y_cov: Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the
        interval_coverage MatrixTable.
    :param norm_cov: Mean coverage level used to define high coverage intervals on the normalization autosome
        (`normalization_contig`). This field must be in the interval_coverage MT.
    :param prop_samples_x: Proportion samples at specified coverage `x_cov` to determine high coverage intervals on
        chromosome X.
    :param prop_samples_y: Proportion samples at specified coverage `y_cov` to determine high coverage intervals on
        chromosome Y.
    :param prop_samples_norm: Proportion samples at specified coverage `norm_cov` to determine high coverage intervals
        on the normalization chromosome specified by `normalization_contig`.
    :param freq_ht: Optional Table to use for f-stat allele frequency cutoff. The input VDS is filtered to sites in
        this Table prior to running Hail's `impute_sex` module, and alternate allele frequency is used from this Table
        with a `min_af` cutoff.
    :param min_af: Minimum alternate allele frequency to be used in f-stat calculations. Default is 0.001.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY
        are above cutoff. Default is -1.0.
    :return: Table with imputed ploidies.
    """
    if (high_cov_per_platform or high_cov_all_platforms) and platform_ht is None:
        raise ValueError(
            "'platform_ht' must be defined if 'high_cov_per_platform' or 'high_cov_all_platforms' is True!"
        )

    if high_cov_per_platform and high_cov_all_platforms:
        raise ValueError(
            "Only one of 'high_cov_per_platform' or 'high_cov_all_platforms' can be True!"
        )

    if high_cov_intervals and high_cov_all_platforms:
        logger.warning(
            "Both 'high_cov_intervals' and 'high_cov_all_platforms' are True, high coverage intervals will be defined"
            " as intervals that match high coverage criteria within all platforms."
        )
    add_globals = {}

    def _get_high_coverage_intervals_ht(
        interval_qc_mt: hl.MatrixTable,
        prefix: str = "",
        agg_func: Callable[
            [hl.expr.BooleanExpression], hl.BooleanExpression
        ] = hl.agg.all,
    ) -> hl.Table:
        """
        Helper function to create a Table filtered to high coverage intervals.

        High coverage intervals are determined using `agg_func`, `x_cov`, `y_cov`, `norm_cov`, `prop_samples_x`,
        `prop_samples_y`, and `prop_samples_norm`.

        :param interval_qc_mt: Input interval QC MatrixTable.
        :param prefix: Prefix of annotations in `interval_qc_mt` that contain the proportion of samples with
            mean DP over coverage cutoffs.
        :param agg_func: Hail aggregation function to determine if an interval coverage meets the
            `cov_*` > `prop_samples_*` criteria.
        :return: Table of high coverage intervals.
        """
        return interval_qc_mt.filter_rows(
            (
                (interval_qc_mt.interval.start.contig == "chrX")
                & agg_func(interval_qc_mt[f"{prefix}fraction_over_{x_cov}x"] > prop_samples_x)
            )
            | (
                (interval_qc_mt.interval.start.contig == "chrY")
                & agg_func(interval_qc_mt[f"{prefix}fraction_over_{y_cov}x"] > prop_samples_y)
            )
            | (
                (interval_qc_mt.interval.start.contig == normalization_contig)
                & agg_func(
                    interval_qc_mt[f"{prefix}fraction_over_{norm_cov}x"] > prop_samples_norm
                )
            )
        ).rows()

    def _annotate_sex(
        vds: hl.vds.VariantDataset, calling_intervals_ht: hl.Table
    ) -> hl.Table:
        """
        Helper function to perform `annotate_sex` using unchanged parameters with changes to the VDS and calling
        intervals.

        :param vds: Input VDS to use for sex ploidy annotation.
        :param calling_intervals_ht: Table including only intervals wanted for sex annotation.
        :return: Table containing sex ploidy estimates for samples in the input VDS.
        """
        ploidy_ht = annotate_sex(
            vds,
            included_intervals=calling_intervals_ht,
            normalization_contig=normalization_contig,
            sites_ht=freq_ht.filter(hl.is_defined(calling_intervals_ht[freq_ht.locus]))
            if freq_ht is not None
            else None,
            aaf_expr="AF",
            gt_expr="LGT",
            f_stat_cutoff=f_stat_cutoff,
            aaf_threshold=min_af,
            variants_only_x_ploidy=variant_depth_only_x_ploidy,
            variants_only_y_ploidy=variant_depth_only_y_ploidy,
            infer_karyotype=False,
        )
        return ploidy_ht

    if high_cov_intervals or high_cov_all_platforms or high_cov_per_platform:
        logger.info(
            "Running sex ploidy estimation using only high coverage intervals: %d percent of samples "
            "with greater than %dx coverage on chrX, %d percent of samples with greater than %dx coverage on chrY, and "
            "%d percent of samples with greater than %dx coverage on %s...",
            int(prop_samples_x * 100),
            x_cov,
            int(prop_samples_y * 100),
            y_cov,
            int(prop_samples_norm * 100),
            norm_cov,
            normalization_contig,
        )

    if high_cov_per_platform:
        logger.info("Running sex ploidy imputation per platform using per platform high coverage intervals...")
        platforms = platform_ht.aggregate(
            hl.agg.collect_as_set(platform_ht.qc_platform)
        )
        per_platform_ploidy_hts = []
        for platform in platforms:
            logger.info(
                "Performing ploidy imputation using high coverage intervals for platform %s...",
                platform,
            )
            ploidy_ht = _annotate_sex(
                hl.vds.filter_samples(
                    vds, platform_ht.filter(platform_ht.qc_platform == platform)
                ),
                _get_high_coverage_intervals_ht(
                    interval_qc_mt.filter_cols(interval_qc_mt.platform == platform),
                    prefix="platform_",
                ),
            )
            per_platform_ploidy_hts.append(ploidy_ht)

        ploidy_ht = per_platform_ploidy_hts[0].union(*per_platform_ploidy_hts[1:])
    elif high_cov_all_platforms:
        logger.info(
            "Running sex ploidy imputation using high coverage intervals across all platforms. Limited to platforms "
            "with at least %s samples...",
            min_platform_size,
        )
        # Excluding small platforms and platform_-1 (platform containing all samples with unassigned platform)
        interval_qc_mt = interval_qc_mt.filter_cols(
            (interval_qc_mt.n_samples >= min_platform_size)
            & (interval_qc_mt.platform != "platform_-1")
        )
        ploidy_ht = _annotate_sex(vds, _get_high_coverage_intervals_ht(interval_qc_mt, prefix="platform_"))
        add_globals["high_cov_all_platforms_min_platform_size"] = min_platform_size
    elif high_cov_intervals:
        logger.info("Running sex ploidy imputation using high coverage intervals across the full sample set...")
        ploidy_ht = _annotate_sex(vds, _get_high_coverage_intervals_ht(interval_qc_mt, agg_func=lambda x: x))
    else:
        logger.info("Running sex ploidy imputation...")
        ploidy_ht = _annotate_sex(vds, interval_qc_mt.rows())

    ploidy_ht = ploidy_ht.annotate_globals(
        high_cov_intervals=high_cov_intervals,
        high_cov_per_platform=high_cov_per_platform,
        high_cov_all_platforms=high_cov_all_platforms,
        x_cov=x_cov,
        y_cov=y_cov,
        norm_cov=norm_cov,
        prop_samples_x=prop_samples_x,
        prop_samples_y=prop_samples_y,
        prop_samples_norm=prop_samples_norm,
        f_stat_min_af=min_af,
        f_stat_cutoff=f_stat_cutoff,
        **add_globals,
    )

    return ploidy_ht


def infer_sex_karyotype_from_ploidy(
    ploidy_ht: hl.Table,
    per_platform: bool = False,
    f_stat_cutoff: float = -1.0,
) -> hl.Table:
    """
    Create a Table with X_karyotype, Y_karyotype, and sex_karyotype.

    :param ploidy_ht: Table with chromosome X and chromosome Y ploidies, and f-stat.
    :param per_platform: Whether the sex karyotype ploidy cutoff inference should be applied per platform.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY
        are above cutoff.
    :return: Table of imputed sex karyotypes.
    """
    logger.info("Running sex karyotype inference")
    if per_platform:
        platforms = ploidy_ht.aggregate(hl.agg.collect_as_set(ploidy_ht.platform))
        per_platform_karyotype_hts = []
        x_ploidy_cutoffs = {}
        y_ploidy_cutoffs = {}

        for platform in platforms:
            logger.info(
                "Performing sex karyotype inference for platform %s...",
                platform,
            )
            karyotype_ht = infer_sex_karyotype(
                ploidy_ht.filter(ploidy_ht.platform == platform),
                f_stat_cutoff,
            )
            per_platform_karyotype_hts.append(karyotype_ht)
            x_ploidy_cutoffs[platform] = karyotype_ht.index_globals().x_ploidy_cutoffs
            y_ploidy_cutoffs[platform] = karyotype_ht.index_globals().y_ploidy_cutoffs

        karyotype_ht = per_platform_karyotype_hts[0].union(
            *per_platform_karyotype_hts[1:]
        )
        karyotype_ht = karyotype_ht.annotate_globals(
            x_ploidy_cutoffs=hl.struct(**x_ploidy_cutoffs),
            y_ploidy_cutoffs=hl.struct(**y_ploidy_cutoffs),
        )
    else:
        karyotype_ht = infer_sex_karyotype(ploidy_ht, f_stat_cutoff)

    return karyotype_ht


def main(args):
    hl.init(
        log="/sex_inference.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default
    hl._set_flags(use_new_shuffle="1")

    test = args.test
    calling_interval_name = args.calling_interval_name
    calling_interval_padding = args.calling_interval_padding
    normalization_contig = args.normalization_contig
    per_platform = args.per_platform
    overwrite = args.overwrite

    try:
        if args.determine_fstat_sites:
            logger.info("Determining sites to use for f-stat computations...")
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            ht = determine_fstat_sites(
                vds,
                approx_af_and_no_callrate=args.approx_af_and_no_callrate,
                min_af=args.min_af,
                min_callrate=args.min_callrate,
            )
            ht.naive_coalesce(args.fstat_n_partitions).write(
                get_checkpoint_path("test_f_stat_sites")
                if test
                else f_stat_sites.path,
                overwrite=overwrite,
            )

        if args.sex_imputation_interval_coverage:
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            calling_intervals_ht = calling_intervals(
                calling_interval_name, calling_interval_padding
            ).ht()
            coverage_mt = generate_sex_imputation_interval_coverage_mt(
                vds,
                calling_intervals_ht,
                contigs=["chrX", "chrY", normalization_contig],
            )
            coverage_mt = coverage_mt.annotate_globals(
                calling_interval_name=calling_interval_name,
                calling_interval_padding=calling_interval_padding,
                normalization_contig=normalization_contig,
            )
            coverage_mt.write(
                get_checkpoint_path("test_sex_imputation_cov", mt=True)
                if test
                else sex_imputation_coverage.path,
                overwrite=overwrite,
            )

        if args.sex_imputation_interval_qc:
            if test:
                coverage_mt = hl.read_matrix_table(
                    get_checkpoint_path("test_sex_imputation_cov", mt=True)
                )
            else:
                coverage_mt = sex_imputation_coverage.mt()

            platform_ht = load_platform_ht(
                test,
                coverage_mt.calling_interval_name,
                coverage_mt.calling_interval_padding,
            )
            platform_mt = generate_sex_imputation_interval_qc_mt(
                coverage_mt,
                platform_ht,
                mean_dp_thresholds=args.mean_dp_thresholds,
            )
            platform_mt.naive_coalesce(args.interval_qc_n_partitions).write(
                get_checkpoint_path("test_sex_imputation_cov.per_platform", mt=True)
                if test
                else sex_imputation_platform_coverage.path,
                overwrite=overwrite,
            )

        if args.impute_sex_ploidy:
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            platform_ht = load_platform_ht(
                test, calling_interval_name, calling_interval_padding
            )
            if args.f_stat_ukb_var:
                # The UK Biobank f-stat table contains only variants that were high callrate (0.99) and common
                # (AF >0.001) within the UK Biobank 200K Regeneron exome dataset and it includes the UK Biobank 200K
                # allele frequency information that can be used in the hl.impute_sex f-stat computation allele frequency
                # cutoff (args.min-af)
                freq_ht = ukb_f_stat.ht()
            else:
                freq_ht = (
                    hl.read_table(get_checkpoint_path("test_f_stat_sites"))
                    if test
                    else f_stat_sites.ht()
                )

            ploidy_ht_path = (
                get_checkpoint_path(f"ploidy_imputation") if test else ploidy.path
            )
            # Added because without this impute_sex_chromosome_ploidy will still run even with overwrite=False
            if overwrite or not file_exists(ploidy_ht_path):
                interval_qc_mt = (
                    hl.read_matrix_table(
                        get_checkpoint_path(
                            "test_sex_imputation_cov.per_platform", mt=True
                        )
                    )
                    if test
                    else sex_imputation_platform_coverage.mt()
                )

                ploidy_ht = compute_sex_ploidy(
                    vds,
                    interval_qc_mt,
                    high_cov_intervals=args.high_cov_intervals,
                    high_cov_per_platform=args.high_cov_per_platform,
                    high_cov_all_platforms=args.high_cov_all_platforms,
                    platform_ht=platform_ht,
                    min_platform_size=args.min_platform_size,
                    normalization_contig=normalization_contig,
                    variant_depth_only_x_ploidy=args.variant_depth_only_x_ploidy,
                    variant_depth_only_y_ploidy=args.variant_depth_only_y_ploidy,
                    x_cov=args.x_cov,
                    y_cov=args.y_cov,
                    norm_cov=args.norm_cov,
                    prop_samples_x=args.prop_samples_x,
                    prop_samples_y=args.prop_samples_y,
                    prop_samples_norm=args.prop_samples_norm,
                    freq_ht=freq_ht,
                    min_af=args.min_af,
                    f_stat_cutoff=args.f_stat_cutoff,
                )

                ploidy_ht = ploidy_ht.annotate(platform=platform_ht[ploidy_ht.key].qc_platform)
                ploidy_ht = ploidy_ht.annotate_globals(f_stat_ukb_var=args.f_stat_ukb_var)
                logger.info("Writing ploidy Table...")
                ploidy_ht.write(ploidy_ht_path, overwrite=True)
            else:
                logger.warning("File exists and overwrite is not set!")

        if args.annotate_sex_karyotype:
            ploidy_ht = (
                hl.read_table(get_checkpoint_path(f"ploidy_imputation"))
                if test
                else ploidy.ht()
            )
            karyotype_ht = infer_sex_karyotype_from_ploidy(
                ploidy_ht,
                per_platform=per_platform,
                f_stat_cutoff=args.f_stat_cutoff,
            )
            sex_ht = ploidy_ht.annotate(**karyotype_ht[ploidy_ht.key])
            sex_ht = sex_ht.annotate_globals(**karyotype_ht.index_globals())

            logger.info("Writing sex HT with karyotype annotation...")
            sex_ht.write(
                get_checkpoint_path("sex") if test else sex.path,
                overwrite=overwrite,
            )

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("sex_inference"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Test the pipeline using the gnomAD v4 test dataset.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    fstat_args = parser.add_argument_group(
        "Determine f-stat sites",
        "Arguments used for determining sites to use for f-stat calculations.",
    )
    fstat_args.add_argument(
        "--determine-fstat-sites",
        help=(
            "Create Table of common (> value specified by '--min-af'), bi-allelic SNPs on chromosome X for f-stat "
            "calculations. Additionally filter to high callrate (> value specified by '--min-callrate') variants "
            "if '--approx-af-and-no-callrate' is not used. NOTE: This requires a densify of chrX!"
        ),
        action="store_true",
    )
    fstat_args.add_argument(
        "--min-callrate", help="Minimum variant callrate.", default=0.99, type=float
    )
    fstat_args.add_argument(
        "--approx-af-and-no-callrate",
        help=(
            "Whether to approximate allele frequency with AC/(n_samples * 2) and use no callrate cutoff for "
            "determination of f-stat sites."
        ),
        action="store_true",
    )
    fstat_args.add_argument(
        "--fstat-n-partitions",
        help="Number of desired partitions for the f-stat sites output Table.",
        default=1000,
        type=int,
    )

    sex_coverage_args = parser.add_argument_group(
        "Sex imputation interval coverage",
        "Arguments used for computing interval coverage for sex imputation.",
    )
    sex_coverage_args.add_argument(
        "--sex-imputation-interval-coverage",
        help=(
            "Create a MatrixTable of interval-by-sample coverage on a specified list of contigs with PAR regions "
            "excluded."
        ),
        action="store_true",
    )
    sex_coverage_args.add_argument(
        "--normalization-contig",
        help="Which autosomal chromosome to use for normalizing the coverage of chromosomes X and Y.",
        type=str,
        default="chr20",
    )
    sex_coverage_args.add_argument(
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'. Only "
            "used if '--test' is set."
        ),
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    sex_coverage_args.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or 50 bp. Only used if '--test' is "
            "set."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )

    sex_interval_qc_args = parser.add_argument_group(
        "Sex chromosome interval QC",
        "Arguments used for making an interval QC HT from the sex imputation interval coverage MT.",
    )
    sex_interval_qc_args.add_argument(
        "--sex-imputation-interval-qc",
        help=(
            "Create a Table of the fraction of samples per interval and per platform with mean DP over thresholds "
            "specified by '--mean-dp-thresholds'."
        ),
        action="store_true",
    )
    sex_interval_qc_args.add_argument(
        "--mean-dp-thresholds",
        help=(
            "List of mean DP cutoffs to determining the fraction of samples with mean coverage >= the cutoff for each "
            "interval."
        ),
        type=int,
        nargs="+",
        default=[5, 10, 15, 20, 25],
    )
    sex_interval_qc_args.add_argument(
        "--interval-qc-n-partitions",
        help="Number of desired partitions for the sex imputation interval QC output Table.",
        default=500,
        type=int,
    )

    sex_ploidy_args = parser.add_argument_group(
        "Impute sex ploidy", "Arguments used for imputing sex chromosome ploidy."
    )
    sex_ploidy_args.add_argument(
        "--impute-sex-ploidy",
        help="Run sex chromosome ploidy imputation.",
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--f-stat-ukb-var",
        help=(
            "Whether to use UK Biobank high callrate (0.99) and common variants (UKB allele frequency > value specified"
            " by '--min-af') for f-stat computation instead of the sites determined by '--determine-fstat-sites'."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--min-af",
        help="Minimum variant allele frequency to retain variant.",
        default=0.001,
        type=float,
    )
    sex_ploidy_args.add_argument(
        "--f-stat-cutoff",
        help=(
            "Cutoff for f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY "
            "are above cutoff."
        ),
        type=float,
        default=-1.0,
    )
    sex_ploidy_high_cov_method_parser = sex_ploidy_args.add_mutually_exclusive_group(
        required=False
    )
    sex_ploidy_high_cov_method_parser.add_argument(
        "--high-cov-intervals",
        help=(
            "Whether to filter to high coverage intervals for the sex ploidy imputation. Can't be used at the same "
            "time as '--high-cov-per-platform' or '--high-cov-all-platforms'."
        ),
        action="store_true",
    )
    sex_ploidy_high_cov_method_parser.add_argument(
        "--high-cov-per-platform",
        help=(
            "Whether to filter to per platform high coverage intervals for the sex ploidy imputation. Can't be used "
             "at the same time as '--high-cov-intervals' or '--high-cov-all-platforms'."
        ),
        action="store_true",
    )
    sex_ploidy_high_cov_method_parser.add_argument(
        "--high-cov-all-platforms",
        help=(
            "Whether to filter to high coverage intervals for the sex ploidy imputation. Use only intervals that are "
            "considered high coverage across all platforms.  Can't be used at the same time as '--high-cov-intervals' "
            "or '--high-cov-per-platform'"
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--min-platform-size",
        help=(
            "Required size of a platform to be considered in '--high-cov-all-platforms'. Only platforms that "
            "have # of samples > 'min_platform_size' are used to determine intervals that have a high coverage across "
            "all platforms."
        ),
        type=int,
        default=100,
    )
    sex_ploidy_args.add_argument(
        "--variant-depth-only-x-ploidy",
        help=(
            "Whether to use depth of variant data for the x ploidy estimation instead of the default behavior that "
            "will use reference blocks."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--variant-depth-only-y-ploidy",
        help=(
            "Whether to use depth of variant data for the y ploidy estimation instead of the default behavior that "
            "will use reference blocks."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--x-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on chromosome X. This field must be in the "
            "interval_coverage MT!"
        ),
        type=int,
        default=10,
    )
    sex_ploidy_args.add_argument(
        "--y-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on chromosome Y. This field must be in the "
            "interval_coverage MT!"
        ),
        type=int,
        default=5,
    )
    sex_ploidy_args.add_argument(
        "--norm-cov",
        help=(
            "Mean coverage level used to define high coverage intervals on the normalization autosome. This field must "
            "be in the interval_coverage MT!"
        ),
        type=int,
        default=20,
    )
    sex_ploidy_args.add_argument(
        "--prop-samples-x",
        help="Proportion samples at specified coverage '--x-cov' to determine high coverage intervals on chromosome X.",
        type=float,
        default=0.80,
    )
    sex_ploidy_args.add_argument(
        "--prop-samples-y",
        help="Proportion samples at specified coverage '--y-cov' to determine high coverage intervals on chromosome Y.",
        type=float,
        default=0.35,
    )
    sex_ploidy_args.add_argument(
        "--prop-samples-norm",
        help=(
            "Proportion samples at specified coverage '--norm-cov' to determine high coverage intervals on the "
            "normalization chromosome specified by '--normalization-contig'."
        ),
        type=float,
        default=0.85,
    )
    sex_ploidy_args.add_argument(
        "--calling-interval-name",
        help=(
            "Name of calling intervals to use for interval coverage. One of: 'ukb', 'broad', or 'intersection'. Only "
            "used if '--test' is set and final coverage MT and/or platform assignment HT are not already written."
        ),
        type=str,
        choices=["ukb", "broad", "intersection"],
        default="intersection",
    )
    sex_ploidy_args.add_argument(
        "--calling-interval-padding",
        help=(
            "Number of base pair padding to use on the calling intervals. One of 0 or 50 bp. Only used if '--test' is "
            "set and final coverage MT and/or platform assignment HT are not already written."
        ),
        type=int,
        choices=[0, 50],
        default=50,
    )
    sex_karyotype_args = parser.add_argument_group(
        "Annotate sex karyotype", "Arguments used for annotating sex karyotype."
    )
    sex_karyotype_args.add_argument(
        "--annotate-sex-karyotype",
        help="Run sex karyotype inference.",
        action="store_true",
    )
    sex_karyotype_args.add_argument(
        "--per-platform",
        help="Whether to run the karyotype inference per platform.",
        action="store_true",
    )

    args = parser.parse_args()
    main(args)

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
