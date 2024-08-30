"""Script to impute chromosomal sex karyotype annotation."""

import argparse
import json
import logging
from typing import Dict, Optional, Tuple, Union

import hail as hl
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.pipeline import annotate_sex, infer_sex_karyotype
from gnomad.sample_qc.sex import get_sex_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    get_logging_path,
    ukb_f_stat,
)
from gnomad_qc.v4.resources.sample_qc import (
    f_stat_sites,
    get_ploidy_cutoff_json_path,
    hard_filtered_samples_no_sex,
    interval_coverage,
    platform,
    ploidy,
    sex,
    sex_chr_coverage,
    sex_imputation_interval_qc,
)
from gnomad_qc.v4.sample_qc.interval_qc import (
    compute_interval_qc,
    get_high_qual_cutoff_dict,
    get_interval_qc_pass,
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

    This Table is designed to be used as a variant filter in sex imputation for f-stat
    computation.

    .. warning::

        By default `approx_af_and_no_callrate` is False and the final Table will be
        filtered to high callrate (> value specified by `min_callrate`) variants. This
        requires a densify of chrX!"

    .. note::

        If `approx_af_and_no_callrate` is True, allele frequency is approximated with
        AC/(n_samples * 2) and no callrate filter is used.

    :param vds: Input VariantDataset.
    :param approx_af_and_no_callrate: Whether to approximate allele frequency with
        AC/(n_samples * 2) and use no callrate cutoff to filter sites.
    :param min_af: Minimum alternate allele frequency cutoff used to filter sites.
    :param min_callrate: Minimum callrate cutoff used to filter sites.
    :return: Table of chromosome X sites to be used for f-stat computation.
    """
    vds = hl.vds.filter_chromosomes(vds, keep=["chrX"])
    vd = vds.variant_data
    vd = vd.filter_rows(
        (hl.len(vd.alleles) == 2) & hl.is_snp(vd.alleles[0], vd.alleles[1])
    )
    vd = vd.transmute_entries(GT=hl.vds.lgt_to_gt(vd.LGT, vd.LA))

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

        If `test` is True and the test platform assignment Table does not exist, the
        function will load the final platform assignment Table instead if it already
        exists.

    :param test: Whether a test platform assignment Table should be loaded.
    :param calling_interval_name: Name of calling intervals to use for interval
        coverage. One of: 'ukb', 'broad', or 'intersection'. Only used if `test` is
        True.
    :param calling_interval_padding: Number of base pair padding to use on the calling
        intervals. One of 0 or 50 bp. Only used if `test` is True.
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
                "Test platform file does not exist for calling interval %s and interval"
                " padding %s, using final platform assignment Table instead. To use a"
                " test platform assignment please run platform_inference.py"
                " --assign-platforms with the --test argument and needed"
                " --calling-interval-name/--calling-interval-padding arguments.",
                calling_interval_name,
                calling_interval_padding,
            )
    elif test:
        raise FileNotFoundError(
            "There is no test platform assignment Table written for calling interval"
            f" {calling_interval_name} and interval padding"
            f" {calling_interval_padding} and a final platform assignment Table does"
            " not exist. Please run platform_inference.py --assign-platforms with the"
            " --test argument and needed"
            " --calling-interval-name/--calling-interval-padding arguments."
        )
    else:
        raise FileNotFoundError(
            f"There is no final platform assignment Table written. Please run:"
            f" platform_inference.py --assign-platforms to compute the platform"
            f" assignment Table."
        )

    return ht


def prepare_sex_imputation_coverage_mt(
    normalization_contig: str = "chr20",
    test: bool = False,
    read_if_exists: bool = False,
) -> hl.MatrixTable:
    """
    Prepare the sex imputation coverage MT.

    Filter the full interval coverage MatrixTable to the specified normalization contig
    and hard filtered samples (before sex hard filter) and union it with the sex
    coverage MatrixTable after excluding intervals that overlap PAR regions.

    :param normalization_contig: Which autosomal chromosome to use for normalizing the
        coverage of chromosomes X and Y. Default is 'chr20'.
    :param test: Whether to use gnomAD v4 test dataset. Default is False.
    :param read_if_exists: Whether to use the sex imputation coverage MT if it already
        exists rather than remaking this intermediate temporary file. Default is False.
    :return: Interval coverage MatrixTable for sex imputation.
    """
    logger.info(
        "Loading the full interval coverage MT and filtering to the desired"
        " normalization contig..."
    )
    coverage_mt = interval_coverage.mt()
    coverage_mt = coverage_mt.filter_rows(
        coverage_mt.interval.start.contig == normalization_contig
    )

    logger.info("Removing hard-filtered samples from the full interval coverage MT...")
    coverage_mt = coverage_mt.filter_cols(
        hl.is_missing(hard_filtered_samples_no_sex.ht()[coverage_mt.col_key])
    )
    if test:
        sex_coverage_mt = hl.read_matrix_table(
            get_checkpoint_path("test_sex_imputation_cov", mt=True)
        )
        logger.info(
            "Filtering to columns in both the coverage MT and the sex coverage MT for"
            " testing...",
        )
        coverage_mt = coverage_mt.semi_join_cols(sex_coverage_mt.cols())
        sex_coverage_mt = sex_coverage_mt.semi_join_cols(coverage_mt.cols())

        # The test sex HT columns are in a different order than the full cov MT because
        # of a different loading order.
        # This will reorder the sex test HT for the row_union, but doesn't need to
        # be done on the full dataset run.
        sex_coverage_mt_tmp = sex_coverage_mt.annotate_cols(
            idx_in_cov_mt=coverage_mt.add_col_index()
            .index_cols(sex_coverage_mt.s)
            .col_idx
        )
        sex_coverage_mt_tmp = sex_coverage_mt_tmp.add_col_index()
        idx_in_sex_coverage_mt = sex_coverage_mt_tmp.aggregate_cols(
            hl.map(
                lambda x: x.col_idx,
                hl.sorted(
                    hl.agg.collect(
                        sex_coverage_mt_tmp.col.select("idx_in_cov_mt", "col_idx")
                    ),
                    key=lambda x: x.idx_in_cov_mt,
                ),
            )
        )
        sex_coverage_mt = sex_coverage_mt.choose_cols(idx_in_sex_coverage_mt)
    else:
        sex_coverage_mt = sex_chr_coverage.mt()

    logger.info(
        "Excluding intervals that overlap PAR regions from the sex coverage MT..."
    )
    sex_coverage_mt = sex_coverage_mt.filter_rows(~sex_coverage_mt.overlap_par)
    coverage_mt = coverage_mt.union_rows(sex_coverage_mt.drop("overlap_par"))
    coverage_mt = coverage_mt.annotate_globals(
        normalization_contig=normalization_contig
    )
    coverage_mt = coverage_mt.checkpoint(
        get_checkpoint_path(
            f"temp_sex_imputation_cov.{normalization_contig}{'.test' if test else ''}",
            mt=True,
        ),
        _read_if_exists=read_if_exists,
        overwrite=not read_if_exists,
    )

    return coverage_mt.drop("gq_thresholds")


def compute_sex_ploidy(
    vds: hl.vds.VariantDataset,
    coverage_mt: hl.MatrixTable,
    interval_qc_ht: Optional[hl.Table] = None,
    high_qual_per_platform: bool = False,
    platform_ht: Optional[hl.Table] = None,
    normalization_contig: str = "chr20",
    variant_depth_only_x_ploidy: bool = False,
    variant_depth_only_y_ploidy: bool = False,
    variant_depth_only_ploidy_filter_lcr: bool = True,
    variant_depth_only_ploidy_filter_segdup: bool = True,
    variant_depth_only_ploidy_snv_only: bool = False,
    compute_x_frac_variants_hom_alt=True,
    freq_ht: Optional[hl.Table] = None,
    min_af: float = 0.001,
    f_stat_cutoff: float = -1.0,
) -> hl.Table:
    """
    Impute sex chromosome ploidy, and optionally chrX heterozygosity and fraction homozygous alternate variants on chrX.

    This function imputes sex chromosome ploidy from a VDS and a sex inference specific
    coverage MT (created by `prepare_sex_imputation_coverage_mt`).

    With no additional parameters passed, chrX and chrY ploidy will be imputed using
    Hail's `hail.vds.impute_sex_chromosome_ploidy` method which computes chromosome
    ploidy using reference block DP per calling interval (using intervals in
    `coverage_mt`). This method breaks up the reference blocks at the calling interval
    boundaries, maintaining all reference block END information for the mean DP per
    interval computation.

    There is also the option to impute ploidy using mean variant depth within the
    specified calling intervals instead of using reference block depths. This can be
    defined differently for chrX and chrY using `variant_depth_only_x_ploidy` and
    `variant_depth_only_y_ploidy`.

    If an `interval_qc_ht` Table is supplied, only high quality intervals will be used
    in sex chromosome ploidy imputation. High quality intervals are defined by
    'interval_qc_ht.pass_interval_qc'.

    If `high_qual_per_platform` is True, `interval_qc_ht` and `platform_ht` must be
    supplied, and 'interval_qc_ht.pass_interval_qc' should be a Struct with one
    BooleanExpression per platform.

    :param vds: Input VDS for use in sex inference.
    :param coverage_mt: Input sex inference specific interval coverage MatrixTable.
    :param interval_qc_ht: Optional interval QC Table to use for filtering to high
        quality intervals before sex ploidy imputation.
    :param high_qual_per_platform: Whether to filter to per platform high quality
        intervals for the sex ploidy imputation. Default is False.
    :param platform_ht: Input platform assignment Table. This is only needed if
        `high_qual_per_platform` is True.
    :param normalization_contig: Which autosomal chromosome to use for normalizing the
        coverage of chromosomes X and Y. Default is 'chr20'.
    :param variant_depth_only_x_ploidy: Whether to use depth of variant data within
        calling intervals instead of reference data for chrX ploidy estimation. Default
        will only use reference data.
    :param variant_depth_only_y_ploidy: Whether to use depth of variant data within
        calling intervals instead of reference data for chrY ploidy estimation. Default
        will only use reference data.
    :param variant_depth_only_ploidy_filter_lcr: Whether to filter out variants in LCR
        regions for variants only ploidy estimation and fraction of homozygous
        alternate variants on chromosome X. Default is True.
    :param variant_depth_only_ploidy_filter_segdup: Whether to filter out variants in
        segdup regions for variants only ploidy estimation and fraction of homozygous
        alternate variants on chromosome X. Default is True.
    :param variant_depth_only_ploidy_snv_only: Whether to filter to only single
        nucleotide variants for variants only ploidy estimation and fraction of
        homozygous alternate variants on chromosome X. Default is False.
    :param compute_x_frac_variants_hom_alt: Whether to return an annotation for the
        fraction of homozygous alternate variants on chromosome X. Default is True.
    :param freq_ht: Optional Table to use for f-stat allele frequency cutoff. The input
        VDS is filtered to sites in this Table prior to running Hail's `impute_sex`
        module, and alternate allele frequency is used from this Table with a `min_af`
        cutoff.
    :param min_af: Minimum alternate allele frequency to be used in f-stat calculations.
        Default is 0.001.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX
        samples are below cutoff and XY are above cutoff. Default is -1.0.
    :return: Table with imputed ploidies.
    """
    if high_qual_per_platform and platform_ht is None:
        raise ValueError(
            "'platform_ht' must be defined if 'high_qual_per_platform' is True!"
        )

    if high_qual_per_platform and interval_qc_ht is None:
        raise ValueError(
            "'interval_qc_ht' must be defined if 'high_qual_per_platform' is True!"
        )

    add_globals = {}

    def _annotate_sex(
        vds: hl.vds.VariantDataset,
        coverage_mt: hl.MatrixTable,
        calling_intervals_ht: hl.Table,
    ) -> hl.Table:
        """
        Perform `annotate_sex` using unchanged parameters with changes to the VDS and calling intervals.

        :param vds: Input VDS to use for sex ploidy annotation.
        :param coverage_mt: Input precomputed coverage MT to use for sex ploidy
            annotation.
        :param calling_intervals_ht: Table including only intervals wanted for sex
            annotation.
        :return: Table containing sex ploidy estimates for samples in the input VDS.
        """
        ploidy_ht = annotate_sex(
            vds,
            included_intervals=calling_intervals_ht,
            normalization_contig=normalization_contig,
            sites_ht=(
                freq_ht.filter(hl.is_defined(calling_intervals_ht[freq_ht.locus]))
                if freq_ht is not None
                else None
            ),
            aaf_expr="AF",
            gt_expr="LGT",
            f_stat_cutoff=f_stat_cutoff,
            aaf_threshold=min_af,
            variants_only_x_ploidy=variant_depth_only_x_ploidy,
            variants_only_y_ploidy=variant_depth_only_y_ploidy,
            variants_filter_lcr=variant_depth_only_ploidy_filter_lcr,
            variants_filter_segdup=variant_depth_only_ploidy_filter_segdup,
            variants_snv_only=variant_depth_only_ploidy_snv_only,
            compute_x_frac_variants_hom_alt=compute_x_frac_variants_hom_alt,
            coverage_mt=coverage_mt,
            infer_karyotype=False,
        )
        return ploidy_ht

    if interval_qc_ht is not None:
        high_qual_cutoffs = (
            interval_qc_ht.index_globals().high_qual_interval_parameters.collect()[0]
        )
        logger.info(
            "Running sex ploidy imputation using only high quality intervals: %s...",
            high_qual_cutoffs,
        )
        add_globals["high_qual_interval_parameters"] = high_qual_cutoffs

        if high_qual_per_platform:
            if not isinstance(interval_qc_ht.pass_interval_qc, hl.StructExpression):
                ValueError(
                    "'interval_qc_ht.pass_interval_qc' is not a StructExpression"
                    " containing interval pass per platform!"
                )
            logger.info(
                "Running sex ploidy imputation per platform using per platform high"
                " quality intervals..."
            )
            platforms = platform_ht.aggregate(
                hl.agg.collect_as_set(platform_ht.qc_platform)
            )
            per_platform_ploidy_hts = []
            for platform in platforms:
                logger.info(
                    "Performing ploidy imputation using high quality intervals for"
                    " platform %s...",
                    platform,
                )
                filtered_platform_ht = platform_ht.filter(
                    platform_ht.qc_platform == platform
                )
                ploidy_ht = _annotate_sex(
                    hl.vds.filter_samples(vds, filtered_platform_ht),
                    coverage_mt.semi_join_cols(filtered_platform_ht),
                    interval_qc_ht.filter(interval_qc_ht.pass_interval_qc[platform]),
                )
                per_platform_ploidy_hts.append(ploidy_ht)

            ploidy_ht = per_platform_ploidy_hts[0].union(*per_platform_ploidy_hts[1:])

        else:
            logger.info(
                "Running sex ploidy imputation using only high quality intervals in"
                " 'interval_qc_ht'..."
            )
            ploidy_ht = _annotate_sex(
                vds, coverage_mt, interval_qc_ht.filter(interval_qc_ht.pass_interval_qc)
            )
    else:
        logger.info("Running sex ploidy imputation...")
        ploidy_ht = _annotate_sex(vds, coverage_mt, coverage_mt.rows())

    ploidy_ht = ploidy_ht.annotate_globals(
        f_stat_min_af=min_af,
        f_stat_cutoff=f_stat_cutoff,
        **add_globals,
    )

    return ploidy_ht


def annotate_sex_karyotype_from_ploidy_cutoffs(
    ploidy_ht: hl.Table,
    sex_karyotype_ploidy_cutoffs: Union[
        Dict[str, Dict[str, Dict[str, float]]],
        Dict[str, Dict[str, float]],
    ],
    per_platform: bool = False,
    apply_x_frac_hom_alt_cutoffs: bool = False,
) -> hl.Table:
    """
    Determine sex karyotype annotation based on chromosome X and chromosome Y ploidy estimates and ploidy cutoffs.

    `ploidy_ht` must include the following annotations:

        - chrX_ploidy: chromosome X ploidy estimate
        - chrY_ploidy: chromosome X ploidy estimate

    The expected format of `sex_karyotype_ploidy_cutoffs` is one of:

        - If `per_platform` is False:

        .. code-block::

            {
                "x_ploidy_cutoffs": {
                    "upper_cutoff_X": 2.6,
                    "lower_cutoff_XX": 1.9,
                    "upper_cutoff_XX": 6.8,
                    "lower_cutoff_XXX": 6.6
                },
                "y_ploidy_cutoffs": {
                    "lower_cutoff_Y": 0.2,
                    "upper_cutoff_Y": 1.3,
                    "lower_cutoff_YY": 1.4
                }
            }

        - If `per_platform` is True:

        .. code-block::

            {
                "x_ploidy_cutoffs": {
                    "platform_1": {
                        "upper_cutoff_X": 2.6,
                        "lower_cutoff_XX": 1.9,
                        "upper_cutoff_XX": 6.2,
                        "lower_cutoff_XXX": 6.6
                    },
                    "platform_0": {
                        "upper_cutoff_X": 1.6,
                        "lower_cutoff_XX": 1.5,
                        "upper_cutoff_XX": 3.3,
                        "lower_cutoff_XXX": 3.5
                    },
                    ...
                },
                "y_ploidy_cutoffs": {
                    "platform_1": {
                        "lower_cutoff_Y": 0.2,
                        "upper_cutoff_Y": 1.3,
                        "lower_cutoff_YY": 1.4
                    },
                    "platform_0": {
                        "lower_cutoff_Y": 0.1,
                        "upper_cutoff_Y": 1.2,
                        "lower_cutoff_YY": 1.1
                    },
                    ...
                }
            }

    Returns a Table with the following annotations:

        - X_karyotype: Sample assigned X karyotype.
        - Y_karyotype: Sample assigned Y karyotype.
        - sex_karyotype: Combination of X_karyotype and Y_karyotype.

    :param ploidy_ht: Table with chromosome X and chromosome Y ploidies.
    :param sex_karyotype_ploidy_cutoffs: Dictionary of sex karyotype ploidy cutoffs.
    :param per_platform: Whether the `sex_karyotype_ploidy_cutoffs` should be applied
        per platform.
    :param apply_x_frac_hom_alt_cutoffs: Whether to apply cutoffs for the fraction
        homozygous alternate genotypes (hom-alt/(hom-alt + het)) on chromosome X.
    :return: Sex karyotype Table.
    """
    x_ploidy_cutoffs = sex_karyotype_ploidy_cutoffs["x_ploidy_cutoffs"]
    y_ploidy_cutoffs = sex_karyotype_ploidy_cutoffs["y_ploidy_cutoffs"]

    all_cutoff_values = list(x_ploidy_cutoffs.values()) + list(
        y_ploidy_cutoffs.values()
    )

    if apply_x_frac_hom_alt_cutoffs:
        x_frac_hom_alt_cutoffs = sex_karyotype_ploidy_cutoffs["x_frac_hom_alt_cutoffs"]
        all_cutoff_values += list(x_frac_hom_alt_cutoffs.values())
    else:
        x_frac_hom_alt_cutoffs = None

    ploidy_cutoffs_isdict = all(map(lambda x: isinstance(x, dict), all_cutoff_values))
    ploidy_cutoffs_isstr = all(map(lambda x: isinstance(x, str), all_cutoff_values))

    def _format_ploidy_cutoffs(
        x_ploidy_cutoffs: Dict[str, float],
        y_ploidy_cutoffs: Dict[str, float],
        x_frac_hom_alt_cutoffs: Optional[Dict[str, float]] = None,
    ) -> Tuple[
        Tuple[float, Tuple[float, float], float],
        Tuple[Tuple[float, float], float],
        Union[Tuple[Tuple[float, float], float], None],
    ]:
        """
        Reformat ploidy cutoffs for input to `get_sex_expr`.

        :param x_ploidy_cutoffs: Dictionary of X ploidy cutoffs for karyotype assignment.
        :param y_ploidy_cutoffs: Dictionary of Y ploidy cutoffs for karyotype assignment.
        :return: Tuple of ploidy cutoff tuples: ((x_ploidy_cutoffs), (y_ploidy_cutoffs))
        """
        x_ploidy_cutoffs = (
            x_ploidy_cutoffs["upper_cutoff_X"],
            (x_ploidy_cutoffs["lower_cutoff_XX"], x_ploidy_cutoffs["upper_cutoff_XX"]),
            x_ploidy_cutoffs["lower_cutoff_XXX"],
        )
        y_ploidy_cutoffs = (
            (y_ploidy_cutoffs["lower_cutoff_Y"], y_ploidy_cutoffs["upper_cutoff_Y"]),
            y_ploidy_cutoffs["lower_cutoff_YY"],
        )

        if x_frac_hom_alt_cutoffs:
            x_frac_hom_alt_cutoffs = (
                (
                    x_frac_hom_alt_cutoffs["lower_cutoff_more_than_one_X"],
                    x_frac_hom_alt_cutoffs["upper_cutoff_more_than_one_X"],
                ),
                x_frac_hom_alt_cutoffs["lower_cutoff_single_X"],
            )

        return x_ploidy_cutoffs, y_ploidy_cutoffs, x_frac_hom_alt_cutoffs

    logger.info("Annotating sex karyotype based on input ploidy cutoffs")
    if per_platform:
        if not ploidy_cutoffs_isdict:
            raise ValueError(
                "If running per platform X ploidy and Y ploidy cutoff dictionary values"
                " must be dictionaries (one per platform)"
            )
        platforms = ploidy_ht.aggregate(hl.agg.collect_as_set(ploidy_ht.platform))
        per_platform_karyotype_hts = []

        for platform in platforms:
            platform_ploidy_ht = ploidy_ht.filter(ploidy_ht.platform == platform)
            (
                x_ploidy_platform_cutoffs,
                y_ploidy_platform_cutoffs,
                x_frac_hom_alt_platform_cutoffs,
            ) = _format_ploidy_cutoffs(
                x_ploidy_cutoffs[platform],
                y_ploidy_cutoffs[platform],
                x_frac_hom_alt_cutoffs[platform] if x_frac_hom_alt_cutoffs else None,
            )
            karyotype_ht = platform_ploidy_ht.select(
                **get_sex_expr(
                    platform_ploidy_ht.chrX_ploidy,
                    platform_ploidy_ht.chrY_ploidy,
                    x_ploidy_platform_cutoffs,
                    y_ploidy_platform_cutoffs,
                    chr_x_frac_hom_alt_expr=(
                        platform_ploidy_ht.chrx_frac_hom_alt_adj
                        if x_frac_hom_alt_cutoffs
                        else None
                    ),
                    chr_x_frac_hom_alt_cutoffs=x_frac_hom_alt_platform_cutoffs,
                )
            )
            per_platform_karyotype_hts.append(karyotype_ht)

        karyotype_ht = per_platform_karyotype_hts[0].union(
            *per_platform_karyotype_hts[1:]
        )
    else:
        if not ploidy_cutoffs_isstr:
            raise ValueError(
                "X ploidy and Y ploidy cutoff dictionary values must be strings when"
                " not running per platform!"
            )
        (
            x_ploidy_cutoffs,
            y_ploidy_cutoffs,
            x_frac_hom_alt_cutoffs,
        ) = _format_ploidy_cutoffs(
            x_ploidy_cutoffs, y_ploidy_cutoffs, x_frac_hom_alt_cutoffs
        )
        karyotype_ht = ploidy_ht.select(
            **get_sex_expr(
                ploidy_ht.chrX_ploidy,
                ploidy_ht.chrY_ploidy,
                x_ploidy_cutoffs,
                y_ploidy_cutoffs,
                chr_x_frac_hom_alt_expr=(
                    None
                    if x_frac_hom_alt_cutoffs is None
                    else ploidy_ht.chrx_frac_hom_alt_adj
                ),
                chr_x_frac_hom_alt_cutoffs=x_frac_hom_alt_cutoffs,
            )
        )

    def _to_struct(x):
        if isinstance(x, dict):
            return hl.struct(**{v: _to_struct(x[v]) for v in x})
        else:
            return x

    karyotype_ht = karyotype_ht.annotate_globals(
        **_to_struct(sex_karyotype_ploidy_cutoffs)
    )

    return karyotype_ht


def infer_sex_karyotype_from_ploidy(
    ploidy_ht: hl.Table,
    per_platform: bool = False,
    f_stat_cutoff: float = -1.0,
    use_gmm_for_ploidy_cutoffs: bool = False,
    apply_x_frac_hom_alt_cutoffs: bool = False,
) -> hl.Table:
    """
    Create a Table with X_karyotype, Y_karyotype, and sex_karyotype.

    :param ploidy_ht: Table with chromosome X and chromosome Y ploidies, and f-stat if
        not `use_gmm_for_ploidy_cutoffs`.
    :param per_platform: Whether the sex karyotype ploidy cutoff inference should be
        applied per platform.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX
        samples are below cutoff and XY are above cutoff.
    :param use_gmm_for_ploidy_cutoffs: Use Gaussian mixture model to split samples into
        'XX' and 'XY' instead of f-stat.
    :param apply_x_frac_hom_alt_cutoffs: Whether to apply cutoffs for the fraction
        homozygous alternate genotypes (hom-alt/(hom-alt + het)) on chromosome X.
    :return: Table of imputed sex karyotypes.
    """
    logger.info("Running sex karyotype inference")
    if per_platform:
        platforms = ploidy_ht.aggregate(hl.agg.collect_as_set(ploidy_ht.platform))
        per_platform_karyotype_hts = []
        x_ploidy_cutoffs = {}
        y_ploidy_cutoffs = {}
        x_frac_hom_alt_cutoffs = {}

        for platform in platforms:
            logger.info(
                "Performing sex karyotype inference for platform %s...",
                platform,
            )
            ploidy_platform_ht = ploidy_ht.filter(ploidy_ht.platform == platform)

            if apply_x_frac_hom_alt_cutoffs:
                chr_x_frac_hom_alt_expr = ploidy_platform_ht.chrx_frac_hom_alt_adj
            else:
                chr_x_frac_hom_alt_expr = None

            karyotype_ht = infer_sex_karyotype(
                ploidy_platform_ht,
                f_stat_cutoff,
                use_gmm_for_ploidy_cutoffs,
                chr_x_frac_hom_alt_expr=chr_x_frac_hom_alt_expr,
            )

            karyotype_ht = karyotype_ht.checkpoint(
                get_checkpoint_path(f"karyotype_platform_{platform}"), overwrite=True
            )
            per_platform_karyotype_hts.append(karyotype_ht)
            x_ploidy_cutoffs[platform] = karyotype_ht.index_globals().x_ploidy_cutoffs
            y_ploidy_cutoffs[platform] = karyotype_ht.index_globals().y_ploidy_cutoffs
            if apply_x_frac_hom_alt_cutoffs:
                x_frac_hom_alt_cutoffs[platform] = (
                    karyotype_ht.index_globals().x_frac_hom_alt_cutoffs
                )

        karyotype_ht = per_platform_karyotype_hts[0].union(
            *per_platform_karyotype_hts[1:]
        )
        karyotype_ht = karyotype_ht.annotate_globals(
            x_ploidy_cutoffs=hl.struct(**x_ploidy_cutoffs),
            y_ploidy_cutoffs=hl.struct(**y_ploidy_cutoffs),
        )
        if apply_x_frac_hom_alt_cutoffs:
            karyotype_ht = karyotype_ht.annotate_globals(
                x_frac_hom_alt_cutoffs=hl.struct(**x_frac_hom_alt_cutoffs),
            )
    else:
        if apply_x_frac_hom_alt_cutoffs:
            chr_x_frac_hom_alt_expr = ploidy_ht.chrx_frac_hom_alt_adj
        else:
            chr_x_frac_hom_alt_expr = None

        karyotype_ht = infer_sex_karyotype(
            ploidy_ht,
            f_stat_cutoff,
            use_gmm_for_ploidy_cutoffs,
            chr_x_frac_hom_alt_expr=chr_x_frac_hom_alt_expr,
        )

    return karyotype_ht


def reformat_ploidy_cutoffs_for_json(
    ht: hl.Table,
    per_platform: bool = False,
    include_x_frac_hom_alt_cutoffs: bool = True,
) -> dict:
    """
    Format x_ploidy_cutoffs and y_ploidy_cutoffs global annotations for JSON export.

    :param ht: Table including globals for x_ploidy_cutoffs and y_ploidy_cutoffs.
    :param per_platform: Whether the ploidy global cutoffs are per platform.
    :param include_x_frac_hom_alt_cutoffs: Whether to include cutoffs for the fraction
        homozygous alternate genotypes (hom-alt/(hom-alt + het)) on chromosome X.
    :return: Dictionary of X and Y ploidy cutoffs for JSON export.
    """
    x_ploidy_cutoffs = dict(ht.index_globals().x_ploidy_cutoffs.collect()[0])
    y_ploidy_cutoffs = dict(ht.index_globals().y_ploidy_cutoffs.collect()[0])
    if include_x_frac_hom_alt_cutoffs:
        x_frac_hom_alt_cutoffs = dict(
            ht.index_globals().x_frac_hom_alt_cutoffs.collect()[0]
        )

    if per_platform:
        x_ploidy_cutoffs = {k: dict(v) for k, v in x_ploidy_cutoffs.items()}
        y_ploidy_cutoffs = {k: dict(v) for k, v in y_ploidy_cutoffs.items()}
        if include_x_frac_hom_alt_cutoffs:
            x_frac_hom_alt_cutoffs = {
                k: dict(v) for k, v in x_frac_hom_alt_cutoffs.items()
            }

    cutoffs = {
        "x_ploidy_cutoffs": x_ploidy_cutoffs,
        "y_ploidy_cutoffs": y_ploidy_cutoffs,
    }
    if include_x_frac_hom_alt_cutoffs:
        cutoffs.update({"x_frac_hom_alt_cutoffs": x_frac_hom_alt_cutoffs})

    return cutoffs


def main(args):
    """Impute chromosomal sex karyotype annotation."""
    hl.init(
        log="/sex_inference.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default.
    hl._set_flags(use_new_shuffle="1")

    test = args.test
    normalization_contig = args.normalization_contig
    high_qual_intervals = args.high_qual_intervals
    high_qual_per_platform = args.high_qual_per_platform
    high_qual_all_platforms = args.high_qual_all_platforms
    per_platform = args.per_platform
    overwrite = args.overwrite
    read_sex_cov_if_exists = args.read_sex_imputation_coverage_mt_if_exists
    apply_x_frac_hom_alt_cutoffs = args.apply_x_frac_hom_alt_cutoffs

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
                get_checkpoint_path("test_f_stat_sites") if test else f_stat_sites.path,
                overwrite=overwrite,
            )

        if args.sex_imputation_interval_qc:
            sex_coverage_mt = prepare_sex_imputation_coverage_mt(
                normalization_contig,
                test,
                read_sex_cov_if_exists,
            )
            platform_ht = load_platform_ht(
                test,
                sex_coverage_mt.calling_interval_name.collect()[0],
                sex_coverage_mt.calling_interval_padding.collect()[0],
            )
            ht = compute_interval_qc(
                sex_coverage_mt,
                platform_ht=platform_ht,
                mean_dp_thresholds=args.mean_dp_thresholds,
            )
            ht.naive_coalesce(args.interval_qc_n_partitions).write(
                (
                    get_checkpoint_path("test_sex_chr_interval_qc")
                    if test
                    else sex_imputation_interval_qc.path
                ),
                overwrite=overwrite,
            )

        if args.impute_sex_ploidy:
            vds = get_gnomad_v4_vds(
                remove_hard_filtered_samples=False,
                remove_hard_filtered_samples_no_sex=True,
                test=test,
            )
            if args.f_stat_ukb_var:
                # The UK Biobank f-stat table contains only variants that were high
                # callrate (0.99) and common (AF >0.001) within the UK Biobank 200K
                # Regeneron exome dataset and it includes the UK Biobank 200K allele
                # frequency information that can be used in the hl.impute_sex f-stat
                # computation allele frequency cutoff (args.min-af).
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

            # Added because without this impute_sex_chromosome_ploidy will still run
            # even with overwrite=False.
            if file_exists(ploidy_ht_path) and not overwrite:
                raise DataException(
                    f"{ploidy_ht_path} already exists and the --overwrite option was"
                    " not used!"
                )

            coverage_mt = prepare_sex_imputation_coverage_mt(
                normalization_contig,
                test,
                read_sex_cov_if_exists,
            )
            platform_ht = load_platform_ht(
                test,
                coverage_mt.calling_interval_name.collect()[0],
                coverage_mt.calling_interval_padding.collect()[0],
            )

            high_qual_cutoffs = None
            if args.high_qual_by_mean_fraction_over_dp_0:
                # The same cutoffs are used for x_non_par, and y_non_par within their
                # respective dictionaries and the same annotations are used for
                # autosome_par, x_non_par, and y_non_par within their respective
                # dictionaries.
                high_qual_cutoffs = get_high_qual_cutoff_dict(
                    args.norm_mean_fraction_over_dp_0,
                    *[args.sex_mean_fraction_over_dp_0] * 2,
                    *["mean_fraction_over_dp_0"] * 3,
                )
            elif args.high_qual_by_fraction_samples_over_cov:
                high_qual_cutoffs = get_high_qual_cutoff_dict(
                    args.fraction_samples_norm,
                    args.fraction_samples_x,
                    args.fraction_samples_y,
                    f"fraction_over_{args.norm_cov}x",
                    f"fraction_over_{args.x_cov}x",
                    f"fraction_over_{args.y_cov}x",
                )

            interval_qc_ht = None
            if high_qual_intervals or high_qual_per_platform or high_qual_all_platforms:
                interval_qc_ht = (
                    hl.read_table(get_checkpoint_path("test_sex_chr_interval_qc"))
                    if test
                    else sex_imputation_interval_qc.ht()
                )
                if (
                    interval_qc_ht.normalization_contig.collect()[0]
                    != normalization_contig
                ):
                    raise ValueError(
                        f"Normalization contig {normalization_contig} is not in the sex"
                        " imputation interval QC HT, rerun"
                        " '--sex-imputation-interval-qc' to use a different"
                        " normalization contig than"
                        f" {interval_qc_ht.normalization_contig.collect()[0]}.",
                    )

                interval_qc_ht = get_interval_qc_pass(
                    interval_qc_ht,
                    high_qual_cutoffs,
                    per_platform=high_qual_per_platform,
                    all_platforms=high_qual_all_platforms,
                    min_platform_size=args.min_platform_size,
                )

            ploidy_ht = compute_sex_ploidy(
                vds,
                coverage_mt,
                interval_qc_ht=interval_qc_ht,
                high_qual_per_platform=high_qual_per_platform,
                platform_ht=platform_ht,
                normalization_contig=normalization_contig,
                variant_depth_only_x_ploidy=args.variant_depth_only_x_ploidy,
                variant_depth_only_y_ploidy=args.variant_depth_only_y_ploidy,
                variant_depth_only_ploidy_filter_lcr=not args.omit_variant_depth_ploidy_lcr_filter,
                variant_depth_only_ploidy_filter_segdup=not args.omit_variant_depth_ploidy_segdup_filter,
                variant_depth_only_ploidy_snv_only=args.variant_depth_ploidy_snv_only,
                compute_x_frac_variants_hom_alt=not args.omit_compute_x_frac_variants_hom_alt,
                freq_ht=freq_ht,
                min_af=args.min_af,
                f_stat_cutoff=args.f_stat_cutoff,
            )
            ploidy_ht = ploidy_ht.annotate(
                platform=platform_ht[ploidy_ht.key].qc_platform
            )
            ploidy_ht = ploidy_ht.annotate_globals(f_stat_ukb_var=args.f_stat_ukb_var)
            logger.info("Writing ploidy Table...")
            ploidy_ht.write(ploidy_ht_path, overwrite=overwrite)

        if args.annotate_sex_karyotype:
            # TODO: Add drop of `is_female` in future versions.
            ploidy_ht = (
                hl.read_table(get_checkpoint_path(f"ploidy_imputation"))
                if test
                else ploidy.ht()
            )

            if args.sex_karyotype_cutoffs:
                with hl.hadoop_open(args.sex_karyotype_cutoffs, "r") as d:
                    ploidy_cutoffs = json.load(d)
                karyotype_ht = annotate_sex_karyotype_from_ploidy_cutoffs(
                    ploidy_ht,
                    ploidy_cutoffs,
                    per_platform=per_platform,
                    apply_x_frac_hom_alt_cutoffs=apply_x_frac_hom_alt_cutoffs,
                )
            else:
                karyotype_ht = infer_sex_karyotype_from_ploidy(
                    ploidy_ht,
                    per_platform=per_platform,
                    f_stat_cutoff=args.f_stat_cutoff,
                    use_gmm_for_ploidy_cutoffs=args.use_gmm_for_ploidy_cutoffs,
                    apply_x_frac_hom_alt_cutoffs=apply_x_frac_hom_alt_cutoffs,
                )
            sex_ht = ploidy_ht.annotate(**karyotype_ht[ploidy_ht.key])
            sex_ht = sex_ht.annotate_globals(**karyotype_ht.index_globals())

            logger.info("Writing sex HT with karyotype annotation...")
            sex_ht.write(
                get_checkpoint_path("sex") if test else sex.path,
                overwrite=overwrite,
            )

            ploidy_cutoffs = reformat_ploidy_cutoffs_for_json(
                sex_ht,
                per_platform=per_platform,
                include_x_frac_hom_alt_cutoffs=apply_x_frac_hom_alt_cutoffs,
            )
            cutoff_json_path = get_ploidy_cutoff_json_path(test=test)
            logger.info("Writing ploidy cutoffs dictionary to %s.", cutoff_json_path)
            with hl.hadoop_open(cutoff_json_path, "w") as d:
                d.write(json.dumps(ploidy_cutoffs))
    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("sex_inference"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
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
            "Create Table of common (> value specified by '--min-af'), bi-allelic SNPs"
            " on chromosome X for f-stat calculations. Additionally filter to high"
            " callrate (> value specified by '--min-callrate') variants if"
            " '--approx-af-and-no-callrate' is not used. NOTE: This requires a densify"
            " of chrX!"
        ),
        action="store_true",
    )
    fstat_args.add_argument(
        "--min-callrate", help="Minimum variant callrate.", default=0.99, type=float
    )
    fstat_args.add_argument(
        "--approx-af-and-no-callrate",
        help=(
            "Whether to approximate allele frequency with AC/(n_samples * 2) and use no"
            " callrate cutoff for determination of f-stat sites."
        ),
        action="store_true",
    )
    fstat_args.add_argument(
        "--fstat-n-partitions",
        help="Number of desired partitions for the f-stat sites output Table.",
        default=1000,
        type=int,
    )

    sex_imputation_interval_qc_args = parser.add_argument_group(
        "Sex chromosome interval QC",
        "Arguments used for making an interval QC HT from the sex imputation interval"
        " coverage MT.",
    )
    sex_imputation_interval_qc_args.add_argument(
        "--sex-imputation-interval-qc",
        help=(
            "Create a Table of the fraction of samples per interval and per platform"
            " with mean DP over thresholds specified by '--mean-dp-thresholds'."
        ),
        action="store_true",
    )
    sex_cov_mt_read_if_exists_action = sex_imputation_interval_qc_args.add_argument(
        "--read-sex-imputation-coverage-mt-if-exists",
        help=(
            "Whether to use the sex imputation coverage MT if it already exists rather"
            " than remaking this intermediate temporary file."
        ),
        action="store_true",
    )
    norm_contig_action = sex_imputation_interval_qc_args.add_argument(
        "--normalization-contig",
        help=(
            "Which autosomal chromosome to use for normalizing the coverage of"
            " chromosomes X and Y."
        ),
        type=str,
        default="chr20",
    )
    sex_imputation_interval_qc_args.add_argument(
        "--mean-dp-thresholds",
        help=(
            "List of mean DP cutoffs to determining the fraction of samples with mean"
            " coverage >= the cutoff for each interval."
        ),
        type=int,
        nargs="+",
        default=[5, 10, 15, 20, 25],
    )
    sex_imputation_interval_qc_args.add_argument(
        "--interval-qc-n-partitions",
        help=(
            "Number of desired partitions for the sex imputation interval QC output"
            " Table."
        ),
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
    # Indicate that the --normalization-contig and
    # --read-sex-imputation-coverage-mt-if-exists options are available for the
    # "sex_ploidy_args" argument group as well.
    sex_ploidy_args._group_actions.append(norm_contig_action)
    sex_ploidy_args._group_actions.append(sex_cov_mt_read_if_exists_action)
    sex_ploidy_args.add_argument(
        "--f-stat-ukb-var",
        help=(
            "Whether to use UK Biobank high callrate (0.99) and common variants (UKB"
            " allele frequency > value specified by '--min-af') for f-stat computation"
            " instead of the sites determined by '--determine-fstat-sites'."
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
            "Cutoff for f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX"
            " samples are below cutoff and XY are above cutoff."
        ),
        type=float,
        default=-1.0,
    )
    sex_ploidy_high_qual_method_parser = sex_ploidy_args.add_mutually_exclusive_group(
        required=False
    )
    sex_ploidy_high_qual_method_parser.add_argument(
        "--high-qual-intervals",
        help=(
            "Whether to filter to high quality intervals for the sex ploidy imputation."
            " Can't be used at the same time as '--high-qual-per-platform' or"
            " '--high-qual-all-platforms'."
        ),
        action="store_true",
    )
    sex_ploidy_high_qual_method_parser.add_argument(
        "--high-qual-per-platform",
        help=(
            "Whether to filter to per platform high quality intervals for the sex"
            " ploidy imputation. Can't be used at the same time as"
            " '--high-qual-intervals' or '--high-qual-all-platforms'."
        ),
        action="store_true",
    )
    sex_ploidy_high_qual_method_parser.add_argument(
        "--high-qual-all-platforms",
        help=(
            "Whether to filter to high quality intervals for the sex ploidy imputation."
            " Use only intervals that are considered high quality across all platforms."
            " Can't be used at the same time as '--high-qual-intervals' or"
            " '--high-qual-per-platform'"
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--min-platform-size",
        help=(
            "Required size of a platform to be considered in"
            " '--high-qual-all-platforms'. Only platforms that have # of samples >"
            " 'min_platform_size' are used to determine intervals that have a high"
            " quality across all platforms."
        ),
        type=int,
        default=100,
    )
    sex_ploidy_args.add_argument(
        "--variant-depth-only-x-ploidy",
        help=(
            "Whether to use depth of variant data for the x ploidy estimation instead"
            " of the default behavior that will use reference blocks."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--variant-depth-only-y-ploidy",
        help=(
            "Whether to use depth of variant data for the y ploidy estimation instead"
            " of the default behavior that will use reference blocks."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--omit-variant-depth-ploidy-lcr-filter",
        help=(
            "Whether to omit filtering out variants in LCR regions for the variants"
            " only ploidy estimation and fraction of homozygous alternate variants on"
            " chromosome X."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--omit-variant-depth-ploidy-segdup-filter",
        help=(
            "Whether to omit filtering out variants in segdup regions for the variants"
            " only ploidy estimation and fraction of homozygous alternate variants on"
            " chromosome X."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--variant-depth-ploidy-snv-only",
        help=(
            "Whether to filter to only single nucleotide variants for variants only"
            " ploidy estimation and fraction of homozygous alternate variants on"
            " chromosome X."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--omit-compute-x-frac-variants-hom-alt",
        help=(
            "Whether to omit the computation of the fraction of homozygous alternate"
            " variants on chromosome X."
        ),
        action="store_true",
    )
    sex_ploidy_high_qual_opt_parser = sex_ploidy_args.add_mutually_exclusive_group(
        required=False
    )
    sex_ploidy_high_qual_opt_parser.add_argument(
        "--high-qual-by-mean-fraction-over-dp-0",
        help=(
            "Whether to use the mean fraction of bases over DP 0 to determine high"
            " quality intervals. Can't be set at the same time as"
            " '--high-qual-by-fraction-samples-over-cov'."
        ),
        action="store_true",
    )
    sex_ploidy_high_qual_opt_parser.add_argument(
        "--high-qual-by-fraction-samples-over-cov",
        help=(
            "Whether to determine high quality intervals using the fraction of samples"
            " with a mean interval quality over a specified quality for chrX (--x-cov),"
            " chrY (--y-cov), and the normalization contig (--norm-cov). Can't be set"
            " at the same time as '--high-qual-by-mean-fraction-over-dp-0'."
        ),
        action="store_true",
    )
    sex_ploidy_args.add_argument(
        "--sex-mean-fraction-over-dp-0",
        help=(
            "Mean fraction of bases over DP 0 used to define high quality intervals on"
            " sex chromosomes."
        ),
        type=float,
        default=0.4,
    )
    sex_ploidy_args.add_argument(
        "--norm-mean-fraction-over-dp-0",
        help=(
            "Mean fraction of bases over DP 0 used to define high quality intervals on"
            " the normalization chromosome."
        ),
        type=float,
        default=0.99,
    )
    sex_ploidy_args.add_argument(
        "--x-cov",
        help=(
            "Mean coverage level used to define high quality intervals on chromosome X."
            " Aggregate mean for this coverage level must be in the sex chromosome"
            " interval QC HT (must be value in '--mean-dp-thresholds' list used to"
            " create the QC HT)!"
        ),
        type=int,
        default=10,
    )
    sex_ploidy_args.add_argument(
        "--y-cov",
        help=(
            "Mean coverage level used to define high quality intervals on chromosome Y."
            " Aggregate mean for this coverage level must be in the sex chromosome"
            " interval QC HT (must be value in '--mean-dp-thresholds' list used to"
            " create the QC HT)!"
        ),
        type=int,
        default=5,
    )
    sex_ploidy_args.add_argument(
        "--norm-cov",
        help=(
            "Mean coverage level used to define high quality intervals on the"
            " normalization autosome. Aggregate mean for this coverage level must be in"
            " the sex chromosome interval QC HT (must be value in"
            " '--mean-dp-thresholds' list used to create the QC HT)!"
        ),
        type=int,
        default=20,
    )
    sex_ploidy_args.add_argument(
        "--fraction-samples-x",
        help=(
            "Fraction samples at specified coverage '--x-cov' to determine high quality"
            " intervals on chromosome X."
        ),
        type=float,
        default=0.80,
    )
    sex_ploidy_args.add_argument(
        "--fraction-samples-y",
        help=(
            "Fraction samples at specified coverage '--y-cov' to determine high quality"
            " intervals on chromosome Y."
        ),
        type=float,
        default=0.35,
    )
    sex_ploidy_args.add_argument(
        "--fraction-samples-norm",
        help=(
            "Fraction samples at specified coverage '--norm-cov' to determine high"
            " quality intervals on the normalization chromosome specified by"
            " '--normalization-contig'."
        ),
        type=float,
        default=0.85,
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
        "--use-gmm-for-ploidy-cutoffs",
        help=(
            "Whether to use Gaussian mixture model to roughly split samples into 'XX'"
            " and 'XY' instead of f-stat."
        ),
        action="store_true",
    )
    sex_karyotype_args.add_argument(
        "--apply-x-frac-hom-alt-cutoffs",
        help=(
            "Whether to apply 'XX' and 'XY' cutoffs for the fraction of homozygous"
            " alternate genotypes on chromosome X and use them to infer sex karyotype."
        ),
        action="store_true",
    )
    sex_karyotype_args.add_argument(
        "--per-platform",
        help="Whether to run the karyotype inference per platform.",
        action="store_true",
    )
    sex_karyotype_args.add_argument(
        "--sex-karyotype-cutoffs",
        help=(
            "Optional path to JSON file containing sex karyotype X and Y ploidy cutoffs"
            " to use for karyotype annotation instead of inferring cutoffs. If"
            " '--apply-x-frac-hom-alt-cutoffs' is used, this file mustalso include"
            " cutoffs for the fraction of homozygous alternate genotypes on"
            " chromosome X."
        ),
        type=str,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if (
        args.high_qual_intervals
        or args.high_qual_per_platform
        or args.high_qual_all_platforms
    ) and not (
        args.high_qual_by_mean_fraction_over_dp_0
        or args.high_qual_by_fraction_samples_over_cov
    ):
        parser.error(
            "One of --high-qual-by-mean-fraction-over-dp-0 or"
            " --high-qual-by-fraction-samples-over-cov is required when a high quality"
            " option (--high-qual-intervals, --high-qual-per-platform, or"
            " --high-qual-all-platforms) is specified."
        )

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
