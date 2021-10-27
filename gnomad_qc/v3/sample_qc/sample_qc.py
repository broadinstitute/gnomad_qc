import argparse
import logging
import math
import pickle
from typing import Any, List, Tuple

import hail as hl

from gnomad.assessment.validity_checks import compare_row_counts
from gnomad.resources.grch38.reference_data import (
    clinvar,
    lcr_intervals,
    purcell_5k_intervals,
    telomeres_and_centromeres,
)
from gnomad.sample_qc.ancestry import assign_population_pcs, run_pca_with_relateds
from gnomad.sample_qc.filtering import (
    compute_qc_metrics_residuals,
    compute_stratified_metrics_filter,
    compute_stratified_sample_qc,
)
from gnomad.sample_qc.pipeline import annotate_sex, get_qc_mt
from gnomad.sample_qc.relatedness import (
    compute_related_samples_to_drop,
    DUPLICATE_OR_TWINS,
    get_relationship_expr,
    PARENT_CHILD,
    SIBLINGS,
    UNRELATED,
)
from gnomad.sample_qc.sex import get_ploidy_cutoffs, get_sex_expr
from gnomad.utils.annotations import bi_allelic_expr, get_adj_expr
from gnomad.utils.filtering import (
    add_filters_expr,
    filter_low_conf_regions,
    filter_to_autosomes,
    filter_to_clinvar_pathogenic,
)
from gnomad.utils.sparse_mt import densify_sites, filter_ref_blocks

from gnomad_qc.v2.resources.sample_qc import get_liftover_v2_qc_mt
from gnomad_qc.v3.resources.annotations import freq, get_info, last_END_position
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta, meta_tsv_path, project_meta
from gnomad_qc.v3.resources.sample_qc import (
    ancestry_pca_eigenvalues,
    ancestry_pca_loadings,
    ancestry_pca_scores,
    sample_clinvar_count,
    get_sample_qc,
    hard_filtered_samples,
    pc_relate_pca_scores,
    pca_related_samples_to_drop,
    pca_samples_rankings,
    picard_metrics,
    pop,
    pop_rf_path,
    pop_tsv_path,
    qc,
    regressed_metrics,
    relatedness,
    release_related_samples_to_drop,
    release_samples_rankings,
    sample_inbreeding,
    sex,
    stratified_metrics,
)


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("sample_qc")
logger.setLevel(logging.INFO)

SUBSETS = [
    "non_topmed",
    "controls_and_biobanks",
    "non_neuro",
    "non_v2",
    "non_cancer",
    "tgp",
    "hgdp",
]
"""
Default list of subsets to include in the sample QC meta HT from the project meta HT
"""


def compute_sample_qc() -> hl.Table:
    """
    Perform sample QC on the raw split matrix table using `compute_stratified_sample_qc`.

    :return: Table containing sample QC metrics
    :rtype: hl.Table
    """
    logger.info("Computing sample QC")
    mt = filter_to_autosomes(
        get_gnomad_v3_mt(
            split=True,
            key_by_locus_and_alleles=True,
            remove_hard_filtered_samples=False,
        )
    )
    mt = mt.select_entries("GT")

    # Filter reference blocks
    mt = filter_ref_blocks(mt)

    # Remove centromeres and telomeres incase they were included
    mt = filter_low_conf_regions(
        mt,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False,
        filter_telomeres_and_centromeres=True,
    )

    sample_qc_ht = compute_stratified_sample_qc(
        mt,
        strata={
            "bi_allelic": bi_allelic_expr(mt),
            "multi_allelic": ~bi_allelic_expr(mt),
        },
        tmp_ht_prefix=get_sample_qc().path[:-3],
    )

    # Remove annotations that cannot be computed from the sparse format
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            x: sample_qc_ht[x].drop(
                "n_called", "n_not_called", "n_filtered", "call_rate"
            )
            for x in sample_qc_ht.row_value
        }
    )

    return sample_qc_ht.repartition(100)


def compute_qc_mt(
    min_af: float = 0.0, min_inbreeding_coeff_threshold: float = -0.025
) -> hl.MatrixTable:
    """
    Generate the QC MT to be used in downstream sample QC.
    
    Filter MT to bi-allelic SNVs at gnomAD v2 QC and p5k sites. Compute and return QC MT.

    :param min_af: Minimum variant allele frequency to retain variant in qc matrix table
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to keep
    :return: MatrixTable filtered to variants for sample QC
    """
    # Load v2 and p5k sites for QC
    v2_qc_sites = get_liftover_v2_qc_mt("joint", ld_pruned=True).rows().key_by("locus")
    qc_sites = v2_qc_sites.union(purcell_5k_intervals.ht(), unify=True)

    qc_sites = qc_sites.filter(hl.is_missing(lcr_intervals.ht()[qc_sites.key]))

    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
    mt = mt.select_entries(
        "END", GT=mt.LGT, adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
    )
    mt = densify_sites(mt, qc_sites, hl.read_table(last_END_position.path))

    # Filter MT to bi-allelic SNVs that are found in the v2 QC and p5k HT
    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2)
        & hl.is_snp(mt.alleles[0], mt.alleles[1])
        & (qc_sites[mt.locus].alleles == mt.alleles)
    )
    mt = mt.checkpoint(
        "gs://gnomad-tmp/gnomad_v3_qc_mt_v2_sites_dense.mt", overwrite=True
    )
    mt = mt.naive_coalesce(5000)
    mt = mt.checkpoint(
        "gs://gnomad-tmp/gnomad_v3_qc_mt_v2_sites_dense_repartitioned.mt",
        overwrite=True,
    )
    info_ht = get_info(split=False).ht()
    info_ht = info_ht.annotate(
        info=info_ht.info.select(
            # No need for AS_annotations since it's bi-allelic sites only
            **{x: info_ht.info[x] for x in info_ht.info if not x.startswith("AS_")}
        )
    )
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
    qc_mt = get_qc_mt(
        mt,
        min_af=min_af,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=None,
        ld_r2=None,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False,
    )
    return qc_mt


def compute_hard_filters(
    cov_threshold: int = 15,
    max_n_snp: float = 3.75e6,
    min_n_snp: float = 2.4e6,
    max_n_singleton: float = 1e5,
    max_r_het_hom_var: float = 3.3,
    max_pct_contamination: float = 5.00,
    max_pct_chimera: float = 5.00,
    min_median_insert_size: int = 250,
    include_sex_filter: bool = True,
) -> hl.Table:
    """
    Apply hard filters to samples and return Table with samples and the reason for filtering. 
    
    This function expects a sex table generated by `impute_sex`.

    :param cov_threshold: Filtering threshold to use for chr20 coverage
    :param max_n_snp: Filtering threshold to use for the max number of SNPs
    :param min_n_snp: Filtering threshold to use for the min number of SNPs
    :param max_n_singleton: Filtering threshold to use for the max number of singletons
    :param max_r_het_hom_var: Filtering threshold to use for the max ratio of heterozygotes to alternate homozygotes
    :param max_pct_contamination: Filtering threshold to use for max percent contamination (this is a percent not a
        proportion, e.g. 5% == 5.00, %5 != 0.05)
    :param max_pct_chimera: Filtering threshold to use for max percent chimera (this is a percent not a proportion,
        e.g. 5% == 5.00, %5 != 0.05)
    :param min_median_insert_size: Filtering threshold to use for min median insert size
    :param include_sex_filter: Should sex inference be used in filtering
    :return: Table of hard filtered samples
    :rtype: hl.Table
    """
    ht = get_gnomad_v3_mt(remove_hard_filtered_samples=False).cols()
    hard_filters = dict()

    # Remove samples failing fingerprinting
    # TODO: Add these into hard filtering metadata when incorporating internal smaples Picard metrics
    hard_filters["failed_fingerprinting"] = hl.array(
        ["09C90823", "10C103592", "S5530"]
    ).contains(ht.s)

    # Remove TCGA tumor samples based on TCGA naming convention: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    hard_filters["TCGA_tumor_sample"] = ht.s.startswith("TCGA") & (
        hl.int(hl.str(ht.s).split("-")[3][:2]) < 10
    )

    # Remove low-coverage samples
    # chrom 20 coverage is computed to infer sex and used here
    cov_ht = sex.ht()
    hard_filters["low_coverage"] = cov_ht[ht.key].chr20_mean_dp < cov_threshold

    # Remove extreme raw bi-allelic sample QC outliers
    # These were determined by visual inspection of the metrics
    bi_allelic_qc_struct = get_sample_qc("bi_allelic").ht()[ht.key]
    hard_filters["bad_qc_metrics"] = (
        (bi_allelic_qc_struct.sample_qc.n_snp > max_n_snp)
        | (bi_allelic_qc_struct.sample_qc.n_snp < min_n_snp)
        | (bi_allelic_qc_struct.sample_qc.n_singleton > max_n_singleton)
        | (bi_allelic_qc_struct.sample_qc.r_het_hom_var > max_r_het_hom_var)
    )

    # Remove samples that fail picard metric thresholds
    picard_ht = picard_metrics.ht()[ht.key]
    hard_filters["contamination"] = (
        picard_ht.bam_metrics.freemix > max_pct_contamination
    )
    hard_filters["chimera"] = picard_ht.bam_metrics.pct_chimeras > max_pct_chimera
    hard_filters["insert_size"] = (
        picard_ht.bam_metrics.median_insert_size < min_median_insert_size
    )

    if include_sex_filter:
        # Remove samples with ambiguous sex assignments
        sex_ht = sex.ht()[ht.key]
        hard_filters["ambiguous_sex"] = sex_ht.sex_karyotype == "ambiguous"
        hard_filters["sex_aneuploidy"] = ~hl.set({"ambiguous", "XX", "XY"}).contains(  # pylint: disable=invalid-unary-operand-type
            sex_ht.sex_karyotype
        )

    ht = ht.annotate(hard_filters=add_filters_expr(filters=hard_filters))

    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    ht = ht.annotate_globals(
        hard_filter_cutoffs=hl.struct(
            min_cov=cov_threshold,
            max_n_snp=max_n_snp,
            min_n_snp=min_n_snp,
            max_n_singleton=max_n_singleton,
            max_r_het_hom_var=max_r_het_hom_var,
            max_pct_contamination=max_pct_contamination,
            max_pct_chimera=max_pct_chimera,
            min_median_insert_size=min_median_insert_size,
        ),
    )

    return ht


def compute_sex(aaf_threshold: float = 0.001, f_stat_cutoff: float = 0.5) -> hl.Table:
    """
    Impute sample sex based on X-chromosome heterozygosity and sex chromosome ploidy.

    Allele frequencies from v3.0 are used to prevent the need to densify the current version's sparse MT.

    :param aaf_threshold: Minimum alternate allele frequency to be used in f-stat calculations.
    :param f_stat_cutoff: f-stat to roughly divide 'XX' from 'XY' samples. Assumes XX samples are below cutoff and XY are above cutoff.
    :return: Table with inferred sex annotation
    :rtype: hl.Table
    """
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, remove_hard_filtered_samples=False,
    )

    # Use AF from v3
    freq_ht = freq.versions["3"].ht()
    freq_ht = freq_ht.select(AF=freq_ht.freq[0].AF)
    freq_ht = freq_ht.filter(freq_ht.AF > aaf_threshold)

    sex_ht = annotate_sex(
        mt,
        excluded_intervals=telomeres_and_centromeres.ht(),
        aaf_threshold=aaf_threshold,
        f_stat_cutoff=f_stat_cutoff,
        sites_ht=freq_ht,
        aaf_expr="AF",
        gt_expr="LGT",
    )

    return sex_ht


def reannotate_sex(
    cov_threshold: int,
    x_ploidy_cutoffs: Tuple[float, Tuple[float, float], float],
    y_ploidy_cutoffs: Tuple[Tuple[float, float], float],
):
    """
    Rerun sex karyotyping annotations without re-computing sex imputation metrics.

    :param cov_threshold: Filtering threshold to use for chr20 coverage in `compute_hard_filters`
    :param x_ploidy_cutoffs: Tuple of X chromosome ploidy cutoffs: (upper cutoff for single X, (lower cutoff for double X, upper cutoff for double X), lower cutoff for triple X)
    :param y_ploidy_cutoffs: Tuple of Y chromosome ploidy cutoffs: ((lower cutoff for single Y, upper cutoff for single Y), lower cutoff for double Y)
    :return: Table with sex karyotyping annotations
    :rtype: hl.Table
    """
    # Copy HT to temp location to overwrite annotation
    sex_ht = sex.ht().checkpoint("gs://gnomad-tmp/sex_ht_checkpoint.ht", overwrite=True)
    hard_filter_ht = compute_hard_filters(cov_threshold, include_sex_filter=False)

    # Copy HT to temp location because it uses sex_ht for chr20 coverage
    hard_filter_ht = hard_filter_ht.checkpoint(
        "gs://gnomad-tmp/hardfilter_checkpoint.ht", overwrite=True
    )
    new_x_ploidy_cutoffs, new_y_ploidy_cutoffs = get_ploidy_cutoffs(
        sex_ht.filter(hl.is_missing(hard_filter_ht[sex_ht.key])), f_stat_cutoff=0.5
    )
    x_ploidy_cutoffs = hl.struct(
        upper_x=x_ploidy_cutoffs[0] if x_ploidy_cutoffs[0] else new_x_ploidy_cutoffs[0],
        lower_xx=x_ploidy_cutoffs[1][0]
        if x_ploidy_cutoffs[1][0]
        else new_x_ploidy_cutoffs[1][0],
        upper_xx=x_ploidy_cutoffs[1][1]
        if x_ploidy_cutoffs[1][1]
        else new_x_ploidy_cutoffs[1][1],
        lower_xxx=x_ploidy_cutoffs[2]
        if x_ploidy_cutoffs[2]
        else new_x_ploidy_cutoffs[2],
    )
    y_ploidy_cutoffs = hl.struct(
        lower_y=y_ploidy_cutoffs[0][0]
        if y_ploidy_cutoffs[0][0]
        else new_y_ploidy_cutoffs[0][0],
        upper_y=y_ploidy_cutoffs[0][1]
        if y_ploidy_cutoffs[0][1]
        else new_y_ploidy_cutoffs[0][1],
        lower_yy=y_ploidy_cutoffs[1]
        if y_ploidy_cutoffs[1]
        else new_y_ploidy_cutoffs[1],
    )
    sex_ht = sex_ht.annotate(
        **get_sex_expr(
            sex_ht.chrX_ploidy,
            sex_ht.chrY_ploidy,
            (
                x_ploidy_cutoffs["upper_x"],
                (x_ploidy_cutoffs["lower_xx"], x_ploidy_cutoffs["upper_xx"]),
                x_ploidy_cutoffs["lower_xxx"],
            ),
            (
                (y_ploidy_cutoffs["lower_y"], y_ploidy_cutoffs["upper_y"]),
                y_ploidy_cutoffs["lower_yy"],
            ),
        )
    )
    sex_ht = sex_ht.annotate_globals(
        x_ploidy_cutoffs=x_ploidy_cutoffs, y_ploidy_cutoffs=y_ploidy_cutoffs,
    )

    return sex_ht


def compute_sample_rankings(use_qc_metrics_filters: bool) -> hl.Table:
    """
    Rank samples by sample filtering and chr20 mean coverage.

    Filtered samples are ranked lowest, followed by non-releasable samples, and then by coverage where higher coverage
    is a higher rank.

    :param use_qc_metrics_filters: Should sample QC metric filters be included or only hard filters
    :return: Table with sample rankings
    :rtype: hl.Table
    """
    project_ht = project_meta.ht()
    project_ht = project_ht.select(
        "releasable",
        "exclude",
        chr20_mean_dp=sex.ht()[project_ht.key].chr20_mean_dp,
        filtered=hl.or_else(
            hl.len(hard_filtered_samples.ht()[project_ht.key].hard_filters) > 0, False
        ),
    )

    if use_qc_metrics_filters:
        project_ht = project_ht.annotate(
            filtered=hl.cond(
                project_ht.filtered,
                True,
                hl.or_else(
                    hl.len(regressed_metrics.ht()[project_ht.key].qc_metrics_filters)
                    > 0,
                    False,
                ),
            )
        )

    project_ht = project_ht.order_by(
        project_ht.filtered,
        hl.desc(project_ht.releasable & ~project_ht.exclude),
        hl.desc(project_ht.chr20_mean_dp),
    ).add_index(name="rank")

    return project_ht.key_by("s").select("filtered", "rank")


def run_pca(
    include_unreleasable_samples: bool, n_pcs: int, related_samples_to_drop: hl.Table
) -> Tuple[List[float], hl.Table, hl.Table]:
    """
    Run population PCA using `run_pca_with_relateds`.

    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA
    :param n_pcs: Number of PCs to compute
    :param related_samples_to_drop: Table of related samples to drop from PCA run
    :return: eigenvalues, scores and loadings
    """
    logger.info("Running population PCA")
    qc_mt = qc.mt()

    samples_to_drop = related_samples_to_drop.select()
    if not include_unreleasable_samples:
        logger.info("Excluding unreleasable samples for PCA.")
        samples_to_drop = samples_to_drop.union(
            qc_mt.filter_cols(
                ~project_meta.ht()[qc_mt.col_key].releasable
                | project_meta.ht()[qc_mt.col_key].exclude
            )
            .cols()
            .select()
        )
    else:
        logger.info("Including unreleasable samples for PCA")

    return run_pca_with_relateds(qc_mt, samples_to_drop, n_pcs=n_pcs)


def assign_pops(
    min_prob: float,
    include_unreleasable_samples: bool,
    max_mislabeled_training_samples: int = 50,  # TODO: Think about this parameter and add it to assign_population_pcs. Maybe should be a fraction? fraction per pop?
    n_pcs: int = 16,
    withhold_prop: float = None,
) -> Tuple[hl.Table, Any]:
    """
    Use a random forest model to assign global population labels based on the results from `assign_pca`.

    The training data used is `project_pop` on the project metadata HT where it is defined, otherwise it uses the
    `v2_pop` if defined and not `oth`.

    The method starts by inferring the pop on all samples, then comparing the training data to the inferred data,
    removing the truth outliers and re-training. This is repeated until the number of truth samples is less than
    or equal to `max_mislabeled_training_samples`.

    :param min_prob: Minimum RF probability for pop assignment
    :param include_unreleasable_samples: Should unreleasable samples be included in the PCA
    :param max_mislabeled_training_samples: Maximum number of training samples that can be mislabelled
    :param n_pcs: Number of PCs to use in the RF
    :param withhold_prop: Proportion of training pop samples to withhold from training will keep all samples if `None`
    :return: Table of pop assignments and the rf model
    :rtype: hl.Table
    """
    logger.info("Assigning global population labels")
    pop_pca_scores_ht = ancestry_pca_scores(include_unreleasable_samples).ht()
    project_meta_ht = project_meta.ht()[pop_pca_scores_ht.key]
    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        training_pop=(
            hl.case()
            .when(
                hl.is_defined(project_meta_ht.project_pop), project_meta_ht.project_pop
            )
            .when(project_meta_ht.v2_pop != "oth", project_meta_ht.v2_pop)
            .or_missing()
        )
    )
    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        training_pop_all=pop_pca_scores_ht.training_pop
    )
    if withhold_prop:
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                hl.is_defined(pop_pca_scores_ht.training_pop)
                & hl.rand_bool(1.0 - withhold_prop),
                pop_pca_scores_ht.training_pop,
            )
        )

    logger.info(
        "Running RF using {} training examples".format(
            pop_pca_scores_ht.aggregate(
                hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
            )
        )
    )

    pop_ht, pops_rf_model = assign_population_pcs(
        pop_pca_scores_ht,
        pc_cols=pop_pca_scores_ht.scores[:n_pcs],
        known_col="training_pop",
        min_prob=min_prob,
    )

    n_mislabeled_samples = pop_ht.aggregate(
        hl.agg.count_where(pop_ht.training_pop != pop_ht.pop)
    )
    pop_ht = pop_ht.annotate(
        training_pop_all=pop_pca_scores_ht[pop_ht.key].training_pop_all
    )
    pop_assignment_iter = 1
    while n_mislabeled_samples > max_mislabeled_training_samples:
        pop_assignment_iter += 1
        logger.info(
            f"Found {n_mislabeled_samples} samples labeled differently from their known pop. Re-running without."
        )

        pop_ht = pop_ht[pop_pca_scores_ht.key]
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                (pop_ht.training_pop == pop_ht.pop), pop_pca_scores_ht.training_pop
            ),
            training_pop_all=hl.or_missing(
                hl.is_missing(pop_ht.training_pop)
                | (pop_ht.training_pop == pop_ht.pop),
                pop_pca_scores_ht.training_pop_all,
            ),
        ).persist()

        logger.info(
            "Running RF using {} training examples".format(
                pop_pca_scores_ht.aggregate(
                    hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
                )
            )
        )

        pop_ht, pops_rf_model = assign_population_pcs(
            pop_pca_scores_ht,
            pc_cols=pop_pca_scores_ht.scores[:n_pcs],
            known_col="training_pop",
            min_prob=min_prob,
        )
        pop_ht = pop_ht.annotate(
            training_pop_all=pop_pca_scores_ht[pop_ht.key].training_pop_all
        )

        n_mislabeled_samples = pop_ht.aggregate(
            hl.agg.count_where(pop_ht.training_pop != pop_ht.pop)
        )

    pop_ht = pop_ht.annotate_globals(
        min_prob=min_prob,
        include_unreleasable_samples=include_unreleasable_samples,
        max_mislabeled_training_samples=max_mislabeled_training_samples,
        pop_assignment_iterations=pop_assignment_iter,
        n_pcs=n_pcs,
    )
    if withhold_prop:
        pop_ht = pop_ht.annotate_globals(withhold_prop=withhold_prop)

    return pop_ht, pops_rf_model


def apply_stratified_filters(
    sample_qc_ht: hl.Table, filtering_qc_metrics: List[str]
) -> hl.Table:
    """
    Use population stratified QC metrics to determine what samples are outliers and should be filtered.

    Use `compute_stratified_metrics_filter` to run hl.sample_qc and return the metrics stratified by assigned
    population with a `qc_metrics_filters` annotation indicating if the sample falls a certain number of MAD
    outside the distribution.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4, we modify this for n_singleton to be
        (4.0, 8.0)

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :return: Table with stratified metrics and filters
    :rtype: hl.Table
    """
    logger.info(
        "Computing stratified QC metrics filters using metrics: "
        + ", ".join(filtering_qc_metrics)
    )
    sample_qc_ht = sample_qc_ht.annotate(qc_pop=pop.ht()[sample_qc_ht.key].pop)
    sample_qc_ht = sample_qc_ht.filter(
        hl.is_missing(hard_filtered_samples.ht()[sample_qc_ht.key])
    )
    stratified_metrics_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={
            metric: sample_qc_ht.sample_qc[metric] for metric in filtering_qc_metrics
        },
        strata={"qc_pop": sample_qc_ht.qc_pop},
        metric_threshold={"n_singleton": (4.0, 8.0)},
    )
    return stratified_metrics_ht


def apply_regressed_filters(
    sample_qc_ht: hl.Table,
    filtering_qc_metrics: List[str],
    include_unreleasable_samples: bool,
    n_pcs: int = 16,
) -> hl.Table:
    """
    Compute sample QC metrics residuals after regressing out population PCs and determine what samples are outliers
    that should be filtered.

    .. note::

        Default left and right MAD for `compute_stratified_metrics_filter` is 4. We modify this for n_singleton_residual
        to be (math.inf, 8.0) and for r_het_hom_var_residual to be (math.inf, 4.0). The math.inf is used to prevent a
        lower cutoff for these metrics.

    :param sample_qc_ht: Sample QC HT
    :param filtering_qc_metrics: Specific metrics to compute
    :param include_unreleasable_samples: Should unreleasable samples be included in the regression
    :param n_pcs: Number of population PCs to use in regression
    :return: Table with stratified metrics and filters
    :rtype: hl.Table
    """
    project_ht = project_meta.ht()
    sample_qc_ht = sample_qc_ht.select(
        **sample_qc_ht.sample_qc,
        **ancestry_pca_scores(include_unreleasable_samples).ht()[sample_qc_ht.key],
        releasable=project_ht[sample_qc_ht.key].releasable,
        exclude=project_ht[sample_qc_ht.key].exclude,
    )
    residuals_ht = compute_qc_metrics_residuals(
        ht=sample_qc_ht,
        pc_scores=sample_qc_ht.scores[:n_pcs],
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        regression_sample_inclusion_expr=sample_qc_ht.releasable
        & ~sample_qc_ht.exclude,
    )
    residuals_ht = residuals_ht.filter(
        hl.is_missing(hard_filtered_samples.ht()[residuals_ht.key])
    )
    stratified_metrics_ht = compute_stratified_metrics_filter(
        ht=residuals_ht,
        qc_metrics=dict(residuals_ht.row_value),
        metric_threshold={
            "n_singleton_residual": (math.inf, 8.0),
            "r_het_hom_var_residual": (math.inf, 4.0),
        },
    )

    residuals_ht = residuals_ht.annotate(**stratified_metrics_ht[residuals_ht.key])
    residuals_ht = residuals_ht.annotate_globals(
        **stratified_metrics_ht.index_globals(), n_pcs=n_pcs,
    )

    return residuals_ht


def get_relatedness_set_ht(relatedness_ht: hl.Table) -> hl.Table:
    """
    Parse relatedness Table to get every relationship (except UNRELATED) per sample.

    Return Table keyed by sample with all sample relationships in a set.

    :param Table relatedness_ht: Table with inferred relationship information output by pc_relate.
        Keyed by sample pair (i, j).
    :return: Table keyed by sample (s) with all relationships annotated as a set.
    :rtype: hl.Table
    """
    relatedness_ht = relatedness_ht.filter(relatedness_ht.relationship != UNRELATED)
    relatedness_ht = relatedness_ht.select("relationship", s=relatedness_ht.i.s).union(
        relatedness_ht.select("relationship", s=relatedness_ht.j.s)
    )
    relatedness_ht = relatedness_ht.group_by(relatedness_ht.s).aggregate(
        relationships=hl.agg.collect_as_set(relatedness_ht.relationship)
    )
    return relatedness_ht


def get_relationship_filter_expr(
    hard_filtered_expr: hl.expr.BooleanExpression,
    relationship: str,
    relationship_set: hl.expr.SetExpression,
) -> hl.expr.builders.CaseBuilder:
    """
    Return case statement to populate relatedness filters in sample_filters struct

    :param hl.expr.BooleanExpression hard_filtered_expr: Boolean for whether sample was hard filtered.
    :param str relationship: Relationship to check for. One of DUPLICATE_OR_TWINS, PARENT_CHILD, or SIBLINGS.
    :param hl.expr.SetExpression relationship_set: Set containing all possible relationship strings for sample.
    :return: Case statement used to population sample_filters related filter field.
    :rtype: hl.expr.builders.CaseBuilder
    """
    return (
        hl.case()
        .when(hard_filtered_expr, hl.null(hl.tbool))
        .when(hl.is_defined(relationship_set), relationship_set.contains(relationship))
        .default(False)
    )


def join_tables(
    left_ht: hl.Table,
    left_key: str,
    right_ht: hl.Table,
    right_key: str,
    join_type: str,
    sample_count_match: bool = True,
) -> hl.Table:
    """
    Joins left and right tables using specified keys and join types and returns result.

    Also prints warning if sample counts are not the same.

    :param Table left_ht: Left Table to be joined
    :param str left_key: Key of left Table
    :param Table right_ht: Right Table to be joined
    :param str right_key: Key of right Table
    :param str join_type: Type of join
    :param bool sample_count_match: Are the sample counts expected to match in the tables
    :return: Table with annotations
    :rtype: Table
    """
    if sample_count_match and not compare_row_counts(left_ht, right_ht):
        logger.warning("Sample counts in left and right tables do not match!")

        in_left_not_right = left_ht.anti_join(right_ht)
        if in_left_not_right.count() != 0:
            logger.warning(
                f"The following {in_left_not_right.count()} samples are found in the left HT, but are not found in the right HT"
            )
            in_left_not_right.select().show(n=-1)

        in_right_not_left = right_ht.anti_join(left_ht)
        if in_right_not_left.count() != 0:
            logger.warning(
                f"The following {in_right_not_left.count()} samples are found in the right HT, but are not found in left HT"
            )
            in_right_not_left.select().show(n=-1)

        if join_type != "outer":
            logger.warning(
                "Join type is not an outer join so some samples will be filtered!"
            )

    return left_ht.key_by(left_key).join(right_ht.key_by(right_key), how=join_type)


def generate_metadata(regressed_metrics_outlier: bool = True) -> hl.Table:
    """
    Pull all sample QC information together in one Table.

    :param regressed_metrics_outlier: Should the outlier table be from the residuals instead of pop stratified
    :return: Table annotated with all metadata
    :rtype: hl.Table
    """
    logging_statement = "Reading in {} and joining with meta HT"

    logger.info(
        "Loading metadata file with subset, age, and releasable information to begin creation of the meta HT"
    )
    left_ht = get_gnomad_v3_mt(remove_hard_filtered_samples=False).cols()

    right_ht = project_meta.ht()
    right_ht = right_ht.select(
        project_meta=hl.struct(**right_ht.row.drop(*(SUBSETS + ["s"]))),
        subsets=hl.struct(**{x: right_ht[x] for x in SUBSETS}),
    )
    left_ht = join_tables(left_ht, "s", right_ht.select_globals(), "s", "right")

    logger.info(logging_statement.format("picard metric HT"))
    right_ht = picard_metrics.ht()
    right_ht = right_ht.select("bam_metrics")
    left_ht = join_tables(left_ht, "s", right_ht, "s", "left", sample_count_match=False)

    logger.info(logging_statement.format("sex HT"))
    impute_stats = ["f_stat", "n_called", "expected_homs", "observed_homs"]
    right_ht = sex.ht()
    right_ht = right_ht.transmute(
        impute_sex_stats=hl.struct(**{x: right_ht[x] for x in impute_stats})
    )

    # Create struct for join
    right_ht = right_ht.select(sex_imputation=hl.struct(**right_ht.row.drop("s")))
    right_ht = right_ht.select_globals(sex_imputation_ploidy_cutoffs=right_ht.globals)
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info(logging_statement.format("sample QC HT"))
    right_ht = get_sample_qc("bi_allelic").ht()

    # Remove annotations that cannot be computed from the sparse format
    not_in_sparse = ["n_called", "n_not_called", "n_filtered", "call_rate"]
    right_ht = right_ht.annotate(
        **{x: right_ht[x].drop(*not_in_sparse) for x in right_ht.row_value}
    )
    left_ht = join_tables(left_ht, "s", right_ht, "s", "right")

    logger.info(logging_statement.format("population PCA HT"))
    right_ht = pop.ht()
    right_ht = right_ht.select_globals(
        population_inference_pca_metrics=right_ht.globals
    )
    right_ht = right_ht.select(population_inference=hl.struct(**right_ht.row.drop("s")))
    left_ht = join_tables(
        left_ht, "s", right_ht, "s", "outer", sample_count_match=False
    )

    logger.info(
        "Reading hard filters HT, renaming hard filters struct to sample_filters, and joining with meta HT"
    )
    right_ht = hard_filtered_samples.ht()
    left_ht = join_tables(
        left_ht, "s", right_ht, "s", "outer", sample_count_match=False
    )

    # Change sample_filters to a struct
    ex_right_ht = right_ht.explode(right_ht.hard_filters)
    hard_filters = ex_right_ht.aggregate(
        hl.agg.collect_as_set(ex_right_ht.hard_filters)
    )
    left_ht = left_ht.transmute(
        sample_filters=hl.struct(
            **{
                v: hl.if_else(
                    hl.is_defined(left_ht.hard_filters),
                    left_ht.hard_filters.contains(v),
                    False,
                )
                for v in hard_filters
            },
            hard_filters=left_ht.hard_filters,
            hard_filtered=hl.is_defined(left_ht.hard_filters)
            & (hl.len(left_ht.hard_filters) > 0),
        )
    )

    logger.info(
        "Reading in PCA related samples to drop HT and preparing to annotate meta HT's sample_filter struct with relatedness booleans"
    )
    related_samples_to_drop_ht = pca_related_samples_to_drop.ht()
    release_related_samples_to_drop_ht = release_related_samples_to_drop.ht()
    relatedness_ht = get_relatedness_set_ht(relatedness.ht())
    related_samples_to_drop_ht = related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[related_samples_to_drop_ht.s].relationships
    )
    release_related_samples_to_drop_ht = release_related_samples_to_drop_ht.annotate(
        relationships=relatedness_ht[release_related_samples_to_drop_ht.s].relationships
    )

    # Annotating meta HT with related filter booleans
    # Any sample that is hard filtered will have missing values for these bools
    # Any sample that was filtered for relatedness will have True for sample_filters.related
    # If a filtered related sample had a relationship with a higher degree than second-degree (duplicate, parent-child, sibling),
    # that filter will also be True
    release_else_expr = release_related_samples_to_drop_ht[left_ht.key].relationships
    all_else_expr = related_samples_to_drop_ht[left_ht.key].relationships
    left_ht = left_ht.annotate(
        sample_filters=left_ht.sample_filters.annotate(
            release_related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(release_related_samples_to_drop_ht[left_ht.key]),
            ),
            release_duplicate=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered,
                DUPLICATE_OR_TWINS,
                release_else_expr,
            ),
            release_parent_child=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered, PARENT_CHILD, release_else_expr,
            ),
            release_sibling=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered, SIBLINGS, release_else_expr,
            ),
            all_samples_related=hl.if_else(
                left_ht.sample_filters.hard_filtered,
                hl.null(hl.tbool),
                hl.is_defined(related_samples_to_drop_ht[left_ht.key]),
            ),
            all_samples_duplicate=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered, DUPLICATE_OR_TWINS, all_else_expr,
            ),
            all_samples_parent_child=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered, PARENT_CHILD, all_else_expr,
            ),
            all_samples_sibling=get_relationship_filter_expr(
                left_ht.sample_filters.hard_filtered, SIBLINGS, all_else_expr,
            ),
        )
    )
    left_ht = left_ht.annotate(
        relatedness_inference=hl.struct(
            relationships=relatedness_ht[left_ht.s].relationships,
        )
    )

    logger.info("Adding relatedness globals (cutoffs)")
    left_ht = left_ht.annotate_globals(
        relatedness_inference_cutoffs=hl.struct(**relatedness_ht.index_globals())
    )

    logger.info(logging_statement.format("outlier HT"))
    if regressed_metrics_outlier:
        right_ht = regressed_metrics.ht()
    else:
        right_ht = stratified_metrics.ht()

    right_ht = right_ht.select_globals(
        outlier_detection_metrics=hl.struct(
            **right_ht.index_globals(), used_regressed_metrics=regressed_metrics_outlier
        )
    )

    left_ht = join_tables(left_ht, "s", right_ht, "s", "outer")
    left_ht = left_ht.transmute(
        sample_filters=left_ht.sample_filters.annotate(
            **{x: left_ht[x] for x in left_ht.row if x.startswith("fail_")},
            qc_metrics_filters=left_ht.qc_metrics_filters,
        )
    )
    if regressed_metrics_outlier:
        left_ht = left_ht.transmute(
            sample_qc=left_ht.sample_qc.annotate(
                **{x: left_ht[x] for x in left_ht.row if x.endswith("_residual")},
            )
        )

    logger.info("Annotating high_quality field")
    left_ht = left_ht.annotate(
        high_quality=~left_ht.sample_filters.hard_filtered
        & (hl.len(left_ht.sample_filters.qc_metrics_filters) == 0)
    )

    logger.info("Annotating releasable field")
    left_ht = left_ht.annotate(
        release=left_ht.project_meta.releasable
        & left_ht.high_quality
        & ~left_ht.project_meta.exclude
        & ~left_ht.sample_filters.release_related
    ).persist()

    logger.info(
        "Release sample counts:"
        f"{left_ht.aggregate(hl.struct(release=hl.agg.count_where(left_ht.release)))}"
    )
    left_ht.describe()
    left_ht.summarize()
    logger.info(f"Final count: {left_ht.count()}")
    logger.info("Complete")

    return left_ht


def main(args):
    hl.init(log="/hail.log", default_reference="GRCh38")

    if args.sample_qc:
        compute_sample_qc().write(get_sample_qc().path, overwrite=args.overwrite)

    if args.impute_sex:
        compute_sex().write(sex.path, overwrite=args.overwrite)
    elif args.reannotate_sex:
        reannotate_sex(
            args.min_cov,
            (args.upper_x, (args.lower_xx, args.upper_xx), args.lower_xxx),
            ((args.lower_y, args.upper_y), args.lower_yy),
        ).write(sex.path, overwrite=args.overwrite)

    if args.compute_hard_filters:
        compute_hard_filters(args.min_cov).write(
            hard_filtered_samples.path, overwrite=args.overwrite
        )

    if args.compute_qc_mt:
        compute_qc_mt().write(qc.path, overwrite=args.overwrite)

    if args.run_pc_relate or args.reannotate_relatedness:
        if args.run_pc_relate:
            logger.info("Running PC-Relate")
            logger.warning(
                "PC-relate requires SSDs and doesn't work with preemptible workers!"
            )
            qc_mt = qc.mt()
            eig, scores, _ = hl.hwe_normalized_pca(
                qc_mt.GT, k=10, compute_loadings=False
            )
            scores = scores.checkpoint(
                pc_relate_pca_scores.path,
                overwrite=args.overwrite,
                _read_if_exists=not args.overwrite,
            )
            relatedness_ht = hl.pc_relate(
                qc_mt.GT,
                min_individual_maf=0.01,
                scores_expr=scores[qc_mt.col_key].scores,
                block_size=4096,
                min_kinship=0.05,
                statistics="all",
            )

        else:
            relatedness_ht = relatedness.ht().checkpoint(
                "gs://gnomad-tmp/relatedness_ht_checkpoint.ht", overwrite=True
            )  # Copy HT to temp location to overwrite annotation
        relatedness_ht = relatedness_ht.annotate(
            relationship=get_relationship_expr(
                kin_expr=relatedness_ht.kin,
                ibd0_expr=relatedness_ht.ibd0,
                ibd1_expr=relatedness_ht.ibd1,
                ibd2_expr=relatedness_ht.ibd2,
                first_degree_kin_thresholds=tuple(args.first_degree_kin_thresholds),
                second_degree_min_kin=args.second_degree_kin_cutoff,
                ibd0_0_max=args.ibd0_0_max,
            )
        )
        relatedness_ht = relatedness_ht.annotate_globals(
            min_individual_maf=0.01,
            min_emission_kinship=0.05,
            ibd0_0_max=args.ibd0_0_max,
            second_degree_kin_cutoff=args.second_degree_kin_cutoff,
            first_degree_kin_thresholds=tuple(args.first_degree_kin_thresholds),
        )
        relatedness_ht.write(relatedness.path, args.overwrite)

    if args.run_pca:
        rank_ht = compute_sample_rankings(
            use_qc_metrics_filters=False
        )  # QC metrics filters do not exist at this point
        rank_ht = rank_ht.checkpoint(
            pca_samples_rankings.path,
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        filtered_samples = hl.literal(
            rank_ht.aggregate(
                hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
            )
        )  # TODO: don't localize once hail bug is fixed
        samples_to_drop = compute_related_samples_to_drop(
            relatedness.ht(),
            rank_ht,
            args.second_degree_kin_cutoff,
            filtered_samples=filtered_samples,
        )
        samples_to_drop = samples_to_drop.key_by(s=samples_to_drop.s.s)
        samples_to_drop.checkpoint(
            pca_related_samples_to_drop.path,
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        pop_pca_eigenvalues, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca(
            args.include_unreleasable_samples, args.n_pcs, samples_to_drop
        )
        pop_pca_scores_ht.write(
            ancestry_pca_scores(args.include_unreleasable_samples).path,
            overwrite=args.overwrite,
        )
        pop_pca_loadings_ht.write(
            ancestry_pca_loadings(args.include_unreleasable_samples).path,
            overwrite=args.overwrite,
        )

        pop_pca_eigenvalues_ht = hl.Table.parallelize(
            hl.literal(
                [
                    {"PC": i + 1, "eigenvalue": x}
                    for i, x in enumerate(pop_pca_eigenvalues)
                ],
                "array<struct{PC: int, eigenvalue: float}>",
            )
        )
        pop_pca_eigenvalues_ht.write(
            ancestry_pca_eigenvalues(args.include_unreleasable_samples).path,
            overwrite=args.overwrite,
        )

    if args.assign_pops:
        n_pcs = args.pop_n_pcs
        pop_ht, pops_rf_model = assign_pops(
            args.min_pop_prob,
            args.include_unreleasable_samples,
            n_pcs=n_pcs,
            withhold_prop=args.withhold_prop,
        )
        pop_ht = pop_ht.checkpoint(
            pop.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite
        )
        pop_ht.transmute(
            **{f"PC{i + 1}": pop_ht.pca_scores[i] for i in range(n_pcs)}
        ).export(pop_tsv_path())

        with hl.hadoop_open(pop_rf_path(), "wb") as out:
            pickle.dump(pops_rf_model, out)

    if args.calculate_inbreeding:
        qc_mt = qc.mt()
        pop_ht = pop.ht()
        qc_mt = qc_mt.annotate_cols(pop=pop_ht[qc_mt.col_key].pop)
        qc_mt = qc_mt.annotate_rows(
            call_stats_by_pop=hl.agg.group_by(
                qc_mt.pop, hl.agg.call_stats(qc_mt.GT, qc_mt.alleles)
            )
        )
        inbreeding_ht = (
            qc_mt.annotate_cols(
                inbreeding=hl.agg.inbreeding(
                    qc_mt.GT, qc_mt.call_stats_by_pop[qc_mt.pop].AF[1]
                )
            )
            .cols()
            .select("inbreeding")
        )
        inbreeding_ht.write(sample_inbreeding.path, overwrite=args.overwrite)

    if args.calculate_clinvar:
        mt = get_gnomad_v3_mt(
            split=True, key_by_locus_and_alleles=True, remove_hard_filtered_samples=True
        )
        clinvar_ht = clinvar.ht()
        mt = mt.filter_rows(hl.is_defined(clinvar_ht[mt.row_key]))
        mt = mt.checkpoint("gs://gnomad-tmp/clinvar_variants.mt", overwrite=True)
        mt = mt.annotate_rows(
            clinvar_path=hl.is_defined(
                filter_to_clinvar_pathogenic(clinvar_ht)[mt.row_key]
            )
        )
        clinvar_sample_ht = mt.select_cols(
            n_clinvar=hl.agg.count_where(mt.GT.is_non_ref()),
            n_clinvar_path=hl.agg.count_where(mt.GT.is_non_ref() & mt.clinvar_path),
        ).cols()
        clinvar_sample_ht.write(sample_clinvar_count.path, overwrite=args.overwrite)

    if args.apply_stratified_filters or args.apply_regressed_filters:
        filtering_qc_metrics = args.filtering_qc_metrics.split(",")
        sample_qc_ht = get_sample_qc("bi_allelic").ht()

        if "inbreeding" in filtering_qc_metrics:
            inbreeding_ht = sample_inbreeding.ht()[sample_qc_ht.key]
            sample_qc_ht = sample_qc_ht.annotate(
                sample_qc=sample_qc_ht.sample_qc.annotate(
                    inbreeding=inbreeding_ht.inbreeding.f_stat
                )
            )

        if args.apply_regressed_filters:
            n_pcs = args.regress_n_pcs
            apply_regressed_filters(
                sample_qc_ht,
                filtering_qc_metrics,
                args.include_unreleasable_samples,
                n_pcs,
            ).write(regressed_metrics.path, overwrite=args.overwrite)
        else:
            apply_stratified_filters(sample_qc_ht, filtering_qc_metrics,).write(
                stratified_metrics.path, overwrite=args.overwrite
            )

    if args.compute_related_samples_to_drop:
        rank_ht = compute_sample_rankings(use_qc_metrics_filters=True)
        rank_ht = rank_ht.checkpoint(
            release_samples_rankings.path,
            overwrite=args.overwrite,
            _read_if_exists=not args.overwrite,
        )
        filtered_samples = hl.literal(
            rank_ht.aggregate(
                hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
            )
        )  # TODO: don't localize once hail bug is fixed
        print(filtered_samples)
        relatedness_ht = relatedness.ht()
        relatedness_ht = relatedness_ht.key_by(
            i=relatedness_ht.i.s, j=relatedness_ht.j.s
        )
        samples_to_drop = compute_related_samples_to_drop(
            relatedness_ht,
            rank_ht,
            args.second_degree_kin_cutoff,
            filtered_samples=filtered_samples,
        )
        samples_to_drop.write(
            release_related_samples_to_drop.path, overwrite=args.overwrite
        )

    if args.generate_metadata:
        meta_ht = generate_metadata(args.regressed_metrics_outlier)
        meta_ht.checkpoint(
            meta.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite
        )
        n_pcs = meta_ht.aggregate(
            hl.agg.min(hl.len(meta_ht.population_inference.pca_scores))
        )

        meta_ht = meta_ht.annotate(
            population_inference=meta_ht.population_inference.transmute(
                **{
                    f"PC{i + 1}": meta_ht.population_inference.pca_scores[i]
                    for i in range(n_pcs)
                }
            ),
            hard_filters=hl.or_missing(
                hl.len(meta_ht.sample_filters.hard_filters) > 0,
                hl.delimit(meta_ht.sample_filters.hard_filters),
            ),
            qc_metrics_filters=hl.or_missing(
                hl.len(meta_ht.sample_filters.qc_metrics_filters) > 0,
                hl.delimit(meta_ht.sample_filters.qc_metrics_filters),
            ),
        )
        meta_ht.flatten().export(meta_tsv_path())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--sample_qc", help="Assigns pops from PCA", action="store_true"
    )
    parser.add_argument(
        "--impute_sex",
        help="Runs sex imputation. Also runs sex karyotyping annotation.",
        action="store_true",
    )
    parser.add_argument(
        "--reannotate_sex",
        help="Runs the sex karyotyping annotations again, without re-computing sex imputation metrics.",
        action="store_true",
    )
    parser.add_argument("--upper_x", help="Upper cutoff for single X", type=float)
    parser.add_argument("--lower_xx", help="Lower cutoff for double X", type=float)
    parser.add_argument("--upper_xx", help="Upper cutoff for double X", type=float)
    parser.add_argument("--lower_xxx", help="Lower cutoff for triple X", type=float)
    parser.add_argument("--lower_y", help="Lower cutoff for single Y", type=float)
    parser.add_argument("--upper_y", help="Upper cutoff for single Y", type=float)
    parser.add_argument("--lower_yy", help="Lower cutoff for double Y", type=float)
    parser.add_argument(
        "--compute_hard_filters",
        help="Computes samples to be hard-filtered",
        action="store_true",
    )
    parser.add_argument(
        "--min_cov",
        help="Minimum coverage for inclusion when computing hard-filters",
        default=15,
        type=int,
    )
    parser.add_argument(
        "--compute_qc_mt",
        help="Creates the QC MT based on liftover of v2 QC and Purcell 5k sites",
        action="store_true",
    )
    parser.add_argument("--run_pc_relate", help="Run PC-relate", action="store_true")
    parser.add_argument(
        "--reannotate_relatedness",
        help="Runs the relatedness annotation without re-running pc-relate",
        action="store_true",
    )
    parser.add_argument(
        "--first_degree_kin_thresholds",
        help="First degree kinship threshold for filtering a pair of samples with a first degree relationship. \
        Default = (0.1767767, 0.4); \
        Defaults taken from Bycroft et al. (2018)",
        nargs=2,
        default=(0.1767767, 0.4),
        type=float,
    )
    parser.add_argument(
        "--second_degree_kin_cutoff",
        help="Minimum kinship threshold for filtering a pair of samples with a second degree relationship\
        in PC relate and filtering related individuals. (Default = 0.1) \
        Bycroft et al. (2018) calculates 0.08838835 but from evaluation of the distributions v3 has used 0.1",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--ibd0_0_max",
        help="IBD0 cutoff to determine parent offspring vs full sibling (Default = 0.05) \
        Default is adjusted from theoretical values; parent-offspring should have an IBD0 = 0. \
        Full siblings should have an IBD0 = 0.25.",
        default=0.05,
    )
    parser.add_argument(
        "--compute_related_samples_to_drop",
        help="Flags related samples to drop",
        action="store_true",
    )
    parser.add_argument(
        "--min_related_hard_filter",
        help="Minimum number of relateds to have to get hard-filterd",
        default=50,
        type=int,
    )
    parser.add_argument("--run_pca", help="Compute PCA", action="store_true")
    parser.add_argument(
        "--n_pcs",
        help="Number of PCs to compute for ancestry PCA",
        default=30,
        type=int,
    )
    parser.add_argument(
        "--include_unreleasable_samples",
        help="Includes unreleasable samples for computing PCA",
        action="store_true",
    )
    parser.add_argument(
        "--assign_pops", help="Assigns pops from PCA", action="store_true"
    )
    parser.add_argument(
        "--pop_n_pcs",
        help="Number of PCs to use for ancestry assignment",
        default=16,
        type=int,
    )
    parser.add_argument(
        "--min_pop_prob",
        help="Minimum RF prob for pop assignment",
        default=0.75,
        type=float,
    )
    parser.add_argument(
        "--withhold_prop",
        help="Proportion of training samples to withhold from pop assignment RF training",
        type=float,
    )
    parser.add_argument(
        "--calculate_inbreeding",
        help="Calculate sample level inbreeding",
        action="store_true",
    )
    parser.add_argument(
        "--calculate_clinvar",
        help="Calculate counts of ClinVar and ClinVar P/LP variants per sample",
        action="store_true",
    )
    parser.add_argument(
        "--filtering_qc_metrics",
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
    )  # used in v3 'n_het', 'n_hom_var',
    parser.add_argument(
        "--apply_stratified_filters",
        help="Compute per pop filtering.",
        action="store_true",
    )
    parser.add_argument(
        "--apply_regressed_filters",
        help="Computes qc_metrics adjusted for pop.",
        action="store_true",
    )
    parser.add_argument(
        "--regress_n_pcs",
        help="Number of PCs to use for qc metric regressions",
        default=10,
        type=int,
    )
    parser.add_argument(
        "--generate_metadata", help="Generates the metadata HT.", action="store_true"
    )
    parser.add_argument(
        "--regressed_metrics_outlier",
        help="Should metadata HT use regression outlier model.",
        action="store_true",
    )

    main(parser.parse_args())
