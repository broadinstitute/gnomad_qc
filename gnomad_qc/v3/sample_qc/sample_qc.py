import argparse
import logging
import pickle
from typing import Any, List, Tuple

import hail as hl
from gnomad.resources.grch38 import (lcr_intervals, purcell_5k_intervals,
                                     telomeres_and_centromeres)
from gnomad.sample_qc.ancestry import (assign_population_pcs,
                                       run_pca_with_relateds)
from gnomad.sample_qc.filtering import (compute_qc_metrics_residuals,
                                        compute_stratified_metrics_filter,
                                        compute_stratified_sample_qc)
from gnomad.sample_qc.pipeline import annotate_sex, get_qc_mt
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.sample_qc.sex import get_ploidy_cutoffs, get_sex_expr
from gnomad.utils.annotations import bi_allelic_expr, get_adj_expr
from gnomad.utils.filtering import add_filters_expr, filter_to_autosomes
from gnomad.utils.sparse_mt import densify_sites

from gnomad_qc.v2.resources.sample_qc import get_liftover_v2_qc_mt
from gnomad_qc.v3.resources.annotations import freq, get_info, last_END_position
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta, meta_tsv_path, project_meta
from gnomad_qc.v3.resources.sample_qc import (ancestry_pca_eigenvalues_path,
                                              ancestry_pca_loadings,
                                              ancestry_pca_scores,
                                              get_sample_qc,
                                              hard_filtered_samples,
                                              pc_relate_pca_scores,
                                              pca_related_samples_to_drop,
                                              pca_samples_rankings,
                                              picard_metrics, pop, pop_rf_path,
                                              pop_tsv_path, qc,
                                              regressed_metrics, relatedness,
                                              release_related_samples_to_drop,
                                              release_samples_rankings,
                                              sample_inbreeding, sex,
                                              stratified_metrics)

logger = logging.getLogger("sample_qc")


def compute_sample_qc() -> hl.Table:
    logger.info("Computing sample QC")
    mt = filter_to_autosomes(
        get_gnomad_v3_mt(
            split=True,
            key_by_locus_and_alleles=True,
            remove_hard_filtered_samples=False
        )
    )
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = mt.select_entries('GT')

    sample_qc_ht = compute_stratified_sample_qc(
        mt,
        strata={
            'bi_allelic': bi_allelic_expr(mt),
            'multi_allelic': ~bi_allelic_expr(mt)
        },
        tmp_ht_prefix=get_sample_qc().path[:-3],
    )

    # Remove annotations that cannot be computed from the sparse format
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            x: sample_qc_ht[x].drop('n_called', 'n_not_called', 'n_filtered', 'call_rate')
            for x in sample_qc_ht.row_value
        }
    )
    return sample_qc_ht.repartition(100)


def compute_qc_mt() -> hl.MatrixTable:
    # Load v2 and p5k sites for QC
    v2_qc_sites = get_liftover_v2_qc_mt('joint', ld_pruned=True).rows().key_by('locus')
    qc_sites = v2_qc_sites.union(purcell_5k_intervals.ht(), unify=True)

    qc_sites = qc_sites.filter(
        hl.is_missing(lcr_intervals.ht()[qc_sites.key])
    )

    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
    mt = mt.select_entries(
        'END',
        GT=mt.LGT,
        adj=get_adj_expr(
            mt.LGT,
            mt.GQ,
            mt.DP,
            mt.LAD
        )
    )
    mt = densify_sites(
        mt,
        qc_sites,
        hl.read_table(last_END_position.path)
    )

    mt = mt.filter_rows(
        (hl.len(mt.alleles) == 2) &
        hl.is_snp(mt.alleles[0], mt.alleles[1]) &
        (qc_sites[mt.row_key].alleles == mt.alleles)

    )
    mt = mt.checkpoint('gs://gnomad-tmp/gnomad_v3_qc_mt_v2_sites_dense.mt', overwrite=True)
    mt = mt.naive_coalesce(5000)
    mt = mt.checkpoint('gs://gnomad-tmp/gnomad_v3_qc_mt_v2_sites_dense_repartitioned.mt', overwrite=True)
    info_ht = get_info(split=False).ht()
    info_ht = info_ht.annotate(
        info=info_ht.info.select(
            # No need for AS_annotations since it's bi-allelic sites only
            **{x: info_ht.info[x] for x in info_ht.info if not x.startswith('AS_')}
        )
    )
    mt = mt.annotate_rows(
        info=info_ht[mt.row_key].info
    )
    qc_mt = get_qc_mt(
        mt,
        min_af=0.0,
        min_inbreeding_coeff_threshold=-0.025,
        min_hardy_weinberg_threshold=None,
        ld_r2=None,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False
    )
    return qc_mt


def compute_hard_filters(cov_threshold: int) -> hl.Table:
    ht = get_gnomad_v3_mt(remove_hard_filtered_samples=False).cols()
    hard_filters = dict()

    # Remove samples failing fingerprinting
    # TODO: Add these into hard filtering metadata when incorporating internal smaples Picard metrics
    hard_filters['failed_fingerprinting'] = hl.array(['09C90823', '10C103592', 'S5530']).contains(ht.s)

    # Remove TCGA tumor samples based on TCGA naming convention: https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
    hard_filters['TCGA_tumor_sample'] = (
        (ht.s.startswith('TCGA') &
         (hl.int(hl.str(ht.s).split("-")[3][:2]) < 10))
    )

    # Remove low-coverage samples
    cov_ht = sex.ht()  # chrom 20 coverage is computed to infer sex and used here
    hard_filters['low_coverage'] = (cov_ht[ht.key].chr20_mean_dp < cov_threshold)

    # Remove extreme raw bi-allelic sample QC outliers
    # These were determined by visual inspection of the metrics in gs://gnomad/sample_qc/  v3_genomes_sample_qc.ipynb
    bi_allelic_qc_ht = hl.read_table('gs://gnomad/sample_qc/ht/genomes_v3/sample_qc_bi_allelic.ht')[ht.key]
    hard_filters['bad_qc_metrics'] = (
            (bi_allelic_qc_ht.sample_qc.n_snp > 3.75e6) |
            (bi_allelic_qc_ht.sample_qc.n_snp < 2.4e6) |
            (bi_allelic_qc_ht.sample_qc.n_singleton > 1e5) |
            (bi_allelic_qc_ht.sample_qc.r_het_hom_var > 3.3)
    )

    # Remove samples with ambiguous sex assignments
    sex_ht = sex.ht()[ht.key]
    hard_filters['ambiguous_sex'] = (sex_ht.sex_karyotype == 'Ambiguous')
    hard_filters['sex_aneuploidy'] = ~hl.set({'Ambiguous', 'XX', 'XY'}).contains(sex_ht.sex_karyotype)

    # Remove samples that fail picard metric thresholds, percents are not divided by 100, e.g. 5% == 5.00, %5 != 0.05
    picard_ht = picard_metrics.ht()[ht.key]
    hard_filters['contamination'] = picard_ht.bam_metrics.freemix > 5.00
    hard_filters['chimera'] = picard_ht.bam_metrics.pct_chimeras > 5.00
    hard_filters['coverage'] = picard_ht.bam_metrics.mean_coverage < 15
    hard_filters['insert_size'] = picard_ht.bam_metrics.median_insert_size < 250

    ht = ht.annotate(
        hard_filters=add_filters_expr(
            filters=hard_filters
        )
    )

    ht = ht.filter(hl.len(ht.hard_filters) > 0)
    return ht


def compute_sex(aaf_threshold=0.001, f_stat_cutoff=0.5) -> hl.Table:
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=False,
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


def compute_sample_rankings(use_qc_metrics_filters: bool) -> hl.Table:
    project_ht = project_meta.ht()
    project_ht = project_ht.select(
        'releasable',
        chr20_mean_dp=sex.ht()[project_ht.key].chr20_mean_dp,
        filtered=hl.or_else(hl.len(hard_filtered_samples.ht()[project_ht.key].hard_filters) > 0, False)
    )

    if use_qc_metrics_filters:
        project_ht = project_ht.annotate(
            filtered=hl.cond(
                project_ht.filtered,
                True,
                hl.or_else(
                    hl.len(regressed_metrics.ht()[project_ht.key].qc_metrics_filters) > 0,
                    False
                )
            )
        )

    project_ht = project_ht.order_by(
        project_ht.filtered,
        hl.desc(project_ht.releasable),
        hl.desc(project_ht.chr20_mean_dp)
    ).add_index(name='rank')

    return project_ht.key_by('s').select('filtered', 'rank')


def run_pca(
        include_unreleasable_samples: bool,
        n_pcs: int,
        related_samples_to_drop: hl.Table
) -> Tuple[List[float], hl.Table, hl.Table]:
    logger.info("Running population PCA")
    qc_mt = qc.mt()

    samples_to_drop = related_samples_to_drop.select()
    if not include_unreleasable_samples:
        logger.info("Excluding unreleasable samples for PCA.")
        samples_to_drop = samples_to_drop.union(
            qc_mt.filter_cols(~project_meta.ht()[qc_mt.col_key].releasable).cols().select()
        )
    else:
        logger.info("Including unreleasable samples for PCA")

    return run_pca_with_relateds(
        qc_mt,
        samples_to_drop,
        n_pcs=n_pcs
    )


def assign_pops(
        min_prob: float,
        include_unreleasable_samples: bool,
        max_mislabeled_training_samples: int = 50  # TODO: Think about this parameter and add it to assign_population_pcs. Maybe should be a fraction? fraction per pop?
) -> Tuple[hl.Table, Any]:
    logger.info("Assigning global population labels")
    pop_pca_scores_ht = ancestry_pca_scores(include_unreleasable_samples).ht()
    project_meta_ht = project_meta.ht()[pop_pca_scores_ht.key]
    pop_pca_scores_ht = pop_pca_scores_ht.annotate(
        training_pop=(
            hl.case()
                .when(hl.is_defined(project_meta_ht.project_pop), project_meta_ht.project_pop)
                .when(project_meta_ht.v2_pop != 'oth', project_meta_ht.v2_pop)
                .or_missing()
        )
    )

    logger.info("Running RF using {} training examples".format(
        pop_pca_scores_ht.aggregate(
            hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
        )
    )
    )

    pop_ht, pops_rf_model = assign_population_pcs(
        pop_pca_scores_ht,
        pc_cols=pop_pca_scores_ht.scores,
        known_col='training_pop',
        min_prob=min_prob
    )

    n_mislabeled_samples = pop_ht.aggregate(hl.agg.count_where(pop_ht.training_pop != pop_ht.pop))
    while n_mislabeled_samples > max_mislabeled_training_samples:
        logger.info(f"Found {n_mislabeled_samples} samples labeled differently from their known pop. Re-running without.")

        pop_ht = pop_ht[pop_pca_scores_ht.key]
        pop_pca_scores_ht = pop_pca_scores_ht.annotate(
            training_pop=hl.or_missing(
                (pop_ht.training_pop == pop_ht.pop),
                pop_pca_scores_ht.training_pop
            )
        ).persist()

        logger.info("Running RF using {} training examples".format(
            pop_pca_scores_ht.aggregate(
                hl.agg.count_where(hl.is_defined(pop_pca_scores_ht.training_pop))
            )
        )
        )

        pop_ht, pops_rf_model = assign_population_pcs(
            pop_pca_scores_ht,
            pc_cols=pop_pca_scores_ht.scores,
            known_col='training_pop',
            min_prob=min_prob
        )
        n_mislabeled_samples = pop_ht.aggregate(hl.agg.count_where(pop_ht.training_pop != pop_ht.pop))

    return pop_ht, pops_rf_model


def apply_stratified_filters(filtering_qc_metrics: List[str]) -> hl.Table:
    logger.info("Computing stratified QC metrics filters using metrics: " + ", ".join(filtering_qc_metrics))
    sample_qc_ht = hl.read_table(get_sample_qc('bi_allelic'))
    sample_qc_ht = sample_qc_ht.annotate(
        qc_pop=pop.ht()[sample_qc_ht.key].pop
    )
    stratified_metrics_ht = compute_stratified_metrics_filter(
        sample_qc_ht,
        qc_metrics={metric: sample_qc_ht.sample_qc[metric] for metric in filtering_qc_metrics},
        strata={'qc_pop': sample_qc_ht.qc_pop},
        metric_threshold={'n_singleton': (4.0, 8.0)}
    )
    return stratified_metrics_ht


def apply_regressed_filters(
        filtering_qc_metrics: List[str],
        include_unreleasable_samples: bool,
) -> hl.Table:
    sample_qc_ht = get_sample_qc('bi_allelic').ht()
    sample_qc_ht = sample_qc_ht.select(
        **sample_qc_ht.sample_qc,
        **ancestry_pca_scores(include_unreleasable_samples).ht()[sample_qc_ht.key],
        releasable=project_meta.ht()[sample_qc_ht.key].releasable
    )
    residuals_ht = compute_qc_metrics_residuals(
        ht=sample_qc_ht,
        pc_scores=sample_qc_ht.scores,
        qc_metrics={metric: sample_qc_ht[metric] for metric in filtering_qc_metrics},
        regression_sample_inclusion_expr=sample_qc_ht.releasable
    )

    stratified_metrics_ht = compute_stratified_metrics_filter(
        ht=residuals_ht,
        qc_metrics=dict(residuals_ht.row_value),
        metric_threshold={'n_singleton_residual': (4.0, 8.0)}
    )

    residuals_ht = residuals_ht.annotate(
        **stratified_metrics_ht[sample_qc_ht.key]
    )
    residuals_ht = residuals_ht.annotate_globals(
        **stratified_metrics_ht.index_globals()
    )

    return residuals_ht


def generate_metadata() -> hl.Table:
    meta_ht = project_meta.ht()
    sample_qc_ht = get_sample_qc("bi_allelic").ht()
    hard_filters_ht = hard_filtered_samples.ht()
    regressed_metrics_ht = regressed_metrics.ht()
    pop_ht = pop.ht()
    release_related_samples_to_drop_ht = release_related_samples_to_drop.ht()
    pca_related_samples_to_drop_ht = pca_related_samples_to_drop.ht()
    sex_ht = sex.ht()
    sex_ht = sex_ht.transmute(
        impute_sex_stats=hl.struct(
            f_stat=sex_ht.f_stat,
            n_called=sex_ht.n_called,
            expected_homs=sex_ht.expected_homs,
            observed_homs=sex_ht.observed_homs
        )
    )

    meta_ht = meta_ht.annotate_globals(
        **regressed_metrics_ht.index_globals()
    )

    meta_ht = meta_ht.annotate(
        sample_qc=sample_qc_ht[meta_ht.key].sample_qc,
        **hard_filters_ht[meta_ht.key],
        **regressed_metrics_ht[meta_ht.key],
        **pop_ht[meta_ht.key],
        **sex_ht[meta_ht.key],
        release_related=hl.is_defined(release_related_samples_to_drop_ht[meta_ht.key]),
        all_samples_related=hl.is_defined(pca_related_samples_to_drop_ht[meta_ht.key]),
    )

    meta_ht = meta_ht.annotate(
        hard_filters=hl.or_else(meta_ht.hard_filters, hl.empty_set(hl.tstr)),
        sample_filters=add_filters_expr(
            filters={
                'related': meta_ht.release_related
            },
            current_filters=meta_ht.hard_filters.union(meta_ht.qc_metrics_filters)
        )
    )

    meta_ht = meta_ht.annotate(
        high_quality=(hl.len(meta_ht.hard_filters) == 0) & (hl.len(meta_ht.qc_metrics_filters) == 0),
        release=meta_ht.releasable & (hl.len(meta_ht.sample_filters) == 0)
    )

    return meta_ht


def main(args):
    hl.init(log='/hail.log', default_reference='GRCh38')

    if args.sample_qc:
        compute_sample_qc().write(get_sample_qc().path, overwrite=args.overwrite)

    if args.impute_sex:
        compute_sex().write(sex.path, overwrite=args.overwrite)
    elif args.reannotate_sex:
        sex_ht = sex.ht().checkpoint('gs://gnomad-tmp/sex_ht_checkpoint.ht', overwrite=True)  # Copy HT to temp location to overwrite annotation
        x_ploidy_cutoff, y_ploidy_cutoff = get_ploidy_cutoffs(sex_ht, f_stat_cutoff=0.5)
        sex_ht = sex_ht.annotate(
            **get_sex_expr(
                sex_ht.chrX_ploidy,
                sex_ht.chrY_ploidy,
                x_ploidy_cutoff,
                y_ploidy_cutoff
            )
        )
        sex_ht.write(sex.path, overwrite=args.overwrite)

    if args.compute_hard_filters:
        compute_hard_filters(
            args.min_cov
        ).write(hard_filtered_samples.path, overwrite=args.overwrite)

    if args.compute_qc_mt:
        compute_qc_mt().write(qc.path, overwrite=args.overwrite)

    if args.run_pc_relate:
        logger.info('Running PC-Relate')
        logger.warn("PC-relate requires SSDs and doesn't work with preemptible workers!")
        qc_mt = qc.mt()
        eig, scores, _ = hl.hwe_normalized_pca(qc_mt.GT, k=10, compute_loadings=False)
        scores = scores.checkpoint(pc_relate_pca_scores.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        relatedness_ht = hl.pc_relate(qc_mt.GT, min_individual_maf=0.01, scores_expr=scores[qc_mt.col_key].scores,
                                      block_size=4096, min_kinship=0.05, statistics='all')
        relatedness_ht.write(relatedness.path, args.overwrite)

    if args.run_pca:
        rank_ht = compute_sample_rankings(use_qc_metrics_filters=False)  # QC metrics filters do not exist at this point
        rank_ht = rank_ht.checkpoint(pca_samples_rankings.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        filtered_samples = hl.literal(rank_ht.aggregate(hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))))  # TODO: don't localize once hail bug is fixed
        samples_to_drop = compute_related_samples_to_drop(
            relatedness.ht(),
            rank_ht,
            args.kin_threshold,
            filtered_samples=filtered_samples
        )
        samples_to_drop.checkpoint(pca_related_samples_to_drop.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        pop_pca_eigenvalues, pop_pca_scores_ht, pop_pca_loadings_ht = run_pca(args.include_unreleasable_samples, args.n_pcs, samples_to_drop)
        pop_pca_scores_ht.write(ancestry_pca_scores(args.include_unreleasable_samples).path, overwrite=args.overwrite)
        pop_pca_loadings_ht.write(ancestry_pca_loadings(args.include_unreleasable_samples).path, overwrite=args.overwrite)
        with hl.utils.hadoop_open(ancestry_pca_eigenvalues_path(args.include_unreleasable_samples), mode='w') as f:
            f.write(",".join([str(x) for x in pop_pca_eigenvalues]))

    if args.assign_pops:
        pop_ht, pops_rf_model = assign_pops(args.min_pop_prob, args.include_unreleasable_samples)
        pop_ht = pop_ht.checkpoint(pop.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        pop_ht.transmute(
            **{f'PC{i + 1}': pop_ht.pca_scores[i] for i in range(0, 10)}
        ).export(pop_tsv_path)

        with hl.hadoop_open(pop_rf_path, 'wb') as out:
            pickle.dump(pops_rf_model, out)

    if args.calculate_inbreeding:
        qc_mt = qc.mt()
        pop_ht = pop.ht()
        qc_mt = qc_mt.annotate_cols(pop=pop_ht[qc_mt.col_key].pop)
        qc_mt = qc_mt.annotate_rows(call_stats_by_pop=hl.agg.group_by(qc_mt.pop, hl.agg.call_stats(qc_mt.GT)))
        inbreeding_ht = qc_mt.annotate_cols(
            inbreeding=hl.agg.inbreeding(qc_mt.GT, qc_mt.call_stats_by_pop[qc_mt.pop].AF[1])
        ).cols().select('inbreeding')
        inbreeding_ht.write(sample_inbreeding.path, overwrite=args.overwrite)

    if args.apply_stratified_filters:
        apply_stratified_filters(
            args.filtering_qc_metrics.split(",")
        ).write(stratified_metrics.path, overwrite=args.overwrite)

    if args.apply_regressed_filters:
        apply_regressed_filters(
            args.filtering_qc_metrics.split(","),
            args.include_unreleasable_samples
        ).write(regressed_metrics.path, overwrite=args.overwrite)

    if args.compute_related_samples_to_drop:
        rank_ht = compute_sample_rankings(use_qc_metrics_filters=True)
        rank_ht = rank_ht.checkpoint(release_samples_rankings.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        filtered_samples = hl.literal(rank_ht.aggregate(hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))))  # TODO: don't localize once hail bug is fixed
        print(filtered_samples)
        samples_to_drop = compute_related_samples_to_drop(
            relatedness.ht(),
            rank_ht,
            args.kin_threshold,
            filtered_samples=filtered_samples
        )
        samples_to_drop.write(release_related_samples_to_drop.path, overwrite=args.overwrite)

    if args.generate_metadata:
        meta_ht = generate_metadata()
        meta_ht.checkpoint(meta.path, overwrite=args.overwrite, _read_if_exists=not args.overwrite)
        n_pcs = meta_ht.aggregate(hl.agg.min(hl.len(meta_ht.pca_scores)))
        meta_ht = meta_ht.transmute(
            **{f'PC{i + 1}': meta_ht.pca_scores[i] for i in range(n_pcs)},
            hard_filters=hl.or_missing(hl.len(meta_ht.hard_filters) > 0, hl.delimit(meta_ht.hard_filters)),
            qc_metrics_filters=hl.or_missing(hl.len(meta_ht.qc_metrics_filters) > 0, hl.delimit(meta_ht.qc_metrics_filters))
        )
        meta_ht.flatten().export(meta_tsv_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--overwrite', help='Overwrite all data from this subset (default: False)', action='store_true')
    parser.add_argument('--sample_qc', help='Assigns pops from PCA', action='store_true')
    parser.add_argument('--impute_sex', help='Runs sex imputation. Also runs sex karyotyping annotation.', action='store_true')
    parser.add_argument('--reannotate_sex', help='Runs the sex karyotyping annotations again, without re-computing sex imputation metrics.', action='store_true')
    parser.add_argument('--compute_hard_filters', help='Computes samples to be hard-filtered', action='store_true')
    parser.add_argument('--min_cov', help="Minimum coverage for inclusion when computing har-filters", default=18, type=int)
    parser.add_argument('--compute_samples_ranking', help='Computes global samples ranking based on hard-filters, releasable and coverage.', action='store_true')
    parser.add_argument('--compute_qc_mt', help='Creates the QC MT based on liftover of v2 QC and Purcell 5k sites', action='store_true')
    parser.add_argument('--run_pc_relate', help='Run PC-relate', action='store_true')
    parser.add_argument('--compute_related_samples_to_drop', help='Flags related samples to drop', action='store_true')
    parser.add_argument('--kin_threshold', help='Maximum kin threshold to be considered unrelated', default=0.1, type=float)
    parser.add_argument('--min_related_hard_filter', help='Minimum number of relateds to have to get hard-filterd', default=50, type=int)
    parser.add_argument('--n_pcs', help='Number of PCs to compute for ancestry PCA', default=30, type=int)
    parser.add_argument('--run_pca', help='Compute PCA', action='store_true')
    parser.add_argument('--include_unreleasable_samples', help='Includes unreleasable samples for computing PCA', action='store_true')
    parser.add_argument('--assign_pops', help='Assigns pops from PCA', action='store_true')
    parser.add_argument('--min_pop_prob', help='Minimum RF prob for pop assignment', default=0.9, type=float)
    parser.add_argument('--calculate_inbreeding', help='Calculate sample level inbreeding', action='store_true')
    parser.add_argument('--filtering_qc_metrics', help="List of QC metrics for filtering.", default=",".join([
        'n_snp', 'n_singleton', 'r_ti_tv', 'r_insertion_deletion', 'n_insertion', 'n_deletion', 'r_het_hom_var',
        'n_het', 'n_hom_var', 'n_transition', 'n_transversion'
    ]))
    parser.add_argument('--apply_stratified_filters', help="Compute per pop filtering.", action='store_true')
    parser.add_argument('--apply_regressed_filters', help='Computes qc_metrics adjusted for pop.', action='store_true')
    parser.add_argument('--generate_metadata', help='Generates the metadata HT.', action='store_true')

    main(parser.parse_args())
