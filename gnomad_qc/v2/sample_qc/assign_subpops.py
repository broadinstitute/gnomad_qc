from gnomad.utils.slack import slack_notifications
from gnomad.sample_qc.ancestry import pc_project
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.sample_qc.joint_sample_qc import split_mt_by_relatedness, run_assign_population_pcs
import argparse
import logging
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subpops")
logger.setLevel(logging.INFO)

# 3 letter country codes are ISO ALPHA-3 (https://www.nationsonline.org/oneworld/country_code_list.htm)


def get_known_populations(ht: hl.Table, pop: str):  # TODO: bring data into separate file and load here
    if pop == 'eur' or pop == 'nfe':
        known_pops = hl.literal({
            'ICR1000': 'nwe',  # gb, ICR
            'ICR142': 'nwe',  # gb, ICR
            'C1017': 'seu',  # it, ATVB
            'C1568': 'seu',  # es, Regicor
            'Bulgarian_Trios': 'bgr',  # Bulgarians
            'C533': 'bgr',  # Bulgarians
            'C821': 'bgr',  # Bulgarians
            'C952': 'bgr',  # Bulgarians
            'G89634': 'est',  # Estonians
            'G94980': 'est',  # Estonians
            # 'C1830': 'neu',  # gb, Leicester
            'C1972': 'nwe',  # gb, Tayside region of Scotland
            # 'G94051': 'seu',  # es, "United States and Spain"
            # 'C1708': 'deu',  # Germans?
            'C1508': 'swe',  # Swedes
            'C1509': 'swe',  # Swedes
        })
    elif pop == 'eas':
        known_pops = hl.literal({
            'C1397': 'oea',  # 'twn',  # Taiwanese trios
            'C1443': 'oea',  # 'twn',  # Taiwanese trios
            'C1506': 'oea',  # 'twn',  # Taiwanese trios
            'C1867': 'oea',  # 'twn',  # Taiwanese trios
            'C978': 'oea',  # 'twn',  # Taiwanese trios
            'C774': 'kor',  # Korean T2D project
            'C1982': 'kor',  # Korean
            'C1940': 'oea',  # 'sgp',  # Singapore
            'C1980': 'oea',  # 'hkg',  # Hong Kong
            '1kg_JPT': 'jpn'
        })
    elif pop == 'afr':
        known_pops = hl.literal({
            'C773': 't2d',  # African American T2D
            'C1002': 't2d',  # African American T2D
            'C1567': 'jhs',  # African American JHS
            'C1956': 'biome',  # African American BioMe
            # TODO: Add 1kg populations here
        })
    else:
        raise ValueError('pop must be one of eur, nfe, eas, afr')
    ht = ht.annotate(known_pop=known_pops.get(ht.meta.project_id))
    if pop == 'eur':
        finns = hl.import_table('gs://gnomad/sample_qc/input_meta/source/99percent_finns_plus_AD_IBD_NFID.tsv.bgz', impute=True)
        finns = finns.filter(finns.percent_finnish > 0.99).key_by('sample_name_in_vcf')
        ht = ht.annotate(known_pop=hl.cond(hl.is_defined(finns[ht.s]), 'fin', ht.known_pop))
    return ht


def main(args):
    hl.init(log='/subpops.log')

    if args.population == 'all':
        pcs = list(range(1, 7))
    elif args.population == 'eur':
        pcs = [1, 2, 3]
    elif args.population == 'eas':
        pcs = [1, 2]
    else:
        pcs = [1, 2]

    if not args.skip_filtering:
        pruned_mt = hl.read_matrix_table(qc_mt_path('joint', ld_pruned=True))
        exome_project_table = hl.read_table(qc_ht_path('exomes', 'hard_filters')).select('project_id')
        exome_platform_table = hl.read_table(qc_ht_path('exomes', 'platforms')).select('qc_platform')
        exome_table = exome_project_table.annotate(qc_platform=hl.str(exome_platform_table[exome_project_table.key].qc_platform))
        genome_table = hl.read_table(qc_ht_path('genomes', 'hard_filters')).select('project_id', 'qc_platform')
        joint_table = exome_table.union(genome_table)
        exome_pop_table = hl.read_table(qc_ht_path('exomes', 'pop_platform')).select('pop')
        genome_pop_table = hl.read_table(qc_ht_path('genomes', 'pop_platform')).select('pop')
        pop_table = exome_pop_table.union(genome_pop_table)
        pop_table = pop_table.annotate(project_id=joint_table[pop_table.key].project_id, qc_platform=joint_table[pop_table.key].qc_platform)
        pruned_mt = pruned_mt.annotate_cols(meta=pop_table[pruned_mt.col_key])
        variants, samples = pruned_mt.count()
        logger.info(f'{samples} samples, {variants} variants found in original joint MT')

        if args.population == 'all':
            sample_criteria = True
        elif args.population == 'eur':
            sample_criteria = (pruned_mt.meta.pop == "nfe") | (pruned_mt.meta.pop == "fin")
        elif args.population == 'eas':
            sample_criteria = (pruned_mt.meta.pop == "eas") & (pruned_mt.data_type == "exomes")
        else:
            sample_criteria = pruned_mt.meta.pop == args.population

        pruned_mt = pruned_mt.filter_cols(sample_criteria)
        variants, samples = pruned_mt.count()
        logger.info(f'{samples} samples, {variants} variants found in {args.population} in joint MT')

        pca_mt, related_mt = split_mt_by_relatedness(pruned_mt)

        # Filter variants by callrate on each platform
        pca_platforms_mt = pca_mt.group_cols_by(pca_mt.meta.qc_platform).aggregate(
            missing=hl.agg.count_where(hl.is_missing(pca_mt.GT)),
            total=hl.agg.count()
        )
        # All variants must have a callrate at least .999 in each platform, or no more than 1 missing sample if platform <= 1000 samples
        pca_platforms_mt = pca_platforms_mt.annotate_entries(
            remove_variant=(hl.case()
                            .when(pca_platforms_mt.total > 1000, pca_platforms_mt.missing / pca_platforms_mt.total > 0.001)
                            .default(pca_platforms_mt.missing > 1))
        )
        pca_platforms_mt = pca_platforms_mt.filter_rows(hl.agg.any(pca_platforms_mt.remove_variant), keep=False)
        pca_mt = pca_mt.filter_rows((hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2 > 0.001) & hl.is_defined(pca_platforms_mt.rows()[pca_mt.row_key]))
        variants, samples = pca_mt.count()
        logger.info(f'{samples} samples, {variants} variants found in {args.population} in PCA MT after filtering variants by AF and platform callrate')

        pca_pruned = hl.ld_prune(pca_mt.GT, r2=0.1)
        pca_mt = pca_mt.filter_rows(hl.is_defined(pca_pruned[pca_mt.row_key]))
        related_mt = related_mt.filter_rows(hl.is_defined(pca_mt.rows()[related_mt.row_key]))
        pca_mt.write(f"{qc_temp_data_prefix('joint')}.{args.population}.unrelated.filtered.mt", args.overwrite)
        related_mt.write(f"{qc_temp_data_prefix('joint')}.{args.population}.related.filtered.mt", args.overwrite)

    pca_mt = hl.read_matrix_table(f"{qc_temp_data_prefix('joint')}.{args.population}.unrelated.filtered.mt")
    related_mt = hl.read_matrix_table(f"{qc_temp_data_prefix('joint')}.{args.population}.related.filtered.mt")

    variants, samples = pca_mt.count()
    logger.info(f'{samples} samples after removing relateds, {variants} variants after filtering and LD pruning')

    if not args.skip_pop_pca:
        pca_evals, pca_scores, pca_loadings = hl.hwe_normalized_pca(pca_mt.GT, k=10, compute_loadings=True)
        pca_mt = pca_mt.annotate_rows(pca_af=hl.agg.mean(pca_mt.GT.n_alt_alleles()) / 2)
        pca_loadings = pca_loadings.annotate(pca_af=pca_mt.rows()[pca_loadings.key].pca_af)
        pca_scores.write(ancestry_pca_scores_ht_path(args.population), args.overwrite)
        pca_loadings.write(ancestry_pca_loadings_ht_path(args.population), args.overwrite)

    pca_scores = hl.read_table(ancestry_pca_scores_ht_path(args.population))
    pca_loadings = hl.read_table(ancestry_pca_loadings_ht_path(args.population))
    pca_mt = pca_mt.annotate_cols(scores=pca_scores[pca_mt.col_key].scores)

    variants, samples = related_mt.count()
    logger.info(f'Projecting population PCs for {samples} related samples...')
    related_ht = pc_project(related_mt, pca_loadings)
    related_mt = related_mt.annotate_cols(scores=related_ht[related_mt.col_key].scores)

    logger.info('Assigning population annotations...')
    pop_colnames = ['related', 'known_pop', 'scores']
    # Join MTs, then annotate with known pops, then select out columns we care about, then assign pop pcs
    joint_ht = pca_mt.cols().union(related_mt.cols())
    joint_ht = get_known_populations(joint_ht, args.population)
    joint_ht = joint_ht.select(*pop_colnames,
                               **{f"PC{i + 1}": joint_ht.scores[i] for i in range(10)})
    joint_pca_ht, joint_pca_fit = run_assign_population_pcs(joint_ht,
                                                            f'{qc_temp_data_prefix("joint")}.RF_pop_assignments.{args.population}.txt.bgz',
                                                            f'{qc_temp_data_prefix("joint")}.RF_fit.{args.population}.pkl',
                                                            pcs=pcs)
    joint_pca_ht = joint_pca_ht.annotate(pop=hl.cond(
        joint_pca_ht.pop == 'oth', hl.literal(f'o{args.population[:2]}'), joint_pca_ht.pop))
    joint_ht = joint_ht.select(**{f"subpop_{args.population}_PC{i}": joint_ht[f"PC{i}"] for i in range(1, 11)},
                               subpop=joint_pca_ht[joint_ht.key].pop,
                               known_subpop=joint_pca_ht[joint_ht.key].known_pop)
    joint_ht.write(subpop_ht_path(args.population), args.overwrite)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--population', help='Which super-populations to select (can be one of the major pops, eur, or all)', required=True)
    parser.add_argument('--skip_filtering', help='Skip calculating filtered joint MT', action='store_true')
    parser.add_argument('--skip_pop_pca', help='Skip calculating population PCs on unrelated samples', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
