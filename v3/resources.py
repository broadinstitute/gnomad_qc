import hail as hl

V3_RELEASE_VERSION = 'r3.0'

v3_qc_mt_path = 'gs://gnomad/sample_qc/mt/genomes_v3/gnomad_v3_qc_mt_v2_sites_dense.mt'
v3_pc_relate_pca_scores_ht_path = 'gs://gnomad/sample_qc/mt/genomes_v3/gnomad_v3_qc_mt_v2_sites_pc_scores.ht' # TODO -> move to /ht
v3_relatedness_ht_path = 'gs://gnomad/sample_qc/mt/genomes_v3/gnomad_v3_qc_mt_v2_sites_relatedness.ht' # TODO -> move to /ht
v3_sex_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_sex.ht'

telomeres_and_centromeres_bed_path = 'gs://gnomad-public/resources/grch38/hg38.telomeresAndMergedCentromeres.bed'
info_vcf_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_info.vcf.bgz'
freq_ht_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3.frequencies.ht'
qual_hists_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3.qual_hists.ht'
vep_ht_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_vep.ht'
gnomad_v2_qc_sites_path = 'gs://gnomad-public/resources/grch38/gnomad_v2_qc_sites_b38.ht'
last_END_position_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_last_END_positions.ht'
internal_ht_path = 'gs://gnomad/release/3.0/ht/gnomad.genomes.r3.0.nested.no_subsets.sites.ht'
qual_hists_json_path = 'gs://gnomad/release/3.0/json/gnomad.genomes.r3.json'
coverage_ht_path = 'gs://gnomad-public/release/3.0/coverage/genomes/gnomad.genomes.r3.0.coverage.ht'
coverage_tsv_path = 'gs://gnomad-public/release/3.0/coverage/genomes/gnomad.genomes.r3.0.coverage.summary.tsv.bgz'

pca_related_samples_to_drop_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_related_samples_to_drop_for_pca.ht'
release_related_samples_to_drop_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_related_release_samples_to_drop.ht'
pop_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_pop.ht'
pop_tsv_path = 'gs://gnomad/sample_qc/temp/genomes_v3/gnomad_v3_RF_pop_assignments.txt.gz'
pop_rf_path = 'gs://gnomad/sample_qc/temp/genomes_v3/gnomad_v3_pop.RF_fit.pickle'
hard_filtered_samples_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_hard_filtered_samples.ht'
project_meta_path= 'gs://gnomad/metadata/genomes_v3/09-09-2019_v3_project_meta.txt'
stratified_metrics_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_stratified_metrics.ht'
regressed_metrics_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_regressed_metrics.ht'
meta_ht_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3_metadata_2019-09-27.ht'
meta_tsv_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3_metadata_2019-09-27.tsv.gz'
pca_samples_rankings_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_pca_samples_ranking.ht'
release_samples_rankings_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_release_samples_ranking.ht'
duplicates_ht_path = 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_duplicates.ht'
ped_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3.ped'
raw_ped_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3_raw.ped'
raw_trios_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3_trios_raw.ped'
trios_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3_trios.ped'
ped_mendel_errors_path = 'gs://gnomad/sample_qc/temp/genomes_v3/gnomad_v3_ped_chr20_mendel_errors.ht'

ac_ht_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_qc_ac.ht'
fam_stats_ht_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_qc_fam_stats.ht'
clinvar_ht_path = 'gs://gnomad-public/resources/grch38/clinvar_20190923.ht'
truth_ht_path = 'gs://gnomad-public/resources/grch38/truth_sites.ht'
dbsnp_b38_vcf_path = 'gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.vcf.bgz'
dbsnp_b38_ht_path = 'gs://gnomad-public/resources/grch38/dbsnp_b151_grch38_all_20180418.ht'
syndip_bed_path = 'gs://gnomad-public/truth-sets/hail-0.2/syndip.b38.bed'
syndip_mt_path = 'gs://gnomad-public/truth-sets/hail-0.2/syndip.b38.mt'
v3_syndip_mt_path = 'gs://gnomad-public/truth-sets/hail-0.2/gnomad_v3_syndip.b38.mt'
na12878_mt_path = 'gs://gnomad-public/truth-sets/hail-0.2/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_PGandRTGphasetransfer.mt'
v3_na12878_mt_path = 'gs://gnomad-public/truth-sets/hail-0.2/gnomad_v3_na12878.mt'
na12878_bed_path = 'gs://gnomad-public/truth-sets/source/HG001_GRCh38_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-10X-SOLID_CHROM1-X_v.3.3.2_highconf_nosomaticdel_noCENorHET7.bed'

synthetic_grch38_vep_ht_path = 'gs://gnomad-resources/context/GRCh38/grch38_context_vep.ht'

truth_samples={
    'syndip': {
        's': "CHMI_CHMI3_WGS2",
        'truth_mt': syndip_mt_path,
        'v3_mt': v3_syndip_mt_path,
        'bed': syndip_bed_path
    },
    'na12878': {
        's': "NWD714003",
        'truth_mt': na12878_mt_path,
        'v3_mt': v3_na12878_mt_path,
        'bed': na12878_bed_path
    }
}


def get_qc_mt(remove_hard_filtered_samples: bool = True):
    qc_mt = hl.read_matrix_table(v3_qc_mt_path)
    if remove_hard_filtered_samples:
        hard_filtered_ht = hl.read_table(hard_filtered_samples_ht_path)
        qc_mt = qc_mt.filter_cols(hl.is_missing(hard_filtered_ht[qc_mt.col_key]))
    return qc_mt


def get_ancestry_pca_ht_path(part: str, include_unreleasable_samples: bool = False):
    return 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_pca_{0}{1}.{2}'.format(
        part,
        '_with_unreleasable_samples' if include_unreleasable_samples else '',
        'txt' if part == 'eigenvalues' else 'ht'
    )


def get_project_meta() -> hl.Table:
    return hl.import_table(project_meta_path, impute=True, key='s')


def get_full_mt_path(split: bool = True) -> str:
    split_str = '.split' if split else ''
    return f"gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned{split_str}.mt"


def get_full_mt(
        split: bool = True,
        key_by_locus_and_alleles: bool = False,
        remove_hard_filtered_samples: bool = True,
        release_only: bool = False
) -> hl.MatrixTable:
    mt = hl.read_matrix_table(get_full_mt_path(split))
    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(hl.ir.MatrixKeyRowsBy(mt._mir, ['locus', 'alleles'], is_sorted=True))

    if remove_hard_filtered_samples:
        hard_filtered_ht = hl.read_table(hard_filtered_samples_ht_path)
        mt = mt.filter_cols(hl.is_missing(hard_filtered_ht[mt.col_key]))

    return mt


def sample_qc_path(strat: str = ""):
    if strat:
        strat = f"_{strat}"
    return f'gs://gnomad/sample_qc/ht/genomes_v3/sample_qc{strat}.ht'


def get_lcr_intervals() -> hl.Table:
    return hl.import_locus_intervals('gs://gnomad-public/resources/grch38/LCRFromHengHg38.txt', reference_genome='GRCh38', skip_invalid_intervals=True)


def get_info_ht_path(split: bool = True):
    return 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_info{}.ht'.format(
        '.split' if split else ''
    )

def get_transmitted_singleton_vcf_path(confidence: str):
    return f'gs://gnomad/variant_qc/genomes_v3/transmitted_singletons_{confidence}.vcf'


def filters_ht_path(model_id: str, split: bool = True, finalized: bool = True): # TODO: Make this generic and take model_id as parameter
    return 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/filtering/{0}{1}{2}.ht'.format(
        model_id,
        '.finalized' if finalized else '',
        '.split' if split else ''
    )


def score_ranking_path(model_id: str, binned: bool):
    return 'gs://gnomad/variant_qc/genomes_v3/{0}.{1}.ht'.format(
        model_id,
        'binned' if binned else 'rank'
    )


def binned_concordance_path(model_id: str, truth_sample: str):
    return f'gs://gnomad/variant_qc/genomes_v3/{truth_sample}_{model_id}_binned_concordance.ht'


def release_ht_path(data_type: str = 'genomes', release_tag: str = V3_RELEASE_VERSION, public=True):
    '''
    Fetch filepath for release (variant-only) Hail Tables
    :param str data_type: 'exomes' or 'genomes'
    :param str release_tag: String describing release version for use in output filepaths
    :param bool nested: If True, fetch Table in which variant annotations (e.g., freq, popmax, faf, and age histograms)
        are in array format ("nested"); if False, fetch Table in which nested variant annotations are unfurled
    :param bool with_subsets: If True, fetch Table in which all gnomAD subsets are present; if False, fetch Table containing
        only gnomAD annotations
    :param bool temp: If True, fetch Table in which nested variant annotations are unfurled but listed under 'info' rather
        than at the top level; used for sanity-checking sites
    :return: Filepath for desired Hail Table
    :rtype: str
    '''
    release = release_tag.lstrip('r')
    if public:
        return f'gs://gnomad-public/release/{release}/ht/{data_type}/gnomad.{data_type}.{release_tag}.sites.ht'
    else:
        f'gs://gnomad/release/{release}/ht/gnomad.{data_type}.{release_tag}.sites.ht'

