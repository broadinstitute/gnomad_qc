import hail as hl
from gnomad_hail import logger
from gnomad_hail.resources.resource_utils import (
    MatrixTableResource,
    TableResource,
    PedigreeResource,
    VersionedPedigreeResource
)
from gnomad_hail.resources.grch38 import na12878_giab
V3_RELEASE_VERSION = 'r3.0'

# MatrixTableResources
v3_qc = MatrixTableResource('gs://gnomad/sample_qc/mt/genomes_v3/gnomad_v3_qc_mt_v2_sites_dense.mt')

# TableResources
gnomad_v3_genotypes = MatrixTableResource("gs://gnomad/raw/hail-0.2/mt/genomes_v3/gnomad_genomes_v3.repartitioned.mt")
last_END_position = TableResource('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_last_END_positions.ht')
v3_pc_relate_pca_scores = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_qc_mt_v2_sites_pc_scores.ht')
v3_relatedness = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_qc_mt_v2_sites_relatedness.ht')
v3_sex = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_sex.ht')
freq = TableResource('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3.frequencies.ht')
qual_hist = TableResource('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3.qual_hists.ht')
vep = TableResource('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_vep.ht')
gnomad_v2_qc_sites = TableResource('gs://gnomad-public/resources/grch38/gnomad_v2_qc_sites_b38.ht')
pca_related_samples_to_drop = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_related_samples_to_drop_for_pca.ht')
release_related_samples_to_drop = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_related_release_samples_to_drop.ht')
pop = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_pop.ht')
hard_filtered_samples = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_hard_filtered_samples.ht')
stratified_metrics = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_stratified_metrics.ht')
regressed_metrics = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_regressed_metrics.ht')
project_meta = TableResource(
    import_func=hl.import_table,
    import_args={
        'path': 'gs://gnomad/metadata/genomes_v3/09-09-2019_v3_project_meta.txt',
        'impute': True,
        'key': 's',
        'min_partitions': 100
    }
)

# Text files
info_vcf_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_info.vcf.bgz'
qual_hists_json_path = 'gs://gnomad/release/3.0/json/gnomad.genomes.r3.json'
pop_tsv_path = 'gs://gnomad/sample_qc/temp/genomes_v3/gnomad_v3_RF_pop_assignments.txt.gz'
pop_rf_path = 'gs://gnomad/sample_qc/temp/genomes_v3/gnomad_v3_pop.RF_fit.pickle'
meta_tsv_path = 'gs://gnomad/metadata/genomes_v3/gnomad_v3_metadata_2019-09-27.tsv.gz'

# TODO: Remove if not used after all python files are in
# internal_ht_path = 'gs://gnomad/release/3.0/ht/gnomad.genomes.r3.0.nested.no_subsets.sites.ht'

meta = TableResource('gs://gnomad/metadata/genomes_v3/gnomad_v3_metadata_2019-09-27.ht')
pca_samples_rankings = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_pca_samples_ranking.ht')
release_samples_rankings = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_release_samples_ranking.ht')
duplicates = TableResource('gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_duplicates.ht')


pedigree = VersionedPedigreeResource(
    'final', # TODO: Make sure "final" is the best label once the family scripts are in
    {
        'raw': PedigreeResource('gs://gnomad/metadata/genomes_v3/gnomad_v3_raw.fam'),
        'final': PedigreeResource('gs://gnomad/metadata/genomes_v3/gnomad_v3.fam')
    }
)

trios = VersionedPedigreeResource( # TODO: Should this be merged with Pedigree into a single resource?
    'final', # TODO: Make sure "final" is the best label once the family scripts are in
    {
        'raw': PedigreeResource('gs://gnomad/metadata/genomes_v3/gnomad_v3_trios_raw.fam'),
        'final': PedigreeResource('gs://gnomad/metadata/genomes_v3/gnomad_v3_trios.fam')
    }
)

ped_mendel_errors = TableResource('gs://gnomad/sample_qc/temp/genomes_v3/gnomad_v3_ped_chr20_mendel_errors.ht')
qc_ac = TableResource('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_qc_ac.ht')
fam_stats = TableResource('gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_qc_fam_stats.ht')

# TODO: All this is now in gnomad_hail --> Adapt variant QC
# truth_samples={
#     'syndip': {
#         's': "CHMI_CHMI3_WGS2",
#         'truth_mt': syndip_mt_path,
#         'v3_mt': v3_syndip_mt_path,
#         'bed': syndip_bed_path
#     },
#     'na12878': {
#         's': "NWD714003",
#         'truth_mt': na12878_giab.path,
#         'v3_mt': v3_na12878_mt_path,
#         'bed': na12878_bed_path
#     }
# }


def _get_ancestry_pca_ht_path(part: str, include_unreleasable_samples: bool = False) -> str:
    return 'gs://gnomad/sample_qc/ht/genomes_v3/gnomad_v3_pca_{0}{1}.{2}'.format(
        part,
        '_with_unreleasable_samples' if include_unreleasable_samples else '',
        'txt' if part == 'eigenvalues' else 'ht'
    )


def get_ancestry_pca_loadings(include_unreleasable_samples: bool = False) -> TableResource:
    return TableResource(_get_ancestry_pca_ht_path('loadings', include_unreleasable_samples))


def get_ancestry_pca_scores(include_unreleasable_samples: bool = False) -> TableResource:
    return TableResource(_get_ancestry_pca_ht_path('scores', include_unreleasable_samples))


def get_ancestry_pca_eigenvalues_path(include_unreleasable_samples: bool = False) -> str:
    return _get_ancestry_pca_ht_path('scores', include_unreleasable_samples)


def get_gnomad_v3_mt(
        key_by_locus_and_alleles: bool = False,
        remove_hard_filtered_samples: bool = True,
        release_only: bool = False,
        samples_meta: bool = False
) -> hl.MatrixTable:
    mt = gnomad_v3_genotypes.mt()
    if key_by_locus_and_alleles:
        mt = hl.MatrixTable(hl.ir.MatrixKeyRowsBy(mt._mir, ['locus', 'alleles'], is_sorted=True))

    if remove_hard_filtered_samples:
        mt = mt.filter_cols(hl.is_missing(hard_filtered_samples.ht()[mt.col_key]))

    if samples_meta:
        mt = mt.annotate_cols(meta=meta.ht()[mt.col_key])

        if release_only:
            mt = mt.filter_cols(mt.meta.release)

    elif release_only:
        mt = mt.filter_cols(meta.ht()[mt.col_key].release)

    return mt


def get_sample_qc(strat: str = "") -> TableResource:
    if strat:
        strat = f"_{strat}"
    return TableResource(f'gs://gnomad/sample_qc/ht/genomes_v3/sample_qc{strat}.ht')


def get_info(split: bool = True) -> TableResource:
    path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_info{}.ht'.format(
        '.split' if split else ''
    )
    return TableResource(path)


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

