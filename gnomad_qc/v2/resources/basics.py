import hail as hl
from typing import *
from gnomad.utils.file_utils import file_exists

CURRENT_HAIL_VERSION = "0.2"
CURRENT_RELEASE = "2.1.1"
CURRENT_GENOME_META = "2018-10-11"  # YYYY-MM-DD
CURRENT_EXOME_META = "2018-10-11"
CURRENT_FAM = '2018-04-12'
CURRENT_DUPS = '2017-10-04'

RELEASES = ["2.0.1", "2.0.2", "2.1", "2.1.1", "test"]

SUBPOPS = {'NFE': ['BGR', 'EST', 'NWE', 'SEU', 'SWE', 'ONF'],
           'EAS': ['KOR', 'JPN', 'OEA']
           }
GENOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH']
EXOME_POPS = ['AFR', 'AMR', 'ASJ', 'EAS', 'FIN', 'NFE', 'OTH', 'SAS']
EXAC_POPS = ["AFR", "AMR", "EAS", "FIN", "NFE", "OTH", "SAS"]

pop_names = {
    'oth': 'Other',
    'afr': 'African-American/African',
    'ami': 'Amish',
    'amr': 'Latino',
    'eas': 'East Asian',
    'fin': 'Finnish',
    'eur': 'European',
    'nfe': 'Non-Finnish European',
    'sas': 'South Asian',
    'mde': 'Middle Eastern',
    'asj': 'Ashkenazi Jewish',
    'uniform': 'Uniform',
    'sas_non_consang': 'South Asian (F < 0.05)',
    'consanguineous': 'South Asian (F > 0.05)',
    'exac': 'ExAC',
    'bgr': 'Bulgarian (Eastern European)',
    'deu': 'German',
    'est': 'Estonian',
    'esp': 'Spanish',
    'gbr': 'British',
    'nwe': 'North-Western European',
    'seu': 'Southern European',
    'ita': 'Italian',
    'swe': 'Swedish',
    'chn': 'Chinese',
    'kor': 'Korean',
    'hkg': 'Hong Kong',
    'sgp': 'Singaporean',
    'twn': 'Taiwanese',
    'jpn': 'Japanese',
    'oea': 'Other East Asian',
    'oeu': 'Other European',
    'onf': 'Other Non-Finnish European',
    'unk': 'Unknown'
}

REFERENCES = {'GRCh37', 'GRCh38'}

def get_release_file(file_path: str, version: str = CURRENT_RELEASE) -> str:
    """
    Tries to get the desired file from the corresponding release version on the google cloud.
    If the file is not found for the desired release version, falls back on previous versions.

    :param  str file_path: Desired file path, with {0} as placeholder(s) for the version number
    :param str version: Desired file version
    :return: Path for closest version of the file available
    """
    if file_exists(file_path.format(version)):
        return file_path.format(version)
    else:
        for v in range(RELEASES.index(version) - 1, -1, -1):
            if file_exists(file_path.format(RELEASES[v])):
                print("WARN: Resource {} could not be found for gnomAD release version {}.\n Loading gnomAD release version {} of the file. ({})".format(
                    file_path.format(version),
                    version,
                    RELEASES[v],
                    file_path.format(RELEASES[v])
                ))
                return file_path.format(RELEASES[v])

        print("ERROR: Resource {} could not be found for any release.".format(
            file_path.format(version)
        ))
        return file_path.format(version)


def public_exomes_ht_path(split=True, version=CURRENT_RELEASE):
    if int(version[0]) > 1 and int(version[2]) > 0:
        return get_release_file('gs://gcp-public-data--gnomad/release/{0}/ht/exomes/gnomad.exomes.r{0}.sites.ht')
    else:
        return get_release_file('gs://gnomad-public/release/{{0}}/vds/exomes/gnomad.exomes.r{{0}}.sites{0}.vds'.format(".split" if split else ""))


def public_genomes_ht_path(split=True, version=CURRENT_RELEASE):
    if int(version[0]) > 1 and int(version[2]) > 0:
        return get_release_file('gs://gcp-public-data--gnomad/release/{0}/ht/genomes/gnomad.genomes.r{0}.sites.ht')
    else:
        return get_release_file('gs://gnomad-public/release/{{0}}/vds/genomes/gnomad.genomes.r{{0}}.sites{0}.vds'.format(".split" if split else ""))


def get_gnomad_public_data(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get public gnomAD data as VDS.

    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Chosen VDS
    :rtype: MatrixTable
    """
    return hl.read_table(get_gnomad_public_data_path(data_type, split=split, version=version))


def get_gnomad_data(data_type: str, adj: bool = False, split: bool = True, raw: bool = False,
                    non_refs_only: bool = False, hail_version: str = CURRENT_HAIL_VERSION,
                    meta_version: str = None, meta_root: Optional[str] = 'meta', full_meta: bool = False,
                    fam_version: str = CURRENT_FAM, fam_root: str = None, duplicate_mapping_root: str = None,
                    release_samples: bool = False, release_annotations: bool = None) -> hl.MatrixTable:
    """
    Wrapper function to get gnomAD data as VDS. By default, returns split hardcalls (with adj annotated but not filtered)

    :param str data_type: One of `exomes` or `genomes`
    :param bool adj: Whether the returned data should be filtered to adj genotypes
    :param bool split: Whether the dataset should be split (only applies to raw=False)
    :param bool raw: Whether to return the raw (10T+) data (not recommended: unsplit, and no special consideration on sex chromosomes)
    :param bool non_refs_only: Whether to return the non-ref-genotype only MT (warning: no special consideration on sex chromosomes)
    :param str hail_version: One of the HAIL_VERSIONs
    :param str meta_version: Version of metadata (None for current)
    :param str meta_root: Where to put metadata. Set to None if no metadata is desired.
    :param str full_meta: Whether to add all metadata (warning: large)
    :param str fam_version: Version of metadata (default to current)
    :param str fam_root: Where to put the pedigree information. Set to None if no pedigree information is desired.
    :param str duplicate_mapping_root: Where to put the duplicate genome/exome samples ID mapping (default is None -- do not annotate)
    :param bool release_samples: When set, filters the data to release samples only
    :param str release_annotations: One of the RELEASES to add variant annotations (into va), or None for no data
    :return: gnomAD hardcalls dataset with chosen annotations
    :rtype: MatrixTable
    """
    from gnomad.utils.filtering import filter_to_adj

    if raw and split:
        raise DataException('No split raw data. Use of hardcalls is recommended.')

    if non_refs_only:
        mt = hl.read_matrix_table(get_gnomad_data_path(data_type, split=split, non_refs_only=non_refs_only, hail_version=hail_version))
    else:
        mt = hl.read_matrix_table(get_gnomad_data_path(data_type, hardcalls=not raw, split=split, hail_version=hail_version))

    if adj:
        mt = filter_to_adj(mt)

    if meta_root:
        meta_ht = get_gnomad_meta(data_type, meta_version, full_meta=full_meta)
        mt = mt.annotate_cols(**{meta_root: meta_ht[mt.s]})

    if duplicate_mapping_root:
        dup_ht = hl.import_table(genomes_exomes_duplicate_ids_tsv_path, impute=True,
                                 key='exome_id' if data_type == "exomes" else 'genome_id')
        mt = mt.annotate_cols(**{duplicate_mapping_root: dup_ht[mt.s]})

    if fam_root:
        fam_ht = hl.import_fam(fam_path(data_type, fam_version))
        mt = mt.annotate_cols(**{fam_root: fam_ht[mt.s]})

    if release_samples:
        mt = mt.filter_cols(mt.meta.release)

    if release_annotations:
        sites_ht = get_gnomad_public_data(data_type, split)
        mt = mt.select_rows(**sites_ht[mt.row_key])
        mt = mt.select_globals(**sites_ht.index_globals())

    return mt


def get_gnomad_meta(data_type: str, version: str = None, full_meta: bool = False) -> hl.Table:
    """
    Wrapper function to get gnomAD metadata as Table

    :param str data_type: One of `exomes`, `genomes` or `joint`
    :param str version: Metadata version (None for current)
    :param bool full_meta: Whether to annotate full metadata (rather than just summarized version)
    :return: Metadata Table
    :rtype: Table
    """
    if data_type == 'joint':
        exomes_ht = hl.read_table(get_gnomad_meta_path('exomes', version)).key_by('s', data_type='exomes')
        exomes_ht = exomes_ht.annotate(qc_platform=hl.str(exomes_ht.qc_platform))
        genomes_ht = hl.read_table(get_gnomad_meta_path('genomes', version)).key_by('s', data_type='genomes')
        ht = exomes_ht.union(genomes_ht, unify=True)
    else:
        ht = hl.read_table(get_gnomad_meta_path(data_type, version)).key_by('s')

    if not full_meta:
        columns = ['age', 'sex',
                   'hard_filters', 'perm_filters', 'pop_platform_filters', 'related',
                   'data_type', 'product', 'product_simplified', 'qc_platform',
                   'project_id', 'project_description', 'internal', 'investigator',
                   'known_pop', 'known_subpop', 'pop', 'subpop',
                   'neuro', 'control', 'topmed',
                   'high_quality', 'release']
        if data_type == 'genomes':
            columns.extend(['pcr_free', 'project_name', 'release_2_0_2'])
        elif data_type == 'exomes':
            columns.extend(['diabetes', 'exac_joint', 'tcga'])

        ht = ht.select(*columns)
    return ht


def get_gnomad_public_data_path(data_type, split=True, version=CURRENT_RELEASE):
    """
    Wrapper function to get paths to gnomAD data

    :param str data_type: One of `exomes` or `genomes`
    :param bool split: Whether the dataset should be split
    :param str version: One of the RELEASEs
    :return: Path to chosen VDS
    :rtype: str
    """
    if version not in RELEASES:
        return DataException("Select version as one of: {}".format(','.join(RELEASES)))

    if data_type == 'exomes':
        return public_exomes_ht_path(split, version)
    elif data_type == 'genomes':
        return public_genomes_ht_path(split, version)
    return DataException("Select data_type as one of 'genomes' or 'exomes'")


def get_gnomad_data_path(data_type, hardcalls=False, split=True, non_refs_only=False, hail_version=CURRENT_HAIL_VERSION):
    """
    Wrapper function to get paths to gnomAD data

    :param str data_type: One of `exomes` or `genomes`
    :param bool hardcalls: Whether hardcalls should be returned
    :param bool split: Whether the dataset should be split (applies to hardcalls and non_refs_only)
    :param bool non_refs_only: Whether non-ref-genotype only MT should be returned
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to chosen VDS
    :rtype: str
    """
    if hardcalls and non_refs_only:
        raise DataException('No dataset with hardcalls and non_refs_only')
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")
    if hardcalls:
        return hardcalls_mt_path(data_type, split, hail_version)
    elif non_refs_only:
        return non_refs_only_mt_path(data_type, split)
    else:
        return raw_exomes_mt_path(hail_version) if data_type == 'exomes' else raw_genomes_mt_path(hail_version)


def get_gnomad_meta_path(data_type, version=None):
    """
    Wrapper function to get paths to gnomAD metadata

    :param str data_type: One of `exomes` or `genomes`
    :param str version: String with version (date) for metadata
    :return: Path to chosen metadata file
    :rtype: str
    """
    if data_type == 'exomes':
        if version:
            return metadata_exomes_ht_path(version)
        return metadata_exomes_ht_path()
    elif data_type == 'genomes':
        if version:
            return metadata_genomes_ht_path(version)
        return metadata_genomes_ht_path()
    return DataException("Select data_type as one of 'genomes' or 'exomes'")


def raw_exomes_mt_path(hail_version=CURRENT_HAIL_VERSION):
    """
    Warning: unsplit and no special consideration on sex chromosomes
    """
    return 'gs://gnomad/raw/hail-{0}/mt/exomes/gnomad.exomes.mt'.format(hail_version)


def raw_genomes_mt_path(hail_version=CURRENT_HAIL_VERSION):
    """
    Warning: unsplit and no special consideration on sex chromosomes
    """
    return 'gs://gnomad/raw/hail-{0}/mt/genomes/gnomad.genomes.mt'.format(hail_version)


def raw_exac_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/raw/hail-{0}/mt/exac/exac.mt'.format(hail_version)


def exac_release_sites_ht_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/raw/hail-{}/ht/exac/exac.r1.sites.vep.ht'.format(hail_version)


def hardcalls_mt_path(data_type, split=True, hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad/hardcalls/hail-{0}/mt/{1}/gnomad.{1}{2}.mt'.format(hail_version, data_type,
                                                                           "" if split else ".unsplit")


def non_refs_only_mt_path(data_type, split=True):
    return f'gs://gnomad/non_refs_only/hail-0.2/mt/{data_type}/gnomad.{data_type}{"" if split else ".unsplit"}.mt'


def pbt_phased_trios_mt_path(data_type: str, split: bool = True, trio_matrix: bool = False, hail_version : str = CURRENT_HAIL_VERSION):
    return "gs://gnomad/hardcalls/hail-{0}/mt/{1}/gnomad.{1}.trios.pbt_phased{2}{3}.mt".format(hail_version, data_type,
                                                                           "" if split else ".unsplit",
                                                                           ".trio_matrix" if trio_matrix else '')

def annotations_ht_path(data_type, annotation_type, hail_version=CURRENT_HAIL_VERSION):
    """
    Get sites-level annotations

    :param str data_type: One of "exomes" or "genomes"
    :param str annotation_type: One of "vep", "qc_stats", "family_stats", "frequencies", "rf", "omes_concordance", "NA12878_concordance", "syndip_concordance", "omes_by_platform_concordance"
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to annotations Table
    :rtype: str
    """
    return 'gs://gnomad/annotations/hail-{0}/ht/{1}/gnomad.{1}.{2}.ht'.format(hail_version, data_type,
                                                                              annotation_type)


def sample_annotations_table_path(data_type, annotation_type, hail_version=CURRENT_HAIL_VERSION):
    """
    Get samples-level annotations

    :param str data_type: One of "exomes" or "genomes"
    :param str annotation_type: One of "family_stats", "downsampling", "omes_concordance", "NA12878_concordance", "syndip_concordance"
    :param str hail_version: One of the HAIL_VERSIONs
    :return: Path to annotations VDS
    :rtype: str
    """
    return 'gs://gnomad/annotations/hail-{0}/sample_tables/{1}/gnomad.{1}.{2}.ht'.format(hail_version, data_type,
                                                                                          annotation_type)

gnomad_pca_mt_path = "gs://gnomad-genomes/sampleqc/gnomad.pca.mt"


def gnomad_public_pca_loadings_ht_path(version: str = CURRENT_RELEASE, subpop: str = ""):
    """
    Returns the path for the public gnomAD VDS containing sites and loadings from the PCA

    :param str version: One of the RELEASEs
    :param str subpop: Can be empty ("") -> global, "eas" or "nfe"
    :return: path to gnomAD public PCA loadings HT
    :rtype: str
    """
    return get_release_file("gs://gnomad-public/release/{{0}}/pca/gnomad.r{{0}}.pca_loadings{0}.ht".format(f'.{subpop}' if subpop else ""),
                            version)


def metadata_genomes_tsv_path(version=CURRENT_GENOME_META):
    return 'gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.tsv.bgz'.format(version)


def metadata_exomes_tsv_path(version=CURRENT_EXOME_META):
    return 'gs://gnomad/metadata/exomes/gnomad.exomes.metadata.{0}.tsv.bgz'.format(version)


def metadata_genomes_ht_path(version=CURRENT_GENOME_META):
    return 'gs://gnomad/metadata/genomes/gnomad.genomes.metadata.{0}.ht'.format(version)


def metadata_exomes_ht_path(version=CURRENT_EXOME_META):
    return 'gs://gnomad/metadata/exomes/gnomad.exomes.metadata.{0}.ht'.format(version)


def coverage_mt_path(data_type, grouped: bool = False) -> str:
    """
    Returns the gnomAD coverage MT.
    If `grouped` is true, then each entry contains aggregated information by platform and sex.

    :param data_type: One of 'exomes' or 'genomes'
    :param grouped: Whether to load the grouped or full MT
    :return: Coverage MT
    :rtype: MatrixTable
    """
    return 'gs://gnomad/coverage/hail-0.2/coverage/{0}/mt/gnomad.{0}.coverage{1}.mt'.format(
        data_type,
        '.grouped' if grouped else ''
    )


def coverage_ht_path(data_type) -> str:
    return f'gs://gnomad-public/release/2.1/coverage/{data_type}/gnomad.{data_type}.r2.1.coverage.ht'


def fam_path(data_type: str, version: str = CURRENT_FAM, true_trios: bool = False) -> str:
    """
    Returns the path to gnomAD pedigree file.

    :param str data_type: One of 'exomes' or 'genomes'
    :param str version: Version of the fam file to get
    :param bool true_trios: If set, removes all families with more than one offspring
    :return: Path to fam file
    :rtype: str
    """
    if not true_trios:
        return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.{version}.fam"
    else:
        return f"gs://gnomad/metadata/{data_type}/gnomad.{data_type}.{version}.true_trios.fam"


def genomes_exomes_duplicate_ids_tsv_path(version: str = CURRENT_DUPS) -> str:
    return f"gs://gnomad/metadata/join/gnomad.genomes_exomes.{version}.duplicate_ids.tsv"


def omni_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_omni2.5.b37.mt'.format(hail_version)


def mills_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/Mills_and_1000G_gold_standard.indels.b37.mt'.format(hail_version)


def hapmap_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hapmap_3.3.b37.mt'.format(hail_version)


def kgp_high_conf_snvs_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000G_phase1.snps.high_confidence.b37.mt'.format(hail_version)


def kgp_phase3_genotypes_mt_path(split: bool = True, hail_version=CURRENT_HAIL_VERSION) -> str:
    """
    1000 Genomes Phase 3 with genotypes (b37)
    Imported from: gs://genomics-public-data/1000-genomes-phase-3/vcf-20150220/ALL.chr*.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf
    Samples populations from: gs://gnomad-public/truth-sets/hail-0.2/1000G.GRCh38.20130502.phase3.sequence.index

    :param bool split: Whether to load to split or non-split version
    :param str hail_version: Hail version
    :return: Path to 1000 Genomes MT
    :rtype: str
    """
    return 'gs://gnomad-public/truth-sets/hail-{0}/1000Genomes_phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes{1}.mt'.format(hail_version, '.split' if split else '')


def NA12878_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.mt'.format(hail_version)


def syndip_mt_path(hail_version=CURRENT_HAIL_VERSION):
    return 'gs://gnomad-public/truth-sets/hail-{0}/hybrid.m37m.mt'.format(hail_version)


def cpg_sites_ht_path():
    return 'gs://gnomad-public/resources/methylation/cpg.ht'


REFERENCE_DATA = {
    'GRCh37': {
        'vep_config': 'gs://hail-common/vep/vep/vep85-loftee-gcloud.json',
        'all_possible': 'gs://gnomad-public/papers/2019-flagship-lof/v1.0/context/Homo_sapiens_assembly19.fasta.snps_only.vep_20181129.ht',
        'methylation': 'gs://gnomad-public/resources/methylation/methylation.ht',
    },
    'GRCh38': {
        'vep_config': 'gs://hail-common/vep/vep/vep95-GRCh38-loftee-gcloud.json',
        'methylation': 'gs://gnomad-resources/methylation/hail-0.2/methylation_GRCh38.ht',
    }
}

def methylation_sites_ht_path(ref: str = 'GRCh37'):
    if ref not in REFERENCE_DATA:
        raise DataException("Select reference as one of: {}".format(','.join(REFERENCES)))
    return REFERENCE_DATA[ref]['methylation']


def context_ht_path(ref: str = 'GRCh37'):
    if ref not in ('GRCh37', ):
        raise DataException("Reference must be GRCh37")
    return REFERENCE_DATA[ref]['all_possible']


def vep_config_path(ref: str = 'GRCh37'):
    if ref not in REFERENCES:
        raise DataException("Select reference as one of: {}".format(','.join(REFERENCES)))
    return REFERENCE_DATA[ref]['vep_config']


vep_config = vep_config_path()  # For backwards-compatibility


dbsnp_vcf_path = "gs://gnomad-public/truth-sets/source/All_20180423.vcf.bgz"
dbsnp_ht_path = "gs://gnomad-public/truth-sets/source/All_20180423.ht"

NA12878_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.bed"
NA12878_high_conf_exome_regions_bed_path = "gs://gnomad-public/truth-sets/source/union13callableMQonlymerged_addcert_nouncert_excludesimplerep_excludesegdups_excludedecoy_excludeRepSeqSTRs_noCNVs_v2.18_2mindatasets_5minYesNoRatio.bed"
syndip_high_conf_regions_bed_path = "gs://gnomad-public/truth-sets/source/hybrid.m37m.bed"
clinvar_vcf_path = "gs://gnomad-resources/clinvar/source/clinvar_20181028.vcf.bgz"
clinvar_ht_path = "gs://gnomad-resources/clinvar/hail-0.2/clinvar_20181028.vep.ht"

# Useful intervals
lcr_intervals_path = "gs://gnomad-public/intervals/LCR.GRCh37_compliant.interval_list"  # "gs://gnomad-public/intervals/LCR.interval_list"
decoy_intervals_path = "gs://gnomad-public/intervals/mm-2-merged.GRCh37_compliant.bed"  # "gs://gnomad-public/intervals/mm-2-merged.bed.gz"
purcell5k_intervals_path = "gs://gnomad-public/intervals/purcell5k.interval_list"
segdup_intervals_path = "gs://gnomad-public/intervals/hg19_self_chain_split_both.bed"

# Exome intervals
exomes_high_conf_regions_intervals_path = "gs://gnomad-public/intervals/exomes_high_coverage.auto.interval_list"
exome_calling_intervals_path = 'gs://gnomad-public/intervals/exome_calling_regions.v1.interval_list'
evaluation_intervals_path = 'gs://gnomad-public/intervals/exome_evaluation_regions.v1.noheader.interval_list'
high_coverage_intervals_path = 'gs://gnomad-public/intervals/high_coverage.auto.interval_list'

# Genome intervals
genome_evaluation_intervals_path = "gs://gnomad-public/intervals/hg19-v0-wgs_evaluation_regions.v1.interval_list"  # from Broad GP
genome_evaluation_intervals_path_hg38 = "gs://gnomad-public/intervals/hg38-v0-wgs_evaluation_regions.hg38.interval_list"
# More can be found at gs://broad-references/hg19


# constraint_ht_path = 'gs://gnomad-public/release/2.1/ht/constraint/constraint.ht'
constraint_ht_path = 'gs://gnomad-public/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.lof_metrics.by_gene.ht'


# Sample QC files
def qc_ht_path(data_type: str):
    return 'gs://gnomad/sample_qc/ht/gnomad.{}.high_callrate_common_biallelic_snps.ht'.format(data_type)


def qc_temp_data_prefix(data_type: str):
    return 'gs://gnomad/sample_qc/temp/{0}/gnomad.{0}'.format(data_type)


def qc_meta_path(data_type: str):
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.streamlined_metadata.2018-10-10.txt.bgz'
    else:
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-10-11.txt.bgz'


# liftover data
def get_gnomad_liftover_data_path(data_type, version=CURRENT_RELEASE) -> str:
    """
    Function to get paths to liftover gnomAD Tables
    :param str data_type: One of `exomes` or `genomes`
    :param str version: One of the RELEASEs
    :return: Path to chosen Table
    :rtype: str
    """
    if version not in RELEASES:
        raise DataException("Select version as one of: {}".format(','.join(RELEASES)))

    if data_type != 'exomes' and data_type != 'genomes': 
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")

    return f'gs://gnomad-public/release/{version}/liftover_grch38/ht/{data_type}/gnomad.{data_type}.r{version}.sites.liftover_grch38.ht'


gnomad_sv_mt_path = 'gs://gnomad/sv/gnomAD-SV_v2_rev1.polished.mt'
gnomad_sv_vcf_path = 'gs://gnomad/sv/gnomAD-SV_v2_rev1.polished.vcf.gz'
gnomad_sv_public_sites_vcf_path = 'gs://gnomad-public/papers/2019-sv/gnomad_v2.1_sv.sites.vcf.gz'
gnomad_sv_public_hists_ht_path = 'gs://gnomad-public/papers/2019-sv/gnomad_sv_hists.ht'
gnomad_sv_release_samples_list_path = 'gs://gnomad/sv/gnomAD-SV_v2_rev1_releasable.samples.list'

class DataException(Exception):
    pass
