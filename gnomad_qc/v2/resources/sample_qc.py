from .basics import DataException, get_gnomad_meta
from gnomad.utils.liftover import get_liftover_genome
import hail as hl


def get_liftover_v2_qc_mt(
        data_type: str,
        ld_pruned: bool,
        release_only: bool = False,
        overwrite: bool = False
) -> hl.MatrixTable:
    """
    Returns MatrixTable for sample QC purposes on build 38: can be exomes, genomes, or joint (joint dataset can also be ld_pruned=True)
    Criteria: callrate > 0.99, AF > 0.001, SNPs only, bi-allelics only
    Note: sites where the locus changes chromosome are discarded
    """
    path = qc_mt_path(data_type, ld_pruned, 'GRCh38')
    if not overwrite and hl.hadoop_exists(path):
        grch38_qc_mt = hl.read_matrix_table(path)
    else:
        grch38_qc_mt = hl.read_matrix_table(qc_mt_path(data_type, ld_pruned=ld_pruned))
        get_liftover_genome(grch38_qc_mt)
        grch38_qc_mt = grch38_qc_mt.key_rows_by()
        grch38_qc_mt = grch38_qc_mt.transmute_rows(
            locus=hl.liftover(grch38_qc_mt.locus, 'GRCh38'),
            locus37=grch38_qc_mt.locus
        )
        grch38_qc_mt = grch38_qc_mt.filter_rows(
            grch38_qc_mt.locus.contig == 'chr' + grch38_qc_mt.locus37.contig
        )
        grch38_qc_mt = grch38_qc_mt.key_rows_by(locus=grch38_qc_mt.locus, alleles=grch38_qc_mt.alleles)
        grch38_qc_mt = grch38_qc_mt.checkpoint(path, overwrite=overwrite)

    if release_only:
        meta = get_gnomad_meta(data_type)
        grch38_qc_mt = grch38_qc_mt.filter_cols(meta[grch38_qc_mt.col_key].release)

    return grch38_qc_mt


def qc_mt_path(data_type: str, ld_pruned: bool = False, reference_genome: str = 'GRCh37') -> str:
    """
    Returns MatrixTable path for sample QC purposes: can be exomes, genomes, or joint (joint dataset can also be ld_pruned=True)
    Criteria: callrate > 0.99, AF > 0.001, bi-allelics SNPs only
    """
    if data_type not in ('exomes', 'genomes', 'joint'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes' or 'joint'")
    if ld_pruned and data_type != 'joint':
        raise DataException("ld_pruned = True is only available for 'joint'")
    ref_str = ''
    if reference_genome == 'GRCh38':
        ref_str='.grch38'
    elif reference_genome != 'GRCh37':
        raise DataException('reference_genome must be one of "GRCh37" or "GRCh38"')

    ld_pruned = '.pruned' if ld_pruned else ''
    return f'gs://gnomad/sample_qc/mt/gnomad.{data_type}.high_callrate_common_biallelic_snps{ld_pruned}{ref_str}.mt'


def qc_ht_path(data_type: str, part: str) -> str:
    """
    Interim sample metadata tables generated in the sample qc process.
    Generally not to be used: use tables from basics.py instead (e.g. metadata_*omes_ht_path)
    hard_filters (contains hard filter, permission information, and sex)
    platforms (contains imputed platform information for exomes only)
    pop_platform (contains final related information and population/platform-specific QC filters)
    """
    if data_type not in ('exomes', 'genomes'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")
    if part not in ('hard_filters', 'platforms', 'pop_platform'):
        raise DataException("Select part as one of 'hard_filters', 'platforms', or 'pop_platform'")
    if data_type == 'genomes' and part == 'platforms':
        raise DataException("'platforms' only available for 'genomes'")
    return f'gs://gnomad/sample_qc/ht/gnomad.{data_type}.{part}.ht'


def rank_annotations_path(data_type: str) -> str:
    """
    Path to annotation data for ranking samples for related pruning. The 'joint' dataset is the results of the ranking.
    """
    if data_type not in ('exomes', 'genomes', 'joint'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes' or 'joint'")
    return f'gs://gnomad/sample_qc/tsv/gnomad.{data_type}.rank_list_annotations.txt.bgz'


def qc_temp_data_prefix(data_type: str) -> str:
    """
    Path to directory with intermediate files for sample QC
    """
    if data_type not in ('exomes', 'genomes', 'joint'):
        raise DataException("Select data_type as one of 'genomes' or 'exomes' or 'joint'")
    return f'gs://gnomad/sample_qc/temp/{data_type}/gnomad.{data_type}'


def qc_meta_path(data_type: str) -> str:
    """
    Input metadata file for sample QC
    """
    if data_type == 'exomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.exomes.streamlined_metadata.2018-10-10.txt.bgz'
    elif data_type == 'genomes':
        return 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.streamlined_metadata.2018-10-10.txt.bgz'
    else:
        raise DataException("Select data_type as one of 'genomes' or 'exomes'")


exome_callrate_scores_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.exomes.callrate_pca_scores.ht'
exome_callrate_mt_path = 'gs://gnomad/sample_qc/mt/gnomad.exomes.callrate.mt'

relatedness_ht_path = 'gs://gnomad/sample_qc/ht/gnomad.joint.relatedness.ht'


def ancestry_pca_scores_ht_path(population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'gs://gnomad/sample_qc/ht/gnomad.joint.unrelated.pca_scores{pop}.ht'


def ancestry_pca_loadings_ht_path (population: str = None) -> str:
    pop = f'.{population}' if population else ''
    return f'gs://gnomad/sample_qc/ht/gnomad.joint.unrelated.pca_loadings{pop}.ht'


def subpop_ht_path(population: str) -> str:
    return f'gs://gnomad/sample_qc/ht/gnomad.joint.subpop_assignments.{population}.ht'


known_population_annotations = 'gs://gnomad/sample_qc/input_meta/gnomad.pop_annots.txt'
estonian_batches = 'gs://gnomad/sample_qc/input_meta/gnomad.genomes.estonian_samples.txt'


# Resources for pedigree inference


def dup_pedigree_tsv_path(data_type: str):
    return f'gs://gnomad/sample_qc/fam/gnomad_{data_type}_dup_pedigree.tsv.bgz'


def raw_fam_path(data_type: str) -> str:
    """

    Returns the path for the raw pedigree path.
    By default each trio in the pedigree is unique, with samples picked amongst duplicates in order to maximize the number of samples within a family belonging to the same project.

    :param str data_type: One of 'exomes' or 'genomes'
    :return: Path to pedigree file
    :rtype: str
    """
    return f'gs://gnomad/sample_qc/fam/gnomad_{data_type}_raw.fam'


def fake_fam_path(data_type: str) -> str:
    """

    As part of the pedigree generation, a set of fake trios are generated by sampling 3 random samples.
    This return a path to the pedigree file containing those fake trios.

    :param str data_type: One of 'exomes' or 'genomes'
    :return: Path to pedigree file
    :rtype: str
    """
    return f'gs://gnomad/sample_qc/fam/gnomad_{data_type}_fake_raw.fam'


def sample_qc_mendel_ht_path(data_type: str, part: str) -> str:
    """
    Returns the path to mendel_errors results computed as part of the pedigree inference

    part should be one of:
    - all_errors
    - per_fam
    - per_sample

    :param data_type:
    :param part:
    :return:
    """
    return f'gs://gnomad/sample_qc/fam/gnomad_{data_type}_{part}_mendel.mt'


def merged_pedigrees_ht_path(data_type: str) -> str:
    return f'gs://gnomad/sample_qc/fam/gnomad_{data_type}_merged_ped.mt'


def get_topmed_shared_sites_ht_path(data_type: str) -> str:
    return f'gs://gnomad/sample_qc/topmed/gnomad_{data_type}_topmed_shared_sites.ht'
