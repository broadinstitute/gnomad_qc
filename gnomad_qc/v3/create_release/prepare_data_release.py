from gnomad.utils.annotations import add_variant_type
from gnomad.resources.grch38 import lcr_intervals
from gnomad.resources.grch38.reference_data import dbsnp
from gnomad.resources.grch38.gnomad import release_vcf_path
from gnomad.utils.reference_genome import get_reference_genome
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad_qc.v3.resources import freq, get_filtering_model, get_info, qual_hist, vep, V3_RELEASE_VERSION, release_ht_path
from gnomad.utils.vcf import INFO_DICT
from gnomad.utils.vep import vep_struct_to_csq
import argparse
import copy
import itertools
import logging
import hail as hl
from typing import Dict, List, Optional


logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("prepare_data_release")
logger.setLevel(logging.INFO)

HISTS = ['gq_hist_alt', 'gq_hist_all', 'dp_hist_alt', 'dp_hist_all', 'ab_hist_alt']

GROUPS = ['adj', 'raw']
SEXES = ['male', 'female']
POPS = ['afr', 'ami', 'amr', 'asj', 'eas', 'fin', 'nfe', 'oth', 'sas']
FAF_POPS = ['afr', 'amr', 'eas', 'nfe', 'sas']
NFE_SUBPOPS = ['onf', 'bgr', 'swe', 'nwe', 'seu', 'est']
EAS_SUBPOPS = ['kor', 'oea', 'jpn']

SORT_ORDER = ['popmax', 'group', 'pop', 'subpop', 'sex']

VEP_CSQ = "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|"
"HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|"
"ALLELE_NUM|DISTANCE|STRAND|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|"
"TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|"
"HGVS_OFFSET|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|LoF|LoF_filter|"
"LoF_flags|LoF_info"


def get_problematic_regions_flags_expr(
        ht: hl.Table,
        non_par: bool = True,
        prob_regions: Dict[str, hl.Table] = None
) -> Dict[str, hl.expr.BooleanExpression]:
    """
    Annotate HT with flags for problematic regions (i.e., LCR, decoy, segdup, and nonpar regions)
    :param Table ht: Table to be annotated
    :return: Table containing new annotations ['lcr', 'decoy', 'segdup', 'nonpar']
    :rtype: Table
    """

    prob_flags_expr = {'non_par': (ht.locus.in_x_nonpar() | ht.locus.in_y_nonpar())} if non_par else {}

    if prob_regions is not None:
        prob_flags_expr.update({
            region_name: hl.is_defined(region_table[ht.locus])
            for region_name, region_table in prob_regions.items()
        })

    return prob_flags_expr


def prepare_table_annotations(
        freq_ht: hl.Table,
        filter_ht: Optional[hl.Table],
        vep_ht: hl.Table,
        dbsnp_ht: hl.Table,
        hist_ht: Optional[hl.Table],
        index_dict,
        info_ht: hl.Table,
        include_project_max: bool = False
) -> hl.Table:
    """
    Join Tables with variant annotations for gnomAD release, dropping sensitive annotations and keeping only variants with nonzero AC
    :param Table freq_ht: Table with frequency annotations
    :param Table filter_ht: Table with random forest variant annotations
    :param Table vep_ht: Table with VEP variant annotations
    :param Table dbsnp_ht: Table with updated dbSNP rsid annotations
    :param Table hist_ht: Table with quality histogram annotations
    :param dict index_dict: Dictionary containing index values for each entry in the frequency HT array, keyed by metadata label
    :param Table info_ht: Table containing allele annotations
    :return: Table containing joined annotations
    :rtype: Table
    """
    freq_ht = hl.filter_intervals(freq_ht, [hl.parse_locus_interval('chrM')], keep=False) # TODO Make it backward compatible?

    if not include_project_max and 'project_max' in freq_ht.row_value:
        freq_ht = freq_ht.drop('project_max')

    raw_idx = index_dict['raw']
    # Remove sites with no release samples AC and lowqual sites
    # TODO: Is this the best place for this logic?
    ht = freq_ht.filter((freq_ht.freq[raw_idx].AC <= 0) | info_ht[freq_ht.key].lowqual, keep=False)

    # Dict that contains all annotations for the Table root
    root_annotations = {f: ht[f] for f in ht.row_value if f != 'info'}
    # Dict that contains all annotations for the Table info struct
    info_annotations = {**ht.info} if 'info' in ht.row_value else {}

    # Add problematic region flags (lcr, non-par here)
    # TODO: This is now added to info -- I think this makes more sense...but maybe not?
    info_annotations.update(
        get_problematic_regions_flags_expr(ht, prob_regions={'lcr': lcr_intervals.ht()})
    )

    # Add variant type annotations
    info_annotations.update({
        **add_variant_type(ht.alleles)
    })

    # Add annotations from the info HT
    info_indexed = info_ht[ht.key]
    root_annotations.update(
        dict(
            a_index=info_indexed.a_index,
            was_split=info_indexed.was_split,
            old_locus=info_indexed.old_locus,
            old_alleles=info_indexed.old_alleles,
            qual=info_indexed.info.QUALapprox
        )
    )

    info_annotations.update(
        {f: info_indexed.info[f] for f in info_indexed.info if f not in ['ac_raw', 'ac', 'QUALapprox']} # TODO make it more generic
    )

    # Add histograms if available
    if hist_ht is not None:
        root_annotations.update(
            **hist_ht[ht.key]
        )

    # Add filters and filtering annotations
    indexed_filter_ht = filter_ht[ht.key]
    root_annotations['filters'] = indexed_filter_ht.filters
    info_annotations.update(**indexed_filter_ht.info)

    # Add vep
    root_annotations['vep']=vep_ht[ht.key].vep

    # Add rsid
    root_annotations['rsid'] = dbsnp_ht[ht.key].rsid

    # Return annotated HT
    ht = ht.select(
        **root_annotations,
        info=hl.struct(**info_annotations)
    )

    # Add global annotations
    return ht.annotate_globals(filtering_model=filter_ht.index_globals().filtering_model)


def make_label_combos(label_groups: Dict[str, List[str]]) -> List[str]:
    '''
    Make combinations of all possible labels for a supplied dictionary of label groups
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :return: list of all possible combinations of values for the supplied label groupings
    :rtype: list[str]
    '''
    copy_label_groups = copy.deepcopy(label_groups)
    if len(copy_label_groups) == 1:
        return [item for sublist in copy_label_groups.values() for item in sublist]
    anchor_group = sorted(copy_label_groups.keys(), key=lambda x: SORT_ORDER.index(x))[0]
    anchor_val = copy_label_groups.pop(anchor_group)
    combos = []
    for x,y in itertools.product(anchor_val, make_label_combos(copy_label_groups)):
        combos.append('{0}_{1}'.format(x,y))
    return combos


def make_info_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    '''
    Make Hail expression for variant annotations to be included in VCF INFO
    :param Table ht: Table containing variant annotations to be reformatted for VCF export
    :return: Dictionary containing Hail expressions for relevant INFO annotations
    :rtype: Dict of str: Expression
    '''
    info_dict = {
        'FS': ht.info_FS,
        'InbreedingCoeff': ht.info_InbreedingCoeff,
        'MQ': ht.info_MQ,
        'MQRankSum': ht.info_MQRankSum,
        'QD': ht.info_QD,
        'ReadPosRankSum': ht.info_ReadPosRankSum,
        'SOR': ht.info_SOR,
        'VQSR_POSITIVE_TRAIN_SITE': ht.info_POSITIVE_TRAIN_SITE,
        'VQSR_NEGATIVE_TRAIN_SITE': ht.info_NEGATIVE_TRAIN_SITE,
        'BaseQRankSum': ht.allele_info.BaseQRankSum,
        'ClippingRankSum': ht.allele_info.ClippingRankSum,
        'DP': ht.allele_info.DP,
        'VQSLOD': ht.allele_info.VQSLOD,
        'VQSR_culprit': ht.allele_info.culprit,
        'segdup': ht.segdup,
        'lcr': ht.lcr,
        'decoy': ht.decoy,
        'nonpar': ht.nonpar,
        'rf_positive_label': ht.tp,
        'rf_negative_label': ht.fail_hard_filters,
        'rf_label': ht.rf_label,
        'rf_train': ht.rf_train,
        'rf_tp_probability': ht.rf_probability,
        'transmitted_singleton': ht.transmitted_singleton,
        'variant_type': ht.variant_type,
        'allele_type': ht.allele_type,
        'n_alt_alleles': ht.n_alt_alleles,
        'was_mixed': ht.was_mixed,
        'has_star': ht.has_star,
        'pab_max': ht.pab_max,
    }
    for hist in HISTS:
        hist_dict = {
            f'{hist}_bin_freq': hl.delimit(ht[hist].bin_freq, delimiter="|"),
            f'{hist}_bin_edges': hl.delimit(ht[hist].bin_edges, delimiter="|"),
            f'{hist}_n_smaller': ht[hist].n_smaller,
            f'{hist}_n_larger': ht[hist].n_larger
        }
        info_dict.update(hist_dict)
    return info_dict


def make_filters_sanity_check_expr(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    '''
    Make Hail Expressions to measure % variants filtered under varying conditions of interest
    :param Table ht: Hail Table containing 'filter' annotation to be examined
    :return: Dictionary containing Hail aggregation expressions to measure filter flags
    :rtype: Dict of str: Expression
    '''
    filters_dict = {
        'n': hl.agg.count(),
        'frac_any_filter': hl.agg.count_where(ht.is_filtered)/hl.agg.count(),
        'frac_inbreed_coeff': hl.agg.count_where(ht.filters.contains('inbreeding_coeff'))/hl.agg.count(),
        'frac_ac0': hl.agg.count_where(ht.filters.contains('AC0'))/hl.agg.count(),
        'frac_rf': hl.agg.count_where(ht.filters.contains('rf'))/hl.agg.count(),
        'frac_inbreed_coeff_only': hl.agg.count_where(ht.filters.contains('inbreeding_coeff') & (ht.filters.length() == 1))/hl.agg.count(),
        'frac_ac0_only': hl.agg.count_where(ht.filters.contains('AC0') & (ht.filters.length() == 1))/hl.agg.count(),
        'frac_rf_only': hl.agg.count_where(ht.filters.contains('rf') & (ht.filters.length() == 1))/hl.agg.count()
    }
    return filters_dict


def sample_sum_check(ht: hl.Table, prefix: str, label_groups: Dict[str, List[str]], verbose: bool, subpop=None):
    '''
    Compute afresh the sum of annotations for a specified group of annotations, and compare to the annotated version;
    display results from checking the sum of the specified annotations in the terminal
    :param Table ht: Hail Table containing annotations to be summed
    :param str prefix: Subset of gnomAD
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks
    :param str subpop: Subpop abbreviation, supplied only if subpopulations are included in the annotation groups being checked
    :rtype: None
    '''
    combo_AC = [ht.info[f'{prefix}_AC_{x}'] for x in make_label_combos(label_groups)]
    combo_AN = [ht.info[f'{prefix}_AN_{x}'] for x in make_label_combos(label_groups)]
    combo_nhomalt = [ht.info[f'{prefix}_nhomalt_{x}'] for x in make_label_combos(label_groups)]
    group = label_groups.pop('group')[0]
    alt_groups = "_".join(sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x)))

    annot_dict = {f'sum_AC_{group}_{alt_groups}': hl.sum(combo_AC),
                  f'sum_AN_{group}_{alt_groups}': hl.sum(combo_AN),
                  f'sum_nhomalt_{group}_{alt_groups}': hl.sum(combo_nhomalt)}

    ht = ht.annotate(**annot_dict)

    for subfield in ['AC', 'AN', 'nhomalt']:
        if not subpop:
            generic_field_check(ht, (ht.info[f'{prefix}_{subfield}_{group}'] != ht[f'sum_{subfield}_{group}_{alt_groups}']),
                                f'{prefix}_{subfield}_{group} = sum({subfield}_{group}_{alt_groups})',
                                [f'info.{prefix}_{subfield}_{group}', f'sum_{subfield}_{group}_{alt_groups}'], verbose)
        else:
            generic_field_check(ht, (ht.info[f'{prefix}_{subfield}_{group}_{subpop}'] != ht[f'sum_{subfield}_{group}_{alt_groups}']),
                                f'{prefix}_{subfield}_{group}_{subpop} = sum({subfield}_{group}_{alt_groups})',
                                [f'info.{prefix}_{subfield}_{group}_{subpop}', f'sum_{subfield}_{group}_{alt_groups}'], verbose)


def generic_field_check(ht: hl.Table, cond_expr, check_description, display_fields, verbose):
    '''
    Check a generic logical condition involving annotations in a Hail Table and print the results to terminal
    :param Table ht: Table containing annotations to be checked
    :param Expression cond_expr: logical expression referring to annotations in ht to be checked
    :param str check_description: String describing the condition being checked; is displayed in terminal summary message
    :param list of str display_fields: List of names of ht annotations to be displayed in case of failure (for troubleshooting purposes);
        these fields are also displayed if verbose is True
    :param bool verbose: If True, show top values of annotations being checked, including checks that pass; if False,
        show only top values of annotations that fail checks
    :rtype: None
    '''
    ht_orig = ht
    ht = ht.filter(cond_expr)
    n_fail = ht.count()
    if n_fail > 0:
        logger.info(f'Found {n_fail} sites that fail {check_description} check:')
        ht = ht.flatten()
        ht.select('locus', 'alleles', *display_fields).show()
    else:
        logger.info(f'PASSED {check_description} check')
        if verbose:
            ht_orig = ht_orig.flatten()
            ht_orig.select(*display_fields).show()


def sanity_check_ht(ht: hl.Table, data_type: str, subsets=None, missingness_threshold=0.5, verbose=False):
    '''
    Perform a battery of sanity checks on a specified group of subsets in a Hail Table containing variant annotations;
    includes summaries of % filter status for different partitions of variants; histogram outlier bin checks; checks on
    AC, AN, and AF annotations; checks that subgroup annotation values add up to the supergroup annotation values;
    checks on sex-chromosome annotations; and summaries of % missingness in variant annotations
    :param Table ht: Table containing variant annotations to check
    :param str data_type: 'exomes' or 'genomes'
    :param list of str subsets: List of subsets to be checked
    :param bool verbose: If True, display top values of relevant annotations being checked, regardless of whether check
        conditions are violated; if False, display only top values of relevant annotations if check conditions are violated
    :return: Terminal display of results from the battery of sanity checks
    :rtype: None
    '''

    xchrom = 'X' if get_reference_genome(ht.locus) == 'GRCh37' else 'chrX'
    ychrom = 'Y' if get_reference_genome(ht.locus) == 'GRCh37' else 'chrY'

    n_sites = ht.count()
    contigs = ht.aggregate(hl.agg.collect_as_set(ht.locus.contig))
    logger.info(f'Found {n_sites} sites for contigs {contigs}')
    info_metrics = list(ht.row.info)
    non_info_metrics = list(ht.row)
    non_info_metrics.remove('info')

    logger.info('VARIANT FILTER SUMMARIES:')
    ht = ht.annotate(is_filtered=ht.filters.length() > 0,
                     in_problematic_region=hl.any(lambda x: x, [ht.lcr]))
                     # in_problematic_region=hl.any(lambda x: x, [ht.info.lcr, ht.info.segdup, ht.info.decoy])) # TODO generalize

    ht_filter_check1 = ht.group_by(ht.is_filtered).aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n'))
    ht_filter_check1.show()

    ht_filter_check2 = ht.group_by(ht.info.allele_type).aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n'))
    ht_filter_check2.show()

    ht_filter_check3 = (ht.group_by(ht.info.allele_type, ht.in_problematic_region)
                        .aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n')))
    ht_filter_check3.show(50, 140)

    ht_filter_check4 = (ht.group_by(ht.info.allele_type, ht.in_problematic_region, ht.info.n_alt_alleles)
                        .aggregate(**make_filters_sanity_check_expr(ht)).order_by(hl.desc('n')))
    ht_filter_check4.show(50, 140)

    logger.info('HISTOGRAM CHECKS:')
    for hist in ['gq_hist_alt', 'gq_hist_all', 'dp_hist_alt', 'dp_hist_all', 'ab_hist_alt']:
        # Check subfield == 0
        generic_field_check(ht, (ht.info[f'{hist}_n_smaller'] != 0), f'{hist}_n_smaller == 0',
                            [f'info.{hist}_n_smaller'], verbose)
        if hist not in ['dp_hist_alt', 'dp_hist_all']:  # NOTE: DP hists can have nonzero values in n_larger bin
            generic_field_check(ht, (ht.info[f'{hist}_n_larger'] != 0), f'{hist}_n_larger == 0',
                                [f'info.{hist}_n_larger'], verbose)

    logger.info('RAW AND ADJ CHECKS:')
    for subfield in ['AC', 'AN', 'AF']:
        # Check AC, AN, AF > 0
        generic_field_check(ht, (ht.info[f'gnomad_{subfield}_raw'] <= 0), f'gnomad {subfield}_raw > 0',
                            [f'info.gnomad_{subfield}_raw'], verbose)
        generic_field_check(ht, (ht.info[f'gnomad_{subfield}_adj'] < 0), f'gnomad {subfield}_adj >= 0',
                            [f'info.gnomad_{subfield}_adj', 'filters'], verbose)

    if subsets is not None:
        for subset in subsets:
            for subfield in ['AC', 'AN', 'nhomalt']:
                # Check AC_raw >= AC adj
                generic_field_check(ht, (ht.info[f'{subset}_{subfield}_raw'] < ht.info[f'{subset}_{subfield}_adj']),
                                    f'{subset} {subfield}_raw >= {subfield}_adj',
                                    [f'info.{subset}_{subfield}_raw', f'info.{subset}_{subfield}_adj'], verbose)

        logger.info('SAMPLE SUM CHECKS:')
        for subset in subsets:
            # Check if pops are present
            pop_adjusted = list(set([x[-3:] for x in info_metrics if subset in x]))
            pop_adjusted = [y for y in pop_adjusted if y in POPS]
            missing_pops = [z for z in POPS if z not in pop_adjusted]
            logger.warn(f'Missing {missing_pops} pops in {subset} subset!')

            # NOTE: group entry in dict is required here to be a list of length 1
            sample_sum_check(ht, subset, dict(group=['adj'], pop=pop_adjusted), verbose)
            sample_sum_check(ht, subset, dict(group=['adj'], sex=SEXES), verbose)
            sample_sum_check(ht, subset, dict(group=['adj'], pop=pop_adjusted, sex=SEXES), verbose)

            # Adjust subpops to those found in subset
            nfe_subpop_adjusted = list(set([x[-3:] for x in info_metrics if subset in x and 'nfe_' in x and 'male' not in x]))
            if nfe_subpop_adjusted != []:
                sample_sum_check(ht, subset, dict(group=['adj'], pop=['nfe'], subpop=nfe_subpop_adjusted), verbose, subpop='nfe')
            if data_type == 'exomes':
                eas_subpop_adjusted = list(set([x[-3:] for x in info_metrics if subset in x and 'eas_' in x and 'male' not in x]))
                if eas_subpop_adjusted != []:
                    sample_sum_check(ht, subset, dict(group=['adj'], pop=['eas'], subpop=eas_subpop_adjusted), verbose, subpop='eas')

    logger.info('SEX CHROMOSOME ANNOTATION CHECKS:')
    female_metrics = [x for x in info_metrics if '_female' in x]
    male_metrics = [x for x in info_metrics if '_male' in x]

    if ychrom in contigs:
        logger.info('Check values of female metrics for Y variants are NA:')
        ht_y = hl.filter_intervals(ht, [hl.parse_locus_interval(ychrom)]) # TODO: Make generic
        metrics_values = {}
        for metric in female_metrics:
            metrics_values[metric] = hl.agg.collect_as_set(ht_y.info[metric])
        output = ht_y.aggregate(hl.struct(**metrics_values))
        for metric,values in dict(output).items():
            if values == {None}:
                logger.info(f"PASSED {metric} = {None} check for {ychrom} variants")
            else:
                logger.info(f"FAILED Y check: Found {values} in {metric}")

    logger.info(f'Check values of male nhomalt metrics for {xchrom} nonpar variants are 0:')
    ht_x = hl.filter_intervals(ht, [hl.parse_locus_interval(xchrom)])
    ht_xnonpar = ht_x.filter(ht_x.locus.in_x_nonpar())
    n = ht_xnonpar.count()
    logger.info(f"Found {n} X nonpar sites")  # Lots of values found in male X nonpar sites

    male_metrics = [x for x in male_metrics if 'nhomalt' in x]
    metrics_values = {}
    for metric in male_metrics:
        metrics_values[metric] = hl.agg.collect_as_set(ht_xnonpar.info[metric])
    output = ht_xnonpar.aggregate(hl.struct(**metrics_values))
    for metric, values in dict(output).items():
        if values == {0}:
            logger.info(f"PASSED {metric} = 0 check for {xchrom} nonpar variants")
        else:
            logger.info(f"FAILED {xchrom} nonpar check: Found {values} in {metric}")

    logger.info(f'Check (nhomalt == nhomalt_female) for {xchrom} nonpar variants:')
    female_metrics = [x for x in female_metrics if 'nhomalt' in x]
    for metric in female_metrics:
        standard_field = metric.replace('_female', '')
        generic_field_check(ht_xnonpar, (ht_xnonpar.info[f'{metric}'] != ht_xnonpar.info[f'{standard_field}']),
                            f'{metric} == {standard_field}', [f'info.{metric}', f'info.{standard_field}'], verbose)

    logger.info('MISSINGNESS CHECKS:')
    metrics_frac_missing = {}
    for x in info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht.info[x])) / n_sites
    for x in non_info_metrics:
        metrics_frac_missing[x] = hl.agg.sum(hl.is_missing(ht[x])) / n_sites
    output = ht.aggregate(hl.struct(**metrics_frac_missing))

    n_fail = 0
    for metric, value in dict(output).items():
        if value > missingness_threshold:
            logger.info("FAILED missing check for {}: {}% missing".format(metric, 100 * value))
            n_fail += 1
        else:
            logger.info("Missingness check for {}: {}% missing".format(metric, 100 * value))
    logger.info("{} missing metrics checks failed".format(n_fail))


def make_freq_meta_index_dict(freq_meta):
    '''
    Make dictionary of the entries in the frequency array annotation, where keys are the grouping combinations and the values
    are the 0-based integer indices
    :param list of str freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    index_dict = index_globals(freq_meta, dict(group=GROUPS))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=['nfe'], subpop=NFE_SUBPOPS)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=['eas'], subpop=EAS_SUBPOPS)))
    return index_dict


def index_globals(freq_meta, label_groups):
    '''
    Create a dictionary keyed by the specified label groupings with values describing the corresponding index of each grouping entry
    in the freq_meta array annotation
    :param list of str freq_meta: Ordered list containing string entries describing all the grouping combinations contained in the
        frequency array annotation
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :return: Dictionary keyed by specified label grouping combinations, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    combos = make_label_combos(label_groups)
    index_dict = {}

    for combo in combos:
        combo_fields = combo.split("_")
        for i,v in enumerate(freq_meta):
            if set(v.values()) == set(combo_fields):
                index_dict.update({f"{combo}": i})
    return index_dict


def unfurl_nested_annotations(ht):
    '''
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays, where the values
    of the dictionary are Hail Expressions describing how to access the corresponding values
    :param Table ht: Hail Table containing the nested variant annotation arrays to be unfurled
    :return: Dictionary containing variant annotations and their corresponding values
    :rtype: Dict of str: Expression
    '''
    expr_dict = dict()
    for k, i in hl.eval(ht.globals.freq_index_dict).items(): # FIXME: This was much more complex -- but I think this is correct?
        if k == 'raw':
            k = 'adj_raw'

        combo_dict = {
            k.replace('adj', 'AC'): ht.freq[i].AC,
            k.replace('adj', 'AN'): ht.freq[i].AN,
            k.replace('adj', 'AF'): ht.freq[i].AF,
            k.replace('adj', 'nhomalt'): ht.freq[i].homozygote_count
        }
        expr_dict.update(combo_dict)

    for k, i in hl.eval(ht.globals.faf_index_dict).items():  # NOTE: faf annotations are all done on adj-only groupings
        entry = k.split("_")
        if entry[0] == "non":
            prefix = "_".join(entry[:2])
            combo_fields = ['adj'] + entry[2:]
        else:
            prefix = entry[0]
            combo_fields = ['adj'] + entry[1:]
        combo = "_".join(combo_fields)
        combo_dict = {
            f"{prefix}_faf95_{combo}": hl.or_missing(hl.set(ht.faf[i].meta.values()) == set(combo_fields),
                                                     ht.faf[i].faf95),
            f"{prefix}_faf99_{combo}": hl.or_missing(hl.set(ht.faf[i].meta.values()) == set(combo_fields),
                                                     ht.faf[i].faf99),
        }
        expr_dict.update(combo_dict)

    if 'popmax_index_dict' in ht.globals:
        for prefix, i in hl.eval(ht.globals.popmax_index_dict).items():
            combo_dict = {
                f'{prefix}_popmax': ht.popmax[i].pop,
                f'{prefix}_AC_popmax': ht.popmax[i].AC,
                f'{prefix}_AN_popmax': ht.popmax[i].AN,
                f'{prefix}_AF_popmax': ht.popmax[i].AF,
                f'{prefix}_nhomalt_popmax': ht.popmax[i].homozygote_count,
            }
            expr_dict.update(combo_dict)
            if prefix == 'gnomad':
                age_hist_dict = {
                    f"{prefix}_age_hist_het_bin_freq": hl.delimit(ht.age_hist_het[i].bin_freq, delimiter="|"),
                    f"{prefix}_age_hist_het_bin_edges": hl.delimit(ht.age_hist_het[i].bin_edges, delimiter="|"),
                    f"{prefix}_age_hist_het_n_smaller": ht.age_hist_het[i].n_smaller,
                    f"{prefix}_age_hist_het_n_larger": ht.age_hist_het[i].n_larger,
                    f"{prefix}_age_hist_hom_bin_freq": hl.delimit(ht.age_hist_hom[i].bin_freq, delimiter="|"),
                    f"{prefix}_age_hist_hom_bin_edges": hl.delimit(ht.age_hist_hom[i].bin_edges, delimiter="|"),
                    f"{prefix}_age_hist_hom_n_smaller": ht.age_hist_hom[i].n_smaller,
                    f"{prefix}_age_hist_hom_n_larger": ht.age_hist_hom[i].n_larger
                }
                expr_dict.update(age_hist_dict)
    return expr_dict


def make_combo_header_text(preposition, group_types, combo_fields, prefix, faf=False):
    '''
    Programmatically generate text to populate the VCF header description for a given variant annotation with specific groupings and subset
    :param str preposition: Relevant preposition to precede automatically generated text
    :param list of str group_types: List of grouping types, e.g. "sex" or "pop"
    :param list of str combo_fields: List of the specific values for each grouping type, for which the text is being generated
    :param str prefix: Subset of gnomAD
    :param bool faf: If True, use alternate logic to automatically populate descriptions for filter allele frequency annotations
    :return: String with automatically generated description text for a given set of combo fields
    :rtype: str
    '''
    combo_dict = dict(zip(group_types, combo_fields))

    if not faf:
        header_text = " " + preposition
        if 'sex' in combo_dict.keys():
            header_text = header_text + " " + combo_dict['sex']
        header_text = header_text + " samples"
        if 'subpop' in combo_dict.keys():
            header_text = header_text + f" of {POP_NAMES[combo_dict['subpop']]} ancestry"
            combo_dict.pop('pop')
        if 'pop' in combo_dict.keys():
            header_text = header_text + f" of {POP_NAMES[combo_dict['pop']]} ancestry"
        if prefix != 'gnomad':
            header_text = header_text + f" in the {prefix} subset"
        if 'group' in group_types:
            if combo_dict['group'] == 'raw':
                header_text = header_text + ", before removing low-confidence genotypes"
    else:
        header_text = ""
        if 'pop' in combo_dict.keys():
            header_text = f" of {POP_NAMES[combo_dict['pop']]} ancestry"
        if prefix != 'gnomad':
            header_text = header_text + f" in the {prefix} subset"
    return header_text


def make_info_dict(prefix, label_groups=None, bin_edges=None, faf=False, popmax=False, age_hist_data=None): # TODO: This needs restructuring
    '''
    Generate dictionary of Number and Description attributes to be used in the VCF header
    :param str prefix: Subset of gnomAD
    :param dict label_groups: Dictionary containing an entry for each label group, where key is the name of the grouping,
        e.g. "sex" or "pop", and value is a list of all possible values for that grouping (e.g. ["male", "female"] or ["afr", "nfe", "amr"])
    :param dict bin_edges: Dictionary keyed by annotation type, with values that reflect the bin edges corresponding to the annotation
    :param bool faf: If True, use alternate logic to auto-populate dictionary values associated with filter allele frequency annotations
    :param bool popmax: If True, use alternate logic to auto-populate dictionary values associated with popmax annotations
    :param str age_hist_data: Pipe-delimited string of age histograms, from get_age_distributions (somewhat required if popmax == True)
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes
    :rtype: Dict of str: (Dict of str: str)
    '''
    info_dict = dict()

    if popmax: # TODO: review this condition
        popmax_text = "" if prefix == 'gnomad' else f" in the {prefix} subset"
        popmax_dict = {
            f'{prefix}_popmax': {"Number": "A",
                                 "Description": "Population with maximum AF{}".format(popmax_text)},
            f'{prefix}_AC_popmax': {"Number": "A",
                                    "Description": "Allele count in the population with the maximum AF{}".format(popmax_text)},
            f'{prefix}_AN_popmax': {"Number": "A",
                                    "Description": "Total number of alleles in the population with the maximum AF{}".format(popmax_text)},
            f'{prefix}_AF_popmax': {"Number": "A",
                                    "Description": "Maximum allele frequency across populations (excluding samples of Ashkenazi, Finnish, and indeterminate ancestry){}".format(popmax_text)},
            f'{prefix}_nhomalt_popmax': {"Number": "A",
                                         "Description": "Count of homozygous individuals in the population with the maximum allele frequency{}".format(popmax_text)}
        }
        info_dict.update(popmax_dict)
        if age_hist_data is not None:
            age_hist_dict = {
                f"{prefix}_age_hist_het_bin_freq": {"Number": "A",
                                                    "Description": f"Histogram of ages of heterozygous individuals; bin edges are: {bin_edges[f'{prefix}_het']}; total number of individuals of any genotype bin: {age_hist_data}"},
                f"{prefix}_age_hist_het_n_smaller": {"Number": "A",
                                                     "Description": "Count of age values falling below lowest histogram bin edge for heterozygous individuals"},
                f"{prefix}_age_hist_het_n_larger": {"Number": "A",
                                                    "Description": "Count of age values falling above highest histogram bin edge for heterozygous individuals"},
                f"{prefix}_age_hist_hom_bin_freq": {"Number": "A",
                                                    "Description": f"Histogram of ages of homozygous alternate individuals; bin edges are: {bin_edges[f'{prefix}_hom']}; total number of individuals of any genotype bin: {age_hist_data}"},
                f"{prefix}_age_hist_hom_n_smaller": {"Number": "A",
                                                     "Description": "Count of age values falling below lowest histogram bin edge for homozygous alternate individuals"},
                f"{prefix}_age_hist_hom_n_larger": {"Number": "A",
                                                    "Description": "Count of age values falling above highest histogram bin edge for homozygous alternate individuals"}
            }
            info_dict.update(age_hist_dict)
    else:
        group_types = sorted(label_groups.keys(), key=lambda x: SORT_ORDER.index(x))
        combos = make_label_combos(label_groups)

        for combo in combos:
            combo_fields = combo.split("_")
            if not faf:
                combo_dict = {
                    f"{prefix}_AC_{combo}": {"Number": "A",
                                             "Description": "Alternate allele count{}".format(make_combo_header_text('for', group_types, combo_fields, prefix))},
                    f"{prefix}_AN_{combo}": {"Number": "1",
                                             "Description": "Total number of alleles{}".format(make_combo_header_text('in', group_types, combo_fields, prefix))},
                    f"{prefix}_AF_{combo}": {"Number": "A",
                                             "Description": "Alternate allele frequency{}".format(make_combo_header_text('in', group_types, combo_fields, prefix))},
                    f"{prefix}_nhomalt_{combo}": {"Number": "A",
                                                  "Description": "Count of homozygous individuals{}".format(make_combo_header_text('in', group_types, combo_fields, prefix))}
                }
            else:
                combo_dict = {
                    f"{prefix}_faf95_{combo}": {"Number": "A",
                                                "Description": "Filtering allele frequency (using Poisson 95% CI) for samples{}".format(make_combo_header_text(None, group_types, combo_fields, prefix, faf=True))},
                    f"{prefix}_faf99_{combo}": {"Number": "A",
                                                "Description": "Filtering allele frequency (using Poisson 99% CI) for samples{}".format(make_combo_header_text(None, group_types, combo_fields, prefix, faf=True))}
                }
            info_dict.update(combo_dict)
    return info_dict


def make_hist_bin_edges_expr(ht, make_age_hist: bool = True):
    '''
    Create dictionary containing variant histogram annotations and their associated bin edges, formatted into a string
    separated by pipe delimiters
    :param Table ht: Table containing histogram variant annotations
    :return: Dictionary keyed by histogram annotation name, with corresponding reformatted bin edges for values
    :rtype: Dict of str: str
    '''
    edges_dict = {}
    if make_age_hist:
        edges_dict = {'gnomad_het': '|'.join(map(lambda x: f'{x:.1f}', ht.take(1)[0].age_hist_het[0].bin_edges)),
                      'gnomad_hom': '|'.join(map(lambda x: f'{x:.1f}', ht.take(1)[0].age_hist_hom[0].bin_edges))}
    for hist in HISTS:
        edges_dict[hist] = '|'.join(map(lambda x: f'{x:.2f}', ht.take(1)[0][hist].bin_edges)) if 'ab' in hist else \
            '|'.join(map(lambda x: str(int(x)), ht.take(1)[0][hist].bin_edges))
    return edges_dict


def make_hist_dict(bin_edges):
    '''
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for histogram annotations
    :param dict bin_edges: Dictionary keyed by histogram annotation name, with corresponding string-reformatted bin edges for values
    :return: Dictionary keyed by VCF INFO annotations, where values are Dictionaries of Number and Description attributes
    :rtype: Dict of str: (Dict of str: str)
    '''
    header_hist_dict = {}
    for hist in HISTS:
        edges = bin_edges[hist]
        hist_fields = hist.split("_")
        hist_text = hist_fields[0].upper()
        if hist_fields[2] == "alt":
            hist_text = hist_text + " in heterozygous individuals"
        hist_dict = {
            f'{hist}_bin_freq': {"Number": "A",
                                 "Description": f"Histogram for {hist_text}; bin edges are: {edges}"},
            f'{hist}_n_smaller': {"Number": "A",
                                  "Description": f"Count of {hist_fields[0].upper()} values falling below lowest histogram bin edge"},
            f'{hist}_n_larger': {"Number": "A",
                                  "Description": f"Count of {hist_fields[0].upper()} values falling above highest histogram bin edge"}
        }
        header_hist_dict.update(hist_dict)
    return header_hist_dict


def make_filter_dict(ht): # TODO: Make this generic
    '''
    Generate dictionary of Number and Description attributes to be used in the VCF header, specifically for FILTER annotations
    :param Table ht: Table containing global annotations of the Random Forests SNP and indel cutoffs
    :return: Dictionary keyed by VCF FILTER annotations, where values are Dictionaries of Number and Description attributes
    :rtype: Dict of str: (Dict of str: str)
    '''
    # snp_cutoff = hl.eval(ht.globals.rf.rf_snv_cutoff)
    # indel_cutoff = hl.eval(ht.globals.rf.rf_indel_cutoff)
    filter_dict = {
        'AC0': {'Description': 'Allele count is zero after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)'},
        'InbreedingCoeff': {'Description': 'InbreedingCoeff < -0.3'},
        'AS_VQSR': {'Description': 'Failed Allele-Specific Variant Quality Score Recalibration threshold'},
        # 'RF': {'Description': 'Failed random forest filtering thresholds of {0}, {1} (probabilities of being a true positive variant) for SNPs, indels'.format(snp_cutoff.min_score, indel_cutoff.min_score)},
        'PASS': {'Description': 'Passed all variant filters'}
    }
    return filter_dict


def make_index_dict(ht):
    '''
    Create a look-up Dictionary for entries contained in the frequency annotation array
    :param Table ht: Table containing freq_meta global annotation to be indexed
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    '''
    freq_meta = hl.eval(ht.globals.freq_meta)
    index_dict = make_freq_meta_index_dict(freq_meta)
    return index_dict


def set_female_y_metrics_to_na(ht):
    '''
    Set AC, AN, and nhomalt Y variant annotations for females to NA (instead of 0)
    :param Table ht: Hail Table containing female variant annotations
    :return: Table with reset annotations
    :rtype: Table
    '''
    metrics = list(ht.row.info)
    female_metrics = [x for x in metrics if '_female' in x]
    female_metrics = [x for x in female_metrics if ('nhomalt' in x) or ('AC' in x) or ('AN' in x)]

    female_metrics_dict = {}
    for metric in female_metrics:
        female_metrics_dict.update({f'{metric}': hl.cond(ht.locus.contig == 'Y', hl.null(hl.tint32), ht.info[f'{metric}'])})
    ht = ht.annotate(info=ht.info.annotate(**female_metrics_dict))
    return ht


def get_array_lengths(ht, subsets):
    '''
    Utility function to measure the array lengths of each frequency variant annotation for each subset in a nested Hail Table
    _before_ subsets are concatenated
    :param Table ht: Hail Table containing frequency array annotations for each gnomAD subset
    :param list of str subsets: List of subsets to be measured
    :return: List of the lengths of each frequency array, in the order of the supplied subset list
    :rtype: List of int
    '''
    l = [len(ht.take(1)[0].freq)]
    for subset in subsets:
        l.append(len(ht.take(1)[0][subset].freq))
    return l


def build_new_index_dict(ht, subsets, array_lengths):
    '''
    Create a new index dictionary for the concatenated frequency array annotation (i.e., after frequency
    arrays for each gnomAD subset have been concatenated into a single frequency array annotation)
    :param Table ht: Hail Table containing global frequency index dictionary annotations for each subset
    :param list of str subsets: List of subsets to be combined for the new index dictionary
    :param list of int array_lengths: List of the lengths of the frequency array for each subset to be concatenated
    :return: Dictionary keyed by subset grouping annotations, where values are the corresponding 0-based indices for the
        subset groupings in the new concatenated frequency array annotation
    :rtype: Dict of str: int
    '''
    temp_index_dict = hl.eval(ht.globals.freq_index_dict)
    new_index_dict = {'gnomad_' + k: v for k, v in temp_index_dict.items()}
    n = array_lengths[0]
    i = 0

    for subset in subsets:
        i += 1
        dict_name = subset + '_index_dict'
        temp_index_dict = hl.eval(ht.globals[dict_name])
        new_index_dict.update({subset + '_' + k: v + n for k, v in temp_index_dict.items()})
        n += array_lengths[i]
    return new_index_dict


def concat_array_expr(ht, subsets, field):
    '''
    Create Hail Expression to concatenate specified variant annotations for a list of subsets into a single variant annotation,
    where annotations to be combined are arrays
    :param Table ht: Table containing annotations to be combined
    :param list of str subsets: List of subsets whose variant annotations should be combined
    :param str field: Variant annotation (array) to be concatenated
    :return: Hail Expression to concatenate specified variant annotations across a list of subsets
    :rtype: Expression
    '''
    copy_subsets = copy.deepcopy(subsets)
    subset = copy_subsets.pop()

    if len(copy_subsets) == 0:
        expr = ht[field].extend(ht[subset][field])
    else:
        expr = concat_array_expr(ht, copy_subsets, field).extend(ht[subset][field])
    return expr


def concat_struct_expr(ht, subsets, field):
    '''
    Create Hail Expression to concatenate specified variant annotations for a list of subsets into a single variant annotation,
    where annotations to be combined are structs
    :param Table ht: Table containing annotations to be combined
    :param list of str subsets: List of subsets whose variant annotations should be combined
    :param str field: Variant annotation (Struct) to be concatenated
    :return: Hail Expression to concatenate specified variant annotations across a list of subsets
    :rtype: Expression
    '''
    copy_subsets = copy.deepcopy(subsets)
    subset = copy_subsets.pop()

    if len(copy_subsets) == 0:
        expr = hl.array([ht[field]]).extend(hl.array([ht[subset][field]]))
    else:
        expr = concat_struct_expr(ht, copy_subsets, field).extend(hl.array([ht[subset][field]]))
    return expr


def make_faf_index_dict(ht, subset=None):  # NOTE: label groups: assumes 'adj', plus all pops
    '''
    Create a look-up Dictionary for entries contained in the filter allele frequency annotation array
    :param Table ht: Table containing filter allele frequency annotations to be indexed
    :param bool subset: If True, use alternate logic to access the correct faf annotation
    :return: Dictionary of faf annotation population groupings, where values are the corresponding 0-based indices for the
        groupings in the faf array
    :rtype: Dict of str: int
    '''
    combos = make_label_combos(dict(group=['adj'], pop=POPS))
    combos = combos + ['adj']
    index_dict = {}
    if subset:
        faf_meta = hl.eval(hl.map(lambda f: f.meta, ht.take(1)[0][subset].faf))
    else:
        faf_meta = hl.eval(hl.map(lambda f: f.meta, ht.take(1)[0].faf))

    for combo in combos:
        combo_fields = combo.split("_")
        for i, v in enumerate(faf_meta):
            if set(v.values()) == set(combo_fields):
                index_dict.update({f"{combo}": i})
    return index_dict


def build_faf_index_dict(ht, subsets):
    '''
    Create a new index dictionary for the concatenated filter allele frequency array annotation (i.e., after faf
    arrays for each gnomAD subset have been concatenated into a single faf annotation)
    :param Table ht: Hail Table containing global faf index dictionary annotations for each subset
    :param list of str subsets: List of subsets to be combined for the new faf index dictionary
    :return: Dictionary keyed by subset grouping annotations, where values are the corresponding 0-based indices for the
        subset groupings in the new concatenated faf annotation
    :rtype: Dict of str: int
    '''
    temp_faf_index = make_faf_index_dict(ht)
    n = len(temp_faf_index)
    new_index_dict = {'gnomad_' + k: v for k, v in temp_faf_index.items()}

    for subset in subsets:
        temp_faf_index = make_faf_index_dict(ht, subset=subset)
        new_index_dict.update({subset + '_' + k: v + n for k, v in temp_faf_index.items()})
        n += len(temp_faf_index)
    return new_index_dict


def get_gnomad_ht(filtering_model_id: str, pcsk9: bool):
    freq_ht = freq.ht()
    index_dict = make_index_dict(freq_ht)

    if pcsk9:
        freq_ht = hl.filter_intervals(freq_ht, [hl.parse_locus_interval('chr1:55030000-55070000')])

    freq_ht = freq_ht.filter(freq_ht.freq[1].AC > 0)  # TODO: This should have been done during freq computation...
    freq_ht = freq_ht.repartition(5000, shuffle=False)
    freq_ht = freq_ht.checkpoint('gs://gnomad-tmp/v3_internal.ht', _read_if_exists=True)

    filters_ht = get_filtering_model(model_id=filtering_model_id).ht()
    vep_ht = vep.ht()
    dbsnp_ht = dbsnp.ht()

    hist_ht = qual_hist.ht()
    hist_ht = hist_ht.select('gq_hist_alt', 'gq_hist_all', 'dp_hist_alt', 'dp_hist_all', 'ab_hist_alt')
    info_ht = get_info().ht()
    # allele_ht = hl.read_table(annotations_ht_path(data_type, 'allele_data'))

    logger.info('Adding annotations...')
    ht = prepare_table_annotations(freq_ht, filters_ht, vep_ht, dbsnp_ht, hist_ht, index_dict, info_ht)

    ht = ht.annotate(
        qual=hl.float64(ht.qual),
        info=ht.info.annotate(
            QD=hl.or_missing(hl.is_finite(ht.info.QD), ht.info.QD) # TODO: Check this upstream
        )
    )

    ht = ht.annotate_globals(freq_index_dict=make_index_dict(ht), faf_index_dict=make_faf_index_dict(ht))
                             # age_distribution=age_hist_data)
    return ht


def prepare_release_vcf(
        ht: hl.Table,
        subset_list: List,
        public: bool,
        pcsk9: bool = False,
        age_hist_data=None # TODO: Figure out how to best handle this -- looks like this may already be in the HT globals?
):

    bin_edges = make_hist_bin_edges_expr(ht, make_age_hist=age_hist_data is not None)

    if age_hist_data is not None:
        age_hist_data = '|'.join(str(x) for x in age_hist_data)

    # Make INFO dictionary for VCF
    for subset in subset_list:
        INFO_DICT.update(make_info_dict(subset, bin_edges=bin_edges, popmax=True,
                                        age_hist_data=age_hist_data if subset == 'gnomad' else None))
        INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS)))
        INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS, sex=SEXES)))
        INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS, pop=POPS)))
        INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS, pop=POPS, sex=SEXES)))
        # INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS, pop=['nfe'], subpop=NFE_SUBPOPS)))
        # INFO_DICT.update(make_info_dict(subset, dict(group=GROUPS, pop=['eas'], subpop=EAS_SUBPOPS)))
        INFO_DICT.update(make_info_dict(subset, dict(group=['adj']), faf=True))
        INFO_DICT.update(make_info_dict(subset, dict(group=['adj'], pop=FAF_POPS), faf=True))
    INFO_DICT.update(make_hist_dict(bin_edges))

    # Adjust keys to remove gnomad, adj tags before exporting to VCF
    _info_dict = {i.replace('gnomad_', '').replace('_adj', ''): j for i, j in INFO_DICT.items()}

    # Construct INFO field
    # ht = ht.annotate(info=hl.struct(**make_info_expr(ht))) # Not needed any more
    ht = ht.annotate(info=ht.info.annotate(**unfurl_nested_annotations(ht)))
    ht = set_female_y_metrics_to_na(ht)

    # Select relevant fields for VCF export
    ht = ht.select('info', 'filters', 'rsid', 'qual', 'vep')
    # ht.write(release_ht_path(public=public), overwrite=overwrite)

    # # Move 'info' annotations to top level for browser release
    # ht = hl.read_table(release_ht_path(data_type, nested=False, temp=True))
    # ht = ht.transmute(**ht.info)
    # ht = ht.select_globals('rf')
    # ht.write(release_ht_path(data_type, nested=False), args.overwrite)


    # ht = hl.read_table(release_ht_path(public=public))
    contigs = hl.eval(ht.aggregate(hl.agg.collect_as_set(ht.locus.contig)))

    # Annotate VEP CSQ
    ht = ht.transmute(
        info=ht.info.annotate(
            vep=vep_struct_to_csq(ht.vep, VEP_CSQ)
        )
    )
    _info_dict['vep'] = {'Description': VEP_CSQ}

    # Remove gnomad_ prefix for VCF export
    row_annots = list(ht.row.info)
    new_row_annots = [x.replace('gnomad_', '').replace('adj_', '') for x in row_annots] # FIXME somehow this was _adj before ???
    info_annot_mapping = dict(zip(new_row_annots, [ht.info[f'{x}'] for x in row_annots]))
    ht = ht.transmute(info=hl.struct(**info_annot_mapping))

    # Rearrange INFO field in desired ordering
    # TODO: These fields are never in info in v3 -- I think it's the right behavior -- might want to update v2 code?
    # drop_hists = [x + '_n_smaller' for x in HISTS] + [x + '_bin_edges' for x in HISTS] + [x + '_n_larger' for x in HISTS if 'dp_' not in x] # + ['age_hist_hom_bin_edges', 'age_hist_het_bin_edges']
    first_fields = ['AC', 'AN', 'AF']
    if 'rf_tp_probability' in ht.row_value:
        first_fields.append('rf_tp_probability')
    ht = ht.annotate(
        info=ht.info.select(
            *first_fields,
            *ht.info.drop(*first_fields)
            # *ht.info.drop(*first_fields, *drop_hists)
        )
    )

    if 'SB' in ht.info:
        ht = ht.annotate(
            info=ht.info.annotate(
                SB=ht.info.SB[0].extend(ht.info.SB[1])
            )
        )
    _info_dict['SB'] = {'Number': '4'}

    # ht = ht.annotate(
    #     qual=hl.float64(ht.qual)
    # )


    # Add VEP annotations
    # vep_csq_ht = hl.read_table(annotations_ht_path(data_type, 'vep_csq'))
    # new_info_dict.update({'vep': {'Description': hl.eval(vep_csq_ht.globals.vep_csq_header)}})
    header_dict = {'info': _info_dict,
                   'filter': make_filter_dict(ht)}
    # ht = ht.annotate(info=ht.info.annotate(vep=vep_csq_ht[ht.key].vep))

    # Export VCFs, full and by chromosome
    gnomad_ref = get_reference_genome(ht.locus)
    # ht = ht.key_by(locus=hl.locus(ht.locus.contig, ht.locus.position, reference_genome=gnomad_ref),
    #                alleles=ht.alleles)

    if pcsk9:
        ht = hl.filter_intervals(ht, [hl.parse_locus_interval('chr1:55030000-55070000')])
    else:
        for contig in contigs:
            contig_ht = hl.filter_intervals(ht, [hl.parse_locus_interval(contig, reference_genome=gnomad_ref)])
            mt = hl.MatrixTable.from_rows_table(contig_ht).key_cols_by(s='foo')
            hl.export_vcf(mt, release_vcf_path('genomes', version=V3_RELEASE_VERSION, contig=contig), metadata=header_dict)

    mt = hl.MatrixTable.from_rows_table(ht).key_cols_by(s='foo')
    hl.export_vcf(mt, release_vcf_path('genomes', version=V3_RELEASE_VERSION), metadata=header_dict)

    # # Export genome VCFs containing variants in exome calling intervals only # TODO: Reactivate once we have v4 calling intervals
    # if data_type == 'genomes':
    #     intervals = hl.import_locus_intervals(exome_calling_intervals_path, skip_invalid_intervals=True, reference_genome=gnomad_ref)
    #     coding_mt = mt.filter_rows(hl.is_defined(intervals[mt.locus]), keep=True)
    #     hl.export_vcf(coding_mt, release_vcf_path(data_type, coding_only=True), metadata=header_dict)

def main(args):
    hl.init(log='/release.log', default_reference='GRCh38')

    # data_type = 'genomes' if args.genomes else 'exomes'
    # import_ht_path = release_ht_path(data_type) if args.include_subset_frequencies else release_ht_path(data_type, with_subsets=False)
    # age_hist_data = get_age_distributions(data_type)

    if args.prepare_internal_ht:
        get_gnomad_ht(args.model_id, args.pcsk9).write(
            release_ht_path(public=args.public) if not args.pcsk9 else 'gs://gnomad-tmp/gnomad_internal_pcks9.ht',
            overwrite=args.overwrite
        )


    # if args.add_subset_frequencies:
    #     ht = hl.read_table(release_ht_path(data_type, with_subsets=False))
    #     controls_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies_control_with_consanguineous'))  # FIXME: revert to plain 'frequencies' for v3+
    #     freq_index_dict = make_index_dict(controls_ht)
    #     print(hl.eval(controls_ht.freq_meta))
    #     controls_ht = controls_ht.drop('helloworld')
    #
    #     non_neuro_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies_neuro_with_consanguineous'))  # FIXME: revert to plain 'frequencies' for v3+
    #     non_topmed_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies_topmed_with_consanguineous'))  # FIXME: revert to plain 'frequencies' for v3+
    #     controls_freq_meta = hl.eval(controls_ht.globals.freq_meta)
    #     [x.update({'subset': 'controls'}) for x in controls_freq_meta]
    #     non_neuro_freq_meta = hl.eval(non_neuro_ht.globals.freq_meta)
    #     [x.update({'subset': 'non_neuro'}) for x in non_neuro_freq_meta]
    #     non_topmed_freq_meta = hl.eval(non_topmed_ht.globals.freq_meta)
    #     [x.update({'subset': 'non_topmed'}) for x in non_topmed_freq_meta]
    #
    #     combined_ht = ht.annotate(controls=controls_ht[ht.key], non_neuro=non_neuro_ht[ht.key], non_topmed=non_topmed_ht[ht.key])
    #     combined_ht = combined_ht.annotate_globals(controls_index_dict=make_index_dict(controls_ht),
    #                                                non_neuro_index_dict=make_index_dict(non_neuro_ht),
    #                                                non_topmed_index_dict=make_index_dict(non_topmed_ht),
    #                                                freq_meta=combined_ht.freq_meta.extend(controls_freq_meta).extend(non_neuro_freq_meta).extend(non_topmed_freq_meta))
    #
    #     if data_type == 'exomes':
    #         non_cancer_ht = hl.read_table(annotations_ht_path(data_type, 'frequencies_tcga_with_consanguineous'))
    #         non_cancer_freq_meta = hl.eval(non_cancer_ht.globals.freq_meta)
    #         [x.update({'subset': 'non_cancer'}) for x in non_cancer_freq_meta]
    #         combined_ht = combined_ht.annotate(non_cancer=non_cancer_ht[combined_ht.key])
    #         combined_ht = combined_ht.annotate_globals(non_cancer_index_dict=make_index_dict(non_cancer_ht),
    #                                                    freq_meta=combined_ht.freq_meta.extend(non_cancer_freq_meta))
    #
    #     # Concatenate subset arrays into unified annotations
    #     subsets = ['controls', 'non_cancer', 'non_neuro', 'non_topmed'] if args.exomes else ['controls', 'non_neuro',
    #                                                                                          'non_topmed']
    #     array_lengths = get_array_lengths(combined_ht, subsets)
    #     combined_ht = combined_ht.select_globals('rf', 'freq_meta',
    #                                              freq_index_dict={k.replace("_adj", ""): v for k, v in
    #                                                               build_new_index_dict(combined_ht, subsets, array_lengths).items()},
    #                                              popmax_index_dict={v: i for i, v in
    #                                                                      enumerate(['gnomad'] + subsets)},
    #                                              age_index_dict={v: i for i, v in
    #                                                                 enumerate(['gnomad'] + subsets)},
    #                                              faf_index_dict={k.replace("_adj", ""): v for k, v in
    #                                                              build_faf_index_dict(combined_ht, subsets).items()})
    #     combined_ht = combined_ht.transmute(freq=concat_array_expr(combined_ht, subsets, 'freq'),
    #                            age_hist_hom=concat_struct_expr(combined_ht, subsets, 'age_hist_hom'),
    #                            age_hist_het=concat_struct_expr(combined_ht, subsets, 'age_hist_het'),
    #                            popmax=concat_struct_expr(combined_ht, subsets, 'popmax'),
    #                            faf=concat_array_expr(combined_ht, subsets, 'faf')
    #                            )
    #     combined_ht.write(release_ht_path(data_type), args.overwrite)

    if args.prepare_release_vcf:
        prepare_release_vcf(
            hl.read_table(release_ht_path(public=args.public)), # FIXME: Think whether this is the desired behavior
            subset_list=['gnomad'],
            public=args.public,
            pcsk9=args.pcsk9
        )

    if args.sanity_check_sites:
        # if data_type == 'exomes':
        #     subset_list = ['gnomad', 'controls', 'non_neuro', 'non_cancer', 'non_topmed'] if args.include_subset_frequencies else ['gnomad']
        # else:
        #     subset_list = ['gnomad', 'controls', 'non_neuro', 'non_topmed'] if args.include_subset_frequencies else ['gnomad']

        ht = hl.read_table(release_ht_path(public=False))
        sanity_check_ht(ht, data_type='v3_genomes', missingness_threshold=0.5, verbose=args.verbose)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pcsk9', help='Only do pcsk9.', action='store_true')
    parser.add_argument('--prepare_internal_ht', help='Prepare internal Hail Table', action='store_true')
    parser.add_argument('--model_id', help='Filtering model ID to use', default='vqsr_alleleSpecificTrans')
    parser.add_argument('--add_subset_frequencies', help='Add frequency annotations for gnomAD subsets', action='store_true')
    parser.add_argument('--include_subset_frequencies', help='Include frequency annotations for gnomAD subsets in release', action='store_true')
    parser.add_argument('--prepare_release_vcf', help='Prepare release VCF', action='store_true')
    parser.add_argument('--sanity_check_sites', help='Run sanity checks function', action='store_true')
    parser.add_argument('--public', help="When set, data is output to gnomad-public bucket", action='store_true')
    parser.add_argument('--verbose', help='Run sanity checks function with verbose output', action='store_true')
    parser.add_argument('--overwrite', help='Overwrite data', action='store_true')
    args = parser.parse_args()

    main(args)
