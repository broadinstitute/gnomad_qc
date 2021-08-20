from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import get_gnomad_data
import argparse
import sys
import logging
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("unified_sample_qc_a")
logger.setLevel(logging.INFO)


def annotate_sex(mt: hl.MatrixTable, out_internal_mt_prefix: str,
                 male_threshold: float = 0.8, female_threshold: float = 0.5) -> hl.MatrixTable:
    """
    Imputes sex, exports data, and annotates mt with this data
    NOTE: Evaluated in R (plots) and decided on cutoff of F<0.5 for females and F>0.8 for males (default) for genomes

    :param MatrixTable mt: MT containing samples to be ascertained for sex
    :param str out_internal_mt_prefix: file path prefix for tsv containing samples and sex imputation annotations
    :return: MatrixTable with imputed sex annotations stashed in column annotation 'sex_check'
    :rtype: MatrixTable
    """
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval('X')])
    sex_ht = hl.impute_sex(mt.GT, aaf_threshold=0.05, female_threshold=female_threshold, male_threshold=male_threshold)
    sex_ht.export(out_internal_mt_prefix + '.sex_check.txt.bgz')
    sex_colnames = ['f_stat', 'is_female']
    sex_ht = sex_ht.select(*sex_colnames)
    mt = mt.annotate_cols(**sex_ht[mt.col_key])
    return mt


def make_hard_filters_expr(ht: hl.Table, data_type: str) -> hl.expr.SetExpression:
    """
    NOTE: additional metadata in Kristen's import file is hard-coded

    :param: Table ht: input MT
    :param: str data_type: 'exomes' or 'genomes'
    :return: output MT
    :rtype: SetExpression
    """
    hard_filters = {
        'contamination': ht.freemix > 0.05,
        'callrate': ht.callrate < 0.85,
        'chimera': ht.pct_chimeras > 0.05,
        'ambiguous_sex': ht.ambiguous_sex
    }

    if data_type == 'exomes':
        hard_filters.update({
            'coverage': ht.mean_chr20_coverage == 0,
            'sex_aneuploidy': ht.sex_aneuploidy
        })
    else:
        hard_filters.update({
            'coverage': ht.mean_dp < 15,
            'insert_size': ht.median_insert_size < 250
        })
    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in hard_filters.items()]))


def make_perm_filters_expr(ht: hl.Table, data_type: str) -> hl.expr.SetExpression:
    """
    NOTE: syndip will remain dropped wrt to permissions, but all possible QC measures will still be calculated

    :param Table ht: input MT
    :param str data_type: 'exomes' or 'genomes'
    :return: output MT
    :rtype: SetExpression
    """
    if data_type == 'genomes':
        perm_filters = {
            'not_releasable': ~ht.releasable_2_1
        }
    else:
        perm_filters = {
            'tcga_tumor': ht.tcga_tumor,
            'tcga_barcode': ht.tcga_weird_barcode,
            'tcga_below_30': ht.tcga_below_30,
            'specific_exclusion': ht.specific_exclusion,
            'esp': ht.esp,
            'not_releasable': ht.non_releasable,
            'syndip': ht.syndip
        }
    return hl.set(hl.filter(lambda x: hl.is_defined(x),
                            [hl.or_missing(filter_expr, name) for name, filter_expr in perm_filters.items()]))


def main(args):
    hl.init()
    data_type = "genomes" if args.genomes else "exomes"

    if not args.skip_write_qc_mt:
        logger.info("Importing data...")
        # 1h40 for exomes, 3h20 for genomes
        mt = get_gnomad_data(data_type, raw=True, split=False)  # NOTE: using full calls since hardcalls doesn't exist at this stage
        logger.info("Filtering to bi-allelic, high-callrate, common SNPs for sample QC...")
        mt = mt.filter_rows((hl.len(mt.alleles) == 2) & hl.is_snp(mt.alleles[0], mt.alleles[1]) &
                            (hl.agg.mean(mt.GT.n_alt_alleles()) / 2 > 0.001) &
                            (hl.agg.fraction(hl.is_defined(mt.GT)) > 0.99))
        mt.annotate_cols(callrate=hl.agg.fraction(hl.is_defined(mt.GT))).naive_coalesce(5000).write(qc_mt_path(data_type), overwrite=args.overwrite)
    qc_mt = hl.read_matrix_table(qc_mt_path(data_type))

    logger.info("Importing metadata...")
    meta_ht = hl.import_table(qc_meta_path(data_type), impute=True, types={'age': hl.tfloat64}).key_by('s')
    qc_mt = qc_mt.annotate_cols(**meta_ht[qc_mt.s])

    logger.info("Inferring sex...")
    qc_ht = annotate_sex(qc_mt, qc_temp_data_prefix(data_type), male_threshold=0.8 if args.genomes else 0.6).cols()
    # Flag Klinefelter's individuals and samples with sex aneuploidies
    if args.exomes:
        qc_ht = qc_ht.annotate(ambiguous_sex=((qc_ht.f_stat >= 0.5) & (hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage <= 0.1))) |
                                             (hl.is_missing(qc_ht.f_stat)) |
                                             ((qc_ht.f_stat >= 0.4) & (qc_ht.f_stat <= 0.6) & (hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))),
                               sex_aneuploidy=(qc_ht.f_stat < 0.4) & hl.is_defined(qc_ht.normalized_y_coverage) & (qc_ht.normalized_y_coverage > 0.1))
    else:
        qc_ht = qc_ht.annotate(ambiguous_sex=hl.is_missing(qc_ht.is_female))

    logger.info("Annotating samples failing hard filters...")
    if args.exomes:
        sex_expr = (hl.case()
                    .when(qc_ht.ambiguous_sex, "ambiguous_sex")
                    .when(qc_ht.sex_aneuploidy, "sex_aneuploidy")
                    .when(qc_ht.is_female, "female")
                    .default("male"))
    else:
        sex_expr = (hl.case()
                    .when(qc_ht.ambiguous_sex, "ambiguous_sex")
                    .when(qc_ht.is_female, "female")
                    .default("male"))
    qc_ht = qc_ht.annotate(hard_filters=make_hard_filters_expr(qc_ht, data_type),
                           perm_filters=make_perm_filters_expr(qc_ht, data_type),
                           sex=sex_expr, data_type=data_type).key_by('data_type', 's')
    qc_ht.write(qc_ht_path(data_type, 'hard_filters'), overwrite=args.overwrite)

    # Export annotations to make rank list for relatedness (in final sample QC)
    if args.exomes:
        colnames = ['internal', 'project_id', 'pct_bases_20x', 'perm_filters']
    else:
        colnames = ['pcr_free', 'mean_dp', 'perm_filters']
    rank_ht = qc_ht.filter(hl.len(qc_ht.hard_filters) == 0, keep=True).select(*colnames)
    (rank_ht.annotate(releasable=(hl.len(rank_ht.perm_filters) == 0)).drop('perm_filters')
     .export(rank_annotations_path(data_type)))

    # Check numbers:
    qc_ht = hl.read_table(qc_ht_path(data_type, 'hard_filters'))
    sample_count = qc_ht.count()
    checkpoint1a = qc_ht.aggregate(hl.agg.count_where(hl.len(qc_ht['hard_filters']) == 0))
    checkpoint1b = qc_ht.aggregate(hl.agg.count_where((hl.len(qc_ht['hard_filters']) == 0) & (hl.len(qc_ht.perm_filters) == 0)))
    logger.info('{} samples found before filtering'.format(sample_count))
    logger.info('{} samples found after checkpoint 1a (hard filters)'.format(checkpoint1a))
    logger.info('{} samples found after checkpoint 1b (hard filters + permissions)'.format(checkpoint1b))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--exomes', help='Input MatrixTable contains exomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--genomes', help='Input MatrixTable contains genomes. One of --exomes or --genomes is required.',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--skip_write_qc_mt', help='Skip writing pre-calculated MatrixTable containing high-quality variants', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if int(args.exomes) + int(args.genomes) != 1:
        sys.exit('Error: One and only one of --exomes or --genomes must be specified')

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
