import argparse
import logging
from typing import Union, Dict

import hail as hl

from gnomad.resources.grch38.reference_data import dbsnp, lcr_intervals, seg_dup_intervals
from gnomad.resources.grch38.gnomad import SUBSETS  # FIXME -- not importing
from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq, allele_data
from gnomad_qc.v3.resources import get_info, get_filtering_model, vep, qual_hist, release_ht_path, get_vqsr_filters # freq, V3_RELEASE_VERSION,
from gnomad.utils.file_utils import file_exists


from gnomad.utils.vcf import (
    AS_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
)


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)

# Rename RF probability in RF fields  # TODO: still applies?
RF_FIELDS.append("rf_tp_probability")

# Remove InbreedingCoeff from allele-specific fields (processed separately from other fields)
AS_FIELDS.remove("InbreedingCoeff")
AS_FIELDS.extend(["AS_QUALapprox", "AS_SB_TABLE"])

SITE_FIELDS.remove("BaseQRankSum")
SITE_FIELDS.extend(["SB", "QUALapprox"])


def region_flag_expr(
        t: Union[hl.Table, hl.MatrixTable],
        non_par: bool = True,
        prob_regions: Dict[str, hl.Table] = None
) -> hl.expr.StructExpression:
    """
    Creates `region_flag` struct.
    Struct contains flags for problematic regions (i.e., LCR, decoy, segdup, and nonpar regions).

    .. note::
        No hg38 resources for decoy or self chain available yet.
    :param Table/MatrixTable t: Input Table/MatrixTable.
    :return: `region_flag` struct row annotation.
    :rtype: hl.expr.StructExpression
    """

    prob_flags_expr = {'non_par': (t.locus.in_x_nonpar() | t.locus.in_y_nonpar())} if non_par else {}

    if prob_regions is not None:
        prob_flags_expr.update({
            region_name: hl.is_defined(region_table[t.locus])
            for region_name, region_table in prob_regions.items()
        })

    return hl.struct(**prob_flags_expr)


VQSR_FIELDS = ["VQSLOD", "culprit", "NEGATIVE_TRAIN_SITE", "POSITIVE_TRAIN_SITE"]  # TODO: remove for v3.1


def prepare_annotations(
    freq_ht: hl.Table, filtering_model_id: str, pcsk9: bool) -> hl.Table:

    """
    Load and join all Tables with variant annotations.

    :param Table freq_ht: Table with frequency annotations.
    :param str filtering_model_id:
    :return: Table containing joined annotations.
    :rtype: hl.Table
    """

    if pcsk9:
        freq_ht = hl.filter_intervals(freq_ht, [hl.parse_locus_interval('chr1:55030000-55070000')])

    logger.info("Loading annotation tables...")
    filters_ht = get_filtering_model(model_id=filtering_model_id).ht()  # TODO: replace with rf_ht?
    vep_ht = vep.versions["3"].ht()  # TODO: remove v3 tag
    dbsnp_ht = dbsnp.ht().select("rsid")  # TODO: update resource
    info_ht = get_info().ht()
    # allele_ht = allele_data.ht()  # TODO: update resource
    vqsr_ht = get_vqsr_filters(model_id=filtering_model_id).versions["3"].ht()  # TODO: remove v3 tag

    logger.info("Filtering lowqual variants and assembling 'info' field...")
    info_fields = SITE_FIELDS + AS_FIELDS
    missing_info_fields = set(info_fields).difference(
        set([field for field in info_ht.take(1)[0].info]))
    select_info_fields = set(info_fields).intersection(set([field for field in info_ht.take(1)[0].info]))
    logger.info(f'Following fields not found in info HT: {missing_info_fields}')
    info_ht = info_ht.transmute(info=info_ht.info.select(*select_info_fields))

    logger.info("Adding annotations...")
    # rf_ht = rf_ht.select(*RF_FIELDS)
    vep_ht = vep_ht.transmute(vep=vep_ht.vep.drop("colocated_variants"))
    vqsr_ht = vqsr_ht.transmute(info=vqsr_ht.info.select(*VQSR_FIELDS)).select("info")

    ht = freq_ht.filter(info_ht[freq_ht.key].AS_lowqual, keep=False)
    ht = ht.annotate(
        a_index=info_ht[ht.key].a_index,
        was_split=info_ht[ht.key].was_split,
        rsid=dbsnp_ht[ht.key].rsid,
        filters=filters_ht[freq_ht.key].filters,  #Replace with: filters=rf_ht[ht.key].filters;
        info=info_ht[ht.key].info,
        vep=vep_ht[ht.key].vep,
        vqsr=vqsr_ht[ht.key].info,
        region_flag=region_flag_expr(ht, non_par=False,
                                     prob_regions={'lcr': lcr_intervals.ht(),  # TODO: Make this dictionary a variable?
                                                   'segdup': seg_dup_intervals.ht()}),
        # allele_info=allele_ht[ht.key].allele_data
        # rf=rf_ht[ht.key]
    )
    ht = ht.transmute(
        info=ht.info.annotate(InbreedingCoeff=ht.InbreedingCoeff)
    )

    # Adjust qual hist formatting -- TODO delete when updated in freq code
    qual_hists = [field for field in ht.row_value if "_hist_" in field]
    adj_hists = [field for field in qual_hists if "_adj" in field]
    raw_hists = [field for field in qual_hists if "_adj" not in field]
    ht = ht.annotate(
        qual_hists=hl.Struct(**{i.replace("_adj", ""): ht[i] for i in adj_hists}),
        raw_qual_hists=hl.Struct(**{i: ht[i] for i in raw_hists}),
    )
    ht = ht.drop(*qual_hists)


    # TODO: check of there are infinite QD values -- may need to set these to missing:
    # ht = ht.annotate(
    #     info=ht.info.annotate(
    #         QD=hl.or_missing(hl.is_finite(ht.info.QD), ht.info.QD) # TODO: Check this upstream
    #     )
    # )
    #NOTE: Laurent annotated with add_variant_type ('variant_type', 'n_alt_alleles'); variant histograms

    # ht = ht.annotate_globals(freq_index_dict=make_index_dict(ht), faf_index_dict=make_faf_index_dict(ht), filtering_model=filter_ht.index_globals().filtering_model, rf_globals=rf_ht.index_globals())
    # TODO: Also add age_distribution, age_index_dict, vep_csq_header

    return ht

# TODO: flag non-par regions?

def pre_process_subset(subset: str, test: bool) -> hl.Table:
    if test:
        ht = hl.read_table(f"gs://gnomad-tmp/gnomad_freq/chr20_1_1000000_freq.{subset}.ht")
    else:
        ht = get_freq(subset=subset).ht()

    ht = ht.annotate_globals(freq_meta=[{**x, **{'subset': subset}} for x in hl.eval(ht.freq_meta)])
    return ht


def main(args):

    hl.init(log="/prepare_internal_release.log", default_reference="GRCh38")
    test_freq_ht_path = "gs://gnomad-tmp/gnomad_freq/chr20_1_1000000_freq.concatenated.ht"

    if args.test and file_exists(test_freq_ht_path):
        freq_ht = hl.read_table(test_freq_ht_path)
        global_freq_ht = hl.read_table("gs://gnomad-tmp/gnomad_freq/chr20_1_1000000_freq.ht")

    else:
        logger.info("Concatenating subset frequencies...")
        if args.test:
            global_freq_ht = hl.read_table("gs://gnomad-tmp/gnomad_freq/chr20_1_1000000_freq.ht").select('freq').select_globals('freq_meta')
            SUBSETS = ['v2_non_neuro', 'v2_control', 'v2_non_topmed']
        else:
            global_freq_ht = get_freq.ht().select('freq').select_globals('freq_meta')

        subset_freq_hts = [pre_process_subset(subset, test=args.test) for subset in SUBSETS]

        #TODO: write check for even length of all subset + freq Tables and issue warning if not

        freq_ht = hl.Table.multi_way_zip_join([global_freq_ht] + subset_freq_hts, data_field_name='freq', global_field_name='freq_meta')
        freq_ht = freq_ht.transmute(freq=freq_ht.freq.flatmap(lambda x: x.freq))
        freq_ht = freq_ht.transmute_globals(freq_meta=freq_ht.freq_meta.flatmap(lambda x: x.freq_meta))

        # TODO: Create index dict on concatenated array

        if args.test:
            freq_ht = freq_ht.checkpoint(test_freq_ht_path)

    # Add back in all global frequency annotations not present in concatenated frequencies HT
    row_fields = set([field for field in global_freq_ht.row_value]).difference(set([field for field in freq_ht.row_value]))
    logger.info(f"Adding back the following row annotations onto concatenated frequencies: {row_fields}")
    freq_ht = freq_ht.annotate(**global_freq_ht[freq_ht.key].select(*row_fields))

    global_fields = set([field for field in global_freq_ht.globals]).difference(set([field for field in freq_ht.globals]))
    logger.info(f"Adding back the following global annotations onto concatenated frequencies: {global_fields}")
    freq_ht = freq_ht.annotate_globals(**global_freq_ht.index_globals().select(*global_fields))
    # TODO: drop 'downsamplings' global

    ht = prepare_annotations(freq_ht, args.model_id, pcsk9=False)  # TODO: change test to PCSK9
    logger.info("Removing chrM...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht.describe()
    ht.show()

    # logger.info("Getting age hist data...")
    # age_hist_data = get_age_distributions(get_age_ht(freeze))
    # ht = ht.annotate_globals(age_distribution=age_hist_data)
    # ht = ht.naive_coalesce(args.n_partitions)

    # TODO: add in new variant annotations from William

    ht = ht.checkpoint("gs://gnomad-tmp/release/v3.1/gnomad.genomes.v3.1.sites.ht" if args.test else release_ht_path(public=False), args.overwrite)
    logger.info(f"Final variant count: {ht.count()}")

# TODO: Check that variant QC root path is updated for v3.1
# TODO: Add 'SOR' (AS_SOR?) to `info` struct from VQSR HT
# TODO: Add omni, mills, transmitted_singleton, singleton to `info` struct

if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('--test', help='Runs a test on chr20:1-1000000', action='store_true')
    parser.add_argument('--model_id', help='Filtering model ID to use', default='vqsr_alleleSpecificTrans')

    parser.add_argument(
        "--n_partitions",
        help="Number of desired partitions for output Tables",
        default=10000,
        type=int,
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)

