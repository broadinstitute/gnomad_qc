from gnomad.resources.grch38.gnomad import (
    SUBSETS,
    GROUPS,
    SEXES,
    POPS,
    KG_POPS,
    HGDP_POPS,
    DOWNSAMPLINGS,
)
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.utils.annotations import region_flag_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import make_label_combos, index_globals
from gnomad.utils.vcf import (
    AS_FIELDS,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
)
from gnomad.utils.vep import VEP_CSQ_HEADER, VEP_CSQ_FIELDS

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq, allele_data
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources import (
    get_info,
    get_filtering_model,
    vep,
    qual_hist,
    release_ht_path,
)

import argparse
import logging
from typing import Union, Dict, List

import hail as hl

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)


AS_FIELDS.remove("InbreedingCoeff")
AS_FIELDS.extend(["AS_QUALapprox", "AS_SB_TABLE"])

SITE_FIELDS.remove("BaseQRankSum")
SITE_FIELDS.extend(["SB", "QUALapprox"])

POPS.extend(KG_POPS)
POPS.extend(HGDP_POPS)
POPS.extend(["global"])


def prepare_annotations(
    freq_ht: hl.Table, filtering_model_id: str, test: bool
) -> hl.Table:

    """
    Load and join all Tables with variant annotations.

    :param Table freq_ht: Table with frequency annotations.
    :param str filtering_model_id: ID for final filtering model chosen for release
    :param bool test: If true, filter table to small section of chr20
    :return: Table containing joined annotations.
    :rtype: hl.Table
    """

    if test:
        freq_ht = hl.filter_intervals(
            freq_ht, [hl.parse_locus_interval("chr20:1-1000000")]
        )

    logger.info("Loading annotation tables...")
    filters_ht = get_filtering_model(model_id=filtering_model_id).ht()
    vep_ht = vep.ht()
    dbsnp_ht = dbsnp.ht().select("rsid")
    info_ht = get_info().ht()

    logger.info("Filtering lowqual variants and assembling 'info' field...")
    info_fields = SITE_FIELDS + AS_FIELDS
    missing_info_fields = set(info_fields).difference(
        set([field for field in info_ht.take(1)[0].info])
    )
    select_info_fields = set(info_fields).intersection(
        set([field for field in info_ht.take(1)[0].info])
    )
    logger.info(f"Following fields not found in info HT: {missing_info_fields}")
    info_ht = info_ht.transmute(info=info_ht.info.select(*select_info_fields))
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            AS_SOR=filters_ht[info_ht.key].AS_SOR,
            SOR=filters_ht[info_ht.key].SOR,
            singleton=filters_ht[info_ht.key].singleton,
            transmitted_singleton=filters_ht[info_ht.key].transmitted_singleton,
            omni=filters_ht[info_ht.key].omni,
            mills=filters_ht[info_ht.key].mills,
            monoallelic=filters_ht[info_ht.key].monoallelic,
            **{f"{score_name}": filters_ht[info_ht.key][f"{score_name}"]},
        )
    )

    logger.info("Adding annotations...")
    vep_ht = vep_ht.transmute(vep=vep_ht.vep.drop("colocated_variants"))
    filters_ht = filters_ht.annotate(
        allele_data=hl.struct(
            has_star=filters_ht.has_star,  # NOTE: no v3.1 variants contained star alleles; consider dropping
            variant_type=filters_ht.variant_type,
            allele_type=filters_ht.allele_type,
            n_alt_alleles=filters_ht.n_alt_alleles,
            was_mixed=filters_ht.was_mixed,
        ),
    )

    ht = freq_ht.filter(info_ht[freq_ht.key].AS_lowqual, keep=False)
    ht = ht.annotate(
        a_index=info_ht[ht.key].a_index,
        was_split=info_ht[ht.key].was_split,
        rsid=dbsnp_ht[ht.key].rsid,
        filters=filters_ht[ht.key].filters,
        info=info_ht[ht.key].info,
        vep=vep_ht[ht.key].vep,
        vqsr=filters_ht[ht.key].vqsr,
        region_flag=region_flag_expr(
            ht,
            non_par=False,
            prob_regions={
                "lcr": lcr_intervals.ht(),
                "segdup": seg_dup_intervals.ht(),
            },
        ),
        allele_info=filters_ht[ht.key].allele_data,
    )
    ht = ht.annotate(info=ht.info.annotate(InbreedingCoeff=ht.InbreedingCoeff))
    ht = ht.drop("InbreedingCoeff")

    ht = ht.annotate_globals(
        vep_version=vep_ht.index_globals().version,
        vep_csq_header=VEP_CSQ_HEADER,
        dbsnp_version=dbsnp.default_version,
        filtering_model=filters_ht.index_globals().filtering_model,
    )

    return ht


def pre_process_subset(subset: str, global_ht: hl.Table, test: bool) -> hl.Table:
    """
    Prepares individual frequency tables for subsets in the release by adding identifying tags to subset freq_meta globals, filling in missing frequency fields for loci present only in the global cohort

    :param subset: Subset ID
    :param global_ht: Hail Table containing all variants discovered in the overall release cohort
    :param test:  If True, read smaller test Hail Table created for specified subset
    :return: Table containing subset frequencies with subset-tagged freq_meta global values and missing freq structs filled in
    """
    if test:
        ht = hl.read_table(f"gs://gnomad-tmp/gnomad_freq/chr20_test_freq.{subset}.ht")
    else:
        ht = get_freq(subset=subset).ht()

    ht = ht.annotate_globals(
        freq_meta=[{**x, **{"subset": subset}} for x in hl.eval(ht.freq_meta)]
    )

    # Keep all loci present in the global frequency HT table and fill in any missing frequency fields with an array of structs containing zero/null callstats data
    # This array should be as long as the freq_meta (each array entry corresponds to a freq_meta entry)
    new_ht = ht.join(global_ht.select().select_globals(), how="right")
    null_struct = hl.struct(
        AC=hl.int32(0),
        AF=hl.float64(0),
        AN=hl.null(hl.tint32),
        homozygote_count=hl.int32(0),
    )
    null_y_struct = hl.struct(
        AC=hl.null(hl.tint32),
        AF=hl.null(hl.tfloat64),
        AN=hl.null(hl.tint32),
        homozygote_count=hl.null(hl.int32),
    )
    null_array_struct = hl.map(
        lambda x: null_struct, hl.range(hl.len(hl.eval(ht.freq_meta)))
    )
    new_ht = new_ht.annotate(
        freq=hl.cond(hl.is_missing(new_ht.freq), null_array_struct, new_ht.freq)
    )

    return new_ht


def make_freq_index_dict(ht):
    """
    Create a look-up Dictionary for entries contained in the frequency annotation array
    :param Table ht: Table containing freq_meta global annotation to be indexed
    :return: Dictionary keyed by grouping combinations in the frequency array, with values describing the corresponding index
        of each grouping entry in the frequency array
    :rtype: Dict of str: int
    """
    freq_meta = hl.eval(ht.globals.freq_meta)

    index_dict = index_globals(freq_meta, dict(group=GROUPS))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)))
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES)))
    index_dict.update(
        index_globals(
            freq_meta, dict(downsampling=DOWNSAMPLINGS, group=["adj"], pop=POPS)
        )
    )
    index_dict.update(index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS)))
    index_dict.update(
        index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, pop=POPS))
    )
    index_dict.update(
        index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, sex=SEXES))
    )
    index_dict.update(
        index_globals(
            freq_meta, dict(group=GROUPS, subset=SUBSETS, pop=POPS, sex=SEXES)
        )
    )

    return index_dict


def main(args):

    hl.init(log="/prepare_internal_release.log", default_reference="GRCh38")
    test_freq_ht_path_concat = (
        "gs://gnomad-tmp/gnomad_freq/chr20_test_freq.concatenated.ht"
    )
    test_freq_ht_path_global = "gs://gnomad-tmp/gnomad_freq/chr20_test_freq.ht"

    if args.test and file_exists(test_freq_ht_path_concat):
        freq_ht = hl.read_table(test_freq_ht_path_concat)
        global_freq_ht = hl.read_table(test_freq_ht_path_global)

    else:
        logger.info("Concatenating subset frequencies...")
        if args.test:
            global_freq_ht = (
                hl.read_table(test_freq_ht_path_global)
                .select("freq")
                .select_globals("freq_meta")
            )
        else:
            global_freq_ht = (
                hl.read_table(get_freq().path)
                .select("freq")
                .select_globals("freq_meta")
            )

        subset_freq_hts = [
            pre_process_subset(subset, global_freq_ht, test=args.test)
            for subset in SUBSETS
        ]

        freq_ht = hl.Table.multi_way_zip_join(
            [global_freq_ht] + subset_freq_hts,
            data_field_name="freq",
            global_field_name="freq_meta",
        )
        freq_ht = freq_ht.transmute(freq=freq_ht.freq.flatmap(lambda x: x.freq))
        freq_ht = freq_ht.transmute_globals(
            freq_meta=freq_ht.freq_meta.flatmap(lambda x: x.freq_meta)
        )

        # Create frequency index dictionary on concatenated array (i.e., including all subsets)
        # NOTE: non-standard downsampling values are created in the frequency script corresponding to population totals, so new_downsampling values must be evaluated and included in the final freq_index_dict
        freq_meta = hl.eval(freq_ht.freq_meta)
        downsamplings = hl.filter(
            lambda x: x.keys().contains("downsampling"), freq_meta
        )
        new_downsamplings = hl.set(
            hl.map(lambda x: hl.int32(x["downsampling"]), downsamplings)
        )
        new_downsamplings_values = hl.eval(
            new_downsamplings.difference(hl.set(DOWNSAMPLINGS))
        )
        DOWNSAMPLINGS.extend([x for x in new_downsamplings_values])

        freq_ht = freq_ht.annotate_globals(
            freq_index_dict=make_freq_index_dict(freq_ht)
        )

        if args.test:
            freq_ht = freq_ht.checkpoint(test_freq_ht_path_concat)

        global_freq_ht = hl.read_table(get_freq().path)

    # Add back in all global frequency annotations not present in concatenated frequencies HT
    row_fields = set([field for field in global_freq_ht.row_value]).difference(
        set([field for field in freq_ht.row_value])
    )
    logger.info(
        f"Adding back the following row annotations onto concatenated frequencies: {row_fields}"
    )
    freq_ht = freq_ht.annotate(**global_freq_ht[freq_ht.key].select(*row_fields))

    global_fields = set([field for field in global_freq_ht.globals]).difference(
        set([field for field in freq_ht.globals])
    )
    global_fields.remove("downsamplings")
    logger.info(
        f"Adding back the following global annotations onto concatenated frequencies: {global_fields}"
    )
    freq_ht = freq_ht.annotate_globals(
        **global_freq_ht.index_globals().select(*global_fields)
    )

    logger.info("Preparing release Table annotations...")
    ht = prepare_annotations(freq_ht, args.model_id, test=False)

    logger.info("Removing chrM and sites without filter...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters), keep=True)

    # Add in new annotations for clinical variant interpretation
    in_silico_ht = hl.read_table(
        "gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1/gnomad_genomes_v3.1.analyst_annotations.ht"
    )  # TODO: update path with versioned Table resource path when available
    ht = ht.annotate(**in_silico_ht[ht.key])

    # Splice in fix to set female metrics to NA on Y chr
    female_idx = [
        hl.eval(ht.freq_index_dict.get(x))
        for x in hl.eval(ht.freq_index_dict)
        if "XX" in x
    ]
    freq_idx_range = hl.eval(hl.range(hl.eval(hl.len(hl.eval(ht.freq_meta)))))
    null_y_struct = hl.struct(
        AC=hl.null(hl.tint32),
        AF=hl.null(hl.tfloat64),
        AN=hl.null(hl.tint32),
        homozygote_count=hl.null(hl.tint32),
    )
    ht = ht.annotate(
        freq=hl.if_else(
            (ht.locus.in_y_nonpar() | ht.locus.in_y_par()),
            hl.map(
                lambda x: hl.if_else(
                    hl.array(female_idx).contains(x), null_y_struct, ht.freq[x]
                ),
                freq_idx_range,
            ),
            ht.freq,
        )
    )

    ht = ht.checkpoint(
        "gs://gnomad-tmp/release/v3.1/gnomad.genomes.v3.1.sites.ht"
        if args.test
        else release_ht_path(public=False),
        args.overwrite,
    )
    logger.info(f"Final variant count: {ht.count()}")
    ht.describe()
    ht.show()
    ht.summarize()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test", help="Runs a test on chr20:1-1000000", action="store_true"
    )
    parser.add_argument(
        "--model_id",
        help="Filtering model ID to use",
        default="vqsr_alleleSpecificTrans",
    )
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
