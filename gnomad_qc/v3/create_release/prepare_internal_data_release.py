import argparse
import logging
from typing import Dict, List, Union

import hail as hl

from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    GROUPS,
    HGDP_POPS,
    KG_POPS,
    POPS,
    SEXES,
    SUBSETS,
)
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import null_callstats_expr, region_flag_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    AS_FIELDS,
    index_globals,
    make_label_combos,
    RF_FIELDS,
    SITE_FIELDS,
    VQSR_FIELDS,
)
from gnomad.utils.vep import VEP_CSQ_HEADER, VEP_CSQ_FIELDS
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import (
    allele_data,
    analyst_annotations,
    get_freq,
    get_info,
    qual_hist,
    vep,
)
from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.release import release_ht_path
from gnomad_qc.v3.resources.variant_qc import final_filter


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("data_release")
logger.setLevel(logging.INFO)

# Remove InbreedingCoeff from allele-specific fields (processed separately from other fields)
AS_FIELDS.remove("InbreedingCoeff")
# Remove BaseQRankSum, as we are keeping the allele-specific version of this annotation instead
SITE_FIELDS.remove("BaseQRankSum")

POPS.extend(KG_POPS)
POPS.extend(HGDP_POPS)
POPS.extend(["global"])


def add_release_annotations(
    freq_ht: hl.Table, test: bool
) -> hl.Table:

    """
    Loads and joins all Tables with variant annotations.

    :param Table freq_ht: Table with frequency annotations
    :param bool test: If true, filter table to small section of chr20
    :return: Table containing joined annotations
    :rtype: hl.Table
    """

    if test:
        freq_ht = hl.filter_intervals(
            freq_ht, [hl.parse_locus_interval("chr20:1-1000000")]
        )

    logger.info("Loading annotation tables...")
    filters_ht = final_filter.ht()
    vep_ht = vep.ht()
    dbsnp_ht = dbsnp.ht().select("rsid")
    info_ht = get_info().ht()
    in_silico_ht = analyst_annotations.ht()

    logger.info("Filtering lowqual variants and assembling 'info' field...")
    info_fields = SITE_FIELDS + AS_FIELDS
    missing_info_fields = set(info_fields).difference(info_ht.info.keys())
    logger.info(
        f"The following fields are not found in the info HT: {missing_info_fields}"
    )

    select_info_fields = set(info_fields).intersection(info_ht.info.keys())
    info_ht = info_ht.transmute(info=info_ht.info.select(*select_info_fields))
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    filters = filters_ht[info_ht.key]
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            AS_SOR=filters.AS_SOR,  # NOTE: AS_SOR will be incorporated into the info HT after v3.1, so no need to add this annotation in future releases
            SOR=filters.SOR,
            singleton=filters.singleton,
            transmitted_singleton=filters.transmitted_singleton,
            omni=filters.omni,
            mills=filters.mills,
            monoallelic=filters.monoallelic,
            **{f"{score_name}": filters[f"{score_name}"]},
        )
    )

    logger.info("Adding annotations...")
    filters_ht = filters_ht.select(
        "filters",
        "vqsr",
        allele_info=hl.struct(
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
        info=info_ht[ht.key].info,
        vep=vep_ht[ht.key].vep.drop("colocated_variants"),
        vqsr=filters_ht[ht.key].vqsr,
        region_flag=region_flag_expr(
            ht,
            non_par=False,
            prob_regions={
                "lcr": lcr_intervals.ht(),
                "segdup": seg_dup_intervals.ht(),
            },
        ),
        **filters_ht[ht.key],
        **in_silico_ht[ht.key],
    )
    ht = ht.transmute(info=ht.info.annotate(InbreedingCoeff=ht.InbreedingCoeff))

    ht = ht.annotate_globals(
        vep_version=vep_ht.index_globals().version,
        vep_csq_header=VEP_CSQ_HEADER,
        dbsnp_version=dbsnp.default_version,
        filtering_model=filters_ht.index_globals().filtering_model,
    )

    return ht


def pre_process_subset_freq(subset: str, global_ht: hl.Table, test: bool) -> hl.Table:
    """
    Prepares individual frequency tables for subsets in the release.

    This function will:
    - Add identifying tags to the subset freq_meta globals
    - Fill in missing frequency fields for loci present only in the global cohort

    :param subset: Subset ID
    :param global_ht: Hail Table containing all variants discovered in the overall release cohort
    :param test: If True, read smaller test Hail Table created for specified subset
    :return: Table containing subset frequencies with subset-tagged freq_meta global values and missing freq structs filled in
    """
    if test:
        ht = hl.read_table(f"gs://gnomad-tmp/gnomad_freq/chr20_test_freq.{subset}.ht")
    else:
        ht = get_freq(subset=subset).ht()

    # Keep all loci present in the global frequency HT table and fill in any missing frequency fields with an array of structs containing zero/null callstats data
    # This array should be as long as the freq_meta (each array entry corresponds to a freq_meta entry)
    new_ht = ht.join(global_ht.select().select_globals(), how="right")
    new_ht = new_ht.annotate(
        freq=hl.if_else(
            hl.is_missing(new_ht.freq),
            hl.map(
                lambda x: null_callstats_expr(), hl.range(hl.len(ht.freq_meta))
            ),
            new_ht.freq,
        )
    )

    return new_ht


def make_freq_index_dict(freq_meta: List[Dict[str, str]]) -> Dict[str, int]:
    """
    Creates a look-up Dictionary for entries contained in the frequency annotation array.

    :param List of Dict freq_meta: Global annotation continaing the set of groupings for each element of the freq array (e.g., [{'group': 'adj'}, {'group': 'adj', 'pop': 'nfe'}])
    :return: Dictionary keyed by the grouping combinations found in the frequency array, where values are the corresponding 0-based indices for the groupings in the freq_meta array
    :rtype: Dict of str: int
    """

    index_dict = {
        **index_globals(freq_meta, dict(group=GROUPS)),
        **index_globals(freq_meta, dict(group=GROUPS, pop=POPS)),
        **index_globals(freq_meta, dict(group=GROUPS, sex=SEXES)),
        **index_globals(freq_meta, dict(group=GROUPS, pop=POPS, sex=SEXES)),
        **index_globals(freq_meta, dict(downsampling=DOWNSAMPLINGS, group=["adj"], pop=POPS)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, pop=POPS)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, sex=SEXES)),
        **index_globals(freq_meta, dict(group=GROUPS, subset=SUBSETS, pop=POPS, sex=SEXES)),
    }

    return index_dict


def main(args):

    hl.init(log="/prepare_internal_release.log", default_reference="GRCh38")

    # The concatenated HT contains all subset frequency annotations, plus the overall cohort frequency annotations,
    # concatenated together in a single freq annotation ('freq')

    if args.test:
        test_freq_ht_path_concat = (
            "gs://gnomad-tmp/gnomad_freq/chr20_test_freq.concatenated.ht"
        )
        test_freq_ht_path_global = "gs://gnomad-tmp/gnomad_freq/chr20_test_freq.ht"

        if file_exists(test_freq_ht_path_concat):
            freq_ht = hl.read_table(test_freq_ht_path_concat)

        else:
            logger.info("Concatenating subset frequencies on test data...")
            if file_exists(test_freq_ht_path_global):
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
                global_freq_ht = hl.filter_intervals(
                    global_freq_ht, [hl.parse_locus_interval("chr20:1-1000000")]
                )
                global_freq_ht = global_freq_ht.checkpoint(test_freq_ht_path_global)

            # TODO: Create concatenated subset HT
            freq_ht = freq_ht.checkpoint(test_freq_ht_path_concat)

    else:
        logger.info("Concatenating subset frequencies...")
        global_freq_ht = (
            hl.read_table(get_freq().path).select("freq").select_globals("freq_meta")
        )

        subset_freq_hts = [
            pre_process_subset_freq(subset, global_freq_ht, test=args.test)
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
            freq_index_dict=make_freq_index_dict(freq_meta)
        )

        global_freq_ht = hl.read_table(get_freq().path)

    # Add back in all global frequency annotations not present in concatenated frequencies HT
    row_fields = global_freq_ht.row_value.keys() - freq_ht.row_value.keys()
    logger.info(
        f"Adding back the following row annotations onto concatenated frequencies: {row_fields}"
    )
    freq_ht = freq_ht.annotate(**global_freq_ht[freq_ht.key].select(*row_fields))

    global_fields = global_freq_ht.globals.keys() - freq_ht.globals.keys()
    global_fields.remove("downsamplings")
    logger.info(
        f"Adding back the following global annotations onto concatenated frequencies: {global_fields}"
    )
    freq_ht = freq_ht.annotate_globals(
        **global_freq_ht.index_globals().select(*global_fields)
    )

    logger.info("Preparing release Table annotations...")
    ht = add_release_annotations(freq_ht, test=False)

    logger.info("Removing chrM and sites without filter...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

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
