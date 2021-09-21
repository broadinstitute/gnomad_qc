import argparse
import logging

import hail as hl

from gnomad.resources.grch38.gnomad import (
    POPS,
    POPS_STORED_AS_SUBPOPS,
    SUBSETS,
)
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import missing_callstats_expr, region_flag_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.release import make_freq_index_dict
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import (
    AS_FIELDS,
    SITE_FIELDS,
)
from gnomad.utils.vep import VEP_CSQ_HEADER

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import (
    analyst_annotations,
    get_freq,
    get_info,
    vep,
)
from gnomad_qc.v3.resources.basics import get_checkpoint_path, qc_temp_prefix
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)

# Remove InbreedingCoeff from allele-specific fields (processed separately from other fields)
AS_FIELDS.remove("InbreedingCoeff")

# Add fine-resolution populations specific to 1KG and HGDP to standard gnomAD pops; used to create frequency index dictionary
POPS.extend(POPS_STORED_AS_SUBPOPS)
# Add 'global' tag used to distinguish cohort-wide vs. subset annotations in frequency index dictionary
POPS.extend(["global"])


def add_release_annotations(freq_ht: hl.Table) -> hl.Table:
    """
    Load and join all Tables with variant annotations.

    :param freq_ht: Table with frequency annotations
    :return: Table containing joined annotations
    """
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
        "The following fields are not found in the info HT: %s", missing_info_fields
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
        region_flag=region_flag_expr(
            ht,
            non_par=False,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
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
        inbreeding_coeff_cutoff=filters_ht.index_globals().inbreeding_coeff_cutoff,
    )

    return ht


def pre_process_subset_freq(
    subset: str, global_ht: hl.Table, test: bool = False
) -> hl.Table:
    """
    Prepare subset frequency Table by filling in missing frequency fields for loci present only in the global cohort.

    .. note::

        The resulting final `freq` array will be as long as the subset `freq_meta` global (i.e., one `freq` entry for each `freq_meta` entry)

    :param subset: subset ID
    :param global_ht: Hail Table containing all variants discovered in the overall release cohort
    :param test: If True, filter to small region on chr20
    :return: Table containing subset frequencies with missing freq structs filled in
    """

    # Read in subset HTs
    subset_ht_path = f"gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1.2/gnomad_genomes_v3.1.2.patch.frequencies.{subset}.ht"
    subset_chr20_ht_path = get_checkpoint_path(
        f"chr20_test_gnomad_genomes_v3.1.2.patch.frequencies.{subset}"
    )

    if test:
        if file_exists(subset_chr20_ht_path):
            logger.info(
                "Loading chr20 %s subset frequency data for testing: %s",
                subset,
                subset_chr20_ht_path,
            )
            subset_ht = hl.read_table(subset_chr20_ht_path)

        elif file_exists(subset_ht_path):
            logger.info(
                "Loading %s subset frequency data for testing: %s",
                subset,
                subset_ht_path,
            )
            subset_ht = hl.read_table(subset_ht_path)
            subset_ht = hl.filter_intervals(
                subset_ht, [hl.parse_locus_interval("chr20:1-1000000")]
            )

    elif file_exists(subset_ht_path):
        logger.info("Loading %s subset frequency data: %s", subset, subset_ht_path)
        subset_ht = hl.read_table(subset_ht_path)

    else:
        raise DataException(
            f"Hail Table containing {subset} subset frequencies not found. You may need to run the script generate_freq_data.py to generate frequency annotations first."
        )

    # Fill in missing freq structs
    ht = subset_ht.join(global_ht.select().select_globals(), how="right")
    ht = ht.annotate(
        freq=hl.if_else(
            hl.is_missing(ht.freq),
            hl.map(lambda x: missing_callstats_expr(), hl.range(hl.len(ht.freq_meta))),
            ht.freq,
        )
    )
    ht = ht.select("freq").select_globals("freq_meta")

    return ht


def main(args):
    hl.init(log="/create_release_ht.log", default_reference="GRCh38")
    global_freq_path = "gs://gnomad/annotations/hail-0.2/ht/genomes_v3.1.2/gnomad_genomes_v3.1.2.patch.frequencies.ht"

    # The concatenated HT contains all subset frequency annotations, plus the overall cohort frequency annotations,
    # concatenated together in a single freq annotation ('freq')

    # Load global frequency Table
    if args.test:
        global_freq_chr20_ht_path = get_checkpoint_path(
            "chr20_test_gnomad_genomes_v3.1.2.patch.frequencies"
        )

        if file_exists(global_freq_chr20_ht_path):
            logger.info(
                "Loading chr20 global frequency data for testing: %s",
                global_freq_chr20_ht_path,
            )
            global_freq_ht = (
                hl.read_table(global_freq_chr20_ht_path)
                .select("freq")
                .select_globals("freq_meta")
            )

        elif file_exists(global_freq_path):
            logger.info(
                "Loading global frequency data for testing: %s", global_freq_path
            )
            global_freq_ht = (
                hl.read_table(global_freq_path)
                .select("freq")
                .select_globals("freq_meta")
            )
            global_freq_ht = hl.filter_intervals(
                global_freq_ht, [hl.parse_locus_interval("chr20:1-1000000")]
            )

    elif file_exists(global_freq_path):
        logger.info("Loading global frequency data: %s", global_freq_path)
        global_freq_ht = (
            hl.read_table(global_freq_path).select("freq").select_globals("freq_meta")
        )

    else:
        raise DataException(
            "Hail Table containing global callset frequencies not found. You may need to run the script to generate frequency annotations first."
        )

    # Load subset frequency Table(s)
    if args.test:
        test_subsets = args.test_subsets
        subset_freq_hts = [
            pre_process_subset_freq(subset, global_freq_ht, test=True)
            for subset in test_subsets
        ]

    else:
        subset_freq_hts = [
            pre_process_subset_freq(subset, global_freq_ht) for subset in SUBSETS
        ]

    logger.info("Concatenating subset frequencies...")
    freq_ht = hl.Table.multi_way_zip_join(
        [global_freq_ht.select("freq").select_globals("freq_meta")] + subset_freq_hts,
        data_field_name="freq",
        global_field_name="freq_meta",
    )
    freq_ht = freq_ht.transmute(freq=freq_ht.freq.flatmap(lambda x: x.freq))
    freq_ht = freq_ht.transmute_globals(
        freq_meta=freq_ht.freq_meta.flatmap(lambda x: x.freq_meta)
    )

    # Create frequency index dictionary on concatenated array (i.e., including all subsets)
    # NOTE: non-standard downsampling values are created in the frequency script corresponding to population totals, so
    # callset-specific DOWNSAMPLINGS must be used instead of the generic DOWNSAMPLING values
    global_freq_ht = hl.read_table(get_freq().path)
    freq_ht = freq_ht.annotate_globals(
        freq_index_dict=make_freq_index_dict(
            freq_meta=hl.eval(freq_ht.freq_meta),
            pops=POPS,
            downsamplings=hl.eval(global_freq_ht.downsamplings),
        )
    )

    # Add back in all global frequency annotations not present in concatenated frequencies HT
    row_fields = global_freq_ht.row_value.keys() - freq_ht.row_value.keys()
    logger.info(
        "Adding back the following row annotations onto concatenated frequencies: %s",
        row_fields,
    )
    freq_ht = freq_ht.annotate(**global_freq_ht[freq_ht.key].select(*row_fields))

    global_fields = global_freq_ht.globals.keys() - freq_ht.globals.keys()
    global_fields.remove("downsamplings")
    logger.info(
        "Adding back the following global annotations onto concatenated frequencies: %s",
        global_fields,
    )
    freq_ht = freq_ht.annotate_globals(
        **global_freq_ht.index_globals().select(*global_fields)
    )

    logger.info("Preparing release Table annotations...")
    ht = add_release_annotations(freq_ht)

    logger.info("Removing chrM and sites without filter...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

    ht = ht.checkpoint(
        qc_temp_prefix(version="3.1.2")
        + "release/gnomad.genomes.v3.1.sites.patch.chr20.ht"
        if args.test
        else "gs://gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.1.sites.patch.ht",
        args.overwrite,
    )

    logger.info("Applying v3.1.2 patch to v3.1.1 release HT...")

    # Using v3.1 allele frequency HT
    v3_1_1_release_ht = release_sites().versions["3.1.1"].ht()
    if args.test:
        logger.info("Filtering to two partitions on chr20")
        v3_1_1_release_ht = hl.filter_intervals(
            v3_1_1_release_ht, [hl.parse_locus_interval("chr20:1-1000000")]
        )
        v3_1_1_release_ht = v3_1_1_release_ht._filter_partitions(range(2))

    ht = v3_1_1_release_ht.transmute(
        **{
            x: hl.coalesce(ht[v3_1_1_release_ht.key][x], v3_1_1_release_ht[x])
            for x in [
                "freq",
                "faf",
                "popmax",
                "qual_hists",
                "raw_qual_hists",
            ]
        },
        info=v3_1_1_release_ht.info.annotate(
            InbreedingCoeff=hl.coalesce(
                ht[v3_1_1_release_ht.key].info.InbreedingCoeff,
                v3_1_1_release_ht.info.InbreedingCoeff,
            )
        ),
    )

    logger.info("Writing out v3.1.2 release HT...")
    ht = ht.checkpoint(
        qc_temp_prefix(version="3.1.2") + "release/gnomad.genomes.v3.1.sites.chr20.ht"
        if args.test
        else "gs://gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.1.sites.ht",
        args.overwrite,
    )

    logger.info("Final variant count: %d", ht.count())
    ht.describe()
    ht.show()
    ht.summarize()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--test", help="Runs a test on chr20:1-1000000", action="store_true"
    )
    parser.add_argument(
        "--test_subsets",
        help="Specify subsets on which to run test, e.g. '--test_subsets non_v2 non_topmed'",
        default=SUBSETS,
        nargs="+",
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
