# noqa: D100

import argparse
import logging
from copy import deepcopy
from datetime import datetime
from functools import reduce
from typing import List

import hail as hl
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.utils.annotations import region_flag_expr
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import AS_FIELDS, SITE_FIELDS

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.annotations.insilico_predictors import get_sift_polyphen_from_vep
from gnomad_qc.v4.resources.annotations import (
    get_freq,
    get_info,
    get_insilico_predictors,
    get_vep,
    get_vrs,
)
from gnomad_qc.v4.resources.basics import qc_temp_prefix
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.release import FIELD_DESCRIPTIONS, release_sites
from gnomad_qc.v4.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)

# Remove InbreedingCoeff from allele-specific fields and SOR from SITE
# (processed separately from other fields)
AS_FIELDS = deepcopy(AS_FIELDS)
AS_FIELDS.remove("InbreedingCoeff")
SITE_FIELDS = deepcopy(SITE_FIELDS)
SITE_FIELDS.remove("SOR")
TABLES_FOR_RELEASE = [
    "dbsnp",
    # "filters",
    "freq",
    # "info",
    "region_flags",
    "in_silico",
    "vep",
]

INSILICO_PREDICTORS = ["spliceai", "pangolin", "revel"]  # "cadd", "phylop"


# Putting this in a function so that it is not evaluated until the script is run.
def get_config(release_exists: bool = False):
    """
    Get configuration dictionary.

    Format:
        '<Name of dataset>': {
                'path': 'gs://path/to/hailtable.ht',
                'select': '<Optional list of fields to select or dict of new field name to location of old fieldin the dataset.>',
                'field_name': '<Optional name of root annotation in combined dataset, defaults to name of dataset.>',
                'custom_select': '<Optional function name of custom select function that is needed for more advanced logic>',
                'select_globals': '<Optional list of globals to select or dict of new global field name to old global field name. If not specified, all globals are selected.>
            },

    :param release_exists: Whether the release HT already exists.
    :return: dict of datasets configs.
    """
    config = {
        "dbsnp": {
            "ht": dbsnp.ht(),
            "path": dbsnp.path,
            "select": ["rsid"],
        },
        # "filters": {
        #     "ht": final_filter().ht(),
        #     "path": final_filter().path,
        #     "select": ["filters", "vqsr"],
        #     "custom_select": custom_filters_select,
        #     "select_globals": ["filtering_model", "inbreeding_coeff_cutoff"],
        # },
        "in_silico": {
            "ht": reduce(
                (lambda joined_ht, ht: joined_ht.join(ht, "outer")),
                [
                    get_insilico_predictors(predictor=predictor).ht()
                    for predictor in INSILICO_PREDICTORS
                ],
            ),
            "path": [
                get_insilico_predictors(predictor=predictor).path
                # TODO: Add phylop maybe? see:
                # https://docs.google.com/document/d/1mHfhDJyY22JTFWoMhqfrDAIeKikTxYBec11v7tOiQeo/edit?disco=AAAA6g_ZPBc
                for predictor in INSILICO_PREDICTORS
            ],
            "field_name": "in_silico_predictors",
            "select": [
                # "cadd",
                "revel_max",
                "spliceai_ds_max",
                "pangolin_largest_ds",
                # "phylop"
            ],
            "custom_select": custom_in_silico_select,
            "select_globals": [
                # "cadd_version",
                "revel_version",
                "spliceai_version",
                "pangolin_version",
                # "phylop_version",
            ],
            "global_name": "tool_versions",
        },
        "info": {
            "ht": get_info().ht(),
            "path": get_info().path,
            "select": ["was_split", "a_index"],
            "custom_select": custom_info_select,
        },
        "freq": {
            "ht": get_freq().ht(),
            "path": get_freq().path,
            "select": [
                "freq",
                "faf",
                "grpmax",
                "gen_anc_faf_max",
                "histograms",
            ],
            "select_globals": [
                "freq_meta",
                "freq_index_dict",
                "freq_meta_sample_count",
                "faf_meta",
                "faf_index_dict",
                "age_distribution",
                "downsamplings",
            ],
        },
        "vep": {
            "ht": get_vep().ht(),
            # TODO: drop 100% missing? Can add to the custom vep select function
            "path": get_vep().path,
            "select": ["vep"],
            # Add in "custom_select" that drops insilicos
            "custom_select": custom_vep_select,
            # TODO: Update to have vep_csq_header -- should this be on the vep table
            # itself?
            "select_globals": ["vep_version", "vep_help", "vep_config"],
            "global_name": "vep_globals",
        },
        "region_flags": {
            "ht": get_freq().ht(),
            "path": get_freq().path,
            "custom_select": custom_region_flags_select,
        },
        "release": {
            "path": release_sites().path,
        },
        # TODO: Fill this in once we have the combined freq HT
        "combined_faf": {},
    }

    if release_exists:
        config["release"].update(
            {
                "ht": release_sites().ht(),
                "select": [r for r in release_sites().ht()._row],
                "select_globals": [g for g in release_sites().ht()._global],
            }
        )
    return config


def custom_in_silico_select(ht: hl.Table):
    """
    Get in silico predictors from VEP for release.

    :param ht: hail Table.
    :return: select expression dict
    """
    vep_in_silico = get_sift_polyphen_from_vep(get_vep().ht())
    selects = {
        "sift_max": vep_in_silico[ht.key].sift_max,
        "polyphen_max": vep_in_silico[ht.key].polyphen_max,
    }
    return selects


def custom_region_flags_select(ht: hl.Table):
    """
    Select region flags for release.

    :param ht: hail Table.
    :return: select expression dict
    """
    selects = {}
    selects = {
        "region_flags": region_flag_expr(
            ht,
            non_par=True,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
        )
    }
    return selects


def custom_filters_select(ht: hl.Table):
    """
    Select gnomad filter HT fields for release dataset.

    Extract fields like 'filters', 'vqsr', and generates 'alelle_info' struct.
    :param ht: hail table.
    :return: select expression dict.
    """
    selects = {}
    selects["allele_info"] = hl.struct(
        variant_type=ht.variant_type,
        allele_type=ht.allele_type,
        n_alt_alleles=ht.n_alt_alleles,
        was_mixed=ht.was_mixed,
    )
    return selects


# TODO: This function needs to be updated depending on which fields we want to keep from
#  info, I noticed there a dups within structs, e.g quasi_info has AS_SOR as does AS_info,
# maybe we need to pull these from filters instead?
def custom_info_select(ht: hl.Table):
    """
    Select fields for info hail Table annotation in release.

    The info field requires fields from the freq HT and the filters HT
    so those are pulled in here along with all info HT fields.

    :param ht: hail table
    :return: select expression dict
    """
    freq_ht = get_config().get("freq")["ht"]
    freq_info_dict = {"inbreeding_coeff": freq_ht[ht.key]["inbreeding_coeff"]}

    filters_ht = get_config().get("filters")["ht"]
    filters = filters_ht[ht.key]
    filters_info_fields = [
        "singleton",
        "transmitted_singleton",
        "omni",
        "mills",
        "monoallelic",
        "SOR",
    ]

    vrs_ht = get_vrs().ht()
    vrs_info_fields = {"vrs": vrs_ht[ht.key].vrs}

    filters_info_dict = {field: filters[field] for field in filters_info_fields}
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    filters_info_dict.update({**{f"{score_name}": filters[f"{score_name}"]}})

    info_dict = {field: ht.info[field] for field in SITE_FIELDS + AS_FIELDS}
    info_dict.update(filters_info_dict)
    info_dict.update(freq_info_dict)
    info_dict.update(vrs_info_fields)

    selects = {"info": hl.struct(**info_dict)}

    return selects


def custom_vep_select(ht: hl.Table):
    """
    Select fields for vep hail Table annotation in release.

    :param ht: hail table
    :return: select expression dict
    """
    selects = {
        "vep": ht.vep.annotate(
            transcript_consequences=ht.vep.transcript_consequences.map(
                lambda x: x.drop(
                    "sift_prediction",
                    "sift_score",
                    "polyphen_prediction",
                    "polyphen_score",
                )
            )
        )
    }
    return selects


def get_select_global_fields(ht: hl.Table):
    """
    Generate a dictionary of globals to select by checking the configs of all tables joined.

    :param ht: Final joined HT with globals.
    """
    t_globals = []
    for t in TABLES_FOR_RELEASE:
        config = get_config().get(t)
        if "select_globals" in config:
            select_globals = get_select_fields(config["select_globals"], ht)
            if "global_name" in config:
                global_name = config.get("global_name")
                select_globals = {global_name: hl.struct(**select_globals)}
            t_globals.append(select_globals)

    t_globals = reduce(lambda a, b: dict(a, **b), t_globals)

    return t_globals


def get_select_fields(selects, base_ht: hl.Table):
    """
    Generate a select dict from traversing the base_ht and extracting annotations.

    :param selects: mapping or list of selections
    :param base_ht: base_ht to traverse
    :return: select mapping from annotation name to base_ht annotation
    """
    select_fields = {}
    if selects is not None:
        if isinstance(selects, list):
            select_fields = {selection: base_ht[selection] for selection in selects}
        elif isinstance(selects, dict):
            for key, val in selects.items():
                # Grab the field and continually select it from the hail table.
                ht = base_ht
                for attr in val.split("."):
                    ht = ht[attr]
                select_fields[key] = ht
    return select_fields


def get_ht(dataset: str, _intervals, test, release_exists) -> hl.Table:
    """
    Return the appropriate hail table with selects applied.

    :param dataset: Hail Table to join.
    :param _intervals: Intervals for reading in hail Table.
    :param test: Whether call is for a test run.
    :return: Hail Table with fields to select.
    """
    logger.info("Getting the %s dataset and its selected annotations...", dataset)
    config = get_config(release_exists=release_exists)[dataset]

    # There is no single path for insilico so this impacts its join efficiency but
    # I dont think its worth the effort of restructuring this code to make it work with
    # multiple input paths thus the special treatment here
    if dataset == "in_silico":
        base_ht = config["ht"]
    else:
        ht_path = config["path"]
        logger.info("Reading in %s", dataset)
        base_ht = hl.read_table(ht_path, _intervals=_intervals)

    if test:
        base_ht = hl.filter_intervals(
            base_ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )

    select_fields = get_select_fields(config.get("select"), base_ht)

    if "custom_select" in config:
        custom_select_fn = config["custom_select"]
        select_fields = {**select_fields, **custom_select_fn(base_ht)}

    if "field_name" in config:
        field_name = config.get("field_name")
        select_query = {field_name: hl.struct(**select_fields)}
    else:
        select_query = select_fields

    return base_ht.select(**select_query)


def join_hts(
    base_table: hl.Table,
    tables: List[str],
    new_partition_percent: float,
    test: bool,
    release_exists: bool,
) -> hl.Table:
    """
    Outer join a list of hail tables.

    :param base_table: Dataset to use for interval partitioning.
    :param tables: List of tables to join.
    :param new_partition_percent: Percent of base_table partitions used for final release hail Table.
    :param test: Whether this is for a test run.
    :param release_exists: Whether the release HT already exists.
    :return: Hail Table with datasets joined.
    """
    logger.info(
        "Reading in %s to determine partition intervals for efficient join",
        base_table,
    )
    base_ht_path = get_config()[base_table]["path"]
    base_ht = hl.read_table(base_ht_path)
    if test:
        # Filter to PCSK9 for testing
        base_ht = hl.filter_intervals(
            base_ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )
    partition_intervals = base_ht._calculate_new_partitions(
        base_ht.n_partitions() * new_partition_percent
    )

    # Reorg list so base table is first
    tables.remove(base_table)
    tables.insert(0, base_table)

    logger.info("Joining datasets: %s...", tables)
    hts = [
        get_ht(
            table,
            _intervals=partition_intervals,
            test=test,
            release_exists=release_exists,
        )
        for table in tables
    ]
    # TODO: Check with hail if an intermediate checkpoint be helpful here
    joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "left")), hts)

    # Track the dataset we've added as well as the source path.
    included_dataset = {
        k: v["path"]
        for k, v in get_config(release_exists=release_exists).items()
        if k in tables
    }

    joined_ht = joined_ht.annotate_globals(
        date=datetime.now().isoformat(),
        datasets=included_dataset,
    )
    joined_ht.describe()
    return joined_ht


def main(args):
    """Create release ht."""
    hl.init(
        log="/create_release_ht.log",
        tmp_dir="gs://gnomad-tmp-4day",
        default_reference="GRCh38",
    )

    logger.info("Creating release HT...")
    ht = join_hts(
        args.base_table,
        args.tables_for_join,
        args.new_partition_percent,
        args.test,
        args.release_exists,
    )

    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    # ht = ht.filter(hl.is_defined(ht.filters))

    t_globals = get_select_global_fields(ht)

    # Previously discussed having a gnomad_qc and gnomad_method repo versions but this
    # does not make sense to me since they have been in continuous development and there
    # is not single version for either.
    ht = ht.select_globals(
        **t_globals,
        #    README=FIELD_DESCRIPTIONS,
        date=ht.date,
        datasets=ht.datasets,
        version=args.version,
    )

    # The dbsnp table does not have a global field for dbsnp_versions, same
    # with sift/polyphen and vrs (still need these)
    ht = ht.annotate_globals(
        tool_versions=ht.tool_versions.annotate(
            dbsnp_version="b156",
            sift_version="FILL THIS IS",
            polyphen_version="FILL THIS IS",
            vrs_version=hl.struct(
                **{
                    "vrs_python_version": "FILL THIS IS",
                    "seqrepo_version": "FILL THIS IS",
                },
            ),
        )
    )

    output_path = (
        qc_temp_prefix() + "release/gnomad.exomes.sites.test.ht"
        if args.test
        else release_sites().path
    )
    logger.info("Writing out release HT to %s", output_path)
    ht = ht.checkpoint(
        output_path,
        args.overwrite,
    )

    logger.info("Final variant count: %d", ht.count())
    ht.describe()


#    ht.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--new-partition-percent",
        help=(
            "Percent of start dataset partitions to use for release HT. Default is 1.1"
            " (110%)"
        ),
        default=1.1,
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "-v",
        "--version",
        help="The version of gnomAD.",
        default=CURRENT_RELEASE,
    )
    parser.add_argument(
        "-t",
        "--test",
        help="Runs a test on PCSK9 region, chr1:55039447-55064852",
        action="store_true",
    )
    parser.add_argument(
        "-j",
        "--tables-for-join",
        help="Tables to join for release",
        default=TABLES_FOR_RELEASE,
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "-b",
        "--base-table",
        help="Base table for interval partition calculation.",
        default="freq",
        choices=TABLES_FOR_RELEASE,
    )
    parser.add_argument(
        "--release-exists",
        help="Whether the release HT already exists.",
        action="store_true",
    )
    parser.add_argument(
        "--slack_channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
