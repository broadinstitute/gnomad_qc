# noqa: D100

import argparse
import logging
from datetime import datetime
from functools import reduce

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: Update gnomAD_methods to dict: key version, value subsets
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
from gnomad.utils.vcf import AS_FIELDS, SITE_FIELDS
from gnomad.utils.vep import VEP_CSQ_HEADER

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import (
    analyst_annotations,
    get_freq,
    get_info,
    vep,
)
from gnomad_qc.v3.resources.basics import get_checkpoint_path, qc_temp_prefix
from gnomad_qc.v3.resources.constants import CURRENT_RELEASE
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)

# Remove InbreedingCoeff from allele-specific fields (processed separately
# from other fields)
AS_FIELDS.remove("InbreedingCoeff")

VERSION = "4.0.0"  # passed arg
OUTPUT_TEMPLATE = (
    f"gs://gnomad-tmp/gnomad/v{version}/release/ht/release_sites_v{version}.ht"
)
TABLES_FOR_RELEASE = [
    "dbsnp",
    "final_filters",
    "freq",
    "info",
    "subsets",
    "region_flags",
]  # "in_silico", "vep",]


"""
Configurations of dataset to combine.
Format:
'<Name of dataset>': {
        'path': 'gs://path/to/hailtable.ht',
        'select': '<Optional list of fields to select or dict of new field name to location of old field
            in the reference dataset. If '#' is at the end, we know to select the appropriate biallelic
            using the a_index.>',
        'field_name': '<Optional name of root annotation in combined dataset, defaults to name of dataset.>',
        'custom_select': '<Optional function name of custom select function>',
    },
"""
CONFIG = {
    "dbsnp": {
        "ht": dbsnp.ht(),
        "path": dbsnp.path,
        "select": ["rsid"],
    },
    "filters": {
        "ht": final_filter().ht(),
        "path": final_filter().path,
        "custom_select": "custom_filters_select",
    },
    #   "in_silico": {
    #       "ht": analyst_annotations.ht(),
    #       "select": ["cadd", "revel", "splice_ai", "primate_ai"],
    #   },
    "info": {
        "ht": get_info().ht(),
        "path": get_info().path,
        "custom_select": "custom_info_select",
    },
    "freq": {
        "ht": get_freq(het_nonref_patch=True).ht(),
        "path": get_freq(het_nonref_patch=True).path,
        "select": ["freq", "faf", "popmax", "qual_hists", "raw_qual_hists"],
    },
    "subsets": {
        "ht": get_freq(het_nonref_patch=True).ht(),
        "path": get_freq(het_nonref_patch=True).path,
        "custom_select": "custom_subset_select",
        "field_name": "subsets",
    },
    #    "vep": {
    #        "ht": vep.ht(),
    # TODO: custom select -- drop 100% missing? Module to do this after all
    # annotations added?
    #        "select": ["vep"],
    #    },
    "region_flags": {
        "ht": get_freq(het_nonref_patch=True).ht(),
        "path": get_freq(het_nonref_patch=True).path,
        "custom_select": "custom_region_flags_select",
    },
}

VERSION = "testing"

# Remove InbreedingCoeff from allele-specific fields (processed separately
# from other fields)
AS_FIELDS.remove("InbreedingCoeff")


def custom_region_flags_select(ht):
    """
    Select region flags for release.

    :param ht: hail Table.
    :return: select expression dict
    """
    selects = {}
    selects = {
        "region_flags": region_flag_expr(
            ht,
            non_par=False,  # TODO: Konrad initially said to remove nonpar flag from HT because it's an easy hail operation, do we still want this?
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
        )
    }
    return selects


def custom_filters_select(ht):
    """
    Select gnomad filter HT fields for release dataset.

    Extract fields like 'filters', 'vqsr', and generates 'alelle_info' struct.
    :param ht: hail table.
    :return: select expression dict.
    """
    selects = {}
    selects["filters"] = ht.filters
    selects["vqsr"] = ht.vqsr
    selects["allele_info"] = hl.struct(
        variant_type=ht.variant_type,
        allele_type=ht.allele_type,
        n_alt_alleles=ht.n_alt_alleles,
        was_mixed=ht.was_mixed,
    )
    return selects


def custom_subset_select(ht):
    """
    Select release subset field using freq HT AN value.

    :param ht: hail table
    :return: select expression dict
    """
    selects = {
        subset: hl.if_else(
            ht.freq[ht.freq_index_dict[f"{subset}-adj"]].AN > 0, True, False
        )
        for subset in SUBSETS
    }
    return selects


def custom_info_select(ht):
    """
    Select fields for info hail Table annotation in release.

    :param ht: hail table
    :return: select expression dict
    """
    selects = {}

    filters_ht = final_filter().ht()
    filters = filters_ht[ht.key]
    filters_info_fields = [
        "singleton",
        "transmitted_singleton",
        "omni",
        "mills",
        "monoallelic",
        "SOR",
    ]
    filters_info_dict = {field: filters[field] for field in filters_info_fields}
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    filters_info_dict.update({**{f"{score_name}": filters[f"{score_name}"]}})

    info_dict = {field: ht.info[field] for field in SITE_FIELDS + AS_FIELDS}
    info_dict.update(filters_info_dict)

    selects["info"] = hl.struct(**info_dict)
    selects["was_split"] = ht.was_split
    selects["a_index"] = ht.a_index
    # TODO: InbreedingCoeff
    return selects


def get_select_fields(selects, base_ht):
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


def get_ht(dataset, _intervals, test) -> hl.Table:
    """
    Return the appropriate hail table with selects applied.

    :param dataset: Hail Table to join.
    :param _intervals: Intervals for reading in hail Table.
    :param test: Whether call is for a test run.
    :return: Hail Table with fields to select.
    """
    config = CONFIG[dataset]
    ht_path = config["path"]
    logger.info("Reading in %s", dataset)
    base_ht = hl.read_table(ht_path, _intervals=_intervals)

    if test:
        base_ht = hl.filter_intervals(
            base_ht,
            [hl.parse_locus_interval("chr1:1-1000000", reference_genome="GRCh38")],
        )

    if config.get("filter"):
        base_ht = base_ht.filter(config["filter"](base_ht))

    # 'select' and 'custom_select's to generate dict.
    select_fields = get_select_fields(config.get("select"), base_ht)
    if "custom_select" in config:
        custom_select_fn_str = config["custom_select"]
        select_fields = {**select_fields, **globals()[custom_select_fn_str](base_ht)}

    if config.get("field_name"):
        field_name = config.get("field_name")
        select_query = {field_name: hl.struct(**select_fields)}
    else:
        select_query = select_fields

    logger.info("%s", dataset)  # TODO: Remove after testing
    logger.info("%s", select_fields)  # TODO: Remove after testing
    return base_ht.select(**select_query)


def join_hts(base_table, tables, new_partition_percent, test, version=VERSION):
    """
    Outer join a list of hail tables.

    :param base_table: Dataset to use for interval partitioning.
    :param tables: List of tables to join.
    :param new_partition_percent: Percent of base_table partitions used for final release hail Table.
    :param test: Whether this is for a test run.
    :param version: Version of release.
    """
    logger.info(
        "Reading in %s to determine partition intervals for efficient join",
        base_table,
    )
    base_ht_path = CONFIG[base_table]["path"]
    base_ht = hl.read_table(base_ht_path)
    if test:
        base_ht = hl.filter_intervals(
            base_ht,
            [hl.parse_locus_interval("chr1:1-1000000", reference_genome="GRCh38")],
        )
    partition_intervals = base_ht._calculate_new_partitions(
        base_ht.n_partitions() * new_partition_percent
    )

    hts = [get_ht(table, _intervals=partition_intervals, test=test) for table in tables]
    joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "left")), hts)

    # Track the dataset we've added as well as the source path.
    included_dataset = {k: v["path"] for k, v in CONFIG.items() if k in tables}
    # Add metadata
    joined_ht = joined_ht.annotate_globals(
        date=datetime.now().isoformat(),
        datasets=hl.dict(included_dataset),
        version=version,
        # TODO: README=field definitions
        # TODO: tool versions= dict of tool versions (spliceAI,CADD, REVEL)
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
    ht = join_hts(
        args.base_table,
        args.tables_for_join,
        args.version,
        args.new_partition_percent,
        args.test,
    )

    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

    output_path = OUTPUT_TEMPLATE
    logger.info("Writing out release HT to %s", output_path)
    ht = ht.checkpoint(
        # qc_temp_prefix() + /release/gnomad.genomes.sites.test.ht"
        "gs://gnomad-tmp-4day/mwilson/release/gnomad.genomes.sites.test.ht"
        if args.test
        else release_sites().path,
        args.overwrite,
    )

    logger.info("Final variant count: %d", ht.count())
    ht.describe()
    ht.show()
    ht.write(output_path)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--new-partition-percent",
        help="Percent of start dataset partitions to use for release HT. Default is 1.1 (110%)",
        default=1.1,
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "-v", "--version", help="The version of gnomAD.", default=VERSION, required=True
    )
    parser.add_argument(
        "-t",
        "--test",
        help="Runs a test on the first two partitions of the HT.",
        action="store_true",
    )
    parser.add_argument(
        "-j",
        "--tables-for-join",
        help="Tables to join for release",
        default=TABLES_FOR_RELEASE,
    )
    parser.add_argument(
        "-b",
        "--base-table",
        help="Base table for interval partition calculation.",
        default="freq",
    )

    args = parser.parse_args()
    main(args)
