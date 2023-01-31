# noqa: D100

import argparse
import logging
from datetime import datetime
from functools import reduce

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: Update gnomAD_methods to dict: key version, value subsets
)
from gnomad.resources.grch38.gnomad import POPS, POPS_STORED_AS_SUBPOPS  # SUBPOPS,
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.resources.resource_utils import DataException
from gnomad.utils.annotations import missing_callstats_expr, region_flag_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.release import make_freq_index_dict  # TODO: Update to have subsets
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
AS_FIELDS.remove("AS_SOR")
SITE_FIELDS.remove("SOR")

# Add fine-resolution populations specific to 1KG/TGP and HGDP to gnomAD sub-populationss; used to create frequency index dictionary
# SUBPOPS.extend(POPS_STORED_AS_SUBPOPS)
# Add sub-populations to standard gnomAD pops; used to create frequency index dictionary
# POPS.extend(SUBPOPS)
# Add 'global' tag used to distinguish cohort-wide vs. subset annotations
# in frequency index dictionary
POPS.extend(["global"])

VERSION = "4.0.0"  # passed arg
OUTPUT_TEMPLATE = (
    "gs://gnomad-tmp/gnomad/v{version}/release/ht/release_sites_v{version}.ht"
)

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
    #    "in_silico": {
    #        "ht": analyst_annotations["3.1.2"].ht(),
    #        "select": ["cadd", "revel", "splice_ai", "primate_ai"],
    #    },
    "info": {
        "ht": get_info().ht(),
        "path": get_info().path,
        "select": {field: "info." + field for field in SITE_FIELDS + AS_FIELDS},
        "field_name": "info",
    },
    "freq": {
        "ht": get_freq(het_nonref_patch=True).ht(),
        "path": get_freq(het_nonref_patch=True).path,
        "select": ["freq", "faf", "popmax", "qual_hists", "raw_qual_hists"],
    },
    #    "vep": {
    #        "ht": vep.ht(),
    # TODO: custom select -- drop 100% missing? Module to do this after all
    # annotations added
    #        "select": ["vep"],
    #    },
}

VERSION = "testing"

# Remove InbreedingCoeff from allele-specific fields (processed separately
# from other fields)
AS_FIELDS.remove("InbreedingCoeff")

# Add fine-resolution populations specific to 1KG/TGP and HGDP to standard
# gnomAD pops; used to create frequency index dictionary
POPS.extend(POPS_STORED_AS_SUBPOPS)
# Add 'global' tag used to distinguish cohort-wide vs. subset annotations
# in frequency index dictionary
POPS.extend(["global"])


def custom_filters_select(ht):
    """
    Select gnomad filter HT fields for release dataset.

    Extract fields like 'filters', 'vqsr', and generates 'alelle_info' struct.
    :param ht: hail table
    :return: select expression dict
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


def get_select_fields(selects, base_ht):
    """
    Take in a select config and base_ht and generate a select dict from traversing the base_ht and extracting annotations.

    If '#' is included at the end of a select field, the appropriate biallelic position will be selected (e.g. 'x#' -> x[base_ht.a_index-1].
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
                    # Select from multi-allelic list.
                    if attr.endswith("#"):
                        attr = attr[:-1]
                        ht = ht[attr][base_ht.a_index - 1]
                    else:
                        ht = ht[attr]
                select_fields[key] = ht
    return select_fields


def get_ht(dataset, _intervals, test):
    """Return the appropriate deduped hail table with selects applied."""
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

    logger.info("%s", select_fields)
    logger.info("%s", dataset)
    return base_ht.select(**select_query).distinct()


def join_hts(
    base_dataset,
    datasets,
    new_partition_percent,
    test,
    version=VERSION,
):
    """Get a list of hail tables and combine into an outer join."""
    logger.info(
        "Reading in %s to determine partition intervals for efficient join",
        base_dataset,
    )
    base_ht_path = CONFIG[base_dataset]["path"]
    base_ht = hl.read_table(base_ht_path)
    base_ht = hl.filter_intervals(
        base_ht, [hl.parse_locus_interval("chr1:1-1000000", reference_genome="GRCh38")]
    )
    partition_intervals = base_ht._calculate_new_partitions(
        base_ht.n_partitions() * new_partition_percent
    )

    hts = [
        get_ht(dataset, _intervals=partition_intervals, test=test)
        for dataset in datasets
    ]
    joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "left")), hts)

    # Track the dataset we've added as well as the source path.
    included_dataset = {k: v["path"] for k, v in CONFIG.items() if k in datasets}
    # Add metadata, but also removes previous globals.
    joined_ht = joined_ht.annotate_globals(
        date=datetime.now().isoformat(),
        datasets=hl.dict(included_dataset),
        version=version,
    )
    joined_ht.describe()
    return joined_ht


def main(args):
    """Create release ht."""
    hl.init(
        log="/create_release_ht.log",
        tmp_dir="gs://gnomad-tmp",
        default_reference="GRCh38",
    )
    ht = join_hts(
        [
            "dbsnp",
            "final_filters",
            "freq",
            "info",
            "in_silico",
            "vep",
        ],
        args.version,
        args.new_partition_percent,
        args.test,
    )

    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters))

    output_path = os.path.join(OUTPUT_TEMPLATE.format(version=VERSION))
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

    args = parser.parse_args()
    main(args)
