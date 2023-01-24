# noqa: D100
import argparse
import logging
import os
from datetime import datetime
from functools import reduce

import hail as hl
from gnomad.resources.grch38.gnomad import POPS, POPS_STORED_AS_SUBPOPS, SUBSETS
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
    "1kg": {
        "path": "gs://seqr-reference-data/GRCh37/1kg/1kg.wgs.phase3.20130502.GRCh37_sites.ht",
        "select": {
            "AC": "info.AC#",
            "AF": "info.AF#",
            "AN": "info.AN",
            "POPMAX_AF": "POPMAX_AF",
        },
        "field_name": "g1k",
    },
    "dbnsfp": {
        "path": "gs://seqr-reference-data/GRCh37/dbNSFP/v2.9.3/dbNSFP2.9.3_variant.ht",
        "select": [
            "SIFT_pred",
            "Polyphen2_HVAR_pred",
            "MutationTaster_pred",
            "FATHMM_pred",
            "MetaSVM_pred",
            "REVEL_score",
            "GERP_RS",
            "phastCons100way_vertebrate",
        ],
    },
    "primate_ai": {
        "path": "gs://seqr-reference-data/GRCh38/primate_ai/PrimateAI_scores_v0.2.liftover_grch38.ht",
        "select": {"score": "info.score"},
    },
    "gnomad_exomes": {
        "path": "gs://gcp-public-data--gnomad/release/2.1.1/liftover_grch38/ht/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht",
        "custom_select": "custom_gnomad_select_v2",
    },
}

# Remove InbreedingCoeff from allele-specific fields (processed separately
# from other fields)
AS_FIELDS.remove("InbreedingCoeff")

# Add fine-resolution populations specific to 1KG/TGP and HGDP to standard
# gnomAD pops; used to create frequency index dictionary
POPS.extend(POPS_STORED_AS_SUBPOPS)
# Add 'global' tag used to distinguish cohort-wide vs. subset annotations
# in frequency index dictionary
POPS.extend(["global"])


def custom_gnomad_select_v2(ht):
    """
    Select for public gnomad v2 dataset (which we did not generate).

    Extract fields like 'AF', 'AN', and generates 'hemi'.
    :param ht: hail table
    :return: select expression dict
    """
    selects = {}
    global_idx = hl.eval(ht.globals.freq_index_dict["gnomad"])
    selects["AF"] = ht.freq[global_idx].AF
    selects["AN"] = ht.freq[global_idx].AN
    selects["AC"] = ht.freq[global_idx].AC
    selects["Hom"] = ht.freq[global_idx].homozygote_count

    selects["AF_POPMAX_OR_GLOBAL"] = hl.or_else(
        ht.popmax[ht.globals.popmax_index_dict["gnomad"]].AF, ht.freq[global_idx].AF
    )
    selects["FAF_AF"] = ht.faf[ht.globals.popmax_index_dict["gnomad"]].faf95
    selects["Hemi"] = hl.if_else(
        ht.locus.in_autosome_or_par(),
        0,
        ht.freq[ht.globals.freq_index_dict["gnomad_male"]].AC,
    )
    return selects


def get_select_fields(selects, base_ht):
    """
    Take in a select config and base_ht and generatee a select dict from traversing the base_ht and extracting annotations.

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


def get_ht(dataset, _intervals):
    """Return the appropriate deduped hail table with selects applied."""
    config = CONFIG[dataset]
    ht_path = config["path"]
    logger.info("Reading in %s", dataset)
    base_ht = hl.read_table(ht_path, _intervals=_intervals)

    if config.get("filter"):
        base_ht = base_ht.filter(config["filter"](base_ht))

    # 'select' and 'custom_select's to generate dict.
    select_fields = get_select_fields(config.get("select"), base_ht)
    if "custom_select" in config:
        custom_select_fn_str = config["custom_select"]
        select_fields = {**select_fields, **globals()[custom_select_fn_str](base_ht)}

    field_name = config.get("field_name") or dataset
    select_query = {field_name: hl.struct(**select_fields)}

    logger.info("%s", select_fields)
    return base_ht.select(**select_query).distinct()


def join_hts(base_dataset, datasets, new_partition_percent, coverage_datasets=[]):
    """Get a list of hail tables and combine into an outer join."""
    logger.info("Reading in %s to determine partition intervals", base_dataset)
    base_ht_path = CONFIG[base_dataset]["path"]
    base_ht = hl.read_table(base_ht_path)
    partition_intervals = base_ht.calculate_new_partitions(
        base_ht.n_partitions() * new_partition_percent
    )

    hts = [get_ht(dataset, _intervals=partition_intervals) for dataset in datasets]
    joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "outer")), hts)

    # Annotate coverages.
    for coverage_dataset in coverage_datasets:
        joined_ht = annotate_coverages(joined_ht, coverage_dataset)

    # Track the dataset we've added as well as the source path.
    included_dataset = {
        k: v["path"] for k, v in CONFIG.items() if k in datasets + coverage_datasets
    }
    # Add metadata, but also removes previous globals.
    joined_ht = joined_ht.select_globals(
        date=datetime.now().isoformat(),
        datasets=hl.dict(included_dataset),
        version=VERSION,
    )
    joined_ht.describe()
    return joined_ht


def run(args):
    """Create release ht."""
    hl.init(
        log="/create_release_ht.log",
        tmp_dir="gs://gnomad-tmp",
        default_reference="GRCh38",
    )
    joined_ht = join_hts(
        [
            "cadd",
            "1kg",
            "mpc",
            "eigen",
            "dbnsfp",
            "topmed",
            "primate_ai",
            "splice_ai",
            "exac",
            "gnomad_genomes",
            "gnomad_exomes",
            "geno2mp",
        ],
        ["gnomad_genome_coverage", "gnomad_exome_coverage"],
        args.version,
        args.new_partition_percent,
    )
    output_path = os.path.join(OUTPUT_TEMPLATE.format(version=VERSION))
    logger.info("Writing to %s", output_path)
    joined_ht.write(os.path.join(output_path))


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-v", "--version", help="The version of gnomAD.", default=VERSION, required=True
    )
    parser.add_argument(
        "--new-partition-percent",
        help="Percent of start dataset partitions to use for release HT. Default is 1.1 (110%)",
        default=1.1,
    )
    args = parser.parse_args()

    run(args)
