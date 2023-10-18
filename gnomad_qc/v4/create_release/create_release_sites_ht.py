"""Script to create release sites HT for exomes."""
import argparse
import json
import logging
from copy import deepcopy
from datetime import datetime
from functools import reduce
from typing import Dict, List, Union

import hail as hl
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.utils.annotations import region_flag_expr
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import AS_FIELDS, SITE_FIELDS
from hail.typecheck import anytype, nullable, sequenceof

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.annotations.insilico_predictors import get_sift_polyphen_from_vep
from gnomad_qc.v4.create_release.create_release_utils import (
    DBSNP_VERSION,
    GENCODE_VERSION,
    MANE_SELECT_VERSION,
    POLYPHEN_VERSION,
    SEQREPO_VERSION,
    SIFT_VERSION,
    VRS_PYTHON_VERSION,
    VRS_SCHEMA_VERSION,
    remove_missing_vep_fields,
)
from gnomad_qc.v4.resources.annotations import (
    get_freq,
    get_info,
    get_insilico_predictors,
    get_vep,
    get_vrs,
)
from gnomad_qc.v4.resources.basics import calling_intervals, qc_temp_prefix
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.release import (
    FREQUENCY_README,
    get_combined_faf_release,
    included_datasets_json_path,
    release_sites,
)
from gnomad_qc.v4.resources.sample_qc import interval_qc_pass
from gnomad_qc.v4.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)

# Remove InbreedingCoeff from allele-specific fields because it is processed separately
# from the other fields.
AS_FIELDS = deepcopy(AS_FIELDS)
AS_FIELDS.remove("InbreedingCoeff")
SITE_FIELDS = deepcopy(SITE_FIELDS)
TABLES_FOR_RELEASE = [
    "dbsnp",
    "filters",
    "freq",
    "info",
    "region_flags",
    "in_silico",
    "vep",
    "joint_faf",
]

INSILICO_PREDICTORS = ["cadd", "revel", "spliceai", "pangolin", "phylop"]

FINALIZED_SCHEMA = {
    "globals": [
        "freq_meta",
        "freq_index_dict",
        "freq_meta_sample_count",
        "faf_meta",
        "faf_index_dict",
        "joint_freq_meta",
        "joint_freq_index_dict",
        "joint_freq_meta_sample_count",
        "joint_faf_meta",
        "joint_faf_index_dict",
        "age_distribution",
        "downsamplings",
        "filtering_model",
        "inbreeding_coeff_cutoff",
        "interval_qc_parameters",
        "tool_versions",
        "vrs_versions",
        "vep_globals",
        "frequency_README",
        "date",
        "version",
    ],
    "rows": [
        "freq",
        "grpmax",
        "faf",
        "fafmax",
        "joint_freq",
        "joint_grpmax",
        "joint_faf",
        "joint_fafmax",
        "a_index",
        "was_split",
        "rsid",
        "filters",
        "info",
        "vep",
        "vqsr_results",
        "region_flags",
        "allele_info",
        "histograms",
        "in_silico_predictors",
    ],
}


# Config is added as a function, so it is not evaluated until the function is called.
def get_config(
    release_exists: bool = False,
) -> Dict[str, Dict[str, hl.expr.Expression]]:
    """
    Get configuration dictionary.

    Format:
        '<Name of dataset>': {
                'ht': '<Optional Hail Table for direct annotation extraction. This is not used for the join.>',
                'path': 'gs://path/to/hail_table.ht',
                'select': '<Optional list of fields to select or dict of new field name to location of old field in the dataset.>',
                'field_name': '<Optional name of root annotation in combined dataset, defaults to name of dataset.>',
                'custom_select': '<Optional function name of custom select function that is needed for more advanced logic>',
                'select_globals': '<Optional list of globals to select or dict of new global field name to old global field name. If not specified, all globals are selected.>'
            },

    .. warning::

        The 'in_silico' key's 'ht' logic is handled separately because it is a list of
        HTs. In this list, the phyloP HT is keyed by locus only and thus the 'ht' code
        below sets the join key to 1, which will grab the first key of
        ht.key.dtype.values() e.g. 'locus', when a HT's keys are not {'locus, 'alleles'}.
        All future in_silico predictors should have the keys confirmed to be 'locus'
        with or without 'alleles' before using this logic.

    :param release_exists: Whether the release HT already exists.
    :return: Dict of dataset's configs.
    """
    config = {
        "dbsnp": {
            "ht": dbsnp.ht(),
            "path": dbsnp.path,
            "select": ["rsid"],
        },
        "filters": {
            "ht": final_filter().ht(),
            "path": final_filter().path,
            "select": ["filters"],
            "custom_select": custom_filters_select,
            "select_globals": ["filtering_model", "inbreeding_coeff_cutoff"],
        },
        "in_silico": {
            "ht": reduce(
                (
                    lambda joined_ht, ht: (
                        joined_ht.join(ht, "outer")
                        if set(ht.key) == {"locus", "alleles"}
                        else joined_ht.join(ht, "outer", _join_key=1)
                    )
                ),
                [
                    get_insilico_predictors(predictor=predictor).ht()
                    for predictor in INSILICO_PREDICTORS
                ],
            ),
            "path": [
                get_insilico_predictors(predictor=predictor).path
                for predictor in INSILICO_PREDICTORS
            ],
            "field_name": "in_silico_predictors",
            "select": [
                "cadd",
                "revel_max",
                "spliceai_ds_max",
                "pangolin_largest_ds",
                "phylop",
            ],
            "custom_select": custom_in_silico_select,
            "select_globals": [
                "cadd_version",
                "revel_version",
                "spliceai_version",
                "pangolin_version",
                "phylop_version",
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
                "histograms",
            ],
            "custom_select": custom_freq_select,
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
            "path": get_vep().path,
            "select": ["vep"],
            "custom_select": custom_vep_select,
            "select_globals": [
                "vep_version",
                "vep_help",
                "vep_config",
            ],
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
        "joint_faf": {
            "ht": get_combined_faf_release().ht(),
            "path": get_combined_faf_release().path,
            "select": ["joint_freq", "joint_faf", "joint_fafmax"],
            "custom_select": custom_joint_faf_select,
            "select_globals": [
                "joint_freq_meta",
                "joint_freq_index_dict",
                "joint_freq_meta_sample_count",
                "joint_faf_meta",
                "joint_faf_index_dict",
            ],
        },
    }

    if release_exists:
        config["release"].update(
            {
                "ht": release_sites().ht(),
                "select": [r for r in release_sites().ht().row],
                "select_globals": [g for g in release_sites().ht().globals],
            }
        )
    return config


def custom_joint_faf_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Drop faf95 from 'grpmax'.

    This annotation will be combined with the others from joint_faf's select in the config.
    See note in `custom_freq_select` explaining why this field is removed.

    :param ht: Joint FAF Hail Table.
    :return: Select expression dict.
    """
    selects = {"joint_grpmax": ht.joint_grpmax.drop("faf95")}

    return selects


def custom_freq_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Drop faf95 from both 'gnomad' and 'non_ukb' in 'grpmax' and rename `gen_anc_faf_max` to `fafmax`.

    These annotations will be combined with the others from freq's select in the config.

    .. note::
        - The faf95 field in the grpmax struct is the FAF of the genetic ancestry group with the largest AF (grpmax AF).
        - The FAF fields within the gen_anc_faf_max struct contains the FAFs from the genetic ancestry group(s) with the largest FAFs
        - These values aren't necessarily the same; the group with the highest AF for a variant isn't necessarily the group with the highest FAF for a variant
        - The filtering allele frequencies that are used by the community are the values within the gen_anc_faf_max struct, NOT grpmax FAF, which is why we are dropping grpmax.faf95 and renaming gen_anc_faf_max

    :param ht: Freq Hail Table
    :return: Select expression dict.
    """
    selects = {
        "grpmax": hl.struct(**{k: ht.grpmax[k].drop("faf95") for k in ht.grpmax}),
        "fafmax": ht.gen_anc_faf_max,
    }

    return selects


def custom_in_silico_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Get in silico predictors from VEP for release.

    This function currently selects only SIFT and Polyphen from VEP.

    :param ht: VEP Hail Table.
    :return: Select expression dict.
    """
    vep_in_silico = get_sift_polyphen_from_vep(get_vep().ht())
    selects = {
        "sift_max": vep_in_silico[ht.key].sift_max,
        "polyphen_max": vep_in_silico[ht.key].polyphen_max,
    }
    return selects


def custom_region_flags_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select region flags for release.

    :param ht: Hail Table.
    :return: Select expression dict.
    """
    selects = {
        "region_flags": region_flag_expr(
            ht,
            non_par=True,
            prob_regions={"lcr": lcr_intervals.ht(), "segdup": seg_dup_intervals.ht()},
        )
    }
    selects["region_flags"] = selects["region_flags"].annotate(
        fail_interval_qc=~interval_qc_pass(all_platforms=True)
        .ht()[ht.locus]
        .pass_interval_qc,
        outside_ukb_capture_region=~hl.is_defined(
            calling_intervals(interval_name="ukb", calling_interval_padding=50).ht()[
                ht.locus
            ]
        ),
        outside_broad_capture_region=~hl.is_defined(
            calling_intervals(interval_name="broad", calling_interval_padding=50).ht()[
                ht.locus
            ]
        ),
    )

    return selects


def custom_filters_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select gnomAD filter HT fields for release dataset.

    Extract "results" field and rename based on filtering method.

    :param ht: Filters Hail Table.
    :return: Select expression dict.
    """
    filter_name = hl.eval(ht.filtering_model.filter_name)
    if filter_name == "RF":
        name = "random_forest_results"
    elif filter_name == "AS_VQSR":
        name = "vqsr_results"
    elif filter_name == "IF":
        name = "isolation_forest_results"
    else:
        raise ValueError(f"Filtering method {filter_name} not recognized.")
    selects = {name: ht.results.annotate(**ht.training_info)}

    return selects


def custom_info_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select fields for info Hail Table annotation in release.

    The info field requires fields from the freq HT and the filters HT so those are
    pulled in here along with all info HT fields. It also adds the `allele_info` struct
    to release HT.

    :param ht: Info Hail Table.
    :return: Select expression dict.
    """
    # Create a dict of the fields from the filters HT that we want to add to the info.
    filters_ht = get_config().get("filters")["ht"]

    # For more information on the selected compute_info_method, please see the
    # run_compute_info in gnomad_qc.v4.annotations.generate_variant_qc_annotations.py
    compute_info_method = hl.eval(
        filters_ht.filtering_model_specific_info.compute_info_method
    )
    score_name = hl.eval(filters_ht.filtering_model.score_name)
    filters_ht = filters_ht.transmute(**filters_ht.truth_sets)
    filters = filters_ht[ht.key]
    filters_info_fields = [
        "singleton",
        "transmitted_singleton",
        "sibling_singleton",
        "omni",
        "mills",
        "monoallelic",
        "only_het",
    ]
    filters_info_dict = {field: filters[field] for field in filters_info_fields}
    filters_info_dict.update({**{f"{score_name}": filters[f"{score_name}"]}})

    # Create a dict of the fields from the freq HT that we want to add to the info.
    freq_ht = get_config().get("freq")["ht"]
    freq_info_dict = {"inbreeding_coeff": freq_ht[ht.key]["inbreeding_coeff"]}

    # Create a dict of the fields from the VRS HT that we want to add to the info.
    vrs_ht = get_vrs().ht()
    vrs_info_fields = {"vrs": vrs_ht[ht.key].vrs}

    # Create a dict of the fields from the info HT that we want keep in the info.
    info_struct = hl.struct(**ht.site_info, **ht[f"{compute_info_method}_info"])
    info_dict = {field: info_struct[field] for field in SITE_FIELDS + AS_FIELDS}
    info_dict.update(filters_info_dict)
    info_dict.update(freq_info_dict)
    info_dict.update(vrs_info_fields)

    # Select the info and allele info annotations. We drop nonsplit_alleles from
    # allele_info so that we don't release alleles that are found in non-releasable
    # samples.
    selects = {
        "info": hl.struct(**info_dict),
        "allele_info": ht.allele_info.drop("nonsplit_alleles"),
    }

    return selects


def custom_vep_select(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select fields for VEP hail Table annotation in release.

    :param ht: VEP Hail table
    :return: Select expression dict.
    """
    vep_expr = remove_missing_vep_fields(ht.vep)
    selects = {
        "vep": vep_expr.annotate(
            transcript_consequences=vep_expr.transcript_consequences.map(
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


def get_select_global_fields(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Generate a dictionary of globals to select by checking the config of all tables joined.

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


def get_select_fields(
    selects: Union[List, Dict], base_ht: hl.Table
) -> Dict[str, hl.expr.Expression]:
    """
    Generate a select dict from traversing the `base_ht` and extracting annotations.

    :param selects: Mapping or list of selections.
    :param base_ht: Base Hail Table to traverse.
    :return: select Mapping from annotation name to `base_ht` annotation.
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


def get_ht(
    dataset: str,
    _intervals: nullable(sequenceof(anytype)),
    test: bool,
    release_exists: bool,
) -> hl.Table:
    """
    Return the appropriate Hail table with selects applied.

    :param dataset: Hail Table to join.
    :param _intervals: Intervals for reading in hail Table. Used to optimize join.
    :param test: Whether call is for a test run.
    :param release_exists: Whether the release HT already exists.
    :return: Hail Table with fields to select.
    """
    logger.info("Getting the %s dataset and its selected annotations...", dataset)
    config = get_config(release_exists=release_exists)[dataset]

    # There is no single path for insilico predictors, so we need to handle this case
    # separately.
    if dataset == "in_silico":
        base_ht = config["ht"]
    else:
        ht_path = config["path"]
        logger.info("Reading in %s", dataset)
        base_ht = hl.read_table(ht_path, _intervals=_intervals)

    if test:
        # Keep only PCSK9.
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
    Outer join a list of Hail Tables.

    :param base_table: Dataset to use for interval partitioning.
    :param tables: List of tables to join.
    :param new_partition_percent: Percent of base_table partitions used for final
        release Hail Table.
    :param test: Whether this is for a test run.
    :param release_exists: Whether the release HT already exists.
    :return: Hail Table with datasets joined.
    """
    if base_table not in tables:
        raise ValueError(f"Base table {base_table} must be in tables to join: {tables}")

    logger.info(
        "Reading in %s to determine partition intervals for efficient join",
        base_table,
    )
    base_ht_path = get_config()[base_table]["path"]
    base_ht = hl.read_table(base_ht_path)
    if test:
        # Filter to PCSK9 for testing.
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

    # Reorg list so base table is first.
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
    joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "left")), hts)

    # Track the datasets we've added as well as the source paths.
    # If release HT is included in tables, read in the included datasets json
    # and update the keys to the path for any new tables
    included_datasets = {}
    if "release" in tables:
        with hl.utils.hadoop_open(
            included_datasets_json_path(
                test=test,
                release_version=hl.eval(get_config()["release"]["ht"].version),
            )
        ) as f:
            included_datasets = json.loads(f.read())

    included_datasets.update(
        {t: get_config(release_exists=release_exists)[t]["path"] for t in tables}
    )
    with hl.utils.hadoop_open(
        included_datasets_json_path(test=test, release_version=args.version), "w"
    ) as f:
        f.write(hl.eval(hl.json(included_datasets)))

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

    # Filter out chrM, AS_lowqual sites (these sites are dropped in the final_filters HT
    # so will not have information in `filters`) and AC_raw == 0.
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters) & (ht.freq[1].AC > 0))

    ht = ht.select_globals(**get_select_global_fields(ht))

    # Add additional globals that were not present on the joined HTs.
    ht = ht.annotate_globals(
        filtering_model=ht.filtering_model.drop("model_id"),
        vep_globals=ht.vep_globals.annotate(
            gencode_version=GENCODE_VERSION,
            mane_select_version=MANE_SELECT_VERSION,
        ),
        tool_versions=ht.tool_versions.annotate(
            dbsnp_version=DBSNP_VERSION,
            sift_version=SIFT_VERSION,
            polyphen_version=POLYPHEN_VERSION,
        ),
        vrs_versions=hl.struct(
            **{
                "vrs_schema_version": VRS_SCHEMA_VERSION,
                "vrs_python_version": VRS_PYTHON_VERSION,
                "seqrepo_version": SEQREPO_VERSION,
            },
        ),
        date=datetime.now().isoformat(),
        version=args.version,
        interval_qc_parameters=interval_qc_pass(all_platforms=True)
        .ht()
        .index_globals()
        .high_qual_interval_parameters,
        frequency_README=FREQUENCY_README,
    )

    # Reorder fields to match final schema.
    ht = ht.select(*FINALIZED_SCHEMA["rows"]).select_globals(
        *FINALIZED_SCHEMA["globals"]
    )

    output_path = (
        f"{qc_temp_prefix()}release/gnomad.exomes.sites.test.updated_101723.ht"
        if args.test
        else release_sites().path
    )
    logger.info("Writing out release HT to %s", output_path)
    ht = ht.checkpoint(
        output_path,
        args.overwrite,
    )

    logger.info("Final release HT schema:")
    ht.describe()

    logger.info("Final variant count: %d", ht.count())


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
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    args = parser.parse_args()
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
