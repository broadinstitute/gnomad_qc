"""Script to create release sites HT for v5 genomes."""

import argparse
import json
import logging
from datetime import datetime
from functools import reduce
from typing import Any, Dict, List, Optional, Union

import hail as hl
from gnomad.resources.grch38.reference_data import (
    dbsnp,
    lcr_intervals,
    seg_dup_intervals,
)
from gnomad.utils.annotations import region_flag_expr
from gnomad.utils.vcf import ALLELE_TYPE_FIELDS, AS_FIELDS, SITE_FIELDS
from hail.utils import new_temp_file

from gnomad_qc.v5.annotations.insilico_predictors import get_sift_polyphen_from_vep
from gnomad_qc.v5.create_release.create_release_utils import (
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
from gnomad_qc.v5.resources.annotations import (
    get_freq,
    get_info_ht,
    get_insilico_predictors,
    get_vep,
    get_vrs,
)
from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import (
    CURRENT_RELEASE,
    GNOMAD_TMP_BUCKET,
    WORKSPACE_BUCKET,
)
from gnomad_qc.v5.resources.release import (
    FREQUENCY_README,
    included_datasets_json_path,
    release_sites,
)
from gnomad_qc.v5.resources.variant_qc import final_filter

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_release_ht")
logger.setLevel(logging.INFO)
# TODO: check AS FIELDS AND SITE FIELDS
TABLES_FOR_RELEASE = [
    "dbsnp",
    "filters",
    "freq",
    "info",
    "region_flags",
    "in_silico",
    "vep",
]

INSILICO_PREDICTORS = ["cadd", "revel", "spliceai", "pangolin", "phylop"]

# check release config
FINALIZED_SCHEMA = {
    "globals": [
        "freq_meta",
        "freq_index_dict",
        "freq_meta_sample_count",
        "faf_meta",
        "faf_index_dict",
        "age_distribution",
        "aou_downsamplings",
        "filtering_model",
        "inbreeding_coeff_cutoff",
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
    release_exists: bool = False, test: bool = False, environment: str = "rwb"
) -> Dict[str, Dict[str, hl.expr.Expression]]:
    """
    Get configuration dictionary for specified data type.

    Format:

    .. code-block::

        '<Name of dataset>': {
            'ht': '<Optional Hail Table for direct annotation extraction. This is not used for the join.>',
            'path': 'gs://path/to/hail_table.ht',
            'select': '<Optional list of fields to select or dict of new field name to location of old field in the dataset.>',
            'field_name': '<Optional name of root annotation in combined dataset, defaults to name of dataset.>',
            'custom_select': '<Optional function name of custom select function that is needed for more advanced logic>',
            'select_globals': '<Optional list of globals to select or dict of new global field name to old global field name. If not specified, all globals are selected.>'
            'custom_globals_select': '<Optional function name of custom globalsselect function that is needed for more advanced logic>'

        },

    .. warning::

        The 'in_silico' key's 'ht' logic is handled separately because it is a list of
        HTs. In this list, the phyloP HT is keyed by locus only and thus the 'ht' code
        below sets the join key to 1, which will grab the first key of
        ht.key.dtype.values() e.g. 'locus', when an HT's keys are not {'locus',
        'alleles'}.
        All future in_silico predictors should have the keys confirmed to be 'locus'
        with or without 'alleles' before using this logic.

    :param release_exists: Whether the release HT already exists.
    :param test: Whether or not to use the test versions of the datasets when available.
    :param environment: Environment to use. Default is "rwb".
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
            "custom_globals_select": custom_filters_select_globals,
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
            "ht": get_info_ht(test=test, environment=environment).ht(),
            "path": get_info_ht(test=test, environment=environment).path,
            "select": ["was_split", "a_index"],
            "custom_select": custom_info_select,
        },
        "freq": {
            "ht": get_freq(
                test=test,
                data_type="genomes",
                data_set="merged",
                environment=environment,
            ).ht(),
            "path": get_freq(
                test=test,
                data_type="genomes",
                data_set="merged",
                environment=environment,
            ).path,
            "select": [
                "freq",
                "faf",
                "histograms",
                "grpmax",
                "fafmax",
            ],
            "select_globals": [
                "freq_meta",
                "freq_index_dict",
                "freq_meta_sample_count",
                "faf_meta",
                "faf_index_dict",
                "age_distribution",
                "aou_downsamplings",
            ],
        },
        "vep": {
            "ht": get_vep(data_type="genomes").ht(),
            "path": get_vep(data_type="genomes").path,
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
            "ht": get_freq(
                test=test,
                data_type="genomes",
                data_set="merged",
                environment=environment,
            ).ht(),
            "path": get_freq(
                test=test,
                data_type="genomes",
                data_set="merged",
                environment=environment,
            ).path,
            "custom_select": custom_region_flags_select,
        },
        "release": {
            "path": release_sites(environment=environment).path,
        },
    }

    if release_exists:
        config["release"].update(
            {
                "ht": release_sites(environment=environment).ht(),
                "select": [r for r in release_sites(environment=environment).ht().row],
                "select_globals": [
                    g for g in release_sites(environment=environment).ht().globals
                ],
            }
        )
    return config


def custom_in_silico_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
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


def custom_region_flags_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
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

    return selects


def custom_filters_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
    """
    Select gnomAD filter HT fields for release dataset.

    Extract "results" field and rename based on filtering method.

    :param ht: Filters Hail Table.
    :return: Select expression dict.
    """
    filter_to_name = {
        "RF": "random_forest_results",
        "AS_VQSR": "vqsr_results",
        "IF": "isolation_forest_results",
    }

    filter_name = hl.eval(ht.filtering_model.filter_name)

    if filter_name not in filter_to_name:
        raise ValueError(f"Filtering method '{filter_name}' not recognized.")

    name = filter_to_name[filter_name]
    selects = {name: ht.results.annotate(**ht.training_info)}

    return selects


def custom_filters_select_globals(ht: hl.Table) -> Dict[str, hl.expr.Expression]:
    """
    Select filter HT globals for release dataset.

    :param ht: Filters Hail Table.
    :return: Select expression dict.
    """
    selects = {
        "filtering_model": hl.struct(
            **{
                "filter_name": ht.filtering_model.filter_name,
                "score_name": ht.filtering_model.score_name,
                "snv_cutoff": ht.filtering_model.snv_cutoff.drop("bin_id"),
                "indel_cutoff": ht.filtering_model.indel_cutoff.drop("bin_id"),
                "snv_training_variables": ht.filtering_model.snv_training_variables,
                "indel_training_variables": ht.filtering_model.indel_training_variables,
            }
        ),
        "inbreeding_coeff_cutoff": ht.inbreeding_coeff_cutoff,
    }

    return selects


def custom_info_select(
    ht: hl.Table, data_type: str, config
) -> Dict[str, hl.expr.Expression]:
    """
    Select fields for info Hail Table annotation in release.

    The info field requires fields from the freq HT and the filters HT so those are
    pulled in here along with all info HT fields. It also adds the `allele_info` struct
    to release HT.

    :param ht: Info Hail Table.
    :return: Select expression dict.
    """
    # Create a dict of the fields from the filters HT that we want to add to the info.
    filters_ht = config.get("filters")["ht"]

    score_name = hl.eval(filters_ht.filtering_model.score_name)
    filters_ht = filters_ht.transmute(**filters_ht.truth_sets)
    filters = filters_ht[ht.key]
    filters_info_fields = [
        "singleton",
        "transmitted_singleton",
        # "sibling_singleton",
        "omni",
        "mills",
        "monoallelic",
        "only_het",
    ]
    # if data_type == "genomes":
    #    filters_info_fields.remove("sibling_singleton")
    filters_info_dict = {field: filters[field] for field in filters_info_fields}
    filters_info_dict[score_name] = filters[score_name]

    # Create a dict of the fields from the freq HT that we want to add to the info.
    freq_ht = config.get("freq")["ht"]
    # TODO: will change back to 'inbreeding_coeff' once we have the new freq_ht?
    freq_info_dict = {"InbreedingCoeff": freq_ht[ht.key]["InbreedingCoeff"]}

    # Create a dict of the fields from the VRS HT that we want to add to the info.
    vrs_ht = get_vrs(data_type=data_type).ht()
    vrs_info_fields = {"vrs": vrs_ht[ht.key].vrs}

    # Create a dict of the fields from the info HT that we want keep in the info.
    # v3 info HT has no SOR or AS_SOR fields. They are computed by VQSR, so we can
    # grab them from the filters HT.
    # TODO: Check back on v4 logic, AS_SOR is already on info ht, is SOR needed?
    info_struct = hl.struct(**ht.info, SOR=filters.SOR)

    info_dict = {field: info_struct[field] for field in SITE_FIELDS + AS_FIELDS}
    info_dict.update(filters_info_dict)
    info_dict.update(freq_info_dict)
    info_dict.update(vrs_info_fields)

    # Select the info and allele info annotations. We drop nonsplit_alleles from
    # allele_info so that we don't release alleles that are found in non-releasable
    # samples.
    selects = {"info": hl.struct(**info_dict)}
    # TODO: Check back on v4 logic.
    selects["allele_info"] = ht.allele_info.drop("nonsplit_alleles")

    return selects


def custom_vep_select(ht: hl.Table, **_) -> Dict[str, hl.expr.Expression]:
    """
    Select fields for VEP hail Table annotation in release.

    :param ht: VEP Hail table
    :return: Select expression dict.
    """
    # TODO: Update based on any updates to the VEP fields.
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


def get_select_global_fields(
    ht: hl.Table,
    config: Dict[str, Dict[str, Any]],
    tables_for_join: List[str] = TABLES_FOR_RELEASE,
) -> Dict[str, hl.expr.Expression]:
    """
    Generate a dictionary of globals to select by checking the config of all tables joined.

    .. note::

        This function will place the globals within the select_globals value above
        any globals returned from custom_select_globals. If ordering is important, use
        only custom_select_globals.

    :param ht: Final joined HT with globals.
    :param config: Dictionary with configuration options for each dataset. Expects the
        dataset name (matching values in `tables_for_join`) as the key and a dictionary
        of options as the value, where the options include any of the following keys:
        'select_globals', 'custom_globals_select', 'global_name'.
    :param tables_for_join: List of tables to join into final release HT.
    :return: Select mapping from global annotation name to `ht` annotation.
    """
    t_globals = []
    select_globals = {}
    for t in tables_for_join:
        t_config = config.get(t)
        if "select_globals" in t_config:
            select_globals = get_select_fields(t_config["select_globals"], ht)
        if "custom_globals_select" in t_config:
            custom_globals_select_fn = t_config["custom_globals_select"]
            select_globals = {
                **select_globals,
                **custom_globals_select_fn(ht),
            }
        if "global_name" in t_config:
            global_name = t_config.get("global_name")
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
        else:
            raise ValueError(f"Invalid selects type: {type(selects)}")
    return select_fields


def get_final_ht_fields(
    ht: hl.Table,
    schema: Dict[str, List[str]] = FINALIZED_SCHEMA,
) -> Dict[str, List[str]]:
    """
    Get the final fields for the release HT.

    Create a dictionary of lists of fields that are in the `schema` and are
    present in the HT. If a field is not present in the HT, log a warning.

    :param ht: Hail Table.
    :param schema: Schema for the release HT.
    :return: Dict of final fields for the release HT.
    """
    final_fields = {"rows": [], "globals": []}
    for field in schema["rows"]:
        if field in ht.row:
            final_fields["rows"].append(field)
        else:
            logger.warning(f"Field {field} from schema not found in HT.")
    for field in schema["globals"]:
        if field in ht.globals:
            final_fields["globals"].append(field)
        else:
            logger.warning(f"Global {field} from schema not found in HT.")

    return final_fields


def get_ht(
    dataset: str,
    config: Dict[str, Dict[str, Any]],
    use_config_ht: bool = False,
    checkpoint: bool = False,
    test: bool = False,
    _intervals: Optional[Any] = None,
    base_ht_filter: Optional[hl.Table] = None,
) -> hl.Table:
    """
    Return the appropriate Hail table with selects applied.

    :param dataset: Hail Table to join.
    :param config: Dictionary with configuration options for each dataset. Expects the
        dataset name as the key and a dictionary of options as the value, where the
        options include some of the following: 'ht', 'path', 'select', 'field_name',
        'custom_select'.
    :param use_config_ht: Whether to use the 'ht' in the dataset config instead of
        reading in `path`. Default is False.
    :param checkpoint: Whether to checkpoint the Hail Table. Default is False.
    :param test: Whether call is for a test run. If True, the dataset will be filtered
        to PCSK9. Default is False.
    :param _intervals: Intervals for reading in hail Table. Used to optimize join.
    :param base_ht_filter: Optional Hail Table to filter on before joining. Default is
        None.
    :return: Hail Table with fields to select.
    """
    logger.info("Getting the %s dataset and its selected annotations...", dataset)
    dataset_config = config[dataset]

    if use_config_ht:
        ht = dataset_config["ht"]
    else:
        logger.info("Reading in %s", dataset)
        ht = hl.read_table(dataset_config["path"], _intervals=_intervals)

    if test:
        # Keep only PCSK9.
        ht = hl.filter_intervals(
            ht,
            [
                hl.parse_locus_interval(
                    "chr1:55039447-55064852", reference_genome="GRCh38"
                )
            ],
        )

    select_fields = get_select_fields(dataset_config.get("select"), ht)

    if "custom_select" in dataset_config:
        custom_select_fn = dataset_config["custom_select"]
        select_fields = {
            **select_fields,
            **custom_select_fn(ht, config=config),
        }

    if "field_name" in dataset_config:
        field_name = dataset_config.get("field_name")
        select_query = {field_name: hl.struct(**select_fields)}
    else:
        select_query = select_fields

    ht = ht.select(**select_query)

    if base_ht_filter:
        ht_key = list(ht.key)
        if list(base_ht_filter.key) != ht_key:
            base_ht_filter = base_ht_filter.key_by(*ht_key)
        ht = ht.semi_join(base_ht_filter)

    if checkpoint:
        ht = ht.checkpoint(new_temp_file(f"{dataset}.for_release", "ht"))

    return ht


def join_hts(
    base_table: str,
    tables: List[str],
    config: Dict[str, Dict[str, Any]],
    version: str,
    new_partition_percent: Optional[float] = None,
    new_n_partitions: Optional[float] = None,
    checkpoint_tables: bool = False,
    track_included_datasets: bool = False,
    use_annotate: bool = False,
    test: bool = False,
    environment: str = "rwb",
) -> hl.Table:
    """
    Outer join a list of Hail Tables.

    :param base_table: Dataset to use for interval partitioning.
    :param tables: List of tables to join.
    :param config: Dictionary with configuration options for each dataset. Expects the
        dataset name (matching values in `tables`) as the key and a dictionary of
        options as the value, where the options include some of the following: 'ht',
        'path', 'select', 'field_name', 'custom_select', 'select_globals',
        'custom_globals_select', 'global_name'.
    :param version: Release version.
    :param new_partition_percent: Percent of base_table partitions used for final
        release Hail Table.
    :param new_n_partitions: Number of partitions for final release Hail Table.
    :param checkpoint_tables: Whether to checkpoint the tables before join. Default is
        False.
    :param track_included_datasets: Whether to track the datasets included in the
        release. This is used to update the included_datasets json file. Default is
        False.
    :param use_annotate: Whether to use annotate instead of join. Default is False.
    :param test: Whether this is for a test run. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Hail Table with datasets joined.
    """
    logger.info(
        "Reading in %s to determine partition intervals for efficient join",
        base_table,
    )

    partition_intervals = None
    if new_n_partitions or new_partition_percent:
        base_ht = get_ht(
            dataset=base_table,
            config=config,
            use_config_ht=True,
            test=test,
        )
        new_n_partitions = new_n_partitions or base_ht.n_partitions()
        if new_partition_percent:
            new_n_partitions = int(base_ht.n_partitions() * new_partition_percent)

        partition_intervals = base_ht._calculate_new_partitions(new_n_partitions)

    base_ht = get_ht(
        dataset=base_table,
        config=config,
        checkpoint=checkpoint_tables,
        _intervals=partition_intervals,
        test=test,
    )

    # Remove base_table from tables if it is in the list and add it to the front.
    if base_table in tables:
        tables.remove(base_table)

    logger.info("Joining datasets: %s...", tables)
    hts = [base_ht] + [
        get_ht(
            dataset=table,
            config=config,
            # There is no single path for insilico predictors, so we need to handle
            # this case separately.
            use_config_ht=True if table in ["in_silico", "exomes_an"] else False,
            checkpoint=checkpoint_tables,
            _intervals=partition_intervals,
            test=test,
            base_ht_filter=base_ht.select(),
        )
        for table in tables
    ]

    if use_annotate:
        joined_ht = reduce(
            lambda joined_ht, ht: joined_ht.annotate(
                **ht.index([*[joined_ht[k] for k in list(ht.key)]])
            ),
            hts,
        )
        hts = [g_ht.index_globals() for g_ht in hts]
        global_struct = reduce(lambda g, ann_g: g.annotate(**ann_g), hts)
        joined_ht = joined_ht.select_globals(**global_struct)
    else:
        joined_ht = reduce((lambda joined_ht, ht: joined_ht.join(ht, "left")), hts)

    joined_ht = joined_ht.select_globals(
        **get_select_global_fields(joined_ht, config, [base_table] + tables)
    )

    if track_included_datasets:
        # Track the datasets we've added as well as the source paths.
        # If release HT is included in tables, read in the included datasets json
        # and update the keys to the path for any new tables
        included_datasets = {}
        if "release" in tables:
            with hl.utils.hadoop_open(
                included_datasets_json_path(
                    test=test,
                    release_version=hl.eval(config["release"]["ht"].version),
                    environment=environment,
                )
            ) as f:
                included_datasets = json.loads(f.read())

        included_datasets.update({t: config[t]["path"] for t in tables})

        with hl.utils.hadoop_open(
            included_datasets_json_path(
                test=test, release_version=version, environment=environment
            ),
            "w",
        ) as f:
            f.write(hl.eval(hl.json(included_datasets)))

    return joined_ht


def check_duplicate_rows_in_config_hts(
    config: dict,
) -> None:
    """
    Check all Hail Tables in a config dict for duplicate rows and raises a ValueError if any HT contains duplicate rows.

    :param config: Dictionary with configuration options for each dataset. Expects the
        dataset name as the key and a dictionary of options as the value. Options that include a
        value for 'ht' will be checked for duplicate rows.
    :param logger: Logger object.
    :return: None.
    """
    dup_errors = []

    for name, cfg in config.items():
        if "ht" not in cfg:
            continue

        ht = cfg["ht"]
        total_count = ht.count()
        distinct_count = ht.select().distinct().count()

        if distinct_count != total_count:
            dup_errors.append(
                f"HT {name} has {distinct_count} distinct rows but "
                f"{total_count} total rows."
            )
        else:
            logger.info(f"HT {name} has no duplicate rows.")

    if dup_errors:
        raise ValueError("\n".join(dup_errors))


def add_global_annotations(ht: hl.Table, version: str) -> hl.Table:
    """
    Add global version annotations to the release HT.

    Adds:
      - VEP metadata (gencode + MANE Select versions)
      - Tool versions (dbSNP, SIFT, PolyPhen)
      - VRS-related versioning
      - Release date and version
      - Frequency README
    :param ht: Hail Table.
    :param version: Version of the release.
    :return: Hail Table with global version annotations.
    """
    return ht.annotate_globals(
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
            vrs_schema_version=VRS_SCHEMA_VERSION,
            vrs_python_version=VRS_PYTHON_VERSION,
            seqrepo_version=SEQREPO_VERSION,
        ),
        date=datetime.now().isoformat(),
        version=version,
        frequency_README=FREQUENCY_README.format(""),
    )


def main(args):
    """Create release ht."""
    if args.rwb:
        environment = "rwb"
        hl.init(
            log="/home/jupyter/workspaces/gnomadproduction/create_release_ht.log",
            tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
        )
    else:
        environment = "batch"
        hl.init(
            tmp_dir=f"gs://{GNOMAD_TMP_BUCKET}-4day",
            log="create_release_ht.log",
        )
        # TODO: Add machine configurations for Batch.
    hl.default_reference("GRCh38")
    overwrite = args.overwrite
    test = args.test

    logger.info(
        "Getting config for release HT to check Tables for duplicate variants..."
    )
    config = get_config(
        release_exists=args.release_exists, test=test, environment=environment
    )
    check_duplicate_rows_in_config_hts(config)

    logger.info(f"Creating release HT...")
    ht = join_hts(
        args.base_table,
        args.tables_for_join,
        config,
        args.version,
        new_partition_percent=args.new_partition_percent,
        track_included_datasets=True,
        test=test,
    )

    # TODO: Check if this filter is still needed.
    # Filter out chrM, AS_lowqual sites (these sites are dropped in the final_filters HT
    # so will not have information in `filters`) and AC_raw == 0.
    logger.info("Filtering out chrM, AS_lowqual, and AC_raw == 0 sites...")
    ht = hl.filter_intervals(ht, [hl.parse_locus_interval("chrM")], keep=False)
    ht = ht.filter(hl.is_defined(ht.filters) & (ht.freq[1].AC > 0))

    logger.info("Finalizing the release HT globals...")
    ht = add_global_annotations(ht, version=args.version)

    # Organize the fields in the release HT to match the order of FINALIZED_SCHEMA when
    # the fields are present in the HT.
    final_fields = get_final_ht_fields(ht)
    ht = ht.select(*final_fields["rows"]).select_globals(*final_fields["globals"])

    output_path = (
        f"{qc_temp_prefix(environment=environment)}release/gnomad.genomes.sites.test.{datetime.today().strftime('%Y-%m-%d')}.ht"
        if test
        else release_sites(environment=environment).path
    )

    logger.info(f"Writing out release HT to %s", output_path)
    ht = ht.naive_coalesce(args.n_partitions).checkpoint(
        output_path,
        args.overwrite,
    )

    logger.info("Final release HT schema:")
    ht.describe()

    logger.info("Final variant count: %d", ht.count())


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rwb",
        help="Run the script in RWB environment.",
        action="store_true",
    )
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
        "--n-partitions",
        help="Number of partitions to naive coalesce the release Table to.",
        type=int,
        default=10000,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    main(args)
