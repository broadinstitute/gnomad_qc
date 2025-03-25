"""
Merge gnomAD v4 metadata with AoU metadata to create gnomAD v5 project metadata.

NOTE:
Use the util_functions.ipynb notebook with function restart_kernel_with_gnomad_package
to clone the gnomad_qc repository, update the gnomAD package, and restart the kernel.If
this is not done, the gnomad_qc imports will fail below.
"""

import logging
import os
from functools import reduce
from itertools import combinations
from typing import Dict, Set

import hail as hl
import pandas as pd
import papermill as pm
from IPython.display import Javascript, display

# This block clones the gnomad repos so you can import gnomad_qc and the most recent
# gnomad_methods commits. Change the branch parameters to your desired branch then
# run in the v5 workspace notebook.
pm.execute_notebook(  # noqa
    "utils_restart_kernel_with_gnomad_packages.ipynb",
    "utils_restart_notebook_output.ipynb",
    parameters={
        "GNOMAD_QC_BRANCH": "main",
        "GNOMAD_METHODS": "main",
    },
)
# Restart the kernel -- this needs to be run here, in the open notebook.
display(Javascript("Jupyter.notebook.kernel.restart()"))  # noqa

from gnomad_qc.v4.resources.meta import meta as meta_v4
from gnomad_qc.v5.resources.basics import get_checkpoint_path
from gnomad_qc.v5.resources.meta import project_meta

hl.default_reference(new_default_reference="GRCh38")

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("merge_project_meta")
logger.setLevel(logging.INFO)

AOU_AUXILIARY_DATA_BUCKET = (
    "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux"
)
FINAL_SCHEMA_FIELDS_AND_TYPES = {
    "s": hl.tstr,
    "project_age": hl.tint,
    "project_contamination_freemix": hl.tfloat,
    "project_contamination_charr": hl.tfloat,
    "project_chimeras_rate": hl.tfloat,
    "project_bases_over_20x_coverage": hl.tfloat,
    "project_mean_depth": hl.tfloat,
    "project_median_insert_size": hl.tfloat,
    "project_callrate": hl.tfloat,
    "project_hard_filters": hl.tset(hl.tstr),
    "project_sex_karyotype": hl.tstr,
    "project_gen_anc": hl.tstr,
    "project_outlier_filters": hl.tset(hl.tstr),
    "project_releasable": hl.tbool,
    "project_release": hl.tstr,
    "project_sample_source": hl.tstr,
    "data_type": hl.tstr,
    "project": hl.tstr,
}

EXPR_TO_TYPE = {
    hl.tstr: hl.str,
    hl.tint: hl.int,
    hl.tfloat: hl.float,
    hl.bool: hl.tbool,
}


def get_meta_config() -> Dict[str, Dict[str, hl.expr.Expression]]:
    """
    Get project metadata configuration.

    :return: Project metadata dictionary configuration.
    """
    return {
        "gnomad": {
            "exomes": [
                {
                    "ht": meta_v4(data_type="exomes").ht(),
                    "file_format": "ht",
                    "field_mappings": {
                        "s": "s",
                        "project_age": "project_meta.age",
                        "project_contamination_freemix": "bam_metrics.contam_rate",
                        "project_contamination_charr": (
                            "hard_filter_metrics.mean_AB_snp_biallelic"
                        ),
                        "project_chimeras_rate": "bam_metrics.chimeras_rate",
                        "project_bases_over_20x_coverage": (
                            "bam_metrics.target_bases_20x_rate"
                        ),
                        "project_mean_depth": "hard_filter_metrics.chr20_mean_dp",
                        "project_median_insert_size": "bam_metrics.median_insert_size",
                        "project_callrate": "hard_filter_metrics.sample_qc_mt_callrate",
                        "project_hard_filters": "sample_filters.hard_filters",
                        "project_sex_karyotype": "sex_imputation.sex_karyotype",
                        "project_gen_anc": "population_inference.pop",
                        "project_qc_metrics_filters": (
                            "sample_filters.qc_metrics_filters"
                        ),
                        "project_releasable": "project_meta.releasable",
                        "project_release": "release",
                    },
                }
            ],
            "genomes": [
                {
                    "ht": meta_v4(data_type="genomes").ht(),
                    "file_format": "ht",
                    "field_mappings": {
                        "s": "s",
                        "project_age": "project_meta.age",
                        # NOTE: Some samples have age_alt, which is the mean age of the
                        # age_bin field.
                        "project_age_alt": "project_meta.age_alt",
                        # TODO: Did we recalculate charr on v4 genomes? Determine/look
                        # into cost and then decide if we want to recalculate.
                        "project_contamination_freemix": "bam_metrics.freemix",
                        "project_chimeras_rate": "bam_metrics.pct_chimeras",
                        "project_bases_over_20x_coverage": "bam_metrics.pct_bases_20x",
                        "project_mean_depth": "sex_imputation.chr20_mean_dp",
                        "project_median_insert_size": "bam_metrics.median_insert_size",
                        "project_hard_filters": "sample_filters.hard_filters",
                        "project_sex_karyotype": "sex_imputation.sex_karyotype",
                        "project_gen_anc": "population_inference.pop",
                        "project_outlier_filters": "sample_filters.qc_metrics_filters",
                        "project_releasable": "project_meta.releasable",
                        "project_release": "release",
                    },
                }
            ],
        },
        "aou": {
            "genomes": [
                {
                    "path": f"{AOU_AUXILIARY_DATA_BUCKET}/ancestry/ancestry_preds.tsv",
                    "file_format": "tsv",
                    "field_mappings": {
                        "s": "research_id",
                        "project_gen_anc": "ancestry_pred",
                    },
                },
                {
                    "path": f"{AOU_AUXILIARY_DATA_BUCKET}/qc/genomic_metrics.tsv",
                    "file_format": "tsv",
                    "field_mappings": {
                        "s": "research_id",
                        "project_sex_karyotype": "dragen_sex_ploidy",
                        "project_mean_depth": "mean_coverage",
                        "project_bases_over_20x_coverage": "genome_coverage",
                        # TODO: Do we want dragen estimates or verifybamid? They show
                        # the two metrics correlate in the QC report.all
                        "project_contamination_freemix": "verify_bam_id2_contamination",
                        # TODO: Check within hard filters if there is a difference --
                        # annotate back on during outlier detection.
                        "project_sample_source": "sample_source",
                    },
                },
                {
                    "path": f"{AOU_AUXILIARY_DATA_BUCKET}/qc/flagged_samples.tsv",
                    "file_format": "tsv",
                    "field_mappings": {
                        "s": "s",
                        "project_outlier_filters": "qc_metrics_filters",
                    },
                },
                {
                    "ht": pull_aou_age(),
                    "file_format": "ht",
                    "field_mappings": {
                        "s": "research_id",
                        "project_age": "age",
                    },
                },
            ],
        },
    }


def pull_aou_age() -> hl.Table:
    """
    Pull age information from AoU metadata.

    :return: Table with AoU age information.
    """
    query = f"""
    SELECT
        CAST(person.person_id as string) as research_id,
        CAST(FLOOR(DATE_DIFF(DATE(CURRENT_DATE),DATE(person.birth_datetime), DAY)/365.25) as int) AS age,
    FROM  `{os.getenv('WORKSPACE_CDR')}.person` person
    """
    age_df = pd.read_gbq(query, dialect="standard", progress_bar_type="tqdm_notebook")

    return hl.Table.from_pandas(age_df, key="research_id")


def select_only_final_fields(
    ht: hl.Table, project: hl.str, data_type: hl.str
) -> hl.Table:
    """
    Extract only fields needed for the final metadata table.

    Filter metadata inputs to include only relevant fields since not all inputs have all
    final fields.

    :param ht: Input Hail Table with metadata
    :param project: Project identifier (e.g., "gnomad", "aou")
    :param data_type: Data type ("genomes" or "exomes")
    :return: Filtered Hail Table containing only fields required for the final schema
    """
    # gnomAD v4 genome metadata has some samples with age_alt, which is the mean age of
    # the age_bin field. This creates a single age field.
    if project == "gnomad" and data_type == "genomes":
        ht = ht.annotate(
            project_age=hl.if_else(
                hl.is_defined(ht.project_age), ht.project_age, ht.project_age_alt
            )
        )

    present_fields = {
        field for field in ht.row if field in FINAL_SCHEMA_FIELDS_AND_TYPES.keys()
    }
    return ht.select(*present_fields - set(ht.key.dtype.fields))


def update_ht_to_final_schema(ht: hl.Table) -> hl.Table:
    """
    Update table to final schema.

    If a field exists within the table being updated, check that the type matches, if
    not, cast to the type in the joint field schema. If it does not exist, annotate
    with missing value.

    :param ht: Table to update to final schema.
    :return: Table with final schema.
    """
    final_fields = {
        k: v
        for k, v in FINAL_SCHEMA_FIELDS_AND_TYPES.items()
        if k not in set(ht.key.dtype.fields)
    }
    annotations = {}
    for field, field_type in final_fields.items():
        if field in ht.row:
            if ht[field].dtype != field_type:
                annotations[field] = EXPR_TO_TYPE[field_type](ht[field])
        else:
            annotations[field] = hl.missing(field_type)

    ht = ht.annotate(**annotations)

    return ht.select(*final_fields)


def get_sample_collisions(meta_hts: Dict[str, hl.Table]) -> Set[str]:
    """
    Get sample ID collisions between projects.

    :param meta_hts: Project metadata tables.
    :return: List of sample IDs that exist in more than one project.
    """
    # Get all possible pairs of projects
    project_pairs = list(combinations(meta_hts.keys(), 2))
    sample_collisions = hl.empty_set(hl.tstr)

    for p1, p2 in project_pairs:
        logger.info(f"project 1 is {p1}, project 2 is {p2}")
        collisions_ht = meta_hts[p1].filter(hl.is_defined(meta_hts[p2][meta_hts[p1].s]))
        collisions = collisions_ht.aggregate(hl.agg.collect_as_set(collisions_ht.s))
        sample_collisions = sample_collisions.union(collisions)

    return sample_collisions


def add_project_prefix_to_sample_collisions(
    ht: hl.Table,
    project: str,
    sample_collisions: Set[str],
    sample_id_field: str = "s",
) -> hl.Table:
    """
    Add project prefix to sample IDs that exist in multiple projects.

    :param ht: Table to add project prefix to sample IDs.
    :param project: Project name.
    :param sample_collisions: Set of sample IDs that exist in multiple projects.
    :return: Table with project prefix added to sample IDs.
    """
    ht = ht.key_by()
    ht = ht.annotate(
        **{
            f"{sample_id_field}": hl.if_else(
                hl.literal(sample_collisions).contains(ht["s"]),
                hl.delimit([ht.project, ht[sample_id_field]], "_"),
                ht[sample_id_field],
            )
        }
    )
    return ht.key_by(sample_id_field)


def main():
    """Create v5 project metadata."""
    meta_hts = {}
    for project, project_metadata in get_meta_config().items():
        project_hts = []
        logger.info("Processing the %s project's metadata...", project)
        for data_type in project_metadata.keys():
            data_type_hts = []
            logger.info("Importing %s's %s metadata ...", project, data_type)
            for file in project_metadata[data_type]:

                if file["file_format"] == "ht":
                    ht = file["ht"]

                elif file["file_format"] == "tsv":
                    file_types = {
                        input_field: FINAL_SCHEMA_FIELDS_AND_TYPES[final_field]
                        for final_field, input_field in file["field_mappings"].items()
                    }
                    ht = hl.import_table(file["path"], types=file_types)

                else:
                    raise ValueError(f"File format {file['file_type']} not supported")

                ht = ht.flatten().select_globals()
                ht = ht.select(
                    **{
                        new_field: ht[field]
                        for new_field, field in file["field_mappings"].items()
                    }
                )
                ht = ht.key_by("s")
                logger.info("Selecting out fields for the final schema...")
                ht = select_only_final_fields(ht, project, data_type)
                data_type_hts.append(ht)

            if len(data_type_hts) > 1:
                logger.info(
                    "Combining the %s metadata files from %s...", data_type, project
                )
                data_type_ht = reduce(
                    (lambda joined_ht, ht: joined_ht.join(ht, how="outer")),
                    data_type_hts,
                )
            else:
                data_type_ht = data_type_hts[0]

            data_type_ht = data_type_ht.annotate(data_type=data_type)
            logger.info("Updating the %s metadata to the final schema...", data_type)
            data_type_ht = update_ht_to_final_schema(data_type_ht)
            project_hts.append(data_type_ht.key_by("s"))

        logger.info(
            "Combining the %s project's files into a single metadata file...", project
        )
        project_ht = reduce((lambda joined_ht, ht: joined_ht.union(ht)), project_hts)
        project_ht = project_ht.annotate(project=project)
        logger.info("Updating the %s metadata to the final schema...", project)
        project_ht = update_ht_to_final_schema(project_ht)
        project_ht = project_ht.checkpoint(
            get_checkpoint_path(f"{project}_meta", mt=False), overwrite=True
        )
        meta_hts[project] = project_ht

    sample_collisions = get_sample_collisions(meta_hts)
    logger.info(
        f"Samples that appear in more than one project: {hl.eval(sample_collisions)}"
    )

    for project in get_meta_config().keys():
        meta_hts[project] = add_project_prefix_to_sample_collisions(
            meta_hts[project],
            project=project,
            sample_collisions=sample_collisions,
        )

    logger.info("Combining all project metadata files into a single metadata file...")
    meta_ht = reduce((lambda joined_ht, ht: joined_ht.union(ht)), meta_hts.values())

    logger.info("Writing combined project metadata to %s...", project_meta.path)
    meta_ht.write(project_meta.path, overwrite=True)


if __name__ == "__main__":
    main()
