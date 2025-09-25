"""
Merge gnomAD v4 metadata with AoU metadata to create gnomAD v5 project metadata.

NOTE:
Use the util_functions.ipynb notebook with function restart_kernel_with_gnomad_packages
to clone the gnomad_qc repository, update the gnomAD package, and restart the kernel. If
this is not done, the gnomad_qc imports will fail below.
"""

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
        "GNOMAD_METHODS_BRANCH": "main",
    },
)
# Restart the kernel -- this needs to be run here, in the open notebook.
display(Javascript("Jupyter.notebook.kernel.restart()"))  # noqa

import logging
import os
from functools import reduce
from itertools import combinations
from typing import Dict, Set

import hail as hl
import pandas as pd

from gnomad_qc.v4.resources.meta import meta as meta_v4
from gnomad_qc.v5.resources.basics import (
    add_project_prefix_to_sample_collisions,
    get_checkpoint_path,
)
from gnomad_qc.v5.resources.meta import project_meta, sample_id_collisions

hl.default_reference(new_default_reference="GRCh38")

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("merge_project_meta")
logger.setLevel(logging.INFO)

AOU_AUXILIARY_DATA_BUCKET = (
    "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/aux"
)
FINAL_SCHEMA_FIELDS_AND_TYPES = {
    "s": hl.tstr,
    "age": hl.tint,
    "contamination_freemix": hl.tfloat,
    "contamination_charr": hl.tfloat,
    "chimeras_rate": hl.tfloat,
    "bases_over_20x_coverage": hl.tfloat,
    "mean_depth": hl.tfloat,
    "median_insert_size": hl.tfloat,
    "callrate": hl.tfloat,
    "hard_filters": hl.tset(hl.tstr),
    "sex_karyotype": hl.tstr,
    "gen_anc": hl.tstr,
    "outlier_filters": hl.tset(hl.tstr),
    "releasable": hl.tbool,
    "release": hl.tbool,
    "sample_source": hl.tstr,
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
                        "age": "project_meta.age",
                        "contamination_freemix": "bam_metrics.contam_rate",
                        "contamination_charr": (
                            "hard_filter_metrics.mean_AB_snp_biallelic"
                        ),
                        "chimeras_rate": "bam_metrics.chimeras_rate",
                        "bases_over_20x_coverage": "bam_metrics.target_bases_20x_rate",
                        "mean_depth": "hard_filter_metrics.chr20_mean_dp",
                        "median_insert_size": "bam_metrics.median_insert_size",
                        "callrate": "hard_filter_metrics.sample_qc_mt_callrate",
                        "hard_filters": "sample_filters.hard_filters",
                        "sex_karyotype": "sex_imputation.sex_karyotype",
                        "gen_anc": "population_inference.pop",
                        "outlier_filters": "sample_filters.qc_metrics_filters",
                        "releasable": "project_meta.releasable",
                        "release": "release",
                    },
                }
            ],
            "genomes": [
                {
                    "ht": meta_v4(data_type="genomes").ht(),
                    "file_format": "ht",
                    "field_mappings": {
                        "s": "s",
                        "age": "project_meta.age",
                        # NOTE: Some samples have age_alt, which is the mean age of the
                        # age_bin field.
                        "age_alt": "project_meta.age_alt",
                        # TODO: Did we recalculate charr on v4 genomes? Determine/look
                        # into cost and then decide if we want to recalculate.
                        "contamination_freemix": "bam_metrics.freemix",
                        "chimeras_rate": "bam_metrics.pct_chimeras",
                        "bases_over_20x_coverage": "bam_metrics.pct_bases_20x",
                        "mean_depth": "sex_imputation.chr20_mean_dp",
                        "median_insert_size": "bam_metrics.median_insert_size",
                        "hard_filters": "sample_filters.hard_filters",
                        "sex_karyotype": "sex_imputation.sex_karyotype",
                        "gen_anc": "population_inference.pop",
                        "outlier_filters": "sample_filters.qc_metrics_filters",
                        "releasable": "project_meta.releasable",
                        "release": "release",
                        "research_project_key": "project_meta.research_project_key",
                    },
                }
            ],
        },
        "aou": {
            "genomes": [
                {
                    "path": f"{AOU_AUXILIARY_DATA_BUCKET}/qc/genomic_metrics.tsv",
                    "file_format": "tsv",
                    "field_mappings": {
                        "s": "research_id",
                        "sex_karyotype": "dragen_sex_ploidy",
                        "mean_depth": "mean_coverage",
                        "bases_over_20x_coverage": "genome_coverage",
                        # TODO: Do we want dragen estimates or verifybamid? They show
                        # the two metrics correlate in the QC report
                        "contamination_freemix": "verify_bam_id2_contamination",
                        # TODO: Check within hard filters if there is a difference --
                        # annotate back on during outlier detection.
                        "sample_source": "sample_source",
                    },
                },
                {
                    "path": f"{AOU_AUXILIARY_DATA_BUCKET}/ancestry/ancestry_preds.tsv",
                    "file_format": "tsv",
                    "field_mappings": {
                        "s": "research_id",
                        # AoU table has two columns with genetic ancestry group assignments,
                        # and this column includes the "oth" (remaining) label in the assignments.
                        # https://support.researchallofus.org/hc/en-us/articles/29475228181908-How-the-All-of-Us-Genomic-data-are-organized#h_01GY7QYC017WWMXFB0YTT233G1
                        "gen_anc": "ancestry_pred_other",
                    },
                },
                {
                    "path": f"{AOU_AUXILIARY_DATA_BUCKET}/qc/all_samples.tsv",
                    "file_format": "tsv",
                    "field_mappings": {
                        "s": "s",
                        "outlier_filters": "qc_metrics_filters",
                    },
                },
                {
                    "ht": pull_aou_age(),
                    "file_format": "ht",
                    "field_mappings": {
                        "s": "research_id",
                        "age": "age",
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

    .. note ::

        There were 897 gnomAD samples that withdrew consent during v5 production.
        Because of the timing of the confirmation of the consent change, the releasable
        field is being reannotated after the original creation of the v5 metadata
        The original metadata was used through generating the release relateds to drop
        HT. However, this use was of no consequence to the outcome of sample QC modules
        as the releasable field is only considered during the final `release` annotation
        in sample QC metadata creation.

    :param ht: Input Hail Table with metadata.
    :param project: Project identifier (e.g., "gnomad", "aou").
    :param data_type: Data type ("genomes" or "exomes").
    :return: Filtered Hail Table containing only fields required for the final schema.
    """
    # gnomAD v4 genome metadata has some samples with age_alt, which is the mean age of
    # the age_bin field. This creates a single age field. Also, 897 samples were consented
    # to be in gnomAD v4 but not v5, so we need to set their releasable field to False.
    if project == "gnomad" and data_type == "genomes":
        ht = ht.annotate(age=hl.if_else(hl.is_defined(ht.age), ht.age, ht.age_alt))
        ht = ht.annotate(
            releasable=hl.if_else(
                (ht.research_project_key == "RP-1061")
                | (ht.research_project_key == "RP-1411"),
                False,
                ht.releasable,
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


def get_sample_collisions(meta_hts: Dict[str, hl.Table]) -> hl.Table:
    """
    Get sample ID collisions between projects.

    :param meta_hts: Project metadata tables.
    :return: Table of sample IDs that exist in more than one project.
    """
    # Get all possible pairs of projects
    project_pairs = list(combinations(meta_hts.keys(), 2))

    sample_collisions = hl.Table.parallelize([], schema=hl.tstruct(s=hl.tstr), key="s")

    for p1, p2 in project_pairs:
        logger.info(f"project 1 is {p1}, project 2 is {p2}")
        collisions_ht = meta_hts[p1].filter(hl.is_defined(meta_hts[p2][meta_hts[p1].s]))
        sample_collisions = sample_collisions.union(collisions_ht.select())

    return sample_collisions


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
                # The left join ensures we keep only AoU samples with srWGS data but
                # requires all samples to be in the first table, which they are.
                data_type_ht = reduce(
                    (lambda joined_ht, ht: joined_ht.join(ht, how="left")),
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
            get_checkpoint_path(f"{project}_meta", mt=False, environment="rwb"),
            overwrite=True,
        )
        meta_hts[project] = project_ht.key_by("s")

    logger.info("Determining if there are any sample ID collisions between projects...")
    sample_collisions = get_sample_collisions(meta_hts)
    sample_collisions = sample_collisions.checkpoint(
        sample_id_collisions.path, overwrite=True
    )
    logger.info(
        f"Samples that appear in more than one project: {hl.eval(sample_collisions.aggregate(hl.agg.collect_as_set(sample_collisions.s)))}"
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
