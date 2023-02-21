"""Script containing generic resources."""
import logging
from typing import Dict, List, Optional

import hail as hl
from gnomad.resources.resource_utils import (
    TableResource,
    VariantDatasetResource,
    VersionedTableResource,
    VersionedVariantDatasetResource,
)
from gnomad.utils.file_utils import check_file_exists_raise_error

from gnomad_qc.v4.resources.constants import CURRENT_VERSION
from gnomad_qc.v4.resources.meta import meta

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)


# Note: Unlike previous versions, the v4 resource directory uses a general format of gs://gnomad/v4.0/<module>/<exomes_or_genomes>/. # noqa
def get_gnomad_v4_vds(
    split: bool = False,
    remove_hard_filtered_samples: bool = True,
    remove_hard_filtered_samples_no_sex: bool = False,
    release_only: bool = False,
    test: bool = False,
    n_partitions: int = None,
) -> hl.vds.VariantDataset:
    """
    Get gnomAD v4 data with desired filtering and metadata annotations.

    :param split: Perform split on VDS - Note: this will perform a split on the VDS rather than grab an already split VDS
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after hard filtering is complete)
    :param remove_hard_filtered_samples_no_sex: Whether to remove samples that failed non sex inference hard filters (only relevant after pre-sex imputation hard filtering is complete)
    :param release_only: Whether to filter the VDS to only samples available for release (can only be used if metadata is present)
    :param test: Whether to use the test VDS instead of the full v4 VDS
    :param n_partitions: Optional argument to read the VDS with a specific number of partitions
    :return: gnomAD v4 dataset with chosen annotations and filters
    """
    if remove_hard_filtered_samples and remove_hard_filtered_samples_no_sex:
        raise ValueError(
            "Only one of 'remove_hard_filtered_samples' or"
            " 'remove_hard_filtered_samples_no_sex' can be set to True."
        )

    if test:
        gnomad_v4_resource = gnomad_v4_testset
    else:
        gnomad_v4_resource = gnomad_v4_genotypes

    if n_partitions:
        vds = hl.vds.read_vds(gnomad_v4_resource.path, n_partitions=n_partitions)
    else:
        vds = gnomad_v4_resource.vds()

    # Count current number of samples in the VDS.
    n_samples = vds.variant_data.count_cols()

    # Remove 75 withdrawn UKB samples (samples with withdrawn consents for application
    # 31063 on 02/22/2022).
    ukb_application_map_ht = ukb_application_map.ht()
    withdrawn_ukb_samples = ukb_application_map_ht.filter(
        ukb_application_map_ht.withdraw
    ).s.collect()

    # Remove 43 samples that are known to be on the pharma's sample remove list. See:
    # https://github.com/broadinstitute/ukbb_qc/blob/70c4268ab32e9efa948fe72f3887e1b81d8acb46/ukbb_qc/resources/basics.py#L308.
    dups_ht = ukb_known_dups.ht()
    ids_to_remove = dups_ht.aggregate(hl.agg.collect(dups_ht["Sample Name - ID1"]))

    # Remove 27 fully duplicated IDs (same exact name for 's' in the VDS). See:
    # https://github.com/broadinstitute/ukbb_qc/blob/70c4268ab32e9efa948fe72f3887e1b81d8acb46/ukbb_qc/resources/basics.py#L286.
    # Confirmed that the col_idx of the 27 dup samples in the original UKB MT match
    # the col_idx of the dup UKB samples in the VDS.
    dup_ids = []
    with hl.hadoop_open(ukb_dups_idx_path, "r") as d:
        for line in d:
            line = line.strip().split("\t")
            dup_ids.append(f"{line[0]}_{line[1]}")

    def _remove_ukb_dup_by_index(
        mt: hl.MatrixTable, dup_ids: hl.expr.ArrayExpression
    ) -> hl.MatrixTable:
        """
        Remove UKB samples with exact duplicate names based on column index.

        :param mt: MatrixTable of either the variant data or reference data of a VDS
        :param dup_ids: ArrayExpression of UKB samples to remove in format of <sample_name>_<col_idx>
        :return: MatrixTable of UKB samples with exact duplicate names removed based on column index
        """
        mt = mt.add_col_index()
        mt = mt.annotate_cols(new_s=hl.format("%s_%s", mt.s, mt.col_idx))
        mt = mt.filter_cols(~dup_ids.contains(mt.new_s)).drop("new_s", "col_idx")
        return mt

    dup_ids = hl.literal(dup_ids)
    vd = _remove_ukb_dup_by_index(vds.variant_data, dup_ids)
    rd = _remove_ukb_dup_by_index(vds.reference_data, dup_ids)
    vds = hl.vds.VariantDataset(rd, vd)

    # Filter withdrawn samples from the VDS.
    withdrawn_ids = withdrawn_ukb_samples + ids_to_remove + hl.eval(dup_ids)
    withdrawn_ids = [i for i in withdrawn_ids if i is not None]

    logger.info("Total number of UKB samples to exclude: %d", len(withdrawn_ids))

    with hl.hadoop_open(all_ukb_samples_to_remove, "w") as d:
        for sample in withdrawn_ids:
            d.write(sample + "\n")

    withdrawn_ht = hl.import_table(all_ukb_samples_to_remove, no_header=True).key_by(
        "f0"
    )

    vds = hl.vds.filter_samples(vds, withdrawn_ht, keep=False, remove_dead_alleles=True)

    # Log number of UKB samples removed from the VDS.
    n_samples_after_exclusion = vds.variant_data.count_cols()
    n_samples_removed = n_samples - n_samples_after_exclusion

    logger.info(
        "Total number of UKB samples removed from the VDS: %d", n_samples_removed
    )

    if remove_hard_filtered_samples or remove_hard_filtered_samples_no_sex:
        if test:
            meta_ht = gnomad_v4_testset_meta.ht()
            if remove_hard_filtered_samples_no_sex:
                hard_filter_expr = meta_ht.rand_sampling_meta.hard_filters_no_sex
            else:
                hard_filter_expr = meta_ht.rand_sampling_meta.hard_filters
            meta_ht = meta_ht.filter(hl.len(hard_filter_expr) == 0)
            vds = hl.vds.filter_samples(vds, meta_ht)
        else:
            from gnomad_qc.v4.resources.sample_qc import (
                hard_filtered_samples,
                hard_filtered_samples_no_sex,
            )

            if remove_hard_filtered_samples:
                hard_filter_ht = hard_filtered_samples.versions[CURRENT_VERSION].ht()
            else:
                hard_filter_ht = hard_filtered_samples_no_sex.versions[
                    CURRENT_VERSION
                ].ht()
            vds = hl.vds.filter_samples(vds, hard_filter_ht, keep=False)

    if release_only:
        if test:
            meta_ht = gnomad_v4_testset_meta.ht()
        else:
            meta_ht = meta.versions[CURRENT_VERSION].ht()
        meta_ht = meta_ht.filter(meta_ht.release)
        vds = hl.vds.filter_samples(vds, meta_ht)

    if split:
        vmt = vds.variant_data
        vmt = vmt.annotate_rows(
            n_unsplit_alleles=hl.len(vmt.alleles),
            mixed_site=(hl.len(vmt.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(vmt.alleles[0], a), vmt.alleles[1:])
            & hl.any(lambda a: hl.is_snp(vmt.alleles[0], a), vmt.alleles[1:]),
        )
        vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
        vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    return vds


_gnomad_v4_genotypes = {
    "4.0": VariantDatasetResource("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds")
}

gnomad_v4_genotypes = VersionedVariantDatasetResource(
    CURRENT_VERSION,
    _gnomad_v4_genotypes,
)

# v4 test dataset VDS.
gnomad_v4_testset = VariantDatasetResource(
    "gs://gnomad/v4.0/raw/exomes/testing/gnomad_v4.0_test.vds"
)
gnomad_v4_testset_meta = TableResource(
    "gs://gnomad/v4.0/raw/exomes/testing/gnomad_v4.0_meta.ht"
)


# UKB data resources.
def _ukb_root_path() -> str:
    """
    Retrieve the path to the UKB data directory.

    :return: String representation of the path to the UKB data directory
    """
    return "gs://gnomad/v4.0/ukbb"


# List of samples to exclude from QC due to withdrawn consents.
# Application 26041 is from the 08/09/2021 list and application 31063 is from the
# 02/22/2022 list.
# These were originally imported as CSVs and then keyed by their respective eids.
ukb_excluded = VersionedTableResource(
    default_version="31063_20220222",
    versions={
        "26041_20210809": TableResource(path=f"{_ukb_root_path()}/w26041_20210809.ht"),
        "31063_20220222": TableResource(
            path=f"{_ukb_root_path()}/ukb31063.withdrawn_participants.20220222.ht"
        ),
    },
)

# UKB map of exome IDs to array sample IDs (application ID: 26041).
ukb_array_sample_map = TableResource(f"{_ukb_root_path()}/array_sample_map_freeze_7.ht")

# UKB full mapping file of sample ID and application IDs 26041, 31063, and 48511.
ukb_application_map = TableResource(
    f"{_ukb_root_path()}/ukbb_application_id_mappings.ht"
)


# Samples known to be on the pharma partners' remove lists.
# All 44 samples are marked as "unresolved duplicates" by the pharma partners.
ukb_known_dups = TableResource(f"{_ukb_root_path()}/pharma_known_dups_7.ht")

# Samples with duplicate names in the VDS and their column index.
# 27 samples to remove based on column index.
ukb_dups_idx_path = f"{_ukb_root_path()}/dup_remove_idx_7.tsv"

# Final list of UKB samples to remove (duplicates that were removed will have their
# original index appended to the sample name).
all_ukb_samples_to_remove = f"{_ukb_root_path()}/all_ukbb_samples_to_remove.txt"

# UKB f-stat sites Table with UKB allele frequencies.
ukb_f_stat = TableResource(f"{_ukb_root_path()}/f_stat_sites.ht")


def qc_temp_prefix(version: str = CURRENT_VERSION) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return
    :return: Path to bucket with temporary QC data
    """
    return f"gs://gnomad-tmp/gnomad.exomes.v{version}.qc_data/"


def get_checkpoint_path(
    name: str, version: str = CURRENT_VERSION, mt: bool = False
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable
    :param version: Version of annotation path to return
    :param bool mt: Whether path is for a MatrixTable, default is False
    :return: Output checkpoint path
    """
    return f'{qc_temp_prefix(version)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(name: str, version: str = CURRENT_VERSION) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file
    :param version: Version of annotation path to return
    :return: Output log path
    """
    return f"{qc_temp_prefix(version)}{name}.log"


def add_meta(
    mt: hl.MatrixTable, version: str = CURRENT_VERSION, meta_name: str = "meta"
) -> hl.MatrixTable:
    """
    Add metadata to MT in 'meta_name' column.

    :param mt: MatrixTable to which 'meta_name' annotation should be added
    :param version: Version of metadata ht to use for annotations
    :return: MatrixTable with metadata added in a 'meta' column
    """
    mt = mt.annotate_cols(**{meta_name: meta.versions[version].ht()[mt.col_key]})

    return mt


def calling_intervals(
    interval_name: str, calling_interval_padding: int
) -> TableResource:
    """
    Return path to capture intervals Table.

    :param interval_name: One of 'ukb', 'broad', or 'intersection'
    :param calling_interval_padding: Padding around calling intervals. Available options are 0 or 50
    :return: Calling intervals resource
    """
    if interval_name not in {"ukb", "broad", "intersection"}:
        raise ValueError(
            "Calling interval name must be one of: 'ukb', 'broad', or 'intersection'!"
        )
    if calling_interval_padding not in {0, 50}:
        raise ValueError("Calling interval padding must be one of: 0 or 50 (bp)!")
    if interval_name == "ukb":
        return TableResource(
            f"gs://gnomad/resources/intervals/xgen_plus_spikein.Homo_sapiens_assembly38.targets.pad{calling_interval_padding}.interval_list.ht"
        )
    if interval_name == "broad":
        return TableResource(
            f"gs://gnomad/resources/intervals/hg38_v0_exome_calling_regions.v1.pad{calling_interval_padding}.interval_list.ht"
        )
    if interval_name == "intersection":
        return TableResource(
            f"gs://gnomad/resources/intervals/xgen.pad{calling_interval_padding}.dsp.pad{calling_interval_padding}.intersection.interval_list.ht"
        )


def check_resource_existence(
    input_step_resources: Optional[Dict[str, List]] = None,
    output_step_resources: Optional[Dict[str, List]] = None,
    overwrite: bool = False,
) -> None:
    """
    Check the existence of all specified input and output resources.

    If any of the input resources (`input_step_resources` values) don't exist, an error
    will be raised indicating which input resources are missing.

    If any of the output resources (`output_step_resources` values) already exist and
    the `overwrite` parameter is not set to True, an error will be raised indicating
    which output resources already exist.

    If no parameters are passed to the function, nothing is done.

    :param input_step_resources: A dictionary with keys as pipeline steps that generate
        input files and the value as a list of the input files to check the existence
        of. Default is None.
    :param output_step_resources: A dictionary with keys as pipeline step that generate
        output files and the value as a list of the output files to check the existence
        of. Default is None.
    :param overwrite: The overwrite parameter used when writing the output files.
        Default is False.
    :return: None.
    """
    # Check if the input resources exist
    if input_step_resources:
        for step, input_resources in input_step_resources.items():
            check_file_exists_raise_error(
                [r if isinstance(r, str) else r.path for r in input_resources],
                error_if_not_exists=True,
                error_if_not_exists_msg=(
                    f"Not all input resources exist. Please add {step} to the command "
                    "line. The following files are missing: "
                ),
            )

    # Check if the output resources exist when `overwrite` is False
    if not overwrite and output_step_resources:
        for step, output_resources in output_step_resources.items():
            check_file_exists_raise_error(
                [r if isinstance(r, str) else r.path for r in output_resources],
                error_if_exists=True,
                error_if_exists_msg=(
                    f"Some of the output resources that will be created by {step} "
                    "already exist and the --overwrite argument was not set. Please "
                    f"rerun {step} with --overwrite. The following files already exist:"
                ),
            )


class PipelineStepResourceCollection:
    def __init__(
        self,
        step_name,
        previous_pipeline_steps=[],
        output_resources=None,
        input_resources=None,
        add_input_resources=None,
        pipeline_name=None,
        overwrite=None,
    ):
        self.pipeline_name = pipeline_name
        self.overwrite = overwrite
        self.pipeline_step = step_name
        self.previous_steps = previous_pipeline_steps

        if output_resources is not None:
            for name, resource in output_resources.items():
                setattr(self, name, resource)

            output_resources = {self.pipeline_step: output_resources.values()}

        self.output_resources = output_resources

        self.input_resources = {}
        if input_resources is None:
            for step in previous_pipeline_steps:
                self.input_resources.update(step.output_resources)
        else:
            self.add_input_resources(input_resources)

        if add_input_resources is not None:
            self.add_input_resources(add_input_resources)

    def add_input_resources(self, input_resources):
        self.input_resources.update(input_resources)

    def check_resource_existance(self, overwrite=None):
        if overwrite is None:
            if self.overwrite is None:
                overwrite = False
            else:
                overwrite = self.overwrite

        check_resource_existence(
            input_step_resources=self.input_resources,
            output_step_resources=self.output_resources,
            overwrite=overwrite,
        )


class PipelineResourceCollection:
    def __init__(
        self,
        pipeline_name,
        pipeline_steps={},
        pipeline_resources=None,
        overwrite=None,
    ):
        self.pipeline_name = pipeline_name
        self.pipeline_steps = {}
        self.add_steps(pipeline_steps)
        self.overwrite = overwrite

        if pipeline_resources is not None:
            for name, resource in pipeline_resources.items():
                setattr(self, name, resource)

    def add_steps(self, steps):
        for step_name, step in steps.items():
            step.pipeline_name = self.pipeline_name
            step.overwrite = self.overwrite
            setattr(self, step_name, step)
        self.pipeline_steps.update(steps)

    def check_resource_existance(self, step, overwrite=None):
        if overwrite is None:
            if self.overwrite is None:
                overwrite = False
            else:
                overwrite = self.overwrite

        self.pipeline_steps[step].check_resource_existance(overwrite=overwrite)
