import logging

import hail as hl
from gnomad.resources.resource_utils import (
    TableResource,
    VariantDatasetResource,
    VersionedTableResource,
    VersionedVariantDatasetResource,
)

from gnomad_qc.v4.resources.constants import CURRENT_VERSION
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.sample_qc import hard_filtered_samples

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)


# Note: Unlike previous versions, the v4 resource directory uses a general format of hgs://gnomad/v4.0/<module>/<exomes_or_genomes>/
def get_gnomad_v4_vds(
    split=False,
    remove_hard_filtered_samples: bool = True,
    release_only: bool = False,
) -> hl.vds.VariantDataset:
    """
    Wrapper function to get gnomAD v4 data with desired filtering and metadata annotations.

    :param split: Perform split on MT - Note: this will perform a split on the MT rather than grab an already split MT
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard filters (only relevant after sample QC)
    :param release_only: Whether to filter the MT to only samples available for release (can only be used if metadata is present)
    :return: gnomAD v4 dataset with chosen annotations and filters
    """
    vds = gnomad_v4_genotypes.vds()
    if remove_hard_filtered_samples:
        vds = hl.vds.filter_samples(
            vds, hard_filtered_samples.versions[CURRENT_VERSION].ht(), keep=False
        )

    if release_only:
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

    # Count current number of samples in the VDS
    n_samples = vds.variant_data.count_cols()

    # Remove 76 withdrawn UKB samples (samples with withdrawn consents for application 26041 on 08/09/2021 and application 31063 on 02/22/2022)
    ukb_application_map_ht = ukb_application_map.ht()
    withdrawn_ukb_samples = ukb_application_map_ht.filter(
        ukb_application_map_ht.withdraw
    ).s.collect()

    # Remove 43 samples that are known to be on the pharma's sample remove list
    # See https://github.com/broadinstitute/ukbb_qc/blob/70c4268ab32e9efa948fe72f3887e1b81d8acb46/ukbb_qc/resources/basics.py#L308
    dups_ht = ukb_known_dups.ht()
    ids_to_remove = dups_ht.aggregate(hl.agg.collect(dups_ht["Sample Name - ID1"]))

    # Remove 27 fully duplicated IDs (same exact name for 's' in the VDS)
    # See https://github.com/broadinstitute/ukbb_qc/blob/70c4268ab32e9efa948fe72f3887e1b81d8acb46/ukbb_qc/resources/basics.py#L286
    # Confirmed that the col_idx of the 27 dup samples in the original UKB MT match the col_idx of the dup UKB samples in the VDS
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

    # Filter withdrawn samples from the VDS
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

    # Log number of UKB samples removed from the VDS
    n_samples_after_exclusion = vds.variant_data.count_cols()
    n_samples_removed = n_samples - n_samples_after_exclusion

    logger.info(
        "Total number of UKB samples removed from the VDS: %d", n_samples_removed
    )

    return vds


_gnomad_v4_genotypes = {
    "4.0": VariantDatasetResource("gs://gnomad/v4.0/raw/exomes/gnomad_v4.0.vds"),
}

gnomad_v4_genotypes = VersionedVariantDatasetResource(
    CURRENT_VERSION,
    _gnomad_v4_genotypes,
)

# v4 test dataset VDS
gnomad_v4_testset = VariantDatasetResource(
    "gs://gnomad/v4.0/raw/exomes/testing/gnomad_v4.0_test.vds"
)
gnomad_v4_testset_meta = TableResource(
    "gs://gnomad/v4.0/raw/exomes/testing/gnomad_v4.0_meta.ht"
)


# UKB data resources
def _ukb_root_path() -> str:
    """
    Retrieve the path to the UKB data directory.

    :return: String representation of the path to the UKB data directory
    """
    return "gs://gnomad/v4.0/ukbb"


# List of samples to exclude from QC due to withdrawn consents for application 26041 on 08/09/2021 and application 31063 on 02/22/2022, originally imported at csvs and keyed by their respective eids.
ukb_excluded = VersionedTableResource(
    default_version="31063_20220222",
    versions={
        "26041_20210809": TableResource(path=f"{_ukb_root_path()}/w26041_20210809.ht"),
        "31063_20220222": TableResource(
            path=f"{_ukb_root_path()}/ukb31063.withdrawn_participants.20220222.ht"
        ),
    },
)

# UKB map of exome IDs to array sample IDs (application ID: 26041)
ukb_array_sample_map = TableResource(f"{_ukb_root_path()}/array_sample_map_freeze_7.ht")

# UKB full mapping file of sample ID and application IDs 26041, 31063, and 48511
ukb_application_map = TableResource(
    f"{_ukb_root_path()}/ukbb_application_id_mappings.ht"
)


def generate_ukb_application_map() -> None:
    """
    Generate the mapping file of application IDs 26041, 31063, and 48511, and add in annotatation 'withdraw', which is set to True if the UKB sample needs to be withdrawn.

    :return: None
    """
    # The array sample map has the most complete list of 's' and 'eid_26041' for UKB samples, so use this file rather than the UKB 'meta.ht' (which already has some samples withdrawn from it)
    ukb_array_sample_map = ukb_array_sample_map.ht()
    ukb_array_sample_map = ukb_array_sample_map.select(
        eid_26041=ukb_array_sample_map.ukbb_app_26041_id
    )

    # Add in bridge mapping from eid_26041 to eid_31063
    ukb_map = ukb_array_sample_map.key_by("eid_26041")
    mapping_ht_26041_31063 = hl.import_table(
        "gs://gnomad/v4.0/ukbb/bridge_26041_31063.csv", delimiter=","
    ).key_by("eid_26041")
    ukb_map = ukb_map.join(mapping_ht_26041_31063, how="outer")

    # Add in bridge mapping from eid_31063 to eid_48511
    ukb_map = ukb_map.key_by("eid_31063")
    mapping_ht_48511_31063 = hl.import_table(
        "gs://gnomad/v4.0/ukbb/bridge_48511_31063", delimiter=","
    ).key_by("eid_31063")
    ukb_map = ukb_map.join(mapping_ht_48511_31063, how="outer")

    # Annotate UKB samples that need to be withdrawn (these samples will be removed when loading the VDS with the 'get_gnomad_v4_vds' function)
    withdrawn_31063_samples = hl.literal(
        ukb_excluded.versions["31063_20220222"].ht().eid_31063.collect()
    )
    withdrawn_26041_samples = hl.literal(
        ukb_excluded.versions["26041_20210809"].ht().eid_26041.collect()
    )
    ukb_map = ukb_map.annotate(
        withdraw=(withdrawn_31063_samples.contains(ukb_map.eid_31063))
        | (withdrawn_26041_samples.contains(ukb_map.eid_26041))
    )
    ukb_map = ukb_map.key_by("s")

    total_samples_in_mapping = ukb_map.count()
    total_eid_26041 = ukb_map.aggregate(
        hl.agg.count_where(~hl.is_missing(ukb_map.eid_26041))
    )
    total_eid_31063 = ukb_map.aggregate(
        hl.agg.count_where(~hl.is_missing(ukb_map.eid_31063))
    )
    total_eid_48511 = ukb_map.aggregate(
        hl.agg.count_where(~hl.is_missing(ukb_map.eid_48511))
    )

    logger.info(
        """Found %d UKB samples in mapping file\n
            %d have an eid_26041\n 
            %d have an eid_31063\n
            %d have an eid_48511\n""",
        total_samples_in_mapping,
        total_eid_26041,
        total_eid_31063,
        total_eid_48511,
    )

    ukb_map.write(ukb_application_map.path, overwrite=True)


# Samples known to be on the pharma partners' remove lists.
# All 44 samples are marked as "unresolved duplicates" by the pharma partners.
ukb_known_dups = TableResource(f"{_ukb_root_path()}/pharma_known_dups_7.ht")

# Samples with duplicate names in the VDS and their column index
# 27 samples to remove based on column index
ukb_dups_idx_path = f"{_ukb_root_path()}/dup_remove_idx_7.tsv"

# Final list of UKB samples to remove (duplicates that were removed will have their original index appended to the sample name)
all_ukb_samples_to_remove = f"{_ukb_root_path()}/all_ukbb_samples_to_remove.txt"


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
    mt = mt.annotate_cols(meta_name=meta.versions[version].ht()[mt.col_key])

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
