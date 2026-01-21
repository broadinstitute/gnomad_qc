"""
Complete VEP annotation for gnomAD context HT.

Background:

    This script was created to recover from a VEP 115 run that failed after processing
    ~99% (37797/38029 partitions) of the gnomAD context HT (all possible SNVs). The job
    ran for an extended period before failing at task 32305 due to a VEP JSON parsing
    error. The error was caused by variant chr18:16770181 A>C, where VEP annotated the
    context field with '-nan', resulting in:

        com.fasterxml.jackson.core.JsonParseException: Unexpected character ('n' (code
        110)) in numeric value: expected digit (0-9) to follow minus sign

    After filtering that variant, another variant in the same region also failed with
    the same error. To ensure successful completion, the entire chr18 centromere region
    (chr18:15460900-20861207) is now excluded from VEP processing.

    Rather than rerun VEP on the entire context HT, this script:

        - Reconstructs the partially written HT by updating metadata files.
        - Identifies which variants still need VEP annotation.
        - Filters out the chr18 centromere region to prevent crashes.
        - Runs VEP only on the remaining unannotated variants.
        - Combines all results into a complete VEP-annotated context HT.

    Note: Variants in the chr18 centromere will have missing VEP annotations and
    should be investigated separately.

Pipeline Steps:

    Step 1: Copy partial HT to temp location.
    Step 2: Extract partition metadata from index files and vep_context HT.
    Step 3: Reconstruct partial HT by updating metadata files.
    Step 4: Filter context HT to variants missing VEP (excluding chr18 centromere).
    Step 5: Run VEP on remaining variants (excludes chr18 centromere).
    Step 6: Run VEP on chr18 centromere variants with modified config.
    Step 7: Combine all VEP results and add metadata to final HT.
"""

import argparse
import gzip
import json
import logging
import subprocess
from pathlib import Path

import hail as hl
from gnomad.resources.grch38.reference_data import vep_context
from gnomad.utils.vep import VEP_CONFIG_PATH, get_vep_help

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("complete_vep_annotation")
logger.setLevel(logging.INFO)


LOGGER_HEADER_TEMPLATE = """

===============================================================================
{message}
===============================================================================

"""

# ==============================================================================
# Path Constants
# ==============================================================================

PARTIAL_HT_PATH = "gs://gnomad-tmp-30day/persist_TableJX3FbwpfyF"

# Schema reference HT was created by running VEP 115 on a small region of the context HT
# to get the correct schema and metadata structure:
#
#   from gnomad.resources.grch38.reference_data import vep_context
#   from gnomad.utils.vep import VEP_CONFIG_PATH, get_vep_help
#
#   ht = vep_context.versions["101"].ht()
#   ht = hl.filter_intervals(
#       ht, [hl.parse_locus_interval("chr1:55039447-55064852", reference_genome="GRCh38")]
#   )
#
#   # Drop old VEP annotations
#   if "vep_proc_id" in list(ht.row):
#       ht = ht.drop("vep", "vep_proc_id")
#   else:
#       ht = ht.drop("vep")
#
#   ht = hl.vep(ht, VEP_CONFIG_PATH)
#
#   vep_help = get_vep_help(VEP_CONFIG_PATH)
#   with hl.hadoop_open(VEP_CONFIG_PATH) as vep_config_file:
#       vep_config = vep_config_file.read()
#
#   ht = ht.annotate_globals(
#       version=f"v115", vep_help=vep_help, vep_config=vep_config
#   )
#
#   ht = ht.checkpoint("gs://gnomad-julia/vep115/pcsk9.ht")
#
SCHEMA_REF_HT_PATH = "gs://gnomad-julia/vep115/pcsk9.ht"

TEMP_PARTIAL_HT_PATH = "gs://gnomad-tmp-30day/reconstructed_partial_vep.ht"
FINAL_HT_PATH = "gs://gnomad-tmp-30day/completed_vep_context.ht"
REVEP_NEEDED_HT_PATH = f"{FINAL_HT_PATH}_revep_needed"
REVEP_DONE_HT_PATH = f"{FINAL_HT_PATH}_revep_done"
CENTROMERE_REVEP_DONE_HT_PATH = f"{FINAL_HT_PATH}_centromere_revep_done"

# ==============================================================================
# VEP Version Constants
# ==============================================================================

VEP_VERSION = "v115"
VEP_CONTEXT_VERSION = "101"
VEP_PROC_ID_FIELD = "vep_proc_id"

# ==============================================================================
# Problematic Region Constants
# ==============================================================================

# Problematic region: chr18 centromere that causes VEP to return '-nan' in context
# field. Variants in this region cause JSON parsing errors and must be excluded from
# VEP processing. This region will be investigated separately.
PROBLEMATIC_REGION = "chr18:15460900-20861207"


# ==============================================================================
# Helper Functions - General Utilities
# ==============================================================================


def copy_partial_ht(partial_path: str, output_path: str) -> None:
    """
    Copy entire partial HT directory to output location.

    :param partial_path: Source path of partial HT.
    :param output_path: Destination path for copied HT.
    :return: None
    """
    logger.info("Copying entire partial HT directory...")

    # Copy everything from partial HT to output using gsutil -m cp -r.
    subprocess.run(
        ["gsutil", "-m", "cp", "-r", f"{partial_path}/*", output_path],
        check=True,
    )
    logger.info("✓ Partial HT copied")


def _read_schema_metadata(schema_ref_path: str) -> dict:
    """
    Read schema metadata from reference HT.

    :param schema_ref_path: Path to schema reference HT.
    :return: Dictionary containing the schema rows metadata.
    """
    logger.info("Reading schema reference metadata...")
    schema_rows_meta_path = f"{schema_ref_path}/rows/metadata.json.gz"
    schema_rows_meta = _read_json_metadata(schema_rows_meta_path)
    logger.info("✓ Schema metadata loaded")

    return schema_rows_meta


def _get_context_ht_path(version: str = VEP_CONTEXT_VERSION) -> str:
    """
    Get the vep_context HT path for a specific version.

    :param version: Version of the context HT (default: VEP_CONTEXT_VERSION).
    :return: Path to the context HT.
    """
    return vep_context.versions[version].path


def _read_json_metadata(file_path: str) -> dict:
    """
    Read and parse a JSON metadata file (handles both gzipped and plain JSON).

    :param file_path: Path to metadata file.
    :return: Parsed JSON as dict.
    """
    with hl.hadoop_open(file_path, "rb") as f:
        return json.loads(f.read().decode())


def _get_sorted_partition_indices(index_to_filename: dict[int, str]) -> list[int]:
    """
    Get sorted list of partition indices from index-to-filename mapping.

    :param index_to_filename: Dict mapping partition index to filename.
    :return: Sorted list of partition indices.
    """
    return sorted(index_to_filename.keys())


def _drop_vep_proc_id_if_present(ht: hl.Table) -> hl.Table:
    """
    Drop vep_proc_id field if present in table.

    :param ht: Hail Table.
    :return: Table with vep_proc_id dropped (if it was present).
    """
    if VEP_PROC_ID_FIELD in list(ht.row):
        return ht.drop(VEP_PROC_ID_FIELD)

    return ht


# ==============================================================================
# vep_context HT Operations - Bounds and Counts
# ==============================================================================


def get_context_ht_partition_counts() -> list[int]:
    """
    Get all partition counts from vep_context HT.

    :return: List of partition counts for all partitions in vep_context HT.
    """
    logger.info(f"Reading partition counts from vep_context HT metadata...")

    meta = _read_json_metadata(f"{_get_context_ht_path()}/metadata.json.gz")

    # Extract partition counts for the specific indices.
    if "components" not in meta or "partition_counts" not in meta["components"]:
        raise ValueError("vep_context HT metadata has no partition_counts field")

    all_counts = meta["components"]["partition_counts"]["counts"]
    logger.info(f"    ✓ Extracted {len(all_counts)} partition counts in vep_context HT")

    return all_counts


def get_context_ht_bounds() -> list[dict]:
    """
    Get range bounds from vep_context HT.

    :return: List of bounds where bounds[i] = bound for partition i.
    """
    logger.info(f"Reading partition bounds from vep_context HT...")

    meta = _read_json_metadata(f"{_get_context_ht_path()}/rows/metadata.json.gz")

    # Extract bounds for the specific partition indices.
    if "_jRangeBounds" not in meta:
        raise ValueError("vep_context HT metadata has no _jRangeBounds field")

    all_bounds = meta["_jRangeBounds"]
    logger.info(f"    ✓ Extracted range bounds for {len(all_bounds)} partitions")

    return all_bounds


# ==============================================================================
# Partition File Operations
# ==============================================================================


def _get_index_directories(ht_path: str) -> dict[int, str]:
    """
    Create mapping from partition index to index directory path for the given HT.

    :param ht_path: Path to HT.
    :return: Dictionary mapping partition index to directory path.
    """
    index_dir = f"{ht_path}/index"
    idx_paths = [
        f["path"]
        for f in hl.hadoop_ls(index_dir)
        if f["path"].endswith(".idx") and f["is_dir"]
    ]

    if not idx_paths:
        raise ValueError(f"No .idx directories found in {index_dir}")

    logger.info(f"    Found {len(idx_paths)} index directories to process")

    # Extract partition indices and create mapping from index to directory name.
    idx_to_filename = {}
    for idx_path in idx_paths:
        fname = Path(idx_path).name
        if fname.startswith("part-"):
            idx = int(fname.split("-")[1])
            idx_to_filename[idx] = idx_path

    logger.info(
        f"    Partition index range: {min(idx_to_filename.keys())} to "
        f"{max(idx_to_filename.keys())}"
    )

    return idx_to_filename


def _validate_partition_count(
    partition_idx: int,
    actual_count: int,
    expected_counts: list[int],
) -> None:
    """
    Validate a partition's row count against expected count from vep_context.

    :param partition_idx: Partition index.
    :param actual_count: Actual count from index metadata.
    :param expected_counts: List of expected counts from vep_context, indexed by partition index.
    :raises ValueError: If count doesn't match expected.
    """
    if partition_idx >= len(expected_counts):
        return

    expected_count = expected_counts[partition_idx]
    if actual_count != expected_count:
        raise ValueError(
            f"Partition {partition_idx} count mismatch: "
            f"expected {expected_count:,} rows from vep_context, "
            f"but found {actual_count:,} rows in index "
            f"(diff: {actual_count - expected_count:+,})"
        )


def _read_partition_counts(
    idx_dir_map: dict[int, str],
    expected_counts: list[int],
) -> dict[int, int]:
    """
    Read row counts from all partition index metadata files.

    :param idx_dir_map: Dict mapping partition index to directory path.
    :param expected_counts: List of expected counts for validation.
    :return: Dict mapping partition index to row count.
    """
    logger.info("    Reading partition counts from index metadata...")
    total_partitions = len(idx_dir_map)
    progress_interval = max(1, total_partitions // 20)  # Report ~20 times.

    # Iterate over partition indices and read counts from each .idx directory.
    partition_counts = {}
    for i, (partition_idx, idx_path) in enumerate(idx_dir_map.items(), 1):

        # Show progress periodically.
        if i % progress_interval == 0 or i == 1 or i == total_partitions:
            logger.info(
                f"      Processing partition {i}/{total_partitions} "
                f"({100*i//total_partitions}%)..."
            )

        # Read the metadata from this .idx directory.
        index_meta = _read_json_metadata(f"{idx_path}/metadata.json.gz")

        if "nKeys" not in index_meta:
            raise ValueError(
                f"No nKeys found in partition {partition_idx} ({idx_path})"
            )

        actual_count = index_meta["nKeys"]
        partition_counts[partition_idx] = actual_count

        # Validate against expected count - fail immediately on mismatch.
        _validate_partition_count(partition_idx, actual_count, expected_counts)

    logger.info("    ✓ All partition counts match vep_context HT!")

    return _convert_counts_to_list(partition_counts)


def _convert_counts_to_list(partition_counts: dict[int, int]) -> list[int]:
    """
    Convert partition counts dict to ordered list and log summary.

    :param partition_counts: Dict mapping partition index to row count.
    :return: List of row counts, ordered by partition index.
    """
    # Convert to ordered list.
    partition_indices = _get_sorted_partition_indices(partition_counts)
    counts = [partition_counts[i] for i in partition_indices]

    total_rows = sum(counts)
    logger.info(f"    ✓ Extracted counts for {len(counts)} partitions")
    logger.info(f"    Total rows across all partitions: {total_rows:,}")

    # Show sample of counts (first 10 and last 10).
    if len(counts) > 20:
        sample_counts = counts[:10] + ["..."] + counts[-10:]
        logger.info(f"    Sample counts: {sample_counts}")
    else:
        logger.info(f"    Counts: {counts}")

    return counts


# ==============================================================================
# Partition Info Persistence - Save/Load
# ==============================================================================


def _save_partition_info(
    counts: list[int],
    bounds: list[dict],
    part_file_names: list[str],
    output_path: str,
) -> None:
    """
    Save partition counts, range bounds, and partition file names to a JSON file.

    :param counts: List of partition counts.
    :param bounds: List of range bounds for each partition.
    :param part_file_names: List of partition file names.
    :param output_path: Path to save the info file.
    :return: None
    """
    info_file = f"{output_path}/partition_info.json"
    logger.info(f"Saving partition info to {info_file}")

    info_data = {
        "partition_counts": counts,
        "partition_bounds": bounds,
        "partition_file_names": part_file_names,
        "total_rows": sum(counts),
    }

    with hl.hadoop_open(info_file, "w") as f:
        f.write(json.dumps(info_data, indent=2))

    logger.info("✓ Partition info saved")


def load_partition_info(
    output_path: str,
) -> tuple[list[int], list[dict], list[str]]:
    """
    Load partition counts, range bounds, and partition file names from a JSON file.

    :param output_path: Path to load the info file from.
    :return: Tuple of (counts list, bounds list, partition file names list).
    """
    info_file = f"{output_path}/partition_info.json"

    if not hl.hadoop_exists(info_file):
        raise ValueError(f"Partition info file not found: {info_file}")

    logger.info(f"Loading partition info from {info_file}")

    with hl.hadoop_open(info_file, "r") as f:
        info_data = json.loads(f.read())

    counts = info_data["partition_counts"]
    bounds = info_data.get("partition_bounds", [])
    part_file_names = info_data["partition_file_names"]
    total_rows = info_data["total_rows"]

    logger.info(
        f"✓ Loaded partition info: {len(counts)} partitions, {total_rows:,} total rows"
    )

    return counts, bounds, part_file_names


def extract_partition_metadata_and_save(
    ht_path: str,
) -> tuple[list[int], list[dict], list[str]]:
    """
    Extract partition metadata from index files and vep_context HT, then save for later reuse.

    This function:

        - Reads partition counts from index metadata files (validated against
          vep_context).
        - Gets partition bounds from vep_context HT.
        - Derives partition file names from index directory names.
        - Saves all metadata to partition_info.json.

    :param ht_path: Path to HT.
    :return: Tuple of (counts list, bounds list, partition file names list).
    """
    # Get expected counts from vep_context for validation.
    expected_counts = get_context_ht_partition_counts()

    # Read counts from index (validates against vep_context).
    logger.info("Reading partition counts from index metadata...")
    idx_dir_map = _get_index_directories(ht_path)
    partition_counts = _read_partition_counts(idx_dir_map, expected_counts)
    partition_indices = _get_sorted_partition_indices(idx_dir_map)

    # Get bounds from vep_context HT and filter to only include partitions in the
    # partial HT.
    partition_bounds = get_context_ht_bounds()
    partition_bounds = [partition_bounds[i] for i in partition_indices]

    # Derive partition file names from index directory names.
    # Index dirs have format "part-NNNNN-UUID.idx", partition files are
    # "part-NNNNN-UUID"
    part_file_names = [
        Path(idx_dir_map[i]).name.replace(".idx", "") for i in partition_indices
    ]

    # Save for later reuse.
    _save_partition_info(partition_counts, partition_bounds, part_file_names, ht_path)

    return partition_counts, partition_bounds, part_file_names


# ==============================================================================
# Metadata Operations - Update and Copy
# ==============================================================================


def _update_and_copy_metadata_files(
    schema_ref_path: str,
    output_path: str,
    schema_rows_meta: dict,
    part_file_names: list[str],
    partition_counts: list[int] = None,
    partition_bounds: list[dict] = None,
) -> None:
    """
    Update metadata files using schema from reference HT and write to output HT.

    Reads main metadata from schema reference, updates it with partition counts and
    bounds, and writes both main and rows metadata to output HT using atomic writes.

    :param schema_ref_path: Path to schema reference HT.
    :param output_path: Path to output HT where metadata will be written.
    :param schema_rows_meta: Schema rows metadata dictionary to be updated and written.
    :param part_file_names: List of partition file names.
    :param partition_counts: Optional list of per-partition row counts to update in
        metadata.
    :param partition_bounds: Optional list of per-partition range bounds to update in
        metadata.
    :return: None
    """
    logger.info("Copying metadata files from schema reference...")

    # Update rows metadata with actual partition files.
    schema_rows_meta["_partFiles"] = part_file_names

    # Read and update main metadata with partition counts.
    schema_meta_path = f"{schema_ref_path}/metadata.json.gz"
    main_meta = _read_json_metadata(schema_meta_path)

    # Update partition counts if we have them.
    if partition_counts:
        total_rows = sum(partition_counts)
        logger.info(
            f"  Updating partition counts for {len(partition_counts)} partitions"
        )
        logger.info(f"  Total rows: {total_rows:,}")

        if "components" in main_meta and "partition_counts" in main_meta["components"]:
            main_meta["components"]["partition_counts"]["counts"] = partition_counts

    # Update range bounds if we have them.
    if (
        partition_bounds
        and isinstance(partition_bounds, list)
        and len(partition_bounds) > 0
    ):
        # Filter out None values.
        valid_bounds = [b for b in partition_bounds if b is not None]
        if valid_bounds:
            logger.info(f"  Updating range bounds for {len(valid_bounds)} partitions")
            schema_rows_meta["_jRangeBounds"] = valid_bounds
        else:
            logger.info("  No valid range bounds found, skipping bounds update")
    else:
        logger.info("  No range bounds provided, skipping bounds update")

    # Write updated main metadata to a temp file first.
    logger.info("  Writing main metadata...")
    temp_main_meta = f"{output_path}/metadata.json.gz.tmp"

    with hl.hadoop_open(temp_main_meta, "wb") as out_f:
        out_f.write(gzip.compress(json.dumps(main_meta).encode()))

    # Move temp file to final location (ensures atomic write).
    subprocess.run(
        ["gsutil", "mv", temp_main_meta, f"{output_path}/metadata.json.gz"],
        check=True,
    )
    logger.info("  ✓ Main metadata updated")

    # Write updated rows metadata to a temp file first.
    logger.info("  Writing rows metadata...")
    temp_rows_meta = f"{output_path}/rows/metadata.json.gz.tmp"

    with hl.hadoop_open(temp_rows_meta, "wb") as out_f:
        out_f.write(gzip.compress(json.dumps(schema_rows_meta).encode()))

    # Move temp file to final location (ensures atomic write).
    subprocess.run(
        ["gsutil", "mv", temp_rows_meta, f"{output_path}/rows/metadata.json.gz"],
        check=True,
    )
    logger.info("  ✓ Rows metadata updated")


def _copy_globals_directory(schema_ref_path: str, output_path: str) -> None:
    """
    Copy globals directory from schema reference to output HT.

    :param schema_ref_path: Path to schema reference HT.
    :param output_path: Path to output HT.
    :return: None
    """
    if hl.hadoop_exists(f"{schema_ref_path}/globals"):
        logger.info("  Replacing globals directory...")
        # Remove existing globals if it exists.
        if hl.hadoop_exists(f"{output_path}/globals"):
            subprocess.run(
                ["gsutil", "-m", "rm", "-r", f"{output_path}/globals"],
                check=True,
            )
        # Copy fresh globals directory.
        subprocess.run(
            ["gsutil", "-m", "cp", "-r", f"{schema_ref_path}/globals", output_path],
            check=True,
        )
        logger.info("  ✓ Globals directory replaced")


def _verify_reconstructed_table(output_path: str) -> hl.Table:
    """
    Verify that the reconstructed table can be read and return it.

    :param output_path: Path to reconstructed HT.
    :return: The reconstructed Hail Table.
    """
    logger.info("Verifying reconstructed table...")
    reconstructed_ht = hl.read_table(output_path)
    logger.info(f"✓ Successfully read reconstructed table")

    # Print schema.
    logger.info("Reconstructed table verification:")
    logger.info(f"    Partition count: {reconstructed_ht.n_partitions()}")

    reconstructed_ht.describe()
    reconstructed_ht.show(5)

    return reconstructed_ht


# ==============================================================================
# HT Reconstruction - Main Orchestration
# ==============================================================================


def reconstruct_partial_ht(
    schema_ref_path: str,
    output_path: str,
    partition_counts: list[int],
    partition_bounds: list[dict],
    part_file_names: list[str],
) -> hl.Table:
    """
    Reconstruct the partially written HT by updating metadata from schema reference.

    Assumes the partial HT has already been copied to output_path and partition info
    has been read.

    :param schema_ref_path: Path to schema reference HT.
    :param output_path: Path to copied partial HT (will update metadata in place).
    :param partition_counts: List of per-partition row counts.
    :param partition_bounds: List of per-partition range bounds.
    :param part_file_names: List of partition file names.
    :return: The reconstructed Hail Table.
    """
    # Read schema metadata from reference.
    schema_rows_meta = _read_schema_metadata(schema_ref_path)

    # Update and copy metadata files.
    _update_and_copy_metadata_files(
        schema_ref_path,
        output_path,
        schema_rows_meta,
        part_file_names,
        partition_counts,
        partition_bounds,
    )

    # Copy globals directory.
    _copy_globals_directory(schema_ref_path, output_path)

    logger.info("✓ Metadata and globals updated")

    # Verify the reconstructed table.
    return _verify_reconstructed_table(output_path)


# ==============================================================================
# VEP Annotation - Load, Prepare, Annotate, Combine
# ==============================================================================


def load_context_ht(version: str = VEP_CONTEXT_VERSION) -> hl.Table:
    """
    Load the gnomAD context HT.

    :param version: Version of the context HT to load.
    :return: Context Hail Table.
    """
    logger.info(f"Loading context HT (vep_context version {version})...")
    context_ht = vep_context.versions[version].ht()
    logger.info(f"Context HT loaded: {context_ht.count():,} variants")

    return context_ht


def load_partial_vep_ht(partial_vep_ht_path: str) -> hl.Table:
    """
    Load the partial VEP HT.

    :param partial_vep_ht_path: Path to partial VEP HT.
    :return: Partial VEP Hail Table.
    """
    logger.info(f"Loading partial VEP HT from: {partial_vep_ht_path}")
    partial_vep_ht = hl.read_table(partial_vep_ht_path)
    logger.info(f"Partial VEP HT loaded: {partial_vep_ht.count():,} variants")

    return partial_vep_ht


def prepare_context_ht(ht: hl.Table) -> hl.Table:
    """
    Prepare context HT by dropping existing VEP annotations.

    :param ht: Context Hail Table.
    :return: Prepared context Hail Table.
    """
    logger.info("Preparing context HT...")
    if "vep" in list(ht.row):
        ht = ht.drop("vep")
    ht = _drop_vep_proc_id_if_present(ht)

    return ht


def get_variants_that_need_vep(
    context_ht: hl.Table, partial_vep_ht: hl.Table
) -> hl.Table:
    """
    Get variants that need VEP.

    :param context_ht: Prepared context Hail Table.
    :param partial_vep_ht: Partial VEP Hail Table with VEP annotations on the context
        HT key.
    :return: Hail Table of variants needing VEP.
    """
    logger.info("Identifying variants needing VEP...")
    ht = context_ht.filter(hl.is_missing(partial_vep_ht[context_ht.key]))

    return ht


def filter_problematic_variants(ht: hl.Table) -> hl.Table:
    """
    Filter out variants in regions that cause VEP to fail.

    Specifically filters the chr18 centromere (chr18:15460900-20861207) which contains
    variants that cause VEP to return '-nan' in the context field, leading to JSON
    parsing errors. This region will be investigated separately.

    :param ht: Hail Table to filter.
    :return: Filtered Hail Table excluding chr18 centromere.
    """
    logger.info(f"Filtering out chr18 centromere: {PROBLEMATIC_REGION}...")
    problematic_interval = hl.parse_locus_interval(
        PROBLEMATIC_REGION, reference_genome="GRCh38"
    )

    # Count variants before filtering.
    count_before = ht.count()

    # Filter out variants in the problematic region.
    ht = ht.filter(~problematic_interval.contains(ht.locus))

    # Count variants after filtering.
    count_after = ht.count()
    filtered_count = count_before - count_after

    logger.info(f"Filtered out {filtered_count:,} variants in chr18 centromere")
    logger.info(f"Remaining variants to VEP: {count_after:,}")

    return ht


def filter_to_centromere_variants(ht: hl.Table) -> hl.Table:
    """
    Filter to ONLY variants in the chr18 centromere region.

    This function is used to isolate centromere variants for VEP processing with
    a modified configuration that excludes the context plugin (which causes '-nan'
    errors).

    :param ht: Hail Table to filter.
    :return: Filtered Hail Table containing only chr18 centromere variants.
    """
    logger.info(f"Filtering to chr18 centromere variants: {PROBLEMATIC_REGION}...")
    centromere_interval = hl.parse_locus_interval(
        PROBLEMATIC_REGION, reference_genome="GRCh38"
    )

    # Filter to only centromere variants.
    ht = ht.filter(centromere_interval.contains(ht.locus))

    count = ht.count()
    logger.info(f"Found {count:,} variants in chr18 centromere to VEP")

    return ht


def run_vep_on_remaining(
    ht: hl.Table, vep_config_path: str = VEP_CONFIG_PATH
) -> hl.Table:
    """
    Run VEP on variants that need it.

    :param ht: Hail Table of variants needing VEP.
    :param vep_config_path: Path to VEP config file.
    :return: Hail Table with VEP annotations (or original if revep_count is 0).
    """
    logger.info("Running VEP on remaining variants...")
    revep_ht = hl.vep(ht, vep_config_path)

    return revep_ht


def run_vep_on_centromere(
    ht: hl.Table, vep_config_path: str = VEP_CONFIG_PATH
) -> hl.Table:
    """
    Run VEP on chr18 centromere variants.

    WARNING: This step requires using a VEP init script that does NOT include the
    'context' plugin in the VEP command, as the context plugin causes '-nan' errors
    for variants in centromeric regions.

    :param ht: Hail Table of centromere variants needing VEP.
    :param vep_config_path: Path to VEP config file.
    :return: Hail Table with VEP annotations.
    """
    logger.warning(
        "=" * 80
        + "\n"
        + "WARNING: Running VEP on chr18 centromere variants!\n"
        + "This step REQUIRES a VEP init script WITHOUT the 'context' plugin.\n"
        + "Using the standard VEP config will cause '-nan' errors.\n"
        + "=" * 80
    )
    logger.info("Running VEP on chr18 centromere variants...")
    revep_ht = hl.vep(ht, vep_config_path)

    return revep_ht


def add_vep_metadata(ht: hl.Table, vep_config_path: str) -> hl.Table:
    """
    Add VEP metadata (version, help, config) to global annotations.

    :param ht: Hail Table with VEP annotations.
    :param vep_config_path: Path to VEP config file.
    :return: Hail Table with metadata annotations.
    """
    logger.info("Adding VEP metadata...")

    # Get VEP metadata.
    vep_help = get_vep_help(vep_config_path)
    with hl.hadoop_open(vep_config_path) as vep_config_file:
        vep_config = vep_config_file.read()

    # Annotate globals.
    ht = ht.annotate_globals(
        version=VEP_VERSION, vep_help=vep_help, vep_config=vep_config
    )

    return ht


def combine_vep_results(
    context_ht: hl.Table,
    partial_vep_ht: hl.Table,
    revep_ht: hl.Table,
    centromere_revep_ht: hl.Table,
    vep_config_path: str = VEP_CONFIG_PATH,
) -> hl.Table:
    """
    Combine VEP results from partial HT and newly VEPed variants.

    Annotates context HT with VEP from multiple sources in priority order:

        1. Partial VEP HT (from original run).
        2. ReVEP HT (newly VEPed non-centromere variants).
        3. Centromere ReVEP HT (VEPed centromere variants).

    Variants not covered by any source will have missing VEP annotations.

    :param context_ht: Context Hail Table.
    :param partial_vep_ht: Partial VEP Hail Table with VEP annotations on the context
        HT key.
    :param revep_ht: Hail Table with newly VEPed variants (excluding centromere).
    :param centromere_revep_ht: Hail Table with VEPed centromere variants.
    :param vep_config_path: Path to VEP config file.
    :return: Final combined Hail Table with VEP metadata.
    """
    logger.info("Combining VEP results...")

    # Start with partial VEP, then fill in with revep results.
    vep_expr = hl.coalesce(
        partial_vep_ht[context_ht.key].vep,
        revep_ht[context_ht.key].vep,
        centromere_revep_ht[context_ht.key].vep,
    )

    ht = context_ht.annotate(vep=vep_expr)
    ht = add_vep_metadata(ht, vep_config_path)

    return ht


# ==============================================================================
# Main Pipeline Orchestration
# ==============================================================================


def main(args) -> None:
    """Run the complete VEP 115 recovery and completion annotation pipeline."""
    # Initialize Hail with tmp directory.
    hl.init(tmp_dir="gs://gnomad-tmp-4day/hail-tmp-complete_vep_annotation/")

    # Step 1: Copy partial HT.
    if args.copy_partial_ht:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                Step 1: Copying entire partial HT directory...
                    Source: {args.partial_ht}
                    Destination: {args.temp_partial_ht}
                """
            )
        )
        copy_partial_ht(args.partial_ht, args.temp_partial_ht)
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message="✓ Step 1 complete: Partial HT copied"
            )
        )

    # Step 2: Extract metadata info from index and vep_context.
    if args.extract_metadata_info:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message="Step 2: Extracting metadata info from index and vep_context..."
            )
        )
        partition_counts, partition_bounds, part_file_names = (
            extract_partition_metadata_and_save(args.temp_partial_ht)
        )
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                ✓ Step 2 complete: Metadata extracted
                    Total variants: {sum(partition_counts)}
                    Total partitions: {len(partition_counts)}
                    Output path: {args.temp_partial_ht}
                """
            )
        )

    # Step 3: Reconstruct partial HT.
    if args.reconstruct_partial_ht:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=(
                    "Step 3: Reconstructing partial HT using schema reference:"
                    f" {args.schema_ref_ht}..."
                )
            )
        )
        partition_counts, partition_bounds, part_file_names = load_partition_info(
            args.temp_partial_ht
        )
        reconstruct_partial_ht(
            args.schema_ref_ht,
            args.temp_partial_ht,
            partition_counts,
            partition_bounds,
            part_file_names,
        )
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                ✓ Step 3 complete: Partial HT reconstructed
                    Total variants: {partition_counts}
                    Output path: {args.temp_partial_ht}
                """
            )
        )

    # Step 4: Filter context HT to variants missing VEP.
    if args.filter_context_ht_to_missing_vep:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message="Step 4: Filtering context HT to variants missing VEP..."
            )
        )
        context_ht = prepare_context_ht(load_context_ht())
        partial_vep_ht = load_partial_vep_ht(args.temp_partial_ht)
        ht = get_variants_that_need_vep(context_ht, partial_vep_ht)
        ht = filter_problematic_variants(ht)
        ht = ht.checkpoint(args.revep_needed_ht, overwrite=args.overwrite)
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                ✓ Step 4 complete: Missing VEP variants identified
                    Total variants: {ht.count():,}
                    Output path: {args.revep_needed_ht}
                """
            )
        )

    # Step 5: Run VEP on remaining variants.
    if args.run_vep_on_remaining:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message="Step 5: Running VEP on remaining variants..."
            )
        )
        ht = (
            hl.read_table(args.revep_needed_ht)
            .naive_coalesce(232)
            .checkpoint(
                "gs://gnomad-tmp-4day/revep_needed_coalesced.ht", overwrite=True
            )
        )
        ht = hl.read_table(
            "gs://gnomad-tmp-4day/revep_needed_coalesced.ht",
            _n_partitions=2000,
        )
        ht = run_vep_on_remaining(ht)
        ht = ht.checkpoint(args.revep_done_ht, overwrite=args.overwrite)
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                ✓ Step 5 complete: VEP run on remaining variants
                    Total variants: {ht.count():,}
                    Output path: {args.revep_done_ht}
                """
            )
        )

    # Step 6: Run VEP on chr18 centromere variants.
    if args.run_vep_on_centromere:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message="Step 6: Running VEP on chr18 centromere variants..."
            )
        )
        logger.warning(
            "\n" + "=" * 80 + "\n"
            "IMPORTANT: This step requires a VEP init script WITHOUT the 'context' plugin!\n"
            "Using the standard VEP configuration will cause '-nan' errors.\n"
            "Make sure your cluster was initialized with the appropriate VEP config.\n"
            + "=" * 80
            + "\n"
        )
        context_ht = prepare_context_ht(load_context_ht())
        partial_vep_ht = load_partial_vep_ht(args.temp_partial_ht)
        ht = get_variants_that_need_vep(context_ht, partial_vep_ht)
        ht = filter_to_centromere_variants(ht).cache().naive_coalesce(10)
        ht = run_vep_on_centromere(ht)
        ht = ht.checkpoint(args.centromere_revep_done_ht, overwrite=args.overwrite)
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                ✓ Step 6 complete: VEP run on chr18 centromere variants
                    Total variants: {ht.count():,}
                    Output path: {args.centromere_revep_done_ht}
                """
            )
        )

    # Step 7: Write final HT.
    if args.write_final_ht:
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(message="Step 7: Writing final HT...")
        )
        context_ht = load_context_ht()
        partial_vep_ht = load_partial_vep_ht(args.temp_partial_ht)
        revep_ht = hl.read_table(args.revep_done_ht)
        centromere_revep_ht = hl.read_table(args.centromere_revep_done_ht)

        ht = combine_vep_results(
            context_ht, partial_vep_ht, revep_ht, centromere_revep_ht
        )
        ht = ht.checkpoint(args.output, overwrite=args.overwrite)
        logger.info(
            LOGGER_HEADER_TEMPLATE.format(
                message=f"""
                ✓ Step 7 complete: Final HT written:
                    Total variants: {ht.count():,}
                    Output path: {args.output}
                """
            )
        )

    logger.info(
        LOGGER_HEADER_TEMPLATE.format(
            message=f"✓ Pipeline steps completed successfully!"
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Complete VEP annotation for gnomAD context HT."
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing files when writing checkpoints.",
    )

    # File path arguments.
    parser.add_argument(
        "--partial-ht",
        default=PARTIAL_HT_PATH,
        help=f"Path to partial HT (default: {PARTIAL_HT_PATH}).",
    )
    parser.add_argument(
        "--schema-ref-ht",
        default=SCHEMA_REF_HT_PATH,
        help=f"Path to schema reference HT (default: {SCHEMA_REF_HT_PATH}).",
    )
    parser.add_argument(
        "--temp-partial-ht",
        default=TEMP_PARTIAL_HT_PATH,
        help=f"Path for reconstructed partial HT (default: {TEMP_PARTIAL_HT_PATH}).",
    )
    parser.add_argument(
        "--output",
        default=FINAL_HT_PATH,
        help=f"Path for final output HT (default: {FINAL_HT_PATH}).",
    )
    parser.add_argument(
        "--revep-needed-ht",
        default=REVEP_NEEDED_HT_PATH,
        help=f"Path for variants needing VEP (default: {REVEP_NEEDED_HT_PATH}).",
    )
    parser.add_argument(
        "--revep-done-ht",
        default=REVEP_DONE_HT_PATH,
        help=f"Path for VEPed variants (default: {REVEP_DONE_HT_PATH}).",
    )
    parser.add_argument(
        "--centromere-revep-done-ht",
        default=CENTROMERE_REVEP_DONE_HT_PATH,
        help=f"Path for centromere VEPed variants (default: {CENTROMERE_REVEP_DONE_HT_PATH}).",
    )

    # Pipeline steps arguments.
    parser.add_argument(
        "--copy-partial-ht",
        action="store_true",
        help="Copy the partial HT to the temp reconstructed HT location.",
    )
    parser.add_argument(
        "--extract-metadata-info",
        action="store_true",
        help=(
            "Extract metadata information from the index files and original context HT "
            "and save to a JSON file."
        ),
    )
    parser.add_argument(
        "--reconstruct-partial-ht",
        action="store_true",
        help="Reconstruct the partial HT from the metadata information.",
    )
    parser.add_argument(
        "--filter-context-ht-to-missing-vep",
        action="store_true",
        help="Filter the context HT to only variants that are missing VEP annotations.",
    )
    parser.add_argument(
        "--run-vep-on-remaining",
        action="store_true",
        help="Run VEP on the remaining variants (excludes chr18 centromere).",
    )
    parser.add_argument(
        "--run-vep-on-centromere",
        action="store_true",
        help=(
            "Run VEP on chr18 centromere variants. "
            "WARNING: Requires VEP init script WITHOUT 'context' plugin."
        ),
    )
    parser.add_argument(
        "--write-final-ht",
        action="store_true",
        help="Write the final HT.",
    )

    main(parser.parse_args())
