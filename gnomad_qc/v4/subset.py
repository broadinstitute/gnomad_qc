"""Script to filter the gnomAD v4 VariantDataset to a subset of specified samples."""

import argparse
import logging
from dataclasses import dataclass
from typing import List, Optional, Union

import hail as hl
from gnomad.utils.vcf import FORMAT_DICT, adjust_vcf_incompatible_types
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import (
    get_gnomad_v4_genomes_vds,
    get_gnomad_v4_vds,
    get_logging_path,
)
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)


@dataclass
class ProcessingConfig:
    """Configuration for data processing operations."""

    # Processing options
    split_multi: bool
    pass_only: bool
    add_variant_qc: bool
    variant_qc_annotations: Optional[List[str]]
    data_type: str

    # Output options
    output_vcf: bool = False
    output_vds: bool = False
    export_meta: bool = False
    keep_data_paths: bool = False
    overwrite: bool = False
    output_partitions: Optional[int] = None

    # Data loading options
    test: bool = False
    rep_on_read_partitions: Optional[int] = None

    # File paths
    output_path: str = ""
    subset_samples: Optional[str] = None
    subset_workspaces: Optional[str] = None
    tmp_dir: str = ""

    def __post_init__(self):
        if self.add_variant_qc and not self.split_multi:
            raise ValueError("Variant QC annotations require split_multi=True.")
        if self.pass_only and not self.split_multi:
            raise ValueError("PASS-only filtering requires split_multi=True.")
        if self.output_vcf and not self.split_multi:
            raise ValueError(
                "VCF export without split multi is not supported at this time."
            )

    @classmethod
    def from_args(cls, args) -> "ProcessingConfig":
        """Create ProcessingConfig from argparse arguments."""
        return cls(
            split_multi=args.split_multi,
            add_variant_qc=args.add_variant_qc,
            pass_only=args.pass_only,
            variant_qc_annotations=args.variant_qc_annotations,
            data_type=args.data_type,
            output_vcf=args.output_vcf,
            output_vds=args.output_vds,
            export_meta=args.export_meta,
            keep_data_paths=args.keep_data_paths,
            overwrite=args.overwrite,
            test=args.test,
            rep_on_read_partitions=args.rep_on_read_partitions,
            output_partitions=args.output_partitions,
            output_path=args.output_path,
            subset_samples=args.subset_samples,
            subset_workspaces=args.subset_workspaces,
            tmp_dir=args.tmp_dir,
        )


def get_gnomad_datasets(data_type: str, n_partitions: Optional[int], test: bool):
    """
    Get requested data type's v4 VariantDataset.

    :param data_type: Type of data to subset.
    :param n_partitions: Number of desired partitions for the VDS, repartioned on read.
    :param test: Whether to filter to the first 2 partitions for testing.
    :return: The gnomAD v4 VariantDataset.
    """
    if data_type == "exomes":
        vds = get_gnomad_v4_vds(
            n_partitions=n_partitions, remove_hard_filtered_samples=False
        )
    else:
        vds = get_gnomad_v4_genomes_vds(
            n_partitions=n_partitions, remove_hard_filtered_samples=False
        )

    if test:
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(range(2)),
            vds.variant_data._filter_partitions(range(2)),
        )
    meta_ht = meta(data_type=data_type).ht()

    return vds, meta_ht


def get_subset_ht(
    subset_samples: Optional[str],
    subset_workspaces: Optional[str],
    meta_ht: Optional[hl.Table],
) -> hl.Table:
    """
    Get the subset HT.

    :param subset_samples: Path to a text file with sample IDs for subsetting and a header: 's'.
    :param subset_workspaces: Path to a text file with Terra workspaces that should be included in the subset and a header: 'terra_workspace'.
    :param meta_ht: The meta HT.
    :param vds_cols: The VDS columns.
    :return: The subset HT.
    """
    if subset_workspaces and not meta_ht:
        raise ValueError("A meta HT must be provided if subsetting by workspaces.")

    if subset_workspaces:
        terra_workspaces = hl.literal(subset_ht.terra_workspace.lower().collect())
        subset_ht = meta_ht.filter(
            terra_workspaces.contains(meta_ht.project_meta.terra_workspace)
        ).select()

    if subset_samples:
        subset_ht = hl.import_table(subset_samples)
        subset_ht = subset_ht.key_by("s")

    return subset_ht


def check_subset_ht(subset_ht: hl.Table, vmt_cols: hl.Table):
    """
    Check that the subset HT is valid.

    :param subset_ht: The subset HT.
    :param vmt_cols: The variant data MatrixTable samples.
    """
    requested_sample_count = subset_ht.count()
    logger.info("Subset request has %d samples...", requested_sample_count)
    anti_join_ht = subset_ht.anti_join(vmt_cols)
    anti_join_ht_count = anti_join_ht.count()
    if anti_join_ht_count != 0:
        missing_samples = anti_join_ht.s.collect()
        message = (
            f"Only {requested_sample_count - anti_join_ht_count} out of {requested_sample_count} "
            f"subsetting-table IDs matched IDs in the dataset.\n"
            f"IDs that aren't in the callset: {missing_samples}\n"
        )
        raise ValueError(message)
    else:
        logger.info("All subsetting-table IDs are in the callset...")


def apply_split_multi_logic(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    config: ProcessingConfig,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply split multi logic to either a MatrixTable (VCF path) or VariantDataset (VDS path).

    :param mtds: MatrixTable or VariantDataset to process.
    :param config: Processing configuration.
    :return: MatrixTable or VariantDataset with multi-allelic sites split.
    """

    def _add_split_annotations(mtds):
        """Add annotations needed for split multi processing."""
        ds = mtds if config.output_vcf else mtds.variant_data
        ds = ds.annotate_rows(
            n_unsplit_alleles=hl.len(ds.alleles),
            mixed_site=(hl.len(ds.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(ds.alleles[0], a), ds.alleles[1:])
            & hl.any(lambda a: hl.is_snp(ds.alleles[0], a), ds.alleles[1:]),
        )

        return (
            ds if config.output_vcf else hl.vds.VariantDataset(mtds.reference_data, ds)
        )

    mtds = _add_split_annotations(mtds)

    if config.output_vcf:
        mtds = hl.experimental.sparse_split_multi(mtds)
        mtds = mtds.filter_rows(hl.agg.any(hl.is_defined(mtds.GT)))
        # Used during splitting multiallelics but no longer needed after.
        mtds = mtds.drop("RGQ")
    else:
        mtds = hl.vds.split_multi(
            mtds,
            filter_changed_loci=True,
        )
        mtds = hl.vds.VariantDataset(mtds.reference_data, mtds.variant_data.drop("RGQ"))

    return mtds


def apply_min_rep_logic(
    vds: hl.vds.VariantDataset,
) -> hl.vds.VariantDataset:
    """
    Apply min_rep logic to an unsplit VariantDataset.

    :param vds: VariantDataset to process.
    :param config: Processing configuration.
    :return: VariantDataset with rows keyed by minimum representation.
    """
    ds = vds.variant_data
    ds = ds.key_rows_by(**hl.min_rep(ds.locus, ds.alleles))

    return hl.vds.VariantDataset(vds.reference_data, ds)


def apply_variant_qc_annotations(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    config: ProcessingConfig,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply variant QC annotations and filtering to either a MatrixTable or VariantDataset.

    :param mtds: MatrixTable or VariantDataset to process.
    :param config: Processing configuration.
    :return: MatrixTable or VariantDataset with variant QC annotations.
    """
    ds = mtds if config.output_vcf else mtds.variant_data
    ht = release_sites(data_type=config.data_type).ht()
    selected_variant_qc_annotations = {}

    if config.variant_qc_annotations is None:
        vqc_annotations = list(ht.row_value.keys())
        vqc_annotations.remove("a_index")
        vqc_annotations.remove("was_split")
    else:
        vqc_annotations = config.variant_qc_annotations

    for ann in vqc_annotations:
        if ann not in ht.row:
            raise ValueError(f"{ann} is not a valid variant QC annotation.")
        selected_variant_qc_annotations.update([(ann, ht[ds.row_key][ann])])

    ds = ds.annotate_rows(**selected_variant_qc_annotations)

    return ds if config.output_vcf else hl.vds.VariantDataset(mtds.reference_data, ds)


def filter_to_pass_only(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    config: ProcessingConfig,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Filter to variants that passed variant QC.

    :param mtds: MatrixTable or VariantDataset to filter.
    :param config: Processing configuration.
    :return: MatrixTable or VariantDataset containing only variants that passed variant QC.
    """
    release_ht = release_sites(data_type=config.data_type).ht()
    ds = mtds if config.output_vcf else mtds.variant_data
    ds = ds.filter_rows(release_ht[ds.row_key].filters.length() == 0)

    return ds if config.output_vcf else hl.vds.VariantDataset(mtds.reference_data, ds)


def format_vcf_info_fields(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply VCF-specific formatting to info fields for export.

    :param mt: MatrixTable to format.
    :return: MatrixTable with fields formatted for VCF export.
    """
    # Note: AS_SB_TABLE should be an array of arrays of int32s.
    mt = mt.annotate_rows(
        info=mt.info.annotate(
            AS_SB_TABLE=hl.array([mt.info.AS_SB_TABLE[:2], mt.info.AS_SB_TABLE[2:]]),
        )
    )
    mt = mt.annotate_rows(info=mt.info.annotate(**mt.info.vrs).drop("vrs"))
    ht = adjust_vcf_incompatible_types(
        mt.select_rows("info").rows(), pipe_delimited_annotations=[]
    )
    mt = mt.annotate_rows(info=ht[mt.row_key].info)
    mt = mt.drop("gvcf_info")
    mt = mt.transmute_rows(rsid=hl.str(";").join(mt.rsid))

    return mt


def process_metadata_export(
    meta_ht: hl.Table,
    vmt: hl.MatrixTable,
    config: ProcessingConfig,
) -> hl.Table:
    """
    Process metadata and return the processed metadata Table.

    :param meta_ht: The metadata Table to process.
    :param vmt: The dataset's subsetted  matrixTable containing the samples to keep.
    :param config: Processing configuration.
    :return: The subsetted metadata Table.
    """
    meta_ht = meta_ht.semi_join(vmt.cols())

    if config.keep_data_paths:
        data_to_drop = {"ukb_meta"}
    else:
        data_to_drop = {"ukb_meta", "cram", "gvcf"}

    meta_ht = meta_ht.annotate(project_meta=meta_ht.project_meta.drop(*data_to_drop))

    return meta_ht


def main(args):
    """Filter the gnomAD v4 VariantDataset to a subset of specified samples."""
    logger.info(
        "Running gnomAD v4 subset script. Make sure you are using a n1-highmem-32 "
        "driver, anything less than that will hit an OOM error..."
    )
    config = ProcessingConfig.from_args(args)

    hl.init(
        log="/v4_subset.log",
        default_reference="GRCh38",
        tmp_dir=config.tmp_dir,
    )

    try:
        logger.info("Getting the v4 %s dataset...", config.data_type)
        vds, meta_ht = get_gnomad_datasets(
            data_type=config.data_type,
            n_partitions=config.rep_on_read_partitions,
            test=config.test,
        )

        logger.info("Loading the HT with sample IDs to subset...")
        subset_ht = get_subset_ht(
            config.subset_samples, config.subset_workspaces, meta_ht
        )

        logger.info("Checking that all requested sample IDs are in the dataset...")
        check_subset_ht(subset_ht, vds.variant_data.cols())

        logger.info("Filtering the VDS to requested samples...")
        vds = hl.vds.filter_samples(
            vds, subset_ht, remove_dead_alleles=False if config.split_multi else True
        )

        if config.output_vcf:
            logger.info("Densifying VDS for VCF export...")
            mtds = hl.vds.to_dense_mt(vds)
            mtds = mtds.checkpoint(new_temp_file("subset", "mt"))
        else:
            mtds = vds

        if config.split_multi:
            logger.info("Splitting multi-allelics...")
            mtds = apply_split_multi_logic(mtds, config)
        else:
            logger.info(
                "Applying min_rep to the variant data MT because remove_dead_alleles in"
                " hl.vds.filter_samples may result in variants that do not have the "
                "minimumn representation in an unsplit VDS..."
            )
            mtds = apply_min_rep_logic(mtds, config)

        if config.add_variant_qc:
            logger.info("Adding variant QC annotations...")
            mtds = apply_variant_qc_annotations(mtds, config)

        if config.pass_only:
            logger.info("Filtering to variants that passed variant QC...")
            mtds = filter_to_pass_only(mtds, config)

        if config.output_vds:
            if config.output_partitions:
                logger.info(
                    "Generating VDS with %d partitions...", config.output_partitions
                )
                temp_vds_path = new_temp_file("subset", "vds")
                mtds = mtds.write(temp_vds_path)
                mtds = hl.vds.read_vds(
                    temp_vds_path, n_partitions=config.output_partitions
                )

            logger.info("Writing VDS subset...")
            mtds.write(f"{config.output_path}/subset.vds", overwrite=config.overwrite)

        if config.output_vcf:
            logger.info("Formatting VCF info fields for export...")
            mtds = format_vcf_info_fields(mtds)
            mtds = mtds.naive_coalesce(
                config.output_partitions
                if config.output_partitions
                else mtds.n_partitions()
            )
            if config.split_multi:
                FORMAT_DICT.pop("RGQ")

            logger.info("Exporting VCF subset...")
            hl.export_vcf(
                mtds,
                f"{config.output_path}/subset.vcf.bgz",
                metadata={"format": FORMAT_DICT},
                parallel="header_per_shard",
                tabix=True,
            )

        if config.export_meta:
            logger.info("Subsetting and exporting metadata...")
            meta_ht = process_metadata_export(meta_ht, vds, config.keep_data_paths)
            meta_ht = meta_ht.checkpoint(
                f"{config.output_path}/metadata.ht", overwrite=config.overwrite
            )
            meta_ht.export(f"{config.output_path}/metadata.tsv.bgz")

        logger.info("Subset operation completed successfully!")

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("subset"))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    arg_parser = argparse.ArgumentParser(
        description=(
            "This script subsets gnomAD using a list of samples or terra workspaces."
        )
    )
    arg_parser.add_argument(
        "--test",
        help="Filter to the first 2 partitions for testing.",
        action="store_true",
    )
    subset_id_parser = arg_parser.add_mutually_exclusive_group(required=True)
    subset_id_parser.add_argument(
        "--subset-samples",
        help="Path to a text file with sample IDs for subsetting and a header: s.",
    )
    subset_id_parser.add_argument(
        "--subset-workspaces",
        help=(
            "Path to a text file with Terra workspaces that should be included in the"
            " subset, must use a header of 'terra_workspace'."
        ),
    )
    arg_parser.add_argument(
        "--data-type",
        help="Type of data to subset.",
        default="exomes",
        choices=["exomes", "genomes"],
    )

    # Make output options mutually exclusive
    output_group = arg_parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument(
        "--output-vds", help="Whether to output a subset VDS.", action="store_true"
    )
    output_group.add_argument(
        "--output-vcf", help="Whether to output a subset VCF.", action="store_true"
    )
    arg_parser.add_argument(
        "--split-multi",
        help="Whether to split multi-allelic variants.",
        action="store_true",
    )
    arg_parser.add_argument(
        "--rep-on-read-partitions",
        help=(
            "Number of partitions to pass when reading in the VDS. If passed, and the "
            "--output-partitions is not, this will be the number of output "
            "partitions. By default, there will be no change in partitioning."
        ),
        type=int,
    )
    arg_parser.add_argument(
        "--add-variant-qc",
        help=(
            "Annotate exported file with gnomAD's variant QC annotations. Defaults to"
            " all annotations if a subset of annotations are not specified using the"
            " --variant-qc-annotations arg"
        ),
        action="store_true",
    )
    arg_parser.add_argument(
        "--pass-only",
        help=(
            "Keep only the variants that passed variant QC, i.e. the filter field is"
            " PASS."
        ),
        action="store_true",
    )
    arg_parser.add_argument(
        "--variant-qc-annotations",
        nargs="+",
        type=str,
        help=(
            "Variant QC annotations to add to the output file. Defaults to all"
            " annotations."
        ),
    )
    arg_parser.add_argument(
        "--export-meta",
        help="Pull sample subset metadata and export to a HT and .tsv.",
        action="store_true",
    )
    arg_parser.add_argument(
        "--keep-data-paths",
        help="Keep CRAM and gVCF paths in the project metadata export.",
        action="store_true",
    )
    arg_parser.add_argument(
        "--output-path",
        help=(
            "Output file path for subsetted VDS/VCF/MT, do not include file extension."
        ),
        required=True,
    )
    arg_parser.add_argument(
        "--output-partitions",
        help="Number of desired partitions for the output file.",
        type=int,
    )
    arg_parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    arg_parser.add_argument(
        "--tmp-dir",
        help="Temporary directory for Hail to write files to.",
        default="gs://gnomad-tmp-4day",
    )

    return arg_parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
