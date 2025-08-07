"""Script to filter the gnomAD v4 VariantDataset to a subset of specified samples."""

import argparse
import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Union

import hail as hl
from gnomad.utils.vcf import FORMAT_DICT, adjust_vcf_incompatible_types
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
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
    tmp_dir: str = "gs://gnomad-tmp-4day"

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

    :param subset_samples: Path to a text file with sample IDs for subsetting and a header: s.
    :param subset_workspaces: Path to a text file with Terra workspaces that should be included in the subset, must use a header of 'terra_workspace'.
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


def check_subset_ht(subset_ht: hl.Table, vds_cols: hl.Table):
    """
    Check that the subset HT is valid.

    :param subset_ht: The subset HT.
    :param vds_cols: The VDS columns.
    """
    requested_sample_count = subset_ht.count()
    logger.info("Subset request has %d samples...", requested_sample_count)

    anti_join_ht = subset_ht.anti_join(vds_cols)
    anti_join_ht_count = anti_join_ht.count()
    if anti_join_ht_count != 0:
        missing_samples = anti_join_ht.s.collect()
        message = (
            f"Only {requested_sample_count - anti_join_ht_count} out of {requested_sample_count} "
            f"subsetting-table IDs matched IDs in the variant callset.\n"
            f"IDs that aren't in the callset: {missing_samples}\n"
        )
        raise ValueError(message)
    else:
        logger.info("All subsetting-table IDs are in the callset...")


def apply_split_multi_logic(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    is_vcf: bool = False,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply split multi logic to either a MatrixTable (VCF path) or VariantDataset (VDS path).

    :param mtds: MatrixTable or VariantDataset to process.
    :param is_vcf: Whether this is the VCF path (MatrixTable) or VDS path (VariantDataset).
    :return: Processed MatrixTable or VariantDataset.
    """
    logger.info("Splitting multi-allelics...")

    def _add_split_annotations(ds):
        """Add annotations needed for split multi processing."""
        return ds.annotate_rows(
            n_unsplit_alleles=hl.len(ds.alleles),
            mixed_site=(hl.len(ds.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(ds.alleles[0], a), ds.alleles[1:])
            & hl.any(lambda a: hl.is_snp(ds.alleles[0], a), ds.alleles[1:]),
        )

    if is_vcf:
        # VCF path: process MatrixTable directly
        mtds = _add_split_annotations(mtds)
        mtds = hl.experimental.sparse_split_multi(mtds)
        mtds = mtds.filter_rows(hl.agg.any(hl.is_defined(mtds.GT)))
        # Used during splitting multiallelics but not longer needed after.
        mtds = mtds.drop("RGQ")
    else:
        # VDS path: process VariantData component
        vd = _add_split_annotations(mtds.variant_data)
        mtds = hl.vds.split_multi(
            hl.vds.VariantDataset(mtds.reference_data, vd.drop("RGQ")),
            filter_changed_loci=True,
        )

    return mtds


def apply_min_rep_logic(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    is_vcf: bool = False,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply min_rep logic to either a MatrixTable (VCF path) or VariantDataset (VDS path).

    :param mtds: MatrixTable or VariantDataset to process.
    :param is_vcf: Whether this is the VCF path (MatrixTable) or VDS path (VariantDataset).
    :return: Processed MatrixTable or VariantDataset.
    """
    logger.info(
        "Applying min_rep to the variant data MT because remove_dead_alleles in"
        " hl.vds.filter_samples may result in variants that do not have the minimum"
        " representation..."
    )

    if is_vcf:
        mtds = mtds.key_rows_by(**hl.min_rep(mtds.locus, mtds.alleles))
    else:
        # VDS path
        vd = mtds.variant_data
        mtds = hl.vds.VariantDataset(
            mtds.reference_data,
            vd.key_rows_by(**hl.min_rep(vd.locus, vd.alleles)),
        )

    return mtds


def apply_variant_qc_annotations(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    config: ProcessingConfig,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply variant QC annotations and filtering to either a MatrixTable or VariantDataset.

    :param mtds: MatrixTable or VariantDataset to process.
    :param config: Processing configuration.
    :return: Processed MatrixTable or VariantDataset.
    """

    def _make_variant_qc_annotations_dict(
        key_expr: hl.expr.StructExpression,
        vqc_annotations: Optional[List[str]] = None,
        data_type: str = "exomes",
    ) -> Dict[str, hl.expr.Expression]:
        """
        Make a dictionary of gnomAD release annotation expressions to annotate onto the subsetted data.

        :param key_expr: Key to join annotations on.
        :param vqc_annotations: Optional list of desired annotations from the release HT.
        :param data_type: Type of data to subset. Defaults to exomes.
        :return: Dictionary containing Hail expressions to annotate onto subset.
        """
        ht = release_sites(data_type=data_type).ht()
        selected_variant_qc_annotations = {}
        if vqc_annotations is None:
            vqc_annotations = list(ht.row_value.keys())
            vqc_annotations.remove("a_index")
            vqc_annotations.remove("was_split")

        for ann in vqc_annotations:
            if ann not in ht.row:
                raise ValueError(f"{ann} is not a valid variant QC annotation.")
            selected_variant_qc_annotations.update([(ann, ht[key_expr][ann])])

        return selected_variant_qc_annotations

    ds = mtds if config.output_vcf else mtds.variant_data

    if config.add_variant_qc:
        logger.info("Adding variant QC annotations...")
        ds = ds.annotate_rows(
            **_make_variant_qc_annotations_dict(
                ds.row_key, config.variant_qc_annotations, config.data_type
            )
        )

    return ds if config.output_vcf else hl.vds.VariantDataset(mtds.reference_data, ds)


def filter_to_pass_only(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    data_type: str,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Filter to variants that passed variant QC.

    :param mtds: MatrixTable or VariantDataset to filter.
    :return: Filtered MatrixTable or VariantDataset.
    """
    logger.info("Filtering to variants that passed variant QC...")
    release_ht = release_sites(data_type=data_type).ht()
    if isinstance(mtds, hl.MatrixTable):
        mtds = mtds.filter_rows(release_ht[mtds.row_key].filters.length() == 0)
    else:
        mtds.variant_data = mtds.variant_data.filter_rows(
            release_ht[mtds.variant_data.row_key].filters.length() == 0
        )
        mtds = hl.vds.VariantDataset(mtds.reference_data, mtds.variant_data)

    return mtds


def format_vcf_info_fields(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Apply VCF-specific formatting to info fields for export.

    :param mt: MatrixTable to format.
    :return: Formatted MatrixTable.
    """
    logger.info("Formatting VCF info fields...")
    # Note: AS_SB_TABLE should be an array of arrays of int32s. However, when we split
    # the info fields for v3 and v4, we use the extend method which adds the second
    # array to the first https://github.com/broadinstitute/gnomad_methods/blob/95b42cb013da1e63311245426ab1574208405f56/gnomad/utils/sparse_mt.py#L775C5-L778C10
    # This means that the AS_SB_TABLE will be an array of four int32s, not an array of
    # two arrays of two int32s.
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


def process_vds_output(
    vds: hl.vds.VariantDataset,
    config: ProcessingConfig,
) -> hl.vds.VariantDataset:
    """
    Process VDS output and return the processed VDS.

    :param vds: The VDS to process.
    :param config: Processing configuration.
    :return: The processed VDS.
    """
    if config.split_multi:
        vds = apply_split_multi_logic(vds)

    vds = apply_min_rep_logic(vds)

    if config.add_variant_qc:
        vds = apply_variant_qc_annotations(vds, config)

    if config.pass_only:
        vds = filter_to_pass_only(vds, config.data_type)

    return vds


def process_vcf_output(
    vds: hl.vds.VariantDataset,
    config: ProcessingConfig,
    output_partitions: Optional[int],
) -> hl.MatrixTable:
    """
    Process VCF output and return the processed MatrixTable.

    :param vds: The VDS to process.
    :param config: Processing configuration.
    :param output_partitions: Number of desired partitions for the output file.
    :return: The processed MatrixTable.
    """
    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.checkpoint(new_temp_file("subset", "mt"))

    # mt = hl.read_matrix_table("gs://gnomad-tmp-4day/subset-2j5jnBbOKUoH3EwQSsjuid.mt")

    if config.split_multi:
        mt = apply_split_multi_logic(mt, is_vcf=True)

    mt = apply_min_rep_logic(mt, is_vcf=True)

    if config.add_variant_qc:
        mt = apply_variant_qc_annotations(mt, config)
    if config.pass_only:
        mt = filter_to_pass_only(mt, config.data_type)

    mt = format_vcf_info_fields(mt)

    # Set final partitioning
    mt = mt.naive_coalesce(
        output_partitions if output_partitions else mt.n_partitions()
    )

    return mt


def process_metadata_export(
    meta_ht: hl.Table,
    vds: hl.vds.VariantDataset,
    keep_data_paths: bool,
) -> hl.Table:
    """
    Process metadata and return the processed metadata Table.

    :param meta_ht: The metadata Table to process.
    :param vds: The VDS to process.
    :param keep_data_paths: Whether to keep data paths in the metadata.
    :return: The processed metadata Table.
    """
    logger.info("Exporting metadata...")
    meta_ht = meta_ht.semi_join(vds.variant_data.cols())

    if keep_data_paths:
        data_to_drop = {"ukb_meta"}
    else:
        data_to_drop = {"ukb_meta", "cram", "gvcf"}

    meta_ht = meta_ht.annotate(project_meta=meta_ht.project_meta.drop(*data_to_drop))
    return meta_ht


def main(args):
    """Filter the gnomAD v4 VariantDataset to a subset of specified samples."""
    # NOTE: genomes VDS needs a n1-highmem-32 driver to run. A 16 will hit a
    # out of heap memory error.

    # Create processing configuration from arguments
    config = ProcessingConfig.from_args(args)

    hl.init(
        log="/v4_subset.log",
        default_reference="GRCh38",
        tmp_dir=config.tmp_dir,
        # Proven Spark optimizations for large genomic datasets
        spark_conf={
            # Adaptive Query Execution - widely recommended for large datasets
            "spark.sql.adaptive.enabled": "true",
            "spark.sql.adaptive.coalescePartitions.enabled": "true",
            "spark.sql.adaptive.skewJoin.enabled": "true",
            # Kryo serialization - faster for large genomic objects
            "spark.serializer": "org.apache.spark.serializer.KryoSerializer",
            "spark.kryo.registrationRequired": "false",
            # Conservative partition management
            "spark.sql.shuffle.partitions": "2000",
        },
    )

    logger.info("Getting gnomAD %s dataset...", config.data_type)
    vds, meta_ht = get_gnomad_datasets(
        data_type=config.data_type,
        n_partitions=config.rep_on_read_partitions,
        test=config.test,
    )

    logger.info("Getting subset HT...")
    subset_ht = get_subset_ht(config.subset_samples, config.subset_workspaces, meta_ht)

    logger.info("Checking that all subsetting-table IDs are in the callset...")
    check_subset_ht(subset_ht, vds.variant_data.cols())

    logger.info("Filtering VDS to %d subset samples...", subset_ht.count())
    vds = hl.vds.filter_samples(
        vds, subset_ht, remove_dead_alleles=False if config.split_multi else True
    )

    # After discussions with the hail team, a VCF export should densify the VDS and
    # checkpoint it before splitting multiallelics for efficiency. This means few
    # functions can be reused downstream as dense MTs and VDSs need to be processed
    # differently. Because of this, we have separate mini-pipelines for VDS and VCF
    # exports.
    # TODO Should we drop that horrific multi-allelic site?
    if config.output_vds:
        logger.info("Producing VDS subset...")
        vds = process_vds_output(vds, config)
        vds.write(f"{config.output_path}/subset.vds", overwrite=config.overwrite)

    if config.output_vcf:
        logger.info("Producing VCF subset...")
        mt = process_vcf_output(vds, config, config.output_partitions)

        if config.split_multi:
            FORMAT_DICT.pop("RGQ")
        hl.export_vcf(
            mt,
            f"{config.output_path}/subset.vcf.bgz",
            metadata={"format": FORMAT_DICT},
            parallel="header_per_shard",
            tabix=True,
        )

    if config.export_meta:
        meta_ht = process_metadata_export(meta_ht, vds, config.keep_data_paths)
        meta_ht = meta_ht.checkpoint(
            f"{config.output_path}/metadata.ht", overwrite=config.overwrite
        )
        meta_ht.export(f"{config.output_path}/metadata.tsv.bgz")


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
