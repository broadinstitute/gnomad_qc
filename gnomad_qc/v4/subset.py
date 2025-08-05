"""Script to filter the gnomAD v4 VariantDataset to a subset of specified samples."""

import argparse
import logging
from typing import Dict, List, Optional, Union

import hail as hl
from gnomad.utils.vcf import FORMAT_DICT
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)


def make_variant_qc_annotations_dict(
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
    vcf_export: bool = False,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply split multi logic to either a MatrixTable (VCF path) or VariantDataset (VDS path).

    :param mtds: MatrixTable or VariantDataset to process.
    :param is_vcf: Whether this is the VCF path (MatrixTable) or VDS path (VariantDataset).
    :return: Processed MatrixTable or VariantDataset.
    """
    logger.info("Splitting multi-allelics...")
    if vcf_export:
        mtds = mtds.annotate_rows(
            n_unsplit_alleles=hl.len(mtds.alleles),
            mixed_site=(hl.len(mtds.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(mtds.alleles[0], a), mtds.alleles[1:])
            & hl.any(lambda a: hl.is_snp(mtds.alleles[0], a), mtds.alleles[1:]),
        )
        # Need to use sparse split multi because it is local-aware.
        mtds = hl.experimental.sparse_split_multi(mtds)
        mtds = mtds.filter_rows(hl.agg.any(hl.is_defined(mtds.GT)))
    else:
        # VDS path
        vd = mtds.variant_data
        vd = vd.annotate_rows(
            n_unsplit_alleles=hl.len(vd.alleles),
            mixed_site=(hl.len(vd.alleles) > 2)
            & hl.any(lambda a: hl.is_indel(vd.alleles[0], a), vd.alleles[1:])
            & hl.any(lambda a: hl.is_snp(vd.alleles[0], a), vd.alleles[1:]),
        )
        mtds = hl.vds.split_multi(
            hl.vds.VariantDataset(mtds.reference_data, vd), filter_changed_loci=True
        )

    return mtds


def apply_min_rep_logic(
    mtds: Union[hl.MatrixTable, hl.vds.VariantDataset],
    vcf_export: bool = False,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply min_rep logic to either a MatrixTable (VCF path) or VariantDataset (VDS path).

    :param mtds: MatrixTable or VariantDataset to process.
    :param vcf_export: Whether this is the VCF path (MatrixTable) or VDS path (VariantDataset).
    :return: Processed MatrixTable or VariantDataset.
    """
    logger.info(
        "Applying min_rep to the variant data MT because remove_dead_alleles in"
        " hl.vds.filter_samples may result in variants that do not have the minimum"
        " representation..."
    )

    if vcf_export:
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
    add_variant_qc: bool,
    pass_only: bool,
    variant_qc_annotations: Optional[List[str]],
    data_type: str,
    vcf_export: bool = False,
) -> Union[hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Apply variant QC annotations and filtering to either a MatrixTable or VariantDataset.

    :param mtds: MatrixTable or VariantDataset to process.
    :param add_variant_qc: Whether to add variant QC annotations.
    :param pass_only: Whether to filter to PASS variants only.
    :param variant_qc_annotations: List of variant QC annotations to add.
    :param data_type: Type of data (exomes or genomes).
    :param vcf_export: Whether this is the VCF path (MatrixTable) or VDS path (VariantDataset).
    :return: Processed MatrixTable or VariantDataset.
    """
    if not (add_variant_qc or pass_only):
        return mtds

    # Get the appropriate data structure for processing
    if vcf_export:
        ds = mtds  # MatrixTable
    else:
        ds = mtds.variant_data  # VariantData from VDS

    ht = release_sites(data_type=data_type).ht()

    if pass_only:
        logger.info("Filtering to variants that passed variant QC...")
        ds = ds.filter_rows(ht[ds.row_key].filters.length() == 0)

    if add_variant_qc:
        logger.info("Adding variant QC annotations...")
        ds = ds.annotate_rows(
            **make_variant_qc_annotations_dict(
                ds.row_key, variant_qc_annotations, data_type
            )
        )
    # Return the appropriate structure
    if vcf_export:
        return ds  # MatrixTable
    else:
        return hl.vds.VariantDataset(mtds.reference_data, ds)  # VDS


def process_vds_output(
    vds: hl.vds.VariantDataset,
    split_multi: bool,
    add_variant_qc: bool,
    pass_only: bool,
    variant_qc_annotations: Optional[List[str]],
    data_type: str,
) -> hl.vds.VariantDataset:
    """
    Process VDS output and return the processed VDS.

    :param vds: The VDS to process.
    :param split_multi: Whether to split multi-allelic variants.
    :param add_variant_qc: Whether to add variant QC annotations.
    :param pass_only: Whether to filter to PASS variants only.
    :param variant_qc_annotations: List of variant QC annotations to add.
    :param data_type: Type of data (exomes or genomes).
    :return: The processed VDS.
    """
    vds = apply_split_multi_logic(vds, vcf_export=False)
    vds = apply_min_rep_logic(vds, vcf_export=False)
    vds = apply_variant_qc_annotations(
        vds,
        add_variant_qc,
        pass_only,
        variant_qc_annotations,
        data_type,
        vcf_export=False,
    )
    return vds


def process_vcf_output(
    vds: hl.vds.VariantDataset,
    split_multi: bool,
    add_variant_qc: bool,
    pass_only: bool,
    variant_qc_annotations: Optional[List[str]],
    data_type: str,
    output_partitions: Optional[int],
) -> hl.MatrixTable:
    """
    Process VCF output and return the processed MatrixTable.

    :param vds: The VDS to process.
    :param split_multi: Whether to split multi-allelic variants.
    :param add_variant_qc: Whether to add variant QC annotations.
    :param pass_only: Whether to filter to PASS variants only.
    :param variant_qc_annotations: List of variant QC annotations to add.
    :param data_type: Type of data (exomes or genomes).
    :param output_partitions: Number of desired partitions for the output file.
    :return: The processed MatrixTable.
    """
    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.checkpoint(new_temp_file("subset", "mt"))

    if split_multi:
        mt = apply_split_multi_logic(mt, vcf_export=True)

    mt = apply_min_rep_logic(mt, vcf_export=True)
    mt = apply_variant_qc_annotations(
        mt,
        add_variant_qc,
        pass_only,
        variant_qc_annotations,
        data_type,
        vcf_export=True,
    )
    mt = mt.drop("gvcf_info")
    mt = mt.transmute_rows(rsid=hl.str(";").join(mt.rsid))
    mt = mt.annotate_rows(info=mt.info.annotate(**mt.info.vrs).drop("vrs"))
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
    hl.init(log="/v4_subset.log", default_reference="GRCh38", tmp_dir=args.tmp_dir)
    test = args.test
    data_type = args.data_type
    output_path = args.output_path
    add_variant_qc = args.add_variant_qc
    variant_qc_annotations = args.variant_qc_annotations
    pass_only = args.pass_only
    split_multi = args.split_multi
    rep_on_read_partitions = args.rep_on_read_partitions
    output_partitions = args.output_partitions
    vcf = args.vcf
    vds = args.vds
    overwrite = args.overwrite

    if vcf and not split_multi:
        raise ValueError(
            "VCF export without split multi is not supported at this time."
        )
    if (add_variant_qc or pass_only) and not split_multi:
        raise ValueError(
            "Cannot annotate variant QC annotations or filter to PASS variants on an"
            " unsplit dataset."
        )

    if not vcf and not vds:
        raise ValueError("One of --vcf or --vds must be specified.")

    logger.info("Getting gnomAD %s dataset...", data_type)
    vds, meta_ht = get_gnomad_datasets(
        data_type=data_type, n_partitions=rep_on_read_partitions, test=test
    )

    logger.info("Getting subset HT...")
    subset_ht = get_subset_ht(args.subset_samples, args.subset_workspaces, meta_ht)

    logger.info("Checking that all subsetting-table IDs are in the callset...")
    check_subset_ht(subset_ht, vds.variant_data.cols())

    logger.info("Filtering VDS to subset samples...")
    vds = hl.vds.filter_samples(
        vds, subset_ht, remove_dead_alleles=False if split_multi else True
    )

    # After discussions with the hail team, a VCF export should densify the VDS and
    # checkpoint it before splitting multiallelics for efficiency. This means few
    # functions can be reused downstream as dense MTs and VDSs need to be processed
    # differently. Because of this, we have separate mini-pipelines for VDS and VCF
    # exports.
    #  # TODO: Do we want to allow only VDS or VCF? Not many steps are shared at the moment so there is not much to gain from a single pipeline for both.
    # TODO Should we drop that horrific multi-allelic site?
    if vds:
        vds = process_vds_output(
            vds,
            split_multi,
            add_variant_qc,
            pass_only,
            variant_qc_annotations,
            data_type,
        )
        vds.write(f"{output_path}/subset.vds", overwrite=overwrite)

    if vcf:
        mt = process_vcf_output(
            vds,
            split_multi,
            add_variant_qc,
            pass_only,
            variant_qc_annotations,
            data_type,
            output_partitions,
        )
        hl.export_vcf(
            mt,
            f"{output_path}/subset.vcf.bgz",
            metadata=FORMAT_DICT,
            parallel="header_per_shard",
            tabix=True,
        )

    if args.export_meta:
        meta_ht = process_metadata_export(
            meta_ht,
            vds,
            args.keep_data_paths,
        )
        meta_ht = meta_ht.checkpoint(f"{output_path}/metadata.ht", overwrite=overwrite)
        meta_ht.export(f"{output_path}/metadata.tsv.bgz")


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "This script subsets gnomAD using a list of samples or terra workspaces."
        )
    )
    parser.add_argument(
        "--test",
        help="Filter to the first 2 partitions for testing.",
        action="store_true",
    )
    subset_id_parser = parser.add_mutually_exclusive_group(required=True)
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
    parser.add_argument(
        "--data-type",
        help="Type of data to subset.",
        default="exomes",
        choices=["exomes", "genomes"],
    )
    parser.add_argument(
        "--vds", help="Whether to make a subset VDS.", action="store_true"
    )
    parser.add_argument(
        "--vcf", help="Whether to make a subset VCF.", action="store_true"
    )
    parser.add_argument(
        "--split-multi",
        help="Whether to split multi-allelic variants.",
        action="store_true",
    )
    parser.add_argument(
        "--rep-on-read-partitions",
        help=(
            "Number of partitions to pass when reading in the VDS. If passed, and the "
            "--output-partitions is not, this will be the number of output "
            "partitions. By default, there will be no change in partitioning."
        ),
        type=int,
    )
    parser.add_argument(
        "--add-variant-qc",
        help=(
            "Annotate exported file with gnomAD's variant QC annotations. Defaults to"
            " all annotations if a subset of annotations are not specified using the"
            " --variant-qc-annotations arg"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--pass-only",
        help=(
            "Keep only the variants that passed variant QC, i.e. the filter field is"
            " PASS."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--variant-qc-annotations",
        nargs="+",
        type=str,
        help=(
            "Variant QC annotations to add to the output file. Defaults to all"
            " annotations."
        ),
    )
    parser.add_argument(
        "--export-meta",
        help="Pull sample subset metadata and export to a HT and .tsv.",
        action="store_true",
    )
    parser.add_argument(
        "--keep-data-paths",
        help="Keep CRAM and gVCF paths in the project metadata export.",
        action="store_true",
    )
    parser.add_argument(
        "--output-path",
        help=(
            "Output file path for subsetted VDS/VCF/MT, do not include file extension."
        ),
        required=True,
    )
    parser.add_argument(
        "--output-partitions",
        help="Number of desired partitions for the output file.",
        type=int,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--tmp-dir",
        help="Temporary directory for Hail to write files to.",
        default="gs://gnomad-tmp-4day",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
