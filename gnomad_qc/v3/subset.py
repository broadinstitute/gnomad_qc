"""Script to filter the gnomAD v3 VDS to a subset of specified samples with optional sample QC and variant QC annotations."""

import argparse
import logging
from typing import Dict, List

import hail as hl
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.vcf import adjust_vcf_incompatible_types

from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds, get_logging_path
from gnomad_qc.v3.resources.meta import meta as metadata
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples
from gnomad_qc.v4.subset import HEADER_DICT, SUBSET_CALLSTATS_INFO_DICT

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

INFO_VCF_AS_PIPE_DELIMITED_FIELDS = [
    "AS_QUALapprox",
    "AS_VarDP",
    "AS_MQ_DP",
    "AS_RAW_MQ",
    "AS_SB_TABLE",
]

VARIANT_QC_ANNOTATIONS = [
    "freq",
    "raw_qual_hists",
    "popmax",
    "qual_hists",
    "faf",
    "rsid",
    "filters",
    "info",
    "vep",
    "vqsr",
    "region_flag",
    "allele_info",
    "age_hist_het",
    "age_hist_hom",
    "cadd",
    "revel",
    "splice_ai",
    "primate_ai",
]


def adjust_vcf_incompatible_types(
    ht: hl.Table,
    pipe_delimited_annotations: List[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,
) -> hl.Table:
    """
    Create a Table ready for vcf export.

    In particular, the following conversions are done:
        - All int64 are coerced to int32
        - Fields specified by `pipe_delimited_annotations` are converted from arrays to pipe-delimited strings

    :param ht: Input Table.
    :param pipe_delimited_annotations: List of info fields (they must be fields of the ht.info Struct).
    :return: Table ready for VCF export.
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), "")), "|")

    # Make sure the HT is keyed by locus, alleles
    ht = ht.key_by("locus", "alleles")

    info_type_convert_expr = {}
    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            logger.warning(
                "Coercing field info.%s from int64 to int32 for VCF output. Value will"
                " be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update(
                {f: hl.int32(hl.min(2**31 - 1, ht.info[f]))}
            )
        elif ft == hl.dtype("array<int64>"):
            logger.warning(
                "Coercing field info.%s from array<int64> to array<int32> for VCF"
                " output. Array values will be capped at int32 max value.",
                f,
            )
            info_type_convert_expr.update(
                {f: ht.info[f].map(lambda x: hl.int32(hl.min(2**31 - 1, x)))}
            )

    ht = ht.annotate(info=ht.info.annotate(**info_type_convert_expr))

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, it
    # is handled below
    for f in pipe_delimited_annotations:
        if f in ht.info and f != "AS_SB_TABLE":
            info_expr[f] = "|" + get_pipe_expr(ht.info[f])

    # Flatten SB if it is an array of arrays
    if "SB" in ht.info and not isinstance(ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr["SB"] = ht.info.SB[0].extend(ht.info.SB[1])

    if "AS_SB_TABLE" in ht.info:
        info_expr["AS_SB_TABLE"] = get_pipe_expr(ht.info.AS_SB_TABLE)

    # Annotate with new expression
    ht = ht.annotate(info=ht.info.annotate(**info_expr))

    return ht


def compute_partitions(
    mt, entry_size=3.5, partition_size=128000000, min_partitions=20
) -> int:  # TODO: move to gnomad_methods?
    """
    Compute a very rough estimate for the optimal number of partitions in a MT.

    This uses the hail recommended partition size of 128MB and the rough estimate
    of 3.5B per entry in the gnomAD sparse MT.

    :param mt: MatrixTable that will be resized
    :param entry_size: Average size in bytes that a single entry requires , defaults to 3.5
    :param partition_size: Amount of data per partition. Hail recommends 128MB, defaults to 128000000
    :param min_partitions: Minimum number of partitions for naive_coalesce, defaults to 20
    :return: number of partitions for resizing MT
    """
    rows, columns = mt.count()
    mt_disk_est = rows * columns * entry_size
    n_partitions = hl.eval(hl.max(int(mt_disk_est / partition_size), min_partitions))
    return n_partitions


def make_variant_qc_annotations_dict(
    key_expr: hl.expr.StructExpression,
    vqc_annotations: List[str] = VARIANT_QC_ANNOTATIONS,
) -> Dict[str, hl.expr.Expression]:
    """
    Make a dictionary of gnomAD release annotation expressions to annotate onto the subsetted data.

    :param release_ht: The hail table containing desired annotations.
    :param key_expr: Key to join annotations on.
    :param vqc_annotations: List of desired annotations from the release HT.
    :return: Dictionary containing Hail experssions to annotate onto subset.
    """
    ht = release_sites().ht()
    selected_variant_qc_annotations = {}
    for ann in vqc_annotations:
        if ann in VARIANT_QC_ANNOTATIONS:
            selected_variant_qc_annotations.update([(ann, ht[key_expr][ann])])

    return selected_variant_qc_annotations


def main(args):
    """Filter the gnomAD v3 VDS to a subset of specified samples."""
    hl.init(
        log="/subset.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp/subset_v3",
    )
    hl._set_flags(no_whole_stage_codegen="1")
    hard_filter = args.hard_filter_samples
    release_only = args.release_only
    output_path = args.output_path
    vds_output = args.vds
    dense_mt = args.dense_mt
    split_multi = args.split_multi
    subset_call_stats = args.subset_call_stats
    add_variant_qc = args.add_variant_qc
    variant_qc_annotations = args.variant_qc_annotations
    pass_only = args.pass_only
    n_partitions = args.num_partitions
    vcf = args.vcf
    header_dict = HEADER_DICT

    try:
        if (add_variant_qc or pass_only) and not split_multi:
            raise ValueError(
                "Cannot annotate variant QC annotations or filter to PASS variants on"
                " an unsplit dataset."
            )

        if n_partitions and vds_output:
            raise ValueError(
                "Number of partitions argument is not supported for VDS output at this"
                " time."
            )

        vds = get_gnomad_v3_vds(
            remove_hard_filtered_samples=hard_filter,
            release_only=release_only,
        )

        if args.test:
            vds = hl.vds.VariantDataset(
                vds.reference_data._filter_partitions(range(2)),
                vds.variant_data._filter_partitions(range(2)),
            )

        subset_ht = hl.import_table(args.subset_samples).key_by("s")
        logger.info("Subsetting to %d samples.", subset_ht.count())
        vds = hl.vds.filter_samples(vds, subset_ht, remove_dead_alleles=True)
        logger.info(
            "Final number of samples being kept in the VDS: %d.",
            vds.variant_data.count_cols(),
        )

        if split_multi:
            logger.info("Splitting multiallelics.")
            vd = vds.variant_data
            vd = vd.annotate_rows(
                n_unsplit_alleles=hl.len(vd.alleles),
                mixed_site=(hl.len(vd.alleles) > 2)
                & hl.any(lambda a: hl.is_indel(vd.alleles[0], a), vd.alleles[1:])
                & hl.any(lambda a: hl.is_snp(vd.alleles[0], a), vd.alleles[1:]),
            )
            vds = hl.vds.split_multi(
                hl.vds.VariantDataset(vds.reference_data, vd), filter_changed_loci=True
            )

            if add_variant_qc or pass_only:
                vd = vds.variant_data
                ht = release_sites().ht()
                if pass_only:
                    logger.info("Filtering to variants that passed variant QC.")
                    vd = vd.filter_rows(ht[vd.row_key].filters.length() == 0)
                if add_variant_qc:
                    logger.info("Adding variant QC annotations.")
                    vd = vd.annotate_rows(
                        **make_variant_qc_annotations_dict(
                            vd.row_key, variant_qc_annotations
                        )
                    )
                    logger.info("Dropping VRS as it is VCF incompatible as structured.")
                    vd = vd.annotate_rows(
                        info=vd.info.drop("vrs")
                    )  # VRS is a struct and cannot be exported to VCF, so we drop it here. Should flatten it in the future.
                vds = hl.vds.VariantDataset(vds.reference_data, vd)

        if vcf or dense_mt or subset_call_stats:
            logger.info("Densifying VDS...")
            mt = hl.vds.to_dense_mt(vds)

        if dense_mt:
            partitions = compute_partitions(mt) if not n_partitions else n_partitions
            logger.info(
                "Naive coalescing to %i partitions and writing output MT.", partitions
            )
            mt = mt.naive_coalesce(n_partitions).checkpoint(
                f"{output_path}/subset.mt", overwrite=args.overwrite
            )

        if vds_output:
            vds.write(f"{output_path}/subset.vds", overwrite=args.overwrite)

        if vcf:  # TODO: Do we want to add in sharded vcf delivery using n_partitions?
            logger.info("Adjusting VCF incompatible types...")
            ht = adjust_vcf_incompatible_types(
                mt.select_rows("info").rows(), pipe_delimited_annotations=[]
            )
            mt = mt.annotate_rows(info=ht[mt.row_key].info)
            logger.info("Exporting subset to VCF.")
            mt = mt.drop("gvcf_info")
            hl.export_vcf(
                mt,
                f"{output_path}/subset.vcf.bgz",
                metadata=header_dict,
                parallel="header_per_shard",
                tabix=True,
            )

        if args.add_sample_qc:
            logger.info("Exporting subset's metadata.")
            meta = metadata.ht()
            meta = meta.semi_join(mt.cols())
            if hard_filter:
                meta = meta.filter(hl.is_missing(hard_filtered_samples.ht()[meta.key]))
            if release_only:
                meta = meta.filter(meta.release)
            meta.export(f"{output_path}/metadata.tsv.bgz")

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log("gs://gnomad-tmp/subset_v3/subset_v3.log")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script subsets gnomAD v3 using a list of samples."
    )
    parser.add_argument(
        "-s",
        "--subset-samples",
        help="Path to a text file with sample IDs for subsetting and a header: s.",
        required=True,
    )
    parser.add_argument(
        "--split-multi",
        help="Whether to split multi-allelic variants.",
        action="store_true",
    )
    parser.add_argument(
        "--dense-mt", help="Whether to make a dense MT.", action="store_true"
    )
    parser.add_argument(
        "--vcf", help="Whether to make a subset VCF.", action="store_true"
    )
    parser.add_argument("--vds", help="Whether to make a VDS.", action="store_true")
    parser.add_argument(
        "--output-path",
        help="Output file path for subsetted file, do not include file.",
        required=True,
    )
    parser.add_argument(
        "--num-partitions", help="Number of partitions in output file.", type=int
    )
    parser.add_argument(
        "--hard-filter-samples",
        help="Removes hard filtered samples.",
        action="store_true",
    )
    parser.add_argument(
        "--release-only", help="Keep release samples only.", action="store_true"
    )
    parser.add_argument(
        "--subset-call-stats",
        help="Adds subset callstats, AC, AN, AF",
        action="store_true",
    )
    parser.add_argument(
        "--add-sample-qc",
        help="Export sample QC from metadata file.",
        action="store_true",
    )
    parser.add_argument(
        "--add-variant-qc",
        help=(
            "Annotate exported file with gnomAD's variant QC annotations. Defaults to"
            " all annotations if a subset of annotations are not specified using the"
            " --variant-qc-fields arg"
        ),
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite outputs if they exist.",
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
    parser.add_argument("--test", help="Run on test dataset.", action="store_true")
    parser.add_argument(
        "--variant-qc-annotations",
        nargs="+",
        type=str,
        choices=VARIANT_QC_ANNOTATIONS,
        default=VARIANT_QC_ANNOTATIONS,
        help=(
            "Variant QC annotations to add to the output file. Defaults to all"
            " annotations."
        ),
    )
    args = parser.parse_args()
    main(args)
