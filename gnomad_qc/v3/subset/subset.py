"""Script to filter the gnomAD v3 VDS to a subset of specified samples with optional sample QC and variant QC annotations."""
import argparse
import logging
from typing import Dict, List

import hail as hl
from gnomad.utils.annotations import get_adj_expr

from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds
from gnomad_qc.v3.resources.meta import meta as metadata
from gnomad_qc.v3.resources.release import release_sites
from gnomad_qc.v3.resources.sample_qc import hard_filtered_samples

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

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


def compute_partitions(
    mt, entry_size=3.5, partition_size=128000000, min_partitions=20
) -> int:  # TODO: move to gnomad_methods
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
        else:
            raise ValueError(
                f"{ann} is not an annotation in the release HT. Possible annotations"
                f" are {VARIANT_QC_ANNOTATIONS}"
            )
    return selected_variant_qc_annotations


def main(args):
    """Filter the gnomAD v3 VDS to a subset of specified samples."""
    hl.init(
        log="/subset.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp/subset_v3",
    )

    hard_filter = args.hard_filter_samples
    release_only = args.release_only
    output_path = args.output_path
    vds = args.vds
    dense_mt = args.dense_mt
    split_multi = args.split_multi
    subset_call_stats = args.subset_call_stats
    variant_qc_annotations = args.variant_qc_annotations
    pass_only = args.pass_only
    n_partitions = args.num_partitions

    if variant_qc_annotations and not split_multi:
        logger.warning("Cannot add variant QC annotations to unsplit multiallelics.")

    if n_partitions and vds:
        raise ValueError(
            "Number of partitions argument is not supported for VDS output at this"
            " time."
        )

    vds = get_gnomad_v3_vds(
        key_by_locus_and_alleles=True,
        remove_hard_filtered_samples=hard_filter,
        release_only=release_only,
    )

    if args.test:
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(range(2)),
            vds.variant_data._filter_partitions(range(2)),
        )

    subset_ht = hl.import_table(args.subset_samples)
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

    if dense_mt or subset_call_stats or variant_qc_annotations or pass_only:
        mt = hl.vds.to_dense_mt(vds)

        if variant_qc_annotations or pass_only:
            logger.info("Adding gnomAD variant QC annotations.")
            if not split_multi:
                logger.info(
                    "Only biallelic variants will be annotated with variant QC"
                    " annotations."
                )
            ht = release_sites().ht()
            if pass_only:
                mt = mt.filter_rows(ht[mt.row_key].filters.contains("PASS"))
            if variant_qc_annotations:
                mt = mt.annotate_rows(
                    **make_variant_qc_annotations_dict(
                        mt.row_key, variant_qc_annotations
                    )
                )
                if vds:
                    vd = vds.variant_data
                    vd = vd.annotate_rows(
                        **make_variant_qc_annotations_dict(
                            vd.row_key, variant_qc_annotations
                        )
                    )
                    vds = hl.vds.VariantDataset(vds.reference_data, vd)

        if subset_call_stats:
            logger.info("Adding subset callstats.")
            if not split_multi:
                mt = mt.annotate_entries(
                    GT=hl.experimental.lgt_to_gt(mt.LGT, mt.LA),
                    adj=get_adj_expr(
                        mt.LGT,
                        mt.GQ,
                        mt.DP,
                        mt.LAD,
                    ),
                )
            else:
                mt = mt.annotate_entries(
                    adj=get_adj_expr(
                        mt.GT,
                        mt.GQ,
                        mt.DP,
                        mt.AD,
                    ),
                )

            ht = mt.annotate_rows(
                subset_callstats_raw=hl.agg.call_stats(mt.GT, mt.alleles),
                subset_callstats_adj=hl.agg.filter(
                    mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)
                ),
            ).rows()
            ht = ht.select(
                info=hl.struct(
                    AC_raw=ht.subset_callstats_raw.AC[1:],
                    AN_raw=ht.subset_callstats_raw.AN,
                    AF_raw=ht.subset_callstats_raw.AF[1:],
                    nhomalt_raw=ht.subset_callstats_raw.homozygote_count[1:],
                    AC=ht.subset_callstats_adj.AC[1:],
                    AN=ht.subset_callstats_adj.AN,
                    AF=ht.subset_callstats_adj.AF[1:],
                    nhomalt=ht.subset_callstats_adj.homozygote_count[1:],
                )
            )
            mt = mt.annotate_rows(info=mt.info.annotate(**ht[mt.row_key].info))

            if args.vds:
                vd = vds.variant_data
                vd = vd.annotate_rows(info=vd.info.annotate(**ht[vd.row_key].info))
                vds = hl.vds.VariantDataset(vds.reference_data, vd)

    if dense_mt:
        partitions = compute_partitions(mt) if not n_partitions else n_partitions
        logger.info(
            "Naive coalescing to %i partitions and writing output MT.", partitions
        )
        mt.naive_coalesce(partitions).write(
            f"{output_path}/subset.mt", overwrite=args.overwrite
        )

    if vds:
        vds.write(f"{output_path}/subset.vds", overwrite=args.overwrite)

    if args.add_sample_qc:
        logger.info("Exporting subset's metadata.")
        meta = metadata.ht()
        meta = meta.semi_join(mt.cols())
        if hard_filter:
            meta = meta.filter(hl.is_missing(hard_filtered_samples.ht()[meta.key]))
        if release_only:
            meta = meta.filter(meta.release)
        meta.export(f"{output_path}/metadata.tsv.bgz")


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
    parser.add_argument("--add-sample-qc", help="Export sample QC from metadata file.")
    parser.add_argument(
        "--variant-qc-annotations",
        help="Annotate exported file with these gnomAD's variant QC annotations.",
        default=VARIANT_QC_ANNOTATIONS,
        type=list,
        nargs="+",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite outputs if they exist.",
        action="store_true",
    )
    parser.add_argument("--vds", help="Output subset file in the hail VDS format.")
    parser.add_argument(
        "--pass-only",
        help=(
            "Keep only the variants that passed variant QC, i.e. the filter field is"
            " PASS."
        ),
        action="store_true",
    )
    parser.add_argument("--test", help="Run on test dataset.")
    args = parser.parse_args()
    main(args)
