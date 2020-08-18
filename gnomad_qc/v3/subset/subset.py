import argparse
import logging
import sys

import hail as hl

from gnomad_qc.v3.resources.meta import meta as metadata
from gnomad_qc.v3.resources.annotations import get_info
from gnomad_qc.v3.resources.raw import get_gnomad_v3_mt
from gnomad.resources.grch38.gnomad import GENOME_POPS
from gnomad.utils.filtering import subset_samples_and_variants


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

HEADER_DICT = {
    "format": {
        "GQ": {"Description": "Genotype Quality"},
        "SB": {
            "Description": "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias"
        },
        "AD": {
            "Description": "Allelic depths for the ref and alt alleles in the order listed"
        },
        "PID": {
            "Description": "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group"
        },
        "GT": {"Description": "Genotype"},
        "MIN_DP": {"Description": "Minimum DP observed within the GVCF block"},
        "PGT": {
            "Description": "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another"
        },
        "PL": {
            "Description": "Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification"
        },
        "DP": {"Description": "Approximate read depth"},
        "RGQ": {"Description": ""},
    },
    "info": {
        "AC": {"Description": "Alternate Allele count in gnomAD post-filtering"},
        "AC_raw": {"Description": "Raw alternate allele count in gnomAD"},
        "AS_ReadPosRankSum": {
            "Description": "allele specific Z-score from Wilcoxon rank sum test of each Alt vs. Ref read position bias"
        },
        "AS_MQRankSum": {
            "Description": "Allele-specifc Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"
        },
        "AS_RAW_MQ": {
            "Description": "Allele-specifc raw root mean square of the mapping quality of reads across all samples"
        },
        "AS_QUALapprox": {
            "Description": "Allele-specifc sum of PL[0] values; used to approximate the QUAL score"
        },
        "AS_MQ_DP": {
            "Description": "Allele-specifc depth over variant samples for better MQ calculation"
        },
        "AS_VarDP": {
            "Description": "Allele-specific (informative) depth over variant genotypes -- including ref, RAW format"
        },
        "AS_MQ": {
            "Description": "Allele-specifc root mean square of the mapping quality of reads across all samples"
        },
        "AS_QD": {
            "Description": "Allele-specifc variant call confidence normalized by depth of sample reads supporting a variant"
        },
        "AS_FS": {
            "Description": "Allele-specifc Phred-scaled p-value of Fisher's exact test for strand bias"
        },
        "AS_SB_TABLE": {
            "Description": "Allele-specific forward/reverse read counts for strand bias tests"
        },
        "ReadPosRankSum": {
            "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias"
        },
        "MQRankSum": {
            "Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"
        },
        "RAW_MQ": {
            "Description": "Raw root mean square of the mapping quality of reads across all samples"
        },
        "QUALapprox": {
            "Description": "Sum of PL[0] values; used to approximate the QUAL score"
        },
        "MQ_DP": {
            "Description": "Depth over variant samples for better MQ calculation"
        },
        "VarDP": {"Description": "(informative) depth over variant genotypes"},
        "SB": {
            "Description": "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias."
        },
        "MQ": {
            "Description": "Root mean square of the mapping quality of reads across all samples"
        },
        "QD": {
            "Description": "Variant call confidence normalized by depth of sample reads supporting a variant"
        },
        "FS": {
            "Description": "Phred-scaled p-value of Fisher's exact test for strand bias"
        },
    },
}


def format_info_for_vcf(ht) -> hl.Table:
    """Formats gnomAD MT fields with types that are not condusive to vcf export
    :param ht: Table to reformat
    :return: Reformatted Table
    """
    logger.info("Formatting the info hail table for VCF export")
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            ht = ht.annotate(
                info=ht.info.annotate(**{f: hl.int32(hl.min(2 ** 31 - 1, ht.info[f]))}) # add warning that number being capped
            )
        elif ft == hl.dtype("float64"):
            ht = ht.annotate(
                info=ht.info.annotate(
                    **{f: hl.float32(hl.min(2 ** 31 - 1, ht.info[f]))}
                )
            )
        elif ft == hl.dtype("array<int64>"):
            ht = ht.annotate(
                info=ht.info.annotate(
                    **{f: ht.info[f].map(lambda x: hl.int32(hl.min(2 ** 31 - 1, x)))}
                )
            )
        elif ft == hl.dtype("array<float64>"):
            ht = ht.annotate(
                info=ht.info.annotate(
                    **{f: ht.info[f].map(lambda x: hl.float32(hl.min(2 ** 31 - 1, x)))}
                )
            )
    return ht


def compute_partitions(mt, entry_size=3.5, partition_size=128000000) -> int:
    """Computes a very rough estimate for the optimal number of partitions in a MT
     using a the hail recommended partition size of 128MB and the rough estimate 
     of 3.5B per entry in the gnomAD sparse MT.

    :param mt: MatrixTable that will be resized
    :param entry_size: Average size in bytes that a single entry requires , defaults to 3.5
    :param partition_size: Amount of data per partition. Hail recommends 128MB, defaults to 128000000
    :return: number of partitions for resizing MT  
    """
    rows, columns = mt.count()
    mt_disk_est = rows * columns * entry_size
    n_partitions = (
        int(mt_disk_est / partition_size) if (mt_disk_est / partition_size) > 20 else 20 # make max
    )
    return n_partitions # bring to gnomad_methods -> parameterize


def main(args):
    hl.init(log="/subset.log", default_reference="GRCh38")
    pop = args.pop
    subset = args.subset
    hard_filter = args.hard_filter_samples
    release = args.release
    sparse = args.sparse
    header_dict = HEADER_DICT

    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, remove_hard_filtered_samples=hard_filter, release_only=release
    )
    info_ht = get_info().ht()
    info_ht = format_info_for_vcf(info_ht)

    if "SB" in info_ht.info and not isinstance(
        info_ht.info.SB, hl.expr.ArrayNumericExpression
    ):
        info_ht = info_ht.annotate(
            info=info_ht.info.annotate(SB=hl.flatten(info_ht.info.SB))
        )

    meta = metadata.ht()

    if pop and pop in GENOME_POPS:
        logger.info(f"Subsetting samples to {pop} population")
        mt = mt.filter_cols(meta[mt.col_key].pop == pop)

    if subset:
        gt_expr = "LGT" if sparse else "GT" #assumes matrix table is not split - infer 
        mt = subset_samples_and_variants(mt, subset, sparse=sparse, gt_expr=gt_expr)

    if args.add_meta:
        meta = meta.semi_join(mt.cols())
        meta.export(f"{args.output_path}/metadata.tsv.bgz")

    if sparse:
        mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
        mt = hl.experimental.densify(mt)
    else:
        mt = hl.split_multi_hts(mt)

    mt = mt.filter_rows(hl.len(mt.alleles) > 1)
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    mt = mt.drop(mt.gvcf_info)
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    if args.subset_call_stats:
        mt = mt.annotate_rows(subset_callstats=hl.agg.call_stats(mt.GT,mt.alleles))     
        mt = mt.annotate_rows(info=mt.info.annotate(subset_AC=mt.subset_callstats.AC[1],
                                                    subset_AN=mt.subset_callstats.AN,
                                                    subset_AF=hl.float32(mt.subset_callstats.AF[1])))
        mt = mt.drop('subset_callstats')
        header_dict['info'].update({"subset_AC": {"Description": "Alternate Allele count in subset post-filtering"} , 
                                    "subset_AN": {"Description": "Ref and Alternate allele number in subset post-filtering"},
                                    "subset_AF": {"Description": "Alternate Allele frequency in subset post-filtering"}})
    
    if args.add_gnomad_freqs:
        from gnomad_qc.v3.resources.annotations import freq
        freq = freq.ht()
        mt = mt.annotate_rows(info=mt.info.annotate(gnomAD_AC_adj=freq[mt.row_key].freq[0].AC, 
                                                    gnomAD_AN_adj=freq[mt.row_key].freq[0].AN,
                                                    gnomAD_AF_adj=hl.float32(freq[mt.row_key].freq[0].AF),
                                                    gnomAD_AC_raw=freq[mt.row_key].freq[1].AC, 
                                                    gnomAD_AN_raw=freq[mt.row_key].freq[1].AN,
                                                    gnomAD_AF_raw=hl.float32(freq[mt.row_key].freq[1].AF)))
        header_dict['info'].update({"gnomAD_AC_adj": {"Description": "High quality alternate allele count in gnomAD post-filtering"},
                                    "gnomAD_AN_adj": {"Description": "High quality allele number. The total number of called alleles in gnomAD post-filtering"},
                                    "gnomAD_AF_adj": {"Description": "High quality alternate allele frequency in gnomAD post-filtering"},
                                    "gnomAD_AC_raw": {"Description": "Raw alternate allele count in gnomAD post-filtering"}, 
                                    "gnomAD_AN_raw": {"Description": "Raw allele number. The total number of called alleles in gnomAD post-filtering"},
                                    "gnomAD_AF_raw": {"Description": "Raw alternate allele frequency in gnomAD post-filtering"}})

    partitions = (
        compute_partitions(mt) if not args.num_vcf_shards else args.num_vcf_shards
    )
    logger.info(f"Naive coalescing to {partitions}")

    mt = mt.naive_coalesce(partitions)
    if args.checkpoint_path:
        mt = mt.checkpoint(args.checkpoint_path, overwrite=True)

    logger.info("Exporting VCF")
    hl.export_vcf(
        mt,
        f"{args.output_path}/sharded_vcf.bgz",
        parallel="header_per_shard",
        metadata=HEADER_DICT,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script subsets gnomAD using a list of samples or population"
    )

    parser.add_argument(
        "-s", "--subset", help="Path to subset table of sample IDs with header: s"
    )
    parser.add_argument("--pop", help="Subset from specific population") # add choices for gnomad pops
    parser.add_argument(
        "--add-meta", help="Pull sample subset metadata", action="store_true"
    )
    parser.add_argument(
        "--output-path", help="Output file path for subsetted VCF, do not include file", required=True,
    )
    parser.add_argument(
        "--num-vcf-shards", help="Number of shards in output VCF", type=int # clarify docstring
    )
    parser.add_argument("--checkpoint-path", help="Path to final mt checkpoint")
    parser.add_argument("--hard-filter-samples", help="Removes hard filtered samples", required=True, action="store_true")
    parser.add_argument("--release", help="Keep release samples only", required=True, action="store_true")
    parser.add_argument(
        "--subset-call-stats", help="Adds subset callstats, AC, AN, AF"
    )
    parser.add_argument(
        "--add-gnomad-freqs", help="Adds gnomAD's adj and raw callstats"
    )
    parser.add_argument("--sparse", help="Use if source MT is sparse", action="store_true")
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    args = parser.parse_args()

    if not (args.subset or args.pop):
        sys.exit('Error: At least subset path or population must be specified')

    main(args)
