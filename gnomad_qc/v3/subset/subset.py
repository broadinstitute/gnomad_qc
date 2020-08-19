import argparse
import logging
import sys

import hail as hl

from gnomad_qc.v3.resources.meta import meta as metadata
from gnomad_qc.v3.resources.annotations import freq, get_info
from gnomad_qc.v3.resources.raw import get_gnomad_v3_mt
from gnomad.resources.grch38.gnomad import GENOME_POPS
from gnomad.utils.filtering import subset_samples_and_variants
from gnomad.utils.annotations import annotate_adj

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

HEADER_DICT = {
    "format": {
        "GT": {"Description": "Genotype", "Number": "1", "Type": "String"},
        "AD": {
            "Description": "Allelic depths for the ref and alt alleles in the order listed",
            "Number": "R",
            "Type": "Integer",
        },
        "DP": {
            "Description": "Approximate read depth",
            "Number": "1",
            "Type": "Integer",
        },
        "GQ": {
            "Description": "Phred-scaled confidence that the genotype assignment is correct. Value is the difference between the second lowest PL and the lowest PL (always normalized to 0).",
            "Number": "1",
            "Type": "Integer",
        },
        "MIN_DP": {
            "Description": "Minimum DP observed within the GVCF block",
            "Number": "1",
            "Type": "Integer",
        },
        "PGT": {
            "Description": "Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another",
            "Number": "1",
            "Type": "String",
        },
        "PID": {
            "Description": "Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group",
            "Number": "1",
            "Type": "String",
        },
        "PL": {
            "Description": "Normalized, phred-scaled likelihoods for genotypes as defined in the VCF specification",
            "Number": "G",
            "Type": "Integer",
        },
        "SB": {
            "Description": "Per-sample component statistics which comprise the Fisher's exact test to detect strand bias. Values are: depth of reference allele on forward strand, depth of reference allele on reverse strand, depth of alternate allele on forward strand, depth of alternate allele on reverse strand.",
            "Number": "4",
            "Type": "Integer",
        },
        "RGQ": {"Description": ""},
        },
    "info": {
        "AS_ReadPosRankSum": {
            "Description": "Allele-specific Z-score from Wilcoxon rank sum test of each alternate vs. reference read position bias"
        },
        "AS_MQRankSum": {
            "Description": "Allele-specific Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"
        },
        "AS_RAW_MQ": {
            "Description": "Allele-specific raw root mean square of the mapping quality of reads across all samples"
        },
        "AS_QUALapprox": {
            "Description": "Allele-specific sum of PL[0] values; used to approximate the QUAL score"
        },
        "AS_MQ_DP": {
            "Description": "Allele-specific depth over variant samples for better MQ calculation"
        },
        "AS_VarDP": {
            "Description": "Allele-specific (informative) depth over variant genotypes -- including ref, RAW format"
        },
        "AS_MQ": {
            "Description": "Allele-specific root mean square of the mapping quality of reads across all samples"
        },
        "AS_QD": {
            "Description": "Allele-specific variant call confidence normalized by depth of sample reads supporting a variant"
        },
        "AS_FS": {
            "Description": "Allele-specific Phred-scaled p-value of Fisher's exact test for strand bias"
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
    """
    Formats gnomAD MT fields with types that are not condusive to vcf export

    :param ht: Table to reformat
    :return: Reformatted Table
    """
    logger.info("Formatting the info hail table for VCF export")
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            ht = ht.annotate(
                info=ht.info.annotate(**{f: hl.int32(hl.min(2 ** 31 - 1, ht.info[f]))}) # TODO: add warning that number being capped
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


def compute_partitions(mt, entry_size=3.5, partition_size=128000000, min_partitions=20) -> int: #TODO: move to gnomad_methods
    """
    Computes a very rough estimate for the optimal number of partitions in a MT
     using a the hail recommended partition size of 128MB and the rough estimate 
     of 3.5B per entry in the gnomAD sparse MT.

    :param mt: MatrixTable that will be resized
    :param entry_size: Average size in bytes that a single entry requires , defaults to 3.5
    :param partition_size: Amount of data per partition. Hail recommends 128MB, defaults to 128000000
    :param min_partitions: Minimum number of partitions for naive_coalesce, defaults to 20
    :return: number of partitions for resizing MT  
    """
    rows, columns = mt.count()  
    mt_disk_est = rows * columns * entry_size
    n_partitions = hl.max(int(mt_disk_est / partition_size), min_partitions)
    return n_partitions 


def main(args):
    hl.init(log="/subset.log", default_reference="GRCh38")
    pop = args.pop
    subset = args.subset
    hard_filter = args.hard_filter_samples
    release_only = args.release_only
    header_dict = HEADER_DICT

    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, remove_hard_filtered_samples=hard_filter, release_only=release_only
    )
    info_ht = get_info().ht()
    info_ht = info_ht.annotate(info=info_ht.info.drop('AC', 'AC_raw')) # Note: AC and AC_raw are computed over all gnomAD samples
    info_ht = format_info_for_vcf(info_ht)

    # Flatten SB if it is an array of arrays
    if "SB" in info_ht.info and not isinstance(
        info_ht.info.SB, hl.expr.ArrayNumericExpression
    ):
        info_ht = info_ht.annotate(
            info=info_ht.info.annotate(SB=hl.flatten(info_ht.info.SB))
        )

    meta = metadata.ht()

    if pop:
        logger.info(f"Subsetting samples to {pop} population")
        mt = mt.filter_cols(meta[mt.col_key].pop == pop)

    if subset:
        mt = subset_samples_and_variants(mt, subset, sparse=True, gt_expr="LGT") 

    if args.export_meta:
        meta = meta.semi_join(mt.cols())
        meta.export(f"{args.output_path}/metadata.tsv.bgz")

    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = hl.experimental.densify(mt)

    
    mt = mt.filter_rows(hl.len(mt.alleles) > 1) # Note: This step is sparse-specific, removing monoallelic sites after densifying
    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))

    mt = mt.drop(mt.gvcf_info) # Note: gvcf_info is sparse-specific 
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    if args.subset_call_stats:
        mt = mt.annotate_rows(subset_callstats_raw=hl.agg.call_stats(mt.GT, mt.alleles))     
        mt = mt.annotate_rows(info=mt.info.annotate(AC_raw=mt.subset_callstats_raw.AC[1],
                                                    AN_raw=mt.subset_callstats_raw.AN,
                                                    AF_raw=hl.float32(mt.subset_callstats_raw.AF[1])))
        mt = mt.drop('subset_callstats_raw')
        mt = annotate_adj(mt)
        mt = mt.annotate_rows(subset_callstats_adj=hl.agg.filter(mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)))
        mt = mt.annotate_rows(info=mt.info.annotate(AC=mt.subset_callstats_adj.AC[1],
                                                    AN=mt.subset_callstats_adj.AN,
                                                    AF=hl.float32(mt.subset_callstats_adj.AF[1])))
        mt = mt.drop('subset_callstats_adj')
        header_dict['info'].update({
            "AC_raw": 
                {
                "Number": "A",
                "Description": "Alternate allele count in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "AN_raw": 
                {
                "Number": "1",
                "Description": "Total number of alleles in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "AF_raw": {
                "Number": "A",
                "Description": "Alternate allele frequency in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "AC": {
                "Number": "A",
                "Description": "Alternate allele count in subset after filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "AN": {
                "Number": "1",
                "Description": "Total number of alleles in subset after filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)"
                },
            "AF": {
                "Number": "A",
                "Description": "Alternate allele frequency in subset after filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                }
            }
        )
    
    if args.add_gnomad_freqs:
        freq_ht = freq.ht()
        mt = mt.annotate_rows(info=mt.info.annotate(gnomAD_AC=freq_ht[mt.row_key].freq[0].AC, 
                                                    gnomAD_AN=freq_ht[mt.row_key].freq[0].AN,
                                                    gnomAD_AF=hl.float32(freq_ht[mt.row_key].freq[0].AF),
                                                    gnomAD_AC_raw=freq_ht[mt.row_key].freq[1].AC, 
                                                    gnomAD_AN_raw=freq_ht[mt.row_key].freq[1].AN,
                                                    gnomAD_AF_raw=hl.float32(freq_ht[mt.row_key].freq[1].AF)))
        header_dict['info'].update({
            "gnomAD_AC": {
                "Number": "A",
                "Description": "High quality alternate allele count in gnomAD after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "gnomAD_AN": {
                "Number": "1",
                "Description": "High quality allele number. The total number of called alleles in gnomAD after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "gnomAD_AF": {
                "Number": "A",
                "Description": "High quality alternate allele frequency in gnomAD after filtering out low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "gnomAD_AC_raw": {
                "Number": "A",
                "Description": "Raw alternate allele count in gnomAD before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                }, 
            "gnomAD_AN_raw": {
                "Number": "1",
                "Description": "Raw allele number. The total number of called alleles in gnomAD before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                },
            "gnomAD_AF_raw": {
                "Number": "A",
                "Description": "Raw alternate allele frequency in gnomAD before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
                }
            }
        )

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
        metadata=header_dict,
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script subsets gnomAD using a list of samples or population"
    )

    parser.add_argument(
        "-s", "--subset", help="Path to a text file with sample IDs for subsetting and a header: s"
    )
    parser.add_argument("--pop", help=f"Subset to a specific population", choices=GENOME_POPS) 
    parser.add_argument(
        "--export-meta", help="Pull sample subset metadata and export to a .tsv", action="store_true"
    )
    parser.add_argument(
        "--output-path", help="Output file path for subsetted VCF, do not include file", required=True,
    )
    parser.add_argument("--num-vcf-shards", help="Number of shards in output VCF", type=int)
    parser.add_argument("--checkpoint-path", help="Path to final mt checkpoint")
    parser.add_argument("--hard-filter-samples", help="Removes hard filtered samples", required=True, action="store_true")
    parser.add_argument("--release-only", help="Keep release samples only", required=True, action="store_true")
    parser.add_argument("--subset-call-stats", help="Adds subset callstats, AC, AN, AF")
    parser.add_argument("--add-gnomad-freqs", help="Adds gnomAD's adj and raw callstats")
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
