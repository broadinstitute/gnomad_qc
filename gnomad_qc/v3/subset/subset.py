import argparse
import logging
import hail as hl
from gnomad_qc.v3.resources.meta import meta as metadata
from gnomad_qc.v3.resources.annotations import get_info
from gnomad_qc.v3.resources.raw import get_gnomad_v3_mt
from gnomad.resources.grch38.gnomad import GENOME_POPS
from gnomad.utils.generic import subset_samples_and_variants


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("seqr_sample_qc")
logger.setLevel(logging.INFO)

INFO_DICT = {
    'info': {
        'AC': {"Description": "Alternate Allele count in gnomAD post-filtering"}, 
        'AC_raw': {"Description": "Raw alternate allele count in gnomAD"}, 
        'AS_ReadPosRankSum': {"Description": "allele specific Z-score from Wilcoxon rank sum test of each Alt vs. Ref read position bias"}, 
        'AS_MQRankSum': {"Description": "Allele-specifc Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"}, 
        'AS_RAW_MQ': {"Description": "Allele-specifc raw root mean square of the mapping quality of reads across all samples"}, 
        'AS_QUALapprox': {"Description": "Allele-specifc sum of PL[0] values; used to approximate the QUAL score"}, 
        'AS_MQ_DP': {"Description": "Allele-specifc depth over variant samples for better MQ calculation"}, 
        'AS_VarDP': {"Description": "Allele-specific (informative) depth over variant genotypes -- including ref, RAW format"},  
        'AS_MQ': {"Description": "Allele-specifc root mean square of the mapping quality of reads across all samples"},  
        'AS_QD': {"Description": "Allele-specifc variant call confidence normalized by depth of sample reads supporting a variant"},  
        'AS_FS': {"Description": "Allele-specifc Phred-scaled p-value of Fisher's exact test for strand bias"},  
        'AS_SB_TABLE': {"Description": "Allele-specific forward/reverse read counts for strand bias tests"}, 
        'ReadPosRankSum': {"Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read position bias"},  
        'MQRankSum': {"Description": "Z-score from Wilcoxon rank sum test of alternate vs. reference read mapping qualities"}, 
        'RAW_MQ': {"Description": "Raw root mean square of the mapping quality of reads across all samples"}, 
        'QUALapprox': {"Description": "Sum of PL[0] values; used to approximate the QUAL score"},  
        'MQ_DP': {"Description": "Depth over variant samples for better MQ calculation"},  
        'VarDP': {"Description": "(informative) depth over variant genotypes"},  
        'SB': {"Description": "Per-sample component statistics which comprise the Fisher's Exact Test to detect strand bias."}, 
        'MQ': {"Description": "Root mean square of the mapping quality of reads across all samples"},
        'QD': {"Description": "Variant call confidence normalized by depth of sample reads supporting a variant"},
        'FS': {"Description": "Phred-scaled p-value of Fisher's exact test for strand bias"},
    }
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
                info=ht.info.annotate(
                    **{f: hl.int32(hl.min(2 ** 31 - 1, ht.info[f]))}
                )
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
                    **{
                        f: ht.info[f].map(
                            lambda x: hl.int32(hl.min(2 ** 31 - 1, x))
                        )
                    }
                )
            )
        elif ft == hl.dtype("array<float64>"):
            ht = ht.annotate(
                info=ht.info.annotate(
                    **{
                        f: ht.info[f].map(
                            lambda x: hl.float32(hl.min(2 ** 31 - 1, x))
                        )
                    }
                )
            )
    return ht    


def compute_partitions(
        mt, 
        entry_size = 3.5,
        partition_size=128000000
) -> int:
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
    n_partitions = int(mt_disk_est/partition_size)
    return n_partitions


def main(args):
    hl.init(
        log="/subset.log", tmp_dir="hdfs:///subset.tmp/", default_reference="GRCh38"
    )
    pop = args.pop
    subset = args.subset
    mt = get_gnomad_v3_mt(
        key_by_locus_and_alleles=True, remove_hard_filtered_samples=False
    )
    info_ht = get_info().ht()
    info_ht = format_info_for_vcf(info_ht)

    if 'SB' in info_ht.info and not isinstance(info_ht.info.SB, hl.expr.ArrayNumericExpression):
        info_ht = info_ht.annotate(info=info_ht.info.annotate(SB=hl.flatten(info_ht.info.SB)))

    meta = metadata.ht()

    if pop and pop in GENOME_POPS:
        logger.info("Subsetting samples to {pop} population")
        mt = mt.annotate_cols(pop=meta[mt.col_key].pop)
        mt = mt.filter_cols(mt.pop == pop)
        mt = mt.cols().drop("pop")

    if subset:
        mt = subset_samples_and_variants(mt, subset, sparse=True, gt_expr="LGT")

    if args.pull_meta:
        meta = meta.semi_join(mt.cols())
        meta.export(f"{args.output_path}/metadata.tsv.bgz")

    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = hl.experimental.densify(mt)

    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

    mt = mt.drop(mt.gvcf_info)
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    partitions = args.num_vcf_shards if args.num_vcf_shards else compute_partitions(mt)
    mt = mt.naive_coalesce(partitions)

    hl.export_vcf(mt, f"{args.output_path}/sharded_vcf.bgz", parallel="header_per_shard", metadata=INFO_DICT)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="This script subsets gnomAD using a list of samples or population"
    )

    parser.add_argument(
        "-s", "--subset", help="Path to subset table of sample IDs with header: s"
    )
    parser.add_argument("-pop", help="Subset from specific population")
    parser.add_argument(
        "--pull-meta", help="Pull sample subset metadata", action="store_true"
    )
    parser.add_argument(
        "--output-path", help="Output file path for subsetted VCF", required=True,
    )
    parser.add_argument("--num-vcf-shards", help="Number of shards in output VCF")
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
