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


def format_info_for_vcf(info_ht) -> hl.Table:
    logger.info(f"Formatting the info hail table for VCF export")
    for f, ft in info_ht.info.dtype.items():
        if ft == hl.dtype("int64"):
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{f: hl.int32(hl.min(2 ** 31 - 1, info_ht.info[f]))}
                )
            )
        elif ft == hl.dtype("float64"):
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{f: hl.float32(hl.min(2 ** 31 - 1, info_ht.info[f]))}
                )
            )
        elif ft == hl.dtype("array<int64>"):
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{
                        f: info_ht.info[f].map(
                            lambda x: hl.int32(hl.min(2 ** 31 - 1, x))
                        )
                    }
                )
            )
        elif ft == hl.dtype("array<float64>"):
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{
                        f: info_ht.info[f].map(
                            lambda x: hl.float32(hl.min(2 ** 31 - 1, x))
                        )
                    }
                )
            )
        elif ft == hl.dtype("array<array<int32>>"):
            info_ht = info_ht.annotate(
                info=info_ht.info.annotate(
                    **{f: hl.flatmap(lambda x: x, info_ht.info[f])}
                )
            )
    return info_ht


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
    meta = hl.read_table(metadata)

    if pop and pop in GENOME_POPS:
        logger.info("Subsetting samples to {pop} population")
        mt = mt.annotate_cols(pop=meta[mt.col_key].pop)
        mt = mt.filter_cols(mt.pop == pop)
        mt = mt.cols().drop("pop")

    if subset:
        mt = subset_samples_and_variants(mt, subset, sparse=True)

    if args.pull_meta:
        meta = meta.semi_join(mt.cols())
        meta.export(f"{args.output_path}/metadata.tsv.bgz")

    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = hl.experimental.densify(mt)

    mt = mt.filter_rows(hl.agg.any(mt.GT.is_non_ref()))
    mt = mt.filter_rows(hl.len(mt.alleles) > 1)

    mt = mt.drop(mt.gvcf_info)
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    mt = mt.naive_coalesce(1000)
    hl.export_vcf(mt, f"{args.output_path}/sharded_vcf.bgz", parallel="header_per_shard")


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
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
