import argparse
import logging

from gnomad.utils.annotations import annotate_adj
from gnomad.utils.vcf import adjust_vcf_incompatible_types
import hail as hl

# TODO: include this when we have the sample QC meta HT: from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

# TODO: Can remove when we use the full sample QC meta HT
from gnomad_qc.v4.resources.meta import project_meta

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)

UKB_BATCHES_FOR_INCLUSION = {"100K", "150K", "200K"}

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
        "RGQ": {
            "Description": "Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)",
            "Number": "1",
            "Type": "Integer",
        },
    },
}

SUBSET_CALLSTATS_INFO_DICT = {
    "AC_raw": {
        "Number": "A",
        "Description": "Alternate allele count in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
    },
    "AN_raw": {
        "Number": "1",
        "Description": "Total number of alleles in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
    },
    "AF_raw": {
        "Number": "A",
        "Description": "Alternate allele frequency in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
    },
    "nhomalt_raw": {
        "Number": "A",
        "Description": "Count of homozygous individuals in subset before filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
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
    },
    "nhomalt": {
        "Number": "A",
        "Description": "Count of homozygous individuals in subset after filtering of low-confidence genotypes (GQ < 20; DP < 10; and AB < 0.2 for het calls)",
    },
}


def main(args):
    hl.init(log="/subset.log", default_reference="GRCh38")
    test = args.test
    output_path = args.output_path
    header_dict = HEADER_DICT

    vds = get_gnomad_v4_vds(
        n_partitions=args.n_partitions, remove_hard_filtered_samples=False
    )

    if test:
        vds = hl.vds.variant_dataset.VariantDataset(
            vds.reference_data._filter_partitions(range(2)),
            vds.variant_data._filter_partitions(range(2)),
        )

    meta_ht = project_meta.ht()
    subset_ht = hl.import_table(args.subset_samples or args.subset_workspaces)

    if args.subset_workspaces:
        terra_workspaces = hl.literal(subset_ht.terra_workspace.lower().collect())
        subset_ht = meta_ht.filter(
            terra_workspaces.contains(meta_ht.project_meta.terra_workspace)
        ).select()
        logger.info(
            "Keeping %d samples using the list of terra workspaces.", subset_ht.count()
        )

    if args.include_ukb_200k:
        ukb_subset_ht = meta_ht.filter(
            (meta_ht.project_meta.terra_workspace == "ukbb_wholeexomedataset")
            & hl.literal(UKB_BATCHES_FOR_INCLUSION).contains(
                meta_ht.project_meta.ukb_meta.ukb_batch
            )
            & ~meta_ht.project_meta.ukb_meta.ukb_withdraw
        ).select()
        subset_ht = subset_ht.union(ukb_subset_ht)
        logger.info(
            "Keeping %d samples after inclusion of the UKB 200K subset.",
            subset_ht.count(),
        )

    vds = hl.vds.filter_samples(vds, subset_ht, remove_dead_alleles=True)
    logger.info(
        "Final number of samples being kept in the VDS: %d.",
        vds.variant_data.count_cols(),
    )

    if args.include_ukb_200k:
        # TODO: add option to provide an application linking file as an argument. Default is ATGU ID
        vds = hl.vds.variant_dataset.VariantDataset(
            vds.reference_data.key_cols_by(
                s=hl.coalesce(
                    meta_ht[vds.reference_data.col_key].project_meta.ukb_meta.eid_31063,
                    vds.reference_data.col_key.s,
                )
            ),
            vds.variant_data.key_cols_by(
                s=hl.coalesce(
                    meta_ht[vds.variant_data.col_key].project_meta.ukb_meta.eid_31063,
                    vds.variant_data.col_key.s,
                )
            ),
        )
        meta_ht = meta_ht.key_by(
            s=hl.coalesce(meta_ht.project_meta.ukb_meta.eid_31063, meta_ht.s)
        )

    if args.split_multi:
        logger.info("Splitting multi-allelics")
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

    if args.vcf or args.dense_mt or args.subset_call_stats:
        logger.info("Densifying VDS")
        mt = hl.vds.to_dense_mt(vds)

        if args.subset_call_stats:
            logger.info("Adding subset callstats")
            mt = annotate_adj(mt)
            ht = mt.annotate_rows(
                subset_callstats_raw=hl.agg.call_stats(mt.GT, mt.alleles),
                subset_callstats_adj=hl.agg.filter(mt.adj, hl.agg.call_stats(mt.GT, mt.alleles))
            ).rows()
            ht = ht.select_rows(
                info=ht.info.annotate(
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
            ht = adjust_vcf_incompatible_types(ht)
            mt = mt.annotate_rows(info=ht.info.annotate(**ht[mt.row_key]))

            if args.vds:
                vd = vds.variant_dataset
                vd = vd.annotate_rows(info=ht.info.annotate(**ht[vd.row_key]))
                vds = hl.vds.VariantDataset(vds.reference_data, vd)

        if args.dense_mt:
            mt.write(f"{output_path}/subset.mt", overwrite=args.overwrite)

        # TODO: add num-vcf-shards where no sharding happens if this is not set
        if args.vcf:
            mt = mt.drop("gvcf_info")
            header_dict['info'] = SUBSET_CALLSTATS_INFO_DICT
            hl.export_vcf(
                mt,
                f"{output_path}.bgz",
                metadata=header_dict,
            )

    if args.vds:
        vds.write(f"{output_path}/subset.vds", overwrite=args.overwrite)

    if args.export_meta:
        logger.info("Exporting metadata")
        meta_ht = meta_ht.semi_join(vds.variant_data.cols())
        # TODO: Dropping the whole ukb_meta struct, but should we keep pop and sex inference if allowed?
        if args.keep_data_paths:
            data_to_drop = {"ukb_meta"}
        else:
            data_to_drop = {"ukb_meta", "cram", "gvcf"}

        meta_ht = meta_ht.annotate(
            project_meta=meta_ht.project_meta.drop(*data_to_drop)
        )
        meta_ht.export(f"{output_path}/metadata.tsv.bgz")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script subsets gnomAD using a list of samples or terra workspaces."
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
            "Path to a text file with Terra workspaces that should be included in the subset, must use a header of "
            "'terra_workspace'."
        ),
    )
    parser.add_argument(
        "--include-ukb-200k",
        help="Whether to include the 200K UK Biobank samples.",
        action="store_true",
    )
    parser.add_argument(
        "--vds", help="Whether to make a subset VDS.", action="store_true"
    )
    parser.add_argument(
        "--vcf", help="Whether to make a subset VCF.", action="store_true"
    )
    parser.add_argument(
        "--dense-mt", help="Whether to make a dense MT", action="store_true"
    )
    parser.add_argument(
        "--split-multi",
        help="Whether to split multi-allelic variants.",
        action="store_true",
    )
    parser.add_argument(
        "--n-partitions",
        help=(
            "Number of desired partitions for the subset VDS if --vds and/or MT if --dense-mt is set and/or the number "
            "of shards in the output VCF if --vcf is set. By default, there will be no change in partitioning."
        ),
        type=int,
    )
    parser.add_argument("--subset-call-stats", help="Adds subset callstats, AC, AN, AF, nhomalt", action="store_true")
    parser.add_argument(
        "--export-meta",
        help="Pull sample subset metadata and export to a .tsv.",
        action="store_true",
    )
    parser.add_argument(
        "--keep-data-paths",
        help="Keep CRAM and gVCF paths in the project metadata export.",
        action="store_true",
    )
    parser.add_argument(
        "--output-path",
        help="Output file path for subsetted VDS/VCF/MT, do not include file extension.",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    args = parser.parse_args()

    main(args)
