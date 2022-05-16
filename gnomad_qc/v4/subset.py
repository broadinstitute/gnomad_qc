import argparse
import logging

import hail as hl

# TODO: include this when we have the sample QC meta HT: from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

# TODO: Can remove when we use the full sample QC meta HT
from gnomad_qc.v4.resources.meta import project_meta


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
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
            "Description": "Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)"
        },
    },
}


def main(args):
    hl.init(log="/subset.log", default_reference="GRCh38")
    header_dict = HEADER_DICT
    test = args.test

    vds = get_gnomad_v4_vds(n_partitions=args.n_partitions)

    if test:
        vds = hl.vds.variant_dataset.VariantDataset(
            vds.reference_data._filter_partitions(range(2)),
            vds.variant_data._filter_partitions(range(2)),
        )

    meta_ht = project_meta.ht()
    subset_ht = hl.import_table(args.subset or args.subset_workspaces)

    if args.subset_workspaces:
        terra_workspaces = hl.literal(subset_ht.terra_workspace.lower().collect())
        subset_ht = meta_ht.filter(
            terra_workspaces.contains(meta_ht.project_meta.terra_workspace)
        ).select()
        logger.info("Keeping %d samples using the list of terra workspaces.", subset_ht.count())

    if args.include_ukb_200k:
        ukb_subset_ht = meta_ht.filter(
            (meta_ht.project_meta.terra_workspace == "ukbb_wholeexomedataset")
            & hl.literal(UKB_BATCHES_FOR_INCLUSION).contains(
                meta_ht.project_meta.ukb_meta.ukb_batch
            )
        ).select()
        subset_ht = subset_ht.union(ukb_subset_ht)
        logger.info("Keeping %d samples after inclusion of the UKB 200K subset.", subset_ht.count())

    vds = hl.vds.filter_samples(vds, subset_ht, remove_dead_alleles=True)

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
        meta_ht = meta_ht.key_by(s=hl.coalesce(meta_ht.project_meta.ukb_meta.eid_31063, meta_ht.s))

    if args.vds:
        vds.write(f"{args.output_path}/subset.vds", overwrite=args.overwrite)

    if args.export_meta:
        logger.info("Exporting metadata")
        meta_ht = meta_ht.semi_join(vds.variant_data.cols())
        # TODO: Dropping the whole ukb_meta struct, but should we keep pop and sex inference if allowed?
        meta_ht = meta_ht.annotate(project_meta=meta_ht.project_meta.drop("ukb_meta"))
        meta_ht.export(f"{args.output_path}/metadata.tsv.bgz")

    if args.split_multi:
        logger.info("Splitting multi-allelics and densifying")
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    if args.vcf or args.dense_mt:
        mt = hl.vds.to_dense_mt(vds)

    if args.dense_mt:
        mt.write(f"{args.output_path}/subset.mt", overwrite=args.overwrite)

    # TODO: add num-vcf-shards where no sharding happens if this is not set
    if args.vcf:
        mt = mt.drop("gvcf_info")
        hl.export_vcf(
            mt,
            f"{args.output_path}.bgz",
            metadata=header_dict,
        )


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
        "--subset",
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
    parser.add_argument(
        "--export-meta",
        help="Pull sample subset metadata and export to a .tsv.",
        action="store_true",
    )
    parser.add_argument(
        "--output-path",
        help="Output file path for subsetted VDS/VCF/MT, do not include file.",
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
