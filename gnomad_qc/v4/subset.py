"""Script to filter the gnomAD v4 VariantDataset to a subset of specified samples."""

import argparse
import logging
from typing import Dict, List, Optional

import hail as hl
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("subset")
logger.setLevel(logging.INFO)


HEADER_DICT = {
    "format": {
        "GT": {"Description": "Genotype", "Number": "1", "Type": "String"},
        "AD": {
            "Description": (
                "Allelic depths for the ref and alt alleles in the order listed"
            ),
            "Number": "R",
            "Type": "Integer",
        },
        "DP": {
            "Description": "Approximate read depth",
            "Number": "1",
            "Type": "Integer",
        },
        "GQ": {
            "Description": (
                "Phred-scaled confidence that the genotype assignment is correct. Value"
                " is the difference between the second lowest PL and the lowest PL"
                " (always normalized to 0)."
            ),
            "Number": "1",
            "Type": "Integer",
        },
        "MIN_DP": {
            "Description": "Minimum DP observed within the GVCF block",
            "Number": "1",
            "Type": "Integer",
        },
        "PGT": {
            "Description": (
                "Physical phasing haplotype information, describing how the alternate"
                " alleles are phased in relation to one another"
            ),
            "Number": "1",
            "Type": "String",
        },
        "PID": {
            "Description": (
                "Physical phasing ID information, where each unique ID within a given"
                " sample (but not across samples) connects records within a phasing"
                " group"
            ),
            "Number": "1",
            "Type": "String",
        },
        "PL": {
            "Description": (
                "Normalized, phred-scaled likelihoods for genotypes as defined in the"
                " VCF specification"
            ),
            "Number": "G",
            "Type": "Integer",
        },
        "SB": {
            "Description": (
                "Per-sample component statistics which comprise the Fisher's exact test"
                " to detect strand bias. Values are: depth of reference allele on"
                " forward strand, depth of reference allele on reverse strand, depth of"
                " alternate allele on forward strand, depth of alternate allele on"
                " reverse strand."
            ),
            "Number": "4",
            "Type": "Integer",
        },
        "RGQ": {
            "Description": (
                "Unconditional reference genotype confidence, encoded as a phred"
                " quality -10*log10 p(genotype call is wrong)"
            ),
            "Number": "1",
            "Type": "Integer",
        },
    },
}


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


def main(args):
    """Filter the gnomAD v4 VariantDataset to a subset of specified samples."""
    hl.init(
        log="/v4_subset.log",
        default_reference="GRCh38",
        tmp_dir=args.tmp_dir,
    )
    test = args.test
    data_type = args.data_type
    output_path = args.output_path
    header_dict = HEADER_DICT
    add_variant_qc = args.add_variant_qc
    variant_qc_annotations = args.variant_qc_annotations
    pass_only = args.pass_only
    split_multi = args.split_multi
    n_partitions = args.n_partitions
    vcf = args.vcf
    dense_mt = args.dense_mt
    vds = args.vds

    if vcf and not split_multi:
        raise ValueError(
            "VCF export without split multi is not supported at this time."
        )
    if (add_variant_qc or pass_only) and not split_multi:
        raise ValueError(
            "Cannot annotate variant QC annotations or filter to PASS variants on an"
            " unsplit dataset."
        )

    if not vcf and not dense_mt and not vds:
        raise ValueError(
            "At least one of --vcf, --dense-mt, or --vds must be specified."
        )

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

    if args.subset_samples:
        subset_ht = hl.import_table(args.subset_samples)
        subset_ht = subset_ht.key_by("s")
        logger.info(
            "Imported %d samples from the subset file.",
            subset_ht.count(),
        )

    if args.subset_workspaces:
        terra_workspaces = hl.literal(subset_ht.terra_workspace.lower().collect())
        subset_ht = meta_ht.filter(
            terra_workspaces.contains(meta_ht.project_meta.terra_workspace)
        ).select()
        logger.info(
            "Keeping %d samples using the list of terra workspaces.", subset_ht.count()
        )

    vds = hl.vds.filter_samples(
        vds, subset_ht, remove_dead_alleles=False if split_multi else True
    )
    logger.info(
        "Final number of samples being kept in the VDS: %d.",
        vds.variant_data.count_cols(),
    )

    if vcf or dense_mt:
        logger.info("Densifying VDS...")
        mt = hl.vds.to_dense_mt(vds)
        mt = mt.checkpoint(new_temp_file("subset", "mt"))
        if split_multi:
            logger.info("Splitting multi-allelics...")
            mt = hl.experimental.sparse_split_multi(mt)
            mt = mt.filter_rows(mt.n_alt_alleles == 1)

    else:
        vd = vds.variant_data
        if split_multi:
            logger.info("Splitting multi-allelics...")
            vd = vd.annotate_rows(
                n_unsplit_alleles=hl.len(vd.alleles),
                mixed_site=(hl.len(vd.alleles) > 2)
                & hl.any(lambda a: hl.is_indel(vd.alleles[0], a), vd.alleles[1:])
                & hl.any(lambda a: hl.is_snp(vd.alleles[0], a), vd.alleles[1:]),
            )
            vds = hl.vds.split_multi(
                hl.vds.VariantDataset(vds.reference_data, vd), filter_changed_loci=True
            )
        else:
            logger.info(
                "Applying min_rep to the variant data MT because remove_dead_alleles in"
                " hl.vds.filter_samples may result in variants that do not have the minimum"
                " representation."
            )
            vds = hl.vds.VariantDataset(
                vds.reference_data,
                vd.key_rows_by(**hl.min_rep(vd.locus, vd.alleles)),
            )

    if add_variant_qc or pass_only:
        ds = vds.variant_data if not vcf or dense_mt else mt
        ht = release_sites(data_type=data_type).ht()
        if pass_only:
            logger.info("Filtering to variants that passed variant QC.")
            ds = ds.filter_rows(ht[vd.row_key].filters.length() == 0)
        if add_variant_qc:
            logger.info("Adding variant QC annotations.")
            ds = ds.annotate_rows(
                **make_variant_qc_annotations_dict(
                    ds.row_key, variant_qc_annotations, data_type
                )
            )
        if not (vcf or dense_mt):
            vds = hl.vds.VariantDataset(vds.reference_data, vd)

    if vcf or dense_mt:
        output_partitions = (
            args.output_partitions if args.output_partitions else mt.n_partitions()
        )
        mt = mt.naive_coalesce(output_partitions)

        if dense_mt:
            mt.naive_coalesce(output_partitions).write(
                f"{output_path}/subset.mt", overwrite=args.overwrite
            )

        # TODO: add num-vcf-shards where no sharding happens if this is not set.
        if vcf:
            mt = mt.drop("gvcf_info")
            mt = mt.transmute_rows(rsid=hl.str(";").join(mt.rsid))
            mt = mt.annotate_rows(info=mt.info.annotate(**mt.info.vrs).drop("vrs"))
            hl.export_vcf(
                mt,
                f"{output_path}/subset.vcf.bgz",
                metadata=header_dict,
                parallel="header_per_shard",
                tabix=True,
            )

    if vds:
        vds.write(f"{output_path}/subset.vds", overwrite=args.overwrite)

    if args.export_meta:
        logger.info("Exporting metadata")
        meta_ht = meta_ht.semi_join(vds.variant_data.cols())
        # TODO: Dropping the whole ukb_meta struct, but should we keep pop and sex inference if allowed? # noqa
        if data_type == "exomes":
            if args.keep_data_paths:
                data_to_drop = {"ukb_meta"}
            else:
                data_to_drop = {"ukb_meta", "cram", "gvcf"}

            meta_ht = meta_ht.annotate(
                project_meta=meta_ht.project_meta.drop(*data_to_drop)
            )
        meta_ht = meta_ht.checkpoint(
            f"{output_path}/metadata.ht", overwrite=args.overwrite
        )
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
            "Number of desired partitions for the subset VDS if --vds and/or MT if"
            " --dense-mt is set and/or the number of shards in the output VCF if --vcf"
            " is set. By default, there will be no change in partitioning."
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
