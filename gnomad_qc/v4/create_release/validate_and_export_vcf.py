# noqa: D100

import argparse
import logging
from copy import deepcopy
from pprint import pprint
from typing import Dict, List, Optional, Set

import hail as hl
from gnomad.assessment.validity_checks import (
    check_global_and_row_annot_lengths,
    pprint_global_anns,
    validate_release_t,
)
from gnomad.resources.grch38.gnomad import (
    HGDP_POPS,
    POPS,
    SUBSETS,
    TGP_POP_NAMES,
    TGP_POPS,
)
from gnomad.sample_qc.ancestry import POP_NAMES
from gnomad.utils.filtering import remove_fields_from_constant
from gnomad.utils.vcf import (
    ALLELE_TYPE_FIELDS,
    AS_FIELDS,
    AS_VQSR_FIELDS,
    FAF_POPS,
    HISTS,
    REGION_FLAG_FIELDS,
    SITE_FIELDS,
    VRS_FIELDS_DICT,
)
from gnomad.utils.vep import VEP_CSQ_HEADER, vep_struct_to_csq

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.v4.resources.basics import get_logging_path
from gnomad_qc.v4.resources.release import release_sites, validated_release_ht

# Add new site fields
NEW_SITE_FIELDS = [
    "monoallelic",
    "only_het",
    "transmitted_singleton",
]
SITE_FIELDS = deepcopy(SITE_FIELDS)
SITE_FIELDS.extend(NEW_SITE_FIELDS)
SITE_FIELDS = {
    "exomes": SITE_FIELDS + ["sibling_singleton"],
    "genomes": SITE_FIELDS,
}

# Add updatedVQSR fields
NEW_AS_VQSR_FIELDS = ["negative_train_site", "positive_train_site"]
AS_VQSR_FIELDSS = deepcopy(AS_VQSR_FIELDS)
AS_VQSR_FIELDS.extend(NEW_AS_VQSR_FIELDS)

# Update inbreeding coeff
AS_FIELDS = remove_fields_from_constant(AS_FIELDS, ["InbreedingCoeff"])
AS_FIELDS.append("inbreeding_coeff")

# Drop decoy, still doesn't exist on 38
REGION_FLAG_FIELDS = deepcopy(REGION_FLAG_FIELDS)
REGION_FLAG_FIELDS = remove_fields_from_constant(
    REGION_FLAG_FIELDS, ["decoy", "nonpar"]
)
REGION_FLAG_FIELDS.append("non_par")
REGION_FLAG_FIELDS = {
    "exomes": REGION_FLAG_FIELDS
    + [
        "fail_interval_qc",
        "outside_ukb_capture_region",
        "outside_broad_capture_region",
    ],
    "genomes": REGION_FLAG_FIELDS,
}

# Remove original alleles for containing non-releasable alleles
ALLELE_TYPE_FIELDS = deepcopy(ALLELE_TYPE_FIELDS)
MISSING_ALLELE_TYPE_FIELDS = ["original_alleles"]
ALLELE_TYPE_FIELDS = remove_fields_from_constant(
    ALLELE_TYPE_FIELDS, MISSING_ALLELE_TYPE_FIELDS
)
# Remove original alleles for containing non-releasable alleles
ALLELE_TYPE_FIELDS = {
    "exomes": ALLELE_TYPE_FIELDS,
    "genomes": remove_fields_from_constant(ALLELE_TYPE_FIELDS, ["has_star"]),
}

INSILICO_FIELDS = [
    "cadd",
    "revel_max",
    "spliceai_ds_max",
    "pangolin_largest_ds",
    "phylop",
    "sift_max",
    "polyphen_max",
]

GENOME_SUBSETS_TO_DROP = remove_fields_from_constant(
    deepcopy(SUBSETS["v3"]), ["hgdp", "tgp"]
)
SUBSETS = {
    "exomes": deepcopy(SUBSETS["v4"]),
    "genomes": remove_fields_from_constant(
        deepcopy(SUBSETS["v3"]), GENOME_SUBSETS_TO_DROP
    ),
}

# Exomes and genomes use the same pops for v4
POPS = deepcopy(POPS["v4"])
# Remove unnecessary pop names from POP_NAMES dict
POPS = {
    pop: POP_NAMES[pop] if pop != "remaining" else "Remaining individuals"
    for pop in POPS
}

SAMPLE_SUM_SETS_AND_POPS = {
    "exomes": {"non_ukb": POPS},
    "genomes": {"hgdp": HGDP_POPS, "tgp": TGP_POPS},
}

# Remove unnecessary pop names from FAF_POPS dict
FAF_POPS = {pop: POP_NAMES[pop] for pop in FAF_POPS}

# Row annotaions and their associated global annotations for length comparison
LEN_COMP_GLOBAL_ROWS = {
    "freq": ["freq_meta", "freq_index_dict", "freq_meta_sample_count"],
    "faf": ["faf_meta", "faf_index_dict"],
    "joint_freq": [
        "joint_freq_meta",
        "joint_freq_index_dict",
        "joint_freq_meta_sample_count",
    ],
    "joint_faf": ["joint_faf_meta", "joint_faf_index_dict"],
}

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def get_export_resources(
    overwrite: bool = False,
    data_type: str = "exomes",
    test: Optional[bool] = False,
) -> PipelineResourceCollection:
    """
    Get export resources.

    :param overwrite: Whether to overwrite existing files.
    :param data_type: Data type to get resources for. One of "exomes" or "genomes".
        Default is "exomes".
    :param test: Whether to use test resources.
    :return: Export resources.
    """
    export_pipeline = PipelineResourceCollection(
        pipeline_name="validate_and_export_vcf",
        overwrite=overwrite,
    )
    validate_release_ht = PipelineStepResourceCollection(
        "--validate-release-ht",
        input_resources={
            "create_release_sites_ht.py": {
                "release_ht": hl.read_table(
                    "gs://gnomad-tmp/gnomad.exomes.v4.0.qc_data/release/gnomad.genomes.sites.test.updated_101923.ht"
                )  # release_sites(data_type=data_type) #TODO: Restore once PR is ready for merge
            }
        },
        output_resources={
            "validated_ht": validated_release_ht(test=test, data_type=data_type)
        },
    )
    export_pipeline.add_steps(
        {
            "validate_release_ht": validate_release_ht,
        }
    )
    return export_pipeline


def filter_to_test(ht: hl.Table, num_partitions: int = 2) -> hl.Table:
    """
    Filter Table to `num_partitions` partitions on chr20, chrX, and chrY for testing.

    :param t: Input Table to filter.
    :param num_partitions: Number of partitions to grab from each chromosome.
    :return: Input Table filtered to `num_partitions` on chr20, chrX, and chrY.
    """
    logger.info(
        "Filtering to %d partitions on chr20, chrX, and chrY (for tests only)...",
        num_partitions,
    )
    ht_chr20 = hl.filter_intervals(ht, [hl.parse_locus_interval("chr20")])
    ht_chr20 = ht_chr20._filter_partitions(range(num_partitions))

    ht_chrx = hl.filter_intervals(ht, [hl.parse_locus_interval("chrX")])
    ht_chrx = ht_chrx._filter_partitions(range(num_partitions))

    ht_chry = hl.filter_intervals(ht, [hl.parse_locus_interval("chrY")])
    ht_chry = ht_chry._filter_partitions(range(num_partitions))
    ht = ht_chr20.union(ht_chrx, ht_chry)

    return ht


def unfurl_nested_annotations(
    ht: hl.Table,
    entries_to_remove: Set[str] = None,
    data_type: str = "exomes",
) -> [hl.expr.StructExpression, Set[str]]:
    """
    Create dictionary keyed by the variant annotation labels to be extracted from variant annotation arrays.

    The values of the returned dictionary are Hail Expressions describing how to access the corresponding values.

    :param ht: Table containing the nested variant annotation arrays to be unfurled.
    :param entries_to_remove: Optional Set of frequency entries to remove for vcf_export.
    :param data_type: Data type to unfurl nested annotations for. One of "exomes" or "genomes".
    :return: StructExpression containing variant annotations and their corresponding expressions and updated entries and set of frequency entries to remove
        to remove from the VCF.
    """
    expr_dict = {}

    # Unfurl freq index dict
    # Cycles through each key and index (e.g., k=adj_afr, i=31)
    logger.info("Unfurling freq data...")
    freq_idx = hl.eval(ht.freq_index_dict)
    expr_dict.update(
        {
            f"{f if f != 'homozygote_count' else 'nhomalt'}_{k}": ht.freq[i][f]
            for k, i in freq_idx.items()
            for f in ht.freq[0].keys()
        }
    )

    logger.info("Unfurling joint freq data...")
    joint_freq_idx = hl.eval(ht.joint_freq_index_dict)
    expr_dict.update(
        {
            f"joint_{f if f != 'homozygote_count' else 'nhomalt'}_{k}": ht.joint_freq[
                i
            ][f]
            for k, i in joint_freq_idx.items()
            for f in ht.joint_freq[0].keys()
        }
    )

    logger.info("Adding grpmax data...")
    grpmax_idx = ht.grpmax.gnomad if data_type == "exomes" else ht.grpmax
    grpmax_dict = {"grpmax": grpmax_idx.gen_anc}
    grpmax_dict.update(
        {
            f"{f if f != 'homozygote_count' else 'nhomalt'}_grpmax": grpmax_idx[f]
            for f in [f for f in grpmax_idx._fields if f != "gen_anc"]
        }
    )
    expr_dict.update(grpmax_dict)

    logger.info("Adding joint grpmax data...")
    joint_grpmax_idx = ht.joint_grpmax
    joint_grpmax_dict = {"joint_grpmax": joint_grpmax_idx.gen_anc}
    joint_grpmax_dict.update(
        {
            f"joint_{f if f != 'homozygote_count' else 'nhomalt'}_grpmax": (
                joint_grpmax_idx[f]
            )
            for f in [f for f in joint_grpmax_idx._fields if f != "gen_anc"]
        }
    )
    expr_dict.update(joint_grpmax_dict)

    logger.info("Unfurling faf data...")
    faf_idx = hl.eval(ht.faf_index_dict)
    expr_dict.update(
        {f"{f}_{k}": ht.faf[i][f] for f in ht.faf[0].keys() for k, i in faf_idx.items()}
    )

    logger.info("Unfurling joint faf data...")
    joint_faf_idx = hl.eval(ht.joint_faf_index_dict)
    expr_dict.update(
        {
            f"joint_{f}_{k}": ht.joint_faf[i][f]
            for f in ht.joint_faf[0].keys()
            for k, i in joint_faf_idx.items()
        }
    )

    logger.info("Unfurling joint fafmax data...")
    joint_fafmax_idx = ht.joint_fafmax
    joint_fafmax_dict = {
        f"joint_fafmax_{f}": joint_fafmax_idx[f] for f in joint_fafmax_idx.keys()
    }
    expr_dict.update(joint_fafmax_dict)

    logger.info("Unfurling age hists...")
    hist_idx = ht.histograms.age_hists
    age_hists = ["age_hist_het", "age_hist_hom"]
    age_hist_dict = {
        f"{hist}_{f}": (
            hl.delimit(hist_idx[hist][f], delimiter="|")
            if "bin" in f
            else hist_idx[hist][f]
        )
        for hist in age_hists
        for f in hist_idx[hist].keys()
    }
    expr_dict.update(age_hist_dict)

    return hl.struct(**expr_dict), entries_to_remove


def make_info_expr(
    t: hl.Table,
    hist_prefix: str = "",
    data_type: str = "exomes",
) -> Dict[str, hl.expr.Expression]:
    """
    Make Hail expression for variant annotations to be included in VCF INFO field.

    :param t: Table containing variant annotations to be reformatted for VCF export.
    :param hist_prefix: Prefix to use for histograms.
    :param data_type: Data type to make info expression for. One of "exomes" or
        "genomes". Default is "exomes".
    :return: Dictionary containing Hail expressions for relevant INFO annotations.
    :rtype: Dict[str, hl.expr.Expression]
    """
    vcf_info_dict = {}
    # Add site-level annotations and AS annotations to vcf_info_dict
    for field in SITE_FIELDS[data_type] + AS_FIELDS:
        vcf_info_dict[field] = t["release_ht_info"][f"{field}"]

    for field in AS_VQSR_FIELDS:
        vcf_info_dict[field] = t["vqsr_results"][f"{field}"]

    # Add region_flag and allele_info fields to info dict
    for field in ALLELE_TYPE_FIELDS[data_type]:
        vcf_info_dict[field] = t["allele_info"][f"{field}"]
    for field in REGION_FLAG_FIELDS[data_type]:
        vcf_info_dict[field] = t["region_flag"][f"{field}"]

    # Add underscore to hist_prefix if it isn't empty
    if hist_prefix != "":
        hist_prefix += "_"

    # Histograms to export are:
    # gq_hist_alt, gq_hist_all, dp_hist_alt, dp_hist_all, ab_hist_alt
    # We previously dropped:
    # _n_smaller for all hists
    # _bin_edges for all hists
    # _n_larger for all hists EXCEPT DP hists
    for hist in HISTS:
        hist_dict = {
            f"{hist}_bin_freq": hl.delimit(
                t.histograms.qual_hists[hist].bin_freq, delimiter="|"
            ),
        }
        vcf_info_dict.update(hist_dict)

        if "dp" in hist:
            vcf_info_dict.update(
                {f"{hist}_n_larger": t.histograms.qual_hists[hist].n_larger},
            )

    # Add in silico annotations to info dict
    insilico_idx = t.in_silico_predictors
    for field in INSILICO_FIELDS:
        if field == "cadd":
            vcf_info_dict[f"{field}_raw_score"] = insilico_idx[field]["raw_score"]
            vcf_info_dict[f"{field}_phred"] = insilico_idx[field]["phred"]
        else:
            vcf_info_dict[field] = insilico_idx[field]

    # Add VRS annotations to info dict
    for field in VRS_FIELDS_DICT:
        vcf_info_dict[field] = t["release_ht_info"]["vrs"][field]

    # Add vep annotations to info dict
    vcf_info_dict["vep"] = t["vep"]

    return vcf_info_dict


def prepare_ht_for_validation(
    ht: hl.Table,
    data_type: str = "exomes",
    freq_entries_to_remove: Optional[List[str]] = None,
    vcf_info_reorder: Optional[List[str]] = None,
) -> hl.Table:
    """
    Prepare HT for validity checks and export.

    :param ht: Release Hail Table
    :param data_type: Data type to prepare HT for. One of "exomes" or "genomes".
        Default is "exomes".
    :param freq_entries_to_remove: List of entries to remove from freq
    :param vcf_info_reorder: Order of VCF INFO fields
    :return: Hail Table prepared for validity checks and export
    """
    logger.info(
        "Unfurling nested gnomAD frequency annotations and add to INFO field..."
    )
    info_struct, freq_entries_to_remove = unfurl_nested_annotations(
        ht, entries_to_remove=freq_entries_to_remove, data_type=data_type
    )

    logger.info("Constructing INFO field")
    ht = ht.annotate(
        region_flag=ht.region_flags,
        release_ht_info=ht.info,
        info=info_struct,
        rsid=hl.str(";").join(ht.rsid),
        vep=vep_struct_to_csq(
            ht.vep, has_polyphen_sift=False
        ),  # TODO: Remove polyhen and sift from vep csq header in vcf
    )
    # Add variant annotations to INFO field
    # This adds the following:
    #   region flag for problematic regions
    #   annotations in ht.release_ht_info (site and allele-specific annotations),
    #   info struct (unfurled data obtained above),
    #   dbSNP rsIDs
    #   all VEP annotations
    ht = ht.annotate(info=ht.info.annotate(**make_info_expr(ht, data_type=data_type)))

    if freq_entries_to_remove:
        ht = ht.annotate_globals(
            vep_csq_header=process_vep_csq_header(VEP_CSQ_HEADER),
            freq_entries_to_remove=freq_entries_to_remove,
        )
    else:
        ht = ht.annotate_globals(
            vep_csq_header=process_vep_csq_header(
                VEP_CSQ_HEADER
            ),  # TODO: Process this to drop polyphen and sift
            freq_entries_to_remove=hl.empty_set(hl.tstr),
        )

    # Select relevant fields for VCF export
    ht = ht.select("info", "filters", "rsid")

    if vcf_info_reorder:
        logger.info("Rearranging fields to desired order...")
        ht = ht.annotate(
            info=ht.info.select(*vcf_info_reorder, *ht.info.drop(*vcf_info_reorder))
        )

    return ht


def process_vep_csq_header(vep_csq_header: str = VEP_CSQ_HEADER) -> str:
    """
    Process VEP CSQ header string, delimited by |, to remove polyphen and sift annotations.

    :param vep_csq_header: VEP CSQ header.
    :return: Processed VEP CSQ header.
    """
    logger.info("Processing VEP CSQ header...")
    vep_csq_header = vep_csq_header.split("|")
    vep_csq_header = [f for f in vep_csq_header if f not in ["PolyPhen", "SIFT"]]
    return vep_csq_header


def check_globals_for_retired_terms(ht: hl.Table) -> None:
    """
    Check list of dictionaries to see if the keys in the dictionaries contain either 'pop and 'oth'.

    :param ht: Input Table
    """
    logger.info("Checking globals for retired terms...")
    errors = []

    for field in ht.globals:
        if field.endswith("meta"):
            for d in hl.eval(ht[field]):
                if "pop" in d.keys():
                    errors.append(
                        f"Found retired term 'pop' in global {field} annotation: {d}"
                    )
                if "oth" in d.values():
                    errors.append(
                        f"Found retired term 'oth' in global {field} annotation: {d}"
                    )
        if "index_dict" in field:
            for k in hl.eval(ht[field]).keys():
                if "oth" in k:
                    errors.append(
                        f"Found retired term 'oth' in global {field} annotation: {k}"
                    )

    if len(errors) > 0:
        logger.info("Failed retired term check")
        pprint(errors)
    else:
        logger.info("Passed retired term check: No retired terms found in globals.")


def main(args):  # noqa: D103
    hl.init(
        log="/validate_and_export_vcf.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    # SSA Logs are easier to troubleshoot with and this hits a
    # class too large error without this no whole stage codegen flag
    # after I added joint_fafmax -- just the straw that broke hails back
    hl._set_flags(use_ssa_logs="1", no_whole_stage_codegen="1")

    overwrite = args.overwrite
    test = args.test
    data_type = args.data_type
    resources = get_export_resources(overwrite, data_type, test)

    try:
        if args.validate_release_ht:
            logger.info("Running release HT validation...")
            res = resources.validate_release_ht
            # res.check_resource_existence() #TODO: Restore once PR is ready for merge
            ht = res.release_ht  # .ht()#TODO: Restore once PR is ready for merge

            # if test: #TODO: Restore once PR is ready for merge
            #     logger.info("Filtering to test partitions...")
            #     ht = filter_to_test(ht)

            logger.info(
                "Checking globals for retired terms and checking their associated row"
                " annotation lengths..."
            )
            check_globals_for_retired_terms(ht)
            pprint_global_anns(ht)
            check_global_and_row_annot_lengths(ht, LEN_COMP_GLOBAL_ROWS)

            logger.info("Preparing HT for validity checks and export...")
            ht = prepare_ht_for_validation(ht, data_type=data_type)
            # Note: Checkpoint saves time in validity checks and final export by not
            # needing to run the VCF HT prep on each chromosome -- more needs to happen
            # before ready for export, but this is an intermediate write.
            logger.info("Writing prepared VCF HT for validity checks and export...")
            ht = ht.checkpoint(res.validated_ht.path, overwrite=overwrite)

            validate_release_t(
                ht,
                subsets=SUBSETS[data_type],
                pops=POPS,
                site_gt_check_expr={
                    "monoallelic": ht.info.monoallelic,
                    "only_het": ht.info.only_het,
                },
                verbose=args.verbose,
                delimiter="_",
                sample_sum_sets_and_pops=SAMPLE_SUM_SETS_AND_POPS[data_type],
                variant_filter_field="AS_VQSR",
                problematic_regions=REGION_FLAG_FIELDS[data_type],
                single_filter_count=True,
            )

            ht.describe()

    finally:
        logger.info("Copying log to logging bucket...")
        hl.copy_log(get_logging_path("validity_checks_and_export"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validate-release-ht",
        help="Run release HT validation",
        action="store_true",
    )
    parser.add_argument(
        "--verbose",
        help="Log successes in addition to failures during validation",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data and start from raw inputs",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Use test resources",
        action="store_true",
    )
    parser.add_argument(
        "-d",
        "--data-type",
        help="Data type to run validity checks on.",
        default="exomes",
        choices=["exomes", "genomes"],
    )

    main(parser.parse_args())
