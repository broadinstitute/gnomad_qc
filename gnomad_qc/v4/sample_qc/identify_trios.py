"""Script to identify trios from relatedness data and filter based on Mendel errors and de novos."""
import argparse
import json
import logging
import random
from collections import Counter, defaultdict
from typing import Dict, Optional, Tuple

import hail as hl
from gnomad.sample_qc.relatedness import (
    create_fake_pedigree,
    get_duplicated_samples,
    get_duplicated_samples_ht,
    infer_families,
)
from gnomad.utils.slack import slack_notifications
from gnomad.utils.vcf import SEXES
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.sample_qc import (
    duplicates,
    finalized_outlier_filtering,
    ped_filter_param_json_path,
    ped_mendel_errors,
    pedigree,
    relatedness,
    sample_rankings,
    sex,
    trios,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("identify_trios")
logger.setLevel(logging.INFO)


def families_to_trios(ped: hl.Pedigree, seed: int = 24) -> hl.Pedigree:
    """
    Convert a Pedigree with families to a Pedigree with only one random trio per family.

    :param ped: Pedigree with families.
    :param seed: Random seed for choosing trio to keep from each family. Default is 24.
    :return: Pedigree with only one trio per family.
    """
    trios_per_fam = defaultdict(list)
    for trio in ped.trios:
        trios_per_fam[trio.fam_id].append(trio)

    message = [
        f"Found a total of {len(ped.trios)} trios in {len(trios_per_fam)} families:"
    ]
    for n_trios, n_fam in Counter([len(t) for t in trios_per_fam.values()]).items():
        message.append(f"{n_fam} with {n_trios} trios.")
    logger.info("\n".join(message))

    random.seed(seed)
    return hl.Pedigree(trios=[random.choice(t) for t in trios_per_fam.values()])


def filter_relatedness_ht(ht: hl.Table, filter_ht: hl.Table) -> hl.Table:
    """
    Filter relatedness Table to only include pairs of samples that are both exomes and not QC-filtered.

    :param ht: Relatedness Table.
    :param filter_ht: Outlier filtering Table.
    :return: Filtered relatedness Table.
    """
    ht = ht.filter((ht.i.data_type == "exomes") & (ht.j.data_type == "exomes"))
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)

    # Remove all pairs with a QC-filtered sample
    ht = ht.filter(
        filter_ht[ht.i].outlier_filtered | filter_ht[ht.j].outlier_filtered,
        keep=False,
    )

    return ht


def run_create_fake_pedigree(
    ped: hl.Pedigree, filter_ht: hl.Table, fake_fam_prop: float = 0.1
) -> hl.Pedigree:
    """
    Generate a fake Pedigree with `fake_fam_prop` defining the proportion of the number of trios in `ped` to use.

    :param ped: Pedigree to use for generating fake Pedigree.
    :param filter_ht: Outlier filtering Table.
    :param fake_fam_prop: Proportion of trios in `ped` to use for generating fake
        Pedigree. Default is 0.1.
    :return: Fake Pedigree.
    """
    n_fake_trios = int(fake_fam_prop * len(ped.complete_trios()))
    logger.info("Generating fake Pedigree with %i trios.", n_fake_trios)
    fake_ped = create_fake_pedigree(
        n=n_fake_trios,
        sample_list=list(filter_ht.filter(filter_ht.outlier_filtered).s.collect()),
        real_pedigree=ped,
    )

    return fake_ped


def run_mendel_errors(
    ped: hl.Pedigree, fake_ped: hl.Pedigree, test: bool = False
) -> hl.Table:
    """
    Run Hail's `mendel_errors` on chr20 of the VDS subset to samples in `ped` and `fake_ped`.

    :param ped: Inferred Pedigree.
    :param fake_ped: Fake Pedigree.
    :param test: Whether to run on five partitions of the VDS for testing. Default is
        False.
    :return: Table with Mendel errors on chr20.
    """
    merged_ped = hl.Pedigree(trios=ped.trios + fake_ped.trios)
    ped_samples = [s for t in merged_ped.trios for s in [t.s, t.pat_id, t.mat_id]]

    logger.info(f"Sub-setting VDS to %i samples.", len(ped_samples))
    # Partition VDS on read to a larger number of partitions to speed up computation.
    vds = get_gnomad_v4_vds(split=True, n_partitions=150000, remove_dead_alleles=False)
    if test:
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(range(5)),
            vds.variant_data._filter_partitions(range(5)),
        )
    else:
        vds = hl.vds.filter_chromosomes(vds, keep="chr20")

    vds = hl.vds.filter_samples(vds, ped_samples)
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.select_entries("GT")

    logger.info(f"Running Mendel errors for %s trios.", len(merged_ped.trios))
    mendel_err_ht, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)

    return mendel_err_ht


def filter_ped(
    ped: hl.Pedigree,
    mendel_ht: hl.Table,
    max_mendel_z: Optional[int] = 3,
    max_de_novo_z: Optional[int] = 3,
    stratify_ukb: bool = False,
    max_mendel_n: Optional[int] = None,
    ukb_max_mendel_n: Optional[int] = None,
    max_de_novo_n: Optional[int] = None,
    ukb_max_de_novo_n: Optional[int] = None,
) -> Tuple[hl.Pedigree, Dict[str, Dict[str, int]]]:
    """
    Filter a Pedigree based on Mendel errors and de novo metrics.

    :param ped: Pedigree to filter.
    :param mendel_ht: Table with Mendel errors.
    :param max_mendel_z: Optional maximum z-score for Mendel error metrics. Default is 3.
    :param max_de_novo_z: Optional maximum z-score for de novo metrics. Default is 3.
    :param stratify_ukb: Whether to stratify z-score cutoffs by trio UK Biobank status.
        Default is False.
    :param max_mendel_n: Optional maximum Mendel error count. If specified and
        `ukb_max_mendel_n` is not, `max_mendel_n` will be used for all samples. If both
        `max_mendel_n` and `ukb_max_mendel_n` are specified, max_mendel_n` will be used
        for non-UKB samples and `ukb_max_mendel_n` will be used for UKB samples.
        Default is None.
    :param ukb_max_mendel_n: Optional maximum Mendel error count for trios in UK
        Biobank. If specified, `max_mendel_n` must also be specified for non-UKB
        samples. If not specified, but `max_mendel_n` is, `max_mendel_n` will be used
        for all samples. Default is None.
    :param max_de_novo_n: Optional maximum de novo count. If specified and
        `ukb_max_de_novo_n` is not, `max_de_novo_n` will be used for all samples. If
        both `max_de_novo_n` and `ukb_max_de_novo_n` are specified, max_de_novo_n` will
        be used for non-UKB samples and `ukb_max_de_novo_n` will be used for UKB
        samples. Default is None.
    :param ukb_max_de_novo_n: Optional maximum de novo count for trios in UK Biobank.
        If specified, `max_de_novo_n` must also be specified for non-UKB samples. If not
        specified, but `max_de_novo_n` is, `max_de_novo_n` will be used for all samples.
        Default is None.
    :return: Tuple of filtered Pedigree and dictionary of filtering parameters.
    """
    if ukb_max_mendel_n and not max_mendel_n:
        raise ValueError("`max_mendel_n` must be specified if `ukb_max_mendel_n` is.")
    if ukb_max_de_novo_n and not max_de_novo_n:
        raise ValueError("`max_de_novo_n` must be specified if `ukb_max_de_novo_n` is.")

    cutoffs = {
        "mendel": (
            max_mendel_z,
            {"UKB": ukb_max_mendel_n, "non-UKB": max_mendel_n}
            if ukb_max_mendel_n
            else max_mendel_n,
        ),
        "de_novo": (
            max_de_novo_z,
            {"UKB": ukb_max_de_novo_n, "non-UKB": max_de_novo_n}
            if ukb_max_de_novo_n
            else max_de_novo_n,
        ),
    }

    # Check that only one of `max_mendel_z` or `max_mendel_n` and one of `max_de_novo_z`
    # or `max_de_novo_n` are set. If both are set favor using the max number over max
    # std dev.
    cutoffs_by_method = defaultdict(dict)
    for m, (max_z, max_n) in cutoffs.items():
        if max_n:
            if max_z:
                logger.warning(
                    "Both `max_%s_z` and `max_%s_n` are set. Using `max_%s_n` of %d%s!",
                    *(m,) * 3,
                    max_n["non-UKB"] if isinstance(max_n, dict) else max_n,
                    f" and `ukb_max_{m}_n` of {max_n['UKB']}"
                    if isinstance(max_n, dict)
                    else "",
                )
            cutoffs_by_method["count"][m] = max_n
        elif max_z:
            cutoffs_by_method["stdev"][m] = max_z

    # Filter Mendel errors Table to only errors in inferred families not from the
    # fake Pedigree.
    mendel_ht = mendel_ht.filter(mendel_ht.fam_id.startswith("fake"), keep=False)

    # Aggregate Mendel errors Table by sample to get per sample mendel error and
    # de novo counts.
    mendel_by_s = mendel_ht.group_by(mendel_ht.s, mendel_ht.fam_id).aggregate(
        n_mendel=hl.agg.count(),
        # Code 2 is parents are hom ref, child is het.
        n_de_novo=hl.agg.count_where(mendel_ht.mendel_code == 2),
    )

    mendel_by_s = mendel_by_s.annotate(ukb=mendel_by_s.s.startswith("UKB"))
    mendel_by_s = mendel_by_s.checkpoint(new_temp_file("filter_ped", extension="ht"))

    # Get aggregate stats (need mean and stdev) for each metric with std dev
    # cutoffs set.
    z_stats_expr = {}
    for m in cutoffs_by_method["stdev"]:
        stats_expr = hl.agg.stats(mendel_by_s[f"n_{m}"])
        # If stratifying by UKB status, get stats for UKB and non-UKB separately.
        if stratify_ukb:
            stats_expr = hl.agg.group_by(mendel_by_s.ukb, stats_expr)
        z_stats_expr[m] = stats_expr

    z_stats = mendel_by_s.aggregate(hl.struct(**z_stats_expr))
    # Build filter expression to filter metrics by requested metrics and methods.
    filter_expr = hl.literal(True)
    cutoffs = {}
    for m, max_z in cutoffs_by_method["stdev"].items():
        cutoffs[m] = {"max_z": max_z}
        log_expr = (
            "Filtering %strios with more than %f %s errors (%i standard deviations "
            "from the mean)"
        )
        # If stratifying by UKB status, get cutoffs for UKB and non-UKB separately.
        if stratify_ukb:
            cutoffs[m]["max_n"] = {}
            ukb_filter_expr = hl.literal(False)
            for ukb in [True, False]:
                max_n = z_stats[m][ukb].mean + max_z * z_stats[m][ukb].stdev
                logger.info(
                    log_expr,
                    "UKB " if ukb else "non-UKB ",
                    max_n,
                    m,
                    max_z,
                )
                ukb_filter_expr |= (mendel_by_s.ukb == ukb) & (
                    mendel_by_s[f"n_{m}"] < max_n
                )
                cutoffs[m]["max_n"][f"{'' if ukb else 'non-'}UKB"] = max_n
            filter_expr &= ukb_filter_expr
        else:
            max_n = z_stats[m].mean + max_z * z_stats[m].stdev
            logger.info(log_expr, "", max_n, m, max_z)
            filter_expr &= mendel_by_s[f"n_{m}"] < max_n
            cutoffs[m]["max_n"] = max_n

    for m, max_n in cutoffs_by_method["count"].items():
        log_expr = "Filtering %strios with more than %d %s errors."
        # If stratifying by UKB status, use different cutoffs for UKB and non-UKB.
        if isinstance(max_n, dict):
            ukb_filter_expr = hl.literal(False)
            for ukb in [True, False]:
                ukb_max_n = max_n["UKB" if ukb else "non-UKB"]
                logger.info(log_expr, "UKB " if ukb else "non-UKB ", ukb_max_n, m)
                ukb_filter_expr |= (mendel_by_s.ukb == ukb) & (
                    mendel_by_s[f"n_{m}"] < ukb_max_n
                )
            filter_expr &= ukb_filter_expr
        else:
            logger.info(log_expr, "", max_n, m)
            filter_expr &= mendel_by_s[f"n_{m}"] < max_n
        cutoffs[m] = {"max_n": max_n}

    # Filter inferred Pedigree to only trios that pass filters defined by `filter_expr`.
    trios = mendel_by_s.aggregate(
        hl.agg.filter(filter_expr, hl.agg.collect(mendel_by_s.s))
    )
    logger.info("Found %i trios passing filters.", len(trios))
    return hl.Pedigree([trio for trio in ped.trios if trio.s in trios]), cutoffs


def get_trio_resources(overwrite: bool, test: bool) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the trio identification pipeline.

    :param overwrite: Whether to overwrite existing resources.
    :param test: Whether to use test resources.
    :return: PipelineResourceCollection containing resources for all steps of the
        trio identification pipeline.
    """
    # Adding resources from previous scripts that are used by multiple steps in the
    # trio identification pipeline.
    rel_ht = relatedness()
    filter_ht = finalized_outlier_filtering()
    pipeline_resources = {
        "outlier_filtering.py --create-finalized-outlier-filter": {
            "filter_ht": filter_ht
        },
        "relatedness.py --finalize-relatedness-ht": {"rel_ht": rel_ht},
    }

    # Initialize trio identification pipeline resource collection.
    trio_pipeline = PipelineResourceCollection(
        pipeline_name="identify_trios",
        pipeline_resources=pipeline_resources,
        overwrite=overwrite,
    )

    # Create resource collection for each step of the trio identification pipeline.
    identify_duplicates = PipelineStepResourceCollection(
        "--identify-duplicates",
        output_resources={"dup_ht": duplicates()},
        input_resources={
            "relatedness.py --compute-related-samples-to-drop": {
                "rank_ht": sample_rankings()
            },
        },
    )
    infer_families = PipelineStepResourceCollection(
        "--infer-families",
        output_resources={"raw_ped": pedigree(finalized=False)},
        pipeline_input_steps=[identify_duplicates],
        add_input_resources={
            "sex_inference.py --annotate-sex-karyotype": {"sex_ht": sex},
        },
    )
    create_fake_pedigree = PipelineStepResourceCollection(
        "--create-fake-pedigree",
        output_resources={"fake_ped": pedigree(finalized=False, fake=True)},
        pipeline_input_steps=[infer_families],
    )
    run_mendel_errors = PipelineStepResourceCollection(
        "--run-mendel-errors",
        output_resources={"mendel_err_ht": ped_mendel_errors(test)},
        pipeline_input_steps=[infer_families, create_fake_pedigree],
    )
    finalize_ped = PipelineStepResourceCollection(
        "--finalize-ped",
        output_resources={
            "final_ped": pedigree(test=test),
            "final_trios": trios(test=test),
            "filter_json": ped_filter_param_json_path(test=test),
        },
        pipeline_input_steps=[infer_families, run_mendel_errors],
    )

    # Add all steps to the trio identification pipeline resource collection.
    trio_pipeline.add_steps(
        {
            "identify_duplicates": identify_duplicates,
            "infer_families": infer_families,
            "create_fake_pedigree": create_fake_pedigree,
            "run_mendel_errors": run_mendel_errors,
            "finalize_ped": finalize_ped,
        }
    )

    return trio_pipeline


def main(args):
    """Identify trios and filter based on Mendel errors and de novos."""
    hl.init(
        log="/identify_trios.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    overwrite = args.overwrite
    test = args.test
    trio_resources = get_trio_resources(overwrite, test)
    trio_resources.check_resource_existence()
    filter_ht = trio_resources.filter_ht.ht()
    rel_ht = filter_relatedness_ht(trio_resources.rel_ht.ht(), filter_ht)

    if args.identify_duplicates:
        logger.info("Selecting best duplicate per duplicated sample set")
        res = trio_resources.identify_duplicates
        res.check_resource_existence()
        ht = get_duplicated_samples_ht(
            get_duplicated_samples(rel_ht),
            res.rank_ht.ht(),
        )
        ht.write(res.dup_ht.path, overwrite=args.overwrite)

    if args.infer_families:
        logger.info("Inferring families")
        res = trio_resources.infer_families
        res.check_resource_existence()
        sex_ht = res.sex_ht.ht()
        sex_ht = sex_ht.filter(hl.literal(SEXES).contains(sex_ht.sex_karyotype))
        sex_ht = sex_ht.annotate(is_female=sex_ht.sex_karyotype == "XX")
        ped = infer_families(rel_ht, sex_ht, res.dup_ht.ht())
        ped.write(res.raw_ped.path)

    if args.create_fake_pedigree:
        res = trio_resources.create_fake_pedigree
        res.check_resource_existence()
        run_create_fake_pedigree(
            res.raw_ped.pedigree(),
            filter_ht,
            fake_fam_prop=args.fake_fam_prop,
        ).write(res.fake_ped.path)

    if args.run_mendel_errors:
        res = trio_resources.run_mendel_errors
        res.check_resource_existence()
        run_mendel_errors(
            res.raw_ped.pedigree(), res.fake_ped.pedigree(), test=test
        ).write(res.mendel_err_ht.path, overwrite=args.overwrite)

    if args.finalize_ped:
        res = trio_resources.finalize_ped
        res.check_resource_existence()
        ped, filters = filter_ped(
            res.raw_ped.pedigree(),
            res.mendel_err_ht.ht(),
            args.max_mendel_z,
            args.max_de_novo_z,
            args.stratify_ukb,
            args.max_mendel,
            args.ukb_max_mendel,
            args.max_de_novo,
            args.ukb_max_de_novo,
        )
        ped.write(res.final_ped.path)
        families_to_trios(ped, args.seed).write(res.final_trios.path)

        # Pedigree has no globals like a HT so write the parameters to a JSON file.
        logger.info(
            "Writing finalized pedigree filter dictionary to %s.", res.filter_json
        )
        with hl.hadoop_open(res.filter_json, "w") as d:
            d.write(json.dumps(filters))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Runs mendel errors on only five partitions of the MT.",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False).",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    identify_dup_args = parser.add_argument_group("Duplicate identification")
    identify_dup_args.add_argument(
        "--identify-duplicates",
        help=(
            "Create a table with duplicate samples indicating which one is the best to"
            " use based on the ranking of all samples after sample QC metric outlier "
            "filtering (ranking used to determine related samples to drop for the "
            "release)."
        ),
        action="store_true",
    )

    inter_fam_args = parser.add_argument_group("Pedigree inference")
    inter_fam_args.add_argument(
        "--infer-families",
        help=(
            "Infer families and trios using the relatedness Table and duplicate Table."
        ),
        action="store_true",
    )

    fake_ped_args = parser.add_argument_group("Fake Pedigree creation")
    fake_ped_args.add_argument(
        "--create-fake-pedigree",
        help=(
            "Create a fake Pedigree from unrelated samples in the data for comparison "
            "to the inferred Pedigree."
        ),
        action="store_true",
    )
    fake_ped_args.add_argument(
        "--fake-fam-prop",
        help=(
            "Number of fake trios to generate as a proportion of the total number of"
            " trios found in the data. Default is 0.1."
        ),
        default=0.1,
        type=float,
    )

    mendel_err_args = parser.add_argument_group("Mendel error calculation")
    mendel_err_args.add_argument(
        "--run-mendel-errors",
        help="Calculate mendel errors for the inferred and fake Pedigrees on chr20.",
        action="store_true",
    )
    finalize_ped_args = parser.add_argument_group(
        "Pedigree filtering for final Pedigree generation"
    )
    finalize_ped_args.add_argument(
        "--finalize-ped",
        help=(
            "Create final families/trios ped files by excluding trios where the number "
            "of Mendel errors or de novos are outliers. Outliers can be defined as "
            "trios that have Mendel errors or de novos higher than those specified in "
            "--max-mendel and --max-de-novo respectively. They can also be defined as "
            "trios that have Mendel errors or de novos higher than --max-mendel-z or "
            "--max-de-novo-z standard deviations above the mean across inferred trios."
        ),
        action="store_true",
    )
    finalize_ped_args.add_argument(
        "--max-mendel-z",
        help=(
            "Max number of standard deviations above the mean Mendel errors across "
            "inferred trios to keep a trio. If flag is set, default is 3."
        ),
        nargs="?",
        const=3,
        type=int,
    )
    finalize_ped_args.add_argument(
        "--max-de-novo-z",
        help=(
            "Max number of standard deviations above the mean de novos across inferred "
            "trios to keep a trio. If flag is set, default is 3."
        ),
        nargs="?",
        const=3,
        type=int,
    )
    finalize_ped_args.add_argument(
        "--stratify-ukb",
        help=(
            "Stratify Mendel errors and de novo standard deviations cutoffs to UKB and "
            "non-UKB samples."
        ),
        action="store_true",
    )
    finalize_ped_args.add_argument(
        "--max-mendel",
        help=(
            "Maximum number of raw Mendel errors for real trios. If specified and "
            "--ukb-max-mendel is not, --max-mendel will be used for all samples. If "
            "both --max-mendel and --ukb-max-mendel are specified, --max-mendel will "
            "be used for non-UKB samples and --ukb-max-mendel will be used for UKB "
            "samples."
        ),
        type=int,
    )
    finalize_ped_args.add_argument(
        "--ukb-max-mendel",
        help=(
            "Maximum number of raw Mendel errors for real trios in UKB samples. If "
            "specified, --max-mendel must also be specified for non-UKB samples. If "
            "not specified, but --max-mendel is, --max-mendel will be used for all "
            "samples."
        ),
        type=int,
    )
    finalize_ped_args.add_argument(
        "--max-de-novo",
        help=(
            "Maximum number of raw de novo mutations for real trios. If specified and "
            "--ukb-max-mendel is not, --max-mendel will be used for all samples. If "
            "both --max-mendel and --ukb-max-mendel are specified, --max-mendel will "
            "be used for non-UKB samples and --ukb-max-mendel will be used for UKB "
            "samples."
        ),
        type=int,
    )
    finalize_ped_args.add_argument(
        "--ukb-max-de-novo",
        help=(
            "Maximum number of raw de novo mutations for real trios in UKB samples. If "
            "specified, --max-mendel must also be specified for non-UKB samples. If "
            "not specified, but --max-mendel is, --max-mendel will be used for all "
            "samples."
        ),
        type=int,
    )
    finalize_ped_args.add_argument(
        "--seed",
        help=(
            "Random seed for choosing one random trio per family to keep after "
            "filtering."
        ),
        type=int,
        default=24,
    )

    args = parser.parse_args()

    if args.ukb_max_mendel and not args.max_mendel:
        parser.error(
            "--ukb-max-mendel requires --max-mendel to be specified for non-UKB "
            "samples."
        )
    if args.ukb_max_de_novo and not args.max_de_novo:
        parser.error(
            "--ukb-max-de-novo requires --max-de-novo to be specified for non-UKB "
            "samples."
        )
    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
