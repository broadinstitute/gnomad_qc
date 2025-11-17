"""Script to identify trios from relatedness data and filter based on Mendel errors and de novos."""

import argparse
import json
import logging
from collections import defaultdict
from typing import Dict, Optional, Tuple

import hail as hl
from gnomad.sample_qc.relatedness import (
    create_fake_pedigree,
    filter_to_trios,
    get_duplicated_samples,
    get_duplicated_samples_ht,
    infer_families,
)
from hail.utils.misc import new_temp_file

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v4.sample_qc.identify_trios import families_to_trios
from gnomad_qc.v5.resources.basics import (
    WORKSPACE_BUCKET,
    get_aou_vds,
    get_logging_path,
)
from gnomad_qc.v5.resources.meta import meta
from gnomad_qc.v5.resources.sample_qc import (
    dense_trios,
    duplicates,
    finalized_outlier_filtering,
    ped_filter_param_json_path,
    ped_mendel_errors,
    pedigree,
    relatedness,
    sample_rankings,
    trios,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("identify_trios")
logger.setLevel(logging.INFO)


def filter_relatedness_ht(ht: hl.Table, filter_ht: hl.Table) -> hl.Table:
    """
    Filter relatedness Table to only include pairs of samples that are both AoU genomes and not QC-filtered.

    :param ht: Relatedness Table.
    :param filter_ht: Outlier filtering Table.
    :return: Filtered relatedness Table.
    """
    ht = ht.filter(
        ((ht.i.data_type == "genomes") & (ht.i.project == "aou"))
        & ((ht.j.data_type == "genomes") & (ht.j.project == "aou"))
    )
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)

    # Remove all pairs with a QC-filtered sample
    ht = ht.filter(
        filter_ht[ht.i].outlier_filtered | filter_ht[ht.j].outlier_filtered,
        keep=False,
    )

    return ht


def run_create_fake_pedigree(
    ped: hl.Pedigree,
    filter_ht: hl.Table,
    fake_fam_prop: float = 0.1,
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
    logger.info("Generating fake Pedigree with %i trios...", n_fake_trios)
    fake_ped = create_fake_pedigree(
        n=n_fake_trios,
        sample_list=list(
            filter_ht.filter(filter_ht.outlier_filtered, keep=False).s.collect()
        ),
        real_pedigree=ped,
    )
    return fake_ped


def run_mendel_errors(
    vds: hl.vds.VariantDataset,
    ped: hl.Pedigree,
    fake_ped: hl.Pedigree,
) -> hl.Table:
    """
    Run Hail's `mendel_errors` on chr20 of the VDS subset to samples in `ped` and `fake_ped`.

    :param vds: Input VariantDataset.
    :param ped: Inferred Pedigree.
    :param fake_ped: Fake Pedigree.
    :return: Table with Mendel errors on chr20.
    """
    merged_ped = hl.Pedigree(trios=ped.trios + fake_ped.trios)
    ped_samples = [s for t in merged_ped.trios for s in [t.s, t.pat_id, t.mat_id]]

    logger.info("Sub-setting VDS to %i samples...", len(ped_samples))
    vds = hl.vds.filter_samples(vds, ped_samples)
    vds.variant_data = vds.variant_data.select_entries("LA", "LGT")
    mt = hl.vds.to_dense_mt(vds)

    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)

    logger.info("Running Mendel errors for %s trios.", len(merged_ped.trios))
    mendel_err_ht, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)

    return mendel_err_ht


def filter_ped(
    ped: hl.Pedigree,
    mendel_ht: hl.Table,
    max_mendel_z: Optional[int] = 3,
    max_de_novo_z: Optional[int] = 3,
    max_mendel_n: Optional[int] = None,
    max_de_novo_n: Optional[int] = None,
) -> Tuple[hl.Pedigree, Dict[str, Dict[str, int]]]:
    """
    Filter a Pedigree based on Mendel errors and de novo metrics.

    :param ped: Pedigree to filter.
    :param mendel_ht: Table with Mendel errors.
    :param max_mendel_z: Optional maximum z-score for Mendel error metrics. Default is 3.
    :param max_de_novo_z: Optional maximum z-score for de novo metrics. Default is 3.
    :param max_mendel_n: Optional maximum Mendel error count. Default is None.
    :param max_de_novo_n: Optional maximum de novo count. Default is None.
    :return: Tuple of filtered Pedigree and dictionary of filtering parameters.
    """
    cutoffs = {
        "mendel": (
            max_mendel_z,
            max_mendel_n,
        ),
        "de_novo": (
            max_de_novo_z,
            max_de_novo_n,
        ),
    }

    # Check that only one of `max_mendel_z` or `max_mendel_n` and one of `max_de_novo_z`
    # or `max_de_novo_n` are set. If both are set, favor using the max number over max
    # std dev.
    cutoffs_by_method = defaultdict(dict)
    for m, (max_z, max_n) in cutoffs.items():
        if max_n:
            if max_z:
                logger.warning(
                    "Both `max_%s_z` and `max_%s_n` are set. Using `max_%s_n` of %d!",
                    *(m,) * 3,
                    max_n,
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
    mendel_by_s = mendel_by_s.checkpoint(new_temp_file("filter_ped", extension="ht"))

    # Get aggregate stats (need mean and stdev) for each metric with std dev
    # cutoffs set.
    z_stats_expr = {}
    for m in cutoffs_by_method["stdev"]:
        stats_expr = hl.agg.stats(mendel_by_s[f"n_{m}"])
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
        max_n = z_stats[m].mean + max_z * z_stats[m].stdev
        logger.info(log_expr, "", max_n, m, max_z)
        filter_expr &= mendel_by_s[f"n_{m}"] < max_n
        cutoffs[m]["max_n"] = max_n

    for m, max_n in cutoffs_by_method["count"].items():
        log_expr = "Filtering %strios with more than %d %s errors."
        logger.info(log_expr, "", max_n, m)
        filter_expr &= mendel_by_s[f"n_{m}"] < max_n
        cutoffs[m] = {"max_n": max_n}

    # Filter inferred Pedigree to only trios that pass filters defined by `filter_expr`.
    trios = mendel_by_s.aggregate(
        hl.agg.filter(filter_expr, hl.agg.collect(mendel_by_s.s))
    )
    logger.info("Found %i trios passing filters.", len(trios))
    return hl.Pedigree([trio for trio in ped.trios if trio.s in trios]), cutoffs


def create_dense_trio_mt(
    fam_ht: hl.Table,
    meta_ht: hl.Table,
    test: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
) -> hl.MatrixTable:
    """
    Create a dense MatrixTable for high quality trios.

    :param fam_ht: Table with family information.
    :param meta_ht: Table with metadata information.
    :param test: Whether to filter to chr20 for testing. Default is False.
    :param naive_coalesce_partitions: Optional Number of partitions to coalesce the VDS
        to. Default is None.
    :return: Dense MatrixTable with high quality trios.
    """
    # Filter the metadata table to only high quality AoU samples.
    meta_ht = meta_ht.filter(
        meta_ht.high_quality & (meta_ht.project_meta.project == "aou")
    )
    fam_ht = fam_ht.filter(
        hl.is_defined(meta_ht[fam_ht.id])
        & hl.is_defined(meta_ht[fam_ht.pat_id])
        & hl.is_defined(meta_ht[fam_ht.mat_id])
    )
    meta_ht = filter_to_trios(meta_ht, fam_ht)

    # Get the gnomAD VDS filtered to high quality releasable trios.
    # Using 'entries_to_keep' to keep all entries that are not `gvcf_info` because it
    # is likely not needed, and removal will reduce the size of the dense MatrixTable.
    vds = get_aou_vds(
        filter_samples=meta_ht,
        chrom="chr20" if test else None,
        entries_to_keep=["LA", "LGT", "LAD", "LPGT", "LPL", "DP", "GQ", "SB"],
        naive_coalesce_partitions=naive_coalesce_partitions,
    )
    return hl.vds.to_dense_mt(vds)


def main(args):
    """Identify trios and filter based on Mendel errors and de novos."""
    if args.rwb:
        environment = "rwb"
        hl.init(
            log="/home/jupyter/workspaces/gnomadproduction/identify_trios.log",
            tmp_dir=f"gs://{WORKSPACE_BUCKET}/tmp/4_day",
        )
    else:
        environment = "batch"
        hl.init(
            backend="batch",
            log="/identify_trios.log",
            tmp_dir="gs://gnomad-tmp-4day",
            gcs_requester_pays_configuration="broad-mpg-gnomad",
            default_reference="GRCh38",
            regions=["us-central1"],
            # TODO: Add machine configurations for Batch.
        )
    hl.default_reference("GRCh38")

    overwrite = args.overwrite
    test = args.test

    try:

        rel_ht = relatedness().ht()
        filter_ht = finalized_outlier_filtering().ht()
        rel_ht = filter_relatedness_ht(rel_ht, filter_ht)
        dup_ht_path = duplicates().path
        raw_ped_path = pedigree(finalized=False).path
        fake_ped_path = pedigree(finalized=False, fake=True).path
        mendel_err_ht_path = ped_mendel_errors(test=test).path
        final_ped_path = pedigree(test=test).path
        filter_json_path = ped_filter_param_json_path(test=test)
        trios_path = trios(test=test).path
        dense_trio_mt_path = dense_trios(test=test).path

        if args.identify_duplicates:
            logger.info("Selecting best duplicate per duplicated sample set...")
            check_resource_existence(
                output_step_resources={"duplicates_ht": [dup_ht_path]},
                overwrite=overwrite,
            )
            ht = get_duplicated_samples_ht(
                get_duplicated_samples(rel_ht),
                sample_rankings(release=True).ht(),
            )
            ht.write(dup_ht_path, overwrite=overwrite)

        if args.infer_families:
            logger.info("Inferring families...")
            check_resource_existence(
                output_step_resources={"raw_pedigree": [raw_ped_path]},
                overwrite=overwrite,
            )
            sex_ht = meta().ht()
            # Filter to AoU release genomes;
            # all of these samples are XX or XY.
            sex_ht = sex_ht.filter(
                (sex_ht.project_meta.project == "aou") & sex_ht.release
            )
            # `hl.Trio` requires a boolean column `is_female`.
            sex_ht = sex_ht.annotate(is_female=sex_ht.sex_karyotype == "XX")
            ped = infer_families(rel_ht, sex_ht, hl.read_table(dup_ht_path))
            ped.write(raw_ped_path)

        if args.create_fake_pedigree:
            logger.info("Creating fake Pedigree...")
            check_resource_existence(
                output_step_resources={"fake_pedigree": [fake_ped_path]},
                overwrite=overwrite,
            )
            fake_ped = run_create_fake_pedigree(
                ped, filter_ht, fake_fam_prop=args.fake_fam_prop
            )
            fake_ped.write(fake_ped_path)

        if args.run_mendel_errors:
            logger.info("Running Mendel errors on chr20...")
            check_resource_existence(
                output_step_resources={"mendel_err_ht": [mendel_err_ht_path]},
                overwrite=overwrite,
            )
            vds = get_aou_vds(
                split=False,
                remove_dead_alleles=False,
                filter_partitions=range(5) if test else None,
                chrom="chr20",
            )

            mendel_err_ht = run_mendel_errors(
                vds,
                hl.Pedigree.read(raw_ped_path),
                hl.Pedigree.read(fake_ped_path),
            )
            mendel_err_ht.write(mendel_err_ht_path, overwrite=overwrite)

        if args.finalize_ped:
            logger.info("Finalizing Pedigree...")
            check_resource_existence(
                output_step_resources={
                    "final_pedigree": [final_ped_path],
                    "final_trios": [trios_path],
                    "filter_json": [filter_json_path],
                },
                overwrite=overwrite,
            )
            ped, filters = filter_ped(
                hl.Pedigree.read(raw_ped_path),
                hl.read_table(mendel_err_ht_path),
                max_mendel_z=args.max_mendel_z,
                max_de_novo_z=args.max_de_novo_z,
                max_mendel_n=args.max_mendel,
                max_de_novo_n=args.max_de_novo,
            )
            ped.write(final_ped_path)
            families_to_trios(ped, args.seed).write(trios_path)

            # Pedigree has no globals like a HT so write the parameters to a JSON file.
            logger.info(
                "Writing finalized pedigree filter dictionary to %s...",
                filter_json_path,
            )
            with hl.hadoop_open(filter_json_path, "w") as d:
                d.write(json.dumps(filters))

        if args.create_dense_trio_mt:
            logger.info("Creating dense trio MT...")
            check_resource_existence(
                output_step_resources={"dense_trio_mt": [dense_trio_mt_path]},
                overwrite=overwrite,
            )
            dense_trio_mt = create_dense_trio_mt(
                hl.Pedigree.read(final_ped_path),
                meta().ht(),
                test=test,
                naive_coalesce_partitions=args.naive_coalesce_partitions,
            )
            dense_trio_mt.write(dense_trio_mt_path, overwrite=overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("identify_trios", environment=environment))


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rwb",
        help="Run the script in RWB environment.",
        action="store_true",
    )
    parser.add_argument(
        "--batch",
        help="Run the script in Batch environment.",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite existing data.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Run mendel errors on only five partitions of the MT.",
        action="store_true",
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
        "--max-de-novo",
        help=(
            "Maximum number of raw de novo mutations for real trios. If specified and"
            " --ukb-max-de-novo is not, --max-de-novo will be used for all samples. If"
            " both --max-de-novo and --ukb-max-de-novo are specified, --max-de-novo"
            " will be used for non-UKB samples and --ukb-max-de-novo will be used for"
            " UKB samples."
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

    dense_trio_mt_args = parser.add_argument_group("Create dense trio MT.")
    dense_trio_mt_args.add_argument(
        "--create-dense-trio-mt",
        help=("Create a dense MT for high quality trios."),
        action="store_true",
    )
    dense_trio_mt_args.add_argument(
        "--naive-coalesce-partitions",
        help=("Number of partitions to coalesce the VDS to."),
        type=int,
        default=5000,
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
