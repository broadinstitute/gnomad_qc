"""Add description."""
import argparse
import logging
from collections import Counter, defaultdict
from typing import Optional

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


def families_to_trios(ped: hl.Pedigree) -> hl.Pedigree:
    """
    Add description.

    :param ped:
    :return:
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

    return hl.Pedigree(trios=[t[0] for t in trios_per_fam.values()])


def filter_relatedness_ht(ht: hl.Table, filter_ht: hl.Table) -> hl.Table:
    """
    Add description.

    :return:
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
    ped: hl.Pedigree, filter_ht: hl.Table, fake_fam_prop: float = 0.01
) -> hl.Pedigree:
    """
    Generate a fake pedigree with `fake_fam_prop` of the number of trios in `ped`.

    :return:
    """
    n_fake_trios = int(fake_fam_prop * len(ped.complete_trios()))
    logger.info("Generating fake pedigree with %i trios.", n_fake_trios)
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
    Add description.

    :return:
    """
    merged_ped = hl.Pedigree(trios=ped.trios + fake_ped.trios)
    ped_samples = [s for t in merged_ped.trios for s in [t.s, t.pat_id, t.mat_id]]

    vds = get_gnomad_v4_vds(split=True)
    if test:
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(range(5)),
            vds.variant_data._filter_partitions(range(5)),
        )
    else:
        # TODO: should we filter to only chr20 like we did in v3, or grab the whole
        #  thing since this is exomes (still probably just want autosomes)?
        vds = hl.vds.filter_chromosomes(vds, keep="chr20")

    vds = hl.vds.filter_samples(vds, ped_samples)
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.select_entries("GT").checkpoint(
        new_temp_file("mendel_errors_chr20", extension="mt")
    )

    logger.info(f"Running Mendel errors for {len(ped.trios)} trios.")
    # TODO: save other metrics or only sample mendel errors?
    mendel_errors, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)

    return mendel_errors


def filter_ped(
    raw_ped: hl.Pedigree,
    mendel: hl.Table,
    max_mendel_z: Optional[int] = 3,
    max_de_novo_z: Optional[int] = 3,
    max_mendel: Optional[int] = None,
    max_de_novo: Optional[int] = None,
) -> hl.Pedigree:
    """
    Add description.

    :param raw_ped:
    :param mendel:
    :param max_de_novo_z:
    :param max_mendel_z:
    :param max_de_novo:
    :param max_mendel:
    :return:
    """
    z_metrics = []
    max_metrics = []
    for m in ["mendel", "de_novo"]:
        max_m = locals()[f"max_{m}"]
        if locals()[f"max_{m}_z"]:
            if max_m:
                logger.warning(
                    "Both `max_%s_z` and `max_%s` are set. Using `max_%s` of %i!",
                    *(m,) * 3,
                    max_m,
                )
                max_metrics.append(m)
            else:
                z_metrics.append(m)

    mendel = mendel.filter(mendel.fam_id.startswith("fake"))
    mendel_by_s = (
        mendel.group_by(mendel.s)
        .aggregate(
            fam_id=hl.agg.take(mendel.fam_id, 1)[0],
            n_mendel=hl.agg.count(),
            n_de_novo=hl.agg.count_where(
                mendel.mendel_code == 2
            ),  # Code 2 is parents are hom ref, child is het
        )
        .checkpoint(new_temp_file("filter_ped", extension="ht"))
    )

    z_stats_expr = {}
    for m in z_metrics:
        z_stats_expr[m] = hl.agg.stats(mendel_by_s[f"n_{m}"])

    filter_expr = hl.literal(True)
    z_stats = mendel_by_s.aggregate(**z_stats_expr)
    for m in ["mendel", "de_novo"]:
        if m in z_stats:
            max_m_z = locals()[f"max_{m}_z"]
            max_m = z_stats[m].mean + max_m_z * z_stats[m].stdev
            logger.info(
                "Filtering trios with more than %f %s errors (%i z-score)",
                max_m,
                m,
                max_m_z,
            )
        elif m in max_metrics:
            max_m = locals()[f"max_{m}"]
            logger.info("Filtering trios with more than %f %s errors", max_m, m)
        else:
            continue
        filter_expr &= mendel_by_s[f"n_{m}"] < max_m

    trios = mendel_by_s.aggregate(
        hl.agg.filter(filter_expr, hl.agg.collect(mendel_by_s.s))
    )
    logger.info("Found %i trios passing filters.", len(trios))
    return hl.Pedigree([trio for trio in raw_ped.trios if trio.s in trios])


def get_trio_resources(overwrite: bool, test: bool) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the trio identification pipeline.

    :return: PipelineResourceCollection containing resources for all steps of the
        trio identification pipeline.
    """
    # Adding resources from previous scripts that are used by multiple steps in the
    # trio identification pipeline.
    rel_ht = relatedness()
    filter_ht = finalized_outlier_filtering()
    filter_input = {
        "outlier_filtering.py --create-finalized-outlier-filter": {
            "filter_ht": filter_ht
        },
    }
    rel_input = {
        "relatedness.py --finalize-relatedness-ht": {"rel_ht": rel_ht},
    }

    # Initialize outlier filtering pipeline resource collection.
    trio_pipeline = PipelineResourceCollection(
        pipeline_name="identify_trios",
        pipeline_resources={"filter_ht": filter_ht, "rel_ht": rel_ht},
        overwrite=overwrite,
    )

    # Create resource collection for each step of the outlier filtering pipeline.
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
        output_resources={
            "raw_ped": pedigree(finalized=False),
            "raw_trios": trios(finalized=False),
        },
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
        },
        pipeline_input_steps=[infer_families, run_mendel_errors],
    )

    # Add all steps to the outlier filtering pipeline resource collection.
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
    """Add description."""
    hl.init(
        log="/identify_trios.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    overwrite = args.overwrite
    test = args.test
    trio_resources = get_trio_resources(overwrite, test)
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
        sex_ht = res.sex.ht()
        sex_ht = sex_ht.filter(hl.literal(SEXES).contains(sex_ht.sex_karyotype))
        sex_ht = sex_ht.annotate(is_female=sex_ht.sex_karyotype == "XX")
        ped = infer_families(rel_ht, sex_ht, res.dup_ht.ht())
        ped.write(res.raw_ped.path)
        # TODO: add options, remove all trios with more than one offspring in the
        #  family, or keep a random one from each family?
        families_to_trios(ped).write(res.raw_trios.path)

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
        ped = filter_ped(
            res.raw_ped.pedigree(),
            res.mendel_err_ht.ht(),
            args.max_mendel_z,
            args.max_de_novo_z,
            args.max_mendel,
            args.max_de_novo,
        )
        ped.write(res.final_ped.path)
        families_to_trios(ped).write(res.final_trios.path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test",
        help="Runs mendel errors on five partitions of the MT.",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    identify_dup = parser.add_argument_group("Duplicate identification")
    identify_dup.add_argument(
        "--identify-duplicates",
        help=(
            "Create a table with duplicate samples indicating which one is the best to"
            " use based on the ranking of all samples after sample QC metric outlier "
            "filtering of (ranking used to determine related samples to drop for the "
            "release)."
        ),
        action="store_true",
    )

    inter_fam = parser.add_argument_group("Pedigree inference")
    inter_fam.add_argument(
        "--infer-families",
        help=(
            "Infer families and trios using the relatedness Table and duplicate Table."
        ),
        action="store_true",
    )

    fake_ped = parser.add_argument_group("Fake pedigree creation")
    fake_ped.add_argument(
        "--create-fake-pedigree",
        help=(
            "Create a fake pedigree from unrelated samples in the data for comparison "
            "to the inferred pedigree."
        ),
        action="store_true",
    )
    fake_ped.add_argument(
        "--fake-fam-prop",
        help=(
            "Number of fake trios to generate as a proportion of the total number of"
            " trios found in the data. Default is 0.1."
        ),
        default=0.1,
        type=float,
    )

    mendel_err = parser.add_argument_group("Mendel error calculation")
    mendel_err.add_argument(
        "--run-mendel-errors",
        help="Calculate mendel errors for the inferred and fake pedigrees.",
        action="store_true",
    )
    finalize_ped = parser.add_argument_group(
        "Pedigree filtering for final pedigree generation"
    )
    finalize_ped.add_argument(
        "--finalize-ped",
        help=(
            "Create final families/trios ped files by excluding trios where the number "
            "of Mendel errors or de novos are outliers. Outliers can be defined as "
            "trios that have Mendel errors or de novos higher than those specified in "
            "--max_mendel and --max_de_novo respectively. They can also be defined as "
            "trios that have Mendel errors or de novos higher than --max-mendel-z or "
            "--max-de_novo-z standard deviations above the mean across inferred trios."
        ),
        action="store_true",
    )
    finalize_ped.add_argument(
        "--max-mendel-z",
        help=(
            "Max number of standard deviations above the mean Mendel errors across "
            "inferred trios to keep a trio. Default is 3."
        ),
        default=3,
        type=int,
    )
    finalize_ped.add_argument(
        "--max-de-novo-z",
        help=(
            "Max number of standard deviations above the mean de novos across inferred "
            "trios to keep a trio. Default is 3."
        ),
        default=3,
        type=int,
    )
    finalize_ped.add_argument(
        "--max-mendel",
        help="Maximum number of raw Mendel errors for real trios.",
        type=int,
    )
    finalize_ped.add_argument(
        "--max-de-novo",
        help="Maximum number of raw de novo mutations for real trios.",
        type=int,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
