"""Add description."""
import argparse
import logging
from collections import Counter, defaultdict

import hail as hl
from gnomad.sample_qc.relatedness import (
    create_fake_pedigree,
    get_duplicated_samples,
    get_duplicated_samples_ht,
    infer_families,
)
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v4.resources.constants import CURRENT_RELEASE
from gnomad_qc.v4.resources.sample_qc import (
    duplicates,
    finalized_outlier_filtering,
    joint_qc_meta,
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


def run_mendel_errors() -> hl.Table:
    """
    Add description.

    :return:
    """
    meta_ht = meta.ht()
    ped = pedigree.versions[f"{CURRENT_RELEASE}_raw"].pedigree()
    logger.info(f"Running Mendel errors for {len(ped.trios)} trios.")

    fake_ped = create_fake_pedigree(
        n=100,
        sample_list=list(
            meta_ht.aggregate(
                hl.agg.filter(
                    hl.rand_bool(0.01)
                    & (
                        (hl.len(meta_ht.qc_metrics_filters) == 0)
                        & hl.or_else(hl.len(meta_ht.hard_filters) == 0, False)
                    ),
                    hl.agg.collect_as_set(meta_ht.s),
                )
            )
        ),
        real_pedigree=ped,
    )
    merged_ped = hl.Pedigree(trios=ped.trios + fake_ped.trios)

    ped_samples = hl.literal(
        set(
            [s for trio in merged_ped.trios for s in [trio.s, trio.pat_id, trio.mat_id]]
        )
    )
    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)
    mt = mt.filter_cols(ped_samples.contains(mt.s))
    mt = hl.filter_intervals(
        mt, [hl.parse_locus_interval("chr20", reference_genome="GRCh38")]
    )
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = mt.select_entries("GT", "END")
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mendel_errors, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)
    return mendel_errors


def run_infer_families() -> hl.Pedigree:
    """
    Add description.

    :return:
    """
    logger.info("Inferring families")
    ped = infer_families(get_relatedness_annotated_ht(), sex.ht(), duplicates.ht())

    # Remove all trios containing any QC-filtered sample
    meta_ht = meta.ht()
    filtered_samples = meta_ht.aggregate(
        hl.agg.filter(
            (hl.len(meta_ht.qc_metrics_filters) > 0)
            | hl.or_else(hl.len(meta_ht.hard_filters) > 0, False),
            hl.agg.collect_as_set(meta_ht.s),
        )
    )

    return hl.Pedigree(
        trios=[
            trio
            for trio in ped.trios
            if trio.s not in filtered_samples
            and trio.pat_id not in filtered_samples
            and trio.mat_id not in filtered_samples
        ]
    )


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
    for n_trios, n_fam in Counter(
        [len(trios) for trios in trios_per_fam.values()]
    ).items():
        message.append(f"{n_fam} with {n_trios} trios.")
    logger.info("\n".join(message))

    return hl.Pedigree(trios=[trios[0] for trios in trios_per_fam.values()])


def filter_ped(
    raw_ped: hl.Pedigree, mendel: hl.Table, max_dnm: int, max_mendel: int
) -> hl.Pedigree:
    """
    Add description.

    :param raw_ped:
    :param mendel:
    :param max_dnm:
    :param max_mendel:
    :return:
    """
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
        .persist()
    )

    good_trios = mendel_by_s.aggregate(
        hl.agg.filter(
            (mendel_by_s.n_mendel < max_mendel) & (mendel_by_s.n_de_novo < max_dnm),
            hl.agg.collect(
                mendel_by_s.s,
            ),
        )
    )
    logger.info(f"Found {len(good_trios)} trios passing filters")
    return hl.Pedigree([trio for trio in raw_ped.trios if trio.s in good_trios])


def main(args):
    """Add description."""
    hl.init(
        log="/identify_trios.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    joint = args.joint
    data_type = "joint" if joint else "exomes"

    if args.identify_duplicates:
        logger.info("Selecting best duplicate per duplicated sample set")
        ht = relatedness().ht()
        # Remove all pairs with a QC-filtered sample
        filter_ht = finalized_outlier_filtering().ht()
        if not joint:
            ht = ht.filter((ht.i.data_type == "exomes") & (ht.j.data_type == "exomes"))
        else:
            joint_qc_meta_ht = joint_qc_meta.ht()
            filter_ht = joint_qc_meta_ht.annotate(
                outlier_filtered=filter_ht[joint_qc_meta_ht.key].outlier_filtered
            )
        ht = ht.key_by(i=ht.i.s, j=ht.j.s)
        # Remove all pairs with a QC-filtered sample
        ht = ht.filter(
            filter_ht[ht.i].outlier_filtered | filter_ht[ht.j].outlier_filtered,
            keep=False,
        )
        ht = get_duplicated_samples_ht(
            get_duplicated_samples(ht), sample_rankings().ht()
        )
        ht.write(duplicates(data_type=data_type).path, overwrite=args.overwrite)

    if args.infer_families:
        ped = run_infer_families()
        ped.write(pedigree.versions[f"{CURRENT_RELEASE}_raw"].path)
        raw_trios = families_to_trios(ped)
        raw_trios.write(trios.versions[f"{CURRENT_RELEASE}_raw"].path)

    if args.run_mendel_errors:
        mendel_errors = run_mendel_errors()
        mendel_errors.write(ped_mendel_errors.path, overwrite=args.overwrite)

    if args.finalize_ped:
        final_ped = filter_ped(
            pedigree.versions[f"{CURRENT_RELEASE}_raw"].pedigree,
            ped_mendel_errors.ht(),
            args.max_dnm,
            args.max_mendel,
        )
        final_ped.write(pedigree.path)

        final_trios = families_to_trios(final_ped)
        final_trios.write(trios.path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite all data from this subset (default: False)",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--joint",
        help=(
            "Whether to include both v3 genomes and v4 exomes. Default is only v4 "
            "exomes."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--identify-duplicates",
        help=(
            "Creates a table with duplicate samples indicating which one is the best to"
            " use."
        ),
        action="store_true",
    )
    parser.add_argument("--infer-families", help="Infers families", action="store_true")
    parser.add_argument(
        "--run-mendel-errors", help="Runs mendel errors", action="store_true"
    )
    parser.add_argument(
        "--finalize-ped",
        help=(
            "Creates a final ped file by excluding families where the number of Mendel"
            " errors or de novos are higher than those specified in --max_dnm and"
            " --max_mendel"
        ),
        action="store_true",
    )
    parser.add_argument(
        "--max-dnm",
        help="Maximum number of raw de novo mutations for real trios",
        default=2200,
        type=int,
    )
    parser.add_argument(
        "--max-mendel",
        help="Maximum number of raw Mendel errors for real trios",
        default=3750,
        type=int,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
