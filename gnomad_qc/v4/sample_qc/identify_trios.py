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
from gnomad.utils.vcf import SEXES
from hail.utils.misc import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds
from gnomad_qc.v4.resources.meta import meta
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


# TODO: change to make fake pedigree, save it, then run mendle errors.
def run_mendel_errors(fake_fam_prop: float = 0.01, test: bool = False) -> hl.Table:
    """
    Add description.

    :return:
    """
    meta_ht = meta.ht()
    ped = pedigree(finalized=False).pedigree()

    n_fake_trios = int(fake_fam_prop * len(ped.complete_trios()))
    logger.info("Generating fake pedigree with %i trios.", n_fake_trios)
    fake_ped = create_fake_pedigree(
        n=n_fake_trios,
        sample_list=list(meta_ht.filter(meta_ht.high_quality).s.collect()),
        real_pedigree=ped,
    )
    merged_ped = hl.Pedigree(trios=ped.trios + fake_ped.trios)
    ped_samples = [s for t in merged_ped.trios for s in [t.s, t.pat_id, t.mat_id]]

    vds = get_gnomad_v4_vds(split=True)
    if test:
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(range(5)),
            vds.variant_data._filter_partitions(range(5)),
        )
    else:
        vds = hl.vds.filter_chromosomes(
            vds, keep="chr20"
        )  # TODO: should we filter to only chr20 like we did in v3, or grab the whole thing since this is exomes (still probably just want autosomes)?

    vds = hl.vds.filter_samples(vds, ped_samples)
    mt = hl.vds.to_dense_mt(vds)
    mt = mt.select_entries("GT").checkpoint(
        new_temp_file("mendel_errors_chr20", extension="mt")
    )

    logger.info(f"Running Mendel errors for {len(ped.trios)} trios.")
    mendel_errors, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)

    return mendel_errors


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


def get_filtered_relatedness_ht() -> hl.Table:
    """
    Add description.

    :return:
    """
    ht = relatedness().ht()
    filter_ht = finalized_outlier_filtering().ht()
    ht = ht.filter((ht.i.data_type == "exomes") & (ht.j.data_type == "exomes"))
    ht = ht.key_by(i=ht.i.s, j=ht.j.s)

    # Remove all pairs with a QC-filtered sample
    ht = ht.filter(
        filter_ht[ht.i].outlier_filtered | filter_ht[ht.j].outlier_filtered,
        keep=False,
    )

    return ht


def main(args):
    """Add description."""
    hl.init(
        log="/identify_trios.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )

    if args.identify_duplicates:
        logger.info("Selecting best duplicate per duplicated sample set")
        ht = get_duplicated_samples_ht(
            get_duplicated_samples(get_filtered_relatedness_ht()),
            sample_rankings().ht(),
        )
        ht.write(duplicates().path, overwrite=args.overwrite)

    if args.infer_families:
        # TODO: For v2, project info was used in trio determination, do we want to do that also?
        #  https://github.com/broadinstitute/gnomad_qc/blob/0df3e1be3cf20cdda71c9f2d9a23089e7a4c79de/gnomad_qc/v2/sample_qc/create_fam.py#L157
        logger.info("Inferring families")
        sex_ht = sex.ht()
        # TODO: This is the easiest way to make this work with infer_families, but let
        #  me know if it would be better to change infer_families instead
        sex_ht = sex_ht.filter(hl.literal(SEXES).contains(sex_ht.sex_karyotype))
        sex_ht = sex_ht.annotate(is_female=sex_ht.sex_karyotype == "XX")
        ped = infer_families(
            get_filtered_relatedness_ht(),
            sex_ht,
            duplicates().ht(),
        )
        ped.write(pedigree(data_type=data_type, finalized=False).path)
        raw_trios = families_to_trios(ped)
        raw_trios.write(trios(data_type=data_type, finalized=False).path)

    if args.run_mendel_errors:
        mendel_errors = run_mendel_errors(pedigree(finalized=False).pedigree())
        mendel_errors.write(ped_mendel_errors().path, overwrite=args.overwrite)

    if args.finalize_ped:
        final_ped = filter_ped(
            pedigree(finalized=False).pedigree,
            ped_mendel_errors().ht(),
            args.max_dnm,
            args.max_mendel,
        )
        final_ped.write(pedigree().path)

        final_trios = families_to_trios(final_ped)
        final_trios.write(trios().path)


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

    select_dup_grp = parser.add_argument_group("Select duplicates")
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
