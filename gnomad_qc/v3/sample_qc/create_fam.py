import argparse
import logging
from collections import Counter, defaultdict

import hail as hl
from gnomad.sample_qc.relatedness import (create_fake_pedigree,
                                          get_duplicated_samples,
                                          get_duplicated_samples_ht,
                                          infer_families)

from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.constants import CURRENT_RELEASE
from gnomad_qc.v3.resources.meta import (meta, ped_mendel_errors, pedigree,
                                         trios)
from gnomad_qc.v3.resources.sample_qc import (duplicates,
                                              get_relatedness_annotated_ht,
                                              release_samples_rankings, sex)

logger = logging.getLogger("create_fam")


def run_mendel_errors() -> hl.Table:
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
    mt = hl.filter_intervals(mt, [hl.parse_locus_interval("chr20", reference_genome='GRCh38')])
    mt = hl.experimental.sparse_split_multi(mt, filter_changed_loci=True)
    mt = mt.select_entries("GT", "END")
    mt = hl.experimental.densify(mt)
    mt = mt.filter_rows(hl.len(mt.alleles) == 2)
    mendel_errors, _, _, _ = hl.mendel_errors(mt["GT"], merged_ped)
    return mendel_errors


def run_infer_families() -> hl.Pedigree:
    logger.info("Inferring families")
    ped = infer_families(
        get_relatedness_annotated_ht(), sex.ht(), duplicates.ht()
    )

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
    mendel = mendel.filter(mendel.fam_id.startswith("fake"))
    mendel_by_s = (
        mendel.group_by(mendel.s)
            .aggregate(
            fam_id=hl.agg.take(mendel.fam_id, 1)[0],
            n_mendel=hl.agg.count(),
            n_de_novo=hl.agg.count_where(mendel.mendel_code == 2), # Code 2 is parents are hom ref, child is het
        )
            .persist()
    )

    good_trios = mendel_by_s.aggregate(
        hl.agg.filter(
            (mendel_by_s.n_mendel < max_mendel) & (mendel_by_s.n_de_novo < max_dnm),
            hl.agg.collect(mendel_by_s.s, ),
        )
    )
    logger.info(f"Found {len(good_trios)} trios passing filters")
    return hl.Pedigree([trio for trio in raw_ped.trios if trio.s in good_trios])


def main(args):
    hl.init(default_reference="GRCh38")

    if args.find_dups:
        logger.info("Selecting best duplicate per duplicated sample set")
        dups = get_duplicated_samples(get_relatedness_annotated_ht())
        dups_ht = get_duplicated_samples_ht(dups, release_samples_rankings.ht())
        dups_ht.write(duplicates.path, overwrite=args.overwrite)

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
        "--find_dups",
        help="Creates a table with duplicate samples indicating which one is the best to use.",
        action="store_true",
    )
    parser.add_argument("--infer_families", help="Infers families", action="store_true")
    parser.add_argument(
        "--run_mendel_errors", help="Runs mendel errors", action="store_true"
    )
    parser.add_argument(
        "--finalize_ped",
        help="Creates a final ped file by excluding families where the number of Mendel errors or de novos are higher than those specified in --max_dnm and --max_mendel",
        action="store_true",
    )
    parser.add_argument(
        "--max_dnm",
        help="Maximum number of raw de novo mutations for real trios",
        defaut=2200,
        type=int,
    )
    parser.add_argument(
        "--max_mendel",
        help="Maximum number of raw Mendel errors for real trios",
        defaut=3750,
        type=int,
    )

    main(parser.parse_args())
