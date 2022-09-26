import argparse
import logging
from typing import Union

from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad.utils.annotations import annotate_adj, get_adj_expr
from gnomad.utils.slack import slack_notifications
import hail as hl

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.meta import meta as v3_meta
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_logging_path,
    get_gnomad_v4_vds,
)
from gnomad_qc.v4.resources.meta import project_meta as v4_meta
from gnomad_qc.v4.resources.sample_qc import (
    get_predetermined_qc,
    hard_filtered_samples,
    joint_qc_meta,
    predetermined_qc_sites,
    joint_qc,
    sample_chr20_mean_dp,
)
from gnomad_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("generate_qc_mt")
logger.setLevel(logging.INFO)


def create_filtered_dense_mt(
    mtds: Union[hl.vds.VariantDataset, hl.MatrixTable],
    split: bool = False,
) -> hl.MatrixTable:
    """
    Filter a sparse MatrixTable or VariantDataset to a set of predetermined QC sites and return a dense MatrixTable.

    :param mtds: Input MatrixTable or VariantDataset.
    :param split: Whether `mtds` should have multi-allelics split before filtering variants.
    :return: Filtered and densified MatrixTable.
    """
    is_vds = isinstance(mtds, hl.vds.VariantDataset)

    if is_vds:
        vds = mtds
    else:
        logger.info("Converting MatrixTable to VariantDataset...")
        mtds = mtds.select_entries("END", "LA", "LGT", "GQ", "DP", "LAD")
        vds = hl.vds.VariantDataset.from_merged_representation(mtds)

    if split:
        logger.info("Performing split-multi...")
        vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info("Filtering variants to predetermined QC variants...")
    vds = hl.vds.filter_variants(vds, predetermined_qc_sites.ht())

    logger.info("Densifying data...")
    mt = hl.vds.to_dense_mt(vds)

    mt = mt.select_entries("GT", "GQ", "DP", "AD")

    return mt


def generate_qc_mt(
    v3_mt: hl.MatrixTable,
    v4_mt: hl.MatrixTable,
    bi_allelic_only: bool = False,
    min_af: float = 0.0001,
    min_callrate: float = 0.99,
    min_inbreeding_coeff_threshold: float = -0.8,
    ld_r2: float = 0.1,
    n_partitions: int = 1000,
) -> hl.MatrixTable:
    """
    Generate combined gnomAD v3 and v4 QC MatrixTable for use in relatedness and ancestry inference.

    :param v3_mt: Dense gnomAD v3 MatrixTable filtered to predetermined sites.
    :param v4_mt: Dense gnomAD v4 MatrixTable filtered to predetermined sites.
    :param bi_allelic_only: Whether to filter to bi-allelic variants.
    :param min_af: Minimum variant allele frequency to retain variant in QC MatrixTable.
    :param min_callrate: Minimum variant callrate to retain variant in QC MatrixTable.
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to retain variant in QC MatrixTable.
    :param ld_r2: LD-pruning cutoff.
    :param n_partitions: Number of partitions to repartition the MT to before LD pruning.
    :return: MatrixTable of sites that pass QC filters.
    """
    v3_count = v3_mt.count()
    v4_count = v4_mt.count()
    logger.info(
        "Number of (variants, samples) in the v3.1 MatrixTable: %s...", v3_count
    )
    logger.info("Number of (variants, samples) in the v4 MatrixTable: %s...", v4_count)
    # Remove v4 hard filtered samples
    v4_mt = v4_mt.anti_join_cols(hard_filtered_samples.ht())

    samples_in_both = v4_mt.cols().semi_join(v3_mt.cols())
    n_samples_in_both = samples_in_both.count()
    if n_samples_in_both > 0:
        logger.info("The following samples are found in both MatrixTables:")
        samples_in_both.show(n_samples_in_both)
        raise ValueError("Overlapping sample names were found in the MatrixTables!")

    n_variants_in_both = v4_mt.rows().semi_join(v3_mt.rows()).count()
    logger.info(
        "Number of variants found in both MatrixTables: %d\n"
        "Number of variants found in only the v3.1 MatrixTable: %d\n"
        "Number of variants found in only the v4 MatrixTable: %d",
        n_variants_in_both,
        v3_count[0] - n_variants_in_both,
        v4_count[0] - n_variants_in_both,
    )

    logger.info(
        "Performing a union of the v3.1 and v4 predetermined QC site MatrixTable columns..."
    )
    mt = v3_mt.union_cols(v4_mt)

    logger.info("Annotating with adj and filtering to pre LD-pruned QC sites...")
    mt = annotate_adj(mt)

    cp_path = get_checkpoint_path("combined_pre_ld_prune_qc_mt", mt=True)
    mt = get_qc_mt(
        mt,
        bi_allelic_only=bi_allelic_only,
        min_af=min_af,
        min_callrate=min_callrate,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=None,  # Already filtered from the initial set of QC variants
        apply_hard_filters=False,
        ld_r2=ld_r2,
        checkpoint_path=cp_path,
        filter_lcr=False,  # Already filtered from the initial set of QC variants
        filter_decoy=False,  # Doesn't exist for hg38
        filter_segdup=False,  # Already filtered from the initial set of QC variants
        n_partitions=n_partitions,
    )
    logger.info(
        "Number of pre LD-pruned QC sites in gnomAD v3 + v4: %d...",
        hl.read_matrix_table(cp_path).count_rows(),
    )

    return mt


def generate_qc_meta_ht() -> hl.Table:
    """
    Combine v3 and v4 sample metadata into a single Table for relatedness and population inference.

    :return: Table with v3 and v4 sample metadata
    """
    v3_meta_ht = v3_meta.ht()
    v3_meta_ht = v3_meta_ht.select(
        v2_meta=hl.struct(
            v2_known_pop=v3_meta_ht.project_meta.v2_known_pop,
            v2_known_subpop=v3_meta_ht.project_meta.v2_known_subpop,
            v2_pop=v3_meta_ht.project_meta.v2_pop,
            v2_subpop=v3_meta_ht.project_meta.v2_subpop,
        ),
        v3_meta=hl.struct(
            v3_project_ancestry=v3_meta_ht.project_meta.project_ancestry,
            v3_project_pop=v3_meta_ht.project_meta.project_pop,
            v3_project_subpop=v3_meta_ht.project_meta.project_subpop,
            v3_subpop_description=v3_meta_ht.project_meta.subpop_description,
            v3_subsets=v3_meta_ht.subsets,
            v3_population_inference=v3_meta_ht.population_inference.select(
                "training_pop", "pca_scores", "pop"
            ),
            v3_relatedness_relationships=v3_meta_ht.relatedness_inference.relationships,
            v3_release_related=v3_meta_ht.sample_filters.release_related,
            v3_qc_metrics_outlier_filters=v3_meta_ht.sample_filters.qc_metrics_filters,
            v3_exclude=v3_meta_ht.project_meta.exclude,
            v3_exclude_reason=v3_meta_ht.project_meta.exclude_reason,
            v3_high_quality=v3_meta_ht.high_quality,
            v3_release=v3_meta_ht.release,
        ),
        chr20_mean_dp=v3_meta_ht.sex_imputation.chr20_mean_dp,
        releasable=v3_meta_ht.project_meta.releasable,
        broad_external=v3_meta_ht.project_meta.broad_external,
        hard_filters=v3_meta_ht.sample_filters.hard_filters,
        hard_filtered=v3_meta_ht.sample_filters.hard_filtered,
        data_type="genomes",
    )

    v4_meta_ht = v4_meta.ht()
    chr20_mean_dp_ht = sample_chr20_mean_dp.ht()
    hardfilter_ht = hard_filtered_samples.ht()

    v4_meta_ht = v4_meta_ht.select(
        v2_meta=hl.struct(
            v2_known_pop=v4_meta_ht.project_meta.v2_meta.v2_known_pop,
            v2_known_subpop=v4_meta_ht.project_meta.v2_meta.v2_known_subpop,
            v2_pop=v4_meta_ht.project_meta.v2_meta.v2_pop,
            v2_subpop=v4_meta_ht.project_meta.v2_meta.v2_subpop,
        ),
        v4_race_ethnicity=v4_meta_ht.project_meta.race_ethnicity,
        chr20_mean_dp=chr20_mean_dp_ht[v4_meta_ht.key].chr20_mean_dp,
        releasable=v4_meta_ht.project_meta.releasable,
        broad_external=v4_meta_ht.project_meta.broad_external,
        hard_filters=hardfilter_ht[v4_meta_ht.key].hard_filters,
        hard_filtered=hl.is_defined(hardfilter_ht[v4_meta_ht.key].hard_filters),
        data_type="exomes",
    )

    samples_in_both = v4_meta_ht.semi_join(v3_meta_ht)
    n_samples_in_both = samples_in_both.count()
    if n_samples_in_both > 0:
        logger.info("The following samples are found in both metadata Tables:")
        samples_in_both.show(n_samples_in_both)
        raise ValueError("Overlapping sample names were found in the metadata Tables!")

    return v3_meta_ht.union(v4_meta_ht, unify=True)


def main(args):
    hl.init(
        log="/generate_qc_mt.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # NOTE: remove this flag when the new shuffle method is the default. This is necessary for hail 0.2.97.
    hl._set_flags(use_new_shuffle="1")

    overwrite = args.overwrite
    ld_r2 = args.ld_r2
    test = args.test

    try:
        if args.create_v3_filtered_dense_mt:
            # Note: This command removes hard filtered samples
            predetermined_qc_path = get_predetermined_qc(version="3.1").test_path if test else get_predetermined_qc(version="3.1").path
            check_resource_existence(predetermined_qc_path, not overwrite)
            mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True, test=test)
            mt = create_filtered_dense_mt(mt, split=True)
            mt = mt.checkpoint(predetermined_qc_path, overwrite=overwrite)
            logger.info(
                "Number of predetermined QC variants found in the gnomAD v3 MatrixTable: %d...",
                mt.count_rows(),
            )

        if args.create_v4_filtered_dense_mt:
            # Note: This subset dense MatrixTable was created before the final hard filtering was determined
            # Hard filtering is performed in `generate_qc_mt` before applying variant filters
            predetermined_qc_path = get_predetermined_qc().test_path if test else get_predetermined_qc().path
            check_resource_existence(predetermined_qc_path, not overwrite)
            vds = get_gnomad_v4_vds(
                split=True, remove_hard_filtered_samples=False, test=test
            )
            mt = create_filtered_dense_mt(vds)
            mt = mt.checkpoint(predetermined_qc_path, overwrite=overwrite)
            logger.info(
                "Number of predetermined QC variants found in the gnomAD v4 VDS: %d...",
                mt.count_rows(),
            )

        if args.generate_qc_mt:
            joint_qc_path = joint_qc.test_path if test else joint_qc.path
            check_resource_existence(joint_qc_path, not overwrite)

            v3_mt = get_predetermined_qc(version="3.1").mt(test=test)
            v4_mt = get_predetermined_qc().mt(test=test)
            mt = generate_qc_mt(
                v3_mt,
                v4_mt,
                bi_allelic_only=args.bi_allelic_only,
                min_af=args.min_af,
                min_callrate=args.min_callrate,
                min_inbreeding_coeff_threshold=args.min_inbreeding_coeff_threshold,
                ld_r2=ld_r2,
                n_partitions=args.n_partitions,
            )
            mt.write(joint_qc_path, overwrite=overwrite)

        if args.generate_qc_meta:
            generate_qc_meta_ht().write(joint_qc_meta.path, overwrite=overwrite)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("generate_qc_mt"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the VDS.", action="store_true"
    )
    parser.add_argument(
        "--create-v3-filtered-dense-mt",
        help="Create a dense MatrixTable from the raw gnomAD v3.1 sparse MatrixTable filtered to predetermined QC variants.",
        action="store_true",
    )
    parser.add_argument(
        "--create-v4-filtered-dense-mt",
        help="Create a dense MatrixTable from the raw gnomAD v4 VariantDataset filtered to predetermined QC variants.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-qc-mt",
        help="Create the final merged gnomAD v3 + v4 QC MatrixTable with all specified filters and LD-pruning.",
        action="store_true",
    )
    parser.add_argument(
        "--bi-allelic-only",
        help="Filter to variants that are bi-allelic.",
        action="store_true",
    )
    parser.add_argument(
        "--min-af",
        help="Minimum variant allele frequency to retain variant in QC MatrixTable.",
        default=0.0001,
        type=float,
    )
    parser.add_argument(
        "--min-callrate",
        help="Minimum variant callrate to retain variant in QC MatrixTable.",
        default=0.99,
        type=float,
    )
    parser.add_argument(
        "--min-inbreeding-coeff-threshold",
        help="Minimum site inbreeding coefficient to retain variant in QC MatrixTable.",
        default=-0.8,
        type=float,
    )
    parser.add_argument(
        "--ld-r2",
        help="LD-pruning cutoff.",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--n-partitions",
        help="Desired number of partitions for output QC MatrixTable.",
        default=1000,
        type=int,
    )
    parser.add_argument(
        "--generate-qc-meta",
        help="Create a merged gnomAD v3 + v4 metadata Table for QC purposes.",
        action="store_true",
    )
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite all data from this subset.",
        action="store_true",
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
