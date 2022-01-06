import argparse
import hail as hl
import logging

from gnomad.sample_qc.ancestry import (
    assign_population_pcs,
    POP_COLORS,
    run_pca_with_relateds,
)

from gnomad_qc.v3.resources.meta import meta
from gnomad_qc.v3.resources.sample_qc import (
    hard_filtered_samples,
    pca_related_samples_to_drop,
)
from gnomad.resources.grch38.reference_data import lcr_intervals
from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.file_utils import file_exists
from gnomad.utils.sparse_mt import densify_sites
from gnomad.sample_qc.pipeline import get_qc_mt
from gnomad_qc.v3.resources import release_sites
from gnomad_qc.v3.resources.basics import get_gnomad_v3_mt
from gnomad_qc.v3.resources.annotations import get_info, last_END_position
from gnomad_qc.v3.resources.sample_qc import get_sample_qc_root


logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("subpop_analysis")
logger.setLevel(logging.INFO)


def compute_subpop_qc_mt(
    mt: hl.MatrixTable, min_popmax_af: float = 0.001,
) -> hl.MatrixTable:
    """
    Generate the subpop QC MT to be used for the subpop analyses.
    
    Filter MT to bi-allelic SNVs at sites with popmax allele frequency greater than `min_popmax_af`, densify, and write MT.

    :param mt: Raw MatrixTable to use for the subpop analysis
    :param min_popmax_af: Minimum population max variant allele frequency to retain variant for the subpop QC MatrixTable
    :return: MatrixTable filtered to variants for subpop analysis
    """
    release_ht = hl.read_table(release_sites().path)

    # Filter to biallelic SNVs not in low-confidence regions and with a popmax above min_popmax_af
    qc_sites = release_ht.filter(
        (release_ht.popmax.AF > min_popmax_af)
        & ~release_ht.was_split
        & hl.is_snp(release_ht.alleles[0], release_ht.alleles[1])
        & hl.is_missing(lcr_intervals.ht()[release_ht.locus])
    ).key_by("locus")

    # Densify the MT and include only sites defined in the qc_sites ht
    mt = mt.select_entries(
        "END", GT=mt.LGT, adj=get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD)
    )
    mt = densify_sites(mt, qc_sites, hl.read_table(last_END_position.path))

    return mt


def run_pop_pca(
    mt: hl.MatrixTable,
    pop: str,
    min_af: float = 0.001,
    min_inbreeding_coeff_threshold: float = -0.25,
    remove_hard_filtered_samples: bool = True,
    release=False,
    high_quality=False,
    outliers=None,
    mt_read_if_exists=True,
    ht_read_if_exists=True,
    outlier_string="",
) -> hl.MatrixTable:
    """
    Generate the QC MT per specified population and generate PCA info.

    :param mt: The QC MT output by the 'compute_subpop_qc_mt' function
    :param min_af: Minimum variant allele frequency to retain variant in QC matrix table
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to keep
    :return: MatrixTable filtered to specific population to use for subpop analysis
    """
    # Adjust names for different project subpops that should actually be the same
    meta_ht = meta.ht()
    meta_ht = meta_ht.annotate(
        project_meta=meta_ht.project_meta.annotate(
            project_subpop=(
                hl.switch(meta_ht.project_meta.project_subpop.lower())
                .when("peru", "pel")
                .when("colombian", "clm")
                .default(meta_ht.project_meta.project_subpop)
            )
        )
    )

    # Add info to the MT
    info_ht = get_info(split=False).ht()
    info_ht = info_ht.annotate(
        info=info_ht.info.select(
            # No need for AS_annotations since it's bi-allelic sites only
            **{x: info_ht.info[x] for x in info_ht.info if not x.startswith("AS_")}
        )
    )
    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)

    # Add sample metadata to the QC MT
    mt = mt.annotate_cols(**meta_ht[mt.col_key])

    # Filter the  QC MT to the specified pop
    pop_mt = qc_mt.filter_cols(mt.pop == pop)

    # Apply specified filters
    if remove_hard_filtered_samples:
        pop_mt = pop_mt.filter_cols(hl.is_missing(hard_filtered_ht[pop_mt.col_key]))
    if release:
        pop_mt = pop_mt.filter_cols(pop_mt.release)
    if high_quality:
        pop_mt = pop_mt.filter_cols(pop_mt.high_quality)
    if outliers is not None:
        pop_mt = pop_mt.filter_cols(hl.is_missing(outliers[pop_mt.col_key]))

    pop_mt = pop_mt.filter_rows(hl.agg.any(pop_mt.GT.is_non_ref()))

    # Generate a QC MT for the given pop
    pop_qc_mt = get_qc_mt(
        pop_mt,
        min_af=min_af,
        min_inbreeding_coeff_threshold=min_inbreeding_coeff_threshold,
        min_hardy_weinberg_threshold=None,
        ld_r2=None,
        filter_lcr=False,
        filter_decoy=False,
        filter_segdup=False,
    )

    # Generate PCA data
    release_string = "_release" if release else ""
    high_quality_string = "_high_quality" if high_quality else ""
    ht_path = (
        get_sample_qc_root()
        + f"/subpop_analysis/{pop}/{pop}_scores{release_string}{high_quality_string}{outlier_string}.ht"
    )

    relateds = pca_related_samples_to_drop.ht()

    if not (ht_read_if_exists and file_exists(ht_path)):
        pop_pca_evals, pop_pca_scores, pop_pca_loadings = run_pca_with_relateds(
            pop_qc_mt, relateds
        )

        pop_ht = pop_pca_scores.annotate(**meta_ht[pop_pca_scores.key])
        pop_ht = pop_ht.annotate(
            project_subpop=pop_ht.project_meta.project_subpop,
            subpop_description=pop_ht.project_meta.subpop_description,
            v2_pop=pop_ht.project_meta.v2_pop,
            v2_subpop=pop_ht.project_meta.v2_subpop,
            pop=pop_ht.population_inference.pop,
            project_id=pop_ht.project_meta.project_id,
            project_pop=pop_ht.project_meta.project_pop,
        )
        pop_ht = pop_ht.checkpoint(ht_path, overwrite=True)
    else:
        pop_ht = hl.read_table(ht_path)

    return pop_ht


def main(args):

    # Read in the raw mt
    mt = get_gnomad_v3_mt(key_by_locus_and_alleles=True)

    # Filter to test partitions if specified
    if args.test:
        logger.info("Filtering to the first two partitions of the MT")
        mt = mt._filter_partitions(range(2))

    # Write out the densified MT
    output_path = (
        get_sample_qc_root()
        + "/subpop_analysis/"
        + "gnomad_v3_qc_mt_subpop_analysis.mt"
    )

    if args.make_subpop_mt:
        logger.info("Generating densified MT to use for subpop analysis...")
        mt = compute_subpop_qc_mt(mt, args.min_popmax_af)
        mt = mt.naive_coalesce(5000)
        mt = mt.checkpoint(output_path, overwrite=True,)

    logger.info("Generating PCs for subpops...")
    # TODO: Add code for AMR


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script generate a qc mt to use for the subpop analysis"
    )
    parser.add_argument(
        "--min-popmax-af",
        help="Minimum population max variant allele frequency to retain variant for the subpop QC matrix table",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--make-subpop-mt",
        help="Runs function to create dense MT to use as QC MT for subpop analysis",
        action="store_true",
    )
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the MT", action="store_true"
    )
