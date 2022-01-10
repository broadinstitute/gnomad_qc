import argparse
import hail as hl
import logging

from gnomad.sample_qc.ancestry import run_pca_with_relateds

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
from gnomad_qc.v3.resources.constants import CURRENT_VERSION
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
    Generate the subpop QC MT to be used for all subpop analyses.
    
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

    # Densify the MT and include only sites defined in the qc_sites HT
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
    release: bool = False,
    high_quality: bool = False,
    outliers=None,
    ht_read_if_exists: bool = True,
    outlier_string="",
) -> hl.MatrixTable:
    """
    Generate the QC MT per specified population and generate PCA info.

    :param mt: The QC MT output by the 'compute_subpop_qc_mt' function
    :param pop: Population to which the Matrix Table should be filtered
    :param min_af: Minimum population variant allele frequency to retain variant in QC MT
    :param min_inbreeding_coeff_threshold: Minimum site inbreeding coefficient to retain variant in QC MT
    :param release: Whether or not to filter to only release samples
    :param high_quality: Whether or not to filter to only high quality samples
    :param outliers: Table keyed by column containing outliers to remove
    :param ht_read_if_exists: Whether or not to read an existing Table of PCA data if it exists
    :param outlier_string: String used to describe which outliers(if any) were removed
    :return: Table with sample metadata and PCA scores for the specified population to use for subpop analysis
    """
    meta_ht = meta.ht()
    relateds = pca_related_samples_to_drop.ht()
    hard_filtered_ht = hard_filtered_samples.ht()

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
    pop_mt = mt.filter_cols(mt.pop == pop)

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

    if not (ht_read_if_exists and file_exists(ht_path)):
        pop_pca_evals, pop_pca_scores, pop_pca_loadings = run_pca_with_relateds(
            pop_qc_mt, relateds
        )

        pop_ht = pop_pca_scores.annotate(**meta_ht[pop_pca_scores.key])
        pop_ht = pop_ht.annotate(
            subpop_description=pop_ht.project_meta.subpop_description,
            v2_pop=pop_ht.project_meta.v2_pop,
            v2_subpop=pop_ht.project_meta.v2_subpop,
            pop=pop_ht.population_inference.pop,
            project_id=pop_ht.project_meta.project_id,
            project_pop=pop_ht.project_meta.project_pop,
        )
        pop_ht = pop_ht.checkpoint(ht_path, overwrite=args.overwrite)
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
        + f"gnomad_v{CURRENT_VERSION}_qc_mt_subpop_analysis.mt"
    )

    if args.make_full_subpop_qc_mt:
        logger.info("Generating densified MT to use for all subpop analyses...")
        mt = compute_subpop_qc_mt(mt, args.min_popmax_af)
        mt = mt.naive_coalesce(5000)
        mt = mt.checkpoint(output_path, overwrite=args.overwrite)

    logger.info("Generating PCs for subpops...")
    # TODO: Add code for AMR using currently unused `run_pop_pca` function


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="This script generates a QC MT and PCA scores to use for subpop analyses"
    )
    parser.add_argument(
        "--min-popmax-af",
        help="Minimum population max variant allele frequency to retain a variant for the subpop QC MT",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--pop",
        help="Population to which the subpop QC MT should be filtered when generating the PCA data",
        type=str,
    )
    parser.add_argument(
        "--min-af",
        help="Minimum population variant allele frequency to retain variant in QC MT when generating the PCA data",
        type=float,
        default=0.001,
    )
    parser.add_argument(
        "--min-inbreeding-coeff-threshold",
        help="Minimum site inbreeding coefficient to keep a variant when generating the PCA data",
        type=float,
        default=-0.25,
    )
    parser.add_argument(
        "--outlier-ht-path",
        help="Path to Table keyed by column containing outliers to remove when generating the PCA data",
        type=str,
    )
    parser.add_argument(
        "--outlier-string",
        help="String used to describe which outliers(if any were set by 'outlier-ht-path') were removed when generating the PCA data",
        type=str,
        default="",
    )
    parser.add_argument(
        "--make-full-subpop-qc-mt",
        help="Runs function to create dense MT to use as QC MT for all subpop analyses. This uses --min-popmax-af to determine variants that need to be retained",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help="Runs a test of the code on only two partitions of the raw gnomAD v3 MT",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files", action="store_true"
    )
    parser.add_argument(
        "--release",
        help="Filter to only release samples when generating the PCA data",
        action="store_true",
    )
    parser.add_argument(
        "--high-quality",
        help="Filter to only high-quality samples when generating the PCA data",
        action="store_true",
    )
    parser.add_argument(
        "--remove-hard-filtered-samples",
        help="Remove hard-filtered samples when generating the PCA data",
        action="store_true",
    )
    parser.add_argument(
        "--ht-read-if-exists",
        help="Read an existing Table of PCA data if it exists",
        action="store_true",
    )
