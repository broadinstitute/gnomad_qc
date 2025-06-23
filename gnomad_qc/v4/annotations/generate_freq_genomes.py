"""
Script to create frequencies HT for v4.0 genomes.

This script is written specifically for the v4.0 genomes release to handle the
addition and subtraction of samples from the HGDP + 1KG subset. The HGDP + 1KG
subset has updated sample QC annotations. The addition of new samples from the
densified MT will include new variants, which will require the recomputation of AN
for the v3.1 releases samples, and merging it to call stats of the v3.1.2 release
sites HT and of updated HGDP + 1KG subset. In addition, the code will also get the
other call stats related annotations, including: filtering allele frequencies ('faf'),
'grpmax', 'gen_anc_faf_max' and 'InbreedingCoeff' as done in v4.0 exomes.
"""

import argparse
import logging
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
from gnomad.resources.grch38.gnomad import POPS_TO_REMOVE_FOR_POPMAX, SUBSETS
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.annotations import (
    annotate_adj,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    gen_anc_faf_max_expr,
    generate_freq_group_membership_array,
    merge_freq_arrays,
    missing_callstats_expr,
    pop_max_expr,
    set_female_y_metrics_to_na_expr,
    update_structured_annotations,
)
from gnomad.utils.filtering import filter_arrays_by_meta
from gnomad.utils.release import make_freq_index_dict_from_meta
from gnomad.utils.slack import slack_notifications
from hail.utils import new_temp_file

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.annotations import get_freq
from gnomad_qc.v3.resources.basics import get_gnomad_v3_vds
from gnomad_qc.v3.resources.basics import meta as v3_meta
from gnomad_qc.v3.resources.release import (
    hgdp_tgp_subset,
    hgdp_tgp_subset_annotations,
    release_sites,
)
from gnomad_qc.v3.resources.sample_qc import relatedness as v3_relatedness
from gnomad_qc.v3.utils import hom_alt_depletion_fix
from gnomad_qc.v4.resources.annotations import get_freq as v4_get_freq
from gnomad_qc.v4.resources.annotations import hgdp_tgp_updated_callstats
from gnomad_qc.v4.resources.basics import meta as v4_meta
from gnomad_qc.v4.resources.sample_qc import (
    hgdp_recomputed_freemix,
    hgdp_tgp_duplicated_to_exomes,
    hgdp_tgp_meta_updated,
    hgdp_tgp_pop_outliers,
    hgdp_tgp_populations_updated,
    hgdp_tgp_related_samples_to_drop,
    hgdp_tgp_related_to_nonsubset,
)
from gnomad_qc.v4.resources.sample_qc import relatedness as v4_relatedness

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger(
    "Create v4.0 frequencies HT with updated HGDP/TGP metadata and new annotations"
)
logger.setLevel(logging.INFO)

SUBSETS = SUBSETS["v3"]
"""List of subsets in the v3.1/v4.0 genomes release."""
POP_MAP = {
    "bantusafrica": "bantusouthafrica",
    "biakaPygmy": "biaka",
    "italian": "bergamoitalian",
    "mbutiPygmy": "mbuti",
    "melanesian": "bougainville",
    "mongola": "mongolian",
    "miaozu": "miao",
    "yizu": "yi",
    "oth": "remaining",
}
"""
Map of HGDP populations that need to be renamed in the v4.0 genomes release. Also includes
the renaming of 'oth' to 'remaining'.
"""
JOIN_FREQS = ["release", "pop_diff", "added", "subtracted"]
"""Frequency Tables to join for creation of the v4.0 genomes release sites HT."""
FREQ_GLOBALS = ("freq_meta", "freq_meta_sample_count")
"""Global annotations on the frequency Table."""


def filter_to_test(
    t: Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset],
    gene_on_chrx: bool = False,
) -> Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Filter to a region in PCSK9 or TRPC5 in Table or MatrixTable for testing purposes.

    :param t: Table or MatrixTable to filter.
    :param gene_on_chrx: Whether to filter to TRPC5 (in the non-PAR region), instead of
        PCSK9, for testing chrX.
    :return: Table or MatrixTable filtered to a region in PCSK9 or TRPC5.
    """
    if gene_on_chrx:
        logger.info("Filtering to TRPC5 on chrX in MT for testing purposes...")
        test_locus = "chrX:111776000-111786000"
    else:
        logger.info("Filtering to 10kb in PCSK9 in MT for testing purposes...")
        test_locus = "chr1:55039447-55064852"

    test_interval = [hl.parse_locus_interval(test_locus, reference_genome="GRCh38")]

    if isinstance(t, hl.vds.VariantDataset):
        t = hl.vds.filter_intervals(t, test_interval, split_reference_blocks=True)
    else:
        t = hl.filter_intervals(t, test_interval)

    return t


def get_hgdp_tgp_related_to_nonsubset(
    v3_meta_ht: hl.Table, rel_ht: hl.Table
) -> hl.Table:
    """
    Get the samples in the HGDP + 1KG subset that were related to samples outside the subset in the v3.1 release and were not included in the v3.1 release.

    :param v3_meta_ht: Table with the v3.1 release metadata.
    :param rel_ht: Table with the v3.1 release relatedness, here we use the results
        from pc_relate of the full dataset.
    :return: Table with the samples in the HGDP + 1KG subset that are related to
        samples outside the subset in v3 release.
    """
    # Filter to samples that are either in the HGDP + 1KG subset or in the v3.1 release.
    v3_meta_ht = v3_meta_ht.annotate(
        hgdp_tgp=v3_meta_ht.subsets.tgp | v3_meta_ht.subsets.hgdp
    )
    v3_meta_ht = v3_meta_ht.filter(v3_meta_ht.release | v3_meta_ht.hgdp_tgp)

    rel_ht = rel_ht.annotate(
        **{
            f"{x}_meta": hl.struct(
                v3_release=hl.coalesce(v3_meta_ht[rel_ht[x].s].release, False),
                hgdp_tgp=hl.coalesce(v3_meta_ht[rel_ht[x].s].hgdp_tgp, False),
            )
            for x in ["i", "j"]
        }
    )

    rel_ht = rel_ht.filter(
        # Filter to pairs with 2nd degree or closer relatedness.
        (rel_ht.relationship != "unrelated")
        # Filter to pairs where at least one of the samples was in the v3.1 release.
        & (rel_ht.i_meta.v3_release | rel_ht.j_meta.v3_release)
        # Filter to pairs where one and only one of the samples is in the HGDP +
        # 1KG subset.
        & ((hl.int(rel_ht.i_meta.hgdp_tgp) + hl.int(rel_ht.j_meta.hgdp_tgp)) == 1)
        # Exclude pairs where the HGDP/1KG sample in the pair is also a v3.1 release
        # sample.
        & ~(rel_ht.i_meta.hgdp_tgp & rel_ht.i_meta.v3_release)
        & ~(rel_ht.j_meta.hgdp_tgp & rel_ht.j_meta.v3_release)
    )

    # Filter the relatedness HT twice to get the HGDP/1KG sample in each pair.
    ht1 = rel_ht.filter(rel_ht.i_meta.hgdp_tgp).key_by().select("i")
    ht1 = ht1.select(s=ht1.i.s)
    ht2 = rel_ht.filter(rel_ht.j_meta.hgdp_tgp).key_by().select("j")
    ht2 = ht2.select(s=ht2.j.s)

    # Merge the two HTs and remove the v3.1 prefix from the sample IDs.
    ht = ht1.union(ht2)
    ht = ht.select(s=ht.s.replace("v3.1::", "")).key_by("s").distinct()

    ht = ht.naive_coalesce(1).checkpoint(new_temp_file("related_to_nonsubset_ht"))
    logger.info(
        "%d HGDP/1KG samples are related to samples outside the subset", ht.count()
    )

    return ht.select_globals()


def get_hgdp_tgp_v4_exome_duplicates(
    v3_meta_ht: hl.Table,
    v4_meta_ht: hl.Table,
    rel_ht: hl.Table,
) -> hl.Table:
    """
    Get the samples in the HGDP + 1KG subset that are duplicates of an exome in the v4.0 release.

    The duplicated samples are defined as samples that were HGDP + 1KG subset samples
    in the v3.1 release and are also in the v4.0 exomes release. The duplicated samples
    have to be removed because we will have combined frequencies rom v4.0 exomes and
    genomes.

    :param v3_meta_ht: Table with the v3.1 release metadata.
    :param v4_meta_ht: Table with the v4.0 exomes release metadata.
    :param rel_ht: Table with the v3.1 and v4.0 joint relatedness, it's based on
        cuKING relatedness results.
    :return: Table with the samples in the HGDP + 1KG subset that are duplicates of an
        exome in the v4.0 exomes release.
    """
    # Get HGDP + 1KG subset samples in v3.1 meta.
    v3_meta_ht = v3_meta_ht.filter(v3_meta_ht.subsets.hgdp | v3_meta_ht.subsets.tgp)

    # Get release samples in v4.0 exomes.
    v4_meta_ht = v4_meta_ht.filter(v4_meta_ht.release)

    # Get samples that are in the v3.1 genomes and are also in the v4.0 exomes.
    rel_ht = rel_ht.filter(rel_ht.gnomad_v3_duplicate)

    # Check if the duplicates are included in the v4.0 exomes release and belong to the
    # HGDP + 1KG subset.
    ht = rel_ht.annotate(
        **{
            f"{x}_meta": hl.struct(
                in_v4_exomes_release=hl.is_defined(v4_meta_ht[rel_ht[x].s]),
                in_hgdp_tgp_subset=hl.is_defined(v3_meta_ht[rel_ht[x].s]),
            )
            for x in ["i", "j"]
        }
    )

    ht = ht.filter(
        (ht.i_meta.in_v4_exomes_release & ht.j_meta.in_hgdp_tgp_subset)
        | (ht.j_meta.in_v4_exomes_release & ht.i_meta.in_hgdp_tgp_subset)
    )

    # Filter the relatedness HT twice to get the HGDP/1KG sample in each pair.
    ht1 = ht.filter(ht.i_meta.in_hgdp_tgp_subset).key_by().select("i")
    ht1 = ht1.select(s=ht1.i.s)
    # This table is likely empty because of the ordering when computing relatedness,
    # but included for completeness.
    ht2 = ht.filter(ht.j_meta.in_hgdp_tgp_subset).key_by().select("j")
    ht2 = ht2.select(s=ht2.j.s)

    # Merge the two HTs and remove the v3.1 prefix from the sample IDs.
    ht = ht1.union(ht2)
    ht = ht.select(s=ht.s.replace("v3.1::", "")).key_by("s").distinct()

    ht = ht.naive_coalesce(1).checkpoint(new_temp_file("duplicate_in_v4_exomes"))
    logger.info(
        "%d HGDP/1KG samples are duplicated in the v4.0 exomes release", ht.count()
    )

    return ht.select_globals()


def add_updated_sample_qc_annotations(ht: hl.Table) -> hl.Table:
    """
    Add updated sample QC annotations to the HGDP + 1KG subset.

    .. note::

        The following annotations need to be updated for the v4.0 genomes release based
        on the latest sample QC results of the subset:

            - `hgdp_tgp_meta.subcontinental_pca.outlier`: to apply the updated pop
              outlier filter implemented by Alicia Martin's group.
            - `gnomad_sample_filters.hard_filtered`: to apply the recomputed freemix
              filter for HGDP samples.
            - `gnomad_sample_filters.related_to_nonsubset`: to further filter out
              samples that are related to samples outside the subset but were not
              included in the v3 release.
            - `gnomad_sample_filters.v4_exome_duplicate`: to further filter out the
              samples in the HGDP + 1KG subset that are duplicates of an exome in the
              v4.0 release.
            - `relatedness_inference.related`: to apply the updated relatedness
              inference implemented by Alicia Martin's group.

    :param ht: Table with the HGDP + 1KG subset metadata from the v3.1.2 release.
    :return: Table with updated sample QC annotations.
    """
    # Load all updated sample QC Tables.
    contamination_ht = hgdp_recomputed_freemix.ht()
    pop_outliers_ht = hgdp_tgp_pop_outliers.ht()
    populations_ht = hgdp_tgp_populations_updated.ht()
    relatedness_ht = hgdp_tgp_related_samples_to_drop.ht()
    related_to_nonsubset_ht = hgdp_tgp_related_to_nonsubset.ht()
    duplicated_to_exomes_ht = hgdp_tgp_duplicated_to_exomes.ht()

    # TODO: get the global & subcontinental PCs for the updated HGDP + 1KG subset:
    #  hgdp_tgp_meta.global_pca_scores,
    #  hgdp_tgp_meta.subcontinental_pca.pca_scores,
    #  hgdp_tgp_meta.subcontinental_pca.pca_scores_outliers_removed,
    #  could update the meta HT after MVP release

    # Update annotations impacted by the recomputed freemix.
    freemix = hl.coalesce(
        contamination_ht[ht.s].recomputed_contam_rate * 100, ht.bam_metrics.freemix
    )
    hard_filters = ht.gnomad_sample_filters.hard_filters | hl.if_else(
        freemix > 5.0, hl.set({"contamination"}), hl.empty_set(hl.tstr)
    )
    hard_filtered = hl.len(hard_filters) > 0
    # Update annotations impacted by the updated relatedness inference.
    related_subset = hl.is_defined(relatedness_ht[ht.key])
    # Update the relatedness to nonsubset samples.
    related_nonsubset = hl.is_defined(related_to_nonsubset_ht[ht.key])
    # Get the duplicated samples in the HGDP + 1KG subset.
    duplicated_to_exomes = hl.is_defined(duplicated_to_exomes_ht[ht.key])
    # Update annotations impacted by the updated HGDP + 1KG metadata.
    outlier = hl.is_defined(pop_outliers_ht[ht.key])
    populations = populations_ht[ht.key]
    # Update gnomAD inferred population 'oth' to be 'remaining'.
    gnomad_pop = ht.gnomad_population_inference.pop.replace("oth", "remaining")

    gnomad_release = (
        ~hard_filtered
        & ~outlier
        & ~related_subset
        & ~related_nonsubset
        & ~duplicated_to_exomes
    )

    sample_annotations_to_update = {
        "bam_metrics": {"freemix": freemix},
        "gnomad_population_inference": {"pop": gnomad_pop},
        "gnomad_sample_filters": {
            "hard_filters": hard_filters,
            "hard_filtered": hard_filtered,
            "related_to_nonsubset": related_nonsubset,
            "v4_exome_duplicate": duplicated_to_exomes,
        },
        "gnomad_high_quality": ht.gnomad_high_quality & ~hard_filtered,
        "gnomad_release": gnomad_release,
        "relatedness_inference": {
            "related": related_subset,
        },
        "hgdp_tgp_meta": {
            "subcontinental_pca": {"outlier": outlier},
            "population": populations.population,
            "latitude": populations.latitude,
            "longitude": populations.longitude,
            "gnomad_labeled_subpop": populations.population.lower(),
        },
        "high_quality": ~hard_filtered & ~outlier,
    }
    updated_ht = update_structured_annotations(
        ht,
        sample_annotations_to_update,
        annotation_update_label="sample_annotations_updated",
    )
    updated_ht = updated_ht.annotate()

    updated_ht = updated_ht.checkpoint(new_temp_file("hgdp_tgp_meta_update", "ht"))
    updated_counts = updated_ht.aggregate(
        hl.struct(
            n_hard_filtered=hl.agg.count_where(
                updated_ht.gnomad_sample_filters.hard_filtered
            ),
            n_related_subset=hl.agg.count_where(
                updated_ht.relatedness_inference.related
            ),
            n_related_nonsubset=hl.agg.count_where(
                ~updated_ht.relatedness_inference.related
                & updated_ht.gnomad_sample_filters.related_to_nonsubset
            ),
            n_duplicate_to_exomes=hl.agg.count_where(
                ~updated_ht.relatedness_inference.related
                & ~updated_ht.gnomad_sample_filters.related_to_nonsubset
                & updated_ht.gnomad_sample_filters.v4_exome_duplicate
            ),
            n_outlier=hl.agg.count_where(
                updated_ht.hgdp_tgp_meta.subcontinental_pca.outlier
            ),
            n_release=hl.agg.count_where(updated_ht.gnomad_release),
            n_diff=hl.agg.explode(
                lambda x: hl.agg.counter(x), updated_ht.sample_annotations_updated
            ),
        )
    )
    logger.info(
        """%d sample(s) hard filtered (%d more samples compared to v3.1.2 subset release);
        %d samples are pop outliers (%d samples different compared to v3.1.2 subset release);
        %d samples have different population labels compared to v3.1.2 subset release;
        %d samples related within the subset (%d samples different compared to v3.1.2 subset release);
        %d samples further filtered out due to their relatedness to samples outside the subset;
        %d samples filtered out because they are duplicated in the v4.0 exomes release;
        %d samples will be in the v4.0 release, compared to 3280 in the v3.1.2 release.""",
        updated_counts["n_hard_filtered"],
        updated_counts["n_diff"]["gnomad_sample_filters.hard_filtered"],
        updated_counts["n_outlier"],
        updated_counts["n_diff"]["hgdp_tgp_meta.subcontinental_pca.outlier"],
        updated_counts["n_diff"]["hgdp_tgp_meta.population"],
        updated_counts["n_related_subset"],
        updated_counts["n_diff"]["relatedness_inference.related"],
        updated_counts["n_related_nonsubset"],
        updated_counts["n_duplicate_to_exomes"],
        updated_counts["n_release"],
    )
    return updated_ht


def get_updated_release_samples(ht: hl.Table) -> Tuple[hl.Table, hl.Table, hl.Table]:
    """
    Get the samples in the HGDP + 1KG subset that will be added, subtracted, or have different pop labels in the v4.0 release compared to the v3.1 release.

    Three sets of samples will be obtained:
        - samples that will have different pop labels in the v4.0 genomes release:
          samples in the to-be-split 'Han' and 'Papuan' populations AND their
          'gnomad_release' status hasn't changed.
        - samples that will be added to the v4.0 genomes release: samples where
          'gnomad_release' status has changed and 'gnomad_release' is now True.
        - samples that will be removed from the v4.0 genomes release: samples where
          'gnomad_release' status has changed and 'gnomad_release' is now False.

    :param ht: Table with the HGDP + 1KG subset metadata with updated sample
        annotations.
    :return: Tuple of Tables with samples that have different pop labels in the v4
        release, to be added, and to be subtracted.
    """
    # Filter to all samples with a change in the 'gnomad_release' status.
    release_diff_ht = ht.filter(
        ht.sample_annotations_updated.contains("gnomad_release")
    )
    release_add_ht = release_diff_ht.filter(release_diff_ht.gnomad_release)
    release_subtract_ht = release_diff_ht.filter(~release_diff_ht.gnomad_release)
    logger.info(
        """
        %d samples will be added after the new filters
        %d sample(s) will be removed after the new filters""",
        release_add_ht.count(),
        release_subtract_ht.count(),
    )

    # Filter to all samples with a pop name split and no change in release status.
    split_pops = hl.literal(["PapuanHighlands", "PapuanSepik", "Han", "NorthernHan"])
    pop_diff_ht = ht.filter(
        split_pops.contains(ht.hgdp_tgp_meta.population)
        & ht.gnomad_release
        & ~ht.sample_annotations_updated.contains("gnomad_release")
    )
    logger.info("%d samples in the split Papuan & Han set", pop_diff_ht.count())

    return pop_diff_ht, release_add_ht, release_subtract_ht


def calculate_callstats_for_selected_samples(
    mt: hl.MatrixTable,
    samples_ht: hl.Table,
    subsets: Optional[List[str]] = None,
) -> hl.Table:
    """
    Calculate the call stats for samples in `samples_ht`.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with selected samples and their metadata.
    :param subsets: Optional List of subsets to be included in the frequency metadata
        'freq_meta' annotation returned by `annotate_freq`. e.g. ['hgdp'], ['tgp'],
        ['hgdp', 'tgp'].
    :return: Table with the call stats for the selected samples.
    """
    logger.info("Filtering MT to sample subset and annotating with metadata...")
    mt = mt.filter_cols(hl.is_defined(samples_ht[mt.col_key]))
    mt = mt.annotate_cols(**samples_ht[mt.col_key])

    ht = annotate_freq(
        mt,
        sex_expr=mt.gnomad_sex_imputation.sex_karyotype,
        pop_expr=(
            mt.hgdp_tgp_meta.population.lower()
            if subsets
            else mt.gnomad_population_inference.pop
        ),
        annotate_mt=False,
    )

    if subsets:
        freq_meta = [
            {**x, **{"subset": "|".join(subsets)}} for x in hl.eval(ht.freq_meta)
        ]
        ht = ht.annotate_globals(freq_meta=freq_meta)

    return ht.checkpoint(new_temp_file("hgdp_tgp_subset_freq", "ht"))


def concatenate_annotations(
    hts: List[hl.Table],
    global_field_names: Union[List[str], Tuple[str]] = FREQ_GLOBALS,
    row_field_names: Union[List[str], Tuple[str]] = ("freq",),
) -> hl.Table:
    """
    Concatenate annotations on multiple Tables into a single Table.

    :param hts: List of Tables with annotations to be concatenated.
    :param global_field_names: Global field names to concatenate.
    :param row_field_names: Row field names to concatenate.
    :return: Table with concatenated annotations.
    """
    global_field_names = list(global_field_names)
    row_field_names = list(row_field_names)
    hts = [
        ht.select_globals(*global_field_names).select(*row_field_names) for ht in hts
    ]
    concat_ht = hl.Table.multi_way_zip_join(
        hts,
        data_field_name="ann_array",
        global_field_name="global_array",
    )

    concat_ht = concat_ht.transmute(
        **{f: concat_ht.ann_array.flatmap(lambda x: x[f]) for f in row_field_names}
    )
    concat_ht = concat_ht.transmute_globals(
        **{
            g: concat_ht.global_array.flatmap(lambda x: x[g])
            for g in global_field_names
        },
    )

    return concat_ht


def get_hgdp_tgp_callstats_for_selected_samples(
    mt: hl.MatrixTable, samples_ht: hl.Table, compute_freq_all: bool = True
) -> hl.Table:
    """
    Calculate call stats for samples in `samples_ht` grouped by all stratifications in the v3.1 release call stats.

    This function calculates call stats for:
        - all samples provided in the `samples_ht` (if `compute_freq_all` is True).
        - samples in only the HGDP subset.
        - samples in only the 1KG (tgp) subset.

    The call stats in these three Tables are then concatenate to create a single Table
    with all call stats.

    :param mt: MatrixTable with the HGDP + 1KG subset.
    :param samples_ht: Table with the samples to filter to before computing call stats.
    :param compute_freq_all: Whether to compute the call stats for all samples in
        `samples_ht` instead of only the subset call stats. If False, only the HGDP and
        1KG subset call stats are computed. Default is True.
    :return: Table with the call stats for the selected samples.
    """
    projects = []
    subset_freq_hts = []
    logger.info("Calculating call stats for selected samples...")
    if compute_freq_all:
        projects.append(("all", None))

    projects.extend([("HGDP", ["hgdp"]), ("1000 Genomes", ["tgp"])])
    project_expr = samples_ht.hgdp_tgp_meta.project
    project_n = samples_ht.aggregate(hl.agg.counter(project_expr))

    for project, subset in projects:
        if project == "all":
            subset_ht = samples_ht
        else:
            subset_ht = samples_ht.filter(project_expr == project)

        if project == "all" or project_n.get(project, 0) > 0:
            freq_ht = calculate_callstats_for_selected_samples(mt, subset_ht, subset)
            freq_ht = freq_ht.select("freq").select_globals(*FREQ_GLOBALS)
            subset_freq_hts.append(freq_ht)

    logger.info("Concatenating call stats Tables for selected samples...")
    freq_ht = concatenate_annotations(subset_freq_hts)

    return freq_ht


def filter_freq_arrays(
    ht: hl.Table,
    items_to_filter: Union[List[str], Dict[str, List[str]]],
    keep: bool = True,
    combine_operator: str = "and",
    annotations: Union[List[str], Tuple[str]] = ("freq",),
) -> hl.Table:
    """
    Use `filter_arrays_by_meta` to filter annotations and globals by 'freq_meta' and annotate them back on `ht`.

    Filter `annotations`, 'freq_meta' and 'freq_meta_sample_count' fields to only
    `items_to_filter` by using the 'freq_meta' array field values.

    :param ht: Input Table.
    :param items_to_filter: Items to filter by.
    :param keep: Whether to keep or remove items. Default is True.
    :param combine_operator: Operator ("and" or "or") to use when combining items in
        'items_to_filter'. Default is "and".
    :param annotations: Annotations in 'ht' to filter by `items_to_filter`.
    :return: Table with filtered 'annotations' and 'freq_meta' array fields.
    """
    freq_meta, array_exprs = filter_arrays_by_meta(
        ht.freq_meta,
        {
            **{a: ht[a] for a in annotations},
            "freq_meta_sample_count": ht.index_globals().freq_meta_sample_count,
        },
        items_to_filter=items_to_filter,
        keep=keep,
        combine_operator=combine_operator,
    )

    ht = ht.annotate(**{a: array_exprs[a] for a in annotations})
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_meta_sample_count=array_exprs["freq_meta_sample_count"],
    )

    return ht


def annotate_v3_subsets_sample_count(ht: hl.Table) -> hl.Table:
    """
    Get the freq sample counts from the v3.1 subsets.

    The freq sample counts on the v3.1.2 release sites HT only contains the counts for
    the non subset groupings. This function adds the freq sample counts from the v3.1
    subset frequency groups to the v3.1.2 sample count array to give a sample count
    array matching the v3.1.2 'freq_meta' global annotation.

    :param ht: Table of the 3.1.2 release sites HT, which only contains the freq
       sample counts for non subset groupings.
    :return: Table with the freq sample counts from the v3.1 subset groups added.
    """
    sc = ht.index_globals().freq_sample_count

    for subset in SUBSETS:
        freq_subset_ht = get_freq(
            version="3.1", subset=subset, het_nonref_patch=True
        ).ht()
        sc_subset = freq_subset_ht.index_globals().freq_sample_count
        sc = sc.extend(sc_subset)

    ht = ht.annotate_globals(freq_meta_sample_count=sc)

    return ht


def update_pop_labels(ht: hl.Table, pop_map: Dict[str, Any]) -> hl.Table:
    """
    Update the population labels in the 'freq_meta' annotation on `ht`.

    :param ht: Table with population labels to update.
    :param pop_map: Dictionary mapping old to new population labels.
    :return: Table with updated population labels in 'freq_meta' annotation.
    """
    pop_map = hl.literal(pop_map)
    ht = ht.annotate_globals(
        freq_meta=ht.freq_meta.map(
            lambda d: d.map_values(lambda x: hl.or_else(pop_map.get(x), x))
        )
    )

    return ht


def join_release_ht_with_subsets(
    ht: hl.Table, freq_hts: Dict[str, hl.Table]
) -> hl.Table:
    """
    Join the release HT with the call stats for added and subtracted samples.

    :param ht: gnomAD v3.1.2 release sites Table.
    :param freq_hts: Dictionary with the call stats Tables for added and subtracted
        samples.
    :return: Table with `ht` and all Tables in `freq_hts` joined.
    """
    logger.info("Adding freq sample counts from v3 subsets...")
    ht = annotate_v3_subsets_sample_count(ht)

    logger.info("Removing 'Han' and 'Papuan' pops from freq and freq_meta...")
    pops_to_remove = {"pop": ["han", "papuan"]}
    ht = filter_freq_arrays(ht, pops_to_remove, keep=False, combine_operator="or")

    logger.info("Updating HGDP pop labels and 'oth' -> 'remaining'...")
    ht = update_pop_labels(ht, POP_MAP)

    freq_hts = {"release": ht, **freq_hts}
    for name, freq_ht in freq_hts.items():
        logger.info("There are %i variants in the %s HT...", freq_ht.count(), name)
        freq_hts[name] = freq_ht.select("freq").select_globals(*FREQ_GLOBALS)

    logger.info("Joining release HT with AFs for added and subtracted samples...")
    freq_hts = [freq_hts[x] for x in JOIN_FREQS]
    ht = hl.Table.multi_way_zip_join(freq_hts, "ann_array", "global_array")

    logger.info("Filtering joint HT to rows where at least one group has AC > 0...")
    ht = ht.filter(hl.any(ht.ann_array.map(lambda x: x.freq.any(lambda y: y.AC > 0))))

    return ht


def get_group_membership_ht_for_an(
    ht: hl.Table, only_pop_diff: bool = False
) -> hl.Table:
    """
    Generate a Table with a 'group_membership' array for each sample indicating whether the sample belongs to specific stratification groups.

    `ht` must have the following annotations:
        - 'meta.population_inference.pop': population label.
        - 'meta.sex_imputation.sex_karyotype`: sex label.
        - 'meta.project_meta.project_subpop': subpopulation label.
        - 'meta.subsets': dictionary with subset labels as keys and boolean values.

    :param ht: Table with the sample metadata.
    :param only_pop_diff: Whether to only include the population stratification groups
        for samples in pop_diff. Default is False.
    :return: Table with the group membership for each sample to be used for computing
        allele number (AN) per group.
    """
    if only_pop_diff:
        pop_expr = ht.hgdp_tgp_meta.population.lower()
        subpop_expr = ht.hgdp_tgp_meta.population.lower()
        sex_expr = ht.gnomad_sex_imputation.sex_karyotype
        subsets = ["hgdp"]
    else:
        pop_expr = ht.meta.population_inference.pop
        sex_expr = ht.meta.sex_imputation.sex_karyotype
        subpop_expr = ht.meta.project_meta.project_subpop
        subsets = ["all"] + SUBSETS

    hts = []
    for subset in subsets:
        if subset in ["tgp", "hgdp"]:
            subset_pop_expr = subpop_expr
        else:
            subset_pop_expr = pop_expr

        # Build strata expressions to get group membership for pop, sex, pop-sex.
        strata_expr = [
            {"pop": subset_pop_expr},
            {"sex": sex_expr},
            {"pop": subset_pop_expr, "sex": sex_expr},
        ]
        # For subsets, add strata expressions to get group membership for subset,
        # subset-pop, subset-sex, and subset-pop-sex.
        if subset != "all" and not only_pop_diff:
            subset_expr = hl.or_missing(ht.meta.subsets[subset], subset)
            strata_expr = [{"subset": subset_expr}] + [
                {"subset": subset_expr, **x} for x in strata_expr
            ]

        subset_ht = generate_freq_group_membership_array(
            ht, strata_expr, remove_zero_sample_groups=True
        )
        subset_globals = subset_ht.index_globals()
        if subset != "all" and not only_pop_diff:
            # The group_membership, freq_meta, and freq_meta_sample_count arrays for
            # subsets is ordered [adj, raw, subset adj, ...] and we want it to be
            # [subset adj, subset raw, ...] where subset adj and subset raw have the
            # same sample grouping.
            subset_adj_idx = 2
            subset_ht = subset_ht.annotate(
                group_membership=hl.array(
                    [subset_ht.group_membership[subset_adj_idx]]
                ).extend(subset_ht.group_membership[subset_adj_idx:])
            )
            freq_meta_sample_count = hl.array(
                [subset_globals.freq_meta_sample_count[subset_adj_idx]]
            ).extend(subset_globals.freq_meta_sample_count[subset_adj_idx:])

            # Modify the freq_meta annotation to match the new group_membership and
            # freq_meta_sample_count arrays.
            subset_ht = subset_ht.annotate_globals(
                freq_meta=hl.array(
                    [
                        {"subset": subset, "group": "adj"},
                        {"subset": subset, "group": "raw"},
                    ]
                ).extend(subset_globals.freq_meta[subset_adj_idx + 1 :]),
                freq_meta_sample_count=freq_meta_sample_count,
            )

        if only_pop_diff:
            # We only want:
            # {'group': 'adj', 'pop': 'han', 'subset': 'hgdp'},
            # {'group': 'adj', 'pop': 'northernhan', 'subset': 'hgdp'},
            # {'group': 'adj', 'pop': 'han', 'sex': 'XX', 'subset': 'hgdp'},
            # {'group': 'adj', 'pop': 'han', 'sex': 'XY', 'subset': 'hgdp'},
            # {'group': 'adj', 'pop': 'northernhan', 'sex': 'XX', 'subset': 'hgdp'},
            # {'group': 'adj', 'pop': 'northernhan', 'sex': 'XY', 'subset': 'hgdp'}
            subset_ht = subset_ht.annotate_globals(
                freq_meta=[
                    {**x, **{"subset": "hgdp"}}
                    for x in hl.eval(subset_ht.freq_meta)[2:]
                ],
                freq_meta_sample_count=subset_ht.freq_meta_sample_count[2:],
            )
            subset_ht = subset_ht.annotate(
                group_membership=subset_ht.group_membership[2:]
            )
            subset_ht = filter_freq_arrays(
                subset_ht, ["pop"], annotations=["group_membership"]
            )
        else:
            # Remove 'Han' and 'Papuan' pops from group_membership, freq_meta, and
            # freq_meta_sample_count.
            subset_ht = filter_freq_arrays(
                subset_ht,
                {"pop": ["han", "papuan"]},
                keep=False,
                combine_operator="or",
                annotations=["group_membership"],
            )

        # Keep track of which groups should aggregate raw genotypes.
        subset_ht = subset_ht.annotate_globals(
            raw_group=subset_ht.freq_meta.map(lambda x: x.get("group", "NA") == "raw"),
        )
        hts.append(subset_ht)

    return concatenate_annotations(
        hts,
        row_field_names=["group_membership"],
        global_field_names=list(FREQ_GLOBALS) + ["raw_group"],
    )


def compute_an_by_group_membership(
    vds: hl.vds.VariantDataset,
    group_membership_ht: hl.Table,
    variant_filter_ht: hl.Table,
) -> hl.Table:
    """
    Compute the allele number for new variants in the v4.0 release by call stats group membership.

    :param vds: VariantDataset with all v3.1 samples (including non-release).
    :param group_membership_ht: Table with the group membership for each sample. This
        is generated by `get_group_membership_ht_for_an`.
    :param variant_filter_ht: Table with all variants that need AN to be computed.
    :return: Table with the allele number for new variants in the v4.0 release.
    """
    n_samples = group_membership_ht.count()
    n_groups = len(group_membership_ht.group_membership.take(1)[0])

    vmt = vds.variant_data
    rmt = vds.reference_data.select_cols().select_rows()

    # Confirm that all variants in variant_filter_ht are present in the
    # vds.variant_dataset and throw an error if this is not the case since all samples
    # added to v4.0 genomes should be in this vds.
    vht = vmt.rows()
    vht = hl.split_multi(vht)
    if variant_filter_ht.aggregate(
        hl.agg.any(hl.is_missing(vht[variant_filter_ht.key]))
    ):
        raise ValueError(
            "Not all variants in the variant_filter_ht are found in the "
            "vds.variant_dataset!"
        )

    # Filter the VDS variant data and reference data to only keep samples that were in
    # the v3.1 release. We do it this way instead of using hl.vds.filter_samples because
    # we want to keep all variants that were in the v3.1 VDS not only those found in the
    # release but `hl.vds.filter_samples` includes vmt =
    # vmt.filter_rows(hl.agg.count() > 0) by default.
    release_s = vmt.aggregate_cols(
        hl.agg.filter(
            hl.is_defined(group_membership_ht[vmt.s]), hl.agg.collect_as_set(vmt.s)
        ),
        _localize=False,
    )._persist()
    vmt = vmt.filter_cols(release_s.contains(vmt.s))
    rmt = rmt.filter_cols(release_s.contains(rmt.s))
    rmt = rmt.filter_rows(hl.agg.count() > 0)

    # Only keep necessary entries and col annotations for AN calculation from the VDS
    # so that only the needed entries are densified. AD, DP, GT, and GQ are needed for
    # adding the adj annotation to compute AN for all adj stratification groups.
    # sex_karyotype is needed by adjust_sex_ploidy to adjust relevant genotypes on chrX
    # and chrY non-PAR before computing AN.
    vmt = (
        vmt.select_cols(
            sex_karyotype=vmt.meta.sex_imputation.sex_karyotype,
            group_membership=group_membership_ht[vmt.col_key].group_membership,
        )
        .select_entries("LA", "LAD", "DP", "LGT", "GQ")
        .select_rows()
    )
    vmt = hl.experimental.sparse_split_multi(vmt, filter_changed_loci=True)
    vmt = vmt.semi_join_rows(variant_filter_ht)
    vds = hl.vds.VariantDataset(rmt, vmt)

    # NOTE: We need to write and read the VDS here to avoid a bug in split_multi_hts
    # that causes code 137 memory errors in Hail 0.2.122, this should be fixed in
    # subsequent versions of Hail.
    tmp_vds_path = new_temp_file("variants_for_an_split", "vds")

    vds.write(tmp_vds_path)
    vds = hl.vds.read_vds(tmp_vds_path)

    # NOTE: The correct thing to do here is the following commented out code, but in
    # order to be consistent with v3.1 we will add the adj annotation before adjusting
    # for sex ploidy.
    # Densify the VDS, adjust GT sex ploidy, and annotate entries with adj. Adj
    # annotation must happen after sex ploidy adjustment because haploid GTs have
    # different adj filtering criteria.
    mt = hl.vds.to_dense_mt(vds)
    # mt = adjust_sex_ploidy(mt, mt.sex_karyotype, male_str="XY", female_str="XX")
    mt = annotate_adj(mt)
    mt = adjust_sex_ploidy(mt, mt.sex_karyotype, male_str="XY", female_str="XX")
    mt = mt.select_entries("GT", "adj")

    # Convert MT to HT with a row annotation that is an array of all samples entries
    # for that variant.
    ht = mt.localize_entries("entries", "cols")

    # For each stratification group in group_membership, determine the indices of the
    # samples that belong to that group.
    ht = ht.annotate_globals(
        indices_by_group=hl.range(n_groups).map(
            lambda g_i: hl.range(n_samples).filter(
                lambda s_i: ht.cols[s_i].group_membership[g_i]
            )
        ),
        raw_group=group_membership_ht.index_globals().raw_group,
    )

    # Keep only the row annotations needed for the AN calculation.
    # Use ploidy to determine the number of alleles for each sample at each variant and
    # then sum the ploidy of all samples for each group indicated by indices_by_group
    # to get the AN for each stratification group.
    ht = ht.select(
        adj_array=ht.entries.map(lambda e: e.adj),
        ploidy_array=ht.entries.map(lambda e: e.GT.ploidy),
    )
    agg_expr = hl.map(
        lambda s_indices, raw: s_indices.aggregate(
            lambda i: hl.if_else(
                raw,
                hl.agg.sum(ht.ploidy_array[i]),
                hl.agg.filter(ht.adj_array[i], hl.agg.sum(ht.ploidy_array[i])),
            )
        ),
        ht.indices_by_group,
        ht.raw_group,
    )

    # Turn the array of ANs for each group into a call stats struct with AC, AF, AN,
    # and homozygote_count to be easily combined with the call stats of the updated
    # samples.
    agg_expr = agg_expr.map(
        lambda x: hl.struct(
            AC=0, AF=hl.missing(hl.tfloat64), AN=hl.int32(x), homozygote_count=0
        )
    )
    ht = ht.select(freq=agg_expr).drop("cols")
    ht = ht.select_globals(**group_membership_ht.index_globals())

    return ht


def generate_v4_genomes_callstats(
    ht: hl.Table, an_ht: hl.Table, pop_diff_an_ht: hl.Table
) -> hl.Table:
    """
    Generate the call stats for the v4.0 genomes release.

    Merge call stats from the v3.1 release, the v3.1 AN of new v4.0 variants,
    samples with updated population labels, added samples, and removed samples.

    Also, compute filtering allele frequencies ('faf'), 'grpmax', 'gen_anc_faf_max' and
    'InbreedingCoeff' on the merged call stats.

    :param ht: Table returned by `join_release_ht_with_subsets`.
    :param an_ht: Table with the allele number for new variants in the v4.0 release.
    :param pop_diff_an_ht: Table with the allele number for samples with updated pop
        labels and variants in the v4.0 release but not present in the callstats of this subset.
    :return: Table with the updated call stats for the v4.0 genomes release.
    """
    logger.info("Updating AN Table HGDP pop labels and 'oth' -> 'remaining'...")
    an_ht = update_pop_labels(an_ht, POP_MAP)

    logger.info(
        "Merging call stats from new variants of v3.1 release samples, from pop_diff "
        "samples and added samples..."
    )
    global_array = ht.index_globals().global_array
    farrays = [
        ht.ann_array[0].freq,
        an_ht[ht.key].freq,
        pop_diff_an_ht[ht.key].freq,
    ] + [ht.ann_array[i].freq for i in range(1, 3)]
    fmeta = [
        global_array[0].freq_meta,
        an_ht.index_globals().freq_meta,
        pop_diff_an_ht.index_globals().freq_meta,
    ] + [global_array[i].freq_meta for i in range(1, 3)]
    count_arrays = [
        global_array[0].freq_meta_sample_count,
        [0] * hl.eval(hl.len(an_ht.freq_meta)),
        [0] * hl.eval(hl.len(pop_diff_an_ht.freq_meta)),
    ] + [global_array[i].freq_meta_sample_count for i in range(1, 3)]

    ht = ht.select(farrays=farrays, sub_array=ht.ann_array[3].freq)
    ht = ht.checkpoint(new_temp_file("farrays", "ht"))

    # Merge the call stats from the v3.1 release, v3.1 AN of new v4.0 variants,
    # pop_diff, and added samples.
    freq_expr, freq_meta, sample_counts = merge_freq_arrays(
        [ht.farrays[i] for i in range(len(fmeta))],
        fmeta,
        count_arrays={"counts": count_arrays},
    )
    ht = ht.annotate(freq=freq_expr)

    # Needs checkpoint here to avoid hanging.
    ht = ht.checkpoint(new_temp_file("added", "ht"))

    logger.info("Merging call stats from subtracted samples...")
    # Normally `set_negatives_to_zero` should be False, but we set it to True
    # because there was a bug in the v3.1 frequency code where we used the v3.0
    # freq_ht for the hotfix of depletion of homozygous frequency calculation,
    # once a variant is new in v3.1 (i.e. not present in v3.0 freq_ht) and is from
    # samples with heterozygous non-reference genotypes with high allele balances,
    # their GT was set to missing in the v3.1 freq_ht, hence we don't have
    # non-zero callstats for these variants in v3.1 release while non-zero callstats
    # are present the corrected "subtract" freq_ht.
    freq_expr, freq_meta, sample_counts = merge_freq_arrays(
        [ht.freq, ht.sub_array],
        [freq_meta, global_array[3].freq_meta],
        operation="diff",
        set_negatives_to_zero=True,
        count_arrays={
            "counts": [sample_counts["counts"], global_array[3].freq_meta_sample_count]
        },
    )
    ht = ht.select(freq=freq_expr)
    ht = ht.select_globals(
        freq_meta=freq_meta,
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta)),
        freq_meta_sample_count=sample_counts["counts"],
    )
    ht = ht.checkpoint(new_temp_file("merged", "ht"))

    logger.info(
        "Set Y-variant frequency call stats for female-specific metrics to missing "
        "structs and the allele frequency to missing for any call stats entries with "
        "AN == 0..."
    )
    ht = ht.annotate(
        freq=set_female_y_metrics_to_na_expr(ht).map(
            lambda x: x.annotate(AF=hl.or_missing(x.AN > 0, x.AF))
        )
    )

    # Change the 'pop' keys in the freq_meta array to 'gen_anc'.
    freq_meta = hl.literal(
        [{("gen_anc" if k == "pop" else k): m[k] for k in m} for m in freq_meta]
    )
    ht = ht.annotate_globals(freq_meta=freq_meta)

    return ht


def finalize_v4_genomes_callstats(
    ht: hl.Table,
    v3_sites_ht: hl.Table,
    v3_meta_ht: hl.Table,
    updated_meta_ht: hl.Table,
) -> hl.Table:
    """
    Finalize the call stats for the v4.0 genomes release.

    The following is done to create the final v4.0 genomes call stats:
        - Compute filtering allele frequencies ('faf'), 'grpmax', 'gen_anc_faf_max' and
          'inbreeding_coeff' on the v4.0 genomes call stats.
        - Drop downsamplings from the call stats.
        - Get the quality histograms from the v3.1 release.
        - Get the age distribution for the v4.0 genomes release.

    :param ht: Table with the updated call stats for the v4.0 genomes release.
    :param v3_sites_ht: Table for the v3.1.2 release sites.
    :param v3_meta_ht: Table for the v3.1 release metadata.
    :param updated_meta_ht: Table for the updated metadata.
    :return: Table with the finalized call stats for the v4.0 genomes release.
    """
    logger.info("Drop downsampling call stats for v4.0 genomes release...")
    ht = filter_freq_arrays(ht, ["downsampling"], keep=False)

    # Change the 'gen_anc' keys in the freq_meta array to 'pop' to compute faf.
    freq_meta = hl.literal(
        [
            {("pop" if k == "gen_anc" else k): m[k] for k in m}
            for m in hl.eval(ht.freq_meta)
        ]
    )

    # Compute filtering allele frequency (faf), grpmax, and gen_anc_faf_max.
    faf, faf_meta = faf_expr(
        ht.freq, freq_meta, ht.locus, POPS_TO_REMOVE_FOR_POPMAX["v3"]
    )
    grpmax = pop_max_expr(ht.freq, freq_meta, POPS_TO_REMOVE_FOR_POPMAX["v3"])
    grpmax = grpmax.annotate(
        gen_anc=grpmax.pop,
        faf95=faf[
            hl.literal(faf_meta).index(lambda y: y.values() == ["adj", grpmax.pop])
        ].faf95,
    ).drop("pop")

    logger.info(
        "Annotating 'faf', 'grpmax', 'gen_anc_faf_max' and 'InbreedingCoeff'..."
    )
    ht = ht.annotate(
        faf=faf,
        grpmax=grpmax,
        gen_anc_faf_max=gen_anc_faf_max_expr(faf, faf_meta),
        InbreedingCoeff=bi_allelic_site_inbreeding_expr(callstats_expr=ht.freq[1]),
    )

    logger.info(
        "Annotating globals 'freq_meta', 'freq_meta_sample_count', 'faf_meta', "
        "'freq_index_dict' and 'faf_index_dict'..."
    )
    faf_index_dict = make_freq_index_dict_from_meta(hl.literal(faf_meta))

    # Change the 'pop' keys back to 'gen_anc' in the freq_meta and faf_meta arrays.
    freq_meta, faf_meta = [
        hl.literal([{("gen_anc" if k == "pop" else k): m[k] for k in m} for m in meta])
        for meta in [hl.eval(freq_meta), faf_meta]
    ]
    ht = ht.annotate_globals(
        freq_meta=freq_meta,
        freq_index_dict=make_freq_index_dict_from_meta(hl.literal(freq_meta)),
        faf_meta=faf_meta,
        faf_index_dict=faf_index_dict,
    )

    logger.info("Getting quality histograms from the v3.1.2 release...")
    ht = get_histograms(ht, v3_sites_ht)

    logger.info("Getting age distribution for v4.0 genomes release...")
    ht = ht.annotate_globals(
        age_distribution=get_age_distribution(v3_meta_ht, updated_meta_ht),
    )

    return ht


def get_histograms(ht: hl.Table, v3_sites_ht: hl.Table) -> hl.Table:
    """
    Get histograms for the v4.0 genomes release from the v3.1 release.

    .. note:
       We didn't compute the variant histograms for the v4.0 genomes release because
       we didn't densify the VDS for all the variants. The histograms for the new
       variants will be missing, and the histograms for the variants included
       in v3.1 release will be the same as the v3.1 release.

    :param ht: Table with the call stats for the v4.0 genomes release.
    :param v3_sites_ht: Table for the v3.1 release sites.
    :return: Table with the histograms added for the v4.0 genomes release.
    """
    logger.info("Getting histograms for the v4.0 genomes release...")
    v3_sites = v3_sites_ht[ht.key]
    ht = ht.annotate(
        histograms=hl.struct(
            qual_hists=v3_sites.qual_hists,
            raw_qual_hists=v3_sites.raw_qual_hists,
            age_hists=hl.struct(
                age_hist_het=v3_sites.age_hist_het,
                age_hist_hom=v3_sites.age_hist_hom,
            ),
        )
    )

    return ht


def set_downsampling_freq_missing(ht: hl.Table, new_variants_ht: hl.Table) -> hl.Table:
    """
    Set the downsampling call stats for new variants in `ht` to missing.

    :param ht: Table with the call stats for all variants in the v4.0 genomes release.
    :param new_variants_ht: Table with new variants in the v4.0 genomes release.
    :return: Table with the downsampling call stats for new variants set to missing.
    """
    downsampling_indices = hl.range(hl.len(ht.freq_meta)).filter(
        lambda i: ht.freq_meta[i].contains("downsampling")
    )

    ht = ht.annotate(
        freq=hl.if_else(
            hl.is_defined(new_variants_ht[ht.key]),
            hl.enumerate(ht.freq).map(
                lambda i: hl.if_else(
                    downsampling_indices.contains(i[0]),
                    missing_callstats_expr(),
                    i[1],
                )
            ),
            ht.freq,
        )
    )

    return ht


def get_age_distribution(
    v3_meta_ht: hl.Table, subset_updated_meta_ht: hl.Table
) -> hl.expr.StructExpression:
    """
    Get the updated 'age_distribution' annotation for v4.0 genomes.

    :param v3_meta_ht: Table with the v3.1.2 release sample metadata.
    :param subset_updated_meta_ht: Table with the updated sample metadata for the
       HGDP + 1KG subset.
    :return: Expression with the age distribution for samples in `subset_updated_meta`.
    """
    logger.info("Getting age distribution for v4.0 genomes release...")
    v3_release_count = v3_meta_ht.filter(v3_meta_ht.release).count()
    v4_meta_ht = v3_meta_ht.select(
        age=v3_meta_ht.project_meta.age,
        release=hl.if_else(
            v3_meta_ht.subsets.tgp | v3_meta_ht.subsets.hgdp,
            subset_updated_meta_ht[v3_meta_ht.s.replace("v3.1::", "")].gnomad_release,
            v3_meta_ht.release,
        ),
    )
    v4_meta_ht = v4_meta_ht.filter(v4_meta_ht.release)
    logger.info(
        "There will be %s samples in the v4.0 genome release, compared to %s in v3.1 "
        "genome release.",
        v4_meta_ht.count(),
        v3_release_count,
    )

    age_expr = v4_meta_ht.aggregate(hl.agg.hist(v4_meta_ht.age, 30, 80, 10))

    return age_expr


def get_pop_diff_v3_vds_group_membership(
    v3_vds: hl.vds.VariantDataset,
    meta_ht: hl.Table,
) -> hl.Table:
    """
    Get the group membership for samples in pop_diff for the v3.1 VDS.

    :param v3_vds: VDS with the v3.1 release samples.
    :param meta_ht: Table with the updated HGDP + 1KG sample metadata.
    :return: Table with the group membership for samples in pop_diff.
    """
    v3_sites_meta_ht = v3_vds.variant_data.cols()
    v3_sites_meta_ht = v3_sites_meta_ht.filter(v3_sites_meta_ht.meta.release)

    logger.info(
        "Annotating call stats group membership for pop diff release samples..."
    )
    # Need to make sure we are using the names that are in the vds for the release
    # samples.
    # Note: the samples in the pop_diff HT are samples in the to-be-split 'Han' and
    # 'Papuan' populations AND their 'gnomad_release' status hasn't changed.
    v3_sites_meta_rename_ht = v3_sites_meta_ht.key_by(
        s_no_prefix=v3_sites_meta_ht.s.replace("v3.1::", "")
    )
    pop_diff_sample_ht = get_updated_release_samples(meta_ht)[0]
    pop_diff_sample_ht = pop_diff_sample_ht.key_by(
        s=v3_sites_meta_rename_ht[pop_diff_sample_ht.s].s
    )
    pop_diff_group_membership_ht = get_group_membership_ht_for_an(
        pop_diff_sample_ht, only_pop_diff=True
    )
    pop_diff_group_membership_ht = pop_diff_group_membership_ht.checkpoint(
        new_temp_file("group_membership_pop_diff", "ht")
    )

    return pop_diff_group_membership_ht


def patch_v4_genomes_callstats(freq_ht: hl.Table) -> hl.Table:
    """
    Patch the call stats for some inconsistent variant calls found during validity checks.

    14 variants were found to have inconsistent calls between v3.1 and v4.0 genomes
    release, we determined the reason for each of these by manual inspection and
    this function applies the needed patches to these variants.

    .. note::
        This is a temporary fix until we can recompute the call stats for all the
        samples.

    :param freq_ht: Table with the call stats for the v4.0 genomes release.
    :return: Table with the patched call stats for the v4.0 genomes release.
    """
    # For most of the variants, there is a single Han sample that was 0/1 in v3.1, but
    # from the dense MT in v4.0 frequencies it was given a missing genotype. These
    # variants are not in v3, which is why they get a missing genotype in the v4.0
    # freq. We don’t understand why they were 0/1 in v3.1, we would expect them to
    # have been missing. It should be 0/1 if missing v3 frequency had been handled
    # correctly in the fix for high AB hets.
    freq_groups1 = hl.set(
        [
            {"group": "adj", "subset": "hgdp", "gen_anc": "han"},
            {"group": "adj", "subset": "hgdp", "gen_anc": "han", "sex": "XX"},
        ]
    )
    patch1_ht = hl.Table.parallelize(
        [
            {"locus": hl.parse_locus("chr1:111544780"), "alleles": ["C", "T"]},
            {"locus": hl.parse_locus("chr11:119744028"), "alleles": ["G", "C"]},
            {"locus": hl.parse_locus("chr12:82828211"), "alleles": ["CTTTTTTT", "C"]},
            {
                "locus": hl.parse_locus("chr16:87738638"),
                "alleles": ["TCATCCATCCACACACCAGCATCTCATC", "T"],
            },
            {"locus": hl.parse_locus("chr22:22399516"), "alleles": ["G", "A"]},
        ]
    ).key_by("locus", "alleles")
    patch1_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups1.contains(m),
            hl.struct(
                AC=x.AC + 1,
                AF=(x.AC + 1) / (x.AN + 2),
                AN=x.AN + 2,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    freq_groups2 = hl.set(
        [
            {"group": "adj", "subset": "hgdp", "gen_anc": "han"},
            {"group": "adj", "subset": "hgdp", "gen_anc": "han", "sex": "XY"},
        ]
    )
    patch2_ht = hl.Table.parallelize(
        [
            {
                "locus": hl.parse_locus("chr12:10380059"),
                "alleles": ["GTTTTTTTTTTTTTTT", "G"],
            },
            {
                "locus": hl.parse_locus("chr16:73308279"),
                "alleles": ["T", "G"],
            },
        ]
    ).key_by("locus", "alleles")
    patch2_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups2.contains(m),
            hl.struct(
                AC=x.AC + 1,
                AF=(x.AC + 1) / (x.AN + 2),
                AN=x.AN + 2,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    freq_groups3 = hl.set(
        [
            {"group": "adj", "subset": "hgdp", "gen_anc": "northernhan"},
            {"group": "adj", "subset": "hgdp", "gen_anc": "northernhan", "sex": "XY"},
        ]
    )
    patch3_ht = hl.Table.parallelize(
        [{"locus": hl.parse_locus("chr7:14689208"), "alleles": ["T", "TTTTTTTTA"]}]
    ).key_by("locus", "alleles")
    patch3_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups3.contains(m),
            hl.struct(
                AC=x.AC + 1,
                AF=(x.AC + 1) / (x.AN + 2),
                AN=x.AN + 2,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    freq_groups4 = hl.set(
        [
            {"group": "adj", "subset": "hgdp", "gen_anc": "northernhan"},
            {"group": "adj", "subset": "hgdp", "gen_anc": "northernhan", "sex": "XX"},
        ]
    )
    patch4_ht = hl.Table.parallelize(
        [{"locus": hl.parse_locus("chr9:43045892"), "alleles": ["C", "A"]}]
    ).key_by("locus", "alleles")
    patch4_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups4.contains(m),
            hl.struct(
                AC=x.AC + 1,
                AF=(x.AC + 1) / (x.AN + 2),
                AN=x.AN + 2,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    # chr16:89834177 ["G","GGCC"], ["GCC","G"], and
    # ["GCCTGGATAAGCATAGCCCGTGTGAATCTGTGAACCTGCCTGTGCTCACGGTTGGCCGTCGTAGAAGCA", "G"].
    # Code added to the HGDP + 1KG subset to remove alleles not in the subset caused
    # the duplication of a single site that changed position during the minrep. This
    # results in 2 of the Han samples being given NA instead of 0/0 as their GT. To fix
    # this we are adding 4 (2 samples * 2 alleles) to AN.
    freq_groups5 = hl.set(
        [
            {"group": "adj", "subset": "hgdp", "gen_anc": "northernhan"},
            {"group": "adj", "subset": "hgdp", "gen_anc": "northernhan", "sex": "XY"},
        ]
    )
    patch5_ht = hl.Table.parallelize(
        [
            {"locus": hl.parse_locus("chr16:89834177"), "alleles": ["G", "GGCC"]},
            {"locus": hl.parse_locus("chr16:89834177"), "alleles": ["GCC", "G"]},
            {
                "locus": hl.parse_locus("chr16:89834177"),
                "alleles": [
                    "GCCTGGATAAGCATAGCCCGTGTGAATCTGTGAACCTGCCTGTGCTCACGGTTGGCCGTCGTAGAAGCA",
                    "G",
                ],
            },
        ]
    ).key_by("locus", "alleles")
    patch5_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups5.contains(m),
            hl.struct(
                AC=x.AC,
                AF=x.AC / (x.AN + 4),
                AN=x.AN + 4,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    # Frequency groups that need a patch for variant chr9-93596925-ATTTTTTTTTTTTTT-A.
    # There is a single Han sample causing a discrepancy in the AC and AN count between
    # v3.1.2 and v4.0. This sample's GT was set to missing in the v3.1.2 release HT,
    # but in the computation for v4.0 the sample's GT was set to 0/1. This is because in
    # the original code for v3.1 (before the correction to account for het_non_ref
    # samples), the hl.if_else statement in the correction for high AB hets (now in
    # gnomad_qc.v3.utils.hom_alt_depletion_fix) didn't include the _het_non_ref check.
    # Without that check, the rest of the hl.if_else statement would evaluate to
    # missing because of the missing v3.0 frequency. However, when the _het_non_ref
    # check is added, the presence of a False (~_het_non_ref evaluates to False for a
    # het_non_ref) in the hl.if_else statement causes it to be given a 0/1 instead.
    # To correct for this we set Han sample GT to missing to match v3.1.2 value.
    freq_groups6 = hl.set(
        [
            {"group": "adj", "subset": "hgdp", "gen_anc": "han"},
            {"group": "adj", "subset": "hgdp", "gen_anc": "han", "sex": "XX"},
        ]
    )
    patch6_ht = hl.Table.parallelize(
        [
            {
                "locus": hl.parse_locus("chr9:93596925"),
                "alleles": ["ATTTTTTTTTTTTTT", "A"],
            }
        ]
    ).key_by("locus", "alleles")
    patch6_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups6.contains(m),
            hl.struct(
                AC=x.AC - 1,
                AF=(x.AC - 1) / (x.AN - 2),
                AN=x.AN - 2,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    # Frequency groups that need a patch for variant chr13-20656122-T-TAAAAAAAAGAA.
    # This variant has an AC raw of 1 and AC adj of 2 because one of the TGP samples
    # being removed has a raw GT of 0/1 but adj GT of None because it fails
    # adj. In v3.1 it was None for both and therefore shouldn't be subtracted
    # from raw either.
    freq_groups7 = hl.set([{"group": "raw", "subset": "tgp"}, {"group": "raw"}])
    patch7_ht = hl.Table.parallelize(
        [{"locus": hl.parse_locus("chr13:20656122"), "alleles": ["T", "TAAAAAAAAGAA"]}]
    ).key_by("locus", "alleles")
    patch7_expr = hl.map(
        lambda x, m: hl.if_else(
            freq_groups7.contains(m),
            hl.struct(
                AC=x.AC + 1,
                AF=(x.AC + 1) / (x.AN + 2),
                AN=x.AN + 2,
                homozygote_count=x.homozygote_count,
            ),
            x,
        ),
        freq_ht.freq,
        freq_ht.freq_meta,
    )

    # Patch the call stats for the variants that need it.
    freq_ht = freq_ht.annotate(
        freq=(
            hl.case(missing_false=True)
            .when(hl.is_defined(patch1_ht[freq_ht.key]), patch1_expr)
            .when(hl.is_defined(patch2_ht[freq_ht.key]), patch2_expr)
            .when(hl.is_defined(patch3_ht[freq_ht.key]), patch3_expr)
            .when(hl.is_defined(patch4_ht[freq_ht.key]), patch4_expr)
            .when(hl.is_defined(patch5_ht[freq_ht.key]), patch5_expr)
            .when(hl.is_defined(patch6_ht[freq_ht.key]), patch6_expr)
            .when(hl.is_defined(patch7_ht[freq_ht.key]), patch7_expr)
            .default(freq_ht.freq)
        )
    )
    freq_ht = freq_ht.checkpoint(new_temp_file("patched", "ht"))

    return freq_ht


def get_v4_genomes_release_resources(
    test: bool, overwrite: bool
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed to create the gnomAD v4.0 genomes release.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        gnomAD v4.0 genomes release pipeline.
    """
    # Initialize gnomAD v4.0 genomes release pipeline resource collection.
    hgdp_tgp_res = {
        "meta_ht": hgdp_tgp_subset_annotations(sample=True).versions["3.1.2"],
        "dense_mt": hgdp_tgp_subset(dense=True, public=False).versions["3.1.2"],
    }
    v3_sites_ht = release_sites(public=True).versions["3.1.2"]
    v4_genome_release_pipeline = PipelineResourceCollection(
        pipeline_name="gnomad_v4_genomes_release",
        overwrite=overwrite,
        pipeline_resources={
            "Released HGDP + 1KG resources": hgdp_tgp_res,
            "gnomAD v3.1.2 release sites HT": {"v3_sites_ht": v3_sites_ht},
        },
    )

    # Create resource collection for each step of the v4.0 genomes release pipeline.
    get_related_to_nonsubset = PipelineStepResourceCollection(
        "--get-related-to-nonsubset",
        input_resources={
            "v3.1 metadata": {"v3_meta_ht": v3_meta},
            "v3.1 relatedness": {"v3_relatedness_ht": v3_relatedness},
        },
        output_resources={"related_to_nonsubset_ht": hgdp_tgp_related_to_nonsubset},
    )
    hgdp_tgp_v4_exome_duplicates = PipelineStepResourceCollection(
        "--get-duplicated-to-exomes",
        input_resources={
            "v3.1 metadata": {"v3_meta_ht": v3_meta},
            "v4.0 metadata": {"v4_meta_ht": v4_meta},
            "v3.1 and v4.0 relatedness": {"v4_relatedness_ht": v4_relatedness()},
        },
        output_resources={
            "hgdp_tgp_v4_exome_duplicate_ht": hgdp_tgp_duplicated_to_exomes
        },
    )
    update_annotations = PipelineStepResourceCollection(
        "--update-annotations",
        input_resources={
            "Released HGDP + 1KG sample metadata": {"meta_ht": hgdp_tgp_res["meta_ht"]}
        },
        output_resources={"updated_meta_ht": hgdp_tgp_meta_updated},
    )
    get_callstats_for_updated_samples = PipelineStepResourceCollection(
        "--get-callstats-for-updated-samples",
        pipeline_input_steps=[update_annotations],
        add_input_resources={
            "gnomAD v3.1.2 HGDP + 1KG dense MT": {"dense_mt": hgdp_tgp_res["dense_mt"]}
        },
        output_resources={
            "freq_added_ht": hgdp_tgp_updated_callstats(subset="added", test=test),
            "freq_subtracted_ht": hgdp_tgp_updated_callstats(
                subset="subtracted", test=test
            ),
            "freq_pop_diff_ht": hgdp_tgp_updated_callstats(
                subset="pop_diff", test=test
            ),
        },
    )
    join_callstats_for_update = PipelineStepResourceCollection(
        "--join-callstats-for-update",
        pipeline_input_steps=[get_callstats_for_updated_samples],
        add_input_resources={
            "gnomAD v3.1.2 release sites HT": {"sites_ht": v3_sites_ht}
        },
        output_resources={
            "freq_join_ht": hgdp_tgp_updated_callstats(subset="join", test=test),
        },
    )
    compute_an_for_new_variants = PipelineStepResourceCollection(
        "--compute-allele-number-for-new-variants",
        pipeline_input_steps=[update_annotations, join_callstats_for_update],
        output_resources={
            "v3_release_an_ht": hgdp_tgp_updated_callstats(
                subset="v3_release_an", test=test
            ),
        },
    )
    compute_an_for_pop_diff = PipelineStepResourceCollection(
        "--compute-allele-number-for-pop-diff",
        pipeline_input_steps=[update_annotations, join_callstats_for_update],
        output_resources={
            "v3_pop_diff_an_ht": hgdp_tgp_updated_callstats(
                subset="v3_pop_diff_an", test=test
            )
        },
    )
    update_release_callstats = PipelineStepResourceCollection(
        "--update-release-callstats",
        pipeline_input_steps=[
            update_annotations,
            join_callstats_for_update,
            compute_an_for_new_variants,
            compute_an_for_pop_diff,
        ],
        add_input_resources={
            "v3.1 metadata": {"v3_meta_ht": v3_meta},
        },
        output_resources={
            "pre_v4_freq_ht": hgdp_tgp_updated_callstats(
                subset="pre_validity_check", test=test
            )
        },
    )
    apply_patch_to_freq_ht = PipelineStepResourceCollection(
        "--apply-patch-to-freq-ht",
        pipeline_input_steps=[update_annotations, update_release_callstats],
        add_input_resources={
            "v3.1 metadata": {"v3_meta_ht": v3_meta},
        },
        output_resources={
            "v4_freq_ht": v4_get_freq(data_type="genomes", test=test),
        },
    )

    # Add all steps to the gnomAD v4.0 genomes release pipeline resource collection.
    v4_genome_release_pipeline.add_steps(
        {
            "get_related_to_nonsubset": get_related_to_nonsubset,
            "get_hgdp_tgp_v4_exome_duplicates": hgdp_tgp_v4_exome_duplicates,
            "update_annotations": update_annotations,
            "get_callstats_for_updated_samples": get_callstats_for_updated_samples,
            "join_callstats_for_update": join_callstats_for_update,
            "compute_an_for_new_variants": compute_an_for_new_variants,
            "compute_an_for_pop_diff": compute_an_for_pop_diff,
            "update_release_callstats": update_release_callstats,
            "apply_patch_to_freq_ht": apply_patch_to_freq_ht,
        }
    )

    return v4_genome_release_pipeline


def main(args):
    """
    Create the v4.0 genomes release Table.

    Update call stats to include samples from the HGDP + 1KG subset that were
    unintentionally excluded whole populations within the HGDP + 1KG subset that were
    the most genetically unique and had small sample sizes (more specifically: San,
    Mbuti Pygmy, Biaka Pygmy, Bougainville, and Papuan) compared to other populations
    within the gnomAD v3.1 callset.

    In order to avoid re-calculating the call stats for the whole subset / whole
    release, we calculate the call stats for the samples that will be added
    and subtracted, then merge the call stats with the old call stats in the
    release HT.

    Changes compared to the v3 release:
      - Some small updates to samples that are hard filtered.
      - Use a population PC outlier approach to filter the HGDP + 1KG samples instead
        of the sample QC metric outlier filtering approach used on the full dataset
        that caused some samples to be unintentionally excluded.
      - HGDP + 1KG release samples are determined using relatedness (pc_relate) run on
        only samples within the subset as well as relatedness to the rest of the
        release.
      - Some population labels were updated.
      - The `Han` and `Papuan` populations were split into more specific groupings.
    """
    test = args.test_pcsk9 or args.test_x_gene
    gene_on_chrx = args.test_x_gene
    overwrite = args.overwrite

    hl.init(
        log="/create_v4.0_genomes_freq.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-30day",
    )
    v4_genome_release_resources = get_v4_genomes_release_resources(
        test=test, overwrite=overwrite
    )
    v3_hgdp_tgp_meta_ht = v4_genome_release_resources.meta_ht.ht()
    v3_hgdp_tgp_dense_mt = v4_genome_release_resources.dense_mt.mt()
    v3_sites_ht = v4_genome_release_resources.v3_sites_ht.ht()

    v3_vds = None
    if (
        args.compute_allele_number_for_new_variants
        or args.compute_allele_number_for_pop_diff
    ):
        v3_vds = get_gnomad_v3_vds(split=False, samples_meta=True)

    if test:
        v3_hgdp_tgp_dense_mt = filter_to_test(
            v3_hgdp_tgp_dense_mt, gene_on_chrx=gene_on_chrx
        )
        v3_sites_ht = filter_to_test(v3_sites_ht, gene_on_chrx=gene_on_chrx)
        if v3_vds is not None:
            v3_vds = filter_to_test(v3_vds, gene_on_chrx=gene_on_chrx)

    # NOTE: there was a bug in the v3.1 frequency code where we used the v3.0
    #  freq_ht for the hotfix of depletion of homozygous frequency calculation,
    #  once a variant is new in v3.1 (i.e. not present in v3.0 freq_ht) and is from
    #  samples with heterozygous non-reference genotypes with high allele balances,
    #  their GT was set to missing in the v3.1 freq_ht. Here we are using a GT that was
    #  created in the same manner so that we can correctly merge with the existing
    #  v3.1.2 release sites HT.

    logger.info(
        "Annotating GT with sex ploidy adjusted genotypes and adj calculated on "
        "genotypes without sex ploidy adjustment..."
    )
    v3_hgdp_tgp_dense_mt = v3_hgdp_tgp_dense_mt.annotate_entries(
        GT=v3_hgdp_tgp_dense_mt.hom_alt_fix_investigation.sex_ploidy_adjusted_GT,
        adj=v3_hgdp_tgp_dense_mt.hom_alt_fix_investigation.unadjusted_adj,
    )

    # Temporary hotfix for depletion of homozygous alternate genotypes
    logger.info(
        "Setting het genotypes at sites with >1% AF (using v3.0 frequencies) and >"
        " 0.9 AB to homalt..."
    )
    # Load v3.0 allele frequencies to avoid an extra frequency calculation
    # NOTE: Using previous callset AF works for small incremental changes to a callset, but we will need to revisit for large increments # noqa
    freq_ht = release_sites(public=True).versions["3.0"].ht().select("freq")
    v3_hgdp_tgp_dense_mt = v3_hgdp_tgp_dense_mt.annotate_entries(
        GT=hom_alt_depletion_fix(
            v3_hgdp_tgp_dense_mt.GT,
            het_non_ref_expr=v3_hgdp_tgp_dense_mt._het_non_ref,
            af_expr=freq_ht[v3_hgdp_tgp_dense_mt.row_key].freq[0].AF,
            ab_expr=v3_hgdp_tgp_dense_mt.AD[1] / v3_hgdp_tgp_dense_mt.DP,
            use_v3_1_correction=True,
        ),
    )

    if args.get_related_to_nonsubset:
        res = v4_genome_release_resources.get_related_to_nonsubset
        res.check_resource_existence()
        ht = get_hgdp_tgp_related_to_nonsubset(
            res.v3_meta_ht.ht(), res.v3_relatedness_ht.ht()
        )
        ht.write(res.related_to_nonsubset_ht.path, overwrite=overwrite)

    if args.get_hgdp_tgp_v4_exome_duplicates:
        res = v4_genome_release_resources.get_hgdp_tgp_v4_exome_duplicates
        res.check_resource_existence()
        ht = get_hgdp_tgp_v4_exome_duplicates(
            res.v3_meta_ht.ht(), res.v4_meta_ht.ht(), res.v4_relatedness_ht.ht()
        )
        ht.write(res.hgdp_tgp_v4_exome_duplicate_ht.path, overwrite=overwrite)

    if args.update_annotations:
        res = v4_genome_release_resources.update_annotations
        res.check_resource_existence()
        logger.info("Adding updated sample QC annotations to meta HT...")
        add_updated_sample_qc_annotations(v3_hgdp_tgp_meta_ht).write(
            res.updated_meta_ht.path, overwrite=overwrite
        )

    if args.get_callstats_for_updated_samples:
        res = v4_genome_release_resources.get_callstats_for_updated_samples
        res.check_resource_existence()
        logger.info(
            "Calculating and concatenating call stats for samples to be added and "
            "samples to be subtracted from the v3.1.2 release sites HT..."
        )
        mt = v3_hgdp_tgp_dense_mt.select_globals()
        filtered_hts = get_updated_release_samples(res.updated_meta_ht.ht())
        for ht, name in zip(filtered_hts, JOIN_FREQS[1:]):
            freq_ht = get_hgdp_tgp_callstats_for_selected_samples(
                mt, ht, compute_freq_all=False if name == "pop_diff" else True
            )
            if name == "pop_diff":
                freq_ht = filter_freq_arrays(freq_ht, ["pop"])

            freq_ht.write(getattr(res, f"freq_{name}_ht").path, overwrite=overwrite)

    if args.join_callstats_for_update:
        res = v4_genome_release_resources.join_callstats_for_update
        res.check_resource_existence()
        logger.info("Joining all call stats HTs with the release sites HT...")
        join_release_ht_with_subsets(
            v3_sites_ht,
            {n: getattr(res, f"freq_{n}_ht").ht() for n in JOIN_FREQS[1:]},
        ).write(res.freq_join_ht.path, overwrite=overwrite)

    if args.compute_allele_number_for_new_variants:
        res = v4_genome_release_resources.compute_an_for_new_variants
        res.check_resource_existence()
        logger.info("Determine variants in the release HT that need a recomputed AN...")
        ht = res.freq_join_ht.ht()
        # First freq in ann_array is the release sites freq and the third is the added.
        ht = ht.filter(
            hl.is_missing(ht.ann_array[0])
            & hl.any(ht.ann_array.freq[2].map(lambda x: x.AC > 0))
        )
        ht = ht.checkpoint(new_temp_file("variants_for_an", "ht"))
        logger.info(
            "There are %i variants with an AC > 0 in the added call stats HTs, but "
            "missing in the release HT...",
            ht.count(),
        )

        v3_sites_meta_ht = v3_vds.variant_data.cols()
        v3_sites_meta_ht = v3_sites_meta_ht.filter(v3_sites_meta_ht.meta.release)

        logger.info(
            "Annotating call stats group membership for each v3.1 release sample..."
        )
        group_membership_ht = get_group_membership_ht_for_an(v3_sites_meta_ht)
        group_membership_ht = group_membership_ht.checkpoint(
            new_temp_file("group_membership_all", "ht")
        )

        logger.info(
            "Computing the AN HT of v3.1 samples for new v4.0 genomes variants..."
        )
        ht = compute_an_by_group_membership(v3_vds, group_membership_ht, ht)
        ht.write(res.v3_release_an_ht.path, overwrite=overwrite)

    if args.compute_allele_number_for_pop_diff:
        res = v4_genome_release_resources.compute_an_for_pop_diff
        res.check_resource_existence()
        pop_diff_group_membership_ht = get_pop_diff_v3_vds_group_membership(
            v3_vds, res.updated_meta_ht.ht()
        ).checkpoint(new_temp_file("group_membership_pop_diff", "ht"))

        logger.info(
            "Computing the AN HT of v3.1 pop diff samples for all v4.0 genomes"
            " variants. Freq meta: %s",
            hl.eval(pop_diff_group_membership_ht.freq_meta),
        )
        # To avoid double counting, we need to filter out the variants that are in
        # v4.0 genomes but not already present in pop_diff samples to get the AN.
        join_ht = res.freq_join_ht.ht()
        join_ht = join_ht.filter(hl.is_missing(join_ht.ann_array[1].freq)).select()
        ht = compute_an_by_group_membership(
            v3_vds, pop_diff_group_membership_ht, join_ht
        )
        ht.write(res.v3_pop_diff_an_ht.path, overwrite=overwrite)

    if args.update_release_callstats:
        res = v4_genome_release_resources.update_release_callstats
        res.check_resource_existence()

        logger.info("Merging all call stats HTs for final v4.0 genomes call stats...")
        ht = generate_v4_genomes_callstats(
            res.freq_join_ht.ht(), res.v3_release_an_ht.ht(), res.v3_pop_diff_an_ht.ht()
        )
        ht = finalize_v4_genomes_callstats(
            ht, v3_sites_ht, res.v3_meta_ht.ht(), res.updated_meta_ht.ht()
        )
        ht.write(res.pre_v4_freq_ht.path, overwrite=overwrite)

    if args.apply_patch_to_freq_ht:
        res = v4_genome_release_resources.apply_patch_to_freq_ht
        res.check_resource_existence()
        ht = patch_v4_genomes_callstats(res.pre_v4_freq_ht.ht())
        ht = finalize_v4_genomes_callstats(
            ht, v3_sites_ht, res.v3_meta_ht.ht(), res.updated_meta_ht.ht()
        )
        ht.write(res.v4_freq_ht.path, overwrite=overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser(
        description="This script creates the v4.0 genomes release HT."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )

    test_args = parser.add_mutually_exclusive_group()
    test_args.add_argument(
        "--test-pcsk9",
        help="Test on a subset of variants in PCSK9 gene.",
        action="store_true",
    )
    test_args.add_argument(
        "--test-x-gene",
        help="Test on a subset of variants in TRPC5 on chrX.",
        action="store_true",
    )

    parser.add_argument(
        "--get-related-to-nonsubset",
        help="Get the relatedness to nonsubset samples.",
        action="store_true",
    )
    parser.add_argument(
        "--get-hgdp-tgp-v4-exome-duplicates",
        help="Get the duplicated samples to exomes.",
        action="store_true",
    )
    parser.add_argument(
        "--update-annotations",
        help="Update HGDP + 1KG sample QC annotations.",
        action="store_true",
    )
    parser.add_argument(
        "--get-callstats-for-updated-samples",
        help=(
            "Get the call stats for the genomes with an updated release status for v4.0"
            " compared to v3.1."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--join-callstats-for-update",
        help=(
            "Join the call stats Tables for the updated samples with the release "
            "call stats Table."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-allele-number-for-new-variants",
        help=(
            "Compute the allele number for variants that are in the updated samples "
            "but not in the release HT."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--compute-allele-number-for-pop-diff",
        help=(
            "Compute the allele number of all v4.0 variants on the samples that have "
            "different pops in v4.0."
        ),
        action="store_true",
    )

    parser.add_argument(
        "--update-release-callstats",
        help=(
            "Update the release call stats by merging the call stats for the updated"
            " samples."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--apply-patch-to-freq-ht",
        help="Apply a patch to the final freq HT.",
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
