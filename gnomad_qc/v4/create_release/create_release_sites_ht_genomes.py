"""Script to create release sites HT for v4.0 genomes."""

import argparse
import logging
from typing import Any, Dict, List, Optional, Tuple, Union

import hail as hl
from gnomad.resources.grch38.gnomad import SUBSETS
from gnomad.sample_qc.sex import adjust_sex_ploidy
from gnomad.utils.annotations import (
    annotate_adj,
    annotate_freq,
    generate_freq_group_membership_array,
    update_structured_annotations,
)
from gnomad.utils.filtering import filter_arrays_by_meta
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
    "Create v4.0 genomes release sites HT with updated HGDP/TGP "
    "metadata and new annotations"
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


def remove_missing_vep_fields(vep_expr: hl.StructExpression) -> hl.StructExpression:
    """
    Remove fields from VEP 105 annotations that have been excluded in past releases or are missing in all rows.

    :param vep_expr: StructExpression containing VEP 105 annotations.
    :return: StructExpression containing VEP 105 annotations with missing fields removed.
    """
    vep_expr = vep_expr.drop("colocated_variants", "context")

    vep_expr = vep_expr.annotate(
        transcript_consequences=vep_expr.transcript_consequences.map(
            lambda x: x.drop("minimised", "swissprot", "trembl", "uniparc")
        )
    )

    for consequence in [
        "intergenic_consequences",
        "motif_feature_consequences",
        "regulatory_feature_consequences",
    ]:
        vep_expr = vep_expr.annotate(
            **{consequence: vep_expr[consequence].map(lambda x: x.drop("minimised"))}
        )
    return vep_expr


def filter_to_test(
    t: Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset],
    gene_on_chrx: bool = False,
    partitions: Optional[List[int]] = None,
) -> Union[hl.Table, hl.MatrixTable, hl.vds.VariantDataset]:
    """
    Filter to a region in PCSK9 or PLCXD1 in Table or MatrixTable for testing purposes.

    :param t: Table or MatrixTable to filter.
    :param gene_on_chrx: Whether to filter to PLCXD1, instead of PCSK9, for testing chrX.
    :param partitions: Optional list of partitions to filter to before applying the
        filter to PCSK9.
    :return: Table or MatrixTable filtered to a region in PCSK9 or PLCXD1.
    """
    if gene_on_chrx:
        logger.info("Filtering to PLCXD1 on chrX in MT for testing purposes...")
        test_locus = "chrX:285000-295000"
    else:
        logger.info("Filtering to 10kb in PCSK9 in MT for testing purposes...")
        test_locus = "chr1:55039447-55064852"

    test_interval = [hl.parse_locus_interval(test_locus, reference_genome="GRCh38")]

    if isinstance(t, hl.vds.VariantDataset):
        if partitions is not None:
            t = hl.vds.VariantDataset(
                t.reference_data._filter_partitions(partitions),
                t.variant_data._filter_partitions(partitions),
            )
        return hl.vds.filter_intervals(t, test_interval, split_reference_blocks=True)
    else:
        if partitions is not None:
            t = t._filter_partitions(partitions)
        return hl.filter_intervals(t, test_interval)


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


def get_group_membership_ht_for_an(ht: hl.Table) -> hl.Table:
    """
    Generate a Table with a 'group_membership' array for each sample indicating whether the sample belongs to specific stratification groups.

    `ht` must have the following annotations:
        - 'meta.population_inference.pop': population label.
        - 'meta.sex_imputation.sex_karyotype`: sex label.
        - 'meta.project_meta.project_subpop': subpopulation label.
        - 'meta.subsets': dictionary with subset labels as keys and boolean values.

    :param ht: Table with the sample metadata.
    :return: Table with the group membership for each sample to be used for computing
        allele number (AN) per group.
    """
    pop_expr = ht.meta.population_inference.pop
    sex_expr = ht.meta.sex_imputation.sex_karyotype
    subpop_expr = ht.meta.project_meta.project_subpop

    hts = []
    for subset in ["all"] + SUBSETS:
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
        if subset != "all":
            subset_expr = hl.or_missing(ht.meta.subsets[subset], subset)
            strata_expr = [{"subset": subset_expr}] + [
                {"subset": subset_expr, **x} for x in strata_expr
            ]

        subset_ht = generate_freq_group_membership_array(
            ht, strata_expr, remove_zero_sample_groups=True
        )
        subset_globals = subset_ht.index_globals()
        if subset != "all":
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

    :param vds: VariantDataset with all v3.1 samples (including non-release) and
        multi-allelics split.
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
    if variant_filter_ht.aggregate(
        hl.agg.any(hl.is_missing(vmt.rows()[variant_filter_ht.key]))
    ):
        raise ValueError(
            "Not all variants in the variant_filter_ht are found in the "
            "vds.variant_dataset!"
        )

    # Filter the VDS variant data and reference data to only keep samples that were in
    # the v3.1 release. We do it this way instead of using hl.vds.filter_samples because
    # we want to keep all variants that were in the v3.1 VDS not only those found in the
    # release and filter_samples includes vmt = vmt.filter_rows(hl.agg.count() > 0).
    release_s = vmt.aggregate_cols(
        hl.agg.filter(vmt.meta.release, hl.agg.collect_as_set(vmt.s)), _localize=False
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
        .select_entries("AD", "DP", "GT", "GQ")
        .select_rows()
    )
    vmt = vmt.semi_join_rows(variant_filter_ht)
    vds = hl.vds.VariantDataset(rmt, vmt)

    # Densify the VDS, adjust GT sex ploidy, and annotate entires with adj. Adj
    # annotation must happen after sex ploidy adjustment because haploid GTs have
    # different adj filtering criteria.
    mt = hl.vds.to_dense_mt(vds)
    mt = adjust_sex_ploidy(mt, mt.sex_karyotype, male_str="XY", female_str="XX")
    mt = annotate_adj(mt)
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
        "dense_mt": hgdp_tgp_subset(dense=True, public=True).versions["3.1.2"],
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
        "--compute-an-for-new-variants",
        pipeline_input_steps=[join_callstats_for_update],
        output_resources={
            "v3_release_an_ht": hgdp_tgp_updated_callstats(
                subset="v3_release_an", test=test
            ),
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
        log="/create_release_v4_genomes.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    v4_genome_release_resources = get_v4_genomes_release_resources(
        test=test, overwrite=overwrite
    )
    v3_hgdp_tgp_meta_ht = v4_genome_release_resources.meta_ht.ht()
    v3_hgdp_tgp_dense_mt = v4_genome_release_resources.dense_mt.mt()
    v3_sites_ht = v4_genome_release_resources.v3_sites_ht.ht()

    v3_vds = None
    if args.compute_allele_number_for_new_variants:
        v3_vds = get_gnomad_v3_vds(split=True, samples_meta=True)

    if test:
        v3_hgdp_tgp_dense_mt = filter_to_test(
            v3_hgdp_tgp_dense_mt, gene_on_chrx=gene_on_chrx
        )
        v3_sites_ht = filter_to_test(v3_sites_ht, gene_on_chrx=gene_on_chrx)
        if v3_vds is not None:
            v3_vds = filter_to_test(
                v3_vds,
                gene_on_chrx=gene_on_chrx,
                partitions=[2183, 2184, 108583],
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
        logger.info("Joining all callstats HTs with the release sites HT...")
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

        logger.info("Annotating call stats group membership for each v3.1 release sample...")
        v3_sites_meta_ht = v3_vds.variant_data.cols()
        v3_sites_meta_ht = v3_sites_meta_ht.filter(v3_sites_meta_ht.meta.release)
        group_membership_ht = get_group_membership_ht_for_an(v3_sites_meta_ht)
        group_membership_ht = group_membership_ht.checkpoint(
            new_temp_file("group_membership_all", "ht")
        )

        logger.info(
            "Computing the AN HT of v3.1 samples for new v4.0 genomes variants..."
        )
        ht = compute_an_by_group_membership(v3_vds, group_membership_ht, ht)
        ht.write(res.v3_release_an_ht.path, overwrite=overwrite)


if __name__ == "__main__":
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
        help="Test on a subset of variants in PLCXD1 on chrX.",
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
            "Join the callstats Tables for the updated samples with the release "
            "callstats Table."
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

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
