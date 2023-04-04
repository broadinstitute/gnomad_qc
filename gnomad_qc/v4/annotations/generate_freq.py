"""Script to generate the frequency data annotations across v4 exomes."""  # TODO: Expand on script description once nearing completion.
import argparse
import logging
from typing import (  # TODO: REMOVE WITH ANNOTATE_FREQ
    Any,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    Union,
)

import hail as hl
from gnomad.resources.grch38.gnomad import (
    SUBSETS,  # TODO: subsets will be changed to UKB, non-UKB,
)
from gnomad.resources.grch38.gnomad import (
    DOWNSAMPLINGS,
    POPS,
    POPS_TO_REMOVE_FOR_POPMAX,
)
from gnomad.sample_qc.sex import adjusted_sex_ploidy_expr
from gnomad.utils.annotations import (
    age_hists_expr,
    annotate_adj,
    annotate_freq,
    bi_allelic_site_inbreeding_expr,
    faf_expr,
    get_adj_expr,
    pop_max_expr,
    qual_hist_expr,
    set_female_y_metrics_to_na_expr,
)
from gnomad.utils.file_utils import file_exists
from gnomad.utils.release import make_faf_index_dict, make_freq_index_dict
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import (  # TODO: Need to generate this resource for raw callstats
    get_freq,
)
from gnomad_qc.v4.resources.basics import (
    get_checkpoint_path,
    get_gnomad_v4_vds,
    qc_temp_prefix,
)
from gnomad_qc.v4.resources.meta import meta
from gnomad_qc.v4.resources.release import release_sites
from gnomad_qc.v4.resources.variant_qc import (  # TODO: Confirm we want this in frequency calculations
    SYNDIP,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("gnomAD_frequency")
logger.setLevel(logging.INFO)


# UNTIL PR GOES IN
def annotate_freq(
    mt: hl.MatrixTable,
    sex_expr: Optional[hl.expr.StringExpression] = None,
    pop_expr: Optional[hl.expr.StringExpression] = None,
    subpop_expr: Optional[hl.expr.StringExpression] = None,
    additional_strata_expr: Optional[Dict[str, hl.expr.StringExpression]] = None,
    additional_strata_grouping_expr: Optional[
        Dict[str, hl.expr.StringExpression]
    ] = None,
    downsamplings: Optional[List[int]] = None,
) -> hl.MatrixTable:
    """
    Annotate `mt` with stratified allele frequencies.

    The output Matrix table will include:
        - row annotation `freq` containing the stratified allele frequencies
        - global annotation `freq_meta` with metadata
        - global annotation `freq_sample_count` with sample count information

    .. note::

        Currently this only supports bi-allelic sites.
        The input `mt` needs to have the following entry fields:
        - GT: a CallExpression containing the genotype
        - adj: a BooleanExpression containing whether the genotype is of high quality or not.
        All expressions arguments need to be expression on the input `mt`.

    .. rubric:: `freq` row annotation

    The `freq` row annotation is an Array of Struct, with each Struct containing the following fields:

        - AC: int32
        - AF: float64
        - AN: int32
        - homozygote_count: int32

    Each element of the array corresponds to a stratification of the data, and the metadata about these annotations is
    stored in the globals.

    .. rubric:: Global `freq_meta` metadata annotation

    The global annotation `freq_meta` is added to the input `mt`. It is a list of dict.
    Each element of the list contains metadata on a frequency stratification and the index in the list corresponds
    to the index of that frequency stratification in the `freq` row annotation.

    .. rubric:: Global `freq_sample_count` annotation

    The global annotation `freq_sample_count` is added to the input `mt`. This is a sample count per sample grouping
    defined in the `freq_meta` global annotation.

    .. rubric:: The `downsamplings` parameter

    If the `downsamplings` parameter is used, frequencies will be computed for all samples and by population
    (if `pop_expr` is specified) by downsampling the number of samples without replacement to each of the numbers specified in the
    `downsamplings` array, provided that there are enough samples in the dataset.
    In addition, if `pop_expr` is specified, a downsampling to each of the exact number of samples present in each population is added.
    Note that samples are randomly sampled only once, meaning that the lower downsamplings are subsets of the higher ones.

    :param mt: Input MatrixTable
    :param sex_expr: When specified, frequencies are stratified by sex. If `pop_expr` is also specified, then a pop/sex stratifiction is added.
    :param pop_expr: When specified, frequencies are stratified by population. If `sex_expr` is also specified, then a pop/sex stratifiction is added.
    :param subpop_expr: When specified, frequencies are stratified by sub-continental population. Note that `pop_expr` is required as well when using this option.
    :param additional_strata_expr: When specified, frequencies are stratified by the given additional strata found in the dict. This can e.g. be used to stratify by platform.
    :param additional_strata_grouping_expr: When specified, frequencies are further stratified by groups within the additional_strata_expr. This can e.g. be used to stratify by platform-population.
    :param downsamplings: When specified, frequencies are computed by downsampling the data to the number of samples given in the list. Note that if `pop_expr` is specified, downsamplings by population is also computed.
    :return: MatrixTable with `freq` annotation
    """
    if subpop_expr is not None and pop_expr is None:
        raise NotImplementedError(
            "annotate_freq requires pop_expr when using subpop_expr"
        )

    if additional_strata_grouping_expr is not None and additional_strata_expr is None:
        raise NotImplementedError(
            "annotate_freq requires additional_strata_expr when using"
            " additional_strata_grouping_expr"
        )

    if additional_strata_expr is None:
        additional_strata_expr = {}

    _freq_meta_expr = hl.struct(**additional_strata_expr)
    if additional_strata_grouping_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(**additional_strata_grouping_expr)
    if sex_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(sex=sex_expr)
    if pop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(pop=pop_expr)
    if subpop_expr is not None:
        _freq_meta_expr = _freq_meta_expr.annotate(subpop=subpop_expr)

    # Annotate cols with provided cuts
    mt = mt.annotate_cols(_freq_meta=_freq_meta_expr)

    if additional_strata_grouping_expr is None:
        additional_strata_grouping_expr = {}

    # Get counters for sex, pop and if set subpop and additional strata
    cut_dict = {
        cut: hl.agg.filter(
            hl.is_defined(mt._freq_meta[cut]), hl.agg.counter(mt._freq_meta[cut])
        )
        for cut in mt._freq_meta
        if cut != "subpop"
    }
    if "subpop" in mt._freq_meta:
        cut_dict["subpop"] = hl.agg.filter(
            hl.is_defined(mt._freq_meta.pop) & hl.is_defined(mt._freq_meta.subpop),
            hl.agg.counter(
                hl.struct(subpop=mt._freq_meta.subpop, pop=mt._freq_meta.pop)
            ),
        )

    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    sample_group_filters = []

    # Create downsamplings if needed
    if downsamplings is not None:
        # Add exact pop size downsampling if pops were provided
        if cut_data.get("pop"):
            downsamplings = list(
                set(downsamplings + list(cut_data.get("pop").values()))
            )  # Add the pops values if not in yet
            downsamplings = sorted(
                [x for x in downsamplings if x <= sum(cut_data.get("pop").values())]
            )
        logger.info("Found %d downsamplings: %s", len(downsamplings), downsamplings)

        # Shuffle the samples, then create a global index for downsampling
        # And a pop-index if pops were provided
        downsampling_ht = mt.cols()
        downsampling_ht = downsampling_ht.annotate(r=hl.rand_unif(0, 1))
        downsampling_ht = downsampling_ht.order_by(downsampling_ht.r)
        scan_expr = {"global_idx": hl.scan.count()}
        if cut_data.get("pop"):
            scan_expr["pop_idx"] = hl.scan.counter(downsampling_ht._freq_meta.pop).get(
                downsampling_ht._freq_meta.pop, 0
            )
        downsampling_ht = downsampling_ht.annotate(**scan_expr)
        downsampling_ht = downsampling_ht.key_by("s").select(*scan_expr)
        mt = mt.annotate_cols(downsampling=downsampling_ht[mt.s])
        mt = mt.annotate_globals(downsamplings=downsamplings)

        # Create downsampled sample groups
        sample_group_filters.extend(
            [
                (
                    {"downsampling": str(ds), "pop": "global"},
                    mt.downsampling.global_idx < ds,
                )
                for ds in downsamplings
            ]
        )
        if cut_data.get("pop"):
            sample_group_filters.extend(
                [
                    (
                        {"downsampling": str(ds), "pop": pop},
                        (mt.downsampling.pop_idx < ds) & (mt._freq_meta.pop == pop),
                    )
                    for ds in downsamplings
                    for pop, pop_count in cut_data.get("pop", {}).items()
                    if ds <= pop_count
                ]
            )

    # Add all desired strata, starting with the full set and ending with
    # downsamplings (if any)
    sample_group_filters = (
        [({}, True)]
        + [({"pop": pop}, mt._freq_meta.pop == pop) for pop in cut_data.get("pop", {})]
        + [({"sex": sex}, mt._freq_meta.sex == sex) for sex in cut_data.get("sex", {})]
        + [
            (
                {"pop": pop, "sex": sex},
                (mt._freq_meta.sex == sex) & (mt._freq_meta.pop == pop),
            )
            for sex in cut_data.get("sex", {})
            for pop in cut_data.get("pop", {})
        ]
        + [
            (
                {"subpop": subpop.subpop, "pop": subpop.pop},
                (mt._freq_meta.pop == subpop.pop)
                & (mt._freq_meta.subpop == subpop.subpop),
            )
            for subpop in cut_data.get("subpop", {})
        ]
        + [
            ({strata: str(s_value)}, mt._freq_meta[strata] == s_value)
            for strata in additional_strata_expr
            for s_value in cut_data.get(strata, {})
        ]
        + sample_group_filters
    )

    # Add additional groupings to strata, e.g. strata-pop, strata-sex, strata-pop-sex
    if additional_strata_grouping_expr is not None:
        for strata in additional_strata_grouping_expr:
            if strata not in cut_data.keys():
                raise KeyError(
                    "%s is not an existing annotation and thus cannot be combined with"
                    " additional strata"
                )
        sample_group_filters.extend(
            [
                (
                    {strata: str(s_value), add_strata: str(as_value)},
                    (mt._freq_meta[strata] == s_value)
                    & (mt._freq_meta[add_strata] == as_value),
                )
                for strata in additional_strata_expr
                for s_value in cut_data.get(strata, {})
                for add_strata in additional_strata_grouping_expr
                for as_value in cut_data.get(add_strata, {})
            ]
        )

    freq_sample_count = mt.aggregate_cols(
        [hl.agg.count_where(x[1]) for x in sample_group_filters]
    )

    # Annotate columns with group_membership
    mt = mt.annotate_cols(group_membership=[x[1] for x in sample_group_filters])

    # Create and annotate global expression with meta and sample count information
    freq_meta_expr = [
        dict(**sample_group[0], group="adj") for sample_group in sample_group_filters
    ]
    freq_meta_expr.insert(1, {"group": "raw"})
    freq_sample_count.insert(1, freq_sample_count[0])
    mt = mt.annotate_globals(
        freq_meta=freq_meta_expr,
        freq_sample_count=freq_sample_count,
    )

    # Create frequency expression array from the sample groups
    # Adding sample_group_filters_range_array to reduce memory usage in this array_agg
    mt = mt.annotate_rows(
        sample_group_filters_range_array=hl.range(len(sample_group_filters))
    )
    freq_expr = hl.agg.array_agg(
        lambda i: hl.agg.filter(
            mt.group_membership[i] & mt.adj, hl.agg.call_stats(mt.GT, mt.alleles)
        ),
        mt.sample_group_filters_range_array,
    )

    # Insert raw as the second element of the array
    freq_expr = (
        freq_expr[:1]
        .extend([hl.agg.call_stats(mt.GT, mt.alleles)])
        .extend(freq_expr[1:])
    )

    # Select non-ref allele (assumes bi-allelic)
    freq_expr = freq_expr.map(
        lambda cs: cs.annotate(
            AC=cs.AC[1],
            AF=cs.AF[
                1
            ],  # TODO This is NA in case AC and AN are 0 -- should we set it to 0?
            homozygote_count=cs.homozygote_count[1],
        )
    )

    # Return MT with freq row annotation
    return mt.annotate_rows(freq=freq_expr).drop("_freq_meta")


def annotate_non_ref_het(vds: hl.vds.VariantDataset) -> hl.vds.VariantDataset:
    """
    Annotate non-ref heterozygous calls to prevent incorrect call adjustment post-split.

    :param vds: Hail VDS to annotate non-ref het onto variant data.
    :return: Hail VDS with _non_ref_het annotation.
    """
    vds.variant_data = vds.variant_data.annotate_entries(
        _het_non_ref=vds.variant_data.LGT.is_het_non_ref()
    )
    return vds


def annotate_high_ab_hets_by_group_membership(
    mt: hl.MatrixTable,
    gatk_expr: hl.expr.StringExpression,
    ab_cutoff: hl.float = 0.9,
    gatk_versions_to_fix: hl.expr.SetExpression = hl.set(["4.0.10.1", "4.1.0.0"]),
) -> hl.MatrixTable:
    """
    Annotate high AB hets by group membership.

    :param mt: MatrixTable to annotate high AB hets onto.
    :param gatk_expr: GATK version expression.
    :param ab_cutoff: Allele balance cutoff for hom alt depletion fix. Defaults to 0.9
    :param gatk_versions_to_fix: GATK versions that need the hom alt depletion fix. Defaults to hl.set(["4.0.10.1", "4.1.0.0"])
    :return: _description_
    """
    logger.info("Annotating number of high AB het sites in each freq group...")
    mt = mt.annotate_rows(
        high_ab_hets_by_group_membership=hl.agg.array_agg(
            lambda i: hl.agg.filter(
                mt.group_membership[i]
                & needs_high_ab_het_fix_expr(
                    mt, gatk_expr, ab_cutoff, gatk_versions_to_fix
                ),
                hl.agg.count(),
            ),
            hl.range(hl.len(mt.group_membership)),
        )
    )
    return mt


def needs_high_ab_het_fix_expr(
    mt: hl.MatrixTable,
    gatk_expr: hl.expr.StringExpression,
    ab_cutoff: hl.float = 0.9,
    gatk_versions_to_fix: hl.expr.SetExpression = hl.set(["4.0.10.1", "4.1.0.0"]),
) -> hl.MatrixTable:
    """
    Annotate the MatrixTable rows with the count of high AB het genotypes per frequency sample grouping.

    This arrary allows the AC to be corrected with the homalt depletion fix.

    :param mt: Hail MatrixTable to annotate.
    :param gatk_expr: Expression indicating which version of GATK a sample's gVCF was called with.
    :param ab_cutoff: Alelle balance threshold at which the GT changes to homalt. Default is 0.9.
    :param gatk_versions_to_fix: Set of GATK versions impacted by hom alt depletion. Default is {"4.0.10.1", "4.1.0.0"}.
    :return mt:
    """
    return (
        (mt.AD[1] / mt.DP > ab_cutoff)
        & mt.adj
        & mt.GT.is_het_ref()
        & (gatk_versions_to_fix.contains(gatk_expr))
        # & ~mt._het_non_ref  # Skip adjusting genotypes if sample originally had a het nonref genotype
    )


def annotate_gatk_version(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Annotate MatrixTable's samples with HaplotypeCaller's GATK version. GATK hom alt depletion fix has been in since GATK version 4.1.4.1.

    :param mt: Hail MatrixTable to annotate.
    :return mt:
    """
    meta_ht = meta.ht()
    gatk_ht = hl.read_table(
        "gs://gnomad-mwilson/v4/homalt_depletion/v4_sample_htc_versions.ht"
    ).key_by(
        "s"
    )  # TODO: Discuss if this should be a resource or annotation on metadata, misssing GATK version info for UKB samples
    mt = mt.annotate_cols(
        gatk_version=hl.case()
        .when(
            hl.is_defined(gatk_ht[mt.col_key].gatk_version),
            gatk_ht[mt.col_key].gatk_version,
        )
        .when(hl.is_defined(meta_ht[mt.col_key].project_meta.ukb_meta.ukb_batch), "ukb")
        .or_missing()
    )
    return mt


def hom_alt_depletion_fix(ht: hl.Table, af_threshold: float = 0.01) -> hl.Table:
    """
    Correct frequencies at sites with an AF greater than the af_threshold.

    :param ht: Hail Table containing freq and high_ab_het annotations.
    :param af_threshold: AF threshold at which to correct frequency. Default is 0.01.
    :return: Hail Table with adjusted frequencies.
    """
    # TODO: Adjust freq AC array using high_ab_het array
    return ht


def main(args):  # noqa: D103
    # TODO: Determine if splitting subset freq from whole callset agg
    subsets = args.subsets
    test = args.test
    chrom = args.chrom
    af_threshold = args.af_threshold

    hl.init(
        log=f"/generate_frequency_data{'.' + '_'.join(subsets) if subsets else ''}.log",
        default_reference="GRCh38",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    # TODO: Update to release = True once outlier filtering is complete,
    # possibly sample_meta=True if added
    vds = get_gnomad_v4_vds(test=test)
    meta_ht = meta.ht()

    logger.info("Adding metadata to VDS variant data cols...")
    vds = hl.vds.VariantDataset(
        vds.reference_data,
        vds.variant_data.annotate_cols(meta=meta_ht[vds.variant_data.col_key]),
    )

    if test or chrom:
        if test:
            logger.info("Filtering to DRD2 in test VDS for testing purposes...")
            test_interval = [
                hl.parse_locus_interval("chr11:113409605-113475691")
            ]  # NOTE: test on gene DRD2 for now, will update to chr22 when original PR goes in
            vds = hl.vds.filter_intervals(vds, test_interval)
            variants, samples = vds.variant_data.count()
            logger.info(
                "Test VDS has %s variants in DRD2 in %s samples...", variants, samples
            )
        else:
            logger.info("Filtering to chromosome %s...")
            vds = hl.vds.filter_chromosomes(vds, keep=f"chr{chrom}")

    # logger.info("Annotating non_ref hets pre-split...")
    # vds = annotate_non_ref_het(vds)

    # if args.subsets:
    #     vds = hl.vds.filter_samples(
    #         vds,
    #     )  # TODO: Where is subset information coming meta.project_meta.ukb is defined for non-ukb and ukd, how do we determine non-topmed v3.1 meta did (topmed=ht.s.startswith("NWD"))

    logger.info("Splitting VDS....")
    vds = hl.vds.split_multi(vds, filter_changed_loci=True)

    logger.info("Densifying VDS...")
    mt = hl.vds.to_dense_mt(vds)

    logger.info("Annotating GATK version for frequency groupings....")
    mt = annotate_gatk_version(mt)

    logger.info("Computing adj and sex adjusted genotypes...")
    mt = mt.annotate_entries(
        GT=adjusted_sex_ploidy_expr(
            mt.locus, mt.GT, mt.meta.sex_imputation.sex_karyotype
        ),
        adj=get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD),
    )

    if args.get_freq_and_high_ab:
        logger.info("Annotating frequencies...")
        mt = annotate_freq(
            mt,
            sex_expr=mt.meta.sex_imputation.sex_karyotype,
            pop_expr=mt.meta.population_inference.pop,
            # downsamplings=DOWNSAMPLINGS["v4"],
            additional_strata_expr={"gatk_version": mt.gatk_version},
            additional_strata_grouping_expr={"pop": mt.meta.population_inference.pop},
        )
        mt = annotate_high_ab_hets_by_group_membership(mt, gatk_expr=mt.gatk_version)
        ht = mt.rows()

    if args.adjust_freqs:
        # TODO: This should apply fix given the AF threshold and adjusts the
        # frequencies using the list created above
        ht = hom_alt_depletion_fix(ht, af_threshold)

    logger.info("Checkpointing frequency table...")
    ht = ht.select("freq", "high_ab_hets_by_group_membership")
    ht = ht.checkpoint(
        f"gs://gnomad-tmp-4day/freq/test_freq_aggs{'adjusted' if args.adjust_freqs else ''}2.ht",
        overwrite=True,
    )
    ht.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Runs a test on two partitions of the MT.", action="store_true"
    )
    parser.add_argument(
        "--chrom",
        help="If passed, script will only run on passed chromosome.",
        type=int,
    )
    parser.add_argument(
        "--overwrite", help="Overwrites existing files.", action="store_true"
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--subsets",
        help="Subsets to run frequency calculation on.",
        choices=SUBSETS["v4"],
    )
    parser.add_argument(
        "--get-freq-and-high-ab",
        help="Calculate frequencies and high AB sites per frequency grouping.",
        action="store_true",
    )
    parser.add_argument(
        "--adjust-freqs",
        help=(
            "Adjust each frequency entry to account for homozygous alternate depletion"
            " present in GATK versions released prior to 4.1.4.1."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--af-threshold",
        help=(
            "Threshold at which to adjust site group frequencies at sites for"
            " homozygous alternate depletion present in GATK versions released prior to"
            " 4.1.4.1."
        ),
        type=float,
        default=0.01,
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
