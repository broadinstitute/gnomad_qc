""" 
Script to generate LD scores for gnomAD v4 Genome SNPs & Indels. Example usage:

hailctl dataproc submit cluster_name generate_ld_data.py --test --overwrite --hgdp-subset --pop afr --re-call-stats \
    --generate-ld-mt \
    --generate-ld-pruned-set \
    --generate-ld-matrix \
    --generate-ld-scores
"""

import argparse
import sys

from gnomad.utils.annotations import get_adj_expr
from gnomad.utils.slack import slack_notifications
from hail.linalg import BlockMatrix
from hail.utils import new_temp_file

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v3.resources.release import hgdp_tgp_subset
from gnomad_qc.v4.resources import *

EXPECTED_CONTIGS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]

from typing import Optional

from gnomad.resources.grch38.gnomad import CURRENT_EXOME_RELEASE, CURRENT_GENOME_RELEASE

# Importing paths and Common+Rare Freqs
from gnomad_qc.v4.resources.ld_resources import *


def get_pop_counters(mt, label) -> dict:
    """
    Calculate and return genetic ancestry counts given an mt and column label.

    :param mt: Input MT
    :param label: Col field containing genetic ancestry information
    :return: Dict of pop counts
    """

    cut_dict = {
        f"{label}": hl.agg.filter(
            hl.is_defined(mt.meta.population_inference.pop)
            & (mt.meta.population_inference.pop != "oth"),
            hl.agg.counter(mt.meta.population_inference.pop),
        ),
    }
    cut_data = mt.aggregate_cols(hl.struct(**cut_dict))
    logger.info(f"Counts: {cut_data}")
    return cut_data


def _filter_call_stats(
    ht,
    freq: float = COMMON_FREQ,
    ac_cutoff: int = 1,
) -> hl.Table:
    """
    Filter ht to relevant variants. Must have field popfreq. This code is used for filtering ht and mt
    """

    ht = ht.filter(
        (ht.pop_freq.AC > ac_cutoff)
        & (ht.pop_freq.AN - ht.pop_freq.AC > 1)
        & (ht.pop_freq.AF > freq)
        & ~((ht.pop_freq.AF == 0.5) & (ht.pop_freq.homozygote_count == 0))
        #  the above filter is needed for "het-only" monomorphic sites
        #  as otherwise variance is 0, and so, row_correlation errors out
    )
    return ht


def _filter_ht_for_ld(
    ht,
    freq: float = COMMON_FREQ,
    label="gen_anc",
    pop="eas",
    ac_cutoff: int = 1,
) -> hl.MatrixTable:
    """
    Filter HT to only variants appropriate for LD analyses for a given genetic ancestry.

    :param mt: Input HT to filter
    :param label: Col field containing genetic ancestry information
    :param pop: Given genetic ancestry group to filter to
    :param common_only: Bool of whether to filter to only common variants (>0.5%)
    :param re_call_stats: Bool to re-calculate callstats. Needed for subsets and when sampling cols.
    :return: MT filtered to appropriate variants and samples
    """
    meta_index = (
        hl.enumerate(ht.freq_meta)
        .find(lambda f: (f[1].get(label) == pop))
        .collect()[0][0]
    )

    # No support for SV datasets in gnomAD v4 Production team code
    # Have to re-do callstats when noted frequencies don't 100% reflect what you are calculating on
    # This is for test sets or 1kg_tdpg
    # CANNOT re_call stats without GT
    ht = ht.filter((hl.len(ht.filters) == 0))
    pop_freq = ht.freq[meta_index]
    ht = ht.annotate(pop_freq=pop_freq)

    return _filter_call_stats(ht, freq=freq, ac_cutoff=ac_cutoff)


def filter_mt_for_ld(
    mt,
    label,
    pop,
    freq: float = COMMON_FREQ,
    re_call_stats: bool = False,
    ld_pruned_path: str = None,
    ac_cutoff: int = 1,
    call_rate_cutoff: float = 0.8,
) -> hl.MatrixTable:
    """
    Filter MT to only variants and samples appropriate for LD analyses

    :param mt: Input MT to filter
    :param label: Col field containing genetic ancestry information
    :param pop: Given genetic ancestry group to filter to
    :param common_only: Bool of whether to filter to only common variants (>0.5%)
    :param re_call_stats: Bool to re-calculate callstats. Needed for subsets and when sampling cols.
    :return: MT filtered to appropriate variants and samples
    """
    logger.info(
        f"With label: {label} and pop: {pop} and freq: {freq} and ac_cutoff: {ac_cutoff} and call_rate_cutoff: {call_rate_cutoff}"
    )

    pop_mt = mt.filter_cols(mt.meta.population_inference.pop == pop)
    meta_index = (
        hl.enumerate(pop_mt.freq_meta)
        .find(lambda f: (f[1].get(label) == pop))
        .collect()[0][0]
    )

    pop_count_int = pop_mt.count_cols()
    logger.info(f"Count1: {pop_mt.count()}")

    # If passed, read in pruned variants and remove them from our analyses
    # Do this here to save having to compute on them
    # ... and nothing here is sample-level, so we can count out
    if ld_pruned_path:
        logger.info("Filtering out LD pruned variants...")
        ld_ht = hl.read_table(ld_pruned_path)
        pop_mt = pop_mt.filter_rows(~hl.is_defined(ld_ht[pop_mt.row_key]))

    logger.info(f"Count2: {pop_mt.count()}")

    # No support for SV datasets in gnomAD v4 Production team code
    # Have to re-do callstats when noted frequencies don't 100% reflect what you are calculating on
    # This is for test sets or 1kg_tdpg
    if re_call_stats:
        logger.info("Regenerating call_stats...")
        call_stats = hl.agg.call_stats(pop_mt.GT, pop_mt.alleles)
        call_stats_bind = hl.bind(
            lambda cs: cs.annotate(
                AC=cs.AC[1], AF=cs.AF[1], homozygote_count=cs.homozygote_count[1]
            ),
            call_stats,
        )
        pop_freq = call_stats_bind
        pop_mt = pop_mt.annotate_rows(pop_freq=pop_freq)

        logger.info(f"Count_recallstats: {pop_mt.count()}")

    else:
        pop_freq = pop_mt.freq[meta_index]
        pop_mt = pop_mt.annotate_rows(pop_freq=pop_freq)
        logger.info(f"Count_not_recall_stats: {pop_mt.count()}")

    # Annotate callrate to check
    pop_mt = pop_mt.annotate_rows(callrate=pop_mt.pop_freq.AN / (pop_count_int * 2))

    pop_mt = pop_mt.filter_rows((hl.len(pop_mt.filters) == 0))
    pop_ht_rows = _filter_call_stats(pop_mt.rows(), freq=freq, ac_cutoff=ac_cutoff)
    pop_mt = pop_mt.filter_rows(hl.is_defined(pop_ht_rows[pop_mt.row_key]))

    logger.info(f"Count_final: {pop_mt.count()}")

    return pop_mt


def generate_ld_pruned_set(
    mt: hl.MatrixTable,
    pop_data: dict,
    data_type: str,
    r2: str = "0.2",
    freq: float = COMMON_FREQ,
    radius: int = 1000000,
    overwrite: bool = False,
    re_call_stats: bool = False,
    version: str = None,
    test: bool = False,
    ld_contig: str = None,
    ac_cutoff: int = 1,
    call_rate_cutoff: float = 0.8,
    adj: bool = False,
) -> None:
    """
    Generate and write set of variants uncorrelated with eachother. Wrapper for hl.ld_prune()

    :param mt: Input MT to filter
    :param pop_data: Dict of population counts
    :param data_type: Genetic data type for exomes or genomes. Exomes are currently too sparse, only runs on genomes.
    :param r2: via Hail: Squared correlation threshold (exclusive upper bound). Must be in the range [0.0, 1.0].
    :param radius: Used for bp_window_size, via Hail: Window size in base pairs (inclusive upper bound).
    :param overwrite: Bool to write over previous outputs or not
    :param re_call_stats: Bool to re-calculate callstats. Needed for subsets and when sampling cols.
    :param version: Version of files, either 'hgdp' or None, which ld_pruned_path() populates with most recent gnomAD genomes version.
    :param test: Filter to test versions and/or write to test paths.
    """

    for label, pops in dict(pop_data).items():
        for pop in pops:
            logger.info(f"Filtering for {pop}...")
            pop_mt = filter_mt_for_ld(
                mt=mt,
                label=label,
                pop=pop,
                freq=freq,
                re_call_stats=re_call_stats,
                ld_pruned_path=None,
                ac_cutoff=ac_cutoff,
                call_rate_cutoff=call_rate_cutoff,
            )
            logger.info(f"Count after filter_mt_for_ld() : {pop_mt.count()}")

            pruned_ht = hl.ld_prune(pop_mt.GT, r2=float(r2), bp_window_size=radius)

            logger.info(f"Count of pruned from hl.ld_prune() : {pruned_ht.count()}")

            ht = pop_mt.rows().select("pop_freq")
            ht = ht.filter(hl.is_defined(pruned_ht[ht.key]))

            # Write out pruned set
            ld_path = ld_pruned_path(
                data_type,
                pop,
                r2,
                version=version,
                test=test,
                freq=freq,
                ld_contig=ld_contig,
                adj=adj,
            )
            logger.info(f"Writing pruned ht for {pop} to {ld_path} ...")
            ht.write(ld_path, overwrite)


def generate_ld_matrix(
    mt,
    pop_data,
    data_type,
    radius: int = 1000000,
    freq: float = COMMON_FREQ,
    adj: bool = False,
    overwrite: bool = False,
    re_call_stats: bool = False,
    version: str = None,
    test: bool = False,
    ld_contig: str = None,
    r2: str = "0.2",
    custom_suffix: str = None,
    ac_cutoff: int = 1,
    call_rate_cutoff: float = 0.8,
) -> None:
    """
    Generate and write matrix of LD correlcations as Hail BlockMatrix using hl.ld_matrix. Read by generate_ld_scores_from_ld_matrix()

    :param mt: Input MT to filter and generate LDs for
    :param pop_data: Dict of population counts
    :param data_type: Genetic data type for exomes or genomes. Exomes are currently too sparse, only runs on genomes.
    :param radius: Used for bp_window_size, via Hail: Window size in base pairs (inclusive upper bound).
    :param common_only: Bool of whether to filter to only common variants (>0.5%)
    :param adj: Whether to use adj ("adjusted" or "passing") frequency
    :param overwrite: Bool to write over previous outputs or not
    :param re_call_stats: Bool to re-calculate callstats. Needed for subsets and when sampling cols.
    :param version: Version of files, either 'hgdp' or None, which ld_pruned_path() populates with most recent gnomAD genomes version.
    :param test: Filter to test versions and/or write to test paths.
    """
    # Takes about 4 hours on 20 n1-standard-8 nodes (with SSD - not sure if necessary) per population
    # Total of ~37 hours ($400)

    for label, pops in dict(pop_data).items():
        for pop in pops:

            pop_mt = filter_mt_for_ld(
                mt,
                label,
                pop,
                freq=freq,
                re_call_stats=re_call_stats,
                ac_cutoff=ac_cutoff,
                call_rate_cutoff=call_rate_cutoff,
                ld_pruned_path=ld_pruned_path(
                    data_type,
                    pop,
                    r2=r2,
                    version=version,
                    test=test,
                    freq=freq,
                    ld_contig=ld_contig,
                    adj=adj,
                ),
            )

            pop_mt.rows().select("pop_freq").add_index().write(
                ld_index_path(
                    data_type=data_type,
                    pop=pop,
                    freq=freq,
                    adj=adj,
                    version=version,
                    test=test,
                    ld_contig=ld_contig,
                    custom_suffix=custom_suffix,
                ),
                overwrite,
            )
            ld = hl.ld_matrix(pop_mt.GT.n_alt_alleles(), pop_mt.locus, radius)
            if data_type != "genomes_snv_sv":
                ld = ld.sparsify_triangle()
            ld.write(
                ld_matrix_path(
                    data_type,
                    pop,
                    freq=freq,
                    adj=adj,
                    version=version,
                    test=test,
                    ld_contig=ld_contig,
                    custom_suffix=custom_suffix,
                ),
                overwrite,
            )


def generate_ld_scores_from_ld_matrix(
    pop_data,
    data_type,
    freq: float = COMMON_FREQ,
    call_rate_cutoff=0.8,
    adj: bool = False,
    radius: int = 1000000,
    overwrite: bool = False,
    version: str = None,
    test: bool = False,
    ld_contig: str = None,
    custom_suffix: str = None,
) -> None:
    """
    Generate LD scores from Hail BlockMatrix written by generate_ld_scores_from_ld_matrix().

    :param pop_data: Dict of population counts
    :param data_type: Genetic data type for exomes or genomes. Exomes are currently too sparse, only runs on genomes.
    :param min_frequency: Lowest frequency for variants to calculate LD scores of.
    :param call_rate_cutoff: Lowest call rate for variants to calculate LD scores of.
    :param common_only: Bool of whether to filter to only common variants (>0.5%)
    :param adj: Whether to use adj ("adjusted" or "passing") frequency
    :param radius: Used for bp_window_size, via Hail: Window size in base pairs (inclusive upper bound).
    :param overwrite: Bool to write over previous outputs or not
    :param version: Version of files, either 'hgdp' or None, which ld_pruned_path() populates with most recent gnomAD genomes version.
    :param test: Filter to test versions and/or write to test paths.
    """
    # This function required a decent number of high-mem machines (with an SSD for good measure) to complete the AFR
    # For the rest, on 20 n1-standard-8's, 1h15m to export block matrix, 15
    # mins to compute LD scores per population (~$150 total)

    # logger.info(
    #     f'Min frequency of {min_frequency} is {">=" if min_frequency >= freq_cutoff else "<"} cutoff of {freq_cutoff}...'
    # )

    for label, pops in dict(pop_data).items():
        for pop, n in pops.items():
            ht = hl.read_table(
                ld_index_path(
                    data_type=data_type,
                    pop=pop,
                    freq=freq,
                    adj=adj,
                    version=version,
                    test=test,
                    ld_contig=ld_contig,
                    custom_suffix=custom_suffix,
                )
            )

            ht = ht.annotate(callrate=(ht.pop_freq.AN / n) / 2)

            ht = ht.filter(
                (ht.pop_freq.AF >= freq)
                & (ht.pop_freq.AF <= 1 - freq)
                & (ht.pop_freq.AN / n >= 2 * call_rate_cutoff)
            ).add_index(name="new_idx")

            indices = ht.idx.collect()

            r2 = BlockMatrix.read(
                ld_matrix_path(
                    data_type=data_type,
                    pop=pop,
                    freq=freq,
                    adj=adj,
                    version=version,
                    test=test,
                    ld_contig=ld_contig,
                    custom_suffix=custom_suffix,
                )
            )
            r2 = r2.filter(indices, indices) ** 2
            r2_adj = ((n - 1.0) / (n - 2.0)) * r2 - (1.0 / (n - 2.0))

            out_name = ld_scores_path(
                data_type,
                pop,
                adj,
                version=version,
                test=test,
                freq=freq,
                ld_contig=ld_contig,
                custom_suffix=custom_suffix,
                call_rate_cutoff=call_rate_cutoff,
            )
            compute_and_annotate_ld_score(ht, r2_adj, radius, out_name, overwrite)


def compute_and_annotate_ld_score(ht, r2_adj, radius, out_name, overwrite) -> None:
    """
    Annotate LD scores onto Hail Table containing variants and write to path.

    :param ht: Hail Table of variants and LD indices
    :param r2: via Hail: Squared correlation threshold (exclusive upper bound). Must be in the range [0.0, 1.0].
    :param radius: Used for bp_window_size, via Hail: Window size in base pairs (inclusive upper bound).
    :param out_name: Path to write outputs to. Designed to be output of ld_scores_path()
    :param overwrite: Bool to write over previous outputs or not
    """
    starts_and_stops = hl.linalg.utils.locus_windows(ht.locus, radius, _localize=False)
    r2_adj = r2_adj._sparsify_row_intervals_expr(starts_and_stops, blocks_only=False)

    l2row = r2_adj.sum(axis=0).T
    l2col = r2_adj.sum(axis=1)
    l2 = l2row + l2col - 1
    l2_bm_tmp = new_temp_file()
    l2_tsv_tmp = new_temp_file()

    l2.write(l2_bm_tmp, force_row_major=True)
    BlockMatrix.export(l2_bm_tmp, l2_tsv_tmp)

    ht_scores = hl.import_table(l2_tsv_tmp, no_header=True, impute=True)
    ht_scores = ht_scores.add_index().rename({"f0": "ld_score"})
    ht_scores = ht_scores.key_by("idx")
    ht = ht.annotate(**ht_scores[ht.new_idx]).select_globals()
    ht.filter(hl.is_defined(ht.ld_score)).write(out_name, overwrite)


def _ld_test_mt(
    mt: hl.MatrixTable,
    mt_contig: bool = False,
    sample_cols: bool = False,
) -> hl.MatrixTable:
    if mt_contig:
        logger.info("Filtering to chr22...")
        mt = mt.filter_rows(mt.locus.contig == "chr22")
    if sample_cols:
        logger.info("Downsampling all cols by 0.1...")
        mt = mt.sample_cols(0.1)

    return mt


def generate_ld_mt(
    ld_mt_path: str,
    is_hgdp: bool = False,
    test: bool = False,
    pop: str = "eas",
    adj: bool = False,
    freq: float = COMMON_FREQ,
    mt_contig: str = None,
    do_v2_samples: bool = False,
    hapmap: bool = False,
    overwrite=False,
) -> hl.MatrixTable:

    if is_hgdp:
        logger.info("Reading in HGDP_TGP Subset...")
        mt = hgdp_tgp_subset(dense=True, public=True).mt()
        mt = mt.annotate_globals(
            freq_meta=mt.gnomad_freq_meta, freq_index_dict=mt.gnomad_freq_index_dict
        )
        mt = mt.annotate_rows(freq=mt.gnomad_freq)
        mt = mt.annotate_cols(
            meta=hl.struct(
                population_inference=hl.struct(pop=mt.gnomad_population_inference.pop)
            )
        )

        # Filter to appropriate genetic ancestry and frequency b efore checkpointing
        mt = mt.filter_cols(mt.meta.population_inference.pop == pop)
        mt = mt.filter_rows(mt.freq[mt.freq_index_dict[f"gnomad_{pop}"]].AF > freq)

        if mt_contig:
            mt = mt.filter_rows(mt.locus.contig == mt_contig)

        if test:
            mt = _ld_test_mt(mt, filter_contig=True, sample_cols=True)

        logger.info(
            f"Checkpointing HGDP_TGP Subset on {mt_contig if mt_contig else 'full'}"
        )
        mt = mt.checkpoint(ld_mt_path, overwrite=overwrite)

    else:
        logger.info("Reading in gnomAD release HT and VDS...")

        ht_release = hl.read_table(release.release_ht_path(data_type="genomes"))
        ht_filter = _filter_ht_for_ld(ht_release, freq=freq, label="gen_anc", pop=pop)
        if mt_contig:
            ht_filter = ht_filter.filter(ht_filter.locus.contig == mt_contig)

        ht_filter = ht_filter.select("freq", "filters").checkpoint(new_temp_file())

        filter_samples_ht = meta(version="4.0", data_type="genomes").ht()
        filter_samples_ht = filter_samples_ht.filter(
            filter_samples_ht.population_inference.pop == pop
        )

        if do_v2_samples:
            logger.info("Filtering to only v2 samples...")
            filter_samples_ht = filter_samples_ht.filter(
                (filter_samples_ht.project_meta.v2_pop == pop)
                & (filter_samples_ht.project_meta.v2_release)
                & (filter_samples_ht.project_meta.v2_high_quality)
            )

        naive_coalesce_partitions = 9800 if mt_contig is None else 1000

        if mt_contig == "chr1":
            naive_coalesce_partitions = 2500

        if hapmap:
            logger.info("Filtering to only hapmap variants...")
            ht_hm3 = hl.read_table(
                f"gs://gnomad-tmp/ld/HM3.UKBB.{pop.upper()}.qc.snplist.grch38.ht"
            )
            ht_filter = ht_filter.filter(hl.is_defined(ht_hm3[ht_filter.key]))
            logger.info(
                f"Count after filtering to hapmap variants: {ht_filter.count()}"
            )

        entries_to_keep = ["GT"] if not adj else ["GT", "GQ", "DP", "AD"]

        vds = basics.get_gnomad_v4_genomes_vds(
            test=test,
            annotate_meta=True,
            release_only=True,
            split=True,
            chrom=[mt_contig] if mt_contig else None,
            filter_variant_ht=ht_filter,
            naive_coalesce_partitions=naive_coalesce_partitions,
            entries_to_keep=entries_to_keep,
            filter_samples_ht=filter_samples_ht,
        )

        vds.variant_data = vds.variant_data.checkpoint(new_temp_file())
        vds.reference_data = vds.reference_data.checkpoint(new_temp_file())

        logger.info("From DENSE MT, performance beware...")
        mt = hl.vds.to_dense_mt(vds)

        logger.info(f"Filtering to relevant cols - only meta - and no VDS/MT rows...")
        mt = mt.select_cols(mt.meta)
        mt = mt.select_rows()

        # NOTE: QUESTION IF THIS SHOULD BE RAN OR NOT
        # NOTE: need to checkpoint either here or step below where noted
        logger.info(
            f"Checkpointing unannotated {pop} on {mt_contig if mt_contig else 'full'}"
        )

        if adj:
            logger.info("Filtering to adj sites only...")
            adj_expr = get_adj_expr(mt.GT, mt.GQ, mt.DP, mt.AD)
            mt = mt.filter_entries(adj_expr)
            mt = mt.select_entries("GT")

        # mt = mt.checkpoint(new_temp_file())

        # logger.info('CUSTOM READ FROM PRIOR CHECKPOINT..')
        # mt = hl.read_matrix_table('gs://gnomad-tmp-4day/ld_tmp//XLmYz1uA7NYEY2jgB55Yuc')
        # This part begins potential performance concerns - what is up here?
        # Ideally, would checkpoint both of these beforehand.
        # Is very slow if MT isn't checkpointed
        # Checkpointing said MT is also very slow...
        logger.info("Annotating information from release HT...")
        mt = mt.annotate_rows(
            freq=ht_filter[mt.row_key].freq, filters=ht_filter[mt.row_key].filters
        )
        mt = mt.annotate_globals(**hl.eval(ht_filter.globals))

        logger.info("Checkpointing annotated MT...")
        mt = mt.checkpoint(ld_mt_path, overwrite=overwrite)

    return mt


def main(args):
    hl.init(
        log="/ld_assessment.log",
        tmp_dir="gs://gnomad-tmp-4day/ld_tmp/",
        gcs_requester_pays_configuration="broad-mpg-gnomad",
    )
    hl._set_flags(use_new_shuffle="1")

    # Only run on genomes so far, no choice to run on exomes, too sparse
    data_type = "genomes"
    if args.exomes:
        raise NotImplementedError("LD Code does not work for exomes yet. Too sparse:(")

    # Define command line args used more than once
    test = args.test
    adj = args.adj
    is_hgdp = args.hgdp_subset
    overwrite = args.overwrite
    pop = args.pop
    custom_mt = args.custom_mt
    do_v2_samples = args.do_v2_samples
    hapmap = args.do_hapmap

    # Chromosome filter logic
    mt_contig = args.mt_contig
    ld_contig = args.ld_contig

    if test:
        mt_contig = "chr22"

    if mt_contig:
        ld_contig = mt_contig

    # Get correct AF based on input logic...
    vds_freq = args.vds_af

    ld_freq = args.ld_af
    ld_ac = args.ld_ac

    if vds_freq < ld_freq:
        pass  # this is the desired behavior, for the vds freq to be more permisive
    if ld_freq < vds_freq:
        ld_freq = vds_freq

    if test and not args.re_call_stats:
        logger.info(
            "WARNING: Test subsets without re-called stats may result in None or inaccurate LD scores."
        )

    # Version is populated via ld_resources.py if None
    version = None
    if is_hgdp:
        version = "hgdp"
    elif do_v2_samples:
        version = "v2samples"

    if version:
        version += "hapmap"
    else:
        version = "hapmap"

    ld_mt_path = ld_mt_checkpoint_path(
        data_type="genomes",
        freq=vds_freq,
        pop=pop,
        version=version,
        test=test,
        mt_contig=mt_contig,
        adj=adj,
    )

    if args.generate_ld_mt:
        mt = generate_ld_mt(
            is_hgdp=is_hgdp,
            test=test,
            pop=pop,
            mt_contig=mt_contig,
            freq=vds_freq,
            ld_mt_path=ld_mt_path,
            overwrite=overwrite,
            do_v2_samples=do_v2_samples,
            adj=adj,
            hapmap=hapmap,
        )
    else:
        if not custom_mt:
            mt = hl.read_matrix_table(ld_mt_path)
        else:
            mt = hl.read_matrix_table(custom_mt)

    label = "gen_anc" if not is_hgdp else "pop"
    pop_data = get_pop_counters(
        mt, label=label
    )  # NOTE: could be irrelevant later? tbd if needed

    contig_list = EXPECTED_CONTIGS if not ld_contig else [ld_contig]

    for ld_contig in contig_list:
        # NOTE: This FULLY serializes the code, but hopefully it works
        # Would we save time without so many separate filter_mt_for_ld calls ?
        # The only point of concern is hl.ld_prune() and the size of the matrix
        # Need to filter the mt to only the contig we are working on
        # This is because the LD matrix is too large otherwise
        mt = mt.filter_rows(mt.locus.contig == ld_contig)
        logger.info(f"Filtering to {ld_contig}...")

        # Previously established that THIS needs a checkpoint to have been ran to work effectively
        if args.generate_ld_pruned_set:
            logger.info(
                "Generating LD Pruned Set of uncorrelated variants, using hl.ld_prune()..."
            )
            generate_ld_pruned_set(
                mt=mt,
                pop_data=pop_data,
                data_type=data_type,
                r2=args.r2,
                freq=ld_freq,
                ac_cutoff=ld_ac,
                radius=args.radius,
                overwrite=overwrite,
                re_call_stats=args.re_call_stats,
                version=version,
                test=test,
                ld_contig=ld_contig,
                adj=adj,
            )

        if args.generate_ld_matrix:
            logger.info(
                "Generating Hail BlockMatrix of variants correlations to other variants, using hl.ld_matrix()..."
            )
            generate_ld_matrix(
                mt=mt,
                pop_data=pop_data,
                data_type=data_type,
                radius=args.radius,
                freq=ld_freq,
                ac_cutoff=ld_ac,
                adj=args.adj,
                overwrite=overwrite,
                re_call_stats=args.re_call_stats,
                version=version,
                test=test,
                ld_contig=ld_contig,
                custom_suffix=args.custom_suffix,
            )

        if args.generate_ld_scores:
            logger.info("Generating in LD scores and annotating onto variant HT...")
            generate_ld_scores_from_ld_matrix(
                pop_data=pop_data,
                data_type=data_type,
                freq=ld_freq,
                call_rate_cutoff=args.min_call_rate,
                adj=args.adj,
                version=version,
                overwrite=overwrite,
                test=test,
                ld_contig=ld_contig,
                custom_suffix=args.custom_suffix,
            )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--exomes",
        help="Input MT is exomes. One of --exomes or --genomes is required.",
        action="store_true",
    )
    parser.add_argument(
        "--generate-ld-mt",
        help="Generate pop-specific Hail MatrixTable of variants and calls to calculate LD for. ",
        action="store_true",
    )
    parser.add_argument(
        "--generate-ld-pruned-set",
        help="Calculates LD pruned set of variants",
        action="store_true",
    )
    parser.add_argument(
        "--generate-ld-matrix", help="Calculates LD matrix", action="store_true"
    )
    parser.add_argument(
        "--vds-af",
        help="Minimum genetic ancestry-specific allele frequency to pull from for VDS.",
        type=float,
        default=RARE_FREQ,
    )
    parser.add_argument(
        "--adj",
        help="Calculates LD matrix only on using high-quality genotypes",
        action="store_true",
    )
    parser.add_argument(
        "--generate-ld-scores",
        help="Calculates LD scores from LD matrix",
        action="store_true",
    )
    parser.add_argument(
        "--min-call-rate",
        help="Minimum call rate to compute LD scores (default 0.8)",
        default=0.8,
        type=float,
    )
    parser.add_argument(
        "--ld-ac",
        help="Minimum allele count to compute LD scores (default 1).",
        default=1,
        type=int,
    )
    parser.add_argument(
        "--ld-af",
        help="Minimum allele frequency to compute LD scores (default 0.0005).",
        default=RARE_FREQ,
        type=float,
    )
    parser.add_argument(
        "--r2", help="r-squared to which to prune LD (default 0.2)", default="0.2"
    )
    parser.add_argument(
        "--radius",
        help="Radius at which to calculate LD information (bp; default 1e6)",
        default=1000000,
        type=int,
    )
    parser.add_argument(
        "--test",
        help="Use test dataset for whichever callset requested.",
        action="store_true",
    )
    parser.add_argument(
        "--hgdp-subset", help="Use hgdp dataset for callset.", action="store_true"
    )
    parser.add_argument(
        "--re-call-stats",
        help="Regenerate callstats for LD work. Can be useful for some subsets.",
        action="store_true",
    )
    parser.add_argument(
        "--pop",
        help="Which individual pop to run LD on.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--mt-contig",
        help="Which chromosome to pull VDS from. Applies for reading and writing.",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--ld-contig",
        help="Which chromosome to run LD on. Applies for reading and writing.",
        type=str,
        required=False,
    )
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument(
        "--custom-mt",
        type=str,
        default=None,
        help="Custom MT to run LD on",
    )
    parser.add_argument(
        "--do-v2-samples",
        action="store_true",
        help="Only use samples also in v2. Good for verification!",
    )
    parser.add_argument(
        "--custom-suffix",
        type=str,
        help="Custom string to append to output files",
        required=False,
    )
    parser.add_argument(
        "--do-hapmap",
        action="store_true",
        help="Filter to HapMap3 SNPs when pulling VDS.",
    )

    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
