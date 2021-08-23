from gnomad.utils.slack import slack_notifications
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v2.resources.sample_qc import *
from gnomad_qc.v2.resources import SUBPOPS, metadata_exomes_tsv_path, metadata_exomes_ht_path, metadata_genomes_ht_path, metadata_genomes_tsv_path
import argparse
import datetime
import logging
import hail as hl

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("finalize_sample_qc")
logger.setLevel(logging.INFO)


def add_release_annotations(ht: hl.Table) -> hl.Table:
    """

    :param Table ht: Table containing meta column annotations for the dataset
    :return: Table containing final 'high_quality' and 'release' sample status annotations
    :rtype: Table
    """
    ht = ht.annotate(high_quality=(hl.len(ht.hard_filters) == 0) & (hl.len(ht.pop_platform_filters) == 0))
    return ht.annotate(release=ht.high_quality & (hl.len(ht.perm_filters) == 0) & ~ht.related)


def collapse_small_pops(ht: hl.Table,  min_pop_size: int) -> hl.Table:
    """

    Collapses (sub)populations that are too small for release into others.
    When collapsing subpops, the name for the other category is composed of "o" +  2 first letters of the superpop
    The original RF population assignments are kept in the `rf_pop` and `rf_subpop` columns.

    :param ht: Input Table
    :return: Table with small populations collapsed
    :rtype: Table
    """

    def get_subpop_oth(pop: str):
        for superpop, subpops in SUBPOPS.items():
            if pop.upper() in subpops:
                return "o"+superpop[:2].lower()

        raise ValueError(f"Subpopulation {pop} not found in possible subpopulations.")

    ht = ht.persist()
    pop_counts = ht.aggregate(hl.agg.filter(ht.release, hl.agg.counter(ht.pop)))
    pop_collapse = {pop: "oth" for pop,n in pop_counts.items() if n < min_pop_size}
    pop_collapse = hl.literal(pop_collapse) if pop_collapse else hl.empty_dict(hl.tstr, hl.tstr)

    subpop_counts = ht.aggregate(hl.agg.filter(ht.release, hl.agg.counter(ht.subpop)))
    subpop_collapse = {subpop: get_subpop_oth(subpop) for subpop,n in subpop_counts.items() if n < min_pop_size}
    subpop_collapse = hl.literal(subpop_collapse) if subpop_collapse else hl.empty_dict(hl.tstr, hl.tstr)

    return ht.annotate(
        pop=pop_collapse.get(ht.pop, ht.pop),
        subpop=subpop_collapse.get(ht.subpop, ht.subpop),
        rf_pop=ht.pop,
        rf_subpop=ht.subpop
    )

# TODO check for possible circularity with producing these annotations.
def add_topmed_annotation(ht: hl.Table, data_type: str, cutoff: float):
    topmed_stats = hl.read_table(get_topmed_shared_sites_ht_path(data_type))
    topmed_stats = topmed_stats.annotate(topmed=(topmed_stats.n_topmed_singletons / topmed_stats.n_singletons) > cutoff)
    return ht.annotate(topmed=topmed_stats[ht.key].topmed)


def main(args):
    subpop_nfe_ht = hl.read_table(subpop_ht_path('nfe')).select('subpop', 'known_subpop')
    subpop_eas_ht = hl.read_table(subpop_ht_path('eas')).select('subpop', 'known_subpop')
    subpop_ht = subpop_nfe_ht.union(subpop_eas_ht)

    exome_ht_a = hl.read_table(qc_ht_path('exomes', 'hard_filters'))
    exome_ht_b = hl.read_table(qc_ht_path('exomes', 'platforms'))
    exome_ht_b = exome_ht_b.rename({f'PC{i}': f'callrate_pc{i}' for i in range(1, 11)})
    exome_ht_c = hl.read_table(qc_ht_path('exomes', 'pop_platform'))
    exome_ht_c = exome_ht_c.transmute(**{f'PC{i+1}': exome_ht_c.scores[i] for i in range(10)}).drop('qc_platform')

    genome_ht_a = hl.read_table(qc_ht_path('genomes', 'hard_filters'))
    genome_ht_c = hl.read_table(qc_ht_path('genomes', 'pop_platform'))
    genome_ht_c = genome_ht_c.transmute(**{f'PC{i+1}': genome_ht_c.scores[i] for i in range(10)}).drop('qc_platform')

    super_meta_exome_ht = exome_ht_a.annotate(**exome_ht_b[exome_ht_a.key],
                                              **exome_ht_c[exome_ht_a.key],
                                              **subpop_ht[exome_ht_a.key])
    super_meta_exome_ht = add_release_annotations(super_meta_exome_ht).repartition(10)
    super_meta_exome_ht = super_meta_exome_ht.key_by('s')
    super_meta_exome_ht = collapse_small_pops(super_meta_exome_ht, args.min_pop_size)
    super_meta_exome_ht = add_topmed_annotation(super_meta_exome_ht, 'exomes', args.exomes_topmed_dup_cutoff)
    super_meta_exome_ht.write(metadata_exomes_ht_path(args.metadata_version_label), args.overwrite)
    super_meta_exome_ht.flatten().export(metadata_exomes_tsv_path(args.metadata_version_label))

    super_meta_genome_ht = genome_ht_a.annotate(**genome_ht_c[genome_ht_a.key],
                                                **subpop_ht[genome_ht_a.key])
    super_meta_genome_ht = add_release_annotations(super_meta_genome_ht).repartition(10)
    super_meta_genome_ht = super_meta_genome_ht.key_by('s')
    super_meta_genome_ht = collapse_small_pops(super_meta_genome_ht, args.min_pop_size)
    super_meta_genome_ht = add_topmed_annotation(super_meta_genome_ht, 'genomes', args.genomes_topmed_dup_cutoff)
    super_meta_genome_ht.write(metadata_genomes_ht_path(args.metadata_version_label), args.overwrite)
    super_meta_genome_ht.flatten().export(metadata_genomes_tsv_path(args.metadata_version_label))

    logger.info('{0} exome, {1} genome high-quality samples remaining after sample QC'.format(
        super_meta_exome_ht.aggregate(hl.agg.count_where(super_meta_exome_ht.high_quality)),
        super_meta_genome_ht.aggregate(hl.agg.count_where(super_meta_genome_ht.high_quality))))
    logger.info('{0} exome, {1} genome release samples remaining after sample QC'.format(
        super_meta_exome_ht.aggregate(hl.agg.count_where(super_meta_exome_ht.release)),
        super_meta_genome_ht.aggregate(hl.agg.count_where(super_meta_genome_ht.release))))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--metadata_version_label', help='Label for metadata version (usually YYYY-MM-DD)', default=datetime.datetime.today().strftime('%Y-%m-%d'))
    parser.add_argument('--min_pop_size', help='Minimum population size for release purpose. Smaller populations are collapsed into oth', default=50, type=int)
    parser.add_argument('--exomes_topmed_dup_cutoff', help='Proportion of shared singletons between an individual exome and TopMed to be considered a duplicate.', default=0.6, type=float)
    parser.add_argument('--genomes_topmed_dup_cutoff', help='Proportion of shared singletons between an individual genome and TopMed to be considered a duplicate.', default=0.8, type=float)
    parser.add_argument('--overwrite', help='Overwrite pre-existing data', action='store_true')
    parser.add_argument('--slack_channel', help='Slack channel to post results and notifications to.')

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
