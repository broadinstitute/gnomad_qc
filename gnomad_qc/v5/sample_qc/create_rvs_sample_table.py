"""Script to retrieve gnomAD samples that share rare variants with an external callset."""

import argparse
import logging

import hail as hl
from gnomad.resources.grch38.reference_data import methylation_sites
from gnomad.utils.filtering import filter_to_autosomes

from gnomad_qc.v4.resources.basics import get_gnomad_v4_genomes_vds, get_gnomad_v4_vds
from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("get_shared_rare_variants")


def create_shared_sites_table(data_type: str, external_ht, test):
    """
    Create a table of shared sites between gnomAD and an external dataset.

    :param data_type: Data type (e.g., exomes or genomes).
    :param external_ht: Path to external dataset.
    :param test: Whether to use test gnomAD VDS, useful for testing only.
    """
    methylation_ht = methylation_sites.ht()
    release = release_sites(data_type=data_type, public=True)
    release = release.filter(
        hl.is_snp(release.alleles[0], release.alleles[1])
        & hl.or_else(
            methylation_ht[release.locus].methylation_level < 6, True
        )  # v2 used 0.6 MEAN but GRCh38 changed annotation and cutoffs, https://the-tgg.slack.com/archives/C03PP8Y98LA/p1698193782207179
        & (release.filters.length() == 0)
    )

    if data_type == "exomes":
        gnomad = get_gnomad_v4_vds(plit=True, release_only=True, test=test)
    else:
        gnomad = get_gnomad_v4_genomes_vds(split=True, release_only=True, test=test)
    gnomad = gnomad.variant_data

    mt = filter_to_autosomes(gnomad)
    mt = mt.filter_rows(hl.is_defined(release[mt.row_key]))
    mt = mt.annotate_rows(
        external_ac=external_ht[mt.row_key].call_stats.AC,
        gnomad_ac=release[mt.row_key].freq[0].AC,
    )

    ht = mt.annotate_cols(
        n_singletons=hl.agg.count_where((mt.gnomad_ac == 1) & mt.GT.is_het()),
        n_doubletons=hl.agg.count_where((mt.gnomad_ac < 3) & mt.GT.is_het()),
        n_tripletons=hl.agg.count_where((mt.gnomad_ac < 4) & mt.GT.is_het()),
        n_external_singletons=hl.agg.count_where(
            (mt.gnomad_ac == 1) & mt.GT.is_het() & hl.is_defined(mt.external_ac)
        ),
        n_external_doubletons=hl.agg.count_where(
            (mt.gnomad_ac < 3) & mt.GT.is_het() & hl.is_defined(mt.external_ac)
        ),
        n_external_tripletons=hl.agg.count_where(
            (mt.gnomad_ac < 4) & mt.GT.is_het() & hl.is_defined(mt.external_ac)
        ),
        n_both_singletons=hl.agg.count_where(
            (mt.gnomad_ac == 1) & mt.GT.is_het() & (mt.external_ac == 1)
        ),
    ).cols()
    return ht


def get_samples_to_drop(rvs_ht_path, r_singletons_cutoff):
    """
    Get samples to drop based on shared singletons ratio.

    :param rvs_ht_path: Path to shared rare variants table.
    :param r_singletons_cutoff: Cutoff for ratio of shared singletons. Samples with ratio above this cutoff will be dropped.
    """
    ht = hl.read_table(rvs_ht_path)
    ht = ht.annotate(
        shared_singletons_ratio=ht.n_both_singletons / ht.n_external_singletons
    )
    ht = ht.filter(ht.shared_singletons_ratio > r_singletons_cutoff)
    return ht


def main(args):
    """Retrieve gnomAD samples that share rare variants with external callset."""
    hl.init(
        log="/tmp/get_shared_rare_variants.log",
        tmp_dir=args.temp_path,
        copy_spark_log_on_error=True,
    )
    hl.default_reference("GRCh38")
    ext_freq_ann = args.external_freq_ann
    rvs_ht_path = args.rvs_ht_path
    external_ht = hl.read_table(args.external_input)

    if ext_freq_ann:
        if ext_freq_ann not in external_ht.row:
            raise ValueError(
                "External dataset does not contain frequency annotation"
                f" {ext_freq_ann}."
            )
        else:
            external_ht = external_ht.annotate(
                freq=external_ht[ext_freq_ann]
            )  # TODO: This is half baked, what if array?  what if struct?

    if args.test:
        external_ht = external_ht.filter(external_ht.locus.contig == "chr22")

    if args.create_rvs_table:
        create_shared_sites_table(
            data_type=args.data_type,
            external_ht=external_ht,
            test=args.test,
        ).write(rvs_ht_path, overwrite=args.overwrite)

    if args.get_samples_to_drop:
        get_samples_to_drop(
            rvs_ht_path=rvs_ht_path,
            r_singletons_cutoff=args.r_singletons_cutoff,
        ).write(args.samples_to_drop_ht_path, overwrite=args.overwrite)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-type",
        help="Data type (exomes or genomes).",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--external-input",
        help="Path to external dataset.",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--output-ht-path",
        help="Path to write output Hail Table to.",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--test",
        help="Whether to use test gnomAD VDS, useful for testing only.",
        action="store_true",
    )
    parser.add_argument(
        "--overwrite",
        help="Whether to overwrite existing table.",
        action="store_true",
    )
    parser.add_argument(
        "--temp-path",
        help="Path to write temporary files to.",
        default="gs://gnomad-tmp-4day/get_shared_rare_variants/",
        type=str,
    )
    parser.add_argument(
        "--create-rvs-table",
        help="Whether to create the shared rare variants table.",
        action="store_true",
    )
    parser.add_argument(
        "--rvs-ht-path",
        help="Path to write shared rare variants Hail Table to.",
        type=str,
    )
    parser.add_argument(
        "--external-freq-ann",
        help="External frequency annotation as string.",
        default="freq",
        type=str,
    )
    parser.add_argument(
        "--get-samples-to-drop",
        help="Whether to compute samples to drop.",
        action="store_true",
    )
    parser.add_argument(
        "--samples-to-drop-ht-path",
        help="Path to write samples to drop Hail Table to.",
        type=str,
    )
    parser.add_argument(
        "--r-singletons-cutoff",
        help=(
            "Cutoff for ratio of shared singletons. Samples with ratio above this"
            " cutoff will be dropped."
        ),
        default=0.8,
        type=float,
    )
    args = parser.parse_args()
    main(args)
