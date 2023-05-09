"""Script containing annotation related resources."""
from typing import Optional

from gnomad.resources.grch37.gnomad import EXOME_RELEASES, GENOME_RELEASES
from gnomad.resources.grch38.gnomad import SUBSETS

# TODO: Add a V4_SUBSETS to replace SUBSETS (currently includes tgp and hgdp which are not relevant for v4) # noqa
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v4.resources.constants import CURRENT_VERSION, VERSIONS


def _annotations_root(version: str = CURRENT_VERSION) -> str:
    """
    Get root path to the variant annotation files.

    :param version: Version of annotation path to return
    :return: Root path of the variant annotation files
    """
    return f"gs://gnomad/v{version}/annotations/exomes"


def get_info(split: bool = True) -> VersionedTableResource:
    """
    Get the gnomAD v4 info TableResource.

    :param split: Whether to return the split or multi-allelic version of the resource
    :return: gnomAD v4 info VersionedTableResource
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version)}/gnomad.exomes.v{version}.info{'.split' if split else ''}.ht"
                )
            )
            for version in VERSIONS
        },
    )


def get_vqsr_filters(
    model_id: str,
    split: bool = True,
    finalized: bool = False,
) -> VersionedTableResource:
    """
    Get the specified VQSR filtering annotation resource.

    :param model_id: VQSR filtering model id
    :param split: Split or multi-allelic version of the filtering file
    :param finalized: Whether to return the raw VQSR table or the finalized VQSR table representing determined cutoffs
    :return: VQSR filtering annotation file
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version)}/vqsr/gnomad.exomes.v{version}.{model_id}{'.finalized' if finalized else ''}{'.split' if split else ''}.ht"
            )
            for version in VERSIONS
        },
    )


def info_vcf_path(version: str = CURRENT_VERSION) -> str:
    """
    Path to sites VCF (input information for running VQSR).

    :param version: Version of annotation path to return
    :return: String for the path to the info VCF
    """
    return f"{_annotations_root(version)}/gnomad.exomes.v{version}.info.vcf.bgz"


def get_transmitted_singleton_vcf_path(
    adj: bool = False, version: str = CURRENT_VERSION
) -> str:
    """
    Provide the path to the transmitted singleton VCF used as input to VQSR.

    :param bool adj: Whether to use adj genotypes
    :param version: Version of transmitted singleton VCF path to return
    :return: String for the path to the transmitted singleton VCF
    """
    return (
        f'{_annotations_root(version)}/gnomad.exomes.v{version}.transmitted_singletons.{"adj" if adj else "raw"}.vcf.bgz'
    )


freq = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.frequencies.ht"
        )
        for version in VERSIONS
    },
)

qual_hist = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.qual_hists.ht"
        )
        for version in VERSIONS
    },
)

vep = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.vep.ht"
        )
        for version in VERSIONS
    },
)

qc_ac = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.qc_ac.ht"
        )
        for version in VERSIONS
    },
)

fam_stats = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.qc_fam_stats.ht"
        )
        for version in VERSIONS
    },
)

allele_data = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.qc_allele_data.ht"
        )
        for version in VERSIONS
    },
)


def get_freq(
    version: str = CURRENT_VERSION,
    subset: Optional[str] = None,
    test: bool = False,
    hom_alt_adjustment=False,
    chr: Optional[str] = None,
) -> VersionedTableResource:
    """
    Get the frequency annotation table for a specified release.

    :param version: Version of annotation path to return
    :param subset: One of the official subsets of the specified release (e.g., non_neuro, non_cancer,
        controls_and_biobanks) or a combination of them split by '-'
    :param hom_alt_adjustment: Whether to return the frequency table before the hom alt adjustment
    :param chr: Chromosome to return frequency table for.
    :return: Hail Table containing subset or overall cohort frequency annotations
    """
    if subset is not None:
        for s in subset.split("-"):
            if s not in SUBSETS:
                raise DataException(
                    f"{subset} subset is not one of the following official subsets:"
                    f" {SUBSETS}"
                )

    return VersionedTableResource(
        version,
        {
            version: TableResource(
                f"gs://gnomad-tmp/julia/dense_freq_10k_partitions/gnomad.exomes.v{version}{'.' + chr if chr else ''}.frequencies{'.' + subset if subset else ''}{'.hom_alt_adjusted' if hom_alt_adjustment else ''}.ht"
            )
            for version in VERSIONS
        },
    )


def get_freq_comparison(version1, data_type1, version2, data_type2):
    """
    Get Table resource for a frequency comparison between two gnomAD versions.

    Table contains results from a chi squared test and fishers exact test comparing the variant frequencies of two
    gnomAD versions/data types.

    :param version1: First gnomAD version in the frequency comparison. Table will be in the root annotation path for this gnomAD version.
    :param data_type1: Data type of first version in the frequency comparison. One of "exomes" or "genomes". Table will be in the root annotation path for this gnomAD data type.
    :param version2: Second gnomAD version in the frequency comparison.
    :param data_type2: Data type of first version in the frequency comparison. One of "exomes" or "genomes".
    :return: Hail Table containing results from chi squared test and fishers exact test
    """
    versions = [r + "_exomes" for r in EXOME_RELEASES] + [
        r + "_genomes" for r in GENOME_RELEASES + VERSIONS
    ]
    if (
        f"{version1}_{data_type1}" not in versions
        or f"{version2}_{data_type2}" not in versions
    ):
        raise DataException(
            "One of the versions/datatypes supplied doesn't exist. Possible options"
            f" are: {versions}"
        )

    ht_path = (
        f"gnomad.{data_type1}_v{version1}_{data_type2}_v{version2}.compare_freq.ht"
    )
    if version1 in VERSIONS:
        ht_path = f"{_annotations_root(version1)}/{ht_path}"
    else:
        ht_path = f"gs://gnomad/annotations/hail-0.2/ht/{data_type1}/{ht_path}"

    return TableResource(ht_path)
