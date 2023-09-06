"""Script containing annotation related resources."""
from typing import Optional

from gnomad.resources.grch37.gnomad import EXOME_RELEASES, GENOME_RELEASES
from gnomad.resources.grch38.gnomad import SUBSETS
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v3.resources.basics import get_checkpoint_path
from gnomad_qc.v4.resources.constants import CURRENT_VERSION, VERSIONS

SUBSETS = SUBSETS["v4"]


def _annotations_root(
    version: str = CURRENT_VERSION,
    test: bool = False,
    data_type: str = "exomes",
) -> str:
    """
    Get root path to the variant annotation files.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "exomes".
    :return: Root path of the variant annotation files.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/annotations/{data_type}"
        if test
        else f"gs://gnomad/v{version}/annotations/{data_type}"
    )


def get_info(split: bool = True, test: bool = False) -> VersionedTableResource:
    """
    Get the gnomAD v4 info VersionedTableResource.

    :param split: Whether to return the split or multi-allelic version of the resource.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :return: gnomAD v4 info VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.info{'.split' if split else ''}.ht"
                )
            )
            for version in VERSIONS
        },
    )


def get_vep(
    version: str = CURRENT_VERSION, test: bool = False, data_type: str = "exomes"
) -> str:
    """
    Get the gnomAD v4 VEP annotation VersionedTableResource.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "exomes".
    :return: gnomAD v4 VEP VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.vep.ht"
                )
            )
            for version in VERSIONS
        },
    )


def validate_vep_path(
    version: str = CURRENT_VERSION, test: bool = False, data_type: str = "exomes"
) -> str:
    """
    Get the gnomAD v4 VEP annotation VersionedTableResource for validation counts.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "exomes".
    :return: gnomAD v4 VEP VersionedTableResource containing validity check.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.vep.validate.ht"
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


def info_vcf_path(version: str = CURRENT_VERSION, test: bool = False) -> str:
    """
    Path to sites VCF (input information for running VQSR).

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :return: String for the path to the info VCF.
    """
    return (
        f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.info.vcf.bgz"
    )


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

fam_stats = VersionedTableResource(
    CURRENT_VERSION,
    {
        version: TableResource(
            f"{_annotations_root(version)}/gnomad.exomes.v{version}.qc_fam_stats.ht"
        )
        for version in VERSIONS
    },
)


def get_freq(
    version: str = CURRENT_VERSION,
    test: bool = False,
    hom_alt_adjusted=False,
    chrom: Optional[str] = None,
    intermediate_subset: Optional[str] = None,
    finalized: bool = True,
) -> VersionedTableResource:
    """
    Get the frequency annotation table for a specified release.

    :param version: Version of annotation path to return.
    :param test: Whether to use a tmp path for tests.
    :param hom_alt_adjusted: Whether to return the hom alt adjusted frequency table.
    :param chrom: Chromosome to return frequency table for. Entire Table will be
        returned if not specified.
    :param intermediate_subset: Optional intermediate subset to return temp frequency
        Table for. Entire Table will be returned if not specified.
    :param finalized: Whether to return the finalized frequency table. Default is True.
    :return: Hail Table containing subset or overall cohort frequency annotations.
    """
    ht_name = f"gnomad.exomes.v{version}"
    if not finalized:
        if chrom:
            ht_name += f".{chrom}"
        if not hom_alt_adjusted:
            ht_name += ".pre_hom_alt_adjustment"
    if intermediate_subset:
        ht_name += f".{intermediate_subset}"
        if test:
            ht_name += ".test"
        ht_path = get_checkpoint_path(ht_name)
    else:
        if finalized:
            ht_name += ".final"
        if test:
            ht_name += ".test"
        ht_path = f"{_annotations_root(version, test)}/{ht_name}.ht"

    return VersionedTableResource(
        version, {version: TableResource(ht_path) for version in VERSIONS}
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


def get_insilico_predictors(
    version: str = CURRENT_VERSION,
    predictor: str = "cadd",
) -> str:
    """
    Get the path to the in silico predictors TableResource for a specified release.

    :param version: Version of annotation path to return.
    :param predictor: One of the in silico predictors available in gnomAD v4, including cadd, revel, primate_ai, splice_ai, and pangolin.
    :return: in silico predictor VersionedTableResource for gnomAD v4.
    """
    return VersionedTableResource(
        CURRENT_VERSION,
        {
            version: TableResource(
                path=f"gs://gnomad/v{version}/annotations/in_silico_predictors/gnomad.v{version}.{predictor}.grch38.ht"
            )
            for version in VERSIONS
        },
    )
