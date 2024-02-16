"""Script containing annotation related resources."""
from typing import Optional

from gnomad.resources.grch38.gnomad import SUBSETS
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VariantDatasetResource,
    VersionedTableResource,
)

from gnomad_qc.v4.resources.constants import (
    ALL_SITES_AN_RELEASES,
    ANNOTATION_VERSIONS,
    COMBINED_FAF_RELEASES,
    CURRENT_ALL_SITES_AN_RELEASE,
    CURRENT_ANNOTATION_VERSION,
    CURRENT_COMBINED_FAF_RELEASE,
    CURRENT_HGDP_TGP_RELEASE,
    CURRENT_VERSION,
    HGDP_TGP_RELEASES,
    VERSIONS,
)

SUBSETS = SUBSETS["v4"]


def _annotations_root(
    version: str = CURRENT_ANNOTATION_VERSION,
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
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.info{'.split' if split else ''}.ht"
                )
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_vep(test: bool = False, data_type: str = "exomes") -> VersionedTableResource:
    """
    Get the gnomAD v4 VEP annotation VersionedTableResource.

    :param test: Whether to use a tmp path for analysis of the test Table instead of the full v4 Table.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "exomes".
    :return: gnomAD v4 VEP VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.vep.ht"
                )
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def validate_vep_path(
    test: bool = False, data_type: str = "exomes"
) -> VersionedTableResource:
    """
    Get the gnomAD v4 VEP annotation VersionedTableResource for validation counts.

    :param test: Whether to use a tmp path for analysis of the test VDS instead of the full v4 VDS.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes". Default is "exomes".
    :return: gnomAD v4 VEP VersionedTableResource containing validity check.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.vep.validate.ht"
                )
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_trio_stats(test: bool = False) -> VersionedTableResource:
    """
    Get the gnomAD v4 trio stats VersionedTableResource.

    :param test: Whether to use a tmp path for testing.
    :return: gnomAD v4 trio stats VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.trio_stats.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_sib_stats(test: bool = False) -> VersionedTableResource:
    """
    Get the gnomAD v4 sibling stats VersionedTableResource.

    :param test: Whether to use a tmp path for testing.
    :return: gnomAD v4 sibling stats VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.sib_stats.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_variant_qc_annotations(test: bool = False) -> VersionedTableResource:
    """
    Return the VersionedTableResource to the RF-ready annotated Table.

    Annotations that are included in the Table:

        Features for RF:
            - variant_type
            - allele_type
            - n_alt_alleles
            - has_star
            - AS_QD
            - AS_pab_max
            - AS_MQRankSum
            - AS_SOR
            - AS_ReadPosRankSum

        Training sites (bool):
            - transmitted_singleton
            - sibling_singleton
            - fail_hard_filters - (ht.QD < 2) | (ht.FS > 60) | (ht.MQ < 30)

    :param test: Whether to use a tmp path for testing.
    :return: Table with variant QC annotations.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.variant_qc_annotations.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def info_vcf_path(
    info_method: str = "AS",
    version: str = CURRENT_ANNOTATION_VERSION,
    split: bool = False,
    test: bool = False,
) -> str:
    """
    Path to sites VCF (input information for running VQSR).

    :param info_method: Method for generating info VCF. Must be one of "AS", "quasi",
        or "set_long_AS_missing". Default is "AS".
    :param version: Version of annotation path to return.
    :param split: Whether to return the split or multi-allelic version of the resource.
    :param test: Whether to use a tmp path for analysis of the test VDS instead of the
        full v4 VDS.
    :return: String for the path to the info VCF.
    """
    if info_method not in ["AS", "quasi", "set_long_AS_missing"]:
        raise ValueError(
            f"Invalid info_method: {info_method}. Must be one of 'AS', 'quasi', or "
            "'long_AS_missing_info'."
        )
    return (
        f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.info.{info_method}{'.split' if split else ''}.vcf.bgz"
    )


def get_true_positive_vcf_path(
    version: str = CURRENT_ANNOTATION_VERSION,
    test: bool = False,
    adj: bool = False,
    true_positive_type: str = "transmitted_singleton",
) -> str:
    """
    Provide the path to the transmitted singleton VCF used as input to VQSR.

    :param version: Version of true positive VCF path to return.
    :param test: Whether to use a tmp path for testing.
    :param adj: Whether to use adj genotypes.
    :param true_positive_type: Type of true positive VCF path to return. Should be one
        of "transmitted_singleton", "sibling_singleton", or
        "transmitted_singleton.sibling_singleton". Default is "transmitted_singleton".
    :return: String for the path to the true positive VCF.
    """
    tp_types = [
        "transmitted_singleton",
        "sibling_singleton",
        "transmitted_singleton.sibling_singleton",
    ]
    if true_positive_type not in tp_types:
        raise ValueError(f"true_positive_type must be one of {tp_types}")
    return (
        f'{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.{true_positive_type}.{"adj" if adj else "raw"}.vcf.bgz'
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


def get_downsampling(
    test: bool = False, subset: Optional[str] = None
) -> VersionedTableResource:
    """
    Get the downsampling annotation table.

    :param test: Whether to use a tmp path for tests. Default is False.
    :param subset: Optional subset to return downsampling Table for. Downsampling for
        entire dataset will be returned if not specified.
    :return: Hail Table containing subset or overall dataset downsampling annotations.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.downsampling{f'.{subset}' if subset else ''}.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_freq(
    version: str = CURRENT_VERSION,
    data_type: str = "exomes",
    test: bool = False,
    hom_alt_adjusted=False,
    chrom: Optional[str] = None,
    intermediate_subset: Optional[str] = None,
    finalized: bool = True,
) -> TableResource:
    """
    Get the frequency annotation table for a specified release.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
    :param test: Whether to use a tmp path for tests.
    :param hom_alt_adjusted: Whether to return the hom alt adjusted frequency table.
    :param chrom: Chromosome to return frequency table for. Entire Table will be
        returned if not specified.
    :param intermediate_subset: Optional intermediate subset to return temp frequency
        Table for. Entire Table will be returned if not specified.
    :param finalized: Whether to return the finalized frequency table. Default is True.
    :return: Hail Table containing subset or overall cohort frequency annotations.
    """
    ht_name = f"gnomad.{data_type}.v{version}.frequencies"
    if data_type == "exomes" and not finalized:
        if chrom:
            ht_name += f".{chrom}"
        if not hom_alt_adjusted:
            ht_name += ".pre_hom_alt_adjustment"

    if data_type == "exomes" and intermediate_subset:
        ht_name += f".{intermediate_subset}"
        if test:
            ht_name += ".test"
        ht_path = f"{_annotations_root(version, test)}/temp/{ht_name}.ht"
    else:
        if finalized:
            ht_name += ".final"
        if test:
            ht_name += ".test"
        ht_path = f"{_annotations_root(version, test, data_type)}/{ht_name}.ht"

    return TableResource(ht_path)


def get_all_sites_an_and_qual_hists(test: bool = False) -> VersionedTableResource:
    """
    Get the all sites AN and qual hists TableResource.

    :param test: Whether to use a tmp path for testing.
    :return: Hail Table containing all sites AN and qual hists annotations.
    """
    return VersionedTableResource(
        CURRENT_ALL_SITES_AN_RELEASE["exomes"],
        {
            version: TableResource(
                f"{_annotations_root(version, test=test)}/gnomad.exomes.v{version}.all_sites_an_and_qual_hists.ht"
            )
            for version in ALL_SITES_AN_RELEASES["exomes"]
        },
    )


def get_combined_frequency(
    test: bool = False, filtered: bool = False
) -> VersionedTableResource:
    """
    Get the combined v4 genome and exome frequency annotation VersionedTableResource.

    :param test: Whether to use a tmp path for testing.
    :param filtered: Whether to return the resource for the filtered combined frequency.
    :return: Hail Table containing combined frequency annotations.
    """
    return VersionedTableResource(
        CURRENT_COMBINED_FAF_RELEASE,
        {
            version: TableResource(
                f"{_annotations_root(version, data_type='joint', test=test)}/gnomad.joint.v{version}.frequencies{'.filtered' if filtered else ''}.ht"
            )
            for version in COMBINED_FAF_RELEASES
        },
    )


def get_freq_comparison(
    method: str, test: bool = False, filtered: bool = False
) -> VersionedTableResource:
    """
    Get VersionedTableResource for a frequency comparison between v4 genomes and exomes.

    Table contains results from one of the following comparison methods:
        - 'contingency_table_test': Hail's contingency table test -- chi-squared or
           Fisher’s exact test of independence depending on min allele count.
        - 'cmh_test': Cochran–Mantel–Haenszel test -- stratified test of independence
           for 2x2xK contingency tables.

    :param method: Method used to compare frequencies between v4 genomes and exomes.
        Can be one of `contingency_table_test` or `cmh_test`.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param filtered: Whether to return the filtered frequency comparison Table.
    :return: VersionedTableResource containing results from the specified comparison
        `method`.
    """
    methods = ["contingency_table_test", "cmh_test"]
    if method not in methods:
        raise DataException(
            f"The method: {method} is not an option. Possible options are: {methods}"
        )

    return VersionedTableResource(
        CURRENT_COMBINED_FAF_RELEASE,
        {
            version: TableResource(
                f"{_annotations_root(version, data_type='joint', test=test)}/gnomad.joint.v{version}.compare_frequencies.{method}{'.filtered' if filtered else ''}.ht"
            )
            for version in COMBINED_FAF_RELEASES
        },
    )


def get_insilico_predictors(predictor: str = "cadd") -> VersionedTableResource:
    """
    Get the path to the in silico predictors TableResource for a specified release.

    :param predictor: One of the in silico predictors available in gnomAD v4, including
        cadd, revel, primate_ai, splice_ai, and pangolin.
    :return: in silico predictor VersionedTableResource for gnomAD v4.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                path=f"gs://gnomad/v{version}/annotations/in_silico_predictors/gnomad.v{version}.{predictor}.grch38.ht"
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def get_vrs(
    original_annotations: bool = False,
    test: bool = False,
    data_type: str = "exomes",
) -> VersionedTableResource:
    """
    Get the gnomAD v4 VersionedTableResource containing VRS annotations.

    :param original_annotations: Whether to obtain the original input Table with
           all its annotations in addition to the added on VRS annotations.
           If set to False, obtain a Table with only the VRS annotations.
    :param test: Whether to use a tmp path for analysis of the test Table instead
           of the full v4 Table.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
           Default is "exomes".
    :return: gnomAD v4 VRS VersionedTableResource.
    """
    return VersionedTableResource(
        CURRENT_ANNOTATION_VERSION,
        {
            version: TableResource(
                path=(
                    f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.original_annotations.vrs.ht"
                    if original_annotations
                    else (
                        f"{_annotations_root(version, test, data_type)}/gnomad.{data_type}.v{version}.vrs.ht"
                    )
                )
            )
            for version in ANNOTATION_VERSIONS
        },
    )


def hgdp_tgp_updated_callstats(
    subset: str, test: bool = False
) -> VersionedTableResource:
    """
    Get the HGDP + 1KG/TGP subset updated call stats TableResource.

    :param subset: The subset of the HGDP + 1KG/TGP release to return,
        must be "added", "subtracted", "pop_diff", "join", "v3_release_an",
        "v3_pop_diff_an", or "pre_validity_check".
    :param test: Whether to return the annotation resource for testing purposes.
    :return: MatrixTableResource for specified subset.
    """
    subsets = [
        "added",
        "subtracted",
        "pop_diff",
        "join",
        "v3_release_an",
        "v3_pop_diff_an",
        "pre_validity_check",
    ]
    if subset not in subsets:
        raise ValueError(f"Operation must be one of {subsets}")

    return VersionedTableResource(
        default_version=CURRENT_HGDP_TGP_RELEASE,
        versions={
            release: TableResource(
                f"{_annotations_root(release, test, 'genomes')}/gnomad.genomes.v{release}.hgdp_1kg_subset_updated_callstats_{subset}.ht"
            )
            for release in HGDP_TGP_RELEASES
        },
    )


def get_split_vds(
    version: str = CURRENT_VERSION,
    data_type: str = "exomes",
    test: bool = False,
) -> VariantDatasetResource:
    """
    Get the gnomAD v4 split VDS.

    This is a temporary resource that will be removed once the split VDS is no longer
    needed. Given the uncertainies around frequency calculation runtimes, we cannot
    store it in gnomad-tmp but this needs to be deleted once frequency work is complete.

    :param version: Version of annotation path to return.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
           Default is "exomes".
    :param test: Whether to use a tmp path for analysis of the test Table instead
           of the full v4 Table.
    :return: gnomAD v4 VariantDatasetResource.
    """
    return VariantDatasetResource(
        f"{_annotations_root(version, test, data_type)}/temp/gnomad.{data_type}.v{version}.split_multi.vds"
    )
