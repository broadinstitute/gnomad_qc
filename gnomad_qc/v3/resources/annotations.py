from typing import Optional

from gnomad.resources.grch38.gnomad import SUBSETS
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v3.resources.constants import (
    CURRENT_INSILICO_ANNOTATION_VERSION,
    DEFAULT_CURRENT_VERSION,
    DEFAULT_VERSIONS,
    INSILICO_ANNOTATION_VERSIONS,
)


def _annotations_root(version: str = DEFAULT_CURRENT_VERSION) -> str:
    """
    Get root path to the variant annotation files

    :param version: Version of annotation path to return
    :return: root path of the variant annotation files
    """
    return f"gs://gnomad/annotations/hail-0.2/ht/genomes_v{version}"


def get_info(split: bool = True) -> VersionedTableResource:
    """
    Gets the gnomAD v3 info TableResource

    :param split: Whether to return the split or multi-allelic version of the resource
    :return: gnomAD v3 info VersionedTableResource
    """

    return VersionedTableResource(
        DEFAULT_CURRENT_VERSION,
        {
            release: TableResource(
                path="{}/gnomad_genomes_v{}_info{}.ht".format(
                    _annotations_root(release), release, ".split" if split else ""
                )
            )
            for release in DEFAULT_VERSIONS
        },
    )


def get_vqsr_filters(
    model_id: str, split: bool = True, finalized: bool = False,
) -> VersionedTableResource:
    """
    Gets the specified VQSR filtering annotation resource.

    :param model_id: VQSR filtering model id
    :param split: Split or multi-allelic version of the filtering file
    :param finalized: Whether to return the raw VQSR table or the finalized VQSR table representing determined cutoffs
    :return: VQSR filtering annotation file
    """
    return VersionedTableResource(
        DEFAULT_CURRENT_VERSION,
        {
            release: TableResource(
                "{}/filtering/{}{}{}.ht".format(
                    _annotations_root(release),
                    model_id,
                    ".finalized" if finalized else "",
                    ".split" if split else "",
                )
            )
            for release in DEFAULT_VERSIONS
        },
    )


def info_vcf_path(version: str = DEFAULT_CURRENT_VERSION) -> str:
    """
    Path to sites VCF (input information for running VQSR)

    :param version: Version of annotation path to return
    :return: String for the path to the info VCF
    """
    return f"{_annotations_root(version)}/gnomad_genomes_v{version}_info.vcf.bgz"


def get_transmitted_singleton_vcf_path(
    adj: bool = False, version: str = DEFAULT_CURRENT_VERSION
) -> str:
    """
    Provides the path to the transmitted singleton VCF used as input to VQSR

    :param bool adj: Whether to use adj genotypes
    :param version: Version of transmitted singleton VCF path to return
    :return:
    """
    return f'{_annotations_root(version)}/transmitted_singletons_{"adj" if adj else "raw"}.vcf.bgz'


last_END_position = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_last_END_positions.ht"
        )
        for release in DEFAULT_VERSIONS
    },
)

freq = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}.frequencies.ht"
        )
        for release in DEFAULT_VERSIONS
    },
)

qual_hist = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}.qual_hists.ht"
        )
        for release in DEFAULT_VERSIONS
    },
)

vep = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_vep.ht"
        )
        for release in DEFAULT_VERSIONS
    },
)

qc_ac = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(f"{_annotations_root(release)}/gnomad_genomes_qc_ac.ht")
        for release in DEFAULT_VERSIONS
    },
)

fam_stats = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_qc_fam_stats.ht"
        )
        for release in DEFAULT_VERSIONS
    },
)

allele_data = VersionedTableResource(
    DEFAULT_CURRENT_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_qc_allele_data.ht"
        )
        for release in DEFAULT_VERSIONS
    },
)


def get_freq(
    version: str = DEFAULT_CURRENT_VERSION, subset: Optional[str] = None
) -> VersionedTableResource:
    """
    Get the frequency annotation table for a specified release.

    :param version: Version of annotation path to return
    :param subset: One of the official subsets of the specified release (e.g., non_neuro, non_cancer,
        controls_and_biobanks) or a combination of them split by '-'
    :return: Hail Table containing subset or overall cohort frequency annotations
    """
    if version == "3" and subset:
        raise DataException("Subsets of gnomAD v3 do not exist")

    if subset is not None:
        for s in subset.split("-"):
            if s not in SUBSETS:
                raise DataException(
                    f"{subset} subset is not one of the following official subsets: {SUBSETS}"
                )

    return VersionedTableResource(
        version,
        {
            release: TableResource(
                f"{_annotations_root(release)}/gnomad_genomes_v{release}.frequencies{'.' + subset if subset else ''}.ht"
            )
            for release in DEFAULT_VERSIONS
        },
    )


analyst_annotations = VersionedTableResource(
    CURRENT_INSILICO_ANNOTATION_VERSION,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_in_silico_predictors.ht"
        )
        for release in INSILICO_ANNOTATION_VERSIONS
    },
)
