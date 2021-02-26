from typing import Optional

from gnomad.resources.grch38.gnomad import SUBSETS
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)

from gnomad_qc.v3.resources.constants import CURRENT_RELEASE, RELEASES


def _annotations_root(version: str = CURRENT_RELEASE) -> str:
    """
    Get root path to the variant annotation files

    :param version: Version of annotation path to return
    :return: root path of the variant annotation files
    """
    return f"gs://gnomad/annotations/hail-0.2/ht/genomes_v{version}"


def get_info(split: bool = True) -> VersionedTableResource:
    """
    Gets the gnomAD v3 info TableResource

    :param version: Version of annotation path to return
    :param split: Whether to return the split or multi-allelic version of the resource
    :return: gnomAD v3 info VersionedTableResource
    """

    return VersionedTableResource(
        CURRENT_RELEASE,
        {
            release: TableResource(
                path="{}/gnomad_genomes_v{}_info{}.ht".format(
                    _annotations_root(release), release, ".split" if split else ""
                )
            )
            for release in RELEASES
        },
    )


def get_filters(
    model_id: str,
    split: bool = True,
    finalized: bool = False,
) -> VersionedTableResource:
    """
    Gets the specified VQSR filtering annotation resource.

    :param model_id: VQSR filtering model id
    :param split: Split or multi-allelic version of the filtering file
    :param finalized: Whether to return the raw VQSR table or the finalized VQSR table representing determined cutoffs
    :return: VQSR filtering annotation file
    """
    return VersionedTableResource(
        CURRENT_RELEASE,
        {
            release: TableResource(
                "{}/filtering/{}{}{}.ht".format(
                    _annotations_root(release),
                    model_id,
                    ".finalized" if finalized else "",
                    ".split" if split else "",
                )
            )
            for release in RELEASES
        },
    )


def info_vcf_path(version: str = CURRENT_RELEASE) -> str:
    """
    Path to sites VCF (input information for running VQSR)

    :param version: Version of annotation path to return
    :return: String for the path to the info VCF
    """
    return f"{_annotations_root(version)}/gnomad_genomes_v{version}_info.vcf.bgz"


last_END_position = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_last_END_positions.ht"
        )
        for release in RELEASES
    },
)

freq = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}.frequencies.ht"
        )
        for release in RELEASES
    },
)

qual_hist = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}.qual_hists.ht"
        )
        for release in RELEASES
    },
)

vep = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_vep.ht"
        )
        for release in RELEASES
    },
)

qc_ac = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(f"{_annotations_root(release)}/gnomad_genomes_qc_ac.ht")
        for release in RELEASES
    },
)

fam_stats = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_qc_fam_stats.ht"
        )
        for release in RELEASES
    },
)

allele_data = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_qc_allele_data.ht"
        )
        for release in RELEASES
    },
)

analyst_annotations = VersionedTableResource(
    CURRENT_RELEASE,
    {
        release: TableResource(
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_in_silico_predictors.ht"
        )
        for release in RELEASES
        if release != "3.0"
    },
)


def get_freq(
    version: str = CURRENT_RELEASE, subset: Optional[str] = None
) -> VersionedTableResource:
    """
    Get the frequency annotation table for a specified release.

    :param version: Version of annotation path to return
    :param subset: One of the official subsets of the specified release (e.g., non_neuro, non_cancer, controls_and_biobanks)
    :return: Hail Table containing subset or overall cohort frequency annotations
    """
    if version == "3" and subset:
        raise DataException("Subsets of gnomAD v3 do not exist")

    if subset and subset not in SUBSETS:
        raise DataException(
            f"{subset} subset is not one of the following official subsets: {SUBSETS}"
        )

    return VersionedTableResource(
        version,
        {
            release: TableResource(
                f"{_annotations_root(release)}/gnomad_genomes_v{release}.frequencies{'.' + subset if subset else ''}.ht"
            )
            for release in RELEASES
        },
    )
