from gnomad.resources.resource_utils import TableResource, VersionedTableResource

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


def get_vqsr_filters(
    model_id: str,
    split: bool = True,
    finalized: bool = False,
) -> VersionedTableResource:
    """
    Gets the specified filtering annotation resource.

    :param model_id: Filtering model id
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
            f"{_annotations_root(release)}/gnomad_genomes_v{release}_analyst_annotations.ht"
        )
        for release in RELEASES if release != "3.0"
    },
)
