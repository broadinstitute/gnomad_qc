from gnomad.resources.resource_utils import TableResource
from gnomad_qc.v3.resources import CURRENT_RELEASE


def get_annotations_root(version: str = CURRENT_RELEASE) -> str:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return f"gs://gnomad/annotations/hail-0.2/ht/genomes_v{version}"


def get_info(version: str = CURRENT_RELEASE, split: bool = True) -> TableResource:
    """
    Gets the gnomAD v3 info TableResource

    :param version: Version of sample QC path to return
    :param split: Whether to return the split or multi-allelic version of the resource
    :return: gnomAD v3 info TableResource
    """
    path = "{}/gnomad_genomes_v{}_info{}.ht".format(
        get_annotations_root(version), version, ".split" if split else ""
    )
    return TableResource(path)


def get_filters(
    model_id: str,
    version: str = CURRENT_RELEASE,
    split: bool = True,
    finalized: bool = False,
) -> TableResource:
    """
    Gets the specified filtering annotation resource.

    :param model_id: Filtering model id
    :param version: Version of sample QC path to return
    :param split: Split or multi-allelic version of the filtering file
    :param finalized:
    :return: Filtering annotation file
    """
    path = f"{get_annotations_root(version)}/filtering/{model_id}{'.finalized' if finalized else ''}{'.split' if split else ''}.ht"
    return TableResource(path)


def last_END_position(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}_last_END_positions.ht"
    )


def freq(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}.frequencies.ht"
    )


def qual_hist(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}.qual_hists.ht"
    )


def vep(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}_vep.ht"
    )


def info_vcf_path(version: str = CURRENT_RELEASE) -> str:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return f"{get_annotations_root(version)}/gnomad_genomes_v{version}_info.vcf.bgz"


def qc_ac(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(f"{get_annotations_root(version)}/gnomad_genomes_qc_ac.ht")


def fam_stats(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_qc_fam_stats.ht"
    )


def allele_data(version: str = CURRENT_RELEASE) -> TableResource:
    """

    :param version: Version of sample QC path to return
    :return:
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_qc_allele_data.ht"
    )
