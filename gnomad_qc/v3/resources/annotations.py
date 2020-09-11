from gnomad.resources.resource_utils import TableResource
from gnomad_qc.v3.resources import CURRENT_RELEASE


def get_annotations_root(version: str = CURRENT_RELEASE) -> str:
    """
    Get root path to the variant  annotation files

    :param version: Version of annotation path to return
    :return: root path of the variant annotation files
    """
    return f"gs://gnomad/annotations/hail-0.2/ht/genomes_v{version}"


def get_info(version: str = CURRENT_RELEASE, split: bool = True) -> TableResource:
    """
    Gets the gnomAD v3 info TableResource

    :param version: Version of annotation path to return
    :param split: Whether to return the split or multi-allelic version of the resource
    :return: gnomAD v3 info TableResource
    """
    path = "{}/gnomad_genomes_v{}_info{}.ht".format(
        get_annotations_root(version), version, ".split" if split else ""
    )
    return TableResource(path)


def get_vqsr_filters(
    model_id: str,
    version: str = CURRENT_RELEASE,
    split: bool = True,
    finalized: bool = False,
) -> TableResource:
    """
    Gets the specified filtering annotation resource.

    :param model_id: Filtering model id
    :param version: Version of annotation path to return
    :param split: Split or multi-allelic version of the filtering file
    :param finalized: Whether to return the raw VQSR table or the finalized VQSR table representing determined cutoffs
    :return: VQSR filtering annotation file
    """
    path = f"{get_annotations_root(version)}/filtering/{model_id}{'.finalized' if finalized else ''}{'.split' if split else ''}.ht"
    return TableResource(path)


def last_END_position(version: str = CURRENT_RELEASE) -> TableResource:
    """
     Retrieves TableResource containing the last END positions annotation

    :param version: Version of annotation path to return
    :return: Table with last_END_position annotation
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}_last_END_positions.ht"
    )


def freq(version: str = CURRENT_RELEASE) -> TableResource:
    """
    Retrieves TableResource containing frequency information

    :param version: Version of annotation path to return
    :return: Table with variant frequency annotations
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}.frequencies.ht"
    )


def qual_hist(version: str = CURRENT_RELEASE) -> TableResource:
    """
    Retrieves TableResource containing quality histograms

    :param version: Version of annotation path to return
    :return: Table with quality histogram information
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}.qual_hists.ht"
    )


def vep(version: str = CURRENT_RELEASE) -> TableResource:
    """
    Retrieves TableResource containing information from Variant Effect Predictor (VEP)

    :param version: Version of annotation path to return
    :return: Table with VEP annotations
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_v{version}_vep.ht"
    )


def info_vcf_path(version: str = CURRENT_RELEASE) -> str:
    """
    Path to sites VCF (input information for running VQSR)

    :param version: Version of annotation path to return
    :return: String for the path to the info VCF
    """
    return f"{get_annotations_root(version)}/gnomad_genomes_v{version}_info.vcf.bgz"


def qc_ac(version: str = CURRENT_RELEASE) -> TableResource:
    """
    Retrieves TableResource containing allele counts per variant

    :param version: Version of annotation path to return
    :return: Table with annotations for variant allele counts
    """
    return TableResource(f"{get_annotations_root(version)}/gnomad_genomes_qc_ac.ht")


def fam_stats(version: str = CURRENT_RELEASE) -> TableResource:
    """
    Retrieves TableResource containing information about trios including statistics about transmission and de novo count

    :param version: Version of annotation path to return
    :return: Table with per variant family statistics
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_qc_fam_stats.ht"
    )


def allele_data(version: str = CURRENT_RELEASE) -> TableResource:
    """
    Retrieves TableResource containing annotations about the alleles

    Includes  the following annotations:

        - variant_type
        - allele_type
        - n_alt_alleles
        - was_mixed
        - has_star

    :param version: Version of annotation path to return
    :return: Table with allele annotations
    """
    return TableResource(
        f"{get_annotations_root(version)}/gnomad_genomes_qc_allele_data.ht"
    )
