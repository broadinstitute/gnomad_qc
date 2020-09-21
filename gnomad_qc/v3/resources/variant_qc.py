from gnomad.resources.resource_utils import TableResource

VARIANT_QC_ROOT = 'gs://gnomad/variant_qc/genomes_v3'


def get_transmitted_singleton_vcf_path(confidence: str):
    return f'{VARIANT_QC_ROOT}/transmitted_singletons_{confidence}.vcf'


def score_ranking_path(model_id: str, binned: bool):
    return '{}/{}.{}.ht'.format(
        VARIANT_QC_ROOT,
        model_id,
        'binned' if binned else 'rank'
    )


def binned_concordance_path(model_id: str, truth_sample: str):
    return f'{VARIANT_QC_ROOT}/{truth_sample}_{model_id}_binned_concordance.ht'


def get_filtering_model(model_id: str, split: bool = True, finalized: bool = True) -> TableResource:
    """
       Gets the specified filtering annotation resource.

       :param model_id: Filtering model id
       :param split: Split or multi-allelic version of the filtering file
       :return: Filtering annotation file
       """
    path = '{}/filtering/{}{}{}.ht'.format(
        VARIANT_QC_ROOT,
        model_id,
        '.finalized' if finalized else '',
        '.split' if split else ''
    )
    return TableResource(path)
