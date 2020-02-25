V3_RELEASE_VERSION = 'r3.0'

qual_hists_json_path = 'gs://gnomad/release/3.0/json/gnomad.genomes.r3.json'


# TODO: Remove if not used after all python files are in
# internal_ht_path = 'gs://gnomad/release/3.0/ht/gnomad.genomes.r3.0.nested.no_subsets.sites.ht'


def release_ht_path(data_type: str = 'genomes', release_tag: str = V3_RELEASE_VERSION, public=True):
    '''
    Fetch filepath for release (variant-only) Hail Tables
    :param str data_type: 'exomes' or 'genomes'
    :param str release_tag: String describing release version for use in output filepaths
    :param bool nested: If True, fetch Table in which variant annotations (e.g., freq, popmax, faf, and age histograms)
        are in array format ("nested"); if False, fetch Table in which nested variant annotations are unfurled
    :param bool with_subsets: If True, fetch Table in which all gnomAD subsets are present; if False, fetch Table containing
        only gnomAD annotations
    :param bool temp: If True, fetch Table in which nested variant annotations are unfurled but listed under 'info' rather
        than at the top level; used for sanity-checking sites
    :return: Filepath for desired Hail Table
    :rtype: str
    '''
    release = release_tag.lstrip('r')
    if public:
        return f'gs://gnomad-public/release/{release}/ht/{data_type}/gnomad.{data_type}.{release_tag}.sites.ht'
    else:
        f'gs://gnomad/release/{release}/ht/gnomad.{data_type}.{release_tag}.sites.ht'
