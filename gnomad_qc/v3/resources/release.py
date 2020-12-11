from gnomad_qc.v3.resources.constants import CURRENT_RELEASE

def qual_hists_json_path(release_version: str= CURRENT_RELEASE) -> str:
    """Fetch filepath for qual histograms JSON

    :param release_version: Release version. Defualts to CURRENT RELEASE
    :return: File path for histogram JSON
    :rtype: str
    """
    return f"gs://gnomad/release/{release_version}/json/gnomad.genomes.r{release_version}.json"


# TODO: Remove if not used after all python files are in
# internal_ht_path = 'gs://gnomad/release/3.0/ht/gnomad.genomes.r3.0.nested.no_subsets.sites.ht'


def release_ht_path(
    data_type: str = "genomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables
    :param data_type: 'exomes' or 'genomes'
    :param release_version: release version
    :param public: Whether to return the desired
    :return: File path for desired Hail Table
    :rtype: str
    """
    if public:
        return f"gs://gnomad-public/release/{release_version}/ht/{data_type}/gnomad.{data_type}.r{release_version}.sites.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/gnomad.{data_type}.r{release_version}.sites.ht"
