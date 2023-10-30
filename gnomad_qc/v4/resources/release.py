"""Script containing release related resources."""
import logging
from typing import Optional

from gnomad.resources.grch38.gnomad import coverage, public_release
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v4.resources.constants import (
    COMBINED_FAF_RELEASES,
    COVERAGE_RELEASES,
    CURRENT_COMBINED_FAF_RELEASE,
    CURRENT_COVERAGE_RELEASE,
    CURRENT_RELEASE,
    RELEASES,
)

FREQUENCY_README = """
The 'freq' row annotation is an array that contains allele frequency information. Each element of the array is a struct that contains the alternate allele count (AC), alternate allele frequency (AF), total number of alleles (AN), and number of homozygous alternate individuals (homozygote_count) for a specific sample grouping.

Use the 'freq_index_dict' global annotation to retrieve frequency information for a specific group of samples from the 'freq' array. This global annotation is a dictionary keyed by sample
grouping combinations whose values are the combination's index in the 'freq' array.

The available keys combinations for the 'freq_index_dict' are as follows:

group, e.g. “adj”, “raw”
sex_group, e.g. “XX_adj”
subset_group, e.g. “non_ukb_raw”
gen-anc_group, e.g. “afr_adj”
gen-anc_sex_group, e.g. “ami_XX_adj”{}
subset_gen-anc_group, e.g. “non_ukb_sas_adj”
subset_gen-anc_group, e.g. “non_ukb_XY_adj”
subset_gen-anc_sex_group, e.g. “non_ukb_mid_XX_adj”,

The example below shows how to access the entry of the high quality genotypes
(group: adj) of XX individuals (sex: XX) labeled as AFR (gen_anc: AFR)
in the HT:

    # Use the key 'afr-XX-adj' to retrieve the index of this groups frequency data in 'freq'
    ht = ht.annotate(afr_XX_freq=ht.freq[ht.freq_index_dict['afr-XX-adj']])

The above example will retrieve the entire frequency struct for each variant. To grab a
certain statistic, such as AC, specify the statistic after the value:

    ht = ht.annotate(afr_XX_AC=ht.freq[ht.freq_index_dict['afr-XX-adj']].AC)

This same approach can be applied to the filtering allele frequency (FAF) array, 'faf',
by using the 'faf_index_dict'. For more information, please visit the FAQ page:
https://gnomad.broadinstitute.org/help#technical-details:~:text=How%20do%20I%20access%20the%20gnomAD%20Hail%20Table%20frequency%20annotation%3F
"""

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("release_resources")
logger.setLevel(logging.INFO)


def _release_root(
    version: str = CURRENT_RELEASE,
    test: bool = False,
    data_type: str = "exomes",
    extension: str = "ht",
) -> str:
    """
    Get root path to the release files.

    :param version: Version of release path to return.
    :param test: Whether to use a tmp path for testing.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
        Default is "exomes".
    :param extension: File extension of release file. Default is "ht".
    :return: Root path of the release files.
    """
    return (
        f"gs://gnomad-tmp/gnomad_v{version}_testing/release/{extension}/{data_type}"
        if test
        else f"gs://gnomad/release/{version}/{extension}/{data_type}"
    )


def annotation_hists_params_path(
    release_version: str = CURRENT_RELEASE,
    data_type: str = "exomes",
) -> str:
    """
    Return path to file containing dictionary of parameters for site metric histograms.

    The keys of the dictionary are the names of the site quality metrics
    while the values are: [lower bound, upper bound, number  of bins].
    For example, "InbreedingCoeff": [-0.25, 0.25, 50].

    :param release_version: Release version. Defaults to CURRENT RELEASE.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
        Default is "exomes".
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Path to file with annotation histograms
    """
    return (
        f"{_release_root(version=release_version, data_type=data_type, extension='json')}/gnomad.{data_type}.v{release_version}_annotation_hist_params.json"
    )


def qual_hists_json_path(
    release_version: str = CURRENT_RELEASE,
    data_type: str = "exomes",
    test: bool = False,
) -> str:
    """
    Fetch filepath for qual histograms JSON.

    :param release_version: Release version. Defaults to CURRENT RELEASE
    :param data_type: Data type 'exomes' or 'genomes'. Default is 'exomes'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: File path for histogram JSON
    """
    return (
        f"{_release_root(release_version, test, data_type, extension='json')}/gnomad.{data_type}.v{release_version}_qual_hists.json"
    )


def get_combined_faf_release(
    test: bool = False, filtered: bool = False
) -> VersionedTableResource:
    """
    Retrieve versioned resource for the combined genome + exome FAF release Table.

    :param test: Whether to use a tmp path for testing. Default is False.
    :param filtered: Whether to use the filtered FAF release. Default is False.
    :return: Combined genome + exome FAF release VersionedTableResource.
    """
    return VersionedTableResource(
        default_version=CURRENT_COMBINED_FAF_RELEASE,
        versions={
            release: TableResource(
                path=(
                    f"gs://gnomad{'-tmp' if test else ''}/release/{release}/ht/joint/gnomad.joint.v{release}.faf{'.filtered' if filtered else ''}.ht"
                )
            )
            for release in COMBINED_FAF_RELEASES
        },
    )


def release_ht_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = False,
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables.

    :param data_type: Data type of release resource to return. Should be one of
        'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version. Default is CURRENT_RELEASE.
    :param public: Whether release sites Table path returned is from public or private
        bucket. Default is False.
    :return: File path for desired release Hail Table.
    """
    if public:
        if file_exists(public_release(data_type).versions[release_version].path):
            return public_release(data_type).versions[release_version].path
        else:
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.filtered_combined_faf.ht"
    else:
        return f"gs://gnomad/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.sites.filtered_combined_faf.ht"


def release_sites(
    data_type: str = "exomes", public: bool = False
) -> VersionedTableResource:
    """
    Retrieve versioned resource for sites-only release Table.

    :param data_type: Data type of release resource to return. Should be one of
        'exomes' or 'genomes'. Default is 'exomes'.
    :param public: Whether release sites Table path returned is from public or private
        bucket. Default is False.
    :return: Sites-only release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_ht_path(
                    data_type=data_type,
                    release_version=release,
                    public=public,
                )
            )
            for release in RELEASES
        },
    )


def release_vcf_path(
    release_version: Optional[str] = None,
    contig: Optional[str] = None,
) -> str:
    """
    Fetch bucket for release (sites-only) VCFs.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used.
    :param contig: String containing the name of the desired reference contig. Defaults to the full (all contigs) sites VCF path
        sites VCF path
    :return: Filepath for the desired VCF
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    if contig:
        return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}.sites.{contig}.vcf.bgz"
    else:
        # If contig is None, return path to sharded vcf bucket.
        # NOTE: need to add .bgz or else hail will not bgzip shards.
        return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}.sites.vcf.bgz"


def release_header_path(release_version: Optional[str] = None) -> str:
    """
    Fetch path to pickle file containing VCF header dictionary.

    :param release_version: Release version. When no release_version is supplied CURRENT_RELEASE is used
    :return: Filepath for header dictionary pickle
    """
    if release_version is None:
        release_version = CURRENT_RELEASE

    return f"gs://gnomad/release/{release_version}/vcf/exomes/gnomad.exomes.v{release_version}_header_dict.pickle"


def append_to_vcf_header_path(
    subset: str, release_version: str = CURRENT_RELEASE
) -> str:
    """
    Fetch path to TSV file containing extra fields to append to VCF header.

    Extra fields are VEP and dbSNP versions.

    :param subset: One of the possible release subsets
    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: Filepath for extra fields TSV file
    """
    return (
        f"gs://gnomad/release/{release_version}/vcf/exomes/extra_fields_for_header{f'_{subset}' if subset else ''}.tsv"
    )


def release_coverage_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
    test: bool = False,
) -> str:
    """
    Fetch filepath for coverage release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version.
    :param public: Determines whether release coverage Table is read from public or
        private bucket. Default is public.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: File path for desired coverage Hail Table.
    """
    if public:
        if test:
            raise ValueError("Cannot use test=True with public=True!")
        try:
            cov = coverage(data_type)
            if release_version in cov.versions:
                path = cov.versions[release_version].path
            else:
                path = None
        except DataException:
            path = None
        if path is None:
            logger.warning(
                "No public coverage Table found for data_type %s and release %s. "
                "Using 'gs://gnomad-public-requester-pays' path.",
                data_type,
                release_version,
            )
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/{data_type}/gnomad.{data_type}.v{release_version}.coverage.ht"
        else:
            return path
    else:
        return (
            f"{_release_root(release_version, test=test, data_type=data_type)}/gnomad.{data_type}.v{release_version}.coverage.ht"
        )


def release_coverage_tsv_path(
    data_type: str = "exomes",
    release_version: str = CURRENT_COVERAGE_RELEASE["exomes"],
    test: bool = False,
) -> str:
    """
    Fetch path to coverage TSV file.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param release_version: Release version. Default is CURRENT_COVERAGE_RELEASE["exomes"].
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Coverage TSV path.
    """
    return (
        f"{_release_root(release_version, test=test, data_type=data_type, extension='tsv')}/gnomad.{data_type}.v{release_version}.coverage.tsv.bgz"
    )


def release_coverage(
    data_type: str = "exomes", public: bool = False, test: bool = False
) -> VersionedTableResource:
    """
    Retrieve versioned resource for coverage release Table.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param public: Determines whether release coverage Table is read from public or
        private bucket. Default is private.
    :param test: Whether to use a tmp path for testing. Default is False.
    :return: Coverage release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_COVERAGE_RELEASE[data_type],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    data_type=data_type,
                    release_version=release,
                    public=public,
                    test=test,
                )
            )
            for release in COVERAGE_RELEASES[data_type]
        },
    )


def included_datasets_json_path(
    data_type: str = "exomes",
    test: bool = False,
    release_version: str = CURRENT_RELEASE,
) -> str:
    """
    Fetch filepath for the JSON containing all datasets used in the release.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param release_version: Release version. Defaults to CURRENT RELEASE
    :return: File path for release versions included datasets JSON
    """
    return (
        f"{_release_root(release_version, test=test, data_type=data_type, extension='json')}/gnomad.exomes.v{release_version}.included_datasets.filtered_combined_faf.json"
    )


def validated_release_ht(
    test: bool = False,
    data_type: str = "exomes",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for validated sites-only release Table.

    :param test: Whether to use a tmp path for testing. Default is False.
    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :return: Validated release Table
    """
    return VersionedTableResource(
        CURRENT_RELEASE,
        {
            version: TableResource(
                path=(
                    f"{_release_root(version, data_type=data_type, test=test)}/gnomad.{data_type}.v{version}.validated_release.filtered_combined_faf.ht"
                )
            )
            for version in RELEASES
        },
    )


def get_freq_array_readme(data_type: str = "exomes") -> str:
    """
    Fetch README for freq array for the specified data_type.

    :param data_type: 'exomes' or 'genomes'. Default is 'exomes'.
    :return: README for freq array.
    """
    if data_type == "exomes":
        return FREQUENCY_README.format(
            "\ndownsampling_group_gen-anc, e.g."
            " “200_eas_adj”,\ndownsampling_group_gen-anc, e.g. “non_ukb_218035_eas_adj”"
        )
    else:
        return FREQUENCY_README.format("\n")
