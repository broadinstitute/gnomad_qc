"""Script containing release related resources."""

import logging

from gnomad.resources.grch38.gnomad import all_sites_an, coverage, public_release
from gnomad.resources.resource_utils import (
    DataException,
    TableResource,
    VersionedTableResource,
)
from gnomad.utils.file_utils import file_exists

from gnomad_qc.v5.resources.basics import qc_temp_prefix
from gnomad_qc.v5.resources.constants import (
    ALL_SITES_AN_RELEASES,
    COVERAGE_RELEASES,
    CURRENT_ALL_SITES_AN_RELEASE,
    CURRENT_COVERAGE_RELEASE,
    CURRENT_RELEASE,
    GNOMAD_BUCKET,
    RELEASES,
    WORKSPACE_BUCKET,
)

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("release_resources")
logger.setLevel(logging.INFO)

FREQUENCY_README = """
The 'freq' row annotation is an array that contains allele frequency information. Each element of the array is a struct that contains the alternate allele count (AC), alternate allele frequency (AF), total number of alleles (AN), and number of homozygous alternate individuals (homozygote_count) for a specific sample grouping.

Use the 'freq_index_dict' global annotation to retrieve frequency information for a specific group of samples from the 'freq' array. This global annotation is a dictionary keyed by sample
grouping combinations whose values are the combination's index in the 'freq' array.

The available keys combinations for the 'freq_index_dict' are as follows:

group, e.g. “adj”, “raw”
sex_group, e.g. “XX_adj”
subset_group, e.g. “non_aou_raw”
gen-anc_group, e.g. “afr_adj”
gen-anc_sex_group, e.g. “ami_XX_adj”
subset_gen-anc_group, e.g. “non_aou_sas_adj”
subset_gen-anc_group, e.g. “non_aou_XY_adj”
subset_gen-anc_sex_group, e.g. “non_aou_mid_XX_adj”,
downsampling_group_gen-anc, e.g. “1373_afr_adj”

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
    data_type: str = "genomes",
    extension: str = "ht",
    environment: str = "rwb",
) -> str:
    """
    Get root path to the release files.

    :param version: Version of release path to return.
    :param test: Whether to use a tmp path for testing.
    :param data_type: Data type of annotation resource. e.g. "exomes" or "genomes".
        Default is "exomes".
    :param extension: File extension of release file. Default is "ht".
    :param environment: Environment to use. Default is "rwb". Must be "rwb" for AoU.
    :return: Root path of the release files.
    """
    path_suffix = f"release/{extension}/{data_type}"
    base_bucket = WORKSPACE_BUCKET if environment == "rwb" else GNOMAD_BUCKET
    if test:
        return (
            f"{qc_temp_prefix(version=version, environment=environment)}{path_suffix}"
        )
    return f"gs://{base_bucket}/v{version}/{path_suffix}"


def release_coverage_path(
    release_version: str = CURRENT_RELEASE,
    public: bool = False,
    test: bool = False,
    coverage_type: str = "coverage",
    environment: str = "rwb",
) -> str:
    """
    Fetch filepath for v5 (AoU + gnomAD v4 genomes) all sites coverage or allele number release Table.

    :param release_version: Release version.
    :param public: Determines whether release coverage Table is read from public (True) or
        private (False) bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param coverage_type: 'coverage' or 'allele_number'. Default is 'coverage'.
    :param environment: Environment to use. Default is "rwb".
    :return: File path for desired coverage Hail Table.
    """
    assert coverage_type in [
        "coverage",
        "allele_number",
    ], "coverage_type must be either 'coverage' or 'allele_number'"

    if public:
        if test:
            raise ValueError("Cannot use test=True with public=True!")
        try:
            if coverage_type == "coverage":
                cov = coverage("genomes")
            else:
                cov = all_sites_an("genomes")
            if release_version in cov.versions:
                path = cov.versions[release_version].path
            else:
                path = None
        except DataException:
            path = None
        if path is None:
            raise ValueError(
                f"No public {coverage_type} Table found for genomes and release {release_version}."
            )
        return path
    else:
        return f"{_release_root(release_version, test=test, environment=environment)}/gnomad.genomes.v{release_version}.{coverage_type}.ht"


def release_coverage_tsv_path(
    release_version: str = CURRENT_COVERAGE_RELEASE["genomes"],
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch path to coverage TSV file.

    :param release_version: Release version. Default is CURRENT_COVERAGE_RELEASE["genomes"].
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Coverage TSV path.
    """
    return f"{_release_root(release_version, test=test, extension='tsv', environment=environment)}/gnomad.genomes.v{release_version}.coverage.tsv.bgz"


def release_all_sites_an_tsv_path(
    release_version: str = None,
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch path to all sites AN TSV file.

    :param release_version: Release version. Default is None.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: All sites AN TSV path.
    """
    release_version = (
        release_version
        if release_version is not None
        else CURRENT_ALL_SITES_AN_RELEASE["genomes"]
    )
    return f"{_release_root(release_version, test=test, extension='tsv', environment=environment)}/gnomad.genomes.v{release_version}.allele_number.tsv.bgz"


def release_coverage(
    public: bool = False,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for coverage release Table.

    :param public: Determines whether release coverage Table is read from public (True) or
        private (False) bucket. Default is False.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Coverage release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_COVERAGE_RELEASE["genomes"],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    release_version=release,
                    public=public,
                    test=test,
                    environment=environment,
                )
            )
            for release in COVERAGE_RELEASES["genomes"]
        },
    )


def release_all_sites_an(
    public: bool = False,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for all sites allele number release Table.

    :param public: Determines whether release allele number Table is read from public or
        private bucket. Default is private.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: All sites allele number release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_ALL_SITES_AN_RELEASE["genomes"],
        versions={
            release: TableResource(
                path=release_coverage_path(
                    release_version=release,
                    public=public,
                    test=test,
                    coverage_type="allele_number",
                    environment=environment,
                )
            )
            for release in ALL_SITES_AN_RELEASES["genomes"]
        },
    )


def included_datasets_json_path(
    test: bool = False,
    release_version: str = CURRENT_RELEASE,
    environment: str = "rwb",
) -> str:
    """
    Fetch filepath for the JSON containing all datasets used in the release.

    :param test: Whether to use a tmp path for testing. Default is False.
    :param release_version: Release version. Default is CURRENT RELEASE.
    :param environment: Environment to use. Default is "rwb".
    :return: File path for release versions included datasets JSON
    """
    return f"{_release_root(release_version, test=test, data_type='genomes', extension='json', environment=environment)}/gnomad.genomes.v{release_version}.included_datasets.json"


def release_ht_path(
    release_version: str = CURRENT_RELEASE,
    public: bool = True,
    test: bool = False,
    environment: str = "rwb",
) -> str:
    """
    Fetch filepath for release (variant-only) Hail Tables.

    :param release_version: Release version. Default is CURRENT_RELEASE.
    :param public: Whether release sites Table path returned is from public instead of private
        bucket. Default is True.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: File path for desired release Hail Table.
    """
    if public:
        if file_exists(
            public_release(data_type="genomes").versions[release_version].path
        ):
            return public_release(data_type="genomes").versions[release_version].path
        else:
            return f"gs://gnomad-public-requester-pays/release/{release_version}/ht/genomes/gnomad.genomes.v{release_version}.sites.ht"
    else:
        # Note: This function was not used to write out the v4.0 joint release HT. That
        # is gs://gnomad/release/4.0/ht/joint/gnomad.joint.v4.0.faf.filtered.ht
        return f"{_release_root(version=release_version, test=test, data_type='genomes', environment=environment)}/gnomad.genomes.v{release_version}.sites.ht"


def release_sites(
    public: bool = True,
    test: bool = False,
    environment: str = "rwb",
) -> VersionedTableResource:
    """
    Retrieve versioned resource for sites-only release Table.

    :param public: Whether release sites Table path returned is from public or private
        bucket. Default is True.
    :param test: Whether to use a tmp path for testing. Default is False.
    :param environment: Environment to use. Default is "rwb".
    :return: Sites-only release Table.
    """
    return VersionedTableResource(
        default_version=CURRENT_RELEASE,
        versions={
            release: TableResource(
                path=release_ht_path(
                    release_version=release,
                    public=public,
                    test=test,
                    environment=environment,
                )
            )
            for release in RELEASES
        },
    )
