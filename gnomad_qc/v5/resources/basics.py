"""Script containing generic resources."""

import logging
from typing import List, Optional, Set, Union

import hail as hl
from gnomad.resources.resource_utils import (
    MatrixTableResource,
    VariantDatasetResource,
    VersionedVariantDatasetResource,
)
from gnomad.utils.file_utils import file_exists
from hail.utils import new_temp_file

from gnomad_qc.v4.resources.basics import _split_and_filter_variant_data_for_loading
from gnomad_qc.v5.resources.constants import (
    AOU_GENOMIC_METRICS_PATH,
    AOU_LOW_QUALITY_PATH,
    AOU_WGS_BUCKET,
    CURRENT_AOU_VERSION,
    CURRENT_VERSION,
    GNOMAD_TMP_BUCKET,
    WORKSPACE_BUCKET,
)
from gnomad_qc.v5.resources.meta import (
    failing_metrics_samples,
    low_quality_samples,
    samples_to_exclude,
)

logger = logging.getLogger("basic_resources")
logger.setLevel(logging.INFO)

# v5 DRAGEN TGP test VDS.
dragen_tgp_vds = VariantDatasetResource(
    "gs://gnomad/v5.0/testing/genomes/dragen_tgp_v5.0_test.vds"
)

_aou_genotypes = {
    "8": VariantDatasetResource(
        "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
    )
}

aou_genotypes = VersionedVariantDatasetResource(
    CURRENT_AOU_VERSION,
    _aou_genotypes,
)

# This test VDS was initially created to count the number of singletons in AoU and compare the results to those
# provided in the sample QC metrics. It contains 10 randomly chosen samples based on if they
# `fail_singleton_residual` (5 passed, 5 failed) from the AoU sample QC data and is also used for testing our code.
# 'hl.vds.filter_samples()' is used to filter the VDS to these samples, hence it includes all the variants present in
# these samples on all the contigs as the original VDS.
aou_test_dataset = VariantDatasetResource(
    f"gs://{WORKSPACE_BUCKET}/v5.0/hard_filtering/10sample_for_singleton_test.vds"
)


def qc_temp_prefix(
    version: str = CURRENT_VERSION, environment: str = "dataproc"
) -> str:
    """
    Return path to temporary QC bucket.

    :param version: Version of annotation path to return.
    :param environment: Compute environment, either 'dataproc' or 'rwb'. Defaults to 'dataproc'.
    :return: Path to bucket with temporary QC data.
    """
    if environment == "rwb":
        env_bucket = f"{WORKSPACE_BUCKET}/tmp"
    elif environment == "dataproc":
        env_bucket = GNOMAD_TMP_BUCKET
    else:
        raise ValueError(
            f"Environment {environment} not recognized. Choose 'rwb' or 'dataproc'."
        )

    return f"gs://{env_bucket}/gnomad.genomes.v{version}.qc_data/"


def get_checkpoint_path(
    name: str,
    version: str = CURRENT_VERSION,
    mt: bool = False,
    environment: str = "dataproc",
) -> str:
    """
    Create a checkpoint path for Table or MatrixTable.

    :param str name: Name of intermediate Table/MatrixTable.
    :param version: Version of annotation path to return.
    :param bool mt: Whether path is for a MatrixTable, default is False.
    :param environment: Compute environment, either 'dataproc' or 'rwb'. Defaults to 'dataproc'.
    :return: Output checkpoint path.
    """
    return f'{qc_temp_prefix(version, environment)}{name}.{"mt" if mt else "ht"}'


def get_logging_path(
    name: str, version: str = CURRENT_VERSION, environment: str = "dataproc"
) -> str:
    """
    Create a path for Hail log files.

    :param name: Name of log file.
    :param version: Version of annotation path to return.
    :param environment: Compute environment, either 'dataproc' or 'rwb'. Defaults to 'dataproc'.
    :return: Output log path.
    """
    return f"{qc_temp_prefix(version, environment)}{name}.log"


def get_aou_vds(
    split: bool = False,
    remove_hard_filtered_samples: bool = True,
    filter_samples: Optional[Union[List[str], hl.Table]] = None,
    test: bool = False,
    filter_partitions: Optional[List[int]] = None,
    chrom: Optional[Union[str, List[str], Set[str]]] = None,
    autosomes_only: bool = False,
    sex_chr_only: bool = False,
    filter_variant_ht: Optional[hl.Table] = None,
    filter_intervals: Optional[List[Union[str, hl.tinterval]]] = None,
    split_reference_blocks: bool = True,
    remove_dead_alleles: bool = True,
    entries_to_keep: Optional[List[str]] = None,
    checkpoint_variant_data: bool = False,
    naive_coalesce_partitions: Optional[int] = None,
) -> hl.vds.VariantDataset:
    """
    Load the AOU VDS.

    :param split: Whether to split multi-allelic variants in the VDS. Note: this will perform a split on the VDS
        rather than grab an already split VDS. Default is False.
    :param remove_hard_filtered_samples: Whether to remove samples that failed hard
        filters (only relevant after hard filtering is complete). Default is True.
    :param filter_samples: Optional samples to filter the VDS to. Can be a list of sample IDs or a Table with sample IDs.
    :param test: Whether to load the test VDS instead of the full VDS. The test VDS includes 10 samples selected from the full dataset for testing purposes. Default is False.
    :param filter_partitions: Optional argument to filter the VDS to a list of specific partitions.
    :param chrom: Optional argument to filter the VDS to a specific chromosome(s).
    :param autosomes_only: Whether to include only autosomes. Default is False.
    :param sex_chr_only: Whether to include only sex chromosomes. Default is False.
    :param filter_variant_ht: Optional argument to filter the VDS to a specific set of variants. Only supported when splitting the VDS.
    :param filter_intervals: Optional argument to filter the VDS to specific intervals.
    :param split_reference_blocks: Whether to split the reference data at the edges of the intervals defined by `filter_intervals`. Default is True.
    :param remove_dead_alleles: Whether to remove dead alleles when removing samples. Default is True.
    :param entries_to_keep: Optional list of entries to keep in the variant data. If splitting the VDS, use the global entries (e.g. 'GT') instead of the local entries (e.g. 'LGT') to keep.
    :param checkpoint_variant_data: Whether to checkpoint the variant data MT after splitting and filtering. Default is False.
    :param naive_coalesce_partitions: Optional number of partitions to coalesce the VDS to. Default is None.
    :return: AoU v8 VDS.
    """
    aou_v8_resource = aou_test_dataset if test else aou_genotypes
    vds = aou_v8_resource.vds()

    if isinstance(chrom, str):
        chrom = [chrom]

    if autosomes_only and sex_chr_only:
        raise ValueError(
            "Only one of 'autosomes_only' or 'sex_chr_only' can be set to True."
        )

    # Apply chromosome filtering.
    if sex_chr_only:
        vds = hl.vds.filter_chromosomes(vds, keep=["chrX", "chrY"])
    elif autosomes_only:
        vds = hl.vds.filter_chromosomes(vds, keep_autosomes=True)
    elif chrom and len(chrom) > 0:
        logger.info("Filtering to chromosome(s) %s...", chrom)
        vds = hl.vds.filter_chromosomes(vds, keep=chrom)

    # Count initial number of samples.
    n_samples_before = vds.variant_data.count_cols()

    # Remove samples that should have been excluded from the AoU v8 release
    # and samples with non-XX/XY ploidies.
    hard_filtered_samples_ht = None
    if remove_hard_filtered_samples:
        from gnomad_qc.v5.resources.sample_qc import hard_filtered_samples

        logger.info("Removing hard filtered samples from AoU VDS...")
        hard_filtered_samples_ht = hard_filtered_samples.ht()
    s_to_exclude = list(hl.eval(get_samples_to_exclude(hard_filtered_samples_ht)))
    vds = hl.vds.filter_samples(
        vds, s_to_exclude, keep=False, remove_dead_alleles=remove_dead_alleles
    )

    # Report final sample exclusion count.
    n_samples_after = vds.variant_data.count_cols()
    logger.info("Removed %d samples from VDS.", n_samples_before - n_samples_after)

    # Remove the chr4 site with excessive numbers of alleles (n=22233) to
    # avoid memory issues with `split_multi`.
    logger.info("Dropping excessively multi-allelic site at chr4:12237652...")
    vds = hl.vds.filter_intervals(
        vds,
        [hl.parse_locus_interval("chr4:12237652-12237653", reference_genome="GRCh38")],
        keep=False,
    )

    if filter_samples:
        logger.info(
            "Filtering to %s samples...",
            (
                len(filter_samples)
                if isinstance(filter_samples, list)
                else filter_samples.count()
            ),
        )
        vds = hl.vds.filter_samples(vds, filter_samples)

    if naive_coalesce_partitions:
        vds = hl.vds.VariantDataset(
            vds.reference_data.naive_coalesce(naive_coalesce_partitions),
            vds.variant_data.naive_coalesce(naive_coalesce_partitions),
        )

    # Apply interval filtering.
    if filter_intervals and len(filter_intervals) > 0:
        logger.info("Filtering to %s intervals...", len(filter_intervals))
        if isinstance(filter_intervals[0], str):
            filter_intervals = [
                hl.parse_locus_interval(x, reference_genome="GRCh38")
                for x in filter_intervals
            ]
        vds = hl.vds.filter_intervals(
            vds, filter_intervals, split_reference_blocks=split_reference_blocks
        )

    # Apply partition filtering.
    if filter_partitions and len(filter_partitions) > 0:
        logger.info("Filtering to %s partitions...", len(filter_partitions))
        vds = hl.vds.VariantDataset(
            vds.reference_data._filter_partitions(filter_partitions),
            vds.variant_data._filter_partitions(filter_partitions),
        )

    vmt = vds.variant_data

    if filter_variant_ht is not None and split is False:
        raise ValueError(
            "Filtering to a specific set of variants is only supported when splitting"
            " the VDS."
        )
    if split:
        vmt = _split_and_filter_variant_data_for_loading(
            vmt, filter_variant_ht, entries_to_keep, checkpoint_variant_data
        )
        # Don't need to filter entries again if already done in split function
        entries_filtered_during_split = entries_to_keep is not None
    else:
        entries_filtered_during_split = False

    if entries_to_keep is not None and not entries_filtered_during_split:
        vmt = vmt.select_entries(*entries_to_keep)

    if checkpoint_variant_data:
        vmt = vmt.checkpoint(new_temp_file("vds_loading.variant_data", "mt"))

    vds = hl.vds.VariantDataset(vds.reference_data, vmt)

    return vds


aou_acaf_mt = MatrixTableResource(
    path=f"gs://{AOU_WGS_BUCKET}/acaf_threshold/splitMT/hail.mt"
)
"""
AoU v8 ACAF (Allele Count/Allele Frequency threshold) MatrixTable.

MatrixTable contains only variants with AF > 1% or AC > 100 in any genetic ancestry group.

See https://support.researchallofus.org/hc/en-us/articles/29475228181908-How-the-All-of-Us-Genomic-data-are-organized#01JJK0HH53FX9XQRDQ5HQFZW9B
and https://support.researchallofus.org/hc/en-us/articles/14929793660948-Smaller-Callsets-for-Analyzing-Short-Read-WGS-SNP-Indel-Data-with-Hail-MT-VCF-and-PLINK
for more information.
"""


aou_exome_mt = MatrixTableResource(path=f"gs://{AOU_WGS_BUCKET}/exome/splitMT/hail.mt")
"""
AoU v8 Exome MatrixTable.

MatrixTable contains only variants in exons (with 15 bp padding on either side) as defined by GENCODE v42 basic.

See same links as above (in `acaf_mt`) for more information.
"""


def get_aou_failing_genomic_metrics_samples() -> hl.expr.SetExpression:
    """
    Import AoU genomic metrics and filter to samples that fail specific quality criteria, including low coverage and ambiguous sex ploidy.

    .. note::

        Samples with low mean coverage (<30x), genome coverage (<90% at 20x), All of Us Hereditary Disease Risk gene
        (AoUHDR) coverage (<95% at 20x), or aligned_q30_bases (<8e10) were expected to be excluded from the AoU
        callset. However, some such samples are still present. AoU is preparing to publish a known issue in their
        quality report related to this. This note will be updated with a link once the issue is published.

        In addition, we exclude samples with ambiguous sex ploidy (i.e., not "XX" or "XY") from the callset.

    :return: SetExpression of samples failing coverage filters or with non-XX-XY sex ploidies.
    """
    types = {
        "research_id": hl.tstr,
        "sample_source": hl.tstr,
        "site_id": hl.tstr,
        "sex_at_birth": hl.tstr,
        "dragen_sex_ploidy": hl.tstr,
        "mean_coverage": hl.tfloat,
        "genome_coverage": hl.tfloat,
        "aou_hdr_coverage": hl.tfloat,
        "dragen_contamination": hl.tfloat,
        "aligned_q30_bases": hl.tint64,
        "verify_bam_id2_contamination": hl.tfloat,
        "biosample_collection_date": hl.tstr,
    }
    ht = hl.import_table(f"gs://{AOU_GENOMIC_METRICS_PATH}", types=types).key_by(
        "research_id"
    )

    low_cov_samples = ht.filter(
        (ht.mean_coverage < 30)
        | (ht.genome_coverage < 0.9)
        | (ht.aou_hdr_coverage < 0.95)
        | (ht.aligned_q30_bases < 8e10)  # AoU threshold; no v8 samples failed this.
    ).select()
    logger.info("%s samples with low coverage...", low_cov_samples.count())

    ambiguous_sex_samples = ht.filter(
        (ht.dragen_sex_ploidy != "XX") & (ht.dragen_sex_ploidy != "XY")
        | hl.is_missing(ht.dragen_sex_ploidy)
    ).select()
    logger.info(
        "%s samples with ambiguous sex ploidy...", ambiguous_sex_samples.count()
    )

    ht = low_cov_samples.union(ambiguous_sex_samples).distinct()
    return ht.aggregate(hl.agg.collect_as_set(ht.research_id))


def get_samples_to_exclude(
    filter_samples: Optional[Union[List[str], hl.Table]] = None,
    overwrite: bool = False,
) -> hl.expr.SetExpression:
    """
    Get set of AoU sample IDs to exclude.

    .. note::

        If `filter_samples` is a Hail Table, it must contain a field named 's' with sample IDs.

    :param filter_samples: Optional additional samples to remove. Can be a list of sample IDs or a Table with sample IDs.
    :param overwrite: Whether to overwrite the existing `samples_to_exclude` resource. Default is False.
    :return: SetExpression containing IDs of samples to exclude from v5 analysis.
    """
    if not file_exists(samples_to_exclude.path) or overwrite:

        if not file_exists(low_quality_samples.path):
            # Load samples flagged in AoU Known Issues #1.
            logger.info("Removing 3 known low-quality samples (Known Issues #1)...")
            low_quality_ht = hl.import_table(f"gs://{AOU_LOW_QUALITY_PATH}").key_by(
                "research_id"
            )
            low_quality_sample_ids = low_quality_ht.aggregate(
                hl.agg.collect_as_set(low_quality_ht.research_id)
            )
            low_quality_sample_ids.write(low_quality_samples.path)

        if not file_exists(failing_metrics_samples.path):
            # Load and count samples failing genomic metrics filters.
            failing_genomic_metrics_samples_ht = (
                get_aou_failing_genomic_metrics_samples()
            )
            logger.info(
                "Removing %d samples failing genomic QC (low coverage or ambiguous sex)...",
                failing_genomic_metrics_samples_ht.count(),
            )
            failing_genomic_metrics_samples = (
                failing_genomic_metrics_samples_ht.aggregate(
                    hl.agg.collect_as_set(
                        failing_genomic_metrics_samples_ht.research_id
                    )
                )
            )
            failing_genomic_metrics_samples.write(failing_metrics_samples.path)

        # Union all samples to exclude and write out.
        low_quality_sample_ids = hl.experimental.read_expression(
            low_quality_samples.path
        )
        failing_genomic_metrics_samples = hl.experimental.read_expression(
            failing_genomic_metrics_samples.path
        )
        s_to_exclude = low_quality_sample_ids.union(failing_genomic_metrics_samples)
        hl.experimental.write_expression(
            s_to_exclude, samples_to_exclude.path, overwrite=True
        )

    s_to_exclude = hl.experimental.read_expression(samples_to_exclude.path)
    filter_sample_ids = filter_samples or []
    if isinstance(filter_sample_ids, hl.Table):
        if "s" not in filter_samples.row:
            raise ValueError(
                "Hail Table must contain a field named 's' with sample IDs."
            )
        filter_sample_ids = filter_samples.aggregate(
            hl.agg.collect_as_set(filter_samples.s)
        )
    elif not isinstance(filter_sample_ids, list):
        raise ValueError(
            "`filter_samples` must be a list of sample IDs or a Hail Table with sample IDs."
        )
    return s_to_exclude.union(set(filter_sample_ids))


def add_project_prefix_to_sample_collisions(
    t: Union[hl.Table, hl.MatrixTable],
    sample_collisions: hl.Table,
    project: Optional[str] = None,
    sample_id_field: str = "s",
) -> hl.Table:
    """
    Add project prefix to sample IDs that exist in multiple projects.

    :param t: Table/MatrixTable to add project prefix to sample IDs.
    :param sample_collisions: Table of sample IDs that exist in multiple projects.
    :param project: Optional project name to prepend to sample collisions. If not set, will use 'ht.project' annotation. Default is None.
    :param sample_id_field: Field name for sample IDs in the table.
    :return: Table with project prefix added to sample IDs.
    """
    logger.info(
        "Adding project prefix to sample IDs that exists in multiple projects..."
    )
    collisions = sample_collisions.aggregate(hl.agg.collect_as_set(sample_collisions.s))

    if project:
        prefix_expr = hl.literal(project)
    else:
        try:
            prefix_expr = t["project"]
        except KeyError:
            raise ValueError(
                "No project name provided and 't' does not contain a 'project' field."
            )

    renaming_expr = {
        f"{sample_id_field}": hl.if_else(
            hl.literal(collisions).contains(t[sample_id_field]),
            hl.delimit([prefix_expr, t[sample_id_field]], "_"),
            t[sample_id_field],
        )
    }
    if isinstance(t, hl.Table):
        return t.key_by(**renaming_expr)

    logger.warning("Input is a MatrixTable, rekeying columns...")
    return t.key_cols_by(**renaming_expr)
