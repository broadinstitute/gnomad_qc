"""Script to generate annotations for variant QC on gnomAD v4."""

import argparse
import logging

import hail as hl
from gnomad.utils.annotations import add_variant_type, annotate_allele_info
from gnomad.utils.slack import slack_notifications
from gnomad.utils.sparse_mt import (
    INFO_INT32_SUM_AGG_FIELDS,
    INFO_SUM_AGG_FIELDS,
    default_compute_info,
    get_as_info_expr,
    get_site_info_expr,
    split_info_annotation,
    split_lowqual_annotation,
)
from gnomad.utils.vcf import adjust_vcf_incompatible_types

from gnomad_qc.resource_utils import (
    PipelineResourceCollection,
    PipelineStepResourceCollection,
)
from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.annotations import get_info, info_vcf_path
from gnomad_qc.v4.resources.basics import get_gnomad_v4_vds

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def split_info(info_ht: hl.Table) -> hl.Table:
    """
    Generate an info Table with split multi-allelic sites from the multi-allelic info Table.

    .. note::

        The split of multi-allelic sites before splitting the 'info' annotation is done
        within `annotate_allele_info` in gnomad_methods so that allele info is also
        annotated on the returned Table.

    :param info_ht: Info Table with unsplit multi-allelics.
    :return: Info Table with split multi-allelics.
    """
    info_ht = annotate_allele_info(info_ht)
    info_ht = info_ht.annotate(
        info=info_ht.info.annotate(
            **split_info_annotation(info_ht.info, info_ht.a_index),
        ),
        AS_lowqual=split_lowqual_annotation(info_ht.AS_lowqual, info_ht.a_index),
    )

    return info_ht


def get_variant_qc_annotation_resources(
    test: bool, overwrite: bool
) -> PipelineResourceCollection:
    """
    Get PipelineResourceCollection for all resources needed in the variant QC annotation pipeline.

    :param test: Whether to gather all resources for testing.
    :param overwrite: Whether to overwrite resources if they exist.
    :return: PipelineResourceCollection containing resources for all steps of the
        variant QC annotation pipeline.
    """
    # Initialize variant QC annotation pipeline resource collection.
    ann_pipeline = PipelineResourceCollection(
        pipeline_name="variant_qc_annotation",
        overwrite=overwrite,
    )

    # Create resource collection for each step of the variant QC annotation pipeline.
    compute_info = PipelineStepResourceCollection(
        "--compute-info",
        output_resources={"info_ht": get_info(split=False, test=test)},
    )
    split_info_ann = PipelineStepResourceCollection(
        "--split-info",
        output_resources={"split_info_ht": get_info(test=test)},
        pipeline_input_steps=[compute_info],
    )
    export_info_vcf = PipelineStepResourceCollection(
        "--export-info-vcf",
        output_resources={"info_vcf": info_vcf_path(test=test)},
        pipeline_input_steps=[compute_info],
    )

    # Add all steps to the variant QC annotation pipeline resource collection.
    ann_pipeline.add_steps(
        {
            "compute_info": compute_info,
            "split_info": split_info_ann,
            "export_info_vcf": export_info_vcf,
        }
    )

    return ann_pipeline


def main(args):
    """Generate all variant annotations needed for variant QC."""
    hl.init(
        default_reference="GRCh38",
        log="/variant_qc_annotations.log",
        tmp_dir="gs://gnomad-tmp-4day",
    )
    test_dataset = args.test_dataset
    test_n_partitions = args.test_n_partitions
    test = test_dataset or test_n_partitions
    overwrite = args.overwrite
    resources = get_variant_qc_annotation_resources(test=test, overwrite=overwrite)
    mt = get_gnomad_v4_vds(
        test=test_dataset,
        high_quality_only=True,
        # Keep control/truth samples because they are used in variant QC.
        high_quality_only_keep_controls=True,
        annotate_meta=True,
    ).variant_data

    if test_n_partitions:
        mt = mt._filter_partitions(range(test_n_partitions))

    if args.compute_info:
        # TODO: is there any reason to also compute info per platform?
        res = resources.compute_info
        res.check_resource_existence()
        if test_dataset:
            unrelated_expr = ~mt.meta.rand_sampling_meta.related
        else:
            unrelated_expr = ~mt.meta.sample_filters.relatedness_filters.related
        default_compute_info(
            mt,
            site_annotations=True,
            ac_filter_groups={"release": mt.meta.release, "unrelated": unrelated_expr},
        ).write(res.info_ht.path, overwrite=overwrite)

    if args.split_info:
        res = resources.split_info
        res.check_resource_existence()
        split_info(res.info_ht.ht()).write(res.split_info_ht.path, overwrite=overwrite)

    if args.export_info_vcf:
        res = resources.export_info_vcf
        res.check_resource_existence()
        hl.export_vcf(adjust_vcf_incompatible_types(res.info_ht.ht()), res.info_vcf)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--slack-channel", help="Slack channel to post results and notifications to."
    )
    parser.add_argument("--overwrite", help="Overwrite data", action="store_true")
    parser.add_argument(
        "--test-dataset", help="Use the test dataset as input.", action="store_true"
    )
    parser.add_argument(
        "--test-n-partitions",
        help="Use only 2 partitions of the VDS as input for testing purposes.",
        nargs="?",
        const=2,
        type=int,
    )
    parser.add_argument("--compute-info", help="Compute info HT.", action="store_true")
    parser.add_argument("--split-info", help="Split info HT.", action="store_true")
    parser.add_argument(
        "--export-info-vcf", help="Export info as VCF.", action="store_true"
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
