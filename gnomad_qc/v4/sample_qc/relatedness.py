#!/usr/bin/env python3

import argparse
import logging
import json
import textwrap

import hail as hl

from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.slack_creds import slack_token
from gnomad_qc.v4.resources.basics import get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    cuking_input_path,
    cuking_output_path,
    pca_related_samples_to_drop,
    pca_samples_rankings,
    qc,
    relatedness,
)

from cuKING.mt_to_cuking_inputs import mt_to_cuking_inputs
from cuKING.cuking_outputs_to_ht import cuking_outputs_to_ht

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def main(args):
    test = args.test
    overwrite = args.overwrite

    if args.print_cuking_command:
        if (
            args.prepare_inputs
            or args.create_relatedness_table
            or args.compute_related_samples_to_drop
        ):
            raise ValueError(
                "--print-cuking-command can't be used simultaneously with other run modes"
            )

        print(
            textwrap.dedent(
                f"""\
                # Make sure that $PROJECT_ID is set.
                cd cuKING && \\
                ./cloud_batch_submit.py \\
                    --location=us-central1 \\
                    --project_id=$PROJECT_ID \\
                    --tag_name=$(git describe --tags) \\
                    --input_uri={cuking_input_path(test=test)} \\
                    --output_uri={cuking_output_path(test=test)} \\
                    --requester_pays_project=$PROJECT_ID \\
                    --kin_threshold={args.second_degree_kin_cutoff} \\
                    --split_factor=4 && \\
                cd ..
                """
            )
        )
        return

    hl.init(log="/relatedness.log", default_reference="GRCh38")

    try:
        if args.prepare_inputs:
            parquet_uri = cuking_input_path(test=test)

            if file_exists(parquet_uri) and not overwrite:
                raise DataException(
                    f"{parquet_uri} already exists and the --overwrite option was not used!"
                )

            mt_to_cuking_inputs(
                mt=qc.mt(), parquet_uri=parquet_uri, overwrite=overwrite
            )

        if args.create_relatedness_table:
            if file_exists(relatedness.path) and not overwrite:
                raise DataException(
                    f"{relatedness.path} already exists and the --overwrite option was not used!"
                )

            ht = cuking_outputs_to_ht(parquet_uri=cuking_output_path(test=test))
            ht = ht.repartition(args.relatedness_n_partitions)
            ht.write(relatedness.path, overwrite=overwrite)

        if args.compute_related_samples_to_drop:
            # compute_related_samples_to_drop uses a rank Table as a tie breaker when
            # pruning samples. We determine rankings as follows: filtered samples are
            # ranked lowest, followed by non-releasable samples, followed by presence in
            # v3 (favoring v3 over v4), and then by coverage.
            def annotate_project_meta_ht(
                project_meta_ht, sex_ht, hard_filtered_samples_ht
            ):
                return project_meta_ht.select(
                    project_meta_ht.releasable,
                    project_meta_ht.exclude,
                    chr20_mean_dp=sex_ht[project_meta_ht.key].chr20_mean_dp,
                    filtered=hl.or_else(
                        hl.len(
                            hard_filtered_samples_ht[project_meta_ht.key].hard_filters
                        )
                        > 0,
                        False,
                    ),
                )

            from gnomad_qc.v3.resources.meta import project_meta as project_meta_v3
            from gnomad_qc.v4.resources.meta import project_meta as project_meta_v4
            from gnomad_qc.v3.resources.sample_qc import (
                hard_filtered_samples as hard_filtered_samples_v3,
            )
            from gnomad_qc.v4.resources.sample_qc import (
                hard_filtered_samples as hard_filtered_samples_v4,
            )
            from gnomad_qc.v3.resources.sample_qc import sex as sex_v3
            from gnomad_qc.v4.resources.sample_qc import sex as sex_v4

            annotated_project_meta_v3_ht = annotate_project_meta_ht(
                project_meta_v3.ht(), sex_v3.ht(), hard_filtered_samples_v3.ht()
            )
            annotated_project_meta_v4_ht = annotate_project_meta_ht(
                project_meta_v4.ht(), sex_v4.ht, hard_filtered_samples_v4.ht()
            )

            rank_ht = annotated_project_meta_v3_ht.union(
                annotated_project_meta_v4_ht, unify=True
            )

            rank_ht = rank_ht.order_by(
                rank_ht.filtered,
                hl.desc(rank_ht.releasable & ~rank_ht.exclude),
                ~hl.is_defined(  # Favor v3 samples.
                    annotated_project_meta_v3_ht[rank_ht.s]
                ),
                hl.desc(rank_ht.chr20_mean_dp),
            ).add_index(name="rank")

            rank_ht = rank_ht.key_by(rank_ht.s)
            rank_ht = rank_ht.select(rank_ht.filtered, rank_ht.rank)

            rank_ht.write(pca_samples_rankings.path)
            rank_ht = pca_samples_rankings.ht()

            filtered_samples = hl.literal(
                rank_ht.aggregate(
                    hl.agg.filter(rank_ht.filtered, hl.agg.collect_as_set(rank_ht.s))
                )
            )

            samples_to_drop = compute_related_samples_to_drop(
                relatedness.ht(),
                rank_ht,
                args.second_degree_kin_cutoff,
                filtered_samples=filtered_samples,
            )
            samples_to_drop = samples_to_drop.key_by(samples_to_drop.s)
            samples_to_drop.write(pca_related_samples_to_drop.path)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("relatedness"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Overwrite output files.",
        action="store_true",
    )
    parser.add_argument(
        "--test", help="Use a test MatrixTableResource as input.", action="store_true"
    )
    parser.add_argument(
        "--prepare-inputs",
        help="Converts the dense QC MatrixTable to a Parquet format suitable for cuKING.",
        action="store_true",
    )
    parser.add_argument(
        "--print-cuking-command",
        help="Print the command to submit a Cloud Batch job for running cuKING.",
        action="store_true",
    )
    parser.add_argument(
        "--create-relatedness-table",
        help="Convert the cuKING outputs to a standard Hail Table.",
        action="store_true",
    )
    parser.add_argument(
        "--compute-related-samples-to-drop",
        help="Determine the minimal set of related samples to prune.",
        action="store_true",
    )
    parser.add_argument(
        "--second-degree-kin-cutoff",
        help="Minimum kinship threshold for filtering a pair of samples with a second degree relationship\
        in KING and filtering related individuals. (Default = 0.1) \
        Bycroft et al. (2018) calculates 0.08838835 but from evaluation of the distributions v3 and v4 has used 0.1",
        default=0.1,
        type=float,
    )
    parser.add_argument(
        "--relatedness-n-partitions",
        help="Number of desired partitions for the relatedness Table.",
        default=50,
        type=int,
    )
    parser.add_argument(
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
