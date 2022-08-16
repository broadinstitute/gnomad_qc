#!/usr/bin/env python3

import argparse
import logging
import json
import textwrap

import hail as hl

from gnomad.utils.file_utils import file_exists
from gnomad.utils.slack import slack_notifications

from gnomad_qc.v4.resources.basics import get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    cuking_input_path,
    cuking_output_path,
    get_predetermined_qc,
    pca_related_samples_to_drop,
    pca_samples_rankings,
    relatedness,
)
from gnomad.resources.resource_utils import DataException
from gnomad.sample_qc.relatedness import compute_related_samples_to_drop
from gnomad_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def main(args):
    test = args.test
overwrite = args.overwrite

    if args.print_cuking_command:
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

    hl.init(
        log="/relatedness.log",
        default_reference="GRCh38",
    )

    try:
        if args.prepare_inputs:
            input_path = cuking_input_path(test=test)

            if file_exists(input_path) and not overwrite:
                raise DataException(f"{input_path} already exists and the --overwrite option was not used!")

            mt_v3 = hl.read_matrix_table(
                get_predetermined_qc(version="3.1", test=test).path
            )
            mt_v4 = hl.read_matrix_table(get_predetermined_qc(test=test).path)

            # Check that there are no duplicate sample names.
            sample_list = mt_v3.s.collect() + mt_v4.s.collect()
            if len(sample_list) != len(set(sample_list)):
                raise DataException("duplicate samples found in v3 + v4")

            # Check that the number of sites match.
            row_count_v3 = mt_v3.count_rows()
            row_count_v4 = mt_v4.count_rows()
            if row_count_v3 != row_count_v4:
                raise DataException("number of sites doesn't match between v3 and v4")

            # Combine v3 and v4 dense QC matrices, as we need to compute relatedness
            # across the union of samples.
            mt = mt_v3.union_cols(mt_v4)

            # Compute the field required for KING.
            mt = mt.select_entries(n_alt_alleles=mt.GT.n_alt_alleles())

            # Create a Table based on entries alone. By adding a row and column index,
            # we can drop the locus, alleles, and sample field, as well as skip writing
            # missing entries.
            mt = mt.select_globals()
            mt = mt.add_row_index()
            mt = mt.add_col_index()
            entries = mt.entries()
            entries = entries.key_by()
            entries = entries.select(
                entries.row_idx, entries.col_idx, entries.n_alt_alleles
            )

            # Export one compressed Parquet file per partition.
            logger.info("Writing Parquet tables...")
            entry_write = entries.to_spark().write.option("compression", "zstd")
            if args.overwrite:
                entry_write = entry_write.mode('overwrite')
            entry_write.parquet(input_path)


            # Write metadata that's useful for postprocessing. Map `col_idx` to `s`
            # explicitly so we don't need to rely on a particular order returned by
            # `collect`.
            logger.info("Writing metadata...")
            col_idx_mapping = hl.struct(col_idx=mt.col_idx, s=mt.s).collect()
            col_idx_mapping.sort(key=lambda e: e.col_idx)
            metadata = {
                "num_sites": row_count_v4,
                "samples": [e.s for e in col_idx_mapping],
            }
            with hl.hadoop_open(f"{input_path}/metadata.json", 'w') as f:
                json.dump(metadata, f)

        if args.create_relatedness_table:
            from hail.utils.java import Env

            spark = Env.spark_session()
            df = spark.read.parquet(cuking_output_path(test=test))
            ht = hl.Table.from_spark(df)
            ht = ht.repartition(64)
            ht = ht.key_by(ht.i, ht.j)
            ht.write(relatedness.path)

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
        "--slack-channel",
        help="Slack channel to post results and notifications to.",
    )

    args = parser.parse_args()

    if args.slack_channel:
        with slack_notifications(slack_token, args.slack_channel):
            main(args)
    else:
        main(args)
