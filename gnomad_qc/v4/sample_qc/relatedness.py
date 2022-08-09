#!/usr/bin/env python3

import argparse
import logging
import json
import textwrap

import hail as hl
from hail.utils import hadoop_exists, hadoop_open

from gnomad.utils.slack import slack_notifications

from gnomad_qc.v4.resources.basics import get_logging_path
from gnomad_qc.v4.resources.sample_qc import (
    get_predetermined_qc_sites,
    cuking_input_path,
    cuking_output_path,
)
from gnomad.resources.resource_utils import DataException
from gnomad_qc.slack_creds import slack_token

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("relatedness")
logger.setLevel(logging.INFO)


def main(args):
    test = args.test

    input_path = cuking_input_path(test=test)
    output_path = cuking_output_path(test=test)

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
                    --input_uri={input_path} \\
                    --output_uri={output_path} \\
                    --requester_pays_project=$PROJECT_ID \\
                    --kin_threshold=0.05 \\
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
            if hadoop_exists(input_path):
                raise DataException(f"{input_path} already exists")

            mt_v3 = hl.read_matrix_table(
                get_predetermined_qc_sites(version="3.1", test=test).path
            )
            mt_v4 = hl.read_matrix_table(get_predetermined_qc_sites(test=test).path)

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

            # Create a table based on entries alone. By adding a row and column index,
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
            entries.to_spark().write.option("compression", "zstd").parquet(input_path)

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
            with hadoop_open(f"{input_path}/metadata.json", 'w') as f:
                json.dump(metadata, f)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(get_logging_path("relatedness"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--test", help="Use a test MatrixTableResource as input.", action="store_true"
    )
    parser.add_argument(
        "--prepare-inputs",
        help="Converts the dense QC matrix to a Parquet format suitable for cuKING.",
        action="store_true",
    )
    parser.add_argument(
        "--print-cuking-command",
        help="Print the command to submit a Cloud Batch job for running cuKING.",
        action="store_true",
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
