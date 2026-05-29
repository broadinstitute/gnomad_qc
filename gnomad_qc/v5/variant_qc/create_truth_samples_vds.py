"""
Script to create a VDS of the 8 Genomes-in-a-Bottle (GiaB) truth samples from their gVCFs.

The GiaB gVCFs were sequenced with the same protocol as the AoU v8 data and live in
``truth_samples_aou_gvcfs_bucket``. They are already reblocked, so they are passed
straight into Hail's VDS combiner (no reblocking step).

The combiner needs the per-gVCF paths up front. The truth-sample bucket cannot be
listed and the sample IDs might be sensitive, so neither the paths nor the IDs are stored
in this repo. Instead the script reads a single-column TSV manifest of gVCF paths
(``truth_samples_gvcf_paths``) by known object path.

This is intended to run in the ``batch`` environment (Hail Batch in the AoU authorization
domain), since that is where the AoU truth-sample gVCFs are readable.
"""

import argparse
import logging
from typing import List

import hail as hl

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.resources.basics import _init_hail, get_logging_path, qc_temp_prefix
from gnomad_qc.v5.resources.variant_qc import (
    get_truth_samples_combiner_plan,
    truth_samples_gvcf_paths,
    truth_samples_vds,
)

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("create_truth_samples_vds")
logger.setLevel(logging.INFO)

# Small interval (PCSK9 gene interval) used by --test.
TEST_INTERVAL = "chr1:55039447-55064852"


def read_gvcf_paths(manifest_path: str) -> List[str]:
    """
    Read the list of truth-sample gVCF paths from a single-column TSV manifest.

    Lines that do not start with ``gs://`` (e.g. a header line, blank lines, or
    comments) are ignored, so the manifest can optionally have a header.

    :param manifest_path: GCS path to the single-column TSV of gVCF paths.
    :return: List of gVCF GCS paths.
    """
    with hl.hadoop_open(manifest_path) as f:
        gvcf_paths = [line.strip() for line in f if line.strip().startswith("gs://")]

    if not gvcf_paths:
        raise ValueError(
            f"No gVCF paths (lines starting with 'gs://') found in {manifest_path}."
        )

    logger.info("Read %d gVCF path(s) from %s.", len(gvcf_paths), manifest_path)
    return gvcf_paths


def _verify_test_vds_against_gvcf(vds: hl.vds.VariantDataset, gvcf_path: str) -> None:
    """
    Cross-check the test VDS's variant calls against the source gVCF over TEST_INTERVAL.

    Imports the single source gVCF directly, restricts to :data:`TEST_INTERVAL`, and
    compares the per-call ``(locus, alleles, genotype)`` set against the combined VDS.
    This confirms the combiner preserved the actual variant data, not just that the VDS
    is structurally well-formed.

    gVCF variant records carry a trailing ``<NON_REF>`` symbolic allele (reference blocks
    have only ``<NON_REF>``); it is dropped so alleles line up with the VDS
    ``variant_data`` representation. The VDS stores local ``LGT``/``LA``, converted to a
    global ``GT`` with :func:`hl.vds.lgt_to_gt` for the comparison.

    :param vds: The combined VDS to check.
    :param gvcf_path: Path to the single source gVCF that was combined.
    :return: None.
    """
    interval = hl.parse_locus_interval(TEST_INTERVAL, reference_genome="GRCh38")
    non_ref = "<NON_REF>"

    # Source gVCF: variant calls over the interval.
    gvcf = hl.import_vcf(
        gvcf_path,
        force_bgz=True,
        reference_genome="GRCh38",
        array_elements_required=False,
    )
    gvcf = hl.filter_intervals(gvcf, [interval])
    gvcf = gvcf.filter_rows(hl.any(gvcf.alleles[1:].map(lambda a: a != non_ref)))
    gvcf = gvcf.annotate_rows(_alleles=gvcf.alleles.filter(lambda a: a != non_ref))
    gvcf_calls = {
        (str(e.locus), tuple(e._alleles)): str(e.GT) for e in gvcf.entries().collect()
    }

    # Combined VDS: variant calls over the interval (LGT/LA -> global GT).
    vd = hl.filter_intervals(vds.variant_data, [interval])
    vd = vd.annotate_entries(_gt=hl.vds.lgt_to_gt(vd.LGT, vd.LA))
    vds_calls = {
        (str(e.locus), tuple(e.alleles)): str(e._gt) for e in vd.entries().collect()
    }

    if vds_calls != gvcf_calls:
        only_gvcf = sorted(set(gvcf_calls) - set(vds_calls))[:5]
        only_vds = sorted(set(vds_calls) - set(gvcf_calls))[:5]
        mismatched = [
            (k, gvcf_calls[k], vds_calls[k])
            for k in set(gvcf_calls) & set(vds_calls)
            if gvcf_calls[k] != vds_calls[k]
        ][:5]
        raise ValueError(
            f"Combined VDS does not match the source gVCF over {TEST_INTERVAL}. "
            f"Sites only in gVCF: {only_gvcf}; sites only in VDS: {only_vds}; "
            f"genotype mismatches (site, gvcf, vds): {mismatched}."
        )

    logger.info(
        "Verified %d variant call(s) in the VDS match the source gVCF over %s.",
        len(vds_calls),
        TEST_INTERVAL,
    )


def validate_vds(vds_path: str, test: bool, manifest_path: str) -> None:
    """
    Verify a combined truth-samples VDS looks correct.

    Runs Hail's :meth:`VariantDataset.validate` (the canonical structural check) plus
    cheap sanity counts: non-empty data, expected sample count, and (in test mode) that
    all variant loci fall within :data:`TEST_INTERVAL`. In test mode it additionally
    cross-checks every variant call against the source gVCF via
    :func:`_verify_test_vds_against_gvcf`.

    :param vds_path: Path to the VDS to validate.
    :param test: Whether the VDS was produced by a ``--test`` run (1 gVCF, one interval).
    :param manifest_path: Manifest path, used to determine the expected sample count for
        a full (non-test) run.
    :return: None.
    """
    vds = hl.vds.read_vds(vds_path)
    vds.validate()

    n_samples = vds.n_samples()
    n_var = vds.variant_data.count_rows()
    n_ref = vds.reference_data.count_rows()
    sample_ids = vds.variant_data.s.collect()
    logger.info(
        "VDS %s: %d sample(s) %s, %d variant rows, %d reference rows.",
        vds_path,
        n_samples,
        sample_ids,
        n_var,
        n_ref,
    )

    if n_var == 0 and n_ref == 0:
        raise ValueError(f"VDS {vds_path} has no variant or reference data.")

    if test:
        if n_samples != 1:
            raise ValueError(f"Expected 1 sample in test VDS, found {n_samples}.")
        interval = hl.parse_locus_interval(TEST_INTERVAL, reference_genome="GRCh38")
        n_outside = vds.variant_data.filter_rows(
            ~interval.contains(vds.variant_data.locus)
        ).count_rows()
        if n_outside:
            raise ValueError(
                f"{n_outside} variant rows fall outside {TEST_INTERVAL}; the interval "
                "restriction did not work."
            )
        # Cross-check the combined data against the single source gVCF.
        _verify_test_vds_against_gvcf(vds, read_gvcf_paths(manifest_path)[0])
    else:
        n_expected = len(read_gvcf_paths(manifest_path))
        if n_samples != n_expected:
            raise ValueError(
                f"Expected {n_expected} samples (one per manifest gVCF), found "
                f"{n_samples}."
            )

    logger.info("VDS %s passed validation.", vds_path)


def main(args):
    """Create a VDS of the GiaB truth samples from their gVCFs."""
    # NOTE: Script assumes we will only run this in Batch.
    environment = "batch"
    overwrite = args.overwrite
    _init_hail("create_truth_samples_vds", environment=environment)

    manifest_path = args.manifest_path
    output_path = (
        f"{qc_temp_prefix(environment=environment)}truth_samples.vds"
        if args.test
        else truth_samples_vds.path
    )

    try:
        if args.create_truth_samples_vds:
            check_resource_existence(
                input_step_resources={"manifest": [manifest_path]},
                output_step_resources={"--create-truth-samples-vds": [output_path]},
                overwrite=overwrite,
            )

            gvcf_paths = read_gvcf_paths(manifest_path)

            if args.test:
                # Cheapest possible real run: combine a single gVCF over one small
                # interval, written to a temp path.
                gvcf_paths = gvcf_paths[:1]
                interval_kwargs = {
                    "intervals": [
                        hl.parse_locus_interval(
                            TEST_INTERVAL, reference_genome="GRCh38"
                        )
                    ]
                }
                logger.info("Test mode: combining 1 gVCF over %s.", TEST_INTERVAL)
            else:
                # Genome-optimized partitioning for the full combine.
                interval_kwargs = {"use_genome_default_intervals": True}

            logger.info(
                "Running the VDS combiner on %d gVCF(s) -> %s",
                len(gvcf_paths),
                output_path,
            )
            combiner = hl.vds.new_combiner(
                output_path=output_path,
                temp_path=qc_temp_prefix(environment=environment, days=4),
                # Save the combiner plan so a failed/interrupted run can resume.
                save_path=get_truth_samples_combiner_plan(test=args.test),
                gvcf_paths=gvcf_paths,
                reference_genome="GRCh38",
                # Overwrite any stale saved combiner plan when --overwrite is set.
                force=overwrite,
                **interval_kwargs,
            )
            combiner.run()

        if args.validate_truth_samples_vds:
            check_resource_existence(
                input_step_resources={"--create-truth-samples-vds": [output_path]},
            )
            validate_vds(output_path, args.test, manifest_path)

    finally:
        logger.info("Copying hail log to logging bucket...")
        hl.copy_log(
            get_logging_path("create_truth_samples_vds", environment=environment)
        )


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--overwrite",
        help="Overwrite the output VDS and any stale saved combiner plan.",
        action="store_true",
    )
    parser.add_argument(
        "--test",
        help=(
            "Cheap end-to-end test: combine only the first gVCF in the manifest over "
            f"a single small interval ({TEST_INTERVAL}) and write to a temporary path."
        ),
        action="store_true",
    )
    parser.add_argument(
        "--manifest-path",
        help=(
            "GCS path to the single-column TSV listing each truth-sample gVCF path "
            "(one 'gs://...' path per line)."
        ),
        default=truth_samples_gvcf_paths,
    )
    parser.add_argument(
        "--create-truth-samples-vds",
        help="Run the VDS combiner to create the truth-samples VDS.",
        action="store_true",
    )
    parser.add_argument(
        "--validate-truth-samples-vds",
        help=(
            "Validate the combined VDS (structural invariants + sample/row sanity "
            "counts). Can be run in the same invocation as --create-truth-samples-vds."
        ),
        action="store_true",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
