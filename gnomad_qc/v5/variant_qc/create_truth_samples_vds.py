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
import hailtop.fs as hfs

from gnomad_qc.resource_utils import check_resource_existence
from gnomad_qc.v5.resources.basics import (
    _get_batch_resource_kwargs,
    _init_hail,
    qc_temp_prefix,
)
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

# Small interval (PCSK9 gene) used by --test to keep the combiner run quick and cheap.
_TEST_CONTIG, _TEST_START, _TEST_END = "chr1", 55039447, 55064852
TEST_INTERVAL = f"{_TEST_CONTIG}:{_TEST_START}-{_TEST_END}"


def _test_interval() -> hl.utils.Interval:
    """
    Build the --test value Interval.

    Returns a value ``Interval`` (not a parsed ``IntervalExpression``) so it can be
    passed straight to ``new_combiner``, which reads ``interval.point_type`` at
    construction. Built lazily here rather than at module load because constructing
    ``hl.Locus`` resolves the reference genome, which would initialize Hail on import
    and break the later ``_init_hail`` call.

    :return: Value Interval spanning :data:`TEST_INTERVAL`.
    """
    return hl.Interval(
        hl.Locus(_TEST_CONTIG, _TEST_START, reference_genome="GRCh38"),
        hl.Locus(_TEST_CONTIG, _TEST_END, reference_genome="GRCh38"),
        includes_end=True,
    )


def read_gvcf_paths(manifest_path: str) -> List[str]:
    """
    Read the list of truth-sample gVCF paths from a single-column TSV manifest.

    Lines that do not start with ``gs://`` (e.g. a header line, blank lines, or
    comments) are ignored, so the manifest can optionally have a header.

    :param manifest_path: GCS path to the single-column TSV of gVCF paths.
    :return: List of gVCF GCS paths.
    """
    with hfs.open(manifest_path) as f:
        gvcf_paths = [line.strip() for line in f if line.strip().startswith("gs://")]

    if not gvcf_paths:
        raise ValueError(
            f"No gVCF paths (lines starting with 'gs://') found in {manifest_path}."
        )

    logger.info("Read %d gVCF path(s) from %s.", len(gvcf_paths), manifest_path)
    return gvcf_paths


def _verify_test_vds_against_gvcf(vds: hl.vds.VariantDataset, gvcf_path: str) -> None:
    """
    Cross-check the test VDS against the source gVCF over TEST_INTERVAL.

    Imports the single source gVCF directly, restricts to :data:`TEST_INTERVAL`, and
    compares both halves of the VDS against it, confirming the combiner preserved the
    actual data (not just that the VDS is structurally well-formed):

    - **Variant calls** (``variant_data``): the per-call ``(locus, alleles, genotype)``
      set. gVCF variant records carry a trailing ``<NON_REF>`` symbolic allele (reference
      blocks have only ``<NON_REF>``); it is dropped so alleles line up with the VDS
      representation. The VDS stores local ``LGT``/``LA``, converted to a global ``GT``
      with :func:`hl.vds.lgt_to_gt` for the comparison.
    - **Reference blocks** (``reference_data``): the per-block ``start_locus -> (END,
      GQ)`` map. Only blocks lying strictly inside the interval are compared: the
      combiner clips reference blocks at the import-interval boundaries, so the first/last
      block would legitimately differ from the gVCF's; interior blocks are untouched and
      must match exactly. Assumes the default combiner reference fields (``GQ`` retained).

    :param vds: The combined VDS to check.
    :param gvcf_path: Path to the single source gVCF that was combined.
    :return: None.
    """
    non_ref = "<NON_REF>"
    test_interval = _test_interval()
    start_pos = test_interval.start.position
    end_pos = test_interval.end.position

    # Source gVCF over the interval.
    gvcf = hl.import_vcf(
        gvcf_path,
        force_bgz=True,
        reference_genome="GRCh38",
        array_elements_required=False,
    )
    gvcf = hl.filter_intervals(gvcf, [test_interval])

    # --- Variant calls ---
    # gVCF variant records have at least one real (non-<NON_REF>) ALT allele.
    var = gvcf.filter_rows(hl.any(gvcf.alleles[1:].map(lambda a: a != non_ref)))
    var = var.annotate_rows(_alleles=var.alleles.filter(lambda a: a != non_ref))
    gvcf_calls = {
        (str(e.locus), tuple(e._alleles)): str(e.GT) for e in var.entries().collect()
    }
    vd = hl.filter_intervals(vds.variant_data, [test_interval])
    vd = vd.annotate_entries(_gt=hl.vds.lgt_to_gt(vd.LGT, vd.LA))
    vds_calls = {
        (str(e.locus), tuple(e.alleles)): str(e._gt) for e in vd.entries().collect()
    }

    # --- Reference blocks (interior only; combiner clips at interval boundaries) ---
    # gVCF reference blocks have <NON_REF> as their only ALT and an INFO/END field.
    ref = gvcf.filter_rows(hl.all(gvcf.alleles[1:].map(lambda a: a == non_ref)))
    gvcf_ref = {
        str(e.locus): (e.info.END, e.GQ)
        for e in ref.entries().collect()
        if start_pos < e.locus.position and e.info.END < end_pos
    }
    rd = hl.filter_intervals(vds.reference_data, [test_interval])
    vds_ref = {
        str(e.locus): (e.END, e.GQ)
        for e in rd.entries().collect()
        if start_pos < e.locus.position and e.END < end_pos
    }

    def _diff(expected: dict, actual: dict) -> str:
        only_expected = sorted(set(expected) - set(actual))[:5]
        only_actual = sorted(set(actual) - set(expected))[:5]
        mismatched = [
            (k, expected[k], actual[k])
            for k in set(expected) & set(actual)
            if expected[k] != actual[k]
        ][:5]
        return (
            f"only in gVCF: {only_expected}; only in VDS: {only_actual}; "
            f"mismatches (key, gvcf, vds): {mismatched}"
        )

    errors = []
    if vds_calls != gvcf_calls:
        errors.append(f"variant calls differ — {_diff(gvcf_calls, vds_calls)}")
    if vds_ref != gvcf_ref:
        errors.append(f"reference blocks differ — {_diff(gvcf_ref, vds_ref)}")
    if errors:
        raise ValueError(
            f"Combined VDS does not match the source gVCF over {TEST_INTERVAL}: "
            + "; ".join(errors)
        )

    logger.info(
        "Verified %d variant call(s) and %d interior reference block(s) in the VDS "
        "match the source gVCF over %s.",
        len(vds_calls),
        len(vds_ref),
        TEST_INTERVAL,
    )


def validate_vds(vds_path: str, test: bool, manifest_path: str) -> None:
    """
    Verify a combined truth-samples VDS looks correct.

    Runs Hail's :meth:`VariantDataset.validate` (the canonical structural check) plus
    cheap sanity counts: non-empty data, expected sample count, and (in test mode) that
    all variant loci fall within :data:`TEST_INTERVAL`. In test mode it additionally
    cross-checks the variant calls and reference blocks against the source gVCF via
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
        n_outside = hl.filter_intervals(
            vds.variant_data, [_test_interval()], keep=False
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
    test = args.test
    _init_hail(
        "create_truth_samples_vds",
        environment=environment,
        billing_project=args.gcp_billing_project,
        experimental=args.experimental,
        **_get_batch_resource_kwargs(args),
    )

    manifest_path = args.manifest_path
    output_path = (
        f"{qc_temp_prefix(environment=environment)}truth_samples.vds"
        if test
        else truth_samples_vds.path
    )

    if args.create_truth_samples_vds:
        check_resource_existence(
            input_step_resources={"manifest": [manifest_path]},
            output_step_resources={"--create-truth-samples-vds": [output_path]},
            overwrite=overwrite,
        )

        gvcf_paths = read_gvcf_paths(manifest_path)

        if test:
            # Cheapest possible real run: combine a single gVCF over one small
            # interval, written to a temp path.
            gvcf_paths = gvcf_paths[:1]
            interval_kwargs = {"intervals": [_test_interval()]}
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
            save_path=get_truth_samples_combiner_plan(test=test),
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
        validate_vds(output_path, test, manifest_path)


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

    # Batch/QoB backend configuration. This script always runs in the batch
    # environment, so these tune the Hail Batch the combiner runs on.
    batch_group = parser.add_argument_group(
        "batch configuration",
        "Optional parameters for the batch/QoB backend.",
    )
    batch_group.add_argument(
        "--gcp-billing-project",
        type=str,
        default="broad-mpg-gnomad",
        help="Google Cloud billing project for reading requester pays buckets.",
    )
    batch_group.add_argument(
        "--experimental",
        action="store_true",
        help=(
            "Route the QoB init through `hl.experimental.init` instead of"
            " `hl.init` and attach the QoB driver to an existing Hail Batch."
            " Requires HAIL_BATCH_ID to be set in the env (Hail Batch"
            " injects this inside batch jobs); raises if not. Without this"
            " flag, each invocation creates its own Hail Batch."
        ),
    )
    batch_group.add_argument(
        "--app-name",
        type=str,
        default=None,
        help="Job name for batch/QoB backend.",
    )
    batch_group.add_argument(
        "--driver-cores",
        type=str,
        default=None,
        help=(
            "Number of cores for the driver node. Pass a power of two"
            " between 0.25 and 16 (as a string, e.g. '2' or '0.5')."
        ),
    )
    batch_group.add_argument(
        "--driver-memory",
        type=str,
        default=None,
        help="Memory type for driver node (e.g., 'highmem').",
    )
    batch_group.add_argument(
        "--jvm-heap-size",
        type=str,
        default=None,
        help=(
            "Max JVM heap size (-Xmx) for the in-process QoB driver under"
            " --experimental, e.g. '5g' or '2500m'. Plumbed to"
            " hl.experimental.init(jvm_heap_size=...). Set to ~50-70%% of"
            " container memory; the rest is for native off-heap"
            " (RegionPool), Python, and OS overhead. Ignored without"
            " --experimental."
        ),
    )
    batch_group.add_argument(
        "--worker-cores",
        type=str,
        default=None,
        help=(
            "Number of cores per worker node. Pass a power of two"
            " between 0.25 and 16 (as a string, e.g. '1' or '0.5')."
        ),
    )
    batch_group.add_argument(
        "--worker-memory",
        type=str,
        default=None,
        help="Memory type for worker nodes (e.g., 'highmem').",
    )

    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
