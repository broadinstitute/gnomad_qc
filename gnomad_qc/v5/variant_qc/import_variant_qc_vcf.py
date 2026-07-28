"""Script to load variant QC result VCF into a Hail Table."""

import argparse
import logging
from typing import List, Optional, Tuple, Union

import hail as hl
from gnomad.utils.sparse_mt import split_info_annotation

from gnomad_qc.v5.resources.basics import (
    _check_resource_existence,
    _get_batch_resource_kwargs,
    _init_hail,
)
from gnomad_qc.v5.resources.variant_qc import IF_FEATURES, get_variant_qc_result

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("import_variant_qc_vcf")
logger.setLevel(logging.INFO)


def import_variant_qc_vcf(
    vcf_path: Union[str, List[str]],
    model_id: str,
    num_partitions: int = 5000,
    import_header_path: Optional[str] = None,
    array_elements_required: bool = False,
    is_split: bool = False,
    deduplicate_check: bool = False,
) -> Union[hl.Table, Tuple[hl.Table, hl.Table]]:
    """
    Import variant QC result sites VCF into a HT.

    :param vcf_path: Path(s) to variant QC result sites VCF(s). Can be a Hadoop glob
        pattern or a list of paths.
    :param model_id: Model ID. Must start with ``rf_``, ``vqsr_``, or ``if_``.
    :param num_partitions: Number of partitions for the output HT.
    :param import_header_path: Optional header file to use for import.
    :param array_elements_required: Value passed to hl.import_vcf.
    :param is_split: Whether the VCF is already split.
    :param deduplicate_check: Whether to remove duplicate variants.
    :return: HT (or split, unsplit HTs) containing variant QC results.
    """
    model_type = model_id.split("_")[0]
    if model_type not in ["rf", "vqsr", "if"]:
        raise ValueError(
            f"Model ID must start with 'rf_', 'vqsr_', or 'if_', but got {model_id}"
        )

    logger.info(
        "Importing variant QC annotations for model %s (array_elements_required=%s)...",
        model_id,
        array_elements_required,
    )
    ht = hl.import_vcf(
        vcf_path,
        force_bgz=True,
        reference_genome="GRCh38",
        array_elements_required=array_elements_required,
        header_file=import_header_path,
    ).rows()

    # Read back with `_n_partitions` to repartition off the table index (no shuffle).
    tmp_path = hl.utils.new_temp_file("imported_vcq_result", "ht")
    ht.write(tmp_path)
    ht = hl.read_table(tmp_path, _n_partitions=num_partitions)

    if deduplicate_check:
        dup_ht = ht.group_by(ht.locus, ht.alleles).aggregate(n=hl.agg.count())
        dup_ht = dup_ht.filter(dup_ht.n > 1)
        n_dup = dup_ht.count()
        if n_dup > 0:
            # NOTE: distinct() keeps one arbitrary row per key, so if shards disagree
            # on INFO values the surviving row is non-deterministic.
            logger.warning(
                "Found %s variants with duplicate keys; keeping distinct rows. "
                "Examples: %s",
                n_dup,
                [(e.locus, e.alleles) for e in dup_ht.take(5)],
            )
            ht = ht.distinct()
        else:
            logger.info("No duplicate variants found.")

    unsplit_count = None
    if not is_split:
        if model_type == "vqsr":
            as_vqslod_expr = {"AS_VQSLOD": ht.info.AS_VQSLOD.map(lambda x: hl.float(x))}
        else:
            as_vqslod_expr = {}
        # Reformat the pipe/comma-delimited AS fields when present. They may be absent
        # (e.g. stripped by the isolation forest v4-test reheader, whose v4 encoding
        # breaks import_vcf) and are not required for the variant QC result.
        # TODO: Confirm these are present and pipe-delimited in the v5 info VCF.
        reformat_expr = {}
        if "AS_QUALapprox" in ht.info:
            reformat_expr["AS_QUALapprox"] = ht.info.AS_QUALapprox.split("\|")[1:].map(
                lambda x: hl.int(x)
            )
        if "AS_VarDP" in ht.info:
            reformat_expr["AS_VarDP"] = ht.info.AS_VarDP.split("\|")[1:].map(
                lambda x: hl.int(x)
            )
        if "AS_SB_TABLE" in ht.info:
            reformat_expr["AS_SB_TABLE"] = ht.info.AS_SB_TABLE.split("\|").map(
                lambda x: hl.or_missing(x != "", x.split(",").map(lambda y: hl.int(y)))
            )
        ht = ht.annotate(info=ht.info.annotate(**as_vqslod_expr, **reformat_expr))

        unsplit_ht = ht.checkpoint(hl.utils.new_temp_file("unsplit_vcq_result", "ht"))
        unsplit_count = unsplit_ht.count()
        split_ht = hl.split_multi_hts(unsplit_ht)
        split_ht = split_ht.annotate(
            info=split_ht.info.annotate(
                **split_info_annotation(split_ht.info, split_ht.a_index)
            ),
        )
    else:
        unsplit_ht = None
        split_ht = ht

    split_ht = split_ht.checkpoint(hl.utils.new_temp_file("split_vcq_result", "ht"))
    split_count = split_ht.count()
    if unsplit_count is None:
        logger.info("Found %s split variants in the VCF.", split_count)
    else:
        logger.info(
            "Found %s unsplit and %s split variants in the VCF.",
            unsplit_count,
            split_count,
        )

    if unsplit_ht is None:
        return split_ht
    else:
        return split_ht, unsplit_ht


def main(args):
    """Load a single variant QC result VCF into a Hail Table."""
    environment = args.environment
    _init_hail(
        "import_variant_qc_vcf",
        environment,
        **_get_batch_resource_kwargs(args),
    )

    # An already-split input yields only the split HT; otherwise split and unsplit.
    splits = (True,) if args.is_split else (True, False)
    resources = [
        get_variant_qc_result(
            args.model_id, test=args.test, split=split, environment=environment
        )
        for split in splits
    ]
    # Check up front so a rerun fails before the import instead of after it, and
    # never after the split HT has already been written.
    _check_resource_existence(
        environment=environment,
        output_step_resources={"variant_qc_result": resources},
        overwrite=args.overwrite,
    )

    hts = import_variant_qc_vcf(
        args.vcf_path,
        args.model_id,
        args.n_partitions,
        args.header_path,
        args.array_elements_required,
        args.is_split,
        args.deduplication_check,
    )
    if not isinstance(hts, tuple):
        hts = (hts,)

    for ht, res in zip(hts, resources):
        ht = ht.annotate_globals(
            model_id=args.model_id,
            snp_features=args.snp_features,
            indel_features=args.indel_features,
        )
        ht.write(res.path, overwrite=args.overwrite)


def get_script_argument_parser() -> argparse.ArgumentParser:
    """Get script argument parser."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o",
        "--overwrite",
        help="Whether to overwrite data already present in the output Table.",
        action="store_true",
    )
    parser.add_argument(
        "--test", help="Whether to use a tmp path for output.", action="store_true"
    )
    parser.add_argument(
        "--environment",
        help="Compute environment.",
        default="batch",
        type=str,
        choices=["rwb", "batch", "dataproc"],
    )
    parser.add_argument(
        "--app-name",
        help="Job name for the batch/QoB backend.",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--driver-cores",
        help="Number of driver cores (Batch only).",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--driver-memory",
        help="Driver memory (Batch only).",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--worker-cores",
        help="Number of worker cores (Batch only).",
        default=None,
        type=int,
    )
    parser.add_argument(
        "--worker-memory",
        help="Worker memory (Batch only).",
        default=None,
        type=str,
    )
    parser.add_argument(
        "--vcf-path",
        help="Path to variant QC result VCF. Can be Hadoop glob patterns.",
        required=True,
    )
    parser.add_argument(
        "--model-id", help="Model ID for the result HT.", type=str, required=True
    )
    parser.add_argument(
        "--n-partitions",
        help="Number of partitions for the output Table.",
        default=5000,
        type=int,
    )
    parser.add_argument("--header-path", help="Optional header file to use for import.")
    parser.add_argument(
        "--array-elements-required",
        action="store_true",
        help="Set array_elements_required=True in import_vcf.",
    )
    parser.add_argument(
        "--is-split", action="store_true", help="Whether the VCF is already split."
    )
    parser.add_argument(
        "--deduplication-check",
        action="store_true",
        help="Remove duplicate variants (for overlapping shards).",
    )
    parser.add_argument(
        "--snp-features",
        help="SNP model features (stored as a global).",
        default=IF_FEATURES["snv"],
        type=str,
        nargs="+",
    )
    parser.add_argument(
        "--indel-features",
        help="Indel model features (stored as a global).",
        default=IF_FEATURES["indel"],
        type=str,
        nargs="+",
    )
    return parser


if __name__ == "__main__":
    parser = get_script_argument_parser()
    main(parser.parse_args())
