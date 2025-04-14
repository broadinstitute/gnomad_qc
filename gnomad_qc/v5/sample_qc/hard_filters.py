"""Script to determine samples that fail hard filtering thresholds."""

import logging

import hail as hl
from gnomad.resources.grch38.reference_data import telomeres_and_centromeres
from gnomad.sample_qc.filtering import compute_stratified_sample_qc
from gnomad.utils.annotations import bi_allelic_expr

from gnomad_qc.v5.resources.basics import get_aou_vds

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("hard_filters")
logger.setLevel(logging.INFO)


def compute_aou_sample_qc(
    test: bool = False,
) -> hl.Table:
    """
    Perform sample QC on the VDS.

    :return: Table containing sample QC metrics
    :rtype: hl.Table
    """
    logger.info("Computing sample QC")

    # Filter the first 100 samples for testing or filter the first 10 partitions?
    vds = get_aou_vds(
        autosomes_only=True, split=True, n_partitions=10 if test else None
    )

    # Remove centromeres and telomeres
    vds = hl.vds.filter_intervals(
        vds, intervals=telomeres_and_centromeres.ht(), keep=False
    )

    sample_qc_ht = compute_stratified_sample_qc(
        vds,
        strata={
            "bi_allelic": bi_allelic_expr(vds.variant_data),
            "multi_allelic": ~bi_allelic_expr(vds.variant_data),
        },
        tmp_ht_prefix=get_sample_qc().path[:-3],
    )

    # Remove annotations that cannot be computed from the sparse format
    sample_qc_ht = sample_qc_ht.annotate(
        **{
            x: sample_qc_ht[x].drop(
                "n_called", "n_not_called", "n_filtered", "call_rate"
            )
            for x in sample_qc_ht.row_value
        }
    )

    return sample_qc_ht.repartition(100)
