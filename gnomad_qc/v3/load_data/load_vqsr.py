import logging
from typing import Optional

import hail as hl

from gnomad.utils.sparse_mt import split_info_annotation
from gnomad_qc.v3.resources.annotations import get_filters

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("import_vqsr")
logger.setLevel(logging.INFO)

def import_vqsr(
    vqsr_path: str,
    vqsr_type: str = "AS",
    num_partitions: int = 5000,
    overwrite: bool = False,
    import_header_path: Optional[str] = None,
) -> None:
    """
    Imports vqsr site vcf into a HT
    :param vqsr_path: Path to input vqsr site vcf. This can be specified as Hadoop glob patterns
    :param vqsr_type: One of `AS` (allele specific) or `AS_TS` (allele specific with transmitted singletons)
    :param num_partitions: Number of partitions to use for the VQSR HT
    :param overwrite: Whether to overwrite imported VQSR HT
    :param import_header_path: Optional path to a header file to use for import
    :return: None
    """

    logger.info(f"Importing VQSR annotations for {vqsr_type} VQSR...")
    mt = hl.import_vcf(
        vqsr_path,
        force_bgz=True,
        reference_genome="GRCh38",
        header_file=import_header_path,
    ).repartition(num_partitions)

    ht = mt.rows()

    ht = ht.annotate(
        info=ht.info.annotate(
            AS_VQSLOD=ht.info.AS_VQSLOD.map(lambda x: hl.float(x)),
            AS_QUALapprox=ht.info.AS_QUALapprox.split("\|")[1:].map(
                lambda x: hl.int(x)
            ),
            AS_VarDP=ht.info.AS_VarDP.split("\|")[1:].map(lambda x: hl.int(x)),
            AS_SB_TABLE=ht.info.AS_SB_TABLE.split("\|").map(
                lambda x: x.split(",").map(lambda y: hl.int(y))
            ),
        )
    )

    ht = ht.checkpoint(
        get_filters('vqsr_alleleSpecificTrans', split=False, finalized=False).path,
        overwrite=overwrite,
    )

    unsplit_count = ht.count()
    ht = hl.split_multi_hts(ht)

    ht = ht.annotate(
        info=ht.info.annotate(**split_info_annotation(ht.info, ht.a_index)),
    )

    ht = ht.checkpoint(
        get_filters('vqsr_alleleSpecificTrans', split=True, finalized=False).path,
        overwrite=overwrite,
    )
    split_count = ht.count()
    logger.info(
        f"Found {unsplit_count} unsplit and {split_count} split variants with VQSR annotations"
    )

import_vqsr(
        "gs://gnomad-julia/gnomad_v3.1/VQSR_AS_TS/gnomAD_3.1.filtered.*.vcf.gz",  #TODO: move to gnomad folder and add path to resources
        "AS_TS",
        5000,
        True,
        import_header_path="gs://gnomad-julia/gnomad_v3.1/VQSR_AS/gnomAD3.1.filtered.header",  # TODO: add to resources
    )

# TODO: add agruments
