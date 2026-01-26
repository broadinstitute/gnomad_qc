"""Script containing annotation related resources."""

from typing import Optional

from gnomad.resources.grch38.gnomad import CURRENT_EXOME_RELEASE, CURRENT_GENOME_RELEASE

COMMON_FREQ = 0.005
RARE_FREQ = 0.0005


def ld_matrix_path(
    data_type: str,
    pop: str,
    freq: float = COMMON_FREQ,
    adj: bool = True,
    ld_contig: str = None,
    version: Optional[str] = None,
    test: bool = False,
    custom_suffix: str = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/matrix/gnomad.{data_type}.{version}{f".test" if test else ""}.{"all_contigs" if not ld_contig else f"{ld_contig}"}.{pop}.{f"af{freq}"}.{"adj." if adj else ""}ld{custom_suffix if custom_suffix else ""}.bm'


def ld_index_path(
    data_type: str,
    pop: str,
    freq: float = COMMON_FREQ,
    adj: bool = True,
    ld_contig: str = None,
    version: Optional[str] = None,
    test: bool = False,
    custom_suffix: str = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/index/gnomad.{data_type}.{version}{f".test" if test else ""}.{"all_contigs" if not ld_contig else f"{ld_contig}"}.{pop}.{f"af{freq}"}.{"adj." if adj else ""}ld.variant_indices{custom_suffix if custom_suffix else ""}.ht'


#
def ld_scores_path(
    data_type: str,
    pop: str,
    adj: bool = True,
    ld_contig: str = None,
    freq: float = COMMON_FREQ,
    version: Optional[str] = None,
    test: bool = False,
    custom_suffix: str = None,
    call_rate_cutoff: float = 0.8,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/scores/gnomad.{data_type}.{version}{f".test" if test else ""}.{"all_contigs" if not ld_contig else f"{ld_contig}"}.{pop}.{f"af{freq}"}.{f"callrate{call_rate_cutoff}"}.{"adj." if adj else ""}ld_scores{f"{custom_suffix}" if custom_suffix else ""}.ht'


def ld_mt_checkpoint_path(
    data_type: str,
    freq: float = COMMON_FREQ,
    pop: str = None,
    version: str = CURRENT_GENOME_RELEASE,
    mt_contig: str = None,
    test: bool = False,
    adj: bool = False,
):
    if version is None:
        version = CURRENT_GENOME_RELEASE
    return f'gs://gnomad-tmp-30day/ld/gnomad.{data_type}.{"all_pops" if not pop else f"{pop}"}.{f"af{freq}"}.{"all_contigs" if not mt_contig else f"{mt_contig}"}.{version}{f".adj" if adj else ""}{f".test" if test else ""}.mt'


def ld_pruned_path(
    data_type: str,
    pop: str,
    r2: str,
    freq: float = COMMON_FREQ,
    ld_contig: str = None,
    version: str = CURRENT_GENOME_RELEASE,
    test: bool = False,
    ld_set: bool = False,
    adj: bool = False,
):
    if version is None:
        version = CURRENT_GENOME_RELEASE
    return f'gs://gnomad-tmp-30day/ld/pruned/gnomad.{data_type}.{version}{f".test" if test else ""}.{f"af{freq}"}.{"all_contigs" if not ld_contig else f"{ld_contig}"}.{pop}.ld.{f"pruned_set" if not ld_set else "ld_set"}{f".adj" if adj else ""}.r2_{r2}.ht'
