"""Script containing annotation related resources."""

from typing import Optional

from gnomad.resources.grch38.gnomad import CURRENT_EXOME_RELEASE, CURRENT_GENOME_RELEASE


def ld_matrix_path(
    data_type: str,
    pop: str,
    common_only: bool = True,
    adj: bool = True,
    version: Optional[str] = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/matrix/gnomad.{data_type}.{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.bm'


def ld_index_path(
    data_type: str,
    pop: str,
    common_only: bool = True,
    adj: bool = True,
    version: Optional[str] = None,
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/index/gnomad.{data_type}.{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.variant_indices.ht'


#
def ld_scores_path(
    data_type: str, pop: str, adj: bool = True, version: Optional[str] = None
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/scores/gnomad.{data_type}.{version}.{pop}.{"adj." if adj else ""}ld_scores.ht'


def ld_mt_checkpoint_path(
    data_type: str,
    pops: str = None,
    version: str = CURRENT_GENOME_RELEASE,
    test: bool = False,
):
    if version is None:
        version = CURRENT_GENOME_RELEASE
    return f'gs://gnomad-tmp-30day/ld/gnomad.{data_type}.{"all_pops" if not pops else f"{pops}"}.{version}{f".test" if test else ""}.mt'


def ld_pruned_path(
    data_type: str, pop: str, r2: str, version: str = CURRENT_GENOME_RELEASE
):
    if version is None:
        version = CURRENT_GENOME_RELEASE
    return f"gs://gnomad-tmp-30day/ld/pruned/gnomad.{data_type}.{version}.{pop}.ld.pruned_set.r2_{r2}.ht"
