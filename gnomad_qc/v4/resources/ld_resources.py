"""Script containing annotation related resources."""

from typing import Optional

from gnomad.resources.grch38.gnomad import CURRENT_EXOME_RELEASE, CURRENT_GENOME_RELEASE
from gnomad.resources.resource_utils import (
    GnomadPublicBlockMatrixResource,
    GnomadPublicTableResource,
)

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
    return f'gs://gnomad-tmp-30day/ld/matrix/gnomad.{data_type}.r{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.bm'

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
    return f'gs://gnomad-tmp-30day/ld/index/gnomad.{data_type}.r{version}.{pop}.{"common." if common_only else ""}{"adj." if adj else ""}ld.variant_indices.ht'
# 
def ld_scores_path(
    data_type: str, pop: str, adj: bool = True, version: Optional[str] = None
):
    if version is None:
        version = (
            CURRENT_EXOME_RELEASE if data_type == "exomes" else CURRENT_GENOME_RELEASE
        )
    return f'gs://gnomad-tmp-30day/ld/scores/gnomad.{data_type}.r{version}.{pop}.{"adj." if adj else ""}ld_scores.ht'


def ld_pruned_path(data_type: str, pop: str, r2: str, version: str = CURRENT_RELEASE):
    return f"gs://gnomad-tmp-30day/ld/pruned/gnomad.{data_type}.r{version}.{pop}.ld.pruned_set.r2_{r2}.ht"

def ld_matrix(pop: str) -> GnomadPublicBlockMatrixResource:
    """Get resource for the LD matrix for the given population."""
    return GnomadPublicBlockMatrixResource(path=ld_matrix_path("genomes", pop))


def ld_index(pop: str) -> GnomadPublicTableResource:
    """Get resource for the LD indices for the given population."""
    return GnomadPublicTableResource(path=ld_index_path("genomes", pop))


def ld_scores(pop: str) -> GnomadPublicTableResource:
    """Get resource for the LD scores for the given population."""
    return GnomadPublicTableResource(path=ld_scores_path("genomes", pop))