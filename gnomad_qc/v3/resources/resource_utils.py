"""Centralized migration warnings for v3 resources moved to new GCS locations."""

import functools
import logging
from typing import Literal


@functools.lru_cache(maxsize=None)
def show_v3_migration_warning(target: Literal["annotations", "sample_qc"]) -> None:
    """
    Show GCS migration warning for v3 resources (only once per target per session).

    :param target: Type of resource that was migrated ("annotations" or "sample_qc")
    """
    logging.warning(
        f"Most objects in 'gs://gnomad/{target}' were moved to either "
        f"'gs://gnomad-archive' or 'gs://gnomad-autoclass'. This function now returns "
        f"'gs://gnomad-autoclass/{target}' paths, if no resource is found there, "
        f"please check 'gs://gnomad-archive/{target}'."
    )
