import logging
import os
from abc import ABC
from enum import Enum
from functools import reduce, wraps
from typing import Iterable, Union

import hail as hl
from hail.linalg import BlockMatrix
from gnomad.resources.resource_utils import (
    BaseResource,
    BlockMatrixResource,
    DataException,
    MatrixTableResource,
    PedigreeResource,
    ResourceNotAvailable,
    TableResource,
    VariantDatasetResource,
)

logger = logging.getLogger("gnomad_qc.resource_utils")


GNOMAD_PRIVATE_BUCKETS = ("gnomad", "gnomad-tmp")
"""
Private buckets for gnomAD data.
"""


class GnomadPrivateResourceSource(Enum):
    """Sources for private gnomAD resources."""

    GNOMAD_PRODUCTION = "gnomAD production private bucket"
    GNOMAD_PRODUCTION_TEMP = "gnomAD production temporary private bucket"


def get_default_private_resource_source() -> Union[GnomadPrivateResourceSource, str]:
    """
    Get the default source for private gnomAD resources.

    The default source is determined by...

    - If the ``GNOMAD_DEFAULT_PRIVATE_RESOURCE_SOURCE`` environment variable is set, use the source configured there.
    - Otherwise, use gnomAD production private bucket.

    :returns: Default resource source
    """
    default_source_from_env = os.getenv("GNOMAD_DEFAULT_PRIVATE_RESOURCE_SOURCE", None)
    if default_source_from_env:
        # Convert to a GnomadPrivateResourceSource enum if possible
        try:
            default_source = GnomadPrivateResourceSource(default_source_from_env)
            logger.info(
                "Using configured source for gnomAD private resources: %s",
                default_source.value,
            )
            return default_source
        except ValueError:
            logger.info(
                "Using configured custom source for gnomAD private resources: %s",
                default_source_from_env,
            )
            return default_source_from_env

    return GnomadPrivateResourceSource.GNOMAD_PRODUCTION


class _GnomadPrivateResourceConfiguration:
    """Configuration for private gnomAD resources."""

    _source: Union[GnomadPrivateResourceSource, str, None] = None

    @property
    def source(self) -> Union[GnomadPrivateResourceSource, str]:
        """
        Get the source for private gnomAD resource files.

        This is used to determine which URLs gnomAD resources will be loaded from.
        :returns: Source name or path to root of resources directory
        """
        if self._source is None:
            self._source = get_default_private_resource_source()

        return self._source

    @source.setter
    def source(self, source: Union[GnomadPrivateResourceSource, str]) -> None:
        """
        Set the default source for resource files.

        This is used to determine which URLs gnomAD resources will be loaded from.
        :param source: Source name or path to root of resources directory
        """
        self._source = source


gnomad_private_resource_configuration = _GnomadPrivateResourceConfiguration()


def check_resource_existence(
    resource: Union[str, BaseResource],
    raise_exist_no_overwrite_error: bool = False,
):
    if not isinstance(resource, str):
        resource = resource.path

    paths_to_test = [resource]
    if any(resource.endswith(ext) for ext in (".ht", ".mt", ".bm", ".parquet")):
        paths_to_test = [f"{resource}/_SUCCESS"]

    if resource.endswith(".vds"):
        paths_to_test = [
            f"{resource}/reference_data/_SUCCESS",
            f"{resource}/variant_data/_SUCCESS",
        ]

    exists = all(
        [hl.current_backend().fs.exists(path_to_test) for path_to_test in paths_to_test]
    )

    if exists and raise_exist_no_overwrite_error:
        raise DataException(
            f"{resource} already exists and the --overwrite option was not used!"
        )

    return exists


def set_gnomad_test():
    gnomad_private_resource_configuration.source = (
        GnomadPrivateResourceSource.GNOMAD_PRODUCTION_TEMP
    )


# TODO check for existence before running if overwrite is False?
class GnomadPrivateResource(BaseResource, ABC):
    """Base class for the gnomAD project's private resources."""

    def _get_path(self, test: bool = False) -> str:
        if test:
            resource_source = GnomadPrivateResourceSource.GNOMAD_PRODUCTION_TEMP
        else:
            resource_source = gnomad_private_resource_configuration.source
        if resource_source == GnomadPrivateResourceSource.GNOMAD_PRODUCTION:
            return self._path

        relative_path = reduce(
            lambda path, bucket: path[5 + len(bucket) :]
            if path.startswith(f"gs://{bucket}/")
            else path,
            GNOMAD_PRIVATE_BUCKETS,
            self._path,
        )

        if resource_source == GnomadPrivateResourceSource.GNOMAD_PRODUCTION_TEMP:
            return f"gs://gnomad-tmp{relative_path}"

        return (
            f"{resource_source.rstrip('/')}{relative_path}"  # pylint: disable=no-member
        )

    def _set_path(self, path):
        if not any(
            path.startswith(f"gs://{bucket}/") for bucket in GNOMAD_PRIVATE_BUCKETS
        ):
            raise ValueError(
                f"GnomadPrivateResource requires a path to a file in one of the private gnomAD buckets ({', '.join(GNOMAD_PRIVATE_BUCKETS)})"
            )

        return super()._set_path(path)

    test_path = property(
        fget=lambda self: self._get_path(test=True),
    )


class GnomadPrivateTableResource(TableResource, GnomadPrivateResource):
    """Resource class for a private Hail Table generated by gnomAD production."""

    def ht(self, force_import: bool = False, test: bool = False) -> hl.Table:
        """
        Read and return the Hail Table resource.

        :return: Hail Table resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_args)
        elif test:
            return hl.read_table(self.test_path)
        else:
            return hl.read_table(self.path)


class GnomadPrivateMatrixTableResource(MatrixTableResource, GnomadPrivateResource):
    """Resource class for a private Hail MatrixTable generated by gnomAD production."""

    def mt(self, force_import: bool = False, test: bool = False) -> hl.MatrixTable:
        """
        Read and return the Hail MatrixTable resource.

        :return: Hail MatrixTable resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_args)
        elif test:
            return hl.read_matrix_table(self.test_path)
        else:
            return hl.read_matrix_table(self.path)


class GnomadPrivatePedigreeResource(PedigreeResource, GnomadPrivateResource):
    """Resource class for a private pedigree generated by gnomAD production."""

    def ht(self, test: bool = False) -> hl.Table:
        """
        Read the pedigree into a family HT using hl.import_fam().

        :return: Family table
        """
        if test:
            path = self.test_path
        else:
            path = self.path

        return hl.import_fam(
            path,
            quant_pheno=self.quant_pheno,
            delimiter=self.delimiter,
            missing=self.missing,
        )

    def pedigree(self, test: bool = False) -> hl.Pedigree:
        """
        Read the pedigree into an hl.Pedigree using hl.Pedigree.read().

        :return: pedigree
        """
        if test:
            path = self.test_path
        else:
            path = self.path

        return hl.Pedigree.read(path, delimiter=self.delimiter)


class GnomadPrivateBlockMatrixResource(BlockMatrixResource, GnomadPrivateResource):
    """Resource class for a private Hail BlockMatrix generated by gnomAD production."""

    def bm(self, test: bool = False) -> BlockMatrix:
        """
        Read and return the Hail MatrixTable resource.

        :return: Hail MatrixTable resource
        """
        if test:
            path = self.test_path
        else:
            path = self.path

        return BlockMatrix.read(path)


class GnomadPrivateVariantDatasetResource(
    VariantDatasetResource, GnomadPrivateResource
):
    """Resource class for a private Hail VariantDataset generated by gnomAD production."""

    def vds(
        self, force_import: bool = False, test: bool = False
    ) -> hl.vds.VariantDataset:
        """
        Read and return the Hail VariantDataset resource.

        :return: Hail VariantDataset resource
        """
        if self.path is None or force_import:
            return self.import_func(**self.import_args)
        elif test:
            return hl.vds.read_vds(self.test_path)
        else:
            return hl.vds.read_vds(self.path)
