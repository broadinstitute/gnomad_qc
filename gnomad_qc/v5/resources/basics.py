"""Script containing generic resources."""

from gnomad.resources.resource_utils import VariantDatasetResource

# v5 DRAGEN TGP test VDS.
dragen_tgp_vds = VariantDatasetResource(
    "gs://gnomad/v5.0/testing/genomes/dragen_tgp_v5.0_test.vds"
)
