import hail as hl
from gnomad.resources.resource_utils import (
    TableResource,
    PedigreeResource,
    VersionedPedigreeResource
)

# Samples metadata
META_ROOT = "gs://gnomad/metadata/genomes_v3"
meta = TableResource(f'{META_ROOT}/gnomad_v3_metadata_2019-09-27.ht')
meta_tsv_path = f'{META_ROOT}/gnomad_v3_metadata_2019-09-27.tsv.gz'
project_meta = TableResource(
    import_func=hl.import_table,
    import_args={
        'path': f'{META_ROOT}/09-09-2019_v3_project_meta.txt',
        'impute': True,
        'key': 's',
        'min_partitions': 100
    }
)
pedigree = VersionedPedigreeResource(
    'final',  # TODO: Make sure "final" is the best label once the family scripts are in
    {
        'raw': PedigreeResource(f'{META_ROOT}/gnomad_v3_raw.fam', delimiter="\t"),
        'final': PedigreeResource(f'{META_ROOT}/gnomad_v3.fam', delimiter="\t")
    }
)

trios = VersionedPedigreeResource(  # TODO: Should this be merged with Pedigree into a single resource?
    'final',  # TODO: Make sure "final" is the best label once the family scripts are in
    {
        'raw': PedigreeResource(f'{META_ROOT}/gnomad_v3_trios_raw.fam'),
        'final': PedigreeResource(f'{META_ROOT}/gnomad_v3_trios.fam')
    }
)

ped_mendel_errors = TableResource(f'{META_ROOT}/gnomad_v3_ped_chr20_mendel_errors.ht')
