import hail as hl

sites_ht = hl.read_table(
    "gs://gnomad/release/4.0/ht/exomes/gnomad.exomes.v4.0.sites.filtered_combined_faf.ht"
)
validated_sites_ht = hl.read_table(
    "gs://gnomad/release/4.0/ht/exomes/gnomad.exomes.v4.0.validated_release.filtered_combined_faf.ht"
)

sites_ht.write(
    "gs://gnomad/release/4.0/ht/exomes/gnomad.exomes.v4.0.sites.ht", overwrite=True
)
validated_sites_ht.write(
    "gs://gnomad/release/4.0/ht/exomes/gnomad.exomes.v4.0.validated_release.ht",
    overwrite=True,
)

# Manually transfer the datasets json file:
# gsutil cp
# gs://gnomad/release/4.0/json/exomes/gnomad.exomes.v4.0.included_datasets.filtered_combined_faf.json
# gs://gnomad/release/4.0/json/exomes/gnomad.exomes.v4.0.included_datasets.json
