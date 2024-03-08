"""Temporary script to set XX samples chrY callstats to missing in final freq HT."""

import hail as hl
from gnomad.utils.annotations import set_female_y_metrics_to_na_expr

from gnomad_qc.v4.resources.annotations import get_freq

hl.init(
    log="/fix_chrY_missing_in_freq.log",
    default_reference="GRCh38",
    tmp_dir="gs://gnomad-tmp-4day",
)

ht = get_freq().ht()
ht = ht.annotate(freq=set_female_y_metrics_to_na_expr(ht))

ht = ht.checkpoint(hl.utils.new_temp_file("freq_temp", "ht"))
ht.write(get_freq().path, overwrite=True)
