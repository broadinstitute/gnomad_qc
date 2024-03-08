import logging

import hail as hl

from gnomad_qc.v4.resources.annotations import get_freq
from gnomad_qc.v4.resources.basics import get_checkpoint_path

logging.basicConfig(format="%(levelname)s (%(name)s %(lineno)s): %(message)s")
logger = logging.getLogger("temp_genome_freq_fix")
logger.setLevel(logging.INFO)

logger.info(
    "Setting allele frequency to missing for any call stats entries with AN == 0..."
)
ht = get_freq(version="4.0", data_type="genomes").ht()
ht = ht.checkpoint(get_checkpoint_path("temp_genome_4.0_freq_fix"))

ht = ht.annotate(
    freq=ht.freq.map(lambda x: x.annotate(AF=hl.or_missing(x.AN > 0, x.AF)))
)
ht.write(get_freq(version="4.0", data_type="genomes").path, overwrite=True)
