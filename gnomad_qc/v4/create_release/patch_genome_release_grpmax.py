"""Script to create release sites HT for v4.0 exomes and genomes."""
import logging

import hail as hl
from gnomad.resources.grch38.gnomad import POPS_TO_REMOVE_FOR_POPMAX
from gnomad.utils.annotations import pop_max_expr

from gnomad_qc.v4.resources.release import release_sites

logging.basicConfig(
    format="%(asctime)s (%(name)s %(lineno)s): %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)
logger = logging.getLogger("patch_grpmax_release_ht")
logger.setLevel(logging.INFO)

ht = release_sites(data_type="genomes").ht()
ht = ht.checkpoint("gs://gnomad-tmp/patch_grpmax_genomes_release.ht")

# Compute grpmax.
grpmax = pop_max_expr(
    ht.freq, ht.freq_meta, POPS_TO_REMOVE_FOR_POPMAX, pop_label="gen_anc"
)
grpmax = grpmax.annotate(
    faf95=ht.faf[
        ht.faf_meta.index(lambda y: y.values() == ["adj", grpmax.gen_anc])
    ].faf95,
)

logger.info("Annotating 'grpmax'...")
ht = ht.annotate(grpmax=grpmax)
ht.write(release_sites(data_type="genomes").path, overwrite=True)
