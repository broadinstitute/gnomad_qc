import hail as hl
from v3.resources import get_full_mt


last_END_position_path = 'gs://gnomad/annotations/hail-0.2/ht/genomes_v3/gnomad_genomes_v3_last_END_positions.ht'

# END RESOURCES

mt = get_full_mt(False)
mt = mt.select_entries('END')
t = mt._localize_entries('__entries', '__cols')
t = t.select(
    last_END_position=hl.or_else(
        hl.min(
            hl.scan.array_agg(
                lambda entry: hl.scan._prev_nonnull(
                    hl.or_missing(
                        hl.is_defined(entry.END),
                        hl.tuple([
                            t.locus,
                            entry.END
                        ])
                    )
                ),
                t.__entries
            ).map(
                lambda x: hl.or_missing(
                    (x[1] >= t.locus.position) & (x[0].contig == t.locus.contig),
                    x[0].position
                )
            )
        ),
        t.locus.position
    )
)
t.write(last_END_position_path, overwrite=True)

