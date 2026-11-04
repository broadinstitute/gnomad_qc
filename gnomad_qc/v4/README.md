## Dependency provenance
The v4 scripts are a record of what produced the published release and are not
maintained against current `gnomad_methods` `main`. The published v4.1 exome tables
were produced with `gnomad_methods` before commit `0e38c42e` and Hail 0.2.122.

`gnomad_methods` 0.9.1 ([#849](https://github.com/broadinstitute/gnomad_methods/pull/849))
changed `get_het_ab_adj_expr` to read the non-ref allele depth for phased alt-first
het calls (`1|0`). v4.1 exome adj AC, AF, AN and adj histograms were computed with the
previous reading, which took the reference depth for those calls. Re-running
`annotations/generate_freq.py`, `annotations/compute_coverage.py` or
`annotations/fix_freq_an.py` against 0.9.1 or later will not reproduce the published
numbers. To reproduce them, pin the pre-fix ref:

```
git+https://github.com/broadinstitute/gnomad_methods@0e38c42e
```

## Contents and ordering of scripts
The scripts below are run approximately in the order they are listed.


### Resources
**resources/** - Resources for all gnomAD input and output files at each step of the QC and release pipeline.
