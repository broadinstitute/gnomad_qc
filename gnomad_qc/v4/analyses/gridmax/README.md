# GridMax analyses

## plot_bin_summary.py

Plots per-bin sample composition from the `bin_summary.tsv` produced by `compute_gridmax_freq.py --export-bin-summary`.

**Requirements:** `matplotlib`, `pandas` (gnomad_qc conda env). Reads TSVs from GCS via `gsutil cat` or from a local path.

**Usage:**
```bash
python plot_bin_summary.py \
    --summary gs://gnomad-tmp-4day/julia/gridmax/ub5.bin_summary.tsv \
    --output-dir /path/to/output \
    [--top-n 50]
```

**Outputs** (all PNG, written to `--output-dir`):

| File | Description |
|------|-------------|
| `bin_by_pop.png` | Stacked bar — top N bins by size, colored by gnomAD ancestry |
| `bin_by_data_type.png` | Stacked bar — top N bins by size, colored by exome/genome |
| `pop_composition.png` | Pie — overall ancestry composition of the passing training set |
| `data_type_composition.png` | Pie — overall exome vs. genome split of the passing training set |
| `bin_size_hist.png` | Histogram of passing bin sizes (log y-scale) |

Bins labeled `suppressed` in the TSV (below `--min-bin-size`) are excluded from all plots except the bin size histogram, which omits them from the histogram but prints the suppressed count to stdout.
