# GridMax analyses

## plot_bin_summary.py

Plots per-bin sample composition from the `bin_summary.tsv` produced by `compute_gridmax_freq.py --export-bin-summary`.

**Requirements:** `matplotlib`, `pandas` (gnomad_qc conda env). Reads TSVs from GCS via `gsutil cat` or from a local path.

**Usage:**
```bash
python plot_bin_summary.py \
    --summary gs://gnomad-tmp-4day/julia/gridmax/rb5.bin_summary.tsv \
    --output-dir /path/to/output \
    [--top-n 50] [--min-bin-occupancy 50]
```

**Outputs** (all PNG, written to `--output-dir`):

| File | Description |
|------|-------------|
| `bin_by_pop.png` | Stacked bar — top N bins by size, colored by gnomAD ancestry |
| `bin_by_data_type.png` | Stacked bar — top N bins by size, colored by exome/genome |
| `pop_composition.png` | Horizontal bar — overall ancestry composition of the passing training set |
| `data_type_composition.png` | Horizontal bar — overall exome vs. genome split of the passing training set |
| `bin_occupancy_hist.png` | Histogram of passing bin occupancies (log y-scale) |

Per-level hierarchy plots (`x{factor}_by_pop.png`, `x{factor}_by_data_type.png`, `x{factor}_size_hist.png`) are also written for each coarsening factor present in the summary.

All composition plots and the histogram use **passing bins only** — those whose training-fold sample count is at least `--min-bin-occupancy` (default 50). Bins below the threshold are dropped, not relabeled. The script prints the total and passing training-sample counts and the number of passing groups at each hierarchy level to stdout.
