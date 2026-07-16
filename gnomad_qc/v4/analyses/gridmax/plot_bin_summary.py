"""
Plot per-bin sample composition from the bin summary TSV produced by compute_gridmax_freq.py.

Usage:
    python plot_bin_summary.py \
        --summary gs://gnomad-tmp-4day/julia/gridmax/ub5.bin_summary.tsv \
        --output-dir /path/to/output \
        [--top-n 50]

Requires: matplotlib, pandas (gnomad_qc conda env).
"""

import argparse
import io
import os
import re
import subprocess

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

matplotlib.use("Agg")

GNOMAD_POP_COLORS = {
    "afr": "#941494",
    "ami": "#FFC0CB",
    "amr": "#ED1E24",
    "asj": "#FF7F00",
    "eas": "#108C44",
    "fin": "#002F6C",
    "mid": "#33CC33",
    "nfe": "#6AA5CD",
    "sas": "#FF9912",
    "oth": "#ABB9B9",
    "remaining": "#ABB9B9",
    "suppressed": "#DDDDDD",
}

DATA_TYPE_COLORS = {
    "exome": "#4C72B0",
    "genome": "#DD8452",
}


def read_tsv(path: str) -> pd.DataFrame:
    """Read a TSV from a local path or GCS URI into a DataFrame."""
    if path.startswith("gs://"):
        result = subprocess.run(
            ["gsutil", "cat", path], capture_output=True, check=True
        )
        return pd.read_csv(io.BytesIO(result.stdout), sep="\t")
    return pd.read_csv(path, sep="\t")


def plot_stacked_bar(
    pivot: pd.DataFrame,
    colors: dict,
    title: str,
    ylabel: str,
    out_path: str,
    top_n: int = 50,
) -> None:
    """Plot a stacked bar chart with bins on x-axis, sorted by total size descending."""
    pivot = pivot.copy()
    pivot["_total"] = pivot.sum(axis=1)
    pivot = pivot.sort_values("_total", ascending=False).drop(columns="_total")
    if top_n:
        pivot = pivot.head(top_n)

    fig, ax = plt.subplots(figsize=(max(12, len(pivot) * 0.25), 6))
    bottom = pd.Series(0, index=pivot.index)
    for col in pivot.columns:
        color = colors.get(col, "#888888")
        ax.bar(
            range(len(pivot)),
            pivot[col],
            bottom=bottom,
            color=color,
            label=col,
            width=0.8,
        )
        bottom = bottom + pivot[col]

    ax.set_title(title, fontsize=13)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Bin (sorted by total size)")
    ax.set_xticks([])
    ax.legend(
        title=pivot.columns.name or "",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        fontsize=9,
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_totals_bar(series: pd.Series, colors: dict, title: str, out_path: str) -> None:
    """Plot a horizontal bar chart of sample totals sorted largest-first."""
    series = series.sort_values(ascending=True)  # ascending so largest bar is at top
    total = series.sum()
    fig, ax = plt.subplots(figsize=(9, max(3, len(series) * 0.55)))
    palette = [colors.get(k, "#888888") for k in series.index]
    bars = ax.barh(range(len(series)), series.values, color=palette, edgecolor="white")
    ax.set_yticks(range(len(series)))
    ax.set_yticklabels(series.index)
    ax.set_xlabel("Training samples")
    ax.set_title(title, fontsize=13)
    for val, bar in zip(series.values, bars):
        pct = 100 * val / total if total > 0 else 0
        ax.text(
            val + total * 0.005,
            bar.get_y() + bar.get_height() / 2,
            f"{val:,} ({pct:.1f}%)",
            va="center",
            fontsize=8,
        )
    ax.set_xlim(0, total * 1.3)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_bin_size_hist(total_per_bin: pd.Series, title: str, out_path: str) -> None:
    """Plot a log-scale histogram of bin sizes and save to out_path."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(total_per_bin, bins=40, color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Bin size (total training samples)")
    ax.set_ylabel("Number of bins")
    ax.set_title(title, fontsize=13)
    ax.set_yscale("log")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def _coarse_col(factor: int) -> str:
    return f"bin_id_x{factor}"


def _hierarchy_factors(df: pd.DataFrame) -> list:
    """Return coarsening factors present in df columns, sorted numerically."""
    factors = []
    for col in df.columns:
        m = re.fullmatch(r"bin_id_x(\d+)", col)
        if m:
            factors.append(int(m.group(1)))
    return sorted(factors)


def _plot_level(df: pd.DataFrame, col: str, min_bin_size: int, top_n: int, output_dir: str, label: str) -> int:
    """
    Plot ancestry and data_type composition for one hierarchy level.

    A group passes when its EXOME training-fold count >= min_bin_size (matching the
    pipeline); composition plots still show all data types within passing groups.
    Returns the number of passing groups.
    """
    exome = df[df["data_type"] == "exome"] if "data_type" in df.columns else df
    totals = exome.groupby(col)["n"].sum()
    passing_ids = totals[totals >= min_bin_size].index
    sub = df[df[col].isin(passing_ids)]
    if sub.empty:
        return 0

    pop_pivot = sub.groupby([col, "pop"])["n"].sum().unstack(fill_value=0)
    pop_pivot.columns.name = "pop"
    plot_stacked_bar(
        pop_pivot,
        GNOMAD_POP_COLORS,
        f"{label}: ancestry (top {top_n} groups by size, n≥{min_bin_size})",
        "Training samples",
        os.path.join(output_dir, f"{label}_by_pop.png"),
        top_n=top_n,
    )

    dt_pivot = sub.groupby([col, "data_type"])["n"].sum().unstack(fill_value=0)
    dt_pivot.columns.name = "data_type"
    plot_stacked_bar(
        dt_pivot,
        DATA_TYPE_COLORS,
        f"{label}: data type (top {top_n} groups by size, n≥{min_bin_size})",
        "Training samples",
        os.path.join(output_dir, f"{label}_by_data_type.png"),
        top_n=top_n,
    )

    plot_bin_size_hist(
        totals,
        f"{label}: group size distribution",
        os.path.join(output_dir, f"{label}_size_hist.png"),
    )

    return len(passing_ids)


def main(args):
    """Generate all bin summary plots and print composition summary."""
    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Reading {args.summary}...")
    df = read_tsv(args.summary)

    # Restrict to training-fold samples (eval_fold==0), matching the pipeline:
    # only training-fold samples contribute to the AF.
    if "eval_fold" in df.columns:
        df = df[df["eval_fold"] == 0]
    df = df[df["bin_id"].notna()]

    # Passing bins match the pipeline: a bin passes when its EXOME training-fold
    # count >= min_bin_size. Genomes are in the joint sample key but do not feed
    # the (exome-only) AF or the pipeline's suppression, so they are excluded from
    # the passing/size determination. Composition plots below still show all data
    # types within the passing bins.
    exome = df[df["data_type"] == "exome"] if "data_type" in df.columns else df
    bin_totals = exome.groupby("bin_id")["n"].sum()
    passing_bins = bin_totals[bin_totals >= args.min_bin_size].index
    passing = df[df["bin_id"].isin(passing_bins)]

    # --- 1. Stacked bar: pop composition per bin (passing only) ---
    pop_pivot = passing.groupby(["bin_id", "pop"])["n"].sum().unstack(fill_value=0)
    pop_pivot.columns.name = "pop"
    plot_stacked_bar(
        pop_pivot,
        GNOMAD_POP_COLORS,
        f"L0 (exact): ancestry (top {args.top_n} bins by size, n≥{args.min_bin_size})",
        "Training samples",
        os.path.join(args.output_dir, "bin_by_pop.png"),
        top_n=args.top_n,
    )

    # --- 2. Stacked bar: data_type composition per bin (passing only) ---
    dt_pivot = passing.groupby(["bin_id", "data_type"])["n"].sum().unstack(fill_value=0)
    dt_pivot.columns.name = "data_type"
    plot_stacked_bar(
        dt_pivot,
        DATA_TYPE_COLORS,
        f"L0 (exact): data type (top {args.top_n} bins by size, n≥{args.min_bin_size})",
        "Training samples",
        os.path.join(args.output_dir, "bin_by_data_type.png"),
        top_n=args.top_n,
    )

    # --- 3. Bar: overall pop composition of passing bins ---
    plot_totals_bar(
        passing.groupby("pop")["n"].sum(),
        GNOMAD_POP_COLORS,
        "Training samples by ancestry (passing bins)",
        os.path.join(args.output_dir, "pop_composition.png"),
    )

    # --- 4. Bar: overall data_type composition of passing bins ---
    plot_totals_bar(
        passing.groupby("data_type")["n"].sum(),
        DATA_TYPE_COLORS,
        "Training samples by data type (passing bins)",
        os.path.join(args.output_dir, "data_type_composition.png"),
    )

    # --- 5. Histogram of passing bin sizes ---
    plot_bin_size_hist(
        bin_totals[bin_totals >= args.min_bin_size],
        "L0 (exact): passing bin size distribution",
        os.path.join(args.output_dir, "bin_size_hist.png"),
    )

    # --- 6. Hierarchy levels ---
    hierarchy_counts = {"L0 (exact)": len(passing_bins)}
    for factor in _hierarchy_factors(df):
        col = _coarse_col(factor)
        label = f"x{factor}"
        n_passing = _plot_level(df, col, args.min_bin_size, args.top_n, args.output_dir, label)
        hierarchy_counts[f"/{factor}"] = n_passing

    # --- Print quick summary ---
    print("\nComposition of passing fine bins (training fold, pop x data_type):")
    summary = passing.groupby(["pop", "data_type"])["n"].sum().unstack(fill_value=0)
    print(summary.to_string())

    print(
        f"\nTraining samples: {int(df['n'].sum()):,} total, "
        f"{int(passing['n'].sum()):,} in passing bins "
        f"({len(passing_bins):,} bins ≥ {args.min_bin_size})."
    )

    print("\nGroups passing privacy filter by hierarchy level:")
    for level, count in hierarchy_counts.items():
        print(f"  {level}: {count:,}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--summary",
        required=True,
        help="GCS or local path to bin_summary.tsv produced by compute_gridmax_freq.py.",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Local directory for PNG output files.",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=50,
        help="Show only the top N bins by total size in bar charts (default: 50).",
    )
    parser.add_argument(
        "--min-bin-size",
        type=int,
        default=50,
        help="Minimum sample count for a group to appear in hierarchy plots (default: 50).",
    )
    args = parser.parse_args()
    main(args)
