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
import subprocess
import sys

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


def plot_pie(series: pd.Series, colors: dict, title: str, out_path: str) -> None:
    """Plot a pie chart of the given series and save to out_path."""
    fig, ax = plt.subplots(figsize=(7, 7))
    palette = [colors.get(k, "#888888") for k in series.index]
    wedges, texts, autotexts = ax.pie(
        series,
        labels=series.index,
        colors=palette,
        autopct=lambda p: f"{p:.1f}%" if p > 1.5 else "",
        startangle=90,
        pctdistance=0.8,
    )
    for t in autotexts:
        t.set_fontsize(8)
    ax.set_title(title, fontsize=13)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def plot_bin_size_hist(total_per_bin: pd.Series, title: str, out_path: str) -> None:
    """Plot a log-scale histogram of bin sizes and save to out_path."""
    fig, ax = plt.subplots(figsize=(8, 5))
    passing = total_per_bin[total_per_bin.index != "suppressed"]
    ax.hist(passing, bins=40, color="#4C72B0", edgecolor="white", linewidth=0.5)
    ax.set_xlabel("Bin size (total training samples)")
    ax.set_ylabel("Number of bins")
    ax.set_title(title, fontsize=13)
    ax.set_yscale("log")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_path}")


def main(args):
    """Generate all bin summary plots and print composition summary."""
    import os

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Reading {args.summary}...")
    df = read_tsv(args.summary)

    passing = df[df["bin_id"].notna() & (df["bin_id"] != "suppressed")]
    total_per_bin = df.groupby("bin_id")["n"].sum()

    # --- 1. Stacked bar: pop composition per bin (passing only) ---
    pop_pivot = passing.groupby(["bin_id", "pop"])["n"].sum().unstack(fill_value=0)
    pop_pivot.columns.name = "pop"
    plot_stacked_bar(
        pop_pivot,
        GNOMAD_POP_COLORS,
        f"Training samples per bin by ancestry (top {args.top_n} bins by size)",
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
        f"Training samples per bin by data type (top {args.top_n} bins by size)",
        "Training samples",
        os.path.join(args.output_dir, "bin_by_data_type.png"),
        top_n=args.top_n,
    )

    # --- 3. Pie: overall pop composition of training set ---
    pop_totals = (
        df[df["bin_id"].notna() & (df["bin_id"] != "suppressed")]
        .groupby("pop")["n"]
        .sum()
        .sort_values(ascending=False)
    )
    plot_pie(
        pop_totals,
        GNOMAD_POP_COLORS,
        "Training samples by ancestry (passing bins)",
        os.path.join(args.output_dir, "pop_composition.png"),
    )

    # --- 4. Pie: overall data_type composition ---
    dt_totals = (
        df[df["bin_id"].notna() & (df["bin_id"] != "suppressed")]
        .groupby("data_type")["n"]
        .sum()
    )
    plot_pie(
        dt_totals,
        DATA_TYPE_COLORS,
        "Training samples by data type (passing bins)",
        os.path.join(args.output_dir, "data_type_composition.png"),
    )

    # --- 5. Histogram of bin sizes ---
    plot_bin_size_hist(
        total_per_bin,
        "Distribution of passing bin sizes",
        os.path.join(args.output_dir, "bin_size_hist.png"),
    )

    # --- Print quick summary table ---
    print("\nOverall composition (passing bins):")
    summary = (
        df[df["bin_id"].notna() & (df["bin_id"] != "suppressed")]
        .groupby(["pop", "data_type"])["n"]
        .sum()
        .unstack(fill_value=0)
    )
    print(summary.to_string())

    suppressed = df[df["bin_id"] == "suppressed"]["n"].sum()
    print(f"\nSamples in suppressed bins: {suppressed:,}")
    print(f"Bins passing filter: {len(passing['bin_id'].unique()):,}")


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
    args = parser.parse_args()
    main(args)
