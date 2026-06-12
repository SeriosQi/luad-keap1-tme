#!/usr/bin/env python3
"""
Task 1 Supplement — Predatory Matrix Correlation Heatmap (Python)
=================================================================
Reproduces correlation heatmap from saved R outputs OR from a CSV expression matrix.

Usage:
  python scripts/correlation_heatmap.py \\
    --expr results/tcga_predatory_expr.csv \\
    --group results/predatory_module_scores.csv \\
    --out results/figures/heatmap_python

  # Or use pre-computed correlation matrix from R:
  python scripts/correlation_heatmap.py \\
    --cor-matrix results/correlation_r_KEAP1-MUT.csv \\
    --out results/figures/heatmap_python_KEAP1_MUT
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats


PREDATORY_GENES = ["SLC7A11", "GGT1", "SLC1A5", "ABCC1", "ABCC2", "ABCC3"]


def spearman_correlation_matrix(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Compute Spearman r and p-value matrices."""
    genes = [g for g in PREDATORY_GENES if g in df.columns]
    sub = df[genes]
    n = len(genes)
    r_mat = pd.DataFrame(np.eye(n), index=genes, columns=genes)
    p_mat = pd.DataFrame(np.zeros((n, n)), index=genes, columns=genes)

    for i, g1 in enumerate(genes):
        for j, g2 in enumerate(genes):
            if i < j:
                rho, p = stats.spearmanr(sub[g1], sub[g2])
                r_mat.loc[g1, g2] = r_mat.loc[g2, g1] = rho
                p_mat.loc[g1, g2] = p_mat.loc[g2, g1] = p
    return r_mat, p_mat


def annotate_matrix(r_mat: pd.DataFrame, p_mat: pd.DataFrame, alpha: float = 0.05) -> np.ndarray:
    """Build cell annotation strings: correlation + significance star."""
    annot = r_mat.copy().astype(str)
    for i in r_mat.index:
        for j in r_mat.columns:
            val = r_mat.loc[i, j]
            star = "*" if (i != j and p_mat.loc[i, j] < alpha) else ""
            annot.loc[i, j] = f"{val:.2f}{star}"
    return annot.values


def plot_heatmap(
    r_mat: pd.DataFrame,
    p_mat: pd.DataFrame | None,
    title: str,
    out_prefix: Path,
    dpi: int = 300,
) -> None:
    """Publication-quality clustered correlation heatmap."""
    sns.set_theme(style="white", context="talk", font_scale=0.85)

    fig, ax = plt.subplots(figsize=(8, 7))
    mask = None
    annot = annotate_matrix(r_mat, p_mat) if p_mat is not None else True

    sns.heatmap(
        r_mat.astype(float),
        ax=ax,
        cmap="RdBu_r",
        center=0,
        vmin=-1,
        vmax=1,
        annot=annot if p_mat is not None else True,
        fmt="" if p_mat is not None else ".2f",
        linewidths=0.5,
        linecolor="white",
        square=True,
        cbar_kws={"label": "Spearman rho", "shrink": 0.8},
        annot_kws={"size": 11},
    )
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    plt.xticks(rotation=45, ha="right", style="italic")
    plt.yticks(rotation=0, style="italic")
    plt.tight_layout()

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    fig.savefig(f"{out_prefix}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_prefix}.pdf / .png")


def load_expression(expr_path: Path, group_path: Path | None, group_filter: str | None):
    expr = pd.read_csv(expr_path, index_col=0)
    if group_path and group_filter:
        meta = pd.read_csv(group_path)
        status_col = "keap1_status" if "keap1_status" in meta.columns else meta.columns[-1]
        sample_col = "sample_id" if "sample_id" in meta.columns else meta.columns[0]
        keep = meta.loc[meta[status_col] == group_filter, sample_col].astype(str)
        cols = [c for c in expr.columns if c in set(keep)]
        expr = expr[cols]
    # genes as rows
    if expr.shape[0] > expr.shape[1]:
        pass
    else:
        expr = expr.T
    return np.log2(expr + 1)


def main():
    parser = argparse.ArgumentParser(description="Predatory matrix correlation heatmap")
    parser.add_argument("--expr", type=Path, help="Expression matrix CSV (genes x samples)")
    parser.add_argument("--group", type=Path, help="Sample metadata CSV with keap1_status")
    parser.add_argument("--group-filter", type=str, default="KEAP1-MUT",
                        help="Filter samples by group label")
    parser.add_argument("--cor-matrix", type=Path, help="Pre-computed correlation matrix CSV")
    parser.add_argument("--title", type=str, default="Predatory Matrix Co-expression")
    parser.add_argument("--out", type=Path, required=True, help="Output prefix (no extension)")
    args = parser.parse_args()

    if args.cor_matrix:
        r_mat = pd.read_csv(args.cor_matrix, index_col=0)
        plot_heatmap(r_mat, p_mat=None, title=args.title, out_prefix=args.out)
        return

    if not args.expr:
        parser.error("Provide --expr or --cor-matrix")

    expr = load_expression(args.expr, args.group, args.group_filter)
    genes = [g for g in PREDATORY_GENES if g in expr.index]
    sub = expr.loc[genes].T
    r_mat, p_mat = spearman_correlation_matrix(sub)
    title = f"{args.title}\n({args.group_filter}, n={sub.shape[0]} samples)"
    plot_heatmap(r_mat, p_mat, title=title, out_prefix=args.out)


if __name__ == "__main__":
    main()
