"""
IAD Diagnostic Panel — Multi-Anchor Variant
============================================

Pipeline stage: COMPUTE → SAVE.

Same logic as iad_analysis.py but the anchor is a 3-gene cholangiocyte
identity score (KRT19, KRT7, EPCAM) instead of KRT19 alone.

Why
---
Single-gene anchors are vulnerable to:
  - dropout (a cell with truly high cholangiocyte identity but zero KRT19
    counts due to scRNA-seq sparsity is missed)
  - upstream-of-KRT19 effects (a regulator could drive cholangiocyte
    identity through KRT7 or EPCAM with KRT19 lagging)

A geometric-mean anchor across three canonical biliary markers is more
robust and is standard in cholangiocyte-state literature (e.g.,
Sampaziotis 2017).

Outputs
-------
results/liver_anchor_coexpression.tsv
results/pancreas_anchor_coexpression.tsv
results/iad_diagnostic_panel_multianchor.tsv
figures/comparison/multianchor_top_coexpressed.png
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import (
    CHOLANGIOCYTE_ANCHORS,
    compute_anchor_score,
    ensure_dir,
)
from iad_analysis import (
    HOUSEKEEPING_EXACT, HOUSEKEEPING_PREFIXES, _is_housekeeping,
    PANEL_SIZE,
)

# ── Config ────────────────────────────────────────────────────────────────────
DATA_DIR    = Path("data")
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")

LIVER_H5AD    = DATA_DIR / "liver"    / "liver_cholangiocytes.h5ad"
PANCREAS_H5AD = DATA_DIR / "pancreas" / "pancreas_ductal.h5ad"

TOP_N        = 50
MIN_EXPR_FRAC = 0.01


# ── Anchor-vs-all gene correlation ────────────────────────────────────────────

def correlate_with_anchor(adata: ad.AnnData, label: str) -> pd.Series:
    """Spearman ρ between every gene and the multi-anchor identity score."""
    from scipy.stats import rankdata

    anchor = compute_anchor_score(adata, CHOLANGIOCYTE_ANCHORS)

    X = adata.X
    X_dense = X.toarray() if hasattr(X, "toarray") else np.array(X, dtype=float)
    expr_frac = (X_dense > 0).mean(axis=0)
    keep = expr_frac >= MIN_EXPR_FRAC
    X_filt = X_dense[:, keep]
    gene_names = adata.var_names[keep]

    print(f"  {label}: correlating {X_filt.shape[1]:,} genes with anchor score …")

    a_ranked = rankdata(anchor)
    X_ranked = np.apply_along_axis(rankdata, 0, X_filt)
    a_c = a_ranked - a_ranked.mean()
    X_c = X_ranked - X_ranked.mean(axis=0)

    num = a_c @ X_c
    denom = np.sqrt((a_c ** 2).sum()) * np.sqrt((X_c ** 2).sum(axis=0))
    denom = np.where(denom == 0, 1e-10, denom)
    corr = num / denom

    result = pd.Series(corr, index=gene_names)
    # Drop the anchor genes themselves
    result = result.drop(labels=list(CHOLANGIOCYTE_ANCHORS), errors="ignore")
    return result.sort_values(ascending=False)


# ── Panel build (intersection + filters) ──────────────────────────────────────

def build_panel(
    l_corr: pd.Series, p_corr: pd.Series,
    n: int = PANEL_SIZE,
    min_rho_liver: float = 0.30,
    max_rho_pancreas: float = 0.50,
) -> pd.DataFrame:
    cols = ["gene", "rho_liver", "rho_pancreas", "specificity", "panel_rank"]
    if l_corr.empty or p_corr.empty:
        return pd.DataFrame(columns=cols)

    shared = l_corr.index.intersection(p_corr.index)
    print(f"  Genes scored in both tissues: {len(shared):,}")

    df = pd.DataFrame({
        "rho_liver":    l_corr.reindex(shared),
        "rho_pancreas": p_corr.reindex(shared),
    })
    df["specificity"] = df["rho_liver"] - df["rho_pancreas"]

    keep_thresh = (df["rho_liver"] >= min_rho_liver) & (df["rho_pancreas"] <= max_rho_pancreas)
    keep_clean  = ~df.index.to_series().map(_is_housekeeping)
    df = df[keep_thresh & keep_clean]
    print(f"  After filters: {len(df):,} candidates.")

    df = df.sort_values("specificity", ascending=False)
    df.insert(0, "gene", df.index)
    df.reset_index(drop=True, inplace=True)
    df["panel_rank"] = df.index + 1
    return df.head(n)


# ── Plot ──────────────────────────────────────────────────────────────────────

def plot_top(l: pd.Series, p: pd.Series, out: Path, n: int = 30) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for ax, corr, title in zip(axes, [l, p], ["Liver", "Pancreas"]):
        if corr.empty:
            continue
        top = corr.head(n)
        colors = ["#16a085" if v > 0 else "#c0392b" for v in top.values]
        ax.barh(range(len(top)), top.values[::-1], color=colors[::-1])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top.index[::-1], fontsize=9)
        ax.set_xlabel("Spearman ρ with anchor score")
        ax.set_title(f"{title}: top {n} genes co-expressed with KRT19+KRT7+EPCAM")
        ax.axvline(0, color="black", linewidth=0.7, linestyle="--")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR)

    print("=" * 60)
    print("IAD diagnostic panel — multi-anchor (KRT19+KRT7+EPCAM)")
    print("=" * 60)

    if not LIVER_H5AD.exists() or not PANCREAS_H5AD.exists():
        print("ERROR: liver or pancreas h5ad missing.", file=sys.stderr)
        sys.exit(1)

    print("\n--- Loading ---")
    liver    = sc.read_h5ad(LIVER_H5AD)
    pancreas = sc.read_h5ad(PANCREAS_H5AD)
    print(f"  Liver:    {liver.n_obs:,} × {liver.n_vars:,}")
    print(f"  Pancreas: {pancreas.n_obs:,} × {pancreas.n_vars:,}")

    print("\n--- Anchor-vs-all correlation ---")
    l_corr = correlate_with_anchor(liver,    "Liver")
    p_corr = correlate_with_anchor(pancreas, "Pancreas")

    l_corr.rename("spearman_rho").to_csv(
        RESULTS_DIR / "liver_anchor_coexpression.tsv", sep="\t", header=True,
    )
    p_corr.rename("spearman_rho").to_csv(
        RESULTS_DIR / "pancreas_anchor_coexpression.tsv", sep="\t", header=True,
    )

    print("\n--- Top co-expressed (multi-anchor) ---")
    plot_top(l_corr, p_corr,
             FIGURES_DIR / "comparison" / "multianchor_top_coexpressed.png")

    print("\n--- Building panel ---")
    panel = build_panel(l_corr, p_corr)
    panel.to_csv(
        RESULTS_DIR / "iad_diagnostic_panel_multianchor.tsv",
        sep="\t", index=False,
    )

    print(f"\nMulti-anchor panel ({len(panel)} genes):")
    if not panel.empty:
        for _, row in panel.iterrows():
            print(
                f"  #{int(row['panel_rank']):>2}  {row['gene']:<18}"
                f" ρ_liv={row['rho_liver']:+.3f}"
                f" ρ_pan={row['rho_pancreas']:+.3f}"
                f" spec={row['specificity']:+.3f}"
            )

    print("\nDone. Compare to results/iad_diagnostic_panel.tsv (KRT19-only).")


if __name__ == "__main__":
    main()
