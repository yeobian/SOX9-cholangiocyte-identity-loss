"""
Step 2 — per-cell PROX1 expression vs PSC pseudotime.

Joins the cholangiocyte AnnData (which now has PROX1) onto the PSC-only
trajectory's per-cell pseudotime. Computes a cell-level Spearman ρ and
saves a scatter with a smoothed trend, colored by disease status.

This is the cell-level equivalent of the donor-level ρ = +0.43 we found in
pseudotime_tf_crossplot.py. Cell-level n is much larger, so the test is
much more powerful.

Outputs
-------
results/prox1_pseudotime_per_cell.tsv
figures/prox1/prox1_vs_pseudotime.png

Usage
-----
    python src/prox1_per_cell_pseudotime.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression

CHOL_PATH        = Path("data") / "liver" / "liver_cholangiocytes.h5ad"
TRAJECTORY_TSV   = Path("results") / "trajectory_psc_per_cell.tsv"
RESULTS_TSV      = Path("results") / "prox1_pseudotime_per_cell.tsv"
FIGURE_PATH      = Path("figures") / "prox1" / "prox1_vs_pseudotime.png"

GENE = "PROX1"


def load_pseudotime() -> pd.DataFrame:
    if not TRAJECTORY_TSV.exists():
        print(f"ERROR: {TRAJECTORY_TSV} missing — run trajectory_psc_only.py",
              file=sys.stderr); sys.exit(1)
    df = pd.read_csv(TRAJECTORY_TSV, sep="\t")
    print(f"  trajectory rows: {len(df):,}")
    df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["pseudotime"])
    print(f"  finite pseudotime rows: {len(df):,}")
    df["cell"] = df["cell"].astype(int)
    return df


def attach_prox1(traj: pd.DataFrame, adata) -> pd.DataFrame:
    """The trajectory 'cell' column holds AnnData obs_names (as ints), not
    positional indices. Look up by name."""
    if GENE not in adata.var_names:
        print(f"ERROR: {GENE} not in AnnData var_names.", file=sys.stderr)
        sys.exit(1)
    prox1_all = get_gene_expression(adata, GENE)

    obs_name_to_pos = {str(name): i for i, name in enumerate(adata.obs_names)}
    traj = traj.copy()
    traj["cell_str"] = traj["cell"].astype(str)
    matched = traj["cell_str"].map(obs_name_to_pos)
    n_missing = int(matched.isna().sum())
    if n_missing == len(traj):
        print(f"ERROR: no cells matched obs_names. Are you using the same "
              f"cholangiocyte AnnData that the trajectory was built on?",
              file=sys.stderr); sys.exit(1)
    if n_missing > 0:
        print(f"  {n_missing:,}/{len(traj):,} trajectory cells not found in "
              f"AnnData — likely dropped at QC; will be excluded.")
        traj = traj[~matched.isna()].copy()
        matched = matched[~matched.isna()]
    pos = matched.astype(int).values
    traj["prox1"] = prox1_all[pos]
    if "dataset_id" in adata.obs.columns:
        traj["dataset_id"] = adata.obs["dataset_id"].astype(str).values[pos]
    return traj


def smoothed(x: np.ndarray, y: np.ndarray, n_bins: int = 30) -> tuple:
    """Bin-mean smoothing — simple, no SciPy LOESS dep."""
    order = np.argsort(x)
    x_s = x[order]; y_s = y[order]
    bins = np.linspace(x_s.min(), x_s.max(), n_bins + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    means   = np.full(n_bins, np.nan)
    for i in range(n_bins):
        mask = (x_s >= bins[i]) & (x_s < bins[i + 1])
        if mask.any():
            means[i] = y_s[mask].mean()
    keep = ~np.isnan(means)
    return centers[keep], means[keep]


def plot(df: pd.DataFrame, out: Path, rho_all: float, p_all: float) -> None:
    fig, ax = plt.subplots(figsize=(10, 6))
    diseases = sorted(df["disease"].unique())
    palette  = dict(zip(diseases, ["#2980b9", "#c0392b", "#27ae60",
                                    "#f39c12", "#8e44ad"][:len(diseases)]))
    # Per-disease: scatter + per-disease smooth line
    for d in diseases:
        sub = df[df["disease"] == d]
        ax.scatter(sub["pseudotime"], sub["prox1"],
                   s=4, c=palette[d], alpha=0.25, label=d)
        if len(sub) > 50:
            xc, ym = smoothed(sub["pseudotime"].values,
                              sub["prox1"].values, n_bins=24)
            ax.plot(xc, ym, color=palette[d], linewidth=2.4)
    # Overall trend (black dashed)
    xc, ym = smoothed(df["pseudotime"].values, df["prox1"].values, n_bins=30)
    ax.plot(xc, ym, color="black", linewidth=2, linestyle="--",
            label=f"All cells trend  (ρ={rho_all:+.3f}, p={p_all:.1e})")
    ax.set_xlabel("PSC pseudotime  (healthy → exhausted)")
    ax.set_ylabel(f"{GENE} expression (log-normalized)")
    ax.set_title(f"Per-cell {GENE} along PSC trajectory  "
                 f"({len(df):,} cholangiocytes)")
    ax.legend(loc="upper left", fontsize=8, frameon=False, markerscale=2)
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_TSV.parent)

    print("=" * 60)
    print("Per-cell PROX1 vs PSC pseudotime")
    print("=" * 60)

    print("\n[1/3] Loading PSC pseudotime …")
    traj = load_pseudotime()

    print(f"\n[2/3] Loading cholangiocytes …")
    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    df = attach_prox1(traj, adata)
    print(f"  joined: {len(df):,} cells with pseudotime + PROX1")

    print("\n[3/3] Cell-level Spearman + plot …")
    rho_all, p_all = spearmanr(df["pseudotime"], df["prox1"])
    print(f"  Cell-level ρ(pseudotime, {GENE}) = {rho_all:+.3f}  p = {p_all:.2e}")

    print("\n  Per-disease cell-level Spearman:")
    rows = []
    for d in sorted(df["disease"].unique()):
        sub = df[df["disease"] == d]
        if len(sub) >= 30:
            rho, p = spearmanr(sub["pseudotime"], sub["prox1"])
        else:
            rho, p = float("nan"), float("nan")
        n_pos = int((sub["prox1"] > 0).sum())
        rows.append({"disease": d, "n_cells": len(sub),
                     "n_prox1_pos": n_pos,
                     "pct_prox1_pos": float(n_pos / len(sub) * 100),
                     "rho_cell_level": rho, "p_cell_level": p})
        print(f"    {d:<35} n={len(sub):,}  "
              f"PROX1+={n_pos/len(sub)*100:5.1f}%  "
              f"ρ={rho:+.3f}  p={p:.2e}")

    pd.DataFrame(rows).to_csv(
        Path("results") / "prox1_pseudotime_per_disease.tsv",
        sep="\t", index=False)

    df.to_csv(RESULTS_TSV, sep="\t", index=False)
    print(f"\n  → {RESULTS_TSV}")

    plot(df, FIGURE_PATH, rho_all, p_all)

    print("\n--- Headline ---")
    print(f"  Cell-level ρ(pseudotime, PROX1) = {rho_all:+.3f}  p = {p_all:.2e}")
    if p_all < 1e-10 and rho_all > 0:
        print("  PROX1 RISES along the PSC trajectory at single-cell resolution.")
    elif rho_all > 0 and p_all < 0.05:
        print("  PROX1 trends UP along the PSC trajectory.")
    else:
        print("  No clear cell-level direction.")


if __name__ == "__main__":
    main()
