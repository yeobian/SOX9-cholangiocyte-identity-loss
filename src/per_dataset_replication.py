"""
Per-dataset replication of the headline TF / PROX1 result.

Tests whether the pooled-Census finding (GATA6/SOX9/HNF1B ↓ + PROX1/HNF4A/FOXA1/2 ↑
along disease) is reproducible *within* individual source datasets — i.e.
without any cross-batch pooling that could create the signal artificially.

For each dataset_id with enough cells and donors:
  donor-level Spearman ρ(mean IAD score, mean gene expression) across donors.

The headline genes have an expected direction:
  PROX1, HNF4A, FOXA1, FOXA2          → expected ρ > 0  (rise with disease)
  GATA6, SOX9, HNF1B, ONECUT1         → expected ρ < 0  (drop with disease)

We then ask: in what fraction of datasets does each gene show the predicted
direction? A gene that replicates in ≥ 2 independent datasets in the right
direction is robust to batch effects.

Usage
-----
    python src/per_dataset_replication.py
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
from iad_score import (
    LIVER_H5AD as DEFAULT_CHOL,
    REFERENCE_PERCENTILE,
    _panel_expression_per_cell, _panel_genes_present, load_panel,
)

CHOL_PATH    = DEFAULT_CHOL
RESULTS_TSV  = Path("results") / "per_dataset_replication.tsv"
SUMMARY_TSV  = Path("results") / "per_dataset_replication_summary.tsv"
FIGURE_PATH  = Path("figures") / "comparison" / "per_dataset_replication_heatmap.png"

# Inclusion thresholds
MIN_CELLS_PER_DATASET  = 200
MIN_DONORS_PER_DATASET = 5

# Genes to test, with their expected disease-progression direction
EXPECTED_UP   = ["PROX1", "HNF4A", "FOXA1", "FOXA2"]
EXPECTED_DOWN = ["GATA6", "SOX9", "HNF1B", "ONECUT1"]
ALL_GENES     = EXPECTED_UP + EXPECTED_DOWN


def per_cell_iad_score(adata) -> np.ndarray:
    panel = load_panel()
    genes = _panel_genes_present(adata, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(adata, genes, weights=weights)
    ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    return 1.0 - np.clip(panel_expr / max(ref, 1e-6), 0.0, 1.0)


def aggregate_per_donor(df: pd.DataFrame, gene_cols: list[str]) -> pd.DataFrame:
    agg = {"iad_score": "mean", "n_cells": "size"}
    agg.update({g: "mean" for g in gene_cols})
    df = df.copy()
    df["n_cells"] = 1
    summary = df.groupby("donor").agg(agg)
    return summary


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_TSV.parent)
    ensure_dir(FIGURE_PATH.parent)

    print("=" * 60)
    print("Per-dataset replication — TF / PROX1 along IAD axis")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)

    print(f"\nLoading {CHOL_PATH} …")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if "dataset_id" not in adata.obs.columns:
        print("ERROR: dataset_id missing in obs.", file=sys.stderr); sys.exit(1)
    if "donor_id" not in adata.obs.columns:
        print("ERROR: donor_id missing in obs.", file=sys.stderr); sys.exit(1)

    # Filter the gene list to what's actually present
    present_genes = [g for g in ALL_GENES if g in adata.var_names]
    missing       = [g for g in ALL_GENES if g not in adata.var_names]
    if missing:
        print(f"  Missing genes (will be skipped): {missing}")
    print(f"  Genes available: {present_genes}")

    # Per-cell IAD score + per-cell gene expression
    print("\nComputing per-cell IAD score …")
    iad = per_cell_iad_score(adata)

    print("Building per-cell long table …")
    base = pd.DataFrame({
        "dataset_id": adata.obs["dataset_id"].astype(str).values,
        "donor":      adata.obs["donor_id"].astype(str).values,
        "iad_score":  iad,
    }, index=range(adata.n_obs))
    for g in present_genes:
        base[g] = get_gene_expression(adata, g)

    # Per-dataset filter: enough cells AND enough donors
    counts = base.groupby("dataset_id").agg(
        n_cells=("iad_score", "size"),
        n_donors=("donor", "nunique"),
    )
    keep = counts[(counts["n_cells"] >= MIN_CELLS_PER_DATASET) &
                  (counts["n_donors"] >= MIN_DONORS_PER_DATASET)].index.tolist()
    print(f"\nDatasets passing inclusion (≥{MIN_CELLS_PER_DATASET} cells "
          f"AND ≥{MIN_DONORS_PER_DATASET} donors): {len(keep)}/{len(counts)}")
    if not keep:
        print("ERROR: no datasets pass thresholds. Relax MIN_CELLS / MIN_DONORS.",
              file=sys.stderr); sys.exit(1)

    # Per-dataset Spearman
    rows = []
    for ds in keep:
        sub = base[base["dataset_id"] == ds]
        donor_sum = aggregate_per_donor(sub, present_genes)
        n_donors  = len(donor_sum)
        for g in present_genes:
            x = donor_sum["iad_score"].values
            y = donor_sum[g].values
            ok = ~np.isnan(x) & ~np.isnan(y)
            if ok.sum() < MIN_DONORS_PER_DATASET:
                rho, p = float("nan"), float("nan")
            else:
                rho, p = spearmanr(x[ok], y[ok])
            rows.append({
                "dataset_id":   ds,
                "n_cells":      int(len(sub)),
                "n_donors":     int(n_donors),
                "gene":         g,
                "expected_dir": "+" if g in EXPECTED_UP else "-",
                "rho":          float(rho),
                "p":            float(p),
                "matches_expected": (
                    (g in EXPECTED_UP and rho > 0) or
                    (g in EXPECTED_DOWN and rho < 0)
                    if not np.isnan(rho) else False
                ),
            })
    long = pd.DataFrame(rows)
    long.to_csv(RESULTS_TSV, sep="\t", index=False)
    print(f"\n  → {RESULTS_TSV}  ({len(long)} dataset×gene rows)")

    # Concordance summary per gene
    print("\n--- Concordance per gene (across datasets) ---")
    sums = long.groupby("gene").agg(
        n_datasets=("dataset_id", "size"),
        n_correct_dir=("matches_expected", "sum"),
        n_significant=("p", lambda s: int((s < 0.05).sum())),
        median_rho=("rho", "median"),
    )
    sums["pct_correct_dir"] = sums["n_correct_dir"] / sums["n_datasets"] * 100
    # Reorder gene rows by expected direction
    ordered = [g for g in ALL_GENES if g in sums.index]
    sums = sums.loc[ordered]
    sums["expected"] = ["+" if g in EXPECTED_UP else "-" for g in sums.index]
    sums = sums[["expected", "n_datasets", "n_correct_dir", "pct_correct_dir",
                 "n_significant", "median_rho"]]
    print(sums.to_string(float_format="%.2f"))
    sums.to_csv(SUMMARY_TSV, sep="\t")
    print(f"  → {SUMMARY_TSV}")

    # Heatmap: datasets (rows) × genes (cols), colored by ρ
    pivot = long.pivot(index="dataset_id", columns="gene", values="rho")
    pivot = pivot[[g for g in ALL_GENES if g in pivot.columns]]
    # Truncate dataset IDs for display
    pivot.index = [str(x)[:14] + "…" for x in pivot.index]
    fig, ax = plt.subplots(figsize=(max(6, 0.7 * len(pivot.columns) + 2),
                                    0.4 * len(pivot) + 2))
    vmax = float(np.nanmax(np.abs(pivot.values)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    im = ax.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    # Mark significant cells with a dot
    for (i, ds_short), ds_full in zip(enumerate(pivot.index), keep):
        for j, g in enumerate(pivot.columns):
            row_match = long[(long["dataset_id"] == ds_full) & (long["gene"] == g)]
            if not row_match.empty and row_match.iloc[0]["p"] < 0.05:
                ax.text(j, i, "•", ha="center", va="center",
                        color="black", fontsize=14, fontweight="bold")
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("Spearman ρ (donor IAD score vs gene mean)")
    ax.set_title("Per-dataset replication — donor-level ρ (• = p<0.05)\n"
                 f"Expected: {', '.join(EXPECTED_UP)} > 0  |  "
                 f"{', '.join(EXPECTED_DOWN)} < 0")
    plt.tight_layout()
    plt.savefig(FIGURE_PATH, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {FIGURE_PATH}")

    # Headline
    print("\n--- Headline ---")
    n_ds = len(keep)
    for g in ordered:
        row = sums.loc[g]
        verdict = ("REPLICATES"   if row["pct_correct_dir"] >= 75
                   else "PARTIAL" if row["pct_correct_dir"] >= 50
                   else "FAILS")
        print(f"  {g:<8} {row['expected']}  "
              f"correct direction in {int(row['n_correct_dir'])}/{n_ds} "
              f"datasets ({row['pct_correct_dir']:.0f}%)   "
              f"{int(row['n_significant'])} of those at p<0.05  "
              f"→ {verdict}")


if __name__ == "__main__":
    main()
