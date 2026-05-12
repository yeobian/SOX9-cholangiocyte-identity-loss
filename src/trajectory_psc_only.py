"""
Trajectory restricted to PSC + normal — sharper disease axis.

Why
---
The default trajectory.py mixes PBC, PSC, IFALD, CRC-metastasis, and normal
into one pseudotime ordering. PBC and PSC are biologically different
cholangiopathies; lumping them weakens the signal.

This script slices to {primary sclerosing cholangitis, normal} only and
re-runs PAGA + DPT. Expected outcome: a cleaner Spearman ρ between
donor-level pseudotime and IAD score, because the disease axis is now
mechanistically homogeneous.

Outputs
-------
results/trajectory_psc_per_cell.tsv
results/trajectory_psc_per_donor.tsv
figures/comparison/paga_graph_psc.png
figures/comparison/pseudotime_umap_psc.png
figures/comparison/pseudotime_vs_iad_psc.png

Usage
-----
    python src/trajectory_psc_only.py
"""

import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir
from iad_score import LIVER_H5AD
from trajectory import (
    run_trajectory, plot_paga, plot_pseudotime_umap,
    plot_pseudotime_vs_iad_per_donor,
)
import pandas as pd

RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures") / "comparison"

# Diseases to keep — focused PSC vs normal contrast
KEEP_DISEASES = {"primary sclerosing cholangitis", "normal"}


def slice_to_psc_normal(adata: ad.AnnData) -> ad.AnnData:
    if "disease" not in adata.obs.columns:
        print("  WARNING: no `disease` column — keeping all cells.")
        return adata
    labels = adata.obs["disease"].astype(str)
    mask = labels.isin(KEEP_DISEASES)
    n_keep = int(mask.sum())
    if n_keep == 0:
        print(f"  ERROR: no cells matched {KEEP_DISEASES}", file=sys.stderr)
        sys.exit(1)
    print(f"  Keeping {n_keep:,}/{adata.n_obs:,} cells from "
          f"{sorted(KEEP_DISEASES)}")
    print("  Per-disease counts after slice:")
    print(labels[mask].value_counts().to_string())
    return adata[mask.values].copy()


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR)
    ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Cholangiocyte trajectory — PSC + normal ONLY")
    print("=" * 60)

    if not LIVER_H5AD.exists():
        print(f"ERROR: {LIVER_H5AD} missing", file=sys.stderr)
        sys.exit(1)

    print("\n--- Loading ---")
    adata = sc.read_h5ad(LIVER_H5AD)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    print("\n--- Slicing to PSC + normal ---")
    adata = slice_to_psc_normal(adata)

    print("\n--- Trajectory ---")
    adata = run_trajectory(adata)

    print("\n--- Plotting ---")
    plot_paga(adata, FIGURES_DIR / "paga_graph_psc.png")
    plot_pseudotime_umap(adata, FIGURES_DIR / "pseudotime_umap_psc.png")

    print("\n--- Per-cell + per-donor tables ---")
    per_cell = pd.DataFrame({
        "cell":       adata.obs_names,
        "leiden":     adata.obs["leiden"].astype(str).values,
        "pseudotime": adata.obs["dpt_pseudotime"].values,
        "iad_score":  adata.obs["iad_score"].values,
        "disease":    adata.obs["disease"].astype(str).values,
        "donor":      adata.obs["donor_id"].astype(str).values,
    })
    out_cell = RESULTS_DIR / "trajectory_psc_per_cell.tsv"
    per_cell.to_csv(out_cell, sep="\t", index=False)
    print(f"  → {out_cell}")

    per_donor = (
        per_cell.groupby(["donor", "disease"], observed=True)
        .agg(mean_pseudotime=("pseudotime", "mean"),
             mean_iad_score=("iad_score", "mean"),
             n_cells=("cell", "count"))
        .reset_index()
        .sort_values("mean_pseudotime", ascending=True)
    )
    out_donor = RESULTS_DIR / "trajectory_psc_per_donor.tsv"
    per_donor.to_csv(out_donor, sep="\t", index=False)
    print(f"  → {out_donor}")

    plot_pseudotime_vs_iad_per_donor(
        per_donor, FIGURES_DIR / "pseudotime_vs_iad_psc.png"
    )

    # Console correlation
    if len(per_donor) >= 3:
        from scipy.stats import spearmanr
        clean = per_donor.replace([float("inf"), float("-inf")], pd.NA).dropna(
            subset=["mean_pseudotime", "mean_iad_score"]
        )
        if len(clean) >= 3:
            rho, p = spearmanr(clean["mean_pseudotime"], clean["mean_iad_score"])
            print(f"\nDonor-level Spearman ρ(pseudotime, IAD score) "
                  f"= {rho:+.3f}  (p = {p:.2e})  [PSC + normal only, "
                  f"n_donors = {len(clean)}]")
            print("Compare to the all-disease run: a stronger ρ here means "
                  "PBC/IFALD/CRC are noise.")

    print("\nDone.")


if __name__ == "__main__":
    main()
