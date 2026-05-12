"""
Cholangiocyte Trajectory: Healthy → Reactive → Exhausted
=========================================================

Pipeline stage: COMPUTE → SAVE.

Use diffusion-map pseudotime rooted in the healthiest cholangiocyte cluster
to order cells along a one-dimensional progression. Combined with PAGA
(partition-based graph abstraction), this gives:

  - a coarse cluster-level graph showing how cholangiocyte substates
    connect (mature ↔ reactive ↔ exhausted)
  - a per-cell pseudotime value
  - per-donor mean pseudotime, which can be correlated with the IAD score

Outputs
-------
results/trajectory_per_cell.tsv         (cell, leiden, pseudotime, iad_score)
results/trajectory_per_donor.tsv        (donor, disease, mean_pseudotime, mean_iad_score)
figures/comparison/paga_graph.png
figures/comparison/pseudotime_umap.png
figures/comparison/pseudotime_vs_iad.png
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir
from iad_score import (
    LIVER_H5AD, REFERENCE_PERCENTILE,
    _panel_expression_per_cell, _panel_genes_present, load_panel,
)

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures") / "comparison"
PER_CELL_TSV  = RESULTS_DIR / "trajectory_per_cell.tsv"
PER_DONOR_TSV = RESULTS_DIR / "trajectory_per_donor.tsv"


# ── IAD score per cell (helper, no donor aggregation) ─────────────────────────

def iad_score_per_cell(adata: ad.AnnData) -> np.ndarray:
    panel = load_panel()
    genes = _panel_genes_present(adata, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(adata, genes, weights=weights)
    ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    return 1.0 - np.clip(panel_expr / max(ref, 1e-6), 0.0, 1.0)


# ── Find healthy root cluster ─────────────────────────────────────────────────

TRAJECTORY_LEIDEN_RES = 0.15   # coarse — readable PAGA, ~5-12 clusters typical
MAX_CLUSTERS_FOR_TRAJECTORY = 15
TRAJECTORY_MAX_CELLS = 5_000   # diffmap on >10k cells is laptop-prohibitive

def find_root_cluster(adata: ad.AnnData, iad_score: np.ndarray) -> str:
    """Cluster with the LOWEST mean IAD score = most-healthy = pseudotime root."""
    df = pd.DataFrame({
        "leiden": adata.obs["leiden"].astype(str).values,
        "iad":    iad_score,
    })
    means = df.groupby("leiden")["iad"].mean().sort_values()
    counts = df.groupby("leiden").size()
    print("  Per-cluster mean IAD score (and n_cells):")
    for cl, m in means.items():
        print(f"    cluster {cl}: mean IAD = {m:.3f}   n_cells = {counts[cl]:,}")
    root = means.index[0]
    print(f"  Root cluster (lowest IAD score): {root}  (n={counts[root]:,})")
    return root


# ── Trajectory computation ────────────────────────────────────────────────────

def _stratified_subsample(adata: ad.AnnData, n_target: int) -> ad.AnnData:
    """
    Subsample to ~n_target cells, stratified by donor so every donor contributes
    proportionally. Diffmap on the full 21k cells is laptop-prohibitive; 5k is
    plenty for trajectory inference (Setty et al. 2019 use 1-2k routinely).
    """
    if adata.n_obs <= n_target:
        return adata
    if "donor_id" in adata.obs.columns:
        donor = adata.obs["donor_id"].astype(str)
        # proportional sample per donor, minimum 1 cell each
        rng = np.random.default_rng(42)
        keep_idx = []
        for d, sub in adata.obs.groupby(donor):
            n_d = len(sub)
            n_keep = max(1, int(round(n_d * n_target / adata.n_obs)))
            n_keep = min(n_keep, n_d)
            keep_idx.extend(rng.choice(sub.index.values, size=n_keep, replace=False))
        sub = adata[keep_idx].copy()
    else:
        sc.pp.subsample(adata, n_obs=n_target, random_state=42)
        sub = adata
    print(f"  Subsampled for trajectory: {sub.n_obs:,} cells "
          f"(was {adata.n_obs:,}, target {n_target:,})")
    return sub


def run_trajectory(adata: ad.AnnData) -> ad.AnnData:
    """
    Build neighbor graph (if not present), Leiden, PAGA, and DPT pseudotime.
    Pseudotime is rooted in the cluster with lowest IAD score.
    """
    # Recluster + subsample BEFORE expensive diffmap step
    adata = _stratified_subsample(adata, TRAJECTORY_MAX_CELLS)

    if "X_pca" not in adata.obsm:
        print("  Recomputing PCA on subsample …")
        # scale only HVGs to keep memory in check
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

    print("  Recomputing neighbors on subsample …")
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)

    print(f"  Re-clustering at resolution={TRAJECTORY_LEIDEN_RES} …")
    sc.tl.leiden(
        adata, resolution=TRAJECTORY_LEIDEN_RES, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )
    n_clu = adata.obs["leiden"].nunique()
    print(f"  → {n_clu} clusters.")

    iad = iad_score_per_cell(adata)
    adata.obs["iad_score"] = iad

    root_cluster = find_root_cluster(adata, iad)
    # Pick the actual cell within that cluster with the LOWEST IAD score
    cl_mask = adata.obs["leiden"].astype(str).values == root_cluster
    candidates = np.where(cl_mask)[0]
    iroot = candidates[np.argmin(iad[candidates])]
    adata.uns["iroot"] = int(iroot)
    print(f"  Root cell index: {iroot} (IAD={iad[iroot]:.3f})")

    print("  Computing diffusion map …")
    sc.tl.diffmap(adata)
    sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_diffmap")

    print("  Computing PAGA …")
    sc.tl.paga(adata, groups="leiden")
    sc.pl.paga(adata, plot=False)  # populates layout

    print("  Computing diffusion pseudotime (DPT) …")
    sc.tl.dpt(adata)

    if "X_umap" not in adata.obsm:
        sc.tl.umap(adata)

    return adata


# ── Plots ─────────────────────────────────────────────────────────────────────

def plot_paga(adata: ad.AnnData, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 6))
    sc.pl.paga(adata, color="leiden", ax=ax, show=False, threshold=0.05)
    ax.set_title("PAGA — cholangiocyte cluster connectivity")
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def plot_pseudotime_umap(adata: ad.AnnData, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(adata, color="dpt_pseudotime", ax=axes[0], show=False,
               title="Diffusion pseudotime")
    sc.pl.umap(adata, color="iad_score", ax=axes[1], show=False,
               title="IAD score (per cell)")
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def plot_pseudotime_vs_iad_per_donor(donor_df: pd.DataFrame, out: Path) -> None:
    if donor_df.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    diseases = sorted(donor_df["disease"].unique())
    palette = sns.color_palette("Set2", n_colors=len(diseases))
    cmap = dict(zip(diseases, palette))
    for d in diseases:
        sub = donor_df[donor_df["disease"] == d]
        ax.scatter(sub["mean_pseudotime"], sub["mean_iad_score"],
                   s=np.clip(sub["n_cells"] / 5, 10, 200),
                   c=[cmap[d]], alpha=0.75, label=d)
    ax.set_xlabel("Mean pseudotime per donor")
    ax.set_ylabel("Mean IAD score per donor")
    ax.set_title("Per-donor: pseudotime vs IAD score")
    ax.legend(fontsize=8, frameon=False, loc="best")
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR)
    ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Cholangiocyte trajectory (PAGA + DPT)")
    print("=" * 60)

    if not LIVER_H5AD.exists():
        print(f"ERROR: {LIVER_H5AD} missing", file=sys.stderr); sys.exit(1)

    print("\n--- Loading ---")
    adata = sc.read_h5ad(LIVER_H5AD)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    print("\n--- Trajectory ---")
    adata = run_trajectory(adata)

    print("\n--- Plotting ---")
    plot_paga(adata, FIGURES_DIR / "paga_graph.png")
    plot_pseudotime_umap(adata, FIGURES_DIR / "pseudotime_umap.png")

    print("\n--- Per-cell + per-donor tables ---")
    per_cell = pd.DataFrame({
        "cell":       adata.obs_names,
        "leiden":     adata.obs["leiden"].astype(str).values,
        "pseudotime": adata.obs["dpt_pseudotime"].values,
        "iad_score":  adata.obs["iad_score"].values,
        "disease":    (adata.obs["disease"].astype(str).values
                       if "disease" in adata.obs.columns else "unknown"),
        "donor":      (adata.obs["donor_id"].astype(str).values
                       if "donor_id" in adata.obs.columns else "unknown"),
    })
    per_cell.to_csv(PER_CELL_TSV, sep="\t", index=False)
    print(f"  → {PER_CELL_TSV}")

    per_donor = (
        per_cell.groupby(["donor", "disease"])
        .agg(mean_pseudotime=("pseudotime", "mean"),
             mean_iad_score=("iad_score", "mean"),
             n_cells=("cell", "count"))
        .reset_index()
        .sort_values("mean_pseudotime", ascending=True)
    )
    per_donor.to_csv(PER_DONOR_TSV, sep="\t", index=False)
    print(f"  → {PER_DONOR_TSV}")

    plot_pseudotime_vs_iad_per_donor(
        per_donor, FIGURES_DIR / "pseudotime_vs_iad.png"
    )

    # Console correlation
    if len(per_donor) >= 3:
        from scipy.stats import spearmanr
        rho, p = spearmanr(per_donor["mean_pseudotime"],
                           per_donor["mean_iad_score"])
        print(f"\nDonor-level Spearman ρ(pseudotime, IAD score) "
              f"= {rho:+.3f}  (p = {p:.2e})")
        print("(Positive = donors further along the trajectory have more "
              "cholangiocyte program loss, supporting the\n trajectory "
              "interpretation as healthy → reactive → exhausted.)")

    print("\nDone.")


if __name__ == "__main__":
    main()
