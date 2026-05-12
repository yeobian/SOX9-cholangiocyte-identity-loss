"""
Intra-PSC trajectory — does PROX1 rise WITHIN PSC, or only across pooled donors?

Method
------
1. Slice cholangiocytes to PSC only (no normal cells in the trajectory).
2. Recompute neighbors + UMAP + Leiden inside the PSC subset.
3. Anchor pseudotime in the PSC cluster with the lowest IAD score (the
   "least-affected" PSC cluster) — i.e. the root is intra-disease, not
   normal.
4. For each panel gene + curated regulon, compute Spearman ρ vs PSC
   pseudotime, both per-cell and donor-level.
5. Save TSVs, three single-feature plots (SOX9 / HNF1B / PROX1), and a
   correlation heatmap.

Outputs
-------
results/trajectory_intra_psc_per_cell.tsv
results/trajectory_intra_psc_per_donor.tsv
results/trajectory_intra_psc_correlations.tsv
figures/comparison/intra_psc_pseudotime_vs_SOX9.png
figures/comparison/intra_psc_pseudotime_vs_HNF1B.png
figures/comparison/intra_psc_pseudotime_vs_PROX1.png
figures/comparison/intra_psc_correlation_heatmap.png

Usage
-----
    python src/trajectory_intra_psc.py
"""

import sys
from pathlib import Path

import anndata as ad
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
from tf_activity import REGULONS

CHOL_PATH    = DEFAULT_CHOL
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "comparison"

PSC_LABEL = "primary sclerosing cholangitis"

# Genes to test directly (must be in adata.var_names)
HEADLINE_GENES = ["PROX1", "HNF4A", "FOXA2", "SOX9", "HNF1B"]
# Regulons / pathway scores to test (computed from REGULONS in tf_activity.py)
HEADLINE_REGULONS = [
    "SOX9", "HNF1B", "ONECUT1", "GATA6", "FOXA1", "FOXA2",
    "PROX1", "HNF4A",
    "BILIARY_id", "HEPATOCYTE_id", "YAP_pathway", "WNT_pathway",
]

# Recluster a bit coarser inside PSC so the root cluster is well populated
PSC_LEIDEN_RES = 0.4
N_PCS          = 30
RANDOM_STATE   = 42


def per_cell_iad_score(adata: ad.AnnData) -> np.ndarray:
    panel = load_panel()
    genes = _panel_genes_present(adata, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(adata, genes, weights=weights)
    ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    return 1.0 - np.clip(panel_expr / max(ref, 1e-6), 0.0, 1.0)


def slice_psc(adata: ad.AnnData) -> ad.AnnData:
    if "disease" not in adata.obs.columns:
        print("ERROR: no `disease` column.", file=sys.stderr); sys.exit(1)
    mask = adata.obs["disease"].astype(str) == PSC_LABEL
    n_keep = int(mask.sum())
    if n_keep == 0:
        print(f"ERROR: 0 PSC cells.", file=sys.stderr); sys.exit(1)
    n_donors = adata.obs.loc[mask, "donor_id"].astype(str).nunique() \
        if "donor_id" in adata.obs.columns else "?"
    print(f"  PSC cells: {n_keep:,} from {n_donors} donors")
    return adata[mask.values].copy()


def recluster(adata: ad.AnnData) -> ad.AnnData:
    """Re-fit PCA, neighbors, UMAP and Leiden inside PSC.

    The X_pca that came with the AnnData was computed on the full atlas
    (all diseases + normal pooled). Within PSC alone, those PCs are not
    the directions of greatest variance, and 15-NN on the inherited PCs
    leaves >90 % of PSC cells in disconnected components after diffmap.
    Re-fitting PCA on PSC + a larger neighbourhood keeps the graph
    connected.
    """
    print("  Re-fitting PCA on PSC subset …")
    n_comps = min(N_PCS, adata.n_obs - 1, adata.n_vars - 1)
    sc.tl.pca(adata, n_comps=n_comps, random_state=RANDOM_STATE,
              zero_center=False)
    print(f"  Recomputing neighbors (n_neighbors=30) …")
    sc.pp.neighbors(adata, n_neighbors=30, n_pcs=n_comps)
    sc.tl.umap(adata, random_state=RANDOM_STATE)
    sc.tl.leiden(
        adata, resolution=PSC_LEIDEN_RES, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )
    n_clu = adata.obs["leiden"].nunique()
    print(f"  → {n_clu} clusters at resolution {PSC_LEIDEN_RES}")
    return adata


def run_trajectory(adata: ad.AnnData) -> ad.AnnData:
    iad = per_cell_iad_score(adata)
    adata.obs["iad_score"] = iad

    df = pd.DataFrame({
        "leiden": adata.obs["leiden"].astype(str).values,
        "iad":    iad,
    })
    means = df.groupby("leiden")["iad"].mean().sort_values()
    counts = df.groupby("leiden").size()
    print("  Per-cluster mean IAD (sorted):")
    for cl, m in means.items():
        print(f"    cluster {cl:<3}  mean IAD = {m:.3f}   n = {int(counts[cl])}")

    root_cluster = means.index[0]
    cl_mask = adata.obs["leiden"].astype(str).values == root_cluster
    cands = np.where(cl_mask)[0]
    iroot = cands[np.argmin(iad[cands])]
    adata.uns["iroot"] = int(iroot)
    print(f"  Root PSC cluster = {root_cluster}  "
          f"(n = {int(counts[root_cluster])}, mean IAD = {means.iloc[0]:.3f})")
    print(f"  Root cell index in PSC subset = {iroot} (IAD = {iad[iroot]:.3f})")

    print("  Computing diffusion map + neighbors on diffmap + DPT …")
    sc.tl.diffmap(adata, n_comps=15)
    # Larger n_neighbors on diffmap to keep PSC graph well-connected
    sc.pp.neighbors(adata, n_neighbors=30, use_rep="X_diffmap")
    sc.tl.dpt(adata)

    # Sanity check: pseudotime should be POSITIVELY correlated with IAD score
    # (root was placed at lowest IAD; cells move away from healthy → higher IAD).
    pt = adata.obs["dpt_pseudotime"].values
    iad = adata.obs["iad_score"].values
    finite = np.isfinite(pt)
    rho_check, _ = spearmanr(pt[finite], iad[finite])
    print(f"  Sanity ρ(pseudotime, iad_score) = {rho_check:+.3f}  "
          f"({'OK' if rho_check >= 0 else 'INVERTED — flipping pseudotime'})")
    if rho_check < 0:
        # Flip so that high pseudotime = high IAD = more cholangiocyte loss
        max_pt = np.nanmax(pt[finite])
        new_pt = np.where(finite, max_pt - pt, pt)
        adata.obs["dpt_pseudotime"] = new_pt
        rho2, _ = spearmanr(new_pt[finite], iad[finite])
        print(f"  Pseudotime flipped: new ρ(pt, iad) = {rho2:+.3f}")

    n_finite = int(np.isfinite(adata.obs["dpt_pseudotime"]).sum())
    print(f"  Cells with finite pseudotime: {n_finite:,}/{adata.n_obs:,} "
          f"({n_finite/adata.n_obs*100:.1f}%)")
    return adata


def score_regulons_local(adata: ad.AnnData) -> dict:
    """Per-cell score for each regulon, restricted to genes present in adata."""
    out = {}
    for tf, targets in REGULONS.items():
        if tf not in HEADLINE_REGULONS:
            continue
        present = [g for g in targets if g in adata.var_names]
        if len(present) < 2:
            print(f"  [skip] regulon {tf}: only {len(present)}/{len(targets)} targets present.")
            continue
        key = f"_score_{tf}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          random_state=RANDOM_STATE)
        out[tf] = adata.obs[key].values.copy()
        del adata.obs[key]
        print(f"  scored {tf:<14}  ({len(present)}/{len(targets)} targets)")
    return out


def correlations(pseudotime: np.ndarray, donor: np.ndarray,
                 features: dict) -> pd.DataFrame:
    """Per-cell + per-donor Spearman ρ of pseudotime against each feature."""
    rows = []
    pt_finite = np.isfinite(pseudotime)
    pt = pseudotime[pt_finite]
    donor_f = donor[pt_finite]
    for name, vec in features.items():
        v = vec[pt_finite]
        ok = ~np.isnan(v)
        if ok.sum() < 30:
            cell_rho, cell_p = float("nan"), float("nan")
        else:
            cell_rho, cell_p = spearmanr(pt[ok], v[ok])
        # donor-level: mean per donor, then spearman
        df = pd.DataFrame({"donor": donor_f, "pt": pt, "x": v}).dropna()
        donor_means = df.groupby("donor").agg(pt=("pt", "mean"),
                                              x=("x", "mean"))
        if len(donor_means) >= 4:
            donor_rho, donor_p = spearmanr(donor_means["pt"],
                                           donor_means["x"])
            n_donors = len(donor_means)
        else:
            donor_rho, donor_p, n_donors = float("nan"), float("nan"), len(donor_means)
        rows.append({
            "feature":   name,
            "cell_rho":  cell_rho,
            "cell_p":    cell_p,
            "donor_rho": donor_rho,
            "donor_p":   donor_p,
            "n_cells":   int(ok.sum()),
            "n_donors":  n_donors,
        })
    return pd.DataFrame(rows)


def smooth(x: np.ndarray, y: np.ndarray, n_bins: int = 25) -> tuple:
    order = np.argsort(x)
    x_s = x[order]; y_s = y[order]
    bins = np.linspace(x_s.min(), x_s.max(), n_bins + 1)
    centers = 0.5 * (bins[:-1] + bins[1:])
    means   = np.full(n_bins, np.nan)
    for i in range(n_bins):
        m = (x_s >= bins[i]) & (x_s < bins[i + 1])
        if m.any():
            means[i] = y_s[m].mean()
    keep = ~np.isnan(means)
    return centers[keep], means[keep]


def single_feature_plot(pt: np.ndarray, y: np.ndarray, name: str,
                        rho: float, p: float, out: Path) -> None:
    finite = np.isfinite(pt) & np.isfinite(y)
    pt_f = pt[finite]; y_f = y[finite]
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(pt_f, y_f, s=4, c="#34495e", alpha=0.25)
    if len(pt_f) > 50:
        xc, ym = smooth(pt_f, y_f, n_bins=24)
        ax.plot(xc, ym, color="#c0392b", linewidth=2.4,
                label=f"binned mean  (cell-level ρ={rho:+.3f}, p={p:.1e})")
    ax.set_xlabel("Intra-PSC pseudotime  (least-affected → most-affected)")
    ax.set_ylabel(f"{name} (log-normalized expr or regulon score)")
    ax.set_title(f"PSC-only trajectory — {name}")
    ax.legend(loc="best", fontsize=9, frameon=False)
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def heatmap(corr: pd.DataFrame, out: Path) -> None:
    if corr.empty:
        return
    df = corr.set_index("feature")[["cell_rho", "donor_rho"]]
    fig, ax = plt.subplots(figsize=(5, 0.4 * len(df) + 1.5))
    vmax = float(np.nanmax(np.abs(df.values)))
    if not np.isfinite(vmax) or vmax == 0:
        vmax = 1.0
    im = ax.imshow(df.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["cell-level", "donor-level"])
    ax.set_yticks(range(len(df))); ax.set_yticklabels(df.index, fontsize=9)
    # significance dots
    for i, (feat, row) in enumerate(df.iterrows()):
        cell_p  = corr.set_index("feature").loc[feat, "cell_p"]
        donor_p = corr.set_index("feature").loc[feat, "donor_p"]
        if not np.isnan(cell_p) and cell_p < 0.05:
            ax.text(0, i, "•", ha="center", va="center",
                    color="black", fontsize=14, fontweight="bold")
        if not np.isnan(donor_p) and donor_p < 0.05:
            ax.text(1, i, "•", ha="center", va="center",
                    color="black", fontsize=14, fontweight="bold")
        for j, v in enumerate(row):
            if not np.isnan(v):
                ax.text(j, i, f"{v:+.2f}",
                        ha="center", va="center", fontsize=8,
                        color="white" if abs(v) > 0.4 else "black")
    plt.colorbar(im, ax=ax, fraction=0.05, pad=0.04,
                 label="Spearman ρ vs PSC pseudotime")
    ax.set_title("Intra-PSC trajectory correlations\n(• = p<0.05)")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR)
    ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Intra-PSC trajectory  (PSC cholangiocytes only)")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)

    print("\n[1/6] Loading + slicing …")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  full: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    adata = slice_psc(adata)

    print("\n[2/6] Re-clustering inside PSC …")
    adata = recluster(adata)

    print("\n[3/6] Trajectory (root in least-affected PSC cluster) …")
    adata = run_trajectory(adata)

    print("\n[4/6] Scoring regulons on PSC cells …")
    regulon_scores = score_regulons_local(adata)

    print("\n[5/6] Computing correlations …")
    pt = adata.obs["dpt_pseudotime"].values
    donor = (adata.obs["donor_id"].astype(str).values
             if "donor_id" in adata.obs.columns
             else np.array(["unknown"] * adata.n_obs))

    features: dict = {}
    for g in HEADLINE_GENES:
        if g in adata.var_names:
            features[f"{g} (mRNA)"] = get_gene_expression(adata, g)
    for tf, vec in regulon_scores.items():
        features[f"{tf} (regulon)"] = vec

    corr = correlations(pt, donor, features)
    print("\n--- Spearman ρ (sorted by donor-level ρ) ---")
    print(corr.sort_values("donor_rho").to_string(
        index=False, float_format="%.3f"))
    out_corr = RESULTS_DIR / "trajectory_intra_psc_correlations.tsv"
    corr.to_csv(out_corr, sep="\t", index=False)
    print(f"\n  → {out_corr}")

    # Per-cell + per-donor tables
    pt_finite = np.isfinite(pt)
    per_cell = pd.DataFrame({
        "cell":       adata.obs_names,
        "leiden":     adata.obs["leiden"].astype(str).values,
        "pseudotime": pt,
        "iad_score":  adata.obs["iad_score"].values,
        "donor":      donor,
    })
    for fname, vec in features.items():
        per_cell[fname] = vec
    per_cell = per_cell.loc[pt_finite].copy()
    out_cell = RESULTS_DIR / "trajectory_intra_psc_per_cell.tsv"
    per_cell.to_csv(out_cell, sep="\t", index=False)
    print(f"  → {out_cell}")

    per_donor = (
        per_cell.groupby("donor")
        .agg({"pseudotime": "mean", "iad_score": "mean",
              **{f: "mean" for f in features.keys()}})
        .reset_index()
        .sort_values("pseudotime")
    )
    per_donor["n_cells"] = per_cell.groupby("donor").size().reindex(
        per_donor["donor"]).values
    out_donor = RESULTS_DIR / "trajectory_intra_psc_per_donor.tsv"
    per_donor.to_csv(out_donor, sep="\t", index=False)
    print(f"  → {out_donor}")

    print("\n[6/6] Plots …")
    for g in ("SOX9", "HNF1B", "PROX1"):
        # Prefer mRNA if available else regulon
        mrna_key = f"{g} (mRNA)"
        reg_key  = f"{g} (regulon)"
        target_key = mrna_key if mrna_key in features else (
            reg_key if reg_key in features else None)
        if target_key is None:
            continue
        row = corr[corr["feature"] == target_key].iloc[0]
        single_feature_plot(
            pt, features[target_key], target_key,
            row["cell_rho"], row["cell_p"],
            FIGURES_DIR / f"intra_psc_pseudotime_vs_{g}.png",
        )
    heatmap(corr, FIGURES_DIR / "intra_psc_correlation_heatmap.png")

    print("\n--- Headline ---")
    sox9 = corr[corr["feature"] == "SOX9 (mRNA)"]
    hnf1b = corr[corr["feature"] == "HNF1B (mRNA)"]
    prox1 = corr[corr["feature"] == "PROX1 (mRNA)"]
    if not sox9.empty:
        r = sox9.iloc[0]
        print(f"  SOX9 mRNA  cell ρ={r['cell_rho']:+.3f} p={r['cell_p']:.1e}  "
              f"donor ρ={r['donor_rho']:+.3f} p={r['donor_p']:.1e}")
    if not hnf1b.empty:
        r = hnf1b.iloc[0]
        print(f"  HNF1B mRNA cell ρ={r['cell_rho']:+.3f} p={r['cell_p']:.1e}  "
              f"donor ρ={r['donor_rho']:+.3f} p={r['donor_p']:.1e}")
    if not prox1.empty:
        r = prox1.iloc[0]
        print(f"  PROX1 mRNA cell ρ={r['cell_rho']:+.3f} p={r['cell_p']:.1e}  "
              f"donor ρ={r['donor_rho']:+.3f} p={r['donor_p']:.1e}")


if __name__ == "__main__":
    main()
