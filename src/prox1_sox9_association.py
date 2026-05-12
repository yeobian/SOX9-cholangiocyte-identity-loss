"""
PROX1 ↔ SOX9 association test in cholangiocytes.

Question
--------
Are PROX1+ cholangiocytes the cells that preserve cholangiocyte identity
(SOX9-high, KRT19+, KRT7+, biliary-identity-high), or are they the cells
that have collapsed into a failed-plasticity intermediate state (SOX9-low,
biliary-identity-low)?

Tests
-----
1. Per-cell Spearman ρ(PROX1, SOX9) — all cholangiocytes, then PSC only.
2. Donor-level Spearman ρ(mean PROX1, mean SOX9) — all, then PSC only.
3. PROX1-high vs PROX1-low (top-quartile vs bottom-three-quartiles):
   - SOX9 expression (Mann–Whitney)
   - Biliary identity score
   - KRT19+, KRT7+, EPCAM+ positivity rates
   - SOX9-high (top-quartile) positivity rate
4. Per-dataset stratified Spearman (per-cell PROX1 vs SOX9).

Outputs
-------
results/prox1_sox9_per_cell_corr.tsv
results/prox1_sox9_per_donor.tsv
results/prox1_sox9_marker_summary.tsv
results/prox1_sox9_per_dataset.tsv
figures/prox1/prox1_vs_sox9_scatter.png
figures/prox1/prox1_high_low_sox9_boxplot.png
figures/prox1/prox1_high_low_marker_bars.png
figures/prox1/prox1_sox9_per_dataset_heatmap.png

Usage
-----
    python src/prox1_sox9_association.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import spearmanr, mannwhitneyu

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression
from iad_score import LIVER_H5AD as DEFAULT_CHOL

CHOL_PATH    = DEFAULT_CHOL
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "prox1"

PSC = "primary sclerosing cholangitis"
PROX1 = "PROX1"
SOX9  = "SOX9"
MARKERS = ["KRT19", "KRT7", "EPCAM"]

# Biliary identity from tf_activity.REGULONS
BILIARY_ID_GENES = ["KRT19", "KRT7", "EPCAM", "SOX9", "HNF1B", "AQP1"]

MIN_CELLS_PER_DATASET = 200
MIN_DONORS_PER_DATASET = 5


def load_cells():
    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    expr = {}
    for g in [PROX1, SOX9] + MARKERS:
        if g in adata.var_names:
            expr[g] = get_gene_expression(adata, g)
        else:
            print(f"  WARN: {g} not in var_names — will skip its tests")

    # Biliary identity score on the full cholangiocyte AnnData
    present = [g for g in BILIARY_ID_GENES if g in adata.var_names]
    sc.tl.score_genes(adata, gene_list=present, score_name="bil_id",
                      random_state=42)
    bil_id = adata.obs["bil_id"].values.copy()
    del adata.obs["bil_id"]
    print(f"  scored BILIARY_id on {len(present)}/{len(BILIARY_ID_GENES)} genes")

    df = pd.DataFrame({
        "donor":   adata.obs["donor_id"].astype(str).values,
        "disease": adata.obs["disease"].astype(str).values,
        "dataset": adata.obs["dataset_id"].astype(str).values
                   if "dataset_id" in adata.obs.columns else "unknown",
        **expr,
        "bil_id":  bil_id,
    })
    return df


def per_cell_corr(df: pd.DataFrame, label: str) -> dict:
    rho, p = spearmanr(df[PROX1], df[SOX9])
    print(f"    [{label}]  per-cell ρ(PROX1,SOX9) = {rho:+.3f}  p = {p:.2e}  "
          f"n = {len(df):,}")
    return {"slice": label, "level": "cell",
            "n": int(len(df)), "rho": float(rho), "p": float(p)}


def per_donor_corr(df: pd.DataFrame, label: str) -> dict:
    grp = df.groupby("donor").agg(p=(PROX1, "mean"), s=(SOX9, "mean"),
                                  n=("donor", "size"))
    grp = grp[grp["n"] >= 20]
    if len(grp) < 4:
        return {"slice": label, "level": "donor",
                "n": int(len(grp)), "rho": float("nan"), "p": float("nan")}
    rho, p = spearmanr(grp["p"], grp["s"])
    print(f"    [{label}]  donor  ρ(mean PROX1, mean SOX9) = {rho:+.3f}  "
          f"p = {p:.2e}  n_donors = {len(grp)}")
    return {"slice": label, "level": "donor",
            "n": int(len(grp)), "rho": float(rho), "p": float(p)}


def high_low_contrast(df: pd.DataFrame, label: str,
                      strict_thr: float) -> pd.DataFrame:
    """Top-quartile PROX1 vs bottom-three-quartiles. Compare SOX9, BILIARY_id,
    and marker positivity rates."""
    high_mask = df[PROX1] > strict_thr
    high = df[high_mask]
    low  = df[~high_mask]
    if len(high) < 30 or len(low) < 30:
        return pd.DataFrame()

    rows = []

    # SOX9 expression
    u, p = mannwhitneyu(high[SOX9], low[SOX9], alternative="two-sided")
    rows.append({"slice": label, "metric": "SOX9 expression",
                 "high_mean":   float(high[SOX9].mean()),
                 "low_mean":    float(low[SOX9].mean()),
                 "high_median": float(high[SOX9].median()),
                 "low_median":  float(low[SOX9].median()),
                 "n_high": int(len(high)), "n_low": int(len(low)),
                 "U": float(u), "p": float(p)})

    # BILIARY_id score
    u, p = mannwhitneyu(high["bil_id"], low["bil_id"], alternative="two-sided")
    rows.append({"slice": label, "metric": "BILIARY_id score",
                 "high_mean":   float(high["bil_id"].mean()),
                 "low_mean":    float(low["bil_id"].mean()),
                 "high_median": float(high["bil_id"].median()),
                 "low_median":  float(low["bil_id"].median()),
                 "n_high": int(len(high)), "n_low": int(len(low)),
                 "U": float(u), "p": float(p)})

    # Marker positivity rates (fraction of cells expressing > 0)
    sox9_high_thr = float(np.percentile(df[SOX9], 75))
    for m in [SOX9 + "-high (top quartile)"] + [g + "+" for g in MARKERS]:
        if "high" in m:
            high_pct = float((high[SOX9] > sox9_high_thr).mean() * 100)
            low_pct  = float((low[SOX9]  > sox9_high_thr).mean() * 100)
        else:
            gene = m[:-1]
            if gene not in df.columns: continue
            high_pct = float((high[gene] > 0).mean() * 100)
            low_pct  = float((low[gene]  > 0).mean() * 100)
        rows.append({"slice": label, "metric": f"%{m}",
                     "high_mean":   high_pct,
                     "low_mean":    low_pct,
                     "high_median": float("nan"),
                     "low_median":  float("nan"),
                     "n_high": int(len(high)), "n_low": int(len(low)),
                     "U": float("nan"), "p": float("nan")})

    out = pd.DataFrame(rows)
    print(f"\n  [{label}]  PROX1-high (top quartile, n={len(high):,}) vs "
          f"PROX1-low (n={len(low):,}):")
    for _, r in out.iterrows():
        if not np.isnan(r["p"]):
            print(f"    {r['metric']:<22}  high mean={r['high_mean']:.3f}  "
                  f"low mean={r['low_mean']:.3f}   p={r['p']:.2e}")
        else:
            print(f"    {r['metric']:<22}  high {r['high_mean']:.1f}%   "
                  f"low {r['low_mean']:.1f}%")
    return out


def per_dataset_corr(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for ds, sub in df.groupby("dataset", observed=True):
        if len(sub) < MIN_CELLS_PER_DATASET: continue
        if sub["donor"].nunique() < MIN_DONORS_PER_DATASET: continue
        rho, p = spearmanr(sub[PROX1], sub[SOX9])
        rows.append({"dataset_id": ds, "n_cells": int(len(sub)),
                     "n_donors":   int(sub["donor"].nunique()),
                     "rho_per_cell": float(rho), "p": float(p)})
    return pd.DataFrame(rows).sort_values("p")


def plot_scatter(df: pd.DataFrame, out: Path) -> None:
    keep = ["normal", "primary sclerosing cholangitis"]
    sub = df[df["disease"].isin(keep)].copy()
    if sub.empty: return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True, sharex=True)
    palette = {"normal": "#16a085",
               "primary sclerosing cholangitis": "#c0392b"}
    for ax, d in zip(axes, keep):
        s = sub[sub["disease"] == d]
        ax.scatter(s[PROX1], s[SOX9], s=4, alpha=0.25, c=palette[d])
        rho, p = spearmanr(s[PROX1], s[SOX9])
        ax.set_title(f"{d}\nper-cell ρ={rho:+.3f}, p={p:.1e}, n={len(s):,}",
                     fontsize=10)
        ax.set_xlabel(f"{PROX1} (log-normalized)")
        ax.set_ylabel(f"{SOX9} (log-normalized)")
        ax.axvline(0, color="grey", linewidth=0.4)
        ax.axhline(0, color="grey", linewidth=0.4)
    fig.suptitle("PROX1 vs SOX9 in cholangiocytes — per cell")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_box(df: pd.DataFrame, strict_thr: float, out: Path) -> None:
    keep = ["normal", "primary sclerosing cholangitis"]
    sub = df[df["disease"].isin(keep)].copy()
    if sub.empty: return
    sub["PROX1_group"] = np.where(sub[PROX1] > strict_thr, "PROX1-high", "PROX1-low")
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, metric, ylabel in zip(
        axes, [SOX9, "bil_id"], ["SOX9 expression", "BILIARY_id score"]
    ):
        sns.boxplot(data=sub, x="disease", y=metric, hue="PROX1_group",
                    palette={"PROX1-high": "#8e44ad", "PROX1-low": "#bdc3c7"},
                    ax=ax, fliersize=0)
        ax.set_xticklabels([t.get_text().replace("primary ", "")
                            for t in ax.get_xticklabels()], rotation=15)
        ax.set_ylabel(ylabel); ax.set_xlabel("")
        ax.legend(fontsize=8, frameon=False, loc="best")
    fig.suptitle("PROX1-high vs PROX1-low cholangiocytes")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_marker_bars(df: pd.DataFrame, strict_thr: float, out: Path) -> None:
    keep = ["normal", "primary sclerosing cholangitis"]
    sub = df[df["disease"].isin(keep)].copy()
    if sub.empty: return
    sub["PROX1_group"] = np.where(sub[PROX1] > strict_thr, "high", "low")

    sox9_thr = float(np.percentile(df[SOX9], 75))
    rows = []
    for d in keep:
        for grp in ["high", "low"]:
            cells = sub[(sub["disease"] == d) & (sub["PROX1_group"] == grp)]
            if cells.empty: continue
            rows.append({"disease": d, "group": f"PROX1-{grp}",
                         "n":      len(cells),
                         "%KRT19+":  float((cells["KRT19"] > 0).mean() * 100)
                                       if "KRT19" in cells else np.nan,
                         "%KRT7+":   float((cells["KRT7"]  > 0).mean() * 100)
                                       if "KRT7"  in cells else np.nan,
                         "%EPCAM+":  float((cells["EPCAM"] > 0).mean() * 100)
                                       if "EPCAM" in cells else np.nan,
                         "%SOX9-high (top Q)":
                                     float((cells[SOX9] > sox9_thr).mean()*100),
                        })
    table = pd.DataFrame(rows)

    metrics = ["%KRT19+", "%KRT7+", "%EPCAM+", "%SOX9-high (top Q)"]
    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4),
                             sharey=True)
    palette = {"PROX1-high": "#8e44ad", "PROX1-low": "#bdc3c7"}
    for ax, m in zip(axes, metrics):
        sns.barplot(data=table, x="disease", y=m, hue="group",
                    palette=palette, ax=ax)
        ax.set_xticklabels([t.get_text().replace("primary ", "")
                            for t in ax.get_xticklabels()], rotation=15,
                           fontsize=8)
        ax.set_ylim(0, 100)
        ax.set_ylabel(m if m == metrics[0] else "")
        ax.set_xlabel("")
        if m == metrics[0]:
            ax.legend(fontsize=8, frameon=False)
        else:
            ax.get_legend().remove()
    fig.suptitle("Marker positivity: PROX1-high vs PROX1-low cholangiocytes")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_dataset_heatmap(per_ds: pd.DataFrame, out: Path) -> None:
    if per_ds.empty: return
    fig, ax = plt.subplots(figsize=(5, 0.4 * len(per_ds) + 1.5))
    short = [s[:14] + "…" for s in per_ds["dataset_id"]]
    vmax = float(np.nanmax(np.abs(per_ds["rho_per_cell"])))
    if not np.isfinite(vmax) or vmax == 0: vmax = 1.0
    rho_arr = per_ds["rho_per_cell"].values.reshape(-1, 1)
    im = ax.imshow(rho_arr, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                   aspect="auto")
    ax.set_yticks(range(len(per_ds))); ax.set_yticklabels(short, fontsize=8)
    ax.set_xticks([0]); ax.set_xticklabels(["per-cell ρ"])
    for i, (rho, p) in enumerate(zip(per_ds["rho_per_cell"], per_ds["p"])):
        ax.text(0, i, f"{rho:+.2f}\np={p:.1e}",
                ha="center", va="center", fontsize=8,
                color="white" if abs(rho) > 0.4 else "black")
    plt.colorbar(im, ax=ax, fraction=0.05, pad=0.04,
                 label=f"ρ({PROX1}, {SOX9})")
    ax.set_title("Per-dataset PROX1 ↔ SOX9 (per-cell Spearman)")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("PROX1 ↔ SOX9 association in cholangiocytes")
    print("=" * 60)

    print("\n[1/5] Loading + scoring …")
    df = load_cells()
    print(f"  cells loaded: {len(df):,}")

    strict_thr = float(np.percentile(df[PROX1], 75))
    print(f"\n[2/5] PROX1-high threshold (top quartile) = {strict_thr:.3f}")

    print("\n[3/5] Correlations  (per-cell + per-donor)")
    rows = []
    for label, sub in [("ALL cholangiocytes", df),
                       ("Normal cholangiocytes", df[df["disease"] == "normal"]),
                       ("PSC cholangiocytes",   df[df["disease"] == PSC])]:
        if len(sub) < 30:
            print(f"    [{label}] only {len(sub)} cells — skipping")
            continue
        rows.append(per_cell_corr(sub, label))
        rows.append(per_donor_corr(sub, label))
    pd.DataFrame(rows).to_csv(RESULTS_DIR / "prox1_sox9_per_cell_corr.tsv",
                              sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'prox1_sox9_per_cell_corr.tsv'}")

    # Donor-level table for plot/inspection
    grp = (df.groupby(["donor", "disease"], observed=True)
             .agg(n_cells=(PROX1, "size"),
                  mean_PROX1=(PROX1, "mean"),
                  mean_SOX9=(SOX9, "mean"),
                  mean_BILID=("bil_id", "mean"))
             .reset_index())
    grp.to_csv(RESULTS_DIR / "prox1_sox9_per_donor.tsv", sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'prox1_sox9_per_donor.tsv'}")

    print("\n[4/5] PROX1-high vs PROX1-low contrasts")
    contrasts = []
    for label, sub in [("ALL", df),
                       ("Normal", df[df["disease"] == "normal"]),
                       ("PSC",   df[df["disease"] == PSC])]:
        out = high_low_contrast(sub, label, strict_thr)
        if not out.empty:
            contrasts.append(out)
    if contrasts:
        all_contrasts = pd.concat(contrasts, ignore_index=True)
        all_contrasts.to_csv(RESULTS_DIR / "prox1_sox9_marker_summary.tsv",
                             sep="\t", index=False)
        print(f"\n  → {RESULTS_DIR / 'prox1_sox9_marker_summary.tsv'}")

    print("\n[5/5] Per-dataset (per-cell ρ)")
    per_ds = per_dataset_corr(df)
    print(per_ds.to_string(index=False, float_format="%.3f"))
    per_ds.to_csv(RESULTS_DIR / "prox1_sox9_per_dataset.tsv",
                  sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'prox1_sox9_per_dataset.tsv'}")

    # Plots
    print("\n--- Plots ---")
    plot_scatter(df, FIGURES_DIR / "prox1_vs_sox9_scatter.png")
    plot_box(df, strict_thr,
             FIGURES_DIR / "prox1_high_low_sox9_boxplot.png")
    plot_marker_bars(df, strict_thr,
                     FIGURES_DIR / "prox1_high_low_marker_bars.png")
    plot_dataset_heatmap(per_ds,
                         FIGURES_DIR / "prox1_sox9_per_dataset_heatmap.png")

    # Headline interpretation
    print("\n--- Headline ---")
    rho_all_cell = next((r for r in rows if r["slice"] == "ALL cholangiocytes"
                         and r["level"] == "cell"), {}).get("rho", np.nan)
    rho_psc_cell = next((r for r in rows if r["slice"] == "PSC cholangiocytes"
                         and r["level"] == "cell"), {}).get("rho", np.nan)
    print(f"  per-cell ρ(PROX1,SOX9):  ALL = {rho_all_cell:+.3f}   "
          f"PSC = {rho_psc_cell:+.3f}")
    n_pos = int((per_ds["rho_per_cell"] > 0).sum())
    print(f"  per-dataset: {n_pos}/{len(per_ds)} datasets show ρ > 0")
    if rho_psc_cell > 0.1:
        verdict = ("PROX1+ cells trend SOX9-PRESERVED within PSC — "
                   "consistent with rescue/identity-preserving sub-state.")
    elif rho_psc_cell < -0.1:
        verdict = ("PROX1+ cells trend SOX9-COLLAPSED within PSC — "
                   "consistent with failed-plasticity intermediate.")
    else:
        verdict = ("PROX1 and SOX9 vary largely INDEPENDENTLY in PSC — "
                   "PROX1 likely an independent heterogeneity marker.")
    print(f"  Verdict (within PSC): {verdict}")


if __name__ == "__main__":
    main()
