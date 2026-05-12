"""
PROX1+ cholangiocyte abundance test.

Tests whether the FRACTION of PROX1+ cholangiocytes (not the per-cell
PROX1 level) changes with disease status, with the IAD score, and with
PSC pseudotime — separately from any per-cell expression effect.

Per-cell PROX1 testing already showed: PROX1 expression does not track
within-PSC pseudotime (donor ρ = +0.20, p = 0.58). The remaining open
question is whether PROX1+ cells are differentially abundant in disease,
i.e. whether disease changes the proportion of cholangiocytes occupying
the PROX1+ sub-state, even if the per-cell level among PROX1+ cells is
fixed.

Pipeline (load → group → compute → save → plot)

Outputs
-------
results/prox1_abundance_per_donor.tsv
results/prox1_abundance_summary.tsv
results/prox1_abundance_per_dataset.tsv
figures/prox1/prox1_abundance_by_disease_boxplot.png
figures/prox1/prox1_abundance_vs_iad_scatter.png
figures/prox1/prox1_abundance_vs_psc_pseudotime.png
figures/prox1/prox1_abundance_per_dataset_dotplot.png

Usage
-----
    python src/prox1_abundance_test.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression
from iad_score import (
    LIVER_H5AD as DEFAULT_CHOL,
    REFERENCE_PERCENTILE,
    _panel_expression_per_cell, _panel_genes_present, load_panel,
)

CHOL_PATH        = DEFAULT_CHOL
PSC_DONOR_TSV    = Path("results") / "trajectory_intra_psc_per_donor.tsv"
RESULTS_DIR      = Path("results")
FIGURES_DIR      = Path("figures") / "prox1"

PROX1_THRESHOLD_LOOSE  = 0.0    # > 0  (any detection)
GENE = "PROX1"

# Minimum cholangiocytes per donor to include in donor-level tests
MIN_CELLS_PER_DONOR = 20


def per_cell_iad_score(adata) -> np.ndarray:
    panel = load_panel()
    genes = _panel_genes_present(adata, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(adata, genes, weights=weights)
    ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    return 1.0 - np.clip(panel_expr / max(ref, 1e-6), 0.0, 1.0)


def aggregate_per_donor(df: pd.DataFrame, strict_thr: float) -> pd.DataFrame:
    """Per-donor PROX1+ stats. df has columns: donor, dataset, disease,
    iad, prox1."""
    rows = []
    for (donor, disease), sub in df.groupby(["donor", "disease"], observed=True):
        n = len(sub)
        n_loose = int((sub["prox1"] > PROX1_THRESHOLD_LOOSE).sum())
        n_strict = int((sub["prox1"] > strict_thr).sum())
        rows.append({
            "donor":               donor,
            "disease":             disease,
            "dataset_id":          sub["dataset"].mode().iloc[0]
                                   if not sub["dataset"].mode().empty else "unknown",
            "n_cells":             n,
            "n_prox1_pos_loose":   n_loose,
            "n_prox1_pos_strict":  n_strict,
            "frac_prox1_pos_loose":  n_loose / n,
            "frac_prox1_pos_strict": n_strict / n,
            "mean_prox1_all":      float(sub["prox1"].mean()),
            "median_prox1_all":    float(sub["prox1"].median()),
            "mean_prox1_positive": (
                float(sub.loc[sub["prox1"] > 0, "prox1"].mean())
                if (sub["prox1"] > 0).any() else 0.0
            ),
            "iad_score":           float(sub["iad"].mean()),
        })
    return pd.DataFrame(rows)


def add_psc_pseudotime(donor_table: pd.DataFrame) -> pd.DataFrame:
    """Join per-PSC-donor pseudotime if available."""
    if not PSC_DONOR_TSV.exists():
        print(f"  [info] {PSC_DONOR_TSV} not present — skipping pseudotime join.")
        donor_table["psc_pseudotime"] = np.nan
        return donor_table
    psc = pd.read_csv(PSC_DONOR_TSV, sep="\t")[["donor", "pseudotime"]]
    psc = psc.rename(columns={"pseudotime": "psc_pseudotime"})
    out = donor_table.merge(psc, on="donor", how="left")
    n_join = out["psc_pseudotime"].notna().sum()
    print(f"  Joined PSC pseudotime onto {n_join} donors.")
    return out


def disease_contrast(per_donor: pd.DataFrame, frac_col: str,
                     ref: str, contrast: str) -> dict:
    a = per_donor[per_donor["disease"] == contrast][frac_col]
    b = per_donor[per_donor["disease"] == ref][frac_col]
    if len(a) < 3 or len(b) < 3:
        return {"contrast": contrast, "ref": ref,
                "n_contrast": len(a), "n_ref": len(b),
                "median_contrast": float("nan"),
                "median_ref": float("nan"),
                "U": float("nan"), "p": float("nan")}
    u, p = mannwhitneyu(a, b, alternative="two-sided")
    return {"contrast": contrast, "ref": ref,
            "n_contrast": len(a), "n_ref": len(b),
            "median_contrast": float(a.median()),
            "median_ref":      float(b.median()),
            "U": float(u), "p": float(p)}


def per_dataset_test(df: pd.DataFrame) -> pd.DataFrame:
    """Per-dataset Spearman: donor PROX1+ fraction vs donor IAD score."""
    rows = []
    for ds, sub in df.groupby("dataset_id", observed=True):
        if sub["donor"].nunique() < 5:
            continue
        ok = sub.dropna(subset=["frac_prox1_pos_loose", "iad_score"])
        if len(ok) < 5:
            continue
        rho, p = spearmanr(ok["frac_prox1_pos_loose"], ok["iad_score"])
        rows.append({
            "dataset_id": ds,
            "n_donors":   int(len(ok)),
            "rho":        float(rho),
            "p":          float(p),
        })
    return pd.DataFrame(rows).sort_values("p")


def plot_box(per_donor: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    keep = ["normal", "primary biliary cholangitis",
            "primary sclerosing cholangitis"]
    sub = per_donor[per_donor["disease"].isin(keep)].copy()
    if sub.empty:
        return
    palette = {"normal": "#16a085",
               "primary biliary cholangitis": "#e67e22",
               "primary sclerosing cholangitis": "#c0392b"}
    sns.boxplot(data=sub, x="disease", y="frac_prox1_pos_loose",
                order=keep, palette=palette, ax=ax, fliersize=0)
    sns.stripplot(data=sub, x="disease", y="frac_prox1_pos_loose",
                  order=keep, color="black", size=4, jitter=True,
                  ax=ax, alpha=0.6)
    ax.set_ylabel("PROX1+ fraction (per donor)")
    ax.set_xlabel("")
    ax.set_xticklabels([t.get_text().replace("primary ", "")
                        for t in ax.get_xticklabels()], rotation=15)
    ax.set_title("PROX1+ cholangiocyte abundance by disease (donor-level)")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_scatter_iad(per_donor: pd.DataFrame, out: Path,
                     rho: float, p: float) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    diseases = sorted(per_donor["disease"].unique())
    palette = {d: c for d, c in zip(
        diseases,
        ["#16a085", "#e67e22", "#c0392b", "#8e44ad", "#2980b9"])}
    for d in diseases:
        sub = per_donor[per_donor["disease"] == d]
        ax.scatter(sub["iad_score"], sub["frac_prox1_pos_loose"],
                   s=np.clip(sub["n_cells"] / 4, 10, 250),
                   c=[palette[d]], alpha=0.7, label=d)
    ax.set_xlabel("Donor mean IAD score (cholangiocyte program loss)")
    ax.set_ylabel("PROX1+ fraction (per donor)")
    ax.set_title(f"PROX1+ abundance vs. IAD score per donor   "
                 f"ρ={rho:+.3f}  p={p:.2e}")
    ax.legend(fontsize=8, frameon=False, loc="best")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_scatter_pt(per_donor: pd.DataFrame, out: Path,
                    rho: float, p: float) -> None:
    sub = per_donor[
        (per_donor["disease"] == "primary sclerosing cholangitis") &
        (per_donor["psc_pseudotime"].notna())
    ]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(sub["psc_pseudotime"], sub["frac_prox1_pos_loose"],
               s=np.clip(sub["n_cells"] / 4, 10, 250),
               c="#c0392b", alpha=0.85)
    for _, r in sub.iterrows():
        ax.text(r["psc_pseudotime"], r["frac_prox1_pos_loose"], r["donor"],
                fontsize=7, alpha=0.6, ha="left", va="bottom")
    ax.set_xlabel("Donor mean intra-PSC pseudotime")
    ax.set_ylabel("PROX1+ fraction (per donor)")
    ax.set_title(f"PROX1+ abundance vs. PSC progression (PSC donors only)   "
                 f"ρ={rho:+.3f}  p={p:.2e}  n={len(sub)}")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_dataset_dotplot(per_donor: pd.DataFrame, out: Path) -> None:
    """Per-dataset PROX1+ fraction stratified by disease."""
    sub = per_donor.copy()
    sub["dataset_short"] = sub["dataset_id"].astype(str).str[:14] + "…"
    keep = ["normal", "primary biliary cholangitis",
            "primary sclerosing cholangitis"]
    sub = sub[sub["disease"].isin(keep)]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(max(8, sub["dataset_short"].nunique() * 0.5),
                                    5))
    palette = {"normal": "#16a085",
               "primary biliary cholangitis": "#e67e22",
               "primary sclerosing cholangitis": "#c0392b"}
    for d in keep:
        s = sub[sub["disease"] == d]
        if s.empty: continue
        ax.scatter(s["dataset_short"], s["frac_prox1_pos_loose"],
                   s=np.clip(s["n_cells"] / 4, 10, 250),
                   c=palette[d], alpha=0.75, label=d, edgecolor="white")
    ax.set_ylabel("PROX1+ fraction per donor")
    ax.set_xlabel("Source dataset (truncated)")
    ax.set_title("PROX1+ abundance per donor, stratified by source dataset")
    plt.xticks(rotation=45, ha="right", fontsize=8)
    ax.legend(fontsize=8, frameon=False, loc="best")
    plt.tight_layout()
    ensure_dir(out.parent); plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("PROX1+ cholangiocyte abundance test")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    print(f"\nLoading {CHOL_PATH} …")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if GENE not in adata.var_names:
        print(f"ERROR: {GENE} not in var.", file=sys.stderr); sys.exit(1)
    if "donor_id" not in adata.obs.columns or "disease" not in adata.obs.columns:
        print("ERROR: donor_id/disease missing.", file=sys.stderr); sys.exit(1)

    prox1 = get_gene_expression(adata, GENE)
    iad   = per_cell_iad_score(adata)

    # Strict threshold = top-quartile of PROX1 expression among ALL cholangiocytes
    strict_thr = float(np.percentile(prox1, 75))
    print(f"\nLoose threshold:   PROX1 > 0  "
          f"({(prox1 > 0).sum():,} cells, "
          f"{(prox1 > 0).mean()*100:.1f}%)")
    print(f"Strict threshold:  PROX1 > {strict_thr:.3f} (top quartile)  "
          f"({(prox1 > strict_thr).sum():,} cells)")

    long = pd.DataFrame({
        "donor":   adata.obs["donor_id"].astype(str).values,
        "disease": adata.obs["disease"].astype(str).values,
        "dataset": (adata.obs["dataset_id"].astype(str).values
                    if "dataset_id" in adata.obs.columns else "unknown"),
        "iad":     iad,
        "prox1":   prox1,
    })

    print("\n--- Per-donor PROX1+ stats ---")
    per_donor = aggregate_per_donor(long, strict_thr)
    per_donor = add_psc_pseudotime(per_donor)
    per_donor = per_donor[per_donor["n_cells"] >= MIN_CELLS_PER_DONOR].copy()
    print(f"  Donors after MIN_CELLS_PER_DONOR={MIN_CELLS_PER_DONOR} filter: "
          f"{len(per_donor)}")
    print("  Per-disease donor counts:")
    print(per_donor["disease"].value_counts().to_string())

    out = RESULTS_DIR / "prox1_abundance_per_donor.tsv"
    per_donor.to_csv(out, sep="\t", index=False)
    print(f"\n  → {out}")

    # -- Test 1: PSC vs normal, PBC vs normal -------------------------------
    print("\n--- Test 1: disease vs normal (Mann–Whitney on PROX1+ fraction) ---")
    contrasts = []
    for cond in ("primary sclerosing cholangitis",
                 "primary biliary cholangitis"):
        for col in ("frac_prox1_pos_loose", "frac_prox1_pos_strict"):
            res = disease_contrast(per_donor, col, "normal", cond)
            res["frac_col"] = col
            contrasts.append(res)
    print(pd.DataFrame(contrasts).to_string(index=False, float_format="%.3f"))
    pd.DataFrame(contrasts).to_csv(
        RESULTS_DIR / "prox1_abundance_summary.tsv",
        sep="\t", index=False)

    # -- Test 2: ρ(PROX1+ fraction, IAD score) across donors ---------------
    print("\n--- Test 2: PROX1+ fraction vs IAD score (across all donors) ---")
    ok = per_donor.dropna(subset=["frac_prox1_pos_loose", "iad_score"])
    rho_iad, p_iad = spearmanr(ok["frac_prox1_pos_loose"], ok["iad_score"])
    print(f"  Spearman ρ = {rho_iad:+.3f}  p = {p_iad:.2e}  n_donors = {len(ok)}")

    # -- Test 3: PROX1+ fraction vs PSC pseudotime among PSC donors --------
    print("\n--- Test 3: PROX1+ fraction vs PSC pseudotime (PSC donors only) ---")
    psc = per_donor[
        (per_donor["disease"] == "primary sclerosing cholangitis") &
        (per_donor["psc_pseudotime"].notna())
    ]
    if len(psc) >= 4:
        rho_pt, p_pt = spearmanr(psc["frac_prox1_pos_loose"],
                                 psc["psc_pseudotime"])
        print(f"  Spearman ρ = {rho_pt:+.3f}  p = {p_pt:.2e}  n_donors = {len(psc)}")
    else:
        rho_pt, p_pt = float("nan"), float("nan")
        print(f"  Only {len(psc)} PSC donors with pseudotime — skipping.")

    # -- Test 4: per-dataset stratification ---------------------------------
    print("\n--- Test 4: per-dataset Spearman (PROX1+ frac vs IAD per donor) ---")
    per_ds = per_dataset_test(per_donor)
    print(per_ds.to_string(index=False, float_format="%.3f"))
    per_ds.to_csv(RESULTS_DIR / "prox1_abundance_per_dataset.tsv",
                  sep="\t", index=False)

    # -- Plots --------------------------------------------------------------
    print("\n--- Plots ---")
    plot_box(per_donor, FIGURES_DIR / "prox1_abundance_by_disease_boxplot.png")
    plot_scatter_iad(per_donor,
                     FIGURES_DIR / "prox1_abundance_vs_iad_scatter.png",
                     rho_iad, p_iad)
    plot_scatter_pt(per_donor,
                    FIGURES_DIR / "prox1_abundance_vs_psc_pseudotime.png",
                    rho_pt, p_pt)
    plot_dataset_dotplot(per_donor,
                         FIGURES_DIR / "prox1_abundance_per_dataset_dotplot.png")

    # -- Headline -----------------------------------------------------------
    print("\n--- Headline ---")
    psc_med = per_donor[per_donor["disease"] == "primary sclerosing cholangitis"]\
        ["frac_prox1_pos_loose"].median()
    nor_med = per_donor[per_donor["disease"] == "normal"]\
        ["frac_prox1_pos_loose"].median()
    print(f"  PROX1+ median fraction: normal={nor_med:.3f}  "
          f"PSC={psc_med:.3f}  delta={psc_med - nor_med:+.3f}")
    psc_p = next((c["p"] for c in contrasts
                  if c["frac_col"] == "frac_prox1_pos_loose"
                  and c["contrast"] == "primary sclerosing cholangitis"),
                 float("nan"))
    print(f"  PSC vs normal Mann–Whitney p = {psc_p:.2e}")
    print(f"  ρ(PROX1+ frac, IAD score) across donors = {rho_iad:+.3f}, p = {p_iad:.2e}")
    if not np.isnan(rho_pt):
        print(f"  ρ(PROX1+ frac, PSC pseudotime) within PSC = {rho_pt:+.3f}, p = {p_pt:.2e}")
    n_pos_dir_ds = int((per_ds["rho"] > 0).sum())
    print(f"  Per-dataset: {n_pos_dir_ds}/{len(per_ds)} datasets have ρ>0 "
          f"(positive = more disease → more PROX1+ cholangiocytes)")


if __name__ == "__main__":
    main()
