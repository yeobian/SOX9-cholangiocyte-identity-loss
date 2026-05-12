"""
Validate the IAD score on the data we already have.

Pipeline stage: COMPUTE → SAVE.

Two checks:
  (1) Tissue specificity   — liver cholangiocytes should score LOW
                             (intact biliary program); pancreas ductal cells
                             should score HIGH (program absent by design).
  (2) Disease stratification — if the liver h5ad carries a `disease` column
                               from CELLxGENE Census (PSC, PBC, biliary atresia,
                               normal …), score distributions should rise with
                               cholestatic disease.

Outputs
-------
results/validation_summary.tsv      — per-group score statistics
figures/validation/score_by_tissue.png
figures/validation/score_by_disease.png   (only if disease labels are present)
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
from iad_score import (
    LIVER_H5AD,
    PANCREAS_H5AD,
    PANEL_PATH,
    REFERENCE_PERCENTILE,
    _panel_expression_per_cell,
    _panel_genes_present,
    load_panel,
    score_anndata,
)
from utils import ensure_dir

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS_DIR     = Path("results")
FIGURES_DIR     = Path("figures") / "validation"
SUMMARY_TSV     = RESULTS_DIR / "validation_summary.tsv"
TISSUE_PLOT     = FIGURES_DIR / "score_by_tissue.png"
DISEASE_PLOT    = FIGURES_DIR / "score_by_disease.png"


# ── Tissue-specificity check ──────────────────────────────────────────────────

def calibrate_reference(liver: ad.AnnData, panel: pd.DataFrame) -> float:
    genes = _panel_genes_present(liver, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(liver, genes, weights=weights)
    ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    print(f"  Healthy reference (liver p{REFERENCE_PERCENTILE:g}): {ref:.4f}")
    return ref


def tissue_check(
    liver: ad.AnnData, pancreas: ad.AnnData, panel: pd.DataFrame, ref: float,
) -> pd.DataFrame:
    liver_score    = score_anndata(liver,    panel=panel, healthy_reference=ref)
    pancreas_score = score_anndata(pancreas, panel=panel, healthy_reference=ref)

    rows = pd.DataFrame({
        "tissue":    ["liver"] * len(liver_score) + ["pancreas"] * len(pancreas_score),
        "iad_score": list(liver_score) + list(pancreas_score),
    })
    return rows


def plot_tissue(scores: pd.DataFrame, out: Path) -> None:
    ensure_dir(out.parent)
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.violinplot(
        data=scores, x="tissue", y="iad_score",
        hue="tissue", palette={"liver": "#16a085", "pancreas": "#2980b9"},
        inner="quartile", legend=False, ax=ax,
    )
    ax.set_title("IAD ductal-loss score — tissue specificity check")
    ax.set_ylabel("IAD score (1 = program absent)")
    ax.set_ylim(-0.05, 1.05)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Disease stratification (if labels exist) ──────────────────────────────────

def has_disease_labels(adata: ad.AnnData) -> bool:
    if "disease" not in adata.obs.columns:
        return False
    n_unique = adata.obs["disease"].nunique()
    return n_unique > 1


def disease_check(
    liver: ad.AnnData, panel: pd.DataFrame, ref: float,
) -> pd.DataFrame:
    score = score_anndata(liver, panel=panel, healthy_reference=ref)
    df = pd.DataFrame({
        "disease":   liver.obs["disease"].astype(str).values,
        "iad_score": score.values,
    })
    return df


def plot_disease(df: pd.DataFrame, out: Path) -> None:
    ensure_dir(out.parent)
    order = (
        df.groupby("disease")["iad_score"].median().sort_values().index.tolist()
    )
    fig, ax = plt.subplots(figsize=(max(7, len(order) * 1.2), 5))
    sns.boxplot(
        data=df, x="disease", y="iad_score",
        order=order, hue="disease", palette="rocket_r",
        legend=False, ax=ax, fliersize=2,
    )
    ax.set_title("IAD ductal-loss score by disease label")
    ax.set_ylabel("IAD score (1 = program absent)")
    ax.set_ylim(-0.05, 1.05)
    plt.xticks(rotation=30, ha="right")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Summary ───────────────────────────────────────────────────────────────────

def summarize_tissue(scores: pd.DataFrame) -> pd.DataFrame:
    return (
        scores.groupby("tissue")["iad_score"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
        .rename(columns={
            "mean": "mean_score", "median": "median_score",
            "std": "sd_score",   "count": "n_cells",
        })
        .assign(group_kind="tissue")
        .rename(columns={"tissue": "group"})
    )


def summarize_disease(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("disease")["iad_score"]
        .agg(["mean", "median", "std", "count"])
        .reset_index()
        .rename(columns={
            "mean": "mean_score", "median": "median_score",
            "std": "sd_score",   "count": "n_cells",
        })
        .assign(group_kind="disease")
        .rename(columns={"disease": "group"})
    )


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR)
    ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("IAD score validation")
    print("=" * 60)

    if not LIVER_H5AD.exists():
        print(f"ERROR: {LIVER_H5AD} not found", file=sys.stderr); sys.exit(1)
    if not PANCREAS_H5AD.exists():
        print(f"ERROR: {PANCREAS_H5AD} not found", file=sys.stderr); sys.exit(1)
    if not PANEL_PATH.exists():
        print(f"ERROR: {PANEL_PATH} not found — run iad_analysis.py first", file=sys.stderr)
        sys.exit(1)

    panel = load_panel()

    print("\n--- Loading data ---")
    liver    = sc.read_h5ad(LIVER_H5AD)
    pancreas = sc.read_h5ad(PANCREAS_H5AD)
    print(f"  Liver:    {liver.n_obs:,} cells")
    print(f"  Pancreas: {pancreas.n_obs:,} cells")

    print("\n--- Calibrating healthy reference ---")
    ref = calibrate_reference(liver, panel)

    print("\n--- Tissue-specificity check ---")
    tissue_df = tissue_check(liver, pancreas, panel, ref)
    plot_tissue(tissue_df, TISSUE_PLOT)
    tissue_summary = summarize_tissue(tissue_df)
    print(tissue_summary.to_string(index=False))

    summaries = [tissue_summary]

    if has_disease_labels(liver):
        print("\n--- Disease-stratification check (liver) ---")
        disease_df = disease_check(liver, panel, ref)
        print(f"  Disease labels: {sorted(disease_df['disease'].unique())}")
        plot_disease(disease_df, DISEASE_PLOT)
        disease_summary = summarize_disease(disease_df)
        print(disease_summary.to_string(index=False))
        summaries.append(disease_summary)
    else:
        print("\n[skip] Liver h5ad has no usable 'disease' column "
              "(only one value or column missing).")

    print("\n--- Saving summary ---")
    out = pd.concat(summaries, ignore_index=True)
    out.to_csv(SUMMARY_TSV, sep="\t", index=False)
    print(f"  → {SUMMARY_TSV}")

    print("\nDone.")


if __name__ == "__main__":
    main()
