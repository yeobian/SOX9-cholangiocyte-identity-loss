"""
PROX1 ↔ IAD Link Analysis
==========================

Pipeline stage: COMPUTE → SAVE.

Aim: connect the diagnostic-panel work for IAD to a candidate
regenerative-therapy mechanism.

Logic
-----
1. PROX1 is a regeneration-suppressor transcription factor in adult mammalian
   retina (suppresses Müller-glia → neuron conversion). Knocking it down
   permits regenerative reprogramming.
2. In adult liver, PROX1 is a hepatocyte-lineage TF. Hepatocyte ↔ cholangiocyte
   plasticity is well-documented; hepatocytes can transdifferentiate into
   cholangiocytes after biliary injury.
3. Hypothesis: if PROX1 suppression releases retinal regeneration, the
   reciprocal PROX1/KRT19 axis in liver may be a tractable target for
   biliary regeneration in IAD.

Quantitative evidence we can extract from the data:
  - PROX1 expression frequency / level in cholangiocytes (should be LOW)
  - PROX1 expression frequency / level in hepatocytes (should be HIGH)
  - PROX1+ KRT19+ "transition" cells (candidate progenitors)
  - Inverse correlation between PROX1 score and IAD score per donor

Outputs
-------
data/liver/liver_hepatocytes.h5ad
results/prox1_iad_link.tsv
figures/prox1/prox1_vs_krt19.png
figures/prox1/prox1_score_vs_iad_score.png
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
from utils import ensure_dir, get_gene_expression
from iad_score import (
    LIVER_H5AD as CHOL_H5AD,
    PANEL_PATH,
    REFERENCE_PERCENTILE,
    _panel_expression_per_cell,
    _panel_genes_present,
    load_panel,
)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR     = Path("data")
HEP_RAW      = DATA_DIR / "liver" / "liver_hepatocytes_raw.h5ad"
HEP_PROC     = DATA_DIR / "liver" / "liver_hepatocytes.h5ad"
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "prox1"
LINK_TSV     = RESULTS_DIR / "prox1_iad_link.tsv"

CENSUS_VERSION = "stable"
HEP_TYPES = [
    "hepatocyte",
    "centrilobular region hepatocyte",
    "midzonal region hepatocyte",
    "periportal region hepatocyte",
]
HEP_MAX_CELLS = 5_000

# QC matches preprocess_liver.py
MIN_GENES, MAX_GENES, MAX_PCT_MT, MIN_CELLS = 200, 8_000, 25, 3


# ── Hepatocyte fetch ──────────────────────────────────────────────────────────

def fetch_hepatocytes() -> ad.AnnData:
    """Pull a small hepatocyte slice from Census so we can compare PROX1 levels."""
    if HEP_RAW.exists():
        print(f"  Reusing {HEP_RAW}")
        return sc.read_h5ad(HEP_RAW)

    try:
        import cellxgene_census
    except ImportError:
        print("ERROR: cellxgene-census not installed.", file=sys.stderr)
        sys.exit(1)

    quoted = ", ".join(f"'{t}'" for t in HEP_TYPES)
    obs_filter = f"tissue_general == 'liver' and cell_type in [{quoted}]"
    print(f"  Census filter: {obs_filter}")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            obs_value_filter=obs_filter,
            obs_column_names=["cell_type", "tissue_general", "disease",
                              "donor_id", "dataset_id"],
        )
    print(f"  Fetched {adata.n_obs:,} hepatocytes × {adata.n_vars:,} genes")

    if adata.n_obs > HEP_MAX_CELLS:
        sc.pp.subsample(adata, n_obs=HEP_MAX_CELLS, random_state=42)
        print(f"  Subsampled → {adata.n_obs:,}")

    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()

    ensure_dir(HEP_RAW.parent)
    adata.write_h5ad(HEP_RAW)
    return adata


# ── QC + normalize hepatocytes (lightweight) ──────────────────────────────────

def quick_preprocess(adata: ad.AnnData) -> ad.AnnData:
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    adata = adata[adata.obs.n_genes_by_counts <= MAX_GENES].copy()
    adata = adata[adata.obs.pct_counts_mt <= MAX_PCT_MT].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


# ── PROX1 / KRT19 quantification ──────────────────────────────────────────────

def quantify_marker(adata: ad.AnnData, gene: str) -> dict:
    if gene not in adata.var_names:
        return {"present": False, "pct_pos": 0.0, "mean_expr_pos": 0.0,
                "mean_expr_all": 0.0, "n_cells": adata.n_obs}
    expr = get_gene_expression(adata, gene)
    pct  = float((expr > 0).mean() * 100)
    mean_pos = float(expr[expr > 0].mean()) if pct > 0 else 0.0
    return {
        "present": True,
        "pct_pos": pct,
        "mean_expr_pos": mean_pos,
        "mean_expr_all": float(expr.mean()),
        "n_cells": adata.n_obs,
    }


def transition_cell_count(adata: ad.AnnData, threshold: float = 0.0) -> int:
    """Count cells co-expressing PROX1 and KRT19 above `threshold` (log-norm)."""
    if "PROX1" not in adata.var_names or "KRT19" not in adata.var_names:
        return 0
    p = get_gene_expression(adata, "PROX1")
    k = get_gene_expression(adata, "KRT19")
    return int(((p > threshold) & (k > threshold)).sum())


# ── Joint scatter / heatmap-of-counts ─────────────────────────────────────────

def plot_prox1_vs_krt19(
    chol: ad.AnnData, hep: ad.AnnData, out: Path,
) -> None:
    rows = []
    for adata, label in [(chol, "Cholangiocyte"), (hep, "Hepatocyte")]:
        if "PROX1" not in adata.var_names or "KRT19" not in adata.var_names:
            continue
        p = get_gene_expression(adata, "PROX1")
        k = get_gene_expression(adata, "KRT19")
        for pp, kk in zip(p, k):
            rows.append({"cell_type": label, "PROX1": pp, "KRT19": kk})
    if not rows:
        return
    df = pd.DataFrame(rows)

    fig, ax = plt.subplots(figsize=(7, 6))
    palette = {"Cholangiocyte": "#16a085", "Hepatocyte": "#c0392b"}
    for label, sub in df.groupby("cell_type"):
        ax.scatter(sub["PROX1"], sub["KRT19"], s=4, alpha=0.35,
                   c=palette.get(label, "#888"), label=label)
    ax.set_xlabel("PROX1 (log-normalized)")
    ax.set_ylabel("KRT19 (log-normalized)")
    ax.set_title("PROX1 vs KRT19 in human liver epithelium")
    ax.legend(markerscale=4)
    ax.axhline(0, color="grey", linewidth=0.5)
    ax.axvline(0, color="grey", linewidth=0.5)
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Per-donor PROX1 ↔ IAD score relationship ──────────────────────────────────

def per_donor_prox1_vs_iadscore(chol: ad.AnnData, panel: pd.DataFrame) -> pd.DataFrame:
    """For each cholangiocyte donor, compute mean PROX1 and mean IAD score."""
    if "PROX1" not in chol.var_names:
        print("  PROX1 not in cholangiocyte var_names — donor table omits PROX1 column.")
        prox1_per_cell = np.zeros(chol.n_obs)
    else:
        prox1_per_cell = get_gene_expression(chol, "PROX1")

    panel_genes = _panel_genes_present(chol, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(chol, panel_genes, weights=weights)
    healthy_ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    iad_score = 1.0 - np.clip(panel_expr / max(healthy_ref, 1e-6), 0.0, 1.0)

    if "donor_id" not in chol.obs.columns:
        return pd.DataFrame({
            "cell": chol.obs_names, "prox1": prox1_per_cell,
            "iad_score": iad_score,
        })

    df = pd.DataFrame({
        "donor":     chol.obs["donor_id"].astype(str).values,
        "disease":   chol.obs["disease"].astype(str).values
                     if "disease" in chol.obs.columns
                     else "unknown",
        "prox1":     prox1_per_cell,
        "iad_score": iad_score,
    })
    out = (
        df.groupby(["donor", "disease"])
        .agg(mean_prox1=("prox1", "mean"),
             mean_iad_score=("iad_score", "mean"),
             n_cells=("prox1", "count"))
        .reset_index()
        .sort_values("mean_iad_score", ascending=False)
    )
    return out


def plot_per_donor(df: pd.DataFrame, out: Path) -> None:
    if df.empty or "mean_prox1" not in df.columns:
        return
    fig, ax = plt.subplots(figsize=(7, 5))
    diseases = sorted(df["disease"].unique())
    palette = sns.color_palette("Set2", n_colors=len(diseases))
    cmap = dict(zip(diseases, palette))
    for d in diseases:
        sub = df[df["disease"] == d]
        ax.scatter(sub["mean_prox1"], sub["mean_iad_score"],
                   s=np.clip(sub["n_cells"] / 5, 10, 200),
                   c=[cmap[d]], alpha=0.75, label=d)
    ax.set_xlabel("Mean PROX1 (log-norm) per donor — cholangiocytes")
    ax.set_ylabel("Mean IAD score per donor")
    ax.set_title("Per-donor: PROX1 vs IAD score (cholangiocyte program loss)")
    ax.legend(fontsize=8, frameon=False, loc="best")
    ax.axvline(df["mean_prox1"].median(), color="grey", linestyle="--", linewidth=0.7)
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
    ensure_dir(HEP_PROC.parent)

    print("=" * 60)
    print("PROX1 ↔ IAD link analysis")
    print("=" * 60)

    if not CHOL_H5AD.exists():
        print(f"ERROR: {CHOL_H5AD} not found", file=sys.stderr)
        sys.exit(1)

    print("\n[1/5] Loading cholangiocytes …")
    chol = sc.read_h5ad(CHOL_H5AD)
    print(f"  {chol.n_obs:,} cells × {chol.n_vars:,} genes")

    print("\n[2/5] Fetching hepatocyte reference slice …")
    if HEP_PROC.exists():
        print(f"  Reusing processed {HEP_PROC}")
        hep = sc.read_h5ad(HEP_PROC)
    else:
        hep = fetch_hepatocytes()
        hep = quick_preprocess(hep)
        hep.write_h5ad(HEP_PROC)
        print(f"  Saved → {HEP_PROC}")
    print(f"  Hepatocytes: {hep.n_obs:,} cells × {hep.n_vars:,} genes")

    print("\n[3/5] Quantifying PROX1 and KRT19 in each population …")
    rows = []
    for label, adata in [("cholangiocyte", chol), ("hepatocyte", hep)]:
        for gene in ("PROX1", "KRT19"):
            q = quantify_marker(adata, gene)
            rows.append({"cell_type": label, "gene": gene, **q})
    quant = pd.DataFrame(rows)
    print(quant.to_string(index=False))

    print("\n[4/5] Counting PROX1+ KRT19+ transition cells …")
    chol_trans = transition_cell_count(chol)
    hep_trans  = transition_cell_count(hep)
    print(f"  Cholangiocyte dataset transition cells: {chol_trans:,} / {chol.n_obs:,}")
    print(f"  Hepatocyte    dataset transition cells: {hep_trans:,} / {hep.n_obs:,}")

    print("\n  Generating PROX1 × KRT19 scatter …")
    plot_prox1_vs_krt19(chol, hep, FIGURES_DIR / "prox1_vs_krt19.png")

    print("\n[5/5] Per-donor PROX1 vs IAD score …")
    panel = load_panel()
    donor_df = per_donor_prox1_vs_iadscore(chol, panel)
    plot_per_donor(donor_df, FIGURES_DIR / "prox1_score_vs_iad_score.png")

    print("\nSaving link summary TSV …")
    quant.assign(group="population_quantification") \
         .to_csv(LINK_TSV, sep="\t", index=False)
    donor_path = RESULTS_DIR / "prox1_iad_link_per_donor.tsv"
    donor_df.to_csv(donor_path, sep="\t", index=False)
    print(f"  → {LINK_TSV}")
    print(f"  → {donor_path}")

    # ── Console-friendly conclusion ───────────────────────────────────────────
    print("\n" + "=" * 60)
    print("INTERPRETATION (numerical)")
    print("=" * 60)
    chol_prox1 = quant.query("cell_type=='cholangiocyte' and gene=='PROX1'").iloc[0]
    hep_prox1  = quant.query("cell_type=='hepatocyte'   and gene=='PROX1'").iloc[0]
    chol_krt19 = quant.query("cell_type=='cholangiocyte' and gene=='KRT19'").iloc[0]
    hep_krt19  = quant.query("cell_type=='hepatocyte'   and gene=='KRT19'").iloc[0]
    print(f"PROX1 % positive   — cholangiocyte: {chol_prox1['pct_pos']:5.2f}%   "
          f"hepatocyte: {hep_prox1['pct_pos']:5.2f}%")
    print(f"KRT19 % positive   — cholangiocyte: {chol_krt19['pct_pos']:5.2f}%   "
          f"hepatocyte: {hep_krt19['pct_pos']:5.2f}%")
    print(f"Transition cells   — cholangiocyte set: {chol_trans:,}   "
          f"hepatocyte set: {hep_trans:,}")
    print(
        "\nThe expected pattern is reciprocal: PROX1+/KRT19- in hepatocytes,\n"
        "PROX1-/KRT19+ in cholangiocytes. A small transition population\n"
        "expressing both is the candidate regenerative substrate.\n"
    )


if __name__ == "__main__":
    main()
