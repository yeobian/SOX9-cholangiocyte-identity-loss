"""
IAD Risk-Classification Scoring
================================

Pipeline stage: COMPUTE.

Turn the diagnostic gene panel into a per-sample 0-1 ductal-loss index.

Concept
-------
The panel genes (results/iad_diagnostic_panel.tsv) are markers whose
expression tracks KRT19 in healthy liver cholangiocytes. In an IAD biopsy
those cells are gone, so panel-gene expression collapses. We turn that into
a continuous score:

    panel_expr_per_cell  = mean(log-normalized panel-gene expression)
    healthy_max          = 95th percentile across the healthy reference
    score_per_cell       = 1 - clip(panel_expr / healthy_max, 0, 1)

    score = 1.0  →  cholangiocyte program absent (suspicious for IAD)
    score = 0.0  →  cholangiocyte program intact (healthy biliary epithelium)

Per-donor scores are then the mean across that donor's cells.

Usage
-----
    # CLI sanity check on existing liver + pancreas h5ad files
    python src/iad_score.py

    # Programmatic
    from iad_score import score_anndata
    s = score_anndata(adata)              # returns pd.Series, length n_obs
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression

# ── Paths ─────────────────────────────────────────────────────────────────────
RESULTS_DIR  = Path("results")
DATA_DIR     = Path("data")
PANEL_PATH   = RESULTS_DIR / "iad_diagnostic_panel.tsv"
LIVER_H5AD   = DATA_DIR / "liver"    / "liver_cholangiocytes.h5ad"
PANCREAS_H5AD = DATA_DIR / "pancreas" / "pancreas_ductal.h5ad"

DONOR_COL = "donor_id"   # CELLxGENE Census standard column

# Genes used as healthy reference percentile (top 5% panel-expressing cells in liver).
REFERENCE_PERCENTILE = 95.0


# ── Load panel ────────────────────────────────────────────────────────────────

def load_panel(path: Path = PANEL_PATH) -> pd.DataFrame:
    if not path.exists():
        print(f"ERROR: panel TSV not found at {path}", file=sys.stderr)
        print("Run:  python src/iad_analysis.py", file=sys.stderr)
        sys.exit(1)
    df = pd.read_csv(path, sep="\t")
    required = {"gene", "rho_liver", "rho_pancreas", "specificity"}
    missing  = required - set(df.columns)
    if missing:
        print(f"ERROR: panel missing columns: {missing}", file=sys.stderr)
        sys.exit(1)
    print(f"  Loaded panel: {len(df)} genes from {path}")
    return df


# ── Compute panel expression per cell ─────────────────────────────────────────

def _panel_genes_present(adata: ad.AnnData, panel: pd.DataFrame) -> list[str]:
    present = [g for g in panel["gene"] if g in adata.var_names]
    missing = [g for g in panel["gene"] if g not in adata.var_names]
    if missing:
        print(f"  Panel genes missing from this dataset ({len(missing)}): {missing}")
    print(f"  Panel genes used: {len(present)}")
    return present


def _panel_expression_per_cell(
    adata: ad.AnnData, panel_genes: list[str], weights: pd.Series | None = None,
) -> np.ndarray:
    """Weighted mean expression across panel genes, per cell."""
    if not panel_genes:
        return np.zeros(adata.n_obs)
    cols = []
    for g in panel_genes:
        cols.append(get_gene_expression(adata, g))
    M = np.vstack(cols)                     # shape: n_genes × n_cells
    if weights is None:
        return M.mean(axis=0)
    w = weights.reindex(panel_genes).fillna(0).values
    if w.sum() == 0:
        return M.mean(axis=0)
    return (M * w[:, None]).sum(axis=0) / w.sum()


# ── Scoring ───────────────────────────────────────────────────────────────────

def score_anndata(
    adata: ad.AnnData,
    panel: pd.DataFrame | None = None,
    healthy_reference: float | None = None,
) -> pd.Series:
    """
    Per-cell IAD ductal-loss score in [0, 1].

    Parameters
    ----------
    adata
        Log-normalized AnnData. Must contain at least some panel genes.
    panel
        Panel DataFrame; if None, loaded from results/iad_diagnostic_panel.tsv.
    healthy_reference
        The expression value mapped to score = 0. If None, computed as the
        REFERENCE_PERCENTILE-th percentile of panel expression in this AnnData
        (i.e. the AnnData itself is treated as the healthy reference).
    """
    panel = panel if panel is not None else load_panel()
    genes = _panel_genes_present(adata, panel)
    if not genes:
        raise ValueError("None of the panel genes are present in this dataset.")

    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(adata, genes, weights=weights)

    if healthy_reference is None:
        healthy_reference = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
        print(f"  Self-calibrated healthy reference (p{REFERENCE_PERCENTILE:g}): "
              f"{healthy_reference:.4f}")

    if healthy_reference <= 0:
        print("  WARNING: healthy reference is 0 — all scores will be 1.")
        healthy_reference = 1e-6

    relative = np.clip(panel_expr / healthy_reference, 0.0, 1.0)
    score = 1.0 - relative
    return pd.Series(score, index=adata.obs_names, name="iad_score")


def aggregate_by_donor(
    adata: ad.AnnData, score: pd.Series, donor_col: str = DONOR_COL,
) -> pd.DataFrame:
    if donor_col not in adata.obs.columns:
        print(f"  [donor agg] '{donor_col}' missing — returning per-cell only.")
        return pd.DataFrame({"cell": score.index, "iad_score": score.values})

    df = pd.DataFrame({
        "donor":     adata.obs[donor_col].astype(str).values,
        "iad_score": score.values,
    })
    out = df.groupby("donor")["iad_score"].agg(["mean", "median", "std", "count"])
    out = out.rename(columns={
        "mean": "mean_score",
        "median": "median_score",
        "std": "sd_score",
        "count": "n_cells",
    }).reset_index().sort_values("mean_score", ascending=False)
    return out


# ── CLI sanity check ──────────────────────────────────────────────────────────

def _sanity_check_one(label: str, h5ad_path: Path, panel: pd.DataFrame,
                      healthy_ref: float | None) -> tuple[pd.Series, pd.DataFrame]:
    print(f"\n--- {label} ---")
    if not h5ad_path.exists():
        print(f"  Skipped: {h5ad_path} not found.")
        return pd.Series(dtype=float), pd.DataFrame()
    adata = sc.read_h5ad(h5ad_path)
    print(f"  Loaded {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    score = score_anndata(adata, panel=panel, healthy_reference=healthy_ref)
    print(f"  Cell-level score: mean={score.mean():.3f}  median={score.median():.3f}  "
          f"min={score.min():.3f}  max={score.max():.3f}")
    donor_df = aggregate_by_donor(adata, score)
    if "mean_score" in donor_df.columns:
        print(f"  Donor-level score (top 5):")
        for _, row in donor_df.head(5).iterrows():
            print(f"    {row['donor']:<30} mean={row['mean_score']:.3f}  n={int(row['n_cells'])}")
    return score, donor_df


def main() -> None:
    sc.settings.verbosity = 1
    print("=" * 60)
    print("IAD Risk-Classification Sanity Check")
    print("=" * 60)

    panel = load_panel()
    ensure_dir(RESULTS_DIR)

    # Calibrate the healthy reference on the LIVER dataset (true healthy biliary
    # epithelium), then apply that same reference to pancreas — pancreas should
    # score ~1.0 because it lacks the liver-cholangiocyte program.
    print("\n[1/2] Calibrating healthy reference on liver cholangiocytes …")
    if not LIVER_H5AD.exists():
        print(f"ERROR: {LIVER_H5AD} not found", file=sys.stderr)
        sys.exit(1)
    liver = sc.read_h5ad(LIVER_H5AD)
    print(f"  Liver: {liver.n_obs:,} cells × {liver.n_vars:,} genes")
    liver_genes = _panel_genes_present(liver, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    liver_panel_expr = _panel_expression_per_cell(liver, liver_genes, weights=weights)
    healthy_ref = float(np.percentile(liver_panel_expr, REFERENCE_PERCENTILE))
    print(f"  Healthy reference (liver p{REFERENCE_PERCENTILE:g}): {healthy_ref:.4f}")

    print("\n[2/2] Scoring liver and pancreas with the same reference …")
    liver_score, liver_donor = _sanity_check_one(
        "Liver cholangiocytes (expect LOW scores)",
        LIVER_H5AD, panel, healthy_ref,
    )
    pancreas_score, pancreas_donor = _sanity_check_one(
        "Pancreas ductal (expect HIGH scores — lacks liver program)",
        PANCREAS_H5AD, panel, healthy_ref,
    )

    # Save
    print("\n--- Saving ---")
    if not liver_donor.empty:
        liver_donor.to_csv(RESULTS_DIR / "iad_score_liver_donors.tsv", sep="\t", index=False)
        print(f"  → {RESULTS_DIR / 'iad_score_liver_donors.tsv'}")
    if not pancreas_donor.empty:
        pancreas_donor.to_csv(RESULTS_DIR / "iad_score_pancreas_donors.tsv", sep="\t", index=False)
        print(f"  → {RESULTS_DIR / 'iad_score_pancreas_donors.tsv'}")

    print("\nDone.")
    print("Healthy liver should be near 0; pancreas should be much higher.")
    print("Run:  python src/validate_score.py   for the full validation report.")


if __name__ == "__main__":
    main()
