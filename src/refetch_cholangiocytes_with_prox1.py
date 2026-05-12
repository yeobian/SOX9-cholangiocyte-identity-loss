"""
Re-build the cholangiocyte AnnData while *forcing PROX1* (and a curated panel)
to be retained through HVG selection.

Why
---
The original preprocess_liver.py runs sc.pp.highly_variable_genes(subset=True),
which drops genes with low variance — including PROX1 in cholangiocytes,
where it is expressed in <1% of cells. That makes downstream PROX1 analysis
on cholangiocytes impossible. This script re-runs preprocessing with
subset=False and a manual gene-keep list that *unconditionally* retains the
PROX1 panel.

Pipeline (load → filter → transform → save)
-------------------------------------------
1. Load data/liver/liver_raw.h5ad (already pulled from CELLxGENE Census)
2. Standard QC + log1p
3. HVG with subset=False
4. Keep = HVG ∪ PROX1_PANEL
5. Subset to cholangiocytes (cell_type contains 'cholangiocyte' or 'bile duct')
6. Save to data/liver/liver_cholangiocytes_with_prox1.h5ad

Usage
-----
    python src/refetch_cholangiocytes_with_prox1.py
Prerequisites: data/liver/liver_raw.h5ad must exist.
"""

import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, subset_to_ductal

# ── Config ────────────────────────────────────────────────────────────────────
DATA_DIR = Path("data")
LIVER_RAW = DATA_DIR / "liver" / "liver_raw.h5ad"
OUT_PATH  = DATA_DIR / "liver" / "liver_cholangiocytes_with_prox1.h5ad"

MIN_GENES   = 200
MAX_GENES   = 8_000
MAX_PCT_MT  = 25
MIN_CELLS   = 3

# Genes that MUST be kept regardless of HVG status
PROX1_PANEL = {
    "PROX1",
    "KRT19", "KRT7", "SOX9", "HNF1B", "EPCAM", "CFTR",
    "ALB", "HNF4A", "HNF1A", "CEBPA", "FOXA2", "TTR",
    "JAG1", "NOTCH2", "NOTCH3", "HES1",
    "SPP1", "TACSTD2", "LGR5", "PROM1",
    "CTGF", "CYR61", "ANKRD1",                    # YAP/TAZ targets
    "ABCB11", "ABCB4", "SLC10A1",                 # bile transport
}


def load_raw() -> ad.AnnData:
    if not LIVER_RAW.exists():
        print(f"ERROR: {LIVER_RAW} not found. Run download_data.py first.",
              file=sys.stderr)
        sys.exit(1)
    print(f"Loading {LIVER_RAW} …")
    adata = sc.read_h5ad(LIVER_RAW)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


def run_qc(adata: ad.AnnData) -> ad.AnnData:
    print("\n--- QC ---")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True,
    )
    print(f"  Before: {adata.n_obs:,} cells | {adata.n_vars:,} genes")
    sc.pp.filter_cells(adata, min_genes=MIN_GENES)
    adata = adata[adata.obs.n_genes_by_counts <= MAX_GENES].copy()
    adata = adata[adata.obs.pct_counts_mt <= MAX_PCT_MT].copy()
    sc.pp.filter_genes(adata, min_cells=MIN_CELLS)
    print(f"  After:  {adata.n_obs:,} cells | {adata.n_vars:,} genes")
    return adata


def normalize_with_panel(adata: ad.AnnData) -> ad.AnnData:
    print("\n--- Preprocessing (PROX1-keeping) ---")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG WITHOUT subsetting — we'll add the panel back manually
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=False,
    )
    n_hvg = int(adata.var["highly_variable"].sum())
    n_panel_present = int(adata.var_names.isin(PROX1_PANEL).sum())
    n_panel_in_hvg  = int(
        (adata.var["highly_variable"] & adata.var_names.isin(PROX1_PANEL)).sum()
    )
    print(f"  HVGs flagged:        {n_hvg:,}")
    print(f"  Panel genes present: {n_panel_present:,}/{len(PROX1_PANEL)}")
    print(f"  Panel ∩ HVG:         {n_panel_in_hvg:,}  "
          f"(rest will be force-kept)")

    keep = adata.var["highly_variable"] | adata.var_names.isin(PROX1_PANEL)
    adata = adata[:, keep].copy()
    print(f"  Kept after union: {adata.n_vars:,} genes")
    print(f"  PROX1 in kept genes? {'PROX1' in adata.var_names}")

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(
        adata, resolution=0.5, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )
    return adata


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(OUT_PATH.parent)

    print("=" * 60)
    print("Re-fetch cholangiocytes WITH PROX1 retained through HVG")
    print("=" * 60)

    adata = load_raw()
    adata = run_qc(adata)
    adata = normalize_with_panel(adata)

    print("\n--- Subsetting to cholangiocytes ---")
    chol = subset_to_ductal(adata, label_col="cell_type")
    print(f"  PROX1 in cholangiocyte var_names? {'PROX1' in chol.var_names}")
    if "PROX1" in chol.var_names:
        from utils import get_gene_expression
        prox1 = get_gene_expression(chol, "PROX1")
        pct  = (prox1 > 0).mean() * 100
        mean = prox1[prox1 > 0].mean() if pct > 0 else 0.0
        print(f"  PROX1 % positive in cholangiocytes: {pct:.2f}%  "
              f"mean(positive cells)={mean:.3f}")

    print(f"\nSaving → {OUT_PATH} ({chol.n_obs:,} cells × {chol.n_vars:,} genes)")
    chol.write_h5ad(OUT_PATH)
    print("\nDone. Now point prox1_iad_link.py / iad_analysis.py at this file:")
    print(f"  CHOL_PATH = '{OUT_PATH.as_posix()}'")


if __name__ == "__main__":
    main()
