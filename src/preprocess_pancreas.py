"""
Preprocess human pancreas scRNA-seq from CELLxGENE Census.

Pipeline stage: FILTER → TRANSFORM.

Mirrors preprocess_liver.py but for pancreas. Produces both the full
pancreas AnnData and a ductal-cell-only subset that serves as the
cross-tissue cholangiocyte reference.

Outputs
-------
data/pancreas/pancreas_processed.h5ad
data/pancreas/pancreas_ductal.h5ad

Usage
-----
    python src/preprocess_pancreas.py
"""

import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, plot_qc, subset_to_ductal

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR        = Path("data")
PANCREAS_RAW    = DATA_DIR / "pancreas" / "pancreas_raw.h5ad"
PANCREAS_FULL   = DATA_DIR / "pancreas" / "pancreas_processed.h5ad"
PANCREAS_DUCT   = DATA_DIR / "pancreas" / "pancreas_ductal.h5ad"
FIGURES_DIR     = Path("figures") / "pancreas"

# ── QC thresholds ─────────────────────────────────────────────────────────────
MIN_GENES   = 200
MAX_GENES   = 8_000
MAX_PCT_MT  = 25
MIN_CELLS   = 3


# ── Load ──────────────────────────────────────────────────────────────────────

def load_raw() -> ad.AnnData:
    if not PANCREAS_RAW.exists():
        print(f"ERROR: raw pancreas file not found at {PANCREAS_RAW}", file=sys.stderr)
        print("Run:  python src/download_data.py pancreas", file=sys.stderr)
        sys.exit(1)
    print(f"Loading {PANCREAS_RAW} …")
    adata = sc.read_h5ad(PANCREAS_RAW)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


# ── QC ────────────────────────────────────────────────────────────────────────

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


# ── Transform ─────────────────────────────────────────────────────────────────

def preprocess(adata: ad.AnnData) -> ad.AnnData:
    print("\n--- Preprocessing ---")
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(
        adata, min_mean=0.0125, max_mean=3, min_disp=0.5, subset=True,
    )
    print(f"  Highly variable genes: {adata.n_vars:,}")
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
    sc.tl.umap(adata)
    sc.tl.leiden(
        adata, resolution=0.5, key_added="leiden",
        flavor="igraph", n_iterations=2, directed=False,
    )
    return adata


# ── Save ──────────────────────────────────────────────────────────────────────

def save_outputs(full: ad.AnnData, duct: ad.AnnData) -> None:
    ensure_dir(PANCREAS_FULL.parent)
    print(f"\nSaving full pancreas → {PANCREAS_FULL}")
    full.write_h5ad(PANCREAS_FULL)
    print(f"Saving ductal subset → {PANCREAS_DUCT}  ({duct.n_obs:,} cells)")
    duct.write_h5ad(PANCREAS_DUCT)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ensure_dir(FIGURES_DIR)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor="white")

    print("=== Pancreas Preprocessing (CELLxGENE Census) ===")
    adata = load_raw()
    adata = run_qc(adata)

    print("\nSaving QC figure …")
    plot_qc(adata, FIGURES_DIR / "qc_violin.png", title="Pancreas QC")

    adata = preprocess(adata)

    print("\n--- Subsetting to ductal cells ---")
    duct = subset_to_ductal(adata, label_col="cell_type")

    save_outputs(adata, duct)
    print("\nDone!\nNext: python src/iad_analysis.py")


if __name__ == "__main__":
    main()
