"""
Preprocess human liver scRNA-seq from CELLxGENE Census.

Pipeline stage: FILTER → TRANSFORM.

1. Load raw AnnData written by download_data.py
2. Validate gene symbols + presence of cholangiocyte markers
3. QC filtering (min/max genes, mitochondrial content)
4. Normalize → log1p → HVGs → scale → PCA → UMAP → Leiden
5. Save full liver to data/liver/liver_processed.h5ad
6. Save cholangiocyte-only subset to data/liver/liver_cholangiocytes.h5ad
   (this is the file the downstream IAD analysis runs on)

Usage
-----
    python src/preprocess_liver.py
"""

import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, plot_qc, subset_to_ductal

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR     = Path("data")
LIVER_RAW    = DATA_DIR / "liver" / "liver_raw.h5ad"
LIVER_FULL   = DATA_DIR / "liver" / "liver_processed.h5ad"
LIVER_CHOL   = DATA_DIR / "liver" / "liver_cholangiocytes.h5ad"
FIGURES_DIR  = Path("figures") / "liver"

# ── QC thresholds (looser MAX_GENES because liver has hepatocytes) ────────────
MIN_GENES   = 200
MAX_GENES   = 8_000
MAX_PCT_MT  = 25
MIN_CELLS   = 3

# Cholangiocyte sentinel — refuse to proceed if this gene is missing.
SENTINEL_GENE = "KRT19"


# ── Load ──────────────────────────────────────────────────────────────────────

def load_raw() -> ad.AnnData:
    if not LIVER_RAW.exists():
        print(f"ERROR: raw liver file not found at {LIVER_RAW}", file=sys.stderr)
        print("Run:  python src/download_data.py liver", file=sys.stderr)
        sys.exit(1)
    print(f"Loading {LIVER_RAW} …")
    adata = sc.read_h5ad(LIVER_RAW)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


def validate_markers(adata: ad.AnnData) -> None:
    if SENTINEL_GENE not in adata.var_names:
        print(
            f"ERROR: sentinel gene '{SENTINEL_GENE}' not in var_names.\n"
            f"  First 10: {list(adata.var_names[:10])}",
            file=sys.stderr,
        )
        sys.exit(1)
    print(f"  Sentinel '{SENTINEL_GENE}' present.")


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

def save_outputs(full: ad.AnnData, chol: ad.AnnData) -> None:
    ensure_dir(LIVER_FULL.parent)
    print(f"\nSaving full liver  → {LIVER_FULL}")
    full.write_h5ad(LIVER_FULL)
    print(f"Saving cholangiocytes → {LIVER_CHOL}  ({chol.n_obs:,} cells)")
    chol.write_h5ad(LIVER_CHOL)


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    ensure_dir(FIGURES_DIR)
    sc.settings.verbosity = 2
    sc.settings.set_figure_params(dpi=100, facecolor="white")

    print("=== Liver Preprocessing (CELLxGENE Census) ===")

    adata = load_raw()
    validate_markers(adata)
    adata = run_qc(adata)

    print("\nSaving QC figure …")
    plot_qc(adata, FIGURES_DIR / "qc_violin.png", title="Liver QC")

    adata = preprocess(adata)

    print("\n--- Subsetting to cholangiocytes ---")
    chol = subset_to_ductal(adata, label_col="cell_type")

    save_outputs(adata, chol)
    print("\nDone!\nNext: python src/iad_analysis.py")


if __name__ == "__main__":
    main()
