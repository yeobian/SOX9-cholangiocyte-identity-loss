"""
Shared utility functions for the IAD diagnostic-panel pipeline.

Mirrors the helpers from the Pancreas project so the analysis recipe
(load → filter → transform → compute → save) stays identical.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad


# ── Filesystem ────────────────────────────────────────────────────────────────

def ensure_dir(path) -> Path:
    """Create directory (and parents) if it doesn't exist."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


# ── QC plot ───────────────────────────────────────────────────────────────────

def plot_qc(adata: ad.AnnData, save_path: Path, title: str = "QC") -> None:
    """Save a QC violin plot (n_genes, total_counts, pct_mito)."""
    ensure_dir(save_path.parent)
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    available = [m for m in metrics if m in adata.obs.columns]
    if not available:
        return
    sc.pl.violin(adata, available, jitter=0.4, multi_panel=True, show=False)
    plt.suptitle(title, y=1.01)
    plt.savefig(save_path, bbox_inches="tight", dpi=150)
    plt.close()


# ── Gene access ───────────────────────────────────────────────────────────────

def get_gene_expression(adata: ad.AnnData, gene: str) -> np.ndarray:
    """Extract a single gene's expression vector as a 1-D float array."""
    if gene not in adata.var_names:
        raise ValueError(
            f"Gene '{gene}' not found in dataset. "
            f"First 10 var_names: {list(adata.var_names[:10])}"
        )
    idx = adata.var_names.get_loc(gene)
    X = adata.X
    if hasattr(X, "toarray"):
        return X[:, idx].toarray().flatten().astype(float)
    return np.array(X[:, idx]).flatten().astype(float)


# ── Cholangiocyte / ductal cell selection ────────────────────────────────────

CHOLANGIOCYTE_KEYWORDS = (
    "cholangiocyte",
    "bile duct",
    "biliary",
    "ductal",
    "duct epithelial",
)


def subset_to_ductal(adata: ad.AnnData, label_col: str = "cell_type") -> ad.AnnData:
    """
    Subset an AnnData to ductal / cholangiocyte cells based on cell-type labels.

    Both liver cholangiocytes and pancreas ductal cells are biliary-tree epithelium,
    so we use one keyword set for both tissues. If no matches are found, returns
    the original AnnData unchanged (useful for unlabeled clustered data).
    """
    if label_col not in adata.obs.columns:
        print(f"  [subset_to_ductal] '{label_col}' missing — returning all cells.")
        return adata

    labels = adata.obs[label_col].astype(str).str.lower()
    mask = labels.apply(
        lambda s: any(kw in s for kw in CHOLANGIOCYTE_KEYWORDS)
    )
    n_match = int(mask.sum())
    if n_match == 0:
        print(f"  [subset_to_ductal] no ductal labels matched — returning all cells.")
        return adata

    print(f"  [subset_to_ductal] keeping {n_match:,} ductal cells of {adata.n_obs:,}.")
    return adata[mask.values].copy()


# ── Multi-anchor cholangiocyte score ─────────────────────────────────────────

CHOLANGIOCYTE_ANCHORS = ("KRT19", "KRT7", "EPCAM")


def compute_anchor_score(
    adata: ad.AnnData,
    anchors: tuple[str, ...] = CHOLANGIOCYTE_ANCHORS,
) -> np.ndarray:
    """
    Per-cell cholangiocyte identity score based on multiple anchor genes.

    Uses the geometric mean of (1 + log-normalized expression) across anchor
    genes — more robust than any single gene because dropout in one gene is
    compensated by the others. Returns a 1-D float array of length n_obs.
    """
    present = [g for g in anchors if g in adata.var_names]
    missing = [g for g in anchors if g not in adata.var_names]
    if missing:
        print(f"  [anchor-score] missing anchors: {missing}; using {present}")
    if not present:
        raise ValueError(
            f"None of the anchor genes {anchors} are in the dataset."
        )

    cols = [get_gene_expression(adata, g) for g in present]
    M = np.vstack(cols)                         # n_anchors × n_cells
    # Geometric mean of (1 + expr); subtract 1 to keep zero-anchored
    log_mean = np.log1p(M).mean(axis=0)
    return np.expm1(log_mean)


# ── Co-expression ─────────────────────────────────────────────────────────────

def compute_gene_correlation(
    adata: ad.AnnData,
    target_gene: str,
    method: str = "spearman",
    min_expr_fraction: float = 0.01,
) -> pd.Series:
    """
    Spearman / Pearson correlation of every gene with `target_gene` across cells.
    Returns a Series sorted descending, with the target gene itself dropped.
    """
    from scipy.stats import rankdata

    target_expr = get_gene_expression(adata, target_gene)

    X = adata.X
    X_dense = X.toarray() if hasattr(X, "toarray") else np.array(X, dtype=float)

    expr_frac = (X_dense > 0).mean(axis=0)
    keep = expr_frac >= min_expr_fraction
    X_filt = X_dense[:, keep]
    gene_names = adata.var_names[keep]

    print(f"    Computing {method} correlations for {X_filt.shape[1]:,} genes …")

    if method == "spearman":
        t_ranked = rankdata(target_expr)
        X_ranked = np.apply_along_axis(rankdata, 0, X_filt)
        t_c = t_ranked - t_ranked.mean()
        X_c = X_ranked - X_ranked.mean(axis=0)
    else:
        t_c = target_expr - target_expr.mean()
        X_c = X_filt - X_filt.mean(axis=0)

    num = t_c @ X_c
    denom = np.sqrt((t_c ** 2).sum()) * np.sqrt((X_c ** 2).sum(axis=0))
    denom = np.where(denom == 0, 1e-10, denom)
    corr = num / denom

    result = pd.Series(corr, index=gene_names)
    result = result.drop(labels=[target_gene], errors="ignore")
    return result.sort_values(ascending=False)
