"""
Pseudobulk Differential Expression — Disease vs Normal Cholangiocytes
======================================================================

Pipeline stage: AGGREGATE → COMPUTE → SAVE.

Why pseudobulk
--------------
Per-cell DEG (e.g., scanpy.tl.rank_genes_groups) treats each cell as an
independent sample, which inflates significance because cells from the
same donor are not independent. The reviewer-grade approach for disease-
vs-control comparisons in scRNA-seq is to:

    1. Aggregate raw counts per donor per cell-type cluster (pseudobulk).
    2. Treat each donor as one observation.
    3. Run a standard bulk-RNA differential expression model (DESeq2)
       with study/dataset_id as a covariate to absorb batch effects.

This matches Squair et al. 2021 (Nature Comm.) and is now the default
for scRNA-seq disease comparisons.

Outputs
-------
results/pseudobulk_PSC_vs_normal.tsv
results/pseudobulk_PBC_vs_normal.tsv
results/pseudobulk_summary.tsv
figures/comparison/pseudobulk_volcano_PSC.png
figures/comparison/pseudobulk_volcano_PBC.png
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR    = Path("data")
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures") / "comparison"
LIVER_H5AD  = DATA_DIR / "liver" / "liver_cholangiocytes.h5ad"

# ── Config ────────────────────────────────────────────────────────────────────
DONOR_COL    = "donor_id"
DISEASE_COL  = "disease"
STUDY_COL    = "dataset_id"   # batch/study covariate

# Disease label normalization (Census labels are verbose)
DISEASE_NORMAL = "normal"
DISEASE_PSC    = "primary sclerosing cholangitis"
DISEASE_PBC    = "primary biliary cholangitis"

MIN_CELLS_PER_DONOR = 20      # drop donors with too few cells
MIN_DONORS_PER_GROUP = 3      # need at least this many donors in each group


# ── Aggregation ───────────────────────────────────────────────────────────────

def get_raw_counts(adata: ad.AnnData) -> np.ndarray:
    """Return a 2D array of integer raw counts. Use 'counts' layer if present."""
    if "counts" in adata.layers:
        X = adata.layers["counts"]
    else:
        # If no counts layer, X may already be raw or normalized.
        # Pancreas and liver pipelines stored counts in layers["counts"].
        X = adata.X
    if hasattr(X, "toarray"):
        X = X.toarray()
    return np.asarray(X)


def pseudobulk(
    adata: ad.AnnData,
    donor_col: str = DONOR_COL,
    disease_col: str = DISEASE_COL,
    study_col: str = STUDY_COL,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate counts per donor.

    Returns
    -------
    counts : DataFrame  (rows = donors, columns = genes)
    meta   : DataFrame  (rows = donors, columns = disease, study, n_cells)
    """
    for col in (donor_col, disease_col):
        if col not in adata.obs.columns:
            print(f"ERROR: required column '{col}' missing from .obs",
                  file=sys.stderr)
            sys.exit(1)

    print(f"  Aggregating {adata.n_obs:,} cells across donors …")
    X = get_raw_counts(adata)
    obs = adata.obs.copy()
    obs["__donor"] = obs[donor_col].astype(str)
    donors = sorted(obs["__donor"].unique())

    rows, meta_rows = [], []
    for d in donors:
        idx = obs["__donor"].values == d
        n = int(idx.sum())
        if n < MIN_CELLS_PER_DONOR:
            continue
        rows.append(X[idx].sum(axis=0).astype(np.int64))
        # disease & study assumed constant within donor; take mode
        disease_vals = obs.loc[idx, disease_col].astype(str)
        study_vals   = obs.loc[idx, study_col].astype(str) \
            if study_col in obs.columns else pd.Series(["unknown"] * n)
        meta_rows.append({
            "donor":   d,
            "disease": disease_vals.mode().iloc[0],
            "study":   study_vals.mode().iloc[0],
            "n_cells": n,
        })

    counts = pd.DataFrame(rows, index=[m["donor"] for m in meta_rows],
                          columns=adata.var_names)
    meta = pd.DataFrame(meta_rows).set_index("donor")
    print(f"  Donors retained: {len(counts):,} (min cells = {MIN_CELLS_PER_DONOR})")
    return counts, meta


# ── DESeq2 via pyDESeq2 ───────────────────────────────────────────────────────

def run_deseq(
    counts: pd.DataFrame, meta: pd.DataFrame,
    contrast_label: str, ref_label: str,
) -> pd.DataFrame:
    """
    Run pyDESeq2 with disease as the design variable and study as covariate.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds  import DeseqStats
    except ImportError:
        print("ERROR: pydeseq2 not installed. Run: pip install pydeseq2",
              file=sys.stderr)
        sys.exit(1)

    sub = meta["disease"].isin([contrast_label, ref_label])
    n_contrast = int((meta.loc[sub, "disease"] == contrast_label).sum())
    n_ref      = int((meta.loc[sub, "disease"] == ref_label).sum())
    if n_contrast < MIN_DONORS_PER_GROUP or n_ref < MIN_DONORS_PER_GROUP:
        print(f"  Skipping {contrast_label} vs {ref_label}: "
              f"n_contrast={n_contrast}, n_ref={n_ref} (need ≥ {MIN_DONORS_PER_GROUP}).")
        return pd.DataFrame()

    c = counts.loc[sub.index[sub]]
    m = meta.loc[sub.index[sub]].copy()
    m["disease"] = pd.Categorical(m["disease"], categories=[ref_label, contrast_label])
    # If only one study, drop the covariate (DESeq2 needs >1 level)
    use_study = m["study"].nunique() > 1
    design = "~ study + disease" if use_study else "~ disease"
    print(f"  pyDESeq2 design: {design}  "
          f"({n_contrast} {contrast_label} vs {n_ref} {ref_label})")

    dds = DeseqDataSet(
        counts=c.astype(int),
        metadata=m,
        design=design,
        refit_cooks=True,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(
        dds,
        contrast=("disease", contrast_label, ref_label),
        quiet=True,
    )
    stats.summary()
    res = stats.results_df.copy()
    res["gene"] = res.index
    res = res.reset_index(drop=True)
    res = res.sort_values("padj", na_position="last")
    return res


# ── Volcano plot ──────────────────────────────────────────────────────────────

def plot_volcano(res: pd.DataFrame, out: Path, title: str,
                 n_label: int = 15) -> None:
    if res.empty:
        return
    df = res.dropna(subset=["log2FoldChange", "padj"]).copy()
    df["nlog10p"] = -np.log10(df["padj"].clip(lower=1e-300))

    # Color
    sig = (df["padj"] < 0.05) & (df["log2FoldChange"].abs() >= 1)
    up   = sig & (df["log2FoldChange"] > 0)
    down = sig & (df["log2FoldChange"] < 0)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "nlog10p"],
               s=8, alpha=0.3, c="#aaaaaa")
    ax.scatter(df.loc[up, "log2FoldChange"], df.loc[up, "nlog10p"],
               s=12, alpha=0.7, c="#c0392b", label="up in disease")
    ax.scatter(df.loc[down, "log2FoldChange"], df.loc[down, "nlog10p"],
               s=12, alpha=0.7, c="#2980b9", label="down in disease")

    # Top-N labels by |log2FC| among significant
    top = df[sig].assign(absL2=df.loc[sig, "log2FoldChange"].abs()) \
                 .nlargest(n_label, "absL2")
    for _, row in top.iterrows():
        ax.text(row["log2FoldChange"], row["nlog10p"], row["gene"],
                fontsize=7, alpha=0.85)

    ax.axvline( 1, color="grey", linestyle="--", linewidth=0.6)
    ax.axvline(-1, color="grey", linestyle="--", linewidth=0.6)
    ax.axhline(-np.log10(0.05), color="grey", linestyle="--", linewidth=0.6)
    ax.set_xlabel("log₂ fold change (disease / normal)")
    ax.set_ylabel("−log₁₀(adjusted p)")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=9)
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

    print("=" * 60)
    print("Pseudobulk differential expression — cholangiocytes")
    print("=" * 60)

    if not LIVER_H5AD.exists():
        print(f"ERROR: {LIVER_H5AD} missing", file=sys.stderr); sys.exit(1)

    print("\n--- Loading liver cholangiocytes ---")
    adata = sc.read_h5ad(LIVER_H5AD)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    print("\n--- Pseudobulk aggregation ---")
    counts, meta = pseudobulk(adata)
    print(f"  Counts matrix: {counts.shape}")
    print(f"  Disease distribution:")
    print(meta["disease"].value_counts().to_string())

    print("\n--- DESeq2 PSC vs normal ---")
    psc = run_deseq(counts, meta, contrast_label=DISEASE_PSC, ref_label=DISEASE_NORMAL)
    if not psc.empty:
        psc.to_csv(RESULTS_DIR / "pseudobulk_PSC_vs_normal.tsv", sep="\t", index=False)
        plot_volcano(psc, FIGURES_DIR / "pseudobulk_volcano_PSC.png",
                     "PSC vs normal cholangiocytes")

    print("\n--- DESeq2 PBC vs normal ---")
    pbc = run_deseq(counts, meta, contrast_label=DISEASE_PBC, ref_label=DISEASE_NORMAL)
    if not pbc.empty:
        pbc.to_csv(RESULTS_DIR / "pseudobulk_PBC_vs_normal.tsv", sep="\t", index=False)
        plot_volcano(pbc, FIGURES_DIR / "pseudobulk_volcano_PBC.png",
                     "PBC vs normal cholangiocytes")

    print("\n--- Cross-disease summary ---")
    rows = []
    for label, df in [("PSC_vs_normal", psc), ("PBC_vs_normal", pbc)]:
        if df.empty:
            continue
        sig_up   = ((df["padj"] < 0.05) & (df["log2FoldChange"] >=  1)).sum()
        sig_down = ((df["padj"] < 0.05) & (df["log2FoldChange"] <= -1)).sum()
        rows.append({
            "comparison":  label,
            "n_genes_tested": len(df),
            "n_sig_up":   int(sig_up),
            "n_sig_down": int(sig_down),
        })
    summary = pd.DataFrame(rows)
    if not summary.empty:
        summary.to_csv(RESULTS_DIR / "pseudobulk_summary.tsv", sep="\t", index=False)
        print(summary.to_string(index=False))

    # Highlight overlap between disease-down genes and the IAD panel
    panel_path = RESULTS_DIR / "iad_diagnostic_panel.tsv"
    if panel_path.exists() and not psc.empty:
        panel = pd.read_csv(panel_path, sep="\t")
        psc_down = set(psc.query("padj < 0.05 and log2FoldChange <= -0.5")["gene"])
        overlap = sorted(set(panel["gene"]) & psc_down)
        print(f"\nIAD panel genes also significantly DOWN in PSC: "
              f"{len(overlap)}/{len(panel)}")
        for g in overlap:
            print(f"  - {g}")

    print("\nDone.")


if __name__ == "__main__":
    main()
