"""
Step 1 — PROX1 spread vs concentration across CELLxGENE Census datasets.

Context
-------
After the HVG-keeping refetch, PROX1 is detected in 45.1% of cholangiocytes.
That is much higher than the literature would suggest. Before reading this as
biology we need to know: is the signal spread evenly across the source
datasets pooled by Census, or is it concentrated in one dataset that may have
ambiguous cholangiocyte annotations (e.g. hybrid hepatocyte–cholangiocyte
cells labeled as cholangiocytes)?

Pipeline
--------
load → group → compute → save → plot

Outputs
-------
results/prox1_dataset_spread.tsv
figures/prox1/prox1_pct_pos_by_dataset.png

Usage
-----
    python src/prox1_dataset_spread.py
"""

import sys
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression

CHOL_PATH    = Path("data") / "liver" / "liver_cholangiocytes.h5ad"
RESULTS_TSV  = Path("results") / "prox1_dataset_spread.tsv"
FIGURE_PATH  = Path("figures") / "prox1" / "prox1_pct_pos_by_dataset.png"

GENE = "PROX1"


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_TSV.parent)
    ensure_dir(FIGURE_PATH.parent)

    print("=" * 60)
    print("PROX1 spread across datasets — cholangiocytes")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    print(f"\nLoading {CHOL_PATH} …")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if GENE not in adata.var_names:
        print(f"ERROR: {GENE} not in var_names — refetch with PROX1 first.",
              file=sys.stderr); sys.exit(1)
    if "dataset_id" not in adata.obs.columns:
        print("ERROR: no dataset_id in obs.", file=sys.stderr); sys.exit(1)

    # PROX1 expression per cell (log1p-normalized values)
    prox1 = get_gene_expression(adata, GENE)
    df = pd.DataFrame({
        "dataset_id": adata.obs["dataset_id"].astype(str).values,
        "prox1":      prox1,
    })

    # Per-dataset summary
    grouped = df.groupby("dataset_id")
    summary = pd.DataFrame({
        "n_cholangiocytes": grouped.size(),
        "n_prox1_pos":      grouped["prox1"].apply(lambda s: int((s > 0).sum())),
        "pct_prox1_pos":    grouped["prox1"].apply(
            lambda s: float((s > 0).mean() * 100)),
        "median_prox1_all":      grouped["prox1"].median(),
        "median_prox1_positive": grouped["prox1"].apply(
            lambda s: float(s[s > 0].median()) if (s > 0).any() else 0.0),
    })

    # Share of all PROX1+ cells contributed by each dataset (concentration check)
    total_pos = int((df["prox1"] > 0).sum())
    summary["share_of_prox1_pos_pool"] = (
        summary["n_prox1_pos"] / max(total_pos, 1) * 100
    )
    summary = summary.sort_values("pct_prox1_pos", ascending=False)
    print(f"\nTotal cholangiocytes: {len(df):,}  "
          f"PROX1+ overall: {total_pos:,} ({total_pos/len(df)*100:.1f}%)")
    print(f"Datasets contributing cholangiocytes: {len(summary)}")
    print("\n--- Per-dataset summary ---")
    pretty = summary.copy()
    pretty["dataset_id"] = pretty.index.astype(str).str[:20] + "…"
    pretty = pretty[[
        "n_cholangiocytes", "n_prox1_pos", "pct_prox1_pos",
        "median_prox1_all", "median_prox1_positive",
        "share_of_prox1_pos_pool"
    ]]
    print(pretty.to_string(float_format="%.3f"))

    # Save TSV (use full dataset_id, not truncated)
    summary.to_csv(RESULTS_TSV, sep="\t")
    print(f"\n  → {RESULTS_TSV}")

    # Bar plot
    n = len(summary)
    fig, ax = plt.subplots(figsize=(max(7, 0.35 * n + 3), 5))
    labels = [str(x)[:14] + ("…" if len(str(x)) > 14 else "")
              for x in summary.index]
    bars = ax.bar(range(n), summary["pct_prox1_pos"].values,
                  color="#8e44ad", edgecolor="white")
    # Annotate cell counts above bars
    for i, (pct, n_cells) in enumerate(zip(summary["pct_prox1_pos"],
                                           summary["n_cholangiocytes"])):
        ax.text(i, pct + 1, f"n={int(n_cells)}",
                ha="center", va="bottom", fontsize=8, rotation=0)
    overall = total_pos / len(df) * 100
    ax.axhline(overall, color="black", linestyle="--", linewidth=1,
               label=f"Pool average = {overall:.1f}%")
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel(f"% of cholangiocytes expressing {GENE} (>0)")
    ax.set_title(f"{GENE} positivity per source dataset")
    ax.legend(loc="upper right", fontsize=9, frameon=False)
    ax.set_ylim(0, max(100, summary["pct_prox1_pos"].max() * 1.1 + 10))
    plt.tight_layout()
    plt.savefig(FIGURE_PATH, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {FIGURE_PATH}")

    # Concentration metrics
    shares = summary["share_of_prox1_pos_pool"].sort_values(ascending=False)
    top1 = float(shares.iloc[0]) if len(shares) else 0.0
    top3 = float(shares.head(3).sum()) if len(shares) >= 3 else float(shares.sum())
    pct_within_5pp = float(
        (summary["pct_prox1_pos"].between(overall - 5, overall + 5)).mean() * 100
    )

    print("\n--- Concentration check ---")
    print(f"  Top dataset contributes {top1:.1f}% of all PROX1+ cells.")
    print(f"  Top-3 datasets contribute {top3:.1f}%.")
    print(f"  {pct_within_5pp:.0f}% of datasets have PROX1+ rate within "
          f"±5 percentage points of the {overall:.1f}% pool average.")
    print()
    if top1 > 60 and len(summary) > 3:
        print("  VERDICT: CONCENTRATED. The 45% rate is dominated by a single "
              "dataset — possible annotation / hybrid-cell artifact. "
              "Inspect cells from that dataset before treating PROX1 in "
              "cholangiocytes as biology.")
    elif top3 > 80 and len(summary) > 5:
        print("  VERDICT: PARTIALLY CONCENTRATED. A few datasets account for "
              "most PROX1+ cells. Worth checking those datasets' "
              "annotation provenance.")
    else:
        print("  VERDICT: SPREAD. PROX1 positivity is broadly distributed "
              "across source datasets, supporting a biological interpretation "
              "rather than an annotation artifact.")


if __name__ == "__main__":
    main()
