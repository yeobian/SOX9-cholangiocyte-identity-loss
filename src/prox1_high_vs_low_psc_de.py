"""
Differential expression: PROX1-high vs PROX1-low cholangiocytes within PSC.

Goal
----
Characterise the PROX1-high cholangiocyte sub-state molecularly. Earlier
analyses showed PROX1-high cells are SOX9-low / biliary-identity-low —
consistent with a failed-plasticity intermediate state. This DE test asks
*what else* moves with PROX1 status, so the sub-state can be named beyond
just "PROX1+ SOX9- KRT19-".

Method
------
Donor-controlled DE. The PSC subset has 7/8 donors from a single dataset,
so any naive DE will be polluted by donor effects. We use a per-donor
pseudobulk approach:

  1. Subset to PSC cholangiocytes.
  2. Define PROX1-high = top quartile of PROX1 expression (matches the
     threshold used in prox1_abundance_test and prox1_sox9_association).
  3. For each donor with ≥ 10 PROX1-high cells AND ≥ 10 PROX1-low cells:
       compute mean gene expression in the high group minus mean in the
       low group → per-donor log-fold-change (in log1p space).
  4. Aggregate across donors: report
        - mean logFC across donors
        - fraction of donors with logFC in the same direction
        - one-sample t-test of per-donor logFCs against zero
       This is the donor-controlled DE.
  5. Sanity-check with uncontrolled Wilcoxon (scanpy.tl.rank_genes_groups)
     to compare magnitudes.

The donor-controlled test treats each donor as an independent replicate,
so a gene needs to move in the same direction in *most* donors to be
called.

Outputs
-------
results/prox1_high_vs_low_psc_de.tsv             — full per-gene table
results/prox1_high_vs_low_psc_de_top.tsv         — top 30 up + top 30 down
results/prox1_high_vs_low_psc_per_donor.tsv      — per-donor cell counts
figures/prox1/prox1_high_low_de_volcano.png
figures/prox1/prox1_high_low_de_top_genes.png

Usage
-----
    python src/prox1_high_vs_low_psc_de.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.stats import ttest_1samp, ranksums

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression
from iad_score import LIVER_H5AD as DEFAULT_CHOL

CHOL_PATH    = DEFAULT_CHOL
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "prox1"

PSC_LABEL = "primary sclerosing cholangitis"
PROX1 = "PROX1"

MIN_CELLS_PER_GROUP_PER_DONOR = 10   # min cells in high AND in low per donor
TOP_N = 30

# Gene annotation lookup — informal, for quick interpretation only
ANNOT = {
    # Hepatocyte identity
    "ALB":"hep","HNF4A":"hep","HNF1A":"hep","CEBPA":"hep","FOXA1":"hep","FOXA2":"hep",
    "TTR":"hep","APOB":"hep","APOA1":"hep","APOA2":"hep","APOC3":"hep","APOE":"hep",
    "FABP1":"hep","SERPINA1":"hep","CYP3A4":"hep","CYP2E1":"hep","CYP2A6":"hep",
    "CYP2B6":"hep","G6PC":"hep","PCK1":"hep","ASS1":"hep","HAL":"hep","GLUL":"hep",
    "AFP":"hep","ADH1B":"hep","ADH4":"hep","RBP4":"hep","SERPINA3":"hep","ORM1":"hep",
    "APCS":"hep","HP":"hep","FGA":"hep","FGB":"hep","FGG":"hep","C3":"hep","CRP":"hep",
    "AHSG":"hep","TF":"hep","ITIH3":"hep","ITIH4":"hep","GC":"hep","AGT":"hep",
    # Biliary identity
    "KRT19":"bil","KRT7":"bil","SOX9":"bil","HNF1B":"bil","EPCAM":"bil","CFTR":"bil",
    "CLDN4":"bil","MUC1":"bil","ANXA4":"bil","TFF3":"bil","AQP1":"bil","FXYD2":"bil",
    "PKHD1":"bil","SLC4A2":"bil","UMOD":"bil","GATA6":"bil","ONECUT1":"bil","FOXJ1":"bil",
    # Plasticity / regeneration / ductular reaction
    "SPP1":"plast","TACSTD2":"plast","LGR5":"plast","PROM1":"plast","NCAM1":"plast",
    "YAP1":"plast","WWTR1":"plast","CTGF":"plast","CYR61":"plast","ANKRD1":"plast",
    "JAG1":"plast","NOTCH2":"plast","NOTCH3":"plast","HES1":"plast","HEY1":"plast",
    # Stress / IEG
    "JUN":"stress","FOS":"stress","ATF3":"stress","EGR1":"stress","DUSP1":"stress",
    "JUNB":"stress","JUND":"stress","FOSB":"stress","ZFP36":"stress","HSPA1A":"stress",
    "HSPA1B":"stress","DNAJB1":"stress","NR4A1":"stress","NR4A2":"stress","BTG2":"stress",
    # Proliferation
    "MKI67":"prolif","TOP2A":"prolif","CCNA2":"prolif","CCNB1":"prolif","CDK1":"prolif",
    "FOXM1":"prolif","MYC":"prolif","AURKB":"prolif","BIRC5":"prolif",
}


def slice_psc(adata):
    if "disease" not in adata.obs.columns:
        print("ERROR: no disease column", file=sys.stderr); sys.exit(1)
    mask = adata.obs["disease"].astype(str) == PSC_LABEL
    return adata[mask.values].copy()


def split_high_low(adata, strict_thr):
    prox1 = get_gene_expression(adata, PROX1)
    high = prox1 > strict_thr
    adata.obs["prox1_group"] = pd.Categorical(
        np.where(high, "high", "low"), categories=["high", "low"],
    )
    return adata


def per_donor_pseudobulk(adata):
    """Per donor × group: mean expression of each gene (log1p-space).
    Returns:
       means_long : DataFrame with columns donor, group, gene, mean
       per_donor_counts : DataFrame donor → n_high, n_low
    """
    obs = adata.obs[["donor_id", "prox1_group"]].copy()
    obs.columns = ["donor", "group"]

    X = adata.X
    if not sparse.issparse(X):
        X = sparse.csr_matrix(X)

    genes = adata.var_names.values
    rows = []
    counts = []
    for donor, dgrp in obs.groupby("donor"):
        n_high = int((dgrp["group"] == "high").sum())
        n_low  = int((dgrp["group"] == "low").sum())
        counts.append({"donor": donor, "n_high": n_high, "n_low": n_low})
        if n_high < MIN_CELLS_PER_GROUP_PER_DONOR or \
           n_low  < MIN_CELLS_PER_GROUP_PER_DONOR:
            continue
        for grp_name, grp_idx in [
            ("high", dgrp.index[dgrp["group"] == "high"]),
            ("low",  dgrp.index[dgrp["group"] == "low"]),
        ]:
            pos = adata.obs_names.get_indexer(grp_idx)
            sub = X[pos]
            mean = np.asarray(sub.mean(axis=0)).ravel()
            for g, m in zip(genes, mean):
                rows.append({"donor": donor, "group": grp_name,
                             "gene": g, "mean": float(m)})
    return pd.DataFrame(rows), pd.DataFrame(counts)


def donor_controlled_de(means_long: pd.DataFrame) -> pd.DataFrame:
    """For each gene, compute per-donor logFC (mean_high - mean_low),
    then aggregate across donors."""
    if means_long.empty:
        return pd.DataFrame()
    wide = (
        means_long.pivot_table(index=["donor", "gene"], columns="group",
                                values="mean", aggfunc="first")
                  .reset_index()
    )
    wide["logFC"] = wide["high"] - wide["low"]

    out_rows = []
    for gene, sub in wide.groupby("gene"):
        lfcs = sub["logFC"].dropna().values
        if len(lfcs) < 3:
            continue
        # one-sample t-test against zero
        if len(set(lfcs)) > 1:
            t, p = ttest_1samp(lfcs, 0.0)
        else:
            t, p = 0.0, 1.0
        out_rows.append({
            "gene":              gene,
            "n_donors":          int(len(lfcs)),
            "mean_logFC":        float(np.mean(lfcs)),
            "median_logFC":      float(np.median(lfcs)),
            "frac_pos_donors":   float((lfcs > 0).mean()),
            "frac_neg_donors":   float((lfcs < 0).mean()),
            "t_stat":            float(t),
            "p_one_sample_t":    float(p),
        })
    df = pd.DataFrame(out_rows)
    # BH FDR
    df = df.sort_values("p_one_sample_t")
    df["rank"] = np.arange(1, len(df) + 1)
    df["q_BH"] = (df["p_one_sample_t"] * len(df) / df["rank"]).clip(upper=1.0)
    return df


def uncontrolled_wilcoxon(adata) -> pd.DataFrame:
    """scanpy.tl.rank_genes_groups Wilcoxon — sanity check."""
    sc.tl.rank_genes_groups(
        adata, groupby="prox1_group", groups=["high"], reference="low",
        method="wilcoxon", n_genes=adata.n_vars,
    )
    r = adata.uns["rank_genes_groups"]
    genes = pd.DataFrame(r["names"])["high"].values
    scores = pd.DataFrame(r["scores"])["high"].values
    logfc = pd.DataFrame(r["logfoldchanges"])["high"].values
    pvals = pd.DataFrame(r["pvals"])["high"].values
    qvals = pd.DataFrame(r["pvals_adj"])["high"].values
    return pd.DataFrame({
        "gene":          genes,
        "wilc_score":    scores,
        "wilc_logFC":    logfc,
        "wilc_p":        pvals,
        "wilc_q":        qvals,
    })


def annotate(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["program"] = df["gene"].map(ANNOT).fillna("")
    return df


def plot_volcano(df: pd.DataFrame, out: Path) -> None:
    if df.empty: return
    fig, ax = plt.subplots(figsize=(8, 6))
    x = df["mean_logFC"].values
    y = -np.log10(df["q_BH"].clip(lower=1e-30).values)
    # Default grey
    colors = np.array(["#bdc3c7"] * len(df), dtype=object)
    # Color by program annotation if present
    palette = {"hep": "#27ae60", "bil": "#2980b9",
               "plast": "#8e44ad", "stress": "#e67e22", "prolif": "#c0392b"}
    for prog, c in palette.items():
        mask = (df["program"] == prog).values
        colors[mask] = c
    sig = (df["q_BH"] < 0.05).values
    ax.scatter(x[~sig], y[~sig], s=10, c=colors[~sig], alpha=0.5)
    ax.scatter(x[sig],  y[sig],  s=14, c=colors[sig],  alpha=0.95,
               edgecolor="black", linewidth=0.3)
    # Label top hits in each direction
    top_up   = df[(df["q_BH"] < 0.05) & (df["mean_logFC"] > 0)]\
                 .nlargest(8, "mean_logFC")
    top_down = df[(df["q_BH"] < 0.05) & (df["mean_logFC"] < 0)]\
                 .nsmallest(8, "mean_logFC")
    for _, r in pd.concat([top_up, top_down]).iterrows():
        ax.annotate(r["gene"], (r["mean_logFC"],
                                 -np.log10(max(r["q_BH"], 1e-30))),
                    fontsize=8, alpha=0.8)
    ax.axvline(0, color="grey", linewidth=0.4)
    ax.axhline(-np.log10(0.05), color="grey", linewidth=0.4, linestyle="--")
    ax.set_xlabel("Mean per-donor logFC  (PROX1-high − PROX1-low)")
    ax.set_ylabel("−log10(BH q)")
    ax.set_title("Donor-controlled DE: PROX1-high vs PROX1-low in PSC")
    # legend
    handles = [plt.Line2D([0],[0], marker="o", color="w",
                          markerfacecolor=c, markersize=8, label=k)
               for k,c in palette.items()]
    handles.append(plt.Line2D([0],[0], marker="o", color="w",
                              markerfacecolor="#bdc3c7", markersize=8,
                              label="other"))
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="best")
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_top_genes(df: pd.DataFrame, out: Path, n: int = 20) -> None:
    if df.empty: return
    top_up   = df[df["q_BH"] < 0.05].nlargest(n, "mean_logFC")
    top_down = df[df["q_BH"] < 0.05].nsmallest(n, "mean_logFC")
    fig, axes = plt.subplots(1, 2, figsize=(13, max(4, 0.3 * n)))
    palette = {"hep": "#27ae60", "bil": "#2980b9",
               "plast": "#8e44ad", "stress": "#e67e22", "prolif": "#c0392b",
               "": "#7f8c8d"}
    for ax, sub, title, color_when_pos in [
        (axes[0], top_up,   f"Top {n} UP in PROX1-high (PSC)",   True),
        (axes[1], top_down, f"Top {n} DOWN in PROX1-high (PSC)", False),
    ]:
        if sub.empty:
            ax.set_title(title + "\n(no significant genes)")
            continue
        sub = sub.iloc[::-1]
        colors = [palette.get(p, "#7f8c8d") for p in sub["program"]]
        ax.barh(range(len(sub)), sub["mean_logFC"], color=colors)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels(sub["gene"], fontsize=9)
        ax.axvline(0, color="black", linewidth=0.4)
        ax.set_xlabel("Mean per-donor logFC")
        ax.set_title(title)
    handles = [plt.Line2D([0],[0], marker="s", color="w",
                          markerfacecolor=c, markersize=10, label=k)
               for k,c in palette.items() if k]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, -0.02))
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("PROX1-high vs PROX1-low DE within PSC (donor-controlled)")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    print(f"\n[1/5] Loading {CHOL_PATH}")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  full: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    adata = slice_psc(adata)
    print(f"  PSC subset: {adata.n_obs:,} cells")
    if PROX1 not in adata.var_names:
        print(f"ERROR: {PROX1} missing", file=sys.stderr); sys.exit(1)

    prox1 = get_gene_expression(adata, PROX1)
    strict_thr = float(np.percentile(prox1, 75))
    print(f"\n[2/5] PROX1-high threshold (top quartile within PSC) = {strict_thr:.3f}")
    adata = split_high_low(adata, strict_thr)
    print(f"  PSC cells PROX1-high: {(adata.obs['prox1_group']=='high').sum():,}  "
          f"PROX1-low: {(adata.obs['prox1_group']=='low').sum():,}")

    print(f"\n[3/5] Building per-donor pseudobulk means …")
    means_long, donor_counts = per_donor_pseudobulk(adata)
    n_ok = int(((donor_counts["n_high"] >= MIN_CELLS_PER_GROUP_PER_DONOR) &
                (donor_counts["n_low"]  >= MIN_CELLS_PER_GROUP_PER_DONOR)).sum())
    print(f"  donors with both groups >= {MIN_CELLS_PER_GROUP_PER_DONOR}: "
          f"{n_ok} / {len(donor_counts)}")
    print("  per-donor counts:")
    print(donor_counts.to_string(index=False))
    donor_counts.to_csv(RESULTS_DIR / "prox1_high_vs_low_psc_per_donor.tsv",
                        sep="\t", index=False)

    print(f"\n[4/5] Donor-controlled DE (one-sample t on per-donor logFC) …")
    de = donor_controlled_de(means_long)
    print(f"  genes with ≥3 informative donors: {len(de):,}")

    # Uncontrolled sanity check (Wilcoxon, ignores donor)
    print(f"\n  Uncontrolled Wilcoxon (scanpy.tl.rank_genes_groups) …")
    wilc = uncontrolled_wilcoxon(adata)
    de = de.merge(wilc, on="gene", how="left")
    de = annotate(de)
    de = de.sort_values("mean_logFC", ascending=False)
    de.to_csv(RESULTS_DIR / "prox1_high_vs_low_psc_de.tsv",
              sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'prox1_high_vs_low_psc_de.tsv'}")

    # Top tables
    sig = de[de["q_BH"] < 0.05]
    top_up   = sig.nlargest(TOP_N, "mean_logFC")
    top_down = sig.nsmallest(TOP_N, "mean_logFC")
    top = pd.concat([top_up.assign(direction="UP"),
                     top_down.assign(direction="DOWN")])
    cols = ["direction","gene","program","n_donors","mean_logFC","median_logFC",
            "frac_pos_donors","frac_neg_donors","p_one_sample_t","q_BH",
            "wilc_logFC","wilc_q"]
    top = top[cols]
    top.to_csv(RESULTS_DIR / "prox1_high_vs_low_psc_de_top.tsv",
               sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'prox1_high_vs_low_psc_de_top.tsv'}")

    print(f"\n--- Top {TOP_N} UP in PROX1-high (q<0.05, by mean_logFC) ---")
    print(top_up[["gene","program","mean_logFC","frac_pos_donors","q_BH"]]
          .to_string(index=False, float_format="%.3f"))
    print(f"\n--- Top {TOP_N} DOWN in PROX1-high (q<0.05, by mean_logFC) ---")
    print(top_down[["gene","program","mean_logFC","frac_pos_donors","q_BH"]]
          .to_string(index=False, float_format="%.3f"))

    print(f"\n[5/5] Plots …")
    plot_volcano(de, FIGURES_DIR / "prox1_high_low_de_volcano.png")
    plot_top_genes(de, FIGURES_DIR / "prox1_high_low_de_top_genes.png",
                   n=20)

    # Program-level summary
    if not sig.empty:
        print(f"\n--- Program-level direction counts (q<0.05) ---")
        sig_annot = sig.copy()
        sig_annot["dir"] = np.where(sig_annot["mean_logFC"] > 0, "UP", "DOWN")
        prog = (sig_annot[sig_annot["program"] != ""]
                .groupby(["program", "dir"]).size().unstack(fill_value=0))
        print(prog.to_string())


if __name__ == "__main__":
    main()
