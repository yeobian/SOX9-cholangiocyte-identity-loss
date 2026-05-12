"""
SOX9-low cholangiocytes as a disease-associated endpoint:
reverse-engineering candidate upstream programs.

Framing
-------
SOX9 collapse is the robust within-disease cholangiocyte signal in this
project (survives donor-level, per-cell, per-dataset, sensitivity, and
technology-stratified tests). SOX9 is treated here as a *marker of the
disease endpoint*, not as the cause. We use SOX9-low PSC cholangiocytes
as a probe to nominate which pathways may sit upstream of cholangiocyte
identity collapse.

No causal claim. No therapy claim. Pathways enriched or depleted in
SOX9-low cells are *candidates* that may suggest rescue axes for future
experimental validation.

Pipeline (load → filter → score → contrast → save)
--------------------------------------------------
1. Load cholangiocyte AnnData (PROX1-keeping refetch).
2. Subset to PSC only.
3. Define SOX9 expression tertiles within PSC; drop the middle.
4. For each curated pathway gene set, score per cell with sc.tl.score_genes.
5. Donor-controlled comparison: per donor, mean(high) − mean(low) per
   pathway. One-sample t-test across donors. BH FDR.
6. Stratify by technology (snRNA vs scRNA) and run separately.
7. Run gene-level pseudobulk DE (same machinery) on the full transcriptome
   for genes-lost / genes-gained tables.
8. Receptor score panel: candidate response receptors for each external
   signal (proxy for ligand-receptor analysis).
9. Sensitivity check: repeat with quartile and quintile splits.

Outputs
-------
results/sox9_de_pseudobulk_psc.tsv             gene-level DE
results/sox9_pathway_scores_psc.tsv            pathway scores per cell
results/sox9_pathway_donor_contrast_psc.tsv    pathway high − low per donor
results/sox9_pathway_summary_psc.tsv           cross-donor summary + verdict
results/sox9_receptor_summary_psc.tsv          receptor expression
results/sox9_threshold_sensitivity_psc.tsv     tertile / quartile / quintile
figures/sox9/sox9_pathway_summary.png
figures/sox9/sox9_de_volcano.png
figures/sox9/sox9_receptor_heatmap.png
figures/sox9/sox9_threshold_sensitivity.png

Usage
-----
    python src/sox9_mechanism_discovery.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import sparse
from scipy.stats import ttest_1samp

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression
from iad_score import LIVER_H5AD as DEFAULT_CHOL

CHOL_PATH    = DEFAULT_CHOL
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "sox9"

PSC = "primary sclerosing cholangitis"
SNRNA_DATASET = "1873a18a-66fd-4a4d-8277-a872c93f5b59"
MIN_CELLS_PER_GROUP_PER_DONOR = 10
GENE = "SOX9"


# ── Curated pathway gene sets ──────────────────────────────────────────────────
# Short, literature-anchored lists. Genes that fail HVG filtering are skipped at
# score time. Keep lists short to keep scores interpretable.
PATHWAYS = {
    # ── Cholangiocyte / biliary identity (expected DOWN in SOX9-low) ─────────
    "biliary_identity":       ["SOX9", "KRT19", "KRT7", "HNF1B", "EPCAM",
                                "CFTR", "ANXA4", "MUC1", "AQP1", "TFF3"],
    "cholangiocyte_function": ["CFTR", "AQP1", "SCTR", "SLC4A2", "ABCB4",
                                "ABCB11", "ABCC2"],

    # ── Progenitor / developmental TFs ───────────────────────────────────────
    "progenitor_TFs":         ["GATA4", "GATA6", "ONECUT1", "ONECUT2",
                                "FOXA1", "FOXA2", "HNF4A", "HNF1A"],

    # ── Maintenance signalling (Notch, BMP, Hippo) ───────────────────────────
    "notch_signalling":       ["NOTCH1", "NOTCH2", "NOTCH3", "JAG1", "JAG2",
                                "DLL1", "DLL3", "DLL4", "HES1", "HEY1",
                                "RBPJ", "MAML1"],
    "notch_targets":          ["HES1", "HEY1", "HEYL", "NRARP", "ID2"],
    "yap_taz_targets":        ["CTGF", "CYR61", "ANKRD1", "AMOTL2", "BIRC5",
                                "AXL", "TEAD1", "TEAD4"],
    "bmp_signalling":         ["BMP2", "BMP4", "BMP7", "BMPR1A", "BMPR1B",
                                "BMPR2", "ID1", "ID2", "ID3", "SMAD1",
                                "SMAD5", "SMAD9"],
    "wnt_signalling":         ["WNT2", "WNT5A", "WNT7B", "FZD1", "FZD4",
                                "FZD6", "FZD7", "LRP5", "LRP6", "LEF1",
                                "TCF7", "AXIN2", "LGR5", "RNF43"],

    # ── Disease / stress programs (expected UP in SOX9-low) ──────────────────
    "tgfb_signalling":        ["TGFB1", "TGFB2", "TGFBR1", "TGFBR2", "TGFBR3",
                                "SMAD2", "SMAD3", "SMAD7", "SERPINE1", "PMEPA1"],
    "tnf_signalling":         ["TNF", "TNFRSF1A", "TNFRSF1B", "NFKB1", "NFKB2",
                                "RELA", "RELB", "BIRC3", "TNFAIP3", "ICAM1",
                                "VCAM1"],
    "il6_signalling":         ["IL6", "IL6R", "IL6ST", "STAT3", "SOCS3",
                                "OSM", "OSMR"],
    "ifn_signalling":         ["IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2", "STAT1",
                                "STAT2", "IRF1", "IRF3", "IRF7", "ISG15",
                                "IFI6", "IFITM3", "MX1"],
    "bile_acid_stress":       ["NR1H4", "NR0B2", "FGFR4", "CYP7A1", "ABCB11",
                                "BAAT", "SLC10A1", "FXR_target_FGF19"],
    "er_stress_upr":          ["HSPA5", "ATF4", "ATF6", "DDIT3", "XBP1",
                                "EIF2AK3", "ERN1", "HSPA1A", "HSPA1B",
                                "HERPUD1"],
    "oxidative_stress":       ["HMOX1", "NQO1", "GPX1", "GPX4", "SOD2",
                                "TXN", "TXNRD1", "GCLC", "GCLM", "NFE2L2"],
    "apoptosis":              ["BAX", "BAK1", "BCL2", "BCL2L1", "MCL1",
                                "CASP3", "CASP7", "CASP8", "CASP9", "BBC3",
                                "PMAIP1", "BID"],
    "senescence":             ["CDKN2A", "CDKN1A", "CDKN2B", "TP53", "GLB1",
                                "LMNB1", "SERPINE1", "IL6", "CXCL8", "IGFBP7"],
    "emt":                    ["VIM", "FN1", "SNAI1", "SNAI2", "ZEB1", "ZEB2",
                                "TWIST1", "TWIST2", "CDH2", "CDH1",
                                "SERPINE1", "TGFB1"],
    "fibrosis_program":       ["COL1A1", "COL1A2", "COL3A1", "ACTA2", "TAGLN",
                                "FAP", "PDGFRB", "TIMP1", "TIMP2", "LOX",
                                "LOXL2"],
    "inflammation_general":   ["CXCL8", "CXCL10", "CCL2", "CCL20", "IL1B",
                                "IL18", "S100A8", "S100A9", "NFKB1"],
}

# ── Receptor panel: response receptors per external signal (proxy for L-R) ───
RECEPTORS = {
    "TGFb_receptors":   ["TGFBR1", "TGFBR2", "TGFBR3"],
    "TNF_receptors":    ["TNFRSF1A", "TNFRSF1B", "TNFRSF10A", "TNFRSF10B"],
    "IL6_receptors":    ["IL6R", "IL6ST", "OSMR", "LIFR"],
    "IFN_receptors":    ["IFNAR1", "IFNAR2", "IFNGR1", "IFNGR2"],
    "Notch_receptors":  ["NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4"],
    "Notch_ligands_on_self": ["JAG1", "JAG2", "DLL1", "DLL3", "DLL4"],
    "BMP_receptors":    ["BMPR1A", "BMPR1B", "BMPR2", "ACVR1", "ACVR2A"],
    "WNT_receptors":    ["FZD1", "FZD2", "FZD3", "FZD4", "FZD5", "FZD6",
                          "FZD7", "FZD8", "LRP5", "LRP6"],
    "Hedgehog":         ["PTCH1", "PTCH2", "SMO", "GLI1", "GLI2", "GLI3"],
    "Bile_acid_receptors": ["NR1H4", "NR0B2", "GPBAR1", "VDR"],
}


# ── Loaders / preprocessing ────────────────────────────────────────────────────

def load_psc(adata) -> "ad.AnnData":
    if "disease" not in adata.obs.columns:
        print("ERROR: no `disease` column.", file=sys.stderr); sys.exit(1)
    psc_mask = adata.obs["disease"].astype(str) == PSC
    sub = adata[psc_mask.values].copy()
    sub.obs["tech"] = np.where(
        sub.obs["dataset_id"].astype(str).values == SNRNA_DATASET,
        "snRNA", "scRNA",
    )
    return sub


def split_sox9_rank(adata: "ad.AnnData", lo_pct=33.33, hi_pct=66.67) -> tuple:
    """Rank-based top / bottom split within PSC.

    SOX9 expression is heavily zero-inflated (most cholangiocytes don't
    detect SOX9 mRNA, and the .X matrix is scaled so all zero-count cells
    share the same negative scaled value). Percentile-on-value splits fail
    because the 33rd and 67th percentiles equal the same scaled-zero value.

    We use rank-based splits with tie-breaking by row index, so the bottom
    `lo_pct`% of cells (by SOX9 rank) go to 'low' and the top (100−hi_pct)%
    go to 'high'. This guarantees exactly the requested fraction in each
    group regardless of zero-inflation.
    """
    sox9 = get_gene_expression(adata, GENE)
    n = len(sox9)
    order = np.argsort(sox9, kind="stable")          # ascending by SOX9, ties stable
    rank  = np.empty(n, dtype=int)
    rank[order] = np.arange(n)
    lo_cut = int(np.floor(n * lo_pct / 100))         # bottom lo_pct%  → 'low'
    hi_cut = int(np.ceil(n * hi_pct / 100))          # top (100-hi_pct)% → 'high'
    group = np.full(n, "mid", dtype=object)
    group[rank < lo_cut] = "low"
    group[rank >= hi_cut] = "high"
    adata.obs["sox9_group"] = pd.Categorical(group,
        categories=["high", "mid", "low"])
    adata.obs["sox9_expr"] = sox9
    sub = adata[adata.obs["sox9_group"].isin(["high", "low"])].copy()
    # Effective thresholds for reporting (per-rank-cut value of SOX9):
    sox9_sorted = sox9[order]
    lo_thr = float(sox9_sorted[max(lo_cut - 1, 0)])
    hi_thr = float(sox9_sorted[min(hi_cut, n - 1)])
    return sub, lo_thr, hi_thr


# ── Pseudobulk DE machinery ────────────────────────────────────────────────────

def per_donor_pseudobulk(adata) -> tuple:
    """Returns per-donor mean expression per gene per group, and donor counts."""
    obs = adata.obs[["donor_id", "sox9_group", "tech"]].copy()
    obs.columns = ["donor", "group", "tech"]
    X = adata.X if sparse.issparse(adata.X) else sparse.csr_matrix(adata.X)
    genes = adata.var_names.values

    rows = []
    counts = []
    for donor, dgrp in obs.groupby("donor", observed=True):
        n_high = int((dgrp["group"] == "high").sum())
        n_low  = int((dgrp["group"] == "low").sum())
        tech   = dgrp["tech"].iloc[0]
        counts.append({"donor": donor, "n_high": n_high,
                       "n_low": n_low, "tech": tech})
        if n_high < MIN_CELLS_PER_GROUP_PER_DONOR or \
           n_low  < MIN_CELLS_PER_GROUP_PER_DONOR:
            continue
        for grp_name, grp_idx in [
            ("high", dgrp.index[dgrp["group"] == "high"]),
            ("low",  dgrp.index[dgrp["group"] == "low"]),
        ]:
            pos = adata.obs_names.get_indexer(grp_idx)
            mean = np.asarray(X[pos].mean(axis=0)).ravel()
            for g, m in zip(genes, mean):
                rows.append({"donor": donor, "group": grp_name,
                             "gene": g, "mean": float(m), "tech": tech})
    return pd.DataFrame(rows), pd.DataFrame(counts)


def donor_controlled_de(means_long: pd.DataFrame) -> pd.DataFrame:
    if means_long.empty:
        return pd.DataFrame()
    wide = means_long.pivot_table(
        index=["donor", "gene"], columns="group",
        values="mean", aggfunc="first",
    ).reset_index()
    wide["logFC"] = wide["high"] - wide["low"]    # high − low; positive = lost in SOX9-low
    rows = []
    for gene, sub in wide.groupby("gene"):
        lfcs = sub["logFC"].dropna().values
        if len(lfcs) < 3:
            continue
        if len(set(lfcs)) > 1:
            t, p = ttest_1samp(lfcs, 0.0)
        else:
            t, p = 0.0, 1.0
        rows.append({
            "gene":            gene,
            "n_donors":        int(len(lfcs)),
            "mean_logFC_HminusL": float(np.mean(lfcs)),
            "median_logFC":    float(np.median(lfcs)),
            "frac_pos_donors": float((lfcs > 0).mean()),
            "frac_neg_donors": float((lfcs < 0).mean()),
            "t_stat":          float(t),
            "p_one_sample_t":  float(p),
        })
    df = pd.DataFrame(rows).sort_values("p_one_sample_t")
    df["rank"] = np.arange(1, len(df) + 1)
    df["q_BH"] = (df["p_one_sample_t"] * len(df) / df["rank"]).clip(upper=1.0)
    return df


# ── Pathway scoring + contrast ─────────────────────────────────────────────────

def score_pathways(adata: "ad.AnnData") -> dict:
    """Returns dict {pathway_name: per-cell score array}."""
    out = {}
    for name, genes in PATHWAYS.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            print(f"  [skip] {name}: only {len(present)}/{len(genes)} genes present.")
            continue
        key = f"_score_{name}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          random_state=42)
        out[name] = adata.obs[key].values.copy()
        del adata.obs[key]
        print(f"  scored {name:<26}  ({len(present)}/{len(genes)} genes)")
    return out


def score_receptors(adata: "ad.AnnData") -> dict:
    out = {}
    for name, genes in RECEPTORS.items():
        present = [g for g in genes if g in adata.var_names]
        if len(present) < 2:
            print(f"  [skip] receptor {name}: only {len(present)}/{len(genes)}.")
            continue
        key = f"_rec_{name}"
        sc.tl.score_genes(adata, gene_list=present, score_name=key,
                          random_state=42)
        out[name] = adata.obs[key].values.copy()
        del adata.obs[key]
        print(f"  scored receptors {name:<22}  ({len(present)}/{len(genes)})")
    return out


def pathway_donor_contrast(adata: "ad.AnnData",
                           scores: dict) -> pd.DataFrame:
    """Per donor with both arms, mean(high) − mean(low). Cross-donor t-test."""
    obs = adata.obs[["donor_id", "sox9_group", "tech"]].copy()
    obs.columns = ["donor", "group", "tech"]
    for name, vec in scores.items():
        obs[f"S_{name}"] = vec

    rows = []
    for donor, dgrp in obs.groupby("donor", observed=True):
        hi = dgrp[dgrp["group"] == "high"]
        lo = dgrp[dgrp["group"] == "low"]
        if len(hi) < MIN_CELLS_PER_GROUP_PER_DONOR or \
           len(lo) < MIN_CELLS_PER_GROUP_PER_DONOR:
            continue
        for name in scores:
            col = f"S_{name}"
            rows.append({
                "donor": donor, "tech": dgrp["tech"].iloc[0],
                "pathway": name,
                "delta_HminusL": float(hi[col].mean() - lo[col].mean()),
                "high_mean": float(hi[col].mean()),
                "low_mean":  float(lo[col].mean()),
                "n_high":  int(len(hi)),
                "n_low":   int(len(lo)),
            })
    if not rows:
        return pd.DataFrame(), pd.DataFrame()
    long = pd.DataFrame(rows)
    summary_rows = []
    for pw, sub in long.groupby("pathway"):
        deltas = sub["delta_HminusL"].values
        if len(deltas) < 3 or len(set(deltas)) < 2:
            t, p = float("nan"), float("nan")
        else:
            t, p = ttest_1samp(deltas, 0.0)
        summary_rows.append({
            "pathway":      pw,
            "n_donors":     int(len(deltas)),
            "mean_delta":   float(np.mean(deltas)),
            "median_delta": float(np.median(deltas)),
            "frac_pos_donors": float((deltas > 0).mean()),
            "frac_neg_donors": float((deltas < 0).mean()),
            "t_stat":       float(t),
            "p":            float(p),
        })
    summary = pd.DataFrame(summary_rows).sort_values("p")
    summary["rank"] = np.arange(1, len(summary) + 1)
    summary["q_BH"] = (summary["p"] * len(summary) / summary["rank"]).clip(upper=1.0)
    return long, summary


def category_label(row) -> str:
    """A/B/C interpretation category from a pathway summary row."""
    if np.isnan(row["mean_delta"]) or row["n_donors"] < 4:
        return "underpowered"
    if row["mean_delta"] > 0.05 and row["frac_pos_donors"] >= 5/7:
        return "A_lost_in_SOX9low"
    if row["mean_delta"] < -0.05 and row["frac_neg_donors"] >= 5/7:
        return "B_gained_in_SOX9low"
    return "ambiguous"


# ── Plots ──────────────────────────────────────────────────────────────────────

def plot_pathway_summary(summary: pd.DataFrame, out: Path) -> None:
    if summary.empty:
        return
    df = summary.sort_values("mean_delta").copy()
    fig, ax = plt.subplots(figsize=(8, 0.35 * len(df) + 1.5))
    palette = {
        "A_lost_in_SOX9low":   "#2980b9",  # lost in SOX9-low (HIGHER in SOX9-high)
        "B_gained_in_SOX9low": "#c0392b",  # gained in SOX9-low (LOWER in SOX9-high)
        "ambiguous":           "#7f8c8d",
        "underpowered":        "#bdc3c7",
    }
    df["cat"] = df.apply(category_label, axis=1)
    colors = [palette.get(c, "#7f8c8d") for c in df["cat"]]
    sig = (df["q_BH"] < 0.1).values
    ax.barh(range(len(df)), df["mean_delta"], color=colors,
            edgecolor=["black" if s else "none" for s in sig])
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["pathway"], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_xlabel("Mean per-donor pathway score Δ (SOX9-high − SOX9-low)")
    ax.set_title("Pathway scores: SOX9-high vs SOX9-low PSC cholangiocytes\n"
                 "Positive Δ = pathway lost in SOX9-low.  "
                 "Black edge = q<0.1.")
    handles = [
        plt.Line2D([0],[0], marker="s", color="w",
                   markerfacecolor=palette["A_lost_in_SOX9low"], markersize=10,
                   label="A: lost maintenance (down in SOX9-low)"),
        plt.Line2D([0],[0], marker="s", color="w",
                   markerfacecolor=palette["B_gained_in_SOX9low"], markersize=10,
                   label="B: gained disease/stress (up in SOX9-low)"),
        plt.Line2D([0],[0], marker="s", color="w",
                   markerfacecolor=palette["ambiguous"], markersize=10,
                   label="ambiguous / small effect"),
    ]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="best")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_de_volcano(de: pd.DataFrame, out: Path) -> None:
    if de.empty: return
    fig, ax = plt.subplots(figsize=(8, 6))
    x = de["mean_logFC_HminusL"].values
    y = -np.log10(de["q_BH"].clip(lower=1e-30).values)
    sig = (de["q_BH"] < 0.1).values
    ax.scatter(x[~sig], y[~sig], s=8, c="#bdc3c7", alpha=0.4)
    ax.scatter(x[sig],  y[sig],  s=14, c="#34495e", alpha=0.85,
               edgecolor="black", linewidth=0.3)
    top_lost = de[(de["q_BH"] < 0.1) & (de["mean_logFC_HminusL"] > 0)]\
                 .nlargest(8, "mean_logFC_HminusL")
    top_gain = de[(de["q_BH"] < 0.1) & (de["mean_logFC_HminusL"] < 0)]\
                 .nsmallest(8, "mean_logFC_HminusL")
    for _, r in pd.concat([top_lost, top_gain]).iterrows():
        ax.annotate(r["gene"], (r["mean_logFC_HminusL"],
                                  -np.log10(max(r["q_BH"], 1e-30))),
                    fontsize=8, alpha=0.85)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.axhline(-np.log10(0.1), color="grey", linewidth=0.4, linestyle="--")
    ax.set_xlabel("Mean per-donor logFC (SOX9-high − SOX9-low)")
    ax.set_ylabel("-log10(BH q)")
    ax.set_title("Donor-controlled DE: SOX9-high vs SOX9-low PSC cholangiocytes\n"
                 "Right = lost in SOX9-low.  Left = gained in SOX9-low.")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_receptor_heatmap(long: pd.DataFrame, summary: pd.DataFrame,
                          out: Path) -> None:
    if long.empty or summary.empty: return
    summary = summary.sort_values("mean_delta")
    fig, ax = plt.subplots(figsize=(6, 0.45 * len(summary) + 1.5))
    pivot = (long.pivot_table(index="pathway", columns="donor",
                              values="delta_HminusL", aggfunc="mean")
                 .reindex(summary["pathway"]))
    sns.heatmap(pivot, cmap="RdBu_r", center=0, ax=ax,
                cbar_kws={"label": "high − low pathway score"},
                annot=False, linewidths=0.4)
    ax.set_title("Per-donor receptor expression contrast\n"
                 "Red = up in SOX9-high (lost in SOX9-low). Blue = up in SOX9-low.")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def plot_threshold_sensitivity(sens: pd.DataFrame, out: Path) -> None:
    if sens.empty: return
    pivot = sens.pivot(index="pathway", columns="threshold", values="mean_delta")
    if pivot.empty: return
    pivot = pivot.reindex(pivot.mean(axis=1).sort_values().index)
    fig, ax = plt.subplots(figsize=(6, 0.35 * len(pivot) + 1.5))
    sns.heatmap(pivot, cmap="RdBu_r", center=0, ax=ax, annot=True, fmt=".2f",
                annot_kws={"fontsize": 7},
                cbar_kws={"label": "mean Δ (high − low)"})
    ax.set_title("Pathway Δ across SOX9 split thresholds (tertile/quartile/quintile)")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("SOX9 mechanism discovery — PSC cholangiocytes")
    print("Frame: SOX9-low as disease endpoint, candidate upstream pathways")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    print(f"\n[1/8] Loading {CHOL_PATH}")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    print(f"\n[2/8] Subsetting to PSC …")
    psc = load_psc(adata)
    print(f"  PSC cells: {psc.n_obs:,}   "
          f"snRNA: {(psc.obs['tech']=='snRNA').sum():,}   "
          f"scRNA: {(psc.obs['tech']=='scRNA').sum():,}")

    if GENE not in psc.var_names:
        print(f"ERROR: {GENE} missing from var.", file=sys.stderr); sys.exit(1)

    print(f"\n[3/8] SOX9 rank-based tertile split within PSC (zero-inflated safe) …")
    psc_hl, lo_thr, hi_thr = split_sox9_rank(psc, 33.33, 66.67)
    pct_zero = float((get_gene_expression(psc, GENE) <= 0).mean() * 100)
    print(f"  SOX9 zero-or-below cells in PSC: {pct_zero:.1f}% "
          f"(this is why we use rank-based, not percentile-on-value)")
    print(f"  effective rank cutoffs: low SOX9 ≤ {lo_thr:.3f}  "
          f"high SOX9 ≥ {hi_thr:.3f}")
    print(f"  PSC cells SOX9-high: {(psc_hl.obs['sox9_group']=='high').sum():,}  "
          f"SOX9-low: {(psc_hl.obs['sox9_group']=='low').sum():,}")

    print(f"\n[4/8] Donor pseudobulk + DE (all PSC, tertile split)")
    long_means, donor_counts = per_donor_pseudobulk(psc_hl)
    print(f"  donor cell counts:")
    print(donor_counts.to_string(index=False))
    de = donor_controlled_de(long_means)
    if not de.empty:
        de.to_csv(RESULTS_DIR / "sox9_de_pseudobulk_psc.tsv",
                  sep="\t", index=False)
        print(f"  → {RESULTS_DIR / 'sox9_de_pseudobulk_psc.tsv'}  "
              f"({len(de):,} genes, {(de['q_BH']<0.1).sum()} at q<0.1)")

    print(f"\n[5/8] Scoring pathways on PSC SOX9-high + SOX9-low cells …")
    scores = score_pathways(psc_hl)
    print(f"\n[6/8] Scoring receptor panels …")
    rec_scores = score_receptors(psc_hl)

    print(f"\n[7/8] Per-donor pathway contrasts …")
    pw_long, pw_summary = pathway_donor_contrast(psc_hl, scores)
    rec_long, rec_summary = pathway_donor_contrast(psc_hl, rec_scores)

    if not pw_summary.empty:
        pw_summary["category"] = pw_summary.apply(category_label, axis=1)
        pw_summary = pw_summary.sort_values("mean_delta")
        pw_long.to_csv(RESULTS_DIR / "sox9_pathway_donor_contrast_psc.tsv",
                       sep="\t", index=False)
        pw_summary.to_csv(RESULTS_DIR / "sox9_pathway_summary_psc.tsv",
                          sep="\t", index=False)
        print(f"  → {RESULTS_DIR / 'sox9_pathway_summary_psc.tsv'}")
        print("\n--- Pathway summary (sorted by mean Δ) ---")
        print(pw_summary[["pathway", "n_donors", "mean_delta",
                           "frac_pos_donors", "frac_neg_donors",
                           "p", "q_BH", "category"]]
              .to_string(index=False, float_format="%.3f"))

    if not rec_summary.empty:
        rec_summary["category"] = rec_summary.apply(category_label, axis=1)
        rec_summary = rec_summary.sort_values("mean_delta")
        rec_long.to_csv(RESULTS_DIR / "sox9_receptor_donor_contrast_psc.tsv",
                        sep="\t", index=False)
        rec_summary.to_csv(RESULTS_DIR / "sox9_receptor_summary_psc.tsv",
                           sep="\t", index=False)
        print(f"  → {RESULTS_DIR / 'sox9_receptor_summary_psc.tsv'}")
        print("\n--- Receptor summary (sorted by mean Δ) ---")
        print(rec_summary[["pathway", "n_donors", "mean_delta",
                            "frac_pos_donors", "frac_neg_donors",
                            "p", "q_BH", "category"]]
              .to_string(index=False, float_format="%.3f"))

    # [8/8] Threshold sensitivity
    print(f"\n[8/8] Threshold sensitivity (quartile, quintile) …")
    sens_rows = []
    for thr_name, lo_pct, hi_pct in [("tertile_33_67", 33.33, 66.67),
                                      ("quartile_25_75", 25.0, 75.0),
                                      ("quintile_20_80", 20.0, 80.0)]:
        sub, _, _ = split_sox9_rank(psc, lo_pct, hi_pct)
        sub_scores = scores if thr_name == "tertile_33_67" else score_pathways(sub)
        _, s_sum = pathway_donor_contrast(sub, sub_scores)
        if not s_sum.empty:
            s_sum["threshold"] = thr_name
            sens_rows.append(s_sum[["pathway", "threshold", "mean_delta",
                                     "p", "q_BH"]])
    if sens_rows:
        sens = pd.concat(sens_rows, ignore_index=True)
        sens.to_csv(RESULTS_DIR / "sox9_threshold_sensitivity_psc.tsv",
                    sep="\t", index=False)
        print(f"  → {RESULTS_DIR / 'sox9_threshold_sensitivity_psc.tsv'}")
        plot_threshold_sensitivity(sens,
            FIGURES_DIR / "sox9_threshold_sensitivity.png")

    # Plots
    print(f"\n--- Plots ---")
    if not pw_summary.empty:
        plot_pathway_summary(pw_summary,
                             FIGURES_DIR / "sox9_pathway_summary.png")
    if not de.empty:
        plot_de_volcano(de, FIGURES_DIR / "sox9_de_volcano.png")
    if not rec_long.empty:
        plot_receptor_heatmap(rec_long, rec_summary,
                              FIGURES_DIR / "sox9_receptor_heatmap.png")

    # Headline
    print(f"\n--- Headline ---")
    if not pw_summary.empty:
        A = pw_summary[pw_summary["category"] == "A_lost_in_SOX9low"]
        B = pw_summary[pw_summary["category"] == "B_gained_in_SOX9low"]
        print(f"  A. Lost maintenance (DOWN in SOX9-low): {len(A)} pathways")
        for _, r in A.iterrows():
            print(f"     {r['pathway']:<26}  Δ={r['mean_delta']:+.3f}  "
                  f"q={r['q_BH']:.2g}  ({r['n_donors']} donors, "
                  f"{r['frac_pos_donors']:.0%} agree)")
        print(f"  B. Gained disease/stress (UP in SOX9-low): {len(B)} pathways")
        for _, r in B.iterrows():
            print(f"     {r['pathway']:<26}  Δ={r['mean_delta']:+.3f}  "
                  f"q={r['q_BH']:.2g}  ({r['n_donors']} donors, "
                  f"{r['frac_neg_donors']:.0%} agree)")


if __name__ == "__main__":
    main()
