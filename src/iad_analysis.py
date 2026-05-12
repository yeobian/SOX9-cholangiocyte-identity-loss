"""
IAD Diagnostic Gene Panel — Liver vs. Pancreas Cholangiocytes
==============================================================

Pipeline stage: COMPUTE → SAVE.

Anchor gene: KRT19 (canonical cholangiocyte marker; loss of KRT19+ ducts
defines Idiopathic Adulthood Ductopenia). This mirrors the role of PROX1 in
the original Pancreas project — a single anchor gene whose co-expression
network defines the cell program of interest.

Strategy
--------
1. Load liver-cholangiocyte and pancreas-ductal AnnData files
2. Compute Spearman ρ between every gene and KRT19 in each tissue
3. Compare the two co-expression programs
   - Shared genes  → universal biliary-tree program (less diagnostic)
   - Liver-only    → cholangiocyte-specific program lost in IAD
                     **these are the diagnostic-panel candidates**
4. gProfiler enrichment on each set
5. Write a ranked diagnostic panel TSV with:
       gene, rho_liver, rho_pancreas, specificity, panel_rank

Outputs
-------
figures/comparison/krt19_umap.png
figures/comparison/krt19_violin.png
figures/comparison/top_coexpressed.png
figures/comparison/venn_overlap.png
figures/{liver,pancreas,comparison}/enrichment.png
results/liver_krt19_coexpression.tsv
results/pancreas_krt19_coexpression.tsv
results/gene_comparison.tsv
results/iad_diagnostic_panel.tsv     ← key deliverable
results/enrichment_*.tsv

Usage
-----
    python src/iad_analysis.py
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
from utils import (
    ensure_dir,
    get_gene_expression,
    compute_gene_correlation,
)

# ── Config ────────────────────────────────────────────────────────────────────
DATA_DIR    = Path("data")
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures")

LIVER_H5AD    = DATA_DIR / "liver"    / "liver_cholangiocytes.h5ad"
PANCREAS_H5AD = DATA_DIR / "pancreas" / "pancreas_ductal.h5ad"

TARGET_GENE = "KRT19"
TOP_N       = 50      # genes carried forward for overlap / enrichment
PANEL_SIZE  = 20      # size of the final diagnostic panel


# ── Load ──────────────────────────────────────────────────────────────────────

def load_data() -> tuple[ad.AnnData, ad.AnnData]:
    for path, label in [(LIVER_H5AD, "liver"), (PANCREAS_H5AD, "pancreas")]:
        if not path.exists():
            print(f"ERROR: {label} H5AD not found at {path}", file=sys.stderr)
            print(f"Run:  python src/preprocess_{label}.py", file=sys.stderr)
            sys.exit(1)

    print("Loading liver cholangiocytes …")
    liver = sc.read_h5ad(LIVER_H5AD)
    print(f"  {liver.n_obs:,} cells × {liver.n_vars:,} genes")

    print("Loading pancreas ductal cells …")
    pancreas = sc.read_h5ad(PANCREAS_H5AD)
    print(f"  {pancreas.n_obs:,} cells × {pancreas.n_vars:,} genes")
    return liver, pancreas


# ── Anchor presence ───────────────────────────────────────────────────────────

def report_anchor(adata: ad.AnnData, tissue: str) -> None:
    if TARGET_GENE not in adata.var_names:
        print(f"  {tissue}: {TARGET_GENE} NOT FOUND.")
        return
    expr = get_gene_expression(adata, TARGET_GENE)
    pct  = (expr > 0).mean() * 100
    mean = expr[expr > 0].mean() if pct > 0 else 0.0
    print(f"  {tissue}: {TARGET_GENE} in {pct:.1f}% of cells | "
          f"mean expr (positive cells) = {mean:.3f}")


# ── Figures ───────────────────────────────────────────────────────────────────

def plot_umap(liver: ad.AnnData, pancreas: ad.AnnData, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, adata, title in zip(axes, [liver, pancreas], ["Liver", "Pancreas"]):
        if TARGET_GENE in adata.var_names and "X_umap" in adata.obsm:
            sc.pl.umap(adata, color=TARGET_GENE, ax=ax, show=False,
                       title=f"{title}: {TARGET_GENE}", colorbar_loc="right")
        else:
            ax.text(0.5, 0.5, f"{TARGET_GENE}\nnot available",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def _cell_type_col(adata: ad.AnnData) -> str:
    for col in ("cell_type", "celltype", "Cell_type", "leiden", "cluster"):
        if col in adata.obs.columns:
            return col
    return adata.obs.columns[0]


def plot_violin(liver: ad.AnnData, pancreas: ad.AnnData, out: Path) -> None:
    records = []
    for adata, tissue in [(liver, "Liver"), (pancreas, "Pancreas")]:
        if TARGET_GENE not in adata.var_names:
            continue
        expr = get_gene_expression(adata, TARGET_GENE)
        ct   = adata.obs[_cell_type_col(adata)].astype(str).values
        for e, c in zip(expr, ct):
            records.append({"tissue": tissue, "cell_type": c, TARGET_GENE: e})

    if not records:
        return

    df = pd.DataFrame(records)
    fig, axes = plt.subplots(1, 2, figsize=(18, 6))
    for ax, tissue in zip(axes, ["Liver", "Pancreas"]):
        sub = df[df["tissue"] == tissue]
        if sub.empty:
            ax.set_title(f"{tissue}: no data")
            continue
        order = (
            sub.groupby("cell_type")[TARGET_GENE]
            .mean().sort_values(ascending=False).index
        )
        sns.violinplot(
            data=sub, x="cell_type", y=TARGET_GENE,
            order=order, ax=ax, hue="cell_type",
            palette="Set2", density_norm="width",
            inner="quartile", legend=False,
        )
        ax.set_title(f"{tissue}: {TARGET_GENE} by cell type")
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(list(order), rotation=45, ha="right", fontsize=8)
        ax.set_xlabel("")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def plot_top_coexpr(
    l_corr: pd.Series, p_corr: pd.Series, out: Path, n: int = 30,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(16, 8))
    for ax, corr, title in zip(axes, [l_corr, p_corr], ["Liver", "Pancreas"]):
        if corr.empty:
            ax.set_title(f"{title}: no data")
            continue
        top    = corr.head(n)
        colors = ["#c0392b" if v > 0 else "#2980b9" for v in top.values]
        ax.barh(range(len(top)), top.values[::-1], color=colors[::-1])
        ax.set_yticks(range(len(top)))
        ax.set_yticklabels(top.index[::-1], fontsize=9)
        ax.set_xlabel(f"Spearman ρ with {TARGET_GENE}")
        ax.set_title(f"{title}: top {n} {TARGET_GENE} co-expressed genes")
        ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def plot_overlap(comp: dict, out: Path) -> None:
    cats   = ["Shared", "Liver only", "Pancreas only"]
    counts = [
        len(comp["shared"]),
        len(comp["liver_unique"]),
        len(comp["pancreas_unique"]),
    ]
    colors = ["#8e44ad", "#16a085", "#2980b9"]
    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(cats, counts, color=colors, edgecolor="white", width=0.5)
    for bar, cnt in zip(bars, counts):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3, str(cnt),
            ha="center", va="bottom", fontweight="bold",
        )
    ax.set_ylabel(f"# of top-{TOP_N} {TARGET_GENE} co-expressed genes")
    ax.set_title(f"Overlap of {TARGET_GENE} co-expression programs\n(Liver vs. Pancreas)")
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Compute ───────────────────────────────────────────────────────────────────

def compute_coexpr(adata: ad.AnnData, tissue: str) -> pd.Series:
    if TARGET_GENE not in adata.var_names:
        print(f"  {tissue}: {TARGET_GENE} not in gene list — skipping.")
        return pd.Series(dtype=float)
    print(f"\n{tissue}: computing Spearman ρ with {TARGET_GENE} …")
    return compute_gene_correlation(adata, TARGET_GENE, method="spearman")


def overlap(l_corr: pd.Series, p_corr: pd.Series) -> dict:
    l_top = set(l_corr.head(TOP_N).index)
    p_top = set(p_corr.head(TOP_N).index)
    return {
        "shared":          sorted(l_top & p_top),
        "liver_unique":    sorted(l_top - p_top),
        "pancreas_unique": sorted(p_top - l_top),
    }


# Housekeeping / non-informative gene prefixes to exclude from the panel.
# These dominate any naive co-expression panel because they're highly expressed
# everywhere; they're not useful diagnostic markers.
HOUSEKEEPING_PREFIXES = (
    "RPS",   # cytosolic small ribosomal subunit
    "RPL",   # cytosolic large ribosomal subunit
    "MRPS",  # mitochondrial ribosomal small
    "MRPL",  # mitochondrial ribosomal large
    "MT-",   # mitochondrial-encoded
    "HLA-",  # MHC class I/II — broad expression
    "HIST",  # histones
    "H1-", "H2A", "H2B", "H3-", "H4-",  # histone variants
)
HOUSEKEEPING_EXACT = {
    # Cytoskeletal / structural
    "ACTB", "ACTG1", "B2M", "GAPDH", "TUBA1A", "TUBA1B", "TUBB",
    # Translation machinery / abundant cytosolic
    "EEF1A1", "EEF1A2", "EEF2", "PABPC1", "TPT1", "FTL", "FTH1",
    # Thymosins — universal cytoskeletal regulators
    "TMSB4X", "TMSB10",
    # Ubiquitous lncRNAs
    "MALAT1", "NEAT1", "XIST",
    # Immediate-early / stress-response transcription factors
    "FOS", "FOSB", "FOSL1", "FOSL2",
    "JUN", "JUNB", "JUND",
    "EGR1", "EGR2", "EGR3",
    "ATF3",
    "MAFF", "MAFG", "MAFK",
    "NR4A1", "NR4A2", "NR4A3",
    "IER2", "IER3", "IER5",
    # RNA-binding stress regulators
    "ZFP36", "ZFP36L1", "ZFP36L2",
    # Heat-shock proteins (stress-induced, not lineage-specific)
    "HSPA1A", "HSPA1B", "HSPA6", "HSPA8",
    "HSPB1", "HSP90AA1", "HSP90AB1",
    "DNAJB1",
    # Oxidative-stress / hypoxia leak-through
    "DDIT3", "DDIT4", "GADD45B",
    "BTG1", "BTG2",
}


def _is_housekeeping(gene: str) -> bool:
    if gene in HOUSEKEEPING_EXACT:
        return True
    return any(gene.startswith(p) for p in HOUSEKEEPING_PREFIXES)


def build_diagnostic_panel(
    l_corr: pd.Series,
    p_corr: pd.Series,
    n: int = PANEL_SIZE,
    min_rho_liver: float = 0.30,
    max_rho_pancreas: float = 0.50,
) -> pd.DataFrame:
    """
    Diagnostic-panel candidates = cholangiocyte-specific markers most likely
    lost in IAD biopsies. Panel is built only from genes present in BOTH
    correlation tables (so missing-data doesn't inflate specificity), with
    housekeeping and immediate-early genes filtered out.

    Inclusion rules
    ---------------
    - gene present in both l_corr and p_corr (intersection)
    - rho_liver    >= min_rho_liver        (must track KRT19 in liver)
    - rho_pancreas <= max_rho_pancreas     (must NOT track KRT19 in pancreas)
    - gene is not a housekeeping / IEG marker

    Ranking: specificity = rho_liver - rho_pancreas, descending.
    """
    cols = ["gene", "rho_liver", "rho_pancreas", "specificity", "panel_rank"]
    if l_corr.empty or p_corr.empty:
        print("  build_diagnostic_panel: missing one tissue's correlation — empty panel.")
        return pd.DataFrame(columns=cols)

    shared_genes = l_corr.index.intersection(p_corr.index)
    print(f"  Genes scored in both tissues: {len(shared_genes):,}")

    df = pd.DataFrame({
        "rho_liver":    l_corr.reindex(shared_genes),
        "rho_pancreas": p_corr.reindex(shared_genes),
    })
    df["specificity"] = df["rho_liver"] - df["rho_pancreas"]

    # Apply filters
    keep_thresh = (df["rho_liver"] >= min_rho_liver) & (df["rho_pancreas"] <= max_rho_pancreas)
    keep_clean  = ~df.index.to_series().map(_is_housekeeping)
    df = df[keep_thresh & keep_clean]
    print(
        f"  After thresholds (ρ_liv≥{min_rho_liver}, ρ_pan≤{max_rho_pancreas}) "
        f"and housekeeping filter: {len(df):,} candidates."
    )

    df = df.sort_values("specificity", ascending=False)
    df.insert(0, "gene", df.index)
    df.reset_index(drop=True, inplace=True)
    df["panel_rank"] = df.index + 1
    return df.head(n)


# ── gProfiler enrichment ──────────────────────────────────────────────────────

def run_enrichment(genes: list[str], label: str) -> pd.DataFrame:
    if not genes:
        return pd.DataFrame()
    try:
        from gprofiler import GProfiler
        gp = GProfiler(return_dataframe=True)
        results = gp.profile(
            organism="hsapiens",
            query=genes,
            sources=["GO:BP", "GO:MF", "KEGG", "REAC"],
            significance_threshold_method="fdr",
            user_threshold=0.05,
            no_evidences=False,
        )
    except Exception as exc:
        print(f"  {label}: gProfiler unavailable ({type(exc).__name__}) — skipping.")
        return pd.DataFrame()

    if results.empty:
        print(f"  {label}: no significant terms (FDR < 0.05)")
    else:
        print(f"  {label}: {len(results)} significant terms")
    return results


def plot_enrichment(df: pd.DataFrame, out: Path, title: str, n: int = 15) -> None:
    if df.empty:
        return
    top = df.nsmallest(n, "p_value")[["name", "p_value"]].copy()
    top["-log10p"] = -np.log10(top["p_value"].clip(lower=1e-30))
    top = top.sort_values("-log10p")
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(range(len(top)), top["-log10p"], color="#27ae60", edgecolor="white")
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top["name"], fontsize=9)
    ax.set_xlabel("−log₁₀(p-value)")
    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Save ──────────────────────────────────────────────────────────────────────

def save_results(
    l_corr: pd.Series, p_corr: pd.Series, comp: dict,
    panel: pd.DataFrame,
    e_liv: pd.DataFrame, e_pan: pd.DataFrame, e_shared: pd.DataFrame,
) -> None:
    ensure_dir(RESULTS_DIR)

    if not l_corr.empty:
        l_corr.rename("spearman_rho").to_csv(
            RESULTS_DIR / "liver_krt19_coexpression.tsv", sep="\t", header=True,
        )
    if not p_corr.empty:
        p_corr.rename("spearman_rho").to_csv(
            RESULTS_DIR / "pancreas_krt19_coexpression.tsv", sep="\t", header=True,
        )

    rows = (
        [("shared",          g) for g in comp.get("shared", [])]
        + [("liver_unique",    g) for g in comp.get("liver_unique", [])]
        + [("pancreas_unique", g) for g in comp.get("pancreas_unique", [])]
    )
    if rows:
        pd.DataFrame(rows, columns=["category", "gene"]).to_csv(
            RESULTS_DIR / "gene_comparison.tsv", sep="\t", index=False,
        )

    if not panel.empty:
        panel.to_csv(
            RESULTS_DIR / "iad_diagnostic_panel.tsv", sep="\t", index=False,
        )

    for df, name in [
        (e_liv,    "enrichment_liver"),
        (e_pan,    "enrichment_pancreas"),
        (e_shared, "enrichment_shared"),
    ]:
        if not df.empty:
            df.to_csv(RESULTS_DIR / f"{name}.tsv", sep="\t", index=False)

    print(f"  Results written to {RESULTS_DIR}/")


# ── Console summary ───────────────────────────────────────────────────────────

def print_summary(
    l_corr: pd.Series, p_corr: pd.Series,
    comp: dict, panel: pd.DataFrame,
) -> None:
    print("\n" + "=" * 60)
    print("IAD DIAGNOSTIC PANEL SUMMARY")
    print("=" * 60)

    for corr, tissue in [(l_corr, "Liver"), (p_corr, "Pancreas")]:
        if not corr.empty:
            print(f"\nTop 10 {TARGET_GENE} co-expressed genes — {tissue}:")
            for gene, rho in corr.head(10).items():
                print(f"  {gene:<22} ρ = {rho:+.4f}")

    if comp.get("shared"):
        print(f"\nShared between tissues (top {TOP_N} each, n={len(comp['shared'])}):")
        for g in comp["shared"][:20]:
            print(f"  {g}")

    if not panel.empty:
        print(f"\nProposed diagnostic panel (top {len(panel)} liver-specific genes):")
        for _, row in panel.iterrows():
            print(
                f"  #{int(row['panel_rank']):>2}  {row['gene']:<18}"
                f" ρ_liv={row['rho_liver']:+.3f}"
                f" ρ_pan={row['rho_pancreas']:+.3f}"
                f" spec={row['specificity']:+.3f}"
            )

    print("\nFull results in results/  |  Figures in figures/")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    for subdir in ("liver", "pancreas", "comparison"):
        ensure_dir(FIGURES_DIR / subdir)
    sc.settings.verbosity = 1
    sc.settings.set_figure_params(dpi=150, facecolor="white")

    print("=" * 60)
    print("IAD DIAGNOSTIC GENE PANEL ANALYSIS")
    print("Liver cholangiocytes vs. Pancreas ductal cells")
    print(f"Anchor gene: {TARGET_GENE}")
    print("=" * 60)

    liver, pancreas = load_data()

    print(f"\n--- {TARGET_GENE} expression ---")
    report_anchor(liver,    "Liver")
    report_anchor(pancreas, "Pancreas")

    print("\n--- Generating expression plots ---")
    plot_umap(liver, pancreas,   FIGURES_DIR / "comparison" / "krt19_umap.png")
    plot_violin(liver, pancreas, FIGURES_DIR / "comparison" / "krt19_violin.png")

    print("\n--- Co-expression ---")
    l_corr = compute_coexpr(liver,    "Liver")
    p_corr = compute_coexpr(pancreas, "Pancreas")

    print("\n--- Generating co-expression plots ---")
    plot_top_coexpr(l_corr, p_corr, FIGURES_DIR / "comparison" / "top_coexpressed.png")

    comp = {}
    if not l_corr.empty and not p_corr.empty:
        comp = overlap(l_corr, p_corr)
        plot_overlap(comp, FIGURES_DIR / "comparison" / "venn_overlap.png")

    print("\n--- Building diagnostic panel ---")
    panel = build_diagnostic_panel(l_corr, p_corr)

    print("\n--- Pathway enrichment (gProfiler) ---")
    e_liv    = run_enrichment(list(l_corr.head(TOP_N).index),    "Liver top genes")
    e_pan    = run_enrichment(list(p_corr.head(TOP_N).index),    "Pancreas top genes")
    e_shared = run_enrichment(comp.get("shared", []),            "Shared genes")

    plot_enrichment(e_liv,    FIGURES_DIR / "liver"      / "enrichment.png",
                    f"Liver: {TARGET_GENE} co-expression pathways")
    plot_enrichment(e_pan,    FIGURES_DIR / "pancreas"   / "enrichment.png",
                    f"Pancreas: {TARGET_GENE} co-expression pathways")
    plot_enrichment(e_shared, FIGURES_DIR / "comparison" / "enrichment_shared.png",
                    f"Shared {TARGET_GENE} co-expression pathways")

    print("\n--- Saving results ---")
    save_results(l_corr, p_corr, comp, panel, e_liv, e_pan, e_shared)

    print_summary(l_corr, p_corr, comp, panel)
    print("\nAnalysis complete.")


if __name__ == "__main__":
    main()
