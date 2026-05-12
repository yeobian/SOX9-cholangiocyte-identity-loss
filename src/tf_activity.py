"""
Transcription-Factor Regulon Activity in Disease Cholangiocytes
================================================================

Pipeline stage: COMPUTE → SAVE.

This version uses **hand-curated biliary/hepatocyte regulons** (compiled from
the cholangiocyte and liver-development literature) rather than fetching
CollecTRI online — so it's reproducible, has no external dependency, and
focuses on the regulators actually relevant to ductopenia.

Per cell we compute, for each TF, the mean log-normalized expression of its
known target gene set (an AUCell-style mean-rank approximation, simplified
via scanpy.tl.score_genes which already implements this with proper
control-gene normalization). The result is a per-cell TF activity matrix.

Outputs
-------
results/tf_activity_per_cell.tsv         (sample of cells)
results/tf_activity_per_donor.tsv        (mean per donor)
results/tf_activity_disease_summary.tsv  (mean per disease)
results/tf_activity_top_disease_TFs.tsv  (PSC vs normal delta)
figures/comparison/tf_activity_heatmap.png
figures/comparison/tf_activity_boxplot.png
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
from utils import ensure_dir

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR    = Path("data")
RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures") / "comparison"
LIVER_H5AD  = DATA_DIR / "liver" / "liver_cholangiocytes.h5ad"


# ── Hand-curated regulons ────────────────────────────────────────────────────
# Each TF has a curated target list compiled from cholangiocyte and liver-
# development literature. Targets are direct or strongly supported indirect.
# Lists deliberately overlap (e.g. KRT19 is target of multiple biliary TFs)
# because that's the actual biology.

REGULONS = {
    # ── Biliary lineage TFs (regenerative-target candidates) ────────────────
    "SOX9":     ["KRT19", "KRT7", "EPCAM", "SPP1", "FOXJ1", "CFTR",
                 "MUC1", "TACSTD2", "JAG1"],
    "HNF1B":    ["KRT19", "KRT7", "AQP1", "CFTR", "SLC4A2", "MUC1",
                 "PKHD1", "UMOD"],
    "ONECUT1":  ["HNF1B", "HNF4A", "FOXA2", "SOX9"],
    "HES1":     ["JAG1", "NOTCH2", "HEY1", "ASCL1"],
    "HEY1":     ["NOTCH2", "JAG1", "HES1"],
    "GATA6":    ["HNF1B", "HNF4A", "MUC1", "EPCAM"],
    "FOXA1":    ["ALB", "HNF4A", "HNF1B", "HNF1A"],
    "FOXA2":    ["ALB", "HNF4A", "HNF1B", "AFP"],

    # ── Hepatocyte lineage (reciprocal axis) ─────────────────────────────────
    "HNF4A":    ["ALB", "TTR", "APOB", "APOC3", "APOA2", "HNF1A", "AFP"],
    "PROX1":    ["HNF4A", "ALB", "TGFB1"],

    # ── Stress / IEG (not a regenerative target — included for orientation) ─
    "JUN":      ["FOS", "ATF3", "EGR1", "DUSP1"],
    "FOS":      ["JUN", "ATF3", "EGR1", "DUSP1"],

    # ── Pathway-level scores treated as pseudo-TFs ──────────────────────────
    "NOTCH_pathway":  ["HES1", "HEY1", "JAG1", "NOTCH2", "NOTCH3"],
    "WNT_pathway":    ["AXIN2", "LGR5", "RNF43", "ZNRF3", "WNT2"],
    "YAP_pathway":    ["CTGF", "CYR61", "AMOTL2", "ANKRD1", "TEAD1"],
    "BILIARY_id":     ["KRT19", "KRT7", "EPCAM", "SOX9", "HNF1B", "AQP1"],
    "HEPATOCYTE_id":  ["ALB", "HNF4A", "CYP3A4", "CYP2E1", "APOC3", "FABP1"],
}


# ── Score regulons via scanpy.tl.score_genes ─────────────────────────────────

def score_regulons(adata: ad.AnnData) -> pd.DataFrame:
    """Per-cell activity for each TF / pathway. Returns DataFrame (cell × TF)."""
    activities = {}
    for tf, targets in REGULONS.items():
        present = [g for g in targets if g in adata.var_names]
        if len(present) < 2:
            print(f"  [skip] {tf}: only {len(present)}/{len(targets)} targets present.")
            continue

        score_key = f"_score_{tf}"
        sc.tl.score_genes(
            adata, gene_list=present, score_name=score_key, random_state=42,
        )
        activities[tf] = adata.obs[score_key].values
        # Clean up — keep only the activity, drop the temporary obs column
        del adata.obs[score_key]
        print(f"  {tf:<18} scored on {len(present)}/{len(targets)} targets present")

    return pd.DataFrame(activities, index=adata.obs_names)


# ── Aggregate ─────────────────────────────────────────────────────────────────

def aggregate(activity: pd.DataFrame, obs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    df = activity.copy()
    df["__donor"]   = obs["donor_id"].astype(str).values \
        if "donor_id" in obs.columns else "unknown"
    df["__disease"] = obs["disease"].astype(str).values \
        if "disease" in obs.columns else "unknown"

    tfs = list(activity.columns)

    per_donor = (
        df.groupby(["__donor", "__disease"])[tfs].mean()
          .reset_index()
          .rename(columns={"__donor": "donor", "__disease": "disease"})
    )
    per_disease = (
        df.groupby("__disease")[tfs].mean()
          .reset_index()
          .rename(columns={"__disease": "disease"})
    )
    return per_donor, per_disease


# ── PSC vs normal contrast ───────────────────────────────────────────────────

def disease_contrast(per_disease: pd.DataFrame,
                     contrast: str = "primary sclerosing cholangitis",
                     reference: str = "normal") -> pd.DataFrame:
    if contrast not in per_disease["disease"].values:
        return pd.DataFrame()
    if reference not in per_disease["disease"].values:
        return pd.DataFrame()

    tfs = [c for c in per_disease.columns if c != "disease"]
    c = per_disease.set_index("disease").loc[contrast, tfs]
    r = per_disease.set_index("disease").loc[reference, tfs]
    delta = c - r
    out = pd.DataFrame({
        "TF":               tfs,
        "activity_normal":  r.values.astype(float),
        "activity_disease": c.values.astype(float),
        "delta_activity":   delta.values.astype(float),
    }).sort_values("delta_activity")
    out["direction"] = np.where(out["delta_activity"] > 0, "up", "down")
    return out


# ── Plots ─────────────────────────────────────────────────────────────────────

def heatmap(per_disease: pd.DataFrame, out: Path) -> None:
    if per_disease.empty:
        return
    tfs = [c for c in per_disease.columns if c != "disease"]
    df = per_disease.set_index("disease")[tfs]
    fig, ax = plt.subplots(figsize=(max(10, len(tfs) * 0.5),
                                    max(4, len(df) * 0.6)))
    sns.heatmap(df, cmap="RdBu_r", center=0, ax=ax,
                cbar_kws={"label": "Mean TF / pathway activity"},
                annot=True, fmt=".2f", annot_kws={"fontsize": 7})
    ax.set_title("Hand-curated TF / pathway activity by disease (cholangiocytes)")
    ax.set_xlabel("")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def boxplot_psc_vs_normal(activity: pd.DataFrame, obs: pd.DataFrame,
                          out: Path, regs_to_show: list[str]) -> None:
    if "disease" not in obs.columns:
        return
    df = activity[regs_to_show].copy()
    df["disease"] = obs["disease"].astype(str).values
    df = df[df["disease"].isin(["normal", "primary sclerosing cholangitis",
                                 "primary biliary cholangitis"])]
    if df.empty:
        return
    long = df.melt(id_vars="disease", var_name="TF", value_name="activity")

    fig, ax = plt.subplots(figsize=(max(8, len(regs_to_show) * 0.9), 5))
    palette = {
        "normal": "#16a085",
        "primary biliary cholangitis": "#e67e22",
        "primary sclerosing cholangitis": "#c0392b",
    }
    sns.boxplot(data=long, x="TF", y="activity", hue="disease",
                palette=palette, ax=ax, fliersize=0)
    ax.axhline(0, color="grey", linewidth=0.5, linestyle="--")
    ax.set_title("TF / pathway activity per cell — PSC and PBC vs normal")
    ax.set_ylabel("Activity (scanpy.tl.score_genes)")
    plt.xticks(rotation=30, ha="right")
    plt.legend(fontsize=8, loc="best", frameon=False)
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR)
    ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("TF / pathway activity — hand-curated regulons")
    print("=" * 60)

    if not LIVER_H5AD.exists():
        print(f"ERROR: {LIVER_H5AD} missing", file=sys.stderr); sys.exit(1)

    print("\n--- Loading ---")
    adata = sc.read_h5ad(LIVER_H5AD)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    print("\n--- Scoring regulons ---")
    activity = score_regulons(adata)
    print(f"  Activity matrix: {activity.shape}")

    print("\n--- Aggregating ---")
    per_donor, per_disease = aggregate(activity, adata.obs)
    per_donor.to_csv(RESULTS_DIR / "tf_activity_per_donor.tsv",
                     sep="\t", index=False)
    per_disease.to_csv(RESULTS_DIR / "tf_activity_disease_summary.tsv",
                       sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'tf_activity_per_donor.tsv'}")
    print(f"  → {RESULTS_DIR / 'tf_activity_disease_summary.tsv'}")

    print("\n--- PSC vs normal contrast ---")
    psc = disease_contrast(per_disease,
                           contrast="primary sclerosing cholangitis",
                           reference="normal")
    if not psc.empty:
        psc.to_csv(RESULTS_DIR / "tf_activity_top_disease_TFs.tsv",
                   sep="\t", index=False)
        print(psc.to_string(index=False))

    print("\n--- PBC vs normal contrast ---")
    pbc = disease_contrast(per_disease,
                           contrast="primary biliary cholangitis",
                           reference="normal")
    if not pbc.empty:
        print(pbc.to_string(index=False))

    print("\n--- Plots ---")
    heatmap(per_disease, FIGURES_DIR / "tf_activity_heatmap.png")
    biliary_focus = [t for t in ["SOX9", "HNF1B", "ONECUT1", "HES1", "HEY1",
                                  "GATA6", "FOXA1", "FOXA2", "HNF4A", "PROX1",
                                  "NOTCH_pathway", "WNT_pathway", "YAP_pathway",
                                  "BILIARY_id", "HEPATOCYTE_id"]
                     if t in activity.columns]
    boxplot_psc_vs_normal(activity, adata.obs,
                          FIGURES_DIR / "tf_activity_boxplot.png",
                          biliary_focus)

    print("\nDone.")
    print("Most regenerative-target-relevant TFs are biliary-lineage:")
    print("  SOX9, HNF1B, ONECUT1, HES1, HEY1, GATA6, FOXA1/2.")
    print("Look for which ones drop in PSC vs normal — those are the")
    print("candidate targets for restoring cholangiocyte identity.")


if __name__ == "__main__":
    main()
