"""
Per-cell SOX9 dose-response within PSC donors.

Question
--------
Donor-level Δ (SOX9-high − SOX9-low pathway score) nominated candidate
upstream pathways (oxidative stress, IL6, YAP loss, IFN priming). The
donor-level test has n = 7 PSC donors and is statistically conservative.
A sharper test: does each candidate pathway score correlate with per-cell
SOX9 *within individual donors*? That uses cell-level statistical power
(thousands per donor) while controlling for donor identity automatically
(each donor is one ρ).

Design
------
1. Subset cholangiocytes to PSC.
2. Score curated pathways per cell (same gene sets as
   sox9_mechanism_discovery.py).
3. For each donor with ≥ 30 cells, compute Spearman ρ(SOX9, pathway_score)
   *within that donor*. Donor-level ρ is the test statistic.
4. Aggregate across donors:
     - mean ρ across donors
     - fraction of donors with same direction
     - one-sample t-test of per-donor ρs against zero, BH-corrected

For a pathway truly upstream of SOX9 (or downstream — direction not
inferrable), we expect a consistent within-donor correlation regardless
of donor identity. For donor-effect-only signals, the within-donor ρ
should average to ~0.

Outputs
-------
results/sox9_per_cell_dose_response.tsv            per-donor ρ per pathway
results/sox9_per_cell_dose_response_summary.tsv    cross-donor summary
figures/sox9/sox9_per_cell_dose_response.png       boxplot of per-donor ρ

Usage
-----
    python src/sox9_per_cell_dose_response.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy.stats import spearmanr, ttest_1samp

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression
from iad_score import LIVER_H5AD as DEFAULT_CHOL
from sox9_mechanism_discovery import PATHWAYS, score_pathways, load_psc

CHOL_PATH    = DEFAULT_CHOL
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "sox9"

GENE = "SOX9"
MIN_CELLS_PER_DONOR = 30
DIRECTION_THRESHOLD = 0.05    # |ρ| > this counts as having direction


# ── Per-donor within-donor Spearman ────────────────────────────────────────────

def per_donor_correlations(adata, scores: dict) -> pd.DataFrame:
    sox9 = get_gene_expression(adata, GENE)
    obs = adata.obs[["donor_id", "tech"]].copy()
    obs.columns = ["donor", "tech"]
    obs["SOX9"] = sox9
    for name, vec in scores.items():
        obs[f"S_{name}"] = vec

    rows = []
    for donor, dgrp in obs.groupby("donor", observed=True):
        if len(dgrp) < MIN_CELLS_PER_DONOR:
            continue
        for name in scores:
            col = f"S_{name}"
            v = dgrp[col].values
            s = dgrp["SOX9"].values
            ok = ~np.isnan(v) & ~np.isnan(s)
            if ok.sum() < 30:
                continue
            # ρ(SOX9, pathway) within this donor
            rho, p = spearmanr(s[ok], v[ok])
            rows.append({
                "donor":   donor,
                "tech":    dgrp["tech"].iloc[0],
                "pathway": name,
                "n_cells": int(ok.sum()),
                "rho":     float(rho),
                "p":       float(p),
            })
    return pd.DataFrame(rows)


# ── Aggregate ──────────────────────────────────────────────────────────────────

def cross_donor_summary(per_donor: pd.DataFrame) -> pd.DataFrame:
    if per_donor.empty:
        return pd.DataFrame()
    rows = []
    for pw, sub in per_donor.groupby("pathway"):
        rhos = sub["rho"].dropna().values
        if len(rhos) < 3 or len(set(rhos)) < 2:
            t, p = float("nan"), float("nan")
        else:
            t, p = ttest_1samp(rhos, 0.0)
        rows.append({
            "pathway":           pw,
            "n_donors":          int(len(rhos)),
            "mean_rho_per_cell": float(np.mean(rhos)),
            "median_rho":        float(np.median(rhos)),
            "min_rho":           float(np.min(rhos)),
            "max_rho":           float(np.max(rhos)),
            "frac_pos_donors":   float((rhos > DIRECTION_THRESHOLD).mean()),
            "frac_neg_donors":   float((rhos < -DIRECTION_THRESHOLD).mean()),
            "t_stat":            float(t),
            "p":                 float(p),
        })
    df = pd.DataFrame(rows).sort_values("p")
    df["rank"] = np.arange(1, len(df) + 1)
    df["q_BH"] = (df["p"] * len(df) / df["rank"]).clip(upper=1.0)
    df = df.sort_values("mean_rho_per_cell")
    return df


def label_category(row) -> str:
    """A/B classification using cell-level evidence.

    A: pathway tracks SOX9 positively (lost in SOX9-low) — concordant
       positive within-donor ρ across ≥ 70 % of donors AND mean ρ > 0.05
    B: pathway tracks SOX9 negatively (gained in SOX9-low) — concordant
       negative ρ across ≥ 70 % of donors AND mean ρ < -0.05
    """
    if np.isnan(row["mean_rho_per_cell"]) or row["n_donors"] < 4:
        return "underpowered"
    if (row["mean_rho_per_cell"] > 0.05 and
        row["frac_pos_donors"] >= 5/7):
        return "A_tracks_SOX9_positively"
    if (row["mean_rho_per_cell"] < -0.05 and
        row["frac_neg_donors"] >= 5/7):
        return "B_tracks_SOX9_negatively"
    return "ambiguous"


# ── Plotting ───────────────────────────────────────────────────────────────────

def plot_per_donor_rho(per_donor: pd.DataFrame, summary: pd.DataFrame,
                       out: Path) -> None:
    if per_donor.empty or summary.empty:
        return
    pw_order = summary.sort_values("mean_rho_per_cell")["pathway"].tolist()
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(pw_order) + 2))
    palette = {
        "A_tracks_SOX9_positively": "#2980b9",
        "B_tracks_SOX9_negatively": "#c0392b",
        "ambiguous":                "#7f8c8d",
        "underpowered":             "#bdc3c7",
    }
    summary["cat"] = summary.apply(label_category, axis=1)
    cat_lookup = summary.set_index("pathway")["cat"].to_dict()

    box_data = [per_donor[per_donor["pathway"] == pw]["rho"].values
                for pw in pw_order]
    colors = [palette.get(cat_lookup.get(pw, ""), "#7f8c8d") for pw in pw_order]
    bp = ax.boxplot(box_data, vert=False, patch_artist=True, widths=0.5,
                    flierprops={"marker": "."})
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.7)
        patch.set_edgecolor("black")
    # overlay per-donor points
    for i, pw in enumerate(pw_order):
        vals = per_donor[per_donor["pathway"] == pw]["rho"].values
        jitter = (np.random.default_rng(42).random(len(vals)) - 0.5) * 0.2
        ax.scatter(vals, [i + 1] * len(vals) + jitter,
                   s=14, c="black", alpha=0.6)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_yticks(range(1, len(pw_order) + 1))
    ax.set_yticklabels(pw_order, fontsize=9)
    ax.set_xlabel("Within-donor Spearman ρ (per-cell SOX9 vs pathway score)")
    ax.set_title("Per-cell SOX9 dose-response within PSC donors\n"
                 "(each black point = one donor's within-donor ρ)")
    handles = [
        plt.Line2D([0],[0], marker="s", color="w",
                   markerfacecolor=palette["A_tracks_SOX9_positively"],
                   markersize=10, label="A: tracks SOX9 + (lost in SOX9-low)"),
        plt.Line2D([0],[0], marker="s", color="w",
                   markerfacecolor=palette["B_tracks_SOX9_negatively"],
                   markersize=10, label="B: tracks SOX9 − (gained in SOX9-low)"),
        plt.Line2D([0],[0], marker="s", color="w",
                   markerfacecolor=palette["ambiguous"],
                   markersize=10, label="ambiguous / small effect"),
    ]
    ax.legend(handles=handles, fontsize=8, frameon=False, loc="best")
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Per-cell SOX9 dose-response within PSC donors")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    print(f"\n[1/4] Loading + subsetting to PSC …")
    adata = sc.read_h5ad(CHOL_PATH)
    psc = load_psc(adata)
    print(f"  PSC cells: {psc.n_obs:,}")

    print(f"\n[2/4] Scoring pathways per cell …")
    scores = score_pathways(psc)
    print(f"  {len(scores)} pathways scored")

    print(f"\n[3/4] Within-donor Spearman ρ(SOX9, pathway) per donor …")
    per_donor = per_donor_correlations(psc, scores)
    print(f"  {per_donor['donor'].nunique() if not per_donor.empty else 0} donors "
          f"with ≥ {MIN_CELLS_PER_DONOR} cells")
    per_donor.to_csv(RESULTS_DIR / "sox9_per_cell_dose_response.tsv",
                     sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'sox9_per_cell_dose_response.tsv'}")

    print(f"\n[4/4] Cross-donor summary …")
    summary = cross_donor_summary(per_donor)
    summary["category"] = summary.apply(label_category, axis=1)
    summary.to_csv(RESULTS_DIR / "sox9_per_cell_dose_response_summary.tsv",
                   sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'sox9_per_cell_dose_response_summary.tsv'}")

    print("\n--- Summary (sorted by mean within-donor ρ) ---")
    print(summary[["pathway", "n_donors", "mean_rho_per_cell",
                     "min_rho", "max_rho",
                     "frac_pos_donors", "frac_neg_donors",
                     "p", "q_BH", "category"]]
          .to_string(index=False, float_format="%.3f"))

    # Plot
    print("\n--- Plot ---")
    plot_per_donor_rho(per_donor, summary,
                       FIGURES_DIR / "sox9_per_cell_dose_response.png")

    # Headline
    print("\n--- Headline (cell-level evidence per pathway) ---")
    A = summary[summary["category"] == "A_tracks_SOX9_positively"]
    B = summary[summary["category"] == "B_tracks_SOX9_negatively"]
    print(f"  A. TRACKS SOX9+ (within-donor positive, lost in SOX9-low): "
          f"{len(A)} pathways")
    for _, r in A.iterrows():
        print(f"     {r['pathway']:<26}  mean ρ={r['mean_rho_per_cell']:+.3f}  "
              f"range [{r['min_rho']:+.2f}, {r['max_rho']:+.2f}]   "
              f"q={r['q_BH']:.2g}  ({r['n_donors']} donors)")
    print(f"  B. TRACKS SOX9- (within-donor negative, gained in SOX9-low): "
          f"{len(B)} pathways")
    for _, r in B.iterrows():
        print(f"     {r['pathway']:<26}  mean ρ={r['mean_rho_per_cell']:+.3f}  "
              f"range [{r['min_rho']:+.2f}, {r['max_rho']:+.2f}]   "
              f"q={r['q_BH']:.2g}  ({r['n_donors']} donors)")


if __name__ == "__main__":
    main()
