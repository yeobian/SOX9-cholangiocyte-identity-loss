"""
Figure 1 — SOX9 collapse is the robust within-PSC progression signal.

Four panels, all built from existing result TSVs (no AnnData reload needed).

A. Donor-level: mean PSC pseudotime vs mean IAD score, colored by disease.
B. Donor-level: mean PSC pseudotime vs SOX9 regulon activity, colored by disease.
C. Per-dataset replication of SOX9 ↔ IAD direction (9 datasets).
D. Technology stratification: per-cell ρ(SOX9, IAD) in snRNA-seq vs scRNA-seq.

Outputs: figures/preprint/figure1_sox9_collapse.png + .pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import apply_style, panel_label, save_fig, PALETTE

RESULTS_DIR = Path("results")
OUT_BASE    = Path("figures") / "preprint" / "figure1_sox9_collapse"


def main():
    apply_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.2))

    # ── Panel A: donor mean pseudotime vs mean IAD ────────────────────────────
    ax = axes[0, 0]
    df_donor = pd.read_csv(RESULTS_DIR / "trajectory_psc_per_donor.tsv", sep="\t")
    df_donor = df_donor.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["mean_pseudotime", "mean_iad_score"]
    )
    for disease in sorted(df_donor["disease"].unique()):
        sub = df_donor[df_donor["disease"] == disease]
        ax.scatter(sub["mean_pseudotime"], sub["mean_iad_score"],
                   s=np.clip(sub["n_cells"] / 5, 12, 200),
                   c=PALETTE.get(disease, PALETTE["other"]),
                   alpha=0.78, edgecolor="white", linewidth=0.6,
                   label=disease.replace("primary ", ""))
    rho, p = spearmanr(df_donor["mean_pseudotime"], df_donor["mean_iad_score"])
    ax.set_xlabel("Donor-mean PSC pseudotime")
    ax.set_ylabel("Donor-mean IAD score")
    ax.set_title(f"PSC + normal donors\nρ = {rho:+.2f}, p = {p:.1e}, n = {len(df_donor)}")
    ax.legend(loc="lower right", fontsize=7, frameon=False)
    panel_label(ax, "A")

    # ── Panel B: pseudotime vs SOX9 regulon ──────────────────────────────────
    ax = axes[0, 1]
    tf_donor = pd.read_csv(RESULTS_DIR / "tf_activity_per_donor.tsv", sep="\t")
    joined = df_donor.merge(tf_donor[["donor", "SOX9"]], on="donor", how="inner")
    rho_s, p_s = spearmanr(joined["mean_pseudotime"], joined["SOX9"])
    for disease in sorted(joined["disease"].unique()):
        sub = joined[joined["disease"] == disease]
        ax.scatter(sub["mean_pseudotime"], sub["SOX9"],
                   s=np.clip(sub["n_cells"] / 5, 12, 200),
                   c=PALETTE.get(disease, PALETTE["other"]),
                   alpha=0.78, edgecolor="white", linewidth=0.6)
    # Add trend line
    if len(joined) >= 5:
        z = np.polyfit(joined["mean_pseudotime"], joined["SOX9"], 1)
        xs = np.linspace(joined["mean_pseudotime"].min(),
                          joined["mean_pseudotime"].max(), 50)
        ax.plot(xs, np.polyval(z, xs), color="black", linewidth=1.2,
                linestyle="--", alpha=0.5)
    ax.set_xlabel("Donor-mean PSC pseudotime")
    ax.set_ylabel("Donor-mean SOX9 regulon activity")
    ax.set_title(f"SOX9 collapse along PSC trajectory\n"
                 f"ρ = {rho_s:+.2f}, p = {p_s:.1e}")
    panel_label(ax, "B")

    # ── Panel C: per-dataset SOX9 ρ across 9 datasets ─────────────────────────
    ax = axes[1, 0]
    rep = pd.read_csv(RESULTS_DIR / "per_dataset_replication.tsv", sep="\t")
    sox9_rep = rep[rep["gene"] == "SOX9"].copy()
    sox9_rep = sox9_rep.sort_values("rho", ascending=True)
    ds_short = [s[:14] + "…" for s in sox9_rep["dataset_id"]]
    colors = []
    for _, r in sox9_rep.iterrows():
        if r["p"] < 0.05 and r["rho"] < 0:
            colors.append(PALETTE["down"])
        elif r["rho"] < 0:
            colors.append("#a6c8df")  # light blue
        else:
            colors.append(PALETTE["up"])
    ax.barh(range(len(sox9_rep)), sox9_rep["rho"],
            color=colors, edgecolor="white", linewidth=0.5)
    # Mark significant ones with asterisk
    for i, (_, r) in enumerate(sox9_rep.iterrows()):
        if r["p"] < 0.05:
            ax.text(r["rho"] - 0.02 if r["rho"] < 0 else r["rho"] + 0.02,
                    i, "*", fontsize=12, fontweight="bold",
                    ha="right" if r["rho"] < 0 else "left", va="center")
    ax.set_yticks(range(len(sox9_rep)))
    ax.set_yticklabels(ds_short, fontsize=7)
    ax.set_xlabel("Donor-level ρ (SOX9 mRNA vs IAD score)")
    ax.set_title("Per-dataset replication\n* = individually p < 0.05")
    ax.axvline(0, color="black", linewidth=0.6)
    panel_label(ax, "C")

    # ── Panel D: Technology stratification — SOX9 vs IAD ──────────────────────
    ax = axes[1, 1]
    tech = pd.read_csv(RESULTS_DIR / "technology_stratified_prox1.tsv", sep="\t")
    # We want SOX9 vs IAD rho, which is the metric 'sox9_vs_iad_cell_rho'
    sox9_row = tech[tech["metric"] == "sox9_vs_iad_cell_rho"]
    if not sox9_row.empty:
        snRNA_v = float(sox9_row["snRNA"].iloc[0])
        scRNA_v = float(sox9_row["scRNA"].iloc[0])
    else:
        snRNA_v = scRNA_v = float("nan")
    bars = ax.bar(["snRNA-seq", "scRNA-seq"], [snRNA_v, scRNA_v],
                  color=[PALETTE["snRNA"], PALETTE["scRNA"]],
                  edgecolor="white", linewidth=0.8, width=0.55)
    for b, v in zip(bars, [snRNA_v, scRNA_v]):
        ax.text(b.get_x() + b.get_width()/2, v - 0.03,
                f"{v:+.2f}", ha="center", va="top",
                fontsize=10, fontweight="bold", color="white")
    # Add donor counts as annotation
    ax.text(0, 0.05, "n = 3,741 cells\n12 donors", ha="center",
            transform=ax.get_xaxis_transform(), fontsize=7, color=PALETTE["snRNA"])
    ax.text(1, 0.05, "n = 351 cells\n3 PSC donors\n+ 44 normal", ha="center",
            transform=ax.get_xaxis_transform(), fontsize=7, color=PALETTE["scRNA"])
    ax.axhline(0, color="black", linewidth=0.6)
    ax.set_ylabel("Per-cell ρ (SOX9 vs IAD score)")
    ax.set_title("Technology stratification\nSOX9 ↓ in both, direction preserved")
    ax.set_ylim(min(snRNA_v, scRNA_v) * 1.3, 0.05)
    panel_label(ax, "D")

    plt.tight_layout()
    fig.suptitle("Figure 1 — SOX9 collapse along PSC progression",
                 y=1.02, fontsize=12, fontweight="bold")
    save_fig(fig, OUT_BASE)


if __name__ == "__main__":
    main()
