"""
Figure 3 — Candidate upstream pathways: donor-level signal vs cell-level
dose-response.

Four panels:
A. Donor-level Δ heatmap (pathway × donor) showing SOX9-high − SOX9-low.
B. Cell-level within-donor ρ boxplot for the same pathways.
C. Side-by-side comparison — donor-level mean Δ vs cell-level mean ρ,
   highlighting which donor candidates collapse at cell-level.
D. Summary verdict bar — A_lost / B_gained / ambiguous classification.

Outputs: figures/preprint/figure3_pathways.png + .pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import apply_style, panel_label, save_fig, PALETTE

RESULTS_DIR = Path("results")
OUT_BASE    = Path("figures") / "preprint" / "figure3_pathways"


def main():
    apply_style()
    fig = plt.figure(figsize=(12.5, 8.5))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0],
                            width_ratios=[1.1, 1.0],
                            hspace=0.42, wspace=0.32)
    axA = fig.add_subplot(gs[0, 0])
    axB = fig.add_subplot(gs[0, 1])
    axC = fig.add_subplot(gs[1, 0])
    axD = fig.add_subplot(gs[1, 1])

    # ── Panel A: per-donor Δ heatmap ──────────────────────────────────────────
    donor_long = pd.read_csv(
        RESULTS_DIR / "sox9_pathway_donor_contrast_psc.tsv", sep="\t")
    summary_donor = pd.read_csv(
        RESULTS_DIR / "sox9_pathway_summary_psc.tsv", sep="\t")
    summary_donor = summary_donor.sort_values("mean_delta")
    pw_order_donor = summary_donor["pathway"].tolist()

    pivot = (donor_long.pivot_table(index="pathway", columns="donor",
                                       values="delta_HminusL", aggfunc="mean")
                        .reindex(pw_order_donor))
    vmax = float(np.nanmax(np.abs(pivot.values))) if pivot.size else 0.5
    im = axA.imshow(pivot.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax,
                     aspect="auto")
    axA.set_yticks(range(len(pivot)))
    axA.set_yticklabels(pivot.index, fontsize=8)
    axA.set_xticks(range(len(pivot.columns)))
    axA.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
    axA.set_xlabel("Donor")
    axA.set_title("Donor-level pathway Δ (SOX9-high − SOX9-low)")
    cbar = plt.colorbar(im, ax=axA, fraction=0.035, pad=0.02)
    cbar.set_label("Δ score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    panel_label(axA, "A")

    # ── Panel B: cell-level within-donor ρ boxplot ────────────────────────────
    cell_long = pd.read_csv(
        RESULTS_DIR / "sox9_per_cell_dose_response.tsv", sep="\t")
    cell_summary = pd.read_csv(
        RESULTS_DIR / "sox9_per_cell_dose_response_summary.tsv", sep="\t")
    cell_summary = cell_summary.sort_values("mean_rho_per_cell")
    pw_order_cell = cell_summary["pathway"].tolist()
    box_data = [cell_long[cell_long["pathway"] == pw]["rho"].values
                for pw in pw_order_cell]
    colors = []
    for _, r in cell_summary.iterrows():
        cat = r["category"]
        if cat == "A_tracks_SOX9_positively":
            colors.append(PALETTE["A_tracks_SOX9_positively"])
        elif cat == "B_tracks_SOX9_negatively":
            colors.append(PALETTE["B_tracks_SOX9_negatively"])
        else:
            colors.append(PALETTE["ambiguous"])
    bp = axB.boxplot(box_data, vert=False, patch_artist=True, widths=0.55,
                       flierprops={"marker": ".", "markersize": 4},
                       medianprops={"color": "white", "linewidth": 1.0})
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.78)
        patch.set_edgecolor("black"); patch.set_linewidth(0.4)
    rng = np.random.default_rng(7)
    for i, pw in enumerate(pw_order_cell):
        vals = cell_long[cell_long["pathway"] == pw]["rho"].values
        jitter = (rng.random(len(vals)) - 0.5) * 0.22
        axB.scatter(vals, [i + 1] * len(vals) + jitter, s=10,
                     c="black", alpha=0.55, zorder=3)
    axB.axvline(0, color="black", linewidth=0.5)
    axB.set_yticks(range(1, len(pw_order_cell) + 1))
    axB.set_yticklabels(pw_order_cell, fontsize=8)
    axB.set_xlabel("Within-donor ρ (per-cell SOX9 vs pathway score)")
    axB.set_title("Cell-level dose-response within donors")
    # Mark FDR-significant pathways
    for i, (_, r) in enumerate(cell_summary.iterrows()):
        if r["q_BH"] < 0.05:
            axB.text(axB.get_xlim()[1] * 0.95, i + 1, "**",
                      fontsize=12, fontweight="bold",
                      ha="right", va="center")
    panel_label(axB, "B")

    # ── Panel C: donor-level Δ vs cell-level ρ ────────────────────────────────
    join = summary_donor[["pathway", "mean_delta",
                            "frac_pos_donors", "frac_neg_donors"]].merge(
        cell_summary[["pathway", "mean_rho_per_cell", "q_BH"]],
        on="pathway", how="inner")
    # Color by quadrant
    def quad_color(r):
        a = r["mean_delta"]; b = r["mean_rho_per_cell"]
        if abs(a) < 0.05 and abs(b) < 0.05:
            return PALETTE["ambiguous"]
        if a > 0 and b > 0:   return PALETTE["A_tracks_SOX9_positively"]
        if a < 0 and b < 0:   return PALETTE["B_tracks_SOX9_negatively"]
        if a < 0 and b > 0:   return "#d4a017"  # donor down, cell up
        if a > 0 and b < 0:   return "#7d3c98"
        return PALETTE["ambiguous"]
    join["color"] = join.apply(quad_color, axis=1)
    axC.axhline(0, color="black", linewidth=0.5)
    axC.axvline(0, color="black", linewidth=0.5)
    axC.scatter(join["mean_delta"], join["mean_rho_per_cell"],
                 s=100, c=join["color"], alpha=0.85, edgecolor="white",
                 linewidth=0.6)
    # Label pathways with non-trivial movement in either axis
    for _, r in join.iterrows():
        if abs(r["mean_delta"]) >= 0.08 or abs(r["mean_rho_per_cell"]) >= 0.08:
            axC.annotate(r["pathway"], (r["mean_delta"],
                                          r["mean_rho_per_cell"]),
                          fontsize=7, alpha=0.85,
                          xytext=(4, 2), textcoords="offset points")
    axC.set_xlabel("Donor-level Δ (SOX9-high − SOX9-low)")
    axC.set_ylabel("Cell-level mean ρ (within donors)")
    axC.set_title("Donor-level vs cell-level evidence per pathway\n"
                   "Diagonal alignment = consistent; off-axis = donor-only signal")
    panel_label(axC, "C")

    # ── Panel D: verdict summary bar ──────────────────────────────────────────
    # Count categories from both analyses
    donor_cats = summary_donor["category"].value_counts()
    cell_cats  = cell_summary["category"].value_counts()
    cat_order = ["A_lost_in_SOX9low", "B_gained_in_SOX9low",
                  "ambiguous", "underpowered"]
    cat_order_cell = ["A_tracks_SOX9_positively",
                       "B_tracks_SOX9_negatively",
                       "ambiguous", "underpowered"]
    cat_labels = ["A: lost\n(maintenance)", "B: gained\n(disease/stress)",
                   "ambiguous", "underpowered"]
    donor_counts = [donor_cats.get(c, 0) for c in cat_order]
    cell_counts  = [cell_cats.get(c, 0)  for c in cat_order_cell]
    x = np.arange(len(cat_labels))
    bw = 0.36
    cols = [PALETTE["A_tracks_SOX9_positively"],
             PALETTE["B_tracks_SOX9_negatively"],
             PALETTE["ambiguous"], PALETTE["underpowered"]]
    axD.bar(x - bw/2, donor_counts, width=bw,
             color=cols, edgecolor="black", linewidth=0.5, alpha=0.85,
             label="donor-level Δ")
    axD.bar(x + bw/2, cell_counts, width=bw,
             color=cols, edgecolor="black", linewidth=0.5, alpha=0.40,
             hatch="//", label="cell-level ρ")
    for i, (d, c) in enumerate(zip(donor_counts, cell_counts)):
        axD.text(x[i] - bw/2, d + 0.2, str(d), ha="center", fontsize=8)
        axD.text(x[i] + bw/2, c + 0.2, str(c), ha="center", fontsize=8)
    axD.set_xticks(x)
    axD.set_xticklabels(cat_labels, fontsize=8)
    axD.set_ylabel("# pathways")
    axD.set_title("Most donor-level candidates do not survive cell-level testing")
    axD.legend(fontsize=8, frameon=False, loc="upper right")
    axD.set_ylim(0, max(max(donor_counts), max(cell_counts)) + 2)
    panel_label(axD, "D")

    fig.suptitle("Figure 3 — Candidate upstream pathways: donor-level vs "
                  "cell-level evidence",
                  y=1.00, fontsize=12, fontweight="bold")
    save_fig(fig, OUT_BASE)


if __name__ == "__main__":
    main()
