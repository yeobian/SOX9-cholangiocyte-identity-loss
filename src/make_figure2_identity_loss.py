"""
Figure 2 — SOX9-low is a coordinated biliary-identity-loss state, and PROX1
marks a separate cholangiocyte sub-state.

Three panels:

A. Within-donor Spearman ρ(SOX9, pathway score) across 8 PSC donors,
   highlighting that biliary_identity is the only pathway robustly tracking
   SOX9 at cell level.
B. Donor-level scatter: mean PROX1 vs mean SOX9 across all donors, by
   disease. Shows PROX1 and SOX9 vary largely independently.
C. PROX1-high vs PROX1-low cholangiocyte markers (SOX9, BIL_id, KRT19+,
   KRT7+, EPCAM+, %SOX9-high) across ALL / Normal / PSC slices.

Outputs: figures/preprint/figure2_identity_loss.png + .pdf
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from figure_style import apply_style, panel_label, save_fig, PALETTE

RESULTS_DIR = Path("results")
OUT_BASE    = Path("figures") / "preprint" / "figure2_identity_loss"


def main():
    apply_style()
    fig = plt.figure(figsize=(11.5, 5.2))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 0.9, 1.4])
    axA = fig.add_subplot(gs[0])
    axB = fig.add_subplot(gs[1])
    axC = fig.add_subplot(gs[2])

    # ── Panel A: within-donor ρ(SOX9, pathway) ────────────────────────────────
    per_donor = pd.read_csv(RESULTS_DIR / "sox9_per_cell_dose_response.tsv",
                             sep="\t")
    summary  = pd.read_csv(RESULTS_DIR / "sox9_per_cell_dose_response_summary.tsv",
                             sep="\t")
    # Sort by mean ρ
    summary = summary.sort_values("mean_rho_per_cell")
    pw_order = summary["pathway"].tolist()
    box_data = [per_donor[per_donor["pathway"] == pw]["rho"].values
                for pw in pw_order]
    colors = []
    for _, r in summary.iterrows():
        cat = r["category"]
        if cat == "A_tracks_SOX9_positively":
            colors.append(PALETTE["A_tracks_SOX9_positively"])
        elif cat == "B_tracks_SOX9_negatively":
            colors.append(PALETTE["B_tracks_SOX9_negatively"])
        else:
            colors.append(PALETTE["ambiguous"])

    bp = axA.boxplot(box_data, vert=False, patch_artist=True, widths=0.55,
                      flierprops={"marker": "."}, medianprops={"color": "white",
                                                                "linewidth": 1.2})
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.75)
        patch.set_edgecolor("black"); patch.set_linewidth(0.5)
    # Per-donor points
    rng = np.random.default_rng(7)
    for i, pw in enumerate(pw_order):
        vals = per_donor[per_donor["pathway"] == pw]["rho"].values
        jitter = (rng.random(len(vals)) - 0.5) * 0.22
        axA.scatter(vals, [i + 1] * len(vals) + jitter, s=12,
                    c="black", alpha=0.55, zorder=3)
    axA.axvline(0, color="black", linewidth=0.5, linestyle="-")
    axA.set_yticks(range(1, len(pw_order) + 1))
    axA.set_yticklabels(pw_order, fontsize=7.5)
    axA.set_xlabel("Within-donor ρ (per-cell SOX9 vs pathway score)")
    axA.set_title("Cell-level dose-response\nacross 8 PSC donors")
    # Mark FDR-significant pathways
    for i, (_, r) in enumerate(summary.iterrows()):
        if r["q_BH"] < 0.05:
            axA.text(0.97, (i + 1) / (len(pw_order) + 1),
                     "**", transform=axA.transAxes, fontsize=10,
                     fontweight="bold", ha="right", va="center")
    panel_label(axA, "A")

    # ── Panel B: donor mean PROX1 vs mean SOX9 ────────────────────────────────
    pscxchol = pd.read_csv(RESULTS_DIR / "prox1_sox9_per_donor.tsv", sep="\t")
    pscxchol = pscxchol[pscxchol["n_cells"] >= 20]   # match abundance threshold
    for disease in ["normal", "primary biliary cholangitis",
                     "primary sclerosing cholangitis"]:
        sub = pscxchol[pscxchol["disease"] == disease]
        if sub.empty: continue
        axB.scatter(sub["mean_PROX1"], sub["mean_SOX9"],
                    s=np.clip(sub["n_cells"] / 5, 12, 200),
                    c=PALETTE.get(disease, PALETTE["other"]),
                    alpha=0.78, edgecolor="white", linewidth=0.6,
                    label=disease.replace("primary ", ""))
    axB.axvline(0, color="grey", linewidth=0.4)
    axB.axhline(0, color="grey", linewidth=0.4)
    axB.set_xlabel("Donor mean PROX1 expression")
    axB.set_ylabel("Donor mean SOX9 expression")
    axB.set_title("PROX1 ↔ SOX9 at donor level\n(largely independent)")
    axB.legend(fontsize=7, frameon=False, loc="best")
    panel_label(axB, "B")

    # ── Panel C: PROX1-high vs PROX1-low marker positivity %s ────────────────
    # Show only percentage metrics on a common 0-100 scale, as grouped bars
    # of high vs low pairs for the ALL and PSC slices.
    marker = pd.read_csv(RESULTS_DIR / "prox1_sox9_marker_summary.tsv", sep="\t")
    pct_metrics = [
        ("%KRT19+",                "KRT19+"),
        ("%KRT7+",                 "KRT7+"),
        ("%EPCAM+",                "EPCAM+"),
        ("%SOX9-high (top quartile)", "SOX9-high"),
    ]
    slc_use = ["ALL", "PSC"]
    sub = marker[marker["slice"].isin(slc_use)
                  & marker["metric"].isin([m[0] for m in pct_metrics])].copy()
    if sub.empty:
        axC.text(0.5, 0.5, "no data", ha="center", va="center",
                  transform=axC.transAxes)
    else:
        metric_order  = [m[0] for m in pct_metrics]
        metric_labels = [m[1] for m in pct_metrics]
        x = np.arange(len(metric_order))
        bw = 0.19
        # Within ALL slice, plot high then low; within PSC slice, plot high then low
        slice_group_colors = {
            ("ALL", "high"):    "#8e44ad",
            ("ALL", "low"):     "#bdc3c7",
            ("PSC", "high"):    "#c0392b",
            ("PSC", "low"):     "#e8b5b0",
        }
        offsets = {"ALL_high": -1.5, "ALL_low": -0.5,
                    "PSC_high":  0.5, "PSC_low":  1.5}
        for slc in slc_use:
            ssub = sub[sub["slice"] == slc].set_index("metric").reindex(metric_order)
            for grp in ("high", "low"):
                col = "high_mean" if grp == "high" else "low_mean"
                vals = ssub[col].values
                offset = offsets[f"{slc}_{grp}"] * bw
                axC.bar(x + offset, vals, width=bw,
                        color=slice_group_colors[(slc, grp)],
                        edgecolor="white", linewidth=0.4,
                        label=f"{slc} · PROX1-{grp}")
        axC.set_xticks(x)
        axC.set_xticklabels(metric_labels, fontsize=9)
        axC.set_ylabel("% positive cells")
        axC.set_ylim(0, max(60, sub[["high_mean", "low_mean"]].values.max() + 8))
        axC.set_title("PROX1-high cells have LOWER biliary markers\n"
                       "(ALL cholangiocytes and PSC slice)")
        axC.legend(fontsize=7, frameon=False, loc="upper right",
                    ncol=2, columnspacing=0.6, handletextpad=0.4)
    panel_label(axC, "C")

    plt.tight_layout()
    fig.suptitle("Figure 2 — SOX9 captures coordinated biliary identity; "
                 "PROX1 marks a separate cholangiocyte sub-state",
                 y=1.04, fontsize=12, fontweight="bold")
    save_fig(fig, OUT_BASE)


if __name__ == "__main__":
    main()
