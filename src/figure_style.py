"""
Shared matplotlib style for the three preprint figures.

Pulls all the cross-figure look-and-feel into one place so panel labels,
fonts, colors and figure sizes are consistent across Figures 1–3.
"""

import logging
import string
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# Silence "font not found" warnings — fallback chain handles missing fonts.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


# ── Color palette ─────────────────────────────────────────────────────────────
PALETTE = {
    "normal":                          "#16a085",
    "primary sclerosing cholangitis":  "#c0392b",
    "primary biliary cholangitis":     "#e67e22",
    "intestinal failure-associated":   "#8e44ad",
    "other":                           "#7f8c8d",
    # Categorical
    "snRNA":                           "#8e44ad",
    "scRNA":                           "#16a085",
    "up":                              "#c0392b",
    "down":                            "#2980b9",
    "neutral":                         "#7f8c8d",
    # Tertile
    "high":                            "#2980b9",
    "low":                             "#c0392b",
    # Categories used in result note
    "A_lost_in_SOX9low":               "#2980b9",
    "B_gained_in_SOX9low":             "#c0392b",
    "ambiguous":                       "#7f8c8d",
    "underpowered":                    "#bdc3c7",
    "A_tracks_SOX9_positively":        "#2980b9",
    "B_tracks_SOX9_negatively":        "#c0392b",
}


# ── Apply consistent rcParams ─────────────────────────────────────────────────
def apply_style() -> None:
    mpl.rcParams.update({
        "font.family":         "sans-serif",
        "font.sans-serif":     ["Helvetica Neue", "Helvetica", "Arial",
                                "DejaVu Sans"],
        "font.size":           10,
        "axes.titlesize":      11,
        "axes.labelsize":      10,
        "xtick.labelsize":     9,
        "ytick.labelsize":     9,
        "legend.fontsize":     9,
        "figure.titlesize":    12,
        "axes.spines.top":     False,
        "axes.spines.right":   False,
        "axes.linewidth":      0.9,
        "xtick.major.width":   0.9,
        "ytick.major.width":   0.9,
        "xtick.major.size":    3.5,
        "ytick.major.size":    3.5,
        "axes.grid":           False,
        "savefig.dpi":         300,
        "figure.dpi":          120,
        "pdf.fonttype":        42,    # Editable text in PDF
        "ps.fonttype":         42,
    })


def panel_label(ax, letter: str, xoff: float = -0.18, yoff: float = 1.05) -> None:
    """Add a bold panel label (A, B, C, D) at the upper-left of an axis."""
    ax.text(xoff, yoff, letter,
            transform=ax.transAxes,
            fontsize=13, fontweight="bold",
            ha="left", va="top")


def save_fig(fig, out_base: Path) -> None:
    """Save as both PNG (high-DPI) and PDF (editable text)."""
    out_base = Path(out_base)
    out_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_base.with_suffix(".png"), bbox_inches="tight", dpi=300)
    fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_base.with_suffix('.png')}")
    print(f"  → {out_base.with_suffix('.pdf')}")
