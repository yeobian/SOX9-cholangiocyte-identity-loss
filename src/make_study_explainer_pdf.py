"""
Polished five-page scientific brief — implements the design package in
results/STUDY_BRIEF_DESIGN_PACKAGE.md.

Pages
-----
1. Cover — title, subtitle, disclaimer pill, abstract, table of contents
2. Figure 1 — SOX9 collapse along PSC progression
3. Figure 2 — SOX9 captures a coordinated biliary-identity program
4. Figure 3 — Donor-level vs cell-level evidence for candidate upstream pathways
5. Limitations and next experimental steps

Typography: two-family system (sans-serif for display + chrome, serif for body)
Palette:    6 tokens (primary text, navy accent, secondary chrome, muted brick
            disclaimer, soft fill, priority green)
Layout:     US Letter portrait, uniform header / footer on every page

All output is a single PDF: figures/preprint/study_explainer.pdf
"""

import logging
import sys
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Rectangle, FancyBboxPatch

# Silence "font not found" warnings — the long font fallback lists below
# are intentional; matplotlib picks whichever is locally installed. The PDF
# renders the same way once a fallback is chosen; the warning is just noise.
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)


# ── Design tokens ─────────────────────────────────────────────────────────────

# Colour palette (hex)
PRIMARY    = "#1A1A1A"
NAVY       = "#1C3A5E"
SECONDARY  = "#6B6B73"
BRICK      = "#8B2E2E"
SOFT_FILL  = "#F5F5F7"
PRIORITY   = "#1E6E3F"

# Page geometry (inches)
PAGE_W, PAGE_H = 8.5, 11.0
MARGIN_L, MARGIN_R = 1.00, 1.00
MARGIN_T, MARGIN_B = 0.90, 0.90
HEADER_Y_IN = PAGE_H - 0.55   # header rule sits at y=10.45
FOOTER_Y_IN = 0.60

CONTENT_LEFT_IN  = MARGIN_L
CONTENT_RIGHT_IN = PAGE_W - MARGIN_R
CONTENT_WIDTH_IN = CONTENT_RIGHT_IN - CONTENT_LEFT_IN

STUDY_TAG = "SOX9 AND CHOLANGIOCYTE IDENTITY LOSS IN PSC"

# Font families. Order matters: matplotlib walks the list and stops at the
# first one it can resolve. macOS-bundled fonts come first so the warning
# log stays clean on a default install. The fancier Source / Crimson / Inter
# entries are kept as graceful upgrades for systems that have them.
SANS_FAMILY  = ["Helvetica Neue", "Helvetica", "Arial",
                "Inter", "DejaVu Sans", "sans-serif"]
SERIF_FAMILY = ["Georgia", "Palatino", "Times New Roman",
                "Source Serif Pro", "Crimson Pro", "DejaVu Serif", "serif"]


FIG_DIR = Path("figures") / "preprint"
OUT_PDF = FIG_DIR / "study_explainer.pdf"


# ── Helpers ───────────────────────────────────────────────────────────────────

def in_to_x(x_in: float) -> float:
    return x_in / PAGE_W


def in_to_y(y_in: float) -> float:
    return y_in / PAGE_H


def w_to_x(w_in: float) -> float:
    return w_in / PAGE_W


def h_to_y(h_in: float) -> float:
    return h_in / PAGE_H


def apply_rcparams() -> None:
    mpl.rcParams.update({
        "font.family":         SERIF_FAMILY,
        "font.size":           10.5,
        "pdf.fonttype":        42,
        "ps.fonttype":         42,
        "axes.spines.top":     False,
        "axes.spines.right":   False,
        "axes.spines.left":    False,
        "axes.spines.bottom":  False,
        "savefig.dpi":         300,
    })


def new_page() -> plt.Figure:
    fig = plt.figure(figsize=(PAGE_W, PAGE_H))
    fig.patch.set_facecolor("white")
    return fig


def draw_header(fig: plt.Figure, page_n: int, total: int,
                section: str) -> None:
    """Hairline rule + small-caps chrome line on either side."""
    # Hairline rule
    rule_y = in_to_y(HEADER_Y_IN)
    fig.add_artist(Rectangle(
        (in_to_x(CONTENT_LEFT_IN), rule_y),
        w_to_x(CONTENT_WIDTH_IN), h_to_y(0.008),
        facecolor=NAVY, edgecolor="none"))
    # Left chrome (study tag + section)
    chrome_y = in_to_y(HEADER_Y_IN + 0.20)
    chrome = f"{STUDY_TAG}  ·  {section.upper()}"
    fig.text(in_to_x(CONTENT_LEFT_IN), chrome_y, chrome,
             fontsize=7.5, color=SECONDARY, family=SANS_FAMILY,
             fontweight="regular", ha="left", va="center")
    # Right chrome (page indicator)
    fig.text(in_to_x(CONTENT_RIGHT_IN), chrome_y, f"PAGE {page_n} / {total}",
             fontsize=7.5, color=SECONDARY, family=SANS_FAMILY,
             ha="right", va="center")


def draw_footer(fig: plt.Figure) -> None:
    """Hairline rule + centred italic disclaimer line."""
    rule_y = in_to_y(FOOTER_Y_IN)
    fig.add_artist(Rectangle(
        (in_to_x(CONTENT_LEFT_IN), rule_y),
        w_to_x(CONTENT_WIDTH_IN), h_to_y(0.005),
        facecolor=SECONDARY, edgecolor="none"))
    fig.text(0.5, in_to_y(FOOTER_Y_IN - 0.22),
             "Exploratory single-cell analysis   ·   Not a causal claim   ·   "
             "Not a therapeutic claim",
             fontsize=8, color=SECONDARY, family=SANS_FAMILY,
             fontstyle="italic", ha="center", va="center")


def draw_page_title(fig: plt.Figure, title: str, y_in: float = 9.85,
                    rule_width_in: float = 1.4,
                    max_chars: int = 52) -> None:
    """Page title in navy semibold + thin accent rule below.

    If the title exceeds max_chars, reduce font size so it fits the content
    width without breaking onto two lines.
    """
    base_fontsize = 18
    fontsize = base_fontsize
    if len(title) > max_chars:
        fontsize = max(13, int(base_fontsize * max_chars / len(title)))
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(y_in), title,
             fontsize=fontsize, color=NAVY, family=SANS_FAMILY,
             fontweight="semibold", ha="left", va="center")
    fig.add_artist(Rectangle(
        (in_to_x(CONTENT_LEFT_IN), in_to_y(y_in - 0.22)),
        w_to_x(rule_width_in), h_to_y(0.018),
        facecolor=NAVY, edgecolor="none"))


def draw_callout(fig: plt.Figure, top_in: float, label: str,
                 body: str, wrap_width: int = 92) -> float:
    """Light left-bordered callout (no heavy fill), height-fitted to text.

    Returns the bottom y_in of the callout so the caller can verify it
    stayed above the footer.
    """
    body_wrapped = textwrap.fill(body, width=wrap_width)
    n_lines = body_wrapped.count("\n") + 1
    # Line height at 10.5pt italic serif with 1.42 linespacing
    line_h_in = 0.20
    label_block_in = 0.20         # label area
    body_block_in  = n_lines * line_h_in
    tail_pad_in    = 0.20
    total_h_in     = label_block_in + body_block_in + tail_pad_in
    bottom_in = top_in - total_h_in

    bar_left_in = CONTENT_LEFT_IN
    bar_width_in = 0.04
    fig.add_artist(Rectangle(
        (in_to_x(bar_left_in), in_to_y(bottom_in)),
        w_to_x(bar_width_in), h_to_y(total_h_in),
        facecolor=NAVY, edgecolor="none"))
    # Label
    fig.text(in_to_x(bar_left_in + 0.22), in_to_y(top_in - 0.10),
             label, fontsize=8.5, color=NAVY, family=SANS_FAMILY,
             fontweight="bold", ha="left", va="center")
    # Body (italic serif)
    fig.text(in_to_x(bar_left_in + 0.22), in_to_y(top_in - 0.32),
             body_wrapped, fontsize=10.5, color=PRIMARY,
             family=SERIF_FAMILY, fontstyle="italic",
             ha="left", va="top", linespacing=1.42)
    return bottom_in


def draw_disclaimer_pill(fig: plt.Figure, center_y_in: float) -> None:
    """Rounded-rectangle pill centred horizontally."""
    pill_w_in = 5.4
    pill_h_in = 0.36
    pill_left_in = (PAGE_W - pill_w_in) / 2
    pill_bot_in  = center_y_in - pill_h_in / 2
    box = FancyBboxPatch(
        (in_to_x(pill_left_in), in_to_y(pill_bot_in)),
        w_to_x(pill_w_in), h_to_y(pill_h_in),
        boxstyle="round,pad=0,rounding_size=0.014",
        linewidth=0.7, edgecolor=BRICK, facecolor="white",
        transform=fig.transFigure,
    )
    fig.add_artist(box)
    fig.text(0.5, in_to_y(center_y_in),
             "EXPLORATORY ANALYSIS   ·   NOT A CAUSAL CLAIM   ·   "
             "NOT A THERAPEUTIC CLAIM",
             fontsize=9, color=BRICK, family=SANS_FAMILY,
             fontweight="bold", ha="center", va="center")


def draw_caption_bullets(fig: plt.Figure, top_in: float,
                         bullets: list[tuple[str, str]]) -> float:
    """Bullet caption: each entry = (panel_letter, body).

    Renders the panel letter as a bold sans-serif label followed by a
    serif body, with a hanging indent so wrapped lines align under the
    body, not under the bullet. Returns the y_in below the last bullet.
    """
    cur_y = top_in
    label_indent_in = 0.18
    body_indent_in  = 0.55
    line_h_in = 0.18
    para_gap_in = 0.20

    for label, body in bullets:
        # Bold panel letter
        fig.text(in_to_x(CONTENT_LEFT_IN + label_indent_in), in_to_y(cur_y),
                 label, fontsize=11.5, color=NAVY, family=SANS_FAMILY,
                 fontweight="bold", ha="left", va="top")
        # Body wrapped, hanging indent
        body_lines = textwrap.wrap(body, width=84)
        for i, line in enumerate(body_lines):
            fig.text(in_to_x(CONTENT_LEFT_IN + body_indent_in),
                     in_to_y(cur_y - i * line_h_in),
                     line, fontsize=10.5, color=PRIMARY,
                     family=SERIF_FAMILY, ha="left", va="top",
                     linespacing=1.4)
        cur_y -= line_h_in * max(len(body_lines), 1) + para_gap_in
    return cur_y


# ── Page builders ─────────────────────────────────────────────────────────────

def page_cover(pdf, total: int) -> None:
    fig = new_page()
    draw_header(fig, 1, total, "Cover")

    # Display title — three lines, left-aligned
    title_lines = [
        "SOX9 collapse in cholangiocytes:",
        "the molecular phenotype of identity loss",
        "in primary sclerosing cholangitis",
    ]
    title_top_in = 9.85
    line_step = 0.45
    for i, line in enumerate(title_lines):
        fig.text(in_to_x(CONTENT_LEFT_IN),
                 in_to_y(title_top_in - i * line_step), line,
                 fontsize=22, color=NAVY, family=SANS_FAMILY,
                 fontweight="semibold", ha="left", va="center")
    # Accent rule under title
    fig.add_artist(Rectangle(
        (in_to_x(CONTENT_LEFT_IN),
         in_to_y(title_top_in - len(title_lines) * line_step + 0.08)),
        w_to_x(1.4), h_to_y(0.02),
        facecolor=NAVY, edgecolor="none"))

    # Subtitle (italic serif, single line)
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(7.95),
             "A single-cell transcriptomic exploration of cholangiocyte "
             "identity loss,", fontsize=11.5, color=SECONDARY,
             family=SERIF_FAMILY, fontstyle="italic", ha="left", va="center")
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(7.72),
             "with idiopathic adulthood ductopenia as the motivating disease "
             "context.",
             fontsize=11.5, color=SECONDARY, family=SERIF_FAMILY,
             fontstyle="italic", ha="left", va="center")

    # Disclaimer pill
    draw_disclaimer_pill(fig, center_y_in=7.05)

    # Abstract heading
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(6.30),
             "ABSTRACT", fontsize=9.5, color=NAVY, family=SANS_FAMILY,
             fontweight="bold", ha="left", va="center")

    # Abstract body — narrower column for elegance
    abstract = (
        "Idiopathic adulthood ductopenia is a rare cholestatic syndrome "
        "characterised by progressive loss of intrahepatic bile ducts; "
        "no public single-cell IAD datasets exist. Using primary "
        "sclerosing cholangitis (PSC) as a tractable cholangiopathy "
        "proxy, we analysed 21,038 cholangiocytes pooled from 16 public "
        "datasets — 4,092 from PSC donors. The most consistent molecular "
        "phenotype of within-disease progression is collapse of SOX9, the "
        "core cholangiocyte transcription factor. SOX9-low PSC "
        "cholangiocytes show coordinated loss of the canonical biliary-"
        "identity program. Candidate upstream pathways remain limited: "
        "oxidative stress retains borderline directionality at single-cell "
        "resolution; IL6 / STAT3 — the strongest donor-level signal — does "
        "not survive cell-level testing. We nominate SOX9-low "
        "cholangiocytes as a candidate molecular surrogate for the "
        "cholangiocyte-loss endpoint relevant to ductopenia. No causal "
        "upstream regulator is identified."
    )
    abstract_lines = textwrap.wrap(abstract, width=92)
    line_h = 0.20
    for i, line in enumerate(abstract_lines):
        fig.text(in_to_x(CONTENT_LEFT_IN),
                 in_to_y(6.05 - i * line_h),
                 line, fontsize=10.5, color=PRIMARY,
                 family=SERIF_FAMILY, ha="left", va="top",
                 linespacing=1.45)

    # Table of contents
    toc_top_in = 2.85
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(toc_top_in),
             "CONTENTS", fontsize=9.5, color=NAVY, family=SANS_FAMILY,
             fontweight="bold", ha="left", va="center")
    entries = [
        ("Figure 1", "SOX9 collapse along PSC progression"),
        ("Figure 2", "SOX9 captures a coordinated biliary-identity program"),
        ("Figure 3", "Donor-level vs cell-level evidence for candidate "
                     "upstream pathways"),
        ("",         "Limitations and next experimental steps"),
        ("",         "References and methods"),
    ]
    entry_y = toc_top_in - 0.32
    for label, body in entries:
        fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(entry_y),
                 label, fontsize=10, color=NAVY, family=SANS_FAMILY,
                 fontweight="bold", ha="left", va="center")
        fig.text(in_to_x(CONTENT_LEFT_IN + 1.10), in_to_y(entry_y),
                 body, fontsize=10, color=PRIMARY, family=SERIF_FAMILY,
                 ha="left", va="center")
        # Thin row separator
        fig.add_artist(Rectangle(
            (in_to_x(CONTENT_LEFT_IN), in_to_y(entry_y - 0.13)),
            w_to_x(CONTENT_WIDTH_IN), h_to_y(0.003),
            facecolor="#dcdcdc", edgecolor="none"))
        entry_y -= 0.30

    draw_footer(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_figure(pdf, page_n: int, total: int, section: str,
                title: str, png_path: Path,
                caption_bullets: list[tuple[str, str]],
                takeaway: str) -> None:
    fig = new_page()
    draw_header(fig, page_n, total, section)
    draw_page_title(fig, title)

    # Embedded figure (crop top 5–6% to remove standalone suptitle)
    img = mpimg.imread(str(png_path))
    crop_top = int(img.shape[0] * 0.055)
    img = img[crop_top:, :]
    h, w = img.shape[0], img.shape[1]
    aspect = h / w

    img_top_in = 9.30
    img_bot_in = 5.20    # 4.10" tall ceiling
    avail_h_in = img_top_in - img_bot_in
    img_w_in = min(CONTENT_WIDTH_IN, avail_h_in / aspect)
    img_h_in = img_w_in * aspect
    if img_h_in > avail_h_in:
        img_h_in = avail_h_in
        img_w_in = img_h_in / aspect
    img_left_in = MARGIN_L + (CONTENT_WIDTH_IN - img_w_in) / 2
    img_bottom_in = img_top_in - img_h_in
    ax_img = fig.add_axes([
        in_to_x(img_left_in), in_to_y(img_bottom_in),
        w_to_x(img_w_in), h_to_y(img_h_in),
    ])
    ax_img.imshow(img); ax_img.axis("off")

    # Caption bullets
    caption_top_in = img_bottom_in - 0.40
    last_y = draw_caption_bullets(fig, caption_top_in, caption_bullets)

    # Key takeaway callout — height-fitted, sits ~0.2" below caption
    callout_top_in = max(last_y - 0.10, 2.30)
    draw_callout(fig, callout_top_in, "KEY TAKEAWAY", takeaway)

    draw_footer(fig)
    pdf.savefig(fig)
    plt.close(fig)


REFERENCES = [
    ("Andrews TS, Nakib D, Perciani CT, et al.",
     "Single-cell, single-nucleus, and spatial transcriptomics "
     "characterization of the immunological landscape in the healthy and PSC "
     "human liver.",
     "J Hepatol 2024. doi:10.1016/j.jhep.2023.12.023"),
    ("CZ CELLxGENE Discover.",
     "Curated single-cell data portal.",
     "https://cellxgene.cziscience.com"),
    ("Wolf FA, Angerer P, Theis FJ.",
     "SCANPY: large-scale single-cell gene expression data analysis.",
     "Genome Biol 2018; 19:15. doi:10.1186/s13059-017-1382-0"),
    ("Virshup I, Rybakov S, Theis FJ, Angerer P, Wolf FA.",
     "anndata: access and store annotated data matrices.",
     "J Open Source Softw 2024; 9:4371. doi:10.21105/joss.04371"),
    ("Haghverdi L, Buettner M, Wolf FA, Buettner F, Theis FJ.",
     "Diffusion pseudotime robustly reconstructs lineage branching.",
     "Nat Methods 2016; 13:845–848. doi:10.1038/nmeth.3971"),
    ("Wolf FA, Hamey FK, Plass M, et al.",
     "PAGA: graph abstraction reconciles clustering with trajectory inference "
     "through a topology preserving map of single cells.",
     "Genome Biol 2019; 20:59. doi:10.1186/s13059-019-1663-x"),
    ("Benjamini Y, Hochberg Y.",
     "Controlling the false discovery rate: a practical and powerful approach "
     "to multiple testing.",
     "J R Stat Soc B 1995; 57:289–300."),
    ("Antoniou A, Raynaud P, Cordi S, et al.",
     "Intrahepatic bile ducts develop according to a new mode of "
     "tubulogenesis regulated by the transcription factor SOX9.",
     "Gastroenterology 2009; 136:2325–2333. doi:10.1053/j.gastro.2009.02.051"),
    ("Ludwig J, Wiesner RH, LaRusso NF.",
     "Idiopathic adulthood ductopenia: a cause of chronic cholestatic liver "
     "disease and biliary cirrhosis.",
     "J Hepatol 1988; 7:193–199. doi:10.1016/S0168-8278(88)80482-3"),
]


def page_references(pdf, page_n: int, total: int) -> None:
    fig = new_page()
    draw_header(fig, page_n, total, "References & methods")
    draw_page_title(fig, "References and methods")

    # References heading
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(9.20),
             "REFERENCES", fontsize=11, color=NAVY, family=SANS_FAMILY,
             fontweight="bold", ha="left", va="center")

    # Numbered references
    cur_y = 8.92
    line_h = 0.16
    item_gap = 0.16
    body_indent_in = CONTENT_LEFT_IN + 0.35

    for idx, (authors, title, venue) in enumerate(REFERENCES):
        # Number marker
        marker = f"{idx + 1}."
        fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(cur_y),
                 marker, fontsize=9.5, color=NAVY, family=SANS_FAMILY,
                 fontweight="bold", ha="left", va="top")
        # Authors (serif, regular)
        wrap_w = int((CONTENT_RIGHT_IN - body_indent_in) / 0.062)
        auth_lines = textwrap.wrap(authors, width=wrap_w)
        for i, line in enumerate(auth_lines):
            fig.text(in_to_x(body_indent_in),
                     in_to_y(cur_y - i * line_h),
                     line, fontsize=9, color=PRIMARY,
                     family=SERIF_FAMILY, ha="left", va="top",
                     linespacing=1.4)
        # Title (italic serif)
        title_y = cur_y - len(auth_lines) * line_h
        title_lines = textwrap.wrap(title, width=wrap_w)
        for i, line in enumerate(title_lines):
            fig.text(in_to_x(body_indent_in),
                     in_to_y(title_y - i * line_h),
                     line, fontsize=9, color=PRIMARY,
                     family=SERIF_FAMILY, fontstyle="italic",
                     ha="left", va="top", linespacing=1.4)
        # Venue (regular)
        venue_y = title_y - len(title_lines) * line_h
        venue_lines = textwrap.wrap(venue, width=wrap_w)
        for i, line in enumerate(venue_lines):
            fig.text(in_to_x(body_indent_in),
                     in_to_y(venue_y - i * line_h),
                     line, fontsize=8.5, color=SECONDARY,
                     family=SERIF_FAMILY, ha="left", va="top",
                     linespacing=1.4)
        total_lines = (len(auth_lines) + len(title_lines) +
                        len(venue_lines))
        cur_y -= total_lines * line_h + item_gap

    # Software & data heading
    sd_y = cur_y - 0.10
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(sd_y),
             "SOFTWARE & DATA AVAILABILITY",
             fontsize=11, color=NAVY, family=SANS_FAMILY,
             fontweight="bold", ha="left", va="center")

    sd_body = (
        "All analyses were performed in Python 3.12 using scanpy (≥ 1.11), "
        "anndata (≥ 0.12), pandas, scipy.stats, and matplotlib. Single-cell "
        "data were obtained via the cellxgene-census Python API; the "
        "cholangiocyte subset comprises 21,038 cells pooled from 16 source "
        "datasets registered in CELLxGENE Discover. Source code, "
        "intermediate result tables, and the full analytical audit trail are "
        "available in the project repository under IAD/src/, IAD/results/, "
        "and IAD/figures/preprint/."
    )
    wrap_w = int(CONTENT_WIDTH_IN / 0.062)
    sd_lines = textwrap.wrap(sd_body, width=wrap_w)
    body_top_y = sd_y - 0.30
    for i, line in enumerate(sd_lines):
        fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(body_top_y - i * line_h),
                 line, fontsize=9.5, color=PRIMARY,
                 family=SERIF_FAMILY, ha="left", va="top",
                 linespacing=1.45)

    draw_footer(fig)
    pdf.savefig(fig)
    plt.close(fig)


def page_limits_and_next(pdf, page_n: int, total: int) -> None:
    fig = new_page()
    draw_header(fig, page_n, total, "Limitations & next steps")
    draw_page_title(fig, "Limitations and next experimental steps")

    # Limitations heading (muted brick)
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(9.20),
             "LIMITATIONS", fontsize=11, color=BRICK, family=SANS_FAMILY,
             fontweight="bold", ha="left", va="center")

    limitations = [
        ("Sample size.",
         "The PSC cohort comprises ten donors across public single-cell "
         "data, with seven from a single source dataset. Within-PSC "
         "inference is statistically conservative."),
        ("Technology heterogeneity.",
         "The dominant PSC dataset is single-nucleus RNA-seq. Directional "
         "findings are preserved across technology; magnitudes are "
         "attenuated for cytoplasmic transcripts."),
        ("Multiple-testing correction.",
         "At this donor count, formal genome-wide FDR does not establish "
         "single-gene significance. Inferential weight rests on direction "
         "concordance across independent datasets and within-donor "
         "cell-level dose-response."),
        ("Receptor coverage.",
         "HVG-based gene selection removed most pathway receptor genes "
         "from the source AnnData, limiting receptor-side analysis."),
        ("Tissue-level signalling.",
         "Cell-cell communication analysis requires non-cholangiocyte "
         "liver cell types not present in the current data."),
        ("Disease scope.",
         "No IAD samples are present. All IAD-relevant inference is by "
         "analogy to PSC."),
        ("Direction of effect.",
         "Pathway correlations with SOX9 cannot distinguish drivers from "
         "consequences without experimental perturbation."),
    ]
    cur_y = 8.90
    label_indent = CONTENT_LEFT_IN + 0.20
    body_indent  = CONTENT_LEFT_IN + 0.20
    line_h = 0.16
    item_gap = 0.14
    for idx, (label, body) in enumerate(limitations):
        marker = f"{idx + 1}."
        fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(cur_y),
                 marker, fontsize=10, color=NAVY, family=SANS_FAMILY,
                 fontweight="bold", ha="left", va="top")
        # Bold label inline
        fig.text(in_to_x(label_indent), in_to_y(cur_y),
                 label, fontsize=10, color=PRIMARY, family=SERIF_FAMILY,
                 fontweight="bold", ha="left", va="top")
        # Body wraps with hanging indent under the body_indent column
        label_run_w_in = len(label) * 0.08    # rough char width
        first_line_indent_in = label_indent + label_run_w_in + 0.06
        wrap_w = max(1, int((CONTENT_RIGHT_IN - first_line_indent_in) / 0.072))
        body_lines = textwrap.wrap(body, width=wrap_w)
        if body_lines:
            # First line goes after the label
            fig.text(in_to_x(first_line_indent_in), in_to_y(cur_y),
                     body_lines[0], fontsize=10, color=PRIMARY,
                     family=SERIF_FAMILY, ha="left", va="top",
                     linespacing=1.4)
            # Subsequent lines wrap under body_indent
            wrap_w2 = max(1, int((CONTENT_RIGHT_IN - body_indent) / 0.072))
            rest = textwrap.wrap(" ".join(body_lines[1:]), width=wrap_w2)
            for j, line in enumerate(rest):
                fig.text(in_to_x(body_indent),
                         in_to_y(cur_y - (j + 1) * line_h),
                         line, fontsize=10, color=PRIMARY,
                         family=SERIF_FAMILY, ha="left", va="top",
                         linespacing=1.4)
            n_lines = 1 + len(rest)
        else:
            n_lines = 1
        cur_y -= n_lines * line_h + item_gap

    # Next steps heading (priority green)
    nx_heading_y = cur_y - 0.10
    fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(nx_heading_y),
             "NEXT EXPERIMENTAL STEPS", fontsize=11, color=PRIORITY,
             family=SANS_FAMILY, fontweight="bold", ha="left", va="center")

    next_steps = [
        ("Priority 1 — In vitro perturbation.",
         "Treat human cholangiocyte organoids with sub-lethal oxidative-"
         "stress conditions (paraquat or H2O2) and quantify SOX9, KRT19, "
         "HNF1B, and AQP1 by qPCR and immunofluorescence at 24, 48, and "
         "72 hours. A monotonic SOX9 decline accompanied by coordinated "
         "biliary-marker loss would support the within-cell rho = -0.12 "
         "signal as a real upstream relationship."),
        ("Priority 2 — Spatial validation.",
         "Andrews et al. (J Hepatol 2024) deposited four Visium sections "
         "from one PSC patient on CELLxGENE Discover. Co-localisation of "
         "SOX9-low transcriptional states with vanishing-duct histology "
         "would extend the finding from dissociated single cells to "
         "tissue context."),
        ("Priority 3 — Computational unblock for cell-cell signalling.",
         "Re-fetch the cholangiocyte AnnData with receptor genes "
         "force-kept, and co-fetch surrounding liver cell types from "
         "Census. This unlocks ligand-receptor analysis to identify "
         "candidate cellular sources of the SOX9-low associated signals."),
    ]
    cur_y = nx_heading_y - 0.32
    for label, body in next_steps:
        fig.text(in_to_x(CONTENT_LEFT_IN), in_to_y(cur_y),
                 label, fontsize=10, color=NAVY, family=SANS_FAMILY,
                 fontweight="bold", ha="left", va="top")
        wrap_w = int(CONTENT_WIDTH_IN / 0.072)
        body_lines = textwrap.wrap(body, width=wrap_w)
        for i, line in enumerate(body_lines):
            fig.text(in_to_x(CONTENT_LEFT_IN + 0.20),
                     in_to_y(cur_y - 0.18 - i * line_h),
                     line, fontsize=10, color=PRIMARY,
                     family=SERIF_FAMILY, ha="left", va="top",
                     linespacing=1.4)
        cur_y -= 0.20 + (len(body_lines) * line_h) + 0.12

    draw_footer(fig)
    pdf.savefig(fig)
    plt.close(fig)


# ── Figure captions & takeaways ────────────────────────────────────────────────

CAP_FIG1 = [
    ("A", "Donor-mean PSC pseudotime vs IAD score. ρ = +0.43, "
          "p = 2 × 10⁻⁴, n = 70 donors — the diffusion-pseudotime "
          "trajectory is consistent with PSC progression."),
    ("B", "Same axis vs donor-mean SOX9 regulon activity. "
          "ρ = −0.72, p = 1.5 × 10⁻¹³."),
    ("C", "Per-dataset replication across 9 source datasets. Bars "
          "marked * are individually significant at p < 0.05; every "
          "significant bar points in the predicted negative direction."),
    ("D", "Technology stratification. Per-cell ρ between SOX9 and "
          "IAD score is negative in both snRNA-seq and scRNA-seq. "
          "snRNA magnitude is attenuated because cytoplasmic mRNAs "
          "such as SOX9 are systematically under-detected in nuclei."),
]
TAKE_FIG1 = (
    "SOX9 collapse along PSC progression replicates across donors, "
    "datasets, and sequencing technologies. This is the most consistent "
    "finding across our tests."
)

CAP_FIG2 = [
    ("A", "Within-donor Spearman ρ between per-cell SOX9 and 19 "
          "curated pathway scores across 8 PSC donors. Only "
          "biliary_identity reaches BH-corrected q < 0.05."),
    ("B", "Donor-mean PROX1 vs donor-mean SOX9 by disease. The two "
          "genes vary largely independently; PROX1 appears to be a "
          "sub-state marker rather than a progression marker."),
    ("C", "PROX1-high vs PROX1-low marker positivity (KRT19+, KRT7+, "
          "EPCAM+, SOX9-high) in the all-cholangiocytes pool and "
          "within PSC. PROX1-high cells consistently carry lower "
          "biliary-marker positivity across both slices."),
]
TAKE_FIG2 = (
    "SOX9 expression tracks the canonical biliary-identity program. "
    "PROX1 marks a separate cholangiocyte sub-state with lower biliary "
    "identity but no disease-progression signal of its own."
)

CAP_FIG3 = [
    ("A", "Donor-level pathway score Δ (SOX9-high − SOX9-low) "
          "across 7 informative PSC donors. Red = lost in SOX9-low; "
          "blue = gained. IL6 / STAT3 and oxidative stress dominate "
          "the gained side."),
    ("B", "Within-donor cell-level ρ between SOX9 and each pathway. "
          "Only biliary_identity passes BH q < 0.05; oxidative stress "
          "is the sole candidate upstream pathway with directionally "
          "consistent within-cell coupling, and remains borderline."),
    ("C", "Donor-level Δ vs cell-level mean ρ. Off-diagonal points "
          "carry donor-level signal that does not replicate at "
          "single-cell resolution — IL6 / STAT3 sits in this region."),
    ("D", "Pathway-category counts under donor-level versus cell-level "
          "criteria. Most donor-level candidates collapse into "
          "ambiguous under cell-level testing."),
]
TAKE_FIG3 = (
    "Only oxidative stress retains direction at cell level, and "
    "borderline. IL6 / STAT3, the strongest donor-level signal, does not "
    "survive the cell-level test. The SOX9-low state appears consistent "
    "with loss of cellular engagement rather than an active stress "
    "response."
)


def main() -> None:
    apply_rcparams()
    if not FIG_DIR.exists():
        print(f"ERROR: {FIG_DIR} missing.", file=sys.stderr); sys.exit(1)
    f1 = FIG_DIR / "figure1_sox9_collapse.png"
    f2 = FIG_DIR / "figure2_identity_loss.png"
    f3 = FIG_DIR / "figure3_pathways.png"
    for p in (f1, f2, f3):
        if not p.exists():
            print(f"ERROR: {p} missing.", file=sys.stderr); sys.exit(1)

    total = 6
    print(f"Building {OUT_PDF} ({total} pages, US Letter portrait) …")
    with PdfPages(OUT_PDF) as pdf:
        page_cover(pdf, total)
        page_figure(pdf, 2, total, "Figure 1",
                    "Figure 1 · SOX9 collapse along PSC progression",
                    f1, CAP_FIG1, TAKE_FIG1)
        page_figure(pdf, 3, total, "Figure 2",
                    "Figure 2 · SOX9 captures coordinated biliary identity",
                    f2, CAP_FIG2, TAKE_FIG2)
        page_figure(pdf, 4, total, "Figure 3",
                    "Figure 3 · Donor-level vs cell-level evidence",
                    f3, CAP_FIG3, TAKE_FIG3)
        page_limits_and_next(pdf, 5, total)
        page_references(pdf, 6, total)
    print(f"  → {OUT_PDF}")


if __name__ == "__main__":
    main()
