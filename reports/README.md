# `reports/`

Polished, public-facing outputs.

| File | What |
|---|---|
| `study_explainer.pdf` | **6-page scientific brief.** Cover with abstract + table of contents; Figures 1, 2, 3 each on their own page with caption + key takeaway; limitations + next-steps page; references + methods page. Designed to be sent to an advisor, collaborator, or external reviewer as a self-contained explanation of the study. |
| `WRITEUP_DRAFT.md` | Preprint-style prose draft (~1,500 words). Background, methods (compressed), five results sub-sections, discussion, limitations, planned-figure rationale, and proposed experimental follow-ups. |

## How to regenerate

```bash
# the standalone figures must exist first
python src/make_figure1_sox9_collapse.py
python src/make_figure2_identity_loss.py
python src/make_figure3_pathways.py

# then the brief
python src/make_study_explainer_pdf.py
```

The brief is built from the three figure PNGs in `figures/preprint/` + the polished text in `src/make_study_explainer_pdf.py`. All page elements (typography, colour palette, layout, spacing) are specified in code so the document is fully reproducible.

## Conservative framing

Both documents follow the project's framing rules: exploratory analysis, no causal claims, no therapeutic claims, IAD scoped strictly to motivating disease context. The cover disclaimer pill and the footer line on every page of `study_explainer.pdf` make this explicit.
