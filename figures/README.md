# `figures/`

Publication-quality figures from the study. Each figure is provided as both PNG (300 DPI) and PDF (editable text, Type 42 fonts) — drop directly into a manuscript, preprint, or talk.

## `figures/preprint/`

| File | Description |
|---|---|
| `figure1_sox9_collapse.{png,pdf}` | **Figure 1** — SOX9 collapse along PSC progression. Four panels: pseudotime vs IAD score, pseudotime vs SOX9 regulon activity, per-dataset replication, technology stratification. |
| `figure2_identity_loss.{png,pdf}` | **Figure 2** — SOX9 captures coordinated biliary identity; PROX1 marks a separate sub-state. Three panels: within-donor ρ(SOX9, pathway), PROX1↔SOX9 at donor level, marker positivity comparison. |
| `figure3_pathways.{png,pdf}` | **Figure 3** — Donor-level vs cell-level evidence for candidate upstream pathways. Four panels: per-donor Δ heatmap, within-donor ρ boxplots, donor-vs-cell scatter, category summary. |

All three figures are also embedded in the 6-page brief at `reports/study_explainer.pdf`.

## What's *not* in this folder (excluded by `.gitignore`)

The standalone figure-generation pipeline produces many intermediate plots under `figures/comparison/`, `figures/sox9/`, `figures/prox1/`, `figures/liver/`, etc., as part of the diagnostic analysis. Those are useful for development but cluttered and version-stamped (`_v2_`, `_v3_`, `_preview_` folders). They are excluded from the public repository. The three preprint figures above are the canonical public artefacts.

## Regenerating

```bash
python src/make_figure1_sox9_collapse.py
python src/make_figure2_identity_loss.py
python src/make_figure3_pathways.py
```

All three scripts read summary TSVs from `results/` and produce both PNG + PDF in this folder.
