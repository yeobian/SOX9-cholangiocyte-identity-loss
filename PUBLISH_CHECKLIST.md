# Publishing checklist

Steps to turn `IAD/github/` into a public GitHub repository safely.

---

## Files that **should** go into the repository

| Location | What | Why |
|---|---|---|
| root | `README.md`, `LICENSE`, `requirements.txt`, `.gitignore`, this checklist | standard repository chrome |
| `docs/` | `PROJECT_SUMMARY.md`, `METHODS.md`, `LIMITATIONS.md` | core documentation |
| `references/` | `references.md` | the 9-citation reference list |
| `src/` | all 31 analysis scripts (`.py`) from `IAD/src/` | the work itself |
| `results/` | curated summary TSVs + `RESULT_NOTE.md` audit trail (see list below) | replicable outputs |
| `figures/preprint/` | `figure1_…`, `figure2_…`, `figure3_…` (PNG + PDF) | publication-quality figures |
| `reports/` | `study_explainer.pdf`, `WRITEUP_DRAFT.md` | the 6-page brief + prose draft |
| `notebooks/` | a `README.md` placeholder (and any future demo notebooks) | scaffolding for future tutorials |
| `data/` | only the `README.md` explaining the data-fetch protocol | data itself is too large for git |

### Curated `results/` TSVs to commit

Commit only the **summary** tables; not the per-cell tables (large, mostly recomputable).

```
results/RESULT_NOTE.md
results/WRITEUP_DRAFT.md         # already moved to reports/, but keep optional copy
results/sox9_pathway_summary_psc.tsv
results/sox9_per_cell_dose_response_summary.tsv
results/sox9_receptor_summary_psc.tsv
results/sox9_threshold_sensitivity_psc.tsv
results/per_dataset_replication_summary.tsv
results/technology_stratified_prox1_verdict.tsv
results/sensitivity_summary.tsv
results/trajectory_psc_per_donor.tsv
results/trajectory_intra_psc_correlations.tsv
results/prox1_abundance_summary.tsv
results/prox1_sox9_marker_summary.tsv
results/prox1_high_vs_low_psc_de_top.tsv
results/iad_diagnostic_panel.tsv
```

(Per-cell tables — `*_per_cell.tsv` — are excluded by `.gitignore`; they regenerate from the scripts.)

---

## Files that **should NOT** go into the repository

| What | Where | Why |
|---|---|---|
| `*.h5ad` AnnData files | `data/liver/`, `data/pancreas/`, `data/ramachandran/` | hundreds of MB to GB; fetch on demand instead |
| `GSE*.tar` archives | `data/` | large GEO raw downloads; reproducible via scripts |
| Preview / version-stamped figure folders | `figures/preprint/_preview_/`, `_v2_/`, `_v3_/`, etc. | diagnostic-only, not canonical |
| `figures/comparison/`, `figures/sox9/`, `figures/prox1/`, `figures/liver/`, etc. | various | intermediate plots; the three preprint figures are the canonical public artefacts |
| `.venv/`, `venv/`, `__pycache__/` | various | environments and bytecode |
| `.DS_Store`, `Thumbs.db` | various | OS noise |
| `*.bak` files | `data/liver/` | local backups |
| Per-cell result tables (`*_per_cell.tsv`) | `results/` | regenerable; keep summary tables only |

All of these are caught by the `.gitignore`. The list is here so you can verify before committing.

---

## One-shot bash to populate `IAD/github/` from the working directory

Run this from the **`IAD/`** folder (the parent of `github/`):

```bash
# Make the target subfolders
mkdir -p github/src github/figures/preprint github/reports github/notebooks github/data
mkdir -p github/results

# Code
cp src/*.py github/src/

# Polished figures
cp figures/preprint/figure1_sox9_collapse.{png,pdf} github/figures/preprint/
cp figures/preprint/figure2_identity_loss.{png,pdf}  github/figures/preprint/
cp figures/preprint/figure3_pathways.{png,pdf}        github/figures/preprint/
cp figures/preprint/study_explainer.pdf               github/reports/

# Audit trail + key writeups
cp results/RESULT_NOTE.md     github/results/
cp results/WRITEUP_DRAFT.md   github/reports/

# Curated summary tables (commit these; not the per-cell tables)
for tsv in \
  sox9_pathway_summary_psc.tsv \
  sox9_per_cell_dose_response_summary.tsv \
  sox9_receptor_summary_psc.tsv \
  sox9_threshold_sensitivity_psc.tsv \
  per_dataset_replication_summary.tsv \
  technology_stratified_prox1_verdict.tsv \
  sensitivity_summary.tsv \
  trajectory_psc_per_donor.tsv \
  trajectory_intra_psc_correlations.tsv \
  prox1_abundance_summary.tsv \
  prox1_sox9_marker_summary.tsv \
  prox1_high_vs_low_psc_de_top.tsv \
  iad_diagnostic_panel.tsv; do
  [ -f "results/$tsv" ] && cp "results/$tsv" github/results/
done
```

After running it, inspect the `github/` folder. It should contain ~30 scripts in `src/`, 6 publication-ready figure files in `figures/preprint/`, ~13 result TSVs in `results/`, two PDFs / markdowns in `reports/`, plus the documentation files. No `.h5ad` files, no raw GEO archives, no preview folders.

---

## Pre-push checks

```bash
cd IAD/github

# 1. Sanity: no large data files leaked in
find . -name "*.h5ad" -o -name "*.tar" -o -name "GSE*"      # should print nothing

# 2. Sanity: total repo size under ~50 MB
du -sh .                                                      # expect <50 MB

# 3. Initialise the repo
git init -b main
git add .
git status                                                    # eyeball before commit
git commit -m "Initial commit: SOX9 cholangiocyte single-cell study"

# 4. Add the remote and push
git remote add origin git@github.com:<your-username>/<repo-name>.git
git push -u origin main
```

---

## Recommended repository settings (after pushing)

* Add the README's three-line tagline as the GitHub **About** description.
* Add topics: `single-cell-rna-seq`, `scanpy`, `cholangiocyte`, `liver`, `primary-sclerosing-cholangitis`, `data-science`, `portfolio`.
* Enable Issues (for community questions); disable Wiki and Projects unless you intend to use them.
* Add a `CITATION.cff` if you later post to bioRxiv — it lets GitHub generate a citation button automatically.

---

## A note on patient privacy

All data used in this study is published, de-identified, and distributed through the CELLxGENE Census public infrastructure. No raw sequencing reads, no patient identifiers, and no clinical metadata beyond what is already public are touched anywhere in this repository. The `data/README.md` makes the access protocol explicit.
