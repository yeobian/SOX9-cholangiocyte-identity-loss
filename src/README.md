# `src/`

All analysis scripts. Each script is self-contained and can be run from the repository root with `python src/<script>.py`. Random seeds are fixed throughout (`random_state=42`).

## Recommended run order

```bash
# 1. Acquire data
python src/download_data.py
python src/refetch_cholangiocytes_with_prox1.py

# 2. Build trajectory
python src/preprocess_liver.py
python src/trajectory_psc_only.py
python src/trajectory_intra_psc.py

# 3. Mechanism discovery
python src/sox9_mechanism_discovery.py
python src/sox9_per_cell_dose_response.py

# 4. Robustness checks
python src/per_dataset_replication.py
python src/sensitivity_drop_1873a18a.py
python src/technology_stratified_prox1.py

# 5. Companion PROX1 analyses (kept for completeness; secondary to SOX9 finding)
python src/prox1_dataset_spread.py
python src/prox1_per_cell_pseudotime.py
python src/prox1_abundance_test.py
python src/prox1_sox9_association.py
python src/prox1_high_vs_low_psc_de.py

# 6. Build figures + brief
python src/make_figure1_sox9_collapse.py
python src/make_figure2_identity_loss.py
python src/make_figure3_pathways.py
python src/make_study_explainer_pdf.py
```

## Script roles

| Script | Role |
|---|---|
| `download_data.py` | CELLxGENE Census fetch (liver + optional pancreas) |
| `refetch_cholangiocytes_with_prox1.py` | re-fetch with PROX1 + curated panel force-kept through HVG |
| `preprocess_liver.py` | QC, normalisation, HVG, clustering |
| `trajectory_psc_only.py` | PAGA + DPT on PSC + normal cells |
| `trajectory_intra_psc.py` | PAGA + DPT on PSC cells only |
| `sox9_mechanism_discovery.py` | 19-pathway scoring + donor-level Δ for SOX9-high vs low |
| `sox9_per_cell_dose_response.py` | within-donor ρ(SOX9, pathway score) at single-cell level |
| `per_dataset_replication.py` | replication of donor-level signal in each individual source dataset |
| `sensitivity_drop_1873a18a.py` | re-run headline tests with the dominant cohort excluded |
| `technology_stratified_prox1.py` | snRNA-seq vs scRNA-seq stratification |
| `make_figure{1,2,3}_*.py` | publication-quality figures (PNG + PDF) |
| `make_study_explainer_pdf.py` | the polished 6-page scientific brief |
| `figure_style.py` | shared matplotlib style for the three preprint figures |
| `utils.py`, `liver_modules.py`, `iad_score.py` | shared helpers and gene-set definitions |

## Style

* Pure Python 3.11+, scanpy / anndata / pandas / scipy / matplotlib.
* No notebooks; scripts only.
* Pipeline shape (`load → filter → transform → compute → save`) used consistently.
* Small functions with single responsibilities.
* Inputs validated early; informative error messages with run-instructions.
* Cautious scientific framing in docstrings — no causal language.
