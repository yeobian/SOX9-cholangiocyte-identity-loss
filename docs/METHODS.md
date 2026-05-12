# Methods

A full methodological description of the analyses performed in this study. For the conservative scientific interpretation of the results, see `LIMITATIONS.md`. For the executive summary, see `PROJECT_SUMMARY.md`. For step-by-step run instructions, see the main `README.md`.

## 1. Data acquisition

Single-cell RNA-seq data was fetched from CELLxGENE Census (CZ CELLxGENE Discover, stable release accessed via the `cellxgene-census` Python API). The query selected cells matching:

* `tissue_general == 'liver'`
* `cell_type` containing one of: `cholangiocyte`, `intrahepatic cholangiocyte`, `epithelial cell of bile duct`, `duct epithelial cell`.

The resulting cholangiocyte AnnData contains **21,038 cells** drawn from **16 source datasets** (105 unique donors), with the following disease distribution: 13,169 normal cells (75 donors), 4,092 PSC cells (10 donors), 3,212 PBC cells (3 donors), and 565 cells from other conditions (intestinal failure–associated liver disease, colorectal carcinoma metastases). The dominant PSC contribution (3,741 of 4,092 PSC cells) comes from one dataset (CELLxGENE `1873a18a-66fd-…`), derived from Andrews TS et al., *J Hepatol* 2024.

## 2. Preprocessing

Standard scanpy pipeline:

1. `sc.pp.calculate_qc_metrics` for mitochondrial percentage and gene-count distribution.
2. Cells filtered by `min_genes=200`, `max_genes=8000`, `pct_counts_mt<=25`.
3. Genes filtered by `min_cells=3`.
4. `sc.pp.normalize_total(target_sum=1e4)` followed by `sc.pp.log1p`.
5. Highly variable gene (HVG) selection (`min_mean=0.0125`, `max_mean=3`, `min_disp=0.5`), with `subset=False` so that a curated panel of cholangiocyte / hepatocyte / plasticity / inflammation markers — including SOX9, KRT19, KRT7, EPCAM, HNF1B, AQP1, CFTR, PROX1, HNF4A, FOXA2 — was force-kept regardless of HVG selection.
6. `sc.pp.scale(max_value=10)`.

The resulting working AnnData has 21,038 cells × 3,610 genes.

## 3. Trajectory inference

Diffusion pseudotime was computed in three nested settings:

1. **Pooled trajectory.** All 16 datasets, all disease states.
2. **PSC + normal trajectory.** Restricted to the cleanest disease vs healthy contrast.
3. **Intra-PSC trajectory.** Restricted to PSC cells only, to ask "what changes *within* disease."

In each case we ran `sc.pp.neighbors` (n_neighbors = 15), `sc.tl.umap`, `sc.tl.leiden`, `sc.tl.diffmap`, and `sc.tl.dpt`, rooting pseudotime in the cluster with the lowest IAD score (a panel-derived measure of cholangiocyte program loss). For the intra-PSC trajectory, PCA was re-fit on the PSC subset and a larger neighbour graph (n_neighbors = 30) was used to keep the diffusion-map connected.

## 4. IAD diagnostic-panel score

Per cell, we compute a weighted score over a 20-gene cholangiocyte diagnostic panel (gene weights = specificity values from earlier exploratory analysis; see `results/iad_diagnostic_panel.tsv`). Each cell's score is normalised to the 95th percentile and inverted, so that:

* IAD score = 0 represents cells with the highest panel expression (healthy cholangiocyte phenotype).
* IAD score = 1 represents cells with low panel expression (candidate "exhausted" cholangiocyte phenotype).

The IAD score is used as the disease-progression axis throughout the analysis.

## 5. Statistical tests applied

### 5.1 Donor-level Spearman

Donor-mean values for pseudotime, IAD score, and individual gene / pathway scores were computed, then Spearman ρ measured across donors. Used to control for cell-state composition differences at the donor level.

### 5.2 Per-cell Spearman within disease

Cell-level Spearman ρ between SOX9 expression and each pathway score, computed separately within each donor, then aggregated cross-donor via one-sample t-test of the per-donor ρ values against zero. Provides much higher statistical power than donor-mean tests at the cost of within-donor structure.

### 5.3 Per-donor pseudobulk DE

For each donor with ≥10 cells in both arms (SOX9-high vs SOX9-low, defined by rank-based tertile split), we computed the per-gene log fold-change as mean(SOX9-high) − mean(SOX9-low) in log1p space. Across donors, we ran a one-sample t-test of the per-donor logFCs against zero, with Benjamini–Hochberg FDR correction across all genes tested.

### 5.4 Per-dataset replication

The donor-level test was repeated *inside* each individual source dataset that contributed ≥200 cells and ≥5 donors. Direction concordance across independent datasets is reported alongside the pooled p-value as a robustness measure.

### 5.5 Sensitivity to dropping the dominant dataset

The headline tests were re-run with the dominant source dataset (`1873a18a-66fd-…`, which contributes 7 of 10 PSC donors) excluded. Findings that survive this exclusion are considered multi-cohort robust; findings that collapse are flagged as one-cohort observations.

### 5.6 Technology stratification (snRNA-seq vs scRNA-seq)

Because the dominant PSC dataset is single-nucleus RNA-seq while the rest are single-cell RNA-seq, all headline tests were re-run separately within the snRNA-seq cohort and the scRNA-seq cohort. snRNA-seq systematically over-represents nuclear-localised transcription-factor transcripts and under-represents cytoplasmic mRNAs (KRT19, SOX9 itself). Findings that preserve direction across both technologies — even when magnitudes differ — are considered technology-independent.

## 6. Pathway scoring

Nineteen curated pathway gene sets (cholangiocyte identity, biliary identity, ductular reaction, hepatocyte identity, progenitor TFs, Notch, YAP/TAZ, BMP, WNT, TGF-β, TNF, IL6, IFN, bile-acid stress, ER stress / UPR, oxidative stress, senescence, EMT, fibrosis, general inflammation) were scored per cell with `scanpy.tl.score_genes`. Gene sets are short, literature-curated lists (8–14 genes each) defined in `src/sox9_mechanism_discovery.py`.

Each pathway was tested for:

* **Donor-level Δ**: mean(SOX9-high) − mean(SOX9-low) within each donor, cross-donor t-test.
* **Per-cell within-donor ρ**: Spearman ρ(SOX9, pathway score) across cells, per donor, aggregated by t-test of the per-donor ρ values.
* **Threshold sensitivity**: tertile / quartile / quintile splits of SOX9 expression.

## 7. Software

* Python 3.11+
* `scanpy >= 1.11`
* `anndata >= 0.10`
* `cellxgene-census >= 1.6`
* `pandas`, `numpy`, `scipy.stats`, `scipy.sparse`
* `matplotlib`, `seaborn`
* `leidenalg`, `python-igraph`, `umap-learn`

See `requirements.txt` for exact versions.

## 8. Reproducibility

Every result table in `results/` is regenerable by running the corresponding script in `src/`. The order is documented in the main `README.md`. Random seeds are fixed throughout (`random_state=42`).

## 9. Code organisation

* `src/download_data.py` and `src/refetch_cholangiocytes_with_prox1.py` — data acquisition with HVG panel force-keep.
* `src/preprocess_liver.py` — QC, normalisation, HVG, clustering.
* `src/trajectory_psc_only.py`, `src/trajectory_intra_psc.py` — diffusion pseudotime, two scopes.
* `src/sox9_mechanism_discovery.py` — pathway scoring + donor-level Δ.
* `src/sox9_per_cell_dose_response.py` — within-donor ρ(SOX9, pathway) cell-level test.
* `src/per_dataset_replication.py` — per-dataset Spearman replication.
* `src/sensitivity_drop_1873a18a.py` — dominant-dataset removal sensitivity.
* `src/technology_stratified_prox1.py` — snRNA-seq vs scRNA-seq stratification.
* `src/make_figure{1,2,3}_*.py` and `src/make_study_explainer_pdf.py` — publication figures and the polished 6-page brief.
