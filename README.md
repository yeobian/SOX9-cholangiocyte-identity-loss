# SOX9 collapse in cholangiocytes — a single-cell exploration

*A computational single-cell transcriptomic study of cholangiocyte identity loss in primary sclerosing cholangitis (PSC), motivated by idiopathic adulthood ductopenia (IAD).*

> **Exploratory analysis · Not a causal claim · Not a therapeutic claim.**

---

## Project overview

Idiopathic adulthood ductopenia (IAD) is a rare cholestatic syndrome characterised by progressive loss of intrahepatic bile ducts. Because no public single-cell IAD datasets exist, this project uses primary sclerosing cholangitis (PSC) as a tractable cholangiopathy proxy and asks one question:

> What transcriptional changes mark the cholangiocyte-identity-loss endpoint relevant to ductopenia?

Across 21,038 human cholangiocytes pooled from 16 CELLxGENE Census datasets (4,092 from PSC donors), the most consistent within-disease signal is collapse of **SOX9**, the core cholangiocyte transcription factor. SOX9-low cells show coordinated loss of the canonical biliary identity program (KRT19, KRT7, EPCAM, HNF1B, AQP1, CFTR). Candidate upstream pathway analysis is intentionally cautious — only oxidative stress retains direction at single-cell resolution, and remains borderline.

---

## Why this matters

* **A defined molecular surrogate.** SOX9-low cholangiocytes are nominated as a candidate molecular surrogate for the cholangiocyte-loss endpoint relevant to ductopenia. The surrogate is observable in existing public data, requires no new sample collection, and survives multi-layer testing.
* **A worked example of cautious single-cell inference.** The analysis trail explicitly records which findings survived which tests — and which did not. PROX1 was initially the focal candidate; it was demoted to a sub-state marker after replication checks. The donor-level IL6/STAT3 signal did not survive cell-level testing. This kind of self-critique is rare in single-cell preprints and is a deliberate feature of the project.
* **Reusable analytical infrastructure.** Every step is scripted, donor-controlled where possible, and accompanied by a sensitivity test. The pipeline can be repointed at any cholangiocyte cohort with minor edits.

---

## Dataset

| Source | Cells | Donors | Notes |
|---|---:|---:|---|
| CELLxGENE Census, human liver | 21,038 cholangiocytes | 105 donors | 16 source datasets |
| · normal | 13,169 | 75 | reference baseline |
| · PSC | 4,092 | 10 | disease proxy |
| · PBC | 3,212 | 3 | secondary cholangiopathy |
| · other (IFALD, CRC-metastasis) | 565 | 17 | descriptive only |

The dominant PSC contribution comes from Andrews TS et al., *J Hepatol* 2024 (CELLxGENE collection `0c8a364b-…`), which deposited a paired snRNA-seq + scRNA-seq cohort plus four Visium spatial sections.

**Data is not committed to this repository.** It is fetched on demand from CELLxGENE Census via the provided scripts; see `data/README.md` for the protocol.

---

## Methods

* **Data acquisition.** CELLxGENE Census stable release, accessed via the `cellxgene-census` Python API. Cells selected by tissue (`liver`) and ontology cell-type (`cholangiocyte`, `intrahepatic cholangiocyte`).
* **Preprocessing.** Standard scanpy pipeline (`normalize_total` → `log1p` → highly variable genes), with the SOX9 / KRT19 / EPCAM / HNF1B / AQP1 / CFTR / PROX1 panel force-kept through HVG selection.
* **Trajectory inference.** PAGA + diffusion pseudotime (DPT), rooted in the cluster with the lowest IAD score. Three trajectories were built: pooled, PSC + normal, and intra-PSC only.
* **Statistical core.** Donor-level Spearman ρ between donor-mean pseudotime and donor-mean gene/score; per-cell Spearman ρ within disease; per-donor pseudobulk DE (one-sample t-test on per-donor logFCs across cross-donor distribution); per-dataset replication; sensitivity to dropping the dominant dataset; technology stratification (snRNA-seq vs scRNA-seq).
* **Pathway scoring.** 19 curated pathway gene sets scored per cell with `scanpy.tl.score_genes`. Per-donor Δ (SOX9-high − SOX9-low) and per-cell within-donor ρ(SOX9, pathway) both computed.
* **Multiple-testing control.** Benjamini–Hochberg FDR within each analysis. Direction concordance across independent datasets and within-donor cell-level dose-response carry inferential weight at this n.

See `docs/METHODS.md` for the full methods write-up.

---

## Key findings

1. **SOX9 collapse along PSC progression** is the most consistent within-disease signal. It replicates across donors, datasets, and sequencing technologies (snRNA-seq vs scRNA-seq).
2. **SOX9-low PSC cholangiocytes are a coordinated biliary-identity-loss state.** KRT19, KRT7, EPCAM, HNF1B, AQP1, and CFTR co-vary with SOX9 across all 8 donors tested.
3. **PROX1 is a separate cholangiocyte sub-state**, not a disease-progression marker. PROX1-high cells consistently have lower biliary-marker positivity than PROX1-low cells, in both healthy and PSC slices.
4. **Candidate upstream pathways are limited at single-cell resolution.** Only oxidative stress retains directional consistency within donors, and remains borderline. The IL6/STAT3 signal apparent at donor-mean level does not survive cell-level testing.
5. **The SOX9-low state appears consistent with loss of cellular engagement**, not active stress response — senescence, inflammation, and TNF programs co-vary *positively* with SOX9 at single-cell resolution.

---

## Limitations

A condensed list — see `docs/LIMITATIONS.md` for the full version.

* Small PSC cohort (n = 10 donors), with seven from a single source dataset.
* The dominant PSC dataset is single-nucleus RNA-seq; cytoplasmic transcripts are systematically under-detected, attenuating magnitudes.
* No genome-wide FDR establishes single-gene significance at this n; direction concordance and within-donor dose-response carry the inferential weight.
* HVG filtering removed most pathway receptor genes; receptor-side ligand-receptor analysis is currently blocked.
* No IAD samples are present. All IAD claims are by analogy.
* Cause vs effect cannot be inferred from expression alone.

---

## How to reproduce

```bash
# 1. Clone the repository
git clone https://github.com/<your-username>/<repo-name>.git
cd <repo-name>

# 2. Set up Python environment (3.11+)
python -m venv .venv
source .venv/bin/activate                       # macOS / Linux
# .venv\Scripts\activate                        # Windows
pip install -r requirements.txt

# 3. Fetch the data (CELLxGENE Census — network required, ~1 GB)
python src/download_data.py
python src/refetch_cholangiocytes_with_prox1.py

# 4. Run the analyses in order
python src/preprocess_liver.py
python src/trajectory_psc_only.py
python src/tf_activity.py
python src/sox9_mechanism_discovery.py
python src/sox9_per_cell_dose_response.py
python src/per_dataset_replication.py
python src/sensitivity_drop_1873a18a.py
python src/technology_stratified_prox1.py

# 5. Build the polished figures and the 6-page brief
python src/make_figure1_sox9_collapse.py
python src/make_figure2_identity_loss.py
python src/make_figure3_pathways.py
python src/make_study_explainer_pdf.py
```

All intermediate result tables land in `results/`. All figures land in `figures/`. The final PDF brief lands in `reports/study_explainer.pdf`.

---

## Repository structure

```
.
├── README.md                  # this file
├── LICENSE                    # MIT
├── requirements.txt           # Python dependencies
├── .gitignore                 # standard Python + data exclusions
├── PUBLISH_CHECKLIST.md       # what to upload, what to keep local
│
├── data/                      # data is fetched on demand, NOT committed
│   └── README.md
│
├── src/                       # all analysis scripts
│   └── README.md              # run order + what each script does
│
├── notebooks/                 # optional tutorial / exploration notebooks
│   └── README.md
│
├── results/                   # result TSVs + audit trail
│   ├── README.md
│   └── RESULT_NOTE.md         # full step-by-step analysis log
│
├── figures/
│   └── preprint/              # publication-ready figures (PNG + PDF)
│       ├── figure1_sox9_collapse.{png,pdf}
│       ├── figure2_identity_loss.{png,pdf}
│       └── figure3_pathways.{png,pdf}
│
├── reports/
│   ├── study_explainer.pdf    # 6-page polished scientific brief
│   └── WRITEUP_DRAFT.md       # preprint-style prose draft
│
├── docs/
│   ├── PROJECT_SUMMARY.md     # short executive summary
│   ├── METHODS.md             # full methods write-up
│   └── LIMITATIONS.md         # full limitations write-up
│
└── references/
    └── references.md          # 9-citation reference list
```

---

## Portfolio note

This project is built as a data-science / computational-biology portfolio piece. It demonstrates:

* End-to-end single-cell RNA-seq analysis on public human liver data (Python, scanpy, pandas, scipy, matplotlib).
* Donor-controlled pseudobulk DE, trajectory inference, gene-set scoring, and per-cell dose-response.
* A complete self-critique trail: per-dataset replication, sensitivity to the dominant dataset, technology stratification (snRNA-seq vs scRNA-seq), and demotion of weaker candidate findings.
* Publication-quality figures (matplotlib + Adobe-editable PDF output) and a polished 6-page scientific brief.
* Conservative scientific framing throughout — every claim is bounded by the evidence available at this sample size.

What the project is *not*: a clinical study, a mechanistic claim, or a treatment proposal.

---

## Disclaimer

This is exploratory single-cell transcriptomic analysis on public data. **It is not a causal claim, not a clinical study, and not a therapeutic recommendation.** Findings are observational and require experimental validation before any mechanistic interpretation. No identifiable patient data is used or distributed. References to specific genes (SOX9, PROX1, etc.) and pathways do not imply they are drug targets.

For correspondence: open an issue on this repository.

---

*Last updated: see commit history.*
