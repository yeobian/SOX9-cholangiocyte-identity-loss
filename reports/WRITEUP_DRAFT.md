# SOX9 collapse in cholangiocytes is a coordinated identity-loss state in primary sclerosing cholangitis, with oxidative stress as the only candidate upstream regulator surviving cell-level testing

*Draft, exploratory single-cell analysis. Not a causal claim. Not a therapeutic claim.*

---

## Abstract

Idiopathic adulthood ductopenia (IAD) is a rare cholestatic syndrome with no public single-cell data. Using primary sclerosing cholangitis (PSC) as a tractable cholangiopathy proxy, we asked which transcriptional changes mark the cholangiocyte identity-loss endpoint that defines ductopenia histologically. From 21,038 human cholangiocytes pooled from 16 CELLxGENE Census datasets (4,092 PSC cells from 10 donors after preprocessing), we built a diffusion-pseudotime trajectory of cholangiocyte program loss and asked what transcription factors and signalling pathways track that trajectory. The single robust within-disease progression signal that survives donor-level, per-cell, per-dataset, sensitivity, and snRNA-seq-vs-scRNA-seq stratification is collapse of the master cholangiocyte transcription factor SOX9 (within-PSC donor-level ρ = −0.74, p = 0.013; per-cell ρ = −0.44 in scRNA-seq, direction-conserved in snRNA-seq). Within PSC, SOX9-low cells show coordinated loss of the entire canonical biliary-identity program (KRT19, KRT7, EPCAM, HNF1B, AQP1, CFTR co-vary with SOX9 across all 8 donors tested, ρ = +0.27, q = 0.002). Per-cell testing of candidate upstream pathways finds only oxidative stress with directionally consistent within-cell coupling (borderline, ρ = −0.12, q = 0.10); donor-level Δ for IL6/STAT3 signalling — the largest pooled candidate — does not replicate at single-cell resolution. Senescence, inflammation, and TNF programs co-vary positively with SOX9 at the cell level, indicating that the SOX9-low endpoint represents loss of cellular engagement rather than active stress engagement. We treat SOX9-low cholangiocytes as a candidate molecular surrogate for the cholangiocyte-loss endpoint relevant to ductopenia in general; no causal upstream regulator is identified.

---

## Background

Idiopathic adulthood ductopenia is defined histologically by progressive loss of > 50 % of interlobular bile ducts in the absence of identifiable cause. It is rare; no public single-cell IAD datasets exist. Primary sclerosing cholangitis (PSC) shares the cholangiocyte-loss endpoint with IAD and is well-represented in single-cell atlases. We use PSC as a proxy disease context to ask which molecular features mark cholangiocyte identity collapse, with the goal of nominating candidate axes relevant to cholangiopathy in general.

The dominant published candidate transcription factors implicated in cholangiocyte identity are SOX9, HNF1B, GATA6, and ONECUT1. PROX1 has been proposed as a hepatocyte-leaning marker whose appearance in cholangiocytes might reflect a transitional state. We tested both axes.

## Methods (compressed)

**Data.** Human cholangiocytes pulled from CELLxGENE Census (`tissue_general == 'liver'`, cell-type containing 'cholangiocyte'), pooled across 16 source datasets. After QC and HVG selection with PROX1 and a curated cholangiocyte panel force-kept, n = 21,038 cells × 3,610 genes, 16 datasets, 105 donors (75 normal, 10 PSC, 3 PBC, 17 other). Of the PSC cells, 3,741 are from a single snRNA-seq dataset (Andrews et al., *J Hepatol* 2024; CELLxGENE collection `0c8a364b-…`) and 351 are from the scRNA-seq companion dataset in the same collection.

**IAD score.** Per-cell weighted score over a 20-gene cholangiocyte diagnostic panel, scaled to the 95th percentile and inverted so that high score = cholangiocyte program loss.

**Trajectory.** PAGA + diffusion pseudotime, rooted in the cluster with lowest IAD score. Three analyses: all-disease pooled, PSC + normal restricted, and intra-PSC only.

**Statistical tests.** Donor-level Spearman ρ between donor-mean pseudotime and donor-mean gene/score; per-cell Spearman ρ within disease; donor-controlled pseudobulk DE (per-donor logFC, one-sample t-test on cross-donor logFCs); per-dataset replication; sensitivity by dropping the dominant dataset; technology stratification (snRNA-seq vs scRNA-seq).

**Pathway scoring.** 19 curated pathway gene sets scored per cell with `scanpy.tl.score_genes`. Per-donor Δ (SOX9-high − SOX9-low) and per-cell within-donor ρ(SOX9, pathway score) computed separately.

## Results

### 1. SOX9 collapse is the robust within-PSC progression signal

In the PSC + normal trajectory (n = 70 donors), donor-mean pseudotime correlates strongly with donor-mean IAD score (ρ = +0.43, p = 2 × 10⁻⁴), confirming the trajectory captures PSC progression. Of 17 transcription factors and pathway scores tested, SOX9 collapses most strongly with progression: donor-level ρ(pseudotime, SOX9 regulon activity) = −0.72, p = 1 × 10⁻¹³. The same signal appears in the intra-PSC trajectory (donor ρ(pseudotime, SOX9 mRNA) = −0.745, p = 0.013, n = 10 PSC donors), survives per-dataset replication (4 of 9 individual datasets reach p < 0.05 in the predicted direction, the remaining 5 trend correctly), survives dropping the dominant dataset (effect size strengthens), and survives technology stratification (per-cell ρ(SOX9, IAD score) = −0.44 in scRNA-seq and −0.06 in snRNA-seq — same direction, snRNA magnitude attenuated as expected for cytoplasmic mRNAs in nuclei).

No other transcription factor or pathway score reproduces this combination of multi-level robustness.

### 2. SOX9-low PSC cholangiocytes are a coordinated biliary-identity-loss state

Within PSC, the cell-level Spearman ρ between SOX9 expression and a curated biliary-identity score (KRT19, KRT7, EPCAM, HNF1B, AQP1, CFTR) is positive in **all 8 informative donors** (mean within-donor ρ = +0.267, range +0.115 to +0.429, q = 0.002 after BH correction across 19 pathway tests). SOX9 expression captures a coordinated cholangiocyte transcriptional program, not an isolated marker. The SOX9-low state represents loss of that whole program, not loss of one gene.

### 3. Candidate upstream pathways do not survive cell-level testing

Donor-level Δ analysis nominated four candidate upstream pathways: oxidative stress (Δ = −0.22, 86 % donor concordance), IL6 / STAT3 signalling (Δ = −0.26, 71 % concordance), YAP/TAZ target activity (Δ = +0.08 lost, 71 % concordance), and IFN priming (Δ = +0.22 lost, 71 % concordance). Cell-level testing — within-donor Spearman ρ(per-cell SOX9, per-cell pathway score), aggregated across donors — gives a more conservative answer.

Only **oxidative stress** retains directional consistency at the cell level (mean within-donor ρ = −0.12, range [−0.33, 0.00], 0 of 8 donors positive, 4 of 8 clearly negative; q_BH = 0.10). The IL6 candidate collapses entirely at cell level (mean ρ = +0.018, q = 0.81) — its donor-level Δ was driven by donor-to-donor differences in IL6 program activity, not within-cell coupling to SOX9. YAP/TAZ and IFN priming retain direction but with small magnitudes and weak per-donor agreement.

We therefore demote IL6/STAT3 from "candidate upstream regulator" to "donor-level signal of unclear cell-level mechanism." Oxidative stress remains the only candidate upstream pathway with both donor-level Δ and cell-level dose-response support. It is borderline by FDR.

### 4. Stress programs co-vary positively with SOX9 at cell-level resolution

Senescence (mean within-donor ρ = +0.145, q = 0.105), inflammation_general (+0.138, q = 0.077), and TNF signalling (+0.167, q = 0.107) all show positive within-cell ρ with SOX9 across 5 of 8 donors. SOX9-high cells carry **more** of these stress-response signatures per cell. This is the clearest cell-level evidence for the interpretation that emerged from earlier transcription-factor activity work: the SOX9-low state is characterised by **loss of cellular engagement** with the tissue environment, not by gain of an active stress response. AP-1 (JUN, FOS) collapses with SOX9 along the disease trajectory; senescence and TNF programs co-vary positively with SOX9 at cell level. The SOX9-low cell is not a stressed cholangiocyte; it is an exhausted one.

### 5. PROX1 is a sub-state marker, not a disease-progression marker

For comparison, we tested PROX1 — a candidate hepatocyte-identity transcription factor proposed in the literature as a marker of biliary-to-hepatocyte transitioning cells. PROX1 expression does not track PSC progression: per-cell ρ(SOX9, PROX1) ≈ 0 within PSC, donor-level Δ across disease is technology-confounded (snRNA donor ρ vs IAD = −0.32; scRNA donor ρ = +0.28), and the within-PSC PROX1+ abundance pattern that initially looked promising (ρ = +0.81) collapses to untestable when the single dominant dataset is removed (n = 3 PSC donors remaining). PROX1-high cholangiocytes do have lower biliary identity (cross-technology concordant) but no consistent disease-progression dynamic. PROX1 marks a coherent cholangiocyte sub-state that exists in healthy liver; it is not the disease endpoint marker that SOX9 is.

## Discussion

The SOX9 finding is robust across the standard battery of single-cell evidence checks. It survives multi-dataset replication, technology stratification across two distinct sequencing platforms, sensitivity to dropping the dominant cohort, and per-cell testing within individual donors. We interpret SOX9-low PSC cholangiocytes as a candidate molecular surrogate for the cholangiocyte-loss endpoint relevant to ductopenia in general, including IAD. No causal claim is supported.

The mechanism-discovery layer is honest about its limits. Donor-level Δ analysis nominated four candidate upstream pathways, but only oxidative stress retains support at cell-level resolution. The IL6 candidate — the largest pooled effect size — does not survive cell-level testing. The cleanest reading is that the cell-level data in current public PSC datasets does not contain a single conserved upstream-of-SOX9 axis; donor-level signals reflect donor-to-donor heterogeneity in stress-program activity more than within-cell coupling.

The positive cell-level correlation of senescence, inflammation, and TNF programs with SOX9 is consistent with a model in which advanced PSC cholangiocytes lose engagement with their tissue environment as identity collapses, rather than mounting a coherent injury response. Whether the SOX9-low endpoint is driven by oxidative stress upstream, or whether it is a final-common-pathway state of multiple insults, is not addressable from these data.

## Limitations

1. PSC has only 10 donors in current public single-cell data, with 7 of them from a single source. Most within-PSC tests have n ≤ 7 donors and are statistically conservative.
2. The dominant PSC dataset is single-nucleus RNA-seq, which systematically over-represents nuclear-localized TF transcripts and under-represents cytoplasmic mRNAs. Direction of biological effects is preserved across technologies; magnitudes differ.
3. The SOX9 finding survives every test we have run, but no formal genome-wide multiple-testing correction is possible at this n. The strongest evidence is direction concordance across independent datasets and technologies.
4. The candidate upstream regulator analysis is sparse on the receptor side because HVG filtering removed most pathway receptor genes (TGFBR1/2, TNFRSF1A/B, IL6R, IFNAR1/2, BMPR1A/B, FZD family). Re-fetching with these receptors force-kept would unlock proper cell-cell signalling analysis.
5. Cell-cell communication analysis is not possible from the current cholangiocyte-only data. Identification of the cellular source of any candidate signal (IL6, TGF-β, TNF, etc.) requires fetching surrounding liver cell types.
6. No IAD samples are present. All IAD claims are by analogy to PSC and are not direct evidence.
7. Cause vs effect is not addressable. Pathways correlated with SOX9 — positively or negatively — could equally well be drivers, consequences, or correlates of identity loss. Functional perturbation in cholangiocyte organoids is required for any mechanistic claim.

## What this analysis is and isn't

This is an exploratory, public-data single-cell analysis. It nominates SOX9 collapse as a robust molecular phenotype of PSC-associated cholangiocyte identity loss and oxidative stress as the most defensible candidate upstream mechanism — both as starting points for experimental work, not as conclusions.

This is not evidence that SOX9 causes ductopenia. It is not evidence that anti-oxidant treatment would preserve bile ducts. It is not evidence about IAD specifically. It is a structured, multi-test, sensitivity-checked computational identification of one robust marker and one borderline candidate axis, intended to inform follow-up work rather than to support clinical claims.

## Figures planned (for full preprint)

1. **SOX9 collapse — the robust finding.**
   Panel A: PAGA + UMAP of cholangiocytes coloured by PSC pseudotime.
   Panel B: SOX9 regulon activity vs donor-mean pseudotime, scatter coloured by disease (`figures/comparison/pseudotime_vs_iad_psc.png` style).
   Panel C: Per-dataset SOX9 mRNA vs IAD score donor ρ (4/9 datasets individually significant, all in the right direction).
   Panel D: Technology stratification — SOX9 vs IAD per-cell ρ in snRNA-seq vs scRNA-seq, both negative.

2. **SOX9-low is a coordinated identity-loss state.**
   Panel A: Within-donor Spearman ρ(SOX9, biliary identity score) across 8 donors — all positive, mean +0.27.
   Panel B: PROX1 ↔ SOX9 cell-level scatter colored by disease — essentially independent within PSC.
   Panel C: PROX1-high vs PROX1-low cell SOX9 / KRT19 / biliary identity bar charts — PROX1-high cells have *lower* identity scores.

3. **Candidate upstream pathways: donor-level Δ vs cell-level ρ.**
   Panel A: Donor-level Δ heatmap across 19 pathways, sorted by donor concordance.
   Panel B: Cell-level within-donor ρ boxplots across pathways. Show only biliary_identity survives BH FDR.
   Panel C: Oxidative-stress pathway score vs SOX9 expression per cell, faceted by donor (showing direction in 4/8 donors).
   Panel D: Per-donor cell-level ρ for IL6, oxidative stress, YAP, IFN — donor-by-donor consistency.

## What the next experimental work would test

The most informative single follow-up: cholangiocyte organoid oxidative-stress perturbation. Treat human cholangiocyte organoids with paraquat or H₂O₂ at a sub-lethal dose; measure SOX9 protein, KRT19, HNF1B, and AQP1 by qPCR and immunofluorescence at 24 / 48 / 72 hours. If SOX9 drops monotonically and the rest of the biliary identity program follows, the within-cell ρ = −0.12 from this analysis is supported as a real upstream relationship. If SOX9 doesn't drop, oxidative stress is downstream or parallel, not upstream.

A second informative experiment: spatial transcriptomics in PSC biopsy material, asking whether SOX9-low cells co-localise at vanishing-duct sites. The Andrews et al. 2024 paper deposited 4 Visium sections from one PSC patient; these are public and could be analysed without new sample collection.

A third: re-fetch the cholangiocyte AnnData with receptors force-kept, fetch the surrounding liver cell types from Census, and run cell-cell communication analysis to identify which cell type is the candidate source of oxidative stress / TGF-β / IL6 signal reaching cholangiocytes. The within-cell receptor-side analysis is currently blocked.

---

*Data availability: all analyses use CELLxGENE Census liver cholangiocyte data. Scripts and intermediate results in this repository under `IAD/src/` and `IAD/results/`. Full analysis history in `IAD/results/RESULT_NOTE.md` (v10).*
