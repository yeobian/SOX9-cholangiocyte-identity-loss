# Project summary

## One-sentence summary

Across 21,038 cholangiocytes from 16 public single-cell datasets, SOX9 collapse is the most consistent within-disease transcriptional signal in primary sclerosing cholangitis (PSC); we treat SOX9-low cells as a candidate molecular surrogate for the cholangiocyte-loss endpoint relevant to idiopathic adulthood ductopenia (IAD).

## Motivation

Idiopathic adulthood ductopenia is a rare cholestatic syndrome marked by progressive loss of intrahepatic bile ducts. Because no public single-cell IAD datasets exist, the question of what drives cholangiocyte identity loss has no direct molecular handle. PSC is biologically and histologically the closest cholangiopathy with public single-cell data, and we use it as a proxy.

## Question

What transcriptional changes — gene-level and pathway-level — most consistently mark the cholangiocyte-identity-loss endpoint that defines ductopenia histologically?

## Approach

1. Pool cholangiocytes from all available CELLxGENE Census datasets (16 sources, 105 donors).
2. Construct a diffusion-pseudotime trajectory rooted in the most-healthy cluster.
3. Identify transcription factors and pathway scores whose activity tracks the trajectory.
4. Stress-test every candidate finding with per-dataset replication, technology stratification, sensitivity to dropping the dominant cohort, and per-cell within-donor dose-response.

## Findings

* **SOX9 collapse** is consistent in donor-level Spearman tests within PSC (ρ = −0.74, p = 0.013, n = 10 donors), in 4 of 9 individual source datasets at p < 0.05 (all in the predicted direction), and in both snRNA-seq and scRNA-seq.
* **The SOX9-low state is coordinated**, not single-gene. Within all 8 PSC donors tested, SOX9 expression positively covaries with KRT19, KRT7, EPCAM, HNF1B, AQP1, and CFTR at the cell level (mean ρ = +0.27, q = 0.002).
* **PROX1 marks a separate cholangiocyte sub-state** with lower biliary identity, but does not track disease progression. PROX1+ cells are abundant in healthy liver and have lower biliary-marker positivity than PROX1− cells.
* **Candidate upstream pathways are limited.** Oxidative stress retains direction at single-cell resolution (mean within-donor ρ = −0.12) but only borderline. IL6 / STAT3 — the strongest donor-mean candidate — does not survive cell-level testing.
* **The SOX9-low endpoint appears consistent with loss of cellular engagement** rather than active stress: senescence, TNF, and inflammation programs co-vary positively with SOX9 at single-cell resolution.

## Defensible scope

This is exploratory single-cell analysis on public data. The repository documents, in addition to the findings, the analyses that *failed* — for example PROX1 was the initial focal candidate and was demoted to a sub-state marker after replication. The headline claim is bounded: SOX9-low cholangiocytes are a *candidate molecular surrogate* for the cholangiocyte-loss endpoint, not a proven cause and not a treatment target.

## What this project demonstrates

A complete, donor-controlled, sensitivity-checked single-cell analysis with: clean Python scripts, publication-quality figures, a polished six-page brief, and a full audit trail of how the strongest claim survived and which weaker claims did not.

## What this project does not claim

* That SOX9 *causes* ductopenia.
* That any of the candidate upstream pathways *drive* SOX9 collapse.
* That the SOX9-low state is reversible or treatable.
* That findings extrapolate to clinical IAD without further validation.

See `LIMITATIONS.md` for the full set of caveats.
