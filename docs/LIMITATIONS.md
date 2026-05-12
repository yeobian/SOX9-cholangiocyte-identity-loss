# Limitations

This document inventories the limitations of the analysis. Read alongside `METHODS.md` and the headline `README.md`. The whole study is positioned as **exploratory single-cell analysis** — these limitations are why.

## 1. Sample size

The PSC cohort consists of ten donors across all current public single-cell data, with seven of them coming from a single source dataset (Andrews TS et al., *J Hepatol* 2024). Within-PSC inference is statistically conservative throughout. Most within-PSC tests have n ≤ 7 donors after donor-level filtering, which is too small to support genome-wide multiple-testing correction at the per-gene level.

## 2. Technology heterogeneity

The dominant PSC dataset is single-nucleus RNA-seq (snRNA-seq), which systematically over-represents nuclear-localised transcripts (transcription factors, lncRNAs) and under-represents cytoplasmic mRNAs (cytokeratins, secreted-protein mRNAs). The remaining datasets are single-cell RNA-seq (scRNA-seq). Direction of biological effects is preserved across technologies in the headline tests, but magnitudes are systematically attenuated for cytoplasmic transcripts in snRNA-seq. We stratify all technology-sensitive tests; even so, the cohort imbalance (PSC dominated by snRNA-seq, normal dominated by scRNA-seq) constrains the donor-level inference.

## 3. Multiple-testing correction

At this donor count, formal genome-wide false-discovery-rate correction does not establish single-gene significance for any candidate other than SOX9 itself. Inferential weight rests on direction concordance across independent datasets and within-donor cell-level dose-response, not on per-gene FDR. Pathway-level scoring is more tractable and is the level at which most candidate mechanism claims are made.

## 4. Receptor coverage

Highly-variable-gene selection in the source AnnData removed most pathway receptor genes (TGFBR1/2, TNFRSF1A/B, IL6R, IFNAR1/2, BMPR1A/B, FZD family). Receptor-side analysis is therefore currently blocked. A receptor-keeping refetch is documented as a known unblock in `RESULT_NOTE.md`, but has not been executed in the current public version.

## 5. Tissue-level signalling

Cell-cell communication analysis (ligand-receptor inference between cholangiocytes and surrounding cell types) requires the non-cholangiocyte liver compartment — Kupffer cells, hepatic stellate cells, hepatocytes, endothelial cells — none of which is present in the cholangiocyte-only AnnData used here. Identifying the cellular *source* of any candidate signal reaching SOX9-low cholangiocytes therefore requires a separate broader fetch, not yet implemented.

## 6. Disease scope and IAD inference

No idiopathic adulthood ductopenia (IAD) samples are present in any current public single-cell dataset. All claims about IAD in this analysis are by analogy to PSC. PSC is biologically and histologically related to IAD (both feature small-intrahepatic bile duct loss), but the analogy is structural, not direct. Findings should not be extrapolated to IAD without IAD-specific data and clinical context.

## 7. Direction of effect

Pathway scores correlated with SOX9 — positively or negatively, at donor or cell level — cannot distinguish drivers from consequences from co-correlates. For example, oxidative stress co-varying inversely with SOX9 across cells could mean (a) oxidative stress drives SOX9 loss, (b) SOX9 loss permits accumulation of an oxidative-stress signature, or (c) both are downstream of a third upstream factor. Distinguishing between these requires experimental perturbation in cholangiocyte organoids; the analysis here cannot.

## 8. PROX1 framing

PROX1 was the project's initial focal candidate. After replication checks, it was demoted from "candidate disease-progression marker" to "candidate sub-state marker." PROX1+ cholangiocytes are abundant (~45 %) in healthy liver, are spread across all 16 source datasets, and consistently show lower biliary-marker positivity than PROX1− cells. None of these findings constitute a disease-specific or causal claim about PROX1.

## 9. IL6 / STAT3 demotion

The IL6 / STAT3 pathway showed the largest donor-mean Δ (SOX9-high − SOX9-low) of any candidate upstream pathway tested. However, the within-donor cell-level ρ(SOX9, IL6 score) is essentially zero (mean +0.018, p = 0.81 by t-test of per-donor ρs). The donor-mean signal is therefore consistent with donor-to-donor heterogeneity in IL6 program activity, not with within-cell coupling to SOX9. The IL6 candidate is reported transparently but explicitly demoted in the headline interpretation.

## 10. No experimental validation

This is a computational analysis. No wet-lab data has been generated as part of this project. The candidate SOX9 surrogate, the borderline oxidative-stress association, and the failed IL6 candidate are all observational. Experimental follow-up in human cholangiocyte organoids would be the next step before any mechanistic interpretation; this is documented as Priority 1 in the next-steps section of the brief, but not implemented here.

## 11. Reproducibility window

The CELLxGENE Census stable release used for data acquisition is dated; future Census releases may update cell-type ontologies, add new donors, or retire older datasets. The analysis pipeline is reproducible against the version pinned in `requirements.txt`, but exact cell counts and per-dataset numbers may shift if the user runs against a more recent Census release.

## 12. What is *not* claimed

* SOX9 collapse causes ductopenia.
* Any candidate upstream pathway (oxidative stress, IL6, YAP/TAZ, IFN priming) drives SOX9 loss.
* The SOX9-low state is reversible.
* Treatment of any candidate pathway would prevent or reverse cholangiocyte loss.
* IAD shares the same molecular endpoint as PSC.

These are *not* findings of this project. They are limits the project intentionally observes.
