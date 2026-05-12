# SOX9 collapse is the robust PSC-cholangiocyte progression signal; PROX1 rise is not

*Result note — exploratory single-cell analysis. Not a causal claim about IAD.*

*v3 — updated after intra-PSC-only trajectory. PROX1 rise is downgraded
further: it does not survive within-disease testing. SOX9 collapse is
upgraded: it survives within-disease, per-dataset replication, AND the
pooled disease-vs-normal contrast.*

## Summary

We built three nested trajectory analyses on 21,038 human cholangiocytes pooled from 16 CELLxGENE Census datasets (75 healthy + 10 PSC + 3 PBC + 12 other donors):

1. **All-disease pooled trajectory** — donor-level ρ(pseudotime, IAD score) = +0.17 (p = 0.10), borderline.
2. **PSC + normal trajectory** — donor-level ρ = +0.43 (p = 2 × 10⁻⁴), strong.
3. **PSC-only trajectory** — donor-level ρ = +0.37 (sanity check), 3,913 / 4,092 cells on connected component.

Pulling apart the levels, the durable claim is **SOX9 collapse along PSC progression**: it survives the within-PSC trajectory (cell-level ρ = −0.26, p = 3 × 10⁻⁶², donor-level ρ = **−0.74, p = 0.013, n = 10 PSC donors**), it is the most-replicated gene across individual source datasets, and it is the strongest pooled-regulon signal (ρ = −0.72, p = 1 × 10⁻¹³). FOXA2 mRNA also drops within PSC (donor ρ = −0.86, p = 0.002).

**The PROX1 rise that the pooled PSC + normal trajectory suggested does not survive within-PSC testing**: cell-level ρ = −0.04 (technically p = 0.007 because n = 3,913, but the effect size is essentially zero); donor-level ρ = +0.20, p = 0.58 across 10 PSC donors. The pooled signal was driven by heterogeneity *among normal donors* — within healthy cholangiocytes, PROX1 marks a sub-population that scores "advanced" on pseudotime, but once disease has set in PROX1 stops moving. PROX1 is therefore a **candidate marker of cholangiocyte heterogeneity** rather than a marker of PSC progression.

## Robust finding

**SOX9 collapse is the molecular phenotype of PSC-associated cholangiocyte exhaustion.**

| Test | n | ρ | p | Interpretation |
|---|---:|---:|---:|---|
| Pooled disease-vs-normal regulon | 21,038 cells | −0.72 | 1 × 10⁻¹³ | Disease cells have lower SOX9 target activity |
| Pooled PSC + normal trajectory (donor) | 76 donors | −0.72 (regulon) | 1 × 10⁻¹³ | SOX9 drops along pooled trajectory |
| Per-dataset SOX9 mRNA (donor) | 9 datasets | 5/9 correct direction; **all 4 significant negatives** | — | Direction is right, replicates inside multiple independent datasets |
| **Within PSC, cell-level SOX9 mRNA** | 3,913 cells | **−0.262** | **3 × 10⁻⁶²** | SOX9 keeps dropping as PSC progresses |
| **Within PSC, donor-level SOX9 mRNA** | 10 donors | **−0.745** | **0.013** | Donors further along PSC progression have less SOX9 |

This is the strongest claim the data supports. **SOX9 collapse in cholangiocytes is the candidate molecular surrogate of PSC-associated bile-duct identity loss** — and by analogy, a candidate molecular surrogate to look for in IAD biopsies if such material becomes available.

A secondary durable finding: **FOXA2 mRNA also drops within PSC** (cell ρ = −0.25, p < 10⁻⁵⁰; donor ρ = −0.86, p = 0.002, n = 10 donors). Interestingly, FOXA2 was the gene that *failed* per-dataset replication on the disease-vs-normal axis but recovers strongly on the within-PSC progression axis. The two axes are testing different things (threshold crossing vs. continuous progression) and FOXA2 evidently lives on the latter.

## Weakened finding

**PROX1 rise is not a within-PSC progression signal at the per-cell expression level.** A separate per-cell vs per-abundance test (added in v4 — see Final PROX1 section below) gives the full picture.

| Test | n | ρ | p | Interpretation |
|---|---:|---:|---:|---|
| Pooled PSC + normal trajectory (donor PROX1 regulon) | 76 donors | +0.43 | 9 × 10⁻⁵ | Suggested PROX1 rises with disease |
| Per-dataset replication (donor PROX1 mRNA) | 9 datasets | 6/9 correct direction; 3 significant; largest dataset wrong direction | — | Mixed signal at the per-dataset level |
| Per-cell stratified — normal cells | 2,995 cells | +0.273 | 4 × 10⁻⁵² | Within healthy: PROX1 marks "advanced" sub-population |
| Per-cell stratified — PSC cells | 1,186 cells | +0.016 | 0.57 | Within disease: PROX1 doesn't track pseudotime |
| **Within PSC, cell-level PROX1 mRNA** | 3,913 cells | **−0.043** | 0.007 | Tiny negative — practically zero |
| **Within PSC, donor-level PROX1 mRNA** | 10 donors | **+0.20** | 0.58 | Not significant |
| **Within PSC, donor-level PROX1 regulon** | 10 donors | +0.44 | 0.20 | Not significant |

The original PROX1 rise in the pooled trajectory came from pooling normal-donor heterogeneity into the disease axis. PROX1+ healthy cholangiocytes are a real sub-population (45 % of pooled cholangiocytes are PROX1+, spread across 16 datasets), and within healthy cholangiocytes those PROX1+ cells score "advanced" on the trajectory.

Also weakened by within-PSC testing: HNF4A regulon (donor ρ = +0.08, p = 0.83), HNF1B regulon (donor ρ = +0.06, p = 0.88), HEPATOCYTE_id (donor ρ = +0.37, p = 0.29). The "hepatocyte-leaning rescue program" we initially saw in the pooled trajectory is mostly an artifact of pooling.

## Final PROX1 section: per-cell expression vs sub-population abundance

We separated two distinct PROX1 questions that were collapsed in earlier analyses:

1. **Per-cell PROX1 expression level** — how much PROX1 mRNA per cholangiocyte? Tested above; does not track within-PSC progression.
2. **PROX1+ sub-population abundance** — what FRACTION of a donor's cholangiocytes are PROX1+? Tested here.

These are independent quantities. A donor whose PROX1+ sub-population doubles in size while individual cells keep the same PROX1 level shows up as flat per-cell expression but as elevated abundance.

### Per-donor abundance results (61 donors with ≥ 20 cholangiocytes)

| Test | n | Test statistic | Interpretation |
|---|---:|---:|---|
| **PSC vs normal — Mann–Whitney on PROX1+ fraction (loose: PROX1 > 0)** | PSC = 8, normal = 44 | median PSC = 0.40, normal = 0.45, **p = 0.34** | No abundance difference at the disease level |
| **PSC vs normal — Mann–Whitney on PROX1+ fraction (strict: top quartile)** | PSC = 8, normal = 44 | medians 0.22 vs 0.22, p = 0.94 | No abundance difference (strict) |
| **PROX1+ fraction vs IAD score across all donors** | 61 donors | **ρ = +0.231, p = 0.073** | Borderline positive trend |
| **Within PSC: PROX1+ fraction vs PSC pseudotime** | **8 donors** | **ρ = +0.810, p = 0.015** | Strong positive |
| Per-dataset Spearman (PROX1+ frac vs IAD per donor) | 4 datasets | 3/4 ρ > 0; none individually significant | Direction mostly right, underpowered |

The within-PSC abundance result is robust to leave-one-out (every drop-one ρ between +0.71 and +0.93, every p ≤ 0.07). It is not driven by a single donor.

### The major caveat that has to travel with this number

**Seven of the eight PSC donors come from a single source dataset** (`1873a18a-66fd-…`). That dataset also happens to be the one where per-cell PROX1 went in the *opposite* of the predicted direction in the per-dataset replication test. The within-PSC abundance pattern is therefore essentially a *within-one-dataset finding*, not multi-dataset replication. The two facts together can be reconciled: in dataset 1873a18a, advanced-PSC donors have more PROX1+ cholangiocytes than mild-PSC donors (drives the abundance signal), but PROX1-expression-per-cell is lower across all that dataset's donors than in healthy cholangiocyte datasets generally (drives the wrong-direction per-dataset replication). The two coexist if PROX1+ cells in this particular dataset express PROX1 at a lower mean level than PROX1+ cells elsewhere — a dataset-specific normalization or capture quirk.

### Final PROX1 interpretation

- **PROX1+ fraction is not different between PSC and normal cholangiocytes overall.** Roughly 45 % of cholangiocytes are PROX1+ regardless of donor disease status.
- **Across all 61 donors, PROX1+ fraction shows a borderline positive trend with IAD score** (ρ = +0.23, p = 0.07). Direction is right; not statistically secured.
- **Within PSC donors specifically, donors further along the disease pseudotime have substantially more PROX1+ cholangiocytes** (ρ = +0.81, p = 0.015, n = 8). Robust to leave-one-out, but **largely confined to a single source dataset** and therefore not multi-dataset replication.
- **The per-cell PROX1 expression level tells a different story** from the abundance: per-cell PROX1 does *not* track within-PSC progression, but the *fraction of PROX1+ cells* does.
- **Most defensible reading**: PROX1 marks a cholangiocyte sub-state that exists in healthy liver and may *expand selectively in advanced PSC*, while its per-cell expression level does not change. Whether this generalises beyond the one dataset is the next thing to find out.
- **Not defensible**: claiming PROX1 is a marker of disease onset (it isn't — healthy and PSC have similar PROX1+ fractions overall), or that PROX1 causes ductopenia.

## PROX1 ↔ SOX9 association — does PROX1 mark rescued or collapsed cholangiocytes?

Test added in v5. With the abundance result suggesting PROX1+ cells expand within advanced PSC, the natural next question is whether those PROX1+ cells are the ones that *preserve* cholangiocyte identity (a hepatocyte-leaning rescue program protecting the duct) or the ones that have *collapsed* into a SOX9-low intermediate state (failed plasticity).

### Per-cell + per-donor correlations

| Slice | n cells | per-cell ρ(PROX1, SOX9) | p | n donors | donor ρ | donor p |
|---|---:|---:|---:|---:|---:|---:|
| ALL cholangiocytes | 21,038 | +0.061 | 1 × 10⁻¹⁸ | 61 | +0.018 | 0.89 |
| Normal | 13,169 | +0.040 | 6 × 10⁻⁶ | 44 | +0.055 | 0.72 |
| **PSC** | 4,092 | +0.017 | 0.27 | 8 | **−0.833** | **0.010** |

Per-cell ρ is essentially zero everywhere. The strong donor-level PSC ρ = −0.83 comes from 8 donors of which 7 sit in the `1873a18a-66fd-…` dataset, the same single-dataset caveat that's followed PROX1 throughout this analysis — so it cannot be claimed as cross-dataset replication on its own.

### PROX1-high (top quartile) vs PROX1-low contrasts

| Slice | SOX9 mean (high vs low) | BILIARY_id (high vs low) | %KRT19+ | %KRT7+ | %EPCAM+ | %SOX9-high (top Q) |
|---|---|---|---|---|---|---|
| **ALL** | −0.063 vs +0.021 (p = 4 × 10⁻¹⁴) | −0.164 vs +0.055 (**p = 1 × 10⁻⁹⁰**) | 21.4 vs 38.2 % | 41.5 vs 53.7 % | 29.0 vs 42.5 % | 22.8 vs 25.7 % |
| Normal | +0.098 vs +0.170 (p = 2 × 10⁻⁷) | −0.031 vs +0.233 (p = 6 × 10⁻⁸¹) | 25.8 vs 46.5 % | 44.0 vs 59.4 % | 35.6 vs 52.0 % | 24.4 vs 25.2 % |
| **PSC** | −0.539 vs −0.454 (p = 1 × 10⁻⁷) | −0.529 vs −0.444 (p = 4 × 10⁻⁵) | 9.6 vs 16.7 % | 33.1 vs 37.8 % | 9.6 vs 15.2 % | 6.7 vs 12.9 % |

**Same direction in every slice**: PROX1-high cholangiocytes have **lower** SOX9, **lower** biliary-identity score, and **lower** positivity for every classical cholangiocyte marker (KRT19, KRT7, EPCAM, SOX9-high) than PROX1-low cholangiocytes. This holds in normal liver and in PSC.

### Per-dataset

| Dataset | n cells | n donors | per-cell ρ | p |
|---|---:|---:|---:|---:|
| 1873a18a-66fd-… | 6,853 | 13 | +0.075 | < 10⁻³ |
| 1d89d081-f43d-… | 3,558 | 27 | +0.082 | < 10⁻³ |
| 84cfa5aa-a0d1-… | 1,511 | 18 | +0.109 | < 10⁻³ |
| 120cbe1a-9434-… | 1,409 | 15 | +0.100 | < 10⁻³ |
| e3ed2ba4-edf5-… | 524 | 6 | +0.153 | < 10⁻³ |
| 0c9a8cfb-6649-… | 507 | 23 | +0.115 | 9 × 10⁻³ |
| 02792605-4760-… | 1,414 | 16 | −0.035 | 0.19 |
| e84f2780-51e8-… | 657 | 10 | +0.031 | 0.43 |
| d1cbed97-d88f-… | 669 | 10 | +0.002 | 0.97 |

8 of 9 datasets show ρ > 0 at the per-cell level, all very small (0.03–0.15). Direction is consistent but magnitude is trivial — PROX1 and SOX9 are mostly independent at the single-cell level.

### Final interpretation: PROX1+ cholangiocytes are SOX9-LOW, not SOX9-preserved

The cleanest reading of all four tests together:

- **Per-cell correlation is near zero** — PROX1 and SOX9 do not vary together on a single linear axis in any slice. They are not coupled in the simple "more PROX1 = more SOX9" sense.
- **But the threshold-based comparison consistently shows PROX1-high cells have less of every cholangiocyte-identity marker** — SOX9, KRT19, KRT7, EPCAM, and the composite biliary-identity score. This holds in normal cholangiocytes and in PSC cholangiocytes alike.
- **In PSC at the donor level, donors with more PROX1+ cholangiocytes have less SOX9 overall** (ρ = −0.83), although this rests on 8 donors mostly from one dataset.

The pattern fits **PROX1-high + SOX9-low → "failed plasticity / intermediate state"** more than the rescue alternative. PROX1+ cholangiocytes look like cells that have *partially lost* the canonical cholangiocyte identity program rather than cells that are *protecting it*. That is consistent with the abundance result (PROX1+ cells expand in advanced PSC) — the expanding sub-population is a less identity-preserving state, not a rescued one.

This is **not** a causal claim. We have no evidence that PROX1 expression *causes* SOX9 collapse — PROX1 may simply mark a sub-state in which cholangiocyte identity is less well maintained. Equally consistent: cholangiocytes that drift toward an intermediate state happen to retain or gain PROX1 because PROX1 is a hepatocyte-identity TF in liver and the intermediate state has hepatocyte-leaning features. Direction of any mechanistic link cannot be inferred from this analysis.

What the data does *not* support: PROX1 marking an "identity-rescuing" or "regeneration-competent" cholangiocyte sub-state. If anything, the opposite.

## What defines the PROX1-high sub-state? Donor-controlled DE within PSC

Test added in v6. With the PROX1+ sub-state characterised as SOX9-low / biliary-low and abundance-expanding within advanced PSC, the next question is what other genes co-mark this state. The signature could be informative about the biological nature of the failed-plasticity intermediate.

### Method

Per-donor pseudobulk DE: within each PSC donor with ≥ 10 PROX1-high cells AND ≥ 10 PROX1-low cells, compute mean expression in PROX1-high minus mean expression in PROX1-low (log1p space). Aggregate per-gene logFCs across donors with a one-sample t-test against zero, BH-correct. This controls for donor identity by treating each donor's logFC as an independent replicate. Conservative, particularly here because only 7 of 10 PSC donors have enough cells in both groups and 7/7 of those come from a single source dataset (`1873a18a-66fd-…`).

### Headline: most apparent DE is donor-driven, not biology-driven

After BH correction across 3,610 genes, only **PROX1** itself survives q < 0.05 — confirming the threshold split (sanity check) but no other gene clears global FDR. Wilcoxon (cell-level, donor-blind) flags many more genes, but the divergence between cell-level Wilcoxon q-values and donor-controlled t-test p-values is exactly the artifact donor control is meant to reveal: most of the cell-level signal was donor-effect, not PROX1-status.

### Consistent-direction trends (≥ 6/7 donors agree)

Below the FDR threshold there is still interpretable subtle signal, captured by per-donor direction concordance rather than by p-values.

| Direction | Gene | Mean per-donor logFC | Donor concordance | t-test p | Note |
|---|---|---:|---:|---:|---|
| UP | **GATA4** | +0.26 | 7 / 7 | 0.019 | hepatic-endoderm TF |
| UP | SCTR | +0.22 | 7 / 7 | 0.001 | secretin receptor (functional cholangiocyte) |
| UP | CUX2 | +0.20 | 6 / 7 | 0.18 | developmental TF |
| UP | DPF3 | +0.23 | 6 / 7 | 0.08 | chromatin remodeller |
| UP | MIR99AHG | +0.19 | 7 / 7 | 0.009 | lncRNA host |
| UP | SHISA9 | +0.24 | 7 / 7 | 0.016 | developmental |
| DOWN | **KRT19** | −0.19 | 7 / 7 | 0.087 | classical biliary marker |
| DOWN | IL18 | −0.26 | 7 / 7 | 0.008 | inflammatory cytokine |
| DOWN | WIPF1 | −0.27 | 7 / 7 | 0.005 | immune / cytoskeletal |
| DOWN | ANXA1 | −0.20 | 6 / 7 | 0.082 | inflammation |
| DOWN | LAMB3 | −0.20 | 7 / 7 | 0.101 | basement membrane |
| DOWN | EREG | −0.19 | 6 / 7 | 0.34 | growth factor |
| DOWN | ITGA6 | −0.27 | 6 / 7 | 0.16 | adhesion |

Counts of genes with ≥ 6/7 donor concordance and |mean logFC| > 0.05: **155 UP, 306 DOWN**. Magnitudes are small (most |logFC| < 0.3).

### Biological reading (cautious)

The PROX1-high cholangiocyte sub-state within PSC shows:
- **A subtle developmental / early-hepatic-endoderm signature** (GATA4, CUX2, DPF3, SHISA9) — consistent with partial re-engagement of an early-lineage program rather than a clean hepatocyte transdifferentiation.
- **Reduced classical biliary identity** (KRT19 down) — consistent with SOX9 / BILIARY_id results above.
- **Reduced inflammation and stress engagement** (IL18, WIPF1, ANXA1 down).
- **Reduced basement-membrane / adhesion / growth-factor signaling** (LAMB3, ITGA6, EREG down).
- **Retention of a functional cholangiocyte gene** — SCTR (secretin receptor) is up.

This favours **partial dedifferentiation toward an early-lineage / progenitor-like program**, with reduced inflammatory and adhesion engagement, rather than either a clean intermediate-plasticity state with strong hepatocyte induction OR a frank regenerative phenotype. The strongest single hit beyond PROX1 itself is **GATA4**, a hepatic-endoderm TF — consistent with cells re-accessing an early developmental program in advanced disease.

### Caveats that travel with this DE

1. **Within-one-dataset finding.** 7/7 informative donors come from `1873a18a-66fd-…`. The DE is essentially testing PROX1-high vs PROX1-low cholangiocytes within a single source cohort.
2. **Conservative test.** Only PROX1 itself reaches q < 0.05. The ≥ 6/7 donor concordance criterion is a softer filter; some hits could be noise.
3. **Small magnitudes.** Mean per-donor logFC ≤ 0.3 throughout. The PROX1-high sub-state is not a dramatic transcriptional shift — it is a subtle drift on a quiet background.
4. **HVG-restricted gene set.** The cholangiocyte AnnData has 3,610 genes after HVG selection (PROX1 force-kept). Genes that would otherwise be informative may have been filtered out at preprocessing.

## Sensitivity analysis: drop dataset 1873a18a

Test added in v7. Across every PSC-specific PROX1 claim in this note, dataset `1873a18a-66fd-…` has been the dominant contributor — 7 of 8 informative PSC donors come from it. We removed it and re-ran the headline tests to see which claims actually survive.

### Donor counts collapse for PSC, not for normal

| | With 1873a18a | Without 1873a18a |
|---|---:|---:|
| Total donors | 61 | 55 |
| Normal donors | 44 | 44 |
| **PSC donors** | **8** | **3** |

Removing 1873a18a essentially removes the PSC cohort from this analysis. Every test that depends on multiple PSC donors becomes statistically unusable.

### Survives without 1873a18a

| Claim | with | without | Note |
|---|---:|---:|---|
| **PROX1+ abundance correlates with donor IAD score** | ρ = +0.231, p = 0.073 | ρ = +0.278, **p = 0.040** | **survives and tightens** |
| **PROX1-high vs PROX1-low SOX9 (ALL cholangiocytes)** | Δ = −0.084 | Δ = −0.117 | same sign, larger |
| **PROX1-high vs PROX1-low BILIARY_id (ALL)** | Δ = −0.219 | Δ = −0.286 | same sign, larger |
| PROX1-high vs PROX1-low SOX9 (Normal slice) | Δ = −0.072 | Δ = −0.065 | same sign |
| PROX1-high vs PROX1-low BILIARY_id (Normal) | Δ = −0.264 | Δ = −0.268 | same sign |

These are general (non-PSC-specific) statements about the PROX1-high sub-state and they do not depend on 1873a18a. PROX1-high cholangiocytes really do carry lower cholangiocyte-identity expression than PROX1-low cholangiocytes — across healthy donors, regardless of which dataset they come from. And the donor-level abundance trend with IAD score actually strengthens when the dominant dataset is removed, which is reassuring.

### Untestable without 1873a18a

| Claim | with | without |
|---|---:|---:|
| Donor ρ(mean PROX1, mean SOX9) within PSC | **−0.833**, p = 0.010, n = 8 | n = 3 — indeterminate |
| PROX1-high vs low SOX9 within PSC | Δ = −0.085 | n = 3, < 30 cells per group — indeterminate |
| PROX1-high vs low BILIARY_id within PSC | Δ = −0.085 | indeterminate |
| PROX1+ abundance vs PSC pseudotime (within PSC) — not re-run here | ρ = +0.81 | only 3 donors retain pseudotime; pseudotime trajectory itself would have to be rebuilt |
| Within-PSC DE signature (GATA4, KRT19↓, IL18↓, …) — not re-run here | 7/7 donors agreed | only 3 donors retain — DE underpowered |

These are all **within-PSC** claims. Removing 1873a18a leaves the data unable to test them, so we cannot say they replicate elsewhere; we also cannot say they don't. They are 1873a18a-dependent in the current public-data state.

### Tiny effects that flip sign

| Claim | with | without |
|---|---:|---:|
| Per-cell ρ(PROX1, SOX9) within PSC | +0.017 | −0.016 |
| Donor ρ(mean PROX1, mean SOX9) across ALL | −0.004 | +0.031 |

Both effect sizes are essentially zero in both directions; the flip is meaningless.

### What this means for the PROX1 story

**Restating the PROX1 findings honestly after sensitivity testing:**

- **Survives multi-dataset testing**: PROX1-high cholangiocytes have lower SOX9 and lower biliary identity than PROX1-low cholangiocytes. PROX1+ abundance correlates with donor IAD score across all 55 non-`1873a18a` donors (ρ = +0.28, p = 0.04). These are pooled-cohort observations that are not PSC-specific.
- **One-dataset findings (not multi-dataset)**: every PSC-specific claim — abundance expansion within PSC progression, donor-level negative PROX1–SOX9 correlation within PSC, the GATA4/CUX2 dedifferentiation DE signature — currently rests on 1873a18a alone. They are internally consistent and survive a leave-one-donor-out test *within* that dataset, but the dataset cannot replicate itself.

This does not falsify the PSC-specific claims. It defines the limit of the current evidence. To upgrade them from "one-dataset findings" to "replicated findings", we need additional PSC cholangiocyte datasets in the Census, or a non-Census PSC dataset added to the analysis. Until then, the PSC-specific PROX1 picture is *consistent with the failed-plasticity / partial-dedifferentiation hypothesis* but only *demonstrated in one source cohort*.

The robust SOX9-collapse finding (the publishable claim) is unaffected by this sensitivity test because it was already shown to replicate at the per-dataset level (`per_dataset_replication.tsv`: 4 datasets with individually significant SOX9 negative ρ vs IAD score, all in the predicted direction). SOX9 is multi-dataset; PROX1 is currently single-cohort.

## Technology stratification: snRNA-seq vs scRNA-seq

Test added in v8. We identified `1873a18a-66fd-…` as **single-nucleus** RNA-seq (Andrews et al., *J Hepatol* 2024, [10.1016/j.jhep.2023.12.023](https://doi.org/10.1016/j.jhep.2023.12.023); same paper also contributed the scRNA-seq dataset `02792605-…` and 8 Visium spatial datasets to the same Census collection). snRNA-seq systematically over-represents nuclear-localised TF transcripts (PROX1, GATA4) and under-represents cytoplasmic mRNAs (KRT19, KRT7, SOX9). Since `1873a18a` carried most of the PSC cells in our pooled cholangiocyte AnnData, every PSC-specific PROX1 claim is potentially technology-confounded. This test runs the headline metrics within snRNA only and within scRNA only.

### The cohort confound

| | normal | PBC | PSC | other |
|---|---:|---:|---:|---:|
| **snRNA cohort (1873a18a, 6,853 cells, 12 donors)** | 781 | 2,331 | **3,741** | 0 |
| **scRNA cohort (15 datasets, 14,185 cells, 55 donors)** | 12,388 | 881 | 351 | 565 |

Disease and technology are heavily confounded. Almost all PSC cholangiocytes live in snRNA; almost all normal cholangiocytes live in scRNA. Any pooled donor-level correlation between PROX1 and disease severity is partly testing snRNA-vs-scRNA effects.

### What survives the technology stratification (concordant direction in both)

| Claim | snRNA | scRNA | Verdict |
|---|---:|---:|---|
| **SOX9 vs IAD (per-cell ρ)** | **−0.06** | **−0.436** | **CONCORDANT** — direction preserved, scRNA magnitude larger (cytoplasmic mRNA detection) |
| PROX1-high vs PROX1-low BILIARY_id Δ (ALL) | −0.031 | **−0.281** | same sign, magnitude differs 9× |
| PROX1-high vs PROX1-low BILIARY_id Δ (Normal) | −0.007 | **−0.263** | same sign, scRNA much stronger |
| cell ρ(PROX1, SOX9) ALL | +0.075 | +0.022 | same sign, both small |
| cell ρ(PROX1, SOX9) Normal | +0.084 | +0.045 | same sign, both small |

The **SOX9 collapse finding is the most robust:** direction is preserved across two independent technologies covering largely non-overlapping disease cohorts. snRNA attenuates the magnitude (typical pattern for cytoplasmic mRNAs in nuclei), but doesn't flip the sign. This is now the strongest single piece of evidence for the publishable claim.

The "PROX1-high cells have lower biliary identity" claim *also* survives at the BILIARY_id-score level, with snRNA showing a flat-but-correct direction and scRNA showing a strong effect.

### What flips sign by technology (technology-confounded, not robust biology)

| Claim | snRNA | scRNA | Verdict |
|---|---:|---:|---|
| **PROX1+ abundance vs donor IAD score** | **−0.32** | **+0.28** | **DISCORDANT — opposite signs** |
| PROX1-high vs PROX1-low SOX9 Δ (ALL) | +0.03 | −0.115 | DISCORDANT |
| PROX1-high vs PROX1-low SOX9 Δ (Normal) | +0.034 | −0.064 | DISCORDANT |
| cell ρ(PROX1, SOX9) PSC | +0.021 | −0.016 | DISCORDANT (both ≈ 0) |

The pooled "PROX1+ abundance correlates with IAD score across donors" we reported earlier (ρ = +0.23 with `1873a18a`, +0.28 without) is now revealed to be averaging over two technologies that point in opposite directions. In the snRNA cohort (most PSC donors), donors with worse IAD score have **fewer** PROX1+ cholangiocytes. In the scRNA cohort (mostly healthy donors with mild variation in IAD), donors with worse IAD score have **more** PROX1+ cholangiocytes. We cannot disentangle "PROX1+ abundance vs disease severity" from "PROX1+ abundance vs technology" with the current data.

The "PROX1-high cells have lower SOX9" claim flips sign by technology — present in scRNA, absent or reversed in snRNA. Likely an snRNA artefact: nuclear PROX1 is over-detected while cytoplasmic SOX9 is under-detected, mechanically decoupling the two.

### What cannot be tested (n too small in one cohort)

The within-PSC PROX1 abundance pseudotime claim (ρ = +0.81), the donor-level within-PSC PROX1↔SOX9 (ρ = −0.83), and the GATA4-led DE signature all live in the snRNA cohort. The scRNA cohort has only 3 PSC donors and 351 PSC cells — not enough to replicate any of those tests at matched statistical power. They remain single-cohort, single-technology findings.

### Reframing the PROX1 story (final, for now)

1. **SOX9 collapse along cholangiocyte program loss** is the robust biology. It survives within-PSC trajectory, per-dataset replication, sensitivity to dropping the dominant dataset, AND now technology stratification. It is the publishable claim.
2. **PROX1-high cholangiocytes have somewhat lower biliary identity** — direction-concordant across technologies (snRNA Δ = −0.03 to −0.007, scRNA Δ = −0.26 to −0.28). The effect-size discrepancy is consistent with snRNA's tendency to under-detect cytoplasmic biliary markers, not with a flip.
3. **The PROX1+ abundance ↔ disease severity relationship is technology-confounded.** The pooled donor-level correlation was an artefact of averaging two technologies that disagree in direction. With the current Census data we cannot resolve this — the snRNA cohort holds nearly all the disease signal, and the scRNA cohort has too few PSC donors to compare.
4. **All within-PSC PROX1 fine-structure findings (abundance pseudotime, donor PROX1↔SOX9, GATA4 DE signature) are single-cohort, single-technology.** They are not falsified, but they cannot currently be replicated.

The next move to firm up the PROX1 story would be a within-paper replication: re-run the within-PSC PROX1 tests separately on `02792605` (the scRNA-seq companion from the same paper, same cohort, same patients). That would be the cleanest comparison because the *patients* are matched and only the *technology* differs. The scRNA-seq companion has 351 PSC cells from 3 donors — borderline for donor-level tests, but enough for cell-level tests.

## Mechanism discovery: what lies upstream of SOX9 collapse?

Test added in v9. Framing: SOX9-low cholangiocytes are treated as a *disease endpoint*, not a cause. We use SOX9-high vs SOX9-low PSC cholangiocytes to nominate candidate upstream pathways. No causal claim, no therapy claim — pathways enriched or depleted in SOX9-low are *candidates that may suggest rescue axes for future validation*.

### Method

Within PSC cholangiocytes (n = 4,092 cells from 10 donors; 7 with cells in both arms), rank-based tertile split (zero-inflation-safe, since most cholangiocytes have SOX9 = 0). Per-donor pseudobulk Δ (SOX9-high − SOX9-low) for 19 curated pathway gene sets; one-sample t-test across donors with BH FDR. Rank-based threshold sensitivity (tertile, quartile, quintile). Per-cell pathway scoring with `sc.tl.score_genes`. Three-category interpretation: A lost maintenance / B gained disease-stress / C candidate upstream regulators.

### What survives donor concordance ≥ 71 %

| Direction | Pathway | Mean Δ (high − low) | Donor agree | Category |
|---|---|---:|---:|---|
| **DOWN in SOX9-low** | biliary_identity | +0.109 | 86 % | A (sanity check, SOX9 in panel) |
| **DOWN in SOX9-low** | ifn_signalling | +0.216 | 71 % | A (lost immune priming) |
| **DOWN in SOX9-low** | senescence | +0.087 | 86 % | A (SOX9-high cells still senescence-engaged) |
| **DOWN in SOX9-low** | yap_taz_targets | +0.084 | 71 % | A (lost mechano-signal) |
| **UP in SOX9-low** | il6_signalling | −0.257 | 71 % | B (strongest gain magnitude) |
| **UP in SOX9-low** | oxidative_stress | −0.224 | 86 % | B (highest concordance gain) |

None reach BH q < 0.1 — pseudobulk DE with 7 donors across 19 pathway tests is conservative under multiple-testing correction. Donor concordance and effect-size patterns are the more interpretable filters at this n.

### Category C — Candidate upstream regulators, ranked by evidence

Ranking criterion = donor concordance × |effect size|. All claims are *candidates*; none are causal.

1. **Oxidative stress / Nrf2 axis** — 86 % donor concordance, Δ = −0.22. SOX9-low cholangiocytes have elevated HMOX1, NQO1, GPX1/4, SOD2, TXN, NFE2L2 expression. Candidate rescue: antioxidant or Nrf2-pathway support. Strongest combination of concordance and magnitude in the gained-program tier.
2. **IL6 / STAT3 signaling** — 71 % concordance, Δ = −0.26 (largest gain magnitude). SOX9-low cells express more IL6, IL6R, IL6ST, STAT3, SOCS3, OSM, OSMR. Clinically tractable rescue: IL6R blockade (tocilizumab class).
3. **YAP/TAZ activity loss** — 71 % concordance, Δ = +0.08 (small magnitude). CTGF, CYR61, ANKRD1, AMOTL2, BIRC5 lower in SOX9-low. Mechano-signalling collapse as a candidate maintenance program.
4. **Loss of IFN priming** — 71 % concordance, Δ = +0.22. IFNAR1/2, IFNGR1/2, STAT1/2, IRF1/3/7, ISG15, IFI6, IFITM3, MX1 lower in SOX9-low. Interpretation is ambiguous: could reflect loss of antiviral / immune-priming engagement, or loss of cholangiocyte responsiveness to inflammatory tissue context.

### What didn't change

TGF-β signalling, BMP, WNT, Notch (both signalling and targets), TNF, bile-acid stress, ER stress / UPR, EMT, fibrosis program, inflammation_general — all donor concordance < 71 % or |Δ| < 0.05. Notably:
- **TGF-β** showed no direction (Δ = −0.001) — surprising given how often it's invoked in cholangiopathy literature. May reflect HVG filter losing TGFBR1/2 + downstream effectors (only 3/10 panel genes present).
- **BMP** also no direction (Δ = −0.08, 57 % concordance) — TGF-β family overall is silent at this analysis level.
- **EMT** Δ = −0.012 — SOX9-low cholangiocytes don't look like overtly EMT cells in this data.

### Receptor panel — mostly blocked by HVG filtering

The receptor proxy for ligand-receptor analysis was largely silenced by HVG filtering: TGFBR1/2, TNFRSF1A/B, IL6R, IFNAR1/2, BMPR1A/B/2, FZD1–10 are all absent from the cholangiocyte AnnData (the HVG step that kept PROX1 didn't keep these). Only Notch_receptors (NOTCH1, NOTCH2), Notch_ligands_on_self (JAG1, JAG2), and Hedgehog (PTCH1, SMO) scored.

To unblock this: add receptor list to the force-keep set in `refetch_cholangiocytes_with_prox1.py` (~30 genes) and re-fetch. One-line code addition. Required before proper cell-cell signaling analysis.

### Limits to the inference

- 7–10 PSC donors total; 7 in the test arms after MIN_CELLS_PER_GROUP_PER_DONOR filtering.
- Almost all PSC cells (3,741 / 4,092) are snRNA-seq from `1873a18a`. This finding is essentially within-one-paper, within-one-technology.
- Pathway scores are sensitive to gene-set composition. Where panels overlap (IL6 appears in both `il6_signalling` and `senescence` panels), Δ values for those pathways are partially correlated.
- Donor-controlled pseudobulk is conservative. Some real biology may sit below the donor-concordance threshold.
- No proper cell-cell communication analysis is possible from the current data (cholangiocyte-only AnnData, no surrounding cell types).
- Direction of mechanistic influence cannot be inferred from expression alone. Pathway UP in SOX9-low could mean "this pathway drives identity loss" or "this pathway responds to identity loss."

### Suggested next steps for mechanism discovery

1. **Re-fetch cholangiocyte AnnData with receptor force-keep list.** Add TGFBR1/2/3, TNFRSF1A/B/10A/B, IL6R, IL6ST, OSMR, IFNAR1/2, IFNGR1/2, BMPR1A/B/2, ACVR1/2A, FZD1–10, LRP5/6, NOTCH1–4, JAG1/2, DLL1/3/4 to `PROX1_PANEL`. Re-run this script. Should unlock the receptor analysis.
2. **Fetch full liver cells (not just cholangiocytes) from Census** to enable proper cell-cell communication analysis. Score expression of IL6, TNF, TGF-β1, JAG1, etc. in surrounding cell types (Kupffer, stellate, hepatocytes, endothelial). Ask: which cell types are the candidate source of the IL6 signal that SOX9-low cholangiocytes are responding to?
3. **Test the candidate upstream pathways for cell-level dose–response.** For each candidate (IL6, oxidative stress, YAP loss), correlate per-cell pathway score with per-cell SOX9 within each donor. Does SOX9 drop monotonically as IL6 signalling rises within individual cholangiocytes? That's a stronger test than donor-mean Δ.
4. **Repeat the whole analysis within PBC.** Only 3 PBC donors in current data — likely underpowered, but worth confirming direction.

### One-paragraph publishable interpretation (cautious)

In PSC cholangiocytes, the SOX9-low transcriptional state — the candidate molecular endpoint of cholangiocyte identity loss — is associated with elevated oxidative-stress and IL6/STAT3 signalling programs, and with reduced YAP/TAZ activity and interferon-priming programs, relative to SOX9-high cholangiocytes within the same donors. None of these associations reach formal multiple-testing significance at n = 7 donors, but the cross-donor direction is concordant in 71 – 86 % of donors for each. These pathways nominate candidate upstream mechanisms — oxidative stress, IL6/STAT3, and mechano-signalling — that may contribute to or accompany cholangiocyte identity collapse in PSC. Whether they are drivers or consequences is not addressable from these data; the candidates would require receptor-resolved cell-cell signalling analysis (currently blocked by HVG-filtered receptor genes) and functional validation in cholangiocyte organoids for any mechanistic claim. By analogy, these are candidate axes for understanding cholangiocyte identity loss in ductopenia in general, including IAD — but no IAD samples are in the current data and no IAD claim is supported.

## Cell-level dose-response: which candidates survive the sharper test?

Added in v10. The donor-level Δ analysis (SOX9-high − SOX9-low pathway score) nominated candidate upstream pathways with n = 7 donors. A much sharper test: within each donor, compute Spearman ρ(per-cell SOX9, per-cell pathway score) across that donor's cells (hundreds to thousands). The within-donor ρ uses cell-level statistical power and is automatically donor-controlled. Cross-donor t-test on the per-donor ρs.

### Within-donor ρ(SOX9, pathway), 8 PSC donors, sorted by mean ρ

| Pathway | Mean ρ | Per-donor range | Donor concordance | q_BH | Category |
|---|---:|---|---|---:|---|
| **biliary_identity** | **+0.267** | [+0.12, +0.43] | 8/8 positive | **0.002** | **A — robust** |
| tnf_signalling | +0.167 | [−0.01, +0.35] | 5/8 positive | 0.107 | ambiguous |
| senescence | +0.145 | [+0.02, +0.41] | 5/8 positive | 0.105 | ambiguous |
| inflammation_general | +0.138 | [+0.01, +0.31] | 5/8 positive | 0.077 | ambiguous |
| ifn_signalling | +0.117 | [−0.14, +0.59] | mixed | 0.485 | ambiguous |
| yap_taz_targets | +0.064 | [−0.18, +0.31] | mixed | 0.603 | ambiguous |
| il6_signalling | +0.018 | [−0.17, +0.13] | mixed | 0.805 | ambiguous |
| (others near 0) | | | | | |
| **oxidative_stress** | **−0.120** | [−0.33, +0.00] | 4/8 negative, 0 positive | 0.096 | borderline B |

### What changes from the donor-level Δ analysis

The donor-level Δ analysis nominated IL6/STAT3, oxidative stress, YAP loss, and IFN-priming loss as candidate upstream mechanisms. The cell-level test forces a more conservative reading:

- **IL6 signalling** — the largest donor-level Δ (−0.26) and our #1 ranked candidate — has within-cell ρ = **+0.018**. Not a within-cell dose-response. The donor-level Δ reflected donor-to-donor differences in IL6 program activity, not within-cell coupling. Demoted from #1 candidate to inconclusive.
- **Oxidative stress** — the only donor-level candidate whose direction survives at cell level. ρ = −0.12 (borderline q = 0.10). Present in half of donors. Still the most defensible *single* mechanism candidate.
- **YAP/TAZ targets, IFN priming** — directions consistent with donor-level Δ but magnitudes small and per-donor agreement weak.
- **Biliary identity (+0.27, 8/8 donors, q = 0.002)** — robust at cell level. This is the identity-axis coherence: SOX9 expression captures a coordinated cholangiocyte-identity program that drops together in SOX9-low cells. Important sanity check, not a candidate upstream regulator.

### A new insight: stress programs co-vary *positively* with SOX9 at cell level

Senescence (ρ = +0.145), inflammation (+0.138), and TNF (+0.167) all have positive within-cell ρ with SOX9 in ≥ 5 of 8 donors. SOX9-high cells carry *more* of these stress-response signatures per cell within donors. This is the clearest cell-level evidence yet for the "SOX9-high cells are still stress-engaged; SOX9-low cells have given up" interpretation: as cells lose SOX9 they don't gain a stress program — they lose all programs including the stress response. Consistent with the AP-1 collapse seen earlier in the TF activity work.

### Revised mechanism reading

At cell-level resolution within PSC cholangiocytes:

1. **Coherent biliary-identity collapse** (KRT19, KRT7, EPCAM, AQP1, CFTR, HNF1B all drop with SOX9 within cells, 8/8 donors). The SOX9-low state is a *coordinated identity-loss state*, not a single-gene drop.
2. **No clean, multi-donor candidate upstream pathway** at cell level. Oxidative stress is the borderline survivor.
3. **Stress programs (senescence, inflammation, TNF) actually co-vary with SOX9 positively at cell level** — the SOX9-low state is characterised by *loss of engagement*, not gain of a specific stress signature. Cells that have lost SOX9 have also lost responsiveness.

### The publishable mechanism statement, with caveats

In PSC cholangiocytes, SOX9-low cells represent a coordinated biliary-identity collapse — the rest of the canonical cholangiocyte transcriptional program (KRT19, KRT7, EPCAM, AQP1, HNF1B, CFTR) co-varies with SOX9 across individual cells in all 8 donors tested. Per-cell analysis does not identify any single upstream pathway that consistently couples to SOX9 across donors; oxidative stress provides borderline within-cell evidence (mean ρ = −0.12) but is donor-inconsistent. Donor-level Δ for IL6/STAT3 signalling was the largest in our pooled summary but does not replicate at the per-cell level. Stress and inflammation pathways co-vary positively with SOX9 at the cell level — i.e. SOX9-low cells are not enriched in active stress response — supporting a model in which advanced PSC cholangiocytes lose engagement with their tissue environment rather than mounting a coherent injury response. Whether the SOX9-low endpoint is driven by oxidative stress upstream, or whether it is a final-common-pathway exhaustion state of multiple insults, is not addressable from these data and would require experimental perturbation in cholangiocyte organoids.

## Remaining open hypotheses

**H-residual-1 (defensible).** Within PSC progression, a sub-population of cholangiocytes loses SOX9 and FOXA2 mRNA progressively. This sub-population is the candidate "exhausted cholangiocyte" implicated in failed bile-duct regeneration. Whether the same population exists in IAD is not testable in any current public dataset.

**H-residual-2 (speculative).** WNT pathway shows a Simpson's-paradox-like split: cell-level ρ = −0.20 within PSC (drops), donor-level ρ = +0.64 (rises, p = 0.048). This means within each donor WNT activity drops along pseudotime, but donors at higher mean pseudotime have higher mean WNT. Possible explanation: regenerative WNT signaling is donor-specific and predicts which donors progress further, even though within an individual donor WNT signaling is being lost. This is interesting enough to flag but too speculative to claim.

**H-residual-3 (no longer well supported).** PROX1 as a disease-progression marker. The data does not support this directly. PROX1 may still be a marker of a hepatocyte-leaning cholangiocyte sub-state that exists in healthy liver and is differentially abundant or activatable in disease, but proving that requires separating "subpopulation abundance" from "per-cell expression" — which the present analysis doesn't do.

## Next validation steps (in order of cost)

1. **Compute the abundance of PROX1+ cholangiocytes per donor**, separately from per-cell PROX1 levels. If PSC donors have a higher *fraction* of PROX1+ cholangiocytes (rather than higher mean PROX1 per cholangiocyte), the heterogeneity-marker hypothesis becomes testable. ~30 lines of code.
2. **Investigate the outlier dataset** (`1873a18a…`, the largest, where PROX1 went the wrong way at the per-dataset level). Identify the original Census collection, donor demographics, dissociation protocol. If it's a particular study choice that systematically biases against PROX1+ cholangiocytes, that affects all downstream interpretation.
3. **Pseudobulk DE between SOX9-high and SOX9-low cholangiocytes within PSC.** Cells already have SOX9 expression measured; splitting them by tertile and running DE within PSC donors would identify what else changes when SOX9 collapses. This is the in-silico equivalent of "what does an exhausted cholangiocyte look like."
4. **Spatial validation on PSC biopsy.** Spatial transcriptomics asking whether SOX9-low cells co-localize at vanishing-duct sites. This is the highest-value experiment but requires clinical material and a collaborator. The GW liver transplant program is the natural starting point.
5. **In vitro causal test.** SOX9 knockdown in human cholangiocyte organoids ± inflammatory cytokines; measure identity loss, proliferation, and apoptosis. Smallest experiment that could falsify the SOX9-collapse mechanism.

## Mechanism diagram (in words) — current best version

Healthy cholangiocyte: SOX9+ FOXA2+ HNF1B+ GATA6+ ONECUT1+, biliary identity high, AP-1 responsive. **Embedded heterogeneity**: roughly half of healthy cholangiocytes additionally express PROX1, marking what looks like a hepatocyte-leaning sub-state present in normal liver.

→ PSC progression →

PSC cholangiocyte: **SOX9 collapses progressively** (the robust within-disease signal). FOXA2 drops in parallel. AP-1 stress program also collapses (cells stop responding rather than continuously responding). HNF1B / HNF4A / hepatocyte-identity scores do *not* clearly rise within PSC — the apparent rise in the pooled analysis was a pooling artifact. PROX1 abundance does not track within-disease progression.

The cell that the data describes is **a cholangiocyte losing SOX9-driven identity, not a cholangiocyte transdifferentiating into a hepatocyte**.

## Limitations (carried forward + new)

- 21,038 cells across 16 datasets gives statistical power but also batch confounding. We have not run formal batch correction (Harmony / scVI / Scanorama).
- "Cholangiocyte" depends on the Census ontology label. PROX1+ KRT19+ HNF4A+ cells may be transitioning cells, hybrid cells, or misannotated hepatocytes.
- Spearman on log-normalized scRNA-seq is dominated by tied zeros. Reported p-values are valid as ordinal tests but should not be interpreted as effect-size estimates.
- The PSC subset has 10 donors; donor-level Spearman tests are underpowered (e.g. PROX1 donor ρ = +0.44 in the within-PSC trajectory has p = 0.20 — direction not testable).
- Diffusion pseudotime depends on the kNN graph parameters. The within-PSC trajectory used n_neighbors = 30 (vs 15 in the pooled one) because the smaller PSC subset needed a more connected graph; results are robust to this choice but not unique.
- The largest source dataset shows PROX1 against the predicted direction. Worth investigating before any PROX1 claim is published.
- IAD itself is not represented in this dataset. All IAD claims are by analogy to PSC.

## Reproducibility

Run from `IAD/`:

```
python src/refetch_cholangiocytes_with_prox1.py
python src/trajectory_psc_only.py
python src/tf_activity.py
python src/pseudotime_tf_crossplot.py
python src/prox1_dataset_spread.py
python src/prox1_per_cell_pseudotime.py
python src/per_dataset_replication.py
python src/trajectory_intra_psc.py
python src/prox1_abundance_test.py
python src/prox1_sox9_association.py
python src/prox1_high_vs_low_psc_de.py
python src/sensitivity_drop_1873a18a.py
python src/technology_stratified_prox1.py
python src/sox9_mechanism_discovery.py
python src/sox9_per_cell_dose_response.py # ← new (this update)
```

Outputs (bold = new this update):
- `results/trajectory_psc_per_cell.tsv`, `results/trajectory_psc_per_donor.tsv`
- `results/pseudotime_vs_tf_activity.tsv`
- `results/tf_activity_disease_summary.tsv`, `results/tf_activity_per_donor.tsv`
- `results/prox1_iad_link.tsv`, `results/prox1_iad_link_per_donor.tsv`
- `results/prox1_dataset_spread.tsv`
- `results/prox1_pseudotime_per_cell.tsv`, `results/prox1_pseudotime_per_disease.tsv`
- `results/per_dataset_replication.tsv`, `results/per_dataset_replication_summary.tsv`
- **`results/trajectory_intra_psc_correlations.tsv`**
- **`results/trajectory_intra_psc_per_cell.tsv`**
- **`results/trajectory_intra_psc_per_donor.tsv`**
- **`figures/comparison/intra_psc_pseudotime_vs_SOX9.png`**
- **`figures/comparison/intra_psc_pseudotime_vs_HNF1B.png`**
- **`figures/comparison/intra_psc_pseudotime_vs_PROX1.png`**
- **`figures/comparison/intra_psc_correlation_heatmap.png`**
- **`results/prox1_abundance_per_donor.tsv`**
- **`results/prox1_abundance_summary.tsv`**
- **`results/prox1_abundance_per_dataset.tsv`**
- **`figures/prox1/prox1_abundance_by_disease_boxplot.png`**
- **`figures/prox1/prox1_abundance_vs_iad_scatter.png`**
- **`figures/prox1/prox1_abundance_vs_psc_pseudotime.png`**
- **`figures/prox1/prox1_abundance_per_dataset_dotplot.png`**
- **`results/prox1_sox9_per_cell_corr.tsv`**
- **`results/prox1_sox9_per_donor.tsv`**
- **`results/prox1_sox9_marker_summary.tsv`**
- **`results/prox1_sox9_per_dataset.tsv`**
- **`figures/prox1/prox1_vs_sox9_scatter.png`**
- **`figures/prox1/prox1_high_low_sox9_boxplot.png`**
- **`figures/prox1/prox1_high_low_marker_bars.png`**
- **`figures/prox1/prox1_sox9_per_dataset_heatmap.png`**
- **`results/prox1_high_vs_low_psc_de.tsv`**
- **`results/prox1_high_vs_low_psc_de_top.tsv`**
- **`results/prox1_high_vs_low_psc_per_donor.tsv`**
- **`figures/prox1/prox1_high_low_de_volcano.png`**
- **`figures/prox1/prox1_high_low_de_top_genes.png`**
- **`results/sensitivity_drop_1873a18a.tsv`**
- **`results/sensitivity_summary.tsv`**
- **`figures/prox1/sensitivity_drop_1873a18a.png`**
- **`results/technology_stratified_prox1.tsv`**
- **`results/technology_stratified_prox1_verdict.tsv`**
- **`figures/prox1/technology_stratified_prox1.png`**
- **`results/sox9_de_pseudobulk_psc.tsv`**
- **`results/sox9_pathway_donor_contrast_psc.tsv`**
- **`results/sox9_pathway_summary_psc.tsv`**
- **`results/sox9_receptor_donor_contrast_psc.tsv`**
- **`results/sox9_receptor_summary_psc.tsv`**
- **`results/sox9_threshold_sensitivity_psc.tsv`**
- **`figures/sox9/sox9_pathway_summary.png`**
- **`figures/sox9/sox9_de_volcano.png`**
- **`figures/sox9/sox9_receptor_heatmap.png`**
- **`figures/sox9/sox9_threshold_sensitivity.png`**
- **`results/sox9_per_cell_dose_response.tsv`**
- **`results/sox9_per_cell_dose_response_summary.tsv`**
- **`figures/sox9/sox9_per_cell_dose_response.png`**
