# `results/`

Curated summary tables and the full analytical audit trail.

## Files

| File | What |
|---|---|
| `RESULT_NOTE.md` | step-by-step analytical audit trail (v10) — how each claim was generated, tested, and either kept or demoted |
| `sox9_pathway_summary_psc.tsv` | donor-level Δ + concordance per pathway, SOX9-high vs SOX9-low in PSC |
| `sox9_per_cell_dose_response_summary.tsv` | within-donor cell-level ρ(SOX9, pathway), per pathway with BH FDR |
| `sox9_receptor_summary_psc.tsv` | receptor-panel donor-level Δ |
| `sox9_threshold_sensitivity_psc.tsv` | pathway Δ across tertile / quartile / quintile SOX9 splits |
| `per_dataset_replication_summary.tsv` | concordance across 9 individual source datasets |
| `technology_stratified_prox1_verdict.tsv` | snRNA-seq vs scRNA-seq verdict per metric |
| `sensitivity_summary.tsv` | which claims survive dropping the dominant cohort |
| `trajectory_psc_per_donor.tsv` | donor-level mean pseudotime + IAD score (PSC + normal) |
| `trajectory_intra_psc_correlations.tsv` | within-PSC cell-level + donor-level correlations |
| `prox1_abundance_summary.tsv` | PROX1+ fraction PSC vs normal (Mann–Whitney) |
| `prox1_sox9_marker_summary.tsv` | PROX1-high vs PROX1-low marker positivity |
| `prox1_high_vs_low_psc_de_top.tsv` | top DE genes between PROX1-high and PROX1-low PSC cholangiocytes |
| `iad_diagnostic_panel.tsv` | the 20-gene cholangiocyte panel used to compute the IAD score |

## Per-cell tables (not committed)

Tables ending in `*_per_cell.tsv` are large (≥100 MB combined), regeneratable from the scripts, and excluded by the project `.gitignore`. To regenerate them, re-run the corresponding script in `src/`.

## How summary tables were produced

Each summary TSV is produced by exactly one script in `src/`. The mapping is documented in the headers of each script. To regenerate any single table, run the script that owns it; full regeneration of all results follows the run order documented in `src/README.md`.
