# `notebooks/`

Reserved for future tutorial and exploration notebooks.

## Current status

No notebooks committed yet. The analysis is currently implemented as Python scripts in `src/`, which is intentional: scripts are easier to test, version, and run end-to-end than notebooks.

## Planned

A future addition that would strengthen the portfolio:

* `01_quickstart.ipynb` — load the cholangiocyte AnnData, reproduce the SOX9 vs IAD score donor-level Spearman in 10 lines of code, and show the headline figure inline.
* `02_within_psc_trajectory.ipynb` — walk through the diffusion-pseudotime construction interactively, with intermediate plots.
* `03_pathway_dose_response.ipynb` — interactively reproduce the within-donor cell-level ρ(SOX9, pathway) test for a single pathway of the user's choice.

## How to contribute a notebook

If you want to add one:

1. Put it in this folder with a `NN_short_name.ipynb` prefix (numbered for run order).
2. Clear outputs before committing (`jupyter nbconvert --clear-output --inplace your_notebook.ipynb`).
3. Reference the scripts in `src/` for the actual heavy work rather than duplicating it.
