# `data/`

**This folder is intentionally empty in the repository.** All single-cell data used by the project is fetched on demand from CELLxGENE Census via the scripts in `src/`.

## Why not committed

The cholangiocyte AnnData alone is ~672 MB. Hepatocyte and full-liver intermediates push the working data above 2 GB. None of these files are appropriate for git distribution. The `.gitignore` at the repository root explicitly excludes everything under `data/` except this README.

## How to fetch

From the repository root, after creating the Python environment described in the main `README.md`:

```bash
python src/download_data.py                    # base liver fetch from Census
python src/refetch_cholangiocytes_with_prox1.py  # HVG with PROX1 + panel kept
```

This will create the working AnnData files under `data/liver/`. The fetch takes ~2–5 minutes on a residential connection and uses the CELLxGENE Census stable release.

## Expected layout after fetch

```
data/
├── liver/
│   ├── liver_raw.h5ad                          # initial Census fetch
│   ├── liver_processed.h5ad                    # post-QC + HVG
│   └── liver_cholangiocytes_with_prox1.h5ad    # the working file used downstream
└── pancreas/                                   # optional — only if running cross-tissue analyses
    └── pancreas_raw.h5ad
```

## Privacy note

All data used here is published, de-identified, and distributed through CELLxGENE Census public infrastructure. No raw sequencing reads, no patient identifiers, and no clinical metadata beyond what is already public are touched anywhere in this repository.
