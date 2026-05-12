"""
Stream human liver + pancreas scRNA-seq data from CELLxGENE Census.

Pipeline stage: LOAD.

Liver cells supply the cholangiocyte population whose loss defines
Idiopathic Adulthood Ductopenia (IAD). Pancreas cells supply the
cross-tissue control (pancreatic ductal cells share the cholangiocyte
program but live in a different organ).

Outputs
-------
data/liver/liver_raw.h5ad
data/pancreas/pancreas_raw.h5ad

Usage
-----
    python src/download_data.py            # fetch both
    python src/download_data.py liver      # fetch one tissue
"""

import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir

# ── Config ────────────────────────────────────────────────────────────────────
DATA_DIR        = Path("data")
LIVER_OUT       = DATA_DIR / "liver"    / "liver_raw.h5ad"
PANCREAS_OUT    = DATA_DIR / "pancreas" / "pancreas_raw.h5ad"

CENSUS_VERSION  = "stable"

# Subsample caps keep peak RAM under control on a laptop.
# Liver has ~500k+ cells across all atlases — fetching everything OOM-kills a
# 16 GB Mac. We prefilter to cholangiocytes at query time (see CHOL_TYPES below)
# which collapses the fetch to a few thousand cells, then cap at 50k for safety.
MAX_CELLS_LIVER     = 50_000
MAX_CELLS_PANCREAS  = 100_000

# Cell-type prefilter for liver. Census cell_type values are CL ontology terms;
# these are the canonical cholangiocyte / bile-duct epithelial labels.
CHOL_TYPES = [
    "cholangiocyte",
    "epithelial cell of bile duct",
    "intrahepatic cholangiocyte",
    "extrahepatic cholangiocyte",
    "duct epithelial cell",
]

OBS_COLS = [
    "cell_type",
    "tissue",
    "tissue_general",
    "assay",
    "disease",
    "sex",
    "development_stage",
    "donor_id",
    "dataset_id",
]


# ── Fetch ─────────────────────────────────────────────────────────────────────

def _check_census():
    try:
        import cellxgene_census  # noqa: F401
    except ImportError:
        print("ERROR: cellxgene-census not installed.", file=sys.stderr)
        print("Run:  pip install cellxgene-census", file=sys.stderr)
        sys.exit(1)


def _build_filter(tissue: str, cell_types: list[str] | None) -> str:
    """Build a Census obs_value_filter string."""
    parts = [f"tissue_general == '{tissue}'"]
    if cell_types:
        quoted = ", ".join(f"'{t}'" for t in cell_types)
        parts.append(f"cell_type in [{quoted}]")
    return " and ".join(parts)


def fetch_tissue(
    tissue: str,
    max_cells: int,
    cell_types: list[str] | None = None,
) -> ad.AnnData:
    """
    Stream human cells of `tissue` from CELLxGENE Census.

    If `cell_types` is given, prefilter to those CL labels at query time —
    drastically reduces peak RAM when the tissue (e.g. liver) is large but
    only a sub-population is needed downstream.
    """
    _check_census()
    import cellxgene_census

    obs_filter = _build_filter(tissue, cell_types)
    print(f"\nConnecting to CELLxGENE Census (version={CENSUS_VERSION}) …")
    print(f"  Filter: {obs_filter}")
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        print(f"Fetching human {tissue} cells …")
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            obs_value_filter=obs_filter,
            obs_column_names=OBS_COLS,
        )

    print(f"  Fetched {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    if adata.n_obs == 0:
        print(
            f"WARNING: 0 cells returned. If you used a cell_type prefilter, "
            f"the labels may not exactly match Census ontology. "
            f"Try widening the cell_types list or removing the prefilter.",
            file=sys.stderr,
        )

    if adata.n_obs > max_cells:
        sc.pp.subsample(adata, n_obs=max_cells, random_state=42)
        print(f"  Subsampled → {adata.n_obs:,} cells (cap={max_cells:,})")

    # Promote gene symbols to var_names if available.
    if "feature_name" in adata.var.columns:
        adata.var_names = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()

    return adata


# ── Save ──────────────────────────────────────────────────────────────────────

def save_raw(adata: ad.AnnData, out_path: Path) -> None:
    ensure_dir(out_path.parent)
    print(f"  Saving → {out_path}")
    adata.write_h5ad(out_path)


# ── Main ──────────────────────────────────────────────────────────────────────

# (out_path, cell_cap, cell_type_prefilter_or_None)
PLAN = {
    "liver":    (LIVER_OUT,    MAX_CELLS_LIVER,    CHOL_TYPES),
    "pancreas": (PANCREAS_OUT, MAX_CELLS_PANCREAS, None),
}


def main() -> None:
    args = sys.argv[1:]
    targets = args if args else list(PLAN.keys())

    unknown = [t for t in targets if t not in PLAN]
    if unknown:
        print(f"ERROR: unknown tissue(s): {unknown}", file=sys.stderr)
        print(f"Valid options: {list(PLAN.keys())}", file=sys.stderr)
        sys.exit(1)

    print("=" * 60)
    print("CELLxGENE Census download")
    print("=" * 60)
    print(f"Targets: {targets}")

    for tissue in targets:
        out_path, cap, cell_types = PLAN[tissue]
        if out_path.exists():
            print(f"\n[skip] {tissue}: {out_path} already exists.")
            continue
        adata = fetch_tissue(tissue, cap, cell_types=cell_types)
        save_raw(adata, out_path)

    print("\nDone.\nNext: python src/preprocess_liver.py")


if __name__ == "__main__":
    main()
