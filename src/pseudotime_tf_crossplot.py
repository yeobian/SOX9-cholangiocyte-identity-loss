"""
Cross trajectory pseudotime against TF / pathway activity.

Reads the existing per-donor outputs:
  results/trajectory_per_donor.tsv   (donor, disease, mean_pseudotime, ...)
  results/tf_activity_per_donor.tsv  (donor, disease, SOX9, HNF1B, ...)

Joins on donor and computes Spearman ρ(mean_pseudotime, activity) per TF.
The hypothesis: as donors progress along pseudotime (healthy → exhausted),
cholangiocyte-identity TFs (BILIARY_id, SOX9, HNF1B) should drop monotonically
(negative ρ), and cells should show signs of stress-program collapse
(JUN/FOS also negative).

Outputs
-------
results/pseudotime_vs_tf_activity.tsv      (TF, n_donors, rho, p)
figures/comparison/pseudotime_vs_tf_bar.png
figures/comparison/pseudotime_vs_tf_scatter.png

Usage
-----
    python src/pseudotime_tf_crossplot.py
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir

RESULTS_DIR = Path("results")
FIGURES_DIR = Path("figures") / "comparison"

TRAJ_DONOR_TSV = RESULTS_DIR / "trajectory_per_donor.tsv"
TF_DONOR_TSV   = RESULTS_DIR / "tf_activity_per_donor.tsv"
OUT_TSV        = RESULTS_DIR / "pseudotime_vs_tf_activity.tsv"

# Optional: prefer the PSC-only trajectory if present
PSC_TRAJ_DONOR_TSV = RESULTS_DIR / "trajectory_psc_per_donor.tsv"


def pick_trajectory_donor() -> Path:
    if PSC_TRAJ_DONOR_TSV.exists():
        print(f"  Using PSC-only trajectory: {PSC_TRAJ_DONOR_TSV.name}")
        return PSC_TRAJ_DONOR_TSV
    print(f"  Using all-disease trajectory: {TRAJ_DONOR_TSV.name}")
    return TRAJ_DONOR_TSV


def load_and_join() -> pd.DataFrame:
    traj_path = pick_trajectory_donor()
    if not traj_path.exists():
        print(f"ERROR: {traj_path} not found. Run trajectory.py first.",
              file=sys.stderr); sys.exit(1)
    if not TF_DONOR_TSV.exists():
        print(f"ERROR: {TF_DONOR_TSV} not found. Run tf_activity.py first.",
              file=sys.stderr); sys.exit(1)

    traj = pd.read_csv(traj_path, sep="\t")
    tf   = pd.read_csv(TF_DONOR_TSV, sep="\t")

    # Drop cells with infinite pseudotime (off the connected component)
    traj = traj.replace([np.inf, -np.inf], np.nan)
    traj = traj.dropna(subset=["mean_pseudotime", "mean_iad_score"])
    print(f"  Trajectory donors (finite mean_pseudotime): {len(traj)}")
    print(f"  TF activity donors:                         {len(tf)}")

    # Inner-join on donor (drop the duplicated `disease` column from one side)
    tf_no_disease = tf.drop(columns=["disease"])
    merged = traj.merge(tf_no_disease, on="donor", how="inner")
    print(f"  Joined: {len(merged)} donors")
    return merged


def compute_correlations(df: pd.DataFrame) -> pd.DataFrame:
    tf_cols = [c for c in df.columns
               if c not in {"donor", "disease", "mean_pseudotime",
                            "mean_iad_score", "n_cells"}]
    rows = []
    for c in tf_cols:
        vals = df[c].astype(float).values
        ok   = ~np.isnan(vals) & ~np.isnan(df["mean_pseudotime"].values)
        if ok.sum() < 4:
            rows.append({"TF": c, "n_donors": int(ok.sum()),
                         "rho": float("nan"), "p": float("nan")})
            continue
        rho, p = spearmanr(df["mean_pseudotime"].values[ok], vals[ok])
        rows.append({"TF": c, "n_donors": int(ok.sum()),
                     "rho": float(rho), "p": float(p)})
    out = pd.DataFrame(rows).sort_values("rho")
    return out


def plot_bar(corr: pd.DataFrame, out: Path) -> None:
    if corr.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 0.45 * len(corr) + 1))
    colors = ["#c0392b" if r > 0 else "#2980b9" for r in corr["rho"]]
    sig    = corr["p"] < 0.05
    edge_colors = ["black" if s else "none" for s in sig]
    ax.barh(corr["TF"], corr["rho"], color=colors, edgecolor=edge_colors)
    ax.axvline(0, color="black", linewidth=0.8, linestyle="--")
    ax.set_xlabel("Spearman ρ (mean pseudotime  vs  TF/pathway activity)")
    ax.set_title("Donor-level: TF activity along trajectory\n"
                 "(black border = p<0.05; negative = drops along trajectory)")
    ensure_dir(out.parent)
    plt.tight_layout()
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def plot_scatter(merged: pd.DataFrame, corr: pd.DataFrame, out: Path,
                 top_n: int = 6) -> None:
    """Scatter mean_pseudotime vs each of the top-magnitude TFs, faceted."""
    if merged.empty or corr.empty:
        return
    top = corr.assign(absrho=lambda d: d["rho"].abs()).nlargest(top_n, "absrho")
    n = len(top)
    cols = 3
    rows = int(np.ceil(n / cols))
    fig, axes = plt.subplots(rows, cols, figsize=(4.5 * cols, 3.5 * rows))
    axes = np.atleast_1d(axes).ravel()
    diseases = sorted(merged["disease"].unique())
    cmap = dict(zip(diseases, plt.get_cmap("Set2").colors[:len(diseases)]))
    for ax, (_, row) in zip(axes, top.iterrows()):
        tf = row["TF"]
        for d in diseases:
            sub = merged[merged["disease"] == d]
            ax.scatter(sub["mean_pseudotime"], sub[tf],
                       s=np.clip(sub["n_cells"] / 5, 10, 200),
                       c=[cmap[d]], alpha=0.75, label=d)
        ax.set_title(f"{tf}   ρ={row['rho']:+.2f}  p={row['p']:.2g}",
                     fontsize=9)
        ax.set_xlabel("mean pseudotime")
        ax.set_ylabel("activity")
    for ax in axes[n:]:
        ax.axis("off")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(diseases),
               frameon=False, fontsize=8, bbox_to_anchor=(0.5, 1.02))
    ensure_dir(out.parent)
    plt.tight_layout(rect=(0, 0, 1, 0.97))
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close()
    print(f"  → {out}")


def main() -> None:
    ensure_dir(RESULTS_DIR)
    ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Pseudotime × TF activity (donor level)")
    print("=" * 60)

    merged = load_and_join()
    if merged.empty:
        print("ERROR: no donors after join.", file=sys.stderr); sys.exit(1)

    corr = compute_correlations(merged)
    print("\n--- Spearman ρ per TF/pathway (sorted ascending) ---")
    print(corr.to_string(index=False, float_format="%.3f"))
    corr.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"  → {OUT_TSV}")

    print("\n--- Plots ---")
    plot_bar(corr, FIGURES_DIR / "pseudotime_vs_tf_bar.png")
    plot_scatter(merged, corr, FIGURES_DIR / "pseudotime_vs_tf_scatter.png")

    print("\n--- Headline ---")
    drops = corr[(corr["rho"] < 0) & (corr["p"] < 0.05)]
    rises = corr[(corr["rho"] > 0) & (corr["p"] < 0.05)]
    if not drops.empty:
        print(f"  TFs that DROP along trajectory (p<0.05): "
              f"{', '.join(drops['TF'])}")
    if not rises.empty:
        print(f"  TFs that RISE along trajectory (p<0.05): "
              f"{', '.join(rises['TF'])}")
    if drops.empty and rises.empty:
        print("  No TF correlation reaches p<0.05 — donor n is small.")
    print("\nDone.")


if __name__ == "__main__":
    main()
