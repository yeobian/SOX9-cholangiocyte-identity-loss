"""
Sensitivity analysis: drop dataset 1873a18a, re-run the key PROX1 tests.

Every PROX1 claim so far (within-PSC abundance, donor-level PROX1↔SOX9
correlation, the DE signature) is essentially within-one-dataset — 7 of 8
informative PSC donors come from `1873a18a-66fd-4a4d-8277-a872c93f5b59`.
Before claiming any PROX1 result, we need to know: what survives if that
dataset is removed?

We compare WITH 1873a18a vs WITHOUT 1873a18a, side by side, on:
  1. PROX1+ abundance: PSC vs normal (Mann–Whitney)
  2. PROX1+ abundance vs IAD score across donors (Spearman)
  3. Per-cell ρ(PROX1, SOX9) — ALL, Normal, PSC
  4. Donor ρ(mean PROX1, mean SOX9) — ALL, Normal, PSC
  5. PROX1-high vs PROX1-low contrasts on SOX9 and BILIARY_id (Mann–Whitney)

Outputs
-------
results/sensitivity_drop_1873a18a.tsv         (side-by-side numbers)
results/sensitivity_summary.tsv               (which results survive)
figures/prox1/sensitivity_drop_1873a18a.png   (paired bar chart of effect sizes)

Usage
-----
    python src/sensitivity_drop_1873a18a.py
"""

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from utils import ensure_dir, get_gene_expression
from iad_score import (
    LIVER_H5AD as DEFAULT_CHOL,
    REFERENCE_PERCENTILE,
    _panel_expression_per_cell, _panel_genes_present, load_panel,
)

CHOL_PATH    = DEFAULT_CHOL
RESULTS_DIR  = Path("results")
FIGURES_DIR  = Path("figures") / "prox1"

PSC = "primary sclerosing cholangitis"
TARGET_DS = "1873a18a-66fd-4a4d-8277-a872c93f5b59"
MIN_CELLS_PER_DONOR = 20

BILIARY_ID_GENES = ["KRT19", "KRT7", "EPCAM", "SOX9", "HNF1B", "AQP1"]


def per_cell_iad(adata) -> np.ndarray:
    panel = load_panel()
    genes = _panel_genes_present(adata, panel)
    weights = pd.Series(panel.set_index("gene")["specificity"])
    panel_expr = _panel_expression_per_cell(adata, genes, weights=weights)
    ref = float(np.percentile(panel_expr, REFERENCE_PERCENTILE))
    return 1.0 - np.clip(panel_expr / max(ref, 1e-6), 0.0, 1.0)


def build_long(adata):
    iad   = per_cell_iad(adata)
    prox1 = get_gene_expression(adata, "PROX1")
    sox9  = get_gene_expression(adata, "SOX9")
    sc.tl.score_genes(adata,
        gene_list=[g for g in BILIARY_ID_GENES if g in adata.var_names],
        score_name="bil_id", random_state=42)
    bil = adata.obs["bil_id"].values.copy()
    del adata.obs["bil_id"]
    df = pd.DataFrame({
        "donor":   adata.obs["donor_id"].astype(str).values,
        "disease": adata.obs["disease"].astype(str).values,
        "dataset": adata.obs["dataset_id"].astype(str).values,
        "iad":     iad,
        "PROX1":   prox1,
        "SOX9":    sox9,
        "bil_id":  bil,
    })
    return df


def per_donor(df: pd.DataFrame, strict_thr: float) -> pd.DataFrame:
    rows = []
    for (donor, disease), sub in df.groupby(["donor", "disease"], observed=True):
        n = len(sub)
        if n < MIN_CELLS_PER_DONOR: continue
        rows.append({
            "donor": donor, "disease": disease,
            "dataset": sub["dataset"].mode().iloc[0],
            "n_cells": n,
            "frac_prox1_pos_loose":  float((sub["PROX1"] > 0).mean()),
            "frac_prox1_pos_strict": float((sub["PROX1"] > strict_thr).mean()),
            "mean_PROX1":            float(sub["PROX1"].mean()),
            "mean_SOX9":             float(sub["SOX9"].mean()),
            "iad_score":             float(sub["iad"].mean()),
        })
    return pd.DataFrame(rows)


def run_all(df: pd.DataFrame, label: str, strict_thr: float) -> dict:
    """Run every key test on `df` and return the numbers."""
    out = {"slice": label}

    pd_donor = per_donor(df, strict_thr)
    out["n_donors_total"]  = int(len(pd_donor))
    out["n_donors_normal"] = int((pd_donor["disease"] == "normal").sum())
    out["n_donors_PSC"]    = int((pd_donor["disease"] == PSC).sum())

    # 1. Abundance PSC vs normal (loose threshold)
    a = pd_donor[pd_donor["disease"] == PSC]["frac_prox1_pos_loose"]
    b = pd_donor[pd_donor["disease"] == "normal"]["frac_prox1_pos_loose"]
    if len(a) >= 3 and len(b) >= 3:
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        out["abund_PSC_vs_normal_p"] = float(p)
        out["abund_PSC_median"] = float(a.median())
        out["abund_normal_median"] = float(b.median())
    else:
        out["abund_PSC_vs_normal_p"] = float("nan")
        out["abund_PSC_median"] = float("nan")
        out["abund_normal_median"] = float("nan")

    # 2. Abundance vs IAD across all donors
    ok = pd_donor.dropna(subset=["frac_prox1_pos_loose", "iad_score"])
    if len(ok) >= 4:
        rho, p = spearmanr(ok["frac_prox1_pos_loose"], ok["iad_score"])
        out["abund_vs_IAD_rho"] = float(rho); out["abund_vs_IAD_p"] = float(p)
        out["abund_vs_IAD_n"] = int(len(ok))
    else:
        out["abund_vs_IAD_rho"] = out["abund_vs_IAD_p"] = float("nan")
        out["abund_vs_IAD_n"] = int(len(ok))

    # 3. per-cell ρ(PROX1, SOX9) within each slice
    for slc, sub in [("ALL", df),
                     ("Normal", df[df["disease"] == "normal"]),
                     ("PSC",    df[df["disease"] == PSC])]:
        if len(sub) >= 30:
            rho, p = spearmanr(sub["PROX1"], sub["SOX9"])
            out[f"cell_rho_PROX1_SOX9_{slc}"] = float(rho)
            out[f"cell_p_PROX1_SOX9_{slc}"]   = float(p)
            out[f"cell_n_{slc}"] = int(len(sub))
        else:
            out[f"cell_rho_PROX1_SOX9_{slc}"] = float("nan")
            out[f"cell_p_PROX1_SOX9_{slc}"]   = float("nan")
            out[f"cell_n_{slc}"] = int(len(sub))

    # 4. donor ρ(mean PROX1, mean SOX9) within each slice
    for slc in ["ALL", "Normal", "PSC"]:
        if slc == "ALL":
            sub = pd_donor
        elif slc == "Normal":
            sub = pd_donor[pd_donor["disease"] == "normal"]
        else:
            sub = pd_donor[pd_donor["disease"] == PSC]
        if len(sub) >= 4:
            rho, p = spearmanr(sub["mean_PROX1"], sub["mean_SOX9"])
            out[f"donor_rho_PROX1_SOX9_{slc}"] = float(rho)
            out[f"donor_p_PROX1_SOX9_{slc}"]   = float(p)
            out[f"donor_n_{slc}"] = int(len(sub))
        else:
            out[f"donor_rho_PROX1_SOX9_{slc}"] = float("nan")
            out[f"donor_p_PROX1_SOX9_{slc}"]   = float("nan")
            out[f"donor_n_{slc}"] = int(len(sub))

    # 5. PROX1-high vs PROX1-low on SOX9 + BILIARY_id (within ALL, Normal, PSC)
    for slc, sub in [("ALL", df),
                     ("Normal", df[df["disease"] == "normal"]),
                     ("PSC",    df[df["disease"] == PSC])]:
        high = sub[sub["PROX1"] > strict_thr]
        low  = sub[~(sub["PROX1"] > strict_thr)]
        if len(high) >= 30 and len(low) >= 30:
            u, p1 = mannwhitneyu(high["SOX9"], low["SOX9"],
                                 alternative="two-sided")
            u, p2 = mannwhitneyu(high["bil_id"], low["bil_id"],
                                 alternative="two-sided")
            out[f"highlow_SOX9_diff_{slc}"]   = float(high["SOX9"].mean() - low["SOX9"].mean())
            out[f"highlow_SOX9_p_{slc}"]      = float(p1)
            out[f"highlow_bilid_diff_{slc}"]  = float(high["bil_id"].mean() - low["bil_id"].mean())
            out[f"highlow_bilid_p_{slc}"]     = float(p2)
        else:
            for k in (f"highlow_SOX9_diff_{slc}", f"highlow_SOX9_p_{slc}",
                      f"highlow_bilid_diff_{slc}", f"highlow_bilid_p_{slc}"):
                out[k] = float("nan")
    return out


def summary_compare(with_ds: dict, without_ds: dict) -> pd.DataFrame:
    """Side-by-side comparison with verdicts."""
    rows = []

    def _row(key, label, expect_dir=None, sig_at=0.05):
        v_w  = with_ds.get(key, float("nan"))
        v_wo = without_ds.get(key, float("nan"))
        rows.append({
            "metric":         label,
            "with_1873":      v_w,
            "without_1873":   v_wo,
        })

    # Abundance PSC vs normal
    _row("abund_PSC_median",          "abundance PSC median")
    _row("abund_normal_median",       "abundance normal median")
    _row("abund_PSC_vs_normal_p",     "abundance PSC vs normal — Mann–Whitney p")

    # Abundance vs IAD
    _row("abund_vs_IAD_rho",          "abundance vs IAD score — donor ρ")
    _row("abund_vs_IAD_p",            "abundance vs IAD score — p")

    # Per-cell ρ
    for slc in ("ALL", "Normal", "PSC"):
        _row(f"cell_rho_PROX1_SOX9_{slc}", f"per-cell ρ(PROX1,SOX9) {slc}")
        _row(f"cell_p_PROX1_SOX9_{slc}",   f"per-cell p {slc}")

    # Donor ρ
    for slc in ("ALL", "Normal", "PSC"):
        _row(f"donor_rho_PROX1_SOX9_{slc}", f"donor ρ(mean PROX1, mean SOX9) {slc}")
        _row(f"donor_p_PROX1_SOX9_{slc}",   f"donor p {slc}")

    # high-low contrasts
    for slc in ("ALL", "Normal", "PSC"):
        _row(f"highlow_SOX9_diff_{slc}",  f"PROX1-high − low SOX9 mean diff ({slc})")
        _row(f"highlow_SOX9_p_{slc}",     f"PROX1-high vs low SOX9 p ({slc})")
        _row(f"highlow_bilid_diff_{slc}", f"PROX1-high − low BILIARY_id diff ({slc})")
        _row(f"highlow_bilid_p_{slc}",    f"PROX1-high vs low BILIARY_id p ({slc})")

    return pd.DataFrame(rows)


def survives(with_v: float, without_v: float, kind: str) -> str:
    """Verdict on whether a finding survives without 1873a18a.

    kind = 'p' (smaller is better, threshold 0.05)
    kind = 'rho' (sign should be preserved, magnitude > 0.1)
    kind = 'diff' (sign should be preserved)
    """
    if any(np.isnan([with_v, without_v])):
        return "indeterminate"
    if kind == "p":
        if with_v < 0.05 and without_v < 0.05: return "SURVIVES"
        if with_v < 0.05 and without_v >= 0.05: return "WEAKENED"
        if with_v >= 0.05 and without_v >= 0.05: return "n.s. either way"
        return "GAINED"
    if kind == "rho":
        same_sign = np.sign(with_v) == np.sign(without_v) and with_v != 0
        if same_sign and abs(without_v) >= 0.1: return "SURVIVES (same sign)"
        if same_sign:                            return "SAME SIGN, smaller"
        return "FLIPS SIGN"
    if kind == "diff":
        same_sign = np.sign(with_v) == np.sign(without_v) and with_v != 0
        return "SAME SIGN" if same_sign else "FLIPS SIGN"
    return ""


def make_verdicts(with_ds: dict, without_ds: dict) -> pd.DataFrame:
    rows = []
    rows.append({"claim": "PROX1+ abundance PSC vs normal",
                 "with_value":    with_ds["abund_PSC_vs_normal_p"],
                 "without_value": without_ds["abund_PSC_vs_normal_p"],
                 "verdict": survives(with_ds["abund_PSC_vs_normal_p"],
                                     without_ds["abund_PSC_vs_normal_p"], "p")})
    rows.append({"claim": "PROX1+ abundance vs IAD score across donors",
                 "with_value":    with_ds["abund_vs_IAD_rho"],
                 "without_value": without_ds["abund_vs_IAD_rho"],
                 "verdict": survives(with_ds["abund_vs_IAD_rho"],
                                     without_ds["abund_vs_IAD_rho"], "rho")})
    for slc in ("ALL", "Normal", "PSC"):
        rows.append({
            "claim": f"per-cell ρ(PROX1,SOX9) {slc}",
            "with_value":    with_ds[f"cell_rho_PROX1_SOX9_{slc}"],
            "without_value": without_ds[f"cell_rho_PROX1_SOX9_{slc}"],
            "verdict": survives(with_ds[f"cell_rho_PROX1_SOX9_{slc}"],
                                without_ds[f"cell_rho_PROX1_SOX9_{slc}"], "rho"),
        })
    for slc in ("ALL", "Normal", "PSC"):
        rows.append({
            "claim": f"donor ρ(mean PROX1, mean SOX9) {slc}",
            "with_value":    with_ds[f"donor_rho_PROX1_SOX9_{slc}"],
            "without_value": without_ds[f"donor_rho_PROX1_SOX9_{slc}"],
            "verdict": survives(with_ds[f"donor_rho_PROX1_SOX9_{slc}"],
                                without_ds[f"donor_rho_PROX1_SOX9_{slc}"], "rho"),
        })
    for slc in ("ALL", "Normal", "PSC"):
        rows.append({
            "claim": f"PROX1-high vs low SOX9 sign ({slc})",
            "with_value":    with_ds[f"highlow_SOX9_diff_{slc}"],
            "without_value": without_ds[f"highlow_SOX9_diff_{slc}"],
            "verdict": survives(with_ds[f"highlow_SOX9_diff_{slc}"],
                                without_ds[f"highlow_SOX9_diff_{slc}"], "diff"),
        })
        rows.append({
            "claim": f"PROX1-high vs low BILIARY_id sign ({slc})",
            "with_value":    with_ds[f"highlow_bilid_diff_{slc}"],
            "without_value": without_ds[f"highlow_bilid_diff_{slc}"],
            "verdict": survives(with_ds[f"highlow_bilid_diff_{slc}"],
                                without_ds[f"highlow_bilid_diff_{slc}"], "diff"),
        })
    return pd.DataFrame(rows)


def plot_paired(with_ds: dict, without_ds: dict, out: Path) -> None:
    metrics = [
        ("abund_vs_IAD_rho", "abundance vs IAD ρ"),
        ("cell_rho_PROX1_SOX9_ALL",    "cell ρ(P,S) ALL"),
        ("cell_rho_PROX1_SOX9_Normal", "cell ρ(P,S) Normal"),
        ("cell_rho_PROX1_SOX9_PSC",    "cell ρ(P,S) PSC"),
        ("donor_rho_PROX1_SOX9_ALL",    "donor ρ(P,S) ALL"),
        ("donor_rho_PROX1_SOX9_Normal", "donor ρ(P,S) Normal"),
        ("donor_rho_PROX1_SOX9_PSC",    "donor ρ(P,S) PSC"),
        ("highlow_SOX9_diff_ALL",   "high−low SOX9 (ALL)"),
        ("highlow_SOX9_diff_PSC",   "high−low SOX9 (PSC)"),
        ("highlow_bilid_diff_ALL", "high−low BIL (ALL)"),
        ("highlow_bilid_diff_PSC", "high−low BIL (PSC)"),
    ]
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(metrics) + 2))
    y = np.arange(len(metrics))
    with_v = [with_ds.get(k, np.nan) for k, _ in metrics]
    wo_v   = [without_ds.get(k, np.nan) for k, _ in metrics]
    h = 0.38
    ax.barh(y - h/2, with_v,  height=h, color="#c0392b", label="with 1873a18a")
    ax.barh(y + h/2, wo_v,    height=h, color="#2980b9", label="without 1873a18a")
    ax.set_yticks(y)
    ax.set_yticklabels([m[1] for m in metrics], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_xlabel("effect size")
    ax.set_title("Sensitivity: dropping dataset 1873a18a")
    ax.legend(fontsize=9, frameon=False)
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Sensitivity: drop dataset 1873a18a, re-run PROX1 tests")
    print("=" * 60)
    print(f"\nLoading {CHOL_PATH} …")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    df_all = build_long(adata)
    strict_thr = float(np.percentile(df_all["PROX1"], 75))
    print(f"  PROX1-high threshold (top quartile, ALL data) = {strict_thr:.3f}")

    n_target = int((df_all["dataset"] == TARGET_DS).sum())
    print(f"  Cells in target dataset {TARGET_DS[:14]}…: {n_target:,}")

    print("\n[1/3] Running tests WITH 1873a18a (baseline) …")
    with_ds = run_all(df_all, "with_1873a18a", strict_thr)
    print(f"  donors: total={with_ds['n_donors_total']}  "
          f"normal={with_ds['n_donors_normal']}  PSC={with_ds['n_donors_PSC']}")

    print("\n[2/3] Running tests WITHOUT 1873a18a …")
    df_drop = df_all[df_all["dataset"] != TARGET_DS].copy()
    print(f"  remaining cells: {len(df_drop):,}/{len(df_all):,}")
    without_ds = run_all(df_drop, "without_1873a18a", strict_thr)
    print(f"  donors: total={without_ds['n_donors_total']}  "
          f"normal={without_ds['n_donors_normal']}  PSC={without_ds['n_donors_PSC']}")

    # Build comparison + verdict tables
    table = summary_compare(with_ds, without_ds)
    table.to_csv(RESULTS_DIR / "sensitivity_drop_1873a18a.tsv",
                 sep="\t", index=False)
    print(f"\n  → {RESULTS_DIR / 'sensitivity_drop_1873a18a.tsv'}")

    verdict = make_verdicts(with_ds, without_ds)
    verdict.to_csv(RESULTS_DIR / "sensitivity_summary.tsv",
                   sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'sensitivity_summary.tsv'}")

    print("\n[3/3] Verdicts:")
    print(verdict.to_string(index=False, float_format="%.3f"))

    plot_paired(with_ds, without_ds,
                FIGURES_DIR / "sensitivity_drop_1873a18a.png")

    # Headline
    print("\n--- Headline (which claims survive) ---")
    surv = verdict[verdict["verdict"].str.startswith(("SURVIVES", "SAME SIGN"))]
    drop = verdict[verdict["verdict"].isin(["FLIPS SIGN", "WEAKENED"])]
    nochg = verdict[verdict["verdict"] == "n.s. either way"]
    print(f"  SURVIVES / SAME SIGN: {len(surv)} claims")
    print(f"  FLIPS SIGN or WEAKENED: {len(drop)} claims")
    print(f"  n.s. in both: {len(nochg)} claims")
    if not drop.empty:
        print("\n  Claims that DON'T survive:")
        for _, r in drop.iterrows():
            print(f"    {r['claim']}  "
                  f"with={r['with_value']:+.3f}  without={r['without_value']:+.3f}  "
                  f"→ {r['verdict']}")


if __name__ == "__main__":
    main()
