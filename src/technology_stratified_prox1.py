"""
Technology-stratified PROX1 analysis: snRNA-seq vs scRNA-seq.

Context
-------
Dataset 1873a18a-66fd-4a4d-8277-a872c93f5b59 is single-nucleus RNA-seq
(per the paper title — https://doi.org/10.1016/j.jhep.2023.12.023, Journal
of Hepatology 2024). The Census `assay` field labels it as "10x 3' v3"
(the chemistry), and the cell-vs-nucleus distinction lives in
`suspension_type`, which our Census fetch did not pull. We therefore
stratify by `dataset_id`:

   snRNA-seq cohort:  1873a18a (6,853 cells, ~33% of pooled cholangiocytes)
   scRNA-seq cohort:  every other dataset (~14,185 cells)

snRNA-seq systematically over-represents nuclear-localised TF transcripts
(PROX1, GATA4) and under-represents cytoplasmic mRNAs (KRT19, KRT7). The
sensitivity analysis already showed every PSC-specific PROX1 finding was
driven by 1873a18a. The question this script asks: which PROX1 findings
are technology-driven (snRNA artefact) and which are biology that
replicates across both technologies?

Tests, run separately within each technology cohort:
  1. PROX1+ % among cholangiocytes (loose: PROX1 > 0; strict: top-quartile
     of the technology's own distribution).
  2. Per-cell ρ(PROX1, SOX9) — ALL, Normal, PSC.
  3. Donor ρ(mean PROX1, mean SOX9) — ALL, Normal, PSC.
  4. PROX1-high vs PROX1-low contrasts on SOX9 and BILIARY_id (Δ + MWU p).
  5. Per-donor PROX1+ abundance vs IAD score (ρ).
  6. SOX9 vs IAD at per-cell level (sanity: SOX9 collapse should hold in
     both technologies if it's biology, not snRNA artefact).

Outputs
-------
results/technology_stratified_prox1.tsv          (side-by-side numbers)
results/technology_stratified_prox1_verdict.tsv  (artefact vs biology calls)
figures/prox1/technology_stratified_prox1.png

Usage
-----
    python src/technology_stratified_prox1.py
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
SNRNA_DATASET = "1873a18a-66fd-4a4d-8277-a872c93f5b59"
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
    tech = np.where(
        adata.obs["dataset_id"].astype(str).values == SNRNA_DATASET,
        "snRNA", "scRNA",
    )
    return pd.DataFrame({
        "donor":   adata.obs["donor_id"].astype(str).values,
        "disease": adata.obs["disease"].astype(str).values,
        "dataset": adata.obs["dataset_id"].astype(str).values,
        "tech":    tech,
        "iad":     iad,
        "PROX1":   prox1,
        "SOX9":    sox9,
        "bil_id":  bil,
    })


def per_donor(df: pd.DataFrame, strict_thr: float) -> pd.DataFrame:
    rows = []
    for (donor, disease, tech), sub in df.groupby(
        ["donor", "disease", "tech"], observed=True
    ):
        if len(sub) < MIN_CELLS_PER_DONOR: continue
        rows.append({
            "donor": donor, "disease": disease, "tech": tech,
            "n_cells": len(sub),
            "frac_prox1_pos_loose":  float((sub["PROX1"] > 0).mean()),
            "frac_prox1_pos_strict": float((sub["PROX1"] > strict_thr).mean()),
            "mean_PROX1":            float(sub["PROX1"].mean()),
            "mean_SOX9":             float(sub["SOX9"].mean()),
            "iad_score":             float(sub["iad"].mean()),
        })
    return pd.DataFrame(rows)


def run_tests(df: pd.DataFrame, tech: str) -> dict:
    """Run the PROX1 test suite within a single technology."""
    out = {"tech": tech, "n_cells_total": int(len(df))}

    # PROX1-high threshold uses THIS technology's own quartile
    strict_thr = float(np.percentile(df["PROX1"], 75))
    out["strict_threshold"] = strict_thr

    # 1. PROX1+ %
    out["pct_prox1_loose"]  = float((df["PROX1"] > 0).mean() * 100)
    out["pct_prox1_strict"] = float((df["PROX1"] > strict_thr).mean() * 100)

    # 2. Per-cell ρ(PROX1, SOX9)
    for slc_name, slc in [
        ("ALL", df),
        ("Normal", df[df["disease"] == "normal"]),
        ("PSC",    df[df["disease"] == PSC]),
    ]:
        if len(slc) >= 30:
            rho, p = spearmanr(slc["PROX1"], slc["SOX9"])
            out[f"cell_rho_PS_{slc_name}"] = float(rho)
            out[f"cell_p_PS_{slc_name}"]   = float(p)
            out[f"cell_n_{slc_name}"]      = int(len(slc))
        else:
            out[f"cell_rho_PS_{slc_name}"] = float("nan")
            out[f"cell_p_PS_{slc_name}"]   = float("nan")
            out[f"cell_n_{slc_name}"]      = int(len(slc))

    # 3. Donor ρ(mean PROX1, mean SOX9) + 4. PROX1-high vs low SOX9/bil_id
    pd_donor = per_donor(df, strict_thr)
    out["n_donors_total"]  = int(len(pd_donor))
    out["n_donors_normal"] = int((pd_donor["disease"] == "normal").sum())
    out["n_donors_PSC"]    = int((pd_donor["disease"] == PSC).sum())

    for slc_name, mask in [("ALL", lambda d: d),
                            ("Normal", lambda d: d[d["disease"] == "normal"]),
                            ("PSC",    lambda d: d[d["disease"] == PSC])]:
        sub = mask(pd_donor)
        if len(sub) >= 4:
            rho, p = spearmanr(sub["mean_PROX1"], sub["mean_SOX9"])
            out[f"donor_rho_PS_{slc_name}"] = float(rho)
            out[f"donor_p_PS_{slc_name}"]   = float(p)
            out[f"donor_n_{slc_name}"]      = int(len(sub))
        else:
            out[f"donor_rho_PS_{slc_name}"] = float("nan")
            out[f"donor_p_PS_{slc_name}"]   = float("nan")
            out[f"donor_n_{slc_name}"]      = int(len(sub))

    # 4. PROX1-high vs PROX1-low: SOX9 and BILIARY_id Δ + p
    for slc_name, slc in [("ALL", df),
                          ("Normal", df[df["disease"] == "normal"]),
                          ("PSC",    df[df["disease"] == PSC])]:
        high = slc[slc["PROX1"] > strict_thr]
        low  = slc[~(slc["PROX1"] > strict_thr)]
        if len(high) >= 30 and len(low) >= 30:
            _, p1 = mannwhitneyu(high["SOX9"],  low["SOX9"],  alternative="two-sided")
            _, p2 = mannwhitneyu(high["bil_id"], low["bil_id"], alternative="two-sided")
            out[f"highlow_SOX9_diff_{slc_name}"]   = float(high["SOX9"].mean() - low["SOX9"].mean())
            out[f"highlow_SOX9_p_{slc_name}"]      = float(p1)
            out[f"highlow_bilid_diff_{slc_name}"]  = float(high["bil_id"].mean() - low["bil_id"].mean())
            out[f"highlow_bilid_p_{slc_name}"]     = float(p2)
        else:
            for k in (f"highlow_SOX9_diff_{slc_name}", f"highlow_SOX9_p_{slc_name}",
                      f"highlow_bilid_diff_{slc_name}", f"highlow_bilid_p_{slc_name}"):
                out[k] = float("nan")

    # 5. PROX1+ abundance vs IAD score, donor-level
    ok = pd_donor.dropna(subset=["frac_prox1_pos_loose", "iad_score"])
    if len(ok) >= 4:
        rho, p = spearmanr(ok["frac_prox1_pos_loose"], ok["iad_score"])
        out["abund_vs_IAD_rho"] = float(rho)
        out["abund_vs_IAD_p"]   = float(p)
        out["abund_vs_IAD_n"]   = int(len(ok))
    else:
        out["abund_vs_IAD_rho"] = out["abund_vs_IAD_p"] = float("nan")
        out["abund_vs_IAD_n"] = int(len(ok))

    # 6. SOX9 vs IAD per cell (sanity — SOX9 collapse should hold in both techs)
    rho, p = spearmanr(df["SOX9"], df["iad"])
    out["sox9_vs_iad_cell_rho"] = float(rho)
    out["sox9_vs_iad_cell_p"]   = float(p)

    return out, pd_donor


def verdict(claim: str, sn: dict, sc: dict, key: str,
            kind: str) -> tuple:
    """Compare values between sn and sc for one metric.

    kind = 'rho' -> sign concordance, both |ρ| > 0.05 = concordant biology
           'diff' -> sign concordance
           'pct' -> compare PROX1+ % magnitudes
    """
    v_sn = sn.get(key, float("nan"))
    v_sc = sc.get(key, float("nan"))
    if any(np.isnan([v_sn, v_sc])):
        return claim, v_sn, v_sc, "indeterminate"
    if kind in ("rho", "diff"):
        same_sign = (np.sign(v_sn) == np.sign(v_sc)) and v_sn != 0 and v_sc != 0
        if same_sign and abs(v_sn) >= 0.05 and abs(v_sc) >= 0.05:
            verd = "CONCORDANT (biology)"
        elif same_sign:
            verd = "same sign, small"
        else:
            verd = "DISCORDANT (likely tech)"
        return claim, v_sn, v_sc, verd
    if kind == "pct":
        diff = v_sn - v_sc
        verd = f"sn − sc = {diff:+.1f} pp"
        return claim, v_sn, v_sc, verd
    return claim, v_sn, v_sc, ""


def plot_compare(sn: dict, sc_d: dict, out: Path) -> None:
    metrics = [
        ("pct_prox1_loose",         "PROX1+ % (loose)"),
        ("abund_vs_IAD_rho",        "abund vs IAD ρ"),
        ("sox9_vs_iad_cell_rho",    "SOX9 vs IAD ρ (cell)"),
        ("cell_rho_PS_ALL",         "cell ρ(P,S) ALL"),
        ("cell_rho_PS_Normal",      "cell ρ(P,S) Normal"),
        ("cell_rho_PS_PSC",         "cell ρ(P,S) PSC"),
        ("donor_rho_PS_ALL",        "donor ρ(P,S) ALL"),
        ("donor_rho_PS_Normal",     "donor ρ(P,S) Normal"),
        ("donor_rho_PS_PSC",        "donor ρ(P,S) PSC"),
        ("highlow_SOX9_diff_ALL",   "PROX1-high SOX9 Δ ALL"),
        ("highlow_bilid_diff_ALL",  "PROX1-high BIL_id Δ ALL"),
        ("highlow_SOX9_diff_PSC",   "PROX1-high SOX9 Δ PSC"),
        ("highlow_bilid_diff_PSC",  "PROX1-high BIL_id Δ PSC"),
    ]
    fig, ax = plt.subplots(figsize=(9, 0.45 * len(metrics) + 2))
    y = np.arange(len(metrics))
    sn_v = [sn.get(k, np.nan) for k, _ in metrics]
    sc_v = [sc_d.get(k, np.nan) for k, _ in metrics]
    h = 0.38
    ax.barh(y - h/2, sn_v, height=h, color="#8e44ad", label="snRNA-seq (1873a18a)")
    ax.barh(y + h/2, sc_v, height=h, color="#16a085", label="scRNA-seq (all other datasets)")
    ax.set_yticks(y); ax.set_yticklabels([m[1] for m in metrics], fontsize=8)
    ax.axvline(0, color="black", linewidth=0.4)
    ax.set_title("PROX1 tests stratified by technology")
    ax.legend(loc="lower right", fontsize=9, frameon=False)
    plt.tight_layout()
    ensure_dir(out.parent)
    plt.savefig(out, bbox_inches="tight", dpi=150)
    plt.close(); print(f"  → {out}")


def main() -> None:
    sc.settings.verbosity = 1
    ensure_dir(RESULTS_DIR); ensure_dir(FIGURES_DIR)

    print("=" * 60)
    print("Technology-stratified PROX1 analysis (sn vs scRNA-seq)")
    print("=" * 60)

    if not CHOL_PATH.exists():
        print(f"ERROR: {CHOL_PATH} missing", file=sys.stderr); sys.exit(1)
    print(f"\nLoading {CHOL_PATH}")
    adata = sc.read_h5ad(CHOL_PATH)
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    df = build_long(adata)
    print(f"  snRNA cells: {(df['tech']=='snRNA').sum():,}")
    print(f"  scRNA cells: {(df['tech']=='scRNA').sum():,}")
    print()
    print("  Disease × technology cross-tab:")
    print(pd.crosstab(df["tech"], df["disease"]).to_string())

    print("\n[1/3] Running tests within snRNA-seq …")
    sn_results, sn_donors = run_tests(df[df["tech"] == "snRNA"], "snRNA")
    print(f"  cells={sn_results['n_cells_total']:,}  "
          f"donors total={sn_results['n_donors_total']}  "
          f"normal={sn_results['n_donors_normal']}  "
          f"PSC={sn_results['n_donors_PSC']}")
    print(f"  PROX1+ % (loose) = {sn_results['pct_prox1_loose']:.1f}%")
    print(f"  PROX1+ threshold (top quartile of snRNA) = "
          f"{sn_results['strict_threshold']:.3f}")

    print("\n[2/3] Running tests within scRNA-seq …")
    sc_results, sc_donors = run_tests(df[df["tech"] == "scRNA"], "scRNA")
    print(f"  cells={sc_results['n_cells_total']:,}  "
          f"donors total={sc_results['n_donors_total']}  "
          f"normal={sc_results['n_donors_normal']}  "
          f"PSC={sc_results['n_donors_PSC']}")
    print(f"  PROX1+ % (loose) = {sc_results['pct_prox1_loose']:.1f}%")
    print(f"  PROX1+ threshold (top quartile of scRNA) = "
          f"{sc_results['strict_threshold']:.3f}")

    # Side-by-side numbers
    print("\n[3/3] Compare …")
    common_keys = sorted(set(sn_results.keys()) & set(sc_results.keys()) -
                          {"tech"})
    rows = []
    for k in common_keys:
        rows.append({"metric": k, "snRNA": sn_results[k], "scRNA": sc_results[k]})
    side_by_side = pd.DataFrame(rows)
    side_by_side.to_csv(RESULTS_DIR / "technology_stratified_prox1.tsv",
                        sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'technology_stratified_prox1.tsv'}")

    # Verdict per claim
    verdicts = [
        verdict("PROX1+ % (loose)",                 sn_results, sc_results,
                "pct_prox1_loose", "pct"),
        verdict("PROX1+ abundance vs donor IAD ρ",  sn_results, sc_results,
                "abund_vs_IAD_rho", "rho"),
        verdict("SOX9 vs IAD (cell ρ)",             sn_results, sc_results,
                "sox9_vs_iad_cell_rho", "rho"),
        verdict("cell ρ(PROX1, SOX9) ALL",          sn_results, sc_results,
                "cell_rho_PS_ALL", "rho"),
        verdict("cell ρ(PROX1, SOX9) Normal",       sn_results, sc_results,
                "cell_rho_PS_Normal", "rho"),
        verdict("cell ρ(PROX1, SOX9) PSC",          sn_results, sc_results,
                "cell_rho_PS_PSC", "rho"),
        verdict("donor ρ(PROX1, SOX9) ALL",         sn_results, sc_results,
                "donor_rho_PS_ALL", "rho"),
        verdict("donor ρ(PROX1, SOX9) Normal",      sn_results, sc_results,
                "donor_rho_PS_Normal", "rho"),
        verdict("donor ρ(PROX1, SOX9) PSC",         sn_results, sc_results,
                "donor_rho_PS_PSC", "rho"),
        verdict("PROX1-high SOX9 Δ (ALL)",          sn_results, sc_results,
                "highlow_SOX9_diff_ALL", "diff"),
        verdict("PROX1-high BIL_id Δ (ALL)",        sn_results, sc_results,
                "highlow_bilid_diff_ALL", "diff"),
        verdict("PROX1-high SOX9 Δ (Normal)",       sn_results, sc_results,
                "highlow_SOX9_diff_Normal", "diff"),
        verdict("PROX1-high BIL_id Δ (Normal)",     sn_results, sc_results,
                "highlow_bilid_diff_Normal", "diff"),
        verdict("PROX1-high SOX9 Δ (PSC)",          sn_results, sc_results,
                "highlow_SOX9_diff_PSC", "diff"),
        verdict("PROX1-high BIL_id Δ (PSC)",        sn_results, sc_results,
                "highlow_bilid_diff_PSC", "diff"),
    ]
    v_df = pd.DataFrame(verdicts, columns=["claim", "snRNA", "scRNA", "verdict"])
    v_df.to_csv(RESULTS_DIR / "technology_stratified_prox1_verdict.tsv",
                sep="\t", index=False)
    print(f"  → {RESULTS_DIR / 'technology_stratified_prox1_verdict.tsv'}")

    print("\n--- Side-by-side verdicts ---")
    print(v_df.to_string(index=False, float_format="%.3f"))

    plot_compare(sn_results, sc_results,
                 FIGURES_DIR / "technology_stratified_prox1.png")

    # Headline
    print("\n--- Headline ---")
    concordant = v_df[v_df["verdict"].str.startswith("CONCORDANT")]
    discordant = v_df[v_df["verdict"].str.startswith("DISCORDANT")]
    print(f"  CONCORDANT (sn and sc agree, biology): {len(concordant)} claims")
    if not concordant.empty:
        for _, r in concordant.iterrows():
            print(f"    {r['claim']:<35}  sn={r['snRNA']:+.3f}  sc={r['scRNA']:+.3f}")
    print(f"\n  DISCORDANT (sn and sc disagree, technology-confounded): "
          f"{len(discordant)} claims")
    if not discordant.empty:
        for _, r in discordant.iterrows():
            print(f"    {r['claim']:<35}  sn={r['snRNA']:+.3f}  sc={r['scRNA']:+.3f}")


if __name__ == "__main__":
    main()
