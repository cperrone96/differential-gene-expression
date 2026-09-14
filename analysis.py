"""Reproducible educational RNA-seq benchmark.

This is intentionally not a replacement for a count-aware negative-binomial
workflow such as DESeq2 or edgeR.  It benchmarks a transparent Welch test on
log2 normalized counts against a composition-aware median-of-ratios
normalization, using repeated synthetic experiments with known truth.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).parent))
from generate_data import generate  # noqa: E402


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    """DESeq-style median-of-ratios factors, robust to library composition."""
    positive = counts.where(counts > 0)
    geo = np.exp(np.log(positive).mean(axis=1, skipna=True))
    valid = geo.notna() & np.isfinite(geo) & (geo > 0)
    ratios = counts.loc[valid].div(geo[valid], axis=0)
    factors = ratios.replace([np.inf, -np.inf], np.nan).median(axis=0)
    return factors / factors.median()


def analyze_tables(counts, info, truth, min_cpm=5.0, min_samples=6):
    """Run one experiment and return per-gene results plus summary metrics."""
    lib = counts.sum(axis=0)
    cpm = counts.div(lib, axis=1) * 1e6
    keep = (cpm >= min_cpm).sum(axis=1) >= min_samples
    counts = counts.loc[keep]
    factors = median_ratio_size_factors(counts)
    normalized = np.log2(counts.div(factors, axis=1) + 1)
    control = info.loc[info.group == "control", "sample"].tolist()
    disease = info.loc[info.group == "disease", "sample"].tolist()
    x, y = normalized[control], normalized[disease]
    test = stats.ttest_ind(y, x, axis=1, equal_var=False, nan_policy="omit")
    diff = y.mean(axis=1) - x.mean(axis=1)
    se = np.sqrt(y.var(axis=1, ddof=1) / len(disease) + x.var(axis=1, ddof=1) / len(control))
    df = (y.var(axis=1, ddof=1) / len(disease) + x.var(axis=1, ddof=1) / len(control)) ** 2 / (
        (y.var(axis=1, ddof=1) / len(disease)) ** 2 / (len(disease) - 1)
        + (x.var(axis=1, ddof=1) / len(control)) ** 2 / (len(control) - 1)
    )
    crit = stats.t.ppf(0.975, df)
    result = pd.DataFrame({"log2fc": diff, "pvalue": test.pvalue, "ci_low": diff - crit * se, "ci_high": diff + crit * se})
    result["fdr"] = multipletests(result.pvalue.fillna(1), method="fdr_bh")[1]
    result["significant"] = (result.fdr < 0.05) & (result.log2fc.abs() >= 1)
    result = result.join(truth.set_index("gene"))
    called, real = result.significant, result.is_de.astype(bool)
    tp = int((called & real).sum()); fp = int((called & ~real).sum())
    fn = int((~called & real).sum())
    metrics = {
        "genes_input": int(len(keep)), "genes_retained": int(keep.sum()),
        "true_de": int(real.sum()), "calls": int(called.sum()), "true_positives": tp, "false_positives": fp,
        "empirical_fdr": float(fp / max(int(called.sum()), 1)), "power": float(tp / max(tp + fn, 1)),
        "bias_mean_log2fc": float((result.loc[real, "log2fc"] - result.loc[real, "true_log2fc"]).mean()),
        "coverage_95pct": float(((result.loc[real, "ci_low"] <= result.loc[real, "true_log2fc"]) & (result.loc[real, "true_log2fc"] <= result.loc[real, "ci_high"])).mean()),
        "normalization_factor_min": float(factors.min()), "normalization_factor_max": float(factors.max()),
    }
    return result, metrics, normalized


def run(seed=2024, repeats=20, output_dir="outputs"):
    out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
    counts, info, truth = generate(seed)
    result, one, normalized = analyze_tables(counts, info, truth)
    result.to_csv(out / "differential_expression_results.csv", index_label="gene")
    de = result[result.is_de == 1].copy()
    de["effect_bin"] = pd.cut(de.true_log2fc.abs(), bins=[0, 1, 1.5, 2, 2.5, np.inf], include_lowest=True)
    sensitivity = de.groupby("effect_bin", observed=False).significant.mean().reset_index(name="power")
    sensitivity.to_csv(out / "sensitivity_by_true_effect.csv", index=False)
    fig, ax = plt.subplots(figsize=(7, 4)); ax.plot(sensitivity.effect_bin.astype(str), sensitivity.power, marker="o")
    ax.set(xlabel="True |log2 fold-change| bin", ylabel="Called proportion (power)", ylim=(0, 1.05), title="Sensitivity versus true effect size")
    ax.tick_params(axis="x", rotation=35); fig.tight_layout(); fig.savefig(out / "sensitivity_vs_effect.png", dpi=160); plt.close(fig)
    # A stable heatmap artifact uses the top 40 FDR-ranked genes.
    top = result.sort_values(["fdr", "pvalue"]).head(40).index
    fig, ax = plt.subplots(figsize=(9, 7))
    z = normalized.loc[top].sub(normalized.loc[top].mean(axis=1), axis=0).div(normalized.loc[top].std(axis=1), axis=0)
    im = ax.imshow(z, aspect="auto", cmap="coolwarm", interpolation="nearest")
    ax.set(title="Top genes (row z-scores; composition-aware normalization)", xlabel="Samples", ylabel="Genes")
    ax.set_xticks(range(z.shape[1]), z.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(z.shape[0]), z.index, fontsize=6)
    fig.colorbar(im, ax=ax, label="z-score"); fig.tight_layout(); fig.savefig(out / "top_genes_heatmap.png", dpi=160); plt.close(fig)

    rows = []
    for s in range(seed, seed + repeats):
        c, i, t = generate(s)
        _, m, _ = analyze_tables(c, i, t)
        m["seed"] = s; rows.append(m)
    repeated = pd.DataFrame(rows)
    repeated.to_csv(out / "simulation_metrics.csv", index=False)
    summary = {"seed": seed, "repeats": repeats, "one_run": one,
               "repeated_mean": {k: float(repeated[k].mean()) for k in ["empirical_fdr", "power", "bias_mean_log2fc", "coverage_95pct"]},
               "repeated_fdr_95th_percentile": float(repeated.empirical_fdr.quantile(0.95)),
               "design": {"normalization": "median-of-ratios", "filter": "CPM >= 5 in at least 6 samples", "test": "Welch t-test on log2 normalized counts", "fdr": "Benjamini-Hochberg"}}
    (out / "metrics.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


if __name__ == "__main__":
    p = argparse.ArgumentParser(); p.add_argument("--seed", type=int, default=2024); p.add_argument("--repeats", type=int, default=20); p.add_argument("--output-dir", default="outputs")
    args = p.parse_args(); print(json.dumps(run(args.seed, args.repeats, args.output_dir), indent=2))
