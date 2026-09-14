# Differential gene expression: a reproducible synthetic benchmark

This is an educational simulation, not a biological result. It asks whether a
small case/control RNA-seq experiment can recover known effects under count
noise and a deliberate composition shift. The observation unit is one gene;
the target is the simulated disease-vs-control log2 fold-change. No human or
clinical data are included.

The primary analysis is a transparent Welch test on log2 counts after
median-of-ratios (DESeq-style) normalization. A low-count filter keeps genes
with CPM >= 5 in at least six of twelve samples. Benjamini-Hochberg is applied
to gene-level p-values, with calls requiring FDR < 0.05 and |log2FC| >= 1.
This is intentionally a teaching benchmark, not a replacement for a
negative-binomial GLM workflow such as DESeq2 or edgeR.

## Results (seed 2024; 10-seed check-in)

The checked-in machine-readable results are in `outputs/metrics.json` and
`outputs/simulation_metrics.csv`. In the current run, the single realization
made 129 calls (107 true positives, 22 false positives; empirical FDR 17.1%)
and power was 63.3%. Across twenty deterministic seeds, mean empirical FDR was
16.2%, mean power 65.0%, mean log2FC bias +0.01, and mean 95% interval coverage
94.6%. These results demonstrate why one favorable simulation cannot establish
FDR control; the composition shift and the simple Welch test produce more
false discoveries than the nominal threshold in this design.

`outputs/differential_expression_results.csv` contains estimates, p-values,
BH FDR, 95% Welch intervals, and truth labels. `outputs/sensitivity_by_true_effect.csv`
and `outputs/sensitivity_vs_effect.png` show power by true effect size;
`outputs/top_genes_heatmap.png` is a descriptive QC visualization, not proof
of a biological/global signal.

## Reproduce

```bash
python -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
python generate_data.py --seed 2024
python analysis.py --seed 2024 --repeats 20 --output-dir outputs
pytest -q
```

The generator accepts any integer seed and writes only synthetic counts,
sample metadata, and known truth. Tests verify deterministic generation,
schema/truth invariants, that the low-count filter removes rows, and that
composition-aware normalization plus interval metrics are exercised.

The notebook is a clean-kernel presentation of the same command-line analysis;
the script and CSV/JSON artifacts are canonical for reproducibility.

## Limitations and nonclaims

- Synthetic negative-binomial counts do not validate a real assay or biological
  pathway, and the truth labels are simulation truth rather than annotations.
- The Welch-on-log2-counts test is not a count-aware GLM and should not be used
  for production RNA-seq decisions.
- FDR, power, bias, and coverage are realization-specific summaries; the
  repeated-seed check quantifies this design only and is not a guarantee for
  another sample size, dispersion, or composition.
