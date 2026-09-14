"""
Generate SYNTHETIC RNA-seq count data for this portfolio project.

None of this is real sequencing data. It is fabricated with a fixed random seed
so the analysis is fully reproducible, while realistically mirroring a small
case-control gene-expression experiment:

  * two groups (e.g. disease vs. control), 6 biological replicates each,
  * ~2,000 genes with realistic baseline expression and per-sample library-size
    variation,
  * negative-binomial count noise (over-dispersed, as real RNA-seq is),
  * a known subset of truly differentially-expressed (DE) genes with defined
    log2 fold-changes, so the analysis can be validated against ground truth.

Run:  python generate_data.py
Writes: data/counts.csv         (genes x samples raw counts)
        data/sample_info.csv    (sample -> group)
        data/gene_truth.csv     (per-gene true DE status + log2FC)
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

N_GENES = 2000
N_PER_GROUP = 6
N_DE = 180                      # true differentially-expressed genes (~9%)
DISPERSION = 0.15              # NB dispersion (variance = mu + phi*mu^2)

def nb_counts(mu, rng):
    """Negative-binomial draw with mean mu and fixed dispersion."""
    mu = np.maximum(mu, 1e-6)
    size = 1.0 / DISPERSION                      # NB 'number of failures'
    p = size / (size + mu)
    return rng.negative_binomial(size, p)


def generate(seed=2024, output_dir="data"):
    """Generate one deterministic synthetic experiment and return its tables.

    The returned tables make simulation-based validation possible without
    writing intermediate files.  Counts are negative-binomial and include a
    mild composition shift so normalization is part of the analysis rather
    than an assumption.
    """
    rng = np.random.default_rng(seed)
    genes = [f"GENE{4000 + i}" for i in range(N_GENES)]
    samples = [f"CTRL_{i+1}" for i in range(N_PER_GROUP)] + [f"DIS_{i+1}" for i in range(N_PER_GROUP)]
    group = np.array(["control"] * N_PER_GROUP + ["disease"] * N_PER_GROUP)

    base_mean = np.exp(rng.normal(4.0, 1.6, N_GENES))
    base_mean = np.clip(base_mean, 3, 30_000)
    lib_factor = rng.normal(1.0, 0.15, len(samples)).clip(0.6, 1.5)

    is_de = np.zeros(N_GENES, dtype=bool)
    de_idx = rng.choice(N_GENES, N_DE, replace=False)
    is_de[de_idx] = True
    true_log2fc = np.zeros(N_GENES)
    signs = rng.choice([-1, 1], N_DE)
    mags = rng.uniform(0.8, 3.0, N_DE)
    true_log2fc[de_idx] = signs * mags

    # A small set of highly expressed genes is shifted in disease.  This is
    # deliberate compositional distortion; a library-size-only method can be
    # biased by it, while median-of-ratios is at least a defensible benchmark.
    comp_idx = rng.choice(np.flatnonzero(~is_de), 20, replace=False)
    composition_shift = np.ones(N_GENES)
    composition_shift[comp_idx] = 8.0
    counts = np.zeros((N_GENES, len(samples)), dtype=int)
    for s, (grp, lf) in enumerate(zip(group, lib_factor)):
        fc = np.where(grp == "disease", 2.0 ** true_log2fc, 1.0)
        fc = fc * (composition_shift if grp == "disease" else 1.0)
        counts[:, s] = nb_counts(base_mean * fc * lf, rng)

    return (
        pd.DataFrame(counts, index=genes, columns=samples).rename_axis("gene"),
        pd.DataFrame({"sample": samples, "group": group}),
        pd.DataFrame({"gene": genes, "is_de": is_de.astype(int), "true_log2fc": true_log2fc.round(3)}),
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=2024)
    parser.add_argument("--output-dir", default="data")
    args = parser.parse_args()
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    counts_df, sample_info, gene_truth = generate(args.seed, out)
    counts_df.to_csv(out / "counts.csv")
    sample_info.to_csv(out / "sample_info.csv", index=False)
    gene_truth.to_csv(out / "gene_truth.csv", index=False)
    print(f"counts.csv: {counts_df.shape[0]} genes x {counts_df.shape[1]} samples")
    print(f"true DE genes: {gene_truth.is_de.sum()} ({gene_truth.is_de.mean():.0%}) | median library size {np.median(counts_df.sum(0)):,.0f} counts")


if __name__ == "__main__":
    main()
