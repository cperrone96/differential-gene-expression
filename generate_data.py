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
import numpy as np
import pandas as pd

RNG = np.random.default_rng(2024)

N_GENES = 2000
N_PER_GROUP = 6
N_DE = 180                      # true differentially-expressed genes (~9%)
DISPERSION = 0.15              # NB dispersion (variance = mu + phi*mu^2)

genes = [f"GENE{4000 + i}" for i in range(N_GENES)]
samples = [f"CTRL_{i+1}" for i in range(N_PER_GROUP)] + \
          [f"DIS_{i+1}"  for i in range(N_PER_GROUP)]
group = np.array(["control"] * N_PER_GROUP + ["disease"] * N_PER_GROUP)

# baseline expression level per gene (counts), heavy-tailed like real data
base_mean = np.exp(RNG.normal(4.0, 1.6, N_GENES))         # ~ median 55 counts
base_mean = np.clip(base_mean, 3, 30_000)

# per-sample library-size factor (sequencing depth differences)
lib_factor = RNG.normal(1.0, 0.15, len(samples)).clip(0.6, 1.5)

# ---------------------------------------------------------------------------
# Assign true DE genes and their log2 fold-changes (disease vs control)
# ---------------------------------------------------------------------------
is_de = np.zeros(N_GENES, dtype=bool)
de_idx = RNG.choice(N_GENES, N_DE, replace=False)
is_de[de_idx] = True

true_log2fc = np.zeros(N_GENES)
# effect sizes between |0.8| and |3.0| log2 units, half up / half down
signs = RNG.choice([-1, 1], N_DE)
mags = RNG.uniform(0.8, 3.0, N_DE)
true_log2fc[de_idx] = signs * mags

# ---------------------------------------------------------------------------
# Simulate negative-binomial counts
# ---------------------------------------------------------------------------
def nb_counts(mu):
    """Negative-binomial draw with mean mu and fixed dispersion."""
    mu = np.maximum(mu, 1e-6)
    size = 1.0 / DISPERSION                      # NB 'number of failures'
    p = size / (size + mu)
    return RNG.negative_binomial(size, p)

counts = np.zeros((N_GENES, len(samples)), dtype=int)
for s, (grp, lf) in enumerate(zip(group, lib_factor)):
    fc = np.where(grp == "disease", 2.0 ** true_log2fc, 1.0)
    mu = base_mean * fc * lf
    counts[:, s] = nb_counts(mu)

counts_df = pd.DataFrame(counts, index=genes, columns=samples)
counts_df.index.name = "gene"
counts_df.to_csv("data/counts.csv")

pd.DataFrame({"sample": samples, "group": group}).to_csv(
    "data/sample_info.csv", index=False)

pd.DataFrame({"gene": genes, "is_de": is_de.astype(int),
              "true_log2fc": true_log2fc.round(3)}).to_csv(
    "data/gene_truth.csv", index=False)

print(f"counts.csv: {N_GENES} genes x {len(samples)} samples")
print(f"true DE genes: {N_DE}  ({N_DE/N_GENES:.0%})  "
      f"| median library size {np.median(counts.sum(0)):,.0f} counts")
