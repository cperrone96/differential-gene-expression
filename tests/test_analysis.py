import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
sys.path.insert(0, str(ROOT))
from generate_data import generate  # noqa: E402
from analysis import analyze_tables, median_ratio_size_factors  # noqa: E402


def test_generator_is_deterministic_and_schema_is_stable():
    a = generate(17); b = generate(17)
    assert a[0].equals(b[0]) and a[1].equals(b[1]) and a[2].equals(b[2])
    assert a[0].shape == (2000, 12)
    assert set(a[1].group) == {"control", "disease"}
    assert int(a[2].is_de.sum()) == 180


def test_filter_and_composition_aware_normalization_are_exercised():
    counts, info, truth = generate(2024)
    factors = median_ratio_size_factors(counts)
    assert np.isfinite(factors).all() and not np.allclose(factors, 1)
    results, metrics, _ = analyze_tables(counts, info, truth)
    assert metrics["genes_retained"] < metrics["genes_input"]
    assert {"fdr", "ci_low", "ci_high", "true_log2fc"}.issubset(results.columns)
    assert 0 <= metrics["empirical_fdr"] <= 1
    assert 0 <= metrics["coverage_95pct"] <= 1
