"""Null calibration: does the method call 'signal' where none exists?

Windows are sampled from a real RLM file. For each, synthetic reads are drawn
with the same depth and the same per-CpG methylation levels but with CpGs
*independent* on each molecule (no linkage). The independence test should then
be calibrated (about 5% of windows at p < 0.05, or fewer because Monte Carlo
p-values on discrete data are conservative). For contrast, the original
v1/v2 rule (observed below the 5th percentile of a pooled null at a fixed 128
reads) is applied to the same synthetic data.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from . import io
from .core import N_CPGS, N_PATTERNS, entropy_from_counts, per_cpg_methylated
from .nulls import (empirical_pvalues, independent_key, pooled_key, simulate_independent,
                    simulate_pooled)

_WEIGHTS = (1 << np.arange(N_CPGS - 1, -1, -1)).astype(np.int64)


def _sample_rows(path, n_windows, min_cov, max_cov, rng, chunksize=1_000_000):
    """Reservoir-free approximation: keep a random subset of each chunk."""
    kept = []
    for chunk in io.read_chunks(path, "rlm", chunksize):
        counts = io.epiallele_matrix(chunk)
        n = counts.sum(axis=1)
        ok = np.flatnonzero((n >= min_cov) & (n <= max_cov))
        take = rng.choice(ok, size=min(ok.size, n_windows), replace=False) if ok.size else ok
        kept.append(counts[take])
    allc = np.concatenate(kept) if kept else np.empty((0, N_PATTERNS), dtype=np.int64)
    if allc.shape[0] > n_windows:
        allc = allc[rng.choice(allc.shape[0], n_windows, replace=False)]
    return allc


def calibrate(input_path, n_windows=2000, n_sims=1000, seed=1, min_coverage=10,
              max_coverage=1000, legacy_reads=128):
    if io.detect_format(input_path) != "rlm":
        raise ValueError("calibration needs the full RLM entropy output (epiallele counts)")
    rng = np.random.default_rng(seed)
    counts = _sample_rows(input_path, n_windows, min_coverage, max_coverage, rng)
    rows = []
    for c in counts:
        n = int(c.sum())
        p = per_cpg_methylated(c) / n
        reads = rng.random((n, N_CPGS)) < p                 # independent CpGs
        synth = np.bincount(reads.astype(np.int64) @ _WEIGHTS, minlength=N_PATTERNS)
        m = per_cpg_methylated(synth)
        obs = float(entropy_from_counts(synth))
        sims = simulate_independent(independent_key(n, m), n_sims, None if seed is None
                                    else seed + 7919)
        p_low, p_high = empirical_pvalues(sims, np.array([obs]))
        meth = m.sum() / (N_CPGS * n)
        k_old = int(round(meth * N_CPGS * legacy_reads))
        old = simulate_pooled(pooled_key(legacy_reads, k_old), 200, seed)
        rows.append((n, meth, obs, p_low[0], p_high[0], obs < np.percentile(old, 5),
                     obs > np.percentile(old, 95)))
    df = pd.DataFrame(rows, columns=["coverage", "methylation", "entropy", "p_low", "p_high",
                                     "old_below_band", "old_above_band"])
    df["p_two"] = np.minimum(1, 2 * np.minimum(df["p_low"], df["p_high"]))
    summary = {
        "windows": len(df),
        "median_coverage": float(df["coverage"].median()) if len(df) else np.nan,
        "frac_p_two_lt_0.05 (target <= 0.05)": float((df["p_two"] < 0.05).mean()),
        "frac_p_low_lt_0.05 (target <= 0.05)": float((df["p_low"] < 0.05).mean()),
        "frac_p_high_lt_0.05 (target <= 0.05)": float((df["p_high"] < 0.05).mean()),
        f"old_method_frac_below_band_{legacy_reads}reads (target ~0.05)":
            float(df["old_below_band"].mean()),
        f"old_method_frac_above_band_{legacy_reads}reads (target ~0.05)":
            float(df["old_above_band"].mean()),
    }
    return summary, df
