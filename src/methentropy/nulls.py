"""Depth-matched null models for 4-CpG methylation entropy.

Two Monte Carlo nulls are provided. Both are evaluated at the window's own
read depth N, so finite-sample bias of the plug-in entropy is shared by the
observation and the null.

``independent`` (recommended)
    Conditional permutation null. Each CpG keeps its observed number of
    methylated reads; reads are shuffled independently at every CpG. This
    breaks linkage between CpGs on the same molecule while preserving the
    per-CpG methylation levels exactly. Observed entropy *below* this null
    means CpGs are co-methylated on the same molecules more than chance.

``pooled``
    The original (v1/v2) null: all 4N CpG calls are pooled and shuffled,
    keeping only the total number of methylated calls. It ignores differences
    between the four CpGs, so it is retained for comparison only.

Null distributions depend only on a small key, so results are cached and
every key gets a deterministic random stream derived from (seed, key): the
output does not depend on chunking or the number of workers.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import N_CPGS, N_PATTERNS, entropy_from_counts

_WEIGHTS = (1 << np.arange(N_CPGS - 1, -1, -1)).astype(np.int64)  # CpG1 = MSB


def independent_key(n_reads: int, meth_per_cpg) -> tuple:
    """Canonical key: entropy is invariant to CpG order and to relabelling
    methylated/unmethylated at any single CpG."""
    folded = sorted(min(int(k), int(n_reads) - int(k)) for k in meth_per_cpg)
    return (int(n_reads), *folded)


def pooled_key(n_reads: int, n_meth_sites: int) -> tuple:
    total = int(n_reads) * N_CPGS
    return (int(n_reads), min(int(n_meth_sites), total - int(n_meth_sites)))


def _rng(seed: int | None, kind: int, key: tuple) -> np.random.Generator:
    if seed is None:
        return np.random.default_rng()
    return np.random.default_rng(np.random.SeedSequence([int(seed), kind, *key]))


def _entropy_of_patterns(meth: np.ndarray) -> np.ndarray:
    """meth: bool array (n_sims, n_reads, 4) -> entropy per simulation."""
    n_sims = meth.shape[0]
    codes = meth.astype(np.int64) @ _WEIGHTS  # (n_sims, n_reads)
    offset = codes + N_PATTERNS * np.arange(n_sims)[:, None]
    counts = np.bincount(offset.ravel(), minlength=n_sims * N_PATTERNS)
    return entropy_from_counts(counts.reshape(n_sims, N_PATTERNS))


def _choose(u: np.ndarray, k: int) -> np.ndarray:
    """Boolean mask selecting exactly k uniformly random positions per row of u."""
    if k <= 0:
        return np.zeros(u.shape, dtype=bool)
    if k >= u.shape[1]:
        return np.ones(u.shape, dtype=bool)
    threshold = np.partition(u, k - 1, axis=1)[:, k - 1:k]
    return u <= threshold  # float64 draws: ties have probability ~0


def simulate_independent(key: tuple, n_sims: int, seed: int | None,
                         batch: int = 500) -> np.ndarray:
    n_reads, ks = key[0], key[1:]
    rng = _rng(seed, 1, key)
    out = np.empty(n_sims)
    for start in range(0, n_sims, batch):
        b = min(batch, n_sims - start)
        meth = np.empty((b, n_reads, N_CPGS), dtype=bool)
        for j, k in enumerate(ks):  # independent shuffle at each CpG
            meth[:, :, j] = _choose(rng.random((b, n_reads)), int(k))
        out[start:start + b] = _entropy_of_patterns(meth)
    return out


def simulate_pooled(key: tuple, n_sims: int, seed: int | None,
                    batch: int = 500) -> np.ndarray:
    n_reads, k = key
    rng = _rng(seed, 2, key)
    out = np.empty(n_sims)
    for start in range(0, n_sims, batch):
        b = min(batch, n_sims - start)
        meth = _choose(rng.random((b, n_reads * N_CPGS)), int(k)).reshape(b, n_reads, N_CPGS)
        out[start:start + b] = _entropy_of_patterns(meth)
    return out


@dataclass
class NullSummary:
    mean: float
    sd: float
    q05: float
    q95: float


def summarize(sims: np.ndarray) -> NullSummary:
    q05, q95 = np.percentile(sims, [5, 95])
    return NullSummary(float(sims.mean()), float(sims.std(ddof=1)) if sims.size > 1 else 0.0,
                       float(q05), float(q95))


def empirical_pvalues(sims: np.ndarray, observed: np.ndarray, tol: float = 1e-9):
    """Lower- and upper-tail Monte Carlo p-values with the +1 correction.

    p_low  = P(null <= observed): small when entropy is lower than chance
             (CpGs co-methylated / ordered patterns).
    p_high = P(null >= observed): small when entropy is higher than chance.
    """
    s = np.sort(sims)
    n = s.size
    obs = np.asarray(observed, dtype=np.float64)
    n_le = np.searchsorted(s, obs + tol, side="right")
    n_ge = n - np.searchsorted(s, obs - tol, side="left")
    return (n_le + 1) / (n + 1), (n_ge + 1) / (n + 1)


def evaluate_independent(task):
    """Worker entry point: (key, observed_array, n_sims, seed) -> results."""
    key, observed, n_sims, seed = task
    sims = simulate_independent(key, n_sims, seed)
    p_low, p_high = empirical_pvalues(sims, observed)
    return key, summarize(sims), p_low, p_high


def evaluate_pooled(task):
    key, n_sims, seed = task
    return key, summarize(simulate_pooled(key, n_sims, seed))
