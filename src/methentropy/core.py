"""Core quantities for 4-CpG methylation entropy.

Conventions (identical to RLM, Hetzel et al. 2021):

* An epiallele is the methylation pattern of one read over 4 consecutive CpGs.
* The 16 epialleles are indexed 0..15 as a 4-bit number where CpG 1 is the
  most significant bit and a set bit means methylated. Index 0 = ``gggg``
  (all unmethylated), index 1 = ``gggG`` (only CpG 4 methylated), ...,
  index 15 = ``GGGG``.
* Methylation entropy = (1/4) * sum_i -p_i log2 p_i over observed epialleles.
"""

from __future__ import annotations

from functools import lru_cache
from itertools import combinations

import numpy as np

N_CPGS = 4
N_PATTERNS = 2 ** N_CPGS
EPIALLELE_NAMES = [
    "".join("G" if (i >> (N_CPGS - 1 - j)) & 1 else "g" for j in range(N_CPGS))
    for i in range(N_PATTERNS)
]
# BITS[i, j] = 1 if epiallele i is methylated at CpG j (j = 0 is CpG 1)
BITS = np.array(
    [[(i >> (N_CPGS - 1 - j)) & 1 for j in range(N_CPGS)] for i in range(N_PATTERNS)],
    dtype=np.int64,
)
N_METH = BITS.sum(axis=1)  # methylated CpGs per epiallele


def entropy_from_counts(counts: np.ndarray, n_cpgs: int = N_CPGS) -> np.ndarray:
    """Per-CpG Shannon entropy for count vectors along the last axis.

    ``counts`` may be 1-D (one window) or 2-D (rows = windows or simulations).
    Rows with zero total return 0.
    """
    counts = np.asarray(counts, dtype=np.float64)
    total = counts.sum(axis=-1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = np.where(total > 0, counts / total, 0.0)
        terms = np.where(p > 0, -p * np.log2(p), 0.0)
    return terms.sum(axis=-1) / n_cpgs


def per_cpg_methylated(counts: np.ndarray) -> np.ndarray:
    """Number of methylated reads at each of the 4 CpGs, shape (..., 4)."""
    return np.asarray(counts, dtype=np.int64) @ BITS


def binary_entropy(p: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(p, dtype=np.float64), 0.0, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        h = -(p * np.log2(p) + (1 - p) * np.log2(1 - p))
    return np.nan_to_num(h, nan=0.0)


def _class_entropy(class_counts, n_reads: int, n_cpgs: int) -> float:
    h = 0.0
    for c in class_counts:
        if c > 0:
            q = c / n_reads
            h -= q * np.log2(q)
    return h / n_cpgs


@lru_cache(maxsize=None)
def minimum_entropy(n_reads: int, n_meth_sites: int, n_cpgs: int = N_CPGS) -> float:
    """Exact minimum methylation entropy for N reads and k methylated CpG calls.

    The minimum is reached when all reads carrying j methylated CpGs share one
    pattern, so the problem reduces to choosing how many reads carry
    j = 0..n_cpgs methylated CpGs (``c_j``) subject to ``sum c_j = N`` and
    ``sum j * c_j = k``. Entropy is concave, so the optimum uses few classes;
    we search every allocation over at most three classes, which is exact for
    this problem (verified against full enumeration in the test-suite).
    """
    n, k = int(n_reads), int(n_meth_sites)
    if n <= 0 or not 0 <= k <= n * n_cpgs:
        raise ValueError(f"invalid (n_reads={n}, n_meth_sites={k})")
    if k % n == 0:  # every read carries k/n methylated CpGs
        return 0.0
    best = np.inf
    classes = range(n_cpgs + 1)
    for a, b, c in combinations(classes, 3):
        # c_a + c_b + c_c = n ; a c_a + b c_b + c c_c = k ; free variable c_c
        for cc in range(0, n + 1):
            rem_n, rem_k = n - cc, k - c * cc
            # c_a + c_b = rem_n ; a c_a + b c_b = rem_k  ->  c_b = (rem_k - a rem_n)/(b - a)
            num = rem_k - a * rem_n
            if rem_n < 0 or num < 0 or num % (b - a):
                continue
            cb = num // (b - a)
            ca = rem_n - cb
            if ca < 0 or cb < 0:
                continue
            h = _class_entropy((ca, cb, cc), n, n_cpgs)
            if h < best:
                best = h
    return float(best)


def bimodal_entropy(n_reads: int, n_meth_sites: int, n_cpgs: int = N_CPGS) -> float:
    """Entropy if every read were fully methylated or fully unmethylated.

    This is the two-epiallele ("allelic" / two-population) reference: a
    fraction m of reads is GGGG and 1 - m is gggg, evaluated at the window's
    depth (number of GGGG reads rounded to the nearest integer).
    """
    n = int(n_reads)
    r = int(round(n_meth_sites / n_cpgs))
    r = min(max(r, 0), n)
    return _class_entropy((r, n - r), n, n_cpgs)
