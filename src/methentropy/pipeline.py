"""Streaming genome-wide run: observed entropy vs depth-matched references."""

from __future__ import annotations

import logging
import multiprocessing as mp
import os
import tempfile
from collections import defaultdict
from contextlib import contextmanager

import numpy as np
import pandas as pd

from . import io
from .core import N_CPGS, bimodal_entropy, entropy_from_counts, minimum_entropy, per_cpg_methylated
from .nulls import (evaluate_independent, evaluate_pooled, independent_key, pooled_key,
                    simulate_pooled, summarize, empirical_pvalues)

log = logging.getLogger("methentropy")

RLM_OUTPUT_COLUMNS = [
    "chrom", "start", "end", "coverage", "methylation",
    "meth_cpg1", "meth_cpg2", "meth_cpg3", "meth_cpg4",
    "entropy", "entropy_min", "entropy_bimodal",
    "null_mean", "null_sd", "null_q05", "null_q95",
    "delta", "z", "scaled_entropy", "p_low", "p_high", "p_two", "q_two",
]
POOLED_COLUMNS = ["pooled_mean", "pooled_q05", "pooled_q95"]
LEGACY_OUTPUT_COLUMNS = [
    "chrom", "start", "end", "assumed_coverage", "methylation", "entropy",
    "entropy_min", "entropy_bimodal", "pooled_mean", "pooled_sd", "pooled_q05",
    "pooled_q95", "delta_pooled", "p_low", "p_high", "p_two", "q_two",
]


def bh_qvalues(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (NaN-aware)."""
    p = np.asarray(p, dtype=np.float64)
    q = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    m = pv.size
    if m == 0:
        return q
    order = np.argsort(pv)
    ranked = pv[order] * m / np.arange(1, m + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(ranked, 1.0)
    q[ok] = out
    return q


@contextmanager
def _pool(workers: int):
    if workers <= 1:
        yield None
        return
    ctx = mp.get_context("fork") if hasattr(os, "fork") else mp.get_context()
    pool = ctx.Pool(workers)
    try:
        yield pool
    finally:
        pool.close()
        pool.join()


def _map(pool, func, tasks):
    if pool is None:
        return map(func, tasks)
    return pool.imap_unordered(func, tasks, chunksize=max(1, len(tasks) // (pool._processes * 16)))


def _per_key_reference(n_arr, k_arr, func):
    cache, out = {}, np.empty(n_arr.size)
    for i, (n, k) in enumerate(zip(n_arr.tolist(), k_arr.tolist())):
        key = (n, k)
        if key not in cache:
            cache[key] = func(n, k)
        out[i] = cache[key]
    return out


def _process_rlm_chunk(chunk, pool, n_sims, seed, min_cov, max_cov, with_pooled, stats):
    counts = io.epiallele_matrix(chunk)
    n = counts.sum(axis=1)
    if "coverage" in chunk and not np.array_equal(n, chunk["coverage"].to_numpy()):
        raise ValueError("epiallele counts do not sum to the RLM coverage column; "
                         "is this an RLM entropy file?")
    keep = (n >= min_cov) & (n <= max_cov)
    stats["low_coverage"] += int((n < min_cov).sum())
    stats["high_coverage"] += int((n > max_cov).sum())
    chunk, counts, n = chunk.loc[keep].reset_index(drop=True), counts[keep], n[keep]
    if chunk.empty:
        return None

    m_cpg = per_cpg_methylated(counts)              # (rows, 4)
    k = m_cpg.sum(axis=1)
    observed = entropy_from_counts(counts)
    diff = np.abs(observed - chunk["rlm_entropy"].to_numpy(dtype=float))
    stats["max_entropy_diff"] = max(stats["max_entropy_diff"], float(np.nanmax(diff)))

    res = pd.DataFrame({
        "chrom": chunk["chrom"], "start": chunk["start"], "end": chunk["end"],
        "coverage": n, "methylation": k / (N_CPGS * n),
    })
    for j in range(N_CPGS):
        res[f"meth_cpg{j + 1}"] = m_cpg[:, j] / n
    res["entropy"] = observed
    res["entropy_min"] = _per_key_reference(n, k, minimum_entropy)
    res["entropy_bimodal"] = _per_key_reference(n, k, bimodal_entropy)

    # group rows sharing the same independence-null key
    groups = defaultdict(list)
    for i, key in enumerate(map(independent_key, n.tolist(), m_cpg.tolist())):
        groups[key].append(i)
    tasks = [(key, observed[idx], n_sims, seed) for key, idx in groups.items()]
    stats["null_keys"] += len(tasks)
    cols = {c: np.empty(len(res)) for c in ("null_mean", "null_sd", "null_q05", "null_q95",
                                              "p_low", "p_high")}
    for key, summ, p_low, p_high in _map(pool, evaluate_independent, tasks):
        idx = groups[key]
        cols["null_mean"][idx] = summ.mean
        cols["null_sd"][idx] = summ.sd
        cols["null_q05"][idx] = summ.q05
        cols["null_q95"][idx] = summ.q95
        cols["p_low"][idx] = p_low
        cols["p_high"][idx] = p_high
    for c, v in cols.items():
        res[c] = v
    res["delta"] = res["entropy"] - res["null_mean"]
    with np.errstate(divide="ignore", invalid="ignore"):
        res["z"] = np.where(res["null_sd"] > 0, res["delta"] / res["null_sd"], np.nan)
        denom = res["null_mean"] - res["entropy_min"]
        res["scaled_entropy"] = np.where(denom > 1e-9,
                                         (res["entropy"] - res["entropy_min"]) / denom, np.nan)
    res["p_two"] = np.minimum(1.0, 2 * np.minimum(res["p_low"], res["p_high"]))
    res["q_two"] = np.nan  # filled after the genome-wide pass

    if with_pooled:
        pk = [pooled_key(a, b) for a, b in zip(n.tolist(), k.tolist())]
        unique = list(dict.fromkeys(pk))
        summaries = dict((key, s) for key, s in
                         _map(pool, evaluate_pooled, [(u, n_sims, seed) for u in unique]))
        res["pooled_mean"] = [summaries[x].mean for x in pk]
        res["pooled_q05"] = [summaries[x].q05 for x in pk]
        res["pooled_q95"] = [summaries[x].q95 for x in pk]
    return res


def _evaluate_pooled_with_obs(task):
    key, observed, n_sims, seed = task
    sims = simulate_pooled(key, n_sims, seed)
    p_low, p_high = empirical_pvalues(sims, observed)
    return key, summarize(sims), p_low, p_high


def _process_legacy_chunk(chunk, pool, n_sims, seed, coverage, stats):
    m = chunk["rlm_mean_methylation"].to_numpy(dtype=float)
    obs = chunk["rlm_entropy"].to_numpy(dtype=float)
    keep = np.isfinite(m) & np.isfinite(obs) & (m >= 0) & (m <= 1)
    stats["invalid"] += int((~keep).sum())
    chunk, m, obs = chunk.loc[keep].reset_index(drop=True), m[keep], obs[keep]
    if chunk.empty:
        return None
    n = np.full(m.size, coverage, dtype=np.int64)
    k = np.rint(m * N_CPGS * coverage).astype(np.int64)
    res = pd.DataFrame({"chrom": chunk["chrom"], "start": chunk["start"], "end": chunk["end"],
                        "assumed_coverage": n, "methylation": m, "entropy": obs})
    res["entropy_min"] = _per_key_reference(n, k, minimum_entropy)
    res["entropy_bimodal"] = _per_key_reference(n, k, bimodal_entropy)
    groups = defaultdict(list)
    for i, key in enumerate(pooled_key(a, b) for a, b in zip(n.tolist(), k.tolist())):
        groups[key].append(i)
    cols = {c: np.empty(len(res)) for c in ("pooled_mean", "pooled_sd", "pooled_q05",
                                              "pooled_q95", "p_low", "p_high")}
    tasks = [(key, obs[idx], n_sims, seed) for key, idx in groups.items()]
    stats["null_keys"] += len(tasks)
    for key, s, p_low, p_high in _map(pool, _evaluate_pooled_with_obs, tasks):
        idx = groups[key]
        cols["pooled_mean"][idx], cols["pooled_sd"][idx] = s.mean, s.sd
        cols["pooled_q05"][idx], cols["pooled_q95"][idx] = s.q05, s.q95
        cols["p_low"][idx], cols["p_high"][idx] = p_low, p_high
    for c, v in cols.items():
        res[c] = v
    res["delta_pooled"] = res["entropy"] - res["pooled_mean"]
    res["p_two"] = np.minimum(1.0, 2 * np.minimum(res["p_low"], res["p_high"]))
    res["q_two"] = np.nan
    return res[LEGACY_OUTPUT_COLUMNS]


def run(input_path, output_path, n_sims=1000, seed=1, workers=1, min_coverage=10,
        max_coverage=1000, chunksize=1_000_000, with_pooled=False, legacy_coverage=None,
        input_format="auto"):
    """Process a whole file and write a tab-separated result (``.gz`` allowed)."""
    fmt = io.detect_format(input_path) if input_format == "auto" else input_format
    if fmt == "legacy" and legacy_coverage is None:
        raise ValueError(
            "Input looks like the legacy 5-column format, which has no read depth. "
            "Re-run RLM to get the full entropy output (recommended), or pass "
            "--legacy-coverage N to assume a fixed depth (pooled null only).")
    log.info("input format: %s", fmt)
    if fmt == "legacy":
        log.warning("legacy input: references use an ASSUMED depth of %d reads for every "
                    "window and only the pooled null; results are not depth-matched.",
                    legacy_coverage)

    stats = defaultdict(float)
    stats["max_entropy_diff"] = 0.0
    p_two_all = []
    out_dir = os.path.dirname(os.path.abspath(output_path))
    tmp = tempfile.NamedTemporaryFile("w", delete=False, suffix=".tsv", dir=out_dir,
                                      prefix=".methentropy_tmp_")
    tmp.close()
    header = True
    try:
        with _pool(workers) as pool:
            for i, chunk in enumerate(io.read_chunks(input_path, fmt, chunksize)):
                stats["rows_in"] += len(chunk)
                if fmt == "rlm":
                    res = _process_rlm_chunk(chunk, pool, n_sims, seed, min_coverage,
                                             max_coverage, with_pooled, stats)
                else:
                    res = _process_legacy_chunk(chunk, pool, n_sims, seed,
                                                int(legacy_coverage), stats)
                if res is None:
                    continue
                stats["rows_out"] += len(res)
                p_two_all.append(res["p_two"].to_numpy())
                res.to_csv(tmp.name, sep="\t", index=False, header=header,
                           mode="w" if header else "a", float_format="%.6g")
                header = False
                log.info("chunk %d: %d rows in, %d kept, %d null keys so far",
                         i + 1, len(chunk), len(res), stats["null_keys"])

        if header:  # nothing written
            cols = RLM_OUTPUT_COLUMNS if fmt == "rlm" else LEGACY_OUTPUT_COLUMNS
            pd.DataFrame(columns=cols).to_csv(output_path, sep="\t", index=False)
            return dict(stats)

        q = bh_qvalues(np.concatenate(p_two_all))
        offset, first = 0, True
        for part in pd.read_csv(tmp.name, sep="\t", chunksize=chunksize,
                                dtype={"chrom": str}):
            part["q_two"] = q[offset:offset + len(part)]
            offset += len(part)
            part.to_csv(output_path, sep="\t", index=False, header=first,
                        mode="w" if first else "a", float_format="%.6g",
                        compression="gzip" if str(output_path).endswith(".gz") else None)
            first = False
    finally:
        if os.path.exists(tmp.name):
            os.remove(tmp.name)

    if fmt == "rlm" and stats["max_entropy_diff"] > 1e-4:
        log.warning("recomputed entropy differs from the RLM entropy column by up to %.3g",
                    stats["max_entropy_diff"])
    log.info("done: %s", dict(stats))
    return dict(stats)
