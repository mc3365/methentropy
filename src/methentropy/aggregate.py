"""Summarise per-window results over genomic regions (DMVs, CGIs, promoters...)."""

from __future__ import annotations

import numpy as np
import pandas as pd

AGG_COLUMNS = ["chrom", "start", "end", "name", "n_windows", "mean_coverage",
               "mean_methylation", "mean_entropy", "mean_null", "mean_delta",
               "weighted_delta", "median_scaled_entropy", "frac_low_sig", "frac_high_sig"]


def read_regions(path) -> pd.DataFrame:
    regions = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    if regions.shape[1] < 3:
        raise ValueError("regions file must have at least chrom, start, end")
    regions = regions.iloc[:, :4] if regions.shape[1] >= 4 else regions.iloc[:, :3]
    regions.columns = ["chrom", "start", "end", "name"][:regions.shape[1]]
    if "name" not in regions:
        regions["name"] = (regions["chrom"] + ":" + regions["start"].astype(str) + "-"
                           + regions["end"].astype(str))
    return regions


def aggregate(results_path, regions_path, output_path, q_threshold=0.05):
    """A window belongs to a region if its first CpG (``start``) lies in [start, end).

    Columns: ``mean_delta`` (observed minus independence-null entropy, averaged
    over windows), ``weighted_delta`` (coverage-weighted), ``frac_low_sig`` /
    ``frac_high_sig`` (fraction of windows with q_two < threshold and entropy
    below / above the null).
    """
    res = pd.read_csv(results_path, sep="\t", dtype={"chrom": str})
    delta_col = "delta" if "delta" in res else "delta_pooled"
    null_col = "null_mean" if "null_mean" in res else "pooled_mean"
    cov_col = "coverage" if "coverage" in res else "assumed_coverage"
    scaled = res["scaled_entropy"] if "scaled_entropy" in res else pd.Series(np.nan, index=res.index)
    sig = res["q_two"].to_numpy() < q_threshold
    res = res.assign(_low=(sig & (res[delta_col] < 0)).astype(float),
                     _high=(sig & (res[delta_col] > 0)).astype(float),
                     _scaled=scaled, _wd=res[delta_col] * res[cov_col])
    regions = read_regions(regions_path)
    rows = []
    by_chrom = {c: g.sort_values("start") for c, g in res.groupby("chrom", sort=False)}
    for reg in regions.itertuples(index=False):
        g = by_chrom.get(reg.chrom)
        if g is None:
            rows.append((reg.chrom, reg.start, reg.end, reg.name, 0) + (np.nan,) * 9)
            continue
        pos = g["start"].to_numpy()
        lo, hi = np.searchsorted(pos, reg.start, "left"), np.searchsorted(pos, reg.end, "left")
        w = g.iloc[lo:hi]
        if w.empty:
            rows.append((reg.chrom, reg.start, reg.end, reg.name, 0) + (np.nan,) * 9)
            continue
        rows.append((reg.chrom, reg.start, reg.end, reg.name, len(w),
                     w[cov_col].mean(), w["methylation"].mean(), w["entropy"].mean(),
                     w[null_col].mean(), w[delta_col].mean(),
                     w["_wd"].sum() / w[cov_col].sum(), w["_scaled"].median(),
                     w["_low"].mean(), w["_high"].mean()))
    out = pd.DataFrame(rows, columns=AGG_COLUMNS)
    out.to_csv(output_path, sep="\t", index=False, float_format="%.6g")
    return out
