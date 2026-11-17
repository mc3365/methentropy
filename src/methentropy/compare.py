"""Region-level comparison of conditions with replicates (``methentropy compare``).

For every region (e.g. DMV) and every replicate result file, windows whose
first CpG lies inside the region (optionally extended by ``flank``) are pooled:

* ``delta``  - mean of observed minus null entropy over the windows
* ``scaled`` - (sum entropy - sum minimum) / (sum null - sum minimum): where the
  region sits between the minimum (0) and chance (1); a ratio of sums, so
  low-methylation windows with a near-zero denominator do not dominate
* ``methylation``, ``entropy``, ``null_mean``, ``n_windows``

Conditions are then compared with each other (every condition against the
first, the reference) in three ways:

1. regions as units: per region, replicate means; paired Wilcoxon signed-rank
   test on the per-region difference (other - reference);
2. replicates as units: per replicate, the mean over regions present in both
   conditions; Welch t-test (n = replicates, the conservative test);
3. per region: Welch t-test across replicates, Benjamini-Hochberg over regions
   (low power with 2-3 replicates; reported as counts).
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd

from .pipeline import bh_qvalues

log = logging.getLogger("methentropy")

METRICS = ("delta", "scaled")
METRIC_LABEL = {"delta": "Observed − null entropy",
                "scaled": "Scaled entropy (0 = minimum, 1 = chance)"}


def _stats():
    try:
        from scipy import stats
    except ImportError as err:  # pragma: no cover
        raise ValueError("compare needs scipy: pip install 'methentropy[plot]'") from err
    return stats


def read_named_regions(path, names_path=None, name_column=4):
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    out = pd.DataFrame({"chrom": df[0].astype(str), "start": df[1].astype(np.int64),
                        "end": df[2].astype(np.int64)})
    if df.shape[1] >= 4:
        out["name"] = df[3].astype(str)
    else:
        out["name"] = out["chrom"] + ":" + out["start"].astype(str) + "-" + out["end"].astype(str)
    if names_path:
        if df.shape[1] < name_column:
            raise ValueError(f"{path} has no column {name_column} to match --names against")
        keep = set(pd.read_csv(names_path, sep=r"\s+", header=None, dtype=str)[0])
        out = out.loc[df[name_column - 1].astype(str).isin(keep).to_numpy()]
    out = out[out["end"] > out["start"]].drop_duplicates(["chrom", "start", "end"])
    if out.empty:
        raise ValueError("no regions left after filtering")
    return out.reset_index(drop=True)


def region_summaries(result_path, regions, flank=0, chunksize=2_000_000):
    """Per-region sums over windows of one result file -> DataFrame (index = region)."""
    header = pd.read_csv(result_path, sep="\t", nrows=0).columns
    pick = lambda *c: next((x for x in c if x in header), None)  # noqa: E731
    cols = {"entropy": "entropy", "null": pick("null_mean", "pooled_mean"),
            "minimum": "entropy_min", "methylation": "methylation",
            "delta": pick("delta", "delta_pooled")}
    if None in cols.values() or any(c not in header for c in cols.values()):
        raise ValueError(f"{result_path}: not a methentropy result file")
    keys = list(cols)
    sums = np.zeros((len(regions), len(keys)))
    n = np.zeros(len(regions))
    by_chrom = {c: g for c, g in regions.groupby("chrom")}
    for chunk in pd.read_csv(result_path, sep="\t", usecols=["chrom", "start", *cols.values()],
                             dtype={"chrom": str}, chunksize=chunksize):
        for chrom, w in chunk.groupby("chrom", sort=False):
            regs = by_chrom.get(chrom)
            if regs is None:
                continue
            w = w.sort_values("start")
            pos = w["start"].to_numpy()
            vals = w[[cols[k] for k in keys]].to_numpy(dtype=float)
            ok_rows = np.isfinite(vals).all(axis=1)
            for r in regs.itertuples():
                lo = np.searchsorted(pos, r.start - flank, "left")
                hi = np.searchsorted(pos, r.end + flank, "left")
                if hi > lo:
                    v = vals[lo:hi][ok_rows[lo:hi]]
                    sums[r.Index] += v.sum(axis=0)
                    n[r.Index] += len(v)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = sums / n[:, None]
        k = {key: i for i, key in enumerate(keys)}
        out = pd.DataFrame({"n_windows": n, "methylation": mean[:, k["methylation"]],
                            "entropy": mean[:, k["entropy"]], "null_mean": mean[:, k["null"]],
                            "delta": mean[:, k["delta"]]}, index=regions.index)
        denom = sums[:, k["null"]] - sums[:, k["minimum"]]
        out["scaled"] = np.where(denom > 1e-6 * np.maximum(n, 1),
                                 (sums[:, k["entropy"]] - sums[:, k["minimum"]]) / denom, np.nan)
    out.loc[n == 0, :] = np.nan
    out["n_windows"] = n
    return out


def compare(samples, regions_path, out_prefix, names_path=None, name_column=4, flank=0,
            min_windows=3, title=None, region_label="region", formats=("png", "pdf")):
    """samples: ordered {condition: [result files]}; first condition is the reference."""
    stats = _stats()
    regions = read_named_regions(regions_path, names_path, name_column)
    conds = list(samples)
    if len(conds) < 2:
        raise ValueError("compare needs at least two conditions (-s ... -s ...)")
    log.info("compare: %d %ss, conditions %s", len(regions), region_label, conds)

    reps = {}
    for c in conds:
        reps[c] = []
        for f in samples[c]:
            log.info("  %s: %s", c, f)
            s = region_summaries(f, regions, flank)
            s.loc[s["n_windows"] < min_windows, ["delta", "scaled", "methylation", "entropy",
                                                 "null_mean"]] = np.nan
            reps[c].append(s)

    table = regions.copy()
    for c in conds:
        for m in ("delta", "scaled", "methylation"):
            stack = np.vstack([r[m].to_numpy() for r in reps[c]])
            with np.errstate(all="ignore"), _quiet():
                table[f"{m}__{c}"] = np.nanmean(stack, axis=0)
        table[f"n_windows__{c}"] = np.nanmean(np.vstack([r["n_windows"] for r in reps[c]]), 0)

    ref = conds[0]
    summary_rows = []
    for c in conds[1:]:
        for m in METRICS:
            a, b = table[f"{m}__{ref}"], table[f"{m}__{c}"]
            both = a.notna() & b.notna()
            diff = (b - a)[both]
            table[f"diff_{m}__{c}_vs_{ref}"] = b - a
            # 1) regions as units
            w_p = stats.wilcoxon(diff).pvalue if len(diff) >= 10 and np.any(diff != 0) else np.nan
            # 2) replicates as units (mean over regions present in both conditions)
            ra = np.array([np.nanmean(r[m].to_numpy()[both.to_numpy()]) for r in reps[ref]])
            rb = np.array([np.nanmean(r[m].to_numpy()[both.to_numpy()]) for r in reps[c]])
            t_p = (stats.ttest_ind(rb, ra, equal_var=False).pvalue
                   if len(ra) >= 2 and len(rb) >= 2 else np.nan)
            # 3) per-region replicate tests
            A = np.vstack([r[m].to_numpy() for r in reps[ref]])
            B = np.vstack([r[m].to_numpy() for r in reps[c]])
            with np.errstate(all="ignore"), _quiet():
                p_reg = stats.ttest_ind(B, A, axis=0, equal_var=False, nan_policy="omit").pvalue
            p_reg = np.asarray(p_reg, dtype=float)
            p_reg[~both.to_numpy()] = np.nan
            q_reg = bh_qvalues(p_reg)
            table[f"p_{m}__{c}_vs_{ref}"] = p_reg
            table[f"q_{m}__{c}_vs_{ref}"] = q_reg
            d = (b - a).to_numpy()
            summary_rows.append({
                "comparison": f"{c} vs {ref}", "metric": m, "n_regions": int(both.sum()),
                "ref_median": float(np.nanmedian(a[both])), "other_median": float(np.nanmedian(b[both])),
                "median_diff": float(np.median(diff)) if len(diff) else np.nan,
                "frac_regions_lower": float((diff < 0).mean()) if len(diff) else np.nan,
                "wilcoxon_p_regions": w_p,
                "ref_replicate_means": ";".join(f"{x:.4f}" for x in ra),
                "other_replicate_means": ";".join(f"{x:.4f}" for x in rb),
                "welch_p_replicates": t_p,
                "n_regions_q05_lower": int(np.sum((q_reg < 0.05) & (d < 0))),
                "n_regions_q05_higher": int(np.sum((q_reg < 0.05) & (d > 0))),
                "n_regions_p05_lower": int(np.sum((p_reg < 0.05) & (d < 0))),
                "n_regions_p05_higher": int(np.sum((p_reg < 0.05) & (d > 0))),
            })
    summary = pd.DataFrame(summary_rows)
    table.to_csv(f"{out_prefix}_regions.tsv", sep="\t", index=False, float_format="%.6g")
    summary.to_csv(f"{out_prefix}_summary.tsv", sep="\t", index=False, float_format="%.6g")
    paths = _plot(table, summary, reps, conds, out_prefix, title, region_label, formats)
    return summary, table, paths


class _quiet:
    def __enter__(self):
        import warnings
        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def __exit__(self, *exc):
        self._w.__exit__(*exc)


def _fmt_p(p):
    if not np.isfinite(p):
        return "n/a"
    return f"p = {p:.2g}" if p >= 1e-4 else "p < 1e-4"


def _plot(table, summary, reps, conds, out_prefix, title, region_label, formats):
    from .plotting import BLUE, INK, INK_2, MUTED, SERIES, _mpl, _save
    plt = _mpl()
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2),
                             gridspec_kw={"width_ratios": [1, 1, 1.1]})
    rng = np.random.default_rng(0)
    for ax, m in zip(axes[:2], METRICS):
        data = [table[f"{m}__{c}"].dropna().to_numpy() for c in conds]
        bp = ax.boxplot(data, widths=0.5, showfliers=False, patch_artist=True,
                        medianprops=dict(color=INK, lw=1.4),
                        whiskerprops=dict(color=MUTED), capprops=dict(color=MUTED),
                        boxprops=dict(edgecolor=MUTED, lw=0.8))
        for j, patch in enumerate(bp["boxes"]):
            patch.set_facecolor(SERIES[j % len(SERIES)])
            patch.set_alpha(0.25)
        for j, c in enumerate(conds):  # replicate means (over regions) as points
            vals = [np.nanmean(r[m]) for r in reps[c]]
            ax.scatter(j + 1 + rng.uniform(-0.08, 0.08, len(vals)), vals, s=26, zorder=3,
                       color=SERIES[j % len(SERIES)], edgecolor="white", lw=1)
        if m == "delta":
            ax.axhline(0, color=INK_2, lw=0.9)
        ax.set_xticks(range(1, len(conds) + 1), conds)
        ax.set_ylabel(METRIC_LABEL[m])
        ax.set_title(f"Per-{region_label} mean (box) and replicate means (dots)", loc="left",
                     fontsize=8.5, color=INK_2)
        ylo, yhi = ax.get_ylim()
        span = yhi - ylo
        ax.set_ylim(ylo, yhi + 0.12 * span * (len(conds) - 1))
        for j, c in enumerate(conds[1:], start=2):
            row = summary[(summary.comparison == f"{c} vs {conds[0]}") & (summary.metric == m)]
            if row.empty:
                continue
            y = yhi + 0.1 * span * (j - 2) + 0.02 * span
            ax.plot([1, 1, j, j], [y - 0.015 * span, y, y, y - 0.015 * span], color=MUTED, lw=0.8)
            r = row.iloc[0]
            ax.text((1 + j) / 2, y + 0.01 * span,
                    f"{region_label}s: {_fmt_p(r.wilcoxon_p_regions)} · reps: "
                    f"{_fmt_p(r.welch_p_replicates)}", ha="center", va="bottom",
                    fontsize=7.5, color=INK_2)
    ax = axes[2]
    last = conds[-1]
    x = table[f"methylation__{last}"] - table[f"methylation__{conds[0]}"]
    y = table[f"diff_delta__{last}_vs_{conds[0]}"]
    ok = x.notna() & y.notna()
    ax.scatter(x[ok], y[ok], s=8, color=BLUE, alpha=0.45, lw=0)
    ax.axhline(0, color=INK_2, lw=0.9)
    ax.axvline(0, color=INK_2, lw=0.9)
    ax.set_xlabel(f"Δ methylation ({last} − {conds[0]})")
    ax.set_ylabel(f"Δ (observed − null) ({last} − {conds[0]})")
    ax.set_title(f"Per-{region_label} change; below 0 = more ordered", loc="left", fontsize=8.5,
                 color=INK_2)
    fig.suptitle(title or "Condition comparison", x=0.01, ha="left", fontsize=11, color=INK)
    fig.tight_layout()
    paths = _save(fig, f"{out_prefix}_compare", formats)
    plt.close(fig)
    return paths
