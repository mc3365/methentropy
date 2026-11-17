"""Figures: feature metaplots (the v1 design) and per-sample QC.

Metaplot (``methentropy metaplot``) follows the original
``example_mef_metaplots.Rmd`` design: every region is scaled to a fixed number
of body bins and extended by fixed-size flanking bins (defaults reproduce
EnrichedHeatmap ``extend = 2000, w = 200, target_ratio = 0.7``: 10 + 47 + 10
bins). For every region and bin the windows whose first CpG falls in the bin
are averaged; replicates of one condition are averaged per region and bin; the
plot shows the median across regions with the 25-75% range shaded.

QC (``methentropy qc``) gives four panels for one result file: entropy vs
methylation with the references, the distribution of observed-minus-null,
observed-minus-null by read depth, and the p-value histogram.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass

import numpy as np
import pandas as pd

log = logging.getLogger("methentropy")

# validated categorical slots (light surface), fixed order; neutrals for references
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
SERIES = [BLUE, ORANGE, AQUA, "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK_2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"

# profiled quantity -> (result column(s) to try, label, colour, line style)
PROFILE_SERIES = {
    "entropy": (("entropy",), "Observed", BLUE, "-"),
    "null_mean": (("null_mean", "pooled_mean"), "Stochastic (null)", ORANGE, "-"),
    "bimodal": (("entropy_bimodal",), "Bimodal / allelic", AQUA, "-"),
    "minimum": (("entropy_min",), "Minimum", MUTED, "--"),
    "methylation": (("methylation",), "Methylation", INK, ":"),
}
DELTA_COLUMNS = ("delta", "delta_pooled")


def _mpl():
    try:
        import matplotlib
    except ImportError as err:  # pragma: no cover
        raise ValueError("plotting needs matplotlib: pip install 'methentropy[plot]'") from err
    matplotlib.use("Agg")
    logging.getLogger("fontTools").setLevel(logging.WARNING)
    logging.getLogger("matplotlib").setLevel(logging.WARNING)
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK_2,
        "xtick.color": INK_2, "ytick.color": INK_2, "axes.titlesize": 10,
        "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "legend.frameon": False, "savefig.dpi": 200, "lines.linewidth": 1.6,
        "pdf.fonttype": 42,
    })
    return plt


def _save(fig, prefix, formats):
    paths = []
    for fmt in formats:
        path = f"{prefix}.{fmt}"
        fig.savefig(path, bbox_inches="tight")
        paths.append(path)
    return paths


def _end_labels(ax, items, min_gap_frac=0.06):
    """Direct labels at the right end of lines, nudged apart so none overlap.

    items: list of (text, y). The gap is a fraction of the y-axis span.
    """
    items = [(t, y) for t, y in items if np.isfinite(y)]
    if not items:
        return
    lo, hi = ax.get_ylim()
    gap = (hi - lo) * min_gap_frac
    items.sort(key=lambda it: it[1])
    placed = []
    for t, y in items:
        if placed and y < placed[-1][1] + gap:
            y = placed[-1][1] + gap
        placed.append((t, y))
    overflow = placed[-1][1] - (hi - gap / 2)
    if overflow > 0:
        placed = [(t, y - overflow) for t, y in placed]
    x = ax.get_xlim()[1]
    for t, y in placed:
        ax.annotate(t, (x, y), xytext=(5, 0), textcoords="offset points", va="center",
                    fontsize=7.5, color=INK_2, annotation_clip=False)


def _smooth(y, k):
    if k <= 1:
        return y
    s = pd.Series(y)
    return s.rolling(k, center=True, min_periods=1).mean().to_numpy()


# ----------------------------------------------------------------- regions
def read_regions(path, names_path=None, name_column=4):
    """BED with chrom, start, end; optional subset by names in one column."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    regions = pd.DataFrame({"chrom": df[0].astype(str), "start": df[1].astype(np.int64),
                            "end": df[2].astype(np.int64)})
    if names_path:
        if df.shape[1] < name_column:
            raise ValueError(f"{path} has no column {name_column} to match --names against")
        keep = set(pd.read_csv(names_path, sep=r"\s+", header=None, dtype=str)[0])
        mask = df[name_column - 1].astype(str).isin(keep).to_numpy()
        regions = regions.loc[mask].reset_index(drop=True)
        if regions.empty:
            raise ValueError(f"no region in {path} column {name_column} matches {names_path}")
    regions = regions[regions["end"] > regions["start"]].reset_index(drop=True)
    return regions


@dataclass
class BinLayout:
    flank: int = 2000
    bin_size: int = 200
    body_bins: int = 47

    @property
    def n_flank(self):
        return int(np.ceil(self.flank / self.bin_size))

    @property
    def n_bins(self):
        return 2 * self.n_flank + self.body_bins

    def assign(self, pos, start, end):
        """Bin index (or -1) for window positions relative to one region."""
        nf = self.n_flank
        out = np.full(pos.shape, -1, dtype=np.int64)
        up = (pos >= start - self.flank) & (pos < start)
        out[up] = nf - 1 - ((start - 1 - pos[up]) // self.bin_size)
        body = (pos >= start) & (pos < end)
        out[body] = nf + np.minimum(((pos[body] - start) * self.body_bins) // (end - start),
                                    self.body_bins - 1)
        down = (pos >= end) & (pos < end + self.flank)
        out[down] = nf + self.body_bins + (pos[down] - end) // self.bin_size
        out[(out < 0) | (out >= self.n_bins)] = -1
        return out


def _resolve_columns(header):
    cols = {}
    for key, (candidates, *_rest) in PROFILE_SERIES.items():
        for c in candidates:
            if c in header:
                cols[key] = c
                break
    for c in DELTA_COLUMNS:
        if c in header:
            cols["delta"] = c
            break
    if "entropy" not in cols:
        raise ValueError("result file has no 'entropy' column; is it methentropy output?")
    return cols


def region_bin_means(result_path, regions, layout, chunksize=2_000_000):
    """Mean of each profiled quantity per region and bin for one result file.

    Returns {quantity: array (n_regions, n_bins)} with NaN where no window fell.
    """
    header = pd.read_csv(result_path, sep="\t", nrows=0).columns
    cols = _resolve_columns(header)
    keys = list(cols)
    n_reg, n_bins = len(regions), layout.n_bins
    sums = np.zeros((len(keys), n_reg, n_bins))
    counts = np.zeros((n_reg, n_bins))
    by_chrom = {c: g for c, g in regions.reset_index().groupby("chrom")}
    usecols = ["chrom", "start", *cols.values()]
    for chunk in pd.read_csv(result_path, sep="\t", usecols=usecols, dtype={"chrom": str},
                             chunksize=chunksize):
        for chrom, wins in chunk.groupby("chrom", sort=False):
            regs = by_chrom.get(chrom)
            if regs is None:
                continue
            wins = wins.sort_values("start")
            pos = wins["start"].to_numpy()
            vals = wins[[cols[k] for k in keys]].to_numpy(dtype=float)
            for r in regs.itertuples():
                lo = np.searchsorted(pos, r.start - layout.flank, "left")
                hi = np.searchsorted(pos, r.end + layout.flank, "left")
                if hi <= lo:
                    continue
                b = layout.assign(pos[lo:hi], r.start, r.end)
                ok = b >= 0
                if not ok.any():
                    continue
                bins, v = b[ok], vals[lo:hi][ok]
                counts[r.Index] += np.bincount(bins, minlength=n_bins)
                for j in range(len(keys)):
                    good = np.isfinite(v[:, j])
                    sums[j, r.Index] += np.bincount(bins[good], weights=v[good, j],
                                                    minlength=n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        means = {k: np.where(counts > 0, sums[j] / counts, np.nan) for j, k in enumerate(keys)}
    return means


def profile_summary(per_replicate):
    """Average replicates per region/bin, then median and IQR across regions."""
    keys = set.intersection(*(set(m) for m in per_replicate))
    out = {}
    for k in keys:
        stack = np.stack([m[k] for m in per_replicate])
        with np.errstate(all="ignore"), _quiet():
            reg = np.nanmean(stack, axis=0)
            out[k] = (np.nanmedian(reg, axis=0), np.nanpercentile(reg, 25, axis=0),
                      np.nanpercentile(reg, 75, axis=0),
                      np.sum(np.isfinite(reg), axis=0))
    return out


class _quiet:
    def __enter__(self):
        import warnings
        self._w = warnings.catch_warnings()
        self._w.__enter__()
        warnings.simplefilter("ignore", RuntimeWarning)

    def __exit__(self, *exc):
        self._w.__exit__(*exc)


def metaplot(samples, regions_path, out_prefix, names_path=None, name_column=4,
             title=None, flank=2000, bin_size=200, body_bins=47, smooth=3,
             formats=("png", "pdf"), region_label="region"):
    """samples: ordered {label: [result files]} (replicates averaged)."""
    plt = _mpl()
    layout = BinLayout(flank, bin_size, body_bins)
    regions = read_regions(regions_path, names_path, name_column)
    log.info("metaplot: %d regions, %d bins", len(regions), layout.n_bins)

    summaries, rows = {}, []
    for label, files in samples.items():
        reps = []
        for f in files:
            log.info("  %s: %s", label, f)
            reps.append(region_bin_means(f, regions, layout))
        summaries[label] = profile_summary(reps)
        for k, (med, q25, q75, n) in summaries[label].items():
            for b in range(layout.n_bins):
                rows.append((label, k, b, med[b], q25[b], q75[b], int(n[b]), len(files)))
    table = pd.DataFrame(rows, columns=["sample", "quantity", "bin", "median", "q25", "q75",
                                        "n_regions", "n_replicates"])
    table.to_csv(f"{out_prefix}_profile.tsv", sep="\t", index=False, float_format="%.6g")

    x = np.arange(layout.n_bins) + 0.5
    nf, nb = layout.n_flank, layout.body_bins
    ticks = [0, nf, nf + nb, layout.n_bins]
    kb = f"{flank / 1000:g}kb"
    ticklabels = [f"-{kb}", "start", "end", f"+{kb}"]

    def style_x(ax):
        ax.set_xlim(0, layout.n_bins)
        ax.set_xticks(ticks, ticklabels)
        for t in (nf, nf + nb):
            ax.axvline(t, color=MUTED, lw=0.8, ls=(0, (3, 3)), zorder=0)

    # figure 1: one row per sample; entropy references | observed - null
    n = len(samples)
    fig, axes = plt.subplots(n, 2, figsize=(10, 2.6 * n + 0.6), squeeze=False,
                             gridspec_kw={"width_ratios": [1.6, 1]})
    has_delta = all("delta" in s for s in summaries.values())
    dlim = 0.0
    for label in samples:
        if has_delta:
            _, q25, q75, _ = summaries[label]["delta"]
            dlim = max(dlim, np.nanmax(np.abs(np.r_[q25, q75])))
    dlim = (dlim or 0.1) * 1.15
    for i, label in enumerate(samples):
        s = summaries[label]
        ax = axes[i, 0]
        end_items = []
        for key, (_c, name, color, ls) in PROFILE_SERIES.items():
            if key not in s:
                continue
            med, q25, q75, _ = s[key]
            if key in ("entropy", "null_mean", "bimodal"):
                ax.fill_between(x, _smooth(q25, smooth), _smooth(q75, smooth), color=color,
                                alpha=0.15, lw=0)
            y = _smooth(med, smooth)
            ax.plot(x, y, color=color, ls=ls, label=name)
            end_items.append((name, y[-1]))
        ax.set_ylim(0, 1)
        style_x(ax)
        if i == 0:
            _end_labels(ax, end_items)
        ax.set_ylabel("Entropy / methylation")
        n_reg = int(np.nanmax(s["entropy"][3])) if "entropy" in s else 0
        ax.set_title(f"{label}  (n = {len(samples[label])} replicate"
                     f"{'s' if len(samples[label]) != 1 else ''}, {n_reg} {region_label}s)",
                     loc="left")
        style_x(ax)
        ax2 = axes[i, 1]
        if has_delta:
            med, q25, q75, _ = s["delta"]
            ax2.axhline(0, color=INK_2, lw=0.9)
            ax2.fill_between(x, _smooth(q25, smooth), _smooth(q75, smooth), color=BLUE,
                             alpha=0.15, lw=0)
            ax2.plot(x, _smooth(med, smooth), color=BLUE)
            ax2.set_ylim(-dlim, dlim)
            ax2.set_ylabel("Observed − null")
            ax2.set_title("below 0 = more ordered than chance", loc="left", fontsize=8,
                          color=INK_2)
        else:
            ax2.set_axis_off()
        style_x(ax2)
    fig.tight_layout()
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=len(labels),
               bbox_to_anchor=(0.5, 1.0))
    if title:
        fig.suptitle(title, y=1.085, fontsize=11, color=INK)
    paths = _save(fig, f"{out_prefix}_metaplot", formats)
    plt.close(fig)

    # figure 2: observed - null of all samples on one axis
    if has_delta:
        fig, ax = plt.subplots(figsize=(6.5, 3.4))
        ax.axhline(0, color=INK_2, lw=0.9)
        ends = []
        for j, label in enumerate(samples):
            med, q25, q75, _ = summaries[label]["delta"]
            color = SERIES[j % len(SERIES)] if j < len(SERIES) else MUTED
            ax.fill_between(x, _smooth(q25, smooth), _smooth(q75, smooth), color=color,
                            alpha=0.10, lw=0)
            y = _smooth(med, smooth)
            ax.plot(x, y, color=color, label=label)
            ends.append((label, y[-1]))
        ax.set_ylim(-dlim, dlim)
        ax.set_ylabel("Observed − null entropy")
        style_x(ax)
        _end_labels(ax, ends)
        ax.legend(loc="upper left", ncol=min(4, n))
        ax.set_title(title or "Observed − null entropy", loc="left")
        fig.tight_layout()
        paths += _save(fig, f"{out_prefix}_delta_overlay", formats)
        plt.close(fig)
    return paths, table


# --------------------------------------------------------------------- QC
def qc(result_path, out_prefix, fraction=0.05, max_points=1_500_000, seed=1,
       chunksize=2_000_000, formats=("png",), q_threshold=0.05, title=None):
    """Four-panel QC figure and a one-line-per-metric summary TSV."""
    plt = _mpl()
    header = pd.read_csv(result_path, sep="\t", nrows=0).columns
    cols = _resolve_columns(header)
    want = ["coverage" if "coverage" in header else "assumed_coverage", "methylation",
            *[cols[k] for k in ("entropy", "null_mean", "bimodal", "minimum") if k in cols]]
    for extra in ("delta", "p_two", "q_two"):
        c = cols.get("delta") if extra == "delta" else extra
        if c and c in header:
            want.append(c)
    want = list(dict.fromkeys(want))
    rng = np.random.default_rng(seed)
    parts = []
    n_total = n_low = n_high = 0
    delta_col = cols.get("delta")
    for chunk in pd.read_csv(result_path, sep="\t", usecols=want, chunksize=chunksize):
        n_total += len(chunk)
        if "q_two" in chunk and delta_col:
            sig = chunk["q_two"] < q_threshold
            n_low += int((sig & (chunk[delta_col] < 0)).sum())
            n_high += int((sig & (chunk[delta_col] > 0)).sum())
        parts.append(chunk.sample(frac=fraction, random_state=int(rng.integers(1 << 31))))
    df = pd.concat(parts, ignore_index=True)
    if len(df) > max_points:
        df = df.sample(n=max_points, random_state=seed)
    cov = want[0]

    fig, axes = plt.subplots(2, 2, figsize=(10, 7.5))
    ax = axes[0, 0]
    hb = ax.hexbin(df["methylation"], df[cols["entropy"]], gridsize=60, bins="log",
                   cmap="Blues", mincnt=1, linewidths=0)
    fig.colorbar(hb, ax=ax, label="windows (log)", shrink=0.85)
    edges = np.linspace(0, 1, 26)
    mid = (edges[:-1] + edges[1:]) / 2
    which = np.clip(np.digitize(df["methylation"], edges) - 1, 0, len(mid) - 1)
    for key in ("entropy", "null_mean", "bimodal", "minimum"):
        if key not in cols:
            continue
        _c, name, color, ls = PROFILE_SERIES[key]
        med = df.groupby(which)[cols[key]].median().reindex(range(len(mid)))
        ax.plot(mid, med.to_numpy(), color=color, ls=ls, label=f"{name} (median)")
    ax.set(xlim=(0, 1), ylim=(0, 1), xlabel="Methylation", ylabel="Entropy",
           title="Entropy vs methylation")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=7.5)

    ax = axes[0, 1]
    if delta_col:
        d = df[delta_col].dropna()
        lim = np.nanpercentile(np.abs(d), 99.5) if len(d) else 1
        ax.hist(d, bins=np.linspace(-lim, lim, 81), color=BLUE, alpha=0.85)
        ax.axvline(0, color=INK_2, lw=0.9)
        ax.set(xlabel="Observed − null entropy", ylabel="windows",
               title="Distribution of observed − null")
        if n_total:
            ax.text(0.02, 0.95,
                    f"q < {q_threshold}:\n{100 * n_low / n_total:.1f}% below null\n"
                    f"{100 * n_high / n_total:.1f}% above null",
                    transform=ax.transAxes, va="top", ha="left", fontsize=8, color=INK_2,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.9, pad=2))
    else:
        ax.set_axis_off()

    ax = axes[1, 0]
    if delta_col:
        cov_edges = np.unique(np.percentile(df[cov], np.linspace(0, 100, 16)).round())
        cb = np.clip(np.digitize(df[cov], cov_edges[1:-1]), 0, len(cov_edges) - 2)
        g = df.groupby(cb)[delta_col]
        centers = df.groupby(cb)[cov].median()
        ax.fill_between(centers, g.quantile(0.25), g.quantile(0.75), color=BLUE, alpha=0.15,
                        lw=0)
        ax.plot(centers, g.median(), color=BLUE, marker="o", ms=3)
        ax.axhline(0, color=INK_2, lw=0.9)
        ax.set(xlabel="Reads per window", ylabel="Observed − null (median, IQR)",
               title="Depth dependence (flat = depth-matched)")
    else:
        ax.set_axis_off()

    ax = axes[1, 1]
    if "p_two" in df:
        ax.hist(df["p_two"].dropna(), bins=np.linspace(0, 1, 41), color=ORANGE, alpha=0.85)
        ax.set(xlabel="p (two-sided)", ylabel="windows",
               title="Two-sided p-values")
    else:
        ax.set_axis_off()
    fig.suptitle(title or str(result_path).split("/")[-1], x=0.01, ha="left", fontsize=11,
                 color=INK)
    fig.tight_layout()
    paths = _save(fig, f"{out_prefix}_qc", formats)
    plt.close(fig)
    summary = {"windows": n_total, "plotted": len(df),
               "frac_q_low": n_low / n_total if n_total else np.nan,
               "frac_q_high": n_high / n_total if n_total else np.nan,
               "median_coverage": float(df[cov].median()),
               "median_delta": float(df[delta_col].median()) if delta_col else np.nan}
    pd.Series(summary).to_csv(f"{out_prefix}_qc_summary.tsv", sep="\t", header=False)
    return paths, summary
