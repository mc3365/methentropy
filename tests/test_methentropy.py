import gzip
import os
import json

import numpy as np
import pandas as pd
import pytest

from methentropy import cli
from methentropy.core import (BITS, EPIALLELE_NAMES, _class_entropy, bimodal_entropy,
                              binary_entropy, entropy_from_counts, minimum_entropy,
                              per_cpg_methylated)
from methentropy.nulls import (empirical_pvalues, independent_key, pooled_key,
                               simulate_independent, simulate_pooled)
from methentropy.pipeline import bh_qvalues, run

WEIGHTS = np.array([8, 4, 2, 1])


# ---------------------------------------------------------------- helpers
def synth_window(rng, n, p, linked=False):
    """Epiallele counts for n reads. linked=True: each read all-or-none."""
    if linked:
        state = rng.random(n) < np.mean(p)
        reads = np.repeat(state[:, None], 4, axis=1)
    else:
        reads = rng.random((n, 4)) < np.asarray(p)
    return np.bincount(reads.astype(int) @ WEIGHTS, minlength=16)


def write_rlm(path, windows, compress=False):
    rows = []
    for i, c in enumerate(windows):
        n = int(c.sum())
        meth = per_cpg_methylated(c).sum() / (4 * n)
        rows.append(["chr1", 100 * i, 100 * i + 2, entropy_from_counts(c), 0.0,
                     *c.tolist(), meth, n])
    df = pd.DataFrame(rows)
    opener = gzip.open if compress else open
    with opener(path, "wt") as h:
        h.write("#chr\tstart\tend\tentropy\tepipolymorphism\t" + "\t".join(EPIALLELE_NAMES)
                + "\tmean_methylation\tcoverage\n")
        df.to_csv(h, sep="\t", header=False, index=False)


# ---------------------------------------------------------------- core
def test_epiallele_convention_matches_rlm():
    assert EPIALLELE_NAMES[0] == "gggg" and EPIALLELE_NAMES[1] == "gggG"
    assert EPIALLELE_NAMES[8] == "Gggg" and EPIALLELE_NAMES[15] == "GGGG"
    assert BITS[1].tolist() == [0, 0, 0, 1]


def test_entropy_known_values():
    c = np.zeros(16); c[0] = 10
    assert entropy_from_counts(c) == 0
    c[15] = 10
    assert entropy_from_counts(c) == pytest.approx(0.25)
    assert entropy_from_counts(np.ones(16)) == pytest.approx(1.0)


def _brute_min(n, k):
    best = np.inf
    for c4 in range(n + 1):
        for c3 in range(n + 1 - c4):
            for c2 in range(n + 1 - c4 - c3):
                c1 = k - 4 * c4 - 3 * c3 - 2 * c2
                if 0 <= c1 <= n - c4 - c3 - c2:
                    c0 = n - c4 - c3 - c2 - c1
                    best = min(best, _class_entropy((c0, c1, c2, c3, c4), n, 4))
    return best


@pytest.mark.parametrize("n", [1, 2, 3, 5, 10, 16, 23])
def test_minimum_entropy_exact(n):
    for k in range(4 * n + 1):
        assert minimum_entropy(n, k) == pytest.approx(_brute_min(n, k), abs=1e-12)


def test_minimum_entropy_reproduces_v1_lookup():
    # three entries from the v1 16-read minimum_entropy_lookup.json
    assert minimum_entropy(16, 1) == pytest.approx(0.08432251665425347)
    assert minimum_entropy(16, 16) == 0.0
    assert minimum_entropy(16, 24) == pytest.approx(0.2028195311147832)


def test_bimodal_is_quarter_binary_entropy():
    assert bimodal_entropy(100, 200) == pytest.approx(binary_entropy(0.5) / 4)
    assert bimodal_entropy(10, 0) == 0


# ---------------------------------------------------------------- nulls
def test_keys_are_symmetric():
    assert independent_key(20, [3, 17, 10, 0]) == independent_key(20, [0, 10, 3, 3])
    assert pooled_key(10, 5) == pooled_key(10, 35)


def test_independent_null_preserves_marginals_and_is_reproducible():
    a = simulate_independent((20, 2, 5, 9, 10), 200, seed=3)
    b = simulate_independent((20, 2, 5, 9, 10), 200, seed=3)
    assert np.array_equal(a, b)
    assert np.all((a >= 0) & (a <= np.log2(20) / 4 + 1e-12))


def test_independent_null_matches_large_n_analytic():
    p = np.array([0.9, 0.6, 0.4, 0.1])
    n = 2000
    ks = np.rint(p * n).astype(int)
    sims = simulate_independent(independent_key(n, ks), 50, seed=1)
    assert sims.mean() == pytest.approx(binary_entropy(p).sum() / 4, abs=0.01)


def test_pooled_null_depth_dependence():
    low = simulate_pooled(pooled_key(10, 20), 500, seed=1).mean()
    high = simulate_pooled(pooled_key(128, 256), 500, seed=1).mean()
    assert low < high - 0.2  # finite-depth bias the v1/v2 pipeline ignored


def test_pvalues():
    sims = np.arange(99, dtype=float)
    lo, hi = empirical_pvalues(sims, np.array([-1.0, 98.0]))
    assert lo[0] == pytest.approx(1 / 100) and hi[1] == pytest.approx(2 / 100)


def test_bh():
    q = bh_qvalues(np.array([0.01, 0.04, 0.03, np.nan, 0.5]))
    assert np.isnan(q[3])
    assert q[0] == pytest.approx(0.04) and q[4] == pytest.approx(0.5)


# ---------------------------------------------------------------- pipeline
@pytest.fixture
def rlm_file(tmp_path):
    rng = np.random.default_rng(0)
    windows = [synth_window(rng, n, p) for n, p in
               [(12, [.5] * 4), (30, [.9, .6, .4, .1]), (25, [.2] * 4)]]
    windows += [synth_window(rng, 30, [.5] * 4, linked=True) for _ in range(5)]
    windows += [synth_window(rng, 5, [.5] * 4)]  # below min coverage
    path = tmp_path / "sample_entropy.bed.gz"
    write_rlm(path, windows, compress=True)
    return path


@pytest.mark.parametrize("workers", [1, 2])
def test_run_end_to_end(rlm_file, tmp_path, workers):
    out = tmp_path / f"out{workers}.tsv.gz"
    stats = run(rlm_file, out, n_sims=300, seed=1, workers=workers, chunksize=4)
    df = pd.read_csv(out, sep="\t")
    assert len(df) == 8 and stats["low_coverage"] == 1
    assert df["q_two"].notna().all()
    linked = df.iloc[3:]
    assert (linked["p_low"] < 0.01).all()      # all-or-none reads: strong linkage
    assert (linked["delta"] < 0).all()
    assert (df["entropy"] >= df["entropy_min"] - 1e-12).all()


def test_results_independent_of_workers_and_chunks(rlm_file, tmp_path):
    run(rlm_file, tmp_path / "a.tsv", n_sims=200, seed=5, workers=1, chunksize=100)
    run(rlm_file, tmp_path / "b.tsv", n_sims=200, seed=5, workers=2, chunksize=3)
    a, b = pd.read_csv(tmp_path / "a.tsv", sep="\t"), pd.read_csv(tmp_path / "b.tsv", sep="\t")
    pd.testing.assert_frame_equal(a, b)


def test_legacy_requires_coverage(tmp_path):
    path = tmp_path / "legacy.bed"
    path.write_text("chr1\t1\t3\t0.5\t0.5\nchr1\t5\t7\t0.1\t0.1\n")
    with pytest.raises(ValueError, match="legacy"):
        run(path, tmp_path / "o.tsv")
    run(path, tmp_path / "o.tsv", n_sims=100, legacy_coverage=128)
    df = pd.read_csv(tmp_path / "o.tsv", sep="\t")
    assert len(df) == 2 and "pooled_mean" in df


def test_cli_aggregate_and_lookup(rlm_file, tmp_path):
    out = tmp_path / "w.tsv"
    assert cli.main(["run", "-i", str(rlm_file), "-o", str(out), "--nsims", "100"]) == 0
    regions = tmp_path / "r.bed"
    regions.write_text("chr1\t0\t250\tfirst\nchr1\t250\t1000\tsecond\nchr2\t0\t10\tnone\n")
    agg_path = tmp_path / "agg.tsv"
    cli.main(["aggregate", "-i", str(out), "-r", str(regions), "-o", str(agg_path)])
    agg = pd.read_csv(agg_path, sep="\t")
    assert agg["n_windows"].tolist() == [3, 5, 0]
    cli.main(["lookup", "--nreads", "16", "128", "-o", str(tmp_path / "lk.tsv")])
    lk = pd.read_csv(tmp_path / "lk.tsv", sep="\t")
    assert len(lk) == 65 + 513


def test_calibration_is_controlled(tmp_path):
    rng = np.random.default_rng(11)
    windows = [synth_window(rng, int(n), rng.uniform(0.05, 0.95, 4))
               for n in rng.integers(10, 40, 300)]
    path = tmp_path / "cal.bed"
    write_rlm(path, windows)
    from methentropy.calibrate import calibrate
    summary, _ = calibrate(path, n_windows=300, n_sims=400, seed=2)
    assert summary["frac_p_two_lt_0.05 (target <= 0.05)"] <= 0.08
    # the v1/v2 rule flags far more than 5% of linkage-free windows
    assert summary["old_method_frac_below_band_128reads (target ~0.05)"] > 0.5


# ---------------------------------------------------------------- plotting
def test_bin_layout_matches_enrichedheatmap_defaults():
    from methentropy.plotting import BinLayout
    lay = BinLayout()
    assert (lay.n_flank, lay.n_bins) == (10, 67)
    pos = np.array([-2001, -2000, -1, 0, 999, 1000, 1999, 2000, 3999, 4000])
    b = lay.assign(pos, start=0, end=2000)
    assert b.tolist() == [-1, 0, 9, 10, 33, 33, 56, 57, 66, -1]


def test_metaplot_and_qc_run(rlm_file, tmp_path):
    pytest.importorskip("matplotlib")
    from methentropy.plotting import metaplot, qc
    out = tmp_path / "w.tsv.gz"
    run(rlm_file, out, n_sims=100)
    regions = tmp_path / "r.bed"
    regions.write_text("chr1\t100\t500\tA\textra\nchr1\t300\t700\tB\textra\n")
    names = tmp_path / "keep.txt"
    names.write_text("A\n")
    paths, table = metaplot({"s1": [str(out)], "s2": [str(out), str(out)]}, regions,
                            str(tmp_path / "mp"), names_path=names, body_bins=4,
                            flank=200, bin_size=100, formats=("png",))
    assert len(paths) == 2 and all(os.path.exists(p) for p in paths)
    assert set(table["quantity"]) >= {"entropy", "null_mean", "delta", "minimum", "bimodal"}
    qpaths, summary = qc(out, str(tmp_path / "q"), fraction=1.0)
    assert os.path.exists(qpaths[0]) and summary["windows"] == 8


def test_compare_detects_gain_in_linkage(tmp_path):
    pytest.importorskip("scipy")
    pytest.importorskip("matplotlib")
    from methentropy.compare import compare
    rng = np.random.default_rng(3)
    regions, results = [], {"A": [], "B": []}
    starts = [1000 * i for i in range(30)]
    for i, s in enumerate(starts):
        regions.append(f"chr1\t{s}\t{s + 600}\tr{i}")
    for cond, linked in (("A", False), ("B", True)):
        for rep in range(3):
            wins = []
            for s in starts:
                for k in range(6):
                    p = rng.uniform(0.3, 0.7)
                    wins.append(synth_window(rng, 30, [p] * 4, linked=linked and k % 2 == 0))
            path = tmp_path / f"{cond}{rep}.bed"
            write_rlm(path, wins)
            df = pd.read_csv(path, sep="\t", comment="#", header=None)
            df[1] = [s + 100 * k for s in starts for k in range(6)]
            df[2] = df[1] + 2
            with open(path, "w") as h:
                h.write("#chr\n")
                df.to_csv(h, sep="\t", header=False, index=False)
            out = tmp_path / f"{cond}{rep}.tsv"
            run(path, out, n_sims=200, seed=rep)
            results[cond].append(str(out))
    bed = tmp_path / "r.bed"
    bed.write_text("\n".join(regions) + "\n")
    summary, table, paths = compare(results, bed, str(tmp_path / "cmp"), formats=("png",))
    d = summary[summary.metric == "delta"].iloc[0]
    assert d.n_regions == 30 and d.median_diff < 0 and d.wilcoxon_p_regions < 1e-4
    assert d.welch_p_replicates < 0.05
    assert os.path.exists(paths[0])
