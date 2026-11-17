# Interpreting methentropy output

This page explains every number `methentropy run` reports for a 4-CpG window,
using one real window as an example, and shows where each number appears in
the figures made by `qc`, `metaplot` and `compare`. For the statistical
background see [methods.md](methods.md).

## One window, worked through

A window covered by **20 reads** over 4 CpGs:

| Pattern | Reads |
|---|---|
| `GGGG` (all four CpGs methylated) | 6 |
| `GGGg` (CpG 4 unmethylated) | 2 |
| `gggg` (all unmethylated) | 12 |

CpGs 1-4 are methylated in 8, 8, 8 and 6 reads: 30 methylated calls out of
80, so methylation = 0.375. Running this window through methentropy gives:

| Column | Value | Meaning |
|---|---|---|
| `coverage` | 20 | reads spanning all four CpGs (N) |
| `methylation` | 0.375 | mean methylation; `meth_cpg1..4` give each CpG (0.40, 0.40, 0.40, 0.30) |
| `entropy` (H) | **0.324** | observed entropy of the 16 epiallele frequencies: (1/4) Σ −p log2 p |
| `entropy_min` | 0.203 | the lowest entropy possible for 20 reads carrying 30 methylated calls (exact; see below) |
| `entropy_bimodal` | 0.243 | entropy if every read were fully methylated or fully unmethylated at this level |
| `null_mean` | 0.828 | mean entropy of 1,000 simulated windows with CpGs methylated independently |
| `null_sd` | 0.039 | spread of those 1,000 simulated entropies |
| `null_q05`, `null_q95` | 0.760, 0.887 | 5th and 95th percentiles of the null |
| `delta` (Δ) | **−0.504** | H − null_mean; negative = more ordered than chance |
| `z` | −12.8 | Δ / null_sd: distance from chance in null standard deviations |
| `scaled_entropy` | 0.194 | (H − min) / (null_mean − min): 0 = minimum, 1 = chance |
| `p_low` | 0.001 | P(null ≤ H): small when the window is more ordered than chance |
| `p_high` | 1.0 | P(null ≥ H): small when the window is more disordered than chance |
| `p_two` | 0.002 | two-sided p-value, 2 × min(p_low, p_high), capped at 1 |
| `q_two` | | `p_two` adjusted for all windows in the file (Benjamini-Hochberg) |

Most reads here are all-or-none, so CpGs on the same molecule agree with each
other far more than independent CpGs would: observed entropy is well below the
null.

## The references

### The null and `null_mean`

The **null** is not one number but a distribution. For each window, the reads
are shuffled independently at every CpG, 1,000 times by default. Each CpG
keeps exactly its observed number of methylated reads (here 8, 8, 8 and 6 of
20), so depth and per-CpG methylation are unchanged; only the coupling
between CpGs on the same read is destroyed. The entropy of every shuffled
window is recorded:

* `null_mean`: the average of the 1,000 values, i.e. the entropy expected by
  chance for this depth and these CpG levels;
* `null_sd`: their standard deviation;
* `null_q05`, `null_q95`: their 5th and 95th percentiles.

`null_q05` and `null_q95` correspond to "stochastic min" and "stochastic max"
of the original v1 pipeline. They are the edges of chance, not the theoretical
limits of entropy.

### `entropy_min`

The lowest entropy any arrangement of the observed methylated calls could
give. The calls are distributed over reads that carry 0, 1, 2, 3 or 4
methylated CpGs, all reads in one class sharing one pattern, and the
allocation with the lowest entropy is found exactly. It is 0 when every read
can carry the same pattern (25, 50 or 75 % methylation with N reads), and it
is a mathematical floor rather than a biological state. In the example it is
0.203 (15 reads with two methylated CpGs and 5 unmethylated reads).

Do not confuse it with `null_q05`: in the example the minimum is 0.203 while
the lower edge of chance is 0.760.

### `entropy_bimodal`

The entropy if each read were either fully methylated or fully unmethylated,
in proportion to the methylation level: h(m) / 4, evaluated at the window's
depth (0.25 at 50 % methylation). This is the expectation for two
homogeneous populations of molecules, for example allele-specific
methylation or a mix of two cell types.

## The statistics

* **`delta`** is the main effect size. It is on the entropy scale (0 to 1),
  so it is directly comparable with the curves in the metaplot.
* **`scaled_entropy`** places H between the minimum and chance. Because both
  references move with methylation, it is easier to compare across windows
  or conditions with different methylation levels. It is left empty when
  `null_mean` and `entropy_min` are nearly equal (e.g. fully methylated
  windows).
* **`z`** expresses Δ in units of the null's own spread. It is useful for
  ranking windows; it is not used in the figures.
* **`p_low` / `p_high`** are Monte Carlo p-values with the +1 correction:
  (number of simulated values at least as extreme + 1) / (simulations + 1).
  With 1,000 simulations the smallest possible value is 1/1001 ≈ 0.001; use
  `--nsims 5000` or more if strict genome-wide thresholds are needed.
* **`p_two`** asks whether the window departs from chance in either
  direction.
* **`q_two`** (Benjamini-Hochberg): with millions of windows, about 5 % would
  reach p < 0.05 by chance alone. BH turns each p-value into a q-value such
  that among all windows with q < 0.05, about 5 % are expected to be false
  discoveries. Use q to call individual windows.

Two patterns to expect in real data:

* Windows where every CpG is fully methylated or fully unmethylated have a
  null with no spread (every shuffle gives the same entropy), so p = 1. They
  form the bar at 1 in the p-value histogram of `qc`.
* Region averages of Δ are much smaller than this example. Most windows sit
  near the null and a minority are strongly ordered; a region average of
  −0.03 usually means more ordered windows, not every window shifting a
  little. The fraction of windows with q < 0.05 (`qc`, `aggregate`) and the
  per-region counts of `compare` show this directly.

## Where each number appears

| Column | `qc` | `metaplot` | `compare` | `aggregate` |
|---|---|---|---|---|
| `entropy` | panel A (hexbins and median) | blue line, left | in scaled entropy | `mean_entropy` |
| `null_mean` | panel A, orange | orange line, left | in scaled entropy | `mean_null` |
| `entropy_min` | panel A, grey dashed | grey dashed, left | in scaled entropy | |
| `entropy_bimodal` | panel A, green | green line, left | | |
| `methylation` | panel A, x axis | black dotted, left | panel C, x axis | `mean_methylation` |
| `delta` | panels B and C | right panel and overlay figure | panel A, panel C y axis | `mean_delta`, `weighted_delta` |
| `scaled_entropy` | | | panel B (region-level, see below) | `median_scaled_entropy` |
| `p_two` | panel D | | | |
| `q_two` | panel B text (% windows below / above null) | | | `frac_low_sig`, `frac_high_sig` |
| `z`, `null_sd`, `null_q05/q95`, `p_low`, `p_high` | | | | |

Columns in the last row are in the result TSV only.

### Reading the metaplot

Each line is the median across regions of each region's per-bin average, and
the shading is the 25-75 % range **across regions**. The shading shows how
different regions are from one another; it is not the spread of the null.
Observed (blue) can therefore lie inside the orange band and still be
consistently below the orange line. Whether observed differs from the null is
read from the right-hand panel (Δ), which is computed window by window.

### Reading `compare`

`compare` works on regions, not windows. For each region and replicate it
pools the windows inside the region:

* region `delta` = mean Δ of its windows;
* region `scaled` = (Σ H − Σ min) / (Σ null_mean − Σ min) over its windows,
  a ratio of sums so that near-zero-variance windows do not dominate.

Each condition is compared with the first (reference):

| Summary column | Test | Unit |
|---|---|---|
| `median_diff`, `frac_regions_lower` | per-region difference (other − reference) | region |
| `wilcoxon_p_regions` | paired Wilcoxon signed-rank test | region |
| `ref_replicate_means`, `other_replicate_means`, `welch_p_replicates` | Welch t-test on per-replicate means over regions | replicate |
| `n_regions_q05_lower/higher` | per-region Welch test across replicates, BH over regions | region |

`welch_p_replicates` is the conservative test (n = number of replicates or
samples) and the one to report for a claim about conditions; the region-level
Wilcoxon p-value is very small whenever many regions move in the same
direction. Per-region tests have little power with 2-3 replicates, so few
regions reach q < 0.05 even when the overall shift is clear.

In the figure, boxes show per-region values (replicates averaged), dots show
replicate means, and the brackets give both p-values. Panel C plots each
region's change in Δ against its change in methylation, to show whether
regions that gain more methylation also become more ordered.
