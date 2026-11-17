# Changelog

## 3.2.1
- `docs/interpreting_output.md`: every output column explained on a worked
  window (null, null_mean, minimum, bimodal, delta, scaled entropy, z,
  p_low/p_high/p_two, BH q) and where each appears in the qc, metaplot,
  compare and aggregate outputs; the example is checked by the test-suite.
- Workflow schematic: minimum example now shows identical reads (entropy 0).

## 3.2.0
- `methentropy compare`: per-region means for every replicate (observed - null,
  scaled entropy, methylation); conditions compared against the first with a
  paired Wilcoxon test over regions, a Welch test over replicate means and
  per-region Welch tests with BH; summary/region tables and a figure.
- Workflow schematic (`docs/figures/methentropy_workflow.{svg,png}`, script
  `make_workflow.py`) shown at the top of the README.
- Metaplot profile table: quantity `null` renamed to `null_mean` (the old name
  is read as missing by pandas and some R readers).
- `scipy` added to the `plot` extra.

## 3.1.0
- `methentropy metaplot`: median (IQR) profiles of observed, null, bimodal and
  minimum entropy plus methylation over scaled regions with flanks, the design of
  the v1 `example_mef_metaplots.Rmd` (EnrichedHeatmap extend 2000 / w 200 /
  target_ratio 0.7 = 10 + 47 + 10 bins); replicates averaged per region; second
  panel and overlay figure for observed - null; profile table written as TSV.
- `methentropy qc`: entropy vs methylation with references, observed - null
  distribution, depth dependence and p-value histogram.
- Optional dependency group `plot` (matplotlib).

## 3.0.0
Rewrite of `rlm_simulation` (v1) and `methylation-entropy-simulation` (v2).

- Reads the full RLM entropy output (epiallele counts, coverage).
- References computed at each window's own depth.
- New independence (per-CpG permutation) null; pooled v1/v2 null optional.
- Exact minimum entropy per (N, k); separate bimodal reference.
- delta, z, scaled entropy, Monte Carlo p-values, BH q-values.
- Deterministic per-key seeding; streaming chunks; multiprocessing.
- `calibrate`, `aggregate` and `lookup` subcommands; test-suite.

Fixed relative to v1: `--nreads` was ignored; 16-read minimum combined with
128-read null; chunking failed for small inputs; no seed.
