# Changelog

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
