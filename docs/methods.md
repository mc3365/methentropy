# Methods

For a column-by-column guide with a worked window, see
[interpreting_output.md](interpreting_output.md).

## Quantity

For a 4-mer of consecutive CpGs spanned by N complete reads with epiallele
counts n_1..n_16, RLM reports methylation entropy

    H = (1/4) * sum_i -(n_i/N) log2(n_i/N)

which lies between 0 (all reads identical) and 1 (all 16 patterns equally
frequent). H depends on three things besides biology:

1. **Mean methylation.** Fewer patterns are likely near 0% or 100%.
2. **Per-CpG methylation levels.** For the same mean, unequal levels at the
   four CpGs (e.g. 0.9/0.6/0.4/0.1) give lower entropy than equal levels.
3. **Read depth.** The plug-in estimator is capped at log2(N)/4 and biased
   downward for small N (expected entropy at 50% methylation with independent
   CpGs: 0.71 at N=10, 0.85 at N=20, 0.93 at N=40, 0.98 at N=128).

A biological claim about heterogeneity has to control for all three. v1/v2
controlled (1) only, at a fixed depth of 128 for the null and 16 for the minimum.
Because typical WGBS windows have 10–40 reads, most windows fall below a 128-read
null for reasons of depth alone.

## References computed per window

**Exact minimum (`entropy_min`).** With k methylated calls in total, entropy is
minimal when all reads carrying j methylated CpGs share one pattern. Choosing
class sizes c_0..c_4 with sum c_j = N and sum j*c_j = k, we minimise the entropy
of the class distribution. Allocations over at most three classes are searched
exhaustively; the test-suite checks equality with full enumeration for N ≤ 23.
This is a mathematical bound, not a biological state: it is 0 at 25/50/75%
because every read can carry the same partially methylated pattern.

**Bimodal (`entropy_bimodal`).** Every read is fully methylated or fully
unmethylated, with the fraction of GGGG reads equal to the methylation level
(rounded at depth N): h(m)/4, e.g. 0.25 at 50%. This is the expectation for
allele-specific methylation or a mixture of two homogeneous cell populations.

**Independence null (`null_*`, p-values).** A conditional permutation test.
Each CpG keeps exactly its observed number of methylated reads; the read
labels are shuffled independently at each CpG, 1,000 times by default. This
destroys linkage between CpGs on a molecule while keeping depth and per-CpG
levels fixed. Under this null the observed entropy is exchangeable with the
simulated values, so the Monte Carlo p-value with +1 correction is valid
(conservative for discrete data). Because entropy is invariant to reordering
the four CpGs and to swapping methylated/unmethylated at any CpG, the null
depends only on (N, sorted min(k_j, N - k_j)), which makes caching effective.

For large N the null mean approaches sum_j h(p_j)/4, the entropy of a product
of independent Bernoullis. `delta` then estimates minus the multi-information
(total correlation) among the four CpGs: the amount of entropy removed by
co-methylation.

**Pooled null (`--with-pooled`, legacy).** The v1/v2 model: all 4N calls
shuffled together. It assumes the four CpGs share one level, so it confounds
per-CpG differences with linkage. Kept for comparison, evaluated at depth N.

## Interpretation

* `delta < 0`, small `p_low`: CpGs are co-methylated on molecules beyond what
  their levels imply. Candidate causes: allele-specific methylation, cell-type
  or cell-state mixtures, processive (de)methylation.
* `delta ≈ 0`: heterogeneity is explained by levels and depth.
* `delta > 0`, small `p_high`: more pattern diversity than independent CpGs
  produce. Rare; check for technical problems (mapping, conversion) first.
* `scaled_entropy` places the window between the minimum (0) and chance (1)
  and is comparable across methylation levels and depths.
* For comparisons between samples (e.g. passages), compare `delta` or
  `scaled_entropy`, not raw entropy, and aggregate over many windows or
  features; ideally test with replicates or permutation across samples.

## Calibration

`methentropy calibrate` samples real windows, draws synthetic reads with the
same depth and per-CpG levels but independent CpGs, and reports how often each
method calls significance. Expected: ≤ 5% at p < 0.05 for v3. For the v1/v2
rule (observed below the 5th percentile of the 128-read pooled null), the
fraction is typically ≫ 5%, which quantifies the bias of the earlier pipeline.

## Known limitations

* Bisulfite non-conversion and sequencing errors add random discordance,
  which raises entropy towards the null; true linkage is underestimated.
* Windows are defined by CpG count, not distance. Linkage decays with
  distance, so stratify by CpG density (e.g. CGI vs non-CGI) when comparing.
* RLM requires reads to span all four CpGs; long 4-mers (sparse CpGs) are
  under-represented for short reads.
* Monte Carlo p-values have a floor of 1/(nsims+1); increase `--nsims` if
  genome-wide FDR control at strict thresholds is needed.
* Model-based alternatives (Ising models: informME for WGBS, CpelNano for
  nanopore) estimate per-CpG levels and couplings jointly and can be used as
  orthogonal validation.
