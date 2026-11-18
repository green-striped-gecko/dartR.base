# Review: gl.fdsim (dartR.base)

- Family mode: analysis
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: abfbae1 (`dev_luis`, contains `origin/dev` aa8b0f9)
- Custodian: Arthur Georges (STY5: changes go through this report)
- Datasets: `testset.gl`, `testset.gs`, `platypus.gl`
- Baseline: `tests/testthat/test-gl.fdsim.R` (snapshot captured
  pre-review, 20 expectations pass)
- Callers: `gl.fixed.diff(test = TRUE)` (once per population pair,
  allopatric only); dartr2shiny exposes `sympatric`, `reps` and `delta`
  as inputs.

## Verdicts

**Standards: Needs work**: the structure follows the house order and the
FBM path gives identical results, but inputs are checked in the wrong
order, progress lines are uncoloured, and the header misdescribes the
return value.

**Spec: Rework**: the p-value is too small by up to 18 orders of
magnitude in the cases that decide significance, and `sympatric = TRUE`
computes from mixed-up rows. Both are errors in the statistic, not in
its presentation.

## Findings

**F1 [HIGH, confidence: high] — p-value from a normal curve fitted to
the wrong spread (spec axis)**
`R/gl.fdsim.r:193-209`: each replicate computes, per locus, the
probability `fd[i]` that sampling alone produces a fixed difference, and
stores the sum. That sum is the *expected* count for the replicate, not
a count. `sdexpected` is therefore the spread of the expectation across
replicates, which leaves out the locus-by-locus chance of a false
positive occurring or not. The p-value then comes from
`pnorm(obs, mn, sdev, lower.tail = FALSE)`. A normal curve is a poor fit
for a count whose mean is near zero, and `lower.tail = FALSE` gives
P(X > obs), not P(X >= obs).
Failure scenario (verified, same model, 5000 simulated counts):

| Pair | obs | `gl.fdsim` p | Simulated P(count >= obs) |
|---|---|---|---|
| SEVERN_BELOW vs TENTERFIELD (`platypus.gl`) | 1 | 4.6e-20 | 0.041 |
| SEVERN_BELOW vs TENTERFIELD | 2 | 3.9e-77 | 0.0002 |
| EmsubRopeMata vs EmmacBurnBara (`testset.gl`) | 4 | 0.022 | 0.152 |
| EmsubRopeMata vs EmmacBurnBara | 5 | 0.0006 | 0.029 |

At `obs = 4` the function reports a significant result where false
positives alone reach 4 in 15% of simulations. For the testset pair the
standard deviation of the expectation is 0.80; that of the count is 1.11.
Every `gl.fixed.diff(test = TRUE)` p-value inherits this.
Proposed change: in each replicate, draw the count itself (one
Bernoulli draw per locus with probability `fd[i]`). Report the mean and
standard deviation of the simulated counts, and take the p-value as the
simulated tail, `(sum(count >= obs) + 1) / (reps + 1)`. The +1 stops a
finite simulation from returning p = 0.

**F2 [HIGH, confidence: high] — `sympatric = TRUE` reads mixed rows
(spec axis)**
`R/gl.fdsim.r:129-132`: for sympatric pairs, `rfA` and `rfB` are both set
to the full table from `gl.allele.freq(by = "popxloc")`. That table
alternates the two populations row by row (86 rows for 43 loci). The
loop reads row `i` for locus `i`, so locus 1 gets population 1 at locus 1,
locus 2 gets population 2 at locus 1, and so on. Only the first 43 rows
are used, and both "populations" get the same frequency and sample size.
Failure scenario (verified): `testset.gl`, EmsubRopeMata vs
EmmacBurnBara: expected false positives 0.0016 (allopatric: 2.33),
p = 0. Swapping the population order gives an identical result.
`gl.fixed.diff` points sympatric users to `gl.pval.sympatry()`, which
does not exist in any dartR package.
Proposed change: for sympatric pairs, sample both populations from their
pooled allele frequency (weighted by sample size), each with its own
sample size. This is the null hypothesis that both samples come from one
gene pool. The intended model needs confirmation from the custodian
(confidence in the fix: medium; confidence in the defect: high).

**F3 [MEDIUM, confidence: high] — SilicoDArT sample size doubled
(DAT1)**
`R/gl.fdsim.r:175-176`: the sample size is `nobs * 2` alleles regardless
of data type. For SilicoDArT (presence/absence, ploidy 1) the unit is the
individual, so the chance of a spurious fixed difference is computed for
twice the real sample and comes out too small.
Failure scenario (verified): `gl.fdsim(testset.gs, ...)` runs without a
message and returns expectations based on 2n.
Proposed change: use `nobs * 2` for SNP and `nobs` for SilicoDArT, keyed
on the data type that `utils.check.datatype()` already returns.

**F4 [LOW, confidence: high] — inputs checked in the wrong order and
incompletely (FS5)**
`R/gl.fdsim.r:74-99`:
- the labels are checked before the length, so
  `poppair = "EmsubRopeMata"` fails with "Population B mislabelled";
- `poppair = c(A, A)` is not rejected;
- `reps = 1` returns `sdexpected = NA` and `prob = NA` without a message;
- `obs` is not checked.
Proposed change: check that `poppair` is two distinct labels, then that
each exists, that `reps` is a whole number of at least 2, and that `obs`
is `NULL` or a single non-negative number, all before any work.

**F5 [LOW, confidence: high] — per-locus loop inside the replicate loop
(STY2)**
`R/gl.fdsim.r:157-200`: the double loop runs `reps × nLoc` iterations in
R. `platypus.gl` (1000 loci) takes 0.38 s for 100 replicates, so about
4 s at the default 1000. `gl.fixed.diff(test = TRUE)` calls the function
once per pair of populations: 190 calls for 20 populations. The
exclusion test and the sample sizes do not change between replicates.
Proposed change: compute the locus filter once and draw all loci per
replicate with vectorised `rbinom()`. This changes the random number
sequence, so seeded results change even where the model does not. It is
worth bundling with change 1, which changes the numbers anyway.

**F6 [LOW, confidence: high] — header misdescribes the function
(DOC1, DOC5 (proposed rule))**
`R/gl.fdsim.r:1-50`:
- `@return` promises "square matrices" `[[1]]`–`[[4]]`; the function
  returns a named list of four numbers (`observed`, `mnexpected`,
  `sdexpected`, `prob`);
- the last sentence of `@description` is garbled ("The probability of the
  observed count of fixed differences is greater than the expected number
  of false positives is calculated");
- no `@details` on the model or on what `sympatric` changes;
- the only example uses `sympatric = TRUE`, the broken path (F2);
- `@family distance` groups it with distance functions; `gl.fixed.diff`
  is its closest relative;
- `@examples` sits before `@return`, out of house order.
Proposed change: rewrite the header. Docs only apart from the example.

**F7 [LOW, confidence: high] — progress lines uncoloured and sample
sizes unlabelled (VRB2)**
`R/gl.fdsim.r:103-118, 212-225`: these use bare `cat()`. At
`verbose = 3` the line "Sample sizes: 11 5" prints in factor-level
order, so for `poppair = c("EmsubRopeMata", "EmmacBurnBara")` the 11 is
EmmacBurnBara's, the reverse of the order the user gave.
Proposed change: wrap them in `report()`, and print each sample size
next to its population name.

## Proposed changes

1. Simulate the false-positive count itself and take the p-value as the
   simulated tail, `(sum(count >= obs) + 1) / (reps + 1)`; report the
   mean and standard deviation of the simulated counts (F1).
   **Consequence: `sdexpected` and `prob` change for every call, and so do
   `gl.fixed.diff(test = TRUE)` p-values. p-values become larger (less
   significant), often by many orders of magnitude when few false
   positives are expected. The smallest reportable p becomes
   1 / (reps + 1).**
2. Sympatric pairs sample both populations from the pooled allele
   frequency, each with its own sample size (F2). Also fix the
   `gl.fixed.diff` help reference to `gl.pval.sympatry()`, which does
   not exist.
   **Consequence: `sympatric = TRUE` results change completely; they
   become symmetric in population order.**
3. Sample size is `nobs` for SilicoDArT, `2 * nobs` for SNP (F3).
   **Consequence: SilicoDArT results change; expected false positives
   rise.**
4. Check `poppair` (two distinct labels, both present), `reps` (whole
   number >= 2) and `obs` (`NULL` or one non-negative number) before any
   work (F4).
   **Consequence: `reps = 1` and `poppair = c(A, A)` now stop with an
   error instead of returning `NA` or running.**
5. Vectorise the per-locus draws (F5).
   **Consequence: results with a fixed seed differ from earlier versions,
   even where the model is the same.**
6. Rewrite the roxygen header: `@return`, `@description`, `@details`,
   allopatric example, `@family`, tag order (F6).
7. Colour the progress lines and label the sample sizes (F7).

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run. PLT not
  applicable (no plot). DEP: only `stats`, imported.
- Spec: behaviour against roxygen on `testset.gl`, `testset.gs`,
  `platypus.gl` — run.
- Numerical check against an independent computation: the per-locus
  model re-implemented outside the function; its expectation (2.337) and
  spread (0.796) match `gl.fdsim` (2.33, 0.81), and its simulated counts
  give the right-hand column of the F1 table — run.
- FBM path (DAT6): `gl.gen2fbm(testset.gl)` gives identical output — run.
- Callers: `gl.fixed.diff` (dartR.base), dartr2shiny config and
  templates — grepped. No other dartR package calls it.
- Method source: the Georges et al. (2018) paper that describes the
  approach was not re-read (not available: no access in this session).
  Changes 1 and 2 change the published method, so the custodian should
  confirm them.
- Google Group / GitHub issues: not searched (not available: no browser
  session).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis |  |
| 2 | approved | Luis | pooled-frequency model |
| 3 | approved | Luis |  |
| 4 | approved | Luis |  |
| 5 | approved | Luis |  |
| 6 | approved | Luis |  |
| 7 | approved | Luis |  |

## Outcome

- Change 1 applied: each replicate draws the count; p = simulated tail.
  `testset.gl` obs 4 / 5 / 6: p = 0.147 / 0.029 / 0.0024 (independent
  re-implementation: 0.152 / 0.029 / 0.0022). `platypus.gl`
  SEVERN_BELOW vs TENTERFIELD obs 1: p = 0.041 (independent: 0.041).
  Testset `sdexpected` 1.20 (was 0.81).
- Change 2 applied: pooled, sample-size-weighted frequency; testset
  sympatric expectation 0.005 in both population orders (within
  simulation error). `gl.fixed.diff` help now points to
  `gl.fdsim(sympatric = TRUE)` instead of the nonexistent
  `gl.pval.sympatry()`.
- Change 3 applied: `testset.gs` expected false positives 6.6 (was 2.6
  with 2n).
- Change 4 applied: five bad-argument cases stop with specific messages.
- Change 5 applied: `platypus.gl`, 1000 replicates, 0.14 s (was about
  3.8 s).
- Change 6 applied: header rewritten; `devtools::document()` run;
  `man/gl.fdsim.Rd` and `man/gl.fixed.diff.Rd` regenerated. Example runs.
- Change 7 applied: "Sample sizes: EmsubRopeMata = 5, EmmacBurnBara = 11".
- Snapshot: 23 expectations pass; every diff from the baseline is marked
  `[approved diff change N]`. FBM output is still identical to in-memory
  output.
- Full suite: the same 9 files fail as before this change, plus
  `test-as-dartR.R:40`. That assertion expects `platypus.gl` to lack the
  `@fbm` slot, and the dartR.data build installed at 14:29 today has
  it. It is unrelated to this change.
- NEWS entry added (escalation gate: numerical output changes).
- PR: pending.

```json
{
  "function": "gl.fdsim",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "abfbae1",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DAT1", "status": "approved", "change": 3},
    {"id": "F4", "severity": "LOW", "confidence": "high", "rule": "FS5", "status": "approved", "change": 4},
    {"id": "F5", "severity": "LOW", "confidence": "high", "rule": "STY2", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "DOC1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "VRB2", "status": "approved", "change": 7}
  ],
  "coverage_skipped": ["PLT: no plot", "method paper not re-read", "Google Group: no browser session"],
  "status": "pr-open",
  "pr": null
}
```
