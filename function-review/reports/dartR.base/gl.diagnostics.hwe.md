# Review: gl.diagnostics.hwe + gl.hwe.pop (dartR.base)

Reviewed together as the remaining Hardy-Weinberg functions. The custodian
asked specifically for `gl.diagnostics.hwe` to be made faster. Each
finding names the file it applies to.

## Provenance

- Family mode: analysis (both)
- Date: 2026-09-24
- Reviewer: Claude (Opus 5.5, claude-opus-5-5, Claude Code),
  dartr-function-review v2.0.0
- Package commit: f9f1be8 (origin/dev), worktree branch
  `review-hwe-diagnostics`
- Datasets: bandicoot.gl (96 individuals, 1000 loci, 5 populations),
  bandicoot.gl[, 1:100]
- Baseline: `tests/testthat/test-gl.diagnostics.hwe.R` (10 expectations)
  and `tests/testthat/test-gl.hwe.pop.R` (7), new files; all pass at the
  reviewed state, as do the existing `test-utils.basic.stats.R` and
  `test-utils.jackknife.R`
- Release state: both functions are on `main`/CRAN (last change on main:
  ed858a4). dartr2shiny wires both (config/functions.csv; its benchmark
  column records 2280 for `gl.diagnostics.hwe` against 0.991 for
  `gl.hwe.pop`). dartR.popgen calls `utils.basic.stats` (gl.ld.haplotype).

## Verdicts

**Standards: Needs work** — the structure is conventional, but the function
prints at `verbose = 0`, both functions return `-1` for a missing package,
and the jackknife's cost grows with the square of the number of loci.
**Spec: Rework** — the headline diagnostic, the ratio of the Fis and Fst
standard errors, is computed from leave-one-out values rounded to four
decimals, and the standard errors themselves are too small by a factor of
(number of loci - 1). The barplot compares tests with loci, and the Fisher
test uses the wrong degrees of freedom.

## Findings

**F1 [HIGH, confidence: high] — jackknife cost grows with loci squared (STY2)**
`R/gl.diagnostics.hwe.r:307-311`, `R/utils.jackknife.R:54-94` — with
`stdErr = TRUE` (the default), `utils.jackknife()` copies the object with
one locus dropped and reruns `utils.basic.stats()` once per locus. Every
per-locus statistic in `utils.basic.stats()` depends only on that locus's
genotypes, and the overall Fis and Fst are ratios of means over loci. Each
leave-one-out value is therefore a sum minus one term, which needs no
rerun.
Failure scenario: bandicoot.gl (1000 loci), `n.cores = 1`: 41.5 s against
1.2 s with `stdErr = FALSE`. Time grows with loci squared, so 10,000 loci
take about 70 min on one core. The default `n.cores = "auto"` spreads the
work over a cluster, so the wait is shorter but every core is busy.
Proposed change: compute the leave-one-out Fis and Fst directly from one
unrounded `utils.basic.stats()` table (O(loci)). This was checked against
brute force on 25 random loci of bandicoot.gl: max difference 3.3e-16.
`n.cores` becomes unused and is kept for compatibility.

**F2 [HIGH, confidence: high] — standard errors wrong in scale and ratio (spec: numerical correctness)**
`R/gl.diagnostics.hwe.r:310-313` — (a) the standard error is
`sqrt(var(theta_i) / n)`. The jackknife standard error is
`sqrt((n - 1) / n * sum((theta_i - mean)^2))`, which equals
`(n - 1) * sqrt(var(theta_i) / n)`. De Meeus 2007 uses the same
pseudovalue form. (b) The leave-one-out values are read from
`$overall`, which `utils.basic.stats()` rounds to 4 decimals. Dropping one
of 1000 loci moves Fis and Fst by less than that, so on bandicoot.gl the
1000 values collapse to 15 distinct Fis and 3 distinct Fst values.
Failure scenario: bandicoot.gl: reported SE Fis 4.43e-06 and SE Fst
7.89e-07, against 0.00432 and 0.000801 for the correct jackknife; the
reported ratio is 5.62 against 5.40. The ratio is what De Meeus (2018)
reads to separate null alleles from a Wahlund effect. Rounding makes it
depend on how many loci collapse to the same 4-decimal value, and more
loci make that worse.
Proposed change: use the standard jackknife standard error on the
unrounded leave-one-out values from change 1.

**F3 [MEDIUM, confidence: high] — barplot compares tests with loci (spec: DOC5 proposed rule)**
`R/gl.diagnostics.hwe.r:207-224` — the "0" bar of the observed series is
the number of non-significant *tests* (loci × populations). The null
series and every other observed bar count *loci*. The null series is also
one random `rbinom()` draw rather than the expectation, so the plot
changes between runs.
Failure scenario: bandicoot.gl: observed "0" bar 4677 (tests), against
863 loci never significant and a null expectation of 767 loci. The plot
suggests a large excess of clean loci that does not exist.
Proposed change: count loci significant in 0, 1, ..., k populations, and
draw the null as the expected count, loci × `dbinom(0:k, k, alpha_val)`.

**F4 [MEDIUM, confidence: high] — Fisher test and expected count use loci, not tests (spec: numerical correctness)**
`R/gl.diagnostics.hwe.r:271,281` — `pchisq(..., 2 * nLoc(x))` and
`nExpected = alpha_val * nLoc(x)` assume every locus was tested in every
population. `gl.report.hwe()` skips loci that are monomorphic or
under-sampled in a population, and Fisher's method has 2k degrees of
freedom for k combined tests.
Failure scenario: bandicoot.gl: 953-976 tests per population against 1000
loci, so the degrees of freedom are 2000 instead of 1906-1952, and
`nExpected` is 50 instead of 47.7-48.8. The effect is small here but
large for sparse populations, where many loci are skipped. The p-value is
then conservative, so real excesses are missed.
Proposed change: use the number of tests per population, k, for both.

**F5 [MEDIUM, confidence: high] — gl.hwe.pop crashes on an object without populations (spec)**
`R/gl.hwe.pop.r:49` — the default population is built with
`dim = nLoc(x)` instead of `nInd(x)`.
Failure scenario: `pop(x) <- NULL; gl.hwe.pop(x)` stops with "Vector
length does no match number of individuals", although the documentation
promises a single population "pop1".
Proposed change: use `nInd(x)`.

**F6 [LOW, confidence: high] — prints at verbose 0 (VRB1, VRB3)**
`R/gl.diagnostics.hwe.r:169,315-328` — the singleton-population warning,
the standard-error line and the summary table print at any verbosity.
Failure scenario: `gl.diagnostics.hwe(x, verbose = 0)` prints the table;
in dartr2shiny and scripts this output cannot be silenced.
Proposed change: warning at `verbose >= 2`; table and standard errors at
`verbose >= 1`, as in the `gl.report.*` functions.

**F7 [LOW, confidence: high] — package guards return -1; wrong package guarded (DEP1)**
`R/gl.diagnostics.hwe.r:146-153`, `R/gl.hwe.pop.r:29-37` — both print an
error and `return(-1)`. `gl.diagnostics.hwe` requires ggtern, which it
never uses (it calls `gl.report.hwe(plot.out = FALSE)`, and ggtern is
needed only for ternary plots), and does not guard HardyWeinberg, which it
needs.
Failure scenario: without ggtern, the diagnostics return -1 although
nothing needs ggtern, and a script that goes on to use the result fails
later with an unrelated error.
Proposed change: `stop(error(...))` per DEP1; the diagnostics guard
HardyWeinberg instead of ggtern.

**F8 [LOW, confidence: high] — histogram reference line off (spec)**
`R/gl.diagnostics.hwe.r:179-190` — the uniform expectation is drawn at
`nrow / (bins - 1)`, and the bins span the data range rather than [0, 1].
Failure scenario: bandicoot.gl: line at 253.9 against the uniform
expectation of 241.2 per bin.
Proposed change: bins fixed on [0, 1] (`boundary = 0`,
`binwidth = 1 / bins`) and the line at `nrow / bins`, labelled "Expected
count under HWE".

**F9 [LOW, confidence: high] — documentation disagrees with behaviour (DOC5 proposed rule, DOC6 proposed rule)**
`R/gl.hwe.pop.r` roxygen: `plot.out` "otherwise returns a dataframe" (it
always returns a list); `plot_colors` default given as `gl.colors(2)` and
`gl.colors("dis")` (the code uses `c("gray90", "deeppink")`); "a barplot"
(the plot is a raster). `R/gl.diagnostics.hwe.r` roxygen: non-ASCII
("Meeûs", curly quotes, en dash) and `n.cores`, which change 1 makes
unused.
Failure scenario: users are told the wrong defaults and return type.
Proposed change: correct the text, use ASCII, and describe `n.cores` as
kept for compatibility. Docs only.

## Proposed changes

1. Leave-one-out Fis and Fst computed directly from one unrounded
   per-locus table instead of re-running per locus (F1). `n.cores` is kept
   but no longer used.
   **Consequence: on bandicoot.gl the run goes from 41.5 s to about 1 s;
   the leave-one-out values are no longer rounded, so the Fis/Fst ratio
   changes slightly (5.62 to 5.40).**
2. Standard jackknife standard error (F2).
   **Consequence: the reported StdErr values change scale by a factor of
   (loci - 1), for example SE Fis 4.43e-06 to 0.00432 on bandicoot.gl.**
3. Barplot counts loci in every bar and shows the expected null counts
   (F3). Plot only; deterministic.
4. Fisher test and nExpected use the number of tests per population (F4).
   **Consequence: `hwe_summary$nExpected` and `pvalue` change wherever
   loci were not tested in every population.**
5. `gl.hwe.pop` default population uses `nInd` (F5).
6. Output gated by verbosity (F6).
7. Package guards stop with an error; the diagnostics guard HardyWeinberg,
   not ggtern (F7).
   **Consequence: with a guarded package missing, the functions error
   instead of returning -1, and the diagnostics no longer need ggtern.**
8. Histogram bins on [0, 1] and the line at the uniform expectation (F8).
   Plot only.
9. Documentation fixes (F9). Docs only.

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, PLT, STY, API — run on both files
- Spec: each documented output against an independent computation on
  bandicoot.gl — run (F2, F3, F4, F8)
- Exact leave-one-out against brute force (25 random loci) — run, max
  difference 3.3e-16
- Profiling: `stdErr = TRUE` vs `FALSE`, one `utils.basic.stats()` call
  (0.031 s), one `gl.report.hwe()` call (0.081 s) — run
- `gl.hwe.pop` significance calls against `gl.report.hwe()` p-values on
  bandicoot.gl — run, 100 % agreement
- `gl.hwe.pop` on an object without populations — run, crashes (F5)
- FBM path (DAT6): SKIPPED — both functions densify through
  `as.matrix()`/`seppop()`; no FBM-specific path to test
- Parallel path (`n.cores > 1`): not timed; change 1 removes it
- Google Group / GitHub issues: not searched (not available: no browser
  session); the custodian reports the slowness directly

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | consequence (ratio shifts, n.cores unused) approved explicitly; requested speed-up |
| 2 | approved | Luis | consequence (StdErr scale change) approved explicitly |
| 3 | approved | Luis | |
| 4 | approved | Luis | consequence (nExpected, pvalue change) approved explicitly |
| 5 | approved | Luis | |
| 6 | approved | Luis | |
| 7 | approved | Luis | consequence (error instead of -1) stated in the option |
| 8 | approved | Luis | |
| 9 | approved | Luis | |

## Outcome

- Change 1 (F1): leave-one-out values from one `utils.basic.stats(x,
  rounded = FALSE)` call. bandicoot.gl `stdErr = TRUE`: 41.5 s to 1.1 s;
  200 individuals x 10,000 loci: 2.0 s. To expose unrounded values,
  `utils.basic.stats()` gained `rounded = TRUE` (default output unchanged;
  its rounded branch is verbatim because single-population callers rely
  on the column name "round(Fis, 4)"). `test-utils.basic.stats.R`,
  `test-utils.jackknife.R` and `test-gl.report.fstat.R` (99 expectations)
  pass.
- Change 2 (F2): standard jackknife SE. The baseline StdErr expectation
  was flipped (approved diff); the new test recomputes by brute force
  (drop each locus, unrounded) and matches to 1e-10. bandicoot.gl: SE Fis
  0.004321, SE Fst 0.000801, ratio 5.40.
- Change 3 (F3): bars count loci; null = loci x dbinom. Test: observed
  and null sum to 991 loci; "0" bar = 863.
- Change 4 (F4): Fisher df = 2k and nExpected = alpha x k. The baseline
  nExpected expectation was flipped (approved diff): 50 becomes
  48.80/47.65/48.55/47.75/48.45; ChiSquare unchanged; pvalue still 1 on
  bandicoot.gl.
- Change 5 (F5): `nInd`; new test, an object without populations gives a
  1 x 1000 matrix for "pop1".
- Change 6 (F6): output gated by verbosity. Also silenced, within the same
  approval: three `gl.colors()` start/end messages (the default arguments
  now pass `verbose = 0`, and `gl.report.hwe()` receives `plot_colors`) and
  the `geom_smooth()` formula message. Test: 0 lines on stdout and 0
  messages at `verbose = 0`. The plot is still displayed, as before.
- Change 7 (F7): `stop()` guards; the diagnostics guard HardyWeinberg, not
  ggtern.
- Change 8 (F8): histogram bins on [0, 1]; the line at 4824 / 20 = 241.2
  (test).
- Change 9 (F9): docs; both R files ASCII only; `man/` regenerated
  (unrelated cross-link drift from PR #421 reverted).
- Tests: diagnostics 16 expectations, gl.hwe.pop 9, utils.basic.stats 8,
  all pass. The diagnostics run end to end at `verbose = 3`.
- PR: #425.

```json
{
  "function": "gl.diagnostics.hwe",
  "companion": "gl.hwe.pop",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "f9f1be8",
  "verdict_standards": "needs_work",
  "verdict_spec": "rework",
  "findings": [
    {"id": "F1", "severity": "HIGH", "confidence": "high", "rule": "STY2", "status": "approved", "change": 1},
    {"id": "F2", "severity": "HIGH", "confidence": "high", "rule": "spec", "status": "approved", "change": 2},
    {"id": "F3", "severity": "MEDIUM", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 3},
    {"id": "F4", "severity": "MEDIUM", "confidence": "high", "rule": "spec", "status": "approved", "change": 4},
    {"id": "F5", "severity": "MEDIUM", "confidence": "high", "rule": "spec", "status": "approved", "change": 5},
    {"id": "F6", "severity": "LOW", "confidence": "high", "rule": "VRB1", "status": "approved", "change": 6},
    {"id": "F7", "severity": "LOW", "confidence": "high", "rule": "DEP1", "status": "approved", "change": 7},
    {"id": "F8", "severity": "LOW", "confidence": "high", "rule": "spec", "status": "approved", "change": 8},
    {"id": "F9", "severity": "LOW", "confidence": "high", "rule": "DOC5", "status": "approved", "change": 9}
  ],
  "coverage_skipped": ["DAT6: no FBM path", "parallel path not timed", "forum/issues: no browser session"],
  "status": "pr-open",
  "pr": 425
}
```
