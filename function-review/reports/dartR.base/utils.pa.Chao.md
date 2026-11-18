# Review: utils.pa.Chao (dartR.base)

- Family mode: analysis (internal helper; the Chao engine of
  `gl.report.pa`)
- Date: 2026-09-24
- Reviewer: Claude (claude-opus-5-5), dartr-function-review v2.0.0
- Package commit: 2b29bea (`origin/dev`)
- Custodian: Bernd Gruber (as for `gl.report.pa`); rewritten by Arthur
  Georges in f758c66 (2026-09-01) under the `gl.report.pa` review
- Datasets: `platypus.gl` (three population pairs), hand-built matrices
- Baseline: `tests/testthat/test-utils.pa.Chao.R` (snapshot captured
  pre-review, 7 expectations pass)
- Prior review: `function-review/reports/dartR.base/gl.report.pa.md` F1
  rewrote this helper (pair-only data, category-matched f1/f2, empty-set
  handling). This review does not revisit those fixes; they check out.

## Verdicts

**Standards: Ready**: an internal helper with a clear header comment,
plain matrix arithmetic, no verbosity or dependency concerns, and
identical output for FBM-backed input.

**Spec: Needs work**: the small-sample factor `(n - 1) / n` in the
Chao1 formula uses the number of private alleles for `n` instead of the
number of allele copies sampled. The effect is under 2% on `platypus.gl`,
but pairs with few private alleles get estimates that are too small.

## Findings

**F1 [LOW, confidence: medium] — `n` in `(n - 1) / n` is the number of
private alleles, not the sample size (spec axis)**
`R/utils.pa.Chao.r:22, 36, 38`: `n <- length(alt.private) +
length(ref.private)` counts the private alleles observed (S_obs in the
Chao notation). In the Chao1 estimator (Chao et al. 2017, Eq. 2c, both
the `f2 > 0` and the bias-corrected `f2 = 0` forms), `n` is the sample
size, the number of sampled units. Here the sampled units are allele
copies, so `n` is the number of allele copies in the pooled pair sample.
The confidence is medium because the paper was not re-read in this
session. The notation is standard for Chao1, but the choice of unit for
allele data (allele copies, not individuals) is this review's reading.
Failure scenario (verified):
- two private singletons among 6 sampled allele copies: the code gives
  `(2 - 1) / 2 × 1 = 0.5`, reported by `gl.report.pa` as 0; with
  `n = 6` the estimate is 0.83, reported as 1;
- `platypus.gl` pairs (37–159 private alleles): the code and the
  corrected factor differ by less than 2% (e.g. SEVERN_BELOW vs
  TENTERFIELD 4.93 vs 5.02; SEVERN_ABOVE vs SEVERN_BELOW 62.60 vs
  62.49). After rounding in `gl.report.pa`, values change by at most 1.
Proposed change: set `n` to the number of non-missing allele copies in
the pooled pair sample, averaged over the private loci (`2 ×` the mean
number of genotyped individuals at those loci), and keep S_obs out of
the formula.

## Proposed changes

1. Use the number of sampled allele copies (pooled pair, mean over
   private loci) as `n` in `(n - 1) / n` (F1).
   **Consequence: `gl.report.pa` Chao1/Chao2 values change; noticeably
   only for pairs with few private alleles, by at most 1 after rounding
   on `platypus.gl`.**

## Coverage

- Standards walk: FS, DOC, VRB, DAT, DEP, STY — run. FS1–FS10 and DOC
  apply to exported functions; this is an internal helper with a header
  comment, which suffices.
- Spec: toy matrices with known f1/f2 and three `platypus.gl` pairs,
  against an independent re-implementation — run.
- FBM path (DAT6): `gl.gen2fbm(platypus.gl)` pairs give identical
  output — run.
- Prior fixes (pair-only data, category-matched f1/f2, n = 0 returns 0):
  confirmed by the baseline tests — run.
- Method source: Chao et al. (2017) not re-read (not available in this
  session); F1 rests on the standard Chao1 notation.
- Callers: `gl.report.pa` only (grep of `R/`).

## Approval

| Change | Decision | By | Note |
|---|---|---|---|
| 1 | approved | Luis | |

## Outcome

- Change 1 applied: `n` is the mean number of non-missing allele copies
  in the pooled pair at the private loci. Two private singletons among
  12 allele copies: 0.5 before, 11/12 after.
- `gl.report.pa(platypus.gl)` Chao1/Chao2 after rounding: 62/22, 18/54,
  5/46 (before: 63/22, 18/54, 5/46).
- Snapshot: 7 expectations pass; each diff from the baseline is marked
  `[approved diff change 1]`. `test-gl.report.pa.R`: 13 pass. FBM
  output still identical.
- NEWS entry added (escalation gate: numerical output changes).
- PR: pending.

```json
{
  "function": "utils.pa.Chao",
  "package": "dartR.base",
  "family": "analysis",
  "skill_version": "2.0.0",
  "commit": "2b29bea",
  "verdict_standards": "ready",
  "verdict_spec": "needs_work",
  "findings": [
    {"id": "F1", "severity": "LOW", "confidence": "medium", "rule": "spec", "status": "approved", "change": 1}
  ],
  "coverage_skipped": ["method paper not re-read"],
  "status": "pr-open",
  "pr": null
}
```
